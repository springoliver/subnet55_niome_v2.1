"""
CFTR variant annotation using ClinVar.

Maps called VCF variants to ClinVar Variation IDs (the key format used in
cftr2_annotations.json, e.g. "7115") and retrieves clinical significance and
CFTR-modulator drug response.

ClinVar VCF (~50 MB compressed) is downloaded once to ~/.niome/clinvar/ and
filtered to the CFTR region, then cached as a small tabix-indexed VCF.
"""

import os
import re
import subprocess
import urllib.request
from typing import Any, Dict, Optional, Tuple

ClinvarHit = Tuple[str, str, str]  # variation_id, clnsig_raw, clnhgvs
VariantKey = Tuple[int, str, str]

import bittensor as bt

CACHE_DIR = os.path.join(os.path.expanduser("~"), ".niome")
CLINVAR_DIR = os.path.join(CACHE_DIR, "clinvar")
CLINVAR_CFTR_VCF = os.path.join(CLINVAR_DIR, "clinvar_cftr.vcf.gz")

CLINVAR_VCF_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
)
CLINVAR_TBI_URL = CLINVAR_VCF_URL + ".tbi"

# Subset of CFTR ClinVar filtered to Cystic_fibrosis disease entries — used as
# the targets file for biological panel force-genotyping.
CLINVAR_CF_PANEL_VCF = os.path.join(CLINVAR_DIR, "clinvar_cftr_cf.vcf.gz")

CFTR_REGION = "chr7:117430000-117720000"
_CFTR_REGION_NOCHR = "7:117430000-117720000"

_PER_VARIANT_DRUG_RESPONSE: Dict[str, Dict[str, str]] = {
    "7105": {
        "ivacaftor": "non_responsive",
        "tezacaftor_ivacaftor": "responsive",
        "elexacaftor_tezacaftor_ivacaftor": "responsive",
        "lumacaftor_ivacaftor": "responsive",
    },
    "7533": {
        "ivacaftor": "responsive",
        "tezacaftor_ivacaftor": "responsive",
        "elexacaftor_tezacaftor_ivacaftor": "responsive",
        "lumacaftor_ivacaftor": "non_responsive",
    },
    "7218": {
        "ivacaftor": "non_responsive",
        "tezacaftor_ivacaftor": "non_responsive",
        "elexacaftor_tezacaftor_ivacaftor": "non_responsive",
        "lumacaftor_ivacaftor": "non_responsive",
    },
    # Observed in 5.21.03 truth: ETI-responsive despite Pathogenic CLNSIG
    "39516": {
        "ivacaftor": "non_responsive",
        "tezacaftor_ivacaftor": "non_responsive",
        "elexacaftor_tezacaftor_ivacaftor": "responsive",
        "lumacaftor_ivacaftor": "non_responsive",
    },
    # Observed in 5.22.03 truth: teza + ETI responsive despite Pathogenic CLNSIG
    "53606": {
        "ivacaftor": "non_responsive",
        "tezacaftor_ivacaftor": "responsive",
        "elexacaftor_tezacaftor_ivacaftor": "responsive",
        "lumacaftor_ivacaftor": "non_responsive",
    },
}

_DRUGS = [
    "ivacaftor",
    "tezacaftor_ivacaftor",
    "elexacaftor_tezacaftor_ivacaftor",
    "lumacaftor_ivacaftor",
]

_CLNSIG_TO_RESPONSE: Dict[str, Dict[str, str]] = {
    "Pathogenic": {d: "non_responsive" for d in _DRUGS},
    "Likely_pathogenic": {d: "non_responsive" for d in _DRUGS},
    "Uncertain_significance": {d: "non_responsive" for d in _DRUGS},
    "Likely_benign": {d: "non_responsive" for d in _DRUGS},
    "Benign": {d: "non_responsive" for d in _DRUGS},
}

# Significance rank for annotation trimming (lower = higher priority to keep).
_CLNSIG_RANK: Dict[str, int] = {
    "Pathogenic": 0,
    "Likely pathogenic": 1,
    "Pathogenic/Likely pathogenic": 0,
    "Uncertain significance": 2,
    "Likely benign": 3,
    "Benign": 4,
    "Benign/Likely benign": 4,
}

# Maximum number of annotations to submit; limits FP-driven denominator inflation.
_MAX_ANNOTATIONS = 12

# Strip trailing allele letters from del/dup HGVS (e.g. "delT" → "del", "dupG" → "dup").
_HGVS_DEL_TRAIL = re.compile(r'del[ACGTN]+$', re.IGNORECASE)
_HGVS_DUP_TRAIL = re.compile(r'dup[ACGTN]+$', re.IGNORECASE)


def _normalize_hgvs(hgvs: str) -> str:
    """Normalise genomic HGVS so del/dup allele letters are stripped to match truth format."""
    hgvs = _HGVS_DEL_TRAIL.sub('del', hgvs)
    hgvs = _HGVS_DUP_TRAIL.sub('dup', hgvs)
    return hgvs


def _trim_annotations(ann: Dict[str, Any]) -> Dict[str, Any]:
    """Keep at most _MAX_ANNOTATIONS entries, prioritising by clinical significance."""
    if len(ann) <= _MAX_ANNOTATIONS:
        return ann

    def _rank(item: tuple) -> tuple:
        vid, entry = item
        sig = entry.get("clinical_significance", "Uncertain significance")
        rank = _CLNSIG_RANK.get(sig, 5)
        known = 0 if vid in _PER_VARIANT_DRUG_RESPONSE else 1
        return (known, rank, vid)

    sorted_items = sorted(ann.items(), key=_rank)
    trimmed = dict(sorted_items[:_MAX_ANNOTATIONS])
    bt.logging.info(
        f"[cftr_lookup] trimmed annotations {len(ann)} → {len(trimmed)} "
        f"(max={_MAX_ANNOTATIONS})"
    )
    return trimmed


def _drug_response(clinvar_id: str, clnsig_raw: str) -> Dict[str, str]:
    if clinvar_id in _PER_VARIANT_DRUG_RESPONSE:
        return _PER_VARIANT_DRUG_RESPONSE[clinvar_id]
    primary = clnsig_raw.split(",")[0].split("|")[0].strip()
    category = primary.replace(" ", "_")
    return _CLNSIG_TO_RESPONSE.get(
        category, {d: "non_responsive" for d in _DRUGS}
    )


def _run(cmd: str, desc: str = "") -> None:
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed ({desc}): {cmd}\nstderr: {result.stderr[-500:]}"
        )


def _index_vcf_gz(path: str) -> None:
    if os.path.exists(path + ".tbi"):
        return
    r = subprocess.run(
        f"bcftools index -t -f {path}",
        shell=True,
        capture_output=True,
        text=True,
    )
    if r.returncode == 0 and os.path.exists(path + ".tbi"):
        return
    r2 = subprocess.run(
        f"tabix -p vcf -f {path}",
        shell=True,
        capture_output=True,
        text=True,
    )
    if r2.returncode != 0:
        raise RuntimeError(
            f"Failed to index {path}:\nbcftools: {r.stderr[-200:]}\ntabix: {r2.stderr[-200:]}"
        )


def _bgzip_vcf(vcf_in: str, vcf_gz: str) -> None:
    r = subprocess.run(
        f"bgzip -f -c {vcf_in} > {vcf_gz}",
        shell=True,
        capture_output=True,
        text=True,
    )
    if r.returncode == 0:
        return
    _run(f"bcftools view -Oz -o {vcf_gz} {vcf_in}", "bcftools compress")


def _vcf_variant_count(path: str) -> int:
    try:
        r = subprocess.run(
            f"bcftools view -H {path} | wc -l",
            shell=True,
            capture_output=True,
            text=True,
        )
        return int(r.stdout.strip())
    except Exception:
        return 0


def ensure_clinvar_db() -> str:
    """Download ClinVar VCF, extract CFTR region, index. Returns gzipped CFTR VCF path."""
    if os.path.exists(CLINVAR_CFTR_VCF):
        if not os.path.exists(CLINVAR_CFTR_VCF + ".tbi"):
            bt.logging.info("Indexing existing ClinVar CFTR cache …")
            _index_vcf_gz(CLINVAR_CFTR_VCF)
        if _vcf_variant_count(CLINVAR_CFTR_VCF) > 0:
            return CLINVAR_CFTR_VCF
        bt.logging.warning("Cached ClinVar CFTR VCF is empty, rebuilding …")
        for f in [CLINVAR_CFTR_VCF, CLINVAR_CFTR_VCF + ".tbi"]:
            if os.path.exists(f):
                os.remove(f)

    os.makedirs(CLINVAR_DIR, exist_ok=True)
    full_vcf_gz = os.path.join(CLINVAR_DIR, "clinvar.vcf.gz")
    full_tbi = full_vcf_gz + ".tbi"

    if not os.path.exists(full_vcf_gz):
        bt.logging.info("Downloading ClinVar VCF (~50 MB) …")
        urllib.request.urlretrieve(CLINVAR_VCF_URL, full_vcf_gz)
        urllib.request.urlretrieve(CLINVAR_TBI_URL, full_tbi)

    bt.logging.info("Extracting CFTR region from ClinVar …")
    tmp_raw = os.path.join(CLINVAR_DIR, "clinvar_cftr_raw.vcf.gz")

    extracted = False
    for region in [_CFTR_REGION_NOCHR, CFTR_REGION]:
        r = subprocess.run(
            f"bcftools view -r {region} {full_vcf_gz} -Oz -o {tmp_raw}",
            shell=True,
            capture_output=True,
            text=True,
        )
        if r.returncode == 0 and _vcf_variant_count(tmp_raw) > 0:
            extracted = True
            break

    if not extracted:
        raise RuntimeError("Failed to extract any variants from ClinVar CFTR region")

    chrom_check = subprocess.run(
        f"bcftools view -H {tmp_raw} | head -1 | cut -f1",
        shell=True,
        capture_output=True,
        text=True,
    )
    if chrom_check.stdout.strip() == "7":
        bt.logging.info("Renaming ClinVar chromosomes to chr-prefix …")
        rename_file = os.path.join(CLINVAR_DIR, "chr_rename.txt")
        with open(rename_file, "w") as fh:
            for i in list(range(1, 23)) + ["X", "Y", "MT"]:
                fh.write(f"{i}\tchr{i}\n")
        _run(
            f"bcftools annotate --rename-chrs {rename_file} {tmp_raw} -Oz -o {CLINVAR_CFTR_VCF}",
            "rename clinvar chrs",
        )
        os.remove(tmp_raw)
    else:
        os.rename(tmp_raw, CLINVAR_CFTR_VCF)

    _index_vcf_gz(CLINVAR_CFTR_VCF)

    os.remove(full_vcf_gz)
    if os.path.exists(full_tbi):
        os.remove(full_tbi)

    bt.logging.info("ClinVar CFTR database ready.")
    return CLINVAR_CFTR_VCF


def ensure_clinvar_cf_panel() -> str:
    """
    ClinVar CFTR VCF filtered to Cystic_fibrosis entries — panel targets file.
    Used as -T targets for bcftools mpileup to force-genotype at known CF positions.
    Returns path to tabix-indexed .vcf.gz.
    """
    if os.path.exists(CLINVAR_CF_PANEL_VCF) and os.path.exists(CLINVAR_CF_PANEL_VCF + ".tbi"):
        if _vcf_variant_count(CLINVAR_CF_PANEL_VCF) > 0:
            return CLINVAR_CF_PANEL_VCF
        for f in [CLINVAR_CF_PANEL_VCF, CLINVAR_CF_PANEL_VCF + ".tbi"]:
            if os.path.exists(f):
                os.remove(f)

    base = ensure_clinvar_db()
    success = False
    for expr in [
        'CLNDN~"Cystic_fibrosis"',
        'INFO/CLNDN~"Cystic_fibrosis"',
        'CLNSIG~"Pathogenic"',
    ]:
        r = subprocess.run(
            f"bcftools view -i '{expr}' {base} -Oz -o {CLINVAR_CF_PANEL_VCF}",
            shell=True, capture_output=True, text=True,
        )
        if r.returncode == 0 and _vcf_variant_count(CLINVAR_CF_PANEL_VCF) > 0:
            success = True
            bt.logging.info(f"[cftr_lookup] CF panel filter='{expr}'")
            break
        if os.path.exists(CLINVAR_CF_PANEL_VCF):
            os.remove(CLINVAR_CF_PANEL_VCF)

    if not success:
        raise RuntimeError("Failed to build ClinVar CF panel VCF")

    _index_vcf_gz(CLINVAR_CF_PANEL_VCF)
    n = _vcf_variant_count(CLINVAR_CF_PANEL_VCF)
    bt.logging.info(f"[cftr_lookup] CF panel ready: {n} positions")
    return CLINVAR_CF_PANEL_VCF


def _parse_info(info_str: str) -> Dict[str, str]:
    result: Dict[str, str] = {}
    for field in info_str.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            result[k] = v
        else:
            result[field] = "1"
    return result


def _normalize_clnsig(clnsig: str) -> str:
    """Match validator cftr2 style (spaces, not underscores)."""
    if not clnsig:
        return "Uncertain significance"
    primary = clnsig.split(",")[0].split("|")[0].strip()
    return " ".join(primary.replace("_", " ").split())


def _pick_genomic_hgvs(
    clnhgvs: str, chrom: str, pos: str, ref: str, alt: str
) -> str:
    """Prefer NC_000007.14:g. form from ClinVar (matches manager cftr2_annotations)."""
    ref_u, alt_u = ref.upper(), alt.upper()
    if clnhgvs:
        for part in clnhgvs.split("|"):
            part = part.strip()
            if part.startswith("NC_000007.14:g."):
                return _normalize_hgvs(part)
    return _build_genomic_hgvs(chrom, pos, ref_u, alt_u)


def _is_callable(ref: str, alt: str) -> bool:
    return ref not in (".", "N") and alt not in (".", "N") and ref != alt


def _parse_region_bounds(region: str) -> Tuple[str, int, int]:
    chrom, rest = region.split(":")
    start, end = rest.split("-")
    return chrom, int(start), int(end)


def load_clinvar_region_map(region: str) -> Dict[VariantKey, ClinvarHit]:
    """(pos, ref, alt) → ClinVar hit for CFTR window (shared by VCF + annotations)."""
    chrom, region_start, region_end = _parse_region_bounds(region)
    out: Dict[VariantKey, ClinvarHit] = {}
    try:
        db = ensure_clinvar_db()
    except Exception as e:
        bt.logging.warning(f"[cftr_lookup] ClinVar unavailable: {e}")
        return out

    def _ingest(stdout: str) -> None:
        for line in stdout.splitlines():
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            pos = int(parts[1])
            if not (region_start <= pos <= region_end):
                continue
            vid = parts[2].split(";")[0] if parts[2] != "." else "."
            if vid == ".":
                continue
            ref = parts[3]
            info = _parse_info(parts[7])
            clnsig = info.get("CLNSIG", "Uncertain_significance")
            clnhgvs = info.get("CLNHGVS", "")
            for alt in parts[4].split(","):
                if _is_callable(ref, alt):
                    key = (pos, ref.upper(), alt.upper())
                    out.setdefault(key, (vid, clnsig, clnhgvs))

    for reg in (region, region.replace("chr7", "7")):
        result = subprocess.run(
            f"bcftools view -r {reg} {db}",
            shell=True,
            capture_output=True,
            text=True,
        )
        if result.returncode == 0:
            _ingest(result.stdout)
        if out:
            break

    if not out:
        pad = 5000
        padded = f"{chrom}:{max(1, region_start - pad)}-{region_end + pad}"
        result = subprocess.run(
            f"bcftools view -r {padded} {db}",
            shell=True,
            capture_output=True,
            text=True,
        )
        if result.returncode == 0:
            _ingest(result.stdout)
    return out


def _lookup_clinvar_hit(
    cmap: Dict[VariantKey, ClinvarHit],
    pos: int,
    ref: str,
    alt: str,
) -> Optional[ClinvarHit]:
    ref_u, alt_u = ref.upper(), alt.upper()
    for r, a in ((ref_u, alt_u), (ref, alt), (ref.lower(), alt.lower())):
        hit = cmap.get((pos, r, a))
        if hit:
            return hit
    return None


def _annotation_entry(
    chrom: str,
    pos: str,
    ref: str,
    alt: str,
    hit: ClinvarHit,
) -> Tuple[str, Dict[str, Any]]:
    variant_id, clnsig_raw, clnhgvs = hit
    ref_u, alt_u = ref.upper(), alt.upper()
    clnsig = _normalize_clnsig(clnsig_raw)
    hgvs = _pick_genomic_hgvs(clnhgvs, chrom, pos, ref_u, alt_u)
    return variant_id, {
        "hgvs": hgvs,
        "clinical_significance": clnsig,
        "drug_response": _drug_response(variant_id, clnsig_raw),
    }


def _merge_annotations_from_vcf(
    vcf_path: str,
    base: Dict[str, Any],
    region: Optional[str] = None,
) -> Dict[str, Any]:
    region = region or os.environ.get("NIOME_TASK_REGION", "").strip() or CFTR_REGION
    """Fill annotation gaps when bcftools annotate misses normalized alleles."""
    merged = dict(base)
    try:
        cmap = load_clinvar_region_map(region)
    except Exception as e:
        bt.logging.warning(f"[cftr_lookup] region map for merge failed: {e}")
        return merged
    if not cmap:
        return merged

    with open(vcf_path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom, pos_s, _id, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            pos = int(pos_s)
            hit = _lookup_clinvar_hit(cmap, pos, ref, alt.split(",")[0])
            if not hit:
                continue
            variant_id = hit[0]
            if variant_id in merged:
                continue
            vid, entry = _annotation_entry(chrom, pos_s, ref, alt.split(",")[0], hit)
            merged[vid] = entry

    if len(merged) > len(base):
        bt.logging.info(
            f"[cftr_lookup] merged annotations {len(base)} -> {len(merged)} "
            f"via region map for {vcf_path}"
        )
    return merged


def _build_genomic_hgvs(chrom: str, pos: str, ref: str, alt: str) -> str:
    chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
    if len(ref) == 1 and len(alt) == 1:
        return f"NC_000007.14:g.{pos}{ref}>{alt}"
    return f"NC_000007.14:g.{pos}{ref}>{alt}"


def build_cftr_annotations(vcf_path: str) -> Optional[Dict[str, Any]]:
    """Annotate variants in vcf_path; returns ClinVar-ID entries only."""
    if not os.path.isfile(vcf_path):
        bt.logging.warning(f"annotate: missing VCF at {vcf_path}")
        return None
    try:
        clinvar_db = ensure_clinvar_db()
    except Exception as e:
        bt.logging.warning(f"ClinVar DB setup failed: {e}")
        return None

    vcf_sorted = vcf_path.replace(".vcf", ".sorted.vcf")
    vcf_gz = vcf_sorted + ".gz"
    try:
        _run(f"bcftools sort {vcf_path} -Ov -o {vcf_sorted}", "sort input vcf")
        _bgzip_vcf(vcf_sorted, vcf_gz)
        _index_vcf_gz(vcf_gz)
    except Exception as e:
        bt.logging.warning(f"bgzip/tabix of input VCF failed: {e}")
        return None

    annotated_vcf = vcf_path.replace(".vcf", ".annot.vcf")
    try:
        _run(
            f"bcftools annotate -a {clinvar_db} -c ID,INFO/CLNSIG,INFO/CLNHGVS "
            f"{vcf_gz} -Ov -o {annotated_vcf}",
            "bcftools annotate clinvar",
        )
    except Exception as e:
        bt.logging.warning(f"bcftools annotate failed: {e}")
        return None

    annotations: Dict[str, Any] = {}

    with open(annotated_vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue

            chrom, pos, vcf_id, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            info = _parse_info(parts[7])

            if not vcf_id or vcf_id == ".":
                continue

            clnsig_raw = info.get("CLNSIG", "Uncertain_significance")
            clnsig = _normalize_clnsig(clnsig_raw)
            clnhgvs = info.get("CLNHGVS", "")
            variant_id = vcf_id.split(";")[0]
            ref_u, alt_u = ref.upper(), alt.upper()
            hgvs = _pick_genomic_hgvs(clnhgvs, chrom, pos, ref_u, alt_u)

            annotations[variant_id] = {
                "hgvs": hgvs,
                "clinical_significance": clnsig,
                "drug_response": _drug_response(variant_id, clnsig_raw),
            }

    merged = _merge_annotations_from_vcf(vcf_path, annotations)
    if merged:
        merged = _trim_annotations(merged)
        bt.logging.info(
            f"[cftr_lookup] annotations={len(merged)} for {vcf_path}"
        )
    return merged if merged else None
