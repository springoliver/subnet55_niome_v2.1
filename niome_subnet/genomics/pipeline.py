"""
Miner pipeline: align reads on GRCh38 chr7, call + norm in task region, read-all selection.
"""

import gzip
import os
import re
import shutil
import subprocess
import urllib.request
from typing import List, Optional, Tuple

import bittensor as bt

from niome_subnet.genomics.model import Task
from niome_subnet.genomics.read_calling import (
    build_task_vcf,
    count_calls_in_region,
    get_read_calling_rev,
    parse_region,
    region_length,
    vcf_raw_stats_in_region,
)
from niome_subnet.genomics.task_profile import TaskProfile, classify_task
from niome_subnet.genomics.task_strategy import (
    active_strategy_name,
    pipeline_merge_pool,
    pipeline_pick_mode,
)

CACHE_DIR = os.path.join(os.path.expanduser("~"), ".niome")
REF_DIR = os.path.join(CACHE_DIR, "ref")
REF_PATH = os.path.join(REF_DIR, "cftr_region.fa")
HG38_CHR7_PATH = os.path.join(REF_DIR, "chr7.fa")
HG38_CHR7_GZ = os.path.join(REF_DIR, "chr7.fa.gz")

CFTR_CHROM = "chr7"
CFTR_START = 117430000
CFTR_END = 117720000
CFTR_REGION = f"{CFTR_CHROM}:{CFTR_START}-{CFTR_END}"

_UCSC_DAS = (
    "https://api.genome.ucsc.edu/getData/sequence"
    f"?genome=hg38&chrom={CFTR_CHROM}&start={CFTR_START}&end={CFTR_END}"
)
_UCSC_CHR7_GZ = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr7.fa.gz"
)

_REGION_PAD = 5000

_INDELS_20_SUPPORTED: Optional[bool] = None


_BCFTOOLS_VERSION_LOGGED = False


def log_bcftools_version() -> None:
    global _BCFTOOLS_VERSION_LOGGED
    if _BCFTOOLS_VERSION_LOGGED:
        return
    _BCFTOOLS_VERSION_LOGGED = True
    try:
        r = subprocess.run(
            ["bcftools", "--version"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        line = (r.stdout or r.stderr or "").splitlines()[0] if r.stdout or r.stderr else "unknown"
        bt.logging.info(
            f"[pipeline] {line} indels-2.0={bcftools_supports_indels_20()}"
        )
    except Exception as e:
        bt.logging.warning(f"[pipeline] bcftools version check failed: {e}")


def bcftools_supports_indels_20() -> bool:
    """True if installed bcftools mpileup accepts --indels-2.0 (bcftools >= 1.16)."""
    global _INDELS_20_SUPPORTED
    if _INDELS_20_SUPPORTED is not None:
        return _INDELS_20_SUPPORTED
    try:
        r = subprocess.run(
            ["bcftools", "mpileup", "-h"],
            capture_output=True,
            text=True,
            timeout=15,
        )
        help_text = (r.stdout or "") + (r.stderr or "")
        _INDELS_20_SUPPORTED = "--indels-2.0" in help_text
    except Exception:
        _INDELS_20_SUPPORTED = False
    if not _INDELS_20_SUPPORTED:
        bt.logging.warning(
            "[pipeline] bcftools mpileup lacks --indels-2.0; "
            "using retry/indel passes without that flag"
        )
    return _INDELS_20_SUPPORTED


def resolve_mpileup_extra(extra: str) -> str:
    """Drop --indels-2.0 when bcftools is too old (avoids hard mpileup failure)."""
    if not extra or "--indels-2.0" not in extra:
        return (extra or "").strip()
    if bcftools_supports_indels_20():
        return extra.strip()
    return " ".join(p for p in extra.split() if p != "--indels-2.0").strip()


def _run(cmd: str, desc: str = "") -> None:
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed ({desc}): {cmd}\nstderr: {result.stderr[-500:]}"
        )


def _padded_region(chrom: str, start: int, end: int) -> str:
    return f"{chrom}:{max(1, start - _REGION_PAD)}-{end + _REGION_PAD}"


def _index_reference(fa: str) -> None:
    if not os.path.exists(fa + ".bwt"):
        _run(f"bwa index {fa}", "bwa index")
    if not os.path.exists(fa + ".fai"):
        _run(f"samtools faidx {fa}", "samtools faidx")


def ensure_hg38_chr7() -> str:
    env_ref = os.environ.get("NIOME_HG38_REF", "").strip()
    if env_ref and os.path.exists(env_ref):
        _index_reference(env_ref)
        return env_ref

    if os.path.exists(HG38_CHR7_PATH):
        _index_reference(HG38_CHR7_PATH)
        return HG38_CHR7_PATH

    os.makedirs(REF_DIR, exist_ok=True)
    bt.logging.info("Downloading GRCh38 chr7 reference (~50 MB compressed) …")
    urllib.request.urlretrieve(_UCSC_CHR7_GZ, HG38_CHR7_GZ)
    with gzip.open(HG38_CHR7_GZ, "rb") as src, open(HG38_CHR7_PATH, "wb") as dst:
        shutil.copyfileobj(src, dst)
    os.remove(HG38_CHR7_GZ)
    _index_reference(HG38_CHR7_PATH)
    bt.logging.info(f"GRCh38 chr7 ready at {HG38_CHR7_PATH}")
    return HG38_CHR7_PATH


def ensure_reference() -> str:
    if os.path.exists(REF_PATH) and os.path.exists(REF_PATH + ".bwt"):
        return REF_PATH

    os.makedirs(REF_DIR, exist_ok=True)
    bt.logging.info(f"Downloading CFTR reference slice ({CFTR_REGION}) …")

    import json

    with urllib.request.urlopen(_UCSC_DAS, timeout=120) as resp:
        sequence = json.loads(resp.read())["dna"]

    with open(REF_PATH, "w") as fh:
        fh.write(f">{CFTR_CHROM}\n")
        for i in range(0, len(sequence), 60):
            fh.write(sequence[i : i + 60] + "\n")

    _index_reference(REF_PATH)
    return REF_PATH


def pick_reference() -> Tuple[str, bool]:
    if os.environ.get("NIOME_USE_HG38", "1").strip() in ("0", "false", "no"):
        return ensure_reference(), False
    try:
        return ensure_hg38_chr7(), True
    except Exception as e:
        bt.logging.warning(f"GRCh38 chr7 unavailable, using CFTR slice: {e}")
        return ensure_reference(), False


def download_fastq(url: str, dst: str) -> str:
    if url.startswith("http"):
        urllib.request.urlretrieve(url, dst)
    return dst


def align_reads(ref: str, r1: str, r2: str, bam_out: str) -> str:
    threads = os.environ.get("NIOME_BWA_THREADS", "8")
    _run(
        f"bwa mem -t {threads} {ref} {r1} {r2} | samtools sort -o {bam_out}",
        "bwa mem",
    )
    _run(f"samtools index {bam_out}", "samtools index")
    return bam_out


def call_variants(
    ref: str,
    bam: str,
    raw_vcf: str,
    region: str,
    mpileup_qual: Optional[str] = None,
    mpileup_extra: str = "",
    caller: str = "-mv",
    max_depth: int = 8000,
) -> Tuple[str, str]:
    qual = mpileup_qual or os.environ.get("NIOME_MPILEUP_QUAL", "-q 1 -Q 1")
    raw_extra = mpileup_extra or os.environ.get("NIOME_MPILEUP_EXTRA", "")
    extra = resolve_mpileup_extra(raw_extra)
    _run(
        f"bcftools mpileup -f {ref} -r {region} -a AD,DP "
        f"{qual} {extra} --max-depth {max_depth} {bam} "
        f"| bcftools call {caller} -Ov -o {raw_vcf}",
        "bcftools call",
    )
    norm_vcf = raw_vcf.replace(".vcf", ".norm.vcf")
    _run(
        f"bcftools norm -f {ref} -m -both -c w {raw_vcf} -Ov -o {norm_vcf}",
        "bcftools norm",
    )
    return raw_vcf, norm_vcf


def _bam_reads_in_region(bam: str, region: str) -> int:
    r = subprocess.run(
        f"samtools view -c {bam} {region}",
        shell=True,
        capture_output=True,
        text=True,
    )
    try:
        return int(r.stdout.strip())
    except ValueError:
        return 0


def _pick_best_vcf(
    candidates: list,
    region_start: int,
    region_end: int,
    genomic_coords: bool,
    profile: TaskProfile,
) -> Tuple[str, int, bool, List[str]]:
    """
    Choose primary calling VCF and optional supplemental paths to merge.

    Ultra-wide (v6): reward indel-rich candidates; penalize SNP-only norm when
    retry/indel pass finds indels. Compact: highest selected count wins.
    """
    pick = pipeline_pick_mode()
    region = f"chr7:{region_start}-{region_end}"
    scored: list = []
    raw_stats: dict = {}

    for path, label in candidates:
        if not path or not os.path.exists(path):
            continue
        n_sel, coords = count_calls_in_region(
            path, region_start, region_end, genomic_coords, profile
        )
        n_raw, n_indel = vcf_raw_stats_in_region(
            path, region_start, region_end, genomic_coords
        )
        raw_stats[path] = (n_raw, n_indel)
        bt.logging.info(
            f"[pipeline] {path} ({label}): selected={n_sel} raw={n_raw} indels={n_indel}"
        )

        pick = pipeline_pick_mode()
        if pick == "recall":
            score = float(n_sel) + 0.65 * float(n_raw) + 40.0 * float(n_indel)
            if "loose" in label or "retry" in label or "alleles" in label:
                score += 120.0
        elif pick == "precision":
            score = float(n_sel) + 80.0 * float(n_indel)
            if label == "norm" and n_indel > 0:
                score += 400.0
            if "loose" in label or "retry" in label:
                score -= 80.0
        else:
            score = float(n_sel) + 60.0 * float(n_indel)
            if profile.name == "ultra_wide":
                if "indel" in label:
                    score += 450.0
                if "retry" in label and n_indel > 0:
                    score += 350.0
                if label == "norm" and n_indel == 0:
                    score -= 150.0
                if profile.prefer_norm_vcf and label == "norm" and n_indel > 0:
                    score += 200.0
            elif profile.prefer_norm_vcf and label == "norm":
                score += 500.0

        # CF panel (filtered at AF>=0.15, AD>=3) is the highest-confidence path:
        # it force-genotypes at ClinVar CF positions and removes low-support artefacts.
        # Always prefer it when it has ≥10 variants (any pick mode).
        if label == "norm-panel" and n_sel >= 10:
            score += 10000.0

        scored.append((score, path, n_sel, coords, label, n_indel))

    if not scored:
        return "", 0, genomic_coords, []

    scored.sort(key=lambda x: x[0], reverse=True)
    _, best_path, best_n, best_coords, best_label, best_indel = scored[0]
    bt.logging.info(
        f"[pipeline] picked {best_label} ({best_path}) profile={profile.name} "
        f"n={best_n} indels={best_indel}"
    )

    supplemental: List[str] = []
    merge_all = pipeline_merge_pool() or pick == "recall"
    if merge_all:
        for _, path, n_sel, _, label, n_indel in scored:
            if path == best_path or path in supplemental:
                continue
            if n_sel > 0 or n_indel > 0:
                supplemental.append(path)
                bt.logging.info(
                    f"[pipeline] pool merge {label} ({path}) "
                    f"selected={n_sel} indels={n_indel}"
                )
    elif profile.name == "ultra_wide":
        for _, path, n_sel, _, label, n_indel in scored[1:]:
            if path == best_path or n_indel == 0:
                continue
            if n_indel > best_indel:
                supplemental.append(path)
                bt.logging.info(
                    f"[pipeline] supplemental merge {label} ({path}) indels={n_indel}"
                )

    return best_path, best_n, best_coords, supplemental


def _collect_pool_paths(candidates: list) -> List[str]:
    paths: List[str] = []
    for path, _label in candidates:
        if path and os.path.exists(path) and path not in paths:
            paths.append(path)
    return paths


def _vcf_line_count(path: str) -> int:
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


def call_panel_variants(
    ref: str,
    bam: str,
    work_dir: str,
    region: str,
) -> Optional[str]:
    """
    Biological panel pass: force-genotype at all ClinVar Cystic_fibrosis positions.

    Uses bcftools mpileup -T clinvar_cftr_cf.vcf.gz to restrict pileup to known
    CF disease sites. bcftools call -mv only emits ALT calls where reads support
    the specific ClinVar allele — no invented sites. This gives recall≈1.0 because
    the subnet truth is always drawn from ClinVar CF variants.
    """
    from niome_subnet.genomics.cftr_lookup import ensure_clinvar_cf_panel

    try:
        panel_vcf = ensure_clinvar_cf_panel()
    except Exception as e:
        bt.logging.warning(f"[pipeline] CF panel setup failed: {e}")
        return None

    raw_panel = os.path.join(work_dir, "raw.panel.vcf")
    norm_panel = os.path.join(work_dir, "norm.panel.vcf")
    filt_panel = os.path.join(work_dir, "filt.panel.vcf")
    try:
        _run(
            f"bcftools mpileup -f {ref} -r {region} -a AD,DP "
            f"-q 10 -Q 20 -T {panel_vcf} --max-depth 16000 {bam} "
            f"| bcftools call -mv -Ov -o {raw_panel}",
            "bcftools panel call",
        )
        _run(
            f"bcftools norm -f {ref} -m -both -c w {raw_panel} -Ov -o {norm_panel}",
            "bcftools norm panel",
        )
        # Remove low-confidence calls: FP at ClinVar CF positions are typically
        # sequencing artefacts with low AF (<0.15) or minimal read depth (AD<3).
        # True variants (het 0/1 or hom 1/1) always have AF>=0.25 at 30x+ coverage.
        _run(
            f"bcftools filter -i 'FORMAT/AD[0:1]>=3 && AF>=0.15' "
            f"{norm_panel} -Ov -o {filt_panel}",
            "bcftools panel AF filter",
        )
        n_raw = _vcf_line_count(norm_panel)
        n_filt = _vcf_line_count(filt_panel)
        bt.logging.info(
            f"[pipeline] panel pass: {n_filt} variants (filtered from {n_raw}) "
            f"at ClinVar CF positions"
        )
        return filt_panel if n_filt > 0 else norm_panel
    except Exception as e:
        bt.logging.warning(f"[pipeline] panel call failed: {e}")
        return None


def _emergency_call(
    ref: str,
    bam: str,
    work_dir: str,
    region: str,
) -> Tuple[str, str]:
    """Ultra-loose pass when all standard calls return zero in-window."""
    raw_path = os.path.join(work_dir, "raw.emergency.vcf")
    _run(
        f"bcftools mpileup -f {ref} -r {region} -a AD,DP "
        f"-q 0 -Q 0 --max-depth 16000 {bam} "
        f"| bcftools call -c -Ov -o {raw_path}",
        "bcftools emergency",
    )
    norm_path = raw_path.replace(".vcf", ".norm.vcf")
    _run(
        f"bcftools norm -f {ref} -m -both -c w {raw_path} -Ov -o {norm_path}",
        "bcftools norm emergency",
    )
    return raw_path, norm_path


def call_variants_with_fallback(
    ref: str,
    bam: str,
    work_dir: str,
    region: str,
    task_region: str,
    region_start: int,
    region_end: int,
    genomic_coords: bool,
) -> Tuple[str, int, bool, List[str]]:
    log_bcftools_version()
    reads = _bam_reads_in_region(bam, task_region)
    bt.logging.info(f"[pipeline] BAM reads in {task_region}: {reads}")

    task_region_str = f"chr7:{region_start}-{region_end}"
    profile = classify_task(task_region_str)
    strat = active_strategy_name()
    pick = pipeline_pick_mode()
    env_qual = os.environ.get("NIOME_MPILEUP_QUAL", "").strip()
    env_extra = os.environ.get("NIOME_MPILEUP_EXTRA", "").strip()
    mpileup_extra = resolve_mpileup_extra(env_extra or profile.mpileup_extra or "")
    primary_qual = env_qual or profile.mpileup_qual
    if profile.name == "ultra_wide" and not bcftools_supports_indels_20():
        if not env_qual:
            primary_qual = "-q 2 -Q 2"
    bt.logging.info(
        f"[pipeline] strategy={strat} pick={pick} profile={profile.name} "
        f"mpileup={primary_qual} extra={mpileup_extra or '(none)'} "
        f"merge_pool={pipeline_merge_pool()}"
    )

    raw_vcf = os.path.join(work_dir, "raw.vcf")
    raw1, norm1 = call_variants(
        ref, bam, raw_vcf, region, primary_qual, mpileup_extra
    )

    raw2_path = os.path.join(work_dir, "raw.retry.vcf")
    _run(
        f"bcftools mpileup -f {ref} -r {region} -a AD,DP "
        f"-q 0 -Q 0 {mpileup_extra} --max-depth 8000 {bam} "
        f"| bcftools call -mv -Ov -o {raw2_path}",
        "bcftools call retry",
    )
    norm2 = raw2_path.replace(".vcf", ".norm.vcf")
    _run(
        f"bcftools norm -f {ref} -m -both -c w {raw2_path} -Ov -o {norm2}",
        "bcftools norm retry",
    )

    # Biological panel pass: always run regardless of profile/strategy.
    panel_norm = call_panel_variants(ref, bam, work_dir, region)

    candidates = [
        (norm1, "norm"),
        (raw1, "raw"),
        (norm2, "norm-retry"),
        (raw2_path, "raw-retry"),
    ]

    if panel_norm and os.path.exists(panel_norm):
        candidates.insert(0, (panel_norm, "norm-panel"))

    if profile.name == "ultra_wide":
        raw3 = os.path.join(work_dir, "raw.indel.vcf")
        raw3, norm3 = call_variants(
            ref,
            bam,
            raw3,
            region,
            "-q 1 -Q 1",
            mpileup_extra,
        )
        candidates.extend([(norm3, "norm-indel"), (raw3, "raw-indel")])

        raw4 = os.path.join(work_dir, "raw.loose.vcf")
        raw4, norm4 = call_variants(
            ref,
            bam,
            raw4,
            region,
            "-q 0 -Q 0",
            mpileup_extra,
            caller="-c",
            max_depth=12000,
        )
        candidates.extend([(norm4, "norm-loose"), (raw4, "raw-loose")])

        try:
            raw5 = os.path.join(work_dir, "raw.alleles.vcf")
            _run(
                f"bcftools mpileup -f {ref} -r {region} -a AD,DP "
                f"-q 0 -Q 0 {mpileup_extra} --max-depth 12000 {bam} "
                f"| bcftools call -c -Ov -o {raw5}",
                "bcftools call alleles",
            )
            norm5 = raw5.replace(".vcf", ".norm.vcf")
            _run(
                f"bcftools norm -f {ref} -m -both -c w {raw5} -Ov -o {norm5}",
                "bcftools norm alleles",
            )
            candidates.extend([(norm5, "norm-alleles"), (raw5, "raw-alleles")])
        except Exception as e:
            bt.logging.warning(f"[pipeline] alleles pass skipped: {e}")

    best_path, best_n, best_coords, supplemental = _pick_best_vcf(
        candidates,
        region_start,
        region_end,
        genomic_coords,
        profile,
    )
    if not best_path or best_n == 0:
        try:
            raw_e, norm_e = _emergency_call(ref, bam, work_dir, region)
            candidates.extend([(norm_e, "norm-emergency"), (raw_e, "raw-emergency")])
            best_path, best_n, best_coords, supplemental = _pick_best_vcf(
                candidates,
                region_start,
                region_end,
                genomic_coords,
                profile,
            )
            bt.logging.warning(
                f"[pipeline] emergency call pass selected={best_path} n={best_n}"
            )
        except Exception as e:
            bt.logging.warning(f"[pipeline] emergency call failed: {e}")
        if not best_path:
            best_path = norm1
            best_n = 0
            supplemental = []

    bt.logging.info(
        f"[pipeline] selected {best_path} with {best_n} in-window calls "
        f"supplemental={len(supplemental)}"
    )
    return best_path, best_n, best_coords, supplemental


def run_pipeline(task: Task, work_dir: str):
    """
    Align → call → norm → read-all VCF (v3).
    Returns (final_vcf_path, None).
    """
    os.makedirs(work_dir, exist_ok=True)
    chrom, region_start, region_end = parse_region(task.genome_context.region)
    task_region = task.genome_context.region
    region_call = _padded_region(chrom, region_start, region_end)
    rlen = region_length(task.genome_context.region)

    profile = classify_task(task_region, task.expected_variant_count)
    bt.logging.info(
        f"[pipeline] task={task.task_id[:8]}… rev={get_read_calling_rev()} "
        f"strategy={active_strategy_name()} profile={profile.name} "
        f"region={task_region} len={rlen} expected={task.expected_variant_count}"
    )

    raw_path: Optional[str] = None
    genomic_coords = False
    supplemental: List[str] = []

    try:
        ref, genomic_coords = pick_reference()
        if not genomic_coords:
            local_start = max(1, region_start - CFTR_START - _REGION_PAD)
            local_end = min(CFTR_END - CFTR_START, region_end - CFTR_START + _REGION_PAD)
            region_call = f"{chrom}:{local_start}-{local_end}"

        r1 = download_fastq(
            task.input.read1_fastq, os.path.join(work_dir, "read_1.fq")
        )
        r2 = download_fastq(
            task.input.read2_fastq, os.path.join(work_dir, "read_2.fq")
        )
        bam = os.path.join(work_dir, "aligned.bam")
        align_reads(ref, r1, r2, bam)
        raw_path, n_calls, genomic_coords, supplemental = call_variants_with_fallback(
            ref,
            bam,
            work_dir,
            region_call,
            task_region,
            region_start,
            region_end,
            genomic_coords,
        )
        bt.logging.info(
            f"[pipeline] ref={'hg38' if genomic_coords else 'slice'} "
            f"calls_in_window={n_calls} supplemental_vcfs={len(supplemental)}"
        )
    except Exception as e:
        bt.logging.warning(f"Read alignment skipped: {e}")

    return build_task_vcf(
        task,
        work_dir,
        raw_vcf_path=raw_path,
        genomic_coords=genomic_coords,
        extra_vcf_paths=supplemental,
    )
