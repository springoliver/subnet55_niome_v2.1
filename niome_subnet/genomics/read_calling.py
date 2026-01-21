"""
NIOME miner: read-evidence VCF via NIOME Native selection (proprietary method).

New subnet rules (expected_variant_count hidden / 0):
  - VCF: call from read evidence, adaptive count from tiers (not fixed N)
  - Truth may include non-ClinVar variants and indels
  - cftr_annotations: ClinVar-mapped variants in the submitted VCF only
"""

import os
import subprocess
from typing import Any, Dict, List, Optional, Set, Tuple

import bittensor as bt

from niome_subnet.genomics.cftr_lookup import ensure_clinvar_db
from niome_subnet.genomics.evidence_selection import METHOD_ID, native_select_variants
from niome_subnet.genomics.model import Task
from niome_subnet.genomics.read_types import ReadCall, gt_from_read_call
from niome_subnet.genomics.task_profile import (
    TaskProfile,
    classify_task,
    curriculum_target_for_profile,
    is_simple_snp,
    parse_region,
    region_length,
    ultra_scoring_core,
)

READ_CALLING_REV = METHOD_ID

_CHR7_LENGTH = 159345973
_CFTR_START = 117430000
_MAX_ALLELE_LEN = 64
_MICRO_PAD = 25000


def _is_callable(ref: str, alt: str) -> bool:
    if not ref or not alt or ref == "." or alt == ".":
        return False
    if alt.startswith("<") or alt == "*":
        return len(ref) <= _MAX_ALLELE_LEN
    return len(ref) <= _MAX_ALLELE_LEN and len(alt) <= _MAX_ALLELE_LEN


def _chrom_ok(chrom: str) -> bool:
    return chrom in ("chr7", "7", "Chr7", "CHR7")


def _parse_gt(parts: List[str]) -> str:
    if len(parts) < 10:
        return "0/1"
    fmt = parts[8].split(":")
    sample = parts[9].split(":")
    if "GT" not in fmt:
        return "0/1"
    gt = sample[fmt.index("GT")].split("/")
    if gt[0] not in (".", "") and gt[-1] not in (".", ""):
        return "/".join(gt[:2])
    return "0/1"


def _parse_ad_dp(parts: List[str]) -> Tuple[int, int]:
    if len(parts) < 10:
        return 0, 0
    fmt = parts[8].split(":")
    sample = parts[9].split(":")
    dp = 0
    alt_ad = 0
    if "DP" in fmt:
        try:
            dp = int(sample[fmt.index("DP")])
        except (ValueError, IndexError):
            dp = 0
    if "AD" in fmt:
        try:
            ads = [int(x) for x in sample[fmt.index("AD")].split(",") if x != "."]
            if len(ads) >= 2:
                alt_ad = max(ads[1:])
            elif ads:
                alt_ad = ads[0]
        except (ValueError, IndexError):
            alt_ad = 0
    return alt_ad, dp


def _infer_genomic_coords(vcf_path: str) -> bool:
    if not vcf_path or not os.path.exists(vcf_path):
        return True
    with open(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            try:
                return int(line.split("\t")[1]) > 1_000_000
            except (IndexError, ValueError):
                continue
    return True


def parse_read_calls(
    vcf_path: Optional[str],
    region_start: int,
    region_end: int,
    genomic_coords: bool = True,
    padded: bool = False,
) -> List[ReadCall]:
    """Parse mpileup/call VCF into ReadCall list within task window."""
    if not vcf_path or not os.path.exists(vcf_path):
        return []

    lo = max(1, region_start - _MICRO_PAD) if padded else region_start
    hi = region_end + _MICRO_PAD if padded else region_end

    calls: List[ReadCall] = []
    with open(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            if not _chrom_ok(parts[0]):
                continue
            pos_abs = int(parts[1]) if genomic_coords else int(parts[1]) + _CFTR_START
            if not (lo <= pos_abs <= hi):
                continue
            ref, alt = parts[3].upper(), parts[4].split(",")[0].upper()
            if not _is_callable(ref, alt):
                continue
            filt = parts[6] if len(parts) > 6 else "PASS"
            pass_filter = filt in ("PASS", ".")
            try:
                qual = float(parts[5])
            except (IndexError, ValueError):
                qual = 0.0
            alt_ad, dp = _parse_ad_dp(parts)
            calls.append(
                ReadCall(
                    pos=pos_abs,
                    ref=ref,
                    alt=alt,
                    qual=qual,
                    gt=_parse_gt(parts),
                    alt_ad=alt_ad,
                    dp=dp,
                    pass_filter=pass_filter,
                )
            )

    calls.sort(key=lambda c: (c.pass_filter, c.qual, c.alt_ad, c.dp), reverse=True)
    return calls


def count_indels_in_calls(calls: List[ReadCall]) -> int:
    return sum(1 for c in calls if not is_simple_snp(c.ref, c.alt))


def vcf_raw_stats_in_region(
    vcf_path: Optional[str],
    region_start: int,
    region_end: int,
    genomic_coords: bool = True,
) -> Tuple[int, int]:
    """Return (total_alleles, indel_alleles) in window before native selection."""
    calls = parse_read_calls(vcf_path, region_start, region_end, genomic_coords)
    if not calls and vcf_path:
        calls = parse_read_calls(
            vcf_path, region_start, region_end, genomic_coords, padded=False
        )
    return len(calls), count_indels_in_calls(calls)


def merge_read_call_pools(
    pools: List[List[ReadCall]],
) -> List[ReadCall]:
    """Dedupe by (pos, ref, alt); keep highest QUAL+AD support."""
    best: Dict[Tuple[int, str, str], ReadCall] = {}
    for pool in pools:
        for call in pool:
            key = (call.pos, call.ref.upper(), call.alt.upper())
            prev = best.get(key)
            score = call.qual + call.alt_ad * 4
            if prev is None or score > prev.qual + prev.alt_ad * 4:
                best[key] = call
    out = list(best.values())
    out.sort(key=lambda c: c.pos)
    return out


def count_calls_in_region(
    vcf_path: Optional[str],
    region_start: int,
    region_end: int,
    genomic_coords: Optional[bool] = None,
    profile: Optional[TaskProfile] = None,
) -> Tuple[int, bool]:
    if genomic_coords is None:
        genomic_coords = _infer_genomic_coords(vcf_path or "")
    region = f"chr7:{region_start}-{region_end}"
    if profile is None:
        profile = classify_task(region)
    kept = select_read_variants(
        parse_read_calls(vcf_path, region_start, region_end, genomic_coords),
        region_start,
        region_end,
        profile,
        None,
    )
    return len(kept), genomic_coords


def _load_clinvar_ids(region: str) -> Dict[Tuple[int, str, str], str]:
    """Map (pos, ref, alt) -> ClinVar variation ID for VCF ID column."""
    chrom, region_start, region_end = parse_region(region)
    out: Dict[Tuple[int, str, str], str] = {}
    try:
        db = ensure_clinvar_db()
    except Exception as e:
        bt.logging.warning(f"[read_calling] ClinVar unavailable: {e}")
        return out

    for reg in (region, region.replace("chr7", "7")):
        result = subprocess.run(
            f"bcftools view -r {reg} {db}",
            shell=True,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            continue
        for line in result.stdout.splitlines():
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            pos = int(parts[1])
            if not (region_start <= pos <= region_end):
                continue
            vid = parts[2].split(";")[0] if parts[2] != "." else "."
            ref = parts[3]
            for alt in parts[4].split(","):
                if _is_callable(ref, alt):
                    out[(pos, ref, alt)] = vid
        if out:
            break

    if not out:
        pad = _MICRO_PAD
        padded = f"{chrom}:{max(1, region_start - pad)}-{region_end + pad}"
        result = subprocess.run(
            f"bcftools view -r {padded} {db}",
            shell=True,
            capture_output=True,
            text=True,
        )
        if result.returncode == 0:
            for line in result.stdout.splitlines():
                if line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 5:
                    continue
                pos = int(parts[1])
                if not (region_start <= pos <= region_end):
                    continue
                vid = parts[2].split(";")[0] if parts[2] != "." else "."
                ref = parts[3]
                for alt in parts[4].split(","):
                    if _is_callable(ref, alt):
                        out.setdefault((pos, ref, alt), vid)

    return out


def select_read_variants(
    calls: List[ReadCall],
    region_start: int,
    region_end: int,
    profile: Optional[TaskProfile] = None,
    clinvar_ids: Optional[Dict[Tuple[int, str, str], str]] = None,
) -> List[ReadCall]:
    """NIOME Native: adaptive read-evidence selection (see evidence_selection.py)."""
    if profile is None:
        profile = classify_task(f"chr7:{region_start}-{region_end}")
    return native_select_variants(
        calls, region_start, region_end, profile, clinvar_ids
    )


def format_vcf(calls: List[ReadCall], clinvar_ids: Dict[Tuple[int, str, str], str]) -> str:
    use_dot_id = os.environ.get("NIOME_VCF_DOT_ID", "1").strip() not in ("0", "false", "no")
    lines = [
        "##fileformat=VCFv4.2",
        f"##source=niome_miner_{READ_CALLING_REV}",
        f"##contig=<ID=chr7,length={_CHR7_LENGTH}>",
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele fraction">',
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    ]

    minimal = os.environ.get("NIOME_VCF_MINIMAL", "0").strip() in ("1", "true", "yes")

    for call in calls:
        call.ref = call.ref.upper()
        call.alt = call.alt.upper()
        key = (call.pos, call.ref, call.alt)
        vid = clinvar_ids.get(key, ".")
        if vid == ".":
            vid = clinvar_ids.get(
                (call.pos, call.ref.lower(), call.alt.lower()), "."
            )
        out_id = "." if use_dot_id or not vid or vid == "." else vid
        qual_str = f"{call.qual:.3f}" if call.qual > 0 else "."
        if minimal:
            lines.append(
                f"chr7\t{call.pos}\t{out_id}\t{call.ref}\t{call.alt}\t"
                f"{qual_str}\tPASS\t.\tGT\t{call.gt}"
            )
            continue
        if call.dp > 0 and call.alt_ad > 0:
            ref_ad = max(0, call.dp - call.alt_ad)
            af = call.alt_ad / call.dp
            info = f"DP={call.dp};AF={af:.4f}"
            fmt = "GT:DP:AD"
            sample = f"{call.gt}:{call.dp}:{ref_ad},{call.alt_ad}"
        elif call.dp > 0:
            info = f"DP={call.dp}"
            fmt = "GT:DP"
            sample = f"{call.gt}:{call.dp}"
        else:
            info = "."
            fmt = "GT"
            sample = call.gt
        lines.append(
            f"chr7\t{call.pos}\t{out_id}\t{call.ref}\t{call.alt}\t{qual_str}\tPASS\t"
            f"{info}\t{fmt}\t{sample}"
        )
    return "\n".join(lines) + "\n"


def build_task_vcf(
    task: Task,
    work_dir: str,
    raw_vcf_path: Optional[str] = None,
    genomic_coords: bool = True,
    extra_vcf_paths: Optional[List[str]] = None,
) -> Tuple[str, Optional[Dict[str, Any]]]:
    """
    Build final VCF from read calls. Returns (path, None) — annotations via cftr_lookup.

    extra_vcf_paths: optional supplemental mpileup VCFs (e.g. indel-rich retry) merged
    before native selection.
    """
    region = task.genome_context.region
    _, region_start, region_end = parse_region(region)
    rlen = region_length(region)
    profile = classify_task(region, task.expected_variant_count)

    pools: List[List[ReadCall]] = []
    for vcf in [raw_vcf_path] + list(extra_vcf_paths or []):
        if not vcf:
            continue
        chunk = parse_read_calls(
            vcf, region_start, region_end, genomic_coords, padded=True
        )
        if not chunk:
            chunk = parse_read_calls(
                vcf, region_start, region_end, genomic_coords, padded=False
            )
        if chunk:
            pools.append(chunk)

    calls = merge_read_call_pools(pools)
    if pools:
        bt.logging.info(
            f"[read_calling] merged_pool={len(calls)} from {len(pools)} vcf source(s)"
        )

    clinvar_ids = _load_clinvar_ids(region)
    selected = select_read_variants(
        calls, region_start, region_end, profile, clinvar_ids
    )

    for call in selected:
        key = (call.pos, call.ref, call.alt)
        if key in clinvar_ids:
            call.clinvar_id = clinvar_ids[key]

    vcf_content = format_vcf(selected, clinvar_ids)
    out_path = os.path.join(work_dir, "final.vcf")
    with open(out_path, "w") as fh:
        fh.write(vcf_content)

    n_read_gt = sum(1 for c in selected if c.dp > 0 and c.alt_ad > 0)
    in_core = sum(1 for c in selected if ultra_scoring_core(c.pos))
    n_indel_raw = count_indels_in_calls(calls)
    n_indel_sub = count_indels_in_calls(selected)
    target_n = curriculum_target_for_profile(profile)
    bt.logging.info(
        f"[read_calling] task={task.task_id[:8]}… rev={READ_CALLING_REV} "
        f"profile={profile.name} region={region} len={rlen} "
        f"curriculum_target={target_n} "
        f"raw_calls={len(calls)} indels_raw={n_indel_raw} "
        f"submitted={len(selected)} indels={n_indel_sub} core={in_core} "
        f"read_gt={n_read_gt} snps={sum(1 for c in selected if is_simple_snp(c.ref, c.alt))} "
        f"clinvar_in_vcf={sum(1 for c in selected if c.clinvar_id != '.')}"
    )

    if target_n > 0 and len(calls) < target_n:
        bt.logging.warning(
            f"[read_calling] pool_shortfall: merged_pool={len(calls)} "
            f"< curriculum_target={target_n} (cannot invent sites; need more calling passes)"
        )
    elif target_n > 0 and len(selected) < target_n:
        bt.logging.warning(
            f"[read_calling] selection_shortfall: submitted={len(selected)} "
            f"< curriculum_target={target_n} pool={len(calls)}"
        )

    if len(selected) == 0:
        bt.logging.error(
            f"[read_calling] ZERO variants for {task.task_id} — check BAM/VCF in {work_dir}"
        )

    return out_path, None
