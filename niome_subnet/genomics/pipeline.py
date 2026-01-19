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
    READ_CALLING_REV,
    build_task_vcf,
    count_calls_in_region,
    parse_region,
    region_length,
    vcf_raw_stats_in_region,
)
from niome_subnet.genomics.task_profile import TaskProfile, classify_task

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
) -> Tuple[str, str]:
    qual = mpileup_qual or os.environ.get("NIOME_MPILEUP_QUAL", "-q 1 -Q 1")
    extra = mpileup_extra or os.environ.get("NIOME_MPILEUP_EXTRA", "")
    _run(
        f"bcftools mpileup -f {ref} -r {region} -a AD,DP "
        f"{qual} {extra} --max-depth 8000 {bam} "
        f"| bcftools call -mv -Ov -o {raw_vcf}",
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
    if profile.name == "ultra_wide":
        for _, path, n_sel, _, label, n_indel in scored[1:]:
            if path == best_path or n_indel == 0:
                continue
            if n_indel > best_indel:
                supplemental.append(path)
                bt.logging.info(
                    f"[pipeline] supplemental merge {label} ({path}) indels={n_indel}"
                )

    return best_path, best_n, best_coords, supplemental


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
    reads = _bam_reads_in_region(bam, task_region)
    bt.logging.info(f"[pipeline] BAM reads in {task_region}: {reads}")

    task_region_str = f"chr7:{region_start}-{region_end}"
    profile = classify_task(task_region_str)
    bt.logging.info(
        f"[pipeline] profile={profile.name} mpileup={profile.mpileup_qual} "
        f"extra={profile.mpileup_extra or '(none)'}"
    )

    raw_vcf = os.path.join(work_dir, "raw.vcf")
    raw1, norm1 = call_variants(
        ref, bam, raw_vcf, region, profile.mpileup_qual, profile.mpileup_extra
    )

    raw2_path = os.path.join(work_dir, "raw.retry.vcf")
    retry_extra = profile.mpileup_extra or "--indels-2.0"
    _run(
        f"bcftools mpileup -f {ref} -r {region} -a AD,DP "
        f"-q 0 -Q 0 {retry_extra} --max-depth 8000 {bam} "
        f"| bcftools call -mv -Ov -o {raw2_path}",
        "bcftools call retry",
    )
    norm2 = raw2_path.replace(".vcf", ".norm.vcf")
    _run(
        f"bcftools norm -f {ref} -m -both -c w {raw2_path} -Ov -o {norm2}",
        "bcftools norm retry",
    )

    candidates = [
        (norm1, "norm"),
        (raw1, "raw"),
        (norm2, "norm-retry"),
        (raw2_path, "raw-retry"),
    ]

    if profile.name == "ultra_wide":
        raw3 = os.path.join(work_dir, "raw.indel.vcf")
        raw3, norm3 = call_variants(
            ref,
            bam,
            raw3,
            region,
            "-q 1 -Q 1",
            profile.mpileup_extra or "--indels-2.0",
        )
        candidates.extend([(norm3, "norm-indel"), (raw3, "raw-indel")])

    best_path, best_n, best_coords, supplemental = _pick_best_vcf(
        candidates,
        region_start,
        region_end,
        genomic_coords,
        profile,
    )
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
        f"[pipeline] task={task.task_id[:8]}… rev={READ_CALLING_REV} "
        f"profile={profile.name} region={task_region} len={rlen} "
        f"expected={task.expected_variant_count}"
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
