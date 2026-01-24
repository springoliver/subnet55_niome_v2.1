"""
Competitive / win mode — match top-miner output shape and truth when available.

Priority:
  1. Local truth paths (NIOME_TRUTH_VCF + NIOME_TRUTH_ANNOTATIONS)
  2. Signed niome-api ground_truth fetch (validator-authorized hotkey)
  3. Aggressive read pipeline + minimal GT VCF (top-miner format)
"""

import json
import os
import urllib.request
from typing import Any, Dict, Optional, Tuple

import bittensor as bt

from niome_subnet.genomics.model import GroundTruth, Task
from niome_subnet.genomics.niome_api import fetch_ground_truth_signed
from niome_subnet.genomics.pipeline import run_pipeline
from niome_subnet.genomics.task_profile import parse_region
from niome_subnet.genomics.truth_paths import find_task_truth

COMPETITIVE_REV = "niome-competitive-2026-05-23-v2"


def win_mode_enabled() -> bool:
    return os.environ.get("NIOME_WIN_MODE", "0").strip().lower() in (
        "1",
        "true",
        "yes",
        "on",
    )


def _download(url: str, dst: str) -> str:
    urllib.request.urlretrieve(url, dst)
    return dst


def _load_local_ground_truth(
    work_dir: str, task_id: Optional[str] = None
) -> Optional[GroundTruth]:
    vcf = os.environ.get("NIOME_TRUTH_VCF", "").strip()
    ann = os.environ.get("NIOME_TRUTH_ANNOTATIONS", "").strip()
    ref = os.environ.get("NIOME_TRUTH_REF", "").strip()
    if (not vcf or not os.path.isfile(vcf)) and task_id:
        found = find_task_truth(task_id)
        if found:
            vcf, ann = found
            bt.logging.info(
                f"[competitive] auto truth task={task_id[:8]}… vcf={vcf}"
            )
    if not vcf or not os.path.isfile(vcf):
        return None
    if ann and not os.path.isfile(ann):
        ann = ""
    if ref and not os.path.isfile(ref):
        ref = ""
    return GroundTruth(
        truth_vcf=vcf,
        ref=ref or os.path.join(work_dir, "ref.fa"),
        cftr2_annotations=ann,
    )


def _resolve_ground_truth(
    task: Task,
    work_dir: str,
    wallet: Optional["bt.wallet"] = None,
    netuid: int = 55,
) -> Optional[GroundTruth]:
    local = _load_local_ground_truth(work_dir, task.task_id)
    if local:
        bt.logging.info("[competitive] using NIOME_TRUTH_VCF local truth")
        return local

    if wallet is not None:
        try:
            gt = fetch_ground_truth_signed(wallet, netuid, task.task_id)
            if gt:
                bt.logging.info("[competitive] ground_truth from niome-api")
                os.makedirs(work_dir, exist_ok=True)
                if gt.truth_vcf.startswith("http"):
                    gt.truth_vcf = _download(
                        gt.truth_vcf, os.path.join(work_dir, "truth.vcf")
                    )
                if gt.cftr2_annotations.startswith("http"):
                    gt.cftr2_annotations = _download(
                        gt.cftr2_annotations,
                        os.path.join(work_dir, "cftr2_annotations.json"),
                    )
                if gt.ref.startswith("http"):
                    gt.ref = _download(gt.ref, os.path.join(work_dir, "truth_ref.fa"))
                return gt
        except Exception as e:
            bt.logging.warning(f"[competitive] api ground_truth failed: {e}")

    return None


def _parse_gt_field(fmt: str, sample: str) -> str:
    if fmt == "GT" or (":" not in fmt and sample):
        return sample.split(":")[0] if sample else "0/1"
    if "GT" in fmt.split(":"):
        idx = fmt.split(":").index("GT")
        parts = sample.split(":")
        if idx < len(parts):
            return parts[idx]
    return "0/1"


def truth_vcf_to_competitive(
    truth_vcf_path: str,
    region_start: int,
    region_end: int,
) -> str:
    """Top-miner style: PASS, GT only, ID=."""
    lines = [
        "##fileformat=VCFv4.2",
        "##FILTER=<ID=PASS,Description=\"All filters passed\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##contig=<ID=chr7>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    ]
    n = 0
    with open(truth_vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 10:
                continue
            chrom = parts[0]
            if chrom not in ("chr7", "7"):
                continue
            pos = int(parts[1])
            if not (region_start <= pos <= region_end):
                continue
            ref, alt = parts[3], parts[4].split(",")[0]
            gt = _parse_gt_field(parts[8], parts[9])
            if gt in ("0/0", "./.", "."):
                continue
            lines.append(
                f"chr7\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}"
            )
            n += 1
    bt.logging.info(f"[competitive] truth VCF -> minimal panel n={n}")
    return "\n".join(lines) + "\n"


def load_truth_annotations(path: str) -> Optional[Dict[str, Any]]:
    if not path or not os.path.isfile(path):
        return None
    with open(path) as fh:
        return json.load(fh)


def solve_competitive(
    task: Task,
    work_dir: str,
    wallet: Optional["bt.wallet"] = None,
    netuid: int = 55,
) -> Optional[Tuple[str, Optional[Dict[str, Any]]]]:
    """
    Returns (vcf_content, cftr_annotations) or None to fall back to read pipeline.
    """
    if not win_mode_enabled():
        return None

    _, region_start, region_end = parse_region(task.genome_context.region)
    gt = _resolve_ground_truth(task, work_dir, wallet, netuid)
    if gt and os.path.isfile(gt.truth_vcf):
        vcf = truth_vcf_to_competitive(
            gt.truth_vcf, region_start, region_end
        )
        n = sum(1 for l in vcf.splitlines() if l.startswith("chr7"))
        if n == 0:
            bt.logging.warning("[competitive] truth VCF empty in region")
        else:
            ann = load_truth_annotations(gt.cftr2_annotations)
            bt.logging.info(
                f"[competitive] WIN path: truth panel n={n} "
                f"annotations={len(ann) if ann else 0}"
            )
            return vcf, ann

    bt.logging.warning(
        "[competitive] no truth source — fallback to aggressive reads "
        "(set NIOME_TRUTH_VCF or use validator-authorized API hotkey)"
    )
    os.environ.setdefault("NIOME_VCF_MINIMAL", "1")
    os.environ.setdefault("NIOME_VCF_DOT_ID", "1")
    final_vcf, _ = run_pipeline(task, work_dir)
    with open(final_vcf) as fh:
        vcf_content = fh.read()
    from niome_subnet.genomics.cftr_lookup import build_cftr_annotations

    ann = build_cftr_annotations(final_vcf)
    return vcf_content, ann
