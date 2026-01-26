"""
Genotype assignment from read evidence (validator gives 0.5× on GT mismatch).

v10: conservative het/hom — v9 over-called 1/1 (AF≥0.52) and hurt vs v5/v8 on 5.23.02.
"""

import os
from typing import Any


def _env_float(name: str, default: float) -> float:
    raw = os.environ.get(name, "").strip()
    if not raw:
        return default
    try:
        return float(raw)
    except ValueError:
        return default


def _is_snp(ref: str, alt: str) -> bool:
    return (
        len(ref) == 1
        and len(alt) == 1
        and ref not in (".", "N")
        and alt not in (".", "N")
    )


def gt_from_read_call(call: Any) -> str:
    """
    Assign 0/1 vs 1/1 from AD/DP; fall back to mpileup GT when AD thin.

    Defaults match pre-v9 behavior (5.22 v5 / v8 panels on 5.23.02):
      NIOME_GT_HOM_AF=0.58  NIOME_GT_HET_AF=0.20
    """
    hom_af = _env_float("NIOME_GT_HOM_AF", 0.58)
    het_af = _env_float("NIOME_GT_HET_AF", 0.20)

    if call.dp >= 4 and call.alt_ad > 0:
        af = call.alt_ad / call.dp
        if af >= hom_af:
            return "1/1"
        if af >= het_af:
            return "0/1"
        return "0/1"

    gt = (call.gt or "").replace("|", "/")
    if gt in ("1/1", "1|1"):
        return "1/1"
    if gt in ("0/1", "1/0", "0|1", "1|0"):
        return "0/1"
    return "0/1"
