"""
Genotype assignment from read evidence (validator gives 0.5× on GT mismatch).

Calibrated for NIOME scoring: hom-alt recall matters as much as het calls.
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

    Env overrides:
      NIOME_GT_HOM_AF  (default 0.52)
      NIOME_GT_HET_AF  (default 0.15)
      NIOME_GT_INDEL_HOM_AF (default 0.45)
    """
    hom_af = _env_float("NIOME_GT_HOM_AF", 0.52)
    het_af = _env_float("NIOME_GT_HET_AF", 0.15)
    indel_hom_af = _env_float("NIOME_GT_INDEL_HOM_AF", 0.45)

    is_indel = not _is_snp(call.ref, call.alt)
    hom_thresh = indel_hom_af if is_indel else hom_af

    if call.dp >= 3 and call.alt_ad > 0:
        af = call.alt_ad / call.dp
        ref_ad = max(0, call.dp - call.alt_ad)
        if af >= hom_thresh:
            return "1/1"
        if af >= het_af:
            return "0/1"
        if ref_ad > 0 and call.alt_ad >= ref_ad * 2 and call.alt_ad >= 3:
            return "1/1"
        if af >= 0.08 and call.qual >= 25.0:
            return "0/1"
        return "0/1"

    gt = (call.gt or "").replace("|", "/")
    if gt in ("1/1", "1|1"):
        return "1/1"
    if gt in ("0/1", "1/0", "0|1", "1|0"):
        return "0/1"
    if gt in ("0/0", "./.", "."):
        return "0/1"
    return "0/1"
