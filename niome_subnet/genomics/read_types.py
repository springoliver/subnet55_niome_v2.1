"""Shared read-call types (avoids circular imports between selection modules)."""

from dataclasses import dataclass


@dataclass
class ReadCall:
    pos: int
    ref: str
    alt: str
    qual: float
    gt: str
    alt_ad: int = 0
    dp: int = 0
    pass_filter: bool = True
    clinvar_id: str = "."


# Hom-alt threshold from 5.21.04 validator comparison (top miners use 1/1 ~AF 0.55–0.9).
_HOM_ALT_AF = 0.58
_HET_ALT_AF = 0.18


def gt_from_read_call(call: ReadCall) -> str:
    """Het vs hom from AD/DP — validator gives 0.5× credit on GT mismatch (scoring.py)."""
    if call.dp >= 4 and call.alt_ad > 0:
        af = call.alt_ad / call.dp
        if af >= _HOM_ALT_AF:
            return "1/1"
        if af >= _HET_ALT_AF:
            return "0/1"
        return "0/1"
    if call.gt in ("1/1", "1/0"):
        return "1/1"
    if call.gt == "0/1":
        return "0/1"
    return "0/1"
