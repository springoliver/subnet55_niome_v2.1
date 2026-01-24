"""Shared read-call types (avoids circular imports between selection modules)."""

from dataclasses import dataclass

from niome_subnet.genomics.gt_tuning import gt_from_read_call

__all__ = ["ReadCall", "gt_from_read_call"]


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
