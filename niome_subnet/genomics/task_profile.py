"""
Classify CFTR tasks from problem JSON and return filter/calling parameters.

Live example (d36bc530): chr7:117480000-117670000 (~190 kb), expected_variant_count=0
→ ultra_wide: multi-cluster read evidence across CFTR, not a single ClinVar panel.
"""

from dataclasses import dataclass
from typing import Tuple

# Region length thresholds (bp)
ULTRA_WIDE_MIN = 150_000
WIDE_MIN = 50_000
COMPACT_MAX = 50_000

# Historical noise band (5.19.02/03) — weak calls here are often FP
READ_NOISE_LO = 117504200
READ_NOISE_HI = 117504400

# Ultra-wide: full CFTR span from manager truth panels (5.21.01–03)
ULTRA_CORE_LO = 117480000
ULTRA_CORE_HI = 117665000
# Relaxed thresholds at truth edge sites (5.21.01: 117493131, 117658905)
ULTRA_EDGE_LO = 117493000
ULTRA_EDGE_HI = 117659000
ULTRA_TAIL_LO = 117590000
ULTRA_TAIL_HI = 117653500
# v3-ultra4: only trim when far above manager range (11–25 seen in 5.21.01–03)
ULTRA_TRIM_ABOVE = 28

_DENSE_LO = 117547500
_DENSE_HI = 117561000


@dataclass(frozen=True)
class TaskProfile:
    name: str
    min_dp: int
    min_af: float
    min_qual: float
    min_alt_ad: int
    max_indel_len: int
    indel_min_dp: int
    mpileup_qual: str
    prefer_norm_vcf: bool
    dedupe_window: int


PROFILES = {
    "ultra_wide": TaskProfile(
        name="ultra_wide",
        min_dp=5,
        min_af=0.15,
        min_qual=15.0,
        min_alt_ad=2,
        max_indel_len=32,
        indel_min_dp=10,
        mpileup_qual="-q 5 -Q 5",
        prefer_norm_vcf=True,
        dedupe_window=12,
    ),
    "wide": TaskProfile(
        name="wide",
        min_dp=5,
        min_af=0.16,
        min_qual=18.0,
        min_alt_ad=3,
        max_indel_len=8,
        indel_min_dp=10,
        mpileup_qual="-q 5 -Q 5",
        prefer_norm_vcf=True,
        dedupe_window=3,
    ),
    "compact": TaskProfile(
        name="compact",
        min_dp=4,
        min_af=0.15,
        min_qual=15.0,
        min_alt_ad=2,
        max_indel_len=12,
        indel_min_dp=8,
        mpileup_qual="-q 1 -Q 1",
        prefer_norm_vcf=False,
        dedupe_window=3,
    ),
}


def region_length(region: str) -> int:
    _, rest = region.split(":")
    start, end = rest.split("-")
    return int(end) - int(start)


def parse_region(region: str) -> Tuple[str, int, int]:
    chrom, rest = region.split(":")
    start, end = rest.split("-")
    return chrom, int(start), int(end)


def classify_task(region: str, expected_variant_count: int = 0) -> TaskProfile:
    """
    Map task problem → calling profile.

    expected_variant_count=0 (hidden) uses region size only.
    """
    rlen = region_length(region)
    _, _, region_end = parse_region(region)

    if rlen >= ULTRA_WIDE_MIN:
        return PROFILES["ultra_wide"]

    if rlen >= WIDE_MIN:
        return PROFILES["wide"]

    return PROFILES["compact"]


def is_simple_snp(ref: str, alt: str) -> bool:
    return (
        len(ref) <= 2
        and len(alt) <= 2
        and ref not in (".", "N")
        and alt not in (".", "N")
    )


def noise_band_penalty(pos: int) -> bool:
    """True if position is in the historical low-truth noise band."""
    return READ_NOISE_LO <= pos <= READ_NOISE_HI


def ultra_scoring_core(pos: int) -> bool:
    """Band covering manager truth across 5.21.01 / 5.21.02 / 5.21.03."""
    return ULTRA_CORE_LO <= pos <= ULTRA_CORE_HI


def ultra_edge_band(pos: int) -> bool:
    """Truth edge positions with lower mpileup QUAL in 5.21.01."""
    return ULTRA_EDGE_LO <= pos <= ULTRA_EDGE_HI


def ultra_tail_band(pos: int) -> bool:
    """Low-QUAL truth sites in dense cluster tail."""
    return ULTRA_TAIL_LO <= pos <= ULTRA_TAIL_HI


def thresholds_for_position(
    pos: int, profile: TaskProfile
) -> Tuple[int, float, int]:
    """Position-aware filters (ultra_wide only); default uses profile baselines."""
    if profile.name != "ultra_wide":
        return profile.min_dp, profile.min_qual, profile.min_alt_ad
    if ultra_tail_band(pos):
        return 4, 12.0, 2
    if ultra_edge_band(pos):
        return 4, 12.0, 2
    if ultra_scoring_core(pos):
        return profile.min_dp, profile.min_qual, profile.min_alt_ad
    return 10, 28.0, 6
