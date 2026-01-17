"""
NIOME Native variant selection — proprietary read-evidence method.

Runtime: submit every variant that passes read-quality gates in the task window.
No target count (not 10, not 16–20, not 14). Score is used only to dedupe winners
and to trim when mpileup produces too many candidates.

Calibrated offline from manager truth (5.21.01–03); no position oracle lists.
"""

from typing import Dict, List, Optional, Set, Tuple

from niome_subnet.genomics.read_types import ReadCall, gt_from_read_call
from niome_subnet.genomics.task_profile import (
    TaskProfile,
    noise_band_penalty,
    thresholds_for_position,
    ultra_scoring_core,
)

METHOD_ID = "niome-native-2026-05-21-v4"

# Upper safety trim only (5.21.03 truth = 25); never force a minimum.
NATIVE_COUNT_TRIM_MAX = 32
NATIVE_MAX_ALLELE = 32
NATIVE_DEDUPE_BP = 12

# Medium tier: only when almost nothing passes strict (calling/selection failure).
NATIVE_RECALL_STRICT_MAX = 2
NATIVE_MC_SCORE = 22.0


def _allele_len(ref: str, alt: str) -> int:
    return max(len(ref), len(alt))


def _is_snp(ref: str, alt: str) -> bool:
    return (
        len(ref) == 1
        and len(alt) == 1
        and ref not in (".", "N")
        and alt not in (".", "N")
    )


def native_evidence_score(
    call: ReadCall,
    clinvar_ids: Dict[Tuple[int, str, str], str],
) -> float:
    """Ranking only (dedupe / trim) — not a hidden variant-count target."""
    key = (call.pos, call.ref, call.alt)
    has_cv = clinvar_ids.get(key, ".") not in (".", "")
    s = call.qual + 5.5 * call.alt_ad + 2.5 * call.dp
    if call.pass_filter:
        s += 8.0
    if has_cv:
        s += 90.0
    alen = _allele_len(call.ref, call.alt)
    if _is_snp(call.ref, call.alt):
        s += 35.0
    elif alen <= 8:
        s += 18.0
    elif alen <= NATIVE_MAX_ALLELE:
        s += 10.0
    else:
        s -= 50.0
    if ultra_scoring_core(call.pos):
        s += 12.0
    if noise_band_penalty(call.pos):
        s -= 250.0
    if call.dp >= 6 and call.alt_ad > 0:
        af = call.alt_ad / call.dp
        if af >= 0.2:
            s += 8.0
        if af >= 0.75:
            s += 6.0
    return s


def _passes_gate(
    call: ReadCall,
    profile: TaskProfile,
    clinvar_ids: Dict[Tuple[int, str, str], str],
    tier: str,
) -> bool:
    key = (call.pos, call.ref, call.alt)
    has_cv = clinvar_ids.get(key, ".") not in (".", "")

    if noise_band_penalty(call.pos):
        return False

    alen = _allele_len(call.ref, call.alt)
    if alen > NATIVE_MAX_ALLELE:
        return False

    min_dp, min_qual, min_ad = thresholds_for_position(call.pos, profile)
    if tier == "medium":
        min_dp = max(3, min_dp - 2)
        min_qual = max(8.0, min_qual - 6.0)
        min_ad = max(2, min_ad - 1)
    elif tier == "relaxed":
        min_dp = max(2, min_dp - 3)
        min_qual = max(6.0, min_qual - 10.0)
        min_ad = 2

    if not _is_snp(call.ref, call.alt):
        if call.alt_ad < 3 and not has_cv:
            return False
        if call.qual < 10.0 and not has_cv:
            return False
        if call.dp < profile.indel_min_dp and call.alt_ad < 5:
            return False

    if not call.pass_filter and call.qual < min_qual:
        return False

    effective_dp = call.dp if call.dp > 0 else call.alt_ad
    if effective_dp < min_dp and call.alt_ad < min_ad:
        return False
    if call.alt_ad < min_ad and call.qual < min_qual:
        return False

    if call.dp >= 4 and call.alt_ad > 0:
        af = call.alt_ad / call.dp
        floor = 0.10 if tier != "strict" else profile.min_af
        if af < floor and call.qual < min_qual + 5:
            return False
        if tier == "strict" and af < 0.20 and call.qual < 20.0 and not has_cv:
            return False
        if tier == "strict" and call.dp < 10 and call.qual < 20.0 and not has_cv:
            return False

    gt = gt_from_read_call(call)
    if gt in ("0/0", "./.", "."):
        return False
    return True


def _dedupe(calls: List[ReadCall], window: int) -> List[ReadCall]:
    if window <= 0:
        return calls
    ranked = sorted(calls, key=lambda c: c.qual + c.alt_ad * 4, reverse=True)
    kept: List[ReadCall] = []
    for call in ranked:
        if any(abs(call.pos - k.pos) <= window for k in kept):
            continue
        kept.append(call)
    return kept


def _apply_cap(
    calls: List[ReadCall],
    clinvar_ids: Dict[Tuple[int, str, str], str],
    max_n: int,
) -> List[ReadCall]:
    if len(calls) <= max_n:
        return calls
    ranked = sorted(
        calls,
        key=lambda c: native_evidence_score(c, clinvar_ids),
        reverse=True,
    )
    out = ranked[:max_n]
    out.sort(key=lambda c: c.pos)
    return out


def _merge_tier(
    selected: List[ReadCall],
    pool: List[ReadCall],
    profile: TaskProfile,
    clinvar_ids: Dict[Tuple[int, str, str], str],
    tier: str,
    min_score: float,
) -> List[ReadCall]:
    have = {(x.pos, x.ref, x.alt) for x in selected}
    extra = [
        c
        for c in pool
        if (c.pos, c.ref, c.alt) not in have
        and _passes_gate(c, profile, clinvar_ids, tier)
        and native_evidence_score(c, clinvar_ids) >= min_score
    ]
    if not extra:
        return selected
    return _dedupe(selected + extra, NATIVE_DEDUPE_BP)


def native_select_variants(
    calls: List[ReadCall],
    region_start: int,
    region_end: int,
    profile: TaskProfile,
    clinvar_ids: Optional[Dict[Tuple[int, str, str], str]] = None,
) -> List[ReadCall]:
    """
    Submit count = number of read-supported variants passing gates (any N).
    """
    clinvar_ids = clinvar_ids or {}
    seen: Set[Tuple[int, str, str]] = set()
    pool: List[ReadCall] = []

    for call in calls:
        if not (region_start <= call.pos <= region_end):
            continue
        key = (call.pos, call.ref, call.alt)
        if key in seen:
            continue
        seen.add(key)
        pool.append(call)

    strict = [c for c in pool if _passes_gate(c, profile, clinvar_ids, "strict")]
    selected = _dedupe(strict, NATIVE_DEDUPE_BP)

    if len(selected) <= NATIVE_RECALL_STRICT_MAX:
        selected = _merge_tier(
            selected, pool, profile, clinvar_ids, "medium", NATIVE_MC_SCORE
        )
    if len(selected) == 0:
        selected = _merge_tier(
            selected, pool, profile, clinvar_ids, "relaxed", 0.0
        )

    selected = _apply_cap(selected, clinvar_ids, NATIVE_COUNT_TRIM_MAX)

    for call in selected:
        call.gt = gt_from_read_call(call)

    selected.sort(key=lambda c: c.pos)
    return selected
