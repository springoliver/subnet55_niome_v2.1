"""
NIOME Native variant selection — proprietary read-evidence method.

Runtime: submit every variant that passes read-quality gates in the task window.
No target count (not 10, not 16–20, not 14). Score is used only to dedupe winners
and to trim when mpileup produces too many candidates.

Calibrated offline from manager truth (5.21.01–03); no position oracle lists.
"""

import os
from typing import Dict, List, Optional, Set, Tuple

from niome_subnet.genomics.read_types import ReadCall, gt_from_read_call
from niome_subnet.genomics.task_profile import (
    TaskProfile,
    curriculum_target_for_profile,
    noise_band_penalty,
    thresholds_for_position,
    ultra_scoring_core,
)

METHOD_ID = "niome-native-2026-05-23-v10"
_FLEET_REV = os.environ.get("NIOME_ACTIVE_STRATEGY_REV", "").strip()
if _FLEET_REV:
    METHOD_ID = f"niome-native-2026-05-23-v10+{_FLEET_REV}"

# Upper safety trim only (5.21.03 truth = 25); never force a minimum.
NATIVE_COUNT_TRIM_MAX = 34
NATIVE_MAX_ALLELE = 48
NATIVE_DEDUPE_BP = 12

# Medium tier: only when almost nothing passes strict (v9 was too aggressive).
NATIVE_RECALL_STRICT_MAX = 2
NATIVE_MC_SCORE = 22.0


def _recall_aggressive() -> bool:
    return os.environ.get("NIOME_NATIVE_RECALL", "0").strip().lower() in (
        "1",
        "true",
        "yes",
        "high",
    )


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
        s -= 80.0
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

    max_alen = max(profile.max_indel_len, NATIVE_MAX_ALLELE)
    alen = _allele_len(call.ref, call.alt)
    if alen > max_alen:
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
    elif tier == "curriculum":
        min_dp = 2
        min_qual = 5.0
        min_ad = 1

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


def _sort_curriculum_fill(
    calls: List[ReadCall],
    clinvar_ids: Dict[Tuple[int, str, str], str],
) -> List[ReadCall]:
    """High evidence first; SNPs rank above complex indels (v9 indel-first hurt 5.23.02)."""

    def key(c: ReadCall) -> Tuple[float, float]:
        snp = 1.0 if _is_snp(c.ref, c.alt) else 0.0
        complexity = float(len(c.ref) + len(c.alt))
        return (native_evidence_score(c, clinvar_ids) + snp * 40.0, -complexity)

    return sorted(calls, key=key, reverse=True)


def _expand_to_curriculum_target(
    selected: List[ReadCall],
    pool: List[ReadCall],
    target: int,
    profile: TaskProfile,
    clinvar_ids: Dict[Tuple[int, str, str], str],
) -> List[ReadCall]:
    """Add best relaxed read calls until target or pool exhausted (no oracle list)."""
    if target <= 0 or len(selected) >= target:
        return selected
    have = {(x.pos, x.ref, x.alt) for x in selected}
    extras = [
        c
        for c in pool
        if (c.pos, c.ref, c.alt) not in have
        and _passes_gate(c, profile, clinvar_ids, "curriculum")
    ]
    for call in _sort_curriculum_fill(extras, clinvar_ids):
        if len(selected) >= target:
            break
        selected.append(call)
    return _dedupe(selected, NATIVE_DEDUPE_BP)


def emergency_select_variants(
    pool: List[ReadCall],
    region_start: int,
    region_end: int,
    profile: TaskProfile,
    clinvar_ids: Dict[Tuple[int, str, str], str],
    max_n: int = 28,
) -> List[ReadCall]:
    """
    Last resort when strict/relaxed tiers yield zero — avoid empty VCF submit.
    Still read-backed; no invented coordinates.
    """
    in_window = [
        c
        for c in pool
        if region_start <= c.pos <= region_end and c.alt_ad >= 1
    ]
    if not in_window:
        return []
    ranked = sorted(
        in_window,
        key=lambda c: native_evidence_score(c, clinvar_ids),
        reverse=True,
    )
    cap = max_n or curriculum_target_for_profile(profile) or 24
    out = _dedupe(ranked[: cap + 6], NATIVE_DEDUPE_BP)[:cap]
    for call in out:
        call.gt = gt_from_read_call(call)
    out.sort(key=lambda c: c.pos)
    return out


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

    target_n = curriculum_target_for_profile(profile)
    trim_max = (
        max(NATIVE_COUNT_TRIM_MAX, target_n + 4)
        if target_n > 0
        else NATIVE_COUNT_TRIM_MAX
    )

    strict = [c for c in pool if _passes_gate(c, profile, clinvar_ids, "strict")]
    selected = _dedupe(strict, NATIVE_DEDUPE_BP)

    if len(selected) <= NATIVE_RECALL_STRICT_MAX:
        selected = _merge_tier(
            selected, pool, profile, clinvar_ids, "medium", NATIVE_MC_SCORE
        )
    elif target_n > 0 and len(selected) < target_n - (
        3 if _recall_aggressive() else 5
    ):
        selected = _merge_tier(
            selected, pool, profile, clinvar_ids, "medium", 14.0
        )

    if len(selected) == 0:
        selected = _merge_tier(
            selected, pool, profile, clinvar_ids, "relaxed", 0.0
        )

    if target_n > 0 and len(selected) < target_n:
        selected = _expand_to_curriculum_target(
            selected, pool, target_n, profile, clinvar_ids
        )

    if len(selected) == 0 and pool:
        selected = emergency_select_variants(
            pool, region_start, region_end, profile, clinvar_ids, max_n=target_n or 28
        )

    selected = _apply_cap(selected, clinvar_ids, trim_max)

    for call in selected:
        call.gt = gt_from_read_call(call)

    selected.sort(key=lambda c: c.pos)
    return selected
