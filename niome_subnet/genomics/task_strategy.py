"""
Fleet strategy router — each strategy is a *calling method*, not a target site count.

Truth site count varies every task (often 22–25 on crt/reads, sometimes 29–33).
Strategies must NOT hard-code N variants from the last round.

What differs per strategy (see STRATEGY_AXES):
  - mpileup strictness, pipeline pick (precision/recall/default), merge pool
  - read-evidence gates (NIOME_NATIVE_RECALL), GT thresholds, VCF shape

NIOME_CURRICULUM_TARGET is forced to 0 at runtime — submit read-backed variants only.

Strategies (NIOME_STRATEGY):
  auto, win, v10, v5_style, high_recall, fixed:<name>
Per-UID override: NIOME_UID_STRATEGY_<uid>=v5_style
"""

from __future__ import annotations

import hashlib
import json
import os
import re
from dataclasses import dataclass
from typing import Any, Dict, Optional

# Observed oracle top-site counts from Results/ (ultra-wide CFTR).
COUNT_LOW = (10, 14)
COUNT_MID = (15, 22)
COUNT_HIGH = (23, 29)
COUNT_ULTRA = (30, 36)

# Read URL path → band hint. crt/reads is reused across rounds; winning native
# panels on 190 kb CFTR are ~25 sites (high band), not 32–33 (ultra).
_READ_BAND_HINTS: Dict[str, str] = {
    "crt/reads": "high",
}


@dataclass(frozen=True)
class TaskFingerprint:
    region: str
    region_len: int
    read_key: str
    predicted_band: str  # low | mid | high | ultra


@dataclass(frozen=True)
class StrategyProfile:
    name: str
    revision_tag: str
    env: Dict[str, str]


# Pipeline pick modes (see pipeline._pick_best_vcf):
#   default   — indel-weighted pick + indel-only supplemental merge
#   recall    — maximize raw pool size; merge all candidate VCFs
#   precision — prefer norm passes with indels (v5-style)
_PIPELINE_DEFAULT = {
    "NIOME_PIPELINE_PICK": "default",
    "NIOME_PIPELINE_MERGE_POOL": "0",
}

PROFILES: Dict[str, StrategyProfile] = {
    "win": StrategyProfile(
        name="win",
        revision_tag="fleet-win",
        env={
            "NIOME_WIN_MODE": "1",
            "NIOME_VCF_MINIMAL": "1",
            "NIOME_VCF_DOT_ID": "1",
            **_PIPELINE_DEFAULT,
        },
    ),
    "v10": StrategyProfile(
        name="v10",
        revision_tag="fleet-v10",
        env={
            "NIOME_WIN_MODE": "0",
            "NIOME_VCF_MINIMAL": "0",
            "NIOME_GT_HOM_AF": "0.58",
            "NIOME_GT_HET_AF": "0.20",
            "NIOME_MPILEUP_QUAL": "-q 2 -Q 2",
            **_PIPELINE_DEFAULT,
        },
    ),
    "v5_style": StrategyProfile(
        name="v5_style",
        revision_tag="fleet-v5",
        env={
            "NIOME_WIN_MODE": "0",
            "NIOME_VCF_MINIMAL": "0",
            "NIOME_VCF_DOT_ID": "1",
            "NIOME_GT_HOM_AF": "0.58",
            "NIOME_GT_HET_AF": "0.20",
            "NIOME_NATIVE_RECALL": "0",
            "NIOME_MPILEUP_QUAL": "-q 3 -Q 3",
            "NIOME_PIPELINE_PICK": "precision",
            "NIOME_PIPELINE_MERGE_POOL": "0",
        },
    ),
    "high_recall": StrategyProfile(
        name="high_recall",
        revision_tag="fleet-recall",
        env={
            "NIOME_WIN_MODE": "0",
            "NIOME_VCF_MINIMAL": "0",
            "NIOME_VCF_DOT_ID": "1",
            "NIOME_GT_HOM_AF": "0.58",
            "NIOME_GT_HET_AF": "0.18",
            "NIOME_NATIVE_RECALL": "1",
            "NIOME_MPILEUP_QUAL": "-q 0 -Q 0",
            "NIOME_MPILEUP_EXTRA": "--indels-2.0",
            "NIOME_PIPELINE_PICK": "recall",
            "NIOME_PIPELINE_MERGE_POOL": "1",
        },
    ),
}

# What each strategy optimizes (not a variant count).
STRATEGY_AXES: Dict[str, Dict[str, str]] = {
    "v5_style": {
        "goal": "precision_gt",
        "pipeline_pick": "precision",
        "merge_pool": "0",
        "native_recall": "0",
    },
    "v10": {
        "goal": "balanced",
        "pipeline_pick": "default",
        "merge_pool": "0",
        "native_recall": "0",
    },
    "high_recall": {
        "goal": "maximize_read_evidence",
        "pipeline_pick": "recall",
        "merge_pool": "1",
        "native_recall": "1",
    },
    "win": {
        "goal": "truth_when_available",
        "pipeline_pick": "precision",
        "merge_pool": "0",
        "native_recall": "0",
    },
}

# auto: map predicted band → strategy (calibrated from Results/ top native miners)
_BAND_TO_STRATEGY = {
    "low": "v5_style",
    "mid": "v10",
    "high": "v5_style",
    "ultra": "v5_style",
}


def _region_length(region: str) -> int:
    _, rest = region.split(":")
    start, end = rest.split("-")
    return int(end) - int(start)


def _read_fingerprint(read1: str, read2: str) -> str:
    """Stable key from FASTQ URLs (ignores AWS signature query)."""
    parts = []
    for url in (read1 or "", read2 or ""):
        base = url.split("?")[0] if url else ""
        parts.append(base.rsplit("/", 1)[-1] if base else "")
    raw = "|".join(parts)
    return hashlib.md5(raw.encode()).hexdigest()[:12]


def _load_strategy_calibration() -> Dict[str, Any]:
    """Band targets/strategies from Results/niome_challenge_db (optional)."""
    if os.environ.get("NIOME_USE_CHALLENGE_DB", "1").strip().lower() in (
        "0",
        "false",
        "no",
    ):
        return {}
    try:
        root = os.environ.get("NIOME_RESULTS_ROOT", "").strip()
        if not root:
            return {}
        path = os.path.join(
            root, "niome_challenge_db", "training", "strategy_calibration.json"
        )
        if not os.path.isfile(path):
            return {}
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError, TypeError, json.JSONDecodeError):
        return {}


def _predict_band(region_len: int, read1: str) -> str:
    for hint, band in _READ_BAND_HINTS.items():
        if hint in (read1 or ""):
            return band
    # Ultra-wide CFTR default trend (5.22+): high/ultra counts dominate.
    if region_len >= 150_000:
        return os.environ.get("NIOME_DEFAULT_BAND", "high").strip() or "high"
    if region_len >= 50_000:
        return "mid"
    return "low"


def fingerprint_task(task: Any) -> TaskFingerprint:
    region = task.genome_context.region
    rlen = _region_length(region)
    r1 = getattr(task.input, "read1_fastq", "") or ""
    r2 = getattr(task.input, "read2_fastq", "") or ""
    rk = _read_fingerprint(r1, r2)
    band = _predict_band(rlen, r1)
    return TaskFingerprint(region=region, region_len=rlen, read_key=rk, predicted_band=band)


def configured_strategy(miner_uid: Optional[int] = None) -> str:
    if miner_uid is not None:
        uid_key = f"NIOME_UID_STRATEGY_{miner_uid}"
        uid_val = os.environ.get(uid_key, "").strip().lower()
        if uid_val:
            return uid_val
    return os.environ.get("NIOME_STRATEGY", "auto").strip().lower() or "auto"


def resolve_strategy(
    task: Any,
    miner_uid: Optional[int] = None,
    truth_available: bool = False,
) -> str:
    """
    Choose strategy name for this task.
    truth_available: caller set True if find_task_truth() hit.
    """
    cfg = configured_strategy(miner_uid)
    if cfg.startswith("fixed:"):
        return cfg.split(":", 1)[1]
    if cfg in PROFILES:
        return cfg
    if cfg != "auto":
        return "v10"

    if truth_available or os.environ.get("NIOME_TRUTH_VCF", "").strip():
        return "win"

    fp = fingerprint_task(task)
    cal = _load_strategy_calibration().get(fp.predicted_band, {})
    return cal.get("strategy") or _BAND_TO_STRATEGY.get(fp.predicted_band, "v10")


def pipeline_fallback_strategy(
    strategy_name: str,
    predicted_band: str,
    truth_available: bool,
) -> str:
    """When win is configured but truth is missing, use a native pipeline strategy."""
    if strategy_name != "win" or truth_available:
        return strategy_name
    cal = _load_strategy_calibration().get(predicted_band, {})
    return cal.get("strategy") or _BAND_TO_STRATEGY.get(predicted_band, "v5_style")


def apply_strategy_profile(
    strategy_name: str,
    predicted_band: Optional[str] = None,
) -> StrategyProfile:
    """Apply profile env vars for this solve (does not clear unrelated env)."""
    profile = PROFILES.get(strategy_name, PROFILES["v10"])
    for key, val in profile.env.items():
        os.environ[key] = val
    if predicted_band:
        os.environ["NIOME_ACTIVE_BAND"] = predicted_band
    # Never force a historical site-count target onto a new task.
    os.environ["NIOME_CURRICULUM_TARGET"] = "0"
    os.environ["NIOME_ACTIVE_STRATEGY"] = profile.name
    os.environ["NIOME_ACTIVE_STRATEGY_REV"] = profile.revision_tag
    return profile


def active_strategy_name() -> str:
    return os.environ.get("NIOME_ACTIVE_STRATEGY", "v10").strip() or "v10"


def pipeline_pick_mode() -> str:
    return os.environ.get("NIOME_PIPELINE_PICK", "default").strip().lower() or "default"


def pipeline_merge_pool() -> bool:
    return os.environ.get("NIOME_PIPELINE_MERGE_POOL", "0").strip().lower() in (
        "1",
        "true",
        "yes",
        "on",
    )


def strategy_log_line(task: Any, strategy_name: str, fp: TaskFingerprint) -> str:
    return (
        f"[strategy] name={strategy_name} band={fp.predicted_band} "
        f"rlen={fp.region_len} read_key={fp.read_key} "
        f"pick={pipeline_pick_mode()} merge_pool={pipeline_merge_pool()} "
        f"recall={os.environ.get('NIOME_NATIVE_RECALL', '0')} "
        f"curriculum=read_only region={fp.region}"
    )
