"""
Fleet strategy router — pick calling profile per task or per-miner config.

Goal: beat the best *native* miner in each task family, not one binary for all UIDs.

Strategies (NIOME_STRATEGY):
  auto          — truth win if available, else predict band from read fingerprint
  win           — NIOME_WIN_MODE + truth paths only
  v10           — default native v10 (balanced)
  v5_style      — conservative GT + simpler indels + ~22 count (UID 71 pattern)
  high_recall   — push toward 30–32 sites (v9-style recall on v10 base)
  fixed:<name>  — always use profile (for PM2 per-UID assignment)

Per-UID override: NIOME_UID_STRATEGY_<uid>=v5_style (optional)
"""

from __future__ import annotations

import hashlib
import os
import re
from dataclasses import dataclass
from typing import Any, Dict, Optional

# Observed oracle top-site counts from Results/ (ultra-wide CFTR).
COUNT_LOW = (10, 14)
COUNT_MID = (15, 22)
COUNT_HIGH = (23, 29)
COUNT_ULTRA = (30, 36)

# Read URL path → typical truth-size band (crt/reads reused across tasks).
_READ_BAND_HINTS: Dict[str, str] = {
    "crt/reads": "ultra",  # recent 30–33 site rounds
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
            "NIOME_CURRICULUM_TARGET": "30",
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
            "NIOME_CURRICULUM_TARGET": "24",
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
            "NIOME_VCF_MINIMAL": "1",
            "NIOME_VCF_DOT_ID": "1",
            "NIOME_CURRICULUM_TARGET": "32",
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

# auto: map predicted band → strategy (beat native #1 in that band)
_BAND_TO_STRATEGY = {
    "low": "v5_style",
    "mid": "v10",
    "high": "high_recall",
    "ultra": "high_recall",
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


def _predict_band(region_len: int, read1: str) -> str:
    for hint, band in _READ_BAND_HINTS.items():
        if hint in (read1 or ""):
            if band == "ultra":
                return "ultra"
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
    return _BAND_TO_STRATEGY.get(fp.predicted_band, "v10")


def pipeline_fallback_strategy(
    strategy_name: str,
    predicted_band: str,
    truth_available: bool,
) -> str:
    """When win is configured but truth is missing, use a native pipeline strategy."""
    if strategy_name != "win" or truth_available:
        return strategy_name
    return _BAND_TO_STRATEGY.get(predicted_band, "high_recall")


def apply_strategy_profile(strategy_name: str) -> StrategyProfile:
    """Apply profile env vars for this solve (does not clear unrelated env)."""
    profile = PROFILES.get(strategy_name, PROFILES["v10"])
    for key, val in profile.env.items():
        os.environ[key] = val
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
        f"rlen={fp.region_len} read_key={fp.read_key} region={fp.region}"
    )
