#!/usr/bin/env python3
import os

from tests.conftest import load_genomics_module

_ts = load_genomics_module("task_strategy.py")
PROFILES = _ts.PROFILES
apply_strategy_profile = _ts.apply_strategy_profile
fingerprint_task = _ts.fingerprint_task
resolve_strategy = _ts.resolve_strategy


class _Task:
    class genome_context:
        region = "chr7:117480000-117670000"

    class input:
        read1_fastq = "https://bucket/crt/reads_1.fq?sig=1"
        read2_fastq = "https://bucket/crt/reads_2.fq?sig=2"


def test_fingerprint_ultra_band():
    fp = fingerprint_task(_Task())
    assert fp.region_len == 190_000
    assert fp.predicted_band in ("high", "ultra")


def test_resolve_auto_with_truth():
    os.environ["NIOME_STRATEGY"] = "auto"
    os.environ.pop("NIOME_TRUTH_VCF", None)
    s = resolve_strategy(_Task(), truth_available=True)
    assert s == "win"


def test_uid_override():
    os.environ["NIOME_STRATEGY"] = "v10"
    os.environ["NIOME_UID_STRATEGY_71"] = "v5_style"
    s = resolve_strategy(_Task(), miner_uid=71)
    assert s == "v5_style"


def test_apply_profile_sets_env():
    p = apply_strategy_profile("high_recall")
    assert os.environ.get("NIOME_NATIVE_RECALL") == "1"
    assert p.name == "high_recall"
