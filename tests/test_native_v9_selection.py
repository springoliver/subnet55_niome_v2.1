#!/usr/bin/env python3
from tests.conftest import load_genomics_module

load_genomics_module("gt_tuning.py")
_rt = load_genomics_module("read_types.py")
_tp = load_genomics_module("task_profile.py")
_ev = load_genomics_module("evidence_selection.py")

METHOD_ID = _ev.METHOD_ID
emergency_select_variants = _ev.emergency_select_variants
native_select_variants = _ev.native_select_variants
PROFILES = _tp.PROFILES
ReadCall = _rt.ReadCall


def test_v9_method_id():
    assert "v9" in METHOD_ID


def test_emergency_select_from_pool():
    profile = PROFILES["ultra_wide"]
    pool = [
        ReadCall(
            pos=117500000 + i * 25,
            ref="A",
            alt="G",
            qual=15.0,
            gt="0/1",
            alt_ad=3,
            dp=10,
        )
        for i in range(5)
    ]
    out = emergency_select_variants(
        pool, 117480000, 117670000, profile, {}, max_n=5
    )
    assert len(out) == 5


def test_indel_curriculum_fill():
    profile = PROFILES["ultra_wide"]
    snp = ReadCall(
        pos=117500100,
        ref="A",
        alt="G",
        qual=40,
        gt="0/1",
        alt_ad=10,
        dp=20,
        pass_filter=True,
    )
    indel = ReadCall(
        pos=117500200,
        ref="AT",
        alt="A",
        qual=18,
        gt="0/1",
        alt_ad=4,
        dp=12,
        pass_filter=True,
    )
    pool = [snp] * 15 + [indel] * 8
    selected = native_select_variants(
        pool, 117480000, 117670000, profile, {}
    )
    assert len(selected) >= 2
    assert any(len(c.ref) != len(c.alt) for c in selected)
