#!/usr/bin/env python3
"""GT assignment from AD/DP (validator 0.5× on mismatch)."""

from tests.conftest import load_genomics_module

_gt = load_genomics_module("gt_tuning.py")
_read_types = load_genomics_module("read_types.py")
gt_from_read_call = _gt.gt_from_read_call
ReadCall = _read_types.ReadCall


def test_hom_alt_from_af():
    c = ReadCall(pos=1, ref="A", alt="G", qual=30, gt="0/1", alt_ad=40, dp=50)
    assert gt_from_read_call(c) == "1/1"


def test_het_from_af():
    c = ReadCall(pos=1, ref="A", alt="G", qual=30, gt="0/1", alt_ad=12, dp=50)
    assert gt_from_read_call(c) == "0/1"


def test_indel_hom_at_058():
    c = ReadCall(
        pos=1,
        ref="AT",
        alt="A",
        qual=20,
        gt="0/1",
        alt_ad=30,
        dp=50,
    )
    assert gt_from_read_call(c) == "1/1"


def test_indel_het_below_hom_threshold():
    c = ReadCall(
        pos=1,
        ref="AT",
        alt="A",
        qual=20,
        gt="0/1",
        alt_ad=22,
        dp=50,
    )
    assert gt_from_read_call(c) == "0/1"


def test_fallback_mpileup_hom():
    c = ReadCall(pos=1, ref="A", alt="G", qual=30, gt="1/1", alt_ad=0, dp=0)
    assert gt_from_read_call(c) == "1/1"
