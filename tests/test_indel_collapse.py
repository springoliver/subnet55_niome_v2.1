#!/usr/bin/env python3
from tests.conftest import load_genomics_module

_rc = load_genomics_module("read_types.py")
_rd = load_genomics_module("read_calling.py")
ReadCall = _rc.ReadCall
merge_read_call_pools = _rd.merge_read_call_pools
collapse_indels_at_position = _rd.collapse_indels_at_position


def test_collapse_prefers_simpler_indel():
    simple = ReadCall(
        pos=117550866,
        ref="GA",
        alt="G",
        qual=100,
        gt="0/1",
        alt_ad=10,
        dp=20,
    )
    complex_ = ReadCall(
        pos=117550866,
        ref="GAAA",
        alt="GAA",
        qual=120,
        gt="0/1",
        alt_ad=12,
        dp=22,
    )
    out = collapse_indels_at_position([simple, complex_])
    assert len(out) == 1
    assert out[0].ref == "GA"


def test_merge_pools_collapses_positions():
    pools = [
        [
            ReadCall(117550866, "GA", "G", 100, "0/1", 10, 20),
            ReadCall(117550866, "GAAA", "GAA", 120, "0/1", 12, 22),
        ]
    ]
    merged = merge_read_call_pools(pools)
    assert len(merged) == 1
