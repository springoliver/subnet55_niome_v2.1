"""Challenge DB builder (no bittensor)."""
import json
from pathlib import Path

import sys
import types

ROOT = Path(__file__).resolve().parents[1]
ns_pkg = types.ModuleType("niome_subnet")
ns_pkg.__path__ = [str(ROOT / "niome_subnet")]
ana_pkg = types.ModuleType("niome_subnet.analysis")
ana_pkg.__path__ = [str(ROOT / "niome_subnet" / "analysis")]
sys.modules["niome_subnet"] = ns_pkg
sys.modules["niome_subnet.analysis"] = ana_pkg

from niome_subnet.analysis.challenge_db import build_database, discover_round_dirs
from niome_subnet.analysis.vcf_panel import band_from_count, compare_panels, parse_vcf_text


def test_band_mapping():
    assert band_from_count(12) == "low"
    assert band_from_count(20) == "mid"
    assert band_from_count(27) == "high"
    assert band_from_count(32) == "ultra"


def test_compare_panels():
    truth = parse_vcf_text(
        "chr7\t117480010\t.\tC\tG\t.\tPASS\t.\tGT\t0/1\n", source="truth"
    )
    miner = parse_vcf_text(
        "chr7\t117480010\t.\tC\tG\t.\tPASS\t.\tGT\t1/1\n", source="native"
    )
    c = compare_panels(truth, miner)
    assert c["tp"] == 1
    assert c["gt_match"] == 0


def test_build_db_on_real_results():
    root = Path(__file__).resolve().parents[1] / "Results"
    if not (root / "5.21.01" / "miner.json").is_file():
        return
    rounds = discover_round_dirs(root)
    assert len(rounds) >= 10
    db = build_database(root)
    manifest = json.loads((db / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["n_rounds"] >= 10
    assert manifest["n_with_truth"] >= 1
