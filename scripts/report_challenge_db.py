#!/usr/bin/env python3
"""Print training calibration summary from niome_challenge_db."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _bootstrap():
    import types

    ns_pkg = types.ModuleType("niome_subnet")
    ns_pkg.__path__ = [str(ROOT / "niome_subnet")]
    ana_pkg = types.ModuleType("niome_subnet.analysis")
    ana_pkg.__path__ = [str(ROOT / "niome_subnet" / "analysis")]
    sys.modules["niome_subnet"] = ns_pkg
    sys.modules["niome_subnet.analysis"] = ana_pkg
    from niome_subnet.analysis.challenge_db import load_manifest

    return load_manifest


load_manifest = _bootstrap()

RESULTS = ROOT / "Results"
TRAINING = RESULTS / "niome_challenge_db" / "training"


def main() -> None:
    manifest = load_manifest(RESULTS)
    cal = json.loads((TRAINING / "strategy_calibration.json").read_text(encoding="utf-8"))
    missed = json.loads(
        (TRAINING / "commonly_missed_positions.json").read_text(encoding="utf-8")
    )
    fleet = json.loads(
        (TRAINING / "fleet_score_by_strategy.json").read_text(encoding="utf-8")
    )

    print("=" * 72)
    print("NIOME CHALLENGE DATABASE — full history summary")
    print("=" * 72)
    print(f"Rounds: {manifest['n_rounds']}  |  With manager truth: {manifest['n_with_truth']}")
    print(f"Built: {manifest['built_at']}")
    print("\n--- Strategy calibration by band (all history) ---")
    for band, info in cal.items():
        tc = info.get("truth_site_counts") or {}
        top = info.get("top_site_counts") or {}
        print(
            f"  {band:5s}  strategy={info['strategy']:12s}  "
            f"curriculum_target={info.get('curriculum_target_suggested')}  "
            f"truth_n median={tc.get('median', '?')}  top_n median={top.get('median', '?')}"
        )
    print("\n--- Fleet mean score by configured strategy ---")
    for strat, st in sorted(fleet.items(), key=lambda x: -x[1].get("mean", 0)):
        print(f"  {strat:12s}  n={st['n']:3d}  mean={st['mean']:.4f}  max={st['max']:.4f}")
    print("\n--- Top 15 commonly missed truth positions (fleet) ---")
    for row in missed[:15]:
        print(f"  chr7:{row['pos']}  missed_in_{row['count']}_rounds")


if __name__ == "__main__":
    main()
