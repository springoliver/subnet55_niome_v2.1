#!/usr/bin/env python3
"""
Review fleet vs native/oracle using full challenge DB (all Results/ history).

Run:
  python scripts/build_challenge_db.py
  python tests/analyze_fleet_strategy.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "Results"


def _bootstrap():
    import types

    ns_pkg = types.ModuleType("niome_subnet")
    ns_pkg.__path__ = [str(ROOT / "niome_subnet")]
    ana_pkg = types.ModuleType("niome_subnet.analysis")
    ana_pkg.__path__ = [str(ROOT / "niome_subnet" / "analysis")]
    sys.modules["niome_subnet"] = ns_pkg
    sys.modules["niome_subnet.analysis"] = ana_pkg
USER_UIDS = {139, 71, 50, 99, 209, 235, 226, 124, 9, 97, 217, 36, 92, 225, 38, 155}


def _load_rounds_from_db():
    _bootstrap()
    from niome_subnet.analysis.challenge_db import load_manifest

    manifest = load_manifest(RESULTS)
    rounds = []
    for meta in manifest["rounds"]:
        path = RESULTS / "niome_challenge_db" / "rounds" / f"{meta['key']}.json"
        rounds.append(json.loads(path.read_text(encoding="utf-8")))
    return rounds


def main():
    try:
        rounds = _load_rounds_from_db()
        print("=" * 72)
        print(f"NIOME fleet review — FULL DATABASE ({len(rounds)} rounds)")
        print("=" * 72)
    except FileNotFoundError:
        print("Challenge DB missing. Run: python scripts/build_challenge_db.py")
        sys.exit(1)

    wins = 0
    for rd in rounds:
        band = rd["band"]
        rec = rd["recommended_strategy"]
        top = rd["top_miner"]
        nat = rd.get("native_best")
        nat_score = nat["final_score"] if nat else 0.0
        nat_n = nat["n_sites"] if nat else 0

        your_best = 0.0
        your_uid = None
        for uid, recf in rd.get("fleet", {}).items():
            if int(uid) not in USER_UIDS:
                continue
            s = recf["final_score"]
            if s > your_best:
                your_best = s
                your_uid = int(uid)

        beat = your_best >= nat_score - 0.001 if nat else False
        if beat:
            wins += 1
        gap = nat_score - your_best

        truth_s = f"truth_n={rd['truth_n']}" if rd["has_truth"] else "no truth"
        print(
            f"\n{rd['round_path']:32s}  task={rd['task_id'][:8]}  "
            f"band={band}  -> `{rec}`  {truth_s}"
        )
        print(
            f"  top: UID {top['miner_uid']}  score={top['final_score']:.4f}  "
            f"n={top['n_sites']}  ({top['panel_class']})"
        )
        if nat:
            print(
                f"  native #1: UID {nat['miner_uid']}  score={nat_score:.4f}  n={nat_n}"
            )
        print(
            f"  your best: UID {your_uid}  score={your_best:.4f}  "
            f"beat_native={'YES' if beat else 'no'}  gap={gap:.4f}"
        )
        dups = rd.get("fleet_duplicate_panels") or []
        if dups:
            print(f"  fleet duplicate panels: {len(dups)} groups")

    print("\n" + "-" * 72)
    print(f"Beat native #1 on {wins}/{len(rounds)} rounds (fleet UIDs)")
    cal_path = RESULTS / "niome_challenge_db" / "training" / "strategy_calibration.json"
    if cal_path.is_file():
        cal = json.loads(cal_path.read_text(encoding="utf-8"))
        print("\nHistorical curriculum targets (from all truth rounds):")
        for band, info in cal.items():
            print(
                f"  {band:5s}  strategy={info['strategy']:12s}  "
                f"target={info.get('curriculum_target_suggested')}"
            )


if __name__ == "__main__":
    main()
