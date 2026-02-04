#!/usr/bin/env python3
"""
Compare truth (when present) vs top miner vs fleet — all DB rounds or subset.

Usage:
  python scripts/build_challenge_db.py
  python tests/analyze_truth_vs_top.py
  python tests/analyze_truth_vs_top.py 5.23.04 5.24.01
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
USER_UIDS = {141, 139, 71, 50, 99, 209, 235, 226, 124, 9, 97, 217, 36, 92, 225, 38, 155}


def main():
    _bootstrap()
    from niome_subnet.analysis.challenge_db import load_manifest

    try:
        manifest = load_manifest(RESULTS)
    except FileNotFoundError:
        print("Run: python scripts/build_challenge_db.py")
        sys.exit(1)

    want = set(sys.argv[1:]) if len(sys.argv) > 1 else None
    rounds_meta = manifest["rounds"]
    if want:
        rounds_meta = [
            m
            for m in rounds_meta
            if m["key"] in want
            or m["path"] in want
            or any(w in m["path"] for w in want)
        ]

    for meta in rounds_meta:
        rd = json.loads(
            (RESULTS / "niome_challenge_db" / "rounds" / f"{meta['key']}.json").read_text(
                encoding="utf-8"
            )
        )
        top = rd["top_miner"]
        print("=" * 72)
        print(
            f"{rd['round_path']}  task={rd['task_id'][:8]}…  band={rd['band']}  "
            f"has_truth={rd['has_truth']}"
        )
        if rd["has_truth"]:
            print(f"  TRUTH sites={rd['truth_n']}")
            top_vs = top.get("vs_truth") or {}
            print(
                f"  TOP uid={top['miner_uid']}  score={top['final_score']:.4f}  "
                f"n={top['n_sites']}  jaccard_truth={top_vs.get('jaccard', 0):.3f}"
            )
        else:
            print(
                f"  TOP uid={top['miner_uid']}  score={top['final_score']:.4f}  "
                f"n={top['n_sites']}  (add truth.vcf + rebuild DB)"
            )

        print("  --- FLEET ---")
        for uid in sorted(USER_UIDS):
            key = str(uid)
            if key not in rd.get("fleet", {}):
                continue
            f = rd["fleet"][key]
            line = (
                f"    UID {uid:3d} ({f.get('configured_strategy', '?'):12s}) "
                f"score={f['final_score']:.4f}  sites={f['n_sites']:2d}"
            )
            if f.get("vs_truth"):
                line += f"  jaccard={f['vs_truth']['jaccard']:.3f}"
            print(line)

        dups = rd.get("fleet_duplicate_panels") or []
        if dups:
            print("  --- DUPLICATE FLEET PANELS ---")
            for g in dups:
                print(f"    n={g['n_sites']} UIDs={g['uids']}")


if __name__ == "__main__":
    main()
