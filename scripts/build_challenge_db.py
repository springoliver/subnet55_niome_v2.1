#!/usr/bin/env python3
"""
Build cumulative NIOME challenge database from all Results/ folders.

  python scripts/build_challenge_db.py
  python scripts/build_challenge_db.py --results-root /path/to/Results
  python scripts/build_challenge_db.py --rounds 5.24.01 5.24.02

Output: Results/niome_challenge_db/
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _bootstrap_analysis():
    """Import analysis package without loading bittensor via niome_subnet.__init__."""
    import types

    if "niome_subnet.analysis.challenge_db" in sys.modules:
        from niome_subnet.analysis.challenge_db import build_database, load_manifest

        return build_database, load_manifest

    ns_pkg = types.ModuleType("niome_subnet")
    ns_pkg.__path__ = [str(ROOT / "niome_subnet")]
    ana_pkg = types.ModuleType("niome_subnet.analysis")
    ana_pkg.__path__ = [str(ROOT / "niome_subnet" / "analysis")]
    sys.modules["niome_subnet"] = ns_pkg
    sys.modules["niome_subnet.analysis"] = ana_pkg

    from niome_subnet.analysis.challenge_db import build_database, load_manifest

    return build_database, load_manifest


build_database, load_manifest = _bootstrap_analysis()


def main() -> None:
    ap = argparse.ArgumentParser(description="Build NIOME challenge Results database")
    ap.add_argument(
        "--results-root",
        type=Path,
        default=ROOT / "Results",
        help="Path to Results/ directory",
    )
    ap.add_argument(
        "--rounds",
        nargs="*",
        help="Only ingest these round paths (e.g. 5.24.01 old_version_result/5.19.01)",
    )
    args = ap.parse_args()

    db_dir = build_database(
        results_root=args.results_root,
        only_rounds=args.rounds,
    )
    manifest = load_manifest(args.results_root)
    print(f"Built: {db_dir}")
    print(f"  rounds: {manifest['n_rounds']}  with_truth: {manifest['n_with_truth']}")
    print(f"  training: {db_dir / 'training'}")
    print("\nRounds indexed:")
    for r in manifest["rounds"]:
        truth = f"truth={r['truth_n']}" if r["has_truth"] else "no truth"
        print(
            f"  {r['path']:32s}  band={r['band']:5s}  "
            f"top_n={r['top_n']:2d}  {truth}"
        )


if __name__ == "__main__":
    main()
