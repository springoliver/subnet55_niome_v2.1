#!/usr/bin/env python3
"""
Aggregate fleet performance BY configured strategy across all Results/ rounds.

Use this to decide how many strategies to run and what each should optimize —
not to copy last round's variant count.

  python tests/analyze_strategy_outcomes.py
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "Results" / "niome_challenge_db" / "rounds"
FLEET = {
    50: "high_recall",
    99: "v5_style",
    209: "high_recall",
    235: "win",
    226: "win",
    124: "v5_style",
    9: "v10",
    97: "v5_style",
    217: "v10",
    141: "auto",
    71: "v5_style",
}


def main() -> None:
    by_strat: dict = defaultdict(
        lambda: {
            "final": [],
            "vcf": [],
            "ann": [],
            "sites": [],
            "rounds": [],
        }
    )

    for path in sorted(ROOT.glob("*.json")):
        rd = json.loads(path.read_text(encoding="utf-8"))
        rk = rd.get("round_key", path.stem)
        top_n = (rd.get("top_miner") or {}).get("n_sites")
        for uid_s, rec in (rd.get("fleet") or {}).items():
            uid = int(uid_s)
            if uid not in FLEET:
                continue
            st = rec.get("configured_strategy") or FLEET[uid]
            if st in ("?", "auto"):
                st = FLEET[uid]
            by_strat[st]["final"].append(float(rec["final_score"]))
            by_strat[st]["vcf"].append(float(rec["vcf_score"]))
            by_strat[st]["ann"].append(float(rec["annotation_score"]))
            by_strat[st]["sites"].append(int(rec["n_sites"]))
            by_strat[st]["rounds"].append(
                {
                    "round": rk,
                    "final": rec["final_score"],
                    "sites": rec["n_sites"],
                    "top_n": top_n,
                }
            )

    print("=" * 72)
    print("FLEET OUTCOMES BY STRATEGY (all indexed rounds)")
    print("=" * 72)
    for st in sorted(by_strat.keys()):
        d = by_strat[st]
        n = len(d["final"])
        if not n:
            continue
        finals = d["final"]
        sites = d["sites"]
        print(f"\n  {st}  (n={n})")
        print(
            f"    final  mean={sum(finals)/n:.4f}  max={max(finals):.4f}  "
            f"min={min(finals):.4f}"
        )
        print(
            f"    vcf    mean={sum(d['vcf'])/n:.4f}  ann mean={sum(d['ann'])/n:.4f}"
        )
        print(
            f"    sites  min={min(sites)} max={max(sites)} "
            f"(top miner N varies per round — do not fix fleet to one N)"
        )
        best = max(d["rounds"], key=lambda x: x["final"])
        print(f"    best round: {best['round']} final={best['final']:.4f} sites={best['sites']}")

    print("\n" + "=" * 72)
    print("Takeaway: optimize VCF+GT per strategy axis; site count follows truth each task.")
    print("=" * 72)


if __name__ == "__main__":
    main()
