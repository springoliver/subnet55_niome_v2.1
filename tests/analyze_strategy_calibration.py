#!/usr/bin/env python3
"""Cross-round: which fleet strategy won vs top native (from challenge DB)."""
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DB = ROOT / "Results" / "niome_challenge_db"
FLEET = {141, 71, 50, 99, 209, 235, 226, 124, 9, 97, 217}


def main():
    rounds_dir = DB / "rounds"
    if not rounds_dir.is_dir():
        print("Run: python scripts/build_challenge_db.py")
        return

    by_strat = defaultdict(list)
    gaps = []

    for path in sorted(rounds_dir.glob("*.json")):
        rd = json.loads(path.read_text(encoding="utf-8"))
        top = rd.get("native_best") or rd.get("top_miner")
        if not top:
            continue
        top_s = float(top["final_score"])
        fleet = rd.get("fleet", {})
        best_uid, best_rec = None, None
        for uid, rec in fleet.items():
            if int(uid) not in FLEET:
                continue
            s = float(rec["final_score"])
            if best_rec is None or s > float(best_rec["final_score"]):
                best_uid, best_rec = uid, rec
        if not best_rec:
            continue
        st = best_rec.get("configured_strategy", "?")
        fs = float(best_rec["final_score"])
        by_strat[st].append(fs)
        gaps.append(
            {
                "round": rd["round_key"],
                "band": rd["band"],
                "top_n": top["n_sites"],
                "top": top_s,
                "fleet_best": fs,
                "gap": top_s - fs,
                "strategy": st,
                "uid": best_uid,
            }
        )

    print("Fleet best score by configured strategy (mean / max):")
    for st, scores in sorted(by_strat.items(), key=lambda x: -sum(x[1]) / len(x[1])):
        print(f"  {st:14s}  n={len(scores):3d}  mean={sum(scores)/len(scores):.4f}  max={max(scores):.4f}")

    print("\nLargest native gaps (recent):")
    for g in sorted(gaps, key=lambda x: -x["gap"])[:12]:
        print(
            f"  {g['round']:10s} band={g['band']:5s} top={g['top']:.3f} "
            f"fleet={g['fleet_best']:.3f} gap={g['gap']:.3f} "
            f"uid={g['uid']} {g['strategy']} top_n={g['top_n']}"
        )


if __name__ == "__main__":
    main()
