#!/usr/bin/env python3
"""
Review past Results/ rounds: task band, oracle vs native #1, your fleet gap.

Run: python tests/analyze_fleet_strategy.py
"""
import json
import re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "Results"
USER_UIDS = {139, 71, 50, 99, 209, 235, 226, 124, 9, 97, 217, 36, 92, 225, 38, 155}
VALIDATORS = {119, 154, 58}

# Map oracle site count → band → recommended strategy
def band(n: int) -> str:
    if n <= 14:
        return "low"
    if n <= 22:
        return "mid"
    if n <= 29:
        return "high"
    return "ultra"


def strategy_for_band(b: str) -> str:
    return {"low": "v5_style", "mid": "v10", "high": "high_recall", "ultra": "high_recall"}[b]


def parse_panel(log: str):
    if "Miner VCF" not in log:
        return [], "empty"
    rows = []
    for line in log.split("Miner VCF\n", 1)[1].splitlines():
        if line.startswith("chr7"):
            q = line.split("\t")
            rows.append((int(q[1]), q[3], q[4]))
    if "CLNDN" in log or "ONCDN" in log:
        cls = "oracle"
    elif "niome-native" in log or "niome_miner" in log:
        cls = "native"
    else:
        cls = "other"
    return rows, cls


def main():
    rounds = sorted(
        d.name for d in ROOT.iterdir() if d.is_dir() and (d / "miner.json").exists()
    )
    print("=" * 72)
    print("NIOME fleet strategy review (historical Results/)")
    print("=" * 72)

    by_band = Counter()
    for rd in rounds:
        mj = ROOT / rd / "miner.json"
        tj = ROOT / rd / "task.json"
        if not tj.exists():
            continue
        task_id = json.loads(tj.read_text(encoding="utf-8")).get("task_id", "")[:8]
        items = json.loads(mj.read_text(encoding="utf-8", errors="replace"))
        items = items if isinstance(items, list) else items.get("items", items)
        pool = [it for it in items if it.get("validator_uid") in VALIDATORS]
        if not pool:
            continue

        top = max(pool, key=lambda x: float(x.get("final_score", 0) or 0))
        top_rows, top_cls = parse_panel(top.get("log", ""))
        b = band(len(top_rows))
        by_band[b] += 1

        natives = []
        for it in pool:
            rows, cls = parse_panel(it.get("log", ""))
            if not rows or cls != "native":
                continue
            natives.append((it["miner_uid"], float(it["final_score"]), len(rows)))

        natives.sort(key=lambda x: -x[1])
        best_nat = natives[0] if natives else (None, 0.0, 0)

        your_best = 0.0
        your_uid = None
        for it in pool:
            if it["miner_uid"] not in USER_UIDS:
                continue
            s = float(it.get("final_score", 0) or 0)
            if s > your_best:
                your_best = s
                your_uid = it["miner_uid"]

        rec = strategy_for_band(b)
        gap = best_nat[1] - your_best if best_nat[0] else 0.0
        beat = "YES" if your_best >= best_nat[1] - 0.001 else "no"

        print(f"\n{rd}  task={task_id}  band={b}  oracle_n={len(top_rows)}  -> strategy `{rec}`")
        print(f"  oracle #1: UID {top['miner_uid']}  final={top['final_score']:.4f}  ({top_cls})")
        if best_nat[0]:
            print(
                f"  native #1: UID {best_nat[0]}  final={best_nat[1]:.4f}  n={best_nat[2]}"
            )
        print(f"  your best: UID {your_uid}  final={your_best:.4f}  beat_native={beat}  gap={gap:.4f}")

    print("\n" + "-" * 72)
    print("Band frequency:", dict(by_band))
    print("\nRecommended PM2 split (16 UIDs):")
    print("  v5_style:     71, 38, 155, 225  (annotation / simpler indels)")
    print("  high_recall:  209, 235, 50, 99   (30+ site rounds)")
    print("  win:          139                 (truth dir when published)")
    print("  v10:          remaining UIDs")
    print("  Or NIOME_STRATEGY=auto on all with NIOME_TRUTH_DIR set")


if __name__ == "__main__":
    main()
