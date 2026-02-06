#!/usr/bin/env python3
"""Analyze Results/5.24.03 round."""
import hashlib
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "Results" / "5.24.03"
VALIDATORS = {119, 154, 58}
FLEET = {
    141: "springhot/auto",
    50: "spring01/high_recall",
    99: "spring02/v5_style",
    209: "spring03/high_recall",
    235: "spring04/win",
    226: "spring05/win",
    124: "spring06/v5_style",
    9: "spring07/v10",
    97: "spring08/v5_style",
    217: "spring09/v10",
    71: "v5_style",
}


def panel(log):
    if "Miner VCF" not in log:
        return []
    rows = []
    for line in log.split("Miner VCF\n", 1)[1].splitlines():
        if line.startswith("chr7"):
            q = line.split("\t")
            rows.append((int(q[1]), q[3], q[4], q[9] if len(q) > 9 else ""))
    return rows


def fp(rows):
    key = "|".join(f"{p}:{r}>{a}:{g}" for p, r, a, g in sorted(rows))
    return hashlib.md5(key.encode()).hexdigest()[:12]


def source_rev(log):
    for line in log.splitlines():
        if line.startswith("##source="):
            return line.split("=", 1)[1]
    return "?"


def main():
    task = json.loads((ROOT / "task.json").read_text(encoding="utf-8"))
    tid = task["task_id"]
    items = json.loads((ROOT / "miner.json").read_text(encoding="utf-8", errors="replace"))
    if isinstance(items, dict):
        items = items["items"]
    pool = [x for x in items if x.get("task_id") == tid and x.get("validator_uid") in VALIDATORS]

    by_uid = {}
    for x in pool:
        u = x["miner_uid"]
        s = float(x.get("final_score") or 0)
        if u not in by_uid or s > by_uid[u][0]:
            by_uid[u] = (s, x)

    print("=" * 72)
    print(f"5.24.03  task={tid[:8]}...  validators={len(pool)}  uids={len(by_uid)}")
    print("=" * 72)

    tops = sorted(by_uid.items(), key=lambda t: -t[1][0])[:12]
    print("\nTOP 12:")
    for u, (s, x) in tops:
        rows = panel(x.get("log", ""))
        cls = "oracle" if "CLNDN" in x.get("log", "") else "native?"
        print(
            f"  uid={u:3d}  final={s:.4f}  vcf={float(x.get('vcf_score') or 0):.4f}  "
            f"ann={float(x.get('annotation_score') or 0):.2f}  sites={len(rows):2d}  {cls}"
        )

    print("\nYOUR FLEET:")
    fleet_entries = []
    for u in sorted(FLEET):
        if u not in by_uid:
            print(f"  uid={u:3d}  {FLEET[u]:22s}  MISSING")
            continue
        s, x = by_uid[u]
        rows = panel(x.get("log", ""))
        rev = source_rev(x.get("log", ""))
        fleet_entries.append((u, s, len(rows), fp(rows), FLEET[u], rev, x))
        print(
            f"  uid={u:3d}  {FLEET[u]:22s}  final={s:.4f}  sites={len(rows):2d}  "
            f"vcf={float(x.get('vcf_score') or 0):.4f}  ann={float(x.get('annotation_score') or 0):.2f}  "
            f"rev={rev.split('+')[-1] if '+' in rev else rev[-30:]}"
        )

    top_score = tops[0][1][0] if tops else 0
    best_fleet = max((e[1] for e in fleet_entries), default=0)
    print(f"\n  Top subnet: {top_score:.4f}  |  Fleet best: {best_fleet:.4f}  gap: {top_score - best_fleet:.4f}")

    fp_g = defaultdict(list)
    for u, s, n, f, st, rev, _ in fleet_entries:
        fp_g[(f, n)].append((u, st, s, rev))
    print("\nDUPLICATE PANELS (same POS/REF/ALT/GT -> same validator score):")
    for (f, n), members in sorted(fp_g.items(), key=lambda x: -len(x[1])):
        if len(members) > 1:
            print(f"  n={n} fp={f}: {[(m[0], m[1], m[2]) for m in members]}")
            revs = {m[3] for m in members}
            if len(revs) > 1:
                print(f"    different rev tags but identical scored variants: {revs}")


if __name__ == "__main__":
    main()
