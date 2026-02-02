#!/usr/bin/env python3
"""Compare fleet UIDs: scores, site counts, VCF fingerprints per round."""
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "Results"
USER_UIDS = {139, 71, 50, 99, 209, 235, 226, 124, 9, 97, 217, 36, 92, 225, 38, 155}
VALIDATORS = {119, 154, 58}
FLEET = {
    50: "spring01/high_recall",
    99: "spring02/high_recall",
    209: "spring03/high_recall",
    235: "spring04/high_recall",
    226: "spring05/win",
    124: "spring06/v5_style",
    9: "spring07/v10",
    97: "spring08/v5_style",
    217: "spring09/v10",
    139: "springhot/auto",
}


def panel(log: str):
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


def analyze_round(rd: str):
    tj = json.loads((ROOT / rd / "task.json").read_text(encoding="utf-8"))
    tid = tj["task_id"]
    mj = json.loads((ROOT / rd / "miner.json").read_text(encoding="utf-8", errors="replace"))
    items = mj if isinstance(mj, list) else mj.get("items", mj)
    pool = [
        x
        for x in items
        if x.get("task_id") == tid and x.get("validator_uid") in VALIDATORS
    ]

    print("=" * 72)
    print(f"ROUND {rd}  task={tid[:8]}…  region={tj['genome_context']['region']}")
    print(f"  validator submissions: {len(pool)}")

    by_uid = {}
    for x in pool:
        uid = x["miner_uid"]
        rows = panel(x.get("log", ""))
        by_uid.setdefault(uid, []).append(
            {
                "score": float(x.get("final_score") or 0),
                "vcf": float(x.get("vcf_score") or 0),
                "ann": float(x.get("annotation_score") or 0),
                "n": len(rows),
                "fp": fp(rows),
            }
        )

    tops = sorted(pool, key=lambda x: -float(x.get("final_score") or 0))[:5]
    print("  TOP miners:")
    for t in tops:
        rows = panel(t.get("log", ""))
        cls = "oracle" if "CLNDN" in t.get("log", "") else "native?"
        print(
            f"    uid={t['miner_uid']:3d}  final={float(t['final_score']):.4f}  "
            f"sites={len(rows)}  ({cls})"
        )

    print("  --- YOUR FLEET ---")
    user_entries = []
    for uid in sorted(USER_UIDS):
        if uid not in by_uid:
            continue
        recs = by_uid[uid]
        sc = sum(r["score"] for r in recs) / len(recs)
        n = recs[0]["n"]
        f = recs[0]["fp"]
        strat = FLEET.get(uid, "?")
        user_entries.append((uid, sc, n, f, strat))
        print(
            f"    UID {uid:3d}  {strat:22s}  score={sc:.4f}  "
            f"sites={n:2d}  fp={f}  (n_validators={len(recs)})"
        )

    fp_groups = defaultdict(list)
    for uid, sc, n, f, strat in user_entries:
        fp_groups[(f, n)].append((uid, strat, sc))
    print("  --- DUPLICATE PANELS (your UIDs, same POS/REF/ALT/GT) ---")
    dupes = False
    for (f, n), members in sorted(fp_groups.items(), key=lambda x: -len(x[1])):
        if len(members) > 1:
            dupes = True
            print(
                f"    sites={n} fp={f}: UIDs={[m[0] for m in members]}  "
                f"strategies={sorted({m[1] for m in members})}  "
                f"scores={[round(m[2], 4) for m in members]}"
            )
    if not dupes:
        print("    (none among fleet UIDs)")

    score_groups = defaultdict(list)
    for uid, recs in by_uid.items():
        sc = round(sum(r["score"] for r in recs) / len(recs), 6)
        score_groups[sc].append(uid)
    collisions = [(sc, uids) for sc, uids in score_groups.items() if len(uids) >= 3]
    if collisions:
        print("  --- SCORE COLLISIONS (3+ UIDs, any miner) ---")
        for sc, uids in sorted(collisions, key=lambda x: -len(x[1]))[:10]:
            fleet_hit = [u for u in uids if u in USER_UIDS]
            print(f"    score={sc}: {len(uids)} UIDs total, fleet={fleet_hit}")


def main():
    rounds = sys.argv[1:] or ["5.23.04", "5.24.01"]
    for rd in rounds:
        analyze_round(rd)


if __name__ == "__main__":
    main()
