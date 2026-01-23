#!/usr/bin/env python3
"""Fetch live scores from niome-api (same data W&B table uses)."""
import json
import re
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

API = "https://niome-api.genomes.io/api/miner_scores"
ROOT = Path(__file__).resolve().parents[1]
USER_UIDS = {225, 38, 155, 139, 235, 209, 71}


def parse_vcf(log: str):
    m = re.search(r"Miner VCF\n(.*)", log, re.S)
    if not m:
        return []
    rows = []
    for line in m.group(1).splitlines():
        if line.startswith("chr7"):
            p = line.split("\t")
            rows.append((int(p[1]), p[3], p[4]))
    return rows


def detect_src(log: str) -> str:
    if "niome-native" in log:
        m = re.search(r"niome-native-[\d-]+-v\d+", log)
        return m.group(0) if m else "niome-native"
    if "niome-miner-v2" in log or "bcftoolsVersion" in log:
        return "bcftools-oracle"
    if "##source=niome_miner" in log:
        return "niome_miner_legacy"
    if parse_vcf(log) and "DP=" not in log:
        return "minimal-gt-oracle"
    return "other"


def main():
    task_path = ROOT / "Results" / "live_task" / "task.json"
    task_id = None
    if task_path.is_file():
        task_id = json.loads(task_path.read_text()).get("task_id")
    print("=" * 78)
    print("NIOME LIVE SCORES (niome-api -> W&B table source)")
    print("=" * 78)
    print(f"Filter task_id: {task_id or '(all recent)'}")

    with urllib.request.urlopen(API, timeout=60) as resp:
        data = json.loads(resp.read().decode())

    items = data.get("items", data if isinstance(data, list) else [])
    print(f"Total records in API: {len(items)}")

    # Group by task_id
    by_task = defaultdict(list)
    for it in items:
        by_task[it.get("task_id", "?")].append(it)

    tasks_sorted = sorted(
        by_task.keys(),
        key=lambda t: max(x.get("created_at", "") for x in by_task[t]),
        reverse=True,
    )
    print("\nRecent tasks:")
    for tid in tasks_sorted[:6]:
        rows = by_task[tid]
        best = max(rows, key=lambda x: x.get("final_score", 0))
        print(
            f"  {tid[:8]}…  n={len(rows)}  "
            f"top UID {best.get('miner_uid')} final={best.get('final_score', 0):.4f}  "
            f"at {rows[0].get('created_at', '')[:19]}"
        )

    target = task_id or (tasks_sorted[0] if tasks_sorted else None)
    if not target:
        print("No tasks found")
        return

    rows = by_task[target]
    # best score per uid
    by_uid = {}
    for it in rows:
        uid = it["miner_uid"]
        fin = float(it.get("final_score", 0) or 0)
        if uid not in by_uid or fin > by_uid[uid]["final"]:
            by_uid[uid] = {
                "final": fin,
                "vcf": float(it.get("vcf_score", 0) or 0),
                "ann": float(it.get("annotation_score", 0) or 0),
                "prec": float(it.get("precision", 0) or 0),
                "rec": float(it.get("recall", 0) or 0),
                "weight": float(it.get("weight", 0) or 0),
                "log": it.get("log", ""),
                "val": it.get("validator_uid"),
            }

    ranked = sorted(by_uid.items(), key=lambda x: -x[1]["final"])
    print(f"\n=== TASK {target} ===")
    print(f"Miners: {len(by_uid)}  validator UIDs: {sorted({r['val'] for r in by_uid.values()})}")

    print("\nTOP 15:")
    for i, (uid, s) in enumerate(ranked[:15], 1):
        v = parse_vcf(s["log"])
        src = detect_src(s["log"])
        print(
            f"  {i:2}. UID {uid:3} final={s['final']:.4f} vcf={s['vcf']:.4f} "
            f"ann={s['ann']:.4f} rec={s['rec']:.3f} prec={s['prec']:.3f} "
            f"n={len(v)} wt={s['weight']:.4f} {src[:28]}"
        )

    print("\nYOUR UIDs:")
    for uid in sorted(USER_UIDS):
        if uid not in by_uid:
            print(f"  UID {uid}: (no score yet)")
            continue
        s = by_uid[uid]
        v = parse_vcf(s["log"])
        rank = next(i for i, (u, _) in enumerate(ranked, 1) if u == uid)
        print(
            f"  UID {uid}: rank #{rank}/{len(ranked)} final={s['final']:.4f} "
            f"vcf={s['vcf']:.4f} n={len(v)} {detect_src(s['log'])[:32]}"
        )

    # Panel clusters in top 20
    panels = defaultdict(list)
    for uid, s in ranked[:30]:
        key = tuple(sorted(parse_vcf(s["log"])))
        panels[key].append(uid)
    print("\nTop-panel clusters (identical VCF):")
    for i, (key, uids) in enumerate(
        sorted(panels.items(), key=lambda x: -len(x[1]))[:4], 1
    ):
        print(f"  cluster {i}: {len(key)} sites, {len(uids)} miners -> UIDs {uids[:12]}")

    counts = Counter(len(parse_vcf(s["log"])) for s in by_uid.values())
    print("\nVariant count distribution:")
    for n, c in sorted(counts.items()):
        print(f"  {n} variants: {c} miners")

    oracle = sum(1 for s in by_uid.values() if detect_src(s["log"]) == "minimal-gt-oracle")
    native = sum(1 for s in by_uid.values() if "niome-native" in detect_src(s["log"]))
    print(f"\nMethod mix: minimal-gt-oracle={oracle}  niome-native={native}  other={len(by_uid)-oracle-native}")
    print("=" * 78)
    print("W&B: https://wandb.ai/genomes/niome/table?nw=nwuserjgenome")
    print("=" * 78)


if __name__ == "__main__":
    main()
