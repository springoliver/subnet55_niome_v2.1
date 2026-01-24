#!/usr/bin/env python3
"""Analyze user miner UIDs across miner.json + validator.json."""
import json
import re
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "Results"

USER_UIDS = {139, 71, 50, 99, 209, 235, 226, 124, 9, 97, 217, 36, 92, 225, 38, 155}
ZERO_SUSPECT = {226, 124, 9, 97, 217}
VALIDATORS = {119, 154, 58}


def parse_vcf(log: str):
    m = re.search(r"Miner VCF\n(.*)", log, re.S)
    if not m:
        return [], "no_vcf"
    rows = []
    for line in m.group(1).splitlines():
        if line.startswith("chr7"):
            p = line.split("\t")
            rows.append((int(p[1]), p[3], p[4]))
    return rows, "ok"


def classify(log: str) -> str:
    if not log or "Miner VCF" not in log:
        return "empty/error"
    if "niome-native" in log or "niome_miner_niome-native" in log:
        return "niome-native"
    if "niome-competitive" in log:
        return "competitive"
    if "CLNDN=Cystic_fibrosis" in log:
        return "oracle-clndn"
    if "bcftoolsVersion" in log:
        return "bcftools-read"
    rows, st = parse_vcf(log)
    if st == "no_vcf":
        return "no_vcf"
    if rows and "DP=" not in log and "GT:DP" not in log:
        return "minimal-oracle"
    if rows:
        return "read-other"
    return "empty"


def load_items(path: Path):
    if not path.is_file():
        return []
    data = json.loads(path.read_text(encoding="utf-8", errors="replace"))
    return data.get("items", data if isinstance(data, list) else [])


def best_per_uid(items, task_id=None):
    by = {}
    for it in items:
        if task_id and it.get("task_id") != task_id:
            continue
        uid = it.get("miner_uid")
        if uid is None:
            continue
        fin = float(it.get("final_score", 0) or 0)
        if uid not in by or fin > by[uid]["final"]:
            log = it.get("log", "")
            vcf, _ = parse_vcf(log)
            by[uid] = {
                "final": fin,
                "vcf": float(it.get("vcf_score", 0) or 0),
                "ann": float(it.get("annotation_score", 0) or 0),
                "rec": float(it.get("recall", 0) or 0),
                "prec": float(it.get("precision", 0) or 0),
                "n": len(vcf),
                "cls": classify(log),
                "val": it.get("validator_uid"),
                "rt": it.get("response_time"),
                "task": it.get("task_id", "")[:8],
                "log_head": log[:120].replace("\n", " "),
            }
    return by


def main():
    rounds = []
    for d in sorted(RESULTS.iterdir()):
        if not d.is_dir() or d.name.startswith("old"):
            continue
        mj = d / "miner.json"
        vj = d / "validator.json"
        tj = d / "task.json"
        if mj.is_file() or vj.is_file():
            task_id = None
            if tj.is_file():
                task_id = json.loads(tj.read_text()).get("task_id")
            rounds.append((d.name, task_id, mj, vj))

    print("=" * 88)
    print("USER MINER ANALYSIS")
    print("UIDs:", sorted(USER_UIDS))
    print("=" * 88)

    latest_task = None
    for name, tid, mj, vj in rounds:
        print(f"\n### {name}  task={tid[:8] if tid else '?'}...")
        items_m = load_items(mj) if mj.is_file() else []
        items_v = load_items(vj) if vj.is_file() else []
        bm = best_per_uid(items_m, tid)
        bv = best_per_uid(items_v, tid)  # validator.json may mix tasks

        # validator.json often has leaderboard not filtered by task - take best any
        if not bm and items_m:
            bm = best_per_uid(items_m)
        if not bv and items_v:
            bv = best_per_uid(items_v)

        if tid:
            latest_task = (name, tid, bm, bv)

        print(f"  miner.json rows={len(items_m)}  validator.json rows={len(items_v)}")
        print(f"  {'UID':>4} {'final':>7} {'vcf':>7} {'n':>4} {'class':<16} source")
        for uid in sorted(USER_UIDS):
            s = bm.get(uid) or bv.get(uid)
            src = "miner" if uid in bm else ("validator" if uid in bv else "-")
            if not s:
                mark = " *** MISSING" if uid in ZERO_SUSPECT else ""
                print(f"  {uid:4}   -       -      -   not in files{mark}")
                continue
            z = " *** ZERO" if s["final"] == 0 else ""
            print(
                f"  {uid:4} {s['final']:7.4f} {s['vcf']:7.4f} {s['n']:4} "
                f"{s['cls']:<16} {src}{z}"
            )

    # Latest round detail for zero suspects
    if latest_task:
        name, tid, bm, bv = latest_task
        print(f"\n{'=' * 88}")
        print(f"ZERO / ERROR DETAIL — latest {name} task {tid[:8]}")
        print("=" * 88)
        items_m = load_items(RESULTS / name / "miner.json")
        for uid in sorted(ZERO_SUSPECT):
            entries = [
                it
                for it in items_m
                if it.get("miner_uid") == uid and it.get("task_id") == tid
            ]
            if not entries:
                entries = [it for it in items_m if it.get("miner_uid") == uid]
            print(f"\nUID {uid}: {len(entries)} record(s)")
            for it in sorted(entries, key=lambda x: -float(x.get("final_score", 0)))[:2]:
                log = it.get("log", "")
                print(
                    f"  val={it.get('validator_uid')} final={it.get('final_score')} "
                    f"vcf={it.get('vcf_score')} prec={it.get('precision')} rec={it.get('recall')}"
                )
                print(f"  class={classify(log)}")
                if "Miner VCF" in log:
                    rows, _ = parse_vcf(log)
                    print(f"  variants={len(rows)}")
                else:
                    print(f"  log: {log[:200]}")

    # Top on latest miner.json
    if latest_task:
        name, tid, bm, _ = latest_task
        all_u = best_per_uid(load_items(RESULTS / name / "miner.json"), tid)
        ranked = sorted(all_u.items(), key=lambda x: -x[1]["final"])
        print(f"\n{'=' * 88}")
        print(f"TOP 10 on {name} (all miners)")
        for i, (u, s) in enumerate(ranked[:10], 1):
            tag = " YOU" if u in USER_UIDS else ""
            print(f"  {i:2}. UID {u:3} final={s['final']:.4f} n={s['n']} {s['cls']}{tag}")

    print("=" * 88)


if __name__ == "__main__":
    main()
