#!/usr/bin/env python3
"""Analyze Results/5.22.04 miner.json + truth trend."""
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MINER = ROOT / "Results" / "5.22.04" / "miner.json"
TASK = ROOT / "Results" / "5.22.04" / "task.json"
TRUTH_503 = ROOT / "Results" / "5.22.03" / "real_correct_result" / "truth.vcf"
USER_UIDS = {225, 38, 155}


def parse_vcf_from_log(log: str):
    m = re.search(r"Miner VCF\n(.*)", log, re.DOTALL)
    if not m:
        return []
    rows = []
    for line in m.group(1).splitlines():
        if line.startswith("chr7"):
            p = line.split("\t")
            gt = "."
            if len(p) > 9:
                fmt = p[8].split(":")
                samp = p[9].split(":")
                if "GT" in fmt:
                    gt = samp[fmt.index("GT")]
            rows.append((int(p[1]), p[3], p[4], gt))
    return rows


def parse_truth(path):
    rows = []
    for line in path.read_text().splitlines():
        if line.startswith("chr7"):
            p = line.split("\t")
            rows.append((int(p[1]), p[3], p[4], p[9]))
    return rows


def jaccard(a, b):
    sa = {(x[0], x[1], x[2]) for x in a}
    sb = {(x[0], x[1], x[2]) for x in b}
    if not sa | sb:
        return 0.0
    return len(sa & sb) / len(sa | sb)


def main():
    task = json.loads(TASK.read_text())
    data = json.loads(MINER.read_text())
    items = data.get("items", data if isinstance(data, list) else [])

    print("=" * 72)
    print("5.22.04 ROUND ANALYSIS")
    print("=" * 72)
    print(f"task_id: {task['task_id']}")
    print(f"region:  {task['genome_context']['region']}")
    print(f"Note: Same task_id/FASTQs as live_task (a9a4db4a), NOT a new sample folder name only.")

    # one row per uid (best final per uid)
    by_uid = {}
    for it in items:
        uid = it.get("miner_uid")
        fin = it.get("final_score", 0)
        if uid not in by_uid or fin > by_uid[uid]["final"]:
            by_uid[uid] = {
                "final": fin,
                "vcf": it.get("vcf_score", 0),
                "ann": it.get("annotation_score", 0),
                "recall": it.get("recall", 0),
                "prec": it.get("precision", 0),
                "vid": it.get("validator_uid"),
                "log": it.get("log", ""),
                "rev": "",
            }
            src = re.search(r"niome-native-[\d-]+-v\d+", it.get("log", ""))
            if src:
                by_uid[uid]["rev"] = src.group(0)

    ranked = sorted(by_uid.items(), key=lambda x: -x[1]["final"])

    print(f"\nMiners scored: {len(by_uid)}")
    print("\nTOP 12:")
    for i, (uid, s) in enumerate(ranked[:12], 1):
        v = parse_vcf_from_log(s["log"])
        indels = sum(1 for x in v if len(x[1]) > 2 or len(x[2]) > 2)
        hom = sum(1 for x in v if x[3] == "1/1")
        print(
            f"  {i:2}. UID {uid:3} final={s['final']:.4f} vcf={s['vcf']:.4f} "
            f"ann={s['ann']:.4f} n={len(v)} indels={indels} hom={hom} rev={s['rev'][:30]}"
        )

    print("\nYOUR UIDs:")
    for uid in sorted(USER_UIDS):
        if uid not in by_uid:
            print(f"  UID {uid}: not in miner.json")
            continue
        s = by_uid[uid]
        v = parse_vcf_from_log(s["log"])
        indels = sum(1 for x in v if len(x[1]) > 2 or len(x[2]) > 2)
        hom = sum(1 for x in v if x[3] == "1/1")
        rank = next(i for i, (u, _) in enumerate(ranked, 1) if u == uid)
        print(
            f"  UID {uid}: rank #{rank}/{len(ranked)} final={s['final']:.4f} "
            f"vcf={s['vcf']:.4f} recall={s['recall']:.3f} prec={s['prec']:.3f} "
            f"n={len(v)} indels={indels} 1/1={hom} rev={s['rev']}"
        )

    # variant count clusters
    counts = Counter()
    for uid, s in by_uid.items():
        n = len(parse_vcf_from_log(s["log"]))
        counts[n] += 1
    print("\nVariant count distribution (all miners):")
    for n, c in sorted(counts.items()):
        print(f"  {n} variants: {c} miners")

    # same VCF cluster?
    sigs = defaultdict(list)
    for uid, s in by_uid.items():
        v = parse_vcf_from_log(s["log"])
        key = tuple(sorted((x[0], x[1], x[2], x[3]) for x in v))
        sigs[len(key)].append(uid)
    print("\nPanel sizes (unique position panels):")
    panels = defaultdict(list)
    for uid, s in by_uid.items():
        v = parse_vcf_from_log(s["log"])
        key = tuple(sorted((x[0], x[1], x[2]) for x in v))
        panels[key].append(uid)
    by_size = sorted(panels.items(), key=lambda x: -len(x[0]))
    for i, (key, uids) in enumerate(by_size[:5], 1):
        print(f"  panel {i}: {len(key)} sites, {len(uids)} miners (e.g. UID {uids[:6]})")

    if TRUTH_503.is_file():
        truth = parse_truth(TRUTH_503)
        print(f"\n--- vs 5.22.03 truth (task 7fc3be20, n={len(truth)}) ---")
        print("WARNING: 5.22.04 is task a9a4db4a — truth may NOT match if sample changed.")
        for uid in sorted(USER_UIDS):
            if uid not in by_uid:
                continue
            v = parse_vcf_from_log(by_uid[uid]["log"])
            j = jaccard(v, truth)
            tk = {(a, b, c) for a, b, c, _ in v}
            tr = {(a, b, c) for a, b, c, _ in truth}
            print(
                f"  UID {uid}: jaccard={j:.3f} TP={len(tk & tr)} "
                f"FP={len(tk - tr)} FN={len(tr - tk)}"
            )
        top_uid = ranked[0][0]
        vtop = parse_vcf_from_log(ranked[0][1]["log"])
        print(
            f"  TOP UID {top_uid}: jaccard={jaccard(vtop, truth):.3f} n={len(vtop)}"
        )

    print("\n--- Manager truth count trend (ultra-wide era) ---")
    rounds = [
        ("5.21.01", ROOT / "Results/5.21.01/real_correct_result/truth.vcf"),
        ("5.21.02", ROOT / "Results/5.21.02/real_correct_result/truth.vcf"),
        ("5.21.03", ROOT / "Results/5.21.03/real_correct_result/truth.vcf"),
        ("5.22.03", TRUTH_503),
    ]
    ns = []
    for name, p in rounds:
        if p.is_file():
            n = len(parse_truth(p))
            ns.append(n)
            print(f"  {name}: {n} truth variants")
    if len(ns) >= 2:
        print(f"  deltas: {[ns[i] - ns[i - 1] for i in range(1, len(ns))]}")
        print("  User hypothesis +2/round: next might be ~28 if trend continues")
        print("  BUT truth N is per-sample, not a formula from task index.")

    print("=" * 72)


if __name__ == "__main__":
    main()
