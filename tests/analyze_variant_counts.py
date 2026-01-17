#!/usr/bin/env python3
"""Variant count per task: manager truth, region, top miners from logs."""
import json
import os
import re
from collections import defaultdict

ROOT = os.path.join(os.path.dirname(__file__), "..", "Results")


def count_vcf_lines(path):
    if not os.path.isfile(path):
        return None
    n = 0
    snp = indel = 0
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        p = line.split("\t")
        if len(p) < 5:
            continue
        n += 1
        ref, alt = p[3], p[4].split(",")[0]
        if len(ref) == 1 and len(alt) == 1:
            snp += 1
        else:
            indel += 1
    return {"n": n, "snp": snp, "indel": indel}


def region_bp(region):
    _, rest = region.split(":")
    a, b = rest.split("-")
    return int(b) - int(a)


def parse_vcf_from_log(log):
    if "Miner VCF\n" not in log:
        return []
    text = log.split("Miner VCF\n", 1)[1]
    rows = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        p = line.split("\t")
        if len(p) >= 5:
            rows.append((int(p[1]), p[3], p[4].split(",")[0]))
    return rows


def load_tasks():
    rows = []
    for name in sorted(os.listdir(ROOT)):
        path = os.path.join(ROOT, name, "task.json")
        if not os.path.isfile(path):
            continue
        with open(path) as f:
            t = json.load(f)
        region = t["genome_context"]["region"]
        truth_path = os.path.join(ROOT, name, "real_correct_result", "truth.vcf")
        truth = count_vcf_lines(truth_path)
        rows.append({
            "round": name,
            "task_id": t["task_id"][:8],
            "region": region,
            "bp": region_bp(region),
            "expected": t.get("expected_variant_count", "?"),
            "truth": truth,
        })
    # old_version truth folders
    old = os.path.join(ROOT, "old_version_result")
    if os.path.isdir(old):
        for name in sorted(os.listdir(old)):
            truth_path = os.path.join(old, name, "truth.vcf")
            if not os.path.isfile(truth_path):
                continue
            task_path = os.path.join(old, name, "task.json")
            region = "?"
            expected = "?"
            tid = name
            if os.path.isfile(task_path):
                with open(task_path) as f:
                    t = json.load(f)
                region = t["genome_context"]["region"]
                expected = t.get("expected_variant_count", "?")
                tid = t["task_id"][:8]
            truth = count_vcf_lines(truth_path)
            rows.append({
                "round": f"old/{name}",
                "task_id": tid,
                "region": region,
                "bp": region_bp(region) if ":" in region else None,
                "expected": expected,
                "truth": truth,
            })
    return rows


def load_miner_counts(task_id_prefix=None):
    """Per task_id: top miner variant counts from miner.json files."""
    by_task = defaultdict(list)
    for name in sorted(os.listdir(ROOT)):
        mpath = os.path.join(ROOT, name, "miner.json")
        if not os.path.isfile(mpath):
            continue
        with open(mpath) as f:
            data = json.load(f)
        for it in data.get("items", []):
            tid = it.get("task_id", "")
            if task_id_prefix and not tid.startswith(task_id_prefix):
                continue
            vcf = parse_vcf_from_log(it.get("log", ""))
            if not vcf:
                continue
            by_task[tid].append({
                "round": name,
                "uid": it.get("miner_uid"),
                "final": it.get("final_score", 0),
                "vcf_score": it.get("vcf_score", 0),
                "n": len(vcf),
                "validator": it.get("validator_uid"),
            })
    return by_task


def main():
    tasks = load_tasks()
    print("=" * 72)
    print("MANAGER TRUTH / TASK SHAPE")
    print("=" * 72)
    print(f"{'round':<14} {'task':<10} {'bp':>7} {'expected':>8} {'truth':>6} {'snp':>4} {'indel':>5} region")
    for r in tasks:
        tr = r["truth"]
        if tr:
            print(
                f"{r['round']:<14} {r['task_id']:<10} {r['bp'] or 0:>7} {str(r['expected']):>8} "
                f"{tr['n']:>6} {tr['snp']:>4} {tr['indel']:>5} {r['region']}"
            )
        else:
            print(f"{r['round']:<14} {r['task_id']:<10} {r['bp'] or 0:>7} {str(r['expected']):>8} {'—':>6} {'—':>4} {'—':>5} {r['region']}")

    print("\n" + "=" * 72)
    print("TOP MINERS ON LIVE ROUNDS (best final_score per task, validator 119)")
    print("=" * 72)
    by_task = load_miner_counts()
    for tid, entries in sorted(by_task.items()):
        v119 = [e for e in entries if e["validator"] == 119]
        if not v119:
            continue
        best = max(v119, key=lambda x: x["final"])
        top_n = sorted(v119, key=lambda x: x["final"], reverse=True)[:5]
        counts = [e["n"] for e in v119 if e["final"] > 0.5]
        print(f"\ntask {tid}…  round(s): {sorted({e['round'] for e in v119})}")
        print(f"  best UID {best['uid']}: final={best['final']:.4f} vcf_n={best['n']}")
        print(f"  top-5 counts: {[e['n'] for e in top_n]}")
        if counts:
            print(f"  high-scorer counts (final>0.5): min={min(counts)} max={max(counts)} median={sorted(counts)[len(counts)//2]}")

    # UID 71 on 5.21.04
    print("\n" + "=" * 72)
    print("YOUR UID 71 (5.21.04)")
    print("=" * 72)
    for tid, entries in by_task.items():
        u71 = [e for e in entries if e["uid"] == 71 and "5.21.04" in e["round"]]
        for e in u71:
            print(f"  validator {e['validator']}: n={e['n']} final={e['final']:.4f} vcf={e['vcf_score']:.4f}")

    print("\n" + "=" * 72)
    print("PATTERN SUMMARY")
    print("=" * 72)
    truths = [r["truth"]["n"] for r in tasks if r.get("truth")]
    bps = [r["bp"] for r in tasks if r.get("truth") and r["bp"]]
    if truths:
        print(f"Manager truth counts: {truths}  (min={min(truths)} max={max(truths)})")
    same_region = [r for r in tasks if r.get("region") == "chr7:117480000-117670000"]
    if len(same_region) >= 2:
        print(f"Same ultra-wide region chr7:117480000-117670000: truth N = "
              f"{[r['truth']['n'] for r in same_region if r['truth']]}")
    print("""
No public formula in task JSON (expected_variant_count=0 on v2.1 live tasks).

What we know:
  - Truth count is NOT fixed; observed manager truth: 10–25 for CFTR tasks.
  - Same 190kb region can be 12, 13, or 25 variants (read sample dependent).
  - Validator count_penalty: min(miner_n, truth_n) / max(miner_n, truth_n)
  - Top miners on 5.21.04 clustered ~20 variants; NOT a subnet constant.

Do NOT submit a fixed N. Match read evidence + stay near truth count.
""")


if __name__ == "__main__":
    main()
