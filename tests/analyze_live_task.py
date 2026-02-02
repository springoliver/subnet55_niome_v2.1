#!/usr/bin/env python3
"""
Preview live task: strategy routing, historical DB context, fleet plan.

  python tests/analyze_live_task.py
  python tests/analyze_live_task.py path/to/task.json
"""
from __future__ import annotations

import json
import os
import sys
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "Results"
GENOMICS = ROOT / "niome_subnet" / "genomics"
ANALYSIS = ROOT / "niome_subnet" / "analysis"


def _bootstrap():
    ns = types.ModuleType("niome_subnet")
    ns.__path__ = [str(ROOT / "niome_subnet")]
    ana = types.ModuleType("niome_subnet.analysis")
    ana.__path__ = [str(ANALYSIS)]
    gen = types.ModuleType("niome_subnet.genomics")
    gen.__path__ = [str(GENOMICS)]
    sys.modules["niome_subnet"] = ns
    sys.modules["niome_subnet.analysis"] = ana
    sys.modules["niome_subnet.genomics"] = gen


def _load_genomics(name: str):
    import importlib.util

    p = GENOMICS / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"niome_genomics.{name}", p)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[f"niome_genomics.{name}"] = mod
    spec.loader.exec_module(mod)
    return mod


def _historical_context(band: str) -> dict:
    cal_path = RESULTS / "niome_challenge_db" / "training" / "strategy_calibration.json"
    if not cal_path.is_file():
        return {}
    cal = json.loads(cal_path.read_text(encoding="utf-8"))
    return cal.get(band, {})


def _compare_recent_tasks(task: dict) -> None:
    tid = task["task_id"]
    r1_base = task["input"]["read1_fastq"].split("?")[0]
    print("\n--- Recent same-region tasks (crt/reads pool) ---")
    refs = [
        "5.24.01",
        "5.23.04",
        "5.23.03",
        "5.23.02",
        "5.22.04",
    ]
    for rd in refs:
        p = RESULTS / rd / "task.json"
        if not p.is_file():
            continue
        ref = json.loads(p.read_text(encoding="utf-8"))
        same_reads = ref["input"]["read1_fastq"].split("?")[0] == r1_base
        print(
            f"  {rd:8s}  task={ref['task_id'][:8]}…  "
            f"same_fastq_key={same_reads}  "
            f"{'SAME task_id' if ref['task_id'] == tid else 'NEW truth loci'}"
        )


def main():
    _bootstrap()
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else RESULTS / "live_task" / "task.json"
    task = json.loads(path.read_text(encoding="utf-8"))

    ts = _load_genomics("task_strategy")
    tp = _load_genomics("task_profile")
    from niome_subnet.genomics.truth_paths import find_task_truth

    region = task["genome_context"]["region"]
    _, rest = region.split(":")
    start, end = map(int, rest.split("-"))
    rlen = end - start

    class _In:
        read1_fastq = task["input"]["read1_fastq"]
        read2_fastq = task["input"]["read2_fastq"]

    class _GC:
        region = region

    class _Task:
        task_id = task["task_id"]
        genome_context = _GC()
        input = _In()
        expected_variant_count = task.get("expected_variant_count", 0)

    fp = ts.fingerprint_task(_Task())
    prof = tp.classify_task(region, 0)
    truth_hit = find_task_truth(task["task_id"]) is not None

    strategies = {
        "auto (springhot)": ts.resolve_strategy(_Task(), truth_available=truth_hit),
        "high_recall": "high_recall",
        "v10": "v10",
        "v5_style": "v5_style",
        "win (spring05)": ts.resolve_strategy(_Task(), truth_available=truth_hit)
        if truth_hit
        else ts.pipeline_fallback_strategy("win", fp.predicted_band, False),
    }

    hist = _historical_context(fp.predicted_band)

    print("=" * 72)
    print("LIVE TASK — fleet strategy plan")
    print("=" * 72)
    print(f"task_id:      {task['task_id']}")
    print(f"region:       {region}  ({rlen:,} bp)")
    print(f"read_key:     {fp.read_key}  (crt/reads — presign date may differ)")
    print(f"band:         {fp.predicted_band}")
    print(f"profile:      {prof.name}")
    print(f"truth_file:   {'yes' if truth_hit else 'NO — win UIDs use high_recall pipeline'}")
    if hist:
        print(
            f"history:      strategy={hist.get('strategy')}  "
            f"curriculum_target≈{hist.get('curriculum_target_suggested')}  "
            f"top_n_median={hist.get('top_site_counts', {}).get('median', '?')}"
        )

    print("\n--- Per-strategy env (after apply_strategy_profile) ---")
    for label, strat in strategies.items():
        ts.apply_strategy_profile(strat)
        pick = ts.pipeline_pick_mode()
        merge = ts.pipeline_merge_pool()
        print(
            f"  {label:22s}  pick={pick:10s}  merge_pool={merge}  "
            f"curriculum={os.environ.get('NIOME_CURRICULUM_TARGET', '?')}"
        )

    _compare_recent_tasks(task)

    print("\n--- Rules for this task ---")
    print("  • NEW task_id → new truth positions; do NOT replay 5.23/5.24 VCF panels.")
    print("  • Same crt/reads FASTQ → same BAM evidence; different loci each round.")
    print("  • Target ~29–33 sites (ultra); high_recall merges all mpileup passes.")
    print("  • After round: save miner.json + task.json under Results/5.24.02/")
    print("    then: python scripts/build_challenge_db.py")

    plan_path = path.parent / "strategy_plan.json"
    plan = {
        "task_id": task["task_id"],
        "region": region,
        "read_key": fp.read_key,
        "band": fp.predicted_band,
        "truth_available": truth_hit,
        "historical": hist,
        "strategies": {k: v for k, v in strategies.items()},
        "fleet_pm2": {
            "high_recall": ["spring01", "spring02", "spring03", "spring04", "miner-2", "miner-4"],
            "v10": ["spring07", "spring09", "miner-1"],
            "v5_style": ["spring06", "spring08", "miner-3"],
            "win": ["spring05", "miner-5"],
            "auto": ["springhot"],
        },
    }
    plan_path.write_text(json.dumps(plan, indent=2), encoding="utf-8")
    print(f"\nWrote: {plan_path}")
    print("=" * 72)


if __name__ == "__main__":
    main()
