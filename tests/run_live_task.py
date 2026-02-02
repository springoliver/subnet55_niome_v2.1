#!/usr/bin/env python3
"""
Offline pipeline test for Results/live_task/task.json (all fleet strategies).

  python tests/run_live_task.py
  python tests/run_live_task.py high_recall
"""
from __future__ import annotations

import json
import os
import sys
import tempfile
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "Results"
TASK_JSON = RESULTS / "live_task" / "task.json"
STRATEGIES = ("high_recall", "v10", "v5_style", "win")


def _bootstrap():
    ns = types.ModuleType("niome_subnet")
    ns.__path__ = [str(ROOT / "niome_subnet")]
    gen = types.ModuleType("niome_subnet.genomics")
    gen.__path__ = [str(ROOT / "niome_subnet" / "genomics")]
    sys.modules["niome_subnet"] = ns
    sys.modules["niome_subnet.genomics"] = gen
    os.environ.setdefault("NIOME_RESULTS_ROOT", str(RESULTS))
    os.environ.setdefault("NIOME_USE_CHALLENGE_DB", "1")


def _run_one(strategy: str, task_data: dict) -> dict:
    from niome_subnet.genomics.cftr_lookup import build_cftr_annotations
    from niome_subnet.genomics.model import Task
    from niome_subnet.genomics.pipeline import run_pipeline
    from niome_subnet.genomics.task_strategy import (
        apply_strategy_profile,
        fingerprint_task,
        pipeline_fallback_strategy,
    )
    from niome_subnet.genomics.truth_paths import find_task_truth
    from niome_subnet.genomics.read_calling import get_read_calling_rev

    task = Task(**task_data)
    truth_hit = find_task_truth(task.task_id) is not None
    fp = fingerprint_task(task)
    strat = strategy
    if strat == "win" and not truth_hit:
        strat = pipeline_fallback_strategy("win", fp.predicted_band, truth_hit)
    apply_strategy_profile(strat)
    os.environ["NIOME_ACTIVE_BAND"] = fp.predicted_band

    with tempfile.TemporaryDirectory(prefix=f"niome_live_{strat}_") as work_dir:
        final_vcf, _ = run_pipeline(task, work_dir)
        with open(final_vcf, encoding="utf-8") as fh:
            vcf = fh.read()
        annotations = build_cftr_annotations(final_vcf) or {}

    lines = [ln for ln in vcf.splitlines() if ln and not ln.startswith("#")]
    positions = [int(ln.split("\t")[1]) for ln in lines] if lines else []
    return {
        "strategy": strat,
        "rev": get_read_calling_rev(),
        "n_variants": len(lines),
        "n_annotations": len(annotations),
        "pos_min": min(positions) if positions else None,
        "pos_max": max(positions) if positions else None,
        "truth_hit": truth_hit,
    }


def main():
    _bootstrap()
    only = sys.argv[1:] if len(sys.argv) > 1 else list(STRATEGIES)
    task_data = json.loads(TASK_JSON.read_text(encoding="utf-8"))

    print(f"task_id: {task_data['task_id']}")
    print(f"region:  {task_data['genome_context']['region']}\n")

    results = []
    for strat in only:
        if strat not in STRATEGIES and strat != "win":
            print(f"skip unknown strategy: {strat}")
            continue
        try:
            r = _run_one(strat, task_data)
            results.append(r)
            print(
                f"  {r['strategy']:12s}  rev={r['rev'][-24:]:24s}  "
                f"sites={r['n_variants']:2d}  ann={r['n_annotations']:2d}  "
                f"truth={r['truth_hit']}"
            )
        except Exception as e:
            print(f"  {strat:12s}  FAILED: {e}")

    fps = {}
    for r in results:
        key = r["n_variants"]
        fps.setdefault(key, []).append(r["strategy"])
    if len(results) > 1:
        print("\nSite-count groups (strategies should differ after pipeline fix):")
        for n, strats in sorted(fps.items()):
            print(f"  n={n}: {strats}")


if __name__ == "__main__":
    main()
