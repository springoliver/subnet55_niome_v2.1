#!/usr/bin/env python3
"""Phase-1: analyze a live task.json (problem only, no truth yet)."""
import importlib.util
import json
import os
import sys
import types

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GENOMICS = os.path.join(ROOT, "niome_subnet", "genomics")

_bt = types.SimpleNamespace()
_bt.logging = types.SimpleNamespace(info=lambda *a, **k: None, warning=lambda *a, **k: None)
sys.modules["bittensor"] = _bt
_pkg = types.ModuleType("niome_subnet")
_genomics = types.ModuleType("niome_subnet.genomics")
_genomics.__path__ = [GENOMICS]
sys.modules["niome_subnet"] = _pkg
sys.modules["niome_subnet.genomics"] = _genomics


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


tp = _load("niome_subnet.genomics.task_profile", os.path.join(GENOMICS, "task_profile.py"))
es = _load("niome_subnet.genomics.evidence_selection", os.path.join(GENOMICS, "evidence_selection.py"))


def analyze_task(task: dict) -> dict:
    region = task["genome_context"]["region"]
    rlen = tp.region_length(region)
    prof = tp.classify_task(region, task.get("expected_variant_count", 0))
    return {
        "task_id": task["task_id"],
        "version": task.get("version"),
        "region": region,
        "region_bp": rlen,
        "expected_variant_count": task.get("expected_variant_count"),
        "profile": prof.name,
        "mpileup": prof.mpileup_qual,
        "max_indel_len": prof.max_indel_len,
        "dedupe_window": prof.dedupe_window,
        "method": es.METHOD_ID,
        "native_trim_max": es.NATIVE_COUNT_TRIM_MAX,
        "native_min_submit": "none (read evidence only)",
    }


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "Results", "5.21.04", "task.json")
    with open(path) as f:
        task = json.load(f)
    r = analyze_task(task)
    print("=== Live problem (phase 1) ===\n")
    for k, v in r.items():
        print(f"  {k}: {v}")
    same = r["region"] == "chr7:117480000-117670000" and r["profile"] == "ultra_wide"
    print()
    if same:
        print("VERDICT: Same ultra-wide CFTR as 5.21.01-03 -> use NIOME Native v3.")
        print("ACTION: Deploy v3 (no 10/16-20 count target). Log must show niome-native-2026-05-21-v3.")
    else:
        print("VERDICT: Region/profile differs from 5.21 rounds → review before submit.")
    print(f"\nSaved task path: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
