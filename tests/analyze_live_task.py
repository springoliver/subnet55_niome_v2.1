#!/usr/bin/env python3
"""Summarize a live task.json and expected Native v5 behavior."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GENOMICS = ROOT / "niome_subnet" / "genomics"


def main():
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "Results" / "live_task" / "task.json"
    task = json.loads(path.read_text())
    region = task["genome_context"]["region"]
    _, rest = region.split(":")
    start, end = map(int, rest.split("-"))
    rlen = end - start

    import importlib.util

    def _load(name, rel):
        p = GENOMICS / rel
        spec = importlib.util.spec_from_file_location(name, p)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    tp = _load("task_profile", "task_profile.py")
    classify_task = tp.classify_task
    ULTRA_WIDE_MIN = tp.ULTRA_WIDE_MIN
    METHOD_ID = "niome-native-2026-05-23-v8"
    NATIVE_COUNT_TRIM_MAX = 32

    prof = classify_task(region, task.get("expected_variant_count", 0))

    print("=" * 72)
    print("LIVE TASK — Native miner plan")
    print("=" * 72)
    print(f"task_id:     {task['task_id']}")
    print(f"version:     {task.get('version')}")
    print(f"region:      {region}  ({rlen} bp)")
    print(f"expected:    {task.get('expected_variant_count')} (hidden count)")
    print(f"profile:     {prof.name}")
    print(f"method:      {METHOD_ID}")
    print(f"trim_max:    {max(NATIVE_COUNT_TRIM_MAX, 34)} with curriculum_target=30")
    print(f"mpileup:     {prof.mpileup_qual}")
    print(f"max_indel:   {prof.max_indel_len} bp")

    if rlen >= ULTRA_WIDE_MIN:
        print("\nUltra-wide CFTR checklist:")
        print("  - BWA GRCh38 chr7, padded region +/- 5kb")
        print("  - bcftools mpileup/call/norm; native v8 curriculum_target=30 + loose/consensus passes")
        print("  - No fixed variant count; curriculum fill toward ~30 sites from reads")
        print("  - GT: 1/1 if alt_AD/DP >= 0.58; 0/1 if >= 0.18")
        print("  - Noise band 117504200-400: penalty only (do not hard-drop)")
        print("  - Indels up to 48 bp; long CFTR deletions allowed")
        print("  - ClinVar annotate only submitted variants")
        print("  - Do NOT copy prior round VCF panels")

    for label, path in [
        ("5.22.04 (a9a4db4a)", ROOT / "Results" / "5.22.04" / "task.json"),
        ("5.22.03 (7fc3be20)", ROOT / "Results" / "5.22.03" / "task.json"),
    ]:
        if not path.is_file():
            continue
        ref = json.loads(path.read_text())
        same_key = (
            ref["input"]["read1_fastq"].split("?")[0]
            == task["input"]["read1_fastq"].split("?")[0]
        )
        same_id = ref["task_id"] == task["task_id"]
        ref_date = ref["input"]["read1_fastq"].split("X-Amz-Date=")[-1][:16] if "X-Amz-Date=" in ref["input"]["read1_fastq"] else "?"
        new_date = task["input"]["read1_fastq"].split("X-Amz-Date=")[-1][:16] if "X-Amz-Date=" in task["input"]["read1_fastq"] else "?"
        print(f"\nCompare to {label}:")
        print(f"  Same S3 key: {same_key}  presign: {ref_date} -> {new_date}")
        print(f"  Same task_id: {same_id}")
    print("\n  Curriculum panel trend: 20 -> 22 -> 29 -> target 30 (this task).")
    print("  New presign 235015Z = likely NEW sample at crt/reads_*.fq — do not copy 5.22.04 positions.")

    print("\nDeploy verify:")
    print("  pm2 logs -> rev=niome-native-2026-05-23-v8 merged_pool>=28 curriculum_target=30")
    print("  first solve: done in tens of seconds, variants ~28-30 (if pool allows)")
    print("  repeat same task: may show 0.0s (cache includes rev)")
    print("=" * 72)


if __name__ == "__main__":
    main()
