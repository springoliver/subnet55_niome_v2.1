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
    METHOD_ID = "niome-native-2026-05-23-v6"
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
    print(f"trim_max:    {NATIVE_COUNT_TRIM_MAX}")
    print(f"mpileup:     {prof.mpileup_qual}")
    print(f"max_indel:   {prof.max_indel_len} bp")

    if rlen >= ULTRA_WIDE_MIN:
        print("\nUltra-wide CFTR checklist:")
        print("  - BWA GRCh38 chr7, padded region +/- 5kb")
        print("  - bcftools mpileup/call/norm (--indels-2.0); native_select_variants (v6)")
        print("  - No fixed variant count; expect ~11-29 sites from reads")
        print("  - GT: 1/1 if alt_AD/DP >= 0.58; 0/1 if >= 0.18")
        print("  - Noise band 117504200-400: penalty only (do not hard-drop)")
        print("  - Indels up to 48 bp; long CFTR deletions allowed")
        print("  - ClinVar annotate only submitted variants")
        print("  - Do NOT copy prior round VCF panels")

    ref_task = ROOT / "Results" / "5.22.03" / "task.json"
    if ref_task.is_file():
        ref = json.loads(ref_task.read_text())
        same_key = (
            ref["input"]["read1_fastq"].split("?")[0]
            == task["input"]["read1_fastq"].split("?")[0]
        )
        same_id = ref["task_id"] == task["task_id"]
        print("\nCompare to 5.22.03 (7fc3be20):")
        print(f"  Same S3 object key (crt/reads_*.fq): {same_key}")
        print(f"  Same task_id: {same_id}")
        if same_key and not same_id:
            print("  -> New task_id but same bucket path: reads may be updated.")
            print("     Manager truth for 5.22.03 had 26 variants — use as shape hint only.")

    print("\nDeploy verify:")
    print("  pm2 logs -> rev=niome-native-2026-05-23-v6")
    print("  first solve: done in tens of seconds, variants ~20-28")
    print("  repeat same task: may show 0.0s (cache includes rev)")
    print("=" * 72)


if __name__ == "__main__":
    main()
