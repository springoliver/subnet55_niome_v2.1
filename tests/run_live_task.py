#!/usr/bin/env python3
"""
Run v3 pipeline on Results/live_task/task.json (ultra-wide live problem).

  export PYTHONPATH="$(pwd)"
  python setup_miner.py
  python tests/run_live_task.py
"""

import json
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from niome_subnet.genomics.cftr_lookup import build_cftr_annotations
from niome_subnet.genomics.model import Task
from niome_subnet.genomics.pipeline import run_pipeline
from niome_subnet.genomics.task_profile import classify_task, region_length

TASK_JSON = os.path.join(
    os.path.dirname(__file__), "..", "Results", "live_task", "task.json"
)


def main():
    with open(TASK_JSON) as fh:
        data = json.load(fh)

    task = Task(**data)
    region = task.genome_context.region
    profile = classify_task(region, task.expected_variant_count)

    print(f"task_id: {task.task_id}")
    print(f"region: {region} ({region_length(region):,} bp)")
    print(f"expected_variant_count: {task.expected_variant_count}")
    print(f"profile: {profile.name}")
    print(f"mpileup: {profile.mpileup_qual}")
    print(f"filters: dp>={profile.min_dp} af>={profile.min_af} qual>={profile.min_qual}")

    with tempfile.TemporaryDirectory(prefix="niome_live_") as work_dir:
        final_vcf, _ = run_pipeline(task, work_dir)
        with open(final_vcf) as fh:
            vcf = fh.read()

        annotations = build_cftr_annotations(final_vcf) or {}

        out_vcf = os.path.join(work_dir, "submitted.vcf")
        with open(out_vcf, "w") as fh:
            fh.write(vcf)

    lines = [ln for ln in vcf.splitlines() if ln and not ln.startswith("#")]
    print(f"\nSubmitted variants: {len(lines)}")
    print(f"ClinVar annotations: {len(annotations)}")
    if lines:
        positions = [int(ln.split("\t")[1]) for ln in lines]
        print(f"POS range: {min(positions)} – {max(positions)}")
        print("\nFirst 15 variants:")
        for ln in lines[:15]:
            p = ln.split("\t")
            print(f"  {p[1]}  {p[3]}>{p[4]}  GT={p[9].split(':')[0] if len(p)>9 else '?'}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"FAILED: {e}")
        sys.exit(1)
