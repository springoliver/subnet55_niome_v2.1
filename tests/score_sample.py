#!/usr/bin/env python3
"""
Score v3 miner pipeline against Results/sample (manager-provided fixture).

Requires: bwa, samtools, bcftools on Linux/WSL.

  export PYTHONPATH="$(pwd)"
  python setup_miner.py
  python tests/score_sample.py
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from niome_subnet.genomics.cftr_lookup import build_cftr_annotations
from niome_subnet.genomics.model import (
    GroundTruth,
    MinerSubmission,
    Task,
    TaskGenomeContext,
    TaskInput,
    TaskOutputSpec,
)
from niome_subnet.genomics.pipeline import ensure_hg38_chr7, run_pipeline
from niome_subnet.genomics.scoring import create_mapping_file, score

SAMPLE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "Results", "sample"
)


def load_sample_task() -> tuple[Task, GroundTruth]:
    truth_vcf = os.path.join(SAMPLE_DIR, "truth.vcf")
    ann_path = os.path.join(SAMPLE_DIR, "cftr2_annotations.json")
    r1 = os.path.join(SAMPLE_DIR, "read_1.fq")
    r2 = os.path.join(SAMPLE_DIR, "read_2.fq")

    truth_lines = [
        line
        for line in open(truth_vcf)
        if line.strip() and not line.startswith("#")
    ]
    positions = [int(line.split("\t")[1]) for line in truth_lines]
    pad = 8000
    region = f"chr7:{min(positions) - pad}-{max(positions) + pad}"

    task = Task(
        task_id="sample-v3",
        version="2.0",
        type="cftr_variant_calling",
        input=TaskInput(read1_fastq=r1, read2_fastq=r2),
        output_spec=TaskOutputSpec(
            format="vcf",
            required_fields=[
                "CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE",
            ],
        ),
        genome_context=TaskGenomeContext(
            chromosome="chr7",
            region=region,
            gene="CFTR",
        ),
        expected_variant_count=0,
    )

    ref = ensure_hg38_chr7()
    gt = GroundTruth(
        truth_vcf=truth_vcf,
        ref=ref,
        cftr2_annotations=ann_path,
    )
    return task, gt


def main():
    task, ground_truth = load_sample_task()
    print(f"Region: {task.genome_context.region}")
    print(f"expected_variant_count: {task.expected_variant_count}")

    with tempfile.TemporaryDirectory(prefix="niome_sample_") as work_dir:
        final_vcf, _ = run_pipeline(task, work_dir)
        with open(final_vcf) as fh:
            vcf_content = fh.read()

        annotations = build_cftr_annotations(final_vcf) or {}

        n_vcf = sum(
            1 for ln in vcf_content.splitlines() if ln and not ln.startswith("#")
        )
        print(f"Submitted variants: {n_vcf}")
        print(f"ClinVar annotations: {len(annotations)}")

        import shutil

        os.makedirs("data", exist_ok=True)
        shutil.copy(ground_truth.truth_vcf, "data/truth.vcf")
        shutil.copy(ground_truth.ref, "data/ref.fa")
        shutil.copy(ground_truth.cftr2_annotations, "data/cftr2_annotations.json")
        shutil.copy(task.input.read1_fastq, "data/read_1.fq")
        shutil.copy(task.input.read2_fastq, "data/read_2.fq")

        ground_truth.truth_vcf = "data/truth.vcf"
        ground_truth.ref = "data/ref.fa"
        ground_truth.cftr2_annotations = "data/cftr2_annotations.json"

        bam = create_mapping_file(ground_truth.ref, "data/read_1.fq", "data/read_2.fq")

        miner_score = score(
            MinerSubmission(
                uid=0,
                vcf_content=vcf_content,
                response_time=1.0,
                cftr_annotations=annotations,
            ),
            ground_truth,
            bam,
        )

        print(f"\n{'='*50}")
        print(f"VCF score:         {miner_score.vcf_score:.4f}")
        print(f"Annotation score:  {miner_score.annotation_score:.4f}")
        print(f"Final score:       {miner_score.final_score:.4f}")
        print(f"Precision:         {miner_score.precision:.4f}")
        print(f"Recall:            {miner_score.recall:.4f}")
        print(f"F1:                {miner_score.f1_score:.4f}")
        print(f"{'='*50}")
        return miner_score.final_score


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"FAILED: {e}")
        print("Run on Linux/WSL with bwa, samtools, bcftools installed.")
        sys.exit(1)
