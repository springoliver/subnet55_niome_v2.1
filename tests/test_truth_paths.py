#!/usr/bin/env python3
import json
import os
import tempfile
from pathlib import Path

from tests.conftest import load_genomics_module

_truth = load_genomics_module("truth_paths.py")
find_task_truth = _truth.find_task_truth


def test_find_task_truth_in_truth_dir():
    with tempfile.TemporaryDirectory() as tmp:
        tid = "abc-123"
        base = Path(tmp) / tid
        base.mkdir()
        (base / "truth.vcf").write_text(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
            encoding="utf-8",
        )
        os.environ["NIOME_TRUTH_VCF"] = ""
        os.environ["NIOME_TRUTH_DIR"] = tmp
        found = find_task_truth(tid)
        assert found is not None
        assert found[0].endswith("truth.vcf")


def test_find_task_truth_in_results_root():
    with tempfile.TemporaryDirectory() as tmp:
        round_dir = Path(tmp) / "5.99.01"
        round_dir.mkdir()
        tid = "task-xyz"
        (round_dir / "task.json").write_text(
            json.dumps({"task_id": tid}),
            encoding="utf-8",
        )
        rc = round_dir / "real_correct_result"
        rc.mkdir()
        (rc / "truth.vcf").write_text(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
            encoding="utf-8",
        )
        os.environ.pop("NIOME_TRUTH_VCF", None)
        os.environ.pop("NIOME_TRUTH_DIR", None)
        os.environ["NIOME_RESULTS_ROOT"] = tmp
        found = find_task_truth(tid)
        assert found is not None
        assert "real_correct_result" in found[0]
