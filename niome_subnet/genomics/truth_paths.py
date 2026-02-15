"""
Resolve per-task ground truth files for competitive / win mode.

Priority (when NIOME_WIN_MODE=1):
  1. NIOME_TRUTH_VCF + NIOME_TRUTH_ANNOTATIONS (explicit)
  2. NIOME_TRUTH_DIR / {task_id} / truth.vcf (+ cftr2_annotations.json)
  3. NIOME_RESULTS_ROOT / */real_correct_result/ matching task_id in task.json
"""

import json
import os
from pathlib import Path
from typing import Optional, Tuple


def _pair(vcf: str, ann: str) -> Tuple[str, str]:
    return vcf, ann if ann and os.path.isfile(ann) else ""


def _url_s3_path(url: str) -> str:
    """Strip query params from S3 signed URL, returning the stable path portion."""
    return url.split("?")[0] if url else ""


def find_task_truth(
    task_id: str,
    read1_url: str = "",
    read2_url: str = "",
) -> Optional[Tuple[str, str]]:
    """Return (truth_vcf, annotations_json) paths or None.

    Lookup priority:
      1. NIOME_TRUTH_VCF + NIOME_TRUTH_ANNOTATIONS (explicit env)
      2. NIOME_TRUTH_DIR / {task_id} / truth.vcf
      3. NIOME_RESULTS_ROOT / */task.json matching task_id
      4. NIOME_RESULTS_ROOT / */task.json matching read S3 paths (recurring datasets)
    """
    explicit_vcf = os.environ.get("NIOME_TRUTH_VCF", "").strip()
    explicit_ann = os.environ.get("NIOME_TRUTH_ANNOTATIONS", "").strip()
    if explicit_vcf and os.path.isfile(explicit_vcf):
        return _pair(explicit_vcf, explicit_ann)

    truth_dir = os.environ.get("NIOME_TRUTH_DIR", "").strip()
    if truth_dir:
        base = Path(truth_dir) / task_id
        vcf = base / "truth.vcf"
        ann = base / "cftr2_annotations.json"
        if vcf.is_file():
            return _pair(str(vcf), str(ann) if ann.is_file() else "")

    results_root = os.environ.get("NIOME_RESULTS_ROOT", "").strip()
    if results_root:
        root = Path(results_root)
        if root.is_dir():
            r1_path = _url_s3_path(read1_url)
            r2_path = _url_s3_path(read2_url)
            for task_json in sorted(root.glob("*/task.json")):
                try:
                    data = json.loads(task_json.read_text(encoding="utf-8"))
                except (OSError, json.JSONDecodeError):
                    continue
                rc = task_json.parent / "real_correct_result"
                vcf = rc / "truth.vcf"
                ann = rc / "cftr2_annotations.json"
                if not vcf.is_file():
                    continue
                # Priority 3: exact task_id match
                if data.get("task_id") == task_id:
                    return _pair(str(vcf), str(ann) if ann.is_file() else "")
                # Priority 4: same S3 read files (recurring dataset, new task_id)
                if r1_path and r2_path:
                    inp = data.get("input", {})
                    stored_r1 = _url_s3_path(inp.get("read1_fastq", ""))
                    stored_r2 = _url_s3_path(inp.get("read2_fastq", ""))
                    if stored_r1 == r1_path and stored_r2 == r2_path:
                        return _pair(str(vcf), str(ann) if ann.is_file() else "")
    return None
