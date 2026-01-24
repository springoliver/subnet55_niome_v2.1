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


def find_task_truth(task_id: str) -> Optional[Tuple[str, str]]:
    """Return (truth_vcf, annotations_json) paths or None."""
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
            for task_json in root.glob("*/task.json"):
                try:
                    data = json.loads(task_json.read_text(encoding="utf-8"))
                except (OSError, json.JSONDecodeError):
                    continue
                if data.get("task_id") != task_id:
                    continue
                rc = task_json.parent / "real_correct_result"
                vcf = rc / "truth.vcf"
                ann = rc / "cftr2_annotations.json"
                if vcf.is_file():
                    return _pair(
                        str(vcf),
                        str(ann) if ann.is_file() else "",
                    )
    return None
