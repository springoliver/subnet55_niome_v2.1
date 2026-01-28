#!/usr/bin/env python3
"""
Offline verify strategy env profiles against manager truth (where available).

Scores exported miner VCF from miner.json logs vs truth.vcf — sanity check only.
Full pipeline verify requires BAM + tests/score_sample.py on fixtures.

Run: python tests/verify_fleet_strategies.py
"""
import json
import os
import re
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from niome_subnet.genomics.model import GroundTruth, MinerSubmission
from niome_subnet.genomics.scoring import create_mapping_file, score

RESULTS = ROOT / "Results"
ROUNDS_WITH_TRUTH = ["5.21.01", "5.21.02", "5.21.03", "5.22.03"]


def extract_vcf(log: str) -> str:
    if "Miner VCF\n" not in log:
        return ""
    return log.split("Miner VCF\n", 1)[1]


def score_vcf_text(vcf_body: str, truth_vcf: str, ann_path: str, bam: str = "data/bam.bam"):
    if not vcf_body.strip():
        return None
    os.makedirs("data", exist_ok=True)
    with open("data/miner.vcf", "w") as fh:
        fh.write(vcf_body if vcf_body.startswith("##") else "##fileformat=VCFv4.2\n" + vcf_body)
    sub = MinerSubmission(vcf_content=open("data/miner.vcf").read(), response_time=10.0)
    gt = GroundTruth(truth_vcf=truth_vcf, ref="", cftr2_annotations=ann_path or "")
    create_mapping_file(gt.truth_vcf, gt.cftr2_annotations)
    return score(sub, gt, bam)


def main():
    print("Offline verification: top native VCF in log vs manager truth\n")
    for rd in ROUNDS_WITH_TRUTH:
        rc = RESULTS / rd / "real_correct_result"
        truth = rc / "truth.vcf"
        ann = rc / "cftr2_annotations.json"
        mj = RESULTS / rd / "miner.json"
        if not truth.is_file() or not mj.is_file():
            print(f"{rd}: skip (no truth or miner.json)")
            continue
        items = json.loads(mj.read_text(encoding="utf-8", errors="replace"))
        items = items if isinstance(items, list) else items.get("items", items)
        natives = []
        for it in items:
            log = it.get("log", "")
            if "niome-native" not in log and "niome_miner" not in log:
                continue
            if not re.search(r"^chr7\t", extract_vcf(log), re.M):
                continue
            natives.append(it)
        if not natives:
            print(f"{rd}: no native VCF in logs")
            continue
        best = max(natives, key=lambda x: x.get("final_score", 0))
        vcf = extract_vcf(best.get("log", ""))
        try:
            ms = score_vcf_text(vcf, str(truth), str(ann) if ann.is_file() else "")
            if ms:
                print(
                    f"{rd}: replay UID {best['miner_uid']} "
                    f"validator_final={best['final_score']:.4f} "
                    f"offline_final={ms.final_score:.4f} n_truth_lines ok"
                )
        except Exception as e:
            print(f"{rd}: score failed: {e}")

    print("\nStrategy module:")
    from niome_subnet.genomics.task_strategy import PROFILES, resolve_strategy

    class T:
        class genome_context:
            region = "chr7:117480000-117670000"

        class input:
            read1_fastq = "https://x/crt/reads_1.fq"
            read2_fastq = "https://x/crt/reads_2.fq"

    t = T()
    for name in ("v5_style", "high_recall", "v10", "win"):
        os.environ.pop("NIOME_TRUTH_VCF", None)
        s = resolve_strategy(t, miner_uid=71 if name == "v5_style" else None, truth_available=False)
        if name != "auto":
            os.environ["NIOME_STRATEGY"] = name
            s = resolve_strategy(t)
        print(f"  NIOME_STRATEGY={name} -> {s} profile={PROFILES.get(s, PROFILES['v10']).revision_tag}")


if __name__ == "__main__":
    main()
