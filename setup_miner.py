#!/usr/bin/env python3
"""
Pre-flight setup for NIOME miner (v3 read-all pipeline).

Downloads and caches:
  1. GRCh38 chr7 reference + BWA/samtools indexes
  2. ClinVar GRCh38 VCF filtered to CFTR region

Run once before starting the miner:
    python setup_miner.py
"""

import shutil
import subprocess
import sys

REQUIRED_TOOLS = ["bwa", "samtools", "bcftools"]
OPTIONAL_TOOLS = ["tabix", "bgzip"]


def check_tools() -> bool:
    print("=== Checking required system tools ===")
    all_ok = True
    for tool in REQUIRED_TOOLS:
        path = shutil.which(tool)
        if path:
            try:
                result = subprocess.run(
                    [tool, "--version"],
                    capture_output=True,
                    text=True,
                )
                version_line = (result.stdout or result.stderr).splitlines()[0]
            except Exception:
                version_line = "(version unknown)"
            print(f"  [OK] {tool:12s}  {version_line}")
        else:
            print(f"  [MISSING] {tool}  — install with: sudo apt-get install {tool}")
            all_ok = False
    print("\n=== Optional tools ===")
    for tool in OPTIONAL_TOOLS:
        path = shutil.which(tool)
        print(f"  [{'OK' if path else 'skip'}] {tool}")
    return all_ok


def warm_reference():
    print("\n=== Warming GRCh38 chr7 reference ===")
    from niome_subnet.genomics.pipeline import HG38_CHR7_PATH, ensure_hg38_chr7
    import os

    env_ref = os.environ.get("NIOME_HG38_REF", "").strip()
    if env_ref:
        print(f"  Using NIOME_HG38_REF={env_ref}")
        ref = ensure_hg38_chr7()
    elif os.path.exists(HG38_CHR7_PATH) and os.path.exists(HG38_CHR7_PATH + ".bwt"):
        print(f"  [SKIP] Already cached at {HG38_CHR7_PATH}")
        ref = HG38_CHR7_PATH
    else:
        ref = ensure_hg38_chr7()
    print(f"  [OK] Reference ready at {ref}")

    print("\n=== Warming CFTR slice fallback ===")
    from niome_subnet.genomics.pipeline import REF_PATH, ensure_reference

    if os.path.exists(REF_PATH) and os.path.exists(REF_PATH + ".bwt"):
        print(f"  [SKIP] Already cached at {REF_PATH}")
    else:
        ensure_reference()
        print(f"  [OK] Slice fallback at {REF_PATH}")


def warm_clinvar():
    print("\n=== Warming ClinVar CFTR region cache ===")
    import os

    from niome_subnet.genomics.cftr_lookup import CLINVAR_CFTR_VCF, ensure_clinvar_db

    if os.path.exists(CLINVAR_CFTR_VCF) and os.path.exists(CLINVAR_CFTR_VCF + ".tbi"):
        print(f"  [SKIP] Already cached at {CLINVAR_CFTR_VCF}")
    else:
        db = ensure_clinvar_db()
        print(f"  [OK] ClinVar CFTR DB ready at {db}")


def main():
    tools_ok = check_tools()
    if not tools_ok:
        print(
            "\n[ERROR] Install: sudo apt-get install bwa samtools bcftools\n"
            "        Miner requires Linux/WSL — not supported on native Windows."
        )
        sys.exit(1)

    try:
        warm_reference()
    except Exception as e:
        print(f"\n[ERROR] Reference setup failed: {e}")
        sys.exit(1)

    try:
        warm_clinvar()
    except Exception as e:
        print(f"\n[ERROR] ClinVar setup failed: {e}")
        sys.exit(1)

    print("\n=== Setup complete (v3 read-all miner) ===")
    print("    Run: python neurons/miner.py --netuid 55 --subtensor.network finney ...")
    print("    Test: python tests/score_sample.py")


if __name__ == "__main__":
    main()
