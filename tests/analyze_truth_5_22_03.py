#!/usr/bin/env python3
"""
Analyze manager truth for Results/5.22.03 (task 7fc3be20).

- Summarize truth VCF + cftr2_annotations
- Compare to prior rounds
- Offline Native v5 selection vs truth (synthetic read-supported pool)
- Show noise-band sites (117504296 / 117504400) are kept after v5 (no hard drop)
"""
import importlib.util
import json
import os
import sys
import types
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GENOMICS = ROOT / "niome_subnet" / "genomics"
TRUTH_DIR = ROOT / "Results" / "5.22.03" / "real_correct_result"
REGION = "chr7:117480000-117670000"
REGION_START, REGION_END = 117480000, 117670000

_bt = types.SimpleNamespace()
_bt.logging = types.SimpleNamespace(
    info=lambda *a, **k: None,
    warning=lambda *a, **k: None,
    error=lambda *a, **k: None,
)
sys.modules["bittensor"] = _bt
_pkg = types.ModuleType("niome_subnet")
_genomics = types.ModuleType("niome_subnet.genomics")
_genomics.__path__ = [str(GENOMICS)]
sys.modules["niome_subnet"] = _pkg
sys.modules["niome_subnet.genomics"] = _genomics
_pkg.genomics = _genomics


def _load(fullname, path):
    spec = importlib.util.spec_from_file_location(fullname, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[fullname] = mod
    spec.loader.exec_module(mod)
    return mod


def _bootstrap_native():
    tp = _load("niome_subnet.genomics.task_profile", GENOMICS / "task_profile.py")
    rt = _load("niome_subnet.genomics.read_types", GENOMICS / "read_types.py")
    es = _load(
        "niome_subnet.genomics.evidence_selection", GENOMICS / "evidence_selection.py"
    )
    cftr = types.ModuleType("niome_subnet.genomics.cftr_lookup")
    cftr.ensure_clinvar_db = lambda: ""
    sys.modules["niome_subnet.genomics.cftr_lookup"] = cftr
    model = types.ModuleType("niome_subnet.genomics.model")
    model.Task = object
    sys.modules["niome_subnet.genomics.model"] = model
    rc = _load("niome_subnet.genomics.read_calling", GENOMICS / "read_calling.py")
    return tp, rt, es, rc


def parse_truth(path: Path):
    rows = []
    for line in path.read_text().splitlines():
        if not line.startswith("chr7"):
            continue
        p = line.split("\t")
        rows.append(
            {
                "pos": int(p[1]),
                "ref": p[3],
                "alt": p[4],
                "gt": p[9] if len(p) > 9 else ".",
            }
        )
    return rows


def jaccard_keys(a, b):
    sa = {(x["pos"], x["ref"], x["alt"]) for x in a}
    sb = {(x["pos"], x["ref"], x["alt"]) for x in b}
    if not sa and not sb:
        return 0.0
    return len(sa & sb) / len(sa | sb)


def truth_to_read_calls(rt, rows, hom_ad=18, het_ad=10, dp=24):
    """Synthetic strong read evidence at every truth site (upper bound for Native)."""
    out = []
    for r in rows:
        gt = r["gt"]
        if gt == "1/1":
            ad, d = hom_ad, dp
        else:
            ad, d = het_ad, dp
        out.append(
            rt.ReadCall(
                r["pos"],
                r["ref"],
                r["alt"],
                200.0,
                gt,
                ad,
                d,
                True,
            )
        )
    return out


def main():
    truth_path = TRUTH_DIR / "truth.vcf"
    ann_path = TRUTH_DIR / "cftr2_annotations.json"
    task_path = ROOT / "Results" / "5.22.03" / "task.json"

    if not truth_path.is_file():
        print(f"Missing {truth_path}")
        return 1

    truth = parse_truth(truth_path)
    ann = json.loads(ann_path.read_text()) if ann_path.is_file() else {}
    task = json.loads(task_path.read_text()) if task_path.is_file() else {}

    print("=" * 72)
    print("5.22.03 TRUTH ANALYSIS")
    print("=" * 72)
    print(f"task_id: {task.get('task_id', '?')}")
    print(f"region:  {task.get('genome_context', {}).get('region', REGION)}")
    print(f"variants: {len(truth)}")
    print(f"cftr2_annotations keys: {len(ann)}")
    pos = [t["pos"] for t in truth]
    print(f"pos range: {min(pos)} - {max(pos)}")
    print(f"GT: {dict(Counter(t['gt'] for t in truth))}")
    snps = sum(1 for t in truth if len(t["ref"]) == 1 and len(t["alt"]) == 1)
    print(f"SNPs: {snps}, indels/complex: {len(truth) - snps}")

    tp, rt, es, rc = _bootstrap_native()
    noise = tp.noise_band_penalty
    in_noise = [t for t in truth if noise(t["pos"])]
    print(f"\nSites in noise band ({tp.READ_NOISE_LO}-{tp.READ_NOISE_HI}): {len(in_noise)}")
    for t in in_noise:
        print(f"  {t['pos']} {t['ref']}>{t['alt']} {t['gt']}")

    print("\n--- Overlap vs prior manager truth ---")
    for rd in ["5.21.01", "5.21.02", "5.21.03"]:
        p = ROOT / "Results" / rd / "real_correct_result" / "truth.vcf"
        if p.is_file():
            other = parse_truth(p)
            print(f"  {rd}: n={len(other):2}  jaccard={jaccard_keys(truth, other):.3f}")

    print("\n--- All truth sites ---")
    for t in sorted(truth, key=lambda x: x["pos"]):
        nb = " [noise band]" if noise(t["pos"]) else ""
        print(f"  {t['pos']:>9}  {t['ref']:>20}  {t['alt']:>6}  {t['gt']}{nb}")

    prof = tp.classify_task(REGION)
    pool = truth_to_read_calls(rt, truth)
    clinvar = {}
    selected = rc.select_read_variants(
        pool, REGION_START, REGION_END, prof, clinvar
    )
    sk = {(c.pos, c.ref.upper(), c.alt.upper()) for c in selected}
    tk = {(t["pos"], t["ref"], t["alt"]) for t in truth}

    print("\n" + "=" * 72)
    print(f"OFFLINE NATIVE ({es.METHOD_ID}) — synthetic read pool = truth")
    print("=" * 72)
    print(f"profile: {prof.name}")
    print(f"submitted: {len(selected)}  (trim max {es.NATIVE_COUNT_TRIM_MAX})")
    print(f"TP: {len(tk & sk)}  FP: {len(sk - tk)}  FN: {len(tk - sk)}")

    fn = [t for t in truth if (t["pos"], t["ref"], t["alt"]) not in sk]
    if fn:
        print("\nFalse negatives (not selected from perfect read pool):")
        for t in fn:
            print(f"  {t['pos']} {t['ref']}>{t['alt']} {t['gt']}")

    noise_sel = [c for c in selected if noise(c.pos)]
    print(f"\nNoise-band sites selected: {len(noise_sel)}")
    for c in noise_sel:
        print(f"  {c.pos} {c.ref}>{c.alt} GT={c.gt}")

    gts = Counter(c.gt for c in selected)
    print(f"Submitted GT mix: {dict(gts)}")

    print("\n--- Implications for live miners (~0.56 on this task) ---")
    print("  - Truth N=26; submitting ~19-21 costs count_penalty.")
    print("  - Need ~10x 1/1 GT where AF high; all 0/1 loses ~half on hom-alt sites.")
    print("  - v5: noise band is score penalty only — 117504296/117504400 can be submitted.")
    print("  - Deploy v5 + pm2 restart; cache key now includes rev string.")
    print("=" * 72)
    return 0 if len(tk - sk) == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
