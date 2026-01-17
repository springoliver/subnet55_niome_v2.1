#!/usr/bin/env python3
"""Simulate v3-ultra3 selection vs manager truth + historical FPs."""
import importlib.util
import os
import sys
import types

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GENOMICS = os.path.join(ROOT, "niome_subnet", "genomics")

_bt = types.SimpleNamespace()
_bt.logging = types.SimpleNamespace(info=lambda *a, **k: None, warning=lambda *a, **k: None, error=lambda *a, **k: None)
sys.modules["bittensor"] = _bt

_pkg = types.ModuleType("niome_subnet")
_genomics = types.ModuleType("niome_subnet.genomics")
_genomics.__path__ = [GENOMICS]
sys.modules["niome_subnet"] = _pkg
sys.modules["niome_subnet.genomics"] = _genomics
_pkg.genomics = _genomics


def _load_into(fullname: str, path: str):
    spec = importlib.util.spec_from_file_location(fullname, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[fullname] = mod
    spec.loader.exec_module(mod)
    return mod


tp = _load_into("niome_subnet.genomics.task_profile", os.path.join(GENOMICS, "task_profile.py"))
_genomics.task_profile = tp

cftr_stub = types.ModuleType("niome_subnet.genomics.cftr_lookup")
cftr_stub.ensure_clinvar_db = lambda: ""
sys.modules["niome_subnet.genomics.cftr_lookup"] = cftr_stub

model_stub = types.ModuleType("niome_subnet.genomics.model")
model_stub.Task = object
sys.modules["niome_subnet.genomics.model"] = model_stub

rc = _load_into("niome_subnet.genomics.read_calling", os.path.join(GENOMICS, "read_calling.py"))
_genomics.read_calling = rc

ReadCall = rc.ReadCall
select_read_variants = rc.select_read_variants
ULTRA_TARGET_MAX = tp.ULTRA_TARGET_MAX
classify_task = tp.classify_task

RESULTS = os.path.join(ROOT, "Results")

FP_501 = [
    (117540258, "GC", "G", 222.0, 30, 16),
    (117559614, "TATAG", "T", 172.9, 34, 10),
    (117602847, "C", "CTG", 117.8, 28, 7),
    (117606672, "AG", "A", 52.6, 26, 6),
    (117615581, "ca", "c", 64.8, 21, 5),
    (117617895, "T", "A", 104.2, 28, 7),
    (117642486, "G", "GC", 128.2, 29, 8),
]

FP_502 = [
    (117508505, "T", "A", 131.3, 17, 7),
    (117509127, "C", "CT", 222.0, 34, 17),
    (117540199, "TC", "T", 80.0, 20, 6),
    (117545137, "T", "C", 90.0, 22, 6),
    (117569408, "T", "G", 70.0, 18, 5),
    (117574451, "CT", "C", 60.0, 16, 4),
    (117592649, "G", "GA", 55.0, 15, 4),
    (117611762, "TG", "T", 50.0, 14, 4),
]


def load_truth_calls(path: str) -> list:
    calls = []
    for line in open(path):
        if not line.strip() or line.startswith("#"):
            continue
        p = line.split("\t")
        calls.append(
            ReadCall(
                pos=int(p[1]),
                ref=p[3],
                alt=p[4],
                qual=200.0,
                gt=p[9] if len(p) > 9 else "0/1",
                alt_ad=12,
                dp=24,
                pass_filter=True,
            )
        )
    return calls


def fp_to_calls(fps) -> list:
    return [
        ReadCall(pos=p, ref=r, alt=a, qual=q, gt="0/1", alt_ad=ad, dp=dp, pass_filter=True)
        for p, r, a, q, dp, ad in fps
    ]


def simulate(round_name: str, truth_n: int, fps) -> dict:
    truth_path = os.path.join(RESULTS, round_name, "real_correct_result", "truth.vcf")
    region_start, region_end = 117480000, 117670000
    profile = classify_task("chr7:117480000-117670000")

    truth = load_truth_calls(truth_path)
    selected = select_read_variants(truth + fp_to_calls(fps), region_start, region_end, profile, {})

    t_keys = {(c.pos, c.ref, c.alt) for c in truth}
    s_keys = {(c.pos, c.ref, c.alt) for c in selected}
    return {
        "round": round_name,
        "truth_n": truth_n,
        "submitted": len(selected),
        "tp": len(t_keys & s_keys),
        "fp": len(s_keys - t_keys),
        "fn": len(t_keys - s_keys),
        "in_band": 11 <= len(selected) <= ULTRA_TARGET_MAX,
        "missed": sorted(p for p, _, _ in (t_keys - s_keys)),
    }


def main():
    print(f"v3-ultra3 rev={rc.READ_CALLING_REV}  target count 11-{ULTRA_TARGET_MAX}\n")
    for r in (
        simulate("5.21.01", 13, FP_501),
        simulate("5.21.02", 12, FP_502),
    ):
        print(
            f"{r['round']}: submitted={r['submitted']} (truth={r['truth_n']})  "
            f"TP={r['tp']} FP={r['fp']} FN={r['fn']}  in_band={r['in_band']}"
        )
        if r["missed"]:
            print(f"  missed: {r['missed']}")
    print("\nNote: live rounds also need BAM calling to recover FN sites from reads.")


if __name__ == "__main__":
    main()
