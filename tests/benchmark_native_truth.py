#!/usr/bin/env python3
"""Offline: NIOME Native selection vs manager truth (5.21.01–03)."""
import importlib.util
import os
import sys
import types

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GENOMICS = os.path.join(ROOT, "niome_subnet", "genomics")
RESULTS = os.path.join(ROOT, "Results")

# Offline truth benchmarks: no curriculum fill (live ultra_wide uses target 30).
os.environ["NIOME_CURRICULUM_TARGET"] = "0"

_bt = types.SimpleNamespace()
_bt.logging = types.SimpleNamespace(
    info=lambda *a, **k: None,
    warning=lambda *a, **k: None,
    error=lambda *a, **k: None,
)
sys.modules["bittensor"] = _bt
_pkg = types.ModuleType("niome_subnet")
_genomics = types.ModuleType("niome_subnet.genomics")
_genomics.__path__ = [GENOMICS]
sys.modules["niome_subnet"] = _pkg
sys.modules["niome_subnet.genomics"] = _genomics
_pkg.genomics = _genomics


def _load(fullname, path):
    spec = importlib.util.spec_from_file_location(fullname, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[fullname] = mod
    spec.loader.exec_module(mod)
    return mod


tp = _load("niome_subnet.genomics.task_profile", os.path.join(GENOMICS, "task_profile.py"))
rt = _load("niome_subnet.genomics.read_types", os.path.join(GENOMICS, "read_types.py"))
_load("niome_subnet.genomics.evidence_selection", os.path.join(GENOMICS, "evidence_selection.py"))
cftr = types.ModuleType("niome_subnet.genomics.cftr_lookup")
cftr.ensure_clinvar_db = lambda: ""
sys.modules["niome_subnet.genomics.cftr_lookup"] = cftr
model = types.ModuleType("niome_subnet.genomics.model")
model.Task = object
sys.modules["niome_subnet.genomics.model"] = model

rc = _load("niome_subnet.genomics.read_calling", os.path.join(GENOMICS, "read_calling.py"))

ReadCall = rt.ReadCall
select_read_variants = rc.select_read_variants
classify_task = tp.classify_task
from niome_subnet.genomics.evidence_selection import (  # noqa: E402
    NATIVE_COUNT_TRIM_MAX,
    METHOD_ID,
)

# Typical FPs seen in ultra2 logs (not in manager truth)
FP_501 = [
    (117540258, "GC", "G", 222.0, 30, 16),
    (117559614, "TATAG", "T", 172.9, 34, 10),
    (117602847, "C", "CTG", 117.8, 28, 7),
]
FP_502 = [
    (117508505, "T", "A", 131.3, 17, 7),
    (117509127, "C", "CT", 222.0, 34, 17),
]
FP_503 = [
    (117509041, "G", "T", 97.0, 22, 6),
    (117510583, "AG", "A", 149.0, 26, 9),
    (117504350, "A", "G", 80.0, 8, 12),
]


def truth_calls(path):
    out = []
    for line in open(path):
        if line.strip() and not line.startswith("#"):
            p = line.split("\t")
            out.append(
                ReadCall(
                    int(p[1]),
                    p[3],
                    p[4],
                    200.0,
                    p[9].split(":")[0] if len(p) > 9 else "0/1",
                    12,
                    24,
                    True,
                )
            )
    return out


def fp_calls(fps):
    return [
        ReadCall(p, r, a, q, "0/1", ad, dp, True)
        for p, r, a, q, dp, ad in fps
    ]


def run_round(name, fps):
    truth_path = os.path.join(RESULTS, name, "real_correct_result", "truth.vcf")
    truth = truth_calls(truth_path)
    prof = classify_task("chr7:117480000-117670000")
    sel = select_read_variants(
        truth + fp_calls(fps), 117480000, 117670000, prof, {}
    )
    tk = {(c.pos, c.ref, c.alt) for c in truth}
    sk = {(c.pos, c.ref, c.alt) for c in sel}
    return {
        "round": name,
        "truth_n": len(truth),
        "submitted": len(sel),
        "tp": len(tk & sk),
        "fp": len(sk - tk),
        "fn": len(tk - sk),
    }


def main():
    print(f"method={METHOD_ID}")
    print(f"trim_max={NATIVE_COUNT_TRIM_MAX} (no minimum submit count)\n")
    ok = True
    # 5.22.03: typical ultra-wide FPs near (not overlapping) truth indels
    FP_52203 = [
        (117504350, "A", "G", 80.0, 8, 12),
        (117509041, "G", "T", 97.0, 22, 6),
        (117510583, "AG", "A", 149.0, 26, 9),
    ]

    for name, fps in [
        ("5.21.01", FP_501),
        ("5.21.02", FP_502),
        ("5.21.03", FP_503),
        ("5.22.03", FP_52203),
    ]:
        r = run_round(name, fps)
        band = 1 <= r["submitted"] <= NATIVE_COUNT_TRIM_MAX
        full = r["fn"] == 0
        low_fp = r["fp"] <= max(3, r["truth_n"] // 3)
        print(
            f"{r['round']}: truth={r['truth_n']} submitted={r['submitted']} "
            f"TP={r['tp']} FP={r['fp']} FN={r['fn']}  "
            f"recall_ok={full} under_trim_max={band} fp_ok={low_fp}"
        )
        if not full or not band or not low_fp:
            ok = False
    print("\n" + ("PASS (offline truth+FP pool)" if ok else "FAIL — tune gates or pipeline"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
