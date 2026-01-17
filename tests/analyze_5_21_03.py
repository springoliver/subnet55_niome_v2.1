#!/usr/bin/env python3
import json
import os

ROOT = os.path.join(os.path.dirname(__file__), "..", "Results", "5.21.03")
MY = {139, 235, 209}


def parse_vcf(text):
    rows = []
    for line in text.splitlines():
        if line and not line.startswith("#"):
            p = line.split("\t")
            rows.append((int(p[1]), p[3], p[4], p[9] if len(p) > 9 else "."))
    return rows


def from_log(log):
    if "Miner VCF\n" not in log:
        return []
    return parse_vcf(log.split("Miner VCF\n", 1)[1])


def compare(truth, miner):
    tk = set((a, b, c) for a, b, c, _ in truth)
    mk = set((a, b, c) for a, b, c, _ in miner)
    tp = tk & mk
    t_pos = {a for a, _, _, _ in truth}
    m_pos = {a for a, _, _, _ in miner}
    n_t, n_m = len(truth), len(miner)
    cp = min(n_m, n_t) / max(n_m, n_t) if n_t and n_m else 0
    return {
        "n": n_m,
        "tp": len(tp),
        "fp": len(mk - tk),
        "fn": len(tk - mk),
        "pos_tp": len(t_pos & m_pos),
        "cp": cp,
        "missed": sorted(a for a, _, _ in tk - mk),
        "extra": sorted(a for a, _, _ in mk - tk),
    }


truth = parse_vcf(open(os.path.join(ROOT, "real_correct_result", "truth.vcf")).read())
snps = sum(1 for p, r, a, _ in truth if len(r) == 1 and len(a) == 1)
print("MANAGER TRUTH 5.21.03")
print(f"  variants: {len(truth)}  (SNPs: {snps}, indels/complex: {len(truth)-snps})")
print(f"  positions: {min(p for p,_,_,_ in truth)} - {max(p for p,_,_,_ in truth)}")

with open(os.path.join(ROOT, "miner.json")) as f:
    items = json.load(f)["items"]

for vid in (119, 154):
    by = {it["miner_uid"]: it for it in items if it["validator_uid"] == vid}
    if not by:
        continue
    print(f"\nValidator {vid} — top 8:")
    for it in sorted(by.values(), key=lambda x: x["final_score"], reverse=True)[:8]:
        m = from_log(it.get("log", ""))
        c = compare(truth, m)
        print(
            f"  UID {it['miner_uid']:3d}  final={it['final_score']:.4f}  "
            f"vcf={it['vcf_score']:.4f}  annot={it['annotation_score']:.4f}  "
            f"n={c['n']}  TP={c['tp']}/{len(truth)}  FP={c['fp']}  FN={c['fn']}  "
            f"count_pen~{c['cp']:.3f}  weight={it.get('weight', 0)}"
        )

def pos_near(truth_pos, miner_pos, w=50):
    return sum(1 for tp in truth_pos if any(abs(tp - mp) <= w for mp in miner_pos))


t_pos = {a for a, _, _, _ in truth}
print("\nPosition overlap (±50bp) vs manager truth:")
for label, uid in [("winner", 21), ("you", 139)]:
    it = next((x for x in items if x["miner_uid"] == uid and x["validator_uid"] == 154), None)
    if it:
        m = from_log(it["log"])
        m_pos = {a for a, _, _, _ in m}
        print(f"  {label} UID {uid}: {pos_near(t_pos, m_pos)} / {len(t_pos)} truth positions near a miner call")

print("\nYOUR MINERS vs TRUTH:")
for uid in sorted(MY):
    for vid in (119, 154):
        it = next((x for x in items if x["miner_uid"] == uid and x["validator_uid"] == vid), None)
        if not it:
            continue
        m = from_log(it["log"])
        c = compare(truth, m)
        src = "ultra2" if "ultra2" in it.get("log", "") else "other"
        print(
            f"  UID {uid} val{vid}: final={it['final_score']:.4f}  n={c['n']}  "
            f"TP={c['tp']} FP={c['fp']} FN={c['fn']}  pen~{c['cp']:.3f}  ({src})"
        )
        if c["missed"][:8]:
            print(f"    missed sample: {c['missed'][:8]}...")
        if c["extra"][:8]:
            print(f"    extra sample: {c['extra'][:8]}...")
