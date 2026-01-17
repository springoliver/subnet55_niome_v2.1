#!/usr/bin/env python3
"""Analyze miner rounds vs manager truth (real_correct_result)."""
import json
import os
import re
from collections import defaultdict

ROOT = os.path.join(os.path.dirname(__file__), "..", "Results")
MY_UIDS = {139, 235, 209}
VALIDATORS = (119, 154)


def parse_vcf_lines(text: str) -> list[dict]:
    rows = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        p = line.split("\t")
        if len(p) < 10:
            continue
        gt = "."
        fmt = p[8] if len(p) > 8 else ""
        sample = p[9] if len(p) > 9 else ""
        if fmt == "GT" or (":" not in fmt and sample):
            gt = sample.split(":")[0] if ":" in sample else sample
        elif "GT" in fmt:
            idx = fmt.split(":").index("GT") if "GT" in fmt.split(":") else 0
            gt = sample.split(":")[idx] if ":" in sample else sample
        rows.append({
            "pos": int(p[1]),
            "ref": p[3],
            "alt": p[4],
            "gt": gt,
            "key": (int(p[1]), p[3], p[4]),
        })
    return rows


def load_truth(round_dir: str) -> list[dict]:
    path = os.path.join(round_dir, "real_correct_result", "truth.vcf")
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        return parse_vcf_lines(f.read())


def extract_vcf_from_log(log: str) -> list[dict]:
    if "Miner VCF\n" not in log:
        return []
    return parse_vcf_lines(log.split("Miner VCF\n", 1)[1])


def compare(truth: list[dict], miner: list[dict]) -> dict:
    t_keys = {r["key"] for r in truth}
    m_keys = {r["key"] for r in miner}
    t_pos = {r["pos"] for r in truth}
    m_pos = {r["pos"] for r in miner}
    tp = t_keys & m_keys
    n_t, n_m = len(truth), len(miner)
    count_penalty = min(n_m, n_t) / max(n_m, n_t) if n_t and n_m else 0.0
    prec = len(tp) / n_m if n_m else 0.0
    rec = len(tp) / n_t if n_t else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    snps_t = sum(1 for r in truth if len(r["ref"]) == 1 and len(r["alt"]) == 1)
    return {
        "truth_n": n_t,
        "truth_snps": snps_t,
        "miner_n": n_m,
        "tp": len(tp),
        "fp": len(m_keys - t_keys),
        "fn": len(t_keys - m_keys),
        "pos_tp": len(t_pos & m_pos),
        "pos_fn": len(t_pos - m_pos),
        "pos_fp": len(m_pos - t_pos),
        "precision": prec,
        "recall": rec,
        "f1": f1,
        "count_penalty": count_penalty,
        "missed_pos": sorted(t_pos - m_pos),
        "extra_pos": sorted(m_pos - t_pos),
    }


def pick_validator_entry(by_uid: dict, uid: int) -> dict | None:
    entries = by_uid.get(uid, [])
    if not entries:
        return None
    for vid in VALIDATORS:
        for e in entries:
            if e.get("validator_uid") == vid:
                return e
    return entries[0]


def load_round(round_name: str) -> dict:
    with open(os.path.join(ROOT, round_name, "miner.json")) as f:
        data = json.load(f)
    items = data if isinstance(data, list) else data.get("items", data)

    by_uid: dict[int, list] = defaultdict(list)
    for it in items:
        if it.get("validator_uid") in VALIDATORS:
            by_uid[it["miner_uid"]].append(it)

    truth = load_truth(os.path.join(ROOT, round_name))
    all_scores = [e for es in by_uid.values() for e in es]
    top = sorted(all_scores, key=lambda x: x.get("final_score", 0), reverse=True)[:8]

    out = {
        "round": round_name,
        "has_truth": bool(truth),
        "truth_n": len(truth),
        "truth_snps": sum(
            1 for r in truth if len(r["ref"]) == 1 and len(r["alt"]) == 1
        ),
        "validators": sorted({e["validator_uid"] for e in all_scores}),
        "top8": [],
        "my_miners": {},
    }

    for t in top:
        rows = extract_vcf_from_log(t.get("log", ""))
        cmp = compare(truth, rows) if truth else {}
        out["top8"].append({
            "uid": t["miner_uid"],
            "val": t["validator_uid"],
            "final": t["final_score"],
            "vcf": t["vcf_score"],
            "annot": t["annotation_score"],
            "weight": t.get("weight", 0),
            "n": cmp.get("miner_n", len(rows)),
            "tp": cmp.get("tp", 0),
            "truth_n": cmp.get("truth_n", len(truth)),
        })

    for uid in MY_UIDS:
        e = pick_validator_entry(by_uid, uid)
        if not e:
            continue
        rows = extract_vcf_from_log(e.get("log", ""))
        cmp = compare(truth, rows) if truth else {}
        m = re.search(r"##source=niome_miner[_-](\S+)", e.get("log", ""))
        out["my_miners"][uid] = {
            "validator_uid": e["validator_uid"],
            "validator_final": e["final_score"],
            "validator_vcf": e["vcf_score"],
            "validator_annot": e["annotation_score"],
            "weight": e.get("weight", 0),
            "source": m.group(1) if m else "unknown",
            "vs_truth": cmp,
        }

    return out


def main():
    rounds = ["5.21.01", "5.21.02", "5.21.03"]
    print("=" * 72)
    print("NIOME ROUNDS vs MANAGER TRUTH — UIDs 139, 235, 209")
    print("Deploy target: rev 2026-05-21-v3-ultra4")
    print("=" * 72)

    for rname in rounds:
        if not os.path.isdir(os.path.join(ROOT, rname)):
            continue
        ro = load_round(rname)
        print(f"\n### {ro['round']}  (validators {ro['validators']})")
        if ro["has_truth"]:
            indels = ro["truth_n"] - ro["truth_snps"]
            print(
                f"Manager truth: {ro['truth_n']} variants "
                f"({ro['truth_snps']} SNPs, {indels} indels/complex)"
            )
        print("Top miners vs manager truth (exact REF/ALT):")
        for i, t in enumerate(ro["top8"][:5], 1):
            print(
                f"  {i}. UID {t['uid']:3d} val{t['val']}  final={t['final']:.4f}  "
                f"n={t['n']}  TP={t['tp']}/{t['truth_n']}  weight={t['weight']}"
            )
        for uid in sorted(ro["my_miners"]):
            m = ro["my_miners"][uid]
            c = m["vs_truth"]
            print(
                f"\n  YOUR UID {uid}  ({m['source']})  validator {m['validator_uid']}"
            )
            print(
                f"    live final={m['validator_final']:.4f}  "
                f"weight={m['weight']}"
            )
            if c:
                print(
                    f"    vs manager truth: n={c['miner_n']} truth={c['truth_n']}  "
                    f"TP={c['tp']} FP={c['fp']} FN={c['fn']}  "
                    f"count_pen~{c['count_penalty']:.3f}  F1~{c['f1']:.3f}"
                )
                if c["missed_pos"]:
                    print(f"    missed ({len(c['missed_pos'])}): {c['missed_pos'][:6]}...")
                if c["extra_pos"]:
                    print(f"    extra ({len(c['extra_pos'])}): {c['extra_pos'][:6]}...")


if __name__ == "__main__":
    main()
