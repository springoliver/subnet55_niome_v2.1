#!/usr/bin/env python3
"""5.22.02: your miners vs top — rev, GT, variant count."""
import json
import re
from collections import Counter

PATH = "Results/5.22.02/miner.json"
MY_UIDS = {50, 71, 139, 209, 235}


def parse_vcf(log):
    if "Miner VCF\n" not in log:
        return [], ""
    text = log.split("Miner VCF\n", 1)[1]
    rev = ""
    for line in text.splitlines():
        if line.startswith("##source="):
            rev = line.split("=", 1)[1].strip()
            break
    rows = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        p = line.split("\t")
        if len(p) < 10:
            continue
        gt = "."
        fmt, sample = p[8], p[9]
        if fmt == "GT":
            gt = sample.split(":")[0]
        elif "GT" in fmt:
            i = fmt.split(":").index("GT")
            gt = sample.split(":")[i].split(":")[0]
        af = None
        if "AF" in (p[7] if len(p) > 7 else ""):
            m = re.search(r"AF=([0-9.]+)", p[7])
            if m:
                af = float(m.group(1))
        if "AD" in fmt and ":" in sample:
            fi = fmt.split(":")
            sp = sample.split(":")
            if "AD" in fi and "DP" in fi:
                ads = [int(x) for x in sp[fi.index("AD")].split(",") if x != "."]
                dp = int(sp[fi.index("DP")])
                if len(ads) >= 2 and dp:
                    af = ads[1] / dp
        rows.append({
            "pos": int(p[1]),
            "ref": p[3],
            "alt": p[4].split(",")[0],
            "gt": gt,
            "af": af,
            "key": (int(p[1]), p[3], p[4].split(",")[0]),
        })
    return rows, rev


def main():
    with open(PATH) as f:
        data = json.load(f)

    by_uid = {}
    all_best = []
    for it in data["items"]:
        if it.get("validator_uid") != 154:
            continue
        uid = it["miner_uid"]
        rows, rev = parse_vcf(it.get("log", ""))
        if not rows:
            continue
        entry = {
            "uid": uid,
            "final": it["final_score"],
            "vcf": it["vcf_score"],
            "ann": it["annotation_score"],
            "prec": it["precision"],
            "rec": it["recall"],
            "n": len(rows),
            "rev": rev,
            "gt": Counter(r["gt"] for r in rows),
            "rows": rows,
        }
        all_best.append(entry)
        if uid in MY_UIDS:
            by_uid[uid] = entry

    all_best.sort(key=lambda x: -x["final"])
    top = all_best[0] if all_best else None

    print("5.22.02 task ddc5b57a  region chr7:117480000-117670000\n")
    print("=" * 70)
    print("YOUR MINERS (validator 154)")
    print("=" * 70)
    for uid in sorted(MY_UIDS):
        e = by_uid.get(uid)
        if not e:
            print(f"  UID {uid}: no entry on v154")
            continue
        hom_af = sum(1 for r in e["rows"] if r["af"] and r["af"] >= 0.58 and r["gt"] == "0/1")
        print(
            f"  UID {uid}: final={e['final']:.4f} vcf={e['vcf']:.4f} ann={e['ann']:.4f} "
            f"n={e['n']} prec={e['prec']:.3f} rec={e['rec']:.3f}"
        )
        print(f"    rev: {e['rev']}")
        print(f"    GT: {dict(e['gt'])}  | should-be-1/1 still 0/1 (AF>=0.58): {hom_af}")

    if top:
        print("\n" + "=" * 70)
        print(f"TOP miner UID {top['uid']}: final={top['final']:.4f} n={top['n']}")
        print(f"  rev: {top['rev'][:80]}")
        print(f"  GT: {dict(top['gt'])}")

    if top and by_uid.get(71):
        u = by_uid[71]
        t = top
        uk = {r["key"]: r for r in u["rows"]}
        tk = {r["key"]: r for r in t["rows"]}
        shared = set(uk) & set(tk)
        gt_mm = [(k, tk[k]["gt"], uk[k]["gt"], uk[k]["af"]) for k in shared if tk[k]["gt"] != uk[k]["gt"]]
        print("\n" + "=" * 70)
        print(f"UID 71 vs TOP UID {t['uid']}")
        print(f"  shared variants: {len(shared)}  only71: {len(set(uk)-set(tk))}  only_top: {len(set(tk)-set(uk))}")
        print(f"  GT mismatches on shared: {len(gt_mm)}")
        for k, tg, g, af in sorted(gt_mm)[:10]:
            print(f"    {k[0]} top={tg} uid71={g} af71={af}")

    print("\n" + "=" * 70)
    print("DEPLOY CHECK")
    print("=" * 70)
    print("  Current repo should log: niome-native-2026-05-21-v4")
    print("  Your 5.22.02 logs show v2/v3 only -> v4 NOT deployed on live hosts")


if __name__ == "__main__":
    main()
