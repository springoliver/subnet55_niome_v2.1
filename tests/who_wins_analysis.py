#!/usr/bin/env python3
"""Did read-based / niome-native ever rank #1?"""
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "Results"
VALIDATORS = {119, 154}


def parse_vcf(log: str):
    m = re.search(r"Miner VCF\n(.*)", log, re.S)
    if not m:
        return 0, False, False
    text = m.group(1)
    n = sum(1 for l in text.splitlines() if l.startswith("chr7"))
    has_dp = "DP=" in text or ":DP:" in text or "\tGT:DP" in text
    has_native = "niome-native" in text or "niome_miner_niome-native" in text
    minimal = n > 0 and not has_dp and "bcftoolsVersion" not in text
    return n, has_dp, has_native or ("niome_miner_" in text and "native" in text)


def classify(log: str) -> str:
    if "niome-native" in log or "niome_miner_niome-native" in log:
        return "niome-native"
    if "bcftoolsVersion" in log or "niome-miner-v2" in log:
        return "bcftools-read"
    if "niome_miner_legacy" in log or ("niome_miner" in log and "bcftools" not in log):
        if "DP=" in log:
            return "legacy-read"
    if parse_vcf(log)[0] and "DP=" not in log:
        return "minimal-oracle"
    return "other"


def analyze_round(path: Path):
    data = json.loads(path.read_text(encoding="utf-8", errors="replace"))
    items = data.get("items", data if isinstance(data, list) else [])
    by_uid = {}
    for it in items:
        if it.get("validator_uid") not in VALIDATORS:
            continue
        uid = it["miner_uid"]
        fin = float(it.get("final_score", 0) or 0)
        if uid not in by_uid or fin > by_uid[uid][0]:
            by_uid[uid] = (fin, it.get("log", ""))

    if not by_uid:
        return None

    ranked = sorted(by_uid.items(), key=lambda x: -x[1][0])
    winner_uid, (winner_fin, winner_log) = ranked[0]
    winner_cls = classify(winner_log)

    # best per method class
    by_cls = {}
    for uid, (fin, log) in by_uid.items():
        c = classify(log)
        if c not in by_cls or fin > by_cls[c][0]:
            by_cls[c] = (fin, uid)

    return {
        "round": path.parent.name,
        "winner_uid": winner_uid,
        "winner_fin": winner_fin,
        "winner_cls": winner_cls,
        "best_native": by_cls.get("niome-native"),
        "best_bcftools": by_cls.get("bcftools-read"),
        "best_oracle": by_cls.get("minimal-oracle"),
        "n_miners": len(by_uid),
    }


def main():
    rounds = sorted(RESULTS.glob("*/miner.json"))
    print("=" * 78)
    print("WHO WINS #1? (validators 119, 154)")
    print("=" * 78)

    native_wins = 0
    read_wins = 0
    oracle_wins = 0

    for mp in rounds:
        if "old_version" in str(mp):
            continue
        r = analyze_round(mp)
        if not r:
            continue
        w = r["winner_cls"]
        if w == "niome-native":
            native_wins += 1
        elif w in ("bcftools-read", "legacy-read"):
            read_wins += 1
        elif w == "minimal-oracle":
            oracle_wins += 1

        bn = r["best_native"]
        bo = r["best_oracle"]
        native_str = (
            f"UID {bn[1]} {bn[0]:.4f}" if bn else "none"
        )
        print(
            f"\n{r['round']}: #1 UID {r['winner_uid']} {r['winner_fin']:.4f} [{r['winner_cls']}]"
        )
        print(f"         best native: {native_str}")
        if bo:
            gap = r["winner_fin"] - (bn[0] if bn else 0)
            print(
                f"         best oracle cluster: UID {bo[1]} {bo[0]:.4f}  "
                f"(native gap {r['winner_fin'] - bn[0]:.4f})" if bn else ""
            )

    print("\n" + "=" * 78)
    print(f"Rounds analyzed: {len([p for p in rounds if 'old_version' not in str(p)])}")
    print(f"#1 minimal-oracle: {oracle_wins}")
    print(f"#1 bcftools/legacy-read: {read_wins}")
    print(f"#1 niome-native: {native_wins}")
    print("=" * 78)


if __name__ == "__main__":
    main()
