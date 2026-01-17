#!/usr/bin/env python3
"""Check if a coldkey/hotkey group is always top miner."""
import json
import os
from collections import defaultdict

ROOT = os.path.join(os.path.dirname(__file__), "..", "Results")

COLDKEY = "5FLELdekVzssX6C6i38KnAWAL17cgkCZsTjmgSphBdm2ae1M"
HOTKEYS = {
    "5CJrjjzWouD91ZnbNgHNCSs2cxqiTPMFktZ9kePzVW7D5bRE",
    "5DEx3gdcDJtx6nx8KUBRMiZ1pt6dPbqZXmrTHW5KX5DqngFf",
    "5ENxyJxk7Vised6BMCaEvwrwyngbwbRzMrvtcvm887ZryQ8p",
    "5HLCjx6HnPjKcxLxdbTkF2RTqQo8iY7mK4goX69AUYURy77B",
    "5DtCGawxYpgrr6fjSqn6Co6wQs7s2QTEWJb2Dfi2spyGo82g",
    "5GFAuKL9TbnRELcnnDsYCbn3ba1BPZRfgtXEfVyc5Cq7WUhP",
    "5FLiSE8TPsG5jiujnet4HZZTHBDQ6VjTr3Y8uo7aA7VfBSa8",
}
GROUP = HOTKEYS | {COLDKEY}


def scan_miner_json(path, round_name):
    if not os.path.isfile(path):
        return []
    with open(path) as f:
        data = json.load(f)
    items = data.get("items", data if isinstance(data, list) else [])
    return [(round_name, it) for it in items]


def main():
    all_items = []
    for base, sub in [(ROOT, os.listdir(ROOT)), (os.path.join(ROOT, "old_version_result"), os.listdir(os.path.join(ROOT, "old_version_result")))]:
        if not os.path.isdir(base):
            continue
        for name in sub:
            m = os.path.join(base, name, "miner.json")
            if os.path.isfile(m):
                era = "old" if "old_version" in base else "new"
                all_items.extend(scan_miner_json(m, f"{era}/{name}"))

    # uid <-> hotkey map
    uid_hotkey = defaultdict(set)
    group_uids = set()
    group_scores = []

    for rnd, it in all_items:
        hk = it.get("miner_hotkey", "")
        uid = it.get("miner_uid")
        if hk:
            uid_hotkey[uid].add(hk)
        if hk in HOTKEYS or hk == COLDKEY:
            group_uids.add(uid)
            group_scores.append({
                "round": rnd,
                "uid": uid,
                "hotkey": hk[:12] + "…",
                "validator": it.get("validator_uid"),
                "final": it.get("final_score", 0),
                "vcf": it.get("vcf_score", 0),
                "weight": it.get("weight", 0),
            })

    print("=" * 78)
    print("HOTKEY GROUP -> UID mapping")
    print("=" * 78)
    for uid in sorted(group_uids):
        hks = uid_hotkey[uid]
        in_group = [h for h in hks if h in HOTKEYS]
        print(f"  UID {uid}: {len(in_group)} hotkey(s) in your list")
        for h in sorted(in_group):
            print(f"    {h}")

    print("\n" + "=" * 78)
    print("GROUP MINER SCORES (when hotkey matches)")
    print("=" * 78)
    by_round = defaultdict(list)
    for g in group_scores:
        by_round[g["round"]].append(g)
    for rnd in sorted(by_round.keys()):
        entries = by_round[rnd]
        best = max(entries, key=lambda x: x["final"])
        print(f"\n  {rnd}:")
        for e in sorted(entries, key=lambda x: -x["final"])[:5]:
            print(f"    UID {e['uid']} v{e['validator']}: final={e['final']:.4f} vcf={e['vcf']:.4f}")

    print("\n" + "=" * 78)
    print("TOP MINER PER ROUND (validator 154, else 119, else 58)")
    print("=" * 78)
    by_round_all = defaultdict(list)
    for rnd, it in all_items:
        by_round_all[rnd].append(it)

    top1_group = 0
    top3_group = 0
    rounds_checked = 0
    for rnd in sorted(by_round_all.keys()):
        if "live_task" in rnd:
            continue
        for vid in (154, 119, 58):
            sub = [it for it in by_round_all[rnd] if it.get("validator_uid") == vid]
            if not sub:
                continue
            best = max(sub, key=lambda x: x.get("final_score", 0))
            rounds_checked += 1
            hk = best.get("miner_hotkey", "")
            uid = best.get("miner_uid")
            is_group = hk in HOTKEYS
            rank_note = "TOP-1 *** GROUP ***" if is_group else ""
            # top 3
            ranked = sorted(sub, key=lambda x: -x.get("final_score", 0))[:3]
            in_top3 = any(it.get("miner_hotkey") in HOTKEYS for it in ranked)
            if is_group:
                top1_group += 1
            if in_top3:
                top3_group += 1
            print(
                f"  {rnd} v{vid}: TOP UID {uid} score={best.get('final_score', 0):.4f} "
                f"hotkey={hk[:20]}… {rank_note}"
            )
            if in_top3 and not is_group:
                g_in = [it for it in ranked if it.get("miner_hotkey") in HOTKEYS]
                if g_in:
                    print(f"         group in top3: UID {g_in[0].get('miner_uid')} "
                          f"score={g_in[0].get('final_score', 0):.4f}")
            break

    print("\n" + "=" * 78)
    print("VERDICT")
    print("=" * 78)
    print(f"  Rounds with validator data: {rounds_checked}")
    print(f"  Group was #1 top scorer: {top1_group} times")
    print(f"  Group in top-3: {top3_group} times")
    print(f"  UID 55 hotkey 5DtCGawx… = top on 5.22.02 (score 1.0)")
    print(f"  UID 4 hotkey 5FLiSE8… appears in high scores on several rounds")
    if top1_group < rounds_checked // 2:
        print("  -> NOT always top. Strong on some rounds, zero or low on others.")


if __name__ == "__main__":
    main()
