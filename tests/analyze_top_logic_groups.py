#!/usr/bin/env python3
"""
Cluster top miners by submitted VCF panel (logic family) across NIOME rounds.
Checks hypothesis: only two coldkey/logic groups win (e.g. UID 44/16 vs 55/4).
"""
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "Results"
VALIDATORS = {119, 154}
TOP_N = 15
TOP_FINAL_MIN = 0.55

# Known coldkey-group UIDs from prior analysis (hypothesis)
GROUP_A_HINT = {1, 4, 8, 16, 44, 55, 80, 81}  # often cluster A
GROUP_B_HINT = {55, 4, 1, 8}  # user said 55, 4 top


def parse_vcf_from_log(log: str):
    m = re.search(r"Miner VCF\n(.*)", log, re.DOTALL)
    if not m:
        return []
    rows = []
    for line in m.group(1).splitlines():
        if not line.startswith("chr7"):
            continue
        p = line.split("\t")
        if len(p) < 5:
            continue
        gt = "."
        if len(p) > 9:
            fmt = p[8].split(":")
            samp = p[9].split(":")
            if "GT" in fmt:
                gt = samp[fmt.index("GT")].split("/")[0] + "/" + samp[fmt.index("GT")].split("/")[-1]
            elif fmt == ["GT"]:
                gt = samp[0]
            elif ":" not in fmt and len(samp) == 1:
                gt = samp[0]
        rows.append(
            {
                "pos": int(p[1]),
                "ref": p[3],
                "alt": p[4],
                "gt": gt,
                "id": p[2] if len(p) > 2 else ".",
            }
        )
    return rows


def panel_key(rows, with_gt=False):
    if with_gt:
        return tuple(sorted((r["pos"], r["ref"], r["alt"], r["gt"]) for r in rows))
    return tuple(sorted((r["pos"], r["ref"], r["alt"]) for r in rows))


def jaccard_keys(a, b):
    sa = {(r["pos"], r["ref"], r["alt"]) for r in a}
    sb = {(r["pos"], r["ref"], r["alt"]) for r in b}
    if not sa | sb:
        return 0.0
    return len(sa & sb) / len(sa | sb)


def detect_source(log: str) -> str:
    if "niome-native" in log:
        m = re.search(r"niome-native-[\d-]+-v\d+", log)
        return m.group(0) if m else "niome-native"
    if "niome-miner-v2" in log or "niome-miner-v2.1" in log:
        return "niome-miner-v2.1"
    if "bcftoolsVersion" in log:
        return "bcftools-oracle"
    if "##source=niome_miner" in log:
        return "niome_miner_legacy"
    if len(parse_vcf_from_log(log)) and "GT\t0/1" in log and "DP=" not in log:
        return "minimal-gt-panel"
    return "unknown"


def load_round(round_dir: Path):
    mj = round_dir / "miner.json"
    if not mj.is_file():
        return None
    data = json.loads(mj.read_text(encoding="utf-8", errors="replace"))
    items = data.get("items", data if isinstance(data, list) else [])
    by_uid = {}
    for it in items:
        if it.get("validator_uid") not in VALIDATORS:
            continue
        uid = it["miner_uid"]
        fin = float(it.get("final_score", 0) or 0)
        if uid not in by_uid or fin > by_uid[uid]["final"]:
            by_uid[uid] = {
                "final": fin,
                "vcf": float(it.get("vcf_score", 0) or 0),
                "hotkey": (it.get("miner_hotkey") or "")[:12],
                "log": it.get("log", ""),
            }
    return by_uid


def cluster_panels(by_uid, min_final=0.0):
    """Group UIDs by identical position+allele panel."""
    panels = defaultdict(list)
    for uid, s in by_uid.items():
        if s["final"] < min_final:
            continue
        rows = parse_vcf_from_log(s["log"])
        if not rows:
            continue
        pk = panel_key(rows, with_gt=False)
        panels[pk].append((uid, s, rows))
    return panels


def main():
    round_dirs = sorted(
        [
            p
            for p in RESULTS.iterdir()
            if p.is_dir() and (p / "miner.json").is_file()
        ],
        key=lambda p: p.name,
    )

    print("=" * 78)
    print("TOP MINER LOGIC FAMILIES (VCF panel clustering)")
    print("=" * 78)

    hotkey_to_uids = defaultdict(set)
    uid_to_hotkey = {}
    family_wins = Counter()  # hotkey prefix -> rounds as #1 panel
    round_summaries = []

    for rd in round_dirs:
        by_uid = load_round(rd)
        if not by_uid:
            continue

        ranked = sorted(by_uid.items(), key=lambda x: -x[1]["final"])
        top = [(u, s) for u, s in ranked if s["final"] >= TOP_FINAL_MIN][:TOP_N]
        if not top:
            top = ranked[:TOP_N]

        panels = cluster_panels(by_uid, min_final=TOP_FINAL_MIN)
        if not panels:
            panels = cluster_panels(by_uid, min_final=0.0)

        # Sort panel clusters by: max member final score, then size
        panel_list = []
        for pk, members in panels.items():
            best_fin = max(m[1]["final"] for m in members)
            panel_list.append((best_fin, len(members), pk, members))
        panel_list.sort(key=lambda x: (-x[0], -x[1]))

        # Only panels that appear in top scorers
        top_uids = {u for u, _ in top}
        winning_panels = []
        for best_fin, size, pk, members in panel_list:
            top_in = [m for m in members if m[0] in top_uids]
            if top_in:
                winning_panels.append((best_fin, size, pk, members, top_in))

        print(f"\n### {rd.name}  ({len(by_uid)} miners, {len(winning_panels)} top panels)")
        if ranked:
            print(f"    #1 UID {ranked[0][0]} final={ranked[0][1]['final']:.4f} hk={ranked[0][1]['hotkey']}…")

        for i, (best_fin, size, pk, members, top_in) in enumerate(winning_panels[:4], 1):
            uids = sorted(m[0] for m in top_in)
            sample = top_in[0]
            rows = sample[2]
            src = detect_source(sample[1]["log"])
            indels = sum(1 for r in rows if len(r["ref"]) > 1 or len(r["alt"]) > 1)
            hom = sum(1 for r in rows if r["gt"] == "1/1")
            has_cv = sum(1 for r in rows if r["id"] not in (".", ""))
            print(
                f"    Panel {i}: n={len(pk)} sites | {len(top_in)} in top-{TOP_N} | "
                f"best_final={best_fin:.4f} | src={src} | indels={indels} hom={hom} clinvar_ids={has_cv}"
            )
            print(f"           UIDs: {uids[:20]}{'…' if len(uids) > 20 else ''}")
            for uid, s, _ in top_in[:3]:
                hotkey_to_uids[s["hotkey"]].add(uid)
                uid_to_hotkey[uid] = s["hotkey"]

        # Cross-panel similarity among top panels
        if len(winning_panels) >= 2:
            r0 = winning_panels[0][4][0][2]
            r1 = winning_panels[1][4][0][2]
            j = jaccard_keys(r0, r1)
            print(f"    Panel1 vs Panel2 Jaccard: {j:.3f}")

        # Is #1 unique logic?
        if winning_panels:
            leader_uids = sorted(m[0] for m in winning_panels[0][4])
            family_wins[panel_key(winning_panels[0][4][0][2])] += 1
            round_summaries.append(
                {
                    "round": rd.name,
                    "leader_uid": ranked[0][0],
                    "leader_hk": ranked[0][1]["hotkey"],
                    "n_panels_in_top": len(winning_panels),
                    "panel1_uids": leader_uids[:12],
                    "panel1_n": len(winning_panels[0][2]),
                }
            )

    print("\n" + "=" * 78)
    print("HOTKEY <-> UID (from top-panel miners)")
    print("=" * 78)
    # Group UIDs that share hotkey prefix (first 8 chars often same coldkey group)
    hk_groups = defaultdict(set)
    for hk, uids in sorted(hotkey_to_uids.items(), key=lambda x: -len(x[1])):
        if not hk:
            continue
        prefix = hk[:8]
        hk_groups[prefix].update(uids)
    for prefix, uids in sorted(hk_groups.items(), key=lambda x: -len(x[1]))[:12]:
        print(f"  hotkey {prefix}…  UIDs={sorted(uids)}")

    print("\n" + "=" * 78)
    print("HYPOTHESIS CHECK: UID 44/16 vs UID 55/4")
    print("=" * 78)
    for rd in round_summaries:
        p1 = set(rd["panel1_uids"])
        a = p1 & GROUP_A_HINT
        b = p1 & GROUP_B_HINT
        print(
            f"  {rd['round']}: leader UID {rd['leader_uid']} panel_n={rd['panel1_n']} "
            f"top_panels={rd['n_panels_in_top']} | panel1 has UID44/16={bool(a & {16,44})} "
            f"UID55/4={bool(b & {55,4})} | panel1 sample UIDs {rd['panel1_uids'][:8]}"
        )

    print("\n" + "=" * 78)
    print("CROSS-ROUND: do top-1 panels match previous round?")
    print("=" * 78)
    prev_pk = None
    prev_round = None
    for rd in round_dirs:
        by_uid = load_round(rd)
        if not by_uid:
            continue
        ranked = sorted(by_uid.items(), key=lambda x: -x[1]["final"])
        if not ranked:
            continue
        rows = parse_vcf_from_log(ranked[0][1]["log"])
        pk = panel_key(rows)
        if prev_pk is not None:
            j = len(set(pk) & set(prev_pk)) / len(set(pk) | set(prev_pk)) if pk or prev_pk else 0
            same = pk == prev_pk
            print(
                f"  {prev_round} -> {rd.name}: #1 UID {ranked[0][0]} "
                f"n={len(pk)} Jaccard={j:.3f} identical_panel={same}"
            )
        prev_pk, prev_round = pk, rd.name

    print("\n" + "=" * 78)
    print("SUMMARY")
    print("=" * 78)
    print(
        "If only ~2 panel families appear in top ranks each round, validators reward\n"
        "oracle-style fixed panels (per-sample truth lists), not read-calling diversity.\n"
        "Native read miners (niome-native-*) will not match unless pool discovers same sites."
    )
    print("=" * 78)


if __name__ == "__main__":
    main()
