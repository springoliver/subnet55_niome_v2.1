#!/usr/bin/env python3
import json
import re
from pathlib import Path

p = Path(__file__).parents[1] / "Results" / "5.22.04" / "miner.json"
items = json.loads(p.read_text(encoding="utf-8", errors="replace"))["items"]
uids = {1, 4, 8, 16, 44, 55, 167, 114, 248, 252}
best = {}
for it in items:
    if it.get("validator_uid") not in (119, 154):
        continue
    u = it["miner_uid"]
    if u not in uids:
        continue
    f = float(it.get("final_score", 0) or 0)
    if u in best and f <= best[u][0]:
        continue
    log = it.get("log", "")
    st = set()
    m = re.search(r"Miner VCF\n(.*)", log, re.S)
    if m:
        for line in m.group(1).splitlines():
            if line.startswith("chr7"):
                q = line.split("\t")
                st.add((int(q[1]), q[3], q[4]))
    best[u] = (f, it.get("miner_hotkey", "")[:24], st)

print("5.22.04 key UIDs:")
for u in sorted(best):
    f, hk, st = best[u]
    print(f"  UID {u:3} final={f:.4f} n={len(st)} hotkey={hk}")

import itertools

for a, b in itertools.combinations(sorted(best), 2):
    sa, sb = best[a][2], best[b][2]
    j = len(sa & sb) / len(sa | sb) if sa | sb else 0
    print(f"  jaccard {a}-{b}: {j:.3f} identical={sa == sb}")
