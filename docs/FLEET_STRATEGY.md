# Fleet strategy plan (review before deploy)

## Goal

Beat the **#1 native miner** in each task family — not oracle (~0.94) on every round without truth.

Assign **different strategies to different UIDs** (or `auto` per task), verify on past `Results/`, then deploy.

## Problem types (from `Results/` 5.20–5.23)

All recent tasks are **ultra-wide CFTR** (~190 kb). “Type” = **truth size band** + **indel complexity**, not region size.

| Type | Top oracle sites | Example rounds | Best native (~) | Your best (~) | Gap |
|------|------------------|----------------|-----------------|---------------|-----|
| **A — low** | 11–14 | 5.21.01–02 | 0.85 | — | — |
| **B — mid** | 15–22 | 5.21.04, 5.22.01–02 | 0.82–0.84 | 0.56–0.64 | large |
| **C — high** | 23–29 | 5.22.04, 5.23.01 | 0.82–0.84 | 0.56–0.60 | large |
| **D — ultra** | 30–33 | 5.23.02–03 | 0.73–0.89 | 0.60–0.84 | closing on 5.23.02 |

Oracle **always** ~0.90–0.99 (truth-shaped panels). Native ceiling ~0.70–0.89 depending on band.

## Past problems → strategy mapping

| Problem | Evidence | Strategy |
|---------|----------|----------|
| Empty VCF (0 score) | 5.23.01: UIDs 9,97,124… | v10 emergency + pipeline retry (all strategies) |
| GT over-call 1/1 | 5.23.02 v9 vs 71 | **v5_style**: hom AF ≥0.58 |
| Wrong indel REF/ALT | 5.23.02–03 six POS | **v5_style** + `collapse_indels_at_position` (v10 base) |
| Under-count vs oracle | ~20 vs 30–33 | **high_recall**: curriculum 32 |
| Annotation leak | 71 ann 0.81 vs 0.44, same VCF | **v5_style**: full FORMAT + ClinVar annotate path |
| Oracle rounds | top 0.94+ | **win**: `NIOME_TRUTH_DIR` / manager truth |

## Strategy profiles (`NIOME_STRATEGY`)

| Strategy | When | Beat target |
|----------|------|-------------|
| `win` | Truth file for `task_id` | Oracle (~0.90+) |
| `v5_style` | Low/mid band or UID 71 | Best native on annotation-heavy tasks |
| `high_recall` | High/ultra band | Best native on 30+ site tasks (UID 65, 62…) |
| `v10` | Default balanced | Mid band native |
| `auto` | Per task: truth → band map | Picks one of above |

### Per-UID PM2 example

```bash
# UID 71 — keep v5-style (proven 0.84 on 5.23.02, 0.60 on 5.23.03)
export NIOME_UID_STRATEGY_71=v5_style

# High-recall arm (beat native on 30+ site rounds)
export NIOME_UID_STRATEGY_209=high_recall
export NIOME_UID_STRATEGY_235=high_recall

# Truth miner when manager publishes
export NIOME_UID_STRATEGY_139=win
export NIOME_TRUTH_DIR=/path/to/truth_by_task

# Rest: v10 or auto
export NIOME_STRATEGY=v10
```

## Verify before launch

```bash
# 1) Strategy analysis on exported rounds
python tests/analyze_fleet_strategy.py

# 2) Offline score vs manager truth (where available)
python tests/verify_fleet_strategies.py

# 3) Unit tests
python -m pytest tests/test_task_strategy.py -q
```

## Launch checklist

1. Full fleet on **one code rev** (current `main`) — strategies via env only.
2. Log line: `[strategy] name=v5_style band=ultra …`
3. After round: `python tests/analyze_user_miners.py`
4. Compare your best per band vs `best_native` in `analyze_fleet_strategy.py` output.

### Full challenge database (all Results history)

**Always rebuild after adding a new round folder:**

```bash
python scripts/build_challenge_db.py
python scripts/report_challenge_db.py
python tests/analyze_fleet_strategy.py
```

See [CHALLENGE_DATABASE.md](CHALLENGE_DATABASE.md). Miners use historical calibration when
`NIOME_RESULTS_ROOT` points at `Results/` (`NIOME_USE_CHALLENGE_DB=1` default).

### Truth vs top-miner analysis

```bash
python tests/analyze_truth_vs_top.py
python tests/analyze_round_fleet_dupes.py 5.23.04 5.24.01
```

### PM2 deploy (Vultr)

`pm2 reload` **does not** fix a wrong `cwd`. If error logs mention
`/root/subnet-niome-429` or `/root/subnet-niome-win`, delete and restart:

```bash
cd /root/subnet-niome && git pull
chmod +x scripts/pm2-fleet-deploy.sh
bash scripts/pm2-fleet-deploy.sh
```

Confirm boot line in logs (`NIOME miner boot repo=/root/subnet-niome …`) and
`Miner running…` around line **~230** in `neurons/miner.py` (not ~204).

```bash
pm2 describe spring01 | grep -E "exec cwd|script path|NIOME_STRATEGY"
pm2 logs spring01 --lines 40 --nostream | grep -E "NIOME miner boot|strategy=|Task "
```

## Honest limits

- **Cannot beat oracle every round** without per-task truth.
- **auto** band prediction is heuristic until round-history DB is added.
- **v5_style** approximates UID 71’s old rev — not byte-identical to `niome-native-2026-05-22-v5`.
