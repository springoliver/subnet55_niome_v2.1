# NIOME challenge database (full Results history)

All analysis and strategy calibration should use the **cumulative database**, not only the latest round.

## Build / update

After you add a new folder under `Results/` (with `task.json` + `miner.json`):

```bash
python scripts/build_challenge_db.py
```

Incremental (rebuild listed rounds + merge existing DB):

```bash
python scripts/build_challenge_db.py --rounds 5.24.02
```

## Output layout

```
Results/niome_challenge_db/
  manifest.json              # index of all rounds
  rounds/<round_key>.json    # per-round summary (scores, panels, vs truth)
  training/
    strategy_calibration.json    # band → strategy + site-count stats
    truth_position_hotspots.json # frequent truth positions (truth rounds)
    commonly_missed_positions.json
    fleet_timeline.jsonl         # fleet UID performance over time
    top_miner_patterns.json
    fleet_score_by_strategy.json
```

## What gets ingested

| Source | Included |
|--------|----------|
| `Results/5.*` | Yes, if `miner.json` + `task.json` |
| `Results/old_version_result/5.*` | Yes |
| `real_correct_result/truth.vcf` | Linked per round |
| `correct_answer_*/truth.vcf` | Indexed by `task_id` |
| `validator.json` | Not stored (too large); scores come from `miner.json` |

## Reports

```bash
python scripts/report_challenge_db.py
python tests/analyze_fleet_strategy.py      # uses DB when present
python tests/analyze_truth_vs_top.py          # all DB rounds if no args
```

## Adding truth for new rounds

Place manager truth under either:

- `Results/<round>/real_correct_result/truth.vcf`
- `Results/<round>/real_correct_result/cftr2_annotations.json`

Then rebuild the DB. `spring05` / `win` strategy uses the same paths via `NIOME_TRUTH_DIR`.

## Training use

- **`strategy_calibration.json`**: set `NIOME_CURRICULUM_TARGET` per band from historical median truth counts.
- **`truth_position_hotspots.json`**: indel/SNP-heavy regions to watch (not oracle replay).
- **`commonly_missed_positions.json`**: where native fleet under-calls vs truth across history.
