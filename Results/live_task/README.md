# Live task `f38d9a47-9da4-423b-9215-c32b31ddb7a3`

| Field | Value |
|-------|--------|
| Region | `chr7:117480000-117670000` (190 kb ultra-wide) |
| Reads | `crt/reads_1.fq` + `crt/reads_2.fq` (new presign 2026-05-24) |
| Band | **ultra** (expect ~29-33 truth sites) |
| Round folder | `Results/5.24.02/` |

## Before deploy

```bash
python tests/analyze_live_task.py
```

See `strategy_plan.json` for per-UID strategy routing.

## After the round

1. Save `miner.json` and `validator.json` into `Results/5.24.02/`
2. If manager publishes truth: `Results/5.24.02/real_correct_result/truth.vcf`
3. Rebuild DB: `python scripts/build_challenge_db.py`

## Important

- **New `task_id`** = new truth loci (do not copy VCF from 5.23.04 / 5.24.01).
- Same FASTQ URLs as prior rounds = same read pool, different scored positions each task.
