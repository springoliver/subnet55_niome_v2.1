# Live task `768bdcb5-4269-4c72-8705-d77d729b3db8`

| Field | Value |
|-------|--------|
| Region | `chr7:117480000-117670000` (190 kb) |
| Reads | `crt/reads_1.fq` + `crt/reads_2.fq` (presign 2026-05-24 12:30 UTC) |
| Band | **high** (5.24.02 top native = **25 sites**, ~0.91) |
| Prior round | `Results/5.24.02/` (`f38d9a47…`, same FASTQ key) |

**springhot UID = 141** — deploy `ecosystem.spring.config.js` with latest `main`.

## Preview routing

```bash
python tests/analyze_live_task.py
```

## After the round

1. Save `task.json`, `miner.json`, `validator.json` under e.g. `Results/5.24.03/`
2. If manager publishes truth: `real_correct_result/truth.vcf`
3. `python scripts/build_challenge_db.py`

## Important

- **New `task_id`** → new truth loci; do not replay 5.24.02 VCF.
- Same `crt/reads` library → expect **~25** sites with **v5_style**, not 27-site recall panels.
