# NIOME Native method

**Revision:** `niome-native-2026-05-22-v5`  
**Module:** `niome_subnet/genomics/evidence_selection.py`

This is the proprietary miner method — not copied from top-miner VCFs or fixed oracle position lists. It is calibrated **offline** from manager `real_correct_result` truth panels (5.21.01–03) and runs **only on read evidence** at task time.

## Why not copy top miners

- Top miners often share the same variant set → similar scores → weight split, not exclusive #1.
- Validator on-chain truth may differ from manager files (5.21.03 showed high scores with low overlap to manager VCF).
- Harder tasks need **indel recall** and **adaptive count**, not a fixed 14- or 17-variant panel.

## Manager truth learnings (offline)

| Round   | Truth count | Composition      |
|---------|-------------|------------------|
| 5.21.01 | 13          | SNPs             |
| 5.21.02 | 12          | SNPs             |
| 5.21.03 | 25          | 16 SNP + 9 indel |

Implications encoded in Native (not as runtime position lists):

- **No target count** — not 10, not 16–20, not 14. Count = variants that pass read gates from mpileup.
- v3 removed **score ≥ 42** on strict tier (that silently clustered many tasks at ~16–18 sites).
- Trim only when calls exceed **32** (safety cap; 5.21.03 truth was 25).
- Ultra-wide CFTR (~190 kb): `-q 5 -Q 5` mpileup, `max_indel_len=32`, 12 bp dedupe for complex loci.
- Noise band `117504200–117504400`: **score penalty only** (v5; 5.22.03 truth has real variants at 117504296 / 117504400).
- Position-aware relaxed thresholds at truth edge/tail (via `task_profile.thresholds_for_position`).

## Pipeline (unchanged flow, native selection)

1. Align reads (GRCh38 chr7 or CFTR slice).
2. `bcftools mpileup/call/norm` on padded task region (+ retry at `-q0` for recall).
3. Pick best candidate VCF by **native in-window call count** (`count_calls_in_region`).
4. **NIOME Native** `native_select_variants()` → final VCF.
5. Annotations: ClinVar IDs only for variants present in submitted VCF (`cftr_lookup`).

## Native selection algorithm

1. Pool all calls in task window (dedupe by pos/ref/alt).
2. **Strict tier:** read gates only (QUAL/AD/DP/AF/indels) → dedupe 12 bp.
3. **Recall:** medium tier only if ≤2 strict sites (calling failed); relaxed only if zero.
4. If count &gt; **32**: trim lowest evidence scores.
5. Assign GT from AD/DP (`gt_from_read_call`).

Evidence score uses: QUAL, AD, DP, PASS, SNP/indel length, core band bonus, noise penalty, ClinVar overlay (for ID/score only — not a submit whitelist).

## Scoring alignment (v2.1)

- `count_penalty = min(miner_n, truth_n) / max(...)` — staying near truth count matters.
- VCF: `0.4×F1 + 0.3×P + 0.3×R`, then × count_penalty.
- Final: `0.7×VCF + 0.3×annotation`.

Native targets high recall with controlled FPs so F1 and count_penalty both stay high.

## Deploy checklist

1. Sync `niome_subnet/genomics/` + `neurons/miner.py` to miner hosts.
2. Restart miner (PM2); confirm log: `rev=niome-native-2026-05-22-v5`.
3. Linux: `bwa`, `bcftools`, `samtools` on PATH; `NIOME_USE_HG38=1` (default).
4. Offline check: `python tests/benchmark_native_truth.py`.
5. Full pipeline (needs BAM): run on `Results/*/real_correct_result/read_*.fq` when tools available.

## Environment

| Variable            | Default | Purpose                          |
|---------------------|---------|----------------------------------|
| `NIOME_USE_HG38`    | 1       | GRCh38 chr7 reference            |
| `NIOME_VCF_DOT_ID`  | 1       | Use `.` for non-ClinVar IDs       |
| `NIOME_BWA_THREADS` | 8       | Alignment threads                |

Do **not** set `NIOME_METHOD=ultra4` — production uses Native only via `select_read_variants`.
