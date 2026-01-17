# NIOME SN55 — New-Version Miner Method (v3)

**Status:** Implemented in `subnet-niome` (`read_calling.py`, `pipeline.py`, `task_profile.py`, `neurons/miner.py`).

**Live ultra-wide task** `3ec52dc3` (`chr7:117480000-117670000`, 190 kb, `expected_variant_count=0`):
- Auto profile: `ultra_wide` — stricter mpileup (`-q10`), prefer norm VCF, filter long indels / 117504 noise.
- **v3-ultra4** (`read_calling` rev `2026-05-21-v3-ultra4`): read-first SNPs + indels (max allele 24 bp), trim only if **>28** calls; supports 5.21.03 (25 truth sites). No fixed 14-cap.
- Test: `python tests/run_live_task.py` (fixture in `Results/live_task/task.json`).

This document compares your **sn55-05** work and **old_version_result** rounds with the subnet manager’s **new CFTR rules**, analyzes **Results/sample**, and defines the pipeline to deploy on the current `subnet-niome` repo.

---

## 1. What changed (old → new)

| Topic | Old version (your sn55-05) | New version (manager) |
|--------|---------------------------|------------------------|
| **Variant count** | `expected_variant_count` = 6–21; miner **must submit exactly N** or validator drops response | `expected_variant_count` = **0** (hidden); submit **all variants supported by reads** |
| **Truth composition** | Often inferred by matching **ClinVar top-N** to N | Truth includes **non-ClinVar** variants; do **not** build truth from ClinVar alone |
| **VCF strategy** | `clinvar_priority` / `read_priority` + **canonical panels** (mid7, dense13) tuned to N | **Read-evidence first**; ClinVar only for **GT overlay** and **annotations** |
| **Annotations** | Per submitted ClinVar row (`_annotations_for`) | **ClinVar variants only**; over-reporting penalized (`max(len(truth), len(miner))` denominator) |
| **Final score** | `0.7 × VCF + 0.3 × annotation` | Same |

**Validator note:** Production validators should skip the count gate when `expected_variant_count == 0`. The public repo still has the old check in `niome_subnet/validator/forward.py` (~line 218); miners on mainnet are evaluated against the **deployed** validator, not necessarily this file.

---

## 2. Review of your sn55-05 work

### What worked

- **Full pipeline:** GRCh38 chr7 ref → BWA → `bcftools mpileup/call/norm` (same stack validators use in `scoring.py`).
- **Read-based GT:** AD/DP → `0/1` vs `1/1` (fixed ~0.81 → ~0.99 on dense rounds when sites were right).
- **SNP preference:** Avoided wrong indel normalization (e.g. `117548795 CGG>C` vs truth `117548801 C>T`).
- **Caching + per-task locks** in `neurons/miner.py` (good for repeated validator queries).
- **CFTR2 annotations** via `cftr_lookup.py` + ClinVar IDs.

### What failed (from `Results/LEARNINGS.md` and round JSON)

| Failure mode | Example | Root cause |
|--------------|---------|------------|
| **Wrong cluster** | 5.19.02/03: 117504 noise vs mid/dense truth | Strategy picked region by **N + region_end**, not reads |
| **Forced panel size** | ab03f860: submit 7 ClinVar sites, truth wanted 7 **different** sites | `_CANONICAL_MID7` + `expected=7` |
| **GT format** | `GT 1/1` with `AD 1,18` parsed as broken GT | VCF FORMAT not normalized |
| **Zero score** | Wrong REF/ALT or missing FORMAT | Indels / wrong chrom naming |
| **Annotation = 0** | Many rounds | No ClinVar ID on rows or empty `cftr_annotations` |

**Conclusion:** sn55-05 was optimized to **guess N and ClinVar panels**. The new subnet rewards **complete read-backed VCF** + **sparse, accurate ClinVar annotations**.

---

## 3. Sample data analysis (`Results/sample`)

| Asset | Content |
|-------|---------|
| `truth.vcf` | **10** variant lines (read truth), positions `117514266` … `117652897` |
| `cftr2_annotations.json` | **8** ClinVar IDs (7181, 53312, …) — **not** 10 |

**Overlap:**

- In truth but **not** in annotations: `117514266`, `117652038` → **non-ClinVar (or no CFTR2 row)** — must appear in **VCF**, must **not** appear in `cftr_annotations`.
- In annotations: 8 positions matching 8 truth SNPs (e.g. 7181 → `117530899G>A`).

**Old sn55-05 behavior on this sample:** Would likely pick **7–9 ClinVar pathogenic sites** (compact/sparse heuristic) → **low recall** (missing `117514266`, `117652038`, etc.).

**New method behavior:** Submit **~10 read-called SNPs** in region; annotate **only 8** ClinVar-mapped variants.

---

## 4. Recommended method: “Read-all VCF + ClinVar-only annotations”

### Phase A — Inputs (unchanged)

1. Download `read1_fastq`, `read2_fastq` from task URLs.
2. Reference: **GRCh38 chr7** (`NIOME_HG38_REF` or `~/.niome/ref/chr7.fa`) — same coordinate system as validators.
3. Region: `task.genome_context.region` with **+5 kb pad** for edge variants.

### Phase B — Alignment & calling (core)

```text
bwa mem -t 4 ref R1 R2 | samtools sort -o aln.bam
samtools index aln.bam
bcftools mpileup -f ref -r <padded_region> -a AD,DP -q 1 -Q 1 --max-depth 8000 aln.bam \
  | bcftools call -mv -Ov -o raw.vcf
bcftools norm -f ref -m -both -c w raw.vcf -Ov -o norm.vcf
```

**Do not** trim to `expected_variant_count`. When `expected_variant_count == 0`, treat as **no limit**.

### Phase C — Variant selection (new logic)

Keep **all** normalized variants in the task window that pass:

| Filter | Suggested rule |
|--------|----------------|
| Region | `region_start ≤ POS ≤ region_end` |
| Quality | `FILTER=PASS` or QUAL ≥ 20 |
| Depth | FORMAT/DP ≥ 4 (validators weight depth &lt;10 lower, but still score) |
| Allele | Split multiallelics; one ALT per line |
| Genotype | Must have callable GT from reads (`0/1`, `1/1`, `0/0` excluded if hom-ref) |
| Zygosity | Derive GT from AD: alt_frac &lt; 0.15 → skip (hom-ref); 0.15–0.85 → `0/1`; ≥0.85 → `1/1` |

**Do not:**

- Rank by ClinVar pathogenicity to **cap count**.
- Inject `_CANONICAL_MID7` / `_CANONICAL_DENSE` to fill slots.
- Prefer 117504 noise block unless read support is strong (keep as **low-priority** filter, not panel driver).

**Optional ClinVar overlay:** If a read call matches a ClinVar REF/ALT at same POS, attach ClinVar ID in VCF `ID` field (helps debugging; use `ID=.` if env `NIOME_VCF_DOT_ID=1` for top-miner style).

### Phase D — Annotations (strict)

1. Run `build_cftr_annotations(vcf_path)` on **final VCF only**.
2. Include **only** variants with a ClinVar variation ID (rsID / numeric key as in sample JSON).
3. **Never** dump full CFTR2 database or all ClinVar candidates in region.
4. Schema per variant (matches `score_annotations`):

```json
{
  "<clinvar_id>": {
    "hgvs": "NC_000007.14:g.<pos><ref>><alt>",
    "clinical_significance": "Pathogenic",
    "drug_response": {
      "ivacaftor": "responsive|non_responsive",
      "tezacaftor_ivacaftor": "...",
      "elexacaftor_tezacaftor_ivacaftor": "...",
      "lumacaftor_ivacaftor": "..."
    }
  }
}
```

Scoring uses **truth keys only**; extra miner keys hurt via `denominator = max(len(truth), len(miner))`.

### Phase E — Response

- `synapse.vcf_content` = full VCF string (header + all selected variants).
- `synapse.cftr_annotations` = ClinVar-only dict (may be `{}` if no ClinVar sites called).
- Finish within **60s** (`FORWARD_TIMEOUT`).

---

## 5. Code migration map (sn55-05 → current repo)

| sn55-05 file | Action |
|--------------|--------|
| `niome_subnet/genomics/pipeline.py` | **Port** — remove `expected_n` from strategy hint; always read-call |
| `niome_subnet/genomics/clinvar_strategy.py` | **Refactor** → `read_calling.py`: drop `[:expected]`, panels, `choose_strategy` N-dependence |
| `niome_subnet/genomics/cftr_lookup.py` | **Port as-is** |
| `niome_subnet/genomics/task_strategy.py` | **Deprecate** or use only `region_len` for mpileup strictness, not N |
| `neurons/miner.py` | **Port** sn55-05 miner; remove cache key on `expected_variant_count` |
| `setup_miner.py` | **Port** |
| `tests/score_correct_answers.py` | **Port** — set `expected_variant_count=0` or count from truth for local tests |

### Functions to delete or gate behind `NIOME_LEGACY_N=1`

- `_CANONICAL_MID7`, `_CANONICAL_DENSE` forced panels
- `_select_variants(..., expected)` top-N ClinVar fill
- `_select_read_variants` loop `if len(selected) >= expected: break`
- `choose_strategy(..., expected_n)` when `expected_n == 0`

### Keep

- `_gt_from_read_call` / AD-based zygosity
- `_is_simple_snp` preference (not exclusion of indels in truth)
- `bcftools norm` matching validator
- Per-task asyncio lock + cache (key = `task_id:region` only)

---

## 6. Local validation workflow

On **Linux/WSL** with `bwa`, `samtools`, `bcftools`:

```bash
export PYTHONPATH="$(pwd)"
python setup_miner.py

# Score against sample (add a small test script or extend score_correct_answers):
# - truth: Results/sample/truth.vcf
# - reads: Results/sample/read_1.fq, read_2.fq
# - annotations: Results/sample/cftr2_annotations.json
# - task.expected_variant_count = 0
```

**Targets on sample:**

- VCF recall: match **10/10** positions with correct REF/ALT/GT
- Annotations: **8 keys** matching `cftr2_annotations.json` (not 10)

---

## 7. Deployment checklist

1. Linux server, Python 3.12, tools installed (`setup_miner.py`).
2. Register netuid **55**, open **axon port**.
3. `pm2 start neurons/miner.py --name niome_miner -- ...`
4. Log line: `strategy=read_all rev=2026-05-20-v3` and `submitted=N variants` (N **not** forced to 0 or 7).
5. Monitor W&B / validator scores: VCF should dominate (70%); annotation mistakes less costly if VCF is complete.

---

## 8. Summary

| Layer | Old sn55-05 | New method |
|-------|-------------|------------|
| **Goal** | Guess N + ClinVar panel | Recover **all read-supported** variants |
| **VCF** | Top-N ClinVar/read hybrid | **Full** filtered call set in region |
| **Annotations** | All submitted ClinVar rows | **Only** ClinVar sites in VCF; no over-report |
| **Panels** | mid7, dense13, compact | **Removed** |
| **Sample** | Would under-call (≈7–8) | Target **10 VCF + 8 annotations** |

Your sn55-05 investment (alignment, norm, GT, ClinVar DB, miner shell) remains the foundation; **remove count-gaming** and **read-all** is the winning shift for the new subnet.
