"""Fleet strategy must change pipeline env and pick mode — not forced site counts."""
import os

from tests.conftest import load_genomics_module

_ts = load_genomics_module("task_strategy.py")
apply_strategy_profile = _ts.apply_strategy_profile
pipeline_merge_pool = _ts.pipeline_merge_pool
pipeline_pick_mode = _ts.pipeline_pick_mode
pipeline_fallback_strategy = _ts.pipeline_fallback_strategy
STRATEGY_AXES = _ts.STRATEGY_AXES

_rc = load_genomics_module("read_calling.py")
get_read_calling_rev = _rc.get_read_calling_rev

_cf = load_genomics_module("cftr_lookup.py")
_normalize_clnsig = _cf._normalize_clnsig
_pick_genomic_hgvs = _cf._pick_genomic_hgvs


def _clear_niome_env():
    for key in list(os.environ):
        if key.startswith("NIOME_"):
            os.environ.pop(key, None)


def test_high_recall_pipeline_differs_from_v10():
    _clear_niome_env()
    apply_strategy_profile("high_recall", predicted_band="high")
    assert pipeline_pick_mode() == "recall"
    assert pipeline_merge_pool() is True
    assert os.environ.get("NIOME_CURRICULUM_TARGET") == "0"
    assert os.environ.get("NIOME_MPILEUP_QUAL") == "-q 0 -Q 0"

    _clear_niome_env()
    apply_strategy_profile("v10", predicted_band="high")
    assert pipeline_pick_mode() == "default"
    assert pipeline_merge_pool() is False
    assert os.environ.get("NIOME_CURRICULUM_TARGET") == "0"


def test_v5_style_precision_pick():
    _clear_niome_env()
    apply_strategy_profile("v5_style", predicted_band="high")
    assert pipeline_pick_mode() == "precision"
    assert os.environ.get("NIOME_NATIVE_RECALL") == "0"
    assert os.environ.get("NIOME_CURRICULUM_TARGET") == "0"


def test_strategy_axes_differ_by_method_not_count():
    assert STRATEGY_AXES["v5_style"]["pipeline_pick"] == "precision"
    assert STRATEGY_AXES["high_recall"]["pipeline_pick"] == "recall"
    assert STRATEGY_AXES["high_recall"]["merge_pool"] == "1"
    assert STRATEGY_AXES["v5_style"]["merge_pool"] == "0"


def test_read_calling_rev_runtime():
    _clear_niome_env()
    apply_strategy_profile("high_recall")
    assert get_read_calling_rev().endswith("+fleet-recall")

    apply_strategy_profile("v10")
    assert get_read_calling_rev().endswith("+fleet-v10")


def test_win_fallback_without_truth():
    _clear_niome_env()
    name = pipeline_fallback_strategy("win", "ultra", False)
    assert name == "v5_style"


def test_crt_reads_band_is_high_not_ultra():
    _clear_niome_env()
    from tests.conftest import load_genomics_module

    ts = load_genomics_module("task_strategy.py")

    class _Ctx:
        region = "chr7:117480000-117670000"

    class _In:
        read1_fastq = "https://bucket/crt/reads_1.fq"
        read2_fastq = "https://bucket/crt/reads_2.fq"

    class _Task:
        genome_context = _Ctx()
        input = _In()

    fp = ts.fingerprint_task(_Task())
    assert fp.predicted_band == "high"
    strat = ts.resolve_strategy(_Task(), miner_uid=141, truth_available=False)
    assert strat == "v5_style"


def test_clnsig_normalize_matches_cftr2():
    assert (
        _normalize_clnsig("Conflicting_classifications_of_pathogenicity")
        == "Conflicting classifications of pathogenicity"
    )
    assert _normalize_clnsig("Likely_pathogenic") == "Likely pathogenic"


def test_pick_genomic_hgvs_prefers_nc():
    hgvs = _pick_genomic_hgvs(
        "NM_123|c.1A>G|NC_000007.14:g.117499010A>T",
        "chr7",
        "117499010",
        "A",
        "T",
    )
    assert hgvs == "NC_000007.14:g.117499010A>T"
