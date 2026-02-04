"""Fleet strategy must change pipeline env and pick mode."""
import os

from tests.conftest import load_genomics_module

_ts = load_genomics_module("task_strategy.py")
apply_strategy_profile = _ts.apply_strategy_profile
pipeline_merge_pool = _ts.pipeline_merge_pool
pipeline_pick_mode = _ts.pipeline_pick_mode
pipeline_fallback_strategy = _ts.pipeline_fallback_strategy

_rc = load_genomics_module("read_calling.py")
get_read_calling_rev = _rc.get_read_calling_rev


def _clear_niome_env():
    for key in list(os.environ):
        if key.startswith("NIOME_"):
            os.environ.pop(key, None)


def test_high_recall_pipeline_differs_from_v10():
    _clear_niome_env()
    apply_strategy_profile("high_recall", predicted_band="high")
    assert pipeline_pick_mode() == "recall"
    assert pipeline_merge_pool() is False
    assert os.environ.get("NIOME_MPILEUP_QUAL") == "-q 0 -Q 0"

    _clear_niome_env()
    apply_strategy_profile("v10")
    assert pipeline_pick_mode() == "default"
    assert pipeline_merge_pool() is False
    assert os.environ.get("NIOME_MPILEUP_QUAL") == "-q 2 -Q 2"


def test_v5_style_precision_pick():
    _clear_niome_env()
    apply_strategy_profile("v5_style")
    assert pipeline_pick_mode() == "precision"
    assert os.environ.get("NIOME_NATIVE_RECALL") == "0"


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
    TaskFingerprint = ts.TaskFingerprint

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
