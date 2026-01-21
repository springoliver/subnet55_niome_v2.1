#!/usr/bin/env python3
"""bcftools --indels-2.0 compatibility shim."""
import importlib.util
import sys
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PIPE = ROOT / "niome_subnet" / "genomics" / "pipeline.py"

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def _load_pipeline():
    if "bittensor" not in sys.modules:
        log = types.SimpleNamespace(info=lambda *a, **k: None, warning=lambda *a, **k: None)
        sys.modules["bittensor"] = types.SimpleNamespace(logging=log)
    spec = importlib.util.spec_from_file_location("niome_pipeline_test", PIPE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_resolve_strips_when_unsupported(monkeypatch=None):
    mod = _load_pipeline()
    mod._INDELS_20_SUPPORTED = False
    assert mod.resolve_mpileup_extra("--indels-2.0") == ""
    assert mod.resolve_mpileup_extra("-q 5 --indels-2.0") == "-q 5"


def test_resolve_keeps_when_supported():
    mod = _load_pipeline()
    mod._INDELS_20_SUPPORTED = True
    assert mod.resolve_mpileup_extra("--indels-2.0") == "--indels-2.0"


if __name__ == "__main__":
    test_resolve_strips_when_unsupported()
    test_resolve_keeps_when_supported()
    print("ok")
