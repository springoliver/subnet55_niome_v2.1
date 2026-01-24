"""Pytest bootstrap: load genomics modules without full niome_subnet package."""
import importlib.util
import sys
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GENOMICS = ROOT / "niome_subnet" / "genomics"
_STUBBED = False


def _ensure_bt_mock():
    if "bittensor" in sys.modules and hasattr(sys.modules["bittensor"], "Synapse"):
        return
    _log = types.SimpleNamespace(
        info=lambda *a, **k: None,
        warning=lambda *a, **k: None,
        error=lambda *a, **k: None,
        debug=lambda *a, **k: None,
    )

    class _Synapse:
        pass

    sys.modules["bittensor"] = types.SimpleNamespace(logging=_log, Synapse=_Synapse)


def _stub_niome_packages():
    global _STUBBED
    if _STUBBED:
        return
    _STUBBED = True
    pkg = types.ModuleType("niome_subnet")
    pkg.__path__ = [str(ROOT / "niome_subnet")]
    gen = types.ModuleType("niome_subnet.genomics")
    gen.__path__ = [str(GENOMICS)]
    sys.modules.setdefault("niome_subnet", pkg)
    sys.modules.setdefault("niome_subnet.genomics", gen)


def load_genomics_module(filename: str):
    """Load niome_subnet/genomics/<filename> into the stub package namespace."""
    _ensure_bt_mock()
    _stub_niome_packages()
    path = GENOMICS / filename
    qual = f"niome_subnet.genomics.{path.stem}"
    if qual in sys.modules:
        return sys.modules[qual]
    spec = importlib.util.spec_from_file_location(qual, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[qual] = mod
    spec.loader.exec_module(mod)
    return mod
