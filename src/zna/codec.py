"""
Backend selection for ZNA codec operations.

Each backend is a Python module that exposes four functions:

    encode_sequence(seq: str) -> bytes
    reverse_complement(seq: str) -> str
    encode_block(seqs, flags, len_bytes_fmt, npolicy, do_rc_r1, do_rc_r2)
        -> (flags_bytes, lengths_bytes, seqs_bytes)
    decode_block(flags_data, lengths_data, seqs_data, len_bytes, count)
        -> list[(seq, is_paired, is_read1, is_read2)]

Available backends
------------------
- ``"python"``  — pure-Python reference implementation (always available).
- ``"accel"``   — C++ nanobind extension (available when built).

Usage::

    from zna.codec import get_backend

    codec = get_backend()          # best available
    codec = get_backend("python")  # force pure-Python
"""
from __future__ import annotations

import importlib
from types import ModuleType
from typing import Optional


# Registry: name → module path
_BACKEND_MODULES = {
    "python": "zna._pycodec",
    "accel": "zna._accel",
}

# Preference order (highest priority first)
_PREFERENCE = ("accel", "python")

# Cached results
_loaded: dict[str, ModuleType] = {}
_default: Optional[ModuleType] = None
_default_name: Optional[str] = None


# ---------------------------------------------------------------------------
# Required API surface — used for validation only
# ---------------------------------------------------------------------------
_REQUIRED_FUNCTIONS = frozenset({
    "encode_sequence",
    "reverse_complement",
    "encode_block",
    "decode_block",
})


def _validate_backend(mod: ModuleType, name: str) -> None:
    """Raise if *mod* does not expose the required codec functions."""
    missing = _REQUIRED_FUNCTIONS - set(dir(mod))
    if missing:
        raise ImportError(
            f"Backend '{name}' is missing required functions: {sorted(missing)}"
        )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def available_backends() -> list[str]:
    """Return the names of all backends that can be imported."""
    result: list[str] = []
    for name, modpath in _BACKEND_MODULES.items():
        try:
            importlib.import_module(modpath)
            result.append(name)
        except ImportError:
            pass
    return result


def get_backend(name: Optional[str] = None) -> ModuleType:
    """Return a backend module by *name*, or the best available one.

    Parameters
    ----------
    name
        ``"python"``, ``"accel"``, or ``None`` (auto-select).

    Returns
    -------
    ModuleType
        A module object whose top-level functions implement the codec API.

    Raises
    ------
    ImportError
        If the requested backend cannot be loaded.
    """
    global _default, _default_name

    if name is not None:
        return _load(name)

    if _default is not None:
        return _default

    # Auto-select: try each in preference order
    for candidate in _PREFERENCE:
        try:
            mod = _load(candidate)
            _default = mod
            _default_name = candidate
            return mod
        except ImportError:
            continue

    raise ImportError("No ZNA codec backend available")


def get_backend_name(name: Optional[str] = None) -> str:
    """Return the canonical name of the backend that ``get_backend(name)`` resolves to."""
    if name is not None:
        _load(name)  # ensure it's valid
        return name

    get_backend()  # ensure _default_name is set
    assert _default_name is not None
    return _default_name


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _load(name: str) -> ModuleType:
    """Import, validate, cache, and return the backend module *name*."""
    if name in _loaded:
        return _loaded[name]

    modpath = _BACKEND_MODULES.get(name)
    if modpath is None:
        raise ImportError(f"Unknown backend: '{name}'. Choose from {list(_BACKEND_MODULES)}")

    mod = importlib.import_module(modpath)

    # The C++ module exports decode_block_columnar — alias it transparently.
    if not hasattr(mod, "decode_block") and hasattr(mod, "decode_block_columnar"):
        mod.decode_block = mod.decode_block_columnar  # type: ignore[attr-defined]

    _validate_backend(mod, name)
    _loaded[name] = mod
    return mod
