"""Backend selection for the merge kernel.

Deliberately the same shape as :mod:`zna.codec`, which selects between ``zna._pycodec``
and ``zna._accel``: a name-keyed registry, a preference order, and a validated required
API. The codec's two implementations are pinned equal by cross-backend tests, and that
discipline is the reason a data-corrupting optimisation was caught in 0.3.4; the merge
kernel gets the same treatment.

The Python backend is the **reference oracle**, not a fallback. It is never deleted and
never optimised at the cost of clarity — it is what the accelerated backend is defined
to agree with, so it has to stay readable enough to be believed.

Backends implement (``docs/archive/MERGE_ACCURACY_PLAN.md`` §6)::

    POLICY_ABI = 2

    scan(s1, s2rc, len1, len2, match_q, step_q, floor_q, adapter_trimmed)
        -> (shift, score_q, overlap_len, mismatches)

    overlap(s1, s2rc, len1, len2, match_q, step_q, t_table, dfit_table, adapter_trimmed)
        -> (verdict, shift, score_q, overlap_len, mismatches, informative_mismatches)

    process_pair(h1, s1, q1, h2, s2, q2, match_q, step_q, t_table, dfit_table,
                 adapter_trimmed, min_read_length, disagree_q[, npolicy, rng_seed,
                 pair_index, rt_check])
        -> (records, outcome, n_dropped, shift, score_q, overlap_len, mismatches,
            bases_consensus_changed, implausible, npolicy_bases, n_rescued,
            detected_bases, detected_mismatches, readthrough_strong,
            detected_overlap_len)

    merge_chunk(buf1, start1, end1, buf2, start2, end2, match_q, step_q, t_table,
                dfit_table, adapter_trimmed, min_read_length, disagree_q, check_sync,
                base_index[, npolicy, rng_seed, rt_check_pairs])
        -> (blob, consumed1, consumed2, counters, len_hist, olen_hist, insert_hist,
            det_olen_hist, need)

    merge_chunk_records(... as merge_chunk ..., base_index, want_headers[, npolicy,
                        rng_seed, rt_check_pairs])
        -> (seqs, ends, consumed1, consumed2, counters, len_hist, olen_hist,
            insert_hist, det_olen_hist, need)

    split_records(buf, start, max_records) -> (offset, found)

All scoring arguments are integers in the fixed-point scale of :mod:`zna.merge.params`,
and the two tables are ``int64`` buffers (``array('q')``) derived there; no float crosses
this boundary, which is what makes the two implementations comparable for exact equality
rather than approximate agreement. :mod:`zna.merge._pymerge` documents each function.

**The ABI marker.** The argument lists above are not the 0.5.x ones, and a compiled
extension built from an older tree would be called with the wrong arguments -- or,
worse, with the right *number* of them meaning different things. So a backend must
declare ``POLICY_ABI`` equal to :data:`POLICY_ABI`; one that does not (every 0.5.x build
lacks it) is refused exactly as if it had not been built, and everything runs on the
reference backend, loudly, as in the no-extension configuration.
"""
from __future__ import annotations

import importlib
from types import ModuleType
from typing import Optional

_BACKEND_MODULES = {
    "python": "zna.merge._pymerge",
    "accel": "zna.merge._accel",
}

#: Highest priority first.
_PREFERENCE = ("accel", "python")

_REQUIRED_FUNCTIONS = frozenset(
    {"scan", "overlap", "process_pair", "merge_chunk", "merge_chunk_records",
     "split_records"})

#: The backend contract this tree's callers speak. See "The ABI marker" above.
POLICY_ABI = 2

_loaded: dict[str, ModuleType] = {}
_default: Optional[ModuleType] = None
_default_name: Optional[str] = None


def available_merge_backends() -> list[str]:
    """Names of every backend that can be imported AND speaks this tree's ABI."""
    out: list[str] = []
    for name in _BACKEND_MODULES:
        try:
            _load(name)
            out.append(name)
        except ImportError:
            pass
    return out


def get_merge_backend(name: Optional[str] = None) -> ModuleType:
    """Return a backend module by *name*, or the best available one.

    Raises ``ImportError`` if the requested backend cannot be loaded.
    """
    global _default, _default_name
    if name is not None and name != "auto":
        return _load(name)
    if _default is not None:
        return _default
    for candidate in _PREFERENCE:
        try:
            mod = _load(candidate)
        except ImportError:
            continue
        _default, _default_name = mod, candidate
        return mod
    raise ImportError("no merge backend available")


def get_merge_backend_name(name: Optional[str] = None) -> str:
    """Canonical name of the backend ``get_merge_backend(name)`` resolves to."""
    if name is not None and name != "auto":
        _load(name)
        return name
    get_merge_backend()
    assert _default_name is not None
    return _default_name


def _load(name: str) -> ModuleType:
    if name in _loaded:
        return _loaded[name]
    modpath = _BACKEND_MODULES.get(name)
    if modpath is None:
        raise ImportError(
            f"unknown merge backend {name!r}; choose from {sorted(_BACKEND_MODULES)}")
    mod = importlib.import_module(modpath)
    abi = getattr(mod, "POLICY_ABI", None)
    if abi != POLICY_ABI:
        raise ImportError(
            f"merge backend {name!r} implements merge-policy ABI {abi}, but this zna "
            f"needs {POLICY_ABI}: it was built from an older source tree. Rebuild the "
            f"extension (pip install -e . with a C++ toolchain).")
    missing = _REQUIRED_FUNCTIONS - set(dir(mod))
    if missing:
        raise ImportError(
            f"merge backend {name!r} is missing required functions: {sorted(missing)}")
    _loaded[name] = mod
    return mod


# --------------------------------------------------------------------------- #
# the active backend
#
# Held here rather than in each consumer so there is one thing to set and one thing to
# read. Resolved on FIRST USE, not at import: the accelerated backend is an extension
# module and the reference one drags in the consensus table, and neither belongs in the
# import cost of code that only wanted a constant.
# --------------------------------------------------------------------------- #

_active: Optional[ModuleType] = None
_active_name: Optional[str] = None


def use(name: Optional[str] = None) -> str:
    """Select the backend for subsequent calls. Returns its canonical name."""
    global _active, _active_name
    _active = get_merge_backend(name)
    _active_name = get_merge_backend_name(name)
    return _active_name


def active() -> ModuleType:
    """The backend in use, selecting the default if none has been chosen."""
    if _active is None:
        use()
    return _active


def active_name() -> str:
    """Canonical name of the backend in use."""
    if _active is None:
        use()
    assert _active_name is not None
    return _active_name
