"""zna read-merge: overlap-merge paired-end reads for ZNA / LLM pretraining.

Replaces the fastp PE-merge step. Every pair is scored once, over one axis of candidate
fragment lengths, as a log-likelihood ratio in bits against chance alignment; the
best-scoring shift is then read at two thresholds:

* **merge** (``score >= t_merge``) — an overlapping pair becomes one full-fragment
  sequence,
* **trim** (``t_trim <= score < t_merge``) — the overlap is real but not worth risking a
  chimera, so both reads are kept and R2's redundant 3' bases are cut (the overlapping
  bases are not counted twice in pretraining),
* **keep** (``score < t_trim``) — both reads unchanged.

Output is a single **mixed interleaved** FASTQ stream (merged reads as singles with the
pair suffix stripped; unmerged pairs as adjacent ``/1``,``/2`` records) consumed by
``zna encode --interleaved``. See ``docs/READ_MERGE_REDESIGN.md``.

**The exports below are resolved lazily** (PEP 562). Importing them eagerly would make
``import zna.merge.args`` — which ``zna/cli.py`` does on *every* invocation, just to
register the subcommand — pull in the kernel, the extension module and the 64 KiB
posterior table of :mod:`zna.merge.params`. None of that belongs in the startup of
``zna inspect``, which advertises itself as fast enough to catalogue a corpus. For the
same reason ``zna/__init__.py`` does not re-export this package at all.

Accessing any name here (``from zna.merge import MergeParams``) imports what it needs,
once, and caches it in the module globals.
"""
from __future__ import annotations

from importlib import import_module

__all__ = [
    "MergeParams",
    "PairOutcome",
    "process_pair",
    "find_overlap",
    "reverse_complement",
    "score_weights",
    "threshold_bits",
    "SCALE",
    "to_bits",
    "to_q",
]

_LAZY = {
    "MergeParams": ".params",
    "PairOutcome": ".pairs",
    "process_pair": ".pairs",
    "find_overlap": ".overlap",
    "reverse_complement": ".overlap",
    "score_weights": ".params",
    "threshold_bits": ".params",
    "SCALE": ".params",
    "to_bits": ".params",
    "to_q": ".params",
}


def __getattr__(name):
    try:
        where = _LAZY[name]
    except KeyError:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from None
    value = getattr(import_module(where, __name__), name)
    globals()[name] = value          # subsequent lookups skip __getattr__ entirely
    return value


def __dir__():
    return sorted(set(globals()) | set(__all__))
