"""zna read-merge: overlap-merge paired-end reads for ZNA / LLM pretraining.

Replaces the fastp PE-merge step. Every pair is scored once, over one axis of candidate
fragment lengths, as a log-likelihood ratio in bits against chance alignment, and the
best-scoring shift is then either **merged** or the pair is **kept** whole
(``docs/archive/MERGE_ACCURACY_PLAN.md`` §2):

* **merge** — the best shift reaches the pair's floor ``T = log2((len1 + len2 - 1) /
  alpha)`` and its mismatches are plausible as sequencing error at the same ``alpha``
  (a repeat outscoring the true overlap is not) -> one full-fragment sequence;
* **keep** — anything else -> both reads unchanged.

The error rate the score and the plausibility test use is ``--error-rate`` (default
0.01), a documented setting rather than an estimate; each run reports the disagreement
its detected overlaps actually show, and warns when at that disagreement the
plausibility test is expected to refuse more than 0.1% of true overlaps.
``--adapter-trimmed`` declares that no read runs past its molecule, which makes
read-through alignments impossible.

Output is a single **mixed interleaved** FASTQ stream (merged reads as singles with the
pair suffix stripped; unmerged pairs as adjacent ``/1``,``/2`` records) consumed by
``zna encode --interleaved``. See ``docs/METHODS.md``.

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
    "scan_unrestricted",
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
    "scan_unrestricted": ".overlap",
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
