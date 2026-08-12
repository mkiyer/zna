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
"""
from __future__ import annotations

from .pairs import MergeParams, PairOutcome, process_pair
from .overlap import find_overlap, reverse_complement, score_weights, threshold_bits

__all__ = [
    "MergeParams",
    "PairOutcome",
    "process_pair",
    "find_overlap",
    "reverse_complement",
    "score_weights",
    "threshold_bits",
]
