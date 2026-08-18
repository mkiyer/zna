"""The ``zna encode --merge-pairs`` input strategy: merge in process, hand over records.

This module is the seam between the merge kernel and the encoder, and it lives
here rather than in ``zna/cli.py`` because it needs ``_RawStream``,
``MergeParams`` and ``DISAGREE_Q`` — none of which belong in encode's startup
(MERGE_PAIRS_PLAN.md §4 step 4).

What crosses the seam is *records with geometry*, not FASTQ text: each emitted
record arrives as ``(seq, is_paired, is_read1, is_read2, has_start, has_end)``
— or with ``labels`` appended — so the encoder never re-derives pairing from
read names.  The name-collision corruption of the two-step path (two
independently merged molecules sharing a read name, silently re-paired into
one fragment) is structurally impossible here.

The geometry transfer is exactly one bit (plan §0.1): ``has_start`` is True
for every record on this path — every emitted read begins at a true fragment
boundary — and ``has_end`` is True only for a merged record, which spans its
whole fragment.  ``IS_RC`` is never supplied: orientation belongs to the
writer's strand-normalization settings, which compose with this stream
unchanged (plan §0.2).
"""
from __future__ import annotations

from typing import Iterator

from .backend import get_merge_backend
from .cli import _READ_BLOCK, _RawStream, _check_drained, _fold, _new_acc
from .fastqio import InputError
from .params import DISAGREE_Q, MergeParams
from ._pymerge import SLOT_MATE2, SLOT_MERGED

#: slot -> (is_paired, is_read1, is_read2, has_start, has_end).  A merged
#: record is a one-record fragment; mates keep their pair identity and reach a
#: true boundary at base 0 only.
_GEOMETRY_BY_SLOT = (
    (False, False, False, True, True),    # SLOT_MERGED
    (True, True, False, True, False),     # SLOT_MATE1
    (True, False, True, True, False),     # SLOT_MATE2
)


def stream_merge_pairs(files, params: MergeParams, *, check_sync: bool = True,
                       merge_backend: str = "auto", label_defs=(), tag_map=None,
                       allow_empty: bool = False, acc=None,
                       quiet: bool = True) -> Iterator[tuple]:
    """Merge two FASTQ files pair by pair, yielding geometry-carrying tuples.

    Yields 6-tuples, or 7-tuples with ``labels`` last when *label_defs* is
    non-empty.  Labels come from each record's OWN source header — R1's for a
    merged record and for MATE1, R2's for MATE2 — exactly as the two-step
    path's extractor reads them (plan §3.5).  A label whose tag is ``ZN`` is
    special: its value is the record's provenance byte taken directly from the
    pair result, with no ``ZN:i:`` tag round-trip (plan §10.4) — the two-step
    path only sees that tag because ``zna merge`` wrote it into a header this
    path never materializes.

    *acc*, when given, is a ``merge/cli`` accumulator the caller can hand to
    ``_assemble_stats`` afterwards for ``--merge-json`` and the run summary.

    Raises :class:`InputError` — *before* the caller's writer can finish, so
    an aborted encode never certifies — on desync, truncation, unequal read
    counts, or (without *allow_empty*) zero emitted records: a zero-pair input
    otherwise writes a 0-record ``.zna`` and a library vanishes from the
    corpus with every stage green.
    """
    backend = get_merge_backend(
        None if merge_backend == "auto" else merge_backend)
    if acc is None:
        acc = _new_acc()
    kw = (params.match_q, params.step_q, params.t_merge_q, params.t_trim_q,
          params.min_read_length, DISAGREE_Q)
    want_headers = bool(label_defs)
    zn_index = None
    if want_headers and tag_map:
        zn = tag_map.get(b"ZN")
        if zn is not None:
            zn_index = zn[0]
    n_labels = len(label_defs)
    target = max(_READ_BLOCK, 2000 * 1024)
    geometry = _GEOMETRY_BY_SLOT

    if want_headers:
        from ..cli import extract_labels_from_header
    r1 = _RawStream(files[0], 1)
    r2 = _RawStream(files[1], 1)
    emitted = 0
    try:
        while True:
            r1.fill(target)
            r2.fill(target)
            if not r1.avail or not r2.avail:
                break
            seqs, ends, c1, c2, counters, lh, oh, ih = (
                backend.merge_chunk_records(
                    r1.buf, r1.pos, len(r1.buf), r2.buf, r2.pos, len(r2.buf),
                    *kw, check_sync, acc[0][0], want_headers,
                    params.npolicy_code, params.rng_seed))
            if not c1 and not c2:
                break              # neither stream holds a complete record
            # One decode per chunk, not per record: sequences are ASCII, so
            # byte offsets are character offsets.
            text = seqs.decode("ascii")
            for seq_off, seq_len, hdr_off, hdr_len, slot, prov in ends:
                seq = text[seq_off:seq_off + seq_len]
                is_paired, is_r1, is_r2, has_start, has_end = geometry[slot]
                if want_headers:
                    hdr = (r2.buf if slot == SLOT_MATE2
                           else r1.buf)[hdr_off:hdr_off + hdr_len]
                    labels = extract_labels_from_header(
                        hdr, tag_map, n_labels, label_defs)
                    if zn_index is not None:
                        labels = (labels[:zn_index] + (prov,)
                                  + labels[zn_index + 1:])
                    yield (seq, is_paired, is_r1, is_r2,
                           has_start, has_end, labels)
                else:
                    yield (seq, is_paired, is_r1, is_r2, has_start, has_end)
            emitted += counters[4]
            _fold(counters, lh, oh, ih, acc)
            r1.pos += c1
            r2.pos += c2
        # Both streams must run out together, and the check runs BEFORE the
        # caller's writer can close: raising out of this generator aborts the
        # encode, so a desynced input leaves an UNCERTIFIED partial file
        # rather than a complete, plausible-looking one (plan §5).
        _check_drained(backend, r1, "R1", "R2")
        _check_drained(backend, r2, "R2", "R1")
        if emitted == 0 and not allow_empty:
            raise InputError(
                "0 records emitted: the inputs hold no pairs. A 0-record .zna "
                "makes a library vanish from a corpus with every stage green; "
                "pass --allow-empty if that is genuinely intended."
            )
    finally:
        r1.close()
        r2.close()
