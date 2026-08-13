"""Score `zna merge` and `fastp --merge` against simulated ground truth.

    python compare.py --sim-prefix sim --out results/ [--threads 4]

Runs both tools on the *same* input, joins their output to the sidecar written by
`simulate.py` by read ID, and writes `report.md`, `summary.json` and — the point of the
exercise — `zna_errors.tsv`, one row for every pair `zna merge` got wrong.

**What "wrong" means here, precisely.** The mates are always full-length, so a merged
record's length determines the shift the tool inferred (`L = shift + len2`). Therefore
`len(emitted) == frag_len` is *equivalent* to "the tool inferred the true overlap", and
the two are not independent checks. That is what makes contract C2 — a merged record is
its fragment exactly — testable with a length comparison, and what makes a chimera
impossible to hide: a pair with no true overlap cannot produce a record longer than
`2 * readlen - 1`, so it can never accidentally land on the true fragment length.

Base 0 of a merged record is structurally R1's base 0, so C1 cannot fail without a code
defect; it is checked anyway, by asking whether the emitted record actually aligns to
the fragment at offset 0 (`frame_violations`).

**Reconstruction accuracy is scored in three columns, not one.** A merged record can
disagree with the fragment because the merger mis-resolved a disagreement, or simply
because both mates were sequenced wrong there. So each merged record is scored against

* the **R1-wins baseline** — what a merger that always keeps R1 in the overlap would
  have emitted, and
* the **oracle floor** — what the best possible consensus could achieve, which is not
  zero: a position where both mates are wrong is unrecoverable.

Both come from the sidecar's per-base error record, so they are exact rather than
estimated, and the useful number is where the tool sits between them.

Nothing here is part of the zna package or its test suite.
"""
from __future__ import annotations

import argparse
import gzip
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from simulate import rc                              # noqa: E402

try:
    from zna.merge.params import SCALE               # the scan's fixed-point scale
except ImportError:                                  # pragma: no cover
    raise SystemExit("this script needs zna importable (it re-runs zna's own scan)")

MERGED, R1, R2 = 0, 1, 2

#: Row cap per category in the errors TSV. The file is meant to be *read*; the summary
#: carries the true counts, and the header records what was dropped.
DEFAULT_ROW_CAP = 20_000


# --------------------------------------------------------------------------- #
# ground truth
# --------------------------------------------------------------------------- #

class Truth:
    """The sidecar, in columns. `frag` is the fragment as sequenced (R1's frame)."""

    COLUMNS = ("read_id", "chrom", "start", "strand", "frag_len", "true_ovl",
               "read_through", "len1", "len2", "n_err1", "n_err2", "err_positions",
               "fragment")

    def __init__(self, path):
        chrom, start, strand, err, frag = [], [], [], [], []
        flen, ovl, rt, l1, l2, e1, e2 = [], [], [], [], [], [], []
        index = {}
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            if tuple(header) != self.COLUMNS:
                raise SystemExit(
                    f"{path}: unexpected columns\n  got  {header}\n  want {list(self.COLUMNS)}")
            for i, line in enumerate(fh):
                f = line.rstrip("\n").split("\t")
                index[f[0]] = i
                chrom.append(f[1]); start.append(int(f[2])); strand.append(f[3])
                flen.append(int(f[4])); ovl.append(int(f[5])); rt.append(int(f[6]))
                l1.append(int(f[7])); l2.append(int(f[8]))
                e1.append(int(f[9])); e2.append(int(f[10]))
                err.append(f[11]); frag.append(f[12].encode())
        self.n = len(frag)
        if len(index) != self.n:
            raise SystemExit(f"{path}: duplicate read IDs ({self.n - len(index)}); the "
                             f"ID is the join key and has to be unique")
        self.index = index                       # read id -> row
        self.index_id = list(index)              # row -> read id (insertion-ordered)
        self.chrom, self.start, self.strand = chrom, start, strand
        self.frag_len = np.array(flen, dtype=np.int32)
        self.true_ovl = np.array(ovl, dtype=np.int32)
        self.read_through = np.array(rt, dtype=np.int8)
        self.len1 = np.array(l1, dtype=np.int32)
        self.len2 = np.array(l2, dtype=np.int32)
        self.n_err1 = np.array(e1, dtype=np.int32)
        self.n_err2 = np.array(e2, dtype=np.int32)
        self.err = err
        self.frag = frag

    def errors(self, i):
        """``({pos: base}, {pos: base})`` for the two mates, in read coordinates."""
        tok = self.err[i]
        if tok == ".":
            return {}, {}
        a, b = {}, {}
        for t in tok.split(","):
            d = a if t[0] == "1" else b
            d[int(t[2:-1])] = ord(t[-1])
        return a, b

# --------------------------------------------------------------------------- #
# reading tool output
# --------------------------------------------------------------------------- #

def _open(path):
    if str(path).endswith(".gz"):
        pigz = shutil.which("pigz")
        if pigz:
            p = subprocess.Popen([pigz, "-dc", str(path)], stdout=subprocess.PIPE,
                                 bufsize=1 << 20)
            return p.stdout, p
        return gzip.open(path, "rb"), None
    return open(path, "rb", buffering=1 << 20), None


def read_fastq(path):
    """Yield ``(header, seq)``; qualities are read and discarded."""
    fh, proc = _open(path)
    try:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().rstrip(b"\r\n")
            fh.readline()
            fh.readline()
            yield h[1:].rstrip(b"\r\n"), s
    finally:
        fh.close()
        if proc is not None:
            proc.wait()


def zna_records(path):
    """`zna merge`'s one mixed interleaved stream, as ``(kind, base_id, seq)``.

    A merged record has had its `/1`,`/2` suffix stripped and carries a trailing
    `merged_<n1>_<n2>` token; an unmerged one keeps its suffix. Classifying on the
    suffix rather than the token is what makes this agree with how ZNA re-pairs records.
    """
    for h, s in read_fastq(path):
        rid = h.split(None, 1)[0]
        if rid.endswith(b"/1"):
            yield R1, rid[:-2].decode(), s
        elif rid.endswith(b"/2"):
            yield R2, rid[:-2].decode(), s
        else:
            yield MERGED, rid.decode(), s


def fastp_records(merged, un1, un2):
    for kind, path in ((MERGED, merged), (R1, un1), (R2, un2)):
        for h, s in read_fastq(path):
            rid = h.split(None, 1)[0]
            if rid.endswith(b"/1") or rid.endswith(b"/2"):
                rid = rid[:-2]
            yield kind, rid.decode(), s


# --------------------------------------------------------------------------- #
# scoring
# --------------------------------------------------------------------------- #

def hamming(a: bytes, b: bytes) -> int:
    return sum(x != y for x, y in zip(a, b))


def edit_distance(a: bytes, b: bytes) -> int:
    """Levenshtein, vectorised per row.

    The insertion term is a running minimum rather than a loop:
    ``cur[j] = min_k<=j (tmp[k] + (j - k))``, i.e. ``min(tmp[k] - k) + j``, so a row is
    a handful of numpy calls instead of `len(b)` Python steps.
    """
    if a == b:
        return 0
    if not a or not b:
        return max(len(a), len(b))
    barr = np.frombuffer(b, dtype=np.uint8)
    j = np.arange(len(b) + 1, dtype=np.int32)
    prev = j.copy()
    for i in range(1, len(a) + 1):
        tmp = np.empty(len(b) + 1, dtype=np.int32)
        tmp[0] = i
        np.minimum(prev[:-1] + (barr != a[i - 1]), prev[1:] + 1, out=tmp[1:])
        run = np.minimum.accumulate(tmp - j)
        prev = np.minimum(tmp, run + j)
    return int(prev[-1])


def best_offset(emitted: bytes, frag: bytes):
    """``(offset, mismatches)`` of the best ungapped placement of *emitted* in *frag*.

    A wrong-length merge inside a repeat is usually real fragment sequence at the wrong
    place, and the offset says which place — a far more direct diagnosis than an edit
    distance, which cannot tell a shift from a pile of substitutions.
    """
    n, m = len(emitted), len(frag)
    if n == 0 or n > m:
        return 0, hamming(emitted, frag)
    e = np.frombuffer(emitted, dtype=np.uint8)
    f = np.frombuffer(frag, dtype=np.uint8)
    best, best_o = None, 0
    for o in range(m - n + 1):
        d = int(np.count_nonzero(f[o:o + n] != e))
        if best is None or d < best:
            best, best_o = d, o
            if d == 0:
                break
    return best_o, best


def overlap_span(frag_len: int, readlen: int):
    """Fragment positions covered by both mates, as ``[start, end)``."""
    return max(0, frag_len - readlen), min(frag_len, readlen)


def baseline_and_floor(truth, i, readlen):
    """``(r1_wins_errors, oracle_floor_errors)`` for one pair, from the sidecar alone.

    * R1-wins: the overlap comes from R1 unconditionally, which is what a merger with no
      quality model does, and what `zna merge` degenerates to on a constant quality
      string.
    * Oracle floor: the best any consensus can do. It is **not zero** — where both mates
      are wrong the position is unrecoverable, and where only one mate covers a position
      there is nothing to vote against.
    """
    L = int(truth.frag_len[i])
    e1, e2 = truth.errors(i)
    ov0, ov1 = overlap_span(L, readlen)
    cov1 = min(L, readlen)

    base = 0
    for p in e1:
        if p < cov1:
            base += 1
    for p in e2:                                # fragment position of R2 read pos p
        jf = L - 1 - p
        if readlen <= jf < L:
            base += 1

    # fragment positions each mate got wrong, restricted to what it covers
    w1 = {p for p in e1 if p < cov1}
    w2 = {L - 1 - p for p in e2 if 0 <= L - 1 - p < L and L - 1 - p >= ov0}
    floor = 0
    for jf in w1:
        if jf < ov0 or (jf in w2):
            floor += 1
    for jf in w2:
        if jf >= ov1:
            floor += 1
    return base, floor


class ToolScore:
    """Everything measured for one tool, plus the rows for its errors TSV.

    Scoring is two-phase. The first phase streams the tool's output and decides each
    pair's category; the second re-reads the input FASTQ for the pairs that went wrong
    and attaches the evidence, so the diagnosis is taken from the bytes the tool was
    actually shown rather than from a reconstruction.
    """

    def __init__(self, name, truth, readlen, row_cap, max_edit, min_read_length,
                 min_trim_overlap=0):
        self.name = name
        self.truth = truth
        self.readlen = readlen
        self.row_cap = row_cap
        self.max_edit = max_edit
        self.min_read_length = min_read_length
        #: Shortest overlap that can clear --threshold-trim, i.e. the shortest one the
        #: tool is *able* to remove. Misses below it are the specification, not defects.
        self.min_trim_overlap = min_trim_overlap
        n = truth.n
        self.state = np.zeros(n, dtype=np.uint8)          # 1 merged, 2 r1, 4 r2
        self.merged_len = np.zeros(n, dtype=np.int32)
        self.r1_len = np.zeros(n, dtype=np.int32)
        self.r2_len = np.zeros(n, dtype=np.int32)
        self.merged_by_ovl = np.zeros(readlen + 2, dtype=np.int64)
        self.merged_rt_by_len = np.zeros(readlen + 2, dtype=np.int64)
        self.exact = 0
        self.wrong_length = 0
        self.frame_violations = 0
        self.c1_violations = 0
        self.chimeras = 0
        self.len_err = {}
        self.base_err_in_ovl = 0
        self.base_err_out_ovl = 0
        self.baseline_err = 0
        self.floor_err = 0
        self.scored_merges = 0
        self.consensus_hurt = 0          # records left worse than plain R1-wins
        self.consensus_hurt_bases = 0
        self.argmax_checked = 0
        self.argmax_below_truth = 0
        self.argmax_margin = []
        self.unknown_ids = 0
        self.cat_counts = {}
        self.cat_written = {}
        self.pending = []            # (pair index, category, emitted seq or None)
        self.rows = []
        self.identity = []           # (category, matched, span) for every wrong merge
        self._edits = 0

    # -- phase 1: per record ---------------------------------------------- #

    def add(self, kind, rid, seq):
        i = self.truth.index.get(rid)
        if i is None:
            self.unknown_ids += 1
            return
        self._check_c1(i, kind, seq)
        if kind == MERGED:
            self.state[i] |= 1
            self.merged_len[i] = len(seq)
            self._score_merged(i, seq)
        elif kind == R1:
            self.state[i] |= 2
            self.r1_len[i] = len(seq)
        else:
            self.state[i] |= 4
            self.r2_len[i] = len(seq)

    #: Prefix compared for C1, and the mismatches tolerated in it. Chance alignment
    #: gives ~18 of 24; a real 5' end gives 0 or 1, since the expected number of
    #: sequencing errors in 24 bases is ~0.05. Anything in between does not occur.
    C1_PREFIX = 24
    C1_TOLERANCE = 5

    def _check_c1(self, i, kind, seq):
        """Contract C1: base 0 of every emitted read is a true fragment boundary.

        R1 and a merged record start at the fragment's 5' end; R2 starts at its 3' end.
        This is the check the whole ZNA fragment-boundary contract rests on, and it is
        the one that has only ever been made against the tool's own inferences — so it
        is made here over *every* emitted record, including the ones whose geometry went
        wrong, which is exactly where a 5' shift would hide.
        """
        t = self.truth
        frag = t.frag[i]
        L = len(frag)
        k = min(self.C1_PREFIX, L, len(seq))
        if k < 8:
            return                       # too short to distinguish; frag_min makes this
        exp = frag[:k] if kind != R2 else rc(frag[-k:])
        if seq[:k] == exp:
            return
        if sum(a != b for a, b in zip(seq[:k], exp)) > self.C1_TOLERANCE:
            self.c1_violations += 1
            self._note(i, "c1_violation", seq if kind == MERGED else None)

    def _score_merged(self, i, seq):
        t = self.truth
        L = int(t.frag_len[i])
        frag = t.frag[i]
        ovl = int(t.true_ovl[i])
        self.scored_merges += 1
        if ovl == 0:
            self.merged_by_ovl[0] += 1
        elif t.read_through[i]:
            self.merged_rt_by_len[min(L, self.readlen + 1)] += 1
        else:
            self.merged_by_ovl[min(ovl, self.readlen + 1)] += 1

        if len(seq) == L:
            if seq == frag:
                self.exact += 1
            ov0, ov1 = overlap_span(L, self.readlen)
            nin = nout = 0
            for j in range(L):
                if seq[j] != frag[j]:
                    if ov0 <= j < ov1:
                        nin += 1
                    else:
                        nout += 1
            self.base_err_in_ovl += nin
            self.base_err_out_ovl += nout
            base, floor = baseline_and_floor(t, i, self.readlen)
            self.baseline_err += base
            self.floor_err += floor
            # A quality-aware consensus is only worth its complexity if it rarely makes
            # a record WORSE than doing nothing, so count that directly rather than
            # inferring it from the aggregate.
            if nin + nout > base:
                self.consensus_hurt += 1
                self.consensus_hurt_bases += nin + nout - base
            # A frame violation would show as a mismatch rate near chance rather than
            # near the error rate; 25% is far above anything sequencing error can reach
            # and far below the 75% of a wrong placement.
            if nin + nout > 5 + 0.25 * L:
                self.frame_violations += 1
                self._note(i, "frame_violation", seq)
            elif nin + nout > floor:
                self._note(i, "consensus_miss", seq)
            return

        self.wrong_length += 1
        d = len(seq) - L
        self.len_err[d] = self.len_err.get(d, 0) + 1
        if ovl == 0:
            self.chimeras += 1
            self._note(i, "chimera", seq)
        else:
            self._note(i, "wrong_length", seq)

    def _note(self, i, category, seq, extra=None):
        """Record one thing the tool got wrong.

        *extra* is ``(emitted_len, len_err)`` for rows whose sequence is not held — a
        trim is described by how much came off, not by the bases that remain.
        """
        self.cat_counts[category] = self.cat_counts.get(category, 0) + 1
        # Wrong merges are always kept, capped or not: their evidence is what turns a
        # chimera count into an explanation, and there should never be many of them.
        detailed = category in ("chimera", "wrong_length", "chimera_dropped",
                                "wrong_length_dropped", "frame_violation",
                                "false_trim", "r1_shortened")
        if detailed or self.cat_written.get(category, 0) < self.row_cap:
            self.pending.append((i, category, seq, extra))

    # -- phase 1b: the merges that were filtered away ---------------------- #

    def find_dropped_merges(self):
        """Pairs the tool emitted nothing for, which can only be a filtered merge.

        Both mates are full-length here, so an unmerged pair always clears
        `--min-read-length` and always emits two records. A pair with no output was
        therefore merged into a record shorter than the filter and dropped — a wrong
        merge that is invisible in the output file, and one that removes the fragment
        from the corpus entirely. It is counted, not ignored.
        """
        t = self.truth
        gone = np.flatnonzero(self.state == 0)
        n = 0
        for i in gone:
            if t.len1[i] < self.min_read_length or t.len2[i] < self.min_read_length:
                continue                      # not attributable to a merge; leave it
            n += 1
            self.wrong_length += 1
            if t.true_ovl[i] == 0:
                self.chimeras += 1
                self._note(int(i), "chimera_dropped", None)
            else:
                self._note(int(i), "wrong_length_dropped", None)
        self.n_merged_dropped = n
        return n

    # -- phase 2: evidence ------------------------------------------------- #

    def write_rows(self, reads, scan):
        """Attach the scan's own evidence to each error row.

        `scan` is `zna merge`'s overlap kernel, run on the pair the tool was given, so
        `scan_*` says what the shipped scoring rule sees: the shift it picks, the
        evidence in bits, and how well the two reads actually agree there. A chimera
        with 90% identity over 89 bases is the genome repeating, which is a different
        finding from a merger inventing an alignment — and the columns say which.
        """
        t = self.truth
        for i, category, seq, extra in self.pending:
            L = int(t.frag_len[i])
            frag = t.frag[i]
            r1, r2 = reads[i]
            shift, score_q, olen, mism = scan(r1, r2)
            ident = f"{olen - mism}/{olen}" if olen else "."
            if olen:
                self.identity.append((category, olen - mism, olen))
            # The evidence at the shift the truth says is right, for comparison.
            true_bits = "."
            if t.true_ovl[i] > 0 and category.startswith("wrong_length"):
                at = scan.score_at(r1, r2, L - int(t.len2[i]))
                if at is not None:
                    true_bits = round(at[0] / SCALE, 2)
                    self.argmax_checked += 1
                    # The scan's contract is "the best shift at or above the floor, else
                    # nothing", so a miss is only a defect when the true shift clears
                    # the floor. Comparing against a `no overlap found` zero would
                    # otherwise report the floor itself as a search failure.
                    missed = (score_q < at[0]) if olen else (at[0] >= scan.floor_q)
                    if missed:
                        self.argmax_below_truth += 1
                    elif olen:
                        self.argmax_margin.append((score_q - at[0]) / SCALE)
            if self.cat_written.get(category, 0) >= self.row_cap:
                continue
            self.cat_written[category] = self.cat_written.get(category, 0) + 1
            if seq is None:                       # a trim row, or a filtered-away merge
                mm, ed, off, shown = ".", ".", ".", "."
                emitted_len, len_err = extra if extra else (".", ".")
            else:
                emitted_len, len_err = len(seq), len(seq) - L
                if len(seq) == L:
                    off, mm = 0, hamming(seq, frag)
                else:
                    off, mm = best_offset(seq, frag)
                ed = -1
                if self._edits < self.max_edit:
                    ed = edit_distance(seq, frag)
                    self._edits += 1
                shown = seq.decode()
            self.rows.append("\t".join(str(x) for x in (
                t.index_id[i], category, t.chrom[i], t.start[i], t.strand[i], L,
                int(t.true_ovl[i]), int(t.read_through[i]), int(t.n_err1[i]),
                int(t.n_err2[i]), emitted_len, len_err, mm, ed,
                shift, round(score_q / SCALE, 2), olen, ident, true_bits, off, shown,
                frag.decode())))

    # -- roll-up ---------------------------------------------------------- #

    def finish(self):
        merged = (self.state & 1) != 0
        both = (self.state & 6) == 6
        orphan = ((self.state & 6) != 0) & ((self.state & 6) != 6)
        self.n_merged = int(merged.sum())
        self.n_pairs_kept = int(both.sum())
        self.n_orphans = int(orphan.sum())
        self.n_no_output = int((self.state == 0).sum())
        # A pair emitted BOTH as a merged single and as a mate pair would double the
        # molecule. Structurally impossible in either tool; asserted rather than assumed.
        self.n_merged_and_paired = int((merged & ((self.state & 6) != 0)).sum())

        self._score_trims(both)
        return self

    # -- the trim band, contract C4 ---------------------------------------- #

    def _score_trims(self, both):
        """What the *unmerged* pairs cost the corpus.

        This half of the tool matters as much as merging to a model trained on the
        output: a pair that does not merge is still encoded, and if its redundant
        overlap survives, the same physical bases appear twice in the corpus with no
        marker saying so. `zna merge` therefore cuts the overlap off R2's **3'** end,
        leaving R1 and R2 tiling the fragment exactly once.

        Scored over **every** kept pair, in three regimes, because the failure modes are
        opposite and averaging them hides both:

        * `true_ovl == 0` — the mates share nothing, so any trim **deletes real
          sequence**. This is the specificity side, and it is the one an analysis
          restricted to overlapping pairs cannot see at all.
        * `0 < true_ovl`, no read-through — exactly `true_ovl` bases are redundant and
          exactly that many should come off. Under-trimming duplicates, over-trimming
          deletes.
        * read-through — both mates carry the *whole* fragment plus adapter. There is no
          trim that fixes this (the geometry is a negative shift, which the trim branch
          does not take); it has to merge, and a kept read-through pair puts the fragment
          in twice, with adapter.

        The two numbers to carry away are the corpus-level ones: bases duplicated, and
        real bases deleted.
        """
        t = self.truth
        RL = self.readlen
        # The trim is SYMMETRIC: the overlap sits at the 3' end of both mates and is
        # split between them, so "how much was removed" is the sum over the pair.
        # Charging only R2 -- which was right when only R2 was cut -- reads a correct
        # balanced trim as half an under-trim, and every trimmed R1 as a violation.
        cut1 = np.where(both, t.len1 - self.r1_len, 0)
        cut2 = np.where(both, t.len2 - self.r2_len, 0)
        removed = cut1 + cut2
        keep2 = np.where(both, self.r2_len, 0)
        # Neither mate may grow, and nothing may come off a 5' end -- the latter is what
        # the C1 prefix check verifies directly, on every emitted record.
        self.r1_shortened = int((both & (self.r1_len > t.len1)).sum())
        # How unequal the emitted pair ended up: the reason for splitting at all.
        tr = both & (removed > 0)
        self.trimmed_pairs = int(tr.sum())
        self.trim_length_gap = (int(np.abs(self.r1_len - self.r2_len)[tr].sum())
                                if tr.any() else 0)
        self.trim_max_gap = int(np.abs(self.r1_len - self.r2_len)[tr].max()) if tr.any() else 0

        no_ovl = both & (t.true_ovl == 0)
        normal = both & (t.read_through == 0) & (t.true_ovl > 0)
        rthru = both & (t.read_through == 1)
        self.kept_no_overlap = int(no_ovl.sum())
        self.kept_overlapping = int(normal.sum())
        self.kept_read_through = int(rthru.sum())

        # (a) nothing to remove: every removed base is real sequence destroyed
        self.trim_false = int((no_ovl & (removed > 0)).sum())
        self.trim_false_bases = int(removed[no_ovl].sum())

        # (b) the trim band proper
        want = t.true_ovl
        self.trim_exact = int((normal & (removed == want)).sum())
        self.trim_over = int((normal & (removed > want)).sum())
        self.trim_under = int((normal & (removed < want) & (removed > 0)).sum())
        self.trim_none = int((normal & (removed == 0)).sum())
        # Overlaps too short to clear --threshold-trim cannot be removed by design, so
        # split the misses: below the floor is the specification, above it is a miss.
        floor = self.min_trim_overlap
        self.trim_none_below_floor = int((normal & (removed == 0)
                                          & (t.true_ovl < floor)).sum())
        self.trim_none_above_floor = self.trim_none - self.trim_none_below_floor

        # (c) read-through kept whole: the entire fragment is in the corpus twice
        rt_dup = np.where(rthru, np.minimum(t.frag_len, keep2), 0)
        rt_adapter = (np.where(rthru, RL - t.frag_len, 0)
                      + np.where(rthru, np.clip(keep2 - t.frag_len, 0, None), 0))

        # corpus-level totals
        dup = np.where(normal, np.clip(want - removed, 0, None), 0) + rt_dup
        deleted = (np.where(normal, np.clip(removed - want, 0, None), 0)
                   + np.where(no_ovl, removed, 0)
                   + np.where(rthru, np.clip(t.frag_len - keep2, 0, None), 0))
        self.bases_duplicated = int(dup.sum())
        self.bases_deleted = int(deleted.sum())
        # The counterfactual, which is what makes the trade legible: how much duplicated
        # sequence would have reached the corpus had nothing been trimmed. The trim is
        # worth its cost exactly insofar as this exceeds `bases_deleted`.
        self.bases_duplicated_untrimmed = int(
            (np.where(normal, t.true_ovl, 0) + np.where(rthru, t.frag_len, 0)).sum())
        self.bases_adapter_left = int(rt_adapter.sum())
        self.bases_emitted_unmerged = int((self.r1_len[both].sum()
                                           + self.r2_len[both].sum()))

        # the sensitivity curve for the trim, the analogue of the merge one
        self.trim_exact_by_ovl = np.zeros(RL + 2, dtype=np.int64)
        sel = np.flatnonzero(normal & (removed == want))
        if sel.size:
            np.add.at(self.trim_exact_by_ovl, np.minimum(t.true_ovl[sel], RL + 1), 1)

        for i in np.flatnonzero(no_ovl & (removed > 0)):
            self._note(int(i), "false_trim", None, (int(keep2[i]), int(removed[i])))
        for i in np.flatnonzero(normal & (removed > want)):
            self._note(int(i), "over_trim", None,
                       (int(keep2[i]), int(removed[i] - want[i])))
        for i in np.flatnonzero(normal & (removed == 0) & (t.true_ovl >= floor)):
            self._note(int(i), "missed_trim", None, (int(keep2[i]), -int(want[i])))
        for i in np.flatnonzero(both & (self.r1_len > t.len1)):
            self._note(int(i), "read_grew", None,
                       (int(self.r1_len[i]), int(self.r1_len[i] - t.len1[i])))


# --------------------------------------------------------------------------- #
# running the tools
# --------------------------------------------------------------------------- #

def run(cmd, log):
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True)
    dt = time.perf_counter() - t0
    log.write_bytes(proc.stdout + b"\n----- stderr -----\n" + proc.stderr)
    if proc.returncode:
        raise SystemExit(f"{cmd[0]} failed (exit {proc.returncode}); see {log}\n"
                         + proc.stderr.decode(errors="replace")[-2000:])
    return dt


def fetch_reads(in1, in2, wanted, truth):
    """The exact input reads for a set of pair indices, straight from the FASTQ.

    The sidecar reconstructs every fragment-derived base, but not the random filler past
    the adapter on a short read-through fragment — and an error row is exactly where a
    reconstruction should not be trusted. One extra pass over the input is cheap next to
    getting the evidence from the same bytes the tools were given.
    """
    want = {truth.index_id[i].encode(): i for i in wanted}
    out = {i: [None, None] for i in wanted}
    for path, slot in ((in1, 0), (in2, 1)):
        for h, s in read_fastq(path):
            rid = h.split(None, 1)[0]
            if rid.endswith(b"/1") or rid.endswith(b"/2"):
                rid = rid[:-2]
            i = want.get(rid)
            if i is not None:
                out[i][slot] = s
    missing = [i for i, v in out.items() if v[0] is None or v[1] is None]
    if missing:
        raise SystemExit(f"{len(missing)} error pairs not found in the input FASTQs")
    return {i: (a, b) for i, (a, b) in out.items()}


def make_scan(t_merge, t_trim):
    """`zna merge`'s own overlap kernel, at the thresholds the run actually used.

    Returns ``scan(r1, r2)`` and ``score_at(r1, r2, shift)``. The second one is what
    separates a *defective* scan from an *ambiguous* input: score the shift the truth
    says is right, and compare. If the tool's pick ever scores lower than the truth's,
    the argmax or its pruning is broken; if it always scores higher, the tool is
    maximising correctly and the sequence really does align better somewhere else.
    """
    from zna.merge import backend as _backend
    from zna.merge.params import MergeParams
    p = MergeParams(t_merge=t_merge, t_trim=t_trim)
    kern = _backend.active().scan

    def scan(r1, r2):
        s2rc = rc(r2)
        return kern(r1, s2rc, len(r1), len(s2rc), p.match_q, p.step_q, p.t_trim_q)

    def score_at(r1, r2, shift):
        """``(score_q, overlap_len, mismatches)`` at one specific shift, or None."""
        s2rc = rc(r2)
        a0, b0 = max(shift, 0), max(-shift, 0)
        n = min(len(r1) - a0, len(s2rc) - b0)
        if n <= 0:
            return None
        d = sum(x != y for x, y in zip(r1[a0:a0 + n], s2rc[b0:b0 + n]))
        return n * p.match_q - d * p.step_q, n, d

    scan.score_at = score_at
    scan.floor_q = p.t_trim_q
    return scan


# --------------------------------------------------------------------------- #
# the report
# --------------------------------------------------------------------------- #

BINS = [(0, 0), (1, 4), (5, 9), (10, 14), (15, 19), (20, 29), (30, 49),
        (50, 99), (100, 149), (150, 10 ** 9)]


def sensitivity_table(truth, scores, readlen):
    normal = truth.read_through == 0
    lines = ["| true overlap | pairs | " + " | ".join(
        f"{s.name} merged" for s in scores) + " |",
        "|---|---:|" + "---:|" * len(scores)]
    for lo, hi in BINS:
        sel = normal & (truth.true_ovl >= lo) & (truth.true_ovl <= hi)
        n = int(sel.sum())
        if not n:
            continue
        cells = []
        for s in scores:
            m = int(s.merged_by_ovl[lo:min(hi, readlen + 1) + 1].sum())
            cells.append(f"{m:,} ({100.0 * m / n:.2f}%)")
        label = f"{lo}" if lo == hi else (f"{lo}+" if hi > readlen else f"{lo}–{hi}")
        lines.append(f"| {label} | {n:,} | " + " | ".join(cells) + " |")
    return "\n".join(lines)


def read_through_table(truth, scores, readlen):
    rt = truth.read_through == 1
    edges = [(1, 39), (40, 59), (60, 79), (80, 99), (100, 119), (120, 149)]
    lines = ["| fragment length | pairs | " + " | ".join(
        f"{s.name} merged" for s in scores) + " |",
        "|---|---:|" + "---:|" * len(scores)]
    for lo, hi in edges:
        sel = rt & (truth.frag_len >= lo) & (truth.frag_len <= hi)
        n = int(sel.sum())
        if not n:
            continue
        cells = []
        for s in scores:
            m = int(s.merged_rt_by_len[lo:hi + 1].sum())
            cells.append(f"{m:,} ({100.0 * m / n:.2f}%)")
        lines.append(f"| {lo}–{hi} | {n:,} | " + " | ".join(cells) + " |")
    return "\n".join(lines)


def trim_curve_table(truth, scores, readlen):
    """Exact-trim rate by true overlap, over pairs the tool did NOT merge.

    The denominator is per tool on purpose: a pair one tool merged and the other kept is
    not a trim failure for the one that merged it, and pooling them would read as one.
    """
    lines = ["### Exact-trim rate by true overlap, over each tool's own kept pairs", "",
             "| true overlap | " + " | ".join(
                 f"{s.name} kept | {s.name} trimmed exactly" for s in scores) + " |",
             "|---|" + "---:|---:|" * len(scores)]
    normal = truth.read_through == 0
    for lo, hi in [(1, 4), (5, 9), (10, 14), (15, 19), (20, 29), (30, 10 ** 9)]:
        band = normal & (truth.true_ovl >= lo) & (truth.true_ovl <= hi)
        if not band.any():
            continue
        cells = []
        for s in scores:
            kept = int((band & ((s.state & 6) == 6)).sum())
            ok = int(s.trim_exact_by_ovl[lo:min(hi, readlen + 1) + 1].sum())
            cells.append(f"{kept:,}")
            cells.append(f"{ok:,} ({100.0 * ok / kept:.1f}%)" if kept else "–")
        label = f"{lo}–{hi}" if hi <= readlen else f"{lo}+"
        lines.append(f"| {label} | " + " | ".join(cells) + " |")
    return "\n".join(lines)


def evidence_table(scores):
    lines = ["| tool | wrong merges | median identity | min identity | ≥80% identical "
             "| median overlap |", "|---|---:|---:|---:|---:|---:|"]
    for s in scores:
        ident = [(m, sp) for c, m, sp in s.identity
                 if c.startswith("chimera") or c.startswith("wrong_length")]
        if not ident:
            lines.append(f"| {s.name} | 0 | – | – | – | – |")
            continue
        fr = sorted(m / sp for m, sp in ident)
        sp = sorted(x for _, x in ident)
        lines.append(
            f"| {s.name} | {len(fr):,} | {fr[len(fr) // 2]:.1%} | {fr[0]:.1%} | "
            f"{sum(1 for x in fr if x >= 0.80):,} | {sp[len(sp) // 2]} |")
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sim-prefix", required=True,
                    help="the --out-prefix given to simulate.py")
    ap.add_argument("--out", required=True, help="results directory")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--min-read-length", type=int, default=40)
    ap.add_argument("--threshold-merge", type=float, default=None,
                    help="passed to zna merge; default is the tool's own")
    ap.add_argument("--threshold-trim", type=float, default=None,
                    help="passed to zna merge; sweep this to price the trim band")
    ap.add_argument("--zna", default=None, help="zna executable (default: on PATH)")
    ap.add_argument("--fastp", default="/Users/mkiyer/sw/miniforge3/envs/fastp/bin/fastp")
    ap.add_argument("--row-cap", type=int, default=DEFAULT_ROW_CAP,
                    help="max rows per category in the errors TSV")
    ap.add_argument("--max-edit", type=int, default=5000,
                    help="edit distances to compute (they are O(n*m); -1 past this)")
    ap.add_argument("--skip-run", action="store_true",
                    help="score tool output already present in --out")
    args = ap.parse_args(argv)

    pre = Path(args.sim_prefix)
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    meta = json.loads(Path(f"{pre}.json").read_text())
    readlen = meta["read_length"]

    zna_exe = args.zna or shutil.which("zna")
    if not zna_exe:
        raise SystemExit("no `zna` on PATH; pass --zna")
    if not Path(args.fastp).exists():
        raise SystemExit(f"fastp not found at {args.fastp}")

    in1, in2 = f"{pre}_1.fq.gz", f"{pre}_2.fq.gz"
    zna_out = outdir / "zna.fq.gz"
    fp_m = outdir / "fastp.merged.fq.gz"
    fp_1 = outdir / "fastp.un1.fq.gz"
    fp_2 = outdir / "fastp.un2.fq.gz"

    zna_cmd = [zna_exe, "merge", "--in1", in1, "--in2", in2, "--out", str(zna_out),
               "--json", str(outdir / "zna.json"),
               "--min-read-length", str(args.min_read_length),
               "--threads", str(args.threads), "-q"]
    if args.threshold_merge is not None:
        zna_cmd += ["--threshold-merge", str(args.threshold_merge)]
    if args.threshold_trim is not None:
        zna_cmd += ["--threshold-trim", str(args.threshold_trim)]
    # fastp does quality filtering, polyG trimming and adapter trimming by default and
    # `zna merge` does none of them, so they are turned off: what is being compared is
    # merging, not preprocessing. `-A` was checked empirically not to change the merge
    # count on read-through pairs (fastp's PE adapter handling is overlap-based, and a
    # merged pair takes the merge path instead). Base correction in the overlap is NOT
    # opt-in here: fastp reports corrected bases in merge mode without `-c`, so both
    # tools are doing consensus.
    fastp_cmd = [args.fastp, "--in1", in1, "--in2", in2,
                 "--merge", "--merged_out", str(fp_m),
                 "--out1", str(fp_1), "--out2", str(fp_2),
                 "-A", "-Q", "-G", "--dont_eval_duplication",
                 "--length_required", str(args.min_read_length),
                 "--json", str(outdir / "fastp.json"), "--html", "/dev/null",
                 "-w", str(min(16, max(1, args.threads)))]

    t_zna = t_fp = None
    if not args.skip_run:
        print("running zna merge ...", file=sys.stderr)
        t_zna = run(zna_cmd, outdir / "zna.log")
        print("running fastp ...", file=sys.stderr)
        t_fp = run(fastp_cmd, outdir / "fastp.log")

    print("loading truth ...", file=sys.stderr)
    truth = Truth(f"{pre}.truth.tsv")

    # The shortest overlap the trim can act on at all, taken from the run's own
    # parameters rather than recomputed: a miss below it is the specification.
    zj = json.loads((outdir / "zna.json").read_text())
    zp = zj["params"]
    min_trim_overlap = -(-zp["threshold_trim_q"] // zp["match_q"])

    scores = []
    for name, records in (("zna", zna_records(zna_out)),
                          ("fastp", fastp_records(fp_m, fp_1, fp_2))):
        print(f"scoring {name} ...", file=sys.stderr)
        sc = ToolScore(name, truth, readlen, args.row_cap, args.max_edit,
                       args.min_read_length, min_trim_overlap)
        for kind, rid, seq in records:
            sc.add(kind, rid, seq)
        sc.finish()
        sc.find_dropped_merges()
        scores.append(sc)

    wanted = sorted({i for sc in scores for i, _, _, _ in sc.pending})
    if wanted:
        print(f"re-reading {len(wanted):,} error pairs from the input ...",
              file=sys.stderr)
        reads = fetch_reads(in1, in2, wanted, truth)
        scan = make_scan(zp["threshold_merge_bits"], zp["threshold_trim_bits"])
        for sc in scores:
            sc.write_rows(reads, scan)

    for sc, tsec in zip(scores, (t_zna, t_fp)):
        sc.elapsed = tsec
    write_outputs(outdir, meta, truth, scores, readlen, zna_cmd, fastp_cmd, args)
    print(f"wrote {outdir}/report.md, summary.json, zna_errors.tsv", file=sys.stderr)
    return 0


def summarise(sc, truth, readlen):
    n = truth.n
    no_ovl = int((truth.true_ovl == 0).sum())
    mergeable = int((truth.true_ovl >= 15).sum())
    got = int(sc.merged_by_ovl[15:].sum() + sc.merged_rt_by_len[15:].sum())
    ident = [(m, s) for c, m, s in sc.identity
             if c.startswith("chimera") or c.startswith("wrong_length")]
    total_merged = sc.n_merged + sc.n_merged_dropped
    d = {
        "input_pairs": n,
        "merged": total_merged,
        "merged_emitted": sc.n_merged,
        "merged_then_dropped_below_min_length": sc.n_merged_dropped,
        "merged_pct": round(100.0 * total_merged / n, 3),
        "pairs_kept_unmerged": sc.n_pairs_kept,
        "orphans": sc.n_orphans,
        "pairs_with_no_output": sc.n_no_output,
        "merged_and_paired": sc.n_merged_and_paired,
        "unknown_read_ids": sc.unknown_ids,
        "chimeras": sc.chimeras,
        "chimera_rate_at_zero_overlap": round(sc.chimeras / no_ovl, 8) if no_ovl else 0.0,
        "sensitivity_at_overlap_ge_15": round(got / mergeable, 6) if mergeable else 0.0,
        "merged_exact_fragment": sc.exact,
        "merged_exact_pct": round(100.0 * sc.exact / sc.n_merged, 3) if sc.n_merged else 0.0,
        "merged_wrong_length": sc.wrong_length,
        "boundary_violations_c2_wrong_length": sc.wrong_length,
        "boundary_violations_c1_base_zero": sc.c1_violations,
        "boundary_violations_c1_frame": sc.frame_violations,
        "records_checked_for_c1": int(sc.n_merged + 2 * sc.n_pairs_kept + sc.n_orphans),
        "base_errors_in_overlap": sc.base_err_in_ovl,
        "base_errors_outside_overlap": sc.base_err_out_ovl,
        "base_errors_total": sc.base_err_in_ovl + sc.base_err_out_ovl,
        "baseline_r1_wins_errors": sc.baseline_err,
        "oracle_floor_errors": sc.floor_err,
        # Where the tool sits between "no quality model" and "the best possible".
        # 100% means every recoverable error in the overlap was recovered.
        "consensus_recovery_pct": round(
            100.0 * (sc.baseline_err - sc.base_err_in_ovl - sc.base_err_out_ovl)
            / (sc.baseline_err - sc.floor_err), 2)
        if sc.baseline_err > sc.floor_err else None,
        "consensus_made_worse_records": sc.consensus_hurt,
        "consensus_made_worse_bases": sc.consensus_hurt_bases,
        # 0 means the scan is maximising correctly and every wrong merge is a genuinely
        # better-scoring alternative alignment, not a search defect.
        "wrong_shifts_scoring_below_the_truth": sc.argmax_below_truth,
        "wrong_shifts_checked_against_the_truth": sc.argmax_checked,
        "wrong_shift_median_margin_bits": round(
            sorted(sc.argmax_margin)[len(sc.argmax_margin) // 2], 2)
        if sc.argmax_margin else None,
        "correct_length_merges_scored": sc.scored_merges - (sc.wrong_length
                                                            - sc.n_merged_dropped),
        # --- the trim band (contract C4) --------------------------------- #
        "kept_pairs_no_true_overlap": sc.kept_no_overlap,
        "kept_pairs_overlapping": sc.kept_overlapping,
        "kept_pairs_read_through": sc.kept_read_through,
        "min_trimmable_overlap": sc.min_trim_overlap,
        "trim_exact": sc.trim_exact,
        "trim_over": sc.trim_over,
        "trim_under": sc.trim_under,
        "trim_absent": sc.trim_none,
        "trim_absent_below_threshold": sc.trim_none_below_floor,
        "trim_absent_above_threshold": sc.trim_none_above_floor,
        # A trim on a pair that shares nothing destroys real sequence. This is the
        # specificity side of the trim, and it must be zero.
        "false_trims": sc.trim_false,
        "false_trim_bases": sc.trim_false_bases,
        "emitted_read_longer_than_input": sc.r1_shortened,
        "trimmed_pairs": sc.trimmed_pairs,
        # The point of splitting the overlap rather than taking it all off R2: the two
        # emitted mates come out the same length.
        "trimmed_mean_length_gap": round(sc.trim_length_gap / sc.trimmed_pairs, 2)
        if sc.trimmed_pairs else 0.0,
        "trimmed_max_length_gap": sc.trim_max_gap,
        # What the corpus actually receives, in bases.
        "bases_emitted_unmerged": sc.bases_emitted_unmerged,
        "fragment_bases_duplicated": sc.bases_duplicated,
        "fragment_bases_duplicated_if_untrimmed": sc.bases_duplicated_untrimmed,
        "fragment_bases_removed_correctly": sc.bases_duplicated_untrimmed
                                            - sc.bases_duplicated,
        "fragment_bases_deleted": sc.bases_deleted,
        "trim_benefit_ratio": round(
            (sc.bases_duplicated_untrimmed - sc.bases_duplicated) / sc.bases_deleted, 2)
        if sc.bases_deleted else None,
        "adapter_bases_left_in_kept_read_through": sc.bases_adapter_left,
        "duplicated_per_10k_unmerged_bases": round(
            1e4 * sc.bases_duplicated / sc.bases_emitted_unmerged, 2)
        if sc.bases_emitted_unmerged else 0.0,
        "error_categories": sc.cat_counts,
        "error_rows_written": sc.cat_written,
    }
    if ident:
        # What every wrong merge was looking at. If these are near-identical over a long
        # stretch, the genome repeats there and no overlap rule could have known better.
        fr = sorted(m / s for m, s in ident)
        sp = sorted(s for _, s in ident)
        d["wrong_merge_evidence"] = {
            "n": len(fr),
            "identity_min": round(fr[0], 4),
            "identity_median": round(fr[len(fr) // 2], 4),
            "identity_ge_0.80": sum(1 for x in fr if x >= 0.80),
            "overlap_len_min": sp[0],
            "overlap_len_median": sp[len(sp) // 2],
        }
    if sc.elapsed:
        d["elapsed_s"] = round(sc.elapsed, 2)
        d["us_per_pair"] = round(1e6 * sc.elapsed / n, 3)
    return d


def write_outputs(outdir, meta, truth, scores, readlen, zna_cmd, fastp_cmd, args):
    summary = {
        "simulation": meta,
        "invocations": {"zna": zna_cmd, "fastp": fastp_cmd},
        "tools": {sc.name: summarise(sc, truth, readlen) for sc in scores},
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    for sc in scores:
        cap_note = "\n".join(
            f"# {c}: {sc.cat_counts[c]:,} total, {sc.cat_written.get(c, 0):,} written"
            for c in sorted(sc.cat_counts))
        head = ("# every pair this tool got wrong, capped per category\n"
                f"{cap_note}\n"
                "# scan_* is zna's own overlap kernel re-run on this pair: the shift it\n"
                "# picks, the evidence in bits, the overlap length, and how well the two\n"
                "# reads agree there. A '.' in emitted_* means the merged record was\n"
                "# below --min-read-length and never written.\n"
                "read_id\tcategory\tchrom\tstart\tstrand\tfrag_len\ttrue_ovl\t"
                "read_through\tn_err1\tn_err2\temitted_len\tlen_err\tn_mismatch\t"
                "edit_distance\tscan_shift\tscan_score_bits\tscan_olen\t"
                "scan_identity\ttrue_shift_score_bits\tbest_offset\temitted_seq\t"
                "true_fragment\n")
        (outdir / f"{sc.name}_errors.tsv").write_text(head + "\n".join(sc.rows) + "\n")

    n = truth.n
    no_ovl = int((truth.true_ovl == 0).sum())

    def pair(key, fmt="{:,}"):
        # A metric can be undefined for a tool rather than zero -- a ratio whose
        # denominator is the thing that never happened -- and "0" would be a lie.
        def cell(v):
            return "–" if v is None else fmt.format(v)
        return f"| {key} | " + " | ".join(
            cell(summary['tools'][s.name][key]) for s in scores) + " |"

    lines = [
        "# `zna merge` vs fastp on simulated ground truth",
        "",
        f"{n:,} pairs, {meta['read_length']} cycles, fragments uniform on "
        f"[{meta['frag_min']}, {meta['frag_max']}], error rate "
        f"{meta['error_rate_realised']} ({meta['quality_model']} qualities), "
        f"genome `{Path(meta['genome']).name}`, seed {meta['seed']}.",
        "",
        "Uniform fragment lengths populate every geometric regime at equal density, so "
        "**the overall merge rate here is not a library number** — read the per-bin "
        "curve instead.",
        "",
        "## Invocations",
        "",
        "```",
        " ".join(zna_cmd),
        "",
        " ".join(fastp_cmd),
        "```",
        "",
        "fastp's quality filtering, polyG trimming and adapter trimming are off, so what "
        "is compared is merging rather than preprocessing; `--length_required` matches "
        "`--min-read-length`. fastp corrects bases in the overlap in merge mode without "
        "`-c`, so both tools are running a consensus.",
        "",
        "## 1–2. Merge sensitivity and specificity, by true overlap",
        "",
        sensitivity_table(truth, scores, readlen),
        "",
        f"The first row is the specificity test: {no_ovl:,} pairs have **no true "
        "overlap at all**, so every merge there is a chimera. This table counts only "
        "merges that were *emitted*; section 3 adds the ones filtered below "
        "`--min-read-length`, so its chimera total is the larger number.",
        "",
        "### Read-through pairs (fragment shorter than the read)",
        "",
        read_through_table(truth, scores, readlen),
        "",
        "## 3. Reconstruction accuracy of merged records",
        "",
        "| metric | " + " | ".join(s.name for s in scores) + " |",
        "|---|" + "---:|" * len(scores),
        pair("merged"),
        pair("merged_emitted"),
        pair("merged_then_dropped_below_min_length"),
        pair("merged_exact_fragment"),
        pair("merged_exact_pct", "{:.3f}%"),
        pair("merged_wrong_length"),
        pair("chimeras"),
        "",
        "`merged_then_dropped_below_min_length` are merges that produced a record under "
        "`--min-read-length` and were filtered away. They are invisible in the output "
        "file but they are still wrong merges, and they delete the fragment from the "
        "corpus rather than emitting it as a pair — so they are counted here, recovered "
        "from the pairs the tool emitted nothing for.",
        "",
        "### What the wrong merges were looking at",
        "",
        evidence_table(scores),
        "",
        "`zna merge`'s own scan, re-run on each pair that merged wrongly: the overlap it "
        "found and the fraction of those bases that actually agree. Near-identity over a "
        "long stretch means the fragment's two ends are genuinely homologous — a real "
        "repeat, which no overlap rule can distinguish from a real overlap.",
        "",
        "### Is the scan wrong, or is the sequence ambiguous?",
        "",
        "For every merge at the wrong shift where a true overlap exists, the shipped "
        "kernel is also asked to score the **true** shift. If its pick ever scored lower "
        "than the truth's, the argmax or its pruning would be defective.",
        "",
        "| metric | " + " | ".join(s.name for s in scores) + " |",
        "|---|" + "---:|" * len(scores),
        pair("wrong_shifts_checked_against_the_truth"),
        pair("wrong_shifts_scoring_below_the_truth"),
        pair("wrong_shift_median_margin_bits"),
        "",
        "The margin is how much *better* the alternative alignment scored than the true "
        "one. Both columns run **zna's** kernel, so the fastp column reads differently: "
        "it is what zna's rule would have seen on the pairs fastp got wrong, and a "
        "margin of 0 there means zna's kernel would have picked the true shift.",
        "",
        "## 4. Base accuracy inside the overlap",
        "",
        "Counted over merged records of the correct length, against the true fragment. "
        "`R1-wins` is what a merger with no quality model would have emitted and "
        "`oracle floor` is the best any consensus could do — a position where **both** "
        "mates are wrong is unrecoverable.",
        "",
        "| metric | " + " | ".join(s.name for s in scores) + " |",
        "|---|" + "---:|" * len(scores),
        pair("correct_length_merges_scored"),
        pair("base_errors_in_overlap"),
        pair("base_errors_outside_overlap"),
        pair("baseline_r1_wins_errors"),
        pair("oracle_floor_errors"),
        pair("consensus_recovery_pct", "{}%"),
        pair("consensus_made_worse_records"),
        pair("consensus_made_worse_bases"),
        "",
        "`consensus_recovery_pct` is where the tool sits between the two: 0% is "
        "R1-wins, 100% is the oracle. Each tool is scored over the records **it** "
        "merged, so the denominators differ where the merge sets do. "
        "`consensus_made_worse_*` counts records the consensus left *worse* than doing "
        "nothing, which is what a quality-aware rule has to keep rare to be worth its "
        "complexity.",
        "",
        "## 5. Boundary violations (contract C1/C2)",
        "",
        "C1 — *base 0 of every emitted read is a true fragment boundary* — is checked "
        "over **every** emitted record, merged or not, by comparing its first "
        f"{ToolScore.C1_PREFIX} bases with the fragment end it claims to start at. A "
        "real 5' end mismatches in 0 or 1 of them; a shifted one mismatches in ~18. "
        "Checking the wrongly-merged records too is the point: a 5' shift would hide "
        "exactly there.",
        "",
        "C2 — *a merged record is its fragment exactly*. With full-length mates the "
        "emitted length determines the inferred shift, so a length mismatch is exactly "
        "a wrong inference, and a pair with no true overlap can never reach the true "
        "length by accident.",
        "",
        "| metric | " + " | ".join(s.name for s in scores) + " |",
        "|---|" + "---:|" * len(scores),
        pair("records_checked_for_c1"),
        pair("boundary_violations_c1_base_zero"),
        pair("boundary_violations_c1_frame"),
        pair("boundary_violations_c2_wrong_length"),
        "",
        "## 6. Unmerged handling (contract C3)",
        "",
        "| metric | " + " | ".join(s.name for s in scores) + " |",
        "|---|" + "---:|" * len(scores),
        pair("pairs_kept_unmerged"),
        pair("orphans"),
        pair("merged_and_paired"),
        pair("pairs_with_no_output"),
        "",
        "`orphans` is contract C3: one mate emitted without the other, which would be "
        "encoded as a spurious full molecule. `merged_and_paired` would be the same "
        "molecule emitted twice. Both must be zero. `pairs_with_no_output` are the "
        "filtered-away merges of section 3, not lost pairs.",
        "",
        "## 7. The trim band (contract C4), and what the corpus receives",
        "",
        "A pair that does not merge is **still encoded**, so its redundant overlap "
        "matters as much as a merge does: leave it in and the same physical bases enter "
        "the corpus twice with nothing marking them as one molecule; cut too much and "
        "real sequence is destroyed. `zna merge` removes the overlap from R2's **3'** "
        "end, leaving R1 and R2 tiling the fragment exactly once.",
        "",
        "Scored over **every** kept pair, split by regime, because the failure modes are "
        "opposite and an average hides both.",
        "",
        "| metric | " + " | ".join(s.name for s in scores) + " |",
        "|---|" + "---:|" * len(scores),
        pair("kept_pairs_no_true_overlap"),
        pair("kept_pairs_overlapping"),
        pair("kept_pairs_read_through"),
        pair("min_trimmable_overlap"),
        pair("trim_exact"),
        pair("trim_over"),
        pair("trim_under"),
        pair("trim_absent"),
        pair("trim_absent_below_threshold"),
        pair("trim_absent_above_threshold"),
        pair("false_trims"),
        pair("false_trim_bases"),
        pair("emitted_read_longer_than_input"),
        pair("trimmed_pairs"),
        pair("trimmed_mean_length_gap"),
        pair("trimmed_max_length_gap"),
        "",
        "`false_trims` is the specificity side and the one an analysis restricted to "
        "overlapping pairs cannot see: a pair whose mates share nothing, trimmed anyway, "
        "has had real sequence deleted. `trim_absent_below_threshold` is the "
        "specification, not a miss — an overlap shorter than `min_trimmable_overlap` "
        "cannot clear `--threshold-trim` and is left in deliberately.",
        "",
        "`trimmed_*_length_gap` is `|len(R1) − len(R2)|` over trimmed pairs, and is why "
        "the overlap is **split** between the two 3' ends rather than taken entirely off "
        "R2: the emitted mates come out the same length, which is what downstream "
        "aligners and models expect. Taking it all off R2 leaves a gap equal to the "
        "whole overlap on every trimmed pair.",
        "",
        "### Bases, which is what a model actually sees",
        "",
        "| metric | " + " | ".join(s.name for s in scores) + " |",
        "|---|" + "---:|" * len(scores),
        pair("bases_emitted_unmerged"),
        pair("fragment_bases_duplicated_if_untrimmed"),
        pair("fragment_bases_removed_correctly"),
        pair("fragment_bases_duplicated"),
        pair("duplicated_per_10k_unmerged_bases"),
        pair("fragment_bases_deleted"),
        pair("trim_benefit_ratio"),
        pair("adapter_bases_left_in_kept_read_through"),
        "",
        "`fragment_bases_duplicated` is the count of fragment positions covered twice in "
        "an emitted pair — the thing the trim exists to remove. `..._if_untrimmed` is "
        "the counterfactual, which is what makes the trade legible: the trim is worth "
        "its cost insofar as `fragment_bases_removed_correctly` exceeds "
        "`fragment_bases_deleted`, and `trim_benefit_ratio` is that quotient. `fastp` "
        "has no redundant-overlap trim at all (it trims *adapter*, not shared "
        "sequence), so its column is the size of the problem rather than a competing "
        "result.",
        "",
        trim_curve_table(truth, scores, readlen),
        "",
        "## 8. Throughput",
        "",
    ]
    if scores[0].elapsed:
        lines += [
            "| tool | wall s | µs/pair |",
            "|---|---:|---:|",
        ] + [f"| {s.name} | {s.elapsed:.2f} | {1e6 * s.elapsed / n:.3f} |"
             for s in scores] + [
            "",
            "Same input, same thread count. fastp also writes three output files and "
            "computes its own statistics, so this is an end-to-end comparison of the "
            "two commands, not of their merge kernels.",
        ]
    lines += ["", "## Where zna merge was wrong", "",
              "`zna_errors.tsv`, by category:", "",
              "| category | count |", "|---|---:|"]
    for c in sorted(scores[0].cat_counts):
        lines.append(f"| {c} | {scores[0].cat_counts[c]:,} |")
    lines += [
        "",
        "`chimera` — merged a pair with no true overlap. `wrong_length` — merged at the "
        "wrong shift. The `_dropped` variants are the same two, for merges that fell "
        "below `--min-read-length` and were never written. `frame_violation` — the "
        "record does not align to the fragment at offset 0. `consensus_miss` — right "
        "length, but more base errors than the oracle floor, i.e. it had the evidence "
        "to fix a base and did not. `false_trim` — cut sequence off a pair that shares "
        "none. `over_trim` / `missed_trim` — cut more or less than the true overlap "
        "(misses are only counted above the trim threshold). `r1_shortened` — moved a "
        "3' end on R1, which nothing should. On trim rows `emitted_len` is R2's emitted "
        "length and `len_err` is bases removed minus the true overlap.",
        "",
        "The `scan_*` columns are `zna merge`'s own kernel re-run on the pair, so a row "
        "carries the evidence the decision was made on rather than only its outcome.",
        "",
    ]
    (outdir / "report.md").write_text("\n".join(lines))


if __name__ == "__main__":
    raise SystemExit(main())
