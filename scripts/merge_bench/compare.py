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

**An unmerged pair is kept whole** (zna 0.6, ``docs/archive/MERGE_ACCURACY_PLAN.md`` §2): there
is no trim band, so a kept pair is scored for being exactly its two input mates -- no
base removed, none added -- and for what it carries into the corpus (the overlap both
mates still hold, and on a kept read-through, adapter). zna's policy flags pass
straight through: ``--alpha``, ``--error-rate``, and ``--adapter-trimmed``, which is a
DECLARATION -- ``simulate.py`` writes raw adapter read-through, so declaring it there is
false and forbids exactly the merges those pairs need; it exists to price that.

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
    raise SystemExit("this script needs zna importable (it re-runs zna's own decision)")

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

    def __init__(self, name, truth, readlen, row_cap, max_edit, min_read_length):
        self.name = name
        self.truth = truth
        self.readlen = readlen
        self.row_cap = row_cap
        self.max_edit = max_edit
        self.min_read_length = min_read_length
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

        *extra* is ``(emitted_len, len_err)`` for rows whose sequence is not held — an
        altered kept mate is described by how its length moved, not by its bases.
        """
        self.cat_counts[category] = self.cat_counts.get(category, 0) + 1
        # Wrong merges are always kept, capped or not: their evidence is what turns a
        # chimera count into an explanation, and there should never be many of them.
        detailed = category in ("chimera", "wrong_length", "chimera_dropped",
                                "wrong_length_dropped", "frame_violation",
                                "kept_mate_altered")
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

        `scan` is `zna merge`'s overlap decision, run on the pair the tool was given at
        the run's own parameters, so `scan_*` says what the shipped rule sees: its
        verdict (merge / implausible / none), the shift, the evidence in bits, and how
        well the two reads actually agree there. A chimera with 90% identity over 89
        bases is the genome repeating, which is a different finding from a merger
        inventing an alignment — and the columns say which.
        """
        t = self.truth
        for i, category, seq, extra in self.pending:
            L = int(t.frag_len[i])
            frag = t.frag[i]
            r1, r2 = reads[i]
            verdict, shift, score_q, olen, mism = scan(r1, r2)
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
                    # The scan's contract is "the best eligible shift at or above the
                    # pair's floor T(len1 + len2 - 1), else nothing", so a miss is only
                    # a defect when the true shift clears that floor. Comparing against
                    # a `no overlap found` zero would otherwise report the floor itself
                    # as a search failure.
                    missed = ((score_q < at[0]) if olen
                              else (at[0] >= scan.floor_q(len(r1), len(r2))))
                    if missed:
                        self.argmax_below_truth += 1
                    elif olen:
                        self.argmax_margin.append((score_q - at[0]) / SCALE)
            if self.cat_written.get(category, 0) >= self.row_cap:
                continue
            self.cat_written[category] = self.cat_written.get(category, 0) + 1
            if seq is None:                    # a kept-mate row, or a filtered-away merge
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
                int(t.n_err2[i]), emitted_len, len_err, mm, ed, verdict,
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

        self._score_kept(both)
        return self

    # -- kept pairs: whole, and what they carry ---------------------------- #

    def _score_kept(self, both):
        """What the *unmerged* pairs are, and what they cost the corpus.

        zna 0.6 keeps an unmerged pair whole -- khorana trains on one randomly chosen
        mate of it, and wants that mate as read -- so the contract is simple: each kept
        mate is its input, exactly as long, and nothing else. Every simulated read is
        free of no-calls, so the N policy never fires here and any change of length is
        a defect (``kept_mate_altered``; the C1 prefix check covers the 5' ends).

        Scored over **every** kept pair, split by regime, because what a kept pair
        carries differs by regime:

        * `true_ovl == 0` -- the mates share nothing; keeping them is the only answer.
        * `0 < true_ovl`, no read-through -- a merge the tool did not make: the overlap
          stays in both mates (``overlap_bases_in_kept_pairs``). Harmless to a model that
          reads one mate, but it is sensitivity forgone -- see sections 1-2.
        * read-through -- both mates carry the WHOLE fragment plus adapter; a kept
          read-through puts adapter into the corpus (``adapter_bases_in_kept_read_through``),
          which is what an undeclared run's read-through merges exist to prevent.
        """
        t = self.truth
        RL = self.readlen
        no_ovl = both & (t.true_ovl == 0)
        normal = both & (t.read_through == 0) & (t.true_ovl > 0)
        rthru = both & (t.read_through == 1)
        self.kept_no_overlap = int(no_ovl.sum())
        self.kept_overlapping = int(normal.sum())
        self.kept_read_through = int(rthru.sum())

        altered = both & ((self.r1_len != t.len1) | (self.r2_len != t.len2))
        self.kept_altered = int(altered.sum())
        self.kept_grew = int((both & ((self.r1_len > t.len1)
                                      | (self.r2_len > t.len2))).sum())
        self.bases_emitted_unmerged = int(self.r1_len[both].sum()
                                          + self.r2_len[both].sum())
        self.bases_overlap_kept = int(np.where(normal, t.true_ovl, 0).sum())
        # Each mate of a read-through pair reads RL - frag_len bases past its molecule.
        self.bases_adapter_kept = int(np.where(rthru, 2 * (RL - t.frag_len), 0).sum())

        for i in np.flatnonzero(altered):
            d = int(self.r1_len[i] + self.r2_len[i] - t.len1[i] - t.len2[i])
            self._note(int(i), "kept_mate_altered", None,
                       (int(self.r1_len[i] + self.r2_len[i]), d))


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


def make_scan(p):
    """`zna merge`'s own overlap decision, at the parameters the run actually used.

    Returns ``scan(r1, r2) -> (verdict, shift, score_q, overlap_len, mismatches)`` -- the
    authoritative decision (``zna.merge.overlap.find_overlap``: the best shift eligible
    under the run's contract, at the pair's own floor, and the plausibility gate's
    verdict on it) -- and ``score_at(r1, r2, shift)``. The second one is what separates
    a *defective* scan from an *ambiguous* input: score the shift the truth says is
    right, and compare. If the tool's pick ever scores lower than the truth's, the
    argmax or its pruning is broken; if it always scores higher, the tool is maximising
    correctly and the sequence really does align better somewhere else.
    """
    from zna.merge.overlap import find_overlap

    def scan(r1, r2):
        o = find_overlap(r1, rc(r2), p)
        return o.verdict, o.shift, o.score_q, o.overlap_len, o.mismatches

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
    scan.floor_q = lambda len1, len2: p.t_q(len1 + len2 - 1)
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
    ap.add_argument("--alpha", default=None,
                    help="passed to zna merge; default is the tool's own")
    ap.add_argument("--error-rate", default=None,
                    help="passed to zna merge; default is the tool's own")
    ap.add_argument("--adapter-trimmed", action="store_true",
                    help="passed to zna merge. FALSE on simulate.py output, which carries "
                         "raw adapter read-through: use it to price the declaration")
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
    if args.alpha is not None:
        zna_cmd += ["--alpha", args.alpha]
    if args.error_rate is not None:
        zna_cmd += ["--error-rate", args.error_rate]
    if args.adapter_trimmed:
        zna_cmd += ["--adapter-trimmed"]
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

    # The run's own policy, rebuilt for the re-scan and checked against the integers
    # the run reports it used -- so the evidence columns are the decision that was made.
    from zna.merge.params import DEFAULT_ALPHA, DEFAULT_ERROR_RATE, MergeParams
    zj = json.loads((outdir / "zna.json").read_text())
    policy = MergeParams(alpha=args.alpha or DEFAULT_ALPHA,
                         error_rate=args.error_rate or DEFAULT_ERROR_RATE,
                         adapter_trimmed=args.adapter_trimmed,
                         min_read_length=args.min_read_length)
    if (zj["params"]["alpha"], zj["params"]["match_q"], zj["params"]["step_q"],
            zj["adapter_trimmed"]) != (float(policy.alpha_exact), policy.match_q,
                                       policy.step_q, policy.adapter_trimmed):
        raise SystemExit(f"{outdir}/zna.json was not made with this policy "
                         f"(--skip-run with different flags?)")

    scores = []
    for name, records in (("zna", zna_records(zna_out)),
                          ("fastp", fastp_records(fp_m, fp_1, fp_2))):
        print(f"scoring {name} ...", file=sys.stderr)
        sc = ToolScore(name, truth, readlen, args.row_cap, args.max_edit,
                       args.min_read_length)
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
        scan = make_scan(policy)
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
        # --- kept pairs: whole (0.6 has no trim band) -------------------- #
        "kept_pairs_no_true_overlap": sc.kept_no_overlap,
        "kept_pairs_overlapping": sc.kept_overlapping,
        "kept_pairs_read_through": sc.kept_read_through,
        # A kept mate is its input, exactly. Must be zero: the simulation has no
        # no-calls, so not even the N policy may change a length.
        "kept_mate_altered": sc.kept_altered,
        "emitted_read_longer_than_input": sc.kept_grew,
        # What the corpus actually receives from kept pairs, in bases.
        "bases_emitted_unmerged": sc.bases_emitted_unmerged,
        "overlap_bases_in_kept_pairs": sc.bases_overlap_kept,
        "adapter_bases_in_kept_read_through": sc.bases_adapter_kept,
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
                "# scan_* is zna's own overlap decision re-run on this pair at the run's\n"
                "# parameters: its verdict, the shift, the evidence in bits, the overlap\n"
                "# length, and how well the two reads agree there. A '.' in emitted_*\n"
                "# means the merged record was below --min-read-length and never written.\n"
                "read_id\tcategory\tchrom\tstart\tstrand\tfrag_len\ttrue_ovl\t"
                "read_through\tn_err1\tn_err2\temitted_len\tlen_err\tn_mismatch\t"
                "edit_distance\tscan_verdict\tscan_shift\tscan_score_bits\tscan_olen\t"
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
        "## 7. Kept pairs, and what the corpus receives from them",
        "",
        "An unmerged pair is **kept whole**: each emitted mate is its input read, "
        "exactly. `kept_mate_altered` and `emitted_read_longer_than_input` must be zero "
        "(the simulation has no no-calls, so not even the N policy may change a length); "
        "the C1 prefix check of section 5 covers the 5' ends.",
        "",
        "| metric | " + " | ".join(s.name for s in scores) + " |",
        "|---|" + "---:|" * len(scores),
        pair("kept_pairs_no_true_overlap"),
        pair("kept_pairs_overlapping"),
        pair("kept_pairs_read_through"),
        pair("kept_mate_altered"),
        pair("emitted_read_longer_than_input"),
        pair("bases_emitted_unmerged"),
        pair("overlap_bases_in_kept_pairs"),
        pair("adapter_bases_in_kept_read_through"),
        "",
        "`overlap_bases_in_kept_pairs` is the true overlap still held by both mates of "
        "an overlapping pair the tool did not merge -- sensitivity forgone, not a "
        "defect of the kept pair. `adapter_bases_in_kept_read_through` is adapter "
        "reaching the corpus through a read-through pair left unmerged: 0.6 merges "
        "those unless `--adapter-trimmed` is declared, and on this raw input the "
        "declaration is false.",
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
        "to fix a base and did not. `kept_mate_altered` — a kept pair whose mates are "
        "not their input lengths; on those rows `emitted_len` is the pair's emitted "
        "bases and `len_err` their change.",
        "",
        "The `scan_*` columns are `zna merge`'s own decision re-run on the pair at the "
        "run's parameters (verdict, shift, bits, overlap), so a row carries the "
        "evidence the decision was made on rather than only its outcome.",
        "",
    ]
    (outdir / "report.md").write_text("\n".join(lines))


if __name__ == "__main__":
    raise SystemExit(main())
