"""Simulate paired-end reads from a genome with exact ground truth.

Draws fragments of **uniform** length from an indexed FASTA, sequences both ends at a
fixed cycle count, runs into adapter when the fragment is shorter than the read, injects
substitution errors, and writes ``<prefix>_1.fq.gz`` / ``<prefix>_2.fq.gz`` beside a
sidecar TSV that records, for every pair, what the answer is.

    python simulate.py --genome hg38.fa --n-pairs 1000000 --read-length 150 \
                       --frag-min 60 --frag-max 450 --error-rate 0.002 \
                       --quality-model novaseq --out-prefix sim --seed 42

**Uniform fragment lengths are deliberate and are not a library.** The point is to
populate every geometric regime at equal density, so the sensitivity curve is measured
where it bends rather than where a real insert distribution happens to put its mass:

| fragment length      | regime                                                  |
|----------------------|---------------------------------------------------------|
| `> 2 * readlen`      | no overlap at all — any merge is a **chimera**           |
| `= 2 * readlen`      | overlap 0, the boundary case                             |
| `readlen < L < 2*rl` | partial overlap, 1..readlen bases — the sensitivity curve|
| `= readlen`          | full overlap, the mates are exact reverse complements    |
| `< readlen`          | **read-through**: both mates run into adapter            |

Do not quote a merge rate from this simulation as if it were a production number; quote
per-bin sensitivity instead.

**The sidecar is the authority, not the FASTQ comment.** Both tools rewrite headers —
`zna merge` renames a merged record `<id> merged_<n1>_<n2>` and fastp has its own
opinions — so the comment is for eyeballing only. The read ID up to the first whitespace
is the join key: both tools preserve it, and `zna merge` strips only the `/1`,`/2`
suffix, which is what makes it stable.

**Why the error model is written fresh rather than borrowed.** khorana's
`khorana/sim/reads.py` has the shape of this (per-read shuffled SNV/insertion/deletion
passes, so errors compound) but it simulates transcripts and **never runs into adapter**,
which is precisely the regime this benchmark exists to measure. Two deliberate
differences from it: a substitution here always changes the base, where khorana redraws
uniformly from ACGT and is therefore a no-op a quarter of the time; and the error
probability is derived from the quality string rather than being independent of it (see
`--quality-model`).

Nothing here is part of the zna package or its test suite.
"""
from __future__ import annotations

import argparse
import json
import mmap
import sys
import time
from pathlib import Path

import numpy as np

# The two mates' post-fragment sequence, i.e. Illumina TruSeq. These are NOT reverse
# complements of each other, which is what keeps the adapter tails from aligning and
# makes a read-through pair's only true alignment the fragment itself.
ADAPTER1 = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER2 = b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

_COMPLEMENT = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")
_BASES = np.frombuffer(b"ACGT", dtype=np.uint8)

#: A->0 C->1 G->2 T->3; everything else 0. Only reached for ACGT (fragments containing
#: N are rejected, and adapters and padding are drawn from ACGT).
_CODE = np.zeros(256, dtype=np.uint8)
_CODE[np.frombuffer(b"ACGT", dtype=np.uint8)] = np.arange(4, dtype=np.uint8)


def rc(seq: bytes) -> bytes:
    return seq.translate(_COMPLEMENT)[::-1]


# --------------------------------------------------------------------------- #
# the genome
# --------------------------------------------------------------------------- #

class IndexedFasta:
    """Random access to an indexed FASTA, straight off the `.fai` and an mmap.

    pysam would do, but it is not in the benchmark environment and the index format is
    five columns. Sequence comes back UPPERCASE: hg38 from UCSC is soft-masked, and
    `zna merge` upper-cases its input, so leaving the case in would make the truth
    disagree with every tool over repeats for no reason.
    """

    def __init__(self, path):
        self.path = str(path)
        fai = Path(self.path + ".fai")
        if not fai.exists():
            raise SystemExit(f"no FASTA index beside {self.path} (run `samtools faidx`)")
        self.index = {}
        for line in fai.read_text().splitlines():
            if not line.strip():
                continue
            name, length, offset, linebases, linewidth = line.split("\t")[:5]
            self.index[name] = (int(length), int(offset), int(linebases), int(linewidth))
        self._fh = open(self.path, "rb")
        self._mm = mmap.mmap(self._fh.fileno(), 0, access=mmap.ACCESS_READ)

    def fetch(self, name: str, start: int, length: int) -> bytes:
        ln, off, lb, lw = self.index[name]
        end = start + length
        if start < 0 or end > ln:
            raise ValueError(f"{name}:{start}-{end} out of bounds ({ln})")
        b0 = off + (start // lb) * lw + (start % lb)
        b1 = off + ((end - 1) // lb) * lw + ((end - 1) % lb) + 1
        return self._mm[b0:b1].replace(b"\n", b"").replace(b"\r", b"").upper()

    def close(self):
        self._mm.close()
        self._fh.close()


def select_chroms(fa: IndexedFasta, pattern: str | None):
    """Chromosomes to draw from, longest first.

    Default is the primary assembly: `chr1..22,X,Y,M` and their unprefixed spellings.
    Alt/random/Un contigs are excluded because they are duplicate and unusually
    repeat-dense sequence, which would quietly weight the repeat tail of the results.
    """
    if pattern:
        wanted = [c.strip() for c in pattern.split(",") if c.strip()]
        missing = [c for c in wanted if c not in fa.index]
        if missing:
            raise SystemExit(f"not in the FASTA index: {', '.join(missing)}")
        return wanted
    primary = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
    primary |= {str(i) for i in range(1, 23)} | {"X", "Y", "MT", "M"}
    out = [c for c in fa.index if c in primary]
    if not out:
        raise SystemExit(
            f"no primary chromosomes found in {fa.path}.fai; pass --chroms explicitly")
    return out


# --------------------------------------------------------------------------- #
# qualities, and the errors they imply
# --------------------------------------------------------------------------- #
#
# The consensus in `zna merge` is quality-driven: when the mates disagree it keeps the
# better-supported base. With a constant quality string every disagreement is a tie and
# R1 always wins, so the whole consensus is untested and metric 4 of the benchmark plan
# measures nothing. Hence two modes, and the accuracy comparison belongs on `novaseq`.
#
# In `novaseq` mode the quality is drawn FIRST, from a four-bin position-dependent
# profile, and the error indicator is then drawn from it at 10^(-Q/10). That ordering is
# the point: it makes Q monotone in the real error probability, which is the only
# property the posterior needs. It is not a claim to reproduce any particular
# instrument's binning.

_NOVASEQ_Q = np.array([37, 23, 12, 2], dtype=np.int16)


def novaseq_profile(readlen: int) -> np.ndarray:
    """(readlen, 4) bin probabilities, degrading with cycle."""
    t = np.linspace(0.0, 1.0, readlen, dtype=np.float64)[:, None]
    w = np.stack([
        0.960 - 0.140 * t[:, 0],     # Q37
        0.030 + 0.090 * t[:, 0],     # Q23
        0.009 + 0.050 * t[:, 0],     # Q12
        np.full(readlen, 0.001),     # Q2, the no-call bin
    ], axis=1)
    return w / w.sum(axis=1, keepdims=True)


class QualityModel:
    """Draws (quality string, error mask) for a block of reads."""

    def __init__(self, kind: str, readlen: int, error_rate: float):
        self.kind = kind
        self.readlen = readlen
        self.error_rate = error_rate
        if kind == "flat":
            self.scale = 1.0
            self._qchar = np.full(readlen, ord("I"), dtype=np.uint8)   # Q40
            self._p = np.full(readlen, error_rate, dtype=np.float64)
            return
        self._prob = novaseq_profile(readlen)
        self._cum = np.cumsum(self._prob, axis=1)
        perr = 10.0 ** (-_NOVASEQ_Q / 10.0)
        natural = float((self._prob * perr).sum() / readlen)
        # Keep --error-rate meaningful in both modes by scaling the implied per-base
        # probabilities to it. The scale is reported: at 1.0 the quality strings are
        # calibrated, and far from 1.0 they are miscalibrated by exactly that factor
        # (still monotone, which is what the consensus uses).
        self.scale = error_rate / natural if natural > 0 else 1.0
        self._perr = perr

    def draw(self, rng, n: int):
        """Returns ``(qual uint8 (n, readlen), error mask bool (n, readlen))``."""
        if self.kind == "flat":
            qual = np.broadcast_to(self._qchar, (n, self.readlen))
            mask = rng.random((n, self.readlen)) < self._p
            return qual, mask
        u = rng.random((n, self.readlen))
        bins = np.zeros((n, self.readlen), dtype=np.uint8)
        for k in range(3):
            bins += (u >= self._cum[:, k]).astype(np.uint8)
        qual = (_NOVASEQ_Q[bins] + 33).astype(np.uint8)
        p = self._perr[bins] * self.scale
        mask = rng.random((n, self.readlen)) < p
        return qual, mask


# --------------------------------------------------------------------------- #
# the simulation
# --------------------------------------------------------------------------- #

def draw_fragments(rng, fa, chroms, cum, lengths, frag_len):
    """Fragment sequences for one block, N-free, uppercase, plus their coordinates.

    Starts are uniform over the concatenated selected sequence, and any window
    containing an `N` is redrawn — which is uniform over the N-free windows, and keeps
    the tools' N policies out of the measurement. hg38's N runs are megabase-scale, so
    the rejection rate is a few percent and the redraw loop terminates immediately.
    """
    n = len(frag_len)
    chrom_i = np.empty(n, dtype=np.int32)
    start = np.empty(n, dtype=np.int64)
    frags = [None] * n
    todo = np.arange(n)
    rounds = 0
    while todo.size:
        rounds += 1
        if rounds > 200:
            raise SystemExit("could not draw N-free fragments; check --chroms")
        g = rng.integers(0, cum[-1], size=todo.size)
        ci = np.searchsorted(cum, g, side="right")
        st = g - np.where(ci > 0, cum[ci - 1], 0)
        still = []
        for k, idx in enumerate(todo):
            c, s = int(ci[k]), int(st[k])
            flen = int(frag_len[idx])
            # A fragment running off the end of a chromosome is redrawn rather than
            # slid back inside it: clamping would pile every overhanging draw onto the
            # same start, which is a small artifact but a real one.
            if s + flen > lengths[c]:
                still.append(idx)
                continue
            seq = fa.fetch(chroms[c], s, flen)
            if b"N" in seq:
                still.append(idx)
                continue
            chrom_i[idx] = c
            start[idx] = s
            frags[idx] = seq
        todo = np.array(still, dtype=np.int64)
    return chrom_i, start, frags


def build_reads(rng, frags, strand, readlen, adapter1, adapter2):
    """The two raw mates for each fragment, at full cycle length.

    ``r1 = (frag + adapter1)[:readlen]`` and ``r2 = (rc(frag) + adapter2)[:readlen]``,
    padded with random bases when the fragment plus its adapter is still short of the
    cycle count. Real chemistry would continue into the index and P7, which is constant
    across reads; random padding is used instead so that nothing downstream of the
    adapter can align by construction.
    """
    n = len(frags)
    out1 = [None] * n
    out2 = [None] * n
    pad_needed = 0
    for i in range(n):
        f = frags[i]
        if strand[i]:
            f = rc(f)
            frags[i] = f
        a = f + adapter1
        b = rc(f) + adapter2
        if len(a) < readlen:
            pad_needed += 1
        out1[i] = a[:readlen]
        out2[i] = b[:readlen]
    if pad_needed:
        # One draw for every read that needs filler, so the padding is deterministic in
        # the seed and never shared between the two mates.
        for i in range(n):
            if len(out1[i]) < readlen:
                out1[i] += rand_bases(rng, readlen - len(out1[i]))
            if len(out2[i]) < readlen:
                out2[i] += rand_bases(rng, readlen - len(out2[i]))
    return out1, out2


def rand_bases(rng, n: int) -> bytes:
    return _BASES[rng.integers(0, 4, size=n)].tobytes()


def inject_errors(rng, reads, mask, readlen):
    """Apply the error mask to a block of equal-length reads, in place on a copy.

    Returns ``(array (n, readlen) uint8, n_err per read, positions per read)``. A
    substitution always CHANGES the base — drawing uniformly from ACGT, as khorana's
    simulator does, would silently make a quarter of the "errors" no-ops and inflate
    every error count reported here.
    """
    n = len(reads)
    arr = np.frombuffer(b"".join(reads), dtype=np.uint8).reshape(n, readlen).copy()
    rows, cols = np.nonzero(mask)
    new = np.empty(0, dtype=np.uint8)
    if rows.size:
        cur = _CODE[arr[rows, cols]]
        bump = rng.integers(1, 4, size=rows.size, dtype=np.uint8)
        new = _BASES[(cur + bump) & 3]
        arr[rows, cols] = new
    n_err = np.bincount(rows, minlength=n)
    bounds = np.searchsorted(rows, np.arange(n + 1))
    return arr, n_err, cols, new, bounds


def self_check(prefix, fa, readlen, adapter1, adapter2, limit):
    """Re-read the output and prove the sidecar tells the truth.

    Everything the benchmark concludes rests on this file being right, so it is checked
    against the genome rather than against the code that wrote it: for each pair, that
    the recorded fragment is really at `chrom:start` on `strand`, that the geometry
    columns follow from the fragment length, and that both reads rebuild **exactly** from
    fragment + adapter + the recorded substitutions. The random filler past the adapter
    is the one thing not reconstructible, and it is never fragment sequence.
    """
    import gzip
    rows = []
    with open(f"{prefix}.truth.tsv") as fh:
        fh.readline()
        for line in fh:
            rows.append(line.rstrip("\n").split("\t"))
            if len(rows) >= limit:
                break
    checked = 0
    for mate, path, adapter in ((1, f"{prefix}_1.fq.gz", adapter1),
                                (2, f"{prefix}_2.fq.gz", adapter2)):
        with gzip.open(path, "rb") as fh:
            for row in rows:
                h = fh.readline().split()[0]
                seq = fh.readline().strip()
                fh.readline()
                qual = fh.readline().strip()
                rid, chrom, start, strand, flen = row[0], row[1], row[2], row[3], int(row[4])
                frag = row[12].encode()
                assert h == b"@%s/%d" % (rid.encode(), mate), (h, rid)
                assert len(seq) == len(qual) == readlen, (rid, len(seq), len(qual))
                if mate == 1:
                    g = fa.fetch(chrom, int(start), flen)
                    assert frag == (rc(g) if strand == "-" else g), f"{rid}: not at {chrom}:{start}"
                    assert len(frag) == flen, rid
                    rt = 1 if flen < readlen else 0
                    ovl = flen if rt else max(0, min(2 * readlen - flen, readlen))
                    assert int(row[5]) == ovl and int(row[6]) == rt, f"{rid}: geometry"
                base = (frag if mate == 1 else rc(frag)) + adapter
                want = bytearray(base[:readlen])
                if row[11] != ".":
                    for tok in row[11].split(","):
                        if int(tok[0]) != mate:
                            continue
                        p = int(tok[2:-1])
                        assert p < readlen, f"{rid}: error past the read"
                        if p < len(want):
                            # A substitution must actually change the base, or the
                            # recorded error counts overstate what was injected.
                            assert ord(tok[-1]) != base[p], f"{rid}: no-op substitution"
                            want[p] = ord(tok[-1])
                assert seq[:len(want)] == bytes(want), f"{rid}/{mate}: read does not rebuild"
                checked += 1
    print(f"self-check: {checked} reads rebuilt exactly from the genome and the sidecar",
          file=sys.stderr)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", required=True, help="indexed FASTA (.fai beside it)")
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--n-pairs", type=int, default=1_000_000)
    ap.add_argument("--read-length", type=int, default=150)
    ap.add_argument("--frag-min", type=int, default=60, help="inclusive")
    ap.add_argument("--frag-max", type=int, default=450, help="inclusive")
    ap.add_argument("--error-rate", type=float, default=0.002,
                    help="mean per-base substitution rate")
    ap.add_argument("--quality-model", choices=("flat", "novaseq"), default="novaseq",
                    help="flat makes every mismatch a tie, which leaves the quality-aware "
                         "consensus untested; novaseq draws Q first and the error from it")
    ap.add_argument("--adapter1", default=ADAPTER1.decode())
    ap.add_argument("--adapter2", default=ADAPTER2.decode())
    ap.add_argument("--chroms", default=None,
                    help="comma-separated; default is the primary assembly")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--block", type=int, default=50_000, help="fragments per block")
    ap.add_argument("--io-threads", type=int, default=4, help="pigz threads per output")
    ap.add_argument("--self-check", type=int, default=2000, metavar="N",
                    help="re-read N pairs and verify them against the genome and the "
                         "sidecar; 0 to skip")
    args = ap.parse_args(argv)

    if args.frag_min < 1 or args.frag_max < args.frag_min:
        raise SystemExit("need 1 <= --frag-min <= --frag-max")
    if args.read_length < 1:
        raise SystemExit("--read-length must be >= 1")
    if not 0.0 <= args.error_rate < 1.0:
        raise SystemExit("--error-rate must be in [0, 1)")

    try:
        from zna.merge.fastqio import FastqWriter
    except ImportError:                                     # pragma: no cover
        raise SystemExit("this script needs zna importable (it reuses its pigz writer)")

    readlen = args.read_length
    adapter1 = args.adapter1.encode().upper()
    adapter2 = args.adapter2.encode().upper()
    rng = np.random.default_rng(args.seed)
    qm = QualityModel(args.quality_model, readlen, args.error_rate)
    if not 0.5 <= qm.scale <= 2.0:
        print(f"note: quality bins imply a mean error rate {1 / qm.scale:.2f}x the "
              f"requested {args.error_rate}; they are scaled to match, so Q is "
              f"miscalibrated by that factor (still monotone)", file=sys.stderr)

    fa = IndexedFasta(args.genome)
    chroms = select_chroms(fa, args.chroms)
    lengths = np.array([fa.index[c][0] for c in chroms], dtype=np.int64)
    if lengths.min() <= args.frag_max:
        chroms = [c for c, ln in zip(chroms, lengths) if ln > args.frag_max]
        lengths = np.array([fa.index[c][0] for c in chroms], dtype=np.int64)
    cum = np.cumsum(lengths)

    prefix = Path(args.out_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    width = max(7, len(str(max(args.n_pairs - 1, 0))))
    t0 = time.perf_counter()

    regimes = {"no_overlap": 0, "boundary": 0, "partial": 0, "full": 0, "read_through": 0}
    tot_err = 0
    tot_bases = 0

    with FastqWriter(f"{prefix}_1.fq.gz", threads=args.io_threads, level=1) as w1, \
         FastqWriter(f"{prefix}_2.fq.gz", threads=args.io_threads, level=1) as w2, \
         open(f"{prefix}.truth.tsv", "w", buffering=1 << 20) as truth:
        truth.write("read_id\tchrom\tstart\tstrand\tfrag_len\ttrue_ovl\tread_through\t"
                    "len1\tlen2\tn_err1\tn_err2\terr_positions\tfragment\n")
        done = 0
        while done < args.n_pairs:
            n = min(args.block, args.n_pairs - done)
            frag_len = rng.integers(args.frag_min, args.frag_max + 1, size=n)
            chrom_i, start, frags = draw_fragments(
                rng, fa, chroms, cum, lengths, frag_len)
            strand = rng.random(n) < 0.5
            raw1, raw2 = build_reads(rng, frags, strand, readlen, adapter1, adapter2)

            q1, m1 = qm.draw(rng, n)
            q2, m2 = qm.draw(rng, n)
            a1, ne1, cols1, sub1, bnd1 = inject_errors(rng, raw1, m1, readlen)
            a2, ne2, cols2, sub2, bnd2 = inject_errors(rng, raw2, m2, readlen)
            tot_err += int(ne1.sum()) + int(ne2.sum())
            tot_bases += 2 * n * readlen

            s1 = a1.tobytes()
            s2 = a2.tobytes()
            qs1 = np.ascontiguousarray(q1).tobytes()
            qs2 = np.ascontiguousarray(q2).tobytes()

            b1, b2, rows = [], [], []
            for i in range(n):
                L = int(frag_len[i])
                rt = 1 if L < readlen else 0
                ovl = L if rt else max(0, min(2 * readlen - L, readlen))
                if rt:
                    regimes["read_through"] += 1
                elif L == readlen:
                    regimes["full"] += 1
                elif L < 2 * readlen:
                    regimes["partial"] += 1
                elif L == 2 * readlen:
                    regimes["boundary"] += 1
                else:
                    regimes["no_overlap"] += 1

                rid = f"sim{done + i:0{width}d}"
                o = i * readlen
                e1, e2 = int(ne1[i]), int(ne2[i])
                cmt = f" fl={L} ov={ovl} rt={rt} e1={e1} e2={e2}"
                b1.append(b"@%s/1%s\n%b\n+\n%b\n" % (
                    rid.encode(), cmt.encode(), s1[o:o + readlen], qs1[o:o + readlen]))
                b2.append(b"@%s/2%s\n%b\n+\n%b\n" % (
                    rid.encode(), cmt.encode(), s2[o:o + readlen], qs2[o:o + readlen]))

                # `<mate>.<pos><observed base>`, so the sidecar reconstructs the reads
                # exactly (everything but the random filler past the adapter) and the
                # scorer never has to open the FASTQs to know what the tools were shown.
                pos = []
                if e1:
                    pos += [f"1.{p}{chr(b)}" for p, b in
                            zip(cols1[bnd1[i]:bnd1[i + 1]], sub1[bnd1[i]:bnd1[i + 1]])]
                if e2:
                    pos += [f"2.{p}{chr(b)}" for p, b in
                            zip(cols2[bnd2[i]:bnd2[i + 1]], sub2[bnd2[i]:bnd2[i + 1]])]
                rows.append(
                    f"{rid}\t{chroms[chrom_i[i]]}\t{start[i]}\t"
                    f"{'-' if strand[i] else '+'}\t{L}\t{ovl}\t{rt}\t"
                    f"{readlen}\t{readlen}\t{e1}\t{e2}\t"
                    f"{','.join(pos) if pos else '.'}\t{frags[i].decode()}\n")

            w1.write_raw(b"".join(b1))
            w2.write_raw(b"".join(b2))
            truth.write("".join(rows))
            done += n
            print(f"\r  {done:,} / {args.n_pairs:,} pairs", end="", file=sys.stderr)
    print("", file=sys.stderr)

    elapsed = time.perf_counter() - t0
    meta = {
        "tool": "merge_bench/simulate.py",
        "genome": str(Path(args.genome).resolve()),
        "chroms": chroms,
        "n_pairs": args.n_pairs,
        "read_length": readlen,
        "frag_min": args.frag_min,
        "frag_max": args.frag_max,
        "frag_length_distribution": "uniform, inclusive",
        "error_rate_requested": args.error_rate,
        "error_rate_realised": round(tot_err / tot_bases, 6) if tot_bases else 0.0,
        "quality_model": args.quality_model,
        # 1.0 means the quality strings are calibrated: P(error) == 10^(-Q/10).
        "quality_calibration_scale": round(qm.scale, 4),
        "adapter1": adapter1.decode(),
        "adapter2": adapter2.decode(),
        "seed": args.seed,
        "regimes": regimes,
        "elapsed_s": round(elapsed, 1),
        "outputs": {
            "r1": f"{prefix}_1.fq.gz",
            "r2": f"{prefix}_2.fq.gz",
            "truth": f"{prefix}.truth.tsv",
        },
    }
    Path(f"{prefix}.json").write_text(json.dumps(meta, indent=2) + "\n")
    if args.self_check > 0:
        self_check(prefix, fa, readlen, adapter1, adapter2,
                   min(args.self_check, args.n_pairs))
    fa.close()
    print(f"wrote {args.n_pairs:,} pairs to {prefix}_[12].fq.gz in {elapsed:.1f}s "
          f"({', '.join(f'{k}={v:,}' for k, v in regimes.items())})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
