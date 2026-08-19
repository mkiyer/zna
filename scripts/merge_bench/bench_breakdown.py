"""Where does the per-pair time actually go? Cumulative stages, min-of-N.

    python bench_breakdown.py r1.fq.gz r2.fq.gz          # REPS=3 by default

The stages are cumulative from the bottom: each adds one layer to the one above, so a
difference between adjacent rows is that layer's cost. Stages 4-7 run entirely from an
in-memory parse, so they exclude I/O; stages 1-3 are the I/O path on its own; stage 8
is one backend chunk call, which is where the shipped tool actually spends its time
(all per-pair work happens inside one GIL-releasing call per chunk -- see cli.py); and
stage 9 is the whole CLI at several thread counts.
"""
import gzip
import os
import sys
import tempfile
import time

from zna.merge import backend as _backend
from zna.merge import cli
from zna.merge.fastqio import read_pairs
from zna.merge.overlap import backend_name, find_overlap, reverse_complement
from zna.merge.pairs import MergeParams, base_name, process_pair
from zna.merge.params import DISAGREE_Q

IN1, IN2 = sys.argv[1], sys.argv[2]
REPS = int(os.environ.get("REPS", "3"))
P = MergeParams()


def load_pairs():
    """Parse once into memory so the compute stages exclude I/O entirely."""
    return [(r1[0], r1[1], r1[2], r2[0], r2[1], r2[2])
            for r1, r2 in read_pairs(IN1, IN2, 1)]


t0 = time.perf_counter()
PAIRS = load_pairs()
t_parse_once = time.perf_counter() - t0
N = len(PAIRS)


def bench(name, fn, reps=REPS):
    best = float("inf")
    for _ in range(reps):
        t = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t)
    print(f"{name:<46} {best / N * 1e6:8.3f} us/pair   ({best:6.3f} s)")
    return best


print(f"# {N} pairs | backend={backend_name()} | python={sys.version.split()[0]}")
print(f"# reps={REPS}, min-of-N\n")

# ---- I/O stages -----------------------------------------------------------
def s_gzip_only():
    for path in (IN1, IN2):
        with gzip.open(path, "rb") as fh:
            while fh.read(1 << 20):
                pass


bench("1. gzip inflate only (stdlib, both files)", s_gzip_only)


def s_pigz_only():
    import shutil
    import subprocess
    pigz = shutil.which("pigz")
    if pigz is None:
        print("2. pigz -dc -p1 inflate only                     (pigz not on PATH)")
        return
    for path in (IN1, IN2):
        p = subprocess.Popen([pigz, "-dc", "-p", "1", path],
                             stdout=subprocess.PIPE, bufsize=1 << 20)
        while p.stdout.read(1 << 20):
            pass
        p.wait()


bench("2. pigz -dc -p1 inflate only (both files)", s_pigz_only)


def s_readlines():
    for _ in read_pairs(IN1, IN2, 1):
        pass


bench("3.   + FASTQ parse (read_pairs, drain)", s_readlines)

# ---- compute stages (in-memory, no I/O) -----------------------------------
def s_revcomp():
    for _h1, _s1, _q1, _h2, s2, _q2 in PAIRS:
        reverse_complement(s2)


bench("4. revcomp(R2) only, in memory", s_revcomp)


def s_scan():
    for _h1, s1, _q1, _h2, s2, _q2 in PAIRS:
        find_overlap(s1, reverse_complement(s2), P)


bench("5.   + the overlap scan (find_overlap)", s_scan)


def s_pp():
    for h1, s1, q1, h2, s2, q2 in PAIRS:
        process_pair(h1, s1, q1, h2, s2, q2, P)


bench("6. process_pair (scan+consensus+build+name)", s_pp)


def s_syncchk():
    for h1, _s1, _q1, h2, _s2, _q2 in PAIRS:
        if base_name(h1) != base_name(h2):
            raise SystemExit("desync")


bench("7. sync check only (base_name x2)", s_syncchk)

# ---- one backend chunk call, which is what the tool actually runs ---------
#
# `merge_chunk` takes the two RAW FASTQ buffers and does everything -- parse, scan,
# consensus, construction, formatting, histograms -- in one call that releases the GIL.
# Rebuilding those buffers here from the parsed records would measure a different
# program, so read the decompressed bytes once and hand them over as they are.
def _raw_buffers():
    out = []
    for path in (IN1, IN2):
        with gzip.open(path, "rb") as fh:
            out.append(fh.read())
    return out


BUF1, BUF2 = _raw_buffers()
_bk = _backend.active()
_CHUNK_KW = dict(match_q=P.match_q, step_q=P.step_q, t_merge_q=P.t_merge_q,
                 t_trim_q=P.t_trim_q, min_read_length=P.min_read_length,
                 disagree_q=DISAGREE_Q, check_sync=True, base_index=0,
                 npolicy=P.npolicy_code, rng_seed=P.rng_seed)


def s_chunk():
    _bk.merge_chunk(BUF1, 0, len(BUF1), BUF2, 0, len(BUF2), **_CHUNK_KW)


bench("8. merge_chunk, whole file as one chunk", s_chunk)

# ---- whole tool -----------------------------------------------------------
def make_args(out, threads=1, level=1):
    return cli.build_parser().parse_args(
        ["--in1", IN1, "--in2", IN2, "--out", out,
         "--threads", str(threads), "--compress-level", str(level), "-q"])


with tempfile.TemporaryDirectory() as td:
    for th in (1, 2, 4):
        def run(th=th, td=td):
            cli.run(make_args(os.path.join(td, f"o{th}.fq.gz"), threads=th))
        bench(f"9. FULL cli.run --threads {th}", run, reps=2)

    def run_plain(td=td):
        cli.run(make_args(os.path.join(td, "plain.fq")))

    bench("10. FULL cli.run -t1, UNCOMPRESSED output", run_plain, reps=2)

print(f"\n# one-off parse into memory: {t_parse_once:.2f} s "
      f"({t_parse_once / N * 1e6:.3f} us/pair)")
