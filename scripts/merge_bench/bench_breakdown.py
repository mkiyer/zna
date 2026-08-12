"""Where does the per-pair time actually go? Cumulative stages, min-of-N."""
import gzip, io, os, sys, time

from zna.merge import cli
from zna.merge.fastqio import read_pairs, FastqWriter
from zna.merge.pairs import MergeParams, PairOutcome, base_name, process_pair
from zna.merge.overlap import HAVE_NUMBA, find_overlap, reverse_complement

IN1, IN2 = sys.argv[1], sys.argv[2]
REPS = int(os.environ.get("REPS", "3"))
P = MergeParams()

# Warm the JIT before timing anything.
process_pair(b"w", b"ACGT" * 40, b"I" * 160, b"w", b"ACGT" * 40, b"I" * 160, P)

def load_pairs():
    """Parse once into memory so the compute stages exclude I/O entirely."""
    return [(r1[0], r1[1], r1[2], r2[0], r2[1], r2[2])
            for r1, r2 in read_pairs(IN1, IN2, 1)]

t0 = time.perf_counter(); PAIRS = load_pairs(); t_parse_once = time.perf_counter() - t0
N = len(PAIRS)

def bench(name, fn, reps=REPS):
    best = min((lambda: (lambda t: t)( (lambda s: (fn(), time.perf_counter()-s)[1])(time.perf_counter()) ))()
               for _ in range(reps))
    print(f"{name:<46} {best/N*1e6:8.3f} us/pair   ({best:6.3f} s)")
    return best

print(f"# {N} pairs | numba={HAVE_NUMBA} | python={sys.version.split()[0]}")
print(f"# reps={REPS}, min-of-N\n")

# ---- I/O stages -----------------------------------------------------------
def s_gzip_only():
    for path in (IN1, IN2):
        with gzip.open(path, "rb") as fh:
            while fh.read(1 << 20): pass
b_gz = bench("1. gzip inflate only (stdlib, both files)", s_gzip_only)

def s_pigz_only():
    import subprocess, shutil
    pigz = shutil.which("pigz")
    for path in (IN1, IN2):
        p = subprocess.Popen([pigz, "-dc", "-p", "1", path], stdout=subprocess.PIPE, bufsize=1<<20)
        while p.stdout.read(1 << 20): pass
        p.wait()
b_pigz = bench("2. pigz -dc -p1 inflate only (both files)", s_pigz_only)

def s_readlines():
    for r in read_pairs(IN1, IN2, 1): pass
b_read = bench("3.   + FASTQ parse (read_pairs, drain)", s_readlines)

# ---- compute stages (in-memory, no I/O) -----------------------------------
def s_revcomp():
    for _h1, _s1, _q1, _h2, s2, _q2 in PAIRS:
        reverse_complement(s2)
b_rc = bench("4. revcomp(R2) only, in memory", s_revcomp)

def s_scan():
    for _h1, s1, _q1, _h2, s2, _q2 in PAIRS:
        find_overlap(s1, reverse_complement(s2), 8.0, 0.01)
b_scan = bench("5.   + the overlap scan (find_overlap)", s_scan)

def s_pp():
    for h1, s1, q1, h2, s2, q2 in PAIRS:
        process_pair(h1, s1, q1, h2, s2, q2, P)
b_pp = bench("6. process_pair (scan+consensus+build+name)", s_pp)

def s_syncchk():
    for h1, _s1, _q1, h2, _s2, _q2 in PAIRS:
        if base_name(h1) != base_name(h2): raise SystemExit("desync")
b_sync = bench("7. sync check only (base_name x2)", s_syncchk)

def s_chunk():
    cli._init_worker(P, True)
    for i in range(0, N, 2000):
        cli._process_chunk((i, PAIRS[i:i+2000]))
b_chunk = bench("8. _process_chunk (6+7+format+histograms)", s_chunk)

# ---- whole tool -----------------------------------------------------------
import tempfile
def make_args(procs, out, level=1, threads=4):
    return cli.build_parser().parse_args(
        ["--in1", IN1, "--in2", IN2, "--out", out, "--processes", str(procs),
         "--threads", str(threads), "--compress-level", str(level), "-q"])

with tempfile.TemporaryDirectory() as td:
    for procs in (1, 2, 4, 8):
        def run(procs=procs, td=td):
            cli.run(make_args(procs, os.path.join(td, f"o{procs}.fq.gz")))
        bench(f"9. FULL cli.run --processes {procs}", run, reps=2)
    def run_plain(td=td):
        cli.run(make_args(1, os.path.join(td, "plain.fq")))
    bench("10. FULL cli.run -p1, UNCOMPRESSED output", run_plain, reps=2)

print(f"\n# one-off parse into memory: {t_parse_once:.2f} s ({t_parse_once/N*1e6:.3f} us/pair)")
