"""Mine a compact truth panel from a ground-truth substrate, with zna 0.5.3's merger.

The tool that made the panel ``panel_eval.py`` evaluates (``docs/archive/MERGE_ACCURACY_PLAN.md``
§8, layer 2): one streaming pass over a simulated substrate keeps EVERY pair zna 0.5.3
gets wrong -- a wrong merge, a lost fragment, a wrong trim, a detected-but-wrong
alignment that was kept, an N-policy drop -- plus a deterministic stratified sample of
the pairs it gets right, each sampled row carrying the weight that projects it back to
the full substrate. Run over 14.1M pairs from 46 substrates it took 22 s on the compiled
backend and produced 245k rows (62 MB).

**It mines with 0.5.3's semantics, and must run under zna 0.5.3.** The ``stratum``
column (``<class>|<overlap bin>``) records what 0.5.3 did with each pair -- including
its trim band, which 0.6 removed -- and that recorded class is the baseline every later
policy is ledgered against. So this script needs the 0.5.3 API (``MergeParams`` with
``t_merge``/``t_trim``, a ``"trimmed"`` outcome, ``find_overlap`` returning a direction)
and refuses to run against a newer tree. Re-mining is only ever needed for a NEW
substrate; reproduce the environment with a pinned install::

    pip install "zna==0.5.3"        # in a scratch env; never the development one
    python mine_panel.py <substrate_dir> <bench_name> <out.tsv.gz>

Substrate format (the policy study's unified format): ``R1.fq[.gz]``, ``R2.fq[.gz]``
and ``truth.tsv.gz`` -- one row per pair in FASTQ order with columns ``pair name L frag
true_r1 true_r2 source`` (``frag`` in R1 orientation; ``source`` a locus id).

Classes (0.5.3):

======== ============================================================================
M+       merged at the true fragment length ``L``, record emitted
M-       merged at a wrong length, record emitted
LOST     merged at a wrong length, record under ``--min-read-length``: fragment gone
DROP-ok  the true fragment is itself under ``--min-read-length``
T+       trimmed at the true length
T-0      trimmed, but the pair has no true overlap
T-ov     trimmed at a wrong length where a true overlap exists
K        kept, no alignment or the true one
Kx       kept with a WRONG alignment detected
LOSTp    kept, but a mate fell under ``--min-read-length`` (trim3 on an N)
======== ============================================================================

Every M-, LOST, T-0, T-ov, Kx and LOSTp pair is kept with weight 1. Each other stratum
keeps its first 50 pairs plus a 0.2% hash-selected sample (by read name, so it does not
depend on file order) and weights them by stratum total / stratum kept.

Output columns: ``bench pair name stratum weight h1 s1 q1 h2 s2 q2 L frag true_r1
true_r2 source``, plus ``<out>.json`` with per-stratum totals. Ported unchanged in
behaviour from the policy study's ``panel/mine.py``.
"""
import argparse
import gzip
import hashlib
import json
import os
import subprocess
import sys
import time
from collections import Counter

#: Overlap-length bins for the stratum label (true overlap = len1 + len2 - L).
BINS = [(0, 0, "none"), (1, 4, "1-4"), (5, 9, "5-9"), (10, 14, "10-14"),
        (15, 19, "15-19"), (20, 29, "20-29"), (30, 49, "30-49"), (50, 99, "50-99"),
        (100, 10 ** 9, "100+")]

#: Classes kept in full (weight 1): everything 0.5.3 got wrong.
WRONG = ("M-", "LOST", "T-0", "T-ov", "Kx", "LOSTp")


def _require_053():
    """Import the 0.5.3 merge API, or exit explaining why this script needs it."""
    try:
        from zna.merge.overlap import find_overlap, reverse_complement, use_backend
        from zna.merge.pairs import process_pair
        from zna.merge.params import MergeParams
        P = MergeParams()
        P.t_merge_q, P.t_trim_q                        # noqa: B018  (0.5.3 fields)
    except (ImportError, AttributeError):
        sys.exit("mine_panel.py records zna 0.5.3's classes (including its trim band) "
                 "and must run under zna 0.5.3; this environment has a newer merge "
                 "policy. Use a pinned 0.5.3 install -- see the module docstring.")
    return find_overlap, reverse_complement, use_backend, process_pair, P


def fq(path):
    """Yield (header, seq, qual) bytes from a FASTQ, decompressing with pigz."""
    cmd = ["pigz", "-dc", path] if path.endswith(".gz") else ["cat", path]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=1 << 20)
    f = p.stdout
    while True:
        h = f.readline()
        if not h:
            break
        s = f.readline().rstrip(b"\n")
        f.readline()
        q = f.readline().rstrip(b"\n")
        yield h[1:].rstrip(b"\n"), s, q
    p.wait()


def truth_rows(path):
    p = subprocess.Popen(["pigz", "-dc", path], stdout=subprocess.PIPE, bufsize=1 << 20)
    head = p.stdout.readline().decode().rstrip("\n").split("\t")
    for line in p.stdout:
        yield dict(zip(head, line.decode().rstrip("\n").split("\t")))
    p.wait()


def geometry_bin(L, l1, l2):
    if L < max(l1, l2):
        return "readthrough"
    if L == max(l1, l2):
        return "full"
    ov = max(0, l1 + l2 - L)
    for lo, hi, name in BINS:
        if lo <= ov <= hi:
            return name


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("sub", help="substrate directory")
    ap.add_argument("bench", help="bench name written into every row")
    ap.add_argument("out", help="output .tsv.gz (a .json summary is written beside it)")
    a = ap.parse_args()
    find_overlap, reverse_complement, use_backend, process_pair, P = _require_053()
    use_backend("accel")
    t0 = time.monotonic()
    r1 = next(p for p in (f"{a.sub}/R1.fq.gz", f"{a.sub}/R1.fq") if os.path.exists(p))
    r2 = r1.replace("R1.fq", "R2.fq")
    totals, sampled, wrong = Counter(), Counter(), Counter()
    rows = []
    for (h1, s1, q1), (h2, s2, q2), t in zip(fq(r1), fq(r2),
                                             truth_rows(f"{a.sub}/truth.tsv.gz")):
        L = int(t["L"])
        l1, l2 = len(s1), len(s2)
        recs, outcome, _dropped, _score, _olen, _diff = process_pair(
            h1, s1, q1, h2, s2, q2, P)
        d, shift, n, _, _ = find_overlap(s1, reverse_complement(s2), P)
        Lhat = (shift + l2) if d == 1 else (n if d == -1 else None)
        geo = geometry_bin(L, l1, l2)
        if outcome == "merged":
            if Lhat == L:
                cls = "M+" if recs else "DROP-ok"    # a genuine molecule under the floor
            else:
                cls = "M-" if recs else "LOST"
        elif outcome == "trimmed":
            cls = "T+" if Lhat == L else ("T-0" if max(0, l1 + l2 - L) == 0 else "T-ov")
        else:
            cls = "K" if (Lhat is None or Lhat == L) else "Kx"
        if not recs and outcome != "merged":
            cls = "LOSTp" if L >= P.min_read_length else "DROP-ok"
        stratum = f"{cls}|{geo}"
        totals[stratum] += 1
        bad = cls in WRONG
        hsel = int.from_bytes(hashlib.blake2b(h1.split()[0], digest_size=8).digest(),
                              "big")
        if bad:
            wrong[cls] += 1
            keep = True
        else:
            keep = hsel % 1000 < 2 or sampled[stratum] < 50
        if keep:
            sampled[stratum] += 1
            rows.append((stratum, bad, h1, s1, q1, h2, s2, q2, t))
    with gzip.open(a.out, "wt", compresslevel=6) as f:
        f.write("bench\tpair\tname\tstratum\tweight\th1\ts1\tq1\th2\ts2\tq2\tL\tfrag"
                "\ttrue_r1\ttrue_r2\tsource\n")
        for stratum, bad, h1, s1, q1, h2, s2, q2, t in rows:
            w = 1.0 if bad else totals[stratum] / sampled[stratum]
            f.write("\t".join([a.bench, t["pair"], t["name"], stratum, f"{w:.6g}",
                               h1.decode(), s1.decode(), q1.decode(), h2.decode(),
                               s2.decode(), q2.decode(), t["L"], t["frag"],
                               t["true_r1"], t["true_r2"], t.get("source", "")]) + "\n")
    summ = dict(bench=a.bench, pairs=sum(totals.values()), kept_rows=len(rows),
                wrong=dict(wrong), strata=dict(totals),
                seconds=round(time.monotonic() - t0, 1))
    with open(a.out.replace(".tsv.gz", ".json"), "w") as fh:
        json.dump(summ, fh, indent=1)
    print(json.dumps({k: summ[k] for k in ("bench", "pairs", "kept_rows", "wrong",
                                           "seconds")}))


if __name__ == "__main__":
    main()
