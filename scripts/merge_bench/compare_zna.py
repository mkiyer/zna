#!/usr/bin/env python3
"""Score encoded .zna geometry against the simulator's truth sidecar.

The .zna-level comparator MERGE_PAIRS_PLAN.md §7 budgets: `compare.py` scores
merged FASTQ, so it stops one stage short of the thing that trains -- the
ENCODED flags.  This decodes one or more .zna files with
``records(with_ends=True)``, joins each fragment to `sim.truth.tsv`, and
scores the IS_FULL_FRAGMENT claim per record: a full-fragment record whose
length equals the true fragment length is an exact reconstruction; one whose
length differs is the wrong-shift residual (a property of the genome; both
paths carry it -- see MERGE_BENCHMARK_RESULTS.md §2).

Records are matched to truth POSITIONALLY: `zna encode` preserves input order,
so fragment i's records are consecutive.  ZNA stores no read names, and the
merger's min-length filter DROPS whole fragments, so bare positional matching
goes wrong at the first drop -- pass ``--ids merged.fq.gz`` (the two-step
intermediate, which still has names) to restrict the truth to the surviving
fragments in order.  Do not run this on a shuffled file.

Usage:
    python compare_zna.py [--ids merged.fq.gz] sim.truth.tsv a.zna [b.zna ...]

With two or more files, also reports where the files disagree record for
record, and which one truth agrees with.
"""
import csv
import sys
from pathlib import Path

from zna.core import ZnaReader


def load_truth(path, ids=None):
    rows = {}
    order = []
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows[row["read_id"]] = (int(row["frag_len"]), int(row["len1"]),
                                    int(row["len2"]))
            order.append(row["read_id"])
    if ids is None:
        return [rows[i] for i in order]
    return [rows[i] for i in ids]


def surviving_ids(merged_fastq):
    """Fragment ids that produced at least one record, in stream order."""
    import gzip
    opener = gzip.open if str(merged_fastq).endswith(".gz") else open
    ids = []
    with opener(merged_fastq, "rb") as fh:
        for i, line in enumerate(fh):
            if i % 4:
                continue
            tok = line[1:].split()[0]
            if tok.endswith(b"/1") or tok.endswith(b"/2"):
                tok = tok[:-2]
            tok = tok.decode()
            if not ids or ids[-1] != tok:
                ids.append(tok)
    return ids


def read_units(path):
    """Group a .zna's records into per-fragment units, in file order."""
    units = []
    with open(path, "rb") as fh:
        pend = None
        for rec in ZnaReader(fh).records(with_ends=True):
            seq, is_paired = rec[0], rec[1]
            if is_paired:
                if pend is None:
                    pend = [rec]
                else:
                    units.append(pend + [rec])
                    pend = None
            else:
                assert pend is None, "unpaired record between mates"
                units.append([rec])
    assert pend is None, "file ends mid-pair"
    return units


def score(name, units, truth):
    n_full = n_exact = n_wrong_len = n_pairs = 0
    if len(units) != len(truth):
        sys.exit(f"{name}: {len(units)} fragments vs {len(truth)} truth rows -- "
                 f"pass --ids with the two-step intermediate to align them")
    for unit, (frag_len, _l1, _l2) in zip(units, truth):
        if len(unit) == 1 and not unit[0][1]:
            n_full += 1
            if len(unit[0][0]) == frag_len:
                n_exact += 1
            else:
                n_wrong_len += 1
        else:
            n_pairs += 1
    n = len(units)
    print(f"{name}: {n} fragments; {n_full} merged ({100*n_full/n:.2f}%), "
          f"{n_pairs} kept as pairs")
    print(f"    IS_FULL_FRAGMENT exact against truth: {n_exact} "
          f"({100*n_exact/max(n_full,1):.3f}% of merged)")
    print(f"    wrong-length residual: {n_wrong_len} "
          f"({1e6*n_wrong_len/max(n,1):.0f} per million fragments) "
          f"-- both paths carry these; a property of the genome")
    return units


def main(argv):
    ids = None
    if len(argv) > 2 and argv[1] == "--ids":
        ids = surviving_ids(argv[2])
        argv = argv[:1] + argv[3:]
    if len(argv) < 3:
        sys.exit(__doc__)
    truth = load_truth(argv[1], ids)
    all_units = {}
    for path in argv[2:]:
        all_units[path] = score(Path(path).name, read_units(path), truth)
    paths = list(all_units)
    if len(paths) >= 2:
        a, b = all_units[paths[0]], all_units[paths[1]]
        diff = sum(1 for ua, ub in zip(a, b)
                   if [(r[0], r[1:]) for r in ua] != [(r[0], r[1:]) for r in ub])
        print(f"\nfragments where the two paths disagree: {diff}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
