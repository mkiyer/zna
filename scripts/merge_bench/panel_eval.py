"""Evaluate the CURRENT zna merge policy on the truth panel, pair by pair.

The panel (``docs/archive/MERGE_ACCURACY_PLAN.md`` §8, layer 2) is 46 benches mined by
``mine_panel.py`` from 14.1M simulated pairs: every pair zna 0.5.3 got wrong, plus a
weighted stratified sample of the pairs it got right. Each row records its 0.5.3 class
in ``stratum`` (``<class>|<overlap bin>``) and a ``weight`` that projects it back onto
the full substrate (1 for every 0.5.3-wrong pair). So one pass over ~245k rows says,
for the whole 14.1M, what the current policy does differently -- and a change is judged
by the **pair-level ledger**: which pairs moved, from what to what.

    PY=/Users/mkiyer/sw/miniforge3/envs/zna_merge/bin/python
    $PY scripts/merge_bench/panel_eval.py <panel_dir> --workers 12 \\
        --json eval.json --ledger ledger.tsv

Runs whichever merge backend zna selects (``--backend``); on the reference backend the
dev panel takes well under a minute with 12 worker processes. ``--error-rate`` is the
run's ``--error-rate`` for every bench (default 0.01, zna's default); pass another value
to see what setting it does.

**The sealed panel is refused.** ``panel/sealed/`` is holdout-2, to be looked at once, at
final qualification; pass ``--final-qualification`` to open it, and only then.

Per pair
--------

1. **The decision**, under the bench's contract (below), with ``--error-rate``,
   ``alpha`` = 1e-6, ``--min-read-length 40`` and ``--npolicy trim3``, as mined.
2. **The run's two diagnostics**, as zna would report them (``MERGE_ACCURACY_PLAN.md``
   §4) -- but weighted, so they project onto the whole substrate:

   * ``detected_overlap_mismatch_rate``: informative mismatches over informative
     positions of every overlap the scan detected, before the gate
     (``PairResult.detected_*``), each row times its weight.
   * ``expected_refused_true_overlap_fraction`` (``refused``): from that rate and the
     weighted histogram of detected overlap lengths (``PairResult.detected_overlap_len``),
     through zna's own :func:`zna.merge.cli.expected_refused_fraction` -- the share of
     true overlaps the gate is expected to refuse at ``--error-rate``. zna warns above
     0.1%; the ``warn e`` column says whether this bench's run would.
   * the read-through check: the weighted share of rows whose UNRESTRICTED best
     alignment (:func:`zna.merge.overlap.scan_unrestricted`) is a read-through reaching
     ``T``. zna computes it over the first 100,000 pairs; here it is over the whole
     projected substrate, which for a FASTQ sorted by fragment length (chr22) is the
     unbiased value. ``warn rt`` says whether a declared run would warn (``> 0.01``).
3. **The class**, against truth (``L``, the true fragment length):

   ======== ===========================================================================
   M+       merged at the true length, record emitted
   M-       merged at a WRONG length, record emitted -- a false full-fragment claim
   LOST     merged at a wrong length that fell under the length filter: the fragment
            is gone
   DROP-ok  the true fragment is under ``min_read_length``; dropping it is correct
   K        kept as two mates
   LOSTp    kept, but a mate fell under the length filter (trim3 on an N), so the
            fragment was dropped although it is long enough
   ======== ===========================================================================

   The 0.5.3 classes also include T+ / T-0 / T-ov (trimmed right / trimmed with no true
   overlap / trimmed the wrong overlap) and Kx (kept with a wrong alignment detected);
   0.6 has no trim band, so its unmerged pairs are all K.

The contract
------------

``adapter_trimmed`` is True where every read ends at or before its molecule's end,
because the substrate either has no adapters or was adapter-trimmed:
``chr22-*`` (khorana's simulator clips reads to the molecule), ``tx-*-clean``,
``txf-clean`` (the khorana profile), ``tx-*-fastp``, ``txf-fastp``, ``sp3hgfp``,
``pglossfp`` (passed through fastp), ``cutq30``, ``cutramp``, ``cutright`` (cut to the
molecule). Every other bench carries raw adapter read-through and runs without it.
"""
from __future__ import annotations

import argparse
import gzip
import json
import multiprocessing as mp
import os
import sys
import time
from collections import Counter
from pathlib import Path

from fractions import Fraction

from zna.merge import backend as zbackend
from zna.merge import cli as zcli
from zna.merge.overlap import reverse_complement, scan_unrestricted
from zna.merge.pairs import process_pair
from zna.merge.params import DEFAULT_ERROR_RATE, MergeParams, decimal_str

MIN_READ_LENGTH = 40
ALPHA = "1e-6"

#: zna's own warning thresholds (zna.merge.cli.run_warnings): the expected share of
#: true overlaps refused, and the read-through share.
WARN_REFUSED = zcli._WARN_REFUSED
WARN_READTHROUGH = float(zcli._WARN_READTHROUGH)

#: 0.5.3's classes that are errors, and 0.6's. T+ and Kx are "not wrong" and "wrong"
#: respectively in 0.5.3's own accounting (mine_panel.py keeps every Kx).
WRONG_053 = ("M-", "LOST", "T-0", "T-ov", "Kx", "LOSTp")


def adapter_trimmed_for(bench: str) -> bool:
    """The bench's contract; see the module docstring."""
    if bench.startswith("chr22-"):
        return True
    if bench.startswith("tx-") and (bench.endswith("-clean") or bench.endswith("-fastp")):
        return True
    return bench in {"txf-clean", "txf-fastp", "sp3hgfp", "pglossfp", "cutq30",
                     "cutramp", "cutright"}


def locus_of(source: str) -> str:
    """A locus key for concentration: the transcript, or a 100 kb genomic bin."""
    src = source.split("|", 1)[0]
    head = src.split(":", 1)[0]
    if head.startswith("ENST") or head.startswith("ENSG"):
        return head
    parts = src.split(":")
    if len(parts) >= 2 and parts[1].isdigit():
        return f"{parts[0]}:{int(parts[1]) // 100_000}"
    return src


# --------------------------------------------------------------------------- #
# data, loaded once in the parent and inherited by forked workers
# --------------------------------------------------------------------------- #

ROWS: list = []          # (bench_idx, pair, stratum, weight, h1, s1, q1, h2, s2, q2, L, source)
BENCHES: list = []       # bench names, indexed by bench_idx
PARAMS: dict = {}        # bench_idx -> MergeParams, set before the workers fork


def load(paths):
    for path in paths:
        with gzip.open(path, "rt") as fh:
            head = fh.readline().rstrip("\n").split("\t")
            col = {k: i for i, k in enumerate(head)}
            bench_idx = None
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if bench_idx is None:
                    BENCHES.append(f[col["bench"]])
                    bench_idx = len(BENCHES) - 1
                ROWS.append((bench_idx, f[col["pair"]], f[col["stratum"]],
                             float(f[col["weight"]]),
                             f[col["h1"]].encode(), f[col["s1"]].encode(),
                             f[col["q1"]].encode(), f[col["h2"]].encode(),
                             f[col["s2"]].encode(), f[col["q2"]].encode(),
                             int(f[col["L"]]), f[col["source"]]))


def _init(backend_name):
    zbackend.use(backend_name)


def _classify(outcome, records, Lhat, L):
    if outcome == "merged":
        if Lhat == L:
            return "M+" if records else "DROP-ok"
        return "M-" if records else "LOST"
    if not records:
        return "LOSTp" if L >= MIN_READ_LENGTH else "DROP-ok"
    return "K"


def _decide_chunk(bounds):
    """The decision, class and diagnostics of every row in [lo, hi)."""
    lo, hi = bounds
    out = []
    for i in range(lo, hi):
        b, pair, stratum, w, h1, s1, q1, h2, s2, q2, L, src = ROWS[i]
        p = PARAMS[b]
        res = process_pair(h1, s1, q1, h2, s2, q2, p)
        Lhat = res.shift + len(s2) if res.olen else None
        new = _classify(res.outcome, res.records, Lhat, L)
        a = scan_unrestricted(s1, reverse_complement(s2), p)
        rt = bool(a.overlap_len) and a.shift + len(s2) < max(len(s1), len(s2))
        out.append((i, new, Lhat, res.implausible, res.detected_bases,
                    res.detected_mismatches, res.detected_overlap_len, rt))
    return out


def _chunks(n, size):
    return [(lo, min(n, lo + size)) for lo in range(0, n, size)]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("panel", help="panel directory (its *.tsv.gz) or TSV files",
                    nargs="+")
    ap.add_argument("--workers", type=int, default=min(12, os.cpu_count() or 1))
    ap.add_argument("--backend", default="auto", choices=("auto", "accel", "python"))
    ap.add_argument("--error-rate", default=DEFAULT_ERROR_RATE,
                    help="zna's --error-rate, for every bench (default: zna's default)")
    ap.add_argument("--bench", action="append", default=None,
                    help="only these benches (repeatable)")
    ap.add_argument("--json", help="write the full per-bench result here")
    ap.add_argument("--ledger", help="write every pair whose class changed (TSV)")
    ap.add_argument("--final-qualification", action="store_true",
                    help="allow the sealed holdout panel: final qualification, once.")
    a = ap.parse_args(argv)
    try:
        e = MergeParams(alpha=ALPHA, error_rate=a.error_rate).e
    except ValueError as err:
        sys.exit(f"--error-rate: {err}")

    paths = []
    for p in map(Path, a.panel):
        paths.extend(sorted(p.glob("*.tsv.gz")) if p.is_dir() else [p])
    if any("sealed" in q.parts for q in paths) and not a.final_qualification:
        sys.exit("refusing the sealed holdout panel; see --final-qualification")
    if a.bench:
        paths = [q for q in paths if q.name[:-len(".tsv.gz")] in set(a.bench)]
    if not paths:
        sys.exit("no panel TSVs found")

    t0 = time.monotonic()
    load(paths)
    backend = zbackend.get_merge_backend_name(None if a.backend == "auto" else a.backend)
    t_load = time.monotonic() - t0

    for b in range(len(BENCHES)):
        PARAMS[b] = MergeParams(alpha=ALPHA, error_rate=a.error_rate,
                                adapter_trimmed=adapter_trimmed_for(BENCHES[b]),
                                min_read_length=MIN_READ_LENGTH, npolicy="trim3")

    # ---- decisions ------------------------------------------------------------
    ctx = mp.get_context("fork")
    chunks = _chunks(len(ROWS), 1500)
    t1 = time.monotonic()
    results = [None] * len(ROWS)
    with ctx.Pool(a.workers, initializer=_init, initargs=(backend,)) as pool:
        for part in pool.imap_unordered(_decide_chunk, chunks):
            for i, *rest in part:
                results[i] = rest
    t_dec = time.monotonic() - t1

    # ---- aggregate -----------------------------------------------------------
    per = {}
    ledger_rows = []
    for b, name in enumerate(BENCHES):
        per[name] = dict(
            bench=name, adapter_trimmed=adapter_trimmed_for(name),
            error_rate=decimal_str(e),
            det_bases_w=0.0, det_mismatches_w=0.0, readthrough_w=0.0, det_len_w=Counter(),
            rows=0, weight=0.0,
            old=Counter(), old_w=Counter(), new=Counter(), new_w=Counter(),
            ledger=Counter(), ledger_w=Counter(), implausible=0, implausible_w=0.0,
            wrong_loci=Counter())
    for i, r in enumerate(ROWS):
        b, pair, stratum, w, _h1, s1, _q1, _h2, s2, _q2, L, src = r
        new, Lhat, impl, det_n, det_d, det_len, rt = results[i]
        old = stratum.split("|", 1)[0]
        d = per[BENCHES[b]]
        d["rows"] += 1
        d["weight"] += w
        d["det_bases_w"] += w * det_n
        d["det_mismatches_w"] += w * det_d
        if det_len:
            d["det_len_w"][det_len] += w
        d["readthrough_w"] += w * rt
        d["old"][old] += 1
        d["old_w"][old] += w
        d["new"][new] += 1
        d["new_w"][new] += w
        d["ledger"][(old, new)] += 1
        d["ledger_w"][(old, new)] += w
        if impl:
            d["implausible"] += 1
            d["implausible_w"] += w
        if new in ("M-", "LOST"):
            d["wrong_loci"][locus_of(src)] += w
        if old != new:
            ledger_rows.append((BENCHES[b], pair, old, new, w, L,
                                "" if Lhat is None else Lhat, int(impl), src))

    def unmerged_053(c):
        return sum(c[k] for k in ("K", "Kx", "T+", "T-0", "T-ov"))

    summary = []
    for name in BENCHES:
        d = per[name]
        ow, nw = d["old_w"], d["new_w"]
        loci = d["wrong_loci"]
        total_wrong = sum(loci.values())
        top10 = sum(v for _k, v in loci.most_common(10))
        det = d["det_mismatches_w"] / d["det_bases_w"] if d["det_bases_w"] else 0.0
        rt = d["readthrough_w"] / d["weight"] if d["weight"] else 0.0
        # zna's closed form over the weighted histogram, at the weighted rate (as an
        # exact rational of the float, so the same function serves both).
        hist = [0.0] * (max(d["det_len_w"], default=0) + 1)
        for n, w in d["det_len_w"].items():
            hist[n] = w
        p = PARAMS[BENCHES.index(name)]
        p.ensure(len(hist) - 1)
        refused = zcli.expected_refused_fraction(Fraction(det), hist, p.dfit_table)
        summary.append(dict(
            bench=name, contract="trimmed" if d["adapter_trimmed"] else "raw",
            e=d["error_rate"],
            detected_rate=round(det, 6),
            expected_refused=float(format(refused, ".4g")),
            warn_error_rate=Fraction(refused) > WARN_REFUSED,
            readthrough_check=round(rt, 6),
            warn_readthrough=d["adapter_trimmed"] and rt > WARN_READTHROUGH,
            wrong_053=ow["M-"], wrong_06=nw["M-"],
            lost_053=ow["LOST"] + ow["LOSTp"], lost_06=nw["LOST"] + nw["LOSTp"],
            correct_053=ow["M+"], correct_06=nw["M+"],
            kept_053=unmerged_053(ow), kept_06=nw["K"],
            implausible_06=d["implausible_w"],
            wrong_loci=len(loci),
            wrong_top10_share=round(top10 / total_wrong, 3) if total_wrong else 0.0,
            wrong_top_loci=[(k, round(v, 1)) for k, v in loci.most_common(5)]))

    # ---- report --------------------------------------------------------------
    hdr = ("bench", "contract", "e", "M- 053", "M- 06", "lost 053", "lost 06",
           "M+ 053", "M+ 06", "kept 053", "kept 06", "impl 06", "det rate", "refused",
           "warn e", "rt chk", "warn rt")
    print("\t".join(hdr))
    tot = Counter()
    for s in summary:
        vals = [s["bench"], s["contract"], s["e"]]
        for k in ("wrong_053", "wrong_06", "lost_053", "lost_06", "correct_053",
                  "correct_06", "kept_053", "kept_06", "implausible_06"):
            vals.append(f"{s[k]:.0f}")
            tot[k] += s[k]
        vals += [f"{s['detected_rate']:.6f}", f"{s['expected_refused']:.3g}",
                 "yes" if s["warn_error_rate"] else "-",
                 f"{s['readthrough_check']:.4f}",
                 "yes" if s["warn_readthrough"] else "-"]
        print("\t".join(vals))
    print("\t".join(["TOTAL", "", ""] + [f"{tot[k]:.0f}" for k in (
        "wrong_053", "wrong_06", "lost_053", "lost_06", "correct_053", "correct_06",
        "kept_053", "kept_06", "implausible_06")] + ["", "", "", "", ""]))
    print(f"\n{len(ROWS)} rows, {len(BENCHES)} benches, backend {backend}, "
          f"--error-rate {decimal_str(e)}, {a.workers} workers: load {t_load:.1f}s, "
          f"decide {t_dec:.1f}s, total {time.monotonic() - t0:.1f}s", file=sys.stderr)

    if a.json:
        def plain(d):
            return {("->".join(k) if isinstance(k, tuple) else k): v for k, v in d.items()}
        out = dict(backend=backend, alpha=ALPHA, error_rate=decimal_str(e),
                   min_read_length=MIN_READ_LENGTH, summary=summary,
                   benches={n: dict(
                       {k: v for k, v in d.items() if not isinstance(v, Counter)},
                       old=dict(d["old"]), old_w=dict(d["old_w"]),
                       new=dict(d["new"]), new_w=dict(d["new_w"]),
                       ledger=plain(d["ledger"]), ledger_w=plain(d["ledger_w"]),
                       wrong_loci=dict(d["wrong_loci"].most_common(50)))
                       for n, d in per.items()})
        Path(a.json).write_text(json.dumps(out, indent=1))
    if a.ledger:
        with open(a.ledger, "w") as fh:
            fh.write("bench\tpair\tclass_053\tclass_06\tweight\tL\tL_hat\timplausible"
                     "\tsource\n")
            for row in ledger_rows:
                fh.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    main()
