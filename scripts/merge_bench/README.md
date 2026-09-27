# Evidence for `docs/METHODS.md` and `docs/MERGE_BENCHMARK_RESULTS.md`

> "retired-design §N" below refers to `docs/MERGE_CPP_DESIGN.md`, consolidated
> away in commit 158a204; read it at `git show 158a204^:docs/MERGE_CPP_DESIGN.md`.
> The surviving algorithmic content is `docs/METHODS.md` §2.

Every number in those documents comes from these scripts — except 0.6's full-scale
qualification, whose scripts live beside its data in the policy study's evidence
directory (layer 3 below). They are kept so the design
can be re-argued against measurements rather than recollection, and so the same numbers
can be taken on a Linux/x86 box. That happened in 0.5.2 and it mattered exactly where
the note used to warn it would — `popcount` — see docs/ROADMAP.md, "Closed by measurement
in 0.5.2".

Nothing here is part of the zna package or the test suite.

Three groups of scripts: `mine_panel.py` + `panel_eval.py` judge a change to the merge
**policy** pair by pair; `simulate.py` + `compare.py` (and `compare_zna.py`) measure
**accuracy** end to end against ground truth; and the rest measure **speed** and pin
the C++ design, with `asan_scan.cpp` keeping the kernel honest under sanitizers.

---

## How a merge change is tested: three layers, seconds each

The method 0.6 was built and qualified with (`docs/archive/MERGE_ACCURACY_PLAN.md` §8). Each
layer answers a different question, and none substitutes for another.

**1. Exact tests, no data** — `tests/test_merge.py`. Everything the decision uses is
derived, so it can be pinned exactly: golden values of the floor `T(N)` and the gate's
`dfit`, checked against their definitions; the gate's sensitivity as a closed form,
`P(Binom(n, e′) > dfit[n])`, against exact rational enumeration; the fixed-point scale by
exhaustive enumeration over `(n, d)`; both contract ranges; the argmax total order with
deliberately built ties; the review's 19 ground-truth pairs (`tests/data/report_cases/`);
and the two backends compared for exact equality at every level. Seconds, in all three
test configurations (compiled; `--merge-backend=python`; no extension).

**2. The truth panel and its ledger** — `mine_panel.py`, `panel_eval.py` (below). 245k
pairs mined from 14.1M simulated ones over 46 substrates: every pair zna 0.5.3 got wrong
plus a weighted stratified sample of those it got right, each row carrying its 0.5.3
class. A change is judged by the **pair-level ledger** — which pairs changed outcome,
from what to what — not by a headline rate, because a rate can improve while the pairs
underneath trade one error for another. Seconds on the compiled backend, under a minute
on the reference one.

**3. Full-scale confirmation, at release only.** Regenerate the substrates from their
recipes (deterministic seeds; the recipes are in the policy study's evidence directory,
`data-recipes/`, and each records the sha256 its output must match), run the compiled
merger over every pair — seconds per million — classify each pair against truth for both
the old and the new version, then delete the data. This is where rates that the panel's
sampling cannot resolve are settled: one sampled correct-merge row weighs ~0.05% of a
bench, five times a 0.01% criterion. For 0.6 it ran over 4.11M pairs
(`docs/MERGE_BENCHMARK_RESULTS.md` §9); the scripts are in the evidence directory's
`stageC-results/scripts/`.

A **sealed holdout** backs layers 2 and 3: a panel from genes and a genome draw never
used while designing, looked at exactly once, after the policy and its acceptance
criteria are frozen. 0.6 spent holdout-2; a future policy change needs a fresh one.

---

## Accuracy: `simulate.py` and `compare.py`

The head-to-head against fastp on simulated ground truth; results in
[MERGE_BENCHMARK_RESULTS.md](../../docs/MERGE_BENCHMARK_RESULTS.md).

```bash
# 1M pairs from hg38, ~12 s, ~200 MB of FASTQ plus a 300 MB truth sidecar.
# Deterministic in the seed, so regenerate it rather than storing it.
python simulate.py --genome hg38.fa --out-prefix sim --n-pairs 1000000 \
                   --read-length 150 --frag-min 60 --frag-max 450 \
                   --error-rate 0.002 --quality-model novaseq --seed 42

# runs both tools, scores both against the sidecar, ~25 s
python compare.py --sim-prefix sim --out results/ --threads 4
```

Writes `results/report.md`, `results/summary.json`, and `results/{zna,fastp}_errors.tsv`
— one row per pair the tool got wrong, carrying the truth, what was emitted, and the
evidence the decision was made on. **That file is the point of the exercise**; read it.

Four things that are easy to get wrong here, all of them load-bearing:

- **`--quality-model flat` measures nothing about the consensus.** With a constant
  quality string every mismatch is a tie and R1 always wins; measured, both tools
  recover exactly **0.0%** of recoverable overlap errors. Use `novaseq`, which draws the
  quality first and the error from it at `10^(-Q/10)`.
- **Uniform fragment lengths are not a library.** They exist to populate every geometric
  regime at equal density. Quote per-bin sensitivity, never the overall merge rate.
- **The sidecar is the authority, not the FASTQ comment.** Both tools rewrite headers.
  The read ID up to the first whitespace is the join key. The sidecar records the
  substituted base at every injected error, so it reconstructs the reads exactly.
- **Put the working environment's `bin/` on `PATH`** or `shutil.which("pigz")` returns
  None and everything silently falls back to stdlib gzip.

`compare.py` re-runs `zna merge`'s own decision -- verdict, shift, bits -- at the run's
parameters on every pair that merged wrongly, and scores the *true* shift alongside the
one the tool chose. That is what distinguishes a defective search from an ambiguous
input, and it is the difference between "fix the code" and "the genome repeats here".

**Kept pairs are scored for being whole.** 0.6 has no trim band (plan §2): an unmerged
pair is its two input mates, exactly, and `kept_mate_altered` must be zero. What a kept
pair still carries into the corpus is reported rather than scored -- the overlap both
mates hold, and on a kept read-through, adapter.

`--alpha`, `--error-rate` and `--adapter-trimmed` pass straight through to `zna merge`,
and the re-scan is rebuilt from the same values and checked against the integers the
run's JSON reports. `simulate.py` writes raw adapter read-through, so
`--adapter-trimmed` is a FALSE declaration on its output: it forbids exactly the merges
those pairs need (on 50k chr22 pairs, 11,493 read-throughs kept whole with 1.04M adapter
bases, and zna's read-through check warns at 23%). Use it to price the declaration, not
as the default.

```bash
for E in 0.005 0.01 0.02; do
  python compare.py --sim-prefix sim --out sweep_e$E --error-rate $E --threads 4
done
```

---

## The 0.6 truth panel: `mine_panel.py` and `panel_eval.py`

How a change to the merge **policy** is judged (`docs/archive/MERGE_ACCURACY_PLAN.md` §8, layer
2). The panel is 245k pairs mined from 14.1M simulated ones over 46 substrates: every
pair zna 0.5.3 got wrong, plus a weighted stratified sample of the ones it got right,
each row recording its 0.5.3 class. It lives with the policy study's evidence, not in this
repo (`sealed/` is holdout-2, opened once for 0.6's qualification, and `panel_eval.py` refuses it without
`--final-qualification`).

```bash
# evaluate the CURRENT policy, pair by pair, against the recorded 0.5.3 classes.
# ~13 s on the reference backend with 12 workers; seconds on the compiled one.
python panel_eval.py <panel_dir> --workers 12 --json eval.json --ledger ledger.tsv
python panel_eval.py <panel_dir> --error-rate 0.03      # zna's --error-rate; default 0.01
```

It prints, per bench, 0.5.3 against the current policy: wrong merges retained (M-),
lost fragments (LOST + LOSTp), correct merges (M+), kept pairs, pairs refused as
implausible, and the run's diagnostics -- the detected-overlap disagreement rate, the
share of true overlaps the gate is expected to refuse at that rate (`refused`, through
zna's own closed form over the weighted detected-length histogram), and the read-through
share -- with whether zna would warn on each. **Counts are weighted**:
a sampled correct row stands for up to a few hundred pairs, so one such row changing
class moves a total by its weight -- read the ledger (`--ledger`, every pair whose class
changed) before reading a small difference as a trend. The contract (`--adapter-trimmed`)
per bench and the weighted form of the diagnostics are in the script's docstring.

`mine_panel.py` is the tool that made the panel. It records 0.5.3's classes, trim band
included, so it runs only under a pinned zna 0.5.3 and refuses a newer tree; it is only
needed to add a substrate:

```bash
# in a scratch environment -- never the development one
pip install "zna==0.5.3"
python mine_panel.py <substrate_dir> <bench_name> <bench_name>.tsv.gz
```

A substrate is `R1.fq[.gz]`, `R2.fq[.gz]` and `truth.tsv.gz` (one row per pair in FASTQ
order: `pair name L frag true_r1 true_r2 source`, the fragment in R1's orientation and
`source` a locus id for concentration statistics). Add the new bench's contract to
`panel_eval.py`'s docstring and table: `--adapter-trimmed` where every read ends at or
before its molecule, undeclared otherwise.

**Reading a result.** Compare against 0.5.3 bench by bench, then read the ledger for the
pairs that moved. A useful check before trusting a new policy: run `panel_eval.py
--backend accel` and `--backend python` and diff the two ledgers — they must be
identical, and for 0.6 they are.

---

## Speed: reproducing the C++ design

```bash
cd scripts/merge_bench

# 1. a representative library: 2x150, insert ~ N(200,70) in [50,400], 0.5% error.
#    Merges at 88.6% against production's measured 88.8%.
python gen_library.py 200000 r1.fq.gz r2.fq.gz

# 2. is the argmax tie-break a specifiable total order?  (retired-design §5)
#    Builds ties deliberately -- random sequence essentially never ties.
python verify_tiebreak.py r1.fq.gz r2.fq.gz

# 3. scan kernel variants, full pruned scans, packing inside the timed region,
#    every variant checked against the shipped kernel pair by pair  (retired-design §6.1)
python dump_pairs.py r1.fq.gz r2.fq.gz 50000 pairs.bin
c++ -O3 -std=c++17 -o bench_scan bench_scan.cpp && ./bench_scan pairs.bin   # scalar vs 2-bit packed
c++ -O3 -std=c++17 -o bench_simd bench_simd.cpp && ./bench_simd pairs.bin   # ...vs byte-wise SIMD
# on x86 see "Taking these on x86" below before adding -mavx2: it enables POPCNT too,
# which silently folds the reduction fix into what looks like a vector-width result
```

**Steps 2 and 3 measure the SCAN, not the 0.6 decision.** `dump_pairs.py` and
`verify_tiebreak.py` call the shipped backend's `scan` over every shift at the fixed
8-bit floor the C++ benches hard-code, with the `e = 0.01` weights; 0.6 wraps that same
loop in a per-pair floor, the `--adapter-trimmed` contract and the plausibility gate,
none of which changes which kernel variant is faster or whether it finds the same
argmax. Ported to the 0.6 API, `dump_pairs.py` writes a `pairs.bin` byte-identical to
the one the 0.5.x script wrote, and `bench_scan` reports 0 mismatches against it;
`verify_tiebreak.py` finds 144 exercised ties and 0 violations of (max score, max n,
min s), scoring in the kernel's integers rather than the floats it used before.

## Retired in 0.6

Deleted because they measured the 0.5.x merge -- its trim band, its two fixed
thresholds, its kernel API -- and have no meaningful 0.6 form. Each is recoverable at
the 0.5.3 tag (`git show v0.5.3:scripts/merge_bench/<file>`).

- **`bench_breakdown.py`** -- cumulative per-stage timing of the 0.5.x Python path
  (retired-design §1), down to `merge_chunk` called with `t_trim_q`. The finding it made
  -- per-pair work belongs in one GIL-releasing chunk call -- is the shipped design; the
  whole run's rate is `pairs_per_second` in `zna merge --json`.
- **`proto_merge.cpp`** -- the whole 0.5.x path in one C++ file, whose point was to emit
  a file byte-identical to 0.5.x `zna merge` (2.32 µs/pair against 8.34 on 200,000
  pairs; retired-design §3 and §7.4). The compiled backend superseded it in 0.4.0, and
  what it proved is now proved continuously: `tests/test_merge.py` holds the compiled
  backend byte-identical to the reference one, and `panel_eval.py --backend accel` vs
  `python` on the truth panel. Porting it would mean re-implementing the 0.6 policy a
  third time with nothing to verify it against but the other two.

## Sanitizers: `asan_scan.cpp`

Runs the header-only core (`merge_core.hpp`, `fastq_chunk.hpp`) under AddressSanitizer
and UBSan, with every read placed so its last byte is the last byte of its allocation and
every policy table allocated at exactly its capacity, so a one-byte overread or a table
lookup one past the end traps instead of landing in slack. Build and run it after any
kernel change (the command is in the file's header); it reports "… clean" or aborts.

## One number that is easy to get wrong

`bench_scan.cpp` measures **with packing inside the timed loop**. An earlier version
hoisted it, which flattered the packed variant by ~20% and would have hidden the finding
that a table-driven packer costs the entire win (1.027 vs 0.631 µs/pair). The audit made
the same mistake in the other direction — its 4.8x SWAR figure was a no-bail sweep doing
~5x the real pruned scan's work. Measure the whole scan, with everything it needs.

## Reproducing the tie-break and fixed-point proofs

`verify_tiebreak.py` builds ties deliberately — random sequence essentially never ties,
which is why an earlier sweep found none and proved nothing. The complementary argument
is arithmetic: a tie across *different* overlap lengths needs `dn = STEP/gcd(M,STEP)`,
which is ~1.4e8 for the shipped weights, so ties only ever occur at equal `n`. Both are
retired-design §5.

The fixed-point scale is chosen by exhaustive enumeration rather than taste — see design
§4 for the table, which is small enough to live as a test.

## Taking these on x86 — done in 0.5.2, and what it found

The SSE2 path in `bench_simd.cpp` had never been run. It runs, and it showed the shipped
kernel getting **1.33×** over scalar on a Xeon E5-2680 v3 against the **2.29×** these
scripts recorded on NEON. The cause was `__builtin_popcount` compiling to
`callq __popcountdi2@plt` — POPCNT is not baseline x86-64 — twice per 32 bases in the
innermost loop. `merge_core.hpp` now reduces with `psadbw` instead.

**`bench_simd.cpp` still contains that popcount kernel**, deliberately: it is the
`before` in the comparison. Two things follow for anyone re-running it on x86.

- **`-mavx2` also enables POPCNT**, so a bare `./bench_simd` against a `-mavx2` build
  compares *two* changes at once and will credit the vector width for the reduction's
  win. Build with `-mpopcnt` alone to separate them; the numbers are in ROADMAP.
- The bail granularity is per-ISA now (`BAIL_VECTORS`): 48 bases on x86, 64 on aarch64.

```bash
c++ -O3 -std=c++17 -o bench_simd bench_simd.cpp && ./bench_simd pairs.bin  # as shipped
c++ -O3 -std=c++17 -mpopcnt  -o bs_pc bench_simd.cpp && ./bs_pc pairs.bin  # isolate popcnt
c++ -O3 -std=c++17 -mavx2    -o bs_v2 bench_simd.cpp && ./bs_v2 pairs.bin  # + 32 B vectors
```

The aarch64 re-measure this section once scheduled was 0.5.3: `BAIL_VECTORS` went
2 -> 4 there (1.22x on the kernel); see docs/ROADMAP.md.
