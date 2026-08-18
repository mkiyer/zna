# Evidence for `docs/METHODS.md` and `docs/MERGE_BENCHMARK_RESULTS.md`

Every number in those documents comes from these scripts. They are kept so the design
can be re-argued against measurements rather than recollection, and so the same numbers
can be taken on a Linux/x86 box — this session's were all aarch64, which matters for
`popcount`. See docs/ROADMAP.md, "0.4.2 — x86 tuning".

Nothing here is part of the zna package or the test suite.

Two groups of scripts: `simulate.py` + `compare.py` measure **accuracy** against ground
truth, and everything else measures **speed** and pins the C++ design.

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

`compare.py` re-runs `zna merge`'s own kernel on every pair that merged wrongly, and
scores the *true* shift alongside the one the tool chose. That is what distinguishes a
defective search from an ambiguous input, and it is the difference between "fix the
code" and "the genome repeats here".

**Trimming is scored as its own contract, not as a footnote to merging.** Most pairs in
a real library do not merge, and an unmerged pair is still encoded, so whether its
redundant overlap came off matters to the corpus exactly as much as a merge does. Two
traps worth knowing:

- **Score every kept pair, not just the overlapping ones.** A trim applied to a pair
  whose mates share nothing deletes real sequence, and an analysis restricted to pairs
  with a true overlap cannot see it. There are 2,733 of those in 1M pairs.
- **Report the counterfactual.** "28,074 duplicated bases survived" means nothing without
  "283,838 would have, untrimmed". `compare.py` emits both, plus the bases deleted, so
  the trade is legible rather than asserted.

`compare.py --threshold-trim T` passes the value straight through to `zna merge`, which
is how the band gets priced — the merge decision does not move with it, so a sweep
isolates the trim cleanly:

```bash
for T in 8 12 16 20; do
  python compare.py --sim-prefix sim --out sweep_t$T --threshold-trim $T --threads 4
done
```

---

## Speed: reproducing the C++ design

```bash
cd scripts/merge_bench

# 1. a representative library: 2x150, insert ~ N(200,70) in [50,400], 0.5% error.
#    Merges at 88.6% against production's measured 88.8%.
python gen_library.py 200000 r1.fq.gz r2.fq.gz

# 2. where the per-pair time actually goes  (design §1)
python bench_breakdown.py r1.fq.gz r2.fq.gz

# 3. is the argmax tie-break a specifiable total order?  (design §5)
#    Builds ties deliberately -- random sequence essentially never ties.
python verify_tiebreak.py r1.fq.gz r2.fq.gz

# 4. scan kernel variants, full pruned scans, packing inside the timed region,
#    every variant checked against the shipped kernel pair by pair  (design §6.1)
python dump_pairs.py r1.fq.gz r2.fq.gz 50000 pairs.bin
c++ -O3 -std=c++17 -o bench_scan bench_scan.cpp && ./bench_scan pairs.bin   # scalar vs 2-bit packed
c++ -O3 -std=c++17 -o bench_simd bench_simd.cpp && ./bench_simd pairs.bin   # ...vs byte-wise SIMD
# on x86, add -mavx2 to bench_simd to enable the 32-byte variants (SSE2 needs no flag)

# 5. the whole path in C++, which must be BYTE-IDENTICAL to `zna merge`  (design §3)
c++ -O3 -std=c++17 -o proto_merge proto_merge.cpp
./proto_merge r1.fq.gz r2.fq.gz cpp.fq 40
zna merge --in1 r1.fq.gz --in2 r2.fq.gz --out py.fq --min-read-length 40 -q
cmp cpp.fq py.fq        # must be silent
```

## `proto_merge.cpp` is a prototype, not the implementation

It is the artifact that proves the design lands where it claims: on 200,000 real pairs
it emits a byte-identical file to `zna merge` at **2.32 µs/pair against 8.34**. It
carries the byte-wise SIMD kernel the design settled on (`neq16` + a 32-base bail) and
the self-sizing `Scratch` arena of design §7.4, and was re-verified byte-identical after
each of those replaced its predecessor — which is the check that matters when swapping a
kernel or an allocation strategy.

Do **not** read it as a model for the real backend. It is single-threaded, its
`popen("pigz -dc")` reader is a stand-in for the chunk protocol in design §7.3, and it
builds the consensus table in C++ where design §4 says to build it once in Python and
pass it across. It also has no error handling worth the name.

`/tmp/proto_v1.cpp`-style comparisons are how §7.4's arena numbers were taken: the
same program with per-pair `std::vector` `.assign()`/`.resize()` and unconditional
consensus copies runs at 2.71 µs/pair end to end and 1.38 µs/pair of compute, against
the arena's 2.32 and 1.00.

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
design §5.

The fixed-point scale is chosen by exhaustive enumeration rather than taste — see design
§4 for the table, which is small enough to live as a test.

## Taking these on x86

Every number in the design is aarch64/NEON. The SSE2 path in `bench_simd.cpp` is written
but has never been run, the AVX2 variants need `-mavx2`, and the optimal bail
granularity is hardware-dependent (32 bases won here). This is design phase G.
