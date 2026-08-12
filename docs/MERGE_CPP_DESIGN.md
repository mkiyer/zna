# Native acceleration for `zna merge`: design

Status: **proposed**, 2026-08-12, against zna at the read-merge port (0.3.6-unreleased).
Target release: **0.4.0** — see §12.

This is the design for replacing the numba overlap kernel with a C++ backend, keeping a
pure-Python implementation as the **reference oracle** rather than a fallback, the way
`_pycodec` and `_accel` already relate.

Every number below was measured on this repo, on 200,000 simulated 2×150 pairs with
insert ~ N(200, 70) truncated to [50, 400] and 0.5% per-base error. That library merges
at **88.6%** against production's measured 88.8%, and its overlap mismatch rate is
0.0075 against production's 0.0091, so it is representative. Machine: Apple M-series,
Python 3.12, numba 0.66, Apple clang 21, `-O3 -std=c++17`, min-of-N.

Companion reading: [READ_MERGE_REDESIGN.md](READ_MERGE_REDESIGN.md) (the scoring rule —
settled, do not re-derive), [MERGE_TOOL_AUDIT.md](MERGE_TOOL_AUDIT.md) §(c) (ideas
measured and rejected), [READ_MERGE_PORT_TO_ZNA.md](READ_MERGE_PORT_TO_ZNA.md) §4 (the
contract C1–C7), [READ_MERGE_ROADMAP.md](READ_MERGE_ROADMAP.md) (status board).

---

## 1. The measurement that reframes the problem

The roadmap's plan for this work is "numba → C++", and it discusses only the scan. That
is not where the time is. Measured per-pair cost of `zna merge --processes 1`:

| stage | µs/pair | share |
|---|---:|---:|
| gzip read + FASTQ parse (wall; pigz inflates in subprocesses) | 0.90 | 11% |
| `revcomp(R2)` | 0.17 | 2% |
| **overlap scan (numba kernel)** | **2.63** | **32%** |
| **posterior consensus `_consensus_r1_overlap`** | **2.33** | **28%** |
| `_build_merged` + `q2[::-1]` | 0.47 | 6% |
| merged-name construction | 0.37 | 4% |
| R1/R2 sync check (`base_name` ×2) | 0.48 | 6% |
| output formatting + histograms | 0.59 | 7% |
| write (pigz -1) | 0.04 | <1% |
| **total** | **8.34** | |
| for reference, `--processes 8` | 1.97 | |

Two things fall out immediately.

**The consensus is as expensive as the scan.** `_consensus_r1_overlap` is a pure-Python
per-position loop that runs on the 43.5% of pairs whose winning overlap carries at least
one mismatch — 5.35 µs per pair that enters it. It is 28% of the tool and the roadmap
does not mention it. **A port that moves only the kernel leaves more than half of the
per-pair Python cost in place**, and would deliver about 1.4x.

**At production settings the merge computation is already nearly free.** The audit
measured p4 ≈ p8 ≈ p12 and concluded the workers are idle; that reproduces here (1.97
µs/pair at p8). So the honest framing of "extremely fast" is not "make the kernel
faster" — it is:

- **collapse 8 cores into 1**, and
- **remove the FASTQ intermediate entirely** (phase 2, `--merge-pairs`), which is where
  the compounding win actually lives.

Anything that only speeds the kernel is measured against a tool that is already I/O
bound at production settings. Say so up front rather than discovering it after.

---

## 2. The floor

| configuration | µs/pair |
|---|---:|
| gzip inflate floor — 2 × `pigz -dc -p1`, concurrent, to `/dev/null` | **0.80** |
| same with `-p4` each (pigz cannot parallelise *inflate*) | 0.95 |
| C++ prototype: gz in → file out, single thread | 2.79 |
| C++ prototype: gz in → `/dev/null` | 2.69 |
| C++ prototype: **plain** in → `/dev/null` (compute only) | **1.43** |
| `zna merge -p1` (numba) | 8.34 |
| `zna merge -p8` (numba) | 1.97 |

The writer is free (2.79 vs 2.69 — pigz level 1 in a subprocess keeps up). Decompression
is not: **0.80 µs/pair is the floor of the current architecture**, and `-p4` makes it
*worse* because pigz threads only the CRC and read-ahead, not the inflate.

So the target is: single-threaded C++ at ~2.7 µs/pair (**3.0x** over p1, on one core
instead of eight), and with a GIL-releasing chunk kernel plus 2–3 threads, compute drops
under the decompression floor and the tool lands at **~0.9–1.1 µs/pair — about 8x over
p1 and 2x over today's p8, on far fewer cores.** Beyond that, gzip is the wall, and the
answer to gzip is to stop writing the intermediate at all (§11).

---

## 3. The prototype, and what it already proves

A complete single-threaded C++ prototype (parse → scan → consensus → build → format →
write) is committed at
[`scripts/merge_bench/proto_merge.cpp`](../scripts/merge_bench/proto_merge.cpp), with
every benchmark behind the numbers in this document beside it and a
[README](../scripts/merge_bench/README.md) that reproduces them. On the 200,000-pair
library it emits
**a byte-identical file to `zna merge`** — 177,296 merged, 5,391 trimmed, 17,313 kept,
222,704 records, `cmp` clean — at **2.79 µs/pair against 8.34**.

That single equality is the strongest evidence in this document, because it
simultaneously validates:

- fixed-point integer scoring ≡ the float score (§4),
- the argmax tie-break (§5),
- packed 2-bit popcount comparison ≡ byte comparison (§6),
- the SWAR packer and its complement-by-XOR (§6),
- the posterior consensus table,
- `_build_merged` from the fragment span, the merged-name construction, the
  merge/trim/keep decisions and the minimum-length filtering.

It is a prototype, not the design: it is single-threaded, its `popen("pigz -dc")` reader
is a stand-in, and several per-pair buffers are copied naively (§7.4). The design below
is what to build; the prototype is the evidence that it lands where it claims.

---

## 4. Fixed-point scoring, and why `SCALE = 2^24`

The roadmap asks for fixed-point so the argmax is exactly reproducible across compilers.
It does not say what scale, and the scale is a correctness decision, not a taste one.

Define, once, in Python:

```
SCALE   = 1 << 24
M       = round(log2((1-e)/0.25) * SCALE)      # 33_311_170 at e = 0.01
MM      = round(log2(0.75/e)     * SCALE)      # 137_813_406 ... (STEP = M + MM)
STEP    = M + MM
T_merge = round(t_merge * SCALE)
T_trim  = round(t_trim  * SCALE)
```

Score is `n*M - d*STEP` in `int64`. **No float appears anywhere in the kernel**, so the
argmax is bit-identical across compilers, optimisation levels and `-ffast-math`.

**Does quantisation ever flip a decision against today's float behaviour?** This is
exhaustively enumerable: a flip needs the true score within the quantisation error of a
threshold, over integer `(n, d)`. Enumerated for `n ≤ 4000, d ≤ 1200`:

| scale | flips at T=8 | flips at T=28 |
|---|---:|---:|
| 2^16 | 2 | 1 |
| 2^20 | 0 | **1** (at n=2575, d=619) |
| **2^24** | **0** | **0** |
| 2^30 | 0 | 0 |

The closest any reachable `(n, d)` gets to a threshold is 7.2×10⁻⁵ bits (T=28, at
n=2575 d=619) and 3.6×10⁻⁴ bits (T=8, at n=1994 d=481). **2^24 is provably
decision-identical to float over the whole enumerated domain**, and 2^20 — the obvious
"millionths of a bit" choice, and what the prototype first used — is not.

Two consequences:

- **Add `--max-read-length` (default 4000, hard error above).** The audit already wanted
  it as a fail-fast for the O(L²) scan. Here it does more: it makes the enumerated
  domain the *reachable* domain, turning the table above from an empirical range check
  into a guarantee. `n*M ≤ 4000 × 2^25 ≈ 2^37`, nowhere near `int64` overflow.
- **Derive `M`/`STEP` in Python only, and pass the integers across.** `log2` is not
  correctly rounded and libm differs between platforms; computing them on both sides
  invites a 1-ULP disagreement that changes a corpus. One call site, integers over the
  boundary, and the values recorded in the JSON stats so any corpus can be audited.

The same argument applies to `_DISAGREE_Q`, which is built with `pow`/`log10`: **build
it in Python and pass the 94×94 byte table into the backend.** The prototype builds it
in C++ and happens to agree here; that is luck, not a guarantee. 8,836 bytes, once per
run.

---

## 5. The argmax is a total order, and it is specifiable

The roadmap says a differential test "must be order-independent — assert the returned
shift is *an* argmax with matching n/d/score", because ties occur on 0.845% of real
pairs and an order-replicating reference would break on the very rewrite it guards.

That advice protects against the wrong hazard, and following it would give up something
valuable. The shipped scan visits shifts in a defined order — the plateau first
(descending `n`, ascending `s`), then both flanks at decreasing `n`, read-through side
first — and improves only on strict `>`. Combined, that is exactly the total order

> **maximise `score`; among ties maximise `n`; among those minimise `s`.**

`n` is a function of `s`, so this is a total order on shifts, computable by any
implementation in any evaluation order. Verified: an exhaustive unpruned scan compared
against `find_overlap` over 4,000 adversarial synthetic pairs, 3,000 real pairs, and —
the part that matters — **144 constructed tie cases, up to 13-way, zero violations.**
(Ties need two shifts with equal `n` *and* equal `d`; since `M/STEP` is irrational-ish
this essentially only happens inside the plateau, i.e. when `len1 ≠ len2`.)

So: **specify the tie-break, and make the differential tests byte-exact.** The
vectorisation opportunity is *within* a shift (popcount over packed bases), not *across*
shifts, so the comparator stays in the scalar outer loop and costs nothing — the
roadmap's SIMD concern does not arise. What we gain is that a given FASTQ produces a
byte-identical output on every platform, which is the actual requirement for training
data, and that the test suite can `cmp` whole files instead of asserting a weaker
property.

Pruning preserves the order for free: shifts are visited in decreasing `n`, and a shift
that can only *tie* loses the tie-break to the incumbent anyway, so rejecting on
`ceiling <= best` is correct. In integers the bail bound becomes exact —
`dmax = (ceiling - best - 1) / STEP` — which also removes the float truncation in
`dmax = int((ceiling - best) / step)` that the roadmap flagged as the one place a float
decides control flow.

---

## 6. The kernel

### 6.1 What was measured

Full pruned scans — plateau, both flanks, ceiling break, per-shift bail — over 50,000
real pairs, **with packing inside the timed region**, each variant checked against the
shipped kernel's `(s, n, d)` on every pair:

| variant | µs/pair | vs numba | vs scalar C++ |
|---|---:|---:|---:|
| numba, as shipped | 2.633 | 1.00x | — |
| scalar C++, direct transliteration | 1.090 | 2.42x | 1.00x |
| byte-SWAR, 8 bases/word, exact byte semantics | 1.040 | 2.53x | 1.05x |
| 2-bit packed, **table** packer, bail 32 | 1.027 | 2.56x | 1.06x |
| 2-bit packed, SWAR packer, bail 32 | 0.631 | 4.17x | 1.73x |
| 2-bit packed, SWAR packer, bail 64 | 0.608 | 4.33x | 1.79x |
| 2-bit packed, SWAR packer, bail 128 | 0.682 | 3.86x | 1.60x |
| **2-bit packed, SWAR packer, fused flank pair** | **0.563** | **4.68x** | **1.94x** |

Zero mismatches in every row. Four findings drive the design:

1. **Merely being C++ is worth 2.42x** — no numba dispatch (the audit measured a null
   njit kernel at 0.221 µs/pair of pure dispatch, a ceiling on what any numba-side
   change could return), and better codegen.
2. **Packing with a per-base table lookup destroys the entire win** (1.027 vs 0.631).
   The SWAR packer is not an optimisation, it is a precondition.
3. **Bail granularity 64 is optimal; 128 is worse than 64.** This **kills the AVX2/NEON
   idea for this workload**: wider vectors mean a coarser bail, and the scan is
   rejection-dominated (~0.32·n comparisons per rejected shift). Portable 64-bit scalar
   code *is* the optimum — no intrinsics, no runtime ISA dispatch, no `-march` flags.
   That is a large simplification bought with a measurement.
4. **Fusing the two flank shifts at each `n`** gives another 1.08x by overlapping two
   serial bail chains.

Net: **4.68x on the scan** against the shipped numba kernel.

### 6.2 How it works

**Packing.** Bases go 2 bits each, base `i` at bits `2*(i%32)` of word `i/32`, so
advancing `k` bases is a right shift by `2k`. For equality comparison the code map only
has to be a bijection on ACGT, and `(c >> 1) & 3` already is one (A→0, C→1, G→3, T→2) —
no table, so eight bases pack in one SWAR step:

```c
uint64_t p = (x >> 1) & 0x0303030303030303;   // 2 bits per byte
p = (p | (p >> 6))  & 0x000F000F000F000F;     // gather to nibbles
p = (p | (p >> 12)) & 0x000000FF000000FF;
p = (p | (p >> 24)) & 0x000000000000FFFF;
```

Under that map, complementing is **XOR 2** (A0↔T2, C1↔G3), so `revcomp(R2)` packs
directly from `R2` with no intermediate buffer — which also deletes today's
`reverse_complement(s2)` allocation.

**Comparison.** `x = a ^ b` is nonzero in a base's 2-bit field iff the bases differ;
`(x | (x >> 1)) & 0x5555…` folds that to one bit per base, and `popcount` counts them.

**Alignment.** For shift `s`, the R1 offset is `max(s,0)` and the R2rc offset is
`max(-s,0)` — so **exactly one of the two is always zero**, hence word-aligned, and only
one side ever needs the double-word extraction. That halves the realignment cost and is
a property of the scan's geometry, not a trick.

### 6.3 N and IUPAC: two paths, one semantics

Today the kernel compares raw bytes, so `N` vs `N` earns a full match and IUPAC codes
compare as themselves (`reverse_complement` passes them through uncomplemented —
deliberate, trap #8). 2-bit packing cannot represent that: `N` would collapse onto a
real base and start manufacturing evidence. The audit lists exactly this as a reason the
numba SWAR prototype was rejected — "an N-semantics change in a corpus-quality tool".

**Resolution: dispatch, do not approximate.** The SWAR packer computes a purity flag as
a by-product (uppercase, then a SWAR test for one of the four bytes). If both mates are
pure ACGT — 99.83% of real pairs, and 100% of the benchmark library — take the packed
path. Otherwise fall back to the byte-comparison scan, which *is* the reference
algorithm with today's exact semantics.

This is not a compromise: on pure-ACGT input the two paths are provably equivalent
(byte equality ≡ 2-bit code equality under a bijection), and on everything else the
semantics are unchanged, bit for bit. The cost is one branch per pair.

---

## 7. Architecture

### 7.1 Backend split, mirroring the codec

zna already has the pattern this needs, and it should be copied rather than reinvented:
`zna/codec.py` selects between `zna._pycodec` and `zna._accel` by name, with
`get_backend()`, `available_backends()`, and a validated required-function set.

```
src/zna/merge/
    __init__.py     lazy exports (unchanged)
    args.py         argparse only (unchanged)
    params.py   NEW MergeParams + the fixed-point derivation (§4). Pure Python,
                    no backend import. One source of truth for M/STEP/T_*/DISAGREE_Q.
    backend.py  NEW get_merge_backend(name=None) -> module   [mirrors codec.py]
    _pymerge.py NEW the reference backend, pure Python, no numba
    _accel.cpp  NEW the C++ backend  -> zna.merge._accel
    overlap.py      thin public API: find_overlap / score_weights / threshold_bits,
                    dispatching to the selected backend. Keeps the `njit` no-op shim,
                    which is public surface the tests import.
    pairs.py        thin public API: MergeParams, PairOutcome, process_pair, base_name
    fastqio.py      unchanged (stream I/O stays in Python)
    cli.py          orchestration
```

A **second CMake target**, installed to `zna/merge/`, rather than adding to
`_accel.cpp`: merging is a different concern from the codec, a merge build failure must
not break the format library, and `import zna` must keep costing nothing (the port
already established that `zna/__init__.py` does not mention this package).

### 7.2 Three levels, one implementation each side

Both backends expose the same three entry points. The fine-grained ones exist so that
failures localise; the coarse one is what production calls.

```python
# L1 — the kernel. The oracle. Per pair, no allocation, no Python objects.
scan(s1: bytes, s2rc: bytes, w: Weights) -> (shift, score_q, olen, diff)

# L2 — one pair: decisions, consensus, record construction.
process_pair(h1, s1, q1, h2, s2, q2, p: ParamsQ) -> (records, outcome, n_dropped,
                                                     score_q, olen, diff)

# L3 — the production path: a slab of raw FASTQ text in, a formatted blob out.
merge_chunk(buf1: bytes, buf2: bytes, p: ParamsQ, check_sync: bool)
    -> (blob: bytes, consumed1: int, consumed2: int, counters: Counters)
```

`ParamsQ` carries only `int64`s and the consensus table, so **no float ever crosses the
boundary**. L3 is a pure function of its arguments: no RNG, no clock, no globals — which
is what makes it fuzzable, thread-safe and byte-deterministic.

### 7.3 Why L3 is the boundary

L1 alone would deliver ~1.4x (§1). L3 captures the parse, the scan, the consensus, the
construction, the formatting and the histograms — everything in the 8.34 µs except the
gzip and the write. It is also the level at which the GIL can be released for a useful
span of work, which is what makes threads viable.

**Chunking protocol.** Python appends raw reads to a per-stream leftover buffer and
calls `merge_chunk`; C++ consumes whole records only and reports `consumed1`/`consumed2`
separately, so the two streams may carry different leftovers and nothing needs to scan
for record boundaries in Python. Read into each buffer only while it holds fewer than
TARGET bytes, so memory is bounded by `(TARGET + CHUNK) × 2 × threads`.

The roadmap deferred exactly this ("raw-blob IPC") because the audit's prototype had
four defects, one disqualifying: *R1 shorter than R2 returned OK and silently dropped
the trailing R2 records*, which the desync check cannot catch because comparing base
names cannot see records that were never read. The protocol above makes that structural:
`merge_chunk` consumes `min(records1, records2)` pairs, and at EOF the caller asserts
**both** buffers are empty. A non-empty leftover on either side is
`InputError("R1 exhausted before R2")` / vice versa. The other three defects — truncated
record → `IndexError`, R2 short → `IndexError`, CRLF surviving into the sequence — are
parser cases with named tests (§9).

**Errors.** The C++ parser throws a distinct type that a nanobind exception translator
raises as `zna.merge.fastqio.InputError`, so `except InputError` in `cli.py` and the
existing tests keep working unchanged, with the same messages.

### 7.4 What the prototype does naively, and what the real one must not

The prototype's 1.43 µs/pair of compute is ~0.56 scan + ~0.87 everything else, and the
"everything else" is dominated by avoidable copies: it unconditionally copies `s1`/`q1`
into mutable buffers (needed only for the 43.5% of pairs the consensus actually
changes), builds `s2rc` byte-at-a-time, reverses `q2` byte-at-a-time, and appends the
merged quality one `push_back` at a time. Copy-on-write for the consensus, a SWAR
revcomp, and `memcpy`-based record assembly should take compute meaningfully below
1 µs/pair. **That is an estimate, not a measurement** — it is the first thing to check
once L3 exists, and the design does not depend on it.

---

## 8. Threading, and deterministic output

With L3 releasing the GIL (`nb::gil_scoped_release` around everything after the input
pointers are taken; `bytes` are immutable and kept alive by the argument reference), the
worker pool becomes a `ThreadPoolExecutor`. That deletes, outright:

- the `fork` context and its Windows fallback,
- `_init_worker` and the `_W` worker-global dict,
- pickling chunks out and blobs back,
- `BrokenProcessPool` handling and the killed-worker test,
- the **sparse** histograms, which exist only because dense ones were being pickled per
  chunk (the audit measured p8 going 4.6 → 10.1 µs/pair). With threads, dense arrays are
  free and simpler.

Reading (`stream.read`) and writing both release the GIL, so a reader-in-main-thread,
N-workers, write-on-completion loop keeps every core busy. Since compute (~1 µs) sits
near the decompression floor (0.80 µs), **`--threads 2–4` is expected to saturate and
beyond that nothing helps** — worth stating in `--help` rather than letting users set 32.

**Write blobs in submission order.** Chunks are independent, so ordering costs only
slight head-of-line blocking, and it buys something today's tool does not have: `zna
merge` output becomes **byte-deterministic regardless of thread count**. That upgrades
contract C6 from "pairs stay adjacent" to "the file is a pure function of the input and
the parameters", lets the test suite `cmp` whole files across `--threads 1/2/8`, and
makes any future corpus defect bisectable.

`--processes` stays as a **deprecated alias for `--threads`** so the pipeline's existing
rule keeps working, with a one-line deprecation notice and a note in the JSON stats.

---

## 9. Testing

The port's own instruments stay and are the acceptance gate; two of them were
re-verified against deliberate mutants during the port and must be re-verified again
after each phase here.

1. **The oracle.** `_pymerge` is the reference implementation, not a fallback. It is
   never deleted and never optimised at the cost of clarity.
2. **Byte-exact differential, at all three levels.** L1: `(shift, score_q, olen, diff)`
   equal on every pair. L2: identical record tuples. L3: identical blob **and** identical
   counters. Because the tie-break is specified (§5), these are equality assertions, not
   "an argmax" assertions.
3. **Whole-file equality** between backends, and across `--threads 1/2/8`, on the
   200k-pair library — the check that caught nothing subtle in the prototype only
   because the prototype was already right.
4. **`test_block_loop_sees_every_position` first.** The roadmap is right that this is the
   test to write before touching the loop: `k += 7` in the unrolled stride changes 6.34%
   of scores and 0.26% of decisions on real data and passes every other test in the
   suite. It is currently the only test that kills that mutant. The packed loop needs its
   own equivalent — a single mismatch swept across all 32 phases of a word *and* across
   the word boundary, plus the partial-word mask.
5. **Fuzz, in the shape of `tests/test_fuzz_roundtrip.py`.** Random lengths (including 0,
   1, and either mate longer), random content including N and IUPAC, random parameters,
   asserting backend equality. That file exists because two 0.3.4 optimisations silently
   corrupted data while 282 tests passed; the same risk applies here exactly.
6. **Independently-drawn mate lengths, everywhere.** The port's non-negotiable. Any new
   geometry fixture that sets `len1 == len2` is structurally blind to the C2 truncation
   class.
7. **ASAN/UBSAN build in CI.** The packed reader deliberately reads one guard word past
   the data (`getw` touches `P[idx+1]`); an off-by-one there is a heap overread that
   produces plausible output. This is the single most dangerous new code in the design
   and it needs a sanitiser, not a code review.
8. **Mutation checks re-run after each phase**: the C2 regate (`direction == FORWARD &&
   shift > 0`) must fail ≥10 tests including both e2e boundary parametrisations, and
   `k += 7` must fail the block-loop test.

---

## 10. Portability

| item | risk | resolution |
|---|---|---|
| `popcount` on the x86-64 baseline | manylinux2014, no `-march`; `__builtin_popcountll` becomes a libgcc call | 2-way runtime dispatch on `__builtin_cpu_supports("popcnt")` (SSE4.2-era, universal in practice) with a SWAR fallback. **Measure both** — if the fallback is within noise, drop the dispatch. |
| endianness | the SWAR packer `memcpy`s 8 bytes into a `uint64_t`; a big-endian build would map bases to different bit positions than the byte-at-a-time paths | `static_assert` little-endian. Every supported target (x86-64, aarch64, Windows) is LE; `_accel.cpp` already assumes it. |
| unaligned loads | `memcpy` into `uint64_t` | already the idiom used in `_accel.cpp`; compilers emit unaligned loads |
| MSVC | no `__builtin_popcountll` | `__popcnt64` / `<intrin.h>` behind one inline |
| C++ standard | `std::popcount` is C++20; the project is `cxx_std_17` | one portable inline; do not raise the toolchain floor |
| numba | disappears | §12 |

---

## 11. Phasing

Every phase ends green, with byte-identical output to the phase before it, so a
regression bisects to one change.

- **A — fixed point, in Python.** Introduce `params.py`, `SCALE = 2^24`, integer
  scoring, the exact `dmax`, `--max-read-length`, and the tie-break spec written down
  and tested. Still numba. *Acceptance: byte-identical output on the 200k library, and
  the §4 enumeration as a test.* This isolates "did fixed point change anything?" from
  "did C++ change anything?".
- **B — backend seam.** Add `backend.py` and `_pymerge.py`; make `overlap.py`/`pairs.py`
  dispatch. Numba still available behind the Python backend. *Acceptance: suite green on
  both backends.*
- **C — C++ L1.** The scan only, with the packed kernel and the byte fallback.
  *Acceptance: L1 differential on the 200k library + fuzz; block-loop test for the packed
  loop; ASAN clean.* Expected: scan 2.63 → 0.56 µs/pair.
- **D — C++ L2.** Consensus, construction, naming. *Acceptance: identical record tuples.*
  Expected: p1 8.34 → ~4 µs/pair.
- **E — C++ L3 + threads + ordered output.** Parser, formatting, histograms, GIL
  release, `ThreadPoolExecutor`, `--processes` → alias. *Acceptance: whole-file equality
  across backends and thread counts; the four IPC defects each have a named test.*
  Expected: ~2.7 µs/pair at 1 thread, ~0.9–1.1 at 2–4.
- **F — remove numba.** Drop the dependency, the `merge` extra, and the numba-specific
  refusal; `zna merge` now refuses without the *accel backend* unless `--allow-slow`,
  exactly parallel to `is_accelerated()` for the codec. Keep the `njit` no-op shim, which
  the tests import.
- **G — optional, measure first.** Replace the pigz subprocesses with libdeflate or
  ISA-L inflate inside the extension. Only worth it once §2's 0.80 µs/pair floor is
  actually the binding constraint — i.e. after E, and probably not before phase 2 below
  removes the intermediate anyway.

**Then** `zna encode --merge-pairs` (the roadmap's phase 2), which feeds the merged
records straight into the existing `(seq, is_paired, is_read1, is_read2, has_start,
has_end)` ingest path and deletes the FASTQ intermediate, its pigz write, its gzip read
and its full re-parse. That is where the remaining I/O goes, and it is worth more than
anything left in the kernel.

---

## 12. Release

**Do not ship 0.3.6 with the numba extra.** The user's instinct here is right and worth
making explicit: phase F removes numba, the `merge` extra and `--allow-slow`'s current
meaning. Publishing them and withdrawing them one release later is user-visible
dependency churn — a `pip install zna[merge]` and a conda `run_constrained` that stop
meaning anything — bought for nothing. The port is committed and green; hold it.

**Release the whole thing as 0.4.0**, not 0.3.6. It adds a subcommand, adds a second
compiled extension, changes the parallelism model from processes to threads, changes
`--processes` to a deprecated alias, and removes an optional dependency. On a 0.x
project that is minor-version material.

If the C++ work slips and something must ship first, ship the port as 0.3.6 *without*
the numba extra by making the pure-Python path the only Python path and keeping the
refusal — i.e. land phase B early. That is a smaller commitment than publishing an
extra.

---

## 13. Deliberately not doing

Each with the number that closes it. [MERGE_TOOL_AUDIT.md](MERGE_TOOL_AUDIT.md) §(c) has
the rest; check it before proposing an improvement.

| idea | why not |
|---|---|
| **AVX2/NEON inner loop** | Measured here: bail at 128 bases is **0.682 µs/pair against 0.608 at 64**. The scan is rejection-dominated, so wider vectors buy throughput on work that should not have happened. Portable 64-bit scalar is the optimum. |
| **byte-SWAR instead of packing** | 1.040 vs 1.090 µs/pair — **1.05x**. Same bail granularity, and the compare was never the cost. |
| **table-driven packing** | 1.027 vs 0.631 µs/pair. Costs ~0.40 µs/pair and gives back the whole win. |
| **k-mer seed-and-extend to propose shifts** | Audit: implemented and measured, **no gain** over the block loop (1.58x vs 1.57x). |
| **Gating merges on shift ambiguity** | Audit §(c): a microsatellite detector, not an identifiability test. Would re-emit **+4.68% duplicated bases**. |
| **Composition-aware `p_null`** | Audit: measured negative three times; −0.55pp merge rate, and ~2pp on exactly the low-complexity reads it was meant to protect. |
| **Estimating `err_rate` from the data** | Audit: `ê = 0.00866` vs 0.01 assumed; adopting it costs −0.04pp merge rate. |
| **Reproducing the reference's iteration order in tests** | Superseded: §5 shows the order *is* a specifiable total order, so tests can be byte-exact instead of weaker. |

---

## 14. Open questions

1. **Does the popcount runtime dispatch earn its complexity on x86-64 baseline?** Needs
   a Linux/x86 measurement; this session's numbers are all aarch64, where `popcount` is
   a native instruction. If the SWAR fallback is within noise, delete the dispatch.
2. **How much of §7.4's ~0.87 µs/pair of non-scan compute survives a careful
   implementation?** Estimated well under 1 µs total; unmeasured.
3. **Does ordered output cost measurable throughput** at 4 threads, or is head-of-line
   blocking lost in the noise? If it costs more than ~2%, make it `--ordered` (default
   on) rather than unconditional.
4. **Is `--max-read-length 4000` the right default?** It must exceed any real read and
   stay inside the enumerated domain of §4. 4000 satisfies both; a long-read chemistry
   would need the enumeration re-run, not just the flag raised.
