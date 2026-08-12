# Native acceleration for `zna merge`: design

Status: **proposed**, 2026-08-12, against zna at the read-merge port.
Target release: **0.4.0**, numba-free — see §12.

This is the design for replacing the numba overlap kernel with a C++ backend, keeping a
pure-Python implementation as the **reference oracle** rather than a fallback, the way
`_pycodec` and `_accel` already relate.

> **Revision 2, 2026-08-12 — after external review.** The review's central technical
> challenge was correct and **the kernel design changed because of it**: byte-wise
> vector comparison beats 2-bit packing (§6), which also deletes the packed path's
> N/IUPAC dispatch. Revision 1's claim that "AVX2/NEON is the wrong answer" was an
> over-generalisation from a sound measurement, and is retracted in §6.1. The tie-break
> rule is simplified and now proved rather than asserted (§5). §15 records what else was
> adopted, and which review items were factually wrong for this codebase.

Every number below was measured on this repo, on 200,000 simulated 2×150 pairs with
insert ~ N(200, 70) truncated to [50, 400] and 0.5% per-base error. That library merges
at **88.6%** against production's measured 88.8%, and its overlap mismatch rate is
0.0075 against production's 0.0091, so it is representative. Machine: Apple M-series
(aarch64/NEON), Python 3.12, numba 0.66, Apple clang 21, `-O3 -std=c++17`, min-of-N.
**No x86 measurements exist yet** — §10, §14.

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
- **remove the FASTQ intermediate entirely** (`--merge-pairs`), which is where the
  compounding win actually lives.

Anything that only speeds the kernel is measured against a tool that is already I/O
bound at production settings. Say so up front rather than discovering it after.

---

## 2. The floor

| configuration | µs/pair |
|---|---:|
| gzip inflate floor — 2 × `pigz -dc -p1`, concurrent, to `/dev/null` | **0.80** |
| same with `-p4` each (pigz cannot parallelise *inflate*) | 0.95 |
| C++ prototype: gz in → file out, single thread | 2.67 |
| C++ prototype: gz in → `/dev/null` | 2.69 |
| C++ prototype: **plain** in → `/dev/null` (compute only) | **1.43** |
| `zna merge -p1` (numba) | 8.34 |
| `zna merge -p8` (numba) | 1.97 |

The writer is free (pigz level 1 in a subprocess keeps up). Decompression is not:
**0.80 µs/pair is the floor of the current architecture**, and `-p4` makes it *worse*
because pigz threads only the CRC and read-ahead, not the inflate.

So the target is: single-threaded C++ at ~2.7 µs/pair (**3.1x** over p1, on one core
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
library it emits **a byte-identical file to `zna merge`** — 177,296 merged, 5,391
trimmed, 17,313 kept, 222,704 records, `cmp` clean — at **2.67 µs/pair against 8.34**.

That single equality simultaneously validates:

- fixed-point integer scoring ≡ the float score (§4),
- the argmax tie-break (§5),
- the vectorised comparison kernel (§6),
- the posterior consensus table,
- `_build_merged` from the fragment span, the merged-name construction, the
  merge/trim/keep decisions and the minimum-length filtering.

It was re-verified byte-identical after the kernel changed from 2-bit packed to
byte-wise SIMD, which is the check that matters when swapping a kernel.

It is a prototype, not the design: it is single-threaded, its `popen("pigz -dc")` reader
is a stand-in, it builds the consensus table in C++ where §4 says to build it once in
Python, and several per-pair buffers are copied naively (§7.4).

---

## 4. Fixed-point scoring, and why `SCALE = 2^24`

Define, once, in Python:

```
SCALE   = 1 << 24
M       = round(log2((1-e)/0.25) * SCALE)      #  33_311_170 at e = 0.01
MM      = round(log2(0.75/e)     * SCALE)      # 104_502_237
STEP    = M + MM                               # 137_813_407
T_merge = round(t_merge * SCALE)
T_trim  = round(t_trim  * SCALE)
```

Score is `n*M - d*STEP` in **`int64`**, and the per-shift bail bound is the exact
integer `dmax = (ceiling - best - 1) / STEP`. **No float appears anywhere in the
kernel**, so the argmax is bit-identical across compilers, optimisation levels and
`-ffast-math`. That also removes the float truncation in `dmax = int((ceiling - best) /
step)` that the roadmap flagged as the one place a float decides control flow.

`int64` is required, not stylistic: `M ≈ 2^25`, so `n*M` exceeds `int32` at **n = 64**.
A 150 bp overlap needs 33 bits.

There are **no lookup tables** in this scoring rule — the score is two integer
multiplies, because the weights are per-match and per-mismatch constants, not
per-position quality-dependent values. (A per-position quality-weighted score, which
*would* want a `(Q1,Q2)` LUT, is redesign §7 and was measured and rejected three times;
see [MERGE_TOOL_AUDIT.md](MERGE_TOOL_AUDIT.md) §(c).)

**Does quantisation ever flip a decision against today's float behaviour?** This is
exhaustively enumerable: a flip needs the true score within the quantisation error of a
threshold, over integer `(n, d)`. Enumerated for `n ≤ 4000, d ≤ 1200`:

| scale | flips at T=8 | flips at T=28 |
|---|---:|---:|
| 2^16 | 2 | 1 |
| 2^20 | 0 | **1** (at n=2575, d=619) |
| **2^24** | **0** | **0** |
| 2^30 | 0 | 0 |

The closest any reachable `(n, d)` gets to a threshold is 7.2×10⁻⁵ bits (T=28) and
3.6×10⁻⁴ bits (T=8). **2^24 is provably decision-identical to float over the whole
enumerated domain**, and 2^20 — the obvious "millionths of a bit" choice, and what the
prototype first used — is not.

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

## 5. The argmax is a total order — maximise score, then minimise `s`

The roadmap says a differential test "must be order-independent — assert the returned
shift is *an* argmax with matching n/d/score", because ties occur on 0.845% of real
pairs and an order-replicating reference would break on the very rewrite it guards.

That advice protects against the wrong hazard, and following it would give up something
valuable. The canonical rule is:

> **maximise `score`; among ties, minimise `s`** (the leftmost / most read-through
> alignment).

Two facts establish it, one measured and one arithmetic:

- **Measured.** An exhaustive unpruned scan compared against `find_overlap` over 4,000
  adversarial synthetic pairs, 3,000 real pairs, and — the part that matters — **144
  constructed tie cases, up to 13-way, zero violations.** Random sequence essentially
  never ties, so the ties had to be built deliberately; they arise inside the plateau,
  i.e. when `len1 ≠ len2`.
- **Arithmetic.** A tie between shifts with *different* overlap lengths needs
  `Δn·M = Δd·STEP`, whose minimal solution is `Δn = STEP/gcd(M,STEP)`. For the shipped
  weights `gcd(M, STEP) = 1`, so `Δn ≥ 137,813,407` — unreachable under any read length,
  and the same holds for `e` ∈ {0.001, 0.005, 0.01, 0.02, 0.05}. **Ties can only occur
  at equal `n`**, so `n` never discriminates and the rule needs only the two keys above.

The shipped scan already implements exactly this: shifts are visited in decreasing `n`
and ascending `s` (plateau first, then each flank pair read-through side first, which is
always the smaller `s`), and improvement is strict `>`. Pruning preserves it for free —
a shift that can only *tie* loses to the incumbent, so rejecting on `ceiling <= best` is
correct.

So: **specify the tie-break, and make the differential tests byte-exact.** The
vectorisation opportunity is *within* a shift, not *across* shifts, so the comparator
stays in the scalar outer loop and costs nothing — the roadmap's SIMD concern does not
arise. What we gain is that a given FASTQ produces a byte-identical output on every
platform, which is the actual requirement for training data, and that the test suite can
`cmp` whole files instead of asserting a weaker property.

---

## 6. The kernel: byte-wise vector comparison

### 6.1 What was measured, and a retraction

Full pruned scans — plateau, both flanks, ceiling break, per-shift bail — over 50,000
real pairs, **with packing inside the timed region** where applicable, every variant
checked against the shipped kernel's `(s, n, d)` on every pair. Zero mismatches
throughout.

| variant | µs/pair | vs numba | vs scalar C++ |
|---|---:|---:|---:|
| numba, as shipped | 2.633 | 1.00x | — |
| scalar C++, direct transliteration | 1.075 | 2.45x | 1.00x |
| byte-SWAR, 8 bases per 64-bit word | 1.040 | 2.53x | 1.03x |
| 2-bit packed, **table** packer | 1.027 | 2.56x | 1.05x |
| 2-bit packed, SWAR packer, bail 64 | 0.608 | 4.33x | 1.77x |
| 2-bit packed, SWAR packer, fused flanks | 0.535 | 4.92x | 2.01x |
| byte SIMD, 16 B vectors, bail every 16 | 0.534 | 4.93x | 2.01x |
| **byte SIMD, 16 B vectors, bail every 32** | **0.470** | **5.60x** | **2.29x** |
| byte SIMD, 16 B vectors, bail every 64 | 0.515 | 5.11x | 2.09x |
| byte SIMD, 16 B, bail 32, fused flanks | 0.555 | 4.74x | 1.94x |

**Retraction.** Revision 1 of this document concluded that "AVX2/NEON is the wrong
answer for this workload", on the grounds that widening the *packed* kernel's bail
granularity from 64 to 128 bases made it slower. That measurement was sound; the
inference from it was not. Bail granularity within the packed scheme is not the same
variable as vector width in a byte-wise scheme, and the byte-wise scheme has a structural
advantage the packed one cannot have: **a shift is a pointer offset, so there is no
cross-word bit realignment at all.** Measured, it wins by 1.14x over the best packed
variant. The external reviewer who pressed this point was right.

Four findings drive the design:

1. **Merely being C++ is worth 2.45x** — no numba dispatch (the audit measured a null
   njit kernel at 0.221 µs/pair of pure dispatch, a ceiling on what any numba-side
   change could return), and better codegen.
2. **Byte-wise vector comparison beats 2-bit packing** (0.470 vs 0.535), needs no
   packing step at all, and — see §6.3 — is semantically exact on N and IUPAC for free.
3. **Bail granularity 32 bases (two 16-byte vectors) is the optimum**, both narrower (16)
   and wider (64) are worse. The scan is rejection-dominated — roughly 0.32·n
   comparisons per rejected shift — so the bail interval, not raw throughput, is the
   variable that matters.
4. **Fusing the two flank shifts hurts the byte kernel** (0.555 vs 0.470) although it
   helped the packed one, because it doubles register pressure and defeats the
   early-exit ordering. Do not carry that idea across.

Net: **5.6x on the scan** against the shipped numba kernel, with the simplest of the
variants tried.

### 6.2 How it works

```c
static inline int neq16(const uint8_t* a, const uint8_t* b) {   // mismatches in 16 bytes
    uint8x16_t eq = vceqq_u8(vld1q_u8(a), vld1q_u8(b));         // NEON
    return 16 - (int)vaddvq_u8(vandq_u8(eq, vdupq_n_u8(1)));
}
// SSE2:  16 - popcount(_mm_movemask_epi8(_mm_cmpeq_epi8(va, vb)))
```

and the per-shift loop is:

```c
const uint8_t* a = s1   + (s > 0 ?  s : 0);      // a shift is just a pointer offset
const uint8_t* b = s2rc + (s < 0 ? -s : 0);
for (k = 0; k + 32 <= n; k += 32) {
    d += neq16(a + k, b + k) + neq16(a + k + 16, b + k + 16);
    if (d > dmax) return REJECT;                  // bail every 32 bases
}
for (; k < n; ++k) d += (a[k] != b[k]);           // scalar tail
```

Unaligned 16-byte loads are penalty-free on every target of interest, and the loop
never reads past `n` because the vector body requires `k + 32 <= n`. **No packing, no
realignment, no guard word, no lookup table, no auxiliary validity mask.**

Note the scalar tail: an overlap in the trim band is 5–14 bases, so short overlaps are
handled entirely scalar. A 16-byte single-vector step for `16 ≤ rem < 32` is a tuning
detail to measure, not a design decision.

### 6.3 N and IUPAC: nothing to do

Today's kernel compares raw bytes, so `N` vs `N` earns a full match and IUPAC codes
compare as themselves (`reverse_complement` passes them through uncomplemented —
deliberate, trap #8). 2-bit packing cannot represent that: `N` would collapse onto a real
base and start manufacturing evidence. Revision 1 handled this with a purity flag and a
dual path, which the audit's own critique of the numba SWAR prototype demanded.

**Byte-wise comparison removes the problem instead of managing it.** `a[k] != b[k]` on
raw bytes *is* the reference semantics — for ACGT, for N, for IUPAC, for lowercase, for
any byte at all. There is no dispatch, no purity test, no second code path to keep in
sync, and no class of input on which the fast path and the reference disagree. This is
the strongest argument for the byte kernel and it is worth more than the 1.14x.

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
    _pymerge.py NEW the reference backend, pure Python, no numba, integer scoring
    _accel.cpp  NEW the C++ backend  -> zna.merge._accel  (nanobind, as zna already uses)
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

### 7.2 The C++ core is I/O-agnostic; adapters sit above it

This is the review's best structural point, and it is what makes `--merge-pairs` cheap
later. The C++ module is two strict layers:

**Core (no I/O, no Python objects, spans only).** Operates on `const uint8_t*` + length
and writes into caller-provided scratch:

```cpp
struct ScanResult { int32_t shift; int32_t olen; int32_t diff; int64_t score_q; };
ScanResult scan(const uint8_t* s1, int len1, const uint8_t* s2rc, int len2,
                const Weights& w);

struct PairOut { ... spans into scratch ...; Outcome outcome; int n_dropped; };
PairOut process_pair(Read r1, Read r2, const ParamsQ& p, Scratch& scratch);
```

**Adapters.** Two of them, and later a third:

- `merge_chunk(buf1, buf2, …) -> (blob, consumed1, consumed2, counters)` — FASTQ text in,
  FASTQ text out. This is `zna merge`.
- the same core, emitting `(seq, is_paired, is_read1, is_read2, has_start, has_end)`
  records instead of formatting FASTQ — this is `zna encode --merge-pairs`, and it
  reuses the core untouched.

Writing the layering down now costs nothing and prevents the core from growing FASTQ
assumptions that would have to be unpicked later.

### 7.3 Three levels across the Python boundary

Both backends expose the same three entry points. The fine-grained ones exist so that
failures localise; the coarse one is what production calls.

```python
# L1 — the kernel. The oracle. Per pair.
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

**Why L3 is the boundary.** L1 alone would deliver ~1.4x (§1). L3 captures the parse,
the scan, the consensus, the construction, the formatting and the histograms —
everything in the 8.34 µs except the gzip and the write. It is also the level at which
the GIL can be released for a useful span of work.

**Chunking protocol.** Python appends raw reads to a per-stream leftover buffer and calls
`merge_chunk`; C++ consumes whole records only and reports `consumed1`/`consumed2`
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

### 7.4 Allocation

The prototype's 1.43 µs/pair of compute is ~0.47 scan + ~0.96 everything else, and the
"everything else" is dominated by avoidable copies: it unconditionally copies `s1`/`q1`
into mutable buffers (needed only for the 43.5% of pairs the consensus actually
changes), builds `s2rc` byte-at-a-time, reverses `q2` byte-at-a-time, and appends the
merged quality one `push_back` at a time.

The design: **one thread-local `Scratch` arena per worker, sized from
`--max-read-length` at startup and never freed**, holding `s2rc`, the reversed `q2`, the
consensus copies and the record staging buffer. Copy-on-write for the consensus, a
vectorised revcomp, and `memcpy`-based record assembly. No allocation in the per-pair
path at all. **The target is under 1 µs/pair of total compute; that is an estimate, not
a measurement**, and it is the first thing to check once L3 exists.

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
makes any future corpus defect bisectable. Ordered output is unconditional — no flag —
unless it measures worse than ~5%.

Ordering is `for fut in futures_in_order: write(fut.result())` over a bounded submit
window, not a lock-free ring buffer. A ring buffer would mean moving the pool into C++,
which trades testable Python orchestration for a concurrency primitive to get right; and
the head-of-line blocking it avoids is bounded by the *variance* in chunk compute time,
which for fixed-size chunks of Illumina reads is a few percent.

`--processes` stays as a **deprecated alias for `--threads`** so the pipeline's existing
rule keeps working, with a one-line deprecation notice and a note in the JSON stats.

---

## 9. Testing

The port's own instruments stay and are the acceptance gate; two of them were
re-verified against deliberate mutants during the port and must be re-verified again
after each phase here.

1. **The oracle.** `_pymerge` is the reference implementation, not a fallback. It is
   never deleted and never optimised at the cost of clarity. It uses the **same integer
   arithmetic** as C++ (§4) — no float in either scoring loop — so backend equality is a
   property of the spec rather than a coincidence of rounding.
2. **Byte-exact differential, at all three levels.** L1: `(shift, score_q, olen, diff)`
   equal on every pair. L2: identical record tuples. L3: identical blob **and** identical
   counters. Because the tie-break is specified (§5), these are equality assertions, not
   "an argmax" assertions.
3. **Whole-file equality** between backends, and across `--threads 1/2/8`, on the
   200k-pair library.
4. **`test_block_loop_sees_every_position` first.** The roadmap is right that this is the
   test to write before touching the loop: `k += 7` in the unrolled stride changes 6.34%
   of scores and 0.26% of decisions on real data and passes every other test in the
   suite. It is currently the only test that kills that mutant. The vector loop needs its
   own equivalent — a single mismatch swept across all 32 positions of the bail block,
   across the vector boundary at 16, and into the scalar tail.
5. **N and IUPAC parity, explicitly.** §6.3 claims byte comparison needs no special
   handling. Pin it: `N` vs `N` scores as a match, `N` vs `A` as a mismatch, and an
   all-IUPAC pair round-trips identically on both backends.
6. **Fuzz, in the shape of `tests/test_fuzz_roundtrip.py`.** Random lengths (including 0,
   1, and either mate longer), random content including N and IUPAC, random parameters,
   asserting backend equality. That file exists because two 0.3.4 optimisations silently
   corrupted data while 282 tests passed; the same risk applies here exactly.
7. **Independently-drawn mate lengths, everywhere.** The port's non-negotiable. Any new
   geometry fixture that sets `len1 == len2` is structurally blind to the C2 truncation
   class.
8. **ASAN/UBSAN build in CI.** The vector loop reads 16 bytes at a time; the bound
   `k + 32 <= n` is what keeps it inside the record, and an off-by-one there is a heap
   overread that produces plausible output. Sanitisers, not code review.
9. **Mutation checks re-run after each phase**: the C2 regate (`direction == FORWARD &&
   shift > 0`) must fail ≥10 tests including both e2e boundary parametrisations, and
   `k += 7` must fail the block-loop test.

---

## 10. Portability

The tool must be **generally fast on every supported architecture**, with x86-specific
tuning as post-0.4.0 work on a Linux box.

| item | resolution |
|---|---|
| vector compare | 16-byte vectors are the baseline everywhere: NEON on aarch64, **SSE2 on x86-64 — which is part of the x86-64 baseline and needs no `-march` flag and no runtime check**. One `neq16` behind two `#ifdef`s. |
| AVX2 (32 B) | Not baseline. Deliberately **out of scope for 0.4.0**: the 16-byte path is already 5.6x, and §6.1 shows the optimum is set by bail granularity rather than width, so a 32-byte vector may well need bail-every-1 to compete. Measure on x86 first (`bench_simd.cpp` builds the variants), then decide whether a runtime-dispatched AVX2 path earns its complexity. |
| `popcount` | Only needed for the x86 `movemask` result — a 16-bit value. Use a 256-entry table or `__builtin_popcount`; either is fine for 16 bits. The NEON path uses `vaddvq_u8` and needs no popcount at all. This is much smaller than revision 1 implied, because the packed kernel is gone. |
| **`__builtin_popcountll` and SIGILL** | The review warns of `SIGILL` on old hardware. That is not the failure mode: without `-mpopcnt` the compiler emits a software fallback, not the POPCNT instruction. `SIGILL` requires compiling *with* an ISA flag the host lacks — which this design never does. The real risk was only a slow libcall, and it is now moot. |
| endianness | Byte comparison is endian-independent, so the packed kernel's little-endian assumption is gone. `static_assert` retained for the length fields only. |
| unaligned loads | `_mm_loadu_si128` / `vld1q_u8` are penalty-free on all targets of interest. |
| manylinux baseline | **manylinux2014**, not `manylinux_2_28` as the review assumed — 0.3.5 pinned it deliberately so RHEL/CentOS 7-era HPC clusters keep working. SSE2 is available there; AVX2 is not guaranteed. |
| MSVC | `_mm_loadu_si128` and `__popcnt16`/table; no `__builtin_*`. |
| C++ standard | stays `cxx_std_17`; nothing here needs C++20. |
| bindings | **nanobind**, which zna already uses (`pyproject.toml` build-requires, `CMakeLists.txt` `nanobind_add_module`). No change. |

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
  dispatch. *Acceptance: suite green on both backends.*
- **C — C++ L1.** The scan only. *Acceptance: L1 differential on the 200k library + fuzz;
  block-loop test for the vector loop; N/IUPAC parity; ASAN clean.* Expected: scan
  2.63 → 0.47 µs/pair.
- **D — C++ L2.** Consensus, construction, naming, the `Scratch` arena. *Acceptance:
  identical record tuples.* Expected: p1 8.34 → ~4 µs/pair.
- **E — C++ L3 + threads + ordered output.** Parser, formatting, histograms, GIL
  release, `ThreadPoolExecutor`, `--processes` → alias. *Acceptance: whole-file equality
  across backends and thread counts; the four IPC defects each have a named test.*
  Expected: ~2.7 µs/pair at 1 thread, ~0.9–1.1 at 2–4.
- **F — remove numba.** Drop the dependency, the `merge` extra, and the numba-specific
  refusal; `zna merge` now refuses without the *accel backend* unless `--allow-slow`,
  exactly parallel to `is_accelerated()` for the codec. Keep the `njit` no-op shim, which
  the tests import.
- **G — post-0.4.0, on Linux/x86.** Re-run `scripts/merge_bench/` on x86: confirm the
  SSE2 path, evaluate AVX2 with a runtime dispatch, and re-tune the bail granularity,
  which is hardware-dependent. Then, if the gzip floor is still binding, evaluate
  libdeflate or ISA-L inflate inside the extension.

**Then** `zna encode --merge-pairs`, which feeds merged records straight into the
existing `(seq, is_paired, is_read1, is_read2, has_start, has_end)` ingest path via the
second adapter of §7.2, and deletes the FASTQ intermediate, its pigz write, its gzip
read and its full re-parse. That is where the remaining I/O goes, and it is worth more
than anything left in the kernel.

---

## 12. Release

**0.4.0, numba-free** (confirmed). The port is committed and green but the release is
held: phase F removes numba, the `merge` extra and `--allow-slow`'s current meaning, so
publishing them now and withdrawing them one release later is user-visible dependency
churn — a `pip install zna[merge]` and a conda `run_constrained` that stop meaning
anything — bought for nothing.

0.4.0 rather than 0.3.6 because it adds a subcommand, adds a second compiled extension,
changes the parallelism model from processes to threads, changes `--processes` to a
deprecated alias, and removes an optional dependency. On a 0.x project that is
minor-version material. `src/zna/__init__.py` and the CHANGELOG heading are renumbered
when phase F lands, not before.

---

## 13. Deliberately not doing

Each with the number that closes it. [MERGE_TOOL_AUDIT.md](MERGE_TOOL_AUDIT.md) §(c) has
the rest; check it before proposing an improvement.

| idea | why not |
|---|---|
| **2-bit packing + popcount** | 0.535 vs 0.470 µs/pair for byte-wise vectors, and it needs a SWAR packer, cross-word bit realignment, a guard word, and a purity dispatch to keep N/IUPAC semantics. Slower *and* more machinery. |
| **AVX2 in 0.4.0** | Not baseline on manylinux2014, and §6.1 shows bail granularity dominates width. Measure on x86 first (phase G). |
| **Fusing the two flank shifts** | Helps the packed kernel 1.08x, hurts the byte kernel (0.555 vs 0.470). |
| **byte-SWAR in a 64-bit word** | 1.040 vs 1.075 µs/pair — **1.03x**. Same bail granularity, and the compare was never the cost. |
| **Google Highway or another SIMD library** | A new third-party dependency for a project whose runtime deps are `zstandard` and `pyyaml`. One `neq16` behind two `#ifdef`s is ~10 lines. |
| **A lock-free ring buffer for ordered output** | Solves head-of-line blocking that the review's own analysis puts at <3%, by moving the thread pool from Python into C++. §8. |
| **k-mer seed-and-extend to propose shifts** | Audit: implemented and measured, **no gain** over the block loop (1.58x vs 1.57x). |
| **Gating merges on shift ambiguity** | Audit §(c): a microsatellite detector, not an identifiability test. Would re-emit **+4.68% duplicated bases**. |
| **Composition-aware `p_null`** | Audit: measured negative three times; −0.55pp merge rate, and ~2pp on exactly the low-complexity reads it was meant to protect. |
| **Estimating `err_rate` from the data** | Audit: `ê = 0.00866` vs 0.01 assumed; adopting it costs −0.04pp merge rate. |
| **Order-independent differential tests** | Superseded: §5 proves the order *is* a specifiable total order, so tests can be byte-exact instead of weaker. |

---

## 14. Open questions

1. **Everything on x86.** All measurements here are aarch64/NEON. The SSE2 path is
   written but unmeasured, the AVX2 variants are unmeasured, and the optimal bail
   granularity is hardware-dependent — 32 bases won here, and there is no reason to
   assume it wins on Zen or Sapphire Rapids. `scripts/merge_bench/bench_simd.cpp` builds
   all of it. This is phase G and the user's post-release work.
2. **How much of §7.4's ~0.96 µs/pair of non-scan compute survives the `Scratch`
   arena?** Target under 1 µs/pair total; unmeasured.
3. **Does ordered output cost measurable throughput** at 4 threads? Ship it
   unconditionally unless it exceeds ~5%.
4. **Is `--max-read-length 4000` right?** It must exceed any real read and stay inside
   the enumerated domain of §4. A long-read chemistry would need the enumeration re-run,
   not just the flag raised.

---

## 15. External review: what was adopted

A third-party review of revision 1 produced three changes and several corrections.

**Adopted, and they improved the design:**

- **Byte-wise SIMD over 2-bit packing** (§6). The central claim — that a shift becomes a
  pointer offset and no bit realignment is needed — is correct, and measured it wins by
  1.14x while deleting the purity dispatch entirely. Revision 1's contrary conclusion is
  retracted in §6.1.
- **An I/O-agnostic C++ core with adapters above it** (§7.2), so `--merge-pairs` reuses
  the core rather than growing a second copy.
- **Thread-local scratch arenas, no allocation in the per-pair path** (§7.4).
- **Ordered output by default with no flag** unless it measures worse than ~5% (§8) —
  tighter than revision 1's 2% threshold, and the reviewer's variance argument is why.

**Already in revision 1**, and unchanged: the reference oracle using the same integer
spec as C++ (§4, §9.1); no float in either scoring loop; nanobind (zna has used it
since 0.3.x); a documented tie-break.

**Corrections, recorded because they change what to build:**

- **`int32_t` fixed-point would overflow at n = 64** (§4). `int64` throughout.
- **There are no log-odds LUTs in this scoring rule** (§4) — the score is two integer
  multiplies. A `(Q1,Q2)` LUT belongs to the quality-aware variant, which was measured
  and rejected three times.
- **The tie-break needs proof that `n` is not a discriminator**, which the review's
  two-key formulation assumed. §5 supplies it via `gcd(M, STEP) = 1`, and the rule
  simplifies to "maximise score, then minimise `s`" as a result.
- **`__builtin_popcountll` does not risk `SIGILL`** (§10) — without an ISA flag the
  compiler emits a software fallback. The risk was a slow libcall, and the byte kernel
  makes it moot.
- **zna targets manylinux2014, not `manylinux_2_28`** (§10), pinned deliberately in
  0.3.5 to keep RHEL7-era HPC clusters. SSE2 is guaranteed there; AVX2 is not.
- **§14.4 of revision 1 was about the fixed-point proof domain**, not stack vs heap
  allocation. The allocation point was still worth taking, and is now §7.4.
