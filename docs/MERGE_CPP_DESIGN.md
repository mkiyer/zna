# Native acceleration for `zna merge`: design

Status: **implemented and shipped in 0.4.0**, 2026-08-12. What was built matches this
document; §16 records where the measurements moved and what that changed.

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
>
> **Revision 3, 2026-08-12.** `--max-read-length` is **gone** — buffers size themselves
> and the fixed-point spec needs no cap (§4, §7.4); measured, the adaptive arena is 27%
> *faster* than fixed-size per-pair vectors. Everything kept for backwards
> compatibility is **deleted**: `--processes`, `--length-required`, `--allow-slow`, the
> `njit` shim (§7.5). x86/AVX2 tuning moves to 0.4.1 (`docs/FUTURE.md`).

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
| C++ prototype: gz in → file out, single thread | **2.32** |
| C++ prototype: **plain** in → `/dev/null` (compute only) | **1.00** |
| `zna merge -p1` (numba) | 8.34 |
| `zna merge -p8` (numba) | 1.97 |

The writer is free (pigz level 1 in a subprocess keeps up). Decompression is not:
**0.80 µs/pair is the floor of the current architecture**, and `-p4` makes it *worse*
because pigz threads only the CRC and read-ahead, not the inflate.

So the target is: single-threaded C++ at **2.32 µs/pair (3.6x over p1, on one core
instead of eight)**, and since compute is 1.00 µs/pair against a 0.80 µs/pair
decompression floor, **two worker threads are enough to make the tool purely I/O bound
at ~0.85–1.0 µs/pair — about 9x over p1 and 2x over today's p8, on a quarter of the
cores.** Beyond that, gzip is the wall, and the answer to gzip is to stop writing the
intermediate at all (§11).

---

## 3. The prototype, and what it already proves

A complete single-threaded C++ prototype (parse → scan → consensus → build → format →
write) is committed at
[`scripts/merge_bench/proto_merge.cpp`](../scripts/merge_bench/proto_merge.cpp), with
every benchmark behind the numbers in this document beside it and a
[README](../scripts/merge_bench/README.md) that reproduces them. On the 200,000-pair
library it emits **a byte-identical file to `zna merge`** — 177,296 merged, 5,391
trimmed, 17,313 kept, 222,704 records, `cmp` clean — at **2.32 µs/pair against 8.34**.

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
3.6×10⁻⁴ bits (T=8). 2^20 — the obvious "millionths of a bit" choice, and what the
prototype first used — is not good enough; 2^24 is.

**No `--max-read-length`, and no cap of any kind.** Revision 1 proposed one to make the
enumerated domain the reachable domain. It is not needed, for three separate reasons:

- **`int64` cannot overflow on any input that fits in memory.** `n*M` overflows at
  n = 2.8×10¹¹ bases and `d*STEP` at d = 6.7×10¹⁰. A FASTQ record that long cannot be
  read.
- **Fixed point is the specification, not an approximation of float.** It is exact and
  deterministic at every read length; there is no domain in which it becomes
  "unreliable". The enumeration is not a safety guard — it is a *characterisation* of
  where adopting the spec perturbs the behaviour being replaced.
- **And that domain is unreachable anyway.** Extending the enumeration with no cap: at
  2^24 the first `(n, d)` on which fixed point and float disagree is
  **n = 32,830 (T=28)** and n = 34,632 (T=8). That is an *overlap* of 32,830 bases,
  needing reads at least that long — two orders of magnitude beyond any Illumina read,
  and the O(L²) scan makes such input impractical for unrelated reasons (the audit
  measured 40 ms/pair at L = 20,000).

So the tool has **three parameters**, as the redesign intended, and read length is not
one of them. Buffers size themselves (§7.4).

**Derive `M`/`STEP` in Python only, and pass the integers across.** `log2` is not
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

### 7.4 Allocation: a self-sizing arena, measured

**One thread-local `Scratch` arena per worker**, holding `s2rc`, the reversed `q2`, the
consensus copies and the record staging buffer. It **starts at 1024 bases and doubles
whenever a longer read appears**, so nothing needs to know the read length up front:

```cpp
inline void ensure(size_t n) {                 // one predictable, never-taken branch
    if (n <= cap) return;
    size_t c = cap ? cap : 1024;
    while (c < n) c <<= 1;
    s2rc.resize(c); q2r.resize(c); s1b.resize(c); q1b.resize(c);
    seq.resize(2*c); qual.resize(2*c);         // a merged record is at most len1+len2
    cap = c;
}
```

Plus **copy-on-write for the consensus** — 56.5% of pairs have a clean winning overlap
and need no mutable copy of `s1`/`q1` at all — and `memcpy`-based record assembly
instead of per-byte appends.

Measured on the 200k library, both variants byte-identical:

| | per-pair `std::vector` (`.assign`/`.resize`) | **adaptive arena + copy-on-write** |
|---|---:|---:|
| end to end, gz in → file out | 2.71 µs/pair | **2.32** |
| compute only | 1.38 µs/pair | **1.00** |
| arena growths over 200,000 pairs | — | **1** |

So growing on demand is not a speed compromise — it is **27% faster** than the
fixed-size version, because dropping the fixed-size assumption is what made
copy-on-write natural. The capacity check costs nothing measurable: it is taken once
per run, and the branch predictor sees it never taken thereafter. Revision 1's estimate
("under 1 µs/pair, unmeasured") is now a measurement: **1.00 µs/pair**, of which 0.47 is
the scan.

The one real downside is inherent to the algorithm rather than to the arena: the scan is
O(L²), so a 20 kb read costs ~40 ms (audit). The tool logs one INFO line the first time
a read exceeds 1024 bases, naming the quadratic cost, so an accidental long-read FASTQ
is diagnosable instead of looking like a hang. It is not an error and there is no flag.

---

### 7.5 The CLI, with nothing kept for compatibility

zna is alpha and forward-only. Every flag that exists to keep something older working is
deleted rather than deprecated:

| deleted | why it existed | replacement |
|---|---|---|
| `--processes` | the process pool | `--threads` |
| `--length-required` | a fastp alias | `--min-read-length` |
| `--allow-slow` | the numba escape hatch | `--backend python` |
| `--max-read-length` | never shipped; §4 | nothing — buffers self-size |
| the `njit` no-op shim in `overlap.py` | numba absence | deleted with numba; `legacy_scan` in the parity test becomes plain Python |
| sparse per-chunk histograms | pickling cost across processes | dense arrays; threads do not pickle |

`--threads` now means **merge worker threads**, which is what a user expects it to mean;
it previously meant pigz threads. Gzip threading becomes `--io-threads`, used for the
*output* only — the reader is always `-p1` because pigz cannot parallelise inflate, and
§2 measures `-p4` as actively worse. Both changes land together, so nothing silently
reinterprets an existing invocation.

The resulting surface, three scoring parameters and no read-length knob:

```
zna merge --in1 R1.fq.gz --in2 R2.fq.gz --out merged.fq.gz [options]

  --threshold-merge BITS    28.0    score >= this merges the pair
  --threshold-trim  BITS     8.0    ...trims R2's redundant 3' end
  --min-read-length N         40    drop emitted reads shorter than this

  --threads N          min(4,ncpu) merge worker threads; 2-3 saturates (§2)
  --io-threads N                 4  pigz threads for the gzipped output
  --chunk-size N              2000  read pairs per work unit
  --compress-level N             1  output is an intermediate; speed beats ratio
  --backend {auto,accel,python}     auto = accel, hard error if unavailable

  --json FILE   --no-sync-check   --allow-empty   -q/--quiet
```

`--backend` replaces `--allow-slow` and does more: it is the supported way to force the
reference implementation, mirrors `get_backend(name)` for the codec, and gives the test
suite a first-class way to run either side.

**Two cross-repo interfaces change, and neither is "backwards compatibility" — they are
live contracts with other repos, so they need coordinated edits:**

- **The JSON stats.** `"numba": bool` becomes `"backend": "accel"|"python"`, and the
  `params` block gains `score_scale`, `match_q`, `step_q`, `threshold_merge_q`,
  `threshold_trim_q` so a corpus can be audited against the exact integers that produced
  it (§4). hulkrna's `gather/tools/read_merge.py` takes whatever scalars it finds, so the
  only breakage is a vanished `numba` cohort field; the new string key needs the
  `(int, float, str)` filter that R5 already widened. Do it now, while nothing has
  shipped.
- **The merged-read name token** `<id> merged_<n1>_<n2>` stays, because khorana's
  `parse_merged_fastq` requires the last token to start with `merged`. It disappears on
  its own when `--merge-pairs` deletes the FASTQ intermediate.

---

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
| AVX2 (32 B) | Not baseline. **Deferred to 0.4.1** — see `docs/FUTURE.md`. The 16-byte path is already 5.6x, and §6.1 shows the optimum is set by bail granularity rather than width, so a 32-byte vector may well need bail-every-1 to compete. `bench_simd.cpp` builds the variants; measure on x86 before writing any dispatch. |
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
  scoring, the exact `dmax`, and the tie-break spec written down and tested. Still
  numba. *Acceptance: byte-identical output on the 200k library, and the §4 enumeration
  as a test.* This isolates "did fixed point change anything?" from "did C++ change
  anything?".
- **B — backend seam.** Add `backend.py` and `_pymerge.py`; make `overlap.py`/`pairs.py`
  dispatch. *Acceptance: suite green on both backends.*
- **C — C++ L1.** The scan only. *Acceptance: L1 differential on the 200k library + fuzz;
  block-loop test for the vector loop; N/IUPAC parity; ASAN clean.* Expected: scan
  2.63 → 0.47 µs/pair.
- **D — C++ L2.** Consensus, construction, naming, the `Scratch` arena. *Acceptance:
  identical record tuples.* Expected: p1 8.34 → ~4 µs/pair.
- **E — C++ L3 + threads + ordered output.** Parser, formatting, dense histograms, GIL
  release, `ThreadPoolExecutor`, and the §7.5 CLI. *Acceptance: whole-file equality
  across backends and thread counts; the four IPC defects each have a named test.*
  Expected: **2.32 µs/pair at 1 thread, ~0.85–1.0 at 2–3.**
- **F — remove numba, and the rest of the legacy.** Drop the dependency, the `merge`
  extra, the `njit` shim and `legacy_scan`'s use of it; `zna merge` refuses without the
  accel backend unless `--backend python`, exactly parallel to `is_accelerated()` for
  the codec. Release 0.4.0.
- **G — 0.4.1, on Linux/x86.** Re-run `scripts/merge_bench/` on x86: confirm the SSE2
  path, evaluate AVX2 with a runtime dispatch, and re-tune the bail granularity, which
  is hardware-dependent. Then, if the gzip floor is still binding, evaluate libdeflate
  or ISA-L inflate inside the extension. Tracked in `docs/FUTURE.md`.

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
changes the parallelism model from processes to threads, removes an optional dependency,
and deletes `--processes`, `--length-required` and `--allow-slow` outright (§7.5).
zna is alpha and forward-only, so those are deletions rather than deprecations — but
they are still the reason this is a minor bump rather than a patch. `src/zna/__init__.py`
and the CHANGELOG heading are renumbered when phase F lands, not before.

---

## 13. Deliberately not doing

Each with the number that closes it. [MERGE_TOOL_AUDIT.md](MERGE_TOOL_AUDIT.md) §(c) has
the rest; check it before proposing an improvement.

| idea | why not |
|---|---|
| **2-bit packing + popcount** | 0.535 vs 0.470 µs/pair for byte-wise vectors, and it needs a SWAR packer, cross-word bit realignment, a guard word, and a purity dispatch to keep N/IUPAC semantics. Slower *and* more machinery. |
| **AVX2 in 0.4.0** | Not baseline on manylinux2014, and §6.1 shows bail granularity dominates width. Deferred to 0.4.1; measure on x86 first (phase G). |
| **A `--max-read-length` cap** | Three independent reasons it is unnecessary (§4): `int64` cannot overflow on any readable input, fixed point is exact at every length, and the first fixed-point/float disagreement needs a 32,830-base *overlap*. Buffers self-size, and doing so measured 27% faster (§7.4). |
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
   all of it. This is 0.4.1.
2. **Does ordered output cost measurable throughput** at 4 threads? Ship it
   unconditionally unless it exceeds ~5%.
3. **Does the thread count want a smarter default than `min(4, ncpu)`?** Compute is
   1.00 µs/pair against a 0.80 µs/pair decompression floor, so 2 threads should
   saturate; 4 is a margin for slower machines and higher `--compress-level`. Confirm
   once phase E exists.

*Closed since revision 2:* how much non-scan compute survives the arena — **measured at
1.00 µs/pair total, 0.47 of it the scan** (§7.4); and whether `--max-read-length` is
needed — **no** (§4).

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
  allocation. The allocation point was still worth taking, and is now §7.4 — where it
  turned out to be worth 27%, and to remove the need for the cap entirely.


---

## 16. What shipped, and where the measurements moved

All six phases landed in order, each byte-identical to the one before it on the
200,000-pair library. Final state:

| | µs/pair |
|---|---:|
| original numba tool, `--processes 1` | 8.34 |
| original numba tool, `--processes 8` | 1.97 |
| **0.4.0, `--threads 1`** | **2.78** |
| **0.4.0, `--threads 2`** | **1.40** |
| 0.4.0, `--threads 4` | 1.43 |

**6.0x single-threaded, and 1.4x faster than the old 8-process configuration on a
quarter of the cores.** The scan alone went 2.633 → 0.470 µs/pair (5.6x).

The design's predictions held: 2 threads saturate (§8), compute lands near 1 µs/pair
(§7.4 measured 1.00), and the tool ends up I/O bound (§2). With gzip removed from both
ends it runs at **0.42 µs/pair**; decompression costs 1.42 and compression 0.52.

**Four things measurement corrected during implementation.**

1. **Dense histograms across the boundary cost 0.004 µs/pair.** §8 budgeted a sparse
   encoding to replace the one that existed for pickling. Unnecessary — dense lists are
   free once nothing is pickled.
2. **Bigger reads made it slower.** 4 MiB blocks measured 2.48 µs/pair against 256
   KiB's 1.71, because a large blocking read starves the workers while it completes.
   The fix was not smaller reads but *overlapping* them: each stream prefetches its next
   block on its own thread, worth 1.71 → 1.40. §7.3's chunk protocol says nothing about
   read latency and should have.
3. **The 2^20 scale would have been wrong**, as §4 predicted, and the arena was worth
   27% rather than the "no downside" §7.4 claimed — dropping the fixed-size assumption
   is what made copy-on-write natural.
4. **`process_pair` in C++ was worth more than estimated**: p1 went 5.97 → 3.34 µs/pair
   at phase D against an estimate of ~4.

**Two defects the work turned up, both caught by tests rather than review.**

- **The two extensions collided.** Both are imported as `_accel`, so both CMake targets
  emitted `_accel.cpython-*.so` into the same build directory and the merge one
  overwrote the codec's. `zna._accel` became the merge scan and `zna.is_accelerated()`
  silently returned False, while the entire merge suite stayed green. Separate
  `LIBRARY_OUTPUT_DIRECTORY`, pinned by `TestExtensionsAreDistinct`.
- **Random sequence never ties, and I built the differential tests from random
  sequence.** Reversing the two flank shifts in the C++ kernel changes the winner on
  every tied pair and passed all 19 original cross-backend tests. This document already
  said ties must be constructed deliberately (§5) and the tests did not. `tie_fixtures()`
  now builds both kinds — plateau ties from unequal-length periodic mates, and *flank*
  ties from equal-length mates held mutually out of phase, which is the only
  construction that separates the two flanks. A related gap: `TestArgmaxTotalOrder` was
  exercising whichever backend was *selected*, so with the extension built the oracle
  went unchecked, making the suite circular. It is now parametrised over every available
  backend.

**A third defect, found by a question rather than a test.** §4 says there is no
read-length cap and §7.4 says the arena "doubles when a longer read appears". The
compiled backend did neither: `next_record` sized its scratch to 1024 bases *before*
parsing and threw on anything longer, so reads over 1024 bp failed on `accel` and
merged fine on `python` — the two backends disagreeing on a whole class of input, which
is exactly what the oracle arrangement exists to prevent. Nothing in the suite fed
`merge_chunk` a long read. Locating a record is now separate from copying it, so the
arena is sized from the record's real length; the boundary is swept in both directions
by `test_reads_longer_than_the_arena`, and `max_read_length` is reported in the JSON
with a one-line notice past 1024 bp, since the scan is O(L²) and that is what explains
a slow run.

**One mutant survives and should.** `dmax = (ceiling - best)` instead of
`(ceiling - best - 1)` is equivalent: a tie scores `== best` and `v > best` rejects it
either way. Recorded so nobody writes a test to chase it.
