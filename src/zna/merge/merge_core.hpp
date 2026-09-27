/**
 * The overlap scan, with no Python and no I/O in sight.
 *
 * Header-only and dependency-free on purpose (docs/METHODS.md, "Layering"): the same
 * core has to serve `zna merge`'s FASTQ path, `zna encode --merge-pairs` later, and a
 * standalone sanitizer driver that cannot link against Python. Keeping it separate from
 * the bindings is what stops FASTQ assumptions leaking into the kernel.
 *
 * The algorithm is `zna/merge/_pymerge.py`, which is the reference oracle: this must
 * agree with it EXACTLY, not approximately, and `tests/test_merge.py` asserts that pair
 * by pair. Two properties make that achievable:
 *
 *   - **Scores are integers.** `score = n*match_q - d*step_q` in int64, in the
 *     fixed-point scale of `zna/merge/params.py`. No float takes part in a comparison,
 *     a bail bound or the argmax, so the result does not depend on the compiler, the
 *     optimisation level, or `-ffast-math`. The weights are derived once in Python and
 *     passed in as integers -- deriving them here from log2() would invite a 1-ULP
 *     disagreement with the Python side, since log2 is not correctly rounded and libm
 *     differs between platforms.
 *
 *   - **The argmax is a specified total order**, not an artifact of iteration order:
 *     maximise score, then minimise s. The visiting order below (plateau first at
 *     maximal overlap and ascending s; then the flanks at decreasing overlap, the
 *     read-through side -- the smaller s -- first) together with strict `>` realises
 *     exactly that. **Do not reorder these loops.** Swapping the two flanks changes the
 *     winner on every tied pair while leaving random-sequence tests green; it is caught
 *     only by the deliberately-periodic fixtures in `tie_fixtures()`.
 *
 * The inner loop compares **raw bytes**, 16 at a time. Byte comparison is not a
 * compromise for speed: it *is* the reference semantics, for ACGT and equally for N,
 * IUPAC codes and lowercase, so there is no fast path that can disagree with the oracle
 * on any input. A 2-bit packed kernel measured 0.535 us/pair against this one's 0.470
 * and would have needed a packer, cross-word bit realignment and a purity dispatch to
 * keep those semantics. Slower and more machinery.
 *
 * Measured on aarch64/NEON, 50k real pairs, full pruned scans: numba 2.633 us/pair,
 * scalar C++ 1.075, this 0.470 (5.6x), at a 32-base bail.
 *
 * **x86-64 was re-measured in 0.5.2 and the kernel changed under it.** Every figure
 * above is Apple silicon; the first Linux/x86 measurement found the SSE2 path running
 * but reducing each vector compare through a PLT call to libgcc's `__popcountdi2`,
 * which cost it 1.88x and had collapsed the SIMD win to 1.33x over scalar. It now
 * reduces with `psadbw`, in groups of three vectors, with a 16-byte step before the
 * byte tail: 3.14 us/pair -> 1.54 on a Xeon E5-2680 v3, and no ISA flag or runtime
 * dispatch is involved. `neq16x` and `BAIL_VECTORS` below carry the numbers, including
 * why AVX2 is not built.
 */
#ifndef ZNA_MERGE_CORE_HPP
#define ZNA_MERGE_CORE_HPP

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>     // Scratch's arenas. Included explicitly: this header documents
                      // itself as standalone, and it was previously relying on the
                      // translation unit to have pulled <vector> in first.

// 16-byte vectors are baseline on every target we build for: NEON on aarch64, SSE2 on
// x86-64 (part of the base ISA -- no -march flag, no runtime check). Anything wider is
// not baseline, and as of 0.5.2 is measured NOT faster rather than merely "not obviously
// faster": with the reduction done right (`neq16x`), a 32-byte AVX2 kernel came in at
// 1.556 us/pair against baseline SSE2's 1.539. The scan is rejection-dominated, so the
// bail interval and the cost of the horizontal reduction matter; the vector width does
// not. docs/ROADMAP.md has the sweep.
#if defined(__ARM_NEON) || defined(__aarch64__)
#  include <arm_neon.h>
#  define ZNA_MERGE_V16 1
#elif defined(__SSE2__) || defined(_M_X64) || defined(_M_AMD64)
#  include <emmintrin.h>
#  define ZNA_MERGE_V16 1
#endif

namespace zna_merge {

/// Population count of a 16-bit mask, the portable fold.
///
/// **Nothing in the scan calls this any more**, and the reason is the point. It existed
/// for `neq16`, which reduced each vector compare with `pmovmskb` + a popcount; the x86
/// kernel now reduces with `psadbw` (see `neq16x`) and needs no popcount at all. It is
/// kept, exported and tested because the lesson it encodes outlived its caller: a
/// popcount that *requires* POPCNT -- `__builtin_popcount` without `-mpopcnt`, or MSVC's
/// `__popcnt16` -- is not a free primitive on a baseline x86-64 build. GCC turns the
/// first into `callq __popcountdi2@plt`, which is what cost the 0.5.1 kernel 1.88x; the
/// second would emit the instruction itself and fault on any pre-Nehalem CPU. If a
/// future kernel wants a 16-bit popcount, this is the one to use.
///
/// `tests/test_merge.py::TestPopcount` checks it against Python's own bit count over all
/// 65,536 inputs, from whichever platform the suite happens to run on -- which is why it
/// is compiled on every platform rather than only where it would be selected.
inline int popcount16_portable(unsigned x) noexcept {
    x &= 0xFFFFu;
    x -= (x >> 1) & 0x5555u;                     // pairs
    x = (x & 0x3333u) + ((x >> 2) & 0x3333u);    // nibbles
    x = (x + (x >> 4)) & 0x0F0Fu;                // bytes
    return static_cast<int>((x * 0x0101u) >> 8) & 0xFFu;   // sum the two bytes
}

/// Number of differing bytes across `V` consecutive 16-byte windows.
///
/// One horizontal reduction for the whole group, not one per vector. The per-lane
/// equality flags are accumulated in a vector register first -- `V` stays well under
/// 255, so a byte lane cannot overflow -- and folded once at the end. That is the shape
/// the NEON path always wanted (`vaddvq_u8` is the expensive instruction, not
/// `vaddq_u8`) and the shape that lets x86 drop `popcount` entirely.
///
/// **x86 reduces with `psadbw`, not `pmovmskb` + `popcount`.** `psadbw` is baseline
/// SSE2 -- no `-march`, no runtime dispatch, no CPUID -- and is the exact analogue of
/// NEON's `vaddvq_u8`. The kernel it replaces called libgcc's `__popcountdi2` through
/// the PLT twice per 32 bases, because POPCNT is SSE4.2-era and this file is compiled
/// for baseline x86-64. Measured on 50,000 real pairs, full pruned scans, Xeon E5-2680
/// v3 (Haswell), g++ -O3 with no ISA flags: **3.14 us/pair before, 1.54 after (2.04x)**.
/// The same source with `-march=native` and 32-byte AVX2 vectors reaches 1.556 -- i.e.
/// *no better than baseline SSE2 with the right reduction* -- which is why no AVX2 path
/// and no dispatch machinery is built. See docs/ROADMAP.md, 0.4.2.
template <int V>
inline int neq16x(const uint8_t* a, const uint8_t* b) noexcept {
    static_assert(V >= 1 && V <= 255, "byte-lane accumulator would overflow");
#if defined(__ARM_NEON) || defined(__aarch64__)
    const uint8x16_t one = vdupq_n_u8(1);
    uint8x16_t acc = vdupq_n_u8(0);              // per-lane count of EQUAL bytes
    for (int j = 0; j < V; ++j) {
        const uint8x16_t eq = vceqq_u8(vld1q_u8(a + j * 16), vld1q_u8(b + j * 16));
        acc = vaddq_u8(acc, vandq_u8(eq, one));  // vceqq gives 0xFF per equal lane
    }
    return 16 * V - static_cast<int>(vaddvq_u8(acc));
#elif defined(ZNA_MERGE_V16)
    __m128i acc = _mm_setzero_si128();           // per-lane count of EQUAL bytes
    for (int j = 0; j < V; ++j) {
        __m128i va, vb;
        std::memcpy(&va, a + j * 16, 16);
        std::memcpy(&vb, b + j * 16, 16);
        // cmpeq gives 0x00/0xFF; 0xFF is -1 as int8, so SUBTRACTING it adds 1 per equal
        // lane -- two instructions per vector, and no constant to keep live.
        acc = _mm_sub_epi8(acc, _mm_cmpeq_epi8(va, vb));
    }
    // psadbw against zero sums the 16 unsigned byte lanes into two 64-bit halves.
    const __m128i sad = _mm_sad_epu8(acc, _mm_setzero_si128());
    const int eq = _mm_cvtsi128_si32(sad) + static_cast<int>(_mm_extract_epi16(sad, 4));
    return 16 * V - eq;
#else
    int d = 0;
    for (int k = 0; k < 16 * V; ++k) d += (a[k] != b[k]);
    return d;
#endif
}

/// Number of differing bytes in one 16-byte window.
///
/// Kept as its own name because it is the scan's *tail* step and because
/// `tests/test_merge.py` pins `VECTOR_WIDTH` at 16.
inline int neq16(const uint8_t* a, const uint8_t* b) noexcept {
    return neq16x<1>(a, b);
}

/// 16-byte vectors between bail checks in `shift_score`.
///
/// A hardware property, not an algorithm property, so it is set per ISA and each value
/// is the measured optimum on that ISA. Changing it cannot change the scan's answer --
/// the mismatch count is monotone in `k`, so "over budget at some checkpoint" and "over
/// budget at the end" are the same predicate however the checkpoints are spaced -- so
/// this is purely a pruning schedule.
///
///   * **3 (48 bases) on x86-64.** 1.539 us/pair against 1.754 at 2 and 1.542 at 4
///     (g++ 8.5 -O3, 50k real pairs); g++ 13.2 agrees on the ordering.
///   * **4 (64 bases) on aarch64.** Re-measured for 0.5.3 on Apple M3 Max (clang -O3,
///     50k real pairs, `scripts/merge_bench/bench_simd.cpp` rows S1-S8), after the
///     grouped reduction moved the optimum exactly as the x86 round predicted it might:
///     the shipped value 2 was the WORST of {2,3,4,5} -- bail 32 at 0.506 us/pair, 48 at
///     0.439, **64 at 0.416**, 80 at 0.415 -- then a cliff: 96 at 0.528, 128 at 0.675,
///     because a 96+ base group fits a 150 bp shift at most once and the scan falls to
///     the step loop. 80 and 64 are within noise of each other across repeats, and 64
///     stays clear of that cliff down to ~80 bp reads, so 64 it is: **1.22x over the
///     shipped 0.5.2 aarch64 path, 0 mismatches.** The old per-vector kernel at bail 32
///     (the 0.5.1 path) measured 0.497, so the grouped reduction alone was a 5% win
///     there too -- the ROADMAP's "should be a win or a wash" is now a measured win.
#if defined(__ARM_NEON) || defined(__aarch64__)
constexpr int BAIL_VECTORS = 4;
#else
constexpr int BAIL_VECTORS = 3;
#endif

constexpr int64_t REJECT = -(static_cast<int64_t>(1) << 62);

/// Score one candidate shift, or REJECT if it cannot strictly beat `best`.
///
/// `score = n*match_q - d*step_q` is monotone in `d` alone, so the largest mismatch
/// count that could still win is known up front and the loop bails the moment it is
/// exceeded. In integers that bound is exact: `score > best` iff
/// `d <= (ceiling - best - 1) / step_q`.
inline int64_t shift_score(const uint8_t* s1, const uint8_t* s2rc,
                           int s, int n, int64_t match_q, int64_t step_q,
                           int64_t best, int* out_d) noexcept {
    const int64_t ceiling = static_cast<int64_t>(n) * match_q;
    if (ceiling <= best) return REJECT;
    const int64_t dmax = (ceiling - best - 1) / step_q;

    const uint8_t* a = s1   + (s > 0 ?  s : 0);   // a shift is just a pointer offset
    const uint8_t* b = s2rc + (s < 0 ? -s : 0);
    int64_t d = 0;
    int k = 0;
#ifdef ZNA_MERGE_V16
    // Vector groups, bailing between groups; then single vectors; then bytes. Every
    // guard is `k + <step> <= n`, so no loop reads a byte the caller did not supply.
    //
    // The middle loop is what makes a wide bail interval affordable. Without it the
    // group loop's leftover -- up to 16*BAIL_VECTORS-1 bases, 47 of a 150 bp read --
    // falls to the byte loop, and a wider interval pays for its cheaper reduction with
    // a longer scalar tail. With it the tail is under 16 bases whatever the interval,
    // which is worth a further 1.08x on x86 (1.664 us/pair -> 1.539).
    constexpr int GROUP = 16 * BAIL_VECTORS;
    for (; k + GROUP <= n; k += GROUP) {
        d += neq16x<BAIL_VECTORS>(a + k, b + k);
        if (d > dmax) return REJECT;
    }
    for (; k + 16 <= n; k += 16) {
        d += neq16(a + k, b + k);
        if (d > dmax) return REJECT;
    }
#endif
    for (; k < n; ++k) d += (a[k] != b[k]);
    if (d > dmax) return REJECT;
    *out_d = static_cast<int>(d);
    return ceiling - d * step_q;
}

struct ScanResult {
    int shift;
    int64_t score_q;
    int overlap_len;   ///< 0 means nothing reached floor_q
    int mismatches;
};

/// Best-scoring ELIGIBLE shift on the signed single axis, with `floor_q` as the least
/// score that counts. Mirrors `_pymerge.scan`.
///
/// Eligible is every s in [-(len2-1), len1-1] by default. Under `adapter_trimmed` -- the
/// declaration that no read extends past its molecule -- only s >= max(0, len1 - len2),
/// i.e. inferred fragment L = s + len2 >= max(len1, len2). In the loops that is exactly
/// two things: the plateau shrinks to its last shift, and the read-through flank is
/// never visited (every shift on it is < plo <= 0). Restricting eligibility removes
/// shifts from the visiting order without reordering the rest, so the argmax total
/// order holds over the eligible set unchanged.
///
/// **The contract is a template parameter, not a runtime test in the flank loop.** With
/// `if (!adapter_trimmed)` inside the loop the undeclared scan measured 0.494 us/pair on
/// 1M chr22 pairs -- 6% SLOWER than 0.5.3's 0.466 despite a floor 20 bits higher; as two
/// instantiations it is 0.400, and `process_pair` 0.612 -> 0.519 (declared 0.362 ->
/// 0.337, dev panel 0.709 -> 0.609 raw and 0.392 -> 0.363 declared; Apple M3 Max,
/// clang -O3, median of 7).
template <bool ADAPTER_TRIMMED>
inline ScanResult scan_eligible(const uint8_t* s1, int len1, const uint8_t* s2rc,
                                int len2, int64_t match_q, int64_t step_q,
                                int64_t floor_q) noexcept {
    const int nmax = len1 < len2 ? len1 : len2;
    if (nmax <= 0) return {0, 0, 0, 0};

    int64_t best = floor_q - 1;        // a score exactly equal to floor_q must win
    int best_s = 0, best_n = 0, best_d = 0;

    // The plateau: every eligible shift achieving the maximal overlap, ascending.
    const int phi = (len1 >= len2) ? len1 - len2 : 0;
    const int plo = ADAPTER_TRIMMED ? phi : ((len1 >= len2) ? 0 : len1 - len2);
    for (int s = plo; s <= phi; ++s) {
        int d = 0;
        const int64_t v = shift_score(s1, s2rc, s, nmax, match_q, step_q, best, &d);
        if (v > best) { best = v; best_s = s; best_n = nmax; best_d = d; }
    }

    // Then both flanks in lockstep at decreasing overlap length, read-through first.
    for (int n = nmax - 1; n > 0; --n) {
        if (static_cast<int64_t>(n) * match_q <= best) break;
        int d = 0;
        int64_t v;
        if (!ADAPTER_TRIMMED) {
            v = shift_score(s1, s2rc, n - len2, n, match_q, step_q, best, &d);
            if (v > best) { best = v; best_s = n - len2; best_n = n; best_d = d; }
            d = 0;
        }
        v = shift_score(s1, s2rc, len1 - n, n, match_q, step_q, best, &d);
        if (v > best) { best = v; best_s = len1 - n; best_n = n; best_d = d; }
    }

    if (best_n == 0) return {0, 0, 0, 0};
    return {best_s, best, best_n, best_d};
}

inline ScanResult scan(const uint8_t* s1, int len1, const uint8_t* s2rc, int len2,
                       int64_t match_q, int64_t step_q, int64_t floor_q,
                       bool adapter_trimmed = false) noexcept {
    return adapter_trimmed
        ? scan_eligible<true>(s1, len1, s2rc, len2, match_q, step_q, floor_q)
        : scan_eligible<false>(s1, len1, s2rc, len2, match_q, step_q, floor_q);
}

/// The two exact policy tables, as `zna/merge/params.py` derives them (int64, borrowed
/// from the caller -- nothing here copies, owns or derives a table).
///
///   * `t_q[N]`, N = len1 + len2 - 1: the pair's merge floor log2(N / alpha) in fixed
///     point. Covers reads up to `t_len / 2` bases.
///   * `dfit[n]`: the most mismatches a TRUE overlap of n bases shows with probability
///     >= alpha. Covers overlaps up to `dfit_len - 1`.
///
/// No float and no libm reaches a decision: both are integers computed once per run in
/// Python with `decimal`/`Fraction`, so the two backends compare the same numbers.
struct Tables {
    const int64_t* t_q;
    size_t t_len;
    const int64_t* dfit;
    size_t dfit_len;

    /// The longest read both tables cover. Mirrors `_pymerge.table_capacity`.
    int64_t capacity() const noexcept {
        const int64_t cap = static_cast<int64_t>(t_len / 2);
        const int64_t dcap = static_cast<int64_t>(dfit_len) - 1;
        return cap < dcap ? cap : dcap;
    }
};

/// `overlap` verdicts. Mirror `VERDICT_*` in `_pymerge.py`.
constexpr int VERDICT_NONE = 0;
constexpr int VERDICT_MERGE = 1;
constexpr int VERDICT_IMPLAUSIBLE = 2;

/// One pair's overlap decision (`_pymerge._decide`): the scan's winner W, the verdict on
/// it, and W's positions where exactly one mate / both mates read `N`.
struct Decision {
    int verdict;
    ScanResult w;       ///< all zero for VERDICT_NONE
    int one_n;          ///< one-sided N positions in W (all of them mismatches)
    int both_n;         ///< N-against-N positions in W (all of them matches)
};

/// `(one_sided, both)`: positions of the overlap at shift `s` where exactly one mate
/// reads `N`, and where both do. Mirrors `_pymerge._n_positions`.
///
/// Neither says anything about whether the mates agree. A one-sided N is a mismatch in
/// the scan, so the gate discounts it; N against N is a match, so the detected-overlap
/// rate discounts it from the compared bases. Byte `N` only: the parser upper-cases,
/// and IUPAC codes compare as themselves.
///
/// Gated on `memchr` over the two overlap slices, which is what keeps an N-free overlap
/// -- nearly every one -- at the cost of two vectorised searches rather than a counting
/// loop. The reference gates on the whole reads; the count is identical either way,
/// since an N outside the overlap contributes nothing to it.
inline void n_positions(const uint8_t* s1, const uint8_t* s2rc, int s, int n,
                        int* one, int* both) noexcept {
    *one = *both = 0;
    if (n <= 0) return;
    const uint8_t* a = s1 + (s > 0 ? s : 0);
    const uint8_t* b = s2rc + (s < 0 ? -s : 0);
    if (!std::memchr(a, 'N', static_cast<size_t>(n)) &&
        !std::memchr(b, 'N', static_cast<size_t>(n))) {
        return;
    }
    int o = 0, t = 0;
    for (int k = 0; k < n; ++k) {
        const int x = a[k] == 'N', y = b[k] == 'N';
        o += x ^ y;
        t += x & y;
    }
    *one = o;
    *both = t;
}

/// The authoritative overlap decision for one pair (plan §2; `_pymerge._decide`).
///
///   W = the best eligible shift with the pair's own floor T_q[len1 + len2 - 1];
///       nothing reaching it is VERDICT_NONE.
///   d_inf = W's mismatches minus its one-sided N positions;
///       d_inf > dfit[n_W] is VERDICT_IMPLAUSIBLE -- abstain, never re-place: nothing
///       is searched for in W's place.
///   otherwise VERDICT_MERGE.
///
/// The gate is one table lookup on the winner and adds no state to the scan loop. The
/// caller guarantees the tables cover max(len1, len2) (`Tables::capacity`); then
/// N <= 2*cap - 1 < t_len and n_W <= cap <= dfit_len - 1, so both lookups are in range.
inline Decision decide(const uint8_t* s1, int len1, const uint8_t* s2rc, int len2,
                       int64_t match_q, int64_t step_q, const Tables& tab,
                       bool adapter_trimmed) noexcept {
    Decision out{VERDICT_NONE, {0, 0, 0, 0}, 0, 0};
    if (len1 <= 0 || len2 <= 0) return out;
    const ScanResult w = scan(s1, len1, s2rc, len2, match_q, step_q,
                              tab.t_q[len1 + len2 - 1], adapter_trimmed);
    if (w.overlap_len == 0) return out;
    out.w = w;
    n_positions(s1, s2rc, w.shift, w.overlap_len, &out.one_n, &out.both_n);
    out.verdict = (w.mismatches - out.one_n > tab.dfit[w.overlap_len])
                      ? VERDICT_IMPLAUSIBLE : VERDICT_MERGE;
    return out;
}

// ===========================================================================
// Level 2: one pair -- consensus, decision, record construction.
//
// Mirrors `process_pair` in `zna/merge/_pymerge.py` exactly. Everything writes into a
// caller-owned Scratch, so the per-pair path allocates nothing.
// ===========================================================================

/// Complement table. A/C/G/T/N in both cases; **everything else passes through
/// uncomplemented**, which is what `bytes.maketrans` does on the Python side and is
/// deliberate: remapping IUPAC codes to N would change the kernel's N-vs-N semantics.
/// rc(b"RYKMSWBDHVN") == b"NVHDBWSMKYR".
struct ComplementTable {
    uint8_t t[256];
    constexpr ComplementTable() : t() {
        for (int i = 0; i < 256; ++i) t[i] = static_cast<uint8_t>(i);
        t[(uint8_t)'A'] = 'T'; t[(uint8_t)'T'] = 'A';
        t[(uint8_t)'C'] = 'G'; t[(uint8_t)'G'] = 'C';
        t[(uint8_t)'N'] = 'N';
        t[(uint8_t)'a'] = 't'; t[(uint8_t)'t'] = 'a';
        t[(uint8_t)'c'] = 'g'; t[(uint8_t)'g'] = 'c';
        t[(uint8_t)'n'] = 'n';
    }
};
inline const ComplementTable COMPLEMENT{};

inline void revcomp_into(const uint8_t* s, int n, uint8_t* out) noexcept {
    for (int i = 0; i < n; ++i) out[i] = COMPLEMENT.t[s[n - 1 - i]];
}

struct Span {
    const uint8_t* p;
    int n;
};

struct Read {
    Span h, s, q;
};

struct OutRec {
    Span h, s, q;
};

/// Pair outcomes: merged into one record or kept as two. There is no third outcome --
/// 0.6 removed the trim band (MERGE_ACCURACY_PLAN.md §2). Mirrors `MERGED`/`KEPT` in
/// `_pymerge.py`.
enum Outcome { OUTCOME_MERGED = 0, OUTCOME_KEPT = 1 };

struct Params {
    int64_t match_q, step_q;
    Tables tab;                  ///< T_q and dfit, borrowed; see Tables
    bool adapter_trimmed;        ///< the declaration: read-through is impossible
    int min_read_length;
    const uint8_t* disagree_q;   ///< 256*256, built in Python (params.py)
    /// What to do with a no-call the overlap could not rescue from the mate.
    /// 0 = keep it (internal/testing), 1 = trim3, 2 = random. Same vocabulary as
    /// `zna encode --npolicy`, deliberately: one flag, one meaning, both tools.
    int npolicy = 1;
    uint64_t rng_seed = 0;       ///< for NPOLICY_RANDOM; see merge_sub_base
};

constexpr int NPOLICY_KEEP = 0;
constexpr int NPOLICY_TRIM3 = 1;
constexpr int NPOLICY_RANDOM = 2;

/// Per-record provenance bits, emitted as the `ZN:i:<bits>` header tag.
///
/// ZNA does not store headers, so the human-readable tokens beside this tag vanish at
/// encode time. This byte is the only per-record provenance that reaches the corpus:
/// `zna encode --label prov:C:ZN` turns it into a `C` (uint8) column, and an absent tag
/// resolves to 0 -- "nothing happened" -- so declaring the column is opt-in and costs
/// files that do not ask for it nothing.
///
/// **There is deliberately no "merged" bit.** This tag carries what would otherwise be
/// LOST, and "merged" is not lost: it is `merged_<n1>_<n2>` in the FASTQ and
/// `IS_FULL_FRAGMENT` in the corpus. Spending a bit on it would put ` ZN:i:1` on ~82% of
/// emitted records to say something two other places already say. The consequence worth
/// knowing: `IS_FULL_FRAGMENT` is only set when `zna encode` is given
/// `--treat-unpaired-as-merged`, so an encode that omits that flag records neither --
/// which is what asking for it means.
///
/// So an absent tag means "nothing happened to this record that you could not already
/// see", and every set bit is a fact with no other home.
///
/// The vocabulary is shared with `zna encode --merge-pairs` (0.5.0), which computes the
/// same `PairResult` and writes the same bits with no FASTQ in between.
///
/// Bit 1 was PROV_TRIMMED until 0.6 removed the trim band. It is retired rather than
/// reused, so a set bit in any corpus still means one thing.
constexpr int PROV_RESCUED   = 2;   ///< >=1 no-call recovered from the mate
constexpr int PROV_NTRIMMED  = 4;   ///< >=1 base removed by --npolicy trim3
constexpr int PROV_NSUBBED   = 8;   ///< >=1 base substituted by --npolicy random

/// splitmix64's finalizer -- the same function as `zna_mix64` in `_accel.cpp` and
/// `_mix64` in `_pycodec.py`. Substitution is position-derived rather than drawn from a
/// running stream, so it cannot depend on how pairs were batched into chunks.
inline uint64_t merge_mix64(uint64_t x) noexcept {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

/// The base substituted for the no-call at `off` of the read keyed `rec`.
inline uint8_t merge_sub_base(uint64_t seed, uint64_t rec, uint64_t off) noexcept {
    return "ACGT"[merge_mix64(seed + 0xBF58476D1CE4E5B9ULL * (rec + 1)
                                   + 0x94D049BB133111EBULL * (off + 1)) & 3ULL];
}

/// Bytes reserved past the header for the provenance tokens appended below.
///
/// Bound, with 11-digit ints: " ZN:i:15"(8) " trim3_N"(18) " subn_N"(17)
/// " rescued_N"(20) " merged_N_N"(31) = 94. That over-counts -- `trim3_` and `subn_` are
/// mutually exclusive, and only a merged record carries `merged_` -- which is the right
/// direction for a buffer bound. 128 is that with room.
///
/// It is a *constant*: the name buffer is sized from the HEADER plus this, never from
/// the read arena. Sizing it from the arena is what overflowed the heap on any FASTQ
/// whose headers outran its reads -- see `Scratch::ensure_name`.
constexpr size_t NAME_RESERVE = 128;

/// One scratch arena per worker. Starts at 1024 bases and doubles when a longer read
/// turns up, so nothing needs to know the read length up front and the per-pair path
/// never allocates. Measured 27% faster than sizing buffers per pair, because dropping
/// the fixed-size assumption is what makes the copy-on-write below natural.
struct Scratch {
    std::vector<uint8_t> s2rc, q2r, s1b, q1b, s2b, seq, qual, name, name2;
    size_t cap = 0;

    void ensure(size_t n) {
        if (n <= cap) return;                 // one predictable, near-never-taken branch
        size_t c = cap ? cap : 1024;
        while (c < n) c <<= 1;
        s2rc.resize(c); q2r.resize(c); s1b.resize(c); q1b.resize(c);
        s2b.resize(c);                           // R2's copy, for --npolicy random only
        seq.resize(2 * c); qual.resize(2 * c);   // a merged record is at most len1+len2
        cap = c;
    }

    /// Size the merged-name buffer from the HEADER, which is the one thing here that is
    /// not a function of read length.
    ///
    /// `name` holds R1's header with the pair suffix stripped plus a
    /// " merged_<n1>_<n2>" tag, so it needs `header + 64` bytes. It used to be resized
    /// inside `ensure()` as `cap + 64` -- from the READ arena -- which overflowed the
    /// heap on any FASTQ whose headers outran its reads: a 16 KB header against 51 bp
    /// reads writes ~15 KB past the end of a 1088-byte buffer. malloc caught it only
    /// sometimes, so the quiet runs were corrupting whatever followed. Separate buffer,
    /// separate reason, separate function.
    void ensure_name(size_t n) {
        if (n <= name.size()) return;
        size_t c = name.size() ? name.size() : 1088;
        while (c < n) c <<= 1;
        name.resize(c);
    }

    /// The second name buffer, for R2 of an unmerged pair.
    ///
    /// A merged pair emits one record and needs one buffer; a kept pair emits two, and
    /// both mates can carry provenance of their own. Sizing this from R2's
    /// header rather than R1's matters -- the two are usually the same length, but
    /// nothing guarantees it.
    void ensure_name2(size_t n) {
        if (n <= name2.size()) return;
        size_t c = name2.size() ? name2.size() : 1088;
        while (c < n) c <<= 1;
        name2.resize(c);
    }
};

struct PairResult {
    OutRec recs[2];
    int n_recs;
    int outcome;
    int n_dropped;
    /// The alignment the pair was merged from -- all zero when no overlap was admitted,
    /// including when one was refused as implausible. Kept when a MERGE verdict falls
    /// back to KEPT because trim3 broke tiling, exactly as the reference reports it.
    int shift;
    int64_t score_q;
    int overlap_len;
    int mismatches;
    int bases_consensus_changed;
    /// 1 when the best alignment reached T but failed the plausibility gate.
    int implausible;
    /// Bases the N policy removed (trim3) or invented (random). Reported so a library
    /// that loses a lot of sequence to no-calls says so, instead of finishing quietly.
    int npolicy_bases;
    /// No-calls the overlap recovered from the mate, which cost nothing. Always R1's:
    /// the merged record takes the overlap from R1, and only it is built from the
    /// consensus.
    int n_rescued;
    /// The N-policy count split per MATE, for the per-record header tokens;
    /// `npolicy_bases == npolicy_1 + npolicy_2` always.
    int npolicy_1, npolicy_2;
    /// The run's diagnostics (plan §4). None of them affects a decision.
    ///   det_bases / det_mismatches: informative positions and informative mismatches
    ///     of the best alignment whenever it reached T -- MERGE and IMPLAUSIBLE alike,
    ///     i.e. BEFORE the gate, which is what a too-low --error-rate makes wrong.
    ///   rt_strong: with the read-through check on, 1 when the UNRESTRICTED best shift
    ///     reaches T as a read-through (L < max(len1, len2)).
    ///   det_len: that detected overlap's length, the n its gate looked up dfit[n] at
    ///     (0 when nothing reached T). Its histogram is what the run's expected share of
    ///     refused true overlaps is computed over.
    int det_bases, det_mismatches, rt_strong, det_len;
    /// Provenance bits per emitted record, parallel to `recs`. See PROV_* above.
    int prov[2];
};

/// Resolve overlap disagreements into R1 alone, by posterior. Returns bases *changed*.
///
/// R1 alone, because the merged record is the only record built from the overlap: it
/// takes the overlap from R1 and R2 contributes only outside it, so R2's copy is
/// discarded. (0.5.x also wrote R2 on its trim path, where each mate kept part of the
/// overlap; the trim path is gone.) A KEPT pair gets no consensus at all -- nothing about
/// it depends on the alignment being right.
///
/// **N rescue.** An `N` carries no base information, so a real call on the other mate
/// beats it whatever the two quality scores say, and the rescued base keeps the
/// surviving mate's own quality rather than a contested-base derating -- there was no
/// contest. Without this the rescue happened only by luck, because an instrument
/// usually assigns an N a low quality: a high-quality N beat a real base and survived
/// into the corpus. Only `N` is rescued, not the IUPAC codes, which do carry partial
/// information.
inline int consensus_r1(uint8_t* s1, uint8_t* q1, const uint8_t* s2rc,
                        const uint8_t* q2r, int s, int olen,
                        const uint8_t* disagree_q, int* rescued) noexcept {
    const int a0 = s > 0 ? s : 0;      // mirrors the scan's overlap alignment
    const int b0 = s < 0 ? -s : 0;
    int changed = 0;
    for (int i = 0; i < olen; ++i) {
        const int ia = a0 + i, ib = b0 + i;
        if (s1[ia] != s2rc[ib]) {
            const uint8_t qa = q1[ia], qb = q2r[ib];
            const bool a_is_n = s1[ia] == 'N', b_is_n = s2rc[ib] == 'N';
            if (a_is_n && !b_is_n) {            // rescue: a real call beats an N
                s1[ia] = s2rc[ib];
                q1[ia] = qb;
                ++changed;
                ++*rescued;
            } else if (b_is_n) {                // R1's real call stands, uncontested
                continue;
            } else if (qb > qa) {               // R2 is the better-supported call
                s1[ia] = s2rc[ib];
                q1[ia] = disagree_q[(size_t)qb * 256 + qa];
                ++changed;
            } else {                            // R1 wins, but contested: derate it
                q1[ia] = disagree_q[(size_t)qa * 256 + qb];
            }
        }
    }
    return changed;
}

/// One mate after the N policy (`_pymerge._npolicy_mate`): the sequence to emit, its
/// length, and the bases the policy touched -- cut off under trim3, substituted under
/// random (`rec` = 2 * pair_index + mate, which is what makes substitution
/// position-derived), none under keep.
///
/// trim3 cuts at the first N, keeping [0, first_N): 3' only, so base 0 -- a fragment
/// terminus -- is never disturbed, and the quality span is the same pointer, shorter.
/// random writes into `buf` (copying `s` there first unless it already IS `buf`, since
/// memcpy with src == dst is undefined) and leaves the length alone.
struct MateOut {
    const uint8_t* s;
    int n;
    int touched;
};

inline MateOut npolicy_mate(const uint8_t* s, int n, uint8_t* buf, int npolicy,
                            uint64_t seed, uint64_t rec) noexcept {
    if (npolicy == NPOLICY_KEEP || n <= 0) return {s, n, 0};
    const uint8_t* first = static_cast<const uint8_t*>(
        std::memchr(s, 'N', static_cast<size_t>(n)));
    if (!first) return {s, n, 0};
    const int k = static_cast<int>(first - s);
    if (npolicy == NPOLICY_TRIM3) return {s, k, n - k};
    if (s != buf) std::memcpy(buf, s, static_cast<size_t>(n));
    int touched = 0;
    for (int i = k; i < n; ++i) {
        if (buf[i] == 'N') { buf[i] = merge_sub_base(seed, rec, i); ++touched; }
    }
    return {buf, n, touched};
}

/// Build one emitted record's name: the input header, then its provenance tokens.
///
/// **Tags pass through untouched.** Everything already in the header is copied verbatim
/// and the tokens are *appended*; nothing is ever removed or rewritten. That is what lets
/// `zna encode --label` read the same `KEY:T:VALUE` tags off a merged record that it
/// would have read off R1, and it is a contract, not an accident -- `strip_suffix` drops
/// only the two bytes of a `/1`/`/2` pair suffix from the ID token itself.
///
/// Token order is fixed, and `merged_<n1>_<n2>` is appended by the caller AFTER these so
/// it stays the final token -- fastp's convention, which costs nothing to keep.
///
/// The colon-less tokens are invisible to ZNA's tag parser, which requires `KEY:T:VALUE`
/// and skips anything else; `ZN:i:<bits>` is the one that is meant to be read, and it is
/// emitted only when some bit is set, so an absent tag resolves to 0 through the label
/// machinery's own missing-value path.
///
/// `cap` is the buffer's capacity; every write is bounded by it. Callers reserve
/// `header + NAME_RESERVE`.
inline int build_name(uint8_t* nm, size_t cap, const Span& h, bool strip_suffix,
                      int bits, int trim3_n, int subn_n, int rescued_n) noexcept {
    int nl = 0;
    if (strip_suffix) {
        int cut = h.n;                  // id_end: first space or tab, else all of it
        for (int i = 0; i < h.n; ++i) {
            if (h.p[i] == ' ' || h.p[i] == '\t') { cut = i; break; }
        }
        int idlen = cut;
        if (idlen >= 2 && h.p[idlen - 2] == '/' &&
            (h.p[idlen - 1] == '1' || h.p[idlen - 1] == '2')) {
            idlen -= 2;
        }
        std::memcpy(nm, h.p, static_cast<size_t>(idlen));                  nl += idlen;
        std::memcpy(nm + nl, h.p + cut, static_cast<size_t>(h.n - cut));
        nl += h.n - cut;
    } else {
        std::memcpy(nm, h.p, static_cast<size_t>(h.n));                    nl = h.n;
    }
    if (bits) {
        nl += std::snprintf(reinterpret_cast<char*>(nm + nl), cap - nl,
                            " ZN:i:%d", bits);
    }
    if (trim3_n) {
        nl += std::snprintf(reinterpret_cast<char*>(nm + nl), cap - nl,
                            " trim3_%d", trim3_n);
    }
    if (subn_n) {
        nl += std::snprintf(reinterpret_cast<char*>(nm + nl), cap - nl,
                            " subn_%d", subn_n);
    }
    if (rescued_n) {
        nl += std::snprintf(reinterpret_cast<char*>(nm + nl), cap - nl,
                            " rescued_%d", rescued_n);
    }
    return nl;
}

/// The name for one emitted mate of an UNMERGED pair: its header, plus whatever
/// provenance tokens it earned.
///
/// Returns the input span **untouched** when there is nothing to say. That is the common
/// case by a wide margin, and it keeps the record zero-copy -- a pointer straight into
/// the caller's input buffer, with no scratch touched and no bytes moved. Only a record
/// the pipeline actually did something to pays for a name.
///
/// `which` selects the buffer: mate 0 uses `name`, mate 1 uses `name2`. A merged pair
/// emits one record and uses `name` alone.
inline Span name_for(Scratch& sc, int which, const Span& h, int bits, int npolicy,
                     int npolicy_n, int rescued_n) {
    if (!bits) return h;
    const int trim3_n = (npolicy == NPOLICY_RANDOM) ? 0 : npolicy_n;
    const int subn_n  = (npolicy == NPOLICY_RANDOM) ? npolicy_n : 0;
    const size_t want = static_cast<size_t>(h.n) + NAME_RESERVE;
    uint8_t* nm;
    size_t room;
    if (which == 0) {
        sc.ensure_name(want);  nm = sc.name.data();  room = sc.name.size();
    } else {
        sc.ensure_name2(want); nm = sc.name2.data(); room = sc.name2.size();
    }
    const int nl = build_name(nm, room, h, /*strip_suffix=*/false,
                              bits, trim3_n, subn_n, rescued_n);
    return {nm, nl};
}

/// Classify one pair and build its output records. Mirrors `_pymerge._process_pair_ex`.
///
/// The decision is `decide`'s verdict, and there are two outcomes:
///
///   MERGE            -> one full-fragment record (R1 wins ties in the posterior
///                       consensus), unless trim3 has cut the mates so far that they no
///                       longer tile the fragment, and then as below
///   NONE/IMPLAUSIBLE -> both mates, unchanged apart from the N policy -- never the
///                       consensus, which only a merged record uses
///
/// Every construction path takes from the 5' end, and trim3 only ever cuts 3' ends, so
/// base 0 of every emitted read stays a true fragment boundary. A merged record is built
/// from the inferred span `L = s + len2`, so its length is `L` identically for every
/// geometry -- do not reintroduce per-direction case analysis.
///
/// `rt_check` asks for the read-through diagnostic (plan §4); the chunk adapters set it
/// for the input's first pairs by input index. Without the declaration every shift is
/// eligible, so the scan's own winner answers it for free; under it, one extra
/// unrestricted scan runs. The caller guarantees the tables cover max(len1, len2).
inline PairResult process_pair(const Read& r1, const Read& r2,
                               const Params& p, Scratch& sc,
                               int64_t pair_index = 0, bool rt_check = false) {
    const int len1 = r1.s.n, len2 = r2.s.n;
    sc.ensure(static_cast<size_t>(len1 > len2 ? len1 : len2));

    uint8_t* s2rc = sc.s2rc.data();
    revcomp_into(r2.s.p, len2, s2rc);

    const Decision dec = decide(r1.s.p, len1, s2rc, len2, p.match_q, p.step_q, p.tab,
                                p.adapter_trimmed);
    const ScanResult& w = dec.w;

    PairResult out{};
    out.implausible = dec.verdict == VERDICT_IMPLAUSIBLE;
    // Compared positions minus the uninformative ones (all zero for VERDICT_NONE).
    out.det_bases = w.overlap_len - dec.one_n - dec.both_n;
    out.det_mismatches = w.mismatches - dec.one_n;
    out.det_len = w.overlap_len;
    if (rt_check && len1 > 0 && len2 > 0) {
        int rs = w.shift, rn = w.overlap_len;
        if (p.adapter_trimmed) {
            const ScanResult u = scan(r1.s.p, len1, s2rc, len2, p.match_q, p.step_q,
                                      p.tab.t_q[len1 + len2 - 1], false);
            rs = u.shift;
            rn = u.overlap_len;
        }
        out.rt_strong = (rn > 0 && rs + len2 < (len1 > len2 ? len1 : len2)) ? 1 : 0;
    }

    // Nothing is built from a refused alignment, and nothing about it is reported as
    // the pair's alignment.
    const bool merge = dec.verdict == VERDICT_MERGE;
    if (merge) {
        out.shift = w.shift;
        out.score_q = w.score_q;
        out.overlap_len = w.overlap_len;
        out.mismatches = w.mismatches;
    }

    const int lr = p.min_read_length;
    const int L = out.shift + len2;        // the inferred fragment length, if merging

    // The consensus is written only into R1, and only on the merge verdict: the merged
    // record takes the overlap from R1. A pair with no admitted overlap is emitted
    // untouched -- an alignment too suspect to merge on is too suspect to rewrite bases
    // on (measured under 0.5.x: of 3,068 kept pairs with a detected overlap, zero had
    // found the true shift, and writing R1 there turned 1,379 correct bases wrong to fix
    // 78).
    //
    // Copy-on-write: a clean winning overlap (56.5% of real pairs) needs no mutable copy
    // of either read at all.
    const uint8_t* S1 = r1.s.p;
    const uint8_t* Q1 = r1.q.p;
    const bool consensus = merge && w.mismatches > 0;
    if (consensus) {
        uint8_t* q2r = sc.q2r.data();
        for (int i = 0; i < len2; ++i) q2r[i] = r2.q.p[len2 - 1 - i];
        uint8_t* s1b = sc.s1b.data();
        uint8_t* q1b = sc.q1b.data();
        std::memcpy(s1b, r1.s.p, static_cast<size_t>(len1));
        std::memcpy(q1b, r1.q.p, static_cast<size_t>(len1));
        out.bases_consensus_changed =
            consensus_r1(s1b, q1b, s2rc, q2r, w.shift, w.overlap_len, p.disagree_q,
                         &out.n_rescued);
        S1 = s1b;
        Q1 = q1b;
    }

    // ---- the N policy, after the rescue, so a no-call the mate could answer costs
    //      nothing. trim3 is 3' only, so both 5' anchors -- the two fragment termini --
    //      are untouched however short the reads get.
    const uint64_t rec1 = static_cast<uint64_t>(pair_index) * 2;
    MateOut m1 = npolicy_mate(S1, len1, sc.s1b.data(), p.npolicy, p.rng_seed, rec1);
    const MateOut m2 = npolicy_mate(r2.s.p, len2, sc.s2b.data(), p.npolicy, p.rng_seed,
                                    rec1 + 1);
    const uint8_t* Q2 = r2.q.p;

    // ---- merge on GEOMETRY, reusing the evidence -----------------------------------
    //
    // The pair still tiles the fragment iff elen1 + elen2 >= L, and when it does the
    // reconstruction IS the fragment, exactly and N-free. Nothing is re-scored: trimming
    // cuts 3' ends, which is where a normal overlap lives, so a re-scan would refuse
    // merges it had ample evidence for a moment earlier. Only trim3 changes a length, so
    // only trim3 can turn a merge verdict into a kept pair here.
    const bool will_merge = merge && (m1.n + m2.n) >= L;

    if (!will_merge && consensus) {
        // A merge verdict that trim3 cut below tiling: the pair is KEPT, and a kept mate
        // is the input with the N policy applied and nothing else (plan §2, and §8's
        // "kept-mate substitutions are zero by construction"). The consensus -- its
        // substitutions, derated qualities and rescues -- existed only to build the
        // merged record, so R1 is re-derived from the input and none of it is counted.
        // Under 0.5.x the rewritten R1 was emitted here (167 kept pairs on the dev
        // panel's N benches, 10 of their substitutions to a wrong base).
        m1 = npolicy_mate(r1.s.p, len1, sc.s1b.data(), p.npolicy, p.rng_seed, rec1);
        Q1 = r1.q.p;
        out.bases_consensus_changed = 0;
        out.n_rescued = 0;
    }

    out.npolicy_1 = m1.touched;
    out.npolicy_2 = m2.touched;
    out.npolicy_bases = m1.touched + m2.touched;

    // Which policy bit a touched record earns. Set per RECORD, from that record's own
    // count -- a kept pair whose R1 lost bases and whose R2 did not says exactly that.
    const int npolicy_bit = (p.npolicy == NPOLICY_RANDOM) ? PROV_NSUBBED : PROV_NTRIMMED;

    bool paired;
    int n_cand;

    if (will_merge) {
        // `shift` is the offset of revcomp(R2) on the shared axis, so it is tied to R2's
        // length. R2 keeps its 5' anchor at fragment position L-1, so a trimmed mate
        // covers [L - elen2, L) and the offset becomes L - elen2. L itself is unchanged
        // -- that is the whole point. s2rc must match the emitted mate, so it is rebuilt
        // when the N policy changed R2 (a cut or a substitution).
        const int elen1 = m1.n, elen2 = m2.n;
        if (m2.s != r2.s.p || elen2 != len2) revcomp_into(m2.s, elen2, s2rc);
        const int s = L - elen2;
        const int take1 = elen1 < L ? elen1 : L;
        const int take2 = L - take1;
        uint8_t* seq = sc.seq.data();
        uint8_t* qual = sc.qual.data();
        std::memcpy(seq, m1.s, static_cast<size_t>(take1));
        std::memcpy(qual, Q1, static_cast<size_t>(take1));
        if (take2) {
            const int b = take1 - s;
            std::memcpy(seq + take1, s2rc + b, static_cast<size_t>(take2));
            for (int i = 0; i < take2; ++i) {
                qual[take1 + i] = Q2[elen2 - 1 - (b + i)];
            }
        }
        // A merged record is built from BOTH mates, so its provenance is the pair's:
        // the trim3/random counts are the two summed, and the rescues are R1's, which
        // are the only ones that reached the emitted bases (consensus_r1 above).
        out.prov[0] = (out.n_rescued ? PROV_RESCUED : 0)
                    | (out.npolicy_bases ? npolicy_bit : 0);

        // fastp-style merged name: "<id> merged_<n1>_<n2>", pair suffix stripped, tags
        // preserved. ZNA ignores the token itself; keeping it LAST is fastp's convention
        // and costs nothing, so the provenance tokens go before it.
        //
        // Sized from the header, not the read: `build_name` reaches `r1.h.n` and the
        // tokens may add NAME_RESERVE more. See `Scratch::ensure_name`.
        const bool rnd = (p.npolicy == NPOLICY_RANDOM);
        sc.ensure_name(static_cast<size_t>(r1.h.n) + NAME_RESERVE);
        uint8_t* nm = sc.name.data();
        int nl = build_name(nm, sc.name.size(), r1.h, /*strip_suffix=*/true, out.prov[0],
                            /*trim3_n=*/rnd ? 0 : out.npolicy_bases,
                            /*subn_n=*/rnd ? out.npolicy_bases : 0,
                            /*rescued_n=*/out.n_rescued);
        nl += std::snprintf(reinterpret_cast<char*>(nm + nl), sc.name.size() - nl,
                            " merged_%d_%d", take1, take2);

        out.recs[0] = {{nm, nl}, {seq, L}, {qual, L}};
        n_cand = 1;
        paired = false;
        out.outcome = OUTCOME_MERGED;
    } else {
        // No admitted overlap, or a merge the N policy broke: keep both reads, touched
        // by the N policy and nothing else -- so a kept record never carries
        // PROV_RESCUED.
        out.prov[0] = out.npolicy_1 ? npolicy_bit : 0;
        out.prov[1] = out.npolicy_2 ? npolicy_bit : 0;
        out.recs[0] = {name_for(sc, 0, r1.h, out.prov[0], p.npolicy, out.npolicy_1, 0),
                       {m1.s, m1.n}, {Q1, m1.n}};
        out.recs[1] = {name_for(sc, 1, r2.h, out.prov[1], p.npolicy, out.npolicy_2, 0),
                       {m2.s, m2.n}, {Q2, m2.n}};
        n_cand = 2;
        paired = true;
        out.outcome = OUTCOME_KEPT;
    }

    // Pair integrity: an unmerged pair is emitted all-or-nothing, because a lone
    // surviving mate would be encoded as a spurious "single" -- a full molecule with
    // both endpoints. A merged read is a genuine full molecule and is filtered alone.
    if (paired) {
        out.n_recs = (out.recs[0].s.n >= lr && out.recs[1].s.n >= lr) ? 2 : 0;
    } else {
        out.n_recs = (out.recs[0].s.n >= lr) ? 1 : 0;
    }
    out.n_dropped = n_cand - out.n_recs;
    return out;
}

}  // namespace zna_merge

#endif  // ZNA_MERGE_CORE_HPP
