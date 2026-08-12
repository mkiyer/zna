/**
 * The overlap scan, with no Python and no I/O in sight.
 *
 * Header-only and dependency-free on purpose (docs/MERGE_CPP_DESIGN.md §7.2): the same
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
 * Measured, 50k real pairs, full pruned scans: numba 2.633 us/pair, scalar C++ 1.075,
 * this 0.470 (5.6x). Bail granularity of 32 bases beat both 16 and 64; see §6.1.
 */
#ifndef ZNA_MERGE_CORE_HPP
#define ZNA_MERGE_CORE_HPP

#include <cstdint>
#include <cstring>

// 16-byte vectors are baseline on every target we build for: NEON on aarch64, SSE2 on
// x86-64 (part of the base ISA -- no -march flag, no runtime check). Anything wider is
// not baseline and, per §6.1, is not obviously faster either: the scan is
// rejection-dominated, so the bail interval matters more than the vector width.
#if defined(__ARM_NEON) || defined(__aarch64__)
#  include <arm_neon.h>
#  define ZNA_MERGE_V16 1
#elif defined(__SSE2__) || defined(_M_X64) || defined(_M_AMD64)
#  include <emmintrin.h>
#  define ZNA_MERGE_V16 1
#endif

namespace zna_merge {

/// Number of differing bytes in a 16-byte window. Unaligned loads are penalty-free on
/// every target of interest.
inline int neq16(const uint8_t* a, const uint8_t* b) noexcept {
#if defined(__ARM_NEON) || defined(__aarch64__)
    const uint8x16_t eq = vceqq_u8(vld1q_u8(a), vld1q_u8(b));
    // vceqq gives 0xFF per equal lane; mask to 1 and sum, so the total cannot overflow
    // a byte lane (16 lanes, max 16).
    return 16 - static_cast<int>(vaddvq_u8(vandq_u8(eq, vdupq_n_u8(1))));
#elif defined(ZNA_MERGE_V16)
    __m128i va, vb;
    std::memcpy(&va, a, 16);
    std::memcpy(&vb, b, 16);
    const int mask = _mm_movemask_epi8(_mm_cmpeq_epi8(va, vb));
    // popcount of a 16-bit value. __builtin_popcount is cheap even without -mpopcnt:
    // the compiler emits a table or a SWAR sequence, never an illegal instruction --
    // POPCNT is only emitted when the ISA is enabled explicitly.
    return 16 - __builtin_popcount(static_cast<unsigned>(mask) & 0xFFFFu);
#else
    int d = 0;
    for (int k = 0; k < 16; ++k) d += (a[k] != b[k]);
    return d;
#endif
}

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
    // Two vectors between bail checks. Measured optimum: 32 bases beats 16 (too much
    // branching) and 64 (too much work past the point of rejection). The guard is
    // `k + 32 <= n`, so the loop never reads a byte the caller did not supply.
    for (; k + 32 <= n; k += 32) {
        d += neq16(a + k, b + k) + neq16(a + k + 16, b + k + 16);
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

/// Best-scoring shift over s in [-(len2-1), len1-1], on the signed single axis.
inline ScanResult scan(const uint8_t* s1, int len1, const uint8_t* s2rc, int len2,
                       int64_t match_q, int64_t step_q, int64_t floor_q) noexcept {
    const int nmax = len1 < len2 ? len1 : len2;
    if (nmax <= 0) return {0, 0, 0, 0};

    int64_t best = floor_q - 1;        // a score exactly equal to floor_q must win
    int best_s = 0, best_n = 0, best_d = 0;

    // The plateau: every shift achieving the maximal overlap, ascending.
    const int plo = (len1 >= len2) ? 0 : len1 - len2;
    const int phi = (len1 >= len2) ? len1 - len2 : 0;
    for (int s = plo; s <= phi; ++s) {
        int d = 0;
        const int64_t v = shift_score(s1, s2rc, s, nmax, match_q, step_q, best, &d);
        if (v > best) { best = v; best_s = s; best_n = nmax; best_d = d; }
    }

    // Then both flanks in lockstep at decreasing overlap length, read-through first.
    for (int n = nmax - 1; n > 0; --n) {
        if (static_cast<int64_t>(n) * match_q <= best) break;
        int d = 0;
        int64_t v = shift_score(s1, s2rc, n - len2, n, match_q, step_q, best, &d);
        if (v > best) { best = v; best_s = n - len2; best_n = n; best_d = d; }
        d = 0;
        v = shift_score(s1, s2rc, len1 - n, n, match_q, step_q, best, &d);
        if (v > best) { best = v; best_s = len1 - n; best_n = n; best_d = d; }
    }

    if (best_n == 0) return {0, 0, 0, 0};
    return {best_s, best, best_n, best_d};
}

}  // namespace zna_merge

#endif  // ZNA_MERGE_CORE_HPP
