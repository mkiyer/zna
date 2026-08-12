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

enum Outcome { OUTCOME_MERGED = 0, OUTCOME_TRIMMED = 1, OUTCOME_KEPT = 2 };

struct Params {
    int64_t match_q, step_q, t_merge_q, t_trim_q;
    int min_read_length;
    const uint8_t* disagree_q;   ///< 256*256, built in Python (params.py)
};

/// One scratch arena per worker. Starts at 1024 bases and doubles when a longer read
/// turns up, so nothing needs to know the read length up front and the per-pair path
/// never allocates. Measured 27% faster than sizing buffers per pair, because dropping
/// the fixed-size assumption is what makes the copy-on-write below natural.
struct Scratch {
    std::vector<uint8_t> s2rc, q2r, s1b, q1b, seq, qual, name;
    size_t cap = 0;

    void ensure(size_t n) {
        if (n <= cap) return;                 // one predictable, near-never-taken branch
        size_t c = cap ? cap : 1024;
        while (c < n) c <<= 1;
        s2rc.resize(c); q2r.resize(c); s1b.resize(c); q1b.resize(c);
        seq.resize(2 * c); qual.resize(2 * c);   // a merged record is at most len1+len2
        name.resize(c + 64);
        cap = c;
    }
};

struct PairResult {
    OutRec recs[2];
    int n_recs;
    int outcome;
    int n_dropped;
    int64_t score_q;
    int overlap_len;
    int mismatches;
    int bases_consensus_changed;
    int trim_guard_fired;
};

/// Resolve overlap disagreements into R1, by posterior. Returns bases *changed*.
///
/// Only R1 is written because the output's overlap always comes from R1 (merge
/// concatenates R1 through the overlap; trim keeps R1 and discards R2's copy), so
/// editing R2's discarded copy could never reach the corpus.
inline int consensus_r1(uint8_t* s1, uint8_t* q1, const uint8_t* s2rc,
                        const uint8_t* q2r, int s, int olen,
                        const uint8_t* disagree_q) noexcept {
    const int a0 = s > 0 ? s : 0;      // mirrors the scan's overlap alignment
    const int b0 = s < 0 ? -s : 0;
    int changed = 0;
    for (int i = 0; i < olen; ++i) {
        const int ia = a0 + i, ib = b0 + i;
        if (s1[ia] != s2rc[ib]) {
            const uint8_t qa = q1[ia], qb = q2r[ib];
            if (qb > qa) {                      // R2 is the better-supported call
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

/// Classify one pair and build its output records.
///
/// Every construction path takes from the 5' end, and trimming only ever cuts R2's 3'
/// end, so base 0 of every emitted read stays a true fragment boundary. A merged record
/// is built from the inferred span `L = s + len2`, so its length is `L` identically for
/// every geometry -- do not reintroduce per-direction case analysis.
inline PairResult process_pair(const Read& r1, const Read& r2,
                               const Params& p, Scratch& sc) {
    const int len1 = r1.s.n, len2 = r2.s.n;
    sc.ensure(static_cast<size_t>(len1 > len2 ? len1 : len2));

    uint8_t* s2rc = sc.s2rc.data();
    revcomp_into(r2.s.p, len2, s2rc);

    const ScanResult r = scan(r1.s.p, len1, s2rc, len2,
                              p.match_q, p.step_q, p.t_trim_q);

    PairResult out{};
    out.score_q = r.score_q;
    out.overlap_len = r.overlap_len;
    out.mismatches = r.mismatches;

    // Copy-on-write: a clean winning overlap (56.5% of real pairs) needs no mutable
    // copy of R1 at all.
    const uint8_t* S1 = r1.s.p;
    const uint8_t* Q1 = r1.q.p;
    if (r.mismatches > 0) {
        uint8_t* q2r = sc.q2r.data();
        for (int i = 0; i < len2; ++i) q2r[i] = r2.q.p[len2 - 1 - i];
        uint8_t* s1b = sc.s1b.data();
        uint8_t* q1b = sc.q1b.data();
        std::memcpy(s1b, r1.s.p, static_cast<size_t>(len1));
        std::memcpy(q1b, r1.q.p, static_cast<size_t>(len1));
        out.bases_consensus_changed =
            consensus_r1(s1b, q1b, s2rc, q2r, r.shift, r.overlap_len, p.disagree_q);
        S1 = s1b;
        Q1 = q1b;
    }

    const int lr = p.min_read_length;
    const bool detected = r.overlap_len > 0;
    bool paired;
    int n_cand;

    if (detected && r.score_q >= p.t_merge_q) {
        // Enough evidence to assert the fragment: collapse to one read.
        const int s = r.shift;
        const int L = s + len2;
        const int take1 = len1 < L ? len1 : L;
        const int take2 = L - take1;
        uint8_t* seq = sc.seq.data();
        uint8_t* qual = sc.qual.data();
        std::memcpy(seq, S1, static_cast<size_t>(take1));
        std::memcpy(qual, Q1, static_cast<size_t>(take1));
        if (take2) {
            const int b = take1 - s;
            std::memcpy(seq + take1, s2rc + b, static_cast<size_t>(take2));
            for (int i = 0; i < take2; ++i) {
                qual[take1 + i] = r2.q.p[len2 - 1 - (b + i)];
            }
        }
        // fastp-style merged name: "<id> merged_<n1>_<n2>", pair suffix stripped, tags
        // preserved. ZNA ignores it; khorana's parse_merged_fastq requires it.
        uint8_t* nm = sc.name.data();
        int cut = 0;
        {   // id_end: first space or tab, else the whole header
            const uint8_t* h = r1.h.p;
            int sp = -1, tb = -1;
            for (int i = 0; i < r1.h.n; ++i) {
                if (h[i] == ' ' && sp < 0) sp = i;
                if (h[i] == '\t' && tb < 0) tb = i;
                if (sp >= 0 && tb >= 0) break;
            }
            cut = (sp < 0) ? (tb < 0 ? r1.h.n : tb) : (tb < 0 ? sp : (sp < tb ? sp : tb));
        }
        int idlen = cut;
        if (idlen >= 2 && r1.h.p[idlen - 2] == '/' &&
            (r1.h.p[idlen - 1] == '1' || r1.h.p[idlen - 1] == '2')) {
            idlen -= 2;
        }
        int nl = 0;
        std::memcpy(nm, r1.h.p, static_cast<size_t>(idlen));               nl += idlen;
        std::memcpy(nm + nl, r1.h.p + cut, static_cast<size_t>(r1.h.n - cut));
        nl += r1.h.n - cut;
        nl += std::snprintf(reinterpret_cast<char*>(nm + nl), 64,
                            " merged_%d_%d", take1, take2);

        out.recs[0] = {{nm, nl}, {seq, L}, {qual, L}};
        n_cand = 1;
        paired = false;
        out.outcome = OUTCOME_MERGED;
    } else if (detected && r.shift >= 0 && r.score_q >= p.t_trim_q &&
               len2 - r.overlap_len >= lr) {
        // Real overlap, but not enough evidence to risk a chimera: keep both reads and
        // trim the redundant overlap off R2's 3' end so it is not counted twice.
        const int keep2 = len2 - r.overlap_len;
        out.recs[0] = {r1.h, {S1, len1}, {Q1, len1}};
        out.recs[1] = {r2.h, {r2.s.p, keep2}, {r2.q.p, keep2}};
        n_cand = 2;
        paired = true;
        out.outcome = OUTCOME_TRIMMED;
    } else {
        // No detectable overlap, an unmergeable read-through, or a trim blocked by the
        // guard: keep both reads exactly as they are.
        if (detected && r.shift >= 0 && r.score_q >= p.t_trim_q) {
            out.trim_guard_fired = 1;
        }
        out.recs[0] = {r1.h, {S1, len1}, {Q1, len1}};
        out.recs[1] = {r2.h, r2.s, r2.q};
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
