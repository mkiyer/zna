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
 * Measured, 50k real pairs, full pruned scans: numba 2.633 us/pair, scalar C++ 1.075,
 * this 0.470 (5.6x). Bail granularity of 32 bases beat both 16 and 64; see §6.1.
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

/// Population count of a 16-bit mask, the portable fold.
///
/// Always compiled, on every platform, even where `popcount16` below does not call it --
/// otherwise the only build that exercises it is the only build that cannot test it.
/// `tests/test_merge.py::TestPopcount` checks it against Python's own bit count over all
/// 65,536 inputs, from whichever platform the suite happens to run on.
inline int popcount16_portable(unsigned x) noexcept {
    x &= 0xFFFFu;
    x -= (x >> 1) & 0x5555u;                     // pairs
    x = (x & 0x3333u) + ((x >> 2) & 0x3333u);    // nibbles
    x = (x + (x >> 4)) & 0x0F0Fu;                // bytes
    return static_cast<int>((x * 0x0101u) >> 8) & 0xFFu;   // sum the two bytes
}

/// Population count of a 16-bit mask, with **no ISA requirement**.
///
/// `__builtin_popcount` is a GCC/Clang extension. MSVC does not have it, and the first
/// Windows build of this extension failed on exactly that -- C3861, at the `neq16` line
/// below -- because `zna merge` is new in 0.4.0 and no MSVC had ever compiled this file.
///
/// MSVC's `__popcnt16` is **not** the fix. It compiles to the POPCNT instruction, which
/// is SSE4.2-era and not baseline x86-64, so the wheel would build cleanly and then
/// fault with an illegal instruction on an older CPU. A compile error on a machine we
/// control beats a crash on a machine we do not.
inline int popcount16(unsigned x) noexcept {
#if defined(__GNUC__) || defined(__clang__)
    // Cheap even without -mpopcnt: the compiler emits a table or a SWAR sequence, and
    // only emits POPCNT itself when the ISA is enabled explicitly.
    return __builtin_popcount(x & 0xFFFFu);
#else
    return popcount16_portable(x);
#endif
}

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
    return 16 - popcount16(static_cast<unsigned>(mask));
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

inline uint8_t complement_base(uint8_t c) noexcept { return COMPLEMENT.t[c]; }

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
    /// What to do with a no-call the overlap could not rescue from the mate.
    /// 0 = keep it (internal/testing), 1 = trim3, 2 = random. Same vocabulary as
    /// `zna encode --npolicy`, deliberately: one flag, one meaning, both tools.
    int npolicy = 1;
    uint64_t rng_seed = 0;       ///< for NPOLICY_RANDOM; see zna_sub_base
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
constexpr int PROV_TRIMMED   = 1;   ///< from a pair whose redundant overlap was split
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
    std::vector<uint8_t> s2rc, q2r, s1b, q1b, s2b, q2b, seq, qual, name, name2;
    size_t cap = 0;

    void ensure(size_t n) {
        if (n <= cap) return;                 // one predictable, near-never-taken branch
        size_t c = cap ? cap : 1024;
        while (c < n) c <<= 1;
        s2rc.resize(c); q2r.resize(c); s1b.resize(c); q1b.resize(c);
        s2b.resize(c); q2b.resize(c);            // R2 gets consensus written into it too
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
    /// A merged pair emits one record and needs one buffer; a kept or trimmed pair emits
    /// two, and both mates can carry provenance of their own. Sizing this from R2's
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
    int64_t score_q;
    int overlap_len;
    int mismatches;
    int bases_consensus_changed;
    int trim_guard_fired;
    /// Bases the N policy removed (trim3) or invented (random). Reported so a library
    /// that loses a lot of sequence to no-calls says so, instead of finishing quietly.
    int npolicy_bases;
    /// No-calls the overlap recovered from the mate, which cost nothing.
    int n_rescued;
    /// The same two quantities split per MATE, for the per-record header tokens.
    ///
    /// `npolicy_bases == npolicy_1 + npolicy_2` and `n_rescued == n_rescued_1 +
    /// n_rescued_2` always, so the run-level counters keep their exact previous values
    /// and the 15-field counter tuple does not grow. That is deliberate: `_fold` in
    /// merge/cli.py sums a fixed prefix, and a counter added past it reports zero.
    int npolicy_1, npolicy_2;
    int n_rescued_1, n_rescued_2;
    /// Provenance bits per emitted record, parallel to `recs`. See PROV_* above.
    int prov[2];
};

/// Resolve overlap disagreements into R1 alone, by posterior. Returns bases *changed*.
///
/// Used for a MERGED pair, where R1's copy of the overlap is the one emitted -- R2
/// contributes only *outside* it, so its copy is discarded and there is nothing to
/// correct. Every other case that corrects at all corrects both mates
/// (`consensus_pair`); see `process_pair` for the rule.
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

/// Resolve overlap disagreements into **both** mates, by posterior. Returns bases
/// *changed* in R1 (the historical counter; the same decisions are written to R2).
///
/// Used **only** on the trim path, where the pair keeps roughly half the overlap in each
/// read, so R2's copy does reach the corpus. Writing only R1 -- correct while R1 kept
/// the whole overlap -- would emit R2's *uncorrected* bases for its half, giving back
/// part of the consensus exactly where the mates disagreed.
///
/// Restricting it to that branch is not an optimisation. Written on every detected
/// overlap it rewrites bases around a fragment boundary on an inference too weak to
/// merge on -- measured, 237 corrupted 5. ends per million pairs.
///
/// Note `shift >= 0` is NOT what makes this safe, though it was once claimed to be. The
/// real bound is  min written R2 index = max(0, len2 - len1 + shift),  which for
/// unequal-length mates reaches index 0 at non-negative shifts: 9 of those 237 are at
/// `shift >= 0`. What keeps the TRIM branch clear of R2. s 5. end is `trim_is_allowed`,
/// which forces `shift >= min_read_length`.
///
/// `s2` and `q2` are R2 in its OWN orientation, which is how it is emitted; overlap
/// slot `ib` of the reverse-complemented axis is R2 index `len2 - 1 - ib`, and the base
/// stored there is the complement of the resolved call.
inline int consensus_pair(uint8_t* s1, uint8_t* q1, uint8_t* s2rc, uint8_t* q2r,
                          uint8_t* s2, uint8_t* q2, int len2, int s, int olen,
                          const uint8_t* disagree_q,
                          int* rescued1, int* rescued2) noexcept {
    const int a0 = s > 0 ? s : 0;      // mirrors the scan's overlap alignment
    const int b0 = s < 0 ? -s : 0;
    int changed = 0;
    for (int i = 0; i < olen; ++i) {
        const int ia = a0 + i, ib = b0 + i;
        if (s1[ia] != s2rc[ib]) {
            const uint8_t qa = q1[ia], qb = q2r[ib];
            const int j2 = len2 - 1 - ib;                  // same base, R2's frame
            const bool a_is_n = s1[ia] == 'N', b_is_n = s2rc[ib] == 'N';
            if (a_is_n != b_is_n) {             // rescue: a real call beats an N
                // Charged to the mate that was REPAIRED, so each emitted record's
                // `rescued_<n>` token counts only its own recovered no-calls.
                if (a_is_n) {                   // R2 rescues R1
                    ++*rescued1;
                    s1[ia] = s2rc[ib];
                    q1[ia] = qb;
                    ++changed;                  // `changed` counts R1 only, by contract
                } else {                        // R1 rescues R2
                    ++*rescued2;
                    s2rc[ib] = s1[ia];
                    q2r[ib] = qa;
                    s2[j2] = complement_base(s1[ia]);
                    q2[j2] = qa;
                }
            } else if (a_is_n) {                // both N: nothing to rescue from
                continue;
            } else if (qb > qa) {               // R2 is the better-supported call
                const uint8_t nq = disagree_q[(size_t)qb * 256 + qa];
                s1[ia] = s2rc[ib];
                q1[ia] = nq;
                q2r[ib] = nq;
                q2[j2] = nq;
                ++changed;
            } else {                            // R1 wins, but contested: derate it
                const uint8_t nq = disagree_q[(size_t)qa * 256 + qb];
                q1[ia] = nq;
                s2rc[ib] = s1[ia];
                q2r[ib] = nq;
                s2[j2] = complement_base(s1[ia]);
                q2[j2] = nq;
            }
        }
    }
    return changed;
}

/// Emitted lengths for a trimmed pair: as close to equal as the geometry allows.
///
/// The pair must tile the fragment exactly once, so `keep1 + keep2 == L` is forced and
/// the only freedom is where the cut falls. Splitting the overlap down the middle
/// (rather than taking all of it off R2) keeps the two emitted reads the same length,
/// which is what downstream aligners and models expect, and it discards the *last*
/// cycles of both reads -- the lowest-quality bases in the pair -- instead of one
/// read's entire copy.
///
/// `keep1` is clamped into `[L - len2, len1]`, the range in which both reads can supply
/// their share; for equal-length mates that clamp never binds and this is exactly
/// "cut `olen / 2` from each".
inline void balanced_split(int L, int len1, int len2, int& keep1, int& keep2) noexcept {
    int k = (L + 1) / 2;                       // odd overlap: the extra base stays on R1
    const int lo = L - len2;
    if (k < lo) k = lo;
    if (k > len1) k = len1;
    keep1 = k;
    keep2 = L - k;
}

/// May this pair be trimmed at all? Each mate must reach at least `lr` bases past the
/// other's 3' end.
///
/// Two things at once, and both matter:
///
///  1. **The guard.** A trim must never turn a usable pair into a dropped fragment.
///     `balanced_split` clamps `keep1` into `[L - len2, len1]`, so `keep1 >= L - len2`
///     and `keep2 >= L - len1`; requiring both of those to be at least `lr` therefore
///     puts both emitted reads at or above the length filter, by construction.
///
///  2. **A ceiling on the overlap the trim band may act on.** For equal-length mates
///     this is exactly the old rule `len2 - olen >= lr`, and it was doing more work than
///     its name suggested: it refuses a trim whose inferred overlap covers nearly the
///     whole read. Such an overlap is spurious almost by definition -- 145 clean bases
///     score 288 bits and would have merged, so reaching the 8-28 bit trim band at that
///     length takes a pile of mismatches. Dropping the rule when the trim became
///     symmetric admitted exactly those pairs and cost 133 extra false trims, 17,214
///     deleted bases and 9 corrupted 5' ends per million pairs -- all of them fragments
///     with no true overlap at all.
inline bool trim_is_allowed(int L, int len1, int len2, int lr) noexcept {
    return (L - len1) >= lr && (L - len2) >= lr;
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

/// Classify one pair and build its output records.
///
/// Every construction path takes from the 5' end, and trimming only ever cuts 3' ends,
/// so base 0 of every emitted read stays a true fragment boundary. A merged record is
/// built from the inferred span `L = s + len2`, so its length is `L` identically for
/// every geometry -- do not reintroduce per-direction case analysis.
inline PairResult process_pair(const Read& r1, const Read& r2,
                               const Params& p, Scratch& sc,
                               int64_t pair_index = 0) {
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

    const int lr = p.min_read_length;
    const bool detected = r.overlap_len > 0;
    const int L = r.shift + len2;          // inferred fragment length (meaningless if !detected)

    // PROVISIONAL decision, on the full reads. It settles where the consensus is written,
    // and it carries the evidence forward: the score is a statement about this pair's
    // FRAGMENT LENGTH, and trimming interior bases cannot change a fragment's length. The
    // final decision is re-taken after trim3 below -- on GEOMETRY, never by re-scoring.
    // See docs/NPOLICY_PLAN.md D4a.
    const bool prov_merge = detected && r.score_q >= p.t_merge_q;
    const bool prov_band = detected && r.shift >= 0 &&
                           r.score_q >= p.t_trim_q && r.score_q < p.t_merge_q;
    const bool prov_trim = prov_band && trim_is_allowed(L, len1, len2, lr);

    // WHERE the consensus is written: into the records whose CONSTRUCTION depends on the
    // overlap being real, and nowhere else.
    //
    //   merged  -> R1 alone. R1's overlap region becomes the merged record; R2
    //              contributes only outside it, so its copy is discarded.
    //   trimmed -> both. Each mate keeps part of the overlap, so both copies are
    //              emitted and both must carry the same call.
    //   kept    -> neither. Nothing emitted depends on the alignment being right.
    //
    // The kept case is not a symmetry nicety, it is what the measurements say. A
    // detection that lands in KEPT is spurious almost by construction, and for two
    // independent reasons:
    //
    //   * shift >= 0 lands here only because `trim_is_allowed` refused it, which needs
    //     shift < min_read_length and hence an inferred overlap over ~110 bases. A
    //     genuine 110-base overlap scores ~218 bits and would have merged at 28. That is
    //     the trim guard's own argument, and it applies verbatim here: an overlap too
    //     suspect to CUT on is too suspect to REWRITE BASES on.
    //   * shift < 0 (read-through) needs the overlap to equal the fragment length, so
    //     scoring in [8, 28) requires a 5-14 bp fragment.
    //
    // Measured on 1M ground-truth pairs: of 3,068 kept pairs carrying a detected
    // overlap, ZERO found the true shift and 97.3% had no true overlap at all. Wrong
    // emitted bases in the overlap window under the three candidate rules --
    // correct-neither 208, correct-R1-only (the old behaviour) 1,509, correct-both
    // 17,870. The old R1 write turned 1,379 correct bases wrong to fix 78.
    //
    // Note `shift >= 0` is NOT a safe proxy for "the write stays off R2's 5' end". The
    // real bound is  min written R2 index = max(0, len2 - len1 + shift),  which for
    // unequal-length mates reaches index 0 at non-negative shifts: 9 of the recorded 237
    // corrupted 5' ends are at shift >= 0.
    const bool write_r1 = prov_merge || prov_trim;
    const bool write_r2 = prov_trim;

    // Copy-on-write: a clean winning overlap (56.5% of real pairs) needs no mutable
    // copy of either read at all, and only a trim needs R2 copied.
    const uint8_t* S1 = r1.s.p;
    const uint8_t* Q1 = r1.q.p;
    const uint8_t* S2 = r2.s.p;
    const uint8_t* Q2 = r2.q.p;
    if (r.mismatches > 0 && write_r1) {
        uint8_t* q2r = sc.q2r.data();
        for (int i = 0; i < len2; ++i) q2r[i] = r2.q.p[len2 - 1 - i];
        uint8_t* s1b = sc.s1b.data();
        uint8_t* q1b = sc.q1b.data();
        std::memcpy(s1b, r1.s.p, static_cast<size_t>(len1));
        std::memcpy(q1b, r1.q.p, static_cast<size_t>(len1));
        if (write_r2) {
            uint8_t* s2b = sc.s2b.data();
            uint8_t* q2b = sc.q2b.data();
            std::memcpy(s2b, r2.s.p, static_cast<size_t>(len2));
            std::memcpy(q2b, r2.q.p, static_cast<size_t>(len2));
            out.bases_consensus_changed =
                consensus_pair(s1b, q1b, s2rc, q2r, s2b, q2b, len2,
                               r.shift, r.overlap_len, p.disagree_q,
                               &out.n_rescued_1, &out.n_rescued_2);
            S2 = s2b;
            Q2 = q2b;
        } else {
            // The merged path rescues into R1 only -- R2's copy of the overlap is
            // discarded -- so every rescue here is charged to mate 1.
            out.bases_consensus_changed =
                consensus_r1(s1b, q1b, s2rc, q2r, r.shift, r.overlap_len, p.disagree_q,
                             &out.n_rescued_1);
        }
        S1 = s1b;
        Q1 = q1b;
    }

    // ---- trim3: cut each read at its first SURVIVING N -----------------------------
    //
    // After the rescue, so a no-call the mate could answer costs nothing. 3' only, so
    // both 5' anchors -- the two fragment termini -- are untouched however short the
    // reads get.
    int elen1 = len1, elen2 = len2;
    int shift = r.shift;
    if (p.npolicy == NPOLICY_RANDOM) {
        // Substitution does not change a length, so the coverage test below is
        // unaffected and `random` never costs a merge -- unlike trim3.
        bool any = false;
        for (int i = 0; i < elen1 && !any; ++i) any = S1[i] == 'N';
        for (int i = 0; i < elen2 && !any; ++i) any = S2[i] == 'N';
        if (any) {
            uint8_t* s1b = sc.s1b.data();
            uint8_t* s2b = sc.s2b.data();
            (void)0;  // S1/S2 may ALREADY be these buffers -- the consensus writes into them --
            // and memcpy with src == dst is undefined. Copy only when they differ.
            if (S1 != s1b) std::memcpy(s1b, S1, static_cast<size_t>(elen1));
            if (S2 != s2b) std::memcpy(s2b, S2, static_cast<size_t>(elen2));
            const uint64_t k1 = static_cast<uint64_t>(pair_index) * 2;
            for (int i = 0; i < elen1; ++i) {
                if (s1b[i] == 'N') { s1b[i] = merge_sub_base(p.rng_seed, k1, i);
                                     ++out.npolicy_1; }
            }
            for (int i = 0; i < elen2; ++i) {
                if (s2b[i] == 'N') { s2b[i] = merge_sub_base(p.rng_seed, k1 + 1, i);
                                     ++out.npolicy_2; }
            }
            S1 = s1b; S2 = s2b;
            revcomp_into(S2, elen2, s2rc);
        }
    } else if (p.npolicy == NPOLICY_TRIM3) {
        int k1 = elen1, k2 = elen2;
        for (int i = 0; i < elen1; ++i) { if (S1[i] == 'N') { k1 = i; break; } }
        for (int i = 0; i < elen2; ++i) { if (S2[i] == 'N') { k2 = i; break; } }
        if (k1 != elen1 || k2 != elen2) {
            out.npolicy_1 = elen1 - k1;
            out.npolicy_2 = elen2 - k2;
            elen1 = k1;
            elen2 = k2;
            // `shift` is the offset of revcomp(R2) on the shared axis, so it is tied to
            // len2. R2 keeps its 5' anchor at fragment position L-1, so the trimmed mate
            // covers [L - elen2, L) and the offset becomes L - elen2. L is unchanged --
            // that is the whole point.
            shift = L - elen2;
            revcomp_into(S2, elen2, s2rc);       // s2rc must match the trimmed mate
        }
    }

    // ---- the final decision, on GEOMETRY, reusing the original evidence -------------
    //
    // The pair still tiles the fragment iff elen1 + elen2 >= L, and when it does the
    // reconstruction IS the fragment, exactly and N-free. Nothing is re-scored.
    const bool covers = detected && (elen1 + elen2) >= L;
    const bool will_merge = prov_merge && covers;
    const bool will_trim = !will_merge && prov_band && (elen1 + elen2) > L &&
                           trim_is_allowed(L, elen1, elen2, lr);

    // The run-level counters are the per-mate ones summed, so they keep their exact
    // previous values and the counter tuple does not grow. See PairResult.
    out.npolicy_bases = out.npolicy_1 + out.npolicy_2;
    out.n_rescued = out.n_rescued_1 + out.n_rescued_2;

    // Which policy bit a touched record earns. Set per RECORD, from that record's own
    // count -- a kept pair whose R1 lost bases and whose R2 did not says exactly that.
    const int npolicy_bit = (p.npolicy == NPOLICY_RANDOM) ? PROV_NSUBBED : PROV_NTRIMMED;

    bool paired;
    int n_cand;

    if (will_merge) {
        // Enough evidence to assert the fragment: collapse to one read.
        const int s = shift;
        const int take1 = elen1 < L ? elen1 : L;
        const int take2 = L - take1;
        uint8_t* seq = sc.seq.data();
        uint8_t* qual = sc.qual.data();
        std::memcpy(seq, S1, static_cast<size_t>(take1));
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
    } else if (will_trim) {
        // Real overlap, but not enough evidence to risk a chimera: keep both reads and
        // split the redundant overlap between them, so the fragment is tiled exactly
        // once and the two emitted reads stay the same length.
        int keep1, keep2;
        balanced_split(L, elen1, elen2, keep1, keep2);
        out.prov[0] = PROV_TRIMMED | (out.n_rescued_1 ? PROV_RESCUED : 0)
                                   | (out.npolicy_1 ? npolicy_bit : 0);
        out.prov[1] = PROV_TRIMMED | (out.n_rescued_2 ? PROV_RESCUED : 0)
                                   | (out.npolicy_2 ? npolicy_bit : 0);
        out.recs[0] = {name_for(sc, 0, r1.h, out.prov[0], p.npolicy,
                                out.npolicy_1, out.n_rescued_1),
                       {S1, keep1}, {Q1, keep1}};
        out.recs[1] = {name_for(sc, 1, r2.h, out.prov[1], p.npolicy,
                                out.npolicy_2, out.n_rescued_2),
                       {S2, keep2}, {Q2, keep2}};
        n_cand = 2;
        paired = true;
        out.outcome = OUTCOME_TRIMMED;
    } else {
        // No detectable overlap, an unmergeable read-through, or a trim blocked by the
        // guard: keep both reads exactly as they are.
        //
        // No consensus is written on this path, so a kept record can never carry
        // PROV_RESCUED -- only the N policy can have touched it.
        if (prov_band && !prov_trim) {
            out.trim_guard_fired = 1;
        }
        out.prov[0] = out.npolicy_1 ? npolicy_bit : 0;
        out.prov[1] = out.npolicy_2 ? npolicy_bit : 0;
        out.recs[0] = {name_for(sc, 0, r1.h, out.prov[0], p.npolicy, out.npolicy_1, 0),
                       {S1, elen1}, {Q1, elen1}};
        out.recs[1] = {name_for(sc, 1, r2.h, out.prov[1], p.npolicy, out.npolicy_2, 0),
                       {S2, elen2}, {Q2, elen2}};
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
