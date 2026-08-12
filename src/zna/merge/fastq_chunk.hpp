/**
 * FASTQ chunk adapter: raw text in, formatted text out.
 *
 * This is an *adapter* over `merge_core.hpp`, not part of the kernel
 * (docs/MERGE_CPP_DESIGN.md §7.2). The core knows nothing about FASTQ; this file knows
 * nothing about the scoring rule. `zna encode --merge-pairs` will add a second adapter
 * beside this one that emits records instead of text, reusing the same core.
 *
 * It is a pure function of its arguments -- no globals, no clock, no RNG -- which is
 * what makes it fuzzable, safely callable with the GIL released, and byte-deterministic.
 *
 * **The consumption protocol.** `merge_chunk` consumes only WHOLE PAIRS and reports how
 * many bytes it took from each stream separately, so the two buffers may carry different
 * leftovers and the caller never has to scan for record boundaries. A partial record at
 * the end of a buffer is not an error -- it is simply not consumed. At EOF the caller
 * asserts BOTH buffers are empty; a non-empty leftover on either side is the "unequal
 * read counts" error.
 *
 * That protocol exists because of a specific defect. The audit's prototype for this
 * consumed `min(records1, records2)` pairs and returned success, SILENTLY DROPPING the
 * trailing R2 records when R1 ran out first -- which the desync check cannot catch,
 * because comparing base names cannot see records that were never read. Here the drop is
 * structurally impossible: unconsumed bytes remain in the caller's buffer and the caller
 * checks them.
 */
#ifndef ZNA_MERGE_FASTQ_CHUNK_HPP
#define ZNA_MERGE_FASTQ_CHUNK_HPP

#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>

#include "merge_core.hpp"

namespace zna_merge {

/// Malformed or inconsistent input: truncation, desync, unequal read counts.
/// Distinct from a bug so the binding can raise `zna.merge.fastqio.InputError` and the
/// CLI can turn it into a clean message rather than a traceback.
struct InputError : std::runtime_error {
    using std::runtime_error::runtime_error;
};

constexpr int HIST_MAX = 1024;   ///< histograms are clamped to this many bins

struct ChunkStats {
    int64_t n_pairs = 0, merged = 0, trimmed = 0, kept = 0;
    int64_t emitted = 0, dropped = 0, bases_trimmed = 0, frags_short = 0;
    int64_t bases_consensus = 0, trim_guard = 0, sum_olen = 0, sum_diff = 0;
    uint32_t len_hist[HIST_MAX + 1] = {};      ///< every emitted record's length
    uint32_t olen_hist[HIST_MAX + 1] = {};     ///< detected overlap lengths
    uint32_t insert_hist[HIST_MAX + 1] = {};   ///< inferred fragment length, merged only
};

/// Scratch for the two uppercased sequences, beside the core's per-pair arena.
struct ChunkScratch {
    Scratch core;
    std::vector<uint8_t> up1, up2;
    void ensure_reads(size_t n) {
        core.ensure(n);
        if (up1.size() < core.cap) { up1.resize(core.cap); up2.resize(core.cap); }
    }
};

/// ASCII upper-case, matching Python's ``bytes.upper()``: a-z only, everything else
/// passes through. The reader has always upper-cased sequences so that comparison is
/// case-insensitive; quality strings are left alone.
inline void upper_into(const uint8_t* src, int n, uint8_t* dst) noexcept {
    for (int i = 0; i < n; ++i) {
        const uint8_t c = src[i];
        dst[i] = (c >= 'a' && c <= 'z') ? static_cast<uint8_t>(c - 32) : c;
    }
}

inline int strip_cr(const uint8_t* p, int n) noexcept {
    return (n > 0 && p[n - 1] == '\r') ? n - 1 : n;
}

/// Read ID up to the first whitespace, minus any ``/1``/``/2`` suffix.
/// Mirrors ZNA's own pairing rule, so the sync check agrees with how ZNA pairs records.
inline Span base_name(const Span& h) noexcept {
    int cut = h.n;
    for (int i = 0; i < h.n; ++i) {
        if (h.p[i] == ' ' || h.p[i] == '\t') { cut = i; break; }
    }
    for (int i = cut - 1; i >= 0; --i) {
        if (h.p[i] == '/') return {h.p, i};
    }
    return {h.p, cut};
}

/// Parse one record starting at *pos*, advancing it. Returns false when no COMPLETE
/// record remains (which is not an error -- the caller refills and retries).
///
/// *up* receives the upper-cased sequence, so the returned Read is valid until the next
/// call that shares the same buffer.
inline bool next_record(const uint8_t* buf, size_t n, size_t& pos,
                        uint8_t* up, size_t up_cap, Read& out, const char* which) {
    if (pos >= n) return false;
    const uint8_t* base = buf + pos;
    const size_t avail = n - pos;

    const uint8_t* e1 = static_cast<const uint8_t*>(std::memchr(base, '\n', avail));
    if (!e1) return false;
    const uint8_t* e2 = static_cast<const uint8_t*>(
        std::memchr(e1 + 1, '\n', static_cast<size_t>(base + avail - (e1 + 1))));
    if (!e2) return false;
    const uint8_t* e3 = static_cast<const uint8_t*>(
        std::memchr(e2 + 1, '\n', static_cast<size_t>(base + avail - (e2 + 1))));
    if (!e3) return false;
    const uint8_t* e4 = static_cast<const uint8_t*>(
        std::memchr(e3 + 1, '\n', static_cast<size_t>(base + avail - (e3 + 1))));
    if (!e4) return false;

    if (base[0] != '@') {
        throw InputError(std::string("malformed FASTQ header in ") + which);
    }
    const int hl = strip_cr(base + 1, static_cast<int>(e1 - base) - 1);
    const int sl = strip_cr(e1 + 1, static_cast<int>(e2 - e1) - 1);
    const int ql = strip_cr(e3 + 1, static_cast<int>(e4 - e3) - 1);
    if (sl != ql) {
        // A file truncated inside its LAST quality line otherwise looks like a complete
        // record and is emitted malformed with a zero exit status.
        throw InputError(std::string("FASTQ record in ") + which + " has " +
                         std::to_string(sl) + " bases but " + std::to_string(ql) +
                         " quality scores (truncated or malformed)");
    }
    if (static_cast<size_t>(sl) > up_cap) {
        throw InputError("read longer than the scratch buffer");   // caller sizes it
    }
    upper_into(e1 + 1, sl, up);
    out.h = {base + 1, hl};
    out.s = {up, sl};
    out.q = {e3 + 1, ql};
    pos += static_cast<size_t>(e4 - base) + 1;
    return true;
}

/// Byte offset just past *max_records* complete records, and how many were found.
///
/// Lets the caller cut a buffer into whole-record chunks for parallel workers without
/// parsing anything itself. Four newlines per record; a trailing partial record is not
/// counted and its bytes are left for the next chunk.
inline void split_records(const uint8_t* buf, size_t n, int64_t max_records,
                          size_t& out_offset, int64_t& out_count) noexcept {
    size_t pos = 0;
    int64_t found = 0;
    while (max_records <= 0 || found < max_records) {
        size_t p = pos;
        bool complete = true;
        for (int line = 0; line < 4; ++line) {
            const uint8_t* nl = static_cast<const uint8_t*>(
                std::memchr(buf + p, '\n', n - p));
            if (!nl) { complete = false; break; }
            p = static_cast<size_t>(nl - buf) + 1;
        }
        if (!complete) break;
        pos = p;
        ++found;
    }
    out_offset = pos;
    out_count = found;
}

inline void emit(std::string& blob, const OutRec& r) {
    blob += '@';
    blob.append(reinterpret_cast<const char*>(r.h.p), static_cast<size_t>(r.h.n));
    blob += '\n';
    blob.append(reinterpret_cast<const char*>(r.s.p), static_cast<size_t>(r.s.n));
    blob += "\n+\n";
    blob.append(reinterpret_cast<const char*>(r.q.p), static_cast<size_t>(r.q.n));
    blob += '\n';
}

/// Merge every whole pair available in both buffers, appending FASTQ text to *blob*.
///
/// *pos1* / *pos2* are advanced to the number of bytes consumed from each stream.
/// *base_index* is the index of the first pair in the input, carried only so a desync
/// can be reported by absolute pair number rather than "somewhere in this chunk".
inline void merge_chunk(const uint8_t* buf1, size_t n1, size_t& pos1,
                        const uint8_t* buf2, size_t n2, size_t& pos2,
                        const Params& p, bool check_sync, int64_t base_index,
                        ChunkScratch& sc, std::string& blob, ChunkStats& st) {
    Read r1, r2;
    for (;;) {
        // Reserve room before parsing, so `next_record` can bound-check against it.
        sc.ensure_reads(1024);
        size_t try1 = pos1, try2 = pos2;
        if (!next_record(buf1, n1, try1, sc.up1.data(), sc.up1.size(), r1, "R1")) break;
        if (!next_record(buf2, n2, try2, sc.up2.data(), sc.up2.size(), r2, "R2")) break;

        // A read longer than the arena needs a bigger arena and a re-parse; rare enough
        // to be worth the retry rather than a pre-scan of every record.
        if (r1.s.n > static_cast<int>(sc.core.cap) ||
            r2.s.n > static_cast<int>(sc.core.cap)) {
            sc.ensure_reads(static_cast<size_t>(r1.s.n > r2.s.n ? r1.s.n : r2.s.n));
            try1 = pos1; try2 = pos2;
            next_record(buf1, n1, try1, sc.up1.data(), sc.up1.size(), r1, "R1");
            next_record(buf2, n2, try2, sc.up2.data(), sc.up2.size(), r2, "R2");
        }

        if (check_sync) {
            const Span a = base_name(r1.h), b = base_name(r2.h);
            if (a.n != b.n || std::memcmp(a.p, b.p, static_cast<size_t>(a.n)) != 0) {
                throw InputError(
                    "R1/R2 out of sync at pair " +
                    std::to_string(base_index + st.n_pairs + 1) + ": '" +
                    std::string(reinterpret_cast<const char*>(a.p), static_cast<size_t>(a.n)) +
                    "' != '" +
                    std::string(reinterpret_cast<const char*>(b.p), static_cast<size_t>(b.n)) + "'");
            }
        }

        const PairResult r = process_pair(r1, r2, p, sc.core);
        ++st.n_pairs;
        st.dropped += r.n_dropped;
        st.bases_consensus += r.bases_consensus_changed;
        st.trim_guard += r.trim_guard_fired;

        if (r.outcome == OUTCOME_MERGED) {
            ++st.merged;
        } else if (r.outcome == OUTCOME_TRIMMED) {
            ++st.trimmed;
            if (r.n_recs) st.bases_trimmed += r2.s.n - r.recs[1].s.n;
        } else {
            ++st.kept;
        }
        if (!r.n_recs && r.outcome != OUTCOME_MERGED) ++st.frags_short;

        if (r.overlap_len) {
            st.olen_hist[r.overlap_len <= HIST_MAX ? r.overlap_len : HIST_MAX] += 1;
            st.sum_olen += r.overlap_len;
            st.sum_diff += r.mismatches;
        }
        for (int i = 0; i < r.n_recs; ++i) {
            emit(blob, r.recs[i]);
            ++st.emitted;
            const int L = r.recs[i].s.n;
            st.len_hist[L <= HIST_MAX ? L : HIST_MAX] += 1;
            if (r.outcome == OUTCOME_MERGED) {
                // A merged record IS the fragment, so its length is the insert size.
                st.insert_hist[L <= HIST_MAX ? L : HIST_MAX] += 1;
            }
        }
        pos1 = try1;
        pos2 = try2;
    }
}

}  // namespace zna_merge

#endif  // ZNA_MERGE_FASTQ_CHUNK_HPP
