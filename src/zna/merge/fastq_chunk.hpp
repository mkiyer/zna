/**
 * FASTQ chunk adapter: raw text in, formatted text out.
 *
 * This is an *adapter* over `merge_core.hpp`, not part of the kernel
 * (docs/METHODS.md, "Layering"). The core knows nothing about FASTQ; this file knows
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
#include <vector>

#include "merge_core.hpp"

namespace zna_merge {

/// Malformed or inconsistent input: truncation, desync, unequal read counts.
/// Distinct from a bug so the binding can raise `zna.merge.fastqio.InputError` and the
/// CLI can turn it into a clean message rather than a traceback.
struct InputError : std::runtime_error {
    using std::runtime_error::runtime_error;
};

struct ChunkStats {
    int64_t n_pairs = 0, merged = 0, trimmed = 0, kept = 0;
    int64_t emitted = 0, dropped = 0, bases_trimmed = 0, frags_short = 0;
    int64_t bases_consensus = 0, trim_guard = 0, sum_olen = 0, sum_diff = 0;
    /// Bases the N policy removed (trim3) or invented (random), and no-calls the
    /// overlap recovered from the mate at no cost.
    int64_t npolicy_bases = 0, n_rescued = 0;
    /// Longest read seen. Reported rather than capped: the scan is O(L^2), so this is
    /// how an accidental long-read FASTQ becomes diagnosable instead of just slow.
    int max_read_len = 0;
    std::vector<uint32_t> len_hist;      ///< every emitted record's length
    std::vector<uint32_t> olen_hist;     ///< detected overlap lengths
    std::vector<uint32_t> insert_hist;   ///< inferred fragment length, merged only

    /// Size the three histograms from the scratch arena's capacity, so indexing them
    /// needs no bound check and no clamp.
    ///
    /// They used to be `uint32_t[1025]` with every index clamped to the last bin, which
    /// silently aggregated the length and insert distributions for reads over 1024 bp --
    /// and read length is otherwise uncapped. The arena grows by doubling and every
    /// value binned here is bounded by it: an overlap is at most `cap`, and an emitted
    /// record is at most `len1 + len2 <= 2 * cap`. So `2 * cap + 1` bins always suffice
    /// and this resizes O(log L) times per chunk, not per record.
    ///
    /// Not a sparse map: that would put a hash lookup where a dense index is one
    /// instruction, in the per-record path.
    void ensure_bins(size_t cap) {
        const size_t want = 2 * cap + 1;
        if (len_hist.size() >= want) return;
        len_hist.resize(want, 0);
        olen_hist.resize(want, 0);
        insert_hist.resize(want, 0);
    }
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

/// A located record: spans into the input buffer, nothing copied yet.
struct RecordSpans {
    Span h, s, q;
    size_t next_pos;
};

/// Locate one record at *pos* and validate it, without copying anything.
///
/// Returns false when no COMPLETE record remains, which is not an error -- the caller
/// refills and retries. Separated from the upper-casing below so the caller can size
/// its scratch from the record's ACTUAL length: doing both at once meant a read longer
/// than the current arena threw instead of growing it, which made the compiled backend
/// reject reads over 1024 bp while the reference merged them happily.
inline bool locate_record(const uint8_t* buf, size_t n, size_t pos,
                          RecordSpans& out, const char* which) {
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
    out.h = {base + 1, hl};
    out.s = {e1 + 1, sl};                 // still pointing at the input, not upper-cased
    out.q = {e3 + 1, ql};
    out.next_pos = pos + static_cast<size_t>(e4 - base) + 1;
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

/// Which mate of the pair an emitted record came from -- the whole geometry
/// transfer of `zna encode --merge-pairs` (MERGE_PAIRS_PLAN.md §0.1).  Slot,
/// not ZNA flag bits, crosses this boundary: the flag layout stays defined in
/// exactly one place, on the Python side.
enum Slot : uint8_t { SLOT_MERGED = 0, SLOT_MATE1 = 1, SLOT_MATE2 = 2 };

/// One record emitted by `merge_chunk_records`.
///
/// `seq_off`/`seq_len` index into the caller's `seqs` blob.  `hdr_off`/`hdr_len`
/// index into the ORIGINAL input buffer -- buf1 for MERGED and MATE1, buf2 for
/// MATE2, derivable from `slot` -- and `hdr_len` is 0 when headers were not
/// requested.  Pointing at the input rather than at `OutRec::h` is not merely
/// cheaper, it is the only CORRECT option for a merged record: its `OutRec::h`
/// points at `Scratch::name`, a single per-pair buffer the next pair
/// overwrites.  It also hands `zna encode --label` byte-identical tag values to
/// the two-step path, whose extractor reads the same input header.
///
/// `prov` is the record's PROV_* byte straight off `PairResult::prov` -- the
/// direct transfer that lets `--merge-pairs` write the provenance column with
/// no `ZN:i:` tag round-trip (MERGE_PAIRS_PLAN.md §10.4).
struct RecordEnd {
    uint32_t seq_off, seq_len;
    uint32_t hdr_off, hdr_len;
    uint8_t slot;
    uint8_t prov;
};

/// The one inner loop both adapters share: locate, sync-check, merge, tally.
///
/// The loop carries the consumption protocol, the desync check and every
/// statistic; the only thing the two output shapes disagree on is what to do
/// with an emitted record, so that is the only thing the emitter functor is
/// handed.  Two copies of this loop would drift -- the reason this is a
/// template, not a convention.
///
/// `emit(r, i, a, b)` is called once per emitted record: the `PairResult`, the
/// record's index within it, and the located input records for both mates.
template <class EmitFn>
inline void merge_chunk_impl(const uint8_t* buf1, size_t n1, size_t& pos1,
                             const uint8_t* buf2, size_t n2, size_t& pos2,
                             const Params& p, bool check_sync, int64_t base_index,
                             ChunkScratch& sc, ChunkStats& st, EmitFn&& emit_rec) {
    RecordSpans a, b;
    Read r1, r2;
    for (;;) {
        if (!locate_record(buf1, n1, pos1, a, "R1")) break;
        if (!locate_record(buf2, n2, pos2, b, "R2")) break;

        // Size the arena from the records we actually have, then copy into it. No cap
        // and no flag: the arena doubles to whatever the input turns out to need.
        const int longest = a.s.n > b.s.n ? a.s.n : b.s.n;
        sc.ensure_reads(static_cast<size_t>(longest));
        st.ensure_bins(sc.core.cap);          // every bin index below is bounded by cap
        if (longest > st.max_read_len) st.max_read_len = longest;

        upper_into(a.s.p, a.s.n, sc.up1.data());
        upper_into(b.s.p, b.s.n, sc.up2.data());
        r1 = {a.h, {sc.up1.data(), a.s.n}, a.q};
        r2 = {b.h, {sc.up2.data(), b.s.n}, b.q};

        if (check_sync) {
            const Span an = base_name(r1.h), bn = base_name(r2.h);
            if (an.n != bn.n || std::memcmp(an.p, bn.p, static_cast<size_t>(an.n)) != 0) {
                throw InputError(
                    "R1/R2 out of sync at pair " +
                    std::to_string(base_index + st.n_pairs + 1) + ": '" +
                    std::string(reinterpret_cast<const char*>(an.p), static_cast<size_t>(an.n)) +
                    "' != '" +
                    std::string(reinterpret_cast<const char*>(bn.p), static_cast<size_t>(bn.n)) + "'");
            }
        }

        const PairResult r = process_pair(r1, r2, p, sc.core,
                                          base_index + st.n_pairs);
        ++st.n_pairs;
        st.dropped += r.n_dropped;
        st.bases_consensus += r.bases_consensus_changed;
        st.trim_guard += r.trim_guard_fired;
        st.npolicy_bases += r.npolicy_bases;
        st.n_rescued += r.n_rescued;

        if (r.outcome == OUTCOME_MERGED) {
            ++st.merged;
        } else if (r.outcome == OUTCOME_TRIMMED) {
            ++st.trimmed;
            // Both mates are cut now, so charge both: this counter is "redundant bases
            // removed", and it is the overlap length either way.
            if (r.n_recs) {
                st.bases_trimmed += (r1.s.n - r.recs[0].s.n) + (r2.s.n - r.recs[1].s.n);
            }
        } else {
            ++st.kept;
        }
        if (!r.n_recs && r.outcome != OUTCOME_MERGED) ++st.frags_short;

        if (r.overlap_len) {
            st.olen_hist[r.overlap_len] += 1;
            st.sum_olen += r.overlap_len;
            st.sum_diff += r.mismatches;
        }
        for (int i = 0; i < r.n_recs; ++i) {
            emit_rec(r, i, a, b);
            ++st.emitted;
            const int L = r.recs[i].s.n;
            st.len_hist[L] += 1;
            if (r.outcome == OUTCOME_MERGED) {
                // A merged record IS the fragment, so its length is the insert size.
                st.insert_hist[L] += 1;
            }
        }
        pos1 = a.next_pos;
        pos2 = b.next_pos;
    }
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
    merge_chunk_impl(
        buf1, n1, pos1, buf2, n2, pos2, p, check_sync, base_index, sc, st,
        [&blob](const PairResult& r, int i, const RecordSpans&,
                const RecordSpans&) { emit(blob, r.recs[i]); });
}

/// Merge every whole pair available, appending RECORDS instead of FASTQ text.
///
/// The second adapter over the same core, for `zna encode --merge-pairs`: no
/// quality strings (ZNA does not store them), no name construction consumed
/// (the core still builds merged names internally; v1 accepts that snprintf --
/// see MERGE_PAIRS_PLAN.md §2.3), and headers as offsets into the caller's
/// input buffers, only when *want_headers*.
inline void merge_chunk_records(const uint8_t* buf1, size_t n1, size_t& pos1,
                                const uint8_t* buf2, size_t n2, size_t& pos2,
                                const Params& p, bool check_sync,
                                int64_t base_index, bool want_headers,
                                ChunkScratch& sc, std::string& seqs,
                                std::vector<RecordEnd>& ends, ChunkStats& st) {
    merge_chunk_impl(
        buf1, n1, pos1, buf2, n2, pos2, p, check_sync, base_index, sc, st,
        [&](const PairResult& r, int i, const RecordSpans& a,
            const RecordSpans& b) {
            const OutRec& rec = r.recs[i];
            const uint8_t slot =
                (r.outcome == OUTCOME_MERGED) ? SLOT_MERGED
                : (i == 0 ? SLOT_MATE1 : SLOT_MATE2);
            uint32_t hdr_off = 0, hdr_len = 0;
            if (want_headers) {
                // The record's OWN source header (MERGED reads R1's): MATE2
                // from buf2, everything else from buf1.
                const Span& h = (slot == SLOT_MATE2) ? b.h : a.h;
                const uint8_t* base = (slot == SLOT_MATE2) ? buf2 : buf1;
                hdr_off = static_cast<uint32_t>(h.p - base);
                hdr_len = static_cast<uint32_t>(h.n);
            }
            ends.push_back({static_cast<uint32_t>(seqs.size()),
                            static_cast<uint32_t>(rec.s.n),
                            hdr_off, hdr_len, slot,
                            static_cast<uint8_t>(r.prov[i])});
            seqs.append(reinterpret_cast<const char*>(rec.s.p),
                        static_cast<size_t>(rec.s.n));
        });
}

}  // namespace zna_merge

#endif  // ZNA_MERGE_FASTQ_CHUNK_HPP
