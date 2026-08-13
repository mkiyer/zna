/**
 * nanobind bindings for the accelerated overlap scan.
 *
 * Thin on purpose: everything with an opinion lives in merge_core.hpp, which knows
 * nothing about Python and can therefore be compiled into a sanitizer driver or reused
 * by `zna encode --merge-pairs` without dragging the interpreter along. This file is
 * argument checking and tuple building.
 */

#include <nanobind/nanobind.h>

#include <stdexcept>
#include <vector>

#include "merge_core.hpp"
#include "fastq_chunk.hpp"

namespace nb = nanobind;

/// Best-scoring shift, as (shift, score_q, overlap_len, mismatches).
static nb::tuple scan(nb::bytes seq1, nb::bytes seq2rc, int len1, int len2,
                      int64_t match_q, int64_t step_q, int64_t floor_q) {
    if (len1 < 0 || len2 < 0 ||
        static_cast<size_t>(len1) > seq1.size() ||
        static_cast<size_t>(len2) > seq2rc.size()) {
        throw std::invalid_argument("scan(): length argument exceeds the buffer");
    }
    if (match_q <= 0 || step_q <= 0) {
        throw std::invalid_argument("scan(): match_q and step_q must be positive");
    }
    const zna_merge::ScanResult r = zna_merge::scan(
        reinterpret_cast<const uint8_t*>(seq1.c_str()), len1,
        reinterpret_cast<const uint8_t*>(seq2rc.c_str()), len2,
        match_q, step_q, floor_q);
    return nb::make_tuple(r.shift, r.score_q, r.overlap_len, r.mismatches);
}

namespace {
/// One arena per thread: the per-pair path must not allocate. `process_pair` below is
/// called under the GIL for now, but making this thread_local costs nothing and is what
/// phase E's GIL-releasing chunk path needs.
thread_local zna_merge::Scratch g_scratch;
thread_local zna_merge::ChunkScratch g_chunk_scratch;

inline zna_merge::Span span_of(nb::bytes& b) {
    return {reinterpret_cast<const uint8_t*>(b.c_str()), static_cast<int>(b.size())};
}
inline nb::bytes to_bytes(const zna_merge::Span& s) {
    return nb::bytes(reinterpret_cast<const char*>(s.p), static_cast<size_t>(s.n));
}
}  // namespace

/// Classify one pair and build its records. Mirrors _pymerge.process_pair exactly.
static nb::tuple process_pair(nb::bytes h1, nb::bytes s1, nb::bytes q1,
                              nb::bytes h2, nb::bytes s2, nb::bytes q2,
                              int64_t match_q, int64_t step_q,
                              int64_t t_merge_q, int64_t t_trim_q,
                              int min_read_length, nb::bytes disagree_q,
                              int npolicy, uint64_t rng_seed, int64_t pair_index) {
    if (s1.size() != q1.size() || s2.size() != q2.size()) {
        throw std::invalid_argument("process_pair(): sequence and quality differ in length");
    }
    if (disagree_q.size() != 256u * 256u) {
        throw std::invalid_argument("process_pair(): disagree_q must be 256*256 bytes");
    }
    if (match_q <= 0 || step_q <= 0) {
        throw std::invalid_argument("process_pair(): match_q and step_q must be positive");
    }
    const zna_merge::Read r1{span_of(h1), span_of(s1), span_of(q1)};
    const zna_merge::Read r2{span_of(h2), span_of(s2), span_of(q2)};
    const zna_merge::Params p{match_q, step_q, t_merge_q, t_trim_q, min_read_length,
                              reinterpret_cast<const uint8_t*>(disagree_q.c_str()),
                              npolicy, rng_seed};

    const zna_merge::PairResult r =
        zna_merge::process_pair(r1, r2, p, g_scratch, pair_index);

    nb::list recs;
    for (int i = 0; i < r.n_recs; ++i) {
        recs.append(nb::make_tuple(to_bytes(r.recs[i].h),
                                   to_bytes(r.recs[i].s),
                                   to_bytes(r.recs[i].q)));
    }
    return nb::make_tuple(recs, r.outcome, r.n_dropped, r.score_q, r.overlap_len,
                          r.mismatches, r.bases_consensus_changed, r.trim_guard_fired,
                          r.npolicy_bases, r.n_rescued);
}

/// Merge every whole pair in the two buffers; return the formatted FASTQ text.
///
/// The GIL is released for the whole parse-and-merge, which is what makes threads worth
/// having: the inputs are immutable `bytes` kept alive by the arguments, and nothing
/// here touches a Python object until it is reacquired.
static nb::tuple merge_chunk(nb::bytes buf1, int64_t start1, int64_t end1,
                             nb::bytes buf2, int64_t start2, int64_t end2,
                             int64_t match_q, int64_t step_q,
                             int64_t t_merge_q, int64_t t_trim_q,
                             int min_read_length, nb::bytes disagree_q,
                             bool check_sync, int64_t base_index, int npolicy,
                             uint64_t rng_seed) {
    if (start1 < 0 || end1 < start1 || static_cast<size_t>(end1) > buf1.size() ||
        start2 < 0 || end2 < start2 || static_cast<size_t>(end2) > buf2.size()) {
        throw std::invalid_argument("merge_chunk(): bad [start, end) range");
    }
    if (disagree_q.size() != 256u * 256u) {
        throw std::invalid_argument("merge_chunk(): disagree_q must be 256*256 bytes");
    }
    if (match_q <= 0 || step_q <= 0) {
        throw std::invalid_argument("merge_chunk(): match_q and step_q must be positive");
    }
    const auto* b1 = reinterpret_cast<const uint8_t*>(buf1.c_str()) + start1;
    const auto* b2 = reinterpret_cast<const uint8_t*>(buf2.c_str()) + start2;
    const size_t n1 = static_cast<size_t>(end1 - start1);
    const size_t n2 = static_cast<size_t>(end2 - start2);
    const zna_merge::Params p{match_q, step_q, t_merge_q, t_trim_q, min_read_length,
                              reinterpret_cast<const uint8_t*>(disagree_q.c_str()),
                              npolicy, rng_seed};

    std::string blob;
    zna_merge::ChunkStats st;
    size_t pos1 = 0, pos2 = 0;
    {
        nb::gil_scoped_release release;
        blob.reserve(n1 + n2);
        zna_merge::merge_chunk(b1, n1, pos1, b2, n2, pos2, p, check_sync, base_index,
                               g_chunk_scratch, blob, st);
    }

    // Trailing zero bins are dropped, so the list ends at the largest value observed.
    // That is the contract the reference backend meets by construction (it grows its
    // lists to exactly the index it is about to increment), which is what lets the
    // cross-backend tests compare the histograms element for element.
    auto hist = [](const std::vector<uint32_t>& h) {
        size_t n = h.size();
        while (n > 0 && h[n - 1] == 0) --n;
        nb::list out;
        for (size_t i = 0; i < n; ++i) out.append(h[i]);
        return out;
    };
    return nb::make_tuple(
        nb::bytes(blob.data(), blob.size()), pos1, pos2,
        nb::make_tuple(st.n_pairs, st.merged, st.trimmed, st.kept, st.emitted,
                       st.dropped, st.bases_trimmed, st.frags_short,
                       st.bases_consensus, st.trim_guard, st.sum_olen, st.sum_diff,
                       st.max_read_len, st.npolicy_bases, st.n_rescued),
        hist(st.len_hist), hist(st.olen_hist), hist(st.insert_hist));
}

/// (offset, n_records) just past `max_records` complete records.
static nb::tuple split_records(nb::bytes buf, int64_t start, int64_t max_records) {
    if (start < 0 || static_cast<size_t>(start) > buf.size()) {
        throw std::invalid_argument("split_records(): bad start");
    }
    size_t off = 0;
    int64_t count = 0;
    const auto* b = reinterpret_cast<const uint8_t*>(buf.c_str()) + start;
    const size_t n = buf.size() - static_cast<size_t>(start);
    {
        nb::gil_scoped_release release;
        zna_merge::split_records(b, n, max_records, off, count);
    }
    return nb::make_tuple(static_cast<int64_t>(off) + start, count);   // absolute
}

NB_MODULE(_accel, m) {
    // Raise the Python InputError the rest of the tool already catches, rather than a
    // new type nobody has an `except` for.
    nb::register_exception_translator(
        [](const std::exception_ptr &pe, void *) {
            try {
                std::rethrow_exception(pe);
            } catch (const zna_merge::InputError &e) {
                nb::object cls =
                    nb::module_::import_("zna.merge.fastqio").attr("InputError");
                PyErr_SetString(cls.ptr(), e.what());
            }
        });

    m.doc() = "Accelerated overlap scan for zna merge";
    m.def("scan", &scan,
          nb::arg("s1"), nb::arg("s2rc"), nb::arg("len1"), nb::arg("len2"),
          nb::arg("match_q"), nb::arg("step_q"), nb::arg("floor_q"),
          "Best-scoring shift, as (shift, score_q, overlap_len, mismatches).\n"
          "Scores are integers in zna.merge.params' fixed-point scale.");
    m.def("process_pair", &process_pair,
          nb::arg("h1"), nb::arg("s1"), nb::arg("q1"),
          nb::arg("h2"), nb::arg("s2"), nb::arg("q2"),
          nb::arg("match_q"), nb::arg("step_q"),
          nb::arg("t_merge_q"), nb::arg("t_trim_q"),
          nb::arg("min_read_length"), nb::arg("disagree_q"),
          nb::arg("npolicy") = 1, nb::arg("rng_seed") = 0,
          nb::arg("pair_index") = 0,
          "Classify one pair and build its output records.\n"
          "Returns (records, outcome, n_dropped, score_q, overlap_len, mismatches,\n"
          "         bases_consensus_changed, trim_guard_fired).");

    m.def("merge_chunk", &merge_chunk,
          nb::arg("buf1"), nb::arg("start1"), nb::arg("end1"),
          nb::arg("buf2"), nb::arg("start2"), nb::arg("end2"),
          nb::arg("match_q"), nb::arg("step_q"),
          nb::arg("t_merge_q"), nb::arg("t_trim_q"),
          nb::arg("min_read_length"), nb::arg("disagree_q"),
          nb::arg("check_sync"), nb::arg("base_index"),
          nb::arg("npolicy") = 1, nb::arg("rng_seed") = 0,
          "Merge every whole pair in the two buffers.\n"
          "Returns (blob, consumed1, consumed2, counters, len_hist, olen_hist,\n"
          "         insert_hist). Releases the GIL.");

    m.def("split_records", &split_records,
          nb::arg("buf"), nb::arg("start"), nb::arg("max_records"),
          "Absolute byte offset just past max_records complete FASTQ records, and how many\n"
          "were found, as (offset, n_records).");

#ifdef ZNA_MERGE_V16
    m.attr("VECTOR_WIDTH") = 16;
#else
    m.attr("VECTOR_WIDTH") = 0;
#endif
}
