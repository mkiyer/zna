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
                              int min_read_length, nb::bytes disagree_q) {
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
                              reinterpret_cast<const uint8_t*>(disagree_q.c_str())};

    const zna_merge::PairResult r = zna_merge::process_pair(r1, r2, p, g_scratch);

    nb::list recs;
    for (int i = 0; i < r.n_recs; ++i) {
        recs.append(nb::make_tuple(to_bytes(r.recs[i].h),
                                   to_bytes(r.recs[i].s),
                                   to_bytes(r.recs[i].q)));
    }
    return nb::make_tuple(recs, r.outcome, r.n_dropped, r.score_q, r.overlap_len,
                          r.mismatches, r.bases_consensus_changed, r.trim_guard_fired);
}

NB_MODULE(_accel, m) {
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
          "Classify one pair and build its output records.\n"
          "Returns (records, outcome, n_dropped, score_q, overlap_len, mismatches,\n"
          "         bases_consensus_changed, trim_guard_fired).");

#ifdef ZNA_MERGE_V16
    m.attr("VECTOR_WIDTH") = 16;
#else
    m.attr("VECTOR_WIDTH") = 0;
#endif
}
