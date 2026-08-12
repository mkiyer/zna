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

NB_MODULE(_accel, m) {
    m.doc() = "Accelerated overlap scan for zna merge";
    m.def("scan", &scan,
          nb::arg("s1"), nb::arg("s2rc"), nb::arg("len1"), nb::arg("len2"),
          nb::arg("match_q"), nb::arg("step_q"), nb::arg("floor_q"),
          "Best-scoring shift, as (shift, score_q, overlap_len, mismatches).\n"
          "Scores are integers in zna.merge.params' fixed-point scale.");
#ifdef ZNA_MERGE_V16
    m.attr("VECTOR_WIDTH") = 16;
#else
    m.attr("VECTOR_WIDTH") = 0;
#endif
}
