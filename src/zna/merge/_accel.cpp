/**
 * nanobind bindings for the accelerated overlap scan.
 *
 * Thin on purpose: everything with an opinion lives in merge_core.hpp, which knows
 * nothing about Python and can therefore be compiled into a sanitizer driver or reused
 * by `zna encode --merge-pairs` without dragging the interpreter along. This file is
 * argument checking and tuple building.
 */

#include <nanobind/nanobind.h>

#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#include "merge_core.hpp"
#include "fastq_chunk.hpp"

namespace nb = nanobind;

/// The contract `zna.merge.backend` requires (`POLICY_ABI` there and in `_pymerge`).
/// Bumped whenever an argument list or a return tuple changes meaning, so an extension
/// built from an older tree is refused rather than called with the wrong arguments.
constexpr int POLICY_ABI = 2;

namespace {

/// One `int64` policy table, borrowed through the buffer protocol for the call.
///
/// `params.py` hands the tables over as `array('q')`; any C-contiguous buffer of 8-byte
/// signed integers is accepted (format `q`, or `l` where long is 64-bit), and so is a
/// plain byte buffer whose length is a multiple of 8, read as native int64 -- that is
/// what `array.tobytes()` produces. Anything else is a TypeError, not a reinterpretation.
///
/// The view is held until the binding returns, which is after the GIL is reacquired:
/// the exporter cannot resize or free the memory while a worker reads it with the GIL
/// released (an `array` with an export outstanding refuses to resize anyway). A buffer
/// that is not 8-byte aligned -- a `bytes` slice, say -- is copied once rather than read
/// through a misaligned pointer.
///
/// The view lives in a member that releases it, not in the constructor's own cleanup: a
/// constructor that throws never runs its destructor, and every refusal below comes
/// AFTER the export succeeded. Released by hand, a refused table leaked its export and a
/// reference -- 1,000 refused `array('i')` calls left the array at refcount 1,002 and
/// unable to resize.
class I64Table {
    struct View {
        Py_buffer v{};
        bool held = false;
        View() = default;
        View(const View&) = delete;
        View& operator=(const View&) = delete;
        ~View() { if (held) PyBuffer_Release(&v); }
    };

public:
    I64Table(nb::handle obj, const char* fn, const char* what) {
        if (PyObject_GetBuffer(obj.ptr(), &view_.v, PyBUF_C_CONTIGUOUS | PyBUF_FORMAT) != 0) {
            PyErr_Clear();
            throw nb::type_error((std::string(fn) + ": " + what +
                                  " must be an int64 buffer (array('q'))").c_str());
        }
        view_.held = true;
        const Py_buffer& v = view_.v;
        const char* f = v.format ? v.format : "B";
        if (*f == '@' || *f == '=' || *f == '<' || *f == '!' || *f == '>') {
            // '<'/'>'/'!' name a byte order; only native (or little-endian on a
            // little-endian host) is the layout the kernel reads.
            const bool native_le = (*f == '@' || *f == '=');
            const uint16_t probe = 1;
            const bool host_le = *reinterpret_cast<const uint8_t*>(&probe) == 1;
            if (!native_le && !((*f == '<') == host_le)) bad(fn, what);
            ++f;
        }
        const bool word = v.itemsize == 8 &&
            (std::strcmp(f, "q") == 0 ||
             (sizeof(long) == 8 && std::strcmp(f, "l") == 0));
        const bool raw = v.itemsize == 1 &&
            (std::strcmp(f, "B") == 0 || std::strcmp(f, "b") == 0 ||
             std::strcmp(f, "c") == 0) && v.len % 8 == 0;
        if (!word && !raw) bad(fn, what);
        n_ = static_cast<size_t>(v.len) / 8;
        if (reinterpret_cast<uintptr_t>(v.buf) % alignof(int64_t) == 0) {
            p_ = static_cast<const int64_t*>(v.buf);
        } else {
            copy_.resize(n_);
            std::memcpy(copy_.data(), v.buf, n_ * 8);
            p_ = copy_.data();
        }
    }
    I64Table(const I64Table&) = delete;
    I64Table& operator=(const I64Table&) = delete;

    const int64_t* data() const { return p_; }
    size_t size() const { return n_; }

private:
    [[noreturn]] void bad(const char* fn, const char* what) {
        throw nb::type_error((std::string(fn) + ": " + what +
                              " must be a buffer of native int64 (array('q'))").c_str());
    }
    View view_;                        // released on every exit, a throw included
    const int64_t* p_ = nullptr;
    size_t n_ = 0;
    std::vector<int64_t> copy_;
};

/// The bounds that keep the kernel's int64 arithmetic exact, mirrored by
/// `_pymerge.MAX_WEIGHT_Q`/`MAX_FLOOR_Q` (the reference's integers are unbounded, so an
/// input past them would make the two disagree -- and here be signed-overflow UB). With
/// `n < 2^31`, a weight `<= 2^30` keeps `n * match_q` and `d * step_q` under 2^61, and a
/// floor within `+-2^61` keeps `floor - 1` and `ceiling - best - 1` under 2^63. Real
/// values are far inside: `step_q` is 1.4e8 at e = 0.01 and 5.3e8 at the smallest
/// error rate accepted (1e-9), and `T` is ~2^29.
constexpr int64_t MAX_WEIGHT_Q = int64_t{1} << 30;
constexpr int64_t MAX_FLOOR_Q = int64_t{1} << 61;

/// The two tables as the kernel takes them, validated: `T_q` needs at least entries 0
/// and 1 (one base per mate), `dfit` at least entry 0, and every `T_q` entry is a floor
/// within `+-MAX_FLOOR_Q`. Beyond that a table's size is its capacity, and a read past
/// it is either a ValueError (`overlap`, `process_pair`, as in the reference) or the
/// chunk functions' `need`. Mirrors `_pymerge._check_tables`.
inline zna_merge::Tables tables_of(const I64Table& t, const I64Table& d, const char* fn) {
    if (t.size() < 2 || d.size() < 1) {
        throw std::invalid_argument(std::string(fn) +
                                    ": the T and dfit tables are empty");
    }
    for (size_t i = 0; i < t.size(); ++i) {
        if (t.data()[i] < -MAX_FLOOR_Q || t.data()[i] > MAX_FLOOR_Q) {
            throw std::invalid_argument(std::string(fn) + ": T table entry " +
                                        std::to_string(i) + " is out of range");
        }
    }
    return {t.data(), t.size(), d.data(), d.size()};
}

/// Mirrors `_pymerge._check_weights`.
inline void check_weights(int64_t match_q, int64_t step_q, const char* fn) {
    if (match_q <= 0 || step_q <= 0 || match_q > MAX_WEIGHT_Q || step_q > MAX_WEIGHT_Q) {
        throw std::invalid_argument(std::string(fn) +
                                    ": match_q and step_q must be in [1, 2^30]");
    }
}

/// Mirrors `_pymerge._check_npolicy`. The kernel reads any code that is not keep or
/// trim3 as random, and would then tag the read with trim3's provenance bit.
inline void check_npolicy(int npolicy, const char* fn) {
    if (npolicy != zna_merge::NPOLICY_KEEP && npolicy != zna_merge::NPOLICY_TRIM3 &&
        npolicy != zna_merge::NPOLICY_RANDOM) {
        throw std::invalid_argument(std::string(fn) + ": npolicy must be 0, 1 or 2");
    }
}

inline void check_capacity(int len1, int len2, const zna_merge::Tables& tab) {
    if (len1 <= 0 || len2 <= 0) return;            // the reference returns NONE first
    const int longest = len1 > len2 ? len1 : len2;
    if (longest > tab.capacity()) {
        throw std::invalid_argument(
            "read of " + std::to_string(longest) + " bases exceeds the policy tables "
            "(capacity " + std::to_string(tab.capacity()) + ")");
    }
}

}  // namespace

/// Best-scoring eligible shift, as (shift, score_q, overlap_len, mismatches).
static nb::tuple scan(nb::bytes seq1, nb::bytes seq2rc, int len1, int len2,
                      int64_t match_q, int64_t step_q, int64_t floor_q,
                      int adapter_trimmed) {
    if (len1 < 0 || len2 < 0 ||
        static_cast<size_t>(len1) > seq1.size() ||
        static_cast<size_t>(len2) > seq2rc.size()) {
        throw std::invalid_argument("scan(): length argument exceeds the buffer");
    }
    check_weights(match_q, step_q, "scan()");
    if (floor_q < -MAX_FLOOR_Q || floor_q > MAX_FLOOR_Q) {
        throw std::invalid_argument("scan(): floor_q is out of range");
    }
    const zna_merge::ScanResult r = zna_merge::scan(
        reinterpret_cast<const uint8_t*>(seq1.c_str()), len1,
        reinterpret_cast<const uint8_t*>(seq2rc.c_str()), len2,
        match_q, step_q, floor_q, adapter_trimmed != 0);
    return nb::make_tuple(r.shift, r.score_q, r.overlap_len, r.mismatches);
}

/// The authoritative decision, as (verdict, shift, score_q, overlap_len, mismatches,
/// informative). Mirrors `_pymerge.overlap`.
static nb::tuple overlap(nb::bytes seq1, nb::bytes seq2rc, int len1, int len2,
                         int64_t match_q, int64_t step_q, nb::handle t_table,
                         nb::handle dfit_table, int adapter_trimmed) {
    if (len1 < 0 || len2 < 0 ||
        static_cast<size_t>(len1) > seq1.size() ||
        static_cast<size_t>(len2) > seq2rc.size()) {
        throw std::invalid_argument("overlap(): length argument exceeds the buffer");
    }
    check_weights(match_q, step_q, "overlap()");
    const I64Table tt(t_table, "overlap()", "t_table");
    const I64Table dt(dfit_table, "overlap()", "dfit_table");
    const zna_merge::Tables tab = tables_of(tt, dt, "overlap()");
    check_capacity(len1, len2, tab);
    const zna_merge::Decision d = zna_merge::decide(
        reinterpret_cast<const uint8_t*>(seq1.c_str()), len1,
        reinterpret_cast<const uint8_t*>(seq2rc.c_str()), len2,
        match_q, step_q, tab, adapter_trimmed != 0);
    return nb::make_tuple(d.verdict, d.w.shift, d.w.score_q, d.w.overlap_len,
                          d.w.mismatches, d.w.mismatches - d.one_n);
}

namespace {
/// One arena per thread: the per-pair path must not allocate, and the chunk paths run
/// with the GIL released on worker threads.
thread_local zna_merge::Scratch g_scratch;
thread_local zna_merge::ChunkScratch g_chunk_scratch;

inline zna_merge::Span span_of(nb::bytes& b) {
    return {reinterpret_cast<const uint8_t*>(b.c_str()), static_cast<int>(b.size())};
}
inline nb::bytes to_bytes(const zna_merge::Span& s) {
    return nb::bytes(reinterpret_cast<const char*>(s.p), static_cast<size_t>(s.n));
}

inline void check_disagree(const nb::bytes& disagree_q, const char* fn) {
    if (disagree_q.size() != 256u * 256u) {
        throw std::invalid_argument(std::string(fn) + ": disagree_q must be 256*256 bytes");
    }
}

inline zna_merge::Params params_of(int64_t match_q, int64_t step_q,
                                   const zna_merge::Tables& tab, int adapter_trimmed,
                                   int min_read_length, const nb::bytes& disagree_q,
                                   int npolicy, uint64_t rng_seed) {
    return zna_merge::Params{match_q, step_q, tab, adapter_trimmed != 0, min_read_length,
                             reinterpret_cast<const uint8_t*>(disagree_q.c_str()),
                             npolicy, rng_seed};
}

/// The 16 counters, in `_pymerge._Tally.counters()` order.
inline nb::tuple counters_of(const zna_merge::ChunkStats& st) {
    return nb::make_tuple(st.n_pairs, st.merged, st.kept, st.emitted, st.dropped,
                          st.frags_short, st.bases_consensus, st.implausible,
                          st.sum_olen, st.sum_diff, st.max_read_len, st.npolicy_bases,
                          st.n_rescued, st.det_bases, st.det_mismatches, st.rt_strong);
}

/// Trailing zero bins are dropped, so the list ends at the largest value observed.
/// That is the contract the reference backend meets by construction (it grows its lists
/// to exactly the index it is about to increment), which is what lets the cross-backend
/// tests compare the histograms element for element.
inline nb::list hist_of(const std::vector<uint32_t>& h) {
    size_t n = h.size();
    while (n > 0 && h[n - 1] == 0) --n;
    nb::list out;
    for (size_t i = 0; i < n; ++i) out.append(h[i]);
    return out;
}

inline void check_range(const nb::bytes& buf1, int64_t start1, int64_t end1,
                        const nb::bytes& buf2, int64_t start2, int64_t end2,
                        const char* fn) {
    if (start1 < 0 || end1 < start1 || static_cast<size_t>(end1) > buf1.size() ||
        start2 < 0 || end2 < start2 || static_cast<size_t>(end2) > buf2.size()) {
        throw std::invalid_argument(std::string(fn) + ": bad [start, end) range");
    }
}
}  // namespace

/// Classify one pair and build its records. Mirrors _pymerge.process_pair exactly.
static nb::tuple process_pair(nb::bytes h1, nb::bytes s1, nb::bytes q1,
                              nb::bytes h2, nb::bytes s2, nb::bytes q2,
                              int64_t match_q, int64_t step_q,
                              nb::handle t_table, nb::handle dfit_table,
                              int adapter_trimmed, int min_read_length,
                              nb::bytes disagree_q, int npolicy, uint64_t rng_seed,
                              int64_t pair_index, int rt_check) {
    if (s1.size() != q1.size() || s2.size() != q2.size()) {
        throw std::invalid_argument("process_pair(): sequence and quality differ in length");
    }
    check_disagree(disagree_q, "process_pair()");
    check_weights(match_q, step_q, "process_pair()");
    check_npolicy(npolicy, "process_pair()");
    const I64Table tt(t_table, "process_pair()", "t_table");
    const I64Table dt(dfit_table, "process_pair()", "dfit_table");
    const zna_merge::Tables tab = tables_of(tt, dt, "process_pair()");
    const zna_merge::Read r1{span_of(h1), span_of(s1), span_of(q1)};
    const zna_merge::Read r2{span_of(h2), span_of(s2), span_of(q2)};
    check_capacity(r1.s.n, r2.s.n, tab);
    const zna_merge::Params p = params_of(match_q, step_q, tab, adapter_trimmed,
                                          min_read_length, disagree_q, npolicy, rng_seed);

    const zna_merge::PairResult r =
        zna_merge::process_pair(r1, r2, p, g_scratch, pair_index, rt_check != 0);

    nb::list recs;
    for (int i = 0; i < r.n_recs; ++i) {
        recs.append(nb::make_tuple(to_bytes(r.recs[i].h),
                                   to_bytes(r.recs[i].s),
                                   to_bytes(r.recs[i].q)));
    }
    return nb::make_tuple(recs, r.outcome, r.n_dropped, r.shift, r.score_q,
                          r.overlap_len, r.mismatches, r.bases_consensus_changed,
                          r.implausible, r.npolicy_bases, r.n_rescued, r.det_bases,
                          r.det_mismatches, r.rt_strong, r.det_len);
}

/// Merge every whole pair in the two buffers; return the formatted FASTQ text.
///
/// The GIL is released for the whole parse-and-merge, which is what makes threads worth
/// having: the inputs are immutable `bytes` kept alive by the arguments, the tables are
/// held by their buffer views, and nothing here touches a Python object until the GIL
/// is reacquired.
static nb::tuple merge_chunk(nb::bytes buf1, int64_t start1, int64_t end1,
                             nb::bytes buf2, int64_t start2, int64_t end2,
                             int64_t match_q, int64_t step_q,
                             nb::handle t_table, nb::handle dfit_table,
                             int adapter_trimmed, int min_read_length,
                             nb::bytes disagree_q, bool check_sync, int64_t base_index,
                             int npolicy, uint64_t rng_seed, int64_t rt_check_pairs) {
    check_range(buf1, start1, end1, buf2, start2, end2, "merge_chunk()");
    check_disagree(disagree_q, "merge_chunk()");
    check_weights(match_q, step_q, "merge_chunk()");
    check_npolicy(npolicy, "merge_chunk()");
    const I64Table tt(t_table, "merge_chunk()", "t_table");
    const I64Table dt(dfit_table, "merge_chunk()", "dfit_table");
    const zna_merge::Tables tab = tables_of(tt, dt, "merge_chunk()");
    const auto* b1 = reinterpret_cast<const uint8_t*>(buf1.c_str()) + start1;
    const auto* b2 = reinterpret_cast<const uint8_t*>(buf2.c_str()) + start2;
    const size_t n1 = static_cast<size_t>(end1 - start1);
    const size_t n2 = static_cast<size_t>(end2 - start2);
    const zna_merge::Params p = params_of(match_q, step_q, tab, adapter_trimmed,
                                          min_read_length, disagree_q, npolicy, rng_seed);

    std::string blob;
    zna_merge::ChunkStats st;
    size_t pos1 = 0, pos2 = 0;
    {
        nb::gil_scoped_release release;
        blob.reserve(n1 + n2);
        zna_merge::merge_chunk(b1, n1, pos1, b2, n2, pos2, p, check_sync, base_index,
                               rt_check_pairs, g_chunk_scratch, blob, st);
    }
    return nb::make_tuple(
        nb::bytes(blob.data(), blob.size()), pos1, pos2, counters_of(st),
        hist_of(st.len_hist), hist_of(st.olen_hist), hist_of(st.insert_hist),
        hist_of(st.det_olen_hist), st.need);
}

/// Merge every whole pair in the two buffers; return RECORDS instead of text.
///
/// The `zna encode --merge-pairs` adapter.  Mirrors merge_chunk's conventions
/// exactly: consumed counts are RELATIVE to start (the caller does pos += c),
/// while every header offset in `ends` is ABSOLUTE into its bytes object,
/// like split_records' return -- into buf1 for MERGED/MATE1 records, buf2 for
/// MATE2, selected by the slot.
static nb::tuple merge_chunk_records(nb::bytes buf1, int64_t start1, int64_t end1,
                                     nb::bytes buf2, int64_t start2, int64_t end2,
                                     int64_t match_q, int64_t step_q,
                                     nb::handle t_table, nb::handle dfit_table,
                                     int adapter_trimmed, int min_read_length,
                                     nb::bytes disagree_q, bool check_sync,
                                     int64_t base_index, bool want_headers,
                                     int npolicy, uint64_t rng_seed,
                                     int64_t rt_check_pairs) {
    check_range(buf1, start1, end1, buf2, start2, end2, "merge_chunk_records()");
    check_disagree(disagree_q, "merge_chunk_records()");
    check_weights(match_q, step_q, "merge_chunk_records()");
    check_npolicy(npolicy, "merge_chunk_records()");
    const I64Table tt(t_table, "merge_chunk_records()", "t_table");
    const I64Table dt(dfit_table, "merge_chunk_records()", "dfit_table");
    const zna_merge::Tables tab = tables_of(tt, dt, "merge_chunk_records()");
    const auto* b1 = reinterpret_cast<const uint8_t*>(buf1.c_str()) + start1;
    const auto* b2 = reinterpret_cast<const uint8_t*>(buf2.c_str()) + start2;
    const size_t n1 = static_cast<size_t>(end1 - start1);
    const size_t n2 = static_cast<size_t>(end2 - start2);
    const zna_merge::Params p = params_of(match_q, step_q, tab, adapter_trimmed,
                                          min_read_length, disagree_q, npolicy, rng_seed);

    std::string seqs;
    std::vector<zna_merge::RecordEnd> ends;
    zna_merge::ChunkStats st;
    size_t pos1 = 0, pos2 = 0;
    {
        nb::gil_scoped_release release;
        seqs.reserve(n1 / 2);
        zna_merge::merge_chunk_records(b1, n1, pos1, b2, n2, pos2, p, check_sync,
                                       base_index, want_headers, rt_check_pairs,
                                       g_chunk_scratch, seqs, ends, st);
    }

    nb::list ends_out;
    for (const auto& e : ends) {
        // Header offsets become absolute into the caller's bytes object.
        const int64_t habs =
            e.hdr_len == 0 ? 0
            : (e.slot == zna_merge::SLOT_MATE2 ? start2 : start1) + e.hdr_off;
        ends_out.append(nb::make_tuple(e.seq_off, e.seq_len, habs, e.hdr_len,
                                       e.slot, e.prov));
    }
    return nb::make_tuple(
        nb::bytes(seqs.data(), seqs.size()), ends_out, pos1, pos2, counters_of(st),
        hist_of(st.len_hist), hist_of(st.olen_hist), hist_of(st.insert_hist),
        hist_of(st.det_olen_hist), st.need);
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
                // Decoded as latin-1, as the reference decodes read names: the message
                // can quote header bytes, and PyErr_SetString's UTF-8 turned a desync
                // on a name holding 0xff into a UnicodeDecodeError.
                const char* what = e.what();
                PyObject* msg = PyUnicode_DecodeLatin1(
                    what, static_cast<Py_ssize_t>(std::strlen(what)), nullptr);
                if (msg) {
                    PyErr_SetObject(cls.ptr(), msg);
                    Py_DECREF(msg);
                }
            }
        });

    m.doc() = "Accelerated overlap scan for zna merge";
    m.attr("POLICY_ABI") = POLICY_ABI;
    m.def("scan", &scan,
          nb::arg("s1"), nb::arg("s2rc"), nb::arg("len1"), nb::arg("len2"),
          nb::arg("match_q"), nb::arg("step_q"), nb::arg("floor_q"),
          nb::arg("adapter_trimmed") = 0,
          "Best-scoring eligible shift, as (shift, score_q, overlap_len, mismatches).\n"
          "Scores are integers in zna.merge.params' fixed-point scale.");
    m.def("overlap", &overlap,
          nb::arg("s1"), nb::arg("s2rc"), nb::arg("len1"), nb::arg("len2"),
          nb::arg("match_q"), nb::arg("step_q"), nb::arg("t_table"),
          nb::arg("dfit_table"), nb::arg("adapter_trimmed"),
          "The overlap decision, as (verdict, shift, score_q, overlap_len, mismatches,\n"
          "informative_mismatches). Tables are int64 buffers (array('q')).");
    m.def("process_pair", &process_pair,
          nb::arg("h1"), nb::arg("s1"), nb::arg("q1"),
          nb::arg("h2"), nb::arg("s2"), nb::arg("q2"),
          nb::arg("match_q"), nb::arg("step_q"),
          nb::arg("t_table"), nb::arg("dfit_table"), nb::arg("adapter_trimmed"),
          nb::arg("min_read_length"), nb::arg("disagree_q"),
          nb::arg("npolicy") = 1, nb::arg("rng_seed") = 0,
          nb::arg("pair_index") = 0, nb::arg("rt_check") = 0,
          "Classify one pair and build its output records.\n"
          "Returns (records, outcome, n_dropped, shift, score_q, overlap_len,\n"
          "         mismatches, bases_consensus_changed, implausible, npolicy_bases,\n"
          "         n_rescued, detected_bases, detected_mismatches,\n"
          "         readthrough_strong, detected_overlap_len).");

    m.def("merge_chunk", &merge_chunk,
          nb::arg("buf1"), nb::arg("start1"), nb::arg("end1"),
          nb::arg("buf2"), nb::arg("start2"), nb::arg("end2"),
          nb::arg("match_q"), nb::arg("step_q"),
          nb::arg("t_table"), nb::arg("dfit_table"), nb::arg("adapter_trimmed"),
          nb::arg("min_read_length"), nb::arg("disagree_q"),
          nb::arg("check_sync"), nb::arg("base_index"),
          nb::arg("npolicy") = 1, nb::arg("rng_seed") = 0,
          nb::arg("rt_check_pairs") = 0,
          "Merge every whole pair in the two buffers.\n"
          "Returns (blob, consumed1, consumed2, counters, len_hist, olen_hist,\n"
          "         insert_hist, det_olen_hist, need). Releases the GIL.");

    m.def("merge_chunk_records", &merge_chunk_records,
          nb::arg("buf1"), nb::arg("start1"), nb::arg("end1"),
          nb::arg("buf2"), nb::arg("start2"), nb::arg("end2"),
          nb::arg("match_q"), nb::arg("step_q"),
          nb::arg("t_table"), nb::arg("dfit_table"), nb::arg("adapter_trimmed"),
          nb::arg("min_read_length"), nb::arg("disagree_q"),
          nb::arg("check_sync"), nb::arg("base_index"),
          nb::arg("want_headers"), nb::arg("npolicy") = 1, nb::arg("rng_seed") = 0,
          nb::arg("rt_check_pairs") = 0,
          "Merge every whole pair in the two buffers, emitting records.\n"
          "Returns (seqs, ends, consumed1, consumed2, counters, len_hist,\n"
          "         olen_hist, insert_hist, det_olen_hist, need).  Each end is\n"
          "         (seq_off, seq_len,\n"
          "         hdr_off, hdr_len, slot, prov); hdr offsets are ABSOLUTE\n"
          "         into buf1 (MERGED/MATE1) or buf2 (MATE2), consumed counts\n"
          "         RELATIVE to start.  Releases the GIL.");

    m.def("split_records", &split_records,
          nb::arg("buf"), nb::arg("start"), nb::arg("max_records"),
          "Absolute byte offset just past max_records complete FASTQ records, and how many\n"
          "were found, as (offset, n_records).");

    // Test hook, not API -- hence the leading underscore, and hence its absence from
    // backend._REQUIRED_FUNCTIONS. It exists so the SWAR popcount that only MSVC
    // actually calls is still checked on every platform the suite runs on; the
    // alternative is a branch whose only build is the one that cannot test it.
    m.def("_popcount16_portable",
          [](unsigned x) { return zna_merge::popcount16_portable(x); },
          nb::arg("x"),
          "Population count of the low 16 bits, via the portable SWAR fold.");

#ifdef ZNA_MERGE_V16
    m.attr("VECTOR_WIDTH") = 16;
#else
    m.attr("VECTOR_WIDTH") = 0;
#endif
}
