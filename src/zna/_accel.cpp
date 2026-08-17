/**
 * ZNA Accelerated Encode/Decode Functions
 *
 * High-performance C++ implementation of DNA sequence encoding/decoding.
 * Uses nanobind for Python bindings.
 *
 * Design note — why this file talks to the CPython C API directly
 * ---------------------------------------------------------------
 * The obvious nanobind signatures (``std::vector<std::string>`` in,
 * ``std::vector<std::tuple<...>>`` out) copy every base twice on top of the
 * codec's own work: once for nanobind's STL round-trip, once more building the
 * Python objects.  On a 4 MiB block that is millions of bases of pure copying
 * either side of the loop that actually does something.
 *
 * So the hot paths here take the Python objects as they are:
 *
 *   - decode writes bases straight into a fresh ``PyUnicode`` buffer, four at a
 *     time from a 256-entry lookup table, and hands back a hand-built list of
 *     tuples;
 *   - encode reads each sequence's character buffer in place, and folds
 *     reverse-complement into the packing loop so no intermediate string exists;
 *   - emitted record tuples are GC-untracked, since a tuple of a string and
 *     some bools cannot be part of a reference cycle.
 *
 * This costs portability: ``PyUnicode_New``, ``PyUnicode_DATA``,
 * ``PyList_SET_ITEM`` and ``PyTuple_SET_ITEM`` are not in the limited API, and
 * on macOS ``PyUnicode_New`` additionally needs an explicit linker opt-in (see
 * ``CMakeLists.txt``).  The module was never actually stable-ABI --
 * ``PyTuple_SET_ITEM`` predates this rewrite -- and ``CMakeLists.txt`` no longer
 * claims otherwise.
 */

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/pair.h>
#include <nanobind/ndarray.h>

#include <cstdint>
#include <cstring>  // for memcpy
#include <cstdlib>  // for strtol, strtod
#include <stdexcept>
#include <string>
#include <vector>

namespace nb = nanobind;

// Base encoding: A=0, C=1, G=2, T=3
constexpr uint8_t BASE_A = 0;
constexpr uint8_t BASE_C = 1;
constexpr uint8_t BASE_G = 2;
constexpr uint8_t BASE_T = 3;
constexpr uint8_t INVALID = 255;

/// splitmix64's finalizer. The substitution and orientation streams are derived from it
/// rather than carried as mutable state, so neither depends on how records are batched.
/// Integer-only and written identically in `_pycodec.py`, so the two backends cannot
/// drift; see docs/NPOLICY_PLAN.md 8.3.
static inline uint64_t zna_mix64(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

/// Which mate of the pair at global record index `r` is reverse-complemented.
static inline uint64_t zna_rc_coin(uint64_t seed, uint64_t r) {
    return zna_mix64(seed + 0x9E3779B97F4A7C15ULL * (r + 1)) & 1ULL;
}

/// The base substituted for the no-call at `off` of the STORED record at index `r`.
/// A different multiplier from the coin, so the two streams never shift each other.
static inline uint8_t zna_sub_base(uint64_t seed, uint64_t r, uint64_t off) {
    return static_cast<uint8_t>(
        zna_mix64(seed + 0xBF58476D1CE4E5B9ULL * (r + 1)
                       + 0x94D049BB133111EBULL * (off + 1)) & 3ULL);
}

// Lookup table for encoding (char -> 2-bit value)
alignas(64) static uint8_t ENCODE_TABLE[256];

// Lookup table for decoding (2-bit value -> char)
static constexpr char DECODE_CHARS[4] = {'A', 'C', 'G', 'T'};

// Complement table for reverse complement (A<->T, C<->G)
static constexpr char COMPLEMENT_TABLE[4] = {'T', 'G', 'C', 'A'};

/// Packed byte -> its four ASCII bases, so decoding is one 4-byte store per
/// input byte instead of four shifts, four masks and four stores.
alignas(64) static uint32_t DECODE_LUT4[256];

/// ASCII base -> its complement, for reversing an already-decoded buffer in
/// place.  Decoded buffers hold only ACGT; the rest of the table is identity so
/// the same table can serve ``reverse_complement`` on arbitrary input.
alignas(64) static char COMPLEMENT_ASCII[256];

// Initialize lookup tables at module load
static void init_tables() noexcept {
    for (int i = 0; i < 256; i++) {
        ENCODE_TABLE[i] = INVALID;
    }
    ENCODE_TABLE['A'] = BASE_A;
    ENCODE_TABLE['a'] = BASE_A;
    ENCODE_TABLE['C'] = BASE_C;
    ENCODE_TABLE['c'] = BASE_C;
    ENCODE_TABLE['G'] = BASE_G;
    ENCODE_TABLE['g'] = BASE_G;
    ENCODE_TABLE['T'] = BASE_T;
    ENCODE_TABLE['t'] = BASE_T;

    for (int v = 0; v < 256; v++) {
        // memcpy in and out preserves byte order on either endianness.
        const char four[4] = {
            DECODE_CHARS[(v >> 6) & 3], DECODE_CHARS[(v >> 4) & 3],
            DECODE_CHARS[(v >> 2) & 3], DECODE_CHARS[v & 3],
        };
        std::memcpy(&DECODE_LUT4[v], four, 4);
    }

    // Identity everywhere, then the eight bases.  Complements are uppercase
    // because that is what this backend has always returned: it round-trips
    // through the 2-bit table, which is case-folding.  The Python backend
    // preserves case instead; tests/test_fuzz_roundtrip.py pins both.
    for (int i = 0; i < 256; i++) {
        COMPLEMENT_ASCII[i] = static_cast<char>(i);
    }
    COMPLEMENT_ASCII[static_cast<uint8_t>('A')] = 'T';
    COMPLEMENT_ASCII[static_cast<uint8_t>('a')] = 'T';
    COMPLEMENT_ASCII[static_cast<uint8_t>('C')] = 'G';
    COMPLEMENT_ASCII[static_cast<uint8_t>('c')] = 'G';
    COMPLEMENT_ASCII[static_cast<uint8_t>('G')] = 'C';
    COMPLEMENT_ASCII[static_cast<uint8_t>('g')] = 'C';
    COMPLEMENT_ASCII[static_cast<uint8_t>('T')] = 'A';
    COMPLEMENT_ASCII[static_cast<uint8_t>('t')] = 'A';
}

// Static initializer
static struct TableInitializer {
    TableInitializer() { init_tables(); }
} table_initializer;


// --- Helper for unaligned little-endian reads ---
inline uint16_t read_u16_le(const uint8_t* ptr) {
    uint16_t val;
    std::memcpy(&val, ptr, 2);
    return val;
}

inline uint32_t read_u32_le(const uint8_t* ptr) {
    uint32_t val;
    std::memcpy(&val, ptr, 4);
    return val;
}


// ---------------------------------------------------------------------------
// Decode kernels
// ---------------------------------------------------------------------------

/// Expand *seq_len* bases from packed *src* into *out*.
///
/// *out* must have room for *seq_len* bytes.  Full 4-base groups are written
/// with one 4-byte store each; the 1-3 base tail is written byte-wise so the
/// last store can never run past the caller's buffer.
static inline void decode_into(char* out, const uint8_t* src, size_t seq_len) noexcept {
    const size_t full = seq_len >> 2;
    for (size_t i = 0; i < full; ++i) {
        std::memcpy(out + (i << 2), &DECODE_LUT4[src[i]], 4);
    }
    const size_t rem = seq_len & 3;
    if (rem) {
        char tail[4];
        std::memcpy(tail, &DECODE_LUT4[src[full]], 4);
        char* dst = out + (full << 2);
        for (size_t k = 0; k < rem; ++k) {
            dst[k] = tail[k];
        }
    }
}

/// Reverse-complement an ASCII buffer in place.
static inline void revcomp_inplace(char* s, size_t n) noexcept {
    if (n < 2) {
        if (n == 1) s[0] = COMPLEMENT_ASCII[static_cast<uint8_t>(s[0])];
        return;
    }
    size_t i = 0, j = n - 1;
    while (i < j) {
        const char a = COMPLEMENT_ASCII[static_cast<uint8_t>(s[i])];
        const char b = COMPLEMENT_ASCII[static_cast<uint8_t>(s[j])];
        s[i++] = b;
        s[j--] = a;
    }
    if (i == j) {
        s[i] = COMPLEMENT_ASCII[static_cast<uint8_t>(s[i])];
    }
}

/// Build a Python ``str`` holding *seq_len* decoded bases.  Returns a new
/// reference, or nullptr with the error indicator set.
///
/// The bases are written straight into the string object's own storage.
/// ``PyUnicode_New(n, 127)`` gives a compact ASCII object with an n+1 byte
/// buffer, which is exactly what the decoder wants to fill — so decoding is one
/// pass over the data and nothing is copied afterwards.
///
/// The obvious alternative, decoding into a scratch buffer and calling
/// ``PyUnicode_FromStringAndSize``, costs three passes rather than one: that
/// function treats its input as UTF-8, so it scans for non-ASCII and then copies.
/// Measured on 200k 150 bp records it gave up a third of the win (raw
/// decode_block +89.8% instead of +120.3%), which is why the linker flag below
/// is worth carrying.
static inline PyObject* make_seq_object(const uint8_t* src, size_t seq_len,
                                        bool rc) {
    if (seq_len == 0) {
        // PyUnicode_New(0, ...) hands back the interned empty string, which must
        // never be written to.
        return PyUnicode_FromStringAndSize("", 0);
    }
    PyObject* s = PyUnicode_New(static_cast<Py_ssize_t>(seq_len), 127);
    if (s == nullptr) {
        return nullptr;
    }
    char* buf = reinterpret_cast<char*>(PyUnicode_DATA(s));
    decode_into(buf, src, seq_len);
    if (rc) {
        revcomp_inplace(buf, seq_len);
    }
    return s;
}

/// Read record *rec*'s stored sequence length.
static inline size_t read_len(const uint8_t* lengths_ptr, int len_bytes, int rec) noexcept {
    if (len_bytes == 1) return lengths_ptr[rec];
    if (len_bytes == 2) return read_u16_le(lengths_ptr + rec * 2);
    return read_u32_le(lengths_ptr + rec * 4);
}

/// Fill the four flag slots common to every emitted record tuple.
static inline void set_flag_items(PyObject* t, uint8_t flag, Py_ssize_t base,
                                  bool with_rc) noexcept {
    PyTuple_SET_ITEM(t, base + 0, Py_NewRef((flag & 4) ? Py_True : Py_False));  // is_paired
    PyTuple_SET_ITEM(t, base + 1, Py_NewRef((flag & 1) ? Py_True : Py_False));  // is_read1
    PyTuple_SET_ITEM(t, base + 2, Py_NewRef((flag & 2) ? Py_True : Py_False));  // is_read2
    if (with_rc) {
        PyTuple_SET_ITEM(t, base + 3, Py_NewRef((flag & 8) ? Py_True : Py_False));
    }
}


/**
 * Encode DNA sequence to 2-bit packed bytes.
 *
 * Each base (A, C, G, T) is encoded as 2 bits (00, 01, 10, 11).
 * Four bases are packed into each byte.
 * N nucleotides are treated as errors (should be handled before calling this).
 */
nb::bytes encode_sequence(const std::string& seq) {
    const size_t length = seq.size();
    const size_t out_len = (length + 3) / 4;

    std::string out;
    out.resize(out_len);
    const uint8_t* seq_ptr = reinterpret_cast<const uint8_t*>(seq.data());
    uint8_t* out_ptr = reinterpret_cast<uint8_t*>(out.data());

    size_t i = 0;
    size_t idx = 0;
    const size_t full_chunks = length / 4;

    // Process 4 bases at a time
    while (idx < full_chunks) {
        uint8_t b0 = ENCODE_TABLE[seq_ptr[i]];
        uint8_t b1 = ENCODE_TABLE[seq_ptr[i + 1]];
        uint8_t b2 = ENCODE_TABLE[seq_ptr[i + 2]];
        uint8_t b3 = ENCODE_TABLE[seq_ptr[i + 3]];

        // Check for invalid characters
        if (b0 == INVALID || b1 == INVALID || b2 == INVALID || b3 == INVALID) {
            throw std::invalid_argument("Invalid character in sequence");
        }

        out_ptr[idx] = (b0 << 6) | (b1 << 4) | (b2 << 2) | b3;
        idx++;
        i += 4;
    }

    // Handle remaining 1-3 bases
    if (i < length) {
        uint8_t val = 0;
        int shift = 6;
        while (i < length) {
            uint8_t b = ENCODE_TABLE[seq_ptr[i]];
            if (b == INVALID) {
                throw std::invalid_argument("Invalid character in sequence");
            }
            val |= (b << shift);
            i++;
            shift -= 2;
        }
        out_ptr[idx] = val;
    }

    return nb::bytes(out.data(), out.size());
}


/**
 * Decode 2-bit packed bytes to DNA sequence string.
 */
std::string decode_sequence(const uint8_t* data, size_t data_len, size_t seq_len) {
    (void)data_len;
    std::string out;
    out.resize(seq_len);
    if (seq_len) {
        decode_into(out.data(), data, seq_len);
    }
    return out;
}


// Record tuple type: (sequence, is_paired, is_read1, is_read2, is_rc)
using Record = std::tuple<std::string, bool, bool, bool, bool>;

/**
 * Decode all records from a block at once.
 *
 * Legacy single-blob entry point, kept for API compatibility; the columnar
 * decoder below is what the reader actually uses.
 */
std::vector<Record> decode_block_records(
    nb::bytes block_data,
    int len_bytes,
    int count
) {
    const uint8_t* data = reinterpret_cast<const uint8_t*>(block_data.c_str());
    const size_t data_size = block_data.size();

    std::vector<Record> results;
    results.reserve(count);

    size_t offset = 0;

    for (int rec = 0; rec < count; rec++) {
        // Safety Check 1: Ensure we have at least 1 byte for flags
        if (offset + 1 > data_size) {
            throw std::runtime_error("Block truncated: cannot read flags byte");
        }

        // Read flags
        uint8_t flags = data[offset++];
        bool is_read1 = (flags & 1) != 0;
        bool is_read2 = (flags & 2) != 0;
        bool is_paired = (flags & 4) != 0;

        // Safety Check 2: Ensure we have length bytes
        if (offset + static_cast<size_t>(len_bytes) > data_size) {
            throw std::runtime_error("Block truncated: cannot read length header");
        }

        // Read sequence length using optimized unaligned reads
        size_t seq_len = 0;
        if (len_bytes == 1) {
            seq_len = data[offset];
        } else if (len_bytes == 2) {
            seq_len = read_u16_le(data + offset);
        } else {  // len_bytes == 4
            seq_len = read_u32_le(data + offset);
        }
        offset += len_bytes;

        // Calculate encoded length
        size_t enc_len = (seq_len + 3) / 4;

        // Safety Check 3: Ensure we have sequence data
        if (offset + enc_len > data_size) {
            throw std::runtime_error("Block truncated: cannot read sequence data");
        }

        // Decode sequence
        std::string seq = decode_sequence(data + offset, enc_len, seq_len);
        offset += enc_len;

        results.emplace_back(std::move(seq), is_paired, is_read1, is_read2, false);
    }

    return results;
}


/**
 * Compute reverse complement of a DNA sequence.
 *
 * Writes straight into the result string's buffer; A<->T, C<->G.  Bases outside
 * ACGTacgt are mirrored without complementing, and complements come back
 * uppercase (this backend folds case; the Python one does not).
 */
nb::object reverse_complement(nb::object seq_obj) {
    PyObject* src_obj = seq_obj.ptr();
    if (!PyUnicode_Check(src_obj)) {
        throw std::invalid_argument("reverse_complement() requires a str");
    }
    Py_ssize_t slen = 0;
    const char* src = PyUnicode_AsUTF8AndSize(src_obj, &slen);
    if (src == nullptr) {
        throw nb::python_error();
    }
    if (slen == 0) {
        return nb::steal(PyUnicode_FromStringAndSize("", 0));
    }

    std::vector<char> buf(static_cast<size_t>(slen));
    char* dst = buf.data();
    for (Py_ssize_t i = 0; i < slen; ++i) {
        dst[i] = COMPLEMENT_ASCII[static_cast<uint8_t>(src[slen - 1 - i])];
    }
    PyObject* out = PyUnicode_FromStringAndSize(dst, slen);
    if (out == nullptr) {
        throw nb::python_error();
    }
    return nb::steal(out);
}


// ---------------------------------------------------------------------------
// Encoder
// ---------------------------------------------------------------------------

/// One input sequence, borrowed from the Python list that owns it.
struct SeqView {
    const uint8_t* data;
    size_t len;
};

/// Pack *slen* bases into *out*, 2 bits each, four bases per byte.
///
/// Templated on direction so the reverse-complement case is a separate compiled
/// loop rather than a branch on every base.  With ``RC`` the source is read back
/// to front and each base is complemented as it is read, which is why encoding a
/// reverse-complemented record needs no intermediate string at all.
///
/// N-policy is applied *after* complementing, matching the two-step form this
/// replaces (reverse-complement leaves N as N, then the policy substitutes it).
/// The substituted base is therefore not itself complemented, and the random
/// policy draws in output order, so the packed bytes are identical either way.
template <bool RC>
static inline void pack_sequence(uint8_t* out, const uint8_t* seq_bytes, size_t slen,
                                 bool has_npolicy, uint64_t rng_seed,
                                 uint64_t rec_index) {
    auto code_at = [&](size_t k) -> uint8_t {
        uint8_t c;
        if constexpr (RC) {
            c = ENCODE_TABLE[seq_bytes[slen - 1 - k]];
        } else {
            c = ENCODE_TABLE[seq_bytes[k]];
        }
        if (c == INVALID) {
            // Which byte was it? `k` is the STORED offset, so report the source index.
            const size_t src_k = RC ? (slen - 1 - k) : k;
            const uint8_t raw = static_cast<uint8_t>(seq_bytes[src_k]);
            // Only a no-call is substitutable. Every other non-ACGT byte -- the IUPAC
            // codes, punctuation, anything -- is an error under EVERY policy. This used
            // to substitute silently, so `--npolicy A` quietly turned an `R` into an `A`
            // and the compiled backend disagreed with the reference.
            if (raw != 'N' && raw != 'n') {
                char msg[80];
                std::snprintf(msg, sizeof(msg),
                              "invalid character '%c' at offset %zu in sequence",
                              (raw >= 32 && raw < 127) ? static_cast<char>(raw) : '?',
                              src_k);
                throw std::invalid_argument(msg);
            }
            if (!has_npolicy) {
                throw std::invalid_argument("Invalid character in sequence");
            }
            // Position-derived, so the answer does not depend on how records are batched.
            return zna_sub_base(rng_seed, static_cast<uint64_t>(rec_index),
                                static_cast<uint64_t>(k));
        }
        if constexpr (RC) {
            return static_cast<uint8_t>(3 - c);
        } else {
            return c;
        }
    };

    size_t s_idx = 0;
    size_t o = 0;
    while (s_idx + 4 <= slen) {
        const uint8_t b0 = code_at(s_idx);
        const uint8_t b1 = code_at(s_idx + 1);
        const uint8_t b2 = code_at(s_idx + 2);
        const uint8_t b3 = code_at(s_idx + 3);
        out[o++] = static_cast<uint8_t>((b0 << 6) | (b1 << 4) | (b2 << 2) | b3);
        s_idx += 4;
    }
    if (s_idx < slen) {
        uint8_t val = 0;
        int shift = 6;
        while (s_idx < slen) {
            val |= static_cast<uint8_t>(code_at(s_idx) << shift);
            s_idx++;
            shift -= 2;
        }
        out[o] = val;
    }
}

/// Columnar output of one encode pass.
struct EncodeResult {
    std::string flags_out;
    std::string lengths_out;
    std::string seqs_out;
};

/// Borrow every sequence's character buffer from the Python list.
///
/// The buffers stay valid because the list holds a reference to each string for
/// the whole call.  Taking them this way is the point of the exercise: the
/// previous ``std::vector<std::string>`` parameter had nanobind copy every base
/// into a fresh heap allocation before the encoder had done any work at all.
static void collect_seq_views(PyObject* seqs, size_t count,
                              std::vector<SeqView>& views, size_t& total_enc_bytes) {
    views.resize(count);
    total_enc_bytes = 0;
    for (size_t i = 0; i < count; ++i) {
        PyObject* item = PyList_GET_ITEM(seqs, static_cast<Py_ssize_t>(i));  // borrowed
        if (!PyUnicode_Check(item)) {
            throw std::invalid_argument(
                "seqs must be a list of str (got a non-string at index "
                + std::to_string(i) + ")");
        }
        Py_ssize_t slen = 0;
        const char* sdata = PyUnicode_AsUTF8AndSize(item, &slen);
        if (sdata == nullptr) {
            throw nb::python_error();
        }
        views[i].data = reinterpret_cast<const uint8_t*>(sdata);
        views[i].len = static_cast<size_t>(slen);
        total_enc_bytes += (static_cast<size_t>(slen) + 3) / 4;
    }
}

/// Shared encode body for both entry points.
///
/// ``encode_block`` and ``encode_block_labeled`` were hand-copied duplicates of
/// this logic, including the strand-normalization rules; a bug fixed in one
/// would not have reached the other.  They now differ only in what they return.
static void encode_core(PyObject* seqs, const std::vector<uint8_t>& flags,
                        int len_bytes_fmt, const std::string& npolicy,
                        bool do_rc_r1, bool do_rc_r2, bool do_random_rc,
                        uint64_t rng_seed, uint64_t base_index,
                        EncodeResult& out) {
    if (!PyList_Check(seqs)) {
        throw std::invalid_argument("seqs must be a list of str");
    }
    const size_t count = static_cast<size_t>(PyList_GET_SIZE(seqs));
    if (count != flags.size()) {
        throw std::invalid_argument("seqs and flags must have the same length");
    }
    if (len_bytes_fmt != 1 && len_bytes_fmt != 2 && len_bytes_fmt != 4) {
        throw std::invalid_argument("len_bytes_fmt must be 1, 2, or 4");
    }

    constexpr uint8_t IS_RC_BIT = 0x08;
    constexpr uint8_t IS_READ1 = 0x01;
    constexpr uint8_t IS_READ2 = 0x02;
    constexpr uint8_t IS_PAIRED = 0x04;

    std::vector<SeqView> views;
    size_t total_enc_bytes = 0;
    collect_seq_views(seqs, count, views, total_enc_bytes);

    out.flags_out.assign(reinterpret_cast<const char*>(flags.data()), count);
    out.lengths_out.assign(static_cast<size_t>(count) * len_bytes_fmt, '\0');
    // Sized exactly from the input lengths, so packing writes through a raw
    // pointer with no capacity check per byte.
    out.seqs_out.assign(total_enc_bytes, '\0');

    char* flags_out = out.flags_out.data();
    char* lengths_out = out.lengths_out.data();
    uint8_t* seqs_out = reinterpret_cast<uint8_t*>(out.seqs_out.data());

    // The policy set is CLOSED. It used to be a chain of equality tests with a silent
    // `A` default, which meant every unrecognised string -- including "drop", the CLI's
    // own default -- substituted A without a word. Reject instead.
    const bool has_npolicy = npolicy == "random";
    if (!npolicy.empty() && npolicy != "random") {
        throw std::invalid_argument(
            "unknown npolicy '" + npolicy + "'; expected 'random' or none. "
            "'trim3' is applied before the codec, and 'drop' is a record filter, "
            "not an encoding policy.");
    }

    size_t max_len = 0;
    if (len_bytes_fmt == 1) max_len = 255;
    else if (len_bytes_fmt == 2) max_len = 65535;
    else max_len = 4294967295UL;

    size_t seq_off = 0;

    for (size_t i = 0; i < count; ++i) {
        uint8_t flag = static_cast<uint8_t>(flags_out[i]);

        const bool is_read1 = (flag & IS_READ1) != 0;
        const bool is_read2 = (flag & IS_READ2) != 0;
        const bool is_paired = (flag & IS_PAIRED) != 0;
        bool needs_rc = false;

        if (do_random_rc) {
            // Unstranded random normalization: exactly one mate of each pair is flipped.
            if (is_paired && is_read1) {
                // The writer guarantees a paired R1's mate is the next record in this
                // same block, which is what lets the coin be applied to the pair as a
                // unit. A caller reaching the backend directly can still break it, and
                // gets this rather than a silent half-normalized fragment -- or, without
                // the bounds half, a read past the end of the flags buffer.
                if (i + 1 >= count || !(flags[i + 1] & IS_READ2)) {
                    throw std::invalid_argument(
                        "record " + std::to_string(base_index + i) + " is a paired R1 "
                        "whose R2 does not follow it in this block; a fragment's reads "
                        "must be consecutive and in the same block"
                    );
                }
                // Which mate is flipped is derived from the pair's GLOBAL record index,
                // not from a running stream, so the assignment cannot depend on how
                // records were batched into blocks.
                if (zna_rc_coin(rng_seed, base_index + i)) {
                    needs_rc = true;
                    flags_out[i] = static_cast<char>(flag | IS_RC_BIT);
                } else {
                    flags_out[i + 1] = static_cast<char>(flags[i + 1] | IS_RC_BIT);
                }
            } else if (!(is_paired && is_read2)) {
                // A single or merged read -- a one-record fragment.  A paired R2 is
                // skipped here because the branch above already flipped the coin for
                // its fragment; it picks up the result via the IS_RC recovery below.
                if (zna_rc_coin(rng_seed, base_index + i)) {
                    needs_rc = true;
                    flags_out[i] = static_cast<char>(flag | IS_RC_BIT);
                }
            }
            // Check if this is an R2 that was already marked for RC by its R1
            if (!needs_rc && (static_cast<uint8_t>(flags_out[i]) & IS_RC_BIT)) {
                needs_rc = true;
            }
        } else {
            // Deterministic stranded normalization.
            // A single/merged read (neither R1 nor R2) is R1-oriented:
            // normalize it with the read1 rule.
            const bool is_single = !is_read1 && !is_read2;
            needs_rc = (do_rc_r1 && (is_read1 || is_single)) || (do_rc_r2 && is_read2);
            if (needs_rc) {
                flags_out[i] = static_cast<char>(flag | IS_RC_BIT);
            }
        }

        const size_t slen = views[i].len;
        if (slen > max_len) {
            throw std::invalid_argument(
                "Sequence length " + std::to_string(slen) +
                " exceeds maximum " + std::to_string(max_len)
            );
        }

        // --- Lengths (little-endian) ---
        if (len_bytes_fmt == 1) {
            lengths_out[i] = static_cast<char>(slen);
        } else if (len_bytes_fmt == 2) {
            uint16_t s = static_cast<uint16_t>(slen);
            std::memcpy(&lengths_out[i * 2], &s, 2);
        } else {
            uint32_t s = static_cast<uint32_t>(slen);
            std::memcpy(&lengths_out[i * 4], &s, 4);
        }

        // --- Sequence (2-bit packing, reverse-complementing as it reads) ---
        if (needs_rc) {
            pack_sequence<true>(seqs_out + seq_off, views[i].data, slen,
                                has_npolicy, rng_seed, base_index + i);
        } else {
            pack_sequence<false>(seqs_out + seq_off, views[i].data, slen,
                                 has_npolicy, rng_seed, base_index + i);
        }
        seq_off += (slen + 3) / 4;
    }
}


/**
 * BATCH ENCODER: builds the three columnar streams entirely in C++.
 *
 * Input: list of str, flags, len format, N policy, strand norm masks
 * Output: Tuple(flags_bytes, lengths_bytes, sequences_bytes)
 */
nb::tuple encode_block(
    nb::object seqs,
    const std::vector<uint8_t>& flags,
    int len_bytes_fmt,
    const std::string& npolicy,
    bool do_rc_r1,
    bool do_rc_r2,
    bool do_random_rc,
    uint64_t rng_seed,
    uint64_t base_index
) {
    EncodeResult r;
    encode_core(seqs.ptr(), flags, len_bytes_fmt, npolicy,
                do_rc_r1, do_rc_r2, do_random_rc, rng_seed, base_index, r);
    return nb::make_tuple(
        nb::bytes(r.flags_out.data(), r.flags_out.size()),
        nb::bytes(r.lengths_out.data(), r.lengths_out.size()),
        nb::bytes(r.seqs_out.data(), r.seqs_out.size())
    );
}


/**
 * BATCH ENCODER WITH LABEL SUPPORT.
 *
 * When label_col_data and label_col_sizes are provided (non-empty),
 * returns a 4-tuple: (flags_bytes, labels_bytes, lengths_bytes, seqs_bytes).
 * The labels_bytes are the pre-packed columnar label data concatenated
 * together (caller packs via struct before calling).
 *
 * When label args are empty, returns the original 3-tuple.
 */
nb::tuple encode_block_labeled(
    nb::object seqs,
    const std::vector<uint8_t>& flags,
    int len_bytes_fmt,
    const std::string& npolicy,
    bool do_rc_r1,
    bool do_rc_r2,
    bool do_random_rc,
    const std::vector<nb::bytes>& label_col_data,
    const std::vector<int>& label_col_sizes,
    uint64_t rng_seed,
    uint64_t base_index
) {
    (void)label_col_sizes;
    EncodeResult r;
    encode_core(seqs.ptr(), flags, len_bytes_fmt, npolicy,
                do_rc_r1, do_rc_r2, do_random_rc, rng_seed, base_index, r);

    if (!label_col_data.empty()) {
        size_t total_label_bytes = 0;
        for (const auto& col : label_col_data) {
            total_label_bytes += col.size();
        }
        std::string labels_out;
        labels_out.reserve(total_label_bytes);
        for (const auto& col : label_col_data) {
            labels_out.append(col.c_str(), col.size());
        }

        return nb::make_tuple(
            nb::bytes(r.flags_out.data(), r.flags_out.size()),
            nb::bytes(labels_out.data(), labels_out.size()),
            nb::bytes(r.lengths_out.data(), r.lengths_out.size()),
            nb::bytes(r.seqs_out.data(), r.seqs_out.size())
        );
    }

    return nb::make_tuple(
        nb::bytes(r.flags_out.data(), r.flags_out.size()),
        nb::bytes(r.lengths_out.data(), r.lengths_out.size()),
        nb::bytes(r.seqs_out.data(), r.seqs_out.size())
    );
}


// ---------------------------------------------------------------------------
// Columnar decoder
// ---------------------------------------------------------------------------

/**
 * Batch decode columnar block data into records.
 *
 * Returns a list of ``(sequence, is_paired, is_read1, is_read2[, is_rc])``.
 *
 * *with_rc* controls the emitted tuple width.  The reader's default path drops
 * IS_RC immediately, and having the decoder emit the width the caller wants
 * saves rebuilding every tuple in Python just to lose one slot.
 *
 * *restore_strand* undoes the encoder's reverse-complement here rather than one
 * record at a time in Python.  It consumes IS_RC, so it requires ``with_rc``
 * to be false.
 */
nb::object decode_block_columnar(
    nb::bytes flags_data,
    nb::bytes lengths_data,
    nb::bytes seqs_data,
    int len_bytes,
    int count,
    bool with_rc,
    bool restore_strand
) {
    if (with_rc && restore_strand) {
        throw std::invalid_argument(
            "with_rc and restore_strand are mutually exclusive: restore_strand "
            "consumes IS_RC to undo the reverse-complement");
    }
    if (count < 0) {
        throw std::invalid_argument("count must be non-negative");
    }
    if (len_bytes != 1 && len_bytes != 2 && len_bytes != 4) {
        throw std::invalid_argument("len_bytes must be 1, 2, or 4");
    }
    if (flags_data.size() < static_cast<size_t>(count) ||
        lengths_data.size() < static_cast<size_t>(count) * len_bytes) {
        throw std::runtime_error("Block truncated: flags or lengths column too short");
    }

    const uint8_t* flags_ptr = reinterpret_cast<const uint8_t*>(flags_data.c_str());
    const uint8_t* lengths_ptr = reinterpret_cast<const uint8_t*>(lengths_data.c_str());
    const uint8_t* seqs_ptr = reinterpret_cast<const uint8_t*>(seqs_data.c_str());
    const size_t seqs_size = seqs_data.size();

    const Py_ssize_t width = with_rc ? 5 : 4;

    nb::object result = nb::steal(PyList_New(count));
    if (!result.is_valid()) {
        throw nb::python_error();
    }
    PyObject* list = result.ptr();

    size_t seq_offset = 0;
    for (int rec = 0; rec < count; rec++) {
        const uint8_t flag = flags_ptr[rec];
        const size_t seq_len = read_len(lengths_ptr, len_bytes, rec);
        const size_t enc_len = (seq_len + 3) / 4;

        if (seq_offset + enc_len > seqs_size) {
            throw std::runtime_error("Block truncated: cannot read sequence data");
        }

        PyObject* seq = make_seq_object(seqs_ptr + seq_offset, seq_len,
                                        restore_strand && (flag & 0x08));
        if (seq == nullptr) {
            throw nb::python_error();
        }
        seq_offset += enc_len;

        PyObject* t = PyTuple_New(width);
        if (t == nullptr) {
            Py_DECREF(seq);
            throw nb::python_error();
        }
        PyTuple_SET_ITEM(t, 0, seq);  // steals
        set_flag_items(t, flag, 1, with_rc);
        // A tuple of a str and some bools holds nothing that can reference it
        // back, so it can never be part of a cycle and the collector need never
        // walk it.
        PyObject_GC_UnTrack(t);

        PyList_SET_ITEM(list, rec, t);  // steals
    }

    return result;
}


/**
 * Decode a block's sequences only, as a plain list of str.
 *
 * The columnar batch path.  A caller that reads flags off the flags column
 * itself has no use for the per-record tuple, and not building one removes a
 * PyTuple allocation, four bool increfs and a GC untrack per record — all of
 * which is pure overhead once the sequence itself is this cheap to produce.
 */
nb::object decode_block_sequences(
    nb::bytes flags_data,
    nb::bytes lengths_data,
    nb::bytes seqs_data,
    int len_bytes,
    int count,
    bool restore_strand
) {
    if (count < 0) {
        throw std::invalid_argument("count must be non-negative");
    }
    if (len_bytes != 1 && len_bytes != 2 && len_bytes != 4) {
        throw std::invalid_argument("len_bytes must be 1, 2, or 4");
    }
    if (lengths_data.size() < static_cast<size_t>(count) * len_bytes) {
        throw std::runtime_error("Block truncated: lengths column too short");
    }
    // Only consulted when undoing strand normalization.
    if (restore_strand && flags_data.size() < static_cast<size_t>(count)) {
        throw std::runtime_error("Block truncated: flags column too short");
    }

    const uint8_t* flags_ptr = reinterpret_cast<const uint8_t*>(flags_data.c_str());
    const uint8_t* lengths_ptr = reinterpret_cast<const uint8_t*>(lengths_data.c_str());
    const uint8_t* seqs_ptr = reinterpret_cast<const uint8_t*>(seqs_data.c_str());
    const size_t seqs_size = seqs_data.size();

    nb::object result = nb::steal(PyList_New(count));
    if (!result.is_valid()) {
        throw nb::python_error();
    }
    PyObject* list = result.ptr();

    size_t seq_offset = 0;
    for (int rec = 0; rec < count; rec++) {
        const size_t seq_len = read_len(lengths_ptr, len_bytes, rec);
        const size_t enc_len = (seq_len + 3) / 4;
        if (seq_offset + enc_len > seqs_size) {
            throw std::runtime_error("Block truncated: cannot read sequence data");
        }
        PyObject* seq = make_seq_object(
            seqs_ptr + seq_offset, seq_len,
            restore_strand && (flags_ptr[rec] & 0x08));
        if (seq == nullptr) {
            throw nb::python_error();
        }
        seq_offset += enc_len;
        PyList_SET_ITEM(list, rec, seq);  // steals
    }

    return result;
}


// ---------------------------------------------------------------------------
// Label-aware decoder
// ---------------------------------------------------------------------------

/// Build the Python object for one label value.  Returns a new reference.
static inline PyObject* make_label_value(char dtype_code, const uint8_t* p) {
    switch (dtype_code) {
        case 'A':  // uint8 (char)
        case 'C':  // uint8
            return PyLong_FromLong(static_cast<long>(*p));
        case 'c': {  // int8
            int8_t v;
            std::memcpy(&v, p, 1);
            return PyLong_FromLong(static_cast<long>(v));
        }
        case 'S': {  // uint16
            uint16_t v;
            std::memcpy(&v, p, 2);
            return PyLong_FromLong(static_cast<long>(v));
        }
        case 's': {  // int16
            int16_t v;
            std::memcpy(&v, p, 2);
            return PyLong_FromLong(static_cast<long>(v));
        }
        case 'I': {  // uint32
            uint32_t v;
            std::memcpy(&v, p, 4);
            return PyLong_FromUnsignedLong(static_cast<unsigned long>(v));
        }
        case 'i': {  // int32
            int32_t v;
            std::memcpy(&v, p, 4);
            return PyLong_FromLong(static_cast<long>(v));
        }
        case 'f': {  // float32
            float v;
            std::memcpy(&v, p, 4);
            return PyFloat_FromDouble(static_cast<double>(v));
        }
        case 'd': {  // float64
            double v;
            std::memcpy(&v, p, 8);
            return PyFloat_FromDouble(v);
        }
        case 'q': {  // int64
            int64_t v;
            std::memcpy(&v, p, 8);
            return PyLong_FromLongLong(static_cast<long long>(v));
        }
        case 'Q': {  // uint64
            uint64_t v;
            std::memcpy(&v, p, 8);
            return PyLong_FromUnsignedLongLong(static_cast<unsigned long long>(v));
        }
        default:
            throw std::runtime_error(
                std::string("Unknown label dtype code: ") + dtype_code);
    }
}

/**
 * Batch decode columnar block data with label columns.
 *
 * Returns a list of ``(seq, is_paired, is_read1, is_read2[, is_rc], labels)``.
 * *with_rc* and *restore_strand* behave as in :func:`decode_block_columnar`.
 */
nb::object decode_block_labeled(
    nb::bytes flags_data,
    nb::bytes lengths_data,
    nb::bytes seqs_data,
    int len_bytes,
    int count,
    const std::vector<nb::bytes>& label_col_data,
    const std::vector<int>& label_col_sizes,
    const std::string& label_dtype_codes,
    bool with_rc,
    bool restore_strand
) {
    if (with_rc && restore_strand) {
        throw std::invalid_argument(
            "with_rc and restore_strand are mutually exclusive");
    }
    if (count < 0) {
        throw std::invalid_argument("count must be non-negative");
    }
    if (len_bytes != 1 && len_bytes != 2 && len_bytes != 4) {
        throw std::invalid_argument("len_bytes must be 1, 2, or 4");
    }

    const uint8_t* flags_ptr = reinterpret_cast<const uint8_t*>(flags_data.c_str());
    const uint8_t* lengths_ptr = reinterpret_cast<const uint8_t*>(lengths_data.c_str());
    const uint8_t* seqs_ptr = reinterpret_cast<const uint8_t*>(seqs_data.c_str());
    const size_t seqs_size = seqs_data.size();

    const size_t num_labels = label_col_data.size();
    if (label_col_sizes.size() < num_labels ||
        label_dtype_codes.size() < num_labels) {
        throw std::invalid_argument(
            "label_col_sizes and label_dtype_codes must cover every label column");
    }

    std::vector<const uint8_t*> col_ptrs(num_labels);
    for (size_t c = 0; c < num_labels; ++c) {
        col_ptrs[c] = reinterpret_cast<const uint8_t*>(label_col_data[c].c_str());
    }

    const Py_ssize_t width = (with_rc ? 5 : 4) + 1;  // + labels slot
    const Py_ssize_t labels_slot = width - 1;

    nb::object result = nb::steal(PyList_New(count));
    if (!result.is_valid()) {
        throw nb::python_error();
    }
    PyObject* list = result.ptr();

    size_t seq_offset = 0;
    for (int rec = 0; rec < count; rec++) {
        const uint8_t flag = flags_ptr[rec];
        const size_t seq_len = read_len(lengths_ptr, len_bytes, rec);
        const size_t enc_len = (seq_len + 3) / 4;

        if (seq_offset + enc_len > seqs_size) {
            throw std::runtime_error("Block truncated: cannot read sequence data");
        }

        PyObject* seq = make_seq_object(seqs_ptr + seq_offset, seq_len,
                                        restore_strand && (flag & 0x08));
        if (seq == nullptr) {
            throw nb::python_error();
        }
        seq_offset += enc_len;

        PyObject* t = PyTuple_New(width);
        if (t == nullptr) {
            Py_DECREF(seq);
            throw nb::python_error();
        }
        PyTuple_SET_ITEM(t, 0, seq);
        set_flag_items(t, flag, 1, with_rc);

        // The labels tuple is built after the record tuple owns the sequence, so
        // an exception here cannot leak it: `result` owns `t` only once it is in
        // the list, so put the empty labels slot in first.
        PyObject* labels = PyTuple_New(static_cast<Py_ssize_t>(num_labels));
        if (labels == nullptr) {
            Py_DECREF(t);
            throw nb::python_error();
        }
        PyTuple_SET_ITEM(t, labels_slot, labels);
        PyObject_GC_UnTrack(t);
        PyList_SET_ITEM(list, rec, t);  // list owns t (and the seq) from here

        for (size_t c = 0; c < num_labels; ++c) {
            PyObject* v = make_label_value(label_dtype_codes[c],
                                           col_ptrs[c] + rec * label_col_sizes[c]);
            if (v == nullptr) {
                throw nb::python_error();
            }
            PyTuple_SET_ITEM(labels, static_cast<Py_ssize_t>(c), v);
        }
        PyObject_GC_UnTrack(labels);
    }

    return result;
}


// ---------------------------------------------------------------------------
// Fast header label extraction (Optimization 1B)
// ---------------------------------------------------------------------------

// Conversion type codes matching Python side
constexpr int CONV_INT = 0;
constexpr int CONV_FLOAT = 1;
constexpr int CONV_ORD = 2;

/**
 * Fast extraction of SAM-style tags from a FASTQ header line.
 *
 * tag_specs: list of (tag_as_2bytes_string, conv_type) pairs
 * missing_values: pre-packed Python tuple of per-label missing values
 *
 * Returns a Python tuple of label values.
 *
 * The header format is: READNAME<ws>TAG:TYPE:VALUE<ws>TAG:TYPE:VALUE...
 * where <ws> is any whitespace (tab or space).
 */
// ---------------------------------------------------------------------------
// Label value conversion.
//
// These must agree with `extract_labels_from_header` in zna/cli.py EXACTLY -- it is the
// reference, and a label column is not self-describing, so a disagreement is invisible
// until a model trains on it.
//
// The parsing here used to be `v = v * 10 + (c - '0')` with no digit test and no range
// test, which is wrong in three ways at once and silent in all of them:
//
//     NM:i:      -> 0        (reference: ValueError)
//     NM:i:abc   -> 5451     (reference: ValueError)   'a'-'0'=49, 'b'-'0'=50, ...
//     NM:i:3.7   -> 287      (reference: ValueError)   '.'-'0' = -2
//     NM:i:<2^63 -> wrapped  (reference: the exact integer) -- and signed overflow is
//                            undefined behaviour, not merely a wrong answer.
//
// So the fast loop stays for the case it is good at -- a short run of ASCII digits --
// and everything else is handed to CPython's own `int()`/`float()`. That way the value
// AND the exception are the reference's by construction, rather than by a
// reimplementation that has to be kept in step by hand.
// ---------------------------------------------------------------------------

/// `int(bytes)`, including the ValueError it raises. Used for anything the fast path
/// declines: empty, signed-and-long, over-long, or not a plain decimal.
static nb::object int_via_python(const char* p, size_t n) {
    const std::string buf(p, n);
    char* end = nullptr;
    PyObject* o = PyLong_FromString(buf.c_str(), &end, 10);
    if (o == nullptr) throw nb::python_error();
    // PyLong_FromString tolerates trailing junk that `int()` does not.
    if (end != buf.c_str() + buf.size()) {
        Py_DECREF(o);
        PyErr_Format(PyExc_ValueError,
                     "invalid literal for int() with base 10: b'%s'", buf.c_str());
        throw nb::python_error();
    }
    return nb::steal<nb::object>(o);
}

/// An integer label value, parsed as Python's `int()` would.
static nb::object parse_int_value(const char* p, size_t n) {
    // Fast path: optional sign, then 1..18 digits -- which always fits in int64.
    const size_t i = (n > 0 && (p[0] == '-' || p[0] == '+')) ? 1u : 0u;
    if (n > i && n - i <= 18) {
        long long v = 0;
        size_t k = i;
        for (; k < n; ++k) {
            const unsigned d = static_cast<unsigned>(static_cast<unsigned char>(p[k])) - '0';
            if (d > 9) break;
            v = v * 10 + static_cast<long long>(d);
        }
        if (k == n) return nb::int_(p[0] == '-' ? -v : v);
    }
    return int_via_python(p, n);
}

/// A float label value, parsed as Python's `float()` would.
///
/// `strtod` alone will not do: it stops at the first byte it cannot use, so "3.7x"
/// silently became 3.7 and "" became 0.0. It is also locale-sensitive about the decimal
/// point, where `float()` never is -- which this handles for free, because a locale that
/// rejects '.' leaves `end` short of the value and falls through to CPython.
static nb::object parse_float_value(const char* p, size_t n) {
    // `strtod` also accepts C99 hex float literals ("0x10" -> 16.0, whole string
    // consumed) where `float()` raises. That is the one form it over-accepts, so send it
    // to CPython rather than trying to out-guess it.
    const size_t sign = (n > 0 && (p[0] == '-' || p[0] == '+')) ? 1u : 0u;
    const bool hex = n >= sign + 2 && p[sign] == '0'
                     && (p[sign + 1] == 'x' || p[sign + 1] == 'X');
    if (n > 0 && !hex) {
        const std::string buf(p, n);
        char* end = nullptr;
        const double v = std::strtod(buf.c_str(), &end);
        if (end == buf.c_str() + buf.size()) return nb::float_(v);
    }
    PyObject* s = PyUnicode_FromStringAndSize(p, static_cast<Py_ssize_t>(n));
    if (s == nullptr) throw nb::python_error();
    PyObject* o = PyFloat_FromString(s);
    Py_DECREF(s);
    if (o == nullptr) throw nb::python_error();
    return nb::steal<nb::object>(o);
}

nb::tuple extract_labels_fast(
    nb::bytes header,
    const std::vector<std::pair<nb::bytes, int>>& tag_specs,
    int num_labels,
    nb::tuple missing_values
) {
    const char* data = header.c_str();
    const size_t data_len = header.size();

    // Build compact tag lookup with variable-length tag strings
    struct TagSpec {
        const char* tag_str;
        size_t tag_len;
        int label_idx;
        int conv_type;
    };
    std::vector<TagSpec> specs;
    specs.reserve(num_labels);
    for (int i = 0; i < num_labels; ++i) {
        const char* tag_str = tag_specs[i].first.c_str();
        size_t tag_len = tag_specs[i].first.size();
        specs.push_back({tag_str, tag_len, i, tag_specs[i].second});
    }

    // Allocate result — start with missing values
    nb::tuple result = nb::steal<nb::tuple>(
        PyTuple_New(static_cast<Py_ssize_t>(num_labels))
    );
    // Initialize with missing values
    for (int i = 0; i < num_labels; ++i) {
        nb::object val = nb::borrow(missing_values[i]);
        PyTuple_SET_ITEM(result.ptr(), i, val.release().ptr());
    }

    int remaining = num_labels;

    // Skip read name — find first whitespace
    size_t pos = 0;
    while (pos < data_len && data[pos] != ' ' && data[pos] != '\t' && data[pos] != '\n') {
        pos++;
    }

    // Process fields
    while (pos < data_len && remaining > 0) {
        // Skip whitespace
        while (pos < data_len && (data[pos] == ' ' || data[pos] == '\t')) {
            pos++;
        }
        if (pos >= data_len) break;

        // Find end of field
        size_t field_start = pos;
        while (pos < data_len && data[pos] != ' ' && data[pos] != '\t' && data[pos] != '\n') {
            pos++;
        }
        size_t field_len = pos - field_start;

        // Check KEY:TYPE:VALUE format — find first colon to extract key
        if (field_len < 4) continue;  // minimum K:X:V
        const char* field = data + field_start;

        // Find the first colon (end of key)
        size_t colon1 = 0;
        while (colon1 < field_len && field[colon1] != ':') colon1++;
        if (colon1 < 1 || colon1 + 2 >= field_len) continue;
        // TYPE is single char after first colon, second colon follows
        if (field[colon1 + 2] != ':') continue;

        // Look up tag by comparing key portion
        for (const auto& spec : specs) {
            if (spec.tag_len == colon1 && std::memcmp(field, spec.tag_str, colon1) == 0) {
                const char* val_start = field + colon1 + 3;
                size_t val_len = field_len - colon1 - 3;

                nb::object py_val;
                if (spec.conv_type == CONV_INT) {
                    py_val = parse_int_value(val_start, val_len);
                } else if (spec.conv_type == CONV_FLOAT) {
                    py_val = parse_float_value(val_start, val_len);
                } else { // CONV_ORD — a single byte is its ordinal, else an integer
                    if (val_len == 1) {
                        py_val = nb::int_(static_cast<long>(static_cast<uint8_t>(val_start[0])));
                    } else {
                        py_val = parse_int_value(val_start, val_len);
                    }
                }

                // Replace missing value in result tuple
                PyObject* old = PyTuple_GET_ITEM(result.ptr(), spec.label_idx);
                Py_INCREF(py_val.ptr());
                PyTuple_SET_ITEM(result.ptr(), spec.label_idx, py_val.ptr());
                Py_XDECREF(old);
                remaining--;
                break;  // found this tag, move to next field
            }
        }
    }

    return result;
}


NB_MODULE(_accel, m) {
    m.doc() = "ZNA accelerated encode/decode functions";

    m.def("encode_sequence", &encode_sequence,
          nb::arg("seq"),
          "Encode DNA sequence to 2-bit packed bytes.\n\n"
          "Each base (A, C, G, T) is encoded as 2 bits.\n"
          "Four bases are packed into each byte.");

    m.def("encode_block", &encode_block,
          nb::arg("seqs"),
          nb::arg("flags"),
          nb::arg("len_bytes_fmt"),
          nb::arg("npolicy"),
          nb::arg("do_rc_r1"),
          nb::arg("do_rc_r2"),
          nb::arg("do_random_rc") = false,
          nb::arg("rng_seed") = 0,
          nb::arg("base_index") = 0,
          "Batch encode sequences into columnar streams (no labels).");

    m.def("encode_block_labeled", &encode_block_labeled,
          nb::arg("seqs"),
          nb::arg("flags"),
          nb::arg("len_bytes_fmt"),
          nb::arg("npolicy"),
          nb::arg("do_rc_r1"),
          nb::arg("do_rc_r2"),
          nb::arg("do_random_rc"),
          nb::arg("label_col_data"),
          nb::arg("label_col_sizes"), nb::arg("rng_seed") = 0,
          nb::arg("base_index") = 0,
          "Batch encode sequences with pre-packed label columns.\n"
          "Returns 4-tuple (flags, labels, lengths, seqs) when labels present,\n"
          "or 3-tuple when label_col_data is empty.");

    m.def("decode_block_records", &decode_block_records,
          nb::arg("block_data"),
          nb::arg("len_bytes"),
          nb::arg("count"),
          "Decode all records from a block at once.");

    m.def("decode_block_columnar", &decode_block_columnar,
          nb::arg("flags_data"),
          nb::arg("lengths_data"),
          nb::arg("seqs_data"),
          nb::arg("len_bytes"),
          nb::arg("count"),
          nb::arg("with_rc") = true,
          nb::arg("restore_strand") = false,
          "Decode columnar block streams into records (no labels).\n"
          "with_rc=False emits 4-tuples; restore_strand=True undoes the\n"
          "encoder's reverse-complement in C++ and also emits 4-tuples.");

    m.def("decode_block", &decode_block_columnar,
          nb::arg("flags_data"),
          nb::arg("lengths_data"),
          nb::arg("seqs_data"),
          nb::arg("len_bytes"),
          nb::arg("count"),
          nb::arg("with_rc") = true,
          nb::arg("restore_strand") = false,
          "Decode block streams into records (no labels).");

    m.def("decode_block_sequences", &decode_block_sequences,
          nb::arg("flags_data"),
          nb::arg("lengths_data"),
          nb::arg("seqs_data"),
          nb::arg("len_bytes"),
          nb::arg("count"),
          nb::arg("restore_strand") = false,
          "Decode a block's sequences only, as a list of str.\n"
          "For callers that read the flags column themselves.");

    m.def("decode_block_labeled", &decode_block_labeled,
          nb::arg("flags_data"),
          nb::arg("lengths_data"),
          nb::arg("seqs_data"),
          nb::arg("len_bytes"),
          nb::arg("count"),
          nb::arg("label_col_data"),
          nb::arg("label_col_sizes"),
          nb::arg("label_dtype_codes"),
          nb::arg("with_rc") = true,
          nb::arg("restore_strand") = false,
          "Decode columnar block streams with label columns.\n"
          "Returns list of (seq, is_paired, is_read1, is_read2[, is_rc], labels).");

    m.def("extract_labels_fast", &extract_labels_fast,
          nb::arg("header"),
          nb::arg("tag_specs"),
          nb::arg("num_labels"),
          nb::arg("missing_values"),
          "Fast SAM-tag extraction from FASTQ header bytes.\n"
          "tag_specs: list of (tag_bytes, conv_type) pairs.\n"
          "missing_values: tuple of default values per label.");

    m.def("reverse_complement", &reverse_complement,
          nb::arg("seq"),
          "Compute reverse complement of a DNA sequence.");
}
