/**
 * ZNA Accelerated Encode/Decode Functions
 * 
 * High-performance C++ implementation of DNA sequence encoding/decoding.
 * Uses nanobind for Python bindings.
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

// Lookup table for encoding (char -> 2-bit value)
alignas(64) static uint8_t ENCODE_TABLE[256];

// Lookup table for decoding (2-bit value -> char)
static constexpr char DECODE_CHARS[4] = {'A', 'C', 'G', 'T'};

// Complement table for reverse complement (A<->T, C<->G)
static constexpr char COMPLEMENT_TABLE[4] = {'T', 'G', 'C', 'A'};

// Initialize encoding table at module load
static void init_encode_table() noexcept {
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
}

// Static initializer
static struct TableInitializer {
    TableInitializer() { init_encode_table(); }
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
    out.reserve(out_len);
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
    std::string out;
    out.reserve(seq_len);
    out.resize(seq_len);
    char* out_ptr = out.data();
    
    size_t j = 0;
    const size_t full_bytes = seq_len / 4;
    
    // Decode full bytes (4 bases each)
    for (size_t i = 0; i < full_bytes; i++) {
        const uint8_t byte_val = data[i];
        out_ptr[j++] = DECODE_CHARS[(byte_val >> 6) & 0x03];
        out_ptr[j++] = DECODE_CHARS[(byte_val >> 4) & 0x03];
        out_ptr[j++] = DECODE_CHARS[(byte_val >> 2) & 0x03];
        out_ptr[j++] = DECODE_CHARS[byte_val & 0x03];
    }
    
    // Handle remaining bases (if any)
    if (j < seq_len && full_bytes < data_len) {
        const uint8_t byte_val = data[full_bytes];
        const size_t remaining = seq_len - j;
        if (remaining >= 1) out_ptr[j++] = DECODE_CHARS[(byte_val >> 6) & 0x03];
        if (remaining >= 2) out_ptr[j++] = DECODE_CHARS[(byte_val >> 4) & 0x03];
        if (remaining >= 3) out_ptr[j++] = DECODE_CHARS[(byte_val >> 2) & 0x03];
    }
    
    return out;
}


// Record tuple type: (sequence, is_paired, is_read1, is_read2, is_rc)
using Record = std::tuple<std::string, bool, bool, bool, bool>;

/**
 * Decode all records from a block at once.
 * 
 * This processes an entire block of records in C++, avoiding Python
 * overhead for each record.
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
 * Uses optimized in-place reversal and complementation.
 * A<->T, C<->G
 */
std::string reverse_complement(const std::string& seq) {
    const size_t len = seq.size();
    std::string result;
    result.resize(len);
    
    // Work from both ends, complement and reverse in one pass
    const char* src = seq.data();
    char* dst = result.data();
    
    for (size_t i = 0; i < len; i++) {
        uint8_t base = ENCODE_TABLE[static_cast<uint8_t>(src[len - 1 - i])];
        if (base == INVALID) {
            // For unknown bases, just reverse without complement
            dst[i] = src[len - 1 - i];
        } else {
            dst[i] = COMPLEMENT_TABLE[base];
        }
    }
    
    return result;
}


/**
 * BATCH ENCODER: The single most important performance optimization.
 * 
 * Takes a list of sequences and flags, builds the three columnar streams
 * (flags, lengths, sequences) entirely in C++, eliminating all per-record
 * Python overhead.
 *
 * For strand-specific libraries, reverse complement is applied in C++ to
 * designated reads before encoding.
 *
 * Input: List of strings, flags, len format, N policy, strand norm masks
 * Output: Tuple(flags_bytes, lengths_bytes, sequences_bytes)
 */
nb::tuple encode_block(
    const std::vector<std::string>& seqs,
    const std::vector<uint8_t>& flags,
    int len_bytes_fmt,
    const std::string& npolicy,
    bool do_rc_r1,
    bool do_rc_r2,
    bool do_random_rc
) {
    const size_t count = seqs.size();
    if (count != flags.size()) {
        throw std::invalid_argument("seqs and flags must have the same length");
    }
    
    // IS_RC flag bit
    constexpr uint8_t IS_RC_BIT = 0x08;
    constexpr uint8_t IS_READ1 = 0x01;
    constexpr uint8_t IS_READ2 = 0x02;
    constexpr uint8_t IS_PAIRED = 0x04;
    
    // 1. Prepare output buffers
    std::string flags_out;
    flags_out.resize(count);
    // Copy flags so we can mutate them
    for (size_t i = 0; i < count; ++i) {
        flags_out[i] = static_cast<char>(flags[i]);
    }
    
    std::string lengths_out;
    lengths_out.resize(count * len_bytes_fmt);
    
    // Heuristic reserve for sequences (~40 encoded bytes per 150bp read)
    std::vector<uint8_t> seqs_out;
    seqs_out.reserve(count * 40);
    
    // N-handling logic setup
    uint8_t n_replace_val = 0; // Default: A
    bool use_random_n = false;
    bool has_npolicy = !npolicy.empty();
    
    if (npolicy == "C" || npolicy == "c") n_replace_val = BASE_C;
    else if (npolicy == "G" || npolicy == "g") n_replace_val = BASE_G;
    else if (npolicy == "T" || npolicy == "t") n_replace_val = BASE_T;
    else if (npolicy == "random") use_random_n = true;
    // else: A (default) or empty (will throw on N)
    
    // Simple PRNG state for random N replacement (xorshift32)
    uint32_t rng_state = 0xDEADBEEF;
    
    // Max sequence length check
    size_t max_len = 0;
    if (len_bytes_fmt == 1) max_len = 255;
    else if (len_bytes_fmt == 2) max_len = 65535;
    else max_len = 4294967295UL;
    
    // 2. Main loop (entirely in C++)
    for (size_t i = 0; i < count; ++i) {
        // --- A. Flags ---
        uint8_t flag = static_cast<uint8_t>(flags_out[i]);
        
        // Determine if this read needs reverse complement for strand normalization
        bool is_read1 = (flag & IS_READ1) != 0;
        bool is_read2 = (flag & IS_READ2) != 0;
        bool is_paired = (flag & IS_PAIRED) != 0;
        bool needs_rc = false;
        
        if (do_random_rc) {
            // Unstranded random normalization
            if (is_paired && is_read1 && (i + 1 < count) && (flags[i + 1] & IS_PAIRED)) {
                // Paired R1+R2: flip a coin to decide which read to RC
                rng_state ^= rng_state << 13;
                rng_state ^= rng_state >> 17;
                rng_state ^= rng_state << 5;
                bool coin = rng_state & 1;
                if (coin) {
                    // RC R1
                    needs_rc = true;
                    flags_out[i] = static_cast<char>(flag | IS_RC_BIT);
                } else {
                    // RC R2
                    flags_out[i + 1] = static_cast<char>(flags[i + 1] | IS_RC_BIT);
                }
            } else if (!(is_paired && is_read2)) {
                // Unpaired/SE or lone R1 (no following R2): random RC
                rng_state ^= rng_state << 13;
                rng_state ^= rng_state >> 17;
                rng_state ^= rng_state << 5;
                if (rng_state & 1) {
                    needs_rc = true;
                    flags_out[i] = static_cast<char>(flag | IS_RC_BIT);
                }
            }
            // Check if this is an R2 that was already marked for RC by its R1 iteration
            if (!needs_rc && (static_cast<uint8_t>(flags_out[i]) & IS_RC_BIT)) {
                needs_rc = true;
                flag = static_cast<uint8_t>(flags_out[i]);
            }
        } else {
            // Deterministic stranded normalization
            needs_rc = (do_rc_r1 && is_read1) || (do_rc_r2 && is_read2);
            if (needs_rc) {
                flags_out[i] = static_cast<char>(flag | IS_RC_BIT);
            }
        }
        
        // Get sequence (possibly reverse complemented)
        const std::string& orig_seq = seqs[i];
        std::string rc_seq;
        const std::string* seq_ptr;
        
        if (needs_rc) {
            rc_seq = reverse_complement(orig_seq);
            seq_ptr = &rc_seq;
        } else {
            seq_ptr = &orig_seq;
        }
        
        const std::string& seq = *seq_ptr;
        size_t slen = seq.size();
        
        if (slen > max_len) {
            throw std::invalid_argument(
                "Sequence length " + std::to_string(slen) + 
                " exceeds maximum " + std::to_string(max_len)
            );
        }
        
        // --- B. Lengths (little-endian) ---
        if (len_bytes_fmt == 1) {
            lengths_out[i] = static_cast<char>(slen);
        } else if (len_bytes_fmt == 2) {
            uint16_t s = static_cast<uint16_t>(slen);
            std::memcpy(&lengths_out[i * 2], &s, 2);
        } else {
            uint32_t s = static_cast<uint32_t>(slen);
            std::memcpy(&lengths_out[i * 4], &s, 4);
        }
        
        // --- C. Encode sequence (2-bit packing with N handling) ---
        const uint8_t* seq_bytes = reinterpret_cast<const uint8_t*>(seq.data());
        size_t s_idx = 0;
        
        // Fast packing loop: 4 bases at a time
        while (s_idx + 4 <= slen) {
            uint8_t b0 = ENCODE_TABLE[seq_bytes[s_idx]];
            uint8_t b1 = ENCODE_TABLE[seq_bytes[s_idx + 1]];
            uint8_t b2 = ENCODE_TABLE[seq_bytes[s_idx + 2]];
            uint8_t b3 = ENCODE_TABLE[seq_bytes[s_idx + 3]];
            
            // Handle N/invalid characters
            if (b0 == INVALID) {
                if (!has_npolicy) throw std::invalid_argument("Invalid character in sequence at position " + std::to_string(s_idx));
                if (use_random_n) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 17; rng_state ^= rng_state << 5; b0 = rng_state & 3; } 
                else b0 = n_replace_val;
            }
            if (b1 == INVALID) {
                if (!has_npolicy) throw std::invalid_argument("Invalid character in sequence");
                if (use_random_n) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 17; rng_state ^= rng_state << 5; b1 = rng_state & 3; }
                else b1 = n_replace_val;
            }
            if (b2 == INVALID) {
                if (!has_npolicy) throw std::invalid_argument("Invalid character in sequence");
                if (use_random_n) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 17; rng_state ^= rng_state << 5; b2 = rng_state & 3; }
                else b2 = n_replace_val;
            }
            if (b3 == INVALID) {
                if (!has_npolicy) throw std::invalid_argument("Invalid character in sequence");
                if (use_random_n) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 17; rng_state ^= rng_state << 5; b3 = rng_state & 3; }
                else b3 = n_replace_val;
            }
            
            seqs_out.push_back((b0 << 6) | (b1 << 4) | (b2 << 2) | b3);
            s_idx += 4;
        }
        
        // Tail: handle remaining 1-3 bases
        if (s_idx < slen) {
            uint8_t val = 0;
            int shift = 6;
            while (s_idx < slen) {
                uint8_t b = ENCODE_TABLE[seq_bytes[s_idx]];
                if (b == INVALID) {
                    if (!has_npolicy) throw std::invalid_argument("Invalid character in sequence");
                    if (use_random_n) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 17; rng_state ^= rng_state << 5; b = rng_state & 3; }
                    else b = n_replace_val;
                }
                val |= (b << shift);
                s_idx++;
                shift -= 2;
            }
            seqs_out.push_back(val);
        }
    }
    
    // Zero-copy return to Python via bytes
    return nb::make_tuple(
        nb::bytes(flags_out.data(), flags_out.size()),
        nb::bytes(lengths_out.data(), lengths_out.size()),
        nb::bytes(seqs_out.data(), seqs_out.size())
    );
}


/**
 * Batch decode columnar block data into records.
 *
 * Takes the three separate columnar streams (flags, lengths, sequences)
 * and returns a list of (sequence, is_paired, is_read1, is_read2, is_rc) tuples.
 * This avoids per-record Python overhead in the decode loop.
 */
std::vector<Record> decode_block_columnar(
    nb::bytes flags_data,
    nb::bytes lengths_data,
    nb::bytes seqs_data,
    int len_bytes,
    int count
) {
    const uint8_t* flags_ptr = reinterpret_cast<const uint8_t*>(flags_data.c_str());
    const uint8_t* lengths_ptr = reinterpret_cast<const uint8_t*>(lengths_data.c_str());
    const uint8_t* seqs_ptr = reinterpret_cast<const uint8_t*>(seqs_data.c_str());
    const size_t seqs_size = seqs_data.size();
    
    std::vector<Record> results;
    results.reserve(count);
    
    size_t seq_offset = 0;
    
    for (int rec = 0; rec < count; rec++) {
        // Read flags
        uint8_t flag = flags_ptr[rec];
        bool is_read1 = (flag & 1) != 0;
        bool is_read2 = (flag & 2) != 0;
        bool is_paired = (flag & 4) != 0;
        
        // Read sequence length
        size_t seq_len = 0;
        if (len_bytes == 1) {
            seq_len = lengths_ptr[rec];
        } else if (len_bytes == 2) {
            seq_len = read_u16_le(lengths_ptr + rec * 2);
        } else {
            seq_len = read_u32_le(lengths_ptr + rec * 4);
        }
        
        // Calculate encoded length
        size_t enc_len = (seq_len + 3) / 4;
        
        if (seq_offset + enc_len > seqs_size) {
            throw std::runtime_error("Block truncated: cannot read sequence data");
        }
        
        // Decode sequence
        std::string seq = decode_sequence(seqs_ptr + seq_offset, enc_len, seq_len);
        seq_offset += enc_len;
        
        results.emplace_back(std::move(seq), is_paired, is_read1, is_read2,
                             (flag & 0x08) != 0);  // is_rc
    }
    
    return results;
}


// ---------------------------------------------------------------------------
// Label-aware encode_block
// ---------------------------------------------------------------------------

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
    const std::vector<std::string>& seqs,
    const std::vector<uint8_t>& flags,
    int len_bytes_fmt,
    const std::string& npolicy,
    bool do_rc_r1,
    bool do_rc_r2,
    bool do_random_rc,
    const std::vector<nb::bytes>& label_col_data,
    const std::vector<int>& label_col_sizes
) {
    const size_t count = seqs.size();
    if (count != flags.size()) {
        throw std::invalid_argument("seqs and flags must have the same length");
    }

    constexpr uint8_t IS_RC_BIT = 0x08;
    constexpr uint8_t IS_READ1 = 0x01;
    constexpr uint8_t IS_READ2 = 0x02;
    constexpr uint8_t IS_PAIRED = 0x04;

    // 1. Prepare output buffers
    std::string flags_out;
    flags_out.resize(count);
    for (size_t i = 0; i < count; ++i) {
        flags_out[i] = static_cast<char>(flags[i]);
    }

    std::string lengths_out;
    lengths_out.resize(count * len_bytes_fmt);

    std::vector<uint8_t> seqs_out;
    seqs_out.reserve(count * 40);

    // N-handling setup
    uint8_t n_replace_val = 0;
    bool use_random_n = false;
    bool has_npolicy = !npolicy.empty();

    if (npolicy == "C" || npolicy == "c") n_replace_val = BASE_C;
    else if (npolicy == "G" || npolicy == "g") n_replace_val = BASE_G;
    else if (npolicy == "T" || npolicy == "t") n_replace_val = BASE_T;
    else if (npolicy == "random") use_random_n = true;

    uint32_t rng_state = 0xDEADBEEF;

    size_t max_len = 0;
    if (len_bytes_fmt == 1) max_len = 255;
    else if (len_bytes_fmt == 2) max_len = 65535;
    else max_len = 4294967295UL;

    // 2. Main encode loop (identical to original encode_block)
    for (size_t i = 0; i < count; ++i) {
        uint8_t flag = static_cast<uint8_t>(flags_out[i]);

        bool is_read1 = (flag & IS_READ1) != 0;
        bool is_read2 = (flag & IS_READ2) != 0;
        bool is_paired = (flag & IS_PAIRED) != 0;
        bool needs_rc = false;

        if (do_random_rc) {
            if (is_paired && is_read1 && (i + 1 < count) && (flags[i + 1] & IS_PAIRED)) {
                rng_state ^= rng_state << 13;
                rng_state ^= rng_state >> 17;
                rng_state ^= rng_state << 5;
                bool coin = rng_state & 1;
                if (coin) {
                    needs_rc = true;
                    flags_out[i] = static_cast<char>(flag | IS_RC_BIT);
                } else {
                    flags_out[i + 1] = static_cast<char>(flags[i + 1] | IS_RC_BIT);
                }
            } else if (!(is_paired && is_read2)) {
                rng_state ^= rng_state << 13;
                rng_state ^= rng_state >> 17;
                rng_state ^= rng_state << 5;
                if (rng_state & 1) {
                    needs_rc = true;
                    flags_out[i] = static_cast<char>(flag | IS_RC_BIT);
                }
            }
            if (!needs_rc && (static_cast<uint8_t>(flags_out[i]) & IS_RC_BIT)) {
                needs_rc = true;
                flag = static_cast<uint8_t>(flags_out[i]);
            }
        } else {
            needs_rc = (do_rc_r1 && is_read1) || (do_rc_r2 && is_read2);
            if (needs_rc) {
                flags_out[i] = static_cast<char>(flag | IS_RC_BIT);
            }
        }

        const std::string& orig_seq = seqs[i];
        std::string rc_seq;
        const std::string* seq_ptr;

        if (needs_rc) {
            rc_seq = reverse_complement(orig_seq);
            seq_ptr = &rc_seq;
        } else {
            seq_ptr = &orig_seq;
        }

        const std::string& seq = *seq_ptr;
        size_t slen = seq.size();

        if (slen > max_len) {
            throw std::invalid_argument(
                "Sequence length " + std::to_string(slen) +
                " exceeds maximum " + std::to_string(max_len)
            );
        }

        // Lengths
        if (len_bytes_fmt == 1) {
            lengths_out[i] = static_cast<char>(slen);
        } else if (len_bytes_fmt == 2) {
            uint16_t s = static_cast<uint16_t>(slen);
            std::memcpy(&lengths_out[i * 2], &s, 2);
        } else {
            uint32_t s = static_cast<uint32_t>(slen);
            std::memcpy(&lengths_out[i * 4], &s, 4);
        }

        // Encode sequence
        const uint8_t* seq_bytes = reinterpret_cast<const uint8_t*>(seq.data());
        size_t s_idx = 0;

        while (s_idx + 4 <= slen) {
            uint8_t b0 = ENCODE_TABLE[seq_bytes[s_idx]];
            uint8_t b1 = ENCODE_TABLE[seq_bytes[s_idx + 1]];
            uint8_t b2 = ENCODE_TABLE[seq_bytes[s_idx + 2]];
            uint8_t b3 = ENCODE_TABLE[seq_bytes[s_idx + 3]];

            if (b0 == INVALID) {
                if (!has_npolicy) throw std::invalid_argument("Invalid character in sequence at position " + std::to_string(s_idx));
                if (use_random_n) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 17; rng_state ^= rng_state << 5; b0 = rng_state & 3; }
                else b0 = n_replace_val;
            }
            if (b1 == INVALID) {
                if (!has_npolicy) throw std::invalid_argument("Invalid character in sequence");
                if (use_random_n) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 17; rng_state ^= rng_state << 5; b1 = rng_state & 3; }
                else b1 = n_replace_val;
            }
            if (b2 == INVALID) {
                if (!has_npolicy) throw std::invalid_argument("Invalid character in sequence");
                if (use_random_n) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 17; rng_state ^= rng_state << 5; b2 = rng_state & 3; }
                else b2 = n_replace_val;
            }
            if (b3 == INVALID) {
                if (!has_npolicy) throw std::invalid_argument("Invalid character in sequence");
                if (use_random_n) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 17; rng_state ^= rng_state << 5; b3 = rng_state & 3; }
                else b3 = n_replace_val;
            }

            seqs_out.push_back((b0 << 6) | (b1 << 4) | (b2 << 2) | b3);
            s_idx += 4;
        }

        if (s_idx < slen) {
            uint8_t val = 0;
            int shift = 6;
            while (s_idx < slen) {
                uint8_t b = ENCODE_TABLE[seq_bytes[s_idx]];
                if (b == INVALID) {
                    if (!has_npolicy) throw std::invalid_argument("Invalid character in sequence");
                    if (use_random_n) { rng_state ^= rng_state << 13; rng_state ^= rng_state >> 17; rng_state ^= rng_state << 5; b = rng_state & 3; }
                    else b = n_replace_val;
                }
                val |= (b << shift);
                s_idx++;
                shift -= 2;
            }
            seqs_out.push_back(val);
        }
    }

    // 3. Build labels output by concatenating pre-packed column bytes
    bool has_labels = !label_col_data.empty();
    if (has_labels) {
        // Concatenate all label column bytes
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
            nb::bytes(flags_out.data(), flags_out.size()),
            nb::bytes(labels_out.data(), labels_out.size()),
            nb::bytes(lengths_out.data(), lengths_out.size()),
            nb::bytes(seqs_out.data(), seqs_out.size())
        );
    }

    return nb::make_tuple(
        nb::bytes(flags_out.data(), flags_out.size()),
        nb::bytes(lengths_out.data(), lengths_out.size()),
        nb::bytes(seqs_out.data(), seqs_out.size())
    );
}


// ---------------------------------------------------------------------------
// Label-aware decode_block
// ---------------------------------------------------------------------------

// Record with labels: (sequence, is_paired, is_read1, is_read2, is_rc, labels_tuple)
using LabeledRecord = std::tuple<std::string, bool, bool, bool, bool, nb::tuple>;

/**
 * Batch decode columnar block data with label columns.
 *
 * label_col_data: list of bytes objects, one per label column (already sliced)
 * label_col_sizes: bytes-per-value for each label column
 * label_dtype_codes: single-char dtype code for each label column
 *
 * Returns list of (seq, is_paired, is_read1, is_read2, is_rc, (label_values...)) tuples.
 */
std::vector<LabeledRecord> decode_block_labeled(
    nb::bytes flags_data,
    nb::bytes lengths_data,
    nb::bytes seqs_data,
    int len_bytes,
    int count,
    const std::vector<nb::bytes>& label_col_data,
    const std::vector<int>& label_col_sizes,
    const std::string& label_dtype_codes
) {
    const uint8_t* flags_ptr = reinterpret_cast<const uint8_t*>(flags_data.c_str());
    const uint8_t* lengths_ptr = reinterpret_cast<const uint8_t*>(lengths_data.c_str());
    const uint8_t* seqs_ptr = reinterpret_cast<const uint8_t*>(seqs_data.c_str());
    const size_t seqs_size = seqs_data.size();

    const size_t num_labels = label_col_data.size();

    // Pre-parse label column pointers
    std::vector<const uint8_t*> col_ptrs(num_labels);
    for (size_t c = 0; c < num_labels; ++c) {
        col_ptrs[c] = reinterpret_cast<const uint8_t*>(label_col_data[c].c_str());
    }

    std::vector<LabeledRecord> results;
    results.reserve(count);

    size_t seq_offset = 0;

    for (int rec = 0; rec < count; rec++) {
        // Read flags
        uint8_t flag = flags_ptr[rec];
        bool is_read1 = (flag & 1) != 0;
        bool is_read2 = (flag & 2) != 0;
        bool is_paired = (flag & 4) != 0;

        // Read sequence length
        size_t seq_len = 0;
        if (len_bytes == 1) {
            seq_len = lengths_ptr[rec];
        } else if (len_bytes == 2) {
            seq_len = read_u16_le(lengths_ptr + rec * 2);
        } else {
            seq_len = read_u32_le(lengths_ptr + rec * 4);
        }

        size_t enc_len = (seq_len + 3) / 4;
        if (seq_offset + enc_len > seqs_size) {
            throw std::runtime_error("Block truncated: cannot read sequence data");
        }

        std::string seq = decode_sequence(seqs_ptr + seq_offset, enc_len, seq_len);
        seq_offset += enc_len;

        // Build label tuple
        nb::tuple labels = nb::steal<nb::tuple>(PyTuple_New(static_cast<Py_ssize_t>(num_labels)));

        for (size_t c = 0; c < num_labels; ++c) {
            const uint8_t* col_ptr = col_ptrs[c];
            int col_size = label_col_sizes[c];
            const uint8_t* val_ptr = col_ptr + rec * col_size;

            nb::object py_val;
            char dtype_code = label_dtype_codes[c];

            switch (dtype_code) {
                case 'A': // uint8 (char)
                case 'C': // uint8
                    py_val = nb::int_(static_cast<long>(*val_ptr));
                    break;
                case 'c': { // int8
                    int8_t v;
                    std::memcpy(&v, val_ptr, 1);
                    py_val = nb::int_(static_cast<long>(v));
                    break;
                }
                case 'S': { // uint16
                    uint16_t v;
                    std::memcpy(&v, val_ptr, 2);
                    py_val = nb::int_(static_cast<long>(v));
                    break;
                }
                case 's': { // int16
                    int16_t v;
                    std::memcpy(&v, val_ptr, 2);
                    py_val = nb::int_(static_cast<long>(v));
                    break;
                }
                case 'I': { // uint32
                    uint32_t v;
                    std::memcpy(&v, val_ptr, 4);
                    py_val = nb::int_(static_cast<unsigned long>(v));
                    break;
                }
                case 'i': { // int32
                    int32_t v;
                    std::memcpy(&v, val_ptr, 4);
                    py_val = nb::int_(static_cast<long>(v));
                    break;
                }
                case 'f': { // float32
                    float v;
                    std::memcpy(&v, val_ptr, 4);
                    py_val = nb::float_(static_cast<double>(v));
                    break;
                }
                case 'd': { // float64
                    double v;
                    std::memcpy(&v, val_ptr, 8);
                    py_val = nb::float_(v);
                    break;
                }
                case 'q': { // int64
                    int64_t v;
                    std::memcpy(&v, val_ptr, 8);
                    py_val = nb::int_(static_cast<long long>(v));
                    break;
                }
                case 'Q': { // uint64
                    uint64_t v;
                    std::memcpy(&v, val_ptr, 8);
                    py_val = nb::int_(static_cast<unsigned long long>(v));
                    break;
                }
                default:
                    throw std::runtime_error(
                        std::string("Unknown label dtype code: ") + dtype_code
                    );
            }

            // Use PyTuple_SET_ITEM (steals reference)
            PyTuple_SET_ITEM(labels.ptr(), static_cast<Py_ssize_t>(c), py_val.release().ptr());
        }

        results.emplace_back(
            std::move(seq), is_paired, is_read1, is_read2,
            (flag & 0x08) != 0, std::move(labels)
        );
    }

    return results;
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
                    // Parse integer from bytes — no allocation
                    long long v = 0;
                    bool negative = false;
                    size_t vi = 0;
                    if (val_len > 0 && val_start[0] == '-') {
                        negative = true;
                        vi = 1;
                    }
                    for (; vi < val_len; ++vi) {
                        v = v * 10 + (val_start[vi] - '0');
                    }
                    if (negative) v = -v;
                    py_val = nb::int_(v);
                } else if (spec.conv_type == CONV_FLOAT) {
                    // Use strtod for float parsing
                    char buf[64];
                    size_t copy_len = val_len < 63 ? val_len : 63;
                    std::memcpy(buf, val_start, copy_len);
                    buf[copy_len] = '\0';
                    double v = std::strtod(buf, nullptr);
                    py_val = nb::float_(v);
                } else { // CONV_ORD
                    if (val_len == 1) {
                        py_val = nb::int_(static_cast<long>(static_cast<uint8_t>(val_start[0])));
                    } else {
                        // Fallback: parse as int
                        long long v = 0;
                        for (size_t vi = 0; vi < val_len; ++vi) {
                            v = v * 10 + (val_start[vi] - '0');
                        }
                        py_val = nb::int_(v);
                    }
                }

                // Replace missing value in result tuple
                // We need to decref the old value and set the new one
                PyObject* old = PyTuple_GET_ITEM(result.ptr(), spec.label_idx);
                Py_INCREF(py_val.ptr());
                PyTuple_SET_ITEM(result.ptr(), spec.label_idx, py_val.ptr());
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
          nb::arg("label_col_sizes"),
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
          "Decode columnar block streams into records (no labels).");
    
    m.def("decode_block", &decode_block_columnar,
          nb::arg("flags_data"),
          nb::arg("lengths_data"),
          nb::arg("seqs_data"),
          nb::arg("len_bytes"),
          nb::arg("count"),
          "Decode block streams into records (no labels).");

    m.def("decode_block_labeled", &decode_block_labeled,
          nb::arg("flags_data"),
          nb::arg("lengths_data"),
          nb::arg("seqs_data"),
          nb::arg("len_bytes"),
          nb::arg("count"),
          nb::arg("label_col_data"),
          nb::arg("label_col_sizes"),
          nb::arg("label_dtype_codes"),
          "Decode columnar block streams with label columns.\n"
          "Returns list of (seq, is_paired, is_read1, is_read2, is_rc, labels) tuples.");

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
