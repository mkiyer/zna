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
#include <nanobind/ndarray.h>

#include <cstdint>
#include <cstring>  // for memcpy
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


// Record tuple type: (sequence, is_paired, is_read1, is_read2)
using Record = std::tuple<std::string, bool, bool, bool>;

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
        
        results.emplace_back(std::move(seq), is_paired, is_read1, is_read2);
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


NB_MODULE(_accel, m) {
    m.doc() = "ZNA accelerated encode/decode functions";
    
    m.def("encode_sequence", &encode_sequence,
          nb::arg("seq"),
          "Encode DNA sequence to 2-bit packed bytes.\n\n"
          "Each base (A, C, G, T) is encoded as 2 bits.\n"
          "Four bases are packed into each byte.");
    
    m.def("decode_block_records", &decode_block_records,
          nb::arg("block_data"),
          nb::arg("len_bytes"),
          nb::arg("count"),
          "Decode all records from a block at once.\n\n"
          "Processes entire block in C++ to avoid Python overhead.");
    
    m.def("reverse_complement", &reverse_complement,
          nb::arg("seq"),
          "Compute reverse complement of a DNA sequence.\n\n"
          "A<->T, C<->G. Used for strand-specific library normalization.");
}
