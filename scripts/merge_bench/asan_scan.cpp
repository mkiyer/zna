// Sanitizer driver for the overlap scan.
//
// The vector loop reads 16 bytes at a time and is kept inside the record only by the
// `k + 32 <= n` guard. An off-by-one there is a heap overread that produces plausible
// output -- the worst kind of bug for a corpus tool, and not one code review reliably
// catches. So run the kernel under ASAN/UBSAN against buffers whose bounds the
// allocator actually enforces.
//
// Every read here is placed so that its LAST byte is the last byte of its allocation,
// which is what makes a one-byte overread trap instead of landing in slack.
//
//   c++ -std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer \
//       -I../../src/zna/merge -o asan_scan asan_scan.cpp && ./asan_scan

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <random>
#include <vector>

#include "merge_core.hpp"
#include "fastq_chunk.hpp"

namespace {

// The shipped weights (zna/merge/params.py at err_rate 0.01, SCALE 2**24).
constexpr int64_t MATCH_Q = 33311170;
constexpr int64_t STEP_Q = 137813407;
constexpr int64_t FLOOR_Q = 8 * (1 << 24);

long long checks = 0;

/// Exactly-sized heap buffers, so ASAN's redzones sit immediately after the data.
void run(const std::vector<uint8_t>& a, const std::vector<uint8_t>& b) {
    uint8_t* p1 = static_cast<uint8_t*>(std::malloc(a.size() ? a.size() : 1));
    uint8_t* p2 = static_cast<uint8_t*>(std::malloc(b.size() ? b.size() : 1));
    if (!a.empty()) std::memcpy(p1, a.data(), a.size());
    if (!b.empty()) std::memcpy(p2, b.data(), b.size());
    volatile auto r = zna_merge::scan(p1, static_cast<int>(a.size()),
                                      p2, static_cast<int>(b.size()),
                                      MATCH_Q, STEP_Q, FLOOR_Q);
    (void)r;
    std::free(p1);
    std::free(p2);
    ++checks;
}

std::vector<uint8_t> draw(std::mt19937& rng, size_t n, const char* alpha, size_t na) {
    std::vector<uint8_t> v(n);
    for (size_t i = 0; i < n; ++i) v[i] = static_cast<uint8_t>(alpha[rng() % na]);
    return v;
}

}  // namespace

int main() {
    std::mt19937 rng(20260812);

    // 1. every length combination through and past the vector and bail boundaries
    for (size_t l1 = 0; l1 <= 80; ++l1) {
        for (size_t l2 = 0; l2 <= 80; ++l2) {
            run(draw(rng, l1, "ACGT", 4), draw(rng, l2, "ACGT", 4));
        }
    }

    // 2. identical mates, so no shift ever bails early and every scan runs to the end
    for (size_t n = 0; n <= 200; ++n) {
        auto s = draw(rng, n, "ACGT", 4);
        run(s, s);
    }

    // 3. periodic content: the deep-tie case, where the scan visits the most shifts
    for (const char* per : {"CA", "ACGT", "A", "AATT"}) {
        const size_t pl = std::strlen(per);
        for (size_t n = 1; n <= 200; ++n) {
            std::vector<uint8_t> s(n);
            for (size_t i = 0; i < n; ++i) s[i] = static_cast<uint8_t>(per[i % pl]);
            std::vector<uint8_t> t(n);
            for (size_t i = 0; i < n; ++i) t[i] = static_cast<uint8_t>(per[(i + 1) % pl]);
            run(s, t);
        }
    }

    // 4. arbitrary bytes, not just nucleotides -- the kernel promises byte semantics
    for (int i = 0; i < 4000; ++i) {
        const size_t l1 = rng() % 300, l2 = rng() % 300;
        std::vector<uint8_t> a(l1), b(l2);
        for (auto& c : a) c = static_cast<uint8_t>(rng() & 0xFF);
        for (auto& c : b) c = static_cast<uint8_t>(rng() & 0xFF);
        run(a, b);
    }

    // 5. long reads, where the O(L^2) scan visits the most shifts of all
    for (size_t n : {512u, 1024u, 2048u}) {
        run(draw(rng, n, "ACGT", 4), draw(rng, n, "ACGT", 4));
        auto s = draw(rng, n, "ACGT", 4);
        run(s, s);
    }

    // ---- the chunk adapter: parser, formatter, counters -------------------------
    //
    // The parser is where the audit's raw-blob prototype had four defects, three of
    // them out-of-bounds reads on malformed input. Feed it truncations at EVERY byte
    // offset, out of exactly-sized allocations so a one-byte overread traps.
    const uint8_t table[256 * 256] = {};                 // contents irrelevant here
    const zna_merge::Params params{33311170, 137813407, 28 * (1 << 24), 8 * (1 << 24),
                                   40, table};
    std::string good;
    for (int i = 0; i < 6; ++i) {
        auto s = draw(rng, 60 + (size_t)(rng() % 40), "ACGT", 4);
        good += "@read" + std::to_string(i) + "/1 tag\n";
        good.append(reinterpret_cast<const char*>(s.data()), s.size());
        good += "\n+\n";
        good.append(s.size(), 'I');
        good += "\n";
    }
    long long chunks = 0;
    zna_merge::ChunkScratch sc;
    for (size_t cut = 0; cut <= good.size(); ++cut) {
        for (size_t cut2 = 0; cut2 <= good.size(); cut2 += 7) {
            uint8_t* a1 = static_cast<uint8_t*>(std::malloc(cut ? cut : 1));
            uint8_t* a2 = static_cast<uint8_t*>(std::malloc(cut2 ? cut2 : 1));
            std::memcpy(a1, good.data(), cut);
            std::memcpy(a2, good.data(), cut2);
            std::string blob;
            zna_merge::ChunkStats st;
            size_t p1 = 0, p2 = 0;
            try {
                zna_merge::merge_chunk(a1, cut, p1, a2, cut2, p2, params, true, 0,
                                       sc, blob, st);
            } catch (const zna_merge::InputError&) {
                // malformed input is supposed to raise; the point is that it does not
                // read out of bounds on the way
            }
            std::free(a1);
            std::free(a2);
            ++chunks;
        }
    }

    // arbitrary bytes: no structure at all
    for (int i = 0; i < 3000; ++i) {
        const size_t n = rng() % 400;
        std::vector<uint8_t> v(n);
        for (auto& c : v) c = static_cast<uint8_t>(rng() & 0xFF);
        uint8_t* a1 = static_cast<uint8_t*>(std::malloc(n ? n : 1));
        // `v.data()` is nullptr when n == 0, and memcpy's arguments are declared
        // non-null, so UBSAN flags the zero-length copy. Guard it the way `run()`
        // above already does.
        if (n) std::memcpy(a1, v.data(), n);
        std::string blob;
        zna_merge::ChunkStats st;
        size_t p1 = 0, p2 = 0;
        try {
            zna_merge::merge_chunk(a1, n, p1, a1, n, p2, params, false, 0, sc, blob, st);
        } catch (const zna_merge::InputError&) {}
        std::free(a1);
        ++chunks;
    }

    std::printf("asan_scan: %lld scans and %lld chunks clean\n", checks, chunks);
    return 0;
}
