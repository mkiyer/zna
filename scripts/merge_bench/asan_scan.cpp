// Sanitizer driver for the overlap scan.
//
// The vector loop reads 16 bytes at a time and is kept inside the record only by the
// `k + 32 <= n` guard. An off-by-one there is a heap overread that produces plausible
// output -- the worst kind of bug for a corpus tool, and not one code review reliably
// catches. So run the kernel under ASAN/UBSAN against buffers whose bounds the
// allocator actually enforces.
//
// Every read here is placed so that its LAST byte is the last byte of its allocation,
// which is what makes a one-byte overread trap instead of landing in slack. The 0.6
// policy tables are allocated the same way, exactly as long as their capacity, so a
// T or dfit lookup one past the end traps too.
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

/// Capacity of the policy tables: the longest read below. T_q[N] needs N = len1 + len2
/// - 1 < 2 * CAP entries, dfit[n] needs n <= CAP.
constexpr size_t CAP = 2048;

long long checks = 0;

/// Exactly-sized tables (see the header comment). The values are plausible, not exact --
/// a 28-bit floor and dfit ~ n/16 -- since only the indexing is under test here.
struct OwnedTables {
    int64_t* t;
    int64_t* d;
    zna_merge::Tables tab;
    OwnedTables() {
        t = static_cast<int64_t*>(std::malloc(2 * CAP * sizeof(int64_t)));
        d = static_cast<int64_t*>(std::malloc((CAP + 1) * sizeof(int64_t)));
        for (size_t i = 0; i < 2 * CAP; ++i) t[i] = 28 * (int64_t(1) << 24);
        for (size_t i = 0; i <= CAP; ++i) d[i] = static_cast<int64_t>(i / 16);
        tab = {t, 2 * CAP, d, CAP + 1};
    }
    ~OwnedTables() { std::free(t); std::free(d); }
};
const OwnedTables& tables() {
    static const OwnedTables o;
    return o;
}

/// Exactly-sized heap buffers, so ASAN's redzones sit immediately after the data.
void run(const std::vector<uint8_t>& a, const std::vector<uint8_t>& b) {
    uint8_t* p1 = static_cast<uint8_t*>(std::malloc(a.size() ? a.size() : 1));
    uint8_t* p2 = static_cast<uint8_t*>(std::malloc(b.size() ? b.size() : 1));
    if (!a.empty()) std::memcpy(p1, a.data(), a.size());
    if (!b.empty()) std::memcpy(p2, b.data(), b.size());
    for (bool at : {false, true}) {
        volatile auto r = zna_merge::scan(p1, static_cast<int>(a.size()),
                                          p2, static_cast<int>(b.size()),
                                          MATCH_Q, STEP_Q, FLOOR_Q, at);
        (void)r;
        if (a.size() <= CAP && b.size() <= CAP) {
            volatile auto dec = zna_merge::decide(p1, static_cast<int>(a.size()),
                                                  p2, static_cast<int>(b.size()),
                                                  MATCH_Q, STEP_Q, tables().tab, at);
            (void)dec;
        }
    }
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
    static uint8_t table[256 * 256];                     // contents irrelevant here
    const zna_merge::Params params{MATCH_Q, STEP_Q, tables().tab, false, 40, table};
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
                                       100000, sc, blob, st);
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
            zna_merge::merge_chunk(a1, n, p1, a1, n, p2, params, false, 0, 100000, sc,
                                   blob, st);
        } catch (const zna_merge::InputError&) {}
        std::free(a1);
        ++chunks;
    }

    // ---- whole pairs with no-calls, every N policy, both contracts ----------------
    //
    // The copy-on-write buffers, the N policy's in-place substitution, and the path a
    // merge verdict takes back to KEPT when trim3 breaks tiling (R1 re-derived from the
    // input) all write into the scratch arena; exactly-sized reads keep them honest.
    long long pairs = 0;
    zna_merge::Scratch ps;
    for (int i = 0; i < 20000; ++i) {
        const size_t fl = 20 + rng() % 300;
        auto frag = draw(rng, fl, "ACGT", 4);
        const size_t l1 = 1 + rng() % 200, l2 = 1 + rng() % 200;
        std::vector<uint8_t> s1(l1), s2(l2);
        for (size_t k = 0; k < l1; ++k) s1[k] = k < fl ? frag[k] : 'A';
        for (size_t k = 0; k < l2; ++k) {
            const uint8_t c = k < fl ? frag[fl - 1 - k] : 'C';
            s2[k] = c == 'A' ? 'T' : c == 'T' ? 'A' : c == 'C' ? 'G' : 'C';
        }
        for (auto* s : {&s1, &s2}) {
            const int nn = static_cast<int>(rng() % 6);
            for (int k = 0; k < nn; ++k) (*s)[rng() % s->size()] = 'N';
            if (rng() % 4 == 0) (*s)[rng() % s->size()] = static_cast<uint8_t>(rng());
        }
        std::vector<uint8_t> q1(l1, 'I'), q2(l2, '5');
        std::string h1 = "p" + std::to_string(i) + "/1", h2 = "p" + std::to_string(i) + "/2";
        auto own = [](const void* src, size_t n) {
            uint8_t* p = static_cast<uint8_t*>(std::malloc(n ? n : 1));
            if (n) std::memcpy(p, src, n);
            return p;
        };
        uint8_t* bs1 = own(s1.data(), l1); uint8_t* bq1 = own(q1.data(), l1);
        uint8_t* bs2 = own(s2.data(), l2); uint8_t* bq2 = own(q2.data(), l2);
        uint8_t* bh1 = own(h1.data(), h1.size()); uint8_t* bh2 = own(h2.data(), h2.size());
        const zna_merge::Read r1{{bh1, (int)h1.size()}, {bs1, (int)l1}, {bq1, (int)l1}};
        const zna_merge::Read r2{{bh2, (int)h2.size()}, {bs2, (int)l2}, {bq2, (int)l2}};
        for (int npol = 0; npol < 3; ++npol) {
            for (bool at : {false, true}) {
                const zna_merge::Params pp{MATCH_Q, STEP_Q, tables().tab, at, 40, table,
                                           npol, 42};
                volatile auto r = zna_merge::process_pair(r1, r2, pp, ps, i, (i & 1) != 0);
                (void)r;
                ++pairs;
            }
        }
        for (uint8_t* p : {bs1, bq1, bs2, bq2, bh1, bh2}) std::free(p);
    }

    std::printf("asan_scan: %lld scans, %lld chunks and %lld pairs clean\n", checks,
                chunks, pairs);
    return 0;
}
