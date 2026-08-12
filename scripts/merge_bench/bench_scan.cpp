// v2: packing is INSIDE the timed loop (v1 hoisted it, which flattered the packed
// variant), plus three further ideas -- a SWAR packer, a wider bail granularity, and
// fusing the two flank shifts to expose instruction-level parallelism.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <chrono>
#include <vector>
#include <algorithm>

static const int64_t SCALE = 1 << 24;
static int64_t M_W, MM_W, STEP, FLOOR_Q, T_MERGE_Q;

struct Pair { const uint8_t* s1; const uint8_t* s2rc; int len1, len2; int exp_s, exp_n, exp_d; };
struct Res { int s, n, d; int64_t score; };

// --- packers ---------------------------------------------------------------
static uint8_t CODE[256];
static void init_code() {
    std::memset(CODE, 0x80, sizeof CODE);
    CODE['A']=CODE['a']=0; CODE['C']=CODE['c']=1; CODE['G']=CODE['g']=2; CODE['T']=CODE['t']=3;
}

/// table packer: one lookup per base, ORs a poison bit for non-ACGT
static inline bool pack_table(const uint8_t* s, int n, uint64_t* out) {
    const int nw = (n + 31) / 32 + 1;
    std::memset(out, 0, (size_t)nw * 8);
    uint8_t acc = 0;
    for (int i = 0; i < n; ++i) {
        uint8_t c = CODE[s[i]];
        acc |= c;
        out[i >> 5] |= (uint64_t)(c & 3) << (2 * (i & 31));
    }
    return (acc & 0x80) == 0;
}

/// SWAR packer. For equality-only comparison the code assignment just has to be a
/// bijection on ACGT, and (c>>1)&3 already is one (A0 C1 G3 T2). Purity is a separate
/// SWAR test: uppercase, then check the byte is one of the four.
static inline uint64_t bytes_eq(uint64_t x, uint64_t m) {   // 0x80 per byte equal to m
    uint64_t y = x ^ m;
    return ~(((y & 0x7F7F7F7F7F7F7F7FULL) + 0x7F7F7F7F7F7F7F7FULL) | y | 0x7F7F7F7F7F7F7F7FULL);
}
static inline bool pack_swar(const uint8_t* s, int n, uint64_t* out) {
    const int nw = (n + 31) / 32 + 1;
    std::memset(out, 0, (size_t)nw * 8);
    const uint64_t A = 0x4141414141414141ULL, C = 0x4343434343434343ULL;
    const uint64_t G = 0x4747474747474747ULL, T = 0x5454545454545454ULL;
    bool pure = true;
    int i = 0;
    for (; i + 8 <= n; i += 8) {
        uint64_t x; std::memcpy(&x, s + i, 8);
        uint64_t u = x & 0xDFDFDFDFDFDFDFDFULL;                       // uppercase
        uint64_t ok = bytes_eq(u,A) | bytes_eq(u,C) | bytes_eq(u,G) | bytes_eq(u,T);
        if (ok != 0x8080808080808080ULL) pure = false;
        uint64_t p = (x >> 1) & 0x0303030303030303ULL;                // 2 bits per byte
        p = (p | (p >> 6))  & 0x000F000F000F000FULL;                  // gather to nibbles
        p = (p | (p >> 12)) & 0x000000FF000000FFULL;
        p = (p | (p >> 24)) & 0x000000000000FFFFULL;
        out[i >> 5] |= p << (2 * (i & 31));
    }
    for (; i < n; ++i) {
        uint8_t c = s[i], u = c & 0xDF;
        if (u!='A' && u!='C' && u!='G' && u!='T') pure = false;
        out[i >> 5] |= (uint64_t)((c >> 1) & 3) << (2 * (i & 31));
    }
    return pure;
}

// --- scan pieces -----------------------------------------------------------
static inline int64_t sc_scalar(const uint8_t* s1, const uint8_t* s2rc,
                                int s, int n, int64_t best, int* out_d) {
    const int64_t ceiling = (int64_t)n * M_W;
    if (ceiling <= best) return INT64_MIN;
    const int64_t dmax = (ceiling - best - 1) / STEP;
    const int i1 = s > 0 ? s : 0, i2 = s < 0 ? -s : 0;
    int64_t d = 0; int k = 0; const int lim = n - 7;
    while (k < lim) {
        d += (s1[i1+k]!=s2rc[i2+k]) + (s1[i1+k+1]!=s2rc[i2+k+1])
           + (s1[i1+k+2]!=s2rc[i2+k+2]) + (s1[i1+k+3]!=s2rc[i2+k+3])
           + (s1[i1+k+4]!=s2rc[i2+k+4]) + (s1[i1+k+5]!=s2rc[i2+k+5])
           + (s1[i1+k+6]!=s2rc[i2+k+6]) + (s1[i1+k+7]!=s2rc[i2+k+7]);
        if (d > dmax) return INT64_MIN;
        k += 8;
    }
    while (k < n) { d += (s1[i1+k]!=s2rc[i2+k]); ++k; }
    if (d > dmax) return INT64_MIN;
    *out_d = (int)d; return ceiling - d * STEP;
}

static inline uint64_t getw(const uint64_t* P, int base, int k) {
    const int idx = (base >> 5) + k, ph = (base & 31) * 2;
    uint64_t lo = P[idx];
    return ph ? ((lo >> ph) | (P[idx + 1] << (64 - ph))) : lo;
}
static inline int mm32(uint64_t a, uint64_t b) {
    uint64_t x = a ^ b;
    x = (x | (x >> 1)) & 0x5555555555555555ULL;
    return __builtin_popcountll(x);
}

template <int BAIL_WORDS>
static inline int64_t sc_packed(const uint64_t* p1, const uint64_t* p2,
                                int s, int n, int64_t best, int* out_d) {
    const int64_t ceiling = (int64_t)n * M_W;
    if (ceiling <= best) return INT64_MIN;
    const int64_t dmax = (ceiling - best - 1) / STEP;
    const int i1 = s > 0 ? s : 0, i2 = s < 0 ? -s : 0;
    int64_t d = 0; int k = 0; const int nfull = n >> 5;
    while (k < nfull) {
        const int stop = std::min(k + BAIL_WORDS, nfull);
        int acc = 0;
        for (; k < stop; ++k) acc += mm32(getw(p1, i1, k), getw(p2, i2, k));
        d += acc;
        if (d > dmax) return INT64_MIN;
    }
    const int rem = n & 31;
    if (rem) {
        const uint64_t mask = (1ULL << (2 * rem)) - 1;
        d += mm32(getw(p1, i1, nfull) & mask, getw(p2, i2, nfull) & mask);
    }
    if (d > dmax) return INT64_MIN;
    *out_d = (int)d; return ceiling - d * STEP;
}

/// two independent shifts at once: their bail chains are serial dependencies, so
/// interleaving the flank pair should hide latency.
static inline void sc_packed2(const uint64_t* p1, const uint64_t* p2,
                              int sa, int sb, int n, int64_t best,
                              int64_t* va, int* da, int64_t* vb, int* db) {
    const int64_t ceiling = (int64_t)n * M_W;
    *va = *vb = INT64_MIN;
    if (ceiling <= best) return;
    const int64_t dmax = (ceiling - best - 1) / STEP;
    const int a1 = sa > 0 ? sa : 0, a2 = sa < 0 ? -sa : 0;
    const int b1 = sb > 0 ? sb : 0, b2 = sb < 0 ? -sb : 0;
    int64_t x = 0, y = 0; bool lx = false, ly = false;
    const int nfull = n >> 5;
    for (int k = 0; k < nfull; ++k) {
        if (!lx) { x += mm32(getw(p1, a1, k), getw(p2, a2, k)); lx = x > dmax; }
        if (!ly) { y += mm32(getw(p1, b1, k), getw(p2, b2, k)); ly = y > dmax; }
        if (lx && ly) return;
    }
    const int rem = n & 31;
    if (rem) {
        const uint64_t mask = (1ULL << (2 * rem)) - 1;
        if (!lx) x += mm32(getw(p1,a1,nfull)&mask, getw(p2,a2,nfull)&mask);
        if (!ly) y += mm32(getw(p1,b1,nfull)&mask, getw(p2,b2,nfull)&mask);
    }
    if (!lx && x <= dmax) { *va = ceiling - x * STEP; *da = (int)x; }
    if (!ly && y <= dmax) { *vb = ceiling - y * STEP; *db = (int)y; }
}

// --- full scans ------------------------------------------------------------
// KIND 0 scalar | 1 packed(table pack) | 2 packed(swar pack) | 3 packed+fused flanks
template <int KIND, int BAIL>
static Res scan(const Pair& P, uint64_t* b1, uint64_t* b2) {
    const int len1 = P.len1, len2 = P.len2;
    const int nmax = len1 < len2 ? len1 : len2;
    if (nmax <= 0) return {0,0,0,0};
    if (KIND >= 1) {
        bool ok = (KIND == 1) ? (pack_table(P.s1, len1, b1) & pack_table(P.s2rc, len2, b2))
                              : (pack_swar (P.s1, len1, b1) & pack_swar (P.s2rc, len2, b2));
        if (!ok) return scan<0, BAIL>(P, b1, b2);   // non-ACGT -> exact byte semantics
    }
    int64_t best = FLOOR_Q - 1;
    int best_s = 0, best_n = 0, best_d = 0;
    auto one = [&](int s, int n, int* dd) -> int64_t {
        if (KIND == 0) return sc_scalar(P.s1, P.s2rc, s, n, best, dd);
        return sc_packed<BAIL>(b1, b2, s, n, best, dd);
    };
    const int plo = (len1 >= len2) ? 0 : len1 - len2;
    const int phi = (len1 >= len2) ? len1 - len2 : 0;
    for (int s = plo; s <= phi; ++s) {
        int d = 0; int64_t v = one(s, nmax, &d);
        if (v > best) { best = v; best_s = s; best_n = nmax; best_d = d; }
    }
    for (int n = nmax - 1; n > 0; --n) {
        if ((int64_t)n * M_W <= best) break;
        if (KIND == 3) {
            int64_t va, vb; int da = 0, db = 0;
            sc_packed2(b1, b2, n - len2, len1 - n, n, best, &va, &da, &vb, &db);
            if (va > best) { best = va; best_s = n - len2; best_n = n; best_d = da; }
            if (vb > best) { best = vb; best_s = len1 - n; best_n = n; best_d = db; }
        } else {
            int d = 0; int64_t v = one(n - len2, n, &d);
            if (v > best) { best = v; best_s = n - len2; best_n = n; best_d = d; }
            d = 0; v = one(len1 - n, n, &d);
            if (v > best) { best = v; best_s = len1 - n; best_n = n; best_d = d; }
        }
    }
    if (best_n == 0) return {0,0,0,0};
    return {best_s, best_n, best_d, best};
}

// --- driver ----------------------------------------------------------------
static std::vector<Pair> pairs;
static std::vector<std::vector<uint8_t>> store;

template <int KIND, int BAIL>
static void run(const char* name, double* base) {
    const uint32_t N = pairs.size();
    std::vector<uint64_t> b1(64), b2(64);
    int bad = 0, merged = 0;
    for (uint32_t i = 0; i < N; ++i) {
        Res r = scan<KIND, BAIL>(pairs[i], b1.data(), b2.data());
        if (r.s != pairs[i].exp_s || r.n != pairs[i].exp_n || r.d != pairs[i].exp_d) ++bad;
        if (r.score >= T_MERGE_Q) ++merged;
    }
    double best_t = 1e30; volatile int64_t sink = 0;
    for (int rep = 0; rep < 7; ++rep) {
        auto t0 = std::chrono::steady_clock::now();
        for (uint32_t i = 0; i < N; ++i) { Res r = scan<KIND, BAIL>(pairs[i], b1.data(), b2.data()); sink += r.score; }
        auto t1 = std::chrono::steady_clock::now();
        best_t = std::min(best_t, std::chrono::duration<double>(t1 - t0).count());
    }
    (void)sink;
    double us = best_t / N * 1e6;
    if (*base == 0) *base = us;
    std::printf("%-56s %7.4f us/pair  %5.2fx  [%d mismatch, %.1f%% merged]\n",
                name, us, *base / us, bad, 100.0 * merged / N);
}

int main(int argc, char** argv) {
    init_code();
    M_W  = (int64_t)llround(std::log2(0.99/0.25) * SCALE);
    MM_W = (int64_t)llround(std::log2(0.75/0.01) * SCALE);
    STEP = M_W + MM_W; FLOOR_Q = 8*SCALE; T_MERGE_Q = 28*SCALE;

    FILE* f = std::fopen(argv[1], "rb");
    uint32_t N; if (std::fread(&N,4,1,f)!=1) return 1;
    pairs.resize(N); store.resize(N);
    for (uint32_t i = 0; i < N; ++i) {
        uint16_t l1,l2,ol,df; int32_t s;
        std::fread(&l1,2,1,f); std::fread(&l2,2,1,f); std::fread(&s,4,1,f);
        std::fread(&ol,2,1,f); std::fread(&df,2,1,f);
        store[i].resize((size_t)l1+l2+16);
        std::fread(store[i].data(),1,(size_t)l1+l2,f);
        pairs[i] = {store[i].data(), store[i].data()+l1, l1, l2, s, ol, df};
    }
    std::fclose(f);
    std::printf("# %u real pairs; packing is inside the timed region\n\n", N);

    double base = 0;
    run<0,1>("A. scalar block-8 (shipped algorithm, in C++)", &base);
    double dummy = base;
    run<1,1>("B. packed, table packer,  bail every 32 bases", &dummy);
    run<2,1>("C. packed, SWAR packer,   bail every 32 bases", &dummy);
    run<2,2>("D. packed, SWAR packer,   bail every 64 bases", &dummy);
    run<2,4>("E. packed, SWAR packer,   bail every 128 bases", &dummy);
    run<3,1>("F. packed, SWAR packer,   fused flank pair (ILP)", &dummy);
    return 0;
}
