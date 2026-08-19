// Byte-wise SIMD vs 2-bit SWAR for the overlap scan.
//
// A reviewer argued that unaligned 8-bit vector loads plus _mm256_cmpeq_epi8 would beat
// 2-bit packing, because a shift becomes a pointer offset and no cross-word bit
// realignment is needed. That is a real advantage and it deserves a measurement, not an
// argument -- byte comparison would also preserve N/IUPAC semantics exactly, removing
// the packed path's purity dispatch.
//
// Every variant runs the COMPLETE pruned scan (plateau, both flanks, ceiling break,
// per-shift bail) and is checked pair by pair against the shipped kernel's own output.
//
//   c++ -O3 -std=c++17 -o bench_simd bench_simd.cpp && ./bench_simd pairs.bin
//   x86:  add -mavx2 to enable the 32-byte variants (baseline SSE2 needs no flag)

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <chrono>
#include <vector>
#include <algorithm>

#if defined(__ARM_NEON) || defined(__aarch64__)
#  include <arm_neon.h>
#  define HAVE_V16 1
#  define ISA_NAME "NEON"
#elif defined(__SSE2__) || defined(_M_X64)
#  include <immintrin.h>
#  define HAVE_V16 1
#  ifdef __AVX2__
#    define HAVE_V32 1
#    define ISA_NAME "AVX2"
#  else
#    define ISA_NAME "SSE2"
#  endif
#else
#  define ISA_NAME "scalar-only"
#endif

static const int64_t SCALE = 1 << 24;
static int64_t M_W, MM_W, STEP, FLOOR_Q;

// --------------------------------------------------------------- byte compare kernels
#ifdef HAVE_V16
/// mismatching bytes in a 16-byte window
static inline int neq16(const uint8_t* a, const uint8_t* b) {
#if defined(__ARM_NEON) || defined(__aarch64__)
    uint8x16_t eq = vceqq_u8(vld1q_u8(a), vld1q_u8(b));
    return 16 - (int)vaddvq_u8(vandq_u8(eq, vdupq_n_u8(1)));
#else
    __m128i va = _mm_loadu_si128((const __m128i*)a);
    __m128i vb = _mm_loadu_si128((const __m128i*)b);
    return 16 - __builtin_popcount((unsigned)_mm_movemask_epi8(_mm_cmpeq_epi8(va, vb)));
#endif
}
#endif
#ifdef HAVE_V32
static inline int neq32(const uint8_t* a, const uint8_t* b) {
    __m256i va = _mm256_loadu_si256((const __m256i*)a);
    __m256i vb = _mm256_loadu_si256((const __m256i*)b);
    return 32 - __builtin_popcount((unsigned)_mm256_movemask_epi8(_mm256_cmpeq_epi8(va, vb)));
}
#endif

/// The 0.5.2 shipped shape: per-lane equality accumulated in a vector register and
/// reduced ONCE per group (mirrors merge_core.hpp neq16x<V>).  On NEON that trades
/// V x vaddvq_u8 for (V-1) x vaddq_u8 + 1 x vaddvq_u8; on x86 it is what lets psadbw
/// replace pmovmskb+popcount.  Benchmarked here against the per-vector `neq16` above,
/// which is the 0.5.1 "before".
#ifdef HAVE_V16
template <int V>
static inline int neq16x(const uint8_t* a, const uint8_t* b) {
#if defined(__ARM_NEON) || defined(__aarch64__)
    const uint8x16_t one = vdupq_n_u8(1);
    uint8x16_t acc = vdupq_n_u8(0);
    for (int j = 0; j < V; ++j) {
        const uint8x16_t eq = vceqq_u8(vld1q_u8(a + j * 16), vld1q_u8(b + j * 16));
        acc = vaddq_u8(acc, vandq_u8(eq, one));
    }
    return 16 * V - (int)vaddvq_u8(acc);
#else
    __m128i acc = _mm_setzero_si128();
    for (int j = 0; j < V; ++j) {
        __m128i va = _mm_loadu_si128((const __m128i*)(a + j * 16));
        __m128i vb = _mm_loadu_si128((const __m128i*)(b + j * 16));
        acc = _mm_sub_epi8(acc, _mm_cmpeq_epi8(va, vb));
    }
    const __m128i sad = _mm_sad_epu8(acc, _mm_setzero_si128());
    return 16 * V - (_mm_cvtsi128_si32(sad) + (int)_mm_extract_epi16(sad, 4));
#endif
}

/// Shipped 0.5.2 shift_score: grouped reduction, bail between groups, then a 16-byte
/// step loop (STEP16=1) or straight to the byte tail (STEP16=0), then bytes.
template <int BAIL, int STEP16>
static inline int64_t sc_grouped(const uint8_t* s1, const uint8_t* s2rc,
                                 int s, int n, int64_t best, int* out_d) {
    const int64_t ceiling = (int64_t)n * M_W;
    if (ceiling <= best) return INT64_MIN;
    const int64_t dmax = (ceiling - best - 1) / STEP;
    const uint8_t* a = s1   + (s > 0 ?  s : 0);
    const uint8_t* b = s2rc + (s < 0 ? -s : 0);
    int64_t d = 0;
    int k = 0;
    constexpr int GROUP = 16 * BAIL;
    for (; k + GROUP <= n; k += GROUP) {
        d += neq16x<BAIL>(a + k, b + k);
        if (d > dmax) return INT64_MIN;
    }
    if (STEP16) {
        for (; k + 16 <= n; k += 16) {
            d += neq16x<1>(a + k, b + k);
            if (d > dmax) return INT64_MIN;
        }
    }
    for (; k < n; ++k) d += (a[k] != b[k]);
    if (d > dmax) return INT64_MIN;
    *out_d = (int)d;
    return ceiling - d * STEP;
}
#endif

/// VW = bytes per vector op, BAIL = vectors between bail checks
template <int VW, int BAIL>
static inline int64_t sc_simd(const uint8_t* s1, const uint8_t* s2rc,
                              int s, int n, int64_t best, int* out_d) {
    const int64_t ceiling = (int64_t)n * M_W;
    if (ceiling <= best) return INT64_MIN;
    const int64_t dmax = (ceiling - best - 1) / STEP;
    const uint8_t* a = s1   + (s > 0 ?  s : 0);   // a shift is just a pointer offset
    const uint8_t* b = s2rc + (s < 0 ? -s : 0);
    int64_t d = 0;
    int k = 0;
    const int step = VW * BAIL;
    for (; k + step <= n; k += step) {
        int acc = 0;
        for (int j = 0; j < BAIL; ++j) {
#ifdef HAVE_V32
            if (VW == 32) { acc += neq32(a + k + j*32, b + k + j*32); continue; }
#endif
#ifdef HAVE_V16
            if (VW == 16) { acc += neq16(a + k + j*16, b + k + j*16); }
#endif
        }
        d += acc;
        if (d > dmax) return INT64_MIN;
    }
    for (; k < n; ++k) d += (a[k] != b[k]);
    if (d > dmax) return INT64_MIN;
    *out_d = (int)d;
    return ceiling - d * STEP;
}

/// the two flank shifts at one n, interleaved: their bail chains are serial
/// dependencies, so running them together hides latency (worth 1.08x on the packed
/// kernel).
template <int VW, int BAIL>
static inline void sc_simd2(const uint8_t* s1, const uint8_t* s2rc,
                            int sa, int sb, int n, int64_t best,
                            int64_t* va, int* da, int64_t* vb, int* db) {
    const int64_t ceiling = (int64_t)n * M_W;
    *va = *vb = INT64_MIN;
    if (ceiling <= best) return;
    const int64_t dmax = (ceiling - best - 1) / STEP;
    const uint8_t* a1 = s1 + (sa>0?sa:0); const uint8_t* a2 = s2rc + (sa<0?-sa:0);
    const uint8_t* b1 = s1 + (sb>0?sb:0); const uint8_t* b2 = s2rc + (sb<0?-sb:0);
    int64_t x = 0, y = 0; bool lx = false, ly = false;
    int k = 0; const int step = VW * BAIL;
    for (; k + step <= n; k += step) {
        if (!lx) { int t=0; for (int j=0;j<BAIL;++j) t += neq16(a1+k+j*VW, a2+k+j*VW);
                   x += t; lx = x > dmax; }
        if (!ly) { int t=0; for (int j=0;j<BAIL;++j) t += neq16(b1+k+j*VW, b2+k+j*VW);
                   y += t; ly = y > dmax; }
        if (lx && ly) return;
    }
    for (int i = k; i < n; ++i) {
        if (!lx) x += (a1[i] != a2[i]);
        if (!ly) y += (b1[i] != b2[i]);
    }
    if (!lx && x <= dmax) { *va = ceiling - x*STEP; *da = (int)x; }
    if (!ly && y <= dmax) { *vb = ceiling - y*STEP; *db = (int)y; }
}

// --------------------------------------------------------------- scalar + packed refs
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

static inline uint64_t bytes_eq(uint64_t x, uint64_t m) {
    uint64_t y = x ^ m;
    return ~(((y & 0x7F7F7F7F7F7F7F7FULL) + 0x7F7F7F7F7F7F7F7FULL) | y | 0x7F7F7F7F7F7F7F7FULL);
}
static inline bool pack_swar(const uint8_t* s, int n, uint64_t* out) {
    std::memset(out, 0, (size_t)((n + 31) / 32 + 1) * 8);
    const uint64_t A=0x4141414141414141ULL,C=0x4343434343434343ULL;
    const uint64_t G=0x4747474747474747ULL,T=0x5454545454545454ULL;
    bool pure = true; int i = 0;
    for (; i + 8 <= n; i += 8) {
        uint64_t x; std::memcpy(&x, s + i, 8);
        uint64_t u = x & 0xDFDFDFDFDFDFDFDFULL;
        if ((bytes_eq(u,A)|bytes_eq(u,C)|bytes_eq(u,G)|bytes_eq(u,T)) != 0x8080808080808080ULL)
            pure = false;
        uint64_t p = (x >> 1) & 0x0303030303030303ULL;
        p = (p | (p >> 6))  & 0x000F000F000F000FULL;
        p = (p | (p >> 12)) & 0x000000FF000000FFULL;
        p = (p | (p >> 24)) & 0x000000000000FFFFULL;
        out[i >> 5] |= p << (2 * (i & 31));
    }
    for (; i < n; ++i) {
        uint8_t c = s[i], u = c & 0xDF;
        if (u!='A'&&u!='C'&&u!='G'&&u!='T') pure = false;
        out[i >> 5] |= (uint64_t)((c >> 1) & 3) << (2 * (i & 31));
    }
    return pure;
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
static inline void sc_packed2(const uint64_t* p1, const uint64_t* p2,
                              int sa, int sb, int n, int64_t best,
                              int64_t* va, int* da, int64_t* vb, int* db) {
    const int64_t ceiling = (int64_t)n * M_W;
    *va = *vb = INT64_MIN;
    if (ceiling <= best) return;
    const int64_t dmax = (ceiling - best - 1) / STEP;
    const int a1 = sa>0?sa:0, a2 = sa<0?-sa:0, b1 = sb>0?sb:0, b2 = sb<0?-sb:0;
    int64_t x=0,y=0; bool lx=false, ly=false;
    const int nfull = n >> 5;
    for (int k = 0; k < nfull; ++k) {
        if (!lx) { x += mm32(getw(p1,a1,k), getw(p2,a2,k)); lx = x > dmax; }
        if (!ly) { y += mm32(getw(p1,b1,k), getw(p2,b2,k)); ly = y > dmax; }
        if (lx && ly) return;
    }
    const int rem = n & 31;
    if (rem) {
        const uint64_t mask = (1ULL << (2*rem)) - 1;
        if (!lx) x += mm32(getw(p1,a1,nfull)&mask, getw(p2,a2,nfull)&mask);
        if (!ly) y += mm32(getw(p1,b1,nfull)&mask, getw(p2,b2,nfull)&mask);
    }
    if (!lx && x <= dmax) { *va = ceiling - x*STEP; *da = (int)x; }
    if (!ly && y <= dmax) { *vb = ceiling - y*STEP; *db = (int)y; }
}
static inline int64_t sc_packed(const uint64_t* p1, const uint64_t* p2,
                                int s, int n, int64_t best, int* out_d) {
    const int64_t ceiling = (int64_t)n * M_W;
    if (ceiling <= best) return INT64_MIN;
    const int64_t dmax = (ceiling - best - 1) / STEP;
    const int i1 = s>0?s:0, i2 = s<0?-s:0;
    int64_t d = 0; int k = 0; const int nfull = n >> 5;
    while (k < nfull) {
        const int stop = std::min(k + 2, nfull);
        int acc = 0;
        for (; k < stop; ++k) acc += mm32(getw(p1,i1,k), getw(p2,i2,k));
        d += acc;
        if (d > dmax) return INT64_MIN;
    }
    const int rem = n & 31;
    if (rem) {
        const uint64_t mask = (1ULL << (2*rem)) - 1;
        d += mm32(getw(p1,i1,nfull)&mask, getw(p2,i2,nfull)&mask);
    }
    if (d > dmax) return INT64_MIN;
    *out_d = (int)d; return ceiling - d * STEP;
}

// --------------------------------------------------------------- the scan
struct Pair { const uint8_t* s1; const uint8_t* s2rc; int len1, len2; int es, en, ed; };
struct Res { int s, n, d; int64_t score; };

// K: 0 scalar | 1 packed+fused | 2 simd VW=16 | 3 simd VW=32
template <int K, int VW, int BAIL>
static Res scan(const Pair& P, uint64_t* b1, uint64_t* b2) {
    const int len1 = P.len1, len2 = P.len2;
    const int nmax = len1 < len2 ? len1 : len2;
    if (nmax <= 0) return {0,0,0,0};
    bool packed = false;
    if (K == 1) {
        packed = pack_swar(P.s1, len1, b1) & pack_swar(P.s2rc, len2, b2);
        if (!packed) return scan<0,VW,BAIL>(P, b1, b2);
    }
    int64_t best = FLOOR_Q - 1;
    int bs = 0, bn = 0, bd = 0;
    auto one = [&](int s, int n, int* dd) -> int64_t {
        if (K == 0) return sc_scalar(P.s1, P.s2rc, s, n, best, dd);
        if (K == 1) return sc_packed(b1, b2, s, n, best, dd);
#ifdef HAVE_V16
        if (K == 5) return sc_grouped<BAIL,0>(P.s1, P.s2rc, s, n, best, dd);
        if (K == 6) return sc_grouped<BAIL,1>(P.s1, P.s2rc, s, n, best, dd);
#endif
        return sc_simd<VW,BAIL>(P.s1, P.s2rc, s, n, best, dd);
    };
    const int plo = (len1>=len2)?0:len1-len2, phi = (len1>=len2)?len1-len2:0;
    for (int s = plo; s <= phi; ++s) {
        int d = 0; int64_t v = one(s, nmax, &d);
        if (v > best) { best=v; bs=s; bn=nmax; bd=d; }
    }
    for (int n = nmax - 1; n > 0; --n) {
        if ((int64_t)n * M_W <= best) break;
        if (K == 4) {
            int64_t va, vb; int da=0, db=0;
            sc_simd2<VW,BAIL>(P.s1, P.s2rc, n-len2, len1-n, n, best, &va, &da, &vb, &db);
            if (va > best) { best=va; bs=n-len2; bn=n; bd=da; }
            if (vb > best) { best=vb; bs=len1-n; bn=n; bd=db; }
        } else if (K == 1) {
            int64_t va, vb; int da=0, db=0;
            sc_packed2(b1, b2, n-len2, len1-n, n, best, &va, &da, &vb, &db);
            if (va > best) { best=va; bs=n-len2; bn=n; bd=da; }
            if (vb > best) { best=vb; bs=len1-n; bn=n; bd=db; }
        } else {
            int d = 0; int64_t v = one(n-len2, n, &d);
            if (v > best) { best=v; bs=n-len2; bn=n; bd=d; }
            d = 0; v = one(len1-n, n, &d);
            if (v > best) { best=v; bs=len1-n; bn=n; bd=d; }
        }
    }
    if (bn == 0) return {0,0,0,0};
    return {bs, bn, bd, best};
}

static std::vector<Pair> pairs;
static std::vector<std::vector<uint8_t>> store;

template <int K, int VW, int BAIL>
static void run(const char* name, double* base) {
    const uint32_t N = pairs.size();
    std::vector<uint64_t> b1(64), b2(64);
    int bad = 0;
    for (uint32_t i = 0; i < N; ++i) {
        Res r = scan<K,VW,BAIL>(pairs[i], b1.data(), b2.data());
        if (r.s!=pairs[i].es || r.n!=pairs[i].en || r.d!=pairs[i].ed) ++bad;
    }
    double bt = 1e30; volatile int64_t sink = 0;
    for (int rep = 0; rep < 7; ++rep) {
        auto t0 = std::chrono::steady_clock::now();
        for (uint32_t i = 0; i < N; ++i) sink += scan<K,VW,BAIL>(pairs[i], b1.data(), b2.data()).score;
        auto t1 = std::chrono::steady_clock::now();
        bt = std::min(bt, std::chrono::duration<double>(t1-t0).count());
    }
    (void)sink;
    double us = bt / N * 1e6;
    if (*base == 0) *base = us;
    std::printf("%-54s %7.4f us/pair  %5.2fx  [%d mismatch]\n", name, us, *base/us, bad);
}

int main(int argc, char** argv) {
    M_W  = (int64_t)llround(std::log2(0.99/0.25) * SCALE);
    MM_W = (int64_t)llround(std::log2(0.75/0.01) * SCALE);
    STEP = M_W + MM_W; FLOOR_Q = 8*SCALE;

    FILE* f = std::fopen(argv[1], "rb");
    uint32_t N; if (std::fread(&N,4,1,f)!=1) return 1;
    pairs.resize(N); store.resize(N);
    for (uint32_t i = 0; i < N; ++i) {
        uint16_t l1,l2,ol,df; int32_t s;
        std::fread(&l1,2,1,f); std::fread(&l2,2,1,f); std::fread(&s,4,1,f);
        std::fread(&ol,2,1,f); std::fread(&df,2,1,f);
        store[i].resize((size_t)l1+l2+64);
        std::fread(store[i].data(),1,(size_t)l1+l2,f);
        pairs[i] = {store[i].data(), store[i].data()+l1, l1, l2, s, ol, df};
    }
    std::fclose(f);
    std::printf("# %u real pairs | vector ISA: %s | packing inside the timed region\n\n",
                N, ISA_NAME);

    double base = 0;
    run<0,0,0>("A. scalar block-8 (shipped algorithm)", &base);
    double d = base;
#ifdef HAVE_V16
    run<2,16,1>("G. byte SIMD 16B, bail every 16 bases", &d);
    run<2,16,2>("H. byte SIMD 16B, bail every 32 bases", &d);
    run<2,16,4>("I. byte SIMD 16B, bail every 64 bases", &d);
#endif
#ifdef HAVE_V32
    run<3,32,1>("J. byte SIMD 32B, bail every 32 bases", &d);
    run<3,32,2>("K. byte SIMD 32B, bail every 64 bases", &d);
#endif
#ifdef HAVE_V16
    run<4,16,2>("L. byte SIMD 16B, bail 32, FUSED flank pair", &d);
    run<4,16,4>("M. byte SIMD 16B, bail 64, FUSED flank pair", &d);
    std::printf("\n# 0.5.2 grouped reduction (neq16x<V>), no 16-byte step loop\n");
    run<5,16,1>("N1. grouped, bail every 16 bases", &d);
    run<5,16,2>("N2. grouped, bail every 32 bases", &d);
    run<5,16,3>("N3. grouped, bail every 48 bases", &d);
    run<5,16,4>("N4. grouped, bail every 64 bases", &d);
    std::printf("\n# 0.5.2 grouped reduction + 16-byte step loop  (== shipped shift_score)\n");
    run<6,16,1>("S1. grouped+step16, bail 16", &d);
    run<6,16,2>("S2. grouped+step16, bail 32   <- SHIPPED on aarch64", &d);
    run<6,16,3>("S3. grouped+step16, bail 48   <- SHIPPED on x86", &d);
    run<6,16,4>("S4. grouped+step16, bail 64", &d);
    run<6,16,5>("S5. grouped+step16, bail 80", &d);
    run<6,16,6>("S6. grouped+step16, bail 96", &d);
    run<6,16,8>("S8. grouped+step16, bail 128", &d);
#endif
    run<1,0,0>("F. 2-bit packed + popcount, fused flanks", &d);
    return 0;
}
