// Full-path prototype: parse -> scan -> consensus -> build -> format, all in C++.
// Writes an uncompressed mixed interleaved FASTQ that must be BYTE-IDENTICAL to
// `zna merge`'s. That equality validates the fixed-point score, the tie-break, the
// consensus table and the record construction all at once.
//
//   ./proto_merge r1.fq.gz r2.fq.gz out.fq [min_read_length]

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <chrono>
#include <string>
#include <vector>
#include <algorithm>
#if defined(__ARM_NEON) || defined(__aarch64__)
#  include <arm_neon.h>
#else
#  include <immintrin.h>
#endif

// --------------------------------------------------------------------------- params
static const int64_t SCALE = 1 << 24;
static int64_t M_W, MM_W, STEP, FLOOR_Q, T_MERGE_Q;
static int MIN_READ_LEN = 40;

// --------------------------------------------------------------------------- tables
static uint8_t DISAGREE_Q[127][127];
static void init_disagree() {
    for (int qw = 33; qw < 127; ++qw)
        for (int ql = 33; ql < 127; ++ql) {
            double pw = std::pow(10.0, -(qw - 33) / 10.0);
            double pl = std::pow(10.0, -(ql - 33) / 10.0);
            double num = pw * (1.0 - pl);
            double den = num + pl * (1.0 - pw);
            double post = (den <= 0.0) ? 0.5 : num / den;
            post = std::min(std::max(post, 1e-10), 0.9999);
            int q = (int)std::lround(-10.0 * std::log10(post));
            int v = q + 33; v = std::min(std::max(v, 33), 126);
            DISAGREE_Q[qw][ql] = (uint8_t)v;
        }
}
static char COMPL[256];
static void init_compl() {
    for (int i = 0; i < 256; ++i) COMPL[i] = (char)i;
    COMPL['A']='T'; COMPL['C']='G'; COMPL['G']='C'; COMPL['T']='A'; COMPL['N']='N';
    COMPL['a']='t'; COMPL['c']='g'; COMPL['g']='c'; COMPL['t']='a'; COMPL['n']='n';
}

// --------------------------------------------------------------------------- packing
static inline uint64_t bytes_eq(uint64_t x, uint64_t m) {
    uint64_t y = x ^ m;
    return ~(((y & 0x7F7F7F7F7F7F7F7FULL) + 0x7F7F7F7F7F7F7F7FULL) | y | 0x7F7F7F7F7F7F7F7FULL);
}
static inline bool pack_swar(const uint8_t* s, int n, uint64_t* out) {
    const int nw = (n + 31) / 32 + 1;
    std::memset(out, 0, (size_t)nw * 8);
    const uint64_t A=0x4141414141414141ULL, C=0x4343434343434343ULL;
    const uint64_t G=0x4747474747474747ULL, T=0x5454545454545454ULL;
    bool pure = true;
    int i = 0;
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

/// pack revcomp(s2) directly, with no intermediate buffer. Under the (c>>1)&3 map
/// (A0 C1 G3 T2) complementing is XOR 2: A0<->T2 and C1<->G3.
static inline bool pack_swar_rc(const uint8_t* s, int n, uint64_t* out) {
    const int nw = (n + 31) / 32 + 1;
    std::memset(out, 0, (size_t)nw * 8);
    bool pure = true;
    for (int i = 0; i < n; ++i) {
        uint8_t c = s[n - 1 - i], u = c & 0xDF;
        if (u!='A'&&u!='C'&&u!='G'&&u!='T') pure = false;
        out[i >> 5] |= (uint64_t)(((c >> 1) & 3) ^ 2) << (2 * (i & 31));
    }
    return pure;
}

// --------------------------------------------------------------------------- scan
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
static inline int64_t sc_packed(const uint64_t* p1, const uint64_t* p2,
                                int s, int n, int64_t best, int* out_d) {
    const int64_t ceiling = (int64_t)n * M_W;
    if (ceiling <= best) return INT64_MIN;
    const int64_t dmax = (ceiling - best - 1) / STEP;
    const int i1 = s > 0 ? s : 0, i2 = s < 0 ? -s : 0;
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
static inline int64_t sc_simd(const uint8_t* s1, const uint8_t* s2rc,
                              int s, int n, int64_t best, int* out_d) {
    const int64_t ceiling = (int64_t)n * M_W;
    if (ceiling <= best) return INT64_MIN;
    const int64_t dmax = (ceiling - best - 1) / STEP;
    const uint8_t* a = s1   + (s > 0 ?  s : 0);
    const uint8_t* b = s2rc + (s < 0 ? -s : 0);
    int64_t d = 0; int k = 0;
    for (; k + 32 <= n; k += 32) {
        d += neq16(a + k, b + k) + neq16(a + k + 16, b + k + 16);
        if (d > dmax) return INT64_MIN;
    }
    for (; k < n; ++k) d += (a[k] != b[k]);
    if (d > dmax) return INT64_MIN;
    *out_d = (int)d; return ceiling - d * STEP;
}

/// exact byte semantics, for reads carrying non-ACGT
static inline int64_t sc_bytes(const uint8_t* s1, const uint8_t* s2rc,
                               int s, int n, int64_t best, int* out_d) {
    const int64_t ceiling = (int64_t)n * M_W;
    if (ceiling <= best) return INT64_MIN;
    const int64_t dmax = (ceiling - best - 1) / STEP;
    const int i1 = s > 0 ? s : 0, i2 = s < 0 ? -s : 0;
    int64_t d = 0;
    for (int k = 0; k < n; ++k) {
        d += (s1[i1+k] != s2rc[i2+k]);
        if ((k & 7) == 7 && d > dmax) return INT64_MIN;
    }
    if (d > dmax) return INT64_MIN;
    *out_d = (int)d; return ceiling - d * STEP;
}

struct Res { int s, n, d; int64_t score; };

static Res scan(const uint8_t* s1, int len1, const uint8_t* s2rc, int len2,
                const uint64_t* p1, const uint64_t* p2, bool packed) {
    const int nmax = len1 < len2 ? len1 : len2;
    if (nmax <= 0) return {0,0,0,0};
    int64_t best = FLOOR_Q - 1;
    int best_s = 0, best_n = 0, best_d = 0;
    auto one = [&](int s, int n, int* dd) -> int64_t {
        (void)packed; (void)p1; (void)p2;
        return sc_simd(s1, s2rc, s, n, best, dd);
    };
    const int plo = (len1 >= len2) ? 0 : len1 - len2;
    const int phi = (len1 >= len2) ? len1 - len2 : 0;
    for (int s = plo; s <= phi; ++s) {
        int d = 0; int64_t v = one(s, nmax, &d);
        if (v > best) { best=v; best_s=s; best_n=nmax; best_d=d; }
    }
    for (int n = nmax - 1; n > 0; --n) {
        if ((int64_t)n * M_W <= best) break;
        int d = 0; int64_t v = one(n - len2, n, &d);
        if (v > best) { best=v; best_s=n-len2; best_n=n; best_d=d; }
        d = 0; v = one(len1 - n, n, &d);
        if (v > best) { best=v; best_s=len1-n; best_n=n; best_d=d; }
    }
    if (best_n == 0) return {0,0,0,0};
    return {best_s, best_n, best_d, best};
}

// --------------------------------------------------------------------------- reader
struct Stream {
    FILE* p; std::string buf; size_t pos = 0; bool eof = false;
    void open(const char* path) {
        std::string cmd = std::string("pigz -dc -p 1 '") + path + "'";
        p = popen(cmd.c_str(), "r");
    }
    bool fill() {
        if (eof) return false;
        if (pos) { buf.erase(0, pos); pos = 0; }
        size_t old = buf.size();
        buf.resize(old + (1 << 20));
        size_t got = fread(&buf[old], 1, 1 << 20, p);
        buf.resize(old + got);
        if (got == 0) { eof = true; return false; }
        return true;
    }
};
struct Rec { const uint8_t *h, *s, *q; int hl, sl, ql; };

static bool next_rec(Stream& st, Rec& r) {
    for (;;) {
        const char* b = st.buf.data() + st.pos;
        size_t avail = st.buf.size() - st.pos;
        if (avail) {
            const char* e1 = (const char*)memchr(b, '\n', avail);
            if (e1) {
                const char* e2 = (const char*)memchr(e1+1, '\n', avail - (e1+1-b));
                if (e2) {
                    const char* e3 = (const char*)memchr(e2+1, '\n', avail - (e2+1-b));
                    if (e3) {
                        const char* e4 = (const char*)memchr(e3+1, '\n', avail - (e3+1-b));
                        if (e4) {
                            r.h = (const uint8_t*)b+1; r.hl = (int)(e1-b)-1;
                            r.s = (const uint8_t*)e1+1; r.sl = (int)(e2-e1)-1;
                            r.q = (const uint8_t*)e3+1; r.ql = (int)(e4-e3)-1;
                            st.pos += (e4 - b) + 1;
                            if (r.sl != r.ql) { fprintf(stderr, "seq/qual length mismatch\n"); exit(1); }
                            return true;
                        }
                    }
                }
            }
        }
        if (!st.fill()) return false;
    }
}

// --------------------------------------------------------------------------- arena
/// One scratch arena per worker. Starts at 1024 bases and doubles when a longer read
/// turns up, so nothing needs to know the read length up front and the per-pair path
/// never allocates. `ensure` is one predictable branch that is essentially never taken.
struct Scratch {
    std::vector<uint8_t> s2rc, q2r, s1b, q1b, seq, qual;
    size_t cap = 0;
    size_t grows = 0;
    inline void ensure(size_t n) {
        if (n <= cap) return;
        size_t c = cap ? cap : 1024;
        while (c < n) c <<= 1;
        s2rc.resize(c); q2r.resize(c); s1b.resize(c); q1b.resize(c);
        seq.resize(2 * c); qual.resize(2 * c);   // a merged record is at most len1+len2
        cap = c; ++grows;
    }
};

// --------------------------------------------------------------------------- helpers
static inline int id_end(const uint8_t* h, int n) {
    int sp = -1, tb = -1;
    for (int i = 0; i < n; ++i) { if (h[i]==' ' && sp<0) sp=i; if (h[i]=='\t' && tb<0) tb=i; }
    if (sp < 0) return tb < 0 ? n : tb;
    if (tb < 0) return sp;
    return sp < tb ? sp : tb;
}

int main(int argc, char** argv) {
    init_disagree(); init_compl();
    M_W  = (int64_t)llround(std::log2(0.99/0.25) * SCALE);
    MM_W = (int64_t)llround(std::log2(0.75/0.01) * SCALE);
    STEP = M_W + MM_W; FLOOR_Q = 8*SCALE; T_MERGE_Q = 28*SCALE;
    if (argc > 4) MIN_READ_LEN = atoi(argv[4]);

    Stream a, b; a.open(argv[1]); b.open(argv[2]);
    FILE* out = fopen(argv[3], "wb");
    std::string blob; blob.reserve(1 << 22);

    std::vector<uint64_t> p1(64), p2(64);
    std::vector<uint8_t> s2rc, q2r, mseq, mqual, s1buf, q1buf;
    long npairs=0, nmerged=0, ntrim=0, nkept=0, nemit=0, ndrop=0, nfrag=0;

    auto t0 = std::chrono::steady_clock::now();
    Rec r1, r2;
    Scratch sc;
    while (next_rec(a, r1)) {
        if (!next_rec(b, r2)) { fprintf(stderr, "R2 exhausted before R1\n"); return 1; }
        ++npairs;
        const int len1 = r1.sl, len2 = r2.sl;
        sc.ensure((size_t)(len1 > len2 ? len1 : len2));

        uint8_t* s2rc = sc.s2rc.data();
        for (int i = 0; i < len2; ++i) s2rc[i] = (uint8_t)COMPL[r2.s[len2-1-i]];

        Res R = scan(r1.s, len1, s2rc, len2, nullptr, nullptr, false);

        // Copy-on-write: 56.5% of pairs have a clean overlap and need no copy at all.
        const uint8_t* S1 = r1.s;
        const uint8_t* Q1 = r1.q;
        if (R.d > 0) {
            uint8_t* q2r = sc.q2r.data();
            for (int i = 0; i < len2; ++i) q2r[i] = r2.q[len2-1-i];
            uint8_t* s1b = sc.s1b.data(); uint8_t* q1b = sc.q1b.data();
            std::memcpy(s1b, r1.s, len1); std::memcpy(q1b, r1.q, len1);
            const int a0 = (R.s >= 0) ? R.s : 0;
            const int b0 = (R.s >= 0) ? 0 : -R.s;
            for (int i = 0; i < R.n; ++i) {
                const int ia = a0 + i, ib = b0 + i;
                if (s1b[ia] != s2rc[ib]) {
                    const uint8_t qa = q1b[ia], qb = q2r[ib];
                    if (qb > qa) { s1b[ia] = s2rc[ib]; q1b[ia] = DISAGREE_Q[qb][qa]; }
                    else         { q1b[ia] = DISAGREE_Q[qa][qb]; }
                }
            }
            S1 = s1b; Q1 = q1b;
        }

        const bool has = (R.n > 0);
        if (has && R.score >= T_MERGE_Q) {
            const int s = R.s, L = s + len2;
            const int take1 = len1 < L ? len1 : L, take2 = L - take1;
            uint8_t* mseq = sc.seq.data(); uint8_t* mqual = sc.qual.data();
            std::memcpy(mseq, S1, take1); std::memcpy(mqual, Q1, take1);
            if (take2) {
                const int off = take1 - s;
                std::memcpy(mseq + take1, s2rc + off, take2);
                for (int i = 0; i < take2; ++i) mqual[take1+i] = r2.q[len2-1-(off+i)];
            }
            ++nmerged;
            if (L >= MIN_READ_LEN) {
                int cut = id_end(r1.h, r1.hl);
                int idlen = cut;
                if (idlen >= 2 && r1.h[idlen-2]=='/' && (r1.h[idlen-1]=='1'||r1.h[idlen-1]=='2')) idlen -= 2;
                char nm[64];
                int nl = snprintf(nm, sizeof nm, " merged_%d_%d", take1, take2);
                blob += '@';
                blob.append((const char*)r1.h, idlen);
                blob.append((const char*)r1.h + cut, r1.hl - cut);
                blob.append(nm, nl);
                blob += '\n';
                blob.append((const char*)mseq, L); blob += "\n+\n";
                blob.append((const char*)mqual, L); blob += '\n';
                ++nemit;
            } else ++ndrop;
        } else if (has && R.s >= 0 && R.score >= FLOOR_Q && len2 - R.n >= MIN_READ_LEN) {
            const int keep2 = len2 - R.n;
            ++ntrim;
            if (len1 >= MIN_READ_LEN && keep2 >= MIN_READ_LEN) {
                blob += '@'; blob.append((const char*)r1.h, r1.hl); blob += '\n';
                blob.append((const char*)S1, len1); blob += "\n+\n";
                blob.append((const char*)Q1, len1); blob += '\n';
                blob += '@'; blob.append((const char*)r2.h, r2.hl); blob += '\n';
                blob.append((const char*)r2.s, keep2); blob += "\n+\n";
                blob.append((const char*)r2.q, keep2); blob += '\n';
                nemit += 2;
            } else { ndrop += 2; ++nfrag; }
        } else {
            ++nkept;
            if (len1 >= MIN_READ_LEN && len2 >= MIN_READ_LEN) {
                blob += '@'; blob.append((const char*)r1.h, r1.hl); blob += '\n';
                blob.append((const char*)S1, len1); blob += "\n+\n";
                blob.append((const char*)Q1, len1); blob += '\n';
                blob += '@'; blob.append((const char*)r2.h, r2.hl); blob += '\n';
                blob.append((const char*)r2.s, len2); blob += "\n+\n";
                blob.append((const char*)r2.q, len2); blob += '\n';
                nemit += 2;
            } else { ndrop += 2; ++nfrag; }
        }
        if (blob.size() > (3u << 20)) { fwrite(blob.data(), 1, blob.size(), out); blob.clear(); }
    }
    fwrite(blob.data(), 1, blob.size(), out);
    fclose(out); pclose(a.p); pclose(b.p);
    auto t1 = std::chrono::steady_clock::now();
    double sec = std::chrono::duration<double>(t1 - t0).count();
    fprintf(stderr, "pairs %ld  merged %ld (%.1f%%)  trim %ld  kept %ld  emitted %ld  dropped %ld\n",
            npairs, nmerged, 100.0*nmerged/npairs, ntrim, nkept, nemit, ndrop);
    fprintf(stderr, "%.3f s   %.4f us/pair   (arena grew %zu times)\n", sec, sec/npairs*1e6, sc.grows);
    return 0;
}
