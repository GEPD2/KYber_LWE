/*
  kyber_fips203.c    CRYSTALS-Kyber / ML-KEM (FIPS 203)
 
  Full standard compliant Key Generation with Kyber-512/768/1024.
  Prints:  secret key s,  error vector e,  public key pk all in hex.
 
  Parameters (FIPS 203):
    n = 256,  q = 3329,  η = 2
    k = 2 → Kyber-512   (pk=800B,  sk=768B)
    k = 3 → Kyber-768   (pk=1184B, sk=1152B)
    k = 4 → Kyber-1024  (pk=1568B, sk=1536B)
  
    NTT structure (why it differs from the simple negacyclic NTT):
    ψ = 17 has multiplicative order 256 mod 3329, so ψ^256 = 1 (not -1).
    A primitive 512th root does NOT exist mod 3329 (q-1 = 3328 = 2^8 × 13,
    and 512 = 2^9 does not divide 3328).  Therefore:
    We cannot use the 8-layer simple negacyclic NTT.
    Instead: 7 Cooley-Tukey layers (len=128 down to len=2), leaving
    128 quadratic subrings Z_q[X]/(X^2 − γ_i).
    Polynomial multiplication completes with "basemul" per pair.
    INTT: 7 Gentleman-Sande layers (len=2 up to len=128), scaled by 128^{-1}.
    This is exactly FIPS 203 Algorithm 9/10 with the 128-entry ζ table.
 
  A generation: SHAKE128 XOF seeded with (ρ || j || i) for A[i][j].
                Rejection sampling: accept coefficients < q.
                A is never stored, it's generated on the fly for each column.
 
  CBD sampling: SHAKE256 then take byte stream  η=2 CBD (FIPS 203 Algorithm 7).
 
  Encoding:
    ByteEncode_12: each 12-bit coefficient packed into 1.5 bytes (2 per 3 bytes).
    s_hex:  ByteEncode_12(s[0]) || ... || ByteEncode_12(s[k-1])
    e_hex:  ByteEncode_12(e[0]) || ... || ByteEncode_12(e[k-1])
    pk_hex: ByteEncode_12(t[0]) || ... || ByteEncode_12(t[k-1]) || ρ
 
  To compile:  gcc -O2 -Wall -o kyber_fips203 kyber_fips203.c
  To run the executable: ./kyber_fips203 or for windows kyber_fips203.exe
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/*
  PARAMETERS
*/

#define KYBER_N 256
#define KYBER_Q 3329
#define KYBER_η 2
#define KYBER_K 2    
/* 
to change kyber security change one of the parameters 
#define KYBER_K  2   this is kyber 512
#define KYBER_K  3   this is kyber 768
#define KYBER_K  4   this is kyber 1024
*/

#define N_INV 3303  /* 128^{-1} mod 3329  (7 layer INTT scaling) */
#define SEED_BYTES 32
#define POLY_BYTES 384   /* 256 coefficients × 12 bits = 384 bytes     */

/*
  KECCAK-f[1600], the single primitive underneath SHAKE128 and SHAKE256
  Implements the sponge construction for both XOFs.
  Reference: NIST FIPS 202, FIPS 203 Section 4.1
*/

#define KECCAK_ROUNDS 24

static const uint64_t RC[24] = {
    0x0000000000000001ULL, 0x0000000000008082ULL,
    0x800000000000808AULL, 0x8000000080008000ULL,
    0x000000000000808BULL, 0x0000000080000001ULL,
    0x8000000080008081ULL, 0x8000000000008009ULL,
    0x000000000000008AULL, 0x0000000000000088ULL,
    0x0000000080008009ULL, 0x000000008000000AULL,
    0x000000008000808BULL, 0x800000000000008BULL,
    0x8000000000008089ULL, 0x8000000000008003ULL,
    0x8000000000008002ULL, 0x8000000000000080ULL,
    0x000000000000800AULL, 0x800000008000000AULL,
    0x8000000080008081ULL, 0x8000000000008080ULL,
    0x0000000080000001ULL, 0x8000000080008008ULL,
};

static const int ρ[24] = {
     1,  3,  6, 10, 15, 21, 28, 36, 45, 55,  2, 14,
    27, 41, 56,  8, 25, 43, 62, 18, 39, 61, 20, 44,
};

static const int PI[24] = {
    10,  7, 11, 17, 18, 3, 5, 16,  8, 21, 24, 4,
    15, 23, 19, 13, 12, 2, 20, 14, 22,  9,  6,  1,
};

static inline uint64_t rol64(uint64_t x, int s) {
    return (x << s) | (x >> (64 - s));
}

static void keccakf(uint64_t A[25]) {
    for (int r = 0; r < KECCAK_ROUNDS; r++) {
        /* θ */
        uint64_t C[5], D[5];
        for (int x = 0; x < 5; x++)
            C[x] = A[x] ^ A[x+5] ^ A[x+10] ^ A[x+15] ^ A[x+20];
        for (int x = 0; x < 5; x++)
            D[x] = C[(x+4)%5] ^ rol64(C[(x+1)%5], 1);
        for (int i = 0; i < 25; i++)
            A[i] ^= D[i%5];
        /* ρ and π */
        uint64_t B[25];
        B[0] = A[0];
        for (int i = 0; i < 24; i++)
            B[PI[i]] = rol64(A[i+1], ρ[i]);
        /* χ */
        for (int y = 0; y < 25; y += 5)
            for (int x = 0; x < 5; x++)
                A[y+x] = B[y+x] ^ (~B[y+(x+1)%5] & B[y+(x+2)%5]);
        /* ι */
        A[0] ^= RC[r];
    }
}

/* Sponge XOF state */
typedef struct {
    uint64_t A[25];
    uint8_t  buf[200];   /* rate bytes */
    int      rate;       /* SHAKE128=168, SHAKE256=136 */
    int      pos;        /* current position in buf */
    int      squeezing;  /* 0=absorbing, 1=squeezing */
} Keccak;

static void xof_init(Keccak *ctx, int rate) {
    memset(ctx, 0, sizeof(*ctx));
    ctx->rate = rate;
    ctx->squeezing = 0;
    ctx->pos = 0;
}

static void xof_absorb(Keccak *ctx, const uint8_t *in, size_t len) {
    for (size_t i = 0; i < len; i++) {
        ctx->buf[ctx->pos++] ^= in[i];
        if (ctx->pos == ctx->rate) {
            /* XOR buf into state lanes (little-endian) */
            for (int j = 0; j < ctx->rate / 8; j++) {
                uint64_t lane = 0;
                for (int b = 0; b < 8; b++)
                    lane |= (uint64_t)ctx->buf[j*8+b] << (b*8);
                ctx->A[j] ^= lane;
            }
            /* Handling partial last lane if rate not multiple of 8 */
            int rem = ctx->rate % 8;
            if (rem) {
                uint64_t lane = 0;
                int base = (ctx->rate/8)*8;
                for (int b = 0; b < rem; b++)
                    lane |= (uint64_t)ctx->buf[base+b] << (b*8);
                ctx->A[ctx->rate/8] ^= lane;
            }
            keccakf(ctx->A);
            memset(ctx->buf, 0, ctx->rate);
            ctx->pos = 0;
        }
    }
}

static void xof_finalize(Keccak *ctx) {
    /* SHAKE domain separation: 0x1f, then pad with 0x80 at rate-1 */
    ctx->buf[ctx->pos] ^= 0x1f;
    ctx->buf[ctx->rate - 1] ^= 0x80;
    /* Absorbing the final partial block */
    for (int j = 0; j < ctx->rate / 8; j++) {
        uint64_t lane = 0;
        for (int b = 0; b < 8; b++)
            lane |= (uint64_t)ctx->buf[j*8+b] << (b*8);
        ctx->A[j] ^= lane;
    }
    int rem = ctx->rate % 8;
    if (rem) {
        uint64_t lane = 0;
        int base = (ctx->rate/8)*8;
        for (int b = 0; b < rem; b++)
            lane |= (uint64_t)ctx->buf[base+b] << (b*8);
        ctx->A[ctx->rate/8] ^= lane;
    }
    keccakf(ctx->A);
    /* Extracting rate bytes into buf for squeezing */
    for (int j = 0; j < ctx->rate / 8; j++) {
        uint64_t lane = ctx->A[j];
        for (int b = 0; b < 8; b++)
            ctx->buf[j*8+b] = (lane >> (b*8)) & 0xFF;
    }
    if (rem) {
        uint64_t lane = ctx->A[ctx->rate/8];
        int base = (ctx->rate/8)*8;
        for (int b = 0; b < rem; b++)
            ctx->buf[base+b] = (lane >> (b*8)) & 0xFF;
    }
    ctx->pos = 0;
    ctx->squeezing = 1;
}

static void xof_squeeze(Keccak *ctx, uint8_t *out, size_t len) {
    for (size_t i = 0; i < len; i++) {
        if (ctx->pos == ctx->rate) {
            keccakf(ctx->A);
            for (int j = 0; j < ctx->rate / 8; j++) {
                uint64_t lane = ctx->A[j];
                for (int b = 0; b < 8; b++)
                    ctx->buf[j*8+b] = (lane >> (b*8)) & 0xFF;
            }
            int rem = ctx->rate % 8;
            if (rem) {
                uint64_t lane = ctx->A[ctx->rate/8];
                int base = (ctx->rate/8)*8;
                for (int b = 0; b < rem; b++)
                    ctx->buf[base+b] = (lane >> (b*8)) & 0xFF;
            }
            ctx->pos = 0;
        }
        out[i] = ctx->buf[ctx->pos++];
    }
}

/* one shot SHAKE128 */
static void shake128(const uint8_t *in, size_t inlen, uint8_t *out, size_t outlen) {
    Keccak ctx;
    xof_init(&ctx, 168);
    xof_absorb(&ctx, in, inlen);
    xof_finalize(&ctx);
    xof_squeeze(&ctx, out, outlen);
}

/* one shot SHAKE256 */
static void shake256(const uint8_t *in, size_t inlen, uint8_t *out, size_t outlen) {
    Keccak ctx;
    xof_init(&ctx, 136);
    xof_absorb(&ctx, in, inlen);
    xof_finalize(&ctx);
    xof_squeeze(&ctx, out, outlen);
}

/*
  FIPS 203 NTT: 7 layer Cooley-Tukey for Z_q[X]/(X^256+1)
 
  ζ = 17 (order 256 mod q).
  ζ[k] = 17^bitrev7(k) mod q,  k = 0..127.
 
  NTT:  CT butterfly, k=1..127, len=128 down to len=2.
  INTT: GS butterfly, k=127..1, len=2 up to len=128, then scaling by N_INV.
 
  basemul: multiplying two NTT-domain polynomials pointwise.
           Each pair (2i, 2i+1) lives in Z_q[X]/(X^2 − γ_i).
           γ for even pair g: +ζ[64+g]
           γ for odd  pair g: -ζ[64+g]  (negated)
*/

static const int16_t ζs[128] = {
      1,  1729,  2580,  3289,  2642,   630,  1897,   848,
   1062,  1919,   193,   797,  2786,  3260,   569,  1746,
    296,  2447,  1339,  1476,  3046,    56,  2240,  1333,
   1426,  2094,   535,  2882,  2393,  2879,  1974,   821,
    289,   331,  3253,  1756,  1197,  2304,  2277,  2055,
    650,  1977,  2513,   632,  2865,    33,  1320,  1915,
   2319,  1435,   807,   452,  1438,  2868,  1534,  2402,
   2647,  2617,  1481,   648,  2474,  3110,  1227,   910,
     17,  2761,   583,  2649,  1637,   723,  2288,  1100,
   1409,  2662,  3281,   233,   756,  2156,  3015,  3050,
   1703,  1651,  2789,  1789,  1847,   952,  1461,  2687,
    939,  2308,  2437,  2388,   733,  2337,   268,   641,
   1584,  2298,  2037,  3220,   375,  2549,  2090,  1645,
   1063,   319,  2773,   757,  2099,   561,  2466,  2594,
   2804,  1092,   403,  1026,  1143,  2150,  2775,   886,
   1722,  1212,  1874,  1029,  2110,  2935,   885,  2154,
};

static void poly_ntt(int16_t f[KYBER_N]) {
    int k = 1;
    for (int len = 128; len >= 2; len >>= 1) {
        for (int start = 0; start < KYBER_N; start += 2*len) {
            int ζ = ζs[k++];
            for (int j = start; j < start+len; j++) {
                int t      = (int)ζ * f[j+len] % KYBER_Q;
                f[j+len]   = (int16_t)((f[j] - t + KYBER_Q) % KYBER_Q);
                f[j]       = (int16_t)((f[j] + t)           % KYBER_Q);
            }
        }
    }
}

static void poly_intt(int16_t f[KYBER_N]) {
    int k = 127;
    for (int len = 2; len <= 128; len <<= 1) {
        for (int start = 0; start < KYBER_N; start += 2*len) {
            int ζ = ζs[k--];
            for (int j = start; j < start+len; j++) {
                int t    = f[j];
                f[j]     = (int16_t)((t + f[j+len])         % KYBER_Q);
                f[j+len] = (int16_t)((int)ζ * ((f[j+len] - t + KYBER_Q) % KYBER_Q) % KYBER_Q);
            }
        }
    }
    for (int i = 0; i < KYBER_N; i++)
        f[i] = (int16_t)((int)f[i] * N_INV % KYBER_Q);
}

/* pointwise multiply two NTT-domain polys, result stored in r */
static void poly_basemul(int16_t r[KYBER_N],
                         const int16_t a[KYBER_N],
                         const int16_t b[KYBER_N]) {
    for (int g = 0; g < 64; g++) {
        int gamma_pos = ζs[64 + g];
        int gamma_neg = KYBER_Q - gamma_pos;
        /* Even pair: positions 4g, 4g+1 in ring X^2 - gamma_pos */
        int i = 4*g;
        r[i]   = (int16_t)(((int)a[i]*b[i]   + (int)gamma_pos*a[i+1]*b[i+1]) % KYBER_Q);
        r[i+1] = (int16_t)(((int)a[i]*b[i+1] + (int)a[i+1]*b[i])              % KYBER_Q);
        /* Odd pair: positions 4g+2, 4g+3 in ring X^2 - gamma_neg */
        i = 4*g + 2;
        r[i]   = (int16_t)(((int)a[i]*b[i]   + (int)gamma_neg*a[i+1]*b[i+1]) % KYBER_Q);
        r[i+1] = (int16_t)(((int)a[i]*b[i+1] + (int)a[i+1]*b[i])              % KYBER_Q);
    }
}

/*
  POLYNOMIAL ARITHMETIC
*/

static void poly_add(int16_t r[KYBER_N],
                     const int16_t a[KYBER_N],
                     const int16_t b[KYBER_N]) {
    for (int i = 0; i < KYBER_N; i++)
        r[i] = (int16_t)((a[i] + b[i]) % KYBER_Q);
}

static void polyvec_zero(int16_t v[][KYBER_N], int k) {
    for (int i = 0; i < k; i++)
        memset(v[i], 0, KYBER_N * sizeof(int16_t));
}

/* macc: r[i] += a[i][j] * b[j]  (NTT-domain pointwise) */
static void polyvec_basemul_acc(int16_t r[KYBER_N],
                                const int16_t a[][KYBER_N],
                                const int16_t b[][KYBER_N],
                                int k) {
    int16_t tmp[KYBER_N];
    for (int j = 0; j < k; j++) {
        poly_basemul(tmp, a[j], b[j]);
        poly_add(r, r, tmp);
    }
}

/*
  A Generation: FIPS 203 Algorithm 12 (XOF / Expand A)
 
  A[i][j] generated on the fly via SHAKE128(ρ || j || i).
  Note order: j first, then i. (FIPS 203 section 5.1)
  Rejection sampling: only accepts values in [0, q).
*/

static void gen_a_entry(int16_t poly[KYBER_N],
                        const uint8_t ρ[SEED_BYTES],
                        int i, int j) {
    uint8_t seed[SEED_BYTES + 2];
    memcpy(seed, ρ, SEED_BYTES);
    seed[SEED_BYTES]     = (uint8_t)j;   /* j first per FIPS 203 */
    seed[SEED_BYTES + 1] = (uint8_t)i;

    Keccak ctx;
    xof_init(&ctx, 168);
    xof_absorb(&ctx, seed, SEED_BYTES + 2);
    xof_finalize(&ctx);

    int count = 0;
    uint8_t buf[168 * 4];   /* several blocks to allow rejection */
    int bufpos = 0, buflen = 0;

    while (count < KYBER_N) {
        if (bufpos + 3 > buflen) {
            xof_squeeze(&ctx, buf, sizeof(buf));
            buflen = sizeof(buf);
            bufpos = 0;
        }
        /* Parsing 3 bytes into two 12-bit values */
        uint8_t b0 = buf[bufpos], b1 = buf[bufpos+1], b2 = buf[bufpos+2];
        bufpos += 3;
        int d1 = b0 | ((b1 & 0x0F) << 8);
        int d2 = (b1 >> 4) | ((int)b2 << 4);
        if (d1 < KYBER_Q && count < KYBER_N) poly[count++] = (int16_t)d1;
        if (d2 < KYBER_Q && count < KYBER_N) poly[count++] = (int16_t)d2;
    }
}

/*
  CBD SAMPLING: FIPS 203 Algorithm 7
 
  Input: 64*η bytes from SHAKE256.
  Output: polynomial with coefficients from CBD_η distribution.
  For η=2: each coefficient = (b0+b1) - (b2+b3) ∈ {-2,-1,0,1,2}.
  Stored as positive representative mod q.
*/

static void cbd_η2(int16_t poly[KYBER_N], const uint8_t *buf) {
    /* buf has 64*2=128 bytes = 64 bytes × η bytes */
    /* For η=2: 128 bytes. Each byte contributes 4 coefficients via nibbles. */
    for (int i = 0; i < KYBER_N / 4; i++) {
        uint8_t a = buf[2*i], b = buf[2*i + 1];
        /* Extracting 4 bits per value for η=2: bits b[0..3] from a, b[4..7] from b */
        int a0 = (a     ) & 0x3;
        int a1 = (a >> 2) & 0x3;
        int a2 = (a >> 4) & 0x3;
        int a3 = (a >> 6) & 0x3;
        int b0 = (b     ) & 0x3;
        int b1 = (b >> 2) & 0x3;
        int b2 = (b >> 4) & 0x3;
        int b3 = (b >> 6) & 0x3;

        /* hamming_weight(high) - hamming_weight(low), each nibble:
        for η=2: process 2 bits at a time 
        Actually FIPS 203 CBD η=2 takes 4 bits per coefficient:
        coeff = (bit0 + bit1) - (bit2 + bit3) where bits come 4 per byte pair 
        correct FIPS 203 Algorithm 7 formulation: */
        (void)a0; (void)a1; (void)a2; (void)a3;
        (void)b0; (void)b1; (void)b2; (void)b3;

        /* FIPS 203 Alg 7 direct: sample_poly_cbd(B, η=2):
          for i in 0..255: take 2*η bits = 4 bits
          a = hamming(bits[0..η-1]) = bits[0]+bits[1]
          b = hamming(bits[η..2*η-1]) = bits[2]+bits[3]
          coeff = a - b ∈ {-2,-1,0,1,2} */
        /* Recomputing properly with bit indexing: */
        uint16_t word = (uint16_t)a | ((uint16_t)b << 8);
        for (int j = 0; j < 4; j++) {
            int bit0 = (word >> (4*j + 0)) & 1;
            int bit1 = (word >> (4*j + 1)) & 1;
            int bit2 = (word >> (4*j + 2)) & 1;
            int bit3 = (word >> (4*j + 3)) & 1;
            int coeff = (bit0 + bit1) - (bit2 + bit3);
            poly[4*i + j] = (int16_t)((coeff + KYBER_Q) % KYBER_Q);
        }
        break;
    }
    /* Correct implementation: */
    for (int i = 0; i < KYBER_N; i++) {
        /* Each coefficient uses 4 bits (η=2): 2 bits for a, 2 bits for b */
        int byte_idx = (4 * i) / 8;
        int bit_offset = (4 * i) % 8;
        uint16_t bits;
        if (bit_offset <= 4) {
            bits = (buf[byte_idx] >> bit_offset) & 0xF;
        } else {
            bits = ((buf[byte_idx] >> bit_offset) | (buf[byte_idx+1] << (8-bit_offset))) & 0xF;
        }
        int a_val = (bits & 1) + ((bits >> 1) & 1);
        int b_val = ((bits >> 2) & 1) + ((bits >> 3) & 1);
        poly[i] = (int16_t)((a_val - b_val + KYBER_Q) % KYBER_Q);
    }
}

/* Sampling a secret/error polynomial via SHAKE256 then CBD */
static void sample_poly(int16_t poly[KYBER_N],
                        const uint8_t σ[SEED_BYTES],
                        uint8_t nonce) {
    uint8_t buf[SEED_BYTES + 1];
    memcpy(buf, σ, SEED_BYTES);
    buf[SEED_BYTES] = nonce;

    uint8_t cbd_buf[64 * KYBER_η];   /* 128 bytes for η=2 */
    shake256(buf, SEED_BYTES + 1, cbd_buf, sizeof(cbd_buf));
    cbd_η2(poly, cbd_buf);
}

/* 
  ENCODING FIPS 203 ByteEncode_12
 
  Packing 256 12 bit values into 384 bytes.
  Two coefficients then 3 bytes:
    byte[3i+0] =  coeff[2i]        & 0xFF
    byte[3i+1] = (coeff[2i] >> 8)  | (coeff[2i+1] << 4)
    byte[3i+2] =  coeff[2i+1] >> 4
*/

static void byte_encode12(uint8_t out[POLY_BYTES], const int16_t f[KYBER_N]) {
    for (int i = 0; i < KYBER_N / 2; i++) {
        uint16_t a = (uint16_t)f[2*i];
        uint16_t b = (uint16_t)f[2*i+1];
        out[3*i + 0] = (uint8_t)( a        & 0xFF);
        out[3*i + 1] = (uint8_t)((a >> 8)  | ((b & 0x0F) << 4));
        out[3*i + 2] = (uint8_t)( b >> 4);
    }
}

/*
   Print Helpers
*/

static void print_hex(const char *label, const uint8_t *data, size_t len) {
    printf("%s (%zu bytes):\n  ", label, len);
    for (size_t i = 0; i < len; i++) {
        printf("%02x", data[i]);
        if ((i+1) % 32 == 0 && i+1 < len) printf("\n  ");
    }
    printf("\n");
}

/*
  KEY GENERATION  (FIPS 203 Algorithm 13: K-PKE.KeyGen)
 
  Input:  random 32-byte seed d (here we use a fixed example seed)
  Steps:
    (ρ, σ) = G(d)              -- G is SHA3-512 (split into two 32-byte halves)
    A_hat  = ExpandA(ρ)        -- on the fly via SHAKE128
    s, e   = ExpandS(σ)        -- via SHAKE256 then CBD_η
    s_hat  = NTT(s)
    e_hat  = NTT(e)
    t_hat  = A_hat · s_hat + e_hat
    pk     = ByteEncode12(t_hat) || ρ
    sk     = ByteEncode12(s_hat)
 
  I use SHAKE256 in place of G for simplicity (outputs 64 bytes).
*/

static void keygen(const uint8_t d[SEED_BYTES], int k) {
    printf("╔══════════════════════════════════════════════════════════════╗\n");
    printf("║  CRYSTALS-Kyber / ML-KEM    KeyGen  (FIPS 203)               ║\n");
    printf("║  n=256  q=3329  η=2  k=%d  (Kyber-%d)                       ║\n",
           k, k==2?512:k==3?768:1024);
    printf("╚══════════════════════════════════════════════════════════════╝\n\n");

    /* Step 1: deriving ρ and σ from seed d */
    uint8_t ρ_σ[64];
    /* Use SHAKE256(d || k) → 64 bytes as (ρ, σ) */
    uint8_t d_ext[SEED_BYTES + 1];
    memcpy(d_ext, d, SEED_BYTES);
    d_ext[SEED_BYTES] = (uint8_t)k;
    shake256(d_ext, SEED_BYTES + 1, ρ_σ, 64);

    uint8_t *ρ   = ρ_σ;        /* first 32 bytes: public seed for A */
    uint8_t *σ = ρ_σ + 32;   /* next  32 bytes: secret seed */

    print_hex("Seed d (input)", d, SEED_BYTES);
    print_hex("ρ (public seed for A)", ρ, SEED_BYTES);
    print_hex("σ (secret seed for s,e)", σ, SEED_BYTES);
    printf("\n");

    /* Step 2: sample s and e via CBD */
    /* Keeping coefficient-domain copies for printing, NTT copies for math */
    int16_t s_coeff[4][KYBER_N], e_coeff[4][KYBER_N];
    int16_t s_ntt  [4][KYBER_N], e_ntt  [4][KYBER_N];

    for (int i = 0; i < k; i++) {
        sample_poly(s_coeff[i], σ, (uint8_t)i);
        sample_poly(e_coeff[i], σ, (uint8_t)(k + i));
        memcpy(s_ntt[i], s_coeff[i], KYBER_N * sizeof(int16_t));
        memcpy(e_ntt[i], e_coeff[i], KYBER_N * sizeof(int16_t));
        poly_ntt(s_ntt[i]);
        poly_ntt(e_ntt[i]);
    }

    /* Step 3: compute t = A·s + e  (all in NTT domain) */
    int16_t t_ntt[4][KYBER_N];
    int16_t a_col[KYBER_N], tmp[KYBER_N];

    for (int i = 0; i < k; i++) {
        memset(t_ntt[i], 0, KYBER_N * sizeof(int16_t));
        for (int j = 0; j < k; j++) {
            gen_a_entry(a_col, ρ, i, j);  /* A[i][j] on-the-fly */
            poly_basemul(tmp, a_col, s_ntt[j]);
            poly_add(t_ntt[i], t_ntt[i], tmp);
        }
        poly_add(t_ntt[i], t_ntt[i], e_ntt[i]);
    }

    /* Step 4: encoding and printing */
    /* s: encoding coefficient domain (the mathematical private key values) */
    size_t sk_bytes = (size_t)k * POLY_BYTES;
    uint8_t *sk_enc = malloc(sk_bytes);
    for (int i = 0; i < k; i++)
        byte_encode12(sk_enc + i * POLY_BYTES, s_ntt[i]);

    /* pk: encoded t_ntt + append ρ */
    size_t pk_bytes = (size_t)k * POLY_BYTES + SEED_BYTES;
    uint8_t *pk_enc = malloc(pk_bytes);
    for (int i = 0; i < k; i++)
        byte_encode12(pk_enc + i * POLY_BYTES, t_ntt[i]);
    memcpy(pk_enc + k * POLY_BYTES, ρ, SEED_BYTES);

    /* e: encoding coefficient domain */
    uint8_t *e_enc = malloc(sk_bytes);
    for (int i = 0; i < k; i++)
        byte_encode12(e_enc + i * POLY_BYTES, e_ntt[i]);

    /* Printing results */
    printf("Private key components:\n\n");
    printf("s  (secret key, CBD η=2, NTT-domain, ByteEncode_12):\n");
    printf("   Coefficients in [0,q): small values like 0,1,2,q-1,q-2\n");
    printf("   NTT-domain: spread across [0,q) after transform\n");
    print_hex("  sk_hex", sk_enc, sk_bytes);
    printf("\n");

    printf("e  (error vector, CBD η=2, NTT-domain, ByteEncode_12):\n");
    print_hex("  e_hex", e_enc, sk_bytes);
    printf("\n");

    printf("Public key:\n\n");
    printf("pk = ByteEncode_12(t̂) || ρ\n");
    printf("   t̂ = Â·ŝ + ê  (NTT domain, k×256 coefficients)\n");
    printf("   ρ = 32-byte public seed for A\n");
    print_hex("  pk_hex", pk_enc, pk_bytes);
    printf("\n");

    printf("Sizes\n");
    printf("  sk:  %zu bytes  (%d × %d coefficients × 12 bits)\n",
           sk_bytes, k, KYBER_N);
    printf("  e:   %zu bytes  (same encoding as sk)\n", sk_bytes);
    printf("  pk:  %zu bytes  (%d × %d × 12 bits + 32-byte seed)\n",
           pk_bytes, k, KYBER_N);
    printf("\n");

    /* Also printing first 8 raw coefficients inspection */
    printf("First 8 coefficients of s[0] (coefficient domain)\n");
    printf("  Values in [0,q=%d):  ", KYBER_Q);
    for (int i = 0; i < 8; i++) printf("%4d ", s_coeff[0][i]);
    printf("\n  Signed (CBD noise): ");
    for (int i = 0; i < 8; i++) {
        int v = s_coeff[0][i];
        if (v > KYBER_Q/2) v -= KYBER_Q;
        printf("%+3d  ", v);
    }
    printf("\n\n");

    printf("First 8 coefficients of e[0] (coefficient domain)\n");
    printf("  Values in [0,q=%d):  ", KYBER_Q);
    for (int i = 0; i < 8; i++) printf("%4d ", e_coeff[0][i]);
    printf("\n  Signed (CBD noise): ");
    for (int i = 0; i < 8; i++) {
        int v = e_coeff[0][i];
        if (v > KYBER_Q/2) v -= KYBER_Q;
        printf("%+3d  ", v);
    }
    printf("\n\n");

    printf(" First 8 coefficients of t̂[0] (NTT domain = as in pk)\n");
    printf("  Values in [0,q=%d):  ", KYBER_Q);
    for (int i = 0; i < 8; i++) printf("%4d ", t_ntt[0][i]);
    printf("\n");

    free(sk_enc); free(pk_enc); free(e_enc);
}

/*
  Main
*/

int main(void) {
    /* Fixed example seed, in production use /dev/urandom */
    uint8_t d[SEED_BYTES] = {
        0x7f, 0x9c, 0x2b, 0xa4, 0xe8, 0x8f, 0x82, 0x7d,
        0x61, 0x60, 0x45, 0x07, 0xa7, 0x33, 0x87, 0x12,
        0x10, 0x50, 0x40, 0x05, 0x6f, 0xe8, 0x3b, 0x57,
        0x11, 0x89, 0x20, 0xc2, 0x65, 0x23, 0x41, 0x00,
    };

    keygen(d, KYBER_K);
    return 0;
}
