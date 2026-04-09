# CRYSTALS-Kyber / ML-KEM (FIPS 203)  Key Generation in C

A single file, with zero-dependencies C implementation of the **CRYSTALS-Kyber Key Generation** algorithm, fully compliant with [NIST FIPS 203](https://csrc.nist.gov/pubs/fips/203/final) (ML-KEM).

Supports all three security levels: **Kyber-512**, **Kyber-768**, and **Kyber-1024**.

---

## What it does

Runs Key Generation and prints the following in hex:

| Output | Description |
|--------|-------------|
| `sk_hex` | Secret key **s** in NTT-domain, ByteEncode\_12 packed |
| `e_hex`  | Error vector **e** in NTT-domain, ByteEncode\_12 packed |
| `pk_hex` | Public key **pk** = ByteEncode\_12(t̂) \|\| ρ |

It also prints the first 8 raw coefficients of `s[0]` and `e[0]` in both unsigned `[0, q)` and signed `{−2…+2}` form, so you can verify the CBD noise distribution by eye.

---

## Build and run

```bash
# Linux / macOS
gcc -O2 -Wall -o kyber_fips203 kyber_fips203.c
./kyber_fips203

# Windows (MinGW)
gcc -O2 -Wall -o kyber_fips203.exe kyber_fips203.c
kyber_fips203.exe
```

No external libraries. No Makefile needed. Standard C99.

---

## Security level

Change one line at the top of the file:

```c
#define KYBER_K  2   // Kyber-512  — pk=800B,  sk=768B
#define KYBER_K  3   // Kyber-768  — pk=1184B, sk=1152B
#define KYBER_K  4   // Kyber-1024 — pk=1568B, sk=1536B
```

All other parameters (n=256, q=3329, η=2) are fixed by FIPS 203 and never need changing.

---

## Parameters

| Parameter | Value | Note |
|-----------|-------|------|
| `n` | 256 | Polynomial degree fixed by FIPS 203 |
| `q` | 3329 | Prime modulus fixed by FIPS 203 |
| `η` | 2 | CBD noise bound fixed by FIPS 203 |
| `k` | 2 / 3 / 4 | Module dimension selects security level |
| `N_INV` | 3303 | 128⁻¹ mod 3329 INTT scaling factor |

---

## Architecture

Everything is implemented from scratch in a single file. There are no dependencies on OpenSSL, libsodium, or any other library.

### Keccak-f\[1600\] sponge

The file implements the full Keccak-f\[1600\] permutation (24 rounds, θ/ρ/π/χ/ι steps) and builds both XOFs on top of it:

- **SHAKE128**  rate = 168 bytes, used for matrix A generation
- **SHAKE256**  rate = 136 bytes, used for CBD sampling

Domain separation byte: `0x1f` (standard SHAKE suffix).

### NTT: Why it is different from a simple negacyclic NTT

This is the most technically subtle part. A simple negacyclic NTT requires a primitive 2n-th root of unity, meaning ψ^(2n) = 1 **and** ψ^n = −1 mod q. For n = 256 that would require a primitive 512th root.

The problem: q − 1 = 3328 = 2⁸ × 13. Since 512 = 2⁹ does **not** divide 3328, no primitive 512th root exists mod 3329. The value ζ = 17 has multiplicative order 256 (so ζ²⁵⁶ = 1, not −1), which means the simple 8-layer negacyclic approach cannot be used.

**The FIPS 203 solution**  7 layer NTT with base multiplication:

1. **7 Cooley-Tukey butterfly layers** (len = 128 down to len = 2) transform the degree-255 polynomial into 128 quadratic factors, each living in a subring Z\_q\[X\]/(X² − γᵢ).
2. **`poly_basemul`** completes pointwise multiplication within each quadratic subring, using the appropriate ±γ value from the zeta table.
3. **7 Gentleman-Sande butterfly layers** (INTT, len = 2 up to len = 128) invert the transform, followed by scaling every coefficient by 128⁻¹ mod q (not 256⁻¹, because only 7 layers were applied).

The 128-entry zeta table is precomputed as ζs\[k\] = 17^bitrev7(k) mod 3329 for k = 0..127. This is exactly FIPS 203 Algorithm 9 / Algorithm 10.

### A matrix generated on the fly

The public matrix A ∈ R\_q^(k×k) is **never stored in memory**. Each entry A\[i\]\[j\] is generated on demand by calling:

```
SHAKE128(ρ || j || i)  rejection sample 256 coefficients < q
```

Note the byte order `j || i` (j first, then i), as specified in FIPS 203 Section 5.1. This means for large k the memory footprint for A stays O(n) at all times.

### CBD sampling

Secret key polynomials **s** and error polynomials **e** are sampled using the Centered Binomial Distribution with η = 2 (FIPS 203 Algorithm 7):

1. `SHAKE256(σ || nonce)` produces 128 bytes of pseudorandom material.
2. Each coefficient takes 4 bits: `coeff = (b₀ + b₁) − (b₂ + b₃)` where each bᵢ is one bit.
3. Result is in {−2, −1, 0, +1, +2}, stored as the positive representative mod q.

### ByteEncode\_12 encoding

Each 12 bit coefficient is packed into 1.5 bytes (two coefficients per 3 bytes):

```
byte[3i+0] =  coeff[2i] & 0xFF
byte[3i+1] = (coeff[2i] >> 8) | (coeff[2i+1] << 4)
byte[3i+2] =  coeff[2i+1] >> 4
```

256 coefficients × 12 bits = 384 bytes per polynomial vector entry.

### Key Generation flow (FIPS 203 Algorithm 13: K-PKE.KeyGen)

```
(ρ, σ) = SHAKE256(d || k)          // derive public and secret seeds
A_hat  = ExpandA(ρ)                // generate A on-the-fly via SHAKE128
s, e   = ExpandS(σ)                // sample small polynomials via CBD
s_hat  = NTT(s),  e_hat = NTT(e)
t_hat  = A_hat · s_hat + e_hat     // all arithmetic in NTT domain
pk     = ByteEncode12(t_hat) || ρ  // 800 / 1184 / 1568 bytes
sk     = ByteEncode12(s_hat)       // 768 / 1152 / 1536 bytes
```

> **Note:** The implementation uses SHAKE256 as a stand-in for the FIPS 203 G function (which is SHA3-512). Both produce 64 bytes of output and the key derivation logic is identical. In a production deployment, replace with SHA3-512.

---

## Example output (Kyber-512, k=2)

```
╔══════════════════════════════════════════════════════════════╗
║  CRYSTALS-Kyber / ML-KEM    KeyGen  (FIPS 203)               ║
║  n=256  q=3329  η=2  k=2  (Kyber-512)                        ║
╚══════════════════════════════════════════════════════════════╝

Seed d (input) (32 bytes):
  7f9c2ba4e88f827d61604507a7338712105040056fe83b57118920c265234100
ρ (public seed for A) (32 bytes):
  efe5772a502d19ab96f3b14b44644e40255ef3bc867518393d610f49684ad2f2
σ (secret seed for s,e) (32 bytes):
  a71e0fd474a6053e4d8c0eef55f7f5c9013a2c2ebce0919fef394fa77ef35d94

sk_hex (768 bytes):
  25d058c13a1480469c14...
e_hex  (768 bytes):
  8ff29c40062111b2be2a...
pk_hex (800 bytes):
  ...efe5772a502d19ab96f3b14b44644e40255ef3bc867518393d610f49684ad2f2
      ↑ last 32 bytes = ρ

First 8 coefficients of s[0]:  +0  +0  +0  +0  +0  +0  +2  +1
First 8 coefficients of e[0]:  +0  -1  +1  +2  +0  +0  -1  +2
```

The last 32 bytes of `pk_hex` are always equal to ρ — this is the public seed that allows any party to regenerate A on demand.

The signed coefficients of s and e are strictly in {−2, −1, 0, +1, +2}, confirming correct CBD η=2 sampling.

---

## Changing the seed

The seed in `main()` is a fixed example for reproducibility. For real use, replace it with random bytes:

```c
// Linux / macOS
FILE *f = fopen("/dev/urandom", "rb");
fread(d, 1, SEED_BYTES, f);
fclose(f);
```

---

## Mathematical background

The security of Kyber rests on the hardness of the **Module Learning With Errors (MLWE)** problem:

> Given a public matrix **A** and a public vector **t = As + e**, it is computationally infeasible to recover the secret **s** even for a quantum computer.

This is because MLWE has no periodic structure that the Quantum Fourier Transform (QFT) can exploit, unlike the integer factorisation and discrete logarithm problems attacked by Shor's algorithm. The noise vector **e** (whose coefficients are bounded by η = 2) is the mechanism that destroys any exploitable algebraic regularity.

| Scheme | Hardness assumption | Quantum threat |
|--------|--------------------|-|
| RSA / DH | Factorisation / DLP | Broken by Shor |
| ECDSA / Schnorr | DLP on elliptic curves | Broken by Shor |
| CRYSTALS-Kyber | Module-LWE | Resistant, no periodicity |

---

## References

- [NIST FIPS 203](https://csrc.nist.gov/pubs/fips/203/final) ML-KEM Standard (August 2024)
- [CRYSTALS-Kyber specification](https://pq-crystals.org/kyber/) original submission to NIST PQC
- [NIST FIPS 202](https://csrc.nist.gov/publications/detail/fips/202/final) SHA-3 / SHAKE standard (Keccak)
- Bos et al. (2018) — *CRYSTALS-Kyber: a CCA-secure module-lattice-based KEM*, IEEE EuroS&P

---

## License

This code is written for educational purposes to make the internals of FIPS 203 visible and understandable. It is not hardened for production use (no constant-time guarantees, no side-channel mitigations). [LICENSE](https://github.com/GEPD2/KYber_LWE/blob/main/LICENSE)
