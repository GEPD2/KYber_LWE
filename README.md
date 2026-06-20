# CRYSTALS-Kyber / ML-KEM (FIPS 203) in C

A single-file, zero-dependency C implementation of the **full CRYSTALS-Kyber Key Encapsulation Mechanism**, compliant with [NIST FIPS 203](https://csrc.nist.gov/pubs/fips/203/final) (ML-KEM).

It implements the complete IND-CCA2 KEM **KeyGen, Encapsulation, and Decapsulation** with the Fujisaki–Okamoto transform and not just key generation. It supports all three security levels: **ML-KEM-512**, **ML-KEM-768**, and **ML-KEM-1024**.

> **Status:** built from scratch, no OpenSSL / libsodium / liboqs. Passes its built-in self-tests (hash known-answer vectors, NTT correctness, KEM round-trip, implicit rejection) on all three parameter sets. Educational reference, **not** constant-time, not hardened for production.

---

## What it does

Running the program executes two phases:

1. **Self-tests** (abort on any failure) that prove the implementation is correct without needing the network go see [Self-tests](#self-tests).
2. A **deterministic KEM exchange**: KeyGen → Encaps → Decaps, printing the public key, the ciphertext, and both sides' shared secret so you can confirm they match.

| Output | Description |
|--------|-------------|
| `ek` | Encapsulation (public) key = ByteEncode₁₂(t̂) ‖ ρ |
| `ct` | Ciphertext = Compress(u) ‖ Compress(v) |
| `K` / `K'` | The 32-byte shared secret, computed by the Encaps and Decaps sides |

Both shared secrets are printed and compared; a correct build prints `Shared secrets match: YES`.

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

No external libraries. No Makefile needed. Compiles cleanly with `-Wall` (zero warnings).

---

## Security level

Change one line at the top of the file:

```c
#define KYBER_K  2   // ML-KEM-512   — ek=800B,  dk=1632B, ct=768B
#define KYBER_K  3   // ML-KEM-768   — ek=1184B, dk=2400B, ct=1088B
#define KYBER_K  4   // ML-KEM-1024  — ek=1568B, dk=3168B, ct=1568B
```

Everything else (`n`, `q`, the per-set η₁/η₂ and dᵤ/dᵥ) is derived automatically from `KYBER_K`, exactly as in the FIPS 203 parameter table.

---

## Parameters

| Parameter | Value | Note |
|-----------|-------|------|
| `n` | 256 | Polynomial degree, fixed by FIPS 203 |
| `q` | 3329 | Prime modulus, fixed by FIPS 203 |
| `η₁` | **3** for ML-KEM-512, **2** for 768/1024 | KeyGen noise width (set-dependent) |
| `η₂` | 2 | Encryption noise width (all sets) |
| `dᵤ`, `dᵥ` | 10/4, 10/4, 11/5 | Compression bit-widths per set |
| `k` | 2 / 3 / 4 | Module rank, selects the security level |
| `N_INV` | 3303 | 128⁻¹ mod 3329, the INTT scaling factor |

> **Note on η₁:** it is *not* a single fixed value. ML-KEM-512 widens its key-generation noise to η₁ = 3 to compensate for its smaller module rank; 768 and 1024 use η₁ = 2. η₂ = 2 for all sets.

---

## Architecture

Everything is implemented from scratch in a single file. There are no dependencies on OpenSSL, libsodium, or any other library.

### Keccak-f[1600] sponge

The file implements the full Keccak-f[1600] permutation (24 rounds, θ/ρ/π/χ/ι) and builds **five** hash roles on top of it. THey are distinguished only by their rate and domain-separation byte:

| Function | Rate | Domain | Role in ML-KEM |
|----------|------|--------|----------------|
| SHAKE128 | 168 | `0x1F` | XOF — expand the matrix **A** from ρ |
| SHAKE256 | 136 | `0x1F` | PRF — sample CBD noise; also **J** (rejection key) |
| SHA3-256 | 136 | `0x06` | **H** — bind ek into the FO transform |
| SHA3-512 | 72  | `0x06` | **G** — derive (ρ, σ) and (K, r) |

### NTT: why it differs from a simple negacyclic NTT

This is the most technically subtle part. A simple negacyclic NTT needs a primitive 2n-th root of unity (ψ²ⁿ = 1 **and** ψⁿ = −1). For n = 256 that would require a primitive 512th root.

The problem: q − 1 = 3328 = 2⁸ × 13. Since 512 = 2⁹ does **not** divide 3328, no primitive 512th root exists mod 3329. The value ζ = 17 has multiplicative order 256 (so ζ²⁵⁶ = 1, not −1), so the simple 8-layer negacyclic approach cannot be used.

**The FIPS 203 solution is a 7-layer NTT with base multiplication:**

1. **7 Cooley–Tukey butterfly layers** (len = 128 down to 2) transform the polynomial into 128 quadratic factors, each in a subring Zq[X]/(X² − γᵢ).
2. **`poly_basemul`** completes the pointwise product within each quadratic subring, using the appropriate ±γ from the zeta table. *(It uses 64-bit intermediates: γ·a·b can reach q³ ≈ 3.7×10¹⁰, which would overflow a 32-bit int.)*
3. **7 Gentleman–Sande butterfly layers** (INTT) invert the transform, then every coefficient is scaled by **128⁻¹** mod q (not 256⁻¹, because only 7 layers were applied).

The 128-entry zeta table is ζs[k] = 17^bitrev7(k) mod 3329 (k = 0..127), in the **plain (non-Montgomery)** domain to match the plain `% q` arithmetic used throughout.

### Matrix A generated on the fly

The public matrix A ∈ Rq^(k×k) is **never stored whole**. Each entry is generated on demand:

```
A[i][j] = SampleNTT( SHAKE128(ρ ‖ j ‖ i) )   — rejection-sample 12-bit values < q
```

Note the byte order `j ‖ i` (column first, then row), as specified in FIPS 203 §5.1.

### CBD sampling

The secret **s** and error **e** polynomials come from the Centered Binomial Distribution (FIPS 203 Algorithm 8):

1. `SHAKE256(σ ‖ nonce)` produces 64·η bytes of pseudorandom material (128 B for η = 2, 192 B for η = 3).
2. Each coefficient uses 2·η bits: `coeff = (b₀+…+b_{η−1}) − (b_η+…+b_{2η−1})`.
3. The result lands in {−η,…,+η}, stored as its positive representative mod q.

### Compression and ByteEncode

Ciphertext components are **compressed** (lossy, dropping low bits to dᵤ/dᵥ bits per coefficient) and decryption still works because the accumulated noise stays below q/4. Coefficients are then packed with `ByteEncode_d`; the 12-bit key vectors pack two coefficients per three bytes (384 bytes per polynomial).

### K-PKE — the IND-CPA core

```
KeyGen:   t̂ = Â·ŝ + ê                              // s, e small (CBD)
Encrypt:  u = INTT(Âᵀ·r̂) + e₁                        // r, e1, e2 fresh CBD
          v = INTT(t̂ᵀ·r̂) + e₂ + Decompress₁(m)       // message scaled by q/2
Decrypt:  w = v − sᵀu  →  round each coeff to a bit
```

### ML-KEM the FO transform (IND-CCA2)

The KEM wraps K-PKE so a tampered ciphertext cannot leak the secret:

```
KeyGen:  dk = dkPKE ‖ ek ‖ H(ek) ‖ z
Encaps:  (K, r) = G(m ‖ H(ek));   c = Encrypt(ek, m, r);   output (K, c)
Decaps:  m' = Decrypt(s, c);  (K', r') = G(m' ‖ H(ek));
         re-encrypt c' = Encrypt(ek, m', r');
         if c' == c:  K = K'            // honest sender
         else:        K = J(z ‖ c)      // implicit rejection
```

**Implicit rejection** is the key CCA2 defence: on any mismatch, Decaps returns a deterministic-but-unpredictable key `J(z‖c)` rather than failing, so an attacker learns nothing from probing.

---

## Self-tests

Five checks run at startup and abort on failure:

1. **Hash KATs** SHA3-256/512 and SHAKE128/256 against published empty-string vectors. This is the only test that catches a *wrong-but-consistent* hash (one that round-trips with itself but doesn't match the standard).
2. **NTT ∘ INTT = identity** over random polynomials.
3. **`basemul` vs schoolbook** the fast NTT-domain product must equal a direct negacyclic convolution in Rq.
4. **KEM round-trip** Decaps(Encaps) must recover the identical shared secret.
5. **Implicit rejection** a tampered ciphertext must yield a different key.

---

## Sizes

| Set | k | Public key `ek` | Decaps key `dk` | Ciphertext `ct` | Shared secret |
|-----|---|-----------------|-----------------|-----------------|---------------|
| ML-KEM-512  | 2 | 800 B  | 1632 B | 768 B  | 32 B |
| ML-KEM-768  | 3 | 1184 B | 2400 B | 1088 B | 32 B |
| ML-KEM-1024 | 4 | 1568 B | 3168 B | 1568 B | 32 B |

(`dk = dkPKE ‖ ek ‖ H(ek) ‖ z`. The bare K-PKE secret key `s` is 768/1152/1536 B.)

---

## Example output (ML-KEM-512, k = 2)

```
ML-KEM-512  (FIPS 203, k=2 eta1=3 eta2=2 du=10 dv=4)
Self-tests:
  [PASS] hash KATs: SHA3-256/512, SHAKE128/256 (empty-string vectors)
  [PASS] NTT then INTT = identity (50 polys)
  [PASS] INTT(basemul(NTT a,NTT b)) = a*b in R_q (50 trials)
  [PASS] KEM round-trip: Decaps(Encaps) recovers shared secret (20 trials)
  [PASS] implicit rejection: tampered ct -> different key (20 trials)
All self-tests passed.

 ML-KEM-512 full KEM (deterministic demo) 
Sizes: ek=800  dk=1632  ct=768  ss=32

ek (first 192 bytes):
  c78209aebc5eb37bc0203254d4f87736f323b5b35687b43decb6a23be45a1916
  ... (800 bytes total)

ct (first 192 bytes):
  52d7007b92725e71faa4a8fa7b90c2f4d472e28d848906a7bf74086b6b141d13
  ... (768 bytes total)

K  (shared secret, Encaps side) (32 bytes):
  39bbbea637d43c9fd825388d63ba93cfb1d200dd612fd08b4e30a84aeee23c43
K' (shared secret, Decaps side) (32 bytes):
  39bbbea637d43c9fd825388d63ba93cfb1d200dd612fd08b4e30a84aeee23c43

Shared secrets match: YES
```

---

## Changing the seed

The seeds in `main()` are fixed for reproducibility. For real use, draw `d`, `z`, and `m` from a cryptographically secure RNG:

```c
// Linux / macOS
FILE *f = fopen("/dev/urandom", "rb");
fread(d, 1, 32, f);   fread(z, 1, 32, f);   fread(m, 1, 32, f);
fclose(f);
```

---

## Mathematical background

Kyber's security rests on the hardness of the **Module Learning With Errors (MLWE)** problem:

> Given a public matrix **A** and a public vector **t = As + e**, recovering the secret **s** is computationally infeasible, even for a quantum computer.

Shor's algorithm breaks RSA and elliptic-curve cryptography because it efficiently solves the *abelian hidden-subgroup / period-finding* problem behind factorisation and discrete logs. MLWE simply **isn't** such a problem, there is no period for the Quantum Fourier Transform to extract which is why no efficient quantum attack on it is known. The small noise **e** is what keeps the secret hidden.

| Scheme | Hardness assumption | Quantum threat |
|--------|---------------------|----------------|
| RSA / DH | Factorisation / DLP | Broken by Shor |
| ECDSA / Schnorr | Elliptic-curve DLP | Broken by Shor |
| CRYSTALS-Kyber | Module-LWE | No known quantum attack |

---

## Recent corrections

This version fixes several issues present in the earlier key-generation-only release:

- **Keccak-f ρ/π step** corrected to the canonical chained walk. The previous version produced hashes that round-tripped with themselves but did **not** match the SHA-3/SHAKE standard (caught by the new hash KAT). All hash-derived output (ρ, σ, A, keys) was affected.
- **`basemul` overflow** now uses 64-bit intermediates; the previous 32-bit code overflowed (γ·a·b ≈ q³) and silently corrupted products.
- **G = SHA3-512** key derivation now uses the correct function. (Earlier versions used SHAKE256 as a stand-in; this is no longer the case.)
- **η₁ is now set-dependent** 3 for ML-KEM-512, 2 for 768/1024. The earlier code hardcoded η = 2, which was non-compliant for ML-KEM-512.
- **Completed to a full KEM** added K-PKE Encrypt/Decrypt and ML-KEM Encaps/Decaps with the Fujisaki–Okamoto transform and implicit rejection.
- **Built-in self-tests** hash KATs, NTT correctness, KEM round-trip, and implicit rejection now run on every invocation.

**Verification scope:** the implementation is validated by the published SHA-3/SHAKE KAT vectors, internal algebraic consistency, the end-to-end KEM round-trip, and exact FIPS 203 key/ciphertext sizes across all three sets. It has **not** yet been checked against NIST's official ACVP/ML-KEM KAT response files.

---

## References

- [NIST FIPS 203](https://csrc.nist.gov/pubs/fips/203/final) ML-KEM Standard (August 2024)
- [CRYSTALS-Kyber specification](https://pq-crystals.org/kyber/) original NIST PQC submission
- [NIST FIPS 202](https://csrc.nist.gov/publications/detail/fips/202/final) SHA-3 / SHAKE (Keccak)
- Bos et al. (2018) *CRYSTALS-Kyber: a CCA-secure module-lattice-based KEM*, IEEE EuroS&P

---

## License

Written for educational purposes, to make the internals of FIPS 203 visible and understandable. It is **not** hardened for production use (no constant-time guarantees, no side-channel mitigations). See [LICENSE](https://github.com/GEPD2/KYber_LWE/blob/main/LICENSE).
