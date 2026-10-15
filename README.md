# FoxMath 🦊

**Minimalist CLI & library for number theory and cryptography: RSA & Diffie-Hellman cryptosystems, BSGS discrete log solver, arithmetic functions (μ, σ, τ), Jacobi symbol & Solovay–Strassen test, Pollard's ρ factorization, totient & Carmichael functions, Tonelli-Shanks roots, interactive REPL, Weierstrass elliptic curves, Chinese Remainder Theorem, and continued fractions.**

**More about me / other projects:** [abhrankan.duckdns.org](https://abhrankan.duckdns.org)

One binary. Clean output. Zero footguns. Educational by design.

Built with curiosity — for students, tinkerers, and crypto/math enthusiasts.

## Features
- Textbook RSA cryptosystem: key generation, encryption, CRT-accelerated decryption, and Pollard's ρ factoring attack (rsa)
- Diffie-Hellman key exchange: key generation, shared secret derivation, and BSGS discrete log attack (dh)
- Legendre symbol
- Jacobi symbol (jacobi)
- Solovay–Strassen probabilistic primality test (solovay)
- Modular square root via Tonelli-Shanks (sqrtmod)
- Deterministic Miller-Rabin primality test (isprime)
- Prime factorization via Pollard's rho & trial division (factor)
- Discrete logarithm solver via Baby-step Giant-step (dlog)
- Euler's totient φ(n) and Carmichael λ(n) (totient)
- Multiplicative order and primitive roots (order)
- Classical arithmetic functions: Möbius μ(n), divisor count τ(n), divisor sum σ(n)
- Extended Euclidean algorithm (egcd) & Bézout coefficients
- Modular multiplicative inverse (modinv)
- Euler/Hermann Machin-like π approximation
- Interactive Python REPL with all functions preloaded (repl)
- Elliptic curve point addition (over finite fields), including
  named curves (secp256k1) with on-curve result validation
- Elliptic curve scalar multiplication k·P (double-and-add), including
  the point at infinity as a first-class result
- Chinese Remainder Theorem solver (handles non-coprime moduli)
- Continued fraction expansion and convergents
- Challenge/quiz mode — practice problems with scoring, or a
  worksheet mode for self-study (Legendre, Jacobi, CRT, SqrtMod, Factor, Totient, DLog, RSA, DH)
- Symlink-friendly (foxmath-rsa, foxmath-dh, foxmath-legendre, foxmath-jacobi,
  foxmath-solovay, foxmath-repl, foxmath-sqrtmod, foxmath-isprime,
  foxmath-factor, foxmath-totient, foxmath-order, foxmath-dlog,
  foxmath-mobius, foxmath-sigma, foxmath-tau, foxmath-egcd,
  foxmath-modinv, foxmath-pi, foxmath-ecadd, foxmath-ecmul,
  foxmath-crt, foxmath-cf, foxmath-challenge)
- JSON output mode

## Installation

```bash
pip install foxmath
```

Or from source:
```bash
git clone https://github.com/Abhrankan-Chakrabarti/foxmath.git
cd foxmath
pip install -e .
```

Check which version you have installed:
```bash
foxmath --version
```
```
foxmath 0.14.0
```

## Usage Examples

```bash
# Toy RSA: keygen, encrypt, decrypt (with CRT), and Pollard's rho crack
foxmath rsa keygen --bits 32 --seed 42
```
```
RSA 32-bit Keypair:
  Modulus   (n) = 2369625997
  Public  e     = 65537
  Private d     = 18514673
  Prime p       = 48817
  Prime q       = 48541
  Euler phi     = 2369528640
  Carmichael λ  = 197460720
```

```bash
foxmath rsa encrypt 42 17 3233
```
```
Ciphertext c = 2557
```

```bash
foxmath rsa decrypt 2557 413 3233 --p 61 --q 53
```
```
Decrypted message m = 42 (CRT accelerated)
```

```bash
# Factor modulus n via Pollard's rho, recover d, and decrypt ciphertext
foxmath rsa crack 3233 17 2790
```
```
Factored n = 3233:
  p = 61
  q = 53
  Euler phi    = 3120
  Carmichael λ = 780
  Recovered d  = 413
  Decrypted m  = 65
```

```bash
# Diffie-Hellman: keygen, shared secret derivation, and BSGS discrete log attack
foxmath dh keygen 23 5 --priv 4
```
```
Diffie-Hellman Keypair (p=23, g=5):
  Private key (a) = 4
  Public key  (A) = 4
```

```bash
foxmath dh derive 4 8 23
```
```
Shared secret S = 2
```

```bash
# Recover Alice's private key a from public key A = 4 via BSGS
foxmath dh crack 5 4 23
```
```
Discrete logarithm solved via BSGS:
  Recovered private key a = 4
  Verification: 5^4 ≡ 4 (mod 23)
```

```bash
foxmath factor 360
```
```
360 = 2^3 * 3^2 * 5
```

```bash
foxmath isprime 104729
```
```
104729 is prime
```

```bash
foxmath isprime 561
```
```
561 is composite
```

```bash
foxmath pi --terms 200
```
```
π ≈ 3.1415926535897932384626433832795028841971693993752  (200 terms)
```

```bash
foxmath legendre 5 17
```
```
Legendre (5/17) = -1
```

```bash
foxmath jacobi 10 13
```
```
Jacobi (10/13) = 1
```

```bash
foxmath solovay 104729
```
```
104729 is probably prime (Solovay-Strassen, 20 rounds)
```

```bash
foxmath sqrtmod 10 13
```
```
sqrt(10) = +/-7 (mod 13)
```

`sqrtmod` uses the Tonelli-Shanks algorithm to find modular square roots
when they exist. It reports an error when the input is a quadratic
non-residue modulo the given odd prime.

```bash
foxmath egcd 30 50
```
```
gcd(30, 50) = 10 = 30*(2) + 50*(-1)
```

```bash
foxmath modinv 3 11
```
```
3^-1 ≡ 4 (mod 11)
```

```bash
foxmath totient 360
```
```
phi(360) = 96
lambda(360) = 12
```

```bash
foxmath order 3 7
```
```
ord_7(3) = 6 (primitive root)
```

```bash
foxmath dlog 3 6 7
```
```
3^3 ≡ 6 (mod 7)
```

```bash
foxmath mobius 30
```
```
mu(30) = -1
```

```bash
foxmath sigma 12
```
```
sigma(12) = 28
```

```bash
foxmath tau 12
```
```
tau(12) = 6
```

```bash
# Interactive REPL (all functions preloaded)
foxmath repl
```
```
FoxMath 0.14.0 Interactive REPL 🦊
All math functions loaded. Type exit() or Ctrl-Z + Enter to quit.
>>> keys = rsa_keygen(bits=32, seed=42)
>>> rsa_encrypt(42, keys["e"], keys["n"])
2557
>>> dh_derive(4, 8, 23)
2
>>> jacobi_symbol(10, 13)
1
>>> solovay_strassen(104729)
True
```

```bash
# Symlinks (optional)
foxmath-ecadd --x1 1 --y1 2 --x2 3 --y2 4 --p 17
```
```
Result point: (14, 2)
```

```bash
# Named curve (secp256k1) — computes 2G by default
foxmath ecadd --curve secp256k1
```
```
Result point: (89565891926547004231252920425935692360644145829622209833684329913297188986597, 12158399299693830322967808612713398636155367887041628176798871954788371653930)
  ✓ on curve
```

```bash
# Scalar multiplication — same curves, O(log k) double-and-add
foxmath ecmul --curve toy193 --k 7
```
```
Result point: (66, 24)
  ✓ on curve
```

```bash
# toy193: y² = x³ + x + 5 (mod 193), #E = 193 (prime) — small enough to
# hand-verify, structured exactly like secp256k1. G = (1, 59).
# k = 0 or a multiple of the group order lands on the identity:
foxmath ecmul --curve toy193 --k 193
```
```
Result: point at infinity (identity)
```

```bash
# Chinese Remainder Theorem
foxmath crt --r 2 3 2 --m 3 5 7
```
```
x ≡ 23 (mod 105)
```

```bash
# Continued fraction expansion + convergents
foxmath cf --num 355 --den 113
```
```
[3; 7, 16]
  3/1
  22/7
  355/113
```

```bash
# Challenge/quiz mode — worksheet (no prompting)
foxmath challenge --topic legendre --count 2 --seed 1 --reveal
```
```
1. Legendre symbol (10/13)
   Answer: 1
2. Legendre symbol (3/7)
   Answer: -1
```

```bash
# Baby-size discrete log: given G and Q = k·G, recover k.
# Easy at this scale — the same problem at secp256k1's ~10^77 is why ECDSA works.
foxmath challenge --topic ecdlp --count 1 --seed 9 --reveal
```
```
1. On y² = x³ + 1x + 5 (mod 193) with G = (1, 59): find k (2 ≤ k ≤ 20) such that k·G = (60, 171)
   Answer: 16
```

Drop `--reveal` to answer interactively instead, with scoring at the end.
`--topic` accepts `rsa`, `dh`, `legendre`, `jacobi`, `crt`, `ecdlp`, `sqrtmod`, `factor`, `totient`, `dlog`, or `mixed` (default).

```bash
# JSON output
foxmath pi --terms 100 --json
```
```json
{
  "tool": "foxmath",
  "command": "pi",
  "terms": 100,
  "approx": "3.1415926535897932384626433832795028841971693993752"
}
```

## Why FoxMath?
Because math and cryptography are more fun when you can play with them instantly in the terminal — with the same reliability as [FoxPipe](https://github.com/foxhackerzdevs/FoxPipe) and [fox-vault](https://github.com/foxhackerzdevs/fox-vault).

**Simple. Practical. Reliable.**

## License
MIT © Abhrankan Chakrabarti

---

**Star if you find it useful!** Issues and PRs welcome.
