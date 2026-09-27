# FoxMath 🦊

**Minimalist CLI for math & cryptography exploration.**

**More about me / other projects:** [abhrankan.netlify.app](https://abhrankan.netlify.app)

One binary. Clean output. Zero footguns. Educational by design.

Built with curiosity — for students, tinkerers, and crypto/math enthusiasts.

## Features
- Legendre symbol
- Euler/Hermann Machin-like π approximation
- Elliptic curve point addition (over finite fields), including
  named curves (secp256k1) with on-curve result validation
- Elliptic curve scalar multiplication k·P (double-and-add), including
  the point at infinity as a first-class result
- Chinese Remainder Theorem solver (handles non-coprime moduli)
- Continued fraction expansion and convergents
- Challenge/quiz mode — practice problems with scoring, or a
  worksheet mode for self-study
- Symlink-friendly (foxmath-legendre, foxmath-pi, foxmath-ecadd,
  foxmath-ecmul, foxmath-crt, foxmath-cf, foxmath-challenge)
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
foxmath 0.6.2
```

## Usage Examples

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
`--topic` accepts `legendre`, `crt`, `ecdlp`, or `mixed` (default).

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
