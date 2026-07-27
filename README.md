# FoxMath 🦊

**Minimalist CLI for math & cryptography exploration.**

One binary. Clean output. Zero footguns. Educational by design.

Built with curiosity — for students, tinkerers, and crypto/math enthusiasts.

## Features
- Legendre symbol
- Euler/Hermann Machin-like π approximation
- Elliptic curve point addition (over finite fields)
- Chinese Remainder Theorem solver (handles non-coprime moduli)
- Symlink-friendly (foxmath-legendre, foxmath-pi, foxmath-ecadd, foxmath-crt)
- JSON output mode
- Coming soon: continued fractions, challenge mode, more curves

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
# Chinese Remainder Theorem
foxmath crt --r 2 3 2 --m 3 5 7
```
```
x ≡ 23 (mod 105)
```

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
