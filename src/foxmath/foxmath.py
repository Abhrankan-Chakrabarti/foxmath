#!/usr/bin/env python3
"""
FoxMath — Minimalist math & crypto explorer 🦊
One binary. Clean. Educational. Safe.
"""

import argparse
import sys
import json
import random
from collections import Counter
from decimal import Decimal, getcontext
from importlib import metadata
from math import gcd, isqrt
from pathlib import Path

def _package_version() -> str:
    """Installed package version, falling back to 'unknown' when run from a bare checkout."""
    try:
        return metadata.version("foxmath")
    except metadata.PackageNotFoundError:
        return "unknown"

# ------------------- Math Functions -------------------

# First 13 primes. As Miller-Rabin bases they make the test *deterministic*
# for every n below 3,317,044,064,679,887,385,961,981 (~3.3e24); 12 bases
# would not be enough (399165290221 * 798330580441, ~3.2e23, fools them).
_PRIME_TEST_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)

def _is_probable_prime(n: int) -> bool:
    """Miller-Rabin primality test.

    Exact for n < ~3.3e24. Above that it is a strong probable-prime test:
    ordinary composites essentially never pass, but a composite built
    specifically to fool these 13 bases could.
    """
    if n < 2:
        return False
    for sp in _PRIME_TEST_BASES:
        if n == sp:
            return True
        if n % sp == 0:
            return False
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in _PRIME_TEST_BASES:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True

def is_prime(n: int) -> bool:
    """Deterministic Miller-Rabin primality test (exact for n < ~3.3e24)."""
    return _is_probable_prime(n)

def _pollard_brent(n: int, rng: random.Random) -> int:
    """Find a non-trivial factor of composite n using Pollard's rho with Brent's variant."""
    if n % 2 == 0:
        return 2
    if is_prime(n):
        return n

    while True:
        y = rng.randint(1, n - 1)
        c = rng.randint(1, n - 1)
        m = 128
        g, r, q = 1, 1, 1
        ys = y
        while g == 1:
            x = y
            for _ in range(r):
                y = (y * y + c) % n
            k = 0
            while k < r and g == 1:
                ys = y
                limit = min(m, r - k)
                for _ in range(limit):
                    y = (y * y + c) % n
                    q = (q * abs(x - y)) % n
                g = gcd(q, n)
                k += limit
            r *= 2

        if g == n:
            g = 1
            while g == 1:
                ys = (ys * ys + c) % n
                g = gcd(abs(x - ys), n)
            if g == n:
                continue
        return g

def _factor_recursive(n: int, rng: random.Random, factors: list) -> None:
    if n <= 1:
        return
    if is_prime(n):
        factors.append(n)
        return
    d = _pollard_brent(n, rng)
    _factor_recursive(d, rng, factors)
    _factor_recursive(n // d, rng, factors)

def factorize(n: int) -> list:
    """Return sorted list of (prime, power) tuples for integer n >= 1.

    Uses small-prime trial division followed by Pollard's rho (Brent's variant).
    """
    if not isinstance(n, int):
        raise TypeError("n must be an integer")
    if n < 1:
        raise ValueError("n must be a positive integer >= 1")
    if n == 1:
        return []

    factors = []
    # Small prime trial division up to 100
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97):
        while n % p == 0:
            factors.append(p)
            n //= p
        if n == 1:
            break

    if n > 1:
        if is_prime(n):
            factors.append(n)
        else:
            rng = random.Random(n ^ 0x5DEECE66D)
            _factor_recursive(n, rng, factors)

    counts = Counter(factors)
    return sorted(counts.items())

def format_factorization(factors: list) -> str:
    """Format [(p, power), ...] into 'p1^k1 * p2^k2' string."""
    if not factors:
        return "1"
    parts = []
    for p, power in factors:
        if power == 1:
            parts.append(str(p))
        else:
            parts.append(f"{p}^{power}")
    return " * ".join(parts)

def euler_phi(n: int) -> int:
    """Euler's totient function phi(n). Counts positive integers <= n coprime to n."""
    if not isinstance(n, int):
        raise TypeError("n must be an integer")
    if n < 1:
        raise ValueError("n must be a positive integer >= 1")
    if n == 1:
        return 1
    result = 1
    for p, power in factorize(n):
        result *= (p ** (power - 1)) * (p - 1)
    return result

def carmichael_lambda(n: int) -> int:
    """Carmichael function lambda(n). Smallest m such that a^m == 1 (mod n) for all coprime a."""
    if not isinstance(n, int):
        raise TypeError("n must be an integer")
    if n < 1:
        raise ValueError("n must be a positive integer >= 1")
    if n == 1:
        return 1
    lam = 1
    for p, power in factorize(n):
        if p == 2:
            if power == 1:
                pk_lam = 1
            elif power == 2:
                pk_lam = 2
            else:
                pk_lam = 1 << (power - 2)
        else:
            pk_lam = (p ** (power - 1)) * (p - 1)
        lam = (lam * pk_lam) // gcd(lam, pk_lam)
    return lam

def multiplicative_order(a: int, m: int) -> int:
    """Compute the multiplicative order of a modulo m.

    Returns the smallest integer k >= 1 such that a^k == 1 (mod m).
    Raises ValueError if m <= 1 or if gcd(a, m) != 1.
    """
    if not isinstance(a, int) or not isinstance(m, int):
        raise TypeError("a and m must be integers")
    if m <= 1:
        raise ValueError("Modulus m must be greater than 1")
    g = gcd(a, m)
    if g != 1:
        raise ValueError(f"{a} and {m} are not coprime (gcd = {g})")

    a %= m
    if a == 1:
        return 1

    lam = carmichael_lambda(m)
    order = lam
    for p, _ in factorize(lam):
        while order % p == 0 and pow(a, order // p, m) == 1:
            order //= p
    return order

def is_primitive_root(a: int, m: int) -> bool:
    """Return True if a is a primitive root modulo m, False otherwise."""
    if not isinstance(a, int) or not isinstance(m, int):
        raise TypeError("a and m must be integers")
    if m <= 1 or gcd(a, m) != 1:
        return False
    return multiplicative_order(a, m) == euler_phi(m)

def discrete_log(g: int, h: int, p: int) -> int:
    """Solve the discrete logarithm problem g^x ≡ h (mod p) using Baby-step Giant-step.

    Returns the smallest integer x >= 0 such that pow(g, x, p) == h % p.
    Raises ValueError if p <= 1, if gcd(g, p) != 1, or if no solution exists.
    """
    if not isinstance(g, int) or not isinstance(h, int) or not isinstance(p, int):
        raise TypeError("g, h, and p must be integers")
    if p <= 1:
        raise ValueError("Modulus p must be greater than 1")

    g %= p
    h %= p

    if h == 1:
        return 0

    d = gcd(g, p)
    if d != 1:
        raise ValueError(f"Base {g} and modulus {p} must be coprime (gcd = {d})")

    order = multiplicative_order(g, p)
    m = isqrt(order)
    if m * m < order:
        m += 1

    # Baby steps: compute g^j mod p for j in [0, m - 1]
    table = {}
    cur = 1
    for j in range(m):
        if cur not in table:
            table[cur] = j
        cur = (cur * g) % p

    # Giant steps: factor = g^(-m) mod p
    factor = pow(g, -m, p)
    gamma = h
    for i in range((order + m - 1) // m + 1):
        if gamma in table:
            ans = i * m + table[gamma]
            if ans < order:
                return ans
        gamma = (gamma * factor) % p

    raise ValueError(f"No discrete logarithm exists for {g}^x ≡ {h} (mod {p})")

def mobius(n: int) -> int:
    """Möbius function mu(n). Returns (-1)^k for square-free n with k prime factors, 0 otherwise."""
    if not isinstance(n, int):
        raise TypeError("n must be an integer")
    if n < 1:
        raise ValueError("n must be a positive integer >= 1")
    if n == 1:
        return 1
    factors = factorize(n)
    for _, power in factors:
        if power > 1:
            return 0
    return -1 if len(factors) % 2 == 1 else 1

def tau(n: int) -> int:
    """Divisor count function tau(n) (also d(n)). Number of positive divisors of n."""
    if not isinstance(n, int):
        raise TypeError("n must be an integer")
    if n < 1:
        raise ValueError("n must be a positive integer >= 1")
    if n == 1:
        return 1
    total = 1
    for _, power in factorize(n):
        total *= (power + 1)
    return total

def sigma(n: int, k: int = 1) -> int:
    """Divisor sum function sigma_k(n). Sum of k-th powers of positive divisors of n."""
    if not isinstance(n, int) or not isinstance(k, int):
        raise TypeError("n and k must be integers")
    if n < 1:
        raise ValueError("n must be a positive integer >= 1")
    if k < 0:
        raise ValueError("k must be a non-negative integer >= 0")
    if n == 1:
        return 1
    if k == 0:
        return tau(n)
    total = 1
    for p, power in factorize(n):
        pk = p ** k
        term = (pk ** (power + 1) - 1) // (pk - 1)
        total *= term
    return total

def _require_odd_prime(p: int) -> None:
    """Enforce the contract the Legendre symbol and Tonelli-Shanks rely on.

    Euler's criterion is only meaningful for a prime modulus. On a composite
    it returns garbage, and Tonelli-Shanks' search for a non-residue can loop
    forever (mod 9 or mod 561, for example, there is none to find).
    """
    if p < 2 or p % 2 == 0 or not _is_probable_prime(p):
        raise ValueError("p must be an odd prime")

def _legendre_unchecked(a: int, p: int) -> int:
    res = pow(a, (p - 1) // 2, p)
    return 1 if res == 1 else (-1 if res == p - 1 else 0)

def legendre_symbol(a: int, p: int) -> int:
    """Compute Legendre symbol (a/p). p must be an odd prime."""
    _require_odd_prime(p)
    return _legendre_unchecked(a, p)

def jacobi_symbol(a: int, n: int) -> int:
    """Compute the Jacobi symbol (a/n) for integer a and positive odd integer n.

    Generalizes the Legendre symbol to any odd positive integer n via reciprocity
    without factoring n. Returns 1, -1, or 0 (if gcd(a, n) > 1).
    Raises TypeError if a or n is not an integer.
    Raises ValueError if n <= 0 or n is even.
    """
    if not isinstance(a, int) or not isinstance(n, int):
        raise TypeError("a and n must be integers")
    if n <= 0 or n % 2 == 0:
        raise ValueError("n must be a positive odd integer")

    a %= n
    result = 1

    while a != 0:
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5):
                result = -result

        if a % 4 == 3 and n % 4 == 3:
            result = -result

        a, n = n % a, a

    return result if n == 1 else 0

def solovay_strassen(n: int, rounds: int = 20, rng: random.Random = None) -> bool:
    """Probabilistic Solovay-Strassen primality test using Euler's criterion.

    Tests whether n is prime by verifying a^((n-1)/2) ≡ (a/n) (mod n) for random bases a.
    Every composite number has at least n/2 Euler witnesses, giving error probability <= 2^(-rounds).
    Raises TypeError if n or rounds is not an integer.
    Raises ValueError if rounds < 1.
    """
    if not isinstance(n, int):
        raise TypeError("n must be an integer")
    if not isinstance(rounds, int):
        raise TypeError("rounds must be an integer")
    if rounds < 1:
        raise ValueError("rounds must be >= 1")

    if n < 2:
        return False
    if n in (2, 3):
        return True
    if n % 2 == 0:
        return False

    if rng is None:
        rng = random.Random()

    for _ in range(rounds):
        a = rng.randint(2, n - 2) if n > 4 else 2
        x = jacobi_symbol(a, n)
        if x == 0:
            return False
        mod_exp = pow(a, (n - 1) // 2, n)
        if mod_exp != (x % n):
            return False
    return True

def sqrt_mod(a: int, p: int) -> tuple:
    """Return the square root(s) of a modulo an odd prime p."""
    _require_odd_prime(p)

    a %= p
    symbol = _legendre_unchecked(a, p)
    if symbol == 0:
        return (0,)
    if symbol == -1:
        raise ValueError(f"{a} has no square root modulo {p}")

    # For p ≡ 3 (mod 4), the root is a^(p+1)/4.
    if p % 4 == 3:
        root = pow(a, (p + 1) // 4, p)
    else:
        # Tonelli-Shanks: write p - 1 = q * 2^s with q odd.
        q, s = p - 1, 0
        while q % 2 == 0:
            q //= 2
            s += 1

        non_residue = 2
        while _legendre_unchecked(non_residue, p) != -1:
            non_residue += 1

        c = pow(non_residue, q, p)
        root = pow(a, (q + 1) // 2, p)
        t = pow(a, q, p)
        m = s

        while t != 1:
            i = 1
            t_squared = (t * t) % p
            while t_squared != 1:
                t_squared = (t_squared * t_squared) % p
                i += 1
                if i == m:
                    raise ValueError(f"{a} has no square root modulo {p}")
            b = pow(c, 1 << (m - i - 1), p)
            root = (root * b) % p
            t = (t * b * b) % p
            c = (b * b) % p
            m = i

    other = (-root) % p
    return tuple(sorted({root, other}))

def _arctan_series(inv_n: Decimal, terms: int) -> Decimal:
    """arctan(1/n) via its Taylor series, evaluated to `terms` terms. `inv_n` is 1/n."""
    total = Decimal(0)
    for k in range(terms):
        sign = Decimal(-1) ** k
        total += sign * inv_n**(2 * k + 1) / (2 * k + 1)
    return total

def euler_hermann_pi_approx(terms: int = 100) -> Decimal:
    """
    Fast-converging π approximation via the Euler/Hermann Machin-like
    identity: pi/4 = arctan(1/2) + arctan(1/3).
    """
    getcontext().prec = 50
    s = _arctan_series(Decimal(1) / 2, terms) + _arctan_series(Decimal(1) / 3, terms)
    return 4 * s

def ec_point_add(x1, y1, x2, y2, a, p):
    """Simple elliptic curve point addition over finite field y² = x³ + a x + b (mod p)"""
    x1, y1, x2, y2 = x1 % p, y1 % p, x2 % p, y2 % p

    if x1 % p == x2 % p and (y1 + y2) % p == 0:
        raise ValueError(
            "Result is the point at infinity (P + (-P)); this simplified "
            "implementation has no way to represent it."
        )
    if x1 == x2 and y1 == y2:
        # Doubling
        lam = (3 * x1**2 + a) * pow(2 * y1, -1, p) % p
    else:
        lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (lam**2 - x1 - x2) % p
    y3 = (lam * (x1 - x3) - y1) % p
    return x3, y3

# Named curve presets: a, b, p, and the standard generator point G (if defined).
CURVES = {
    "secp256k1": {
        "a": 0,
        "b": 7,
        "p": 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F,
        "gx": 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
        "gy": 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8,
    },
    "toy97": {  # small curve for quick hand-checkable examples: y^2 = x^3 + 2x + 3 (mod 97)
        "a": 2,
        "b": 3,
        "p": 97,
        "gx": 3,
        "gy": 6,
    },
    "toy193": {  # teaching curve: y^2 = x^3 + x + 5 (mod 193); #E = 193 (prime),
                 # so every non-identity point generates the whole group —
                 # hand-checkable discrete-log experiments live here.
        "a": 1,
        "b": 5,
        "p": 193,
        "gx": 1,
        "gy": 59,
    },
}

def is_on_curve(x: int, y: int, a: int, b: int, p: int) -> bool:
    """Check y^2 = x^3 + a*x + b (mod p)."""
    return (y * y - (x ** 3 + a * x + b)) % p == 0

def _ec_add_raw(P, Q, a, p):
    """
    Point addition with an explicit identity: None represents the point at
    infinity. Used by ec_mul; ec_point_add above is the user-facing variant,
    which raises on infinity instead of representing it.
    """
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if P == Q:
        lam = (3 * x1 * x1 + a) * pow(2 * y1, -1, p) % p
    else:
        lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (lam * lam - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)

def ec_mul(k: int, x: int, y: int, a: int, p):
    """
    Scalar multiplication k·P via double-and-add on y² = x³ + a x + b (mod p).
    Returns (x, y), or None for the point at infinity (k ≡ 0 mod ord(P)).

    Double-and-add: O(log k) group operations instead of k — this is exactly
    why the forward direction (compute k·G) is easy while the reverse
    (recover k from k·G, the discrete logarithm) is believed hard for large
    groups. That asymmetry is the engine inside ECDSA and ECDH.
    """
    if k < 0:
        raise ValueError("k must be non-negative")
    result = None  # identity
    addend = (x % p, y % p)
    while k > 0:
        if k & 1:
            result = _ec_add_raw(result, addend, a, p)
        addend = _ec_add_raw(addend, addend, a, p)
        k >>= 1
    return result

def extended_gcd(a: int, b: int) -> tuple:
    """Returns (g, x, y) such that a*x + b*y == g == gcd(a, b), with g >= 0."""
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r != 0:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    if old_r < 0:
        old_r, old_s, old_t = -old_r, -old_s, -old_t
    return old_r, old_s, old_t

_extended_gcd = extended_gcd

def mod_inverse(a: int, m: int) -> int:
    """Compute the modular multiplicative inverse of a modulo m."""
    if m <= 1:
        raise ValueError("Modulus m must be greater than 1")
    try:
        return pow(a, -1, m)
    except ValueError:
        raise ValueError(f"{a} has no modular inverse modulo {m} (gcd({a}, {m}) != 1)")

def _combine_crt(r1: int, m1: int, r2: int, m2: int):
    """Combine x≡r1 (mod m1) and x≡r2 (mod m2) into a single x≡r (mod lcm(m1,m2))."""
    g, p, _ = _extended_gcd(m1, m2)
    if (r2 - r1) % g != 0:
        raise ValueError(
            f"No solution exists: moduli {m1} and {m2} are inconsistent "
            f"for remainders {r1} and {r2}."
        )
    lcm = m1 // g * m2
    x = (r1 + (r2 - r1) // g * p % (m2 // g) * m1) % lcm
    return x, lcm

def crt_solve(remainders, moduli):
    """
    Chinese Remainder Theorem, generalized to non-pairwise-coprime moduli.
    Returns (x, m) such that x is the unique solution mod m = lcm(moduli),
    or raises ValueError if the system is inconsistent.
    """
    if len(remainders) != len(moduli):
        raise ValueError("remainders and moduli must have the same length")
    if not remainders:
        raise ValueError("need at least one congruence")
    if any(m <= 0 for m in moduli):
        raise ValueError("all moduli must be positive")

    r, m = remainders[0] % moduli[0], moduli[0]
    for ri, mi in zip(remainders[1:], moduli[1:]):
        r, m = _combine_crt(r, m, ri % mi, mi)
    return r, m

def continued_fraction(num: int, den: int, max_terms: int = 100):
    """
    Continued fraction expansion [a0; a1, a2, ...] of num/den via the
    Euclidean algorithm. Terminates naturally for any rational input.
    """
    if den == 0:
        raise ValueError("denominator cannot be zero")
    cf = []
    n, d = num, den
    while d != 0 and len(cf) < max_terms:
        a = n // d
        cf.append(a)
        n, d = d, n - a * d
    return cf

def cf_convergents(cf):
    """
    Convergents (p_n, q_n) of a continued fraction, via the standard
    recurrence h_n = a_n*h_(n-1) + h_(n-2), same for k_n, with
    h_(-2)=0, h_(-1)=1, k_(-2)=1, k_(-1)=0.
    """
    convergents = []
    h2, h1 = 0, 1
    k2, k1 = 1, 0
    for a in cf:
        h = a * h1 + h2
        k = a * k1 + k2
        convergents.append((h, k))
        h2, h1 = h1, h
        k2, k1 = k1, k
    return convergents

# ------------------- Challenge / Quiz Mode -------------------

_SMALL_ODD_PRIMES = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
                      53, 59, 61, 67, 71, 73, 79, 83, 89, 97]
_CRT_MODULI_CANDIDATES = [3, 4, 5, 7, 9, 11, 13]

def _generate_legendre_problem(rng: random.Random) -> dict:
    p = rng.choice(_SMALL_ODD_PRIMES)
    a = rng.randint(1, p - 1)
    return {
        "topic": "legendre",
        "question": f"Legendre symbol ({a}/{p})",
        "answer": legendre_symbol(a, p),
    }

def _generate_crt_problem(rng: random.Random) -> dict:
    while True:
        m1, m2 = rng.sample(_CRT_MODULI_CANDIDATES, 2)
        if gcd(m1, m2) == 1:
            break
    r1, r2 = rng.randrange(m1), rng.randrange(m2)
    x, m = crt_solve([r1, r2], [m1, m2])
    return {
        "topic": "crt",
        "question": f"x ≡ {r1} (mod {m1}), x ≡ {r2} (mod {m2}); find x (0 ≤ x < {m})",
        "answer": x,
    }

_ECDLP_CURVE = "toy193"  # prime group order → every point is a generator

def _generate_ecdlp_problem(rng: random.Random) -> dict:
    """
    A baby-size elliptic-curve discrete logarithm: given G and Q = k·G on the
    teaching curve, recover k. With 2 ≤ k ≤ 20 it is solvable by hand through
    repeated addition — which is precisely the lesson: on secp256k1 the same
    search space is ~10^77, and brute force stops being an option.
    """
    c = CURVES[_ECDLP_CURVE]
    k = rng.randint(2, 20)
    qx, qy = ec_mul(k, c["gx"], c["gy"], c["a"], c["p"])
    return {
        "topic": "ecdlp",
        "question": (
            f"On y² = x³ + {c['a']}x + {c['b']} (mod {c['p']}) with G = "
            f"({c['gx']}, {c['gy']}): find k (2 ≤ k ≤ 20) such that k·G = ({qx}, {qy})"
        ),
        "answer": k,
    }

_SQRTMOD_PRIMES = [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]

def _generate_sqrtmod_problem(rng: random.Random) -> dict:
    p = rng.choice(_SQRTMOD_PRIMES)
    r = rng.randint(2, (p - 1) // 2)
    a = (r * r) % p
    return {
        "topic": "sqrtmod",
        "question": f"Find the smaller positive square root of {a} modulo {p} (1 ≤ r ≤ {p // 2})",
        "answer": r,
    }

_FACTOR_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

def _generate_factor_problem(rng: random.Random) -> dict:
    count = rng.choice([2, 2, 2, 3])
    chosen = [rng.choice(_FACTOR_PRIMES) for _ in range(count)]
    n = 1
    for p in chosen:
        n *= p
    min_p = min(chosen)
    return {
        "topic": "factor",
        "question": f"Find the smallest prime factor of {n}",
        "answer": min_p,
    }

_TOTIENT_NUMS = [
    6, 8, 9, 10, 12, 14, 15, 16, 18, 20, 21, 22, 24, 25, 26, 27, 28, 30,
    35, 36, 40, 42, 45, 49, 50, 54, 60, 72, 84, 100,
]

def _generate_totient_problem(rng: random.Random) -> dict:
    n = rng.choice(_TOTIENT_NUMS)
    return {
        "topic": "totient",
        "question": f"Compute Euler's totient phi({n})",
        "answer": euler_phi(n),
    }

_DLOG_PROBLEMS = [
    (2, 13), (3, 17), (2, 19), (2, 29), (3, 31),
    (2, 37), (3, 41), (3, 43), (2, 53), (2, 59),
]

def _generate_dlog_problem(rng: random.Random) -> dict:
    g, p = rng.choice(_DLOG_PROBLEMS)
    order = multiplicative_order(g, p)
    x = rng.randint(1, order - 1)
    h = pow(g, x, p)
    return {
        "topic": "dlog",
        "question": f"Find x (0 ≤ x < {order}) such that {g}^x ≡ {h} (mod {p})",
        "answer": x,
    }

_JACOBI_MODULI = [
    9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 39, 45, 51, 55, 63, 65, 77, 85, 91, 99,
]

def _generate_jacobi_problem(rng: random.Random) -> dict:
    n = rng.choice(_JACOBI_MODULI)
    a = rng.randint(1, n - 1)
    return {
        "topic": "jacobi",
        "question": f"Jacobi symbol ({a}/{n})",
        "answer": jacobi_symbol(a, n),
    }

_TOPIC_GENERATORS = {
    "legendre": _generate_legendre_problem,
    "crt": _generate_crt_problem,
    "ecdlp": _generate_ecdlp_problem,
    "sqrtmod": _generate_sqrtmod_problem,
    "factor": _generate_factor_problem,
    "totient": _generate_totient_problem,
    "dlog": _generate_dlog_problem,
    "jacobi": _generate_jacobi_problem,
}

def generate_problem(topic: str, rng: random.Random) -> dict:
    if topic == "mixed":
        topic = rng.choice(list(_TOPIC_GENERATORS.keys()))
    if topic not in _TOPIC_GENERATORS:
        raise ValueError(f"Unknown challenge topic: {topic}")
    return _TOPIC_GENERATORS[topic](rng)

# Convenience aliases
dlog = discrete_log
egcd = extended_gcd
modinv = mod_inverse
crt = crt_solve

def cmd_repl() -> None:
    """Launch interactive FoxMath REPL."""
    import code

    banner = (
        f"FoxMath {_package_version()} Interactive REPL 🦊\n"
        "All math functions loaded. Type exit() or Ctrl-Z + Enter to quit.\n"
        "Functions: jacobi_symbol, solovay_strassen, discrete_log (dlog),\n"
        "           mobius, sigma, tau, factorize, is_prime, euler_phi, ..."
    )
    print(banner)
    namespace = {
        "__name__": "__console__",
        "__doc__": None,
        "is_prime": is_prime,
        "factorize": factorize,
        "format_factorization": format_factorization,
        "legendre_symbol": legendre_symbol,
        "jacobi_symbol": jacobi_symbol,
        "solovay_strassen": solovay_strassen,
        "sqrt_mod": sqrt_mod,
        "discrete_log": discrete_log,
        "dlog": discrete_log,
        "euler_phi": euler_phi,
        "carmichael_lambda": carmichael_lambda,
        "multiplicative_order": multiplicative_order,
        "is_primitive_root": is_primitive_root,
        "mobius": mobius,
        "tau": tau,
        "sigma": sigma,
        "extended_gcd": extended_gcd,
        "egcd": extended_gcd,
        "mod_inverse": mod_inverse,
        "modinv": mod_inverse,
        "crt_solve": crt_solve,
        "crt": crt_solve,
        "euler_hermann_pi_approx": euler_hermann_pi_approx,
        "ec_point_add": ec_point_add,
        "ec_mul": ec_mul,
        "is_on_curve": is_on_curve,
        "CURVES": CURVES,
        "continued_fraction": continued_fraction,
        "cf_convergents": cf_convergents,
    }
    code.interact(banner="", local=namespace, exitmsg="")

# ------------------- CLI -------------------

def get_command():
    """Support symlinks: foxmath-legendre, foxmath-pi, etc."""
    name = Path(sys.argv[0]).stem.lower()
    if name.startswith("foxmath-"):
        return name.split("-", 1)[1]
    return None

def main():
    if hasattr(sys.stdout, "reconfigure"):
        try:
            sys.stdout.reconfigure(encoding="utf-8")
        except Exception:
            pass
    if hasattr(sys.stderr, "reconfigure"):
        try:
            sys.stderr.reconfigure(encoding="utf-8")
        except Exception:
            pass

    cmd = get_command()

    # Handle --json independently of argparse subparser scoping: a subparser's
    # own default for a same-named flag silently overwrites whatever the
    # parent parser already parsed for it, so --json only ever worked in one
    # position when defined via argparse directly. Strip it out manually
    # instead, so it works before or after the subcommand.
    argv = sys.argv[1:]
    json_output = "--json" in argv
    argv = [a for a in argv if a != "--json"]

    if cmd:
        # Symlink invocation (e.g. `foxmath-legendre 2 7`): argparse's
        # subparsers positional still expects a command-name token next,
        # so without this it tries to match the first real argument
        # (e.g. "2") against {legendre,pi,ecadd} and fails. Inject the
        # detected command so it parses exactly as if typed explicitly.
        argv = [cmd] + argv

    parser = argparse.ArgumentParser(description="FoxMath — Math & Crypto Explorer 🦊")
    parser.add_argument("--version", action="version", version=f"foxmath {_package_version()}")
    subparsers = parser.add_subparsers(dest="command", required=not cmd)

    # Legendre
    leg = subparsers.add_parser("legendre", help="Legendre symbol (a/p)")
    leg.add_argument("a", type=int, help="Integer a")
    leg.add_argument("p", type=int, help="Odd prime p")

    sqrt = subparsers.add_parser("sqrtmod", help="Square root(s) modulo an odd prime")
    sqrt.add_argument("a", type=int, help="Integer a")
    sqrt.add_argument("p", type=int, help="Odd prime p")

    prime_cmd = subparsers.add_parser("isprime", help="Deterministic Miller-Rabin primality test")
    prime_cmd.add_argument("n", type=int, help="Integer to test")

    factor_cmd = subparsers.add_parser("factor", help="Prime factorization into prime powers")
    factor_cmd.add_argument("n", type=int, help="Positive integer to factor (n >= 1)")

    # Pi
    pi_cmd = subparsers.add_parser("pi", help="Euler/Hermann Machin-like π approximation")
    pi_cmd.add_argument("--terms", type=int, default=100)

    # Elliptic curve add
    ec = subparsers.add_parser("ecadd", help="Elliptic curve point addition")
    ec.add_argument("--curve", choices=sorted(CURVES.keys()),
                     help="Named curve (provides a, b, p, and generator G)")
    ec.add_argument("--x1", type=int, help="Defaults to G's x if --curve is given")
    ec.add_argument("--y1", type=int, help="Defaults to G's y if --curve is given")
    ec.add_argument("--x2", type=int, help="Defaults to (x1,y1) if omitted (i.e. doubling)")
    ec.add_argument("--y2", type=int, help="Defaults to (x1,y1) if omitted (i.e. doubling)")
    ec.add_argument("--a", type=int, help="Curve parameter a (overrides --curve if given)")
    ec.add_argument("--p", type=int, help="Prime modulus (overrides --curve if given)")

    # Elliptic curve scalar multiplication
    ecm = subparsers.add_parser("ecmul", help="Elliptic curve scalar multiplication (k·P)")
    ecm.add_argument("--curve", choices=sorted(CURVES.keys()),
                     help="Named curve (provides a, b, p, and generator G)")
    ecm.add_argument("--k", type=int, required=True, help="Scalar multiplier")
    ecm.add_argument("--x", type=int, help="Point x (defaults to G's x if --curve is given)")
    ecm.add_argument("--y", type=int, help="Point y (defaults to G's y if --curve is given)")
    ecm.add_argument("--a", type=int, help="Curve parameter a (overrides --curve if given)")
    ecm.add_argument("--p", type=int, help="Prime modulus (overrides --curve if given)")

    # Chinese Remainder Theorem
    crt = subparsers.add_parser("crt", help="Solve a system of congruences x = r (mod m)")
    crt.add_argument("--r", type=int, nargs="+", required=True, help="Remainders, e.g. --r 2 3 2")
    crt.add_argument("--m", type=int, nargs="+", required=True, help="Moduli, e.g. --m 3 5 7")

    # Continued fractions
    cf_cmd = subparsers.add_parser("cf", help="Continued fraction expansion and convergents")
    cf_cmd.add_argument("--num", type=int, required=True, help="Numerator")
    cf_cmd.add_argument("--den", type=int, required=True, help="Denominator")
    cf_cmd.add_argument("--terms", type=int, default=100, help="Max terms to expand")

    # Extended GCD
    egcd_cmd = subparsers.add_parser("egcd", help="Extended Euclidean algorithm (Bézout's identity)")
    egcd_cmd.add_argument("a", type=int, help="Integer a")
    egcd_cmd.add_argument("b", type=int, help="Integer b")

    # Modular Inverse
    minv_cmd = subparsers.add_parser("modinv", help="Modular multiplicative inverse a^-1 (mod m)")
    minv_cmd.add_argument("a", type=int, help="Integer a")
    minv_cmd.add_argument("m", type=int, help="Modulus m (> 1)")

    # Totient
    totient_cmd = subparsers.add_parser("totient", help="Euler's totient phi(n) and Carmichael lambda(n)")
    totient_cmd.add_argument("n", type=int, help="Positive integer to evaluate (n >= 1)")

    # Multiplicative Order
    order_cmd = subparsers.add_parser("order", help="Multiplicative order of a modulo m")
    order_cmd.add_argument("a", type=int, help="Base integer a")
    order_cmd.add_argument("m", type=int, help="Modulus m (> 1)")

    # Discrete Logarithm
    dlog_cmd = subparsers.add_parser("dlog", help="Discrete logarithm solver: g^x ≡ h (mod p)")
    dlog_cmd.add_argument("g", type=int, help="Base integer g")
    dlog_cmd.add_argument("h", type=int, help="Target value h")
    dlog_cmd.add_argument("p", type=int, help="Modulus p (> 1)")

    # Möbius function
    mob_cmd = subparsers.add_parser("mobius", help="Möbius function mu(n)")
    mob_cmd.add_argument("n", type=int, help="Positive integer to evaluate (n >= 1)")

    # Sum of divisors
    sig_cmd = subparsers.add_parser("sigma", help="Sum of divisors function sigma_k(n)")
    sig_cmd.add_argument("n", type=int, help="Positive integer to evaluate (n >= 1)")
    sig_cmd.add_argument("--k", type=int, default=1, help="Exponent k (default: 1)")

    # Divisor count
    tau_cmd = subparsers.add_parser("tau", help="Divisor count function tau(n)")
    tau_cmd.add_argument("n", type=int, help="Positive integer to evaluate (n >= 1)")

    # Jacobi
    jac_cmd = subparsers.add_parser("jacobi", help="Jacobi symbol (a/n) for positive odd n")
    jac_cmd.add_argument("a", type=int, help="Integer a")
    jac_cmd.add_argument("n", type=int, help="Positive odd integer n")

    # Solovay-Strassen primality test
    sol_cmd = subparsers.add_parser("solovay", help="Solovay-Strassen probabilistic primality test")
    sol_cmd.add_argument("n", type=int, help="Integer to test")
    sol_cmd.add_argument("--rounds", type=int, default=20, help="Number of test rounds (default: 20)")

    # Interactive REPL
    subparsers.add_parser("repl", help="Interactive FoxMath REPL")

    # Challenge / quiz mode
    challenge = subparsers.add_parser("challenge", help="Practice problems (Legendre, Jacobi, CRT, SqrtMod, Factor, Totient, DLog)")
    challenge.add_argument("--topic", choices=["legendre", "crt", "ecdlp", "sqrtmod", "factor", "totient", "dlog", "jacobi", "mixed"], default="mixed")
    challenge.add_argument("--count", type=int, default=5, help="Number of problems")
    challenge.add_argument("--seed", type=int, help="Random seed for reproducible problem sets")
    challenge.add_argument("--reveal", action="store_true",
                            help="Print problems with answers, no prompting (worksheet mode)")

    args = parser.parse_args(argv)

    result = {"tool": "foxmath"}

    try:
        if cmd == "legendre" or args.command == "legendre":
            ls = legendre_symbol(args.a, args.p)
            result.update({"command": "legendre", "a": args.a, "p": args.p, "value": ls})
            if not json_output:
                print(f"Legendre ({args.a}/{args.p}) = {ls}")

        elif cmd == "sqrtmod" or args.command == "sqrtmod":
            roots = sqrt_mod(args.a, args.p)
            result.update({"command": "sqrtmod", "a": args.a, "p": args.p, "roots": list(roots)})
            if not json_output:
                if len(roots) == 1:
                    print(f"sqrt({args.a}) = {roots[0]} (mod {args.p})")
                else:
                    print(f"sqrt({args.a}) = +/-{roots[1]} (mod {args.p})")

        elif cmd == "isprime" or args.command == "isprime":
            prime = is_prime(args.n)
            result.update({"command": "isprime", "n": args.n, "is_prime": prime})
            if not json_output:
                if args.n < 2:
                    print(f"{args.n} is not prime")
                else:
                    verdict = "prime" if prime else "composite"
                    print(f"{args.n} is {verdict}")

        elif cmd == "factor" or args.command == "factor":
            factors = factorize(args.n)
            is_p = is_prime(args.n)
            expr = format_factorization(factors)
            prime_factors = [p for p, _ in factors]
            factors_dict = [{"prime": p, "power": power} for p, power in factors]
            result.update({
                "command": "factor",
                "n": args.n,
                "factors": factors_dict,
                "prime_factors": prime_factors,
                "expression": expr,
                "is_prime": is_p,
            })
            if not json_output:
                if is_p:
                    print(f"{args.n} is prime")
                elif args.n == 1:
                    print("1 has no prime factors")
                else:
                    print(f"{args.n} = {expr}")

        elif cmd == "pi" or args.command == "pi":
            approx = euler_hermann_pi_approx(args.terms)
            result.update({"command": "pi", "terms": args.terms, "approx": str(approx)})
            if not json_output:
                print(f"π ≈ {approx}  ({args.terms} terms)")

        elif cmd == "ecadd" or args.command == "ecadd":
            a, b = args.a, None
            p = args.p
            x1, y1 = args.x1, args.y1

            if args.curve:
                preset = CURVES[args.curve]
                a = preset["a"] if a is None else a
                p = preset["p"] if p is None else p
                b = preset["b"]
                if x1 is None:
                    x1, y1 = preset["gx"], preset["gy"]
            elif a is None:
                a = -3  # historical default for manual (non-curve) mode

            if a is None or p is None or x1 is None or y1 is None:
                raise ValueError(
                    "Provide --curve NAME, or at least --p --x1 --y1 manually "
                    "(--a defaults to -3 if omitted)."
                )

            x2 = x1 if args.x2 is None else args.x2
            y2 = y1 if args.y2 is None else args.y2

            x3, y3 = ec_point_add(x1, y1, x2, y2, a, p)
            result.update({
                "command": "ecadd", "curve": args.curve,
                "p1": [x1, y1], "p2": [x2, y2], "result": (x3, y3),
            })
            if b is not None:
                result["result_on_curve"] = is_on_curve(x3, y3, a, b, p)
            if not json_output:
                print(f"Result point: ({x3}, {y3})")
                if b is not None:
                    status = "✓ on curve" if result["result_on_curve"] else "✗ NOT on curve"
                    print(f"  {status}")

        elif cmd == "ecmul" or args.command == "ecmul":
            a, b = args.a, None
            p = args.p
            x1, y1 = args.x, args.y

            if args.curve:
                preset = CURVES[args.curve]
                a = preset["a"] if a is None else a
                p = preset["p"] if p is None else p
                b = preset["b"]
                if x1 is None or y1 is None:
                    x1, y1 = preset["gx"], preset["gy"]
            elif a is None:
                a = -3  # historical default for manual (non-curve) mode

            if a is None or p is None or x1 is None or y1 is None:
                raise ValueError(
                    "Provide --curve NAME, or at least --p --x --y manually "
                    "(--a defaults to -3 if omitted)."
                )

            R = ec_mul(args.k, x1, y1, a, p)
            result.update({
                "command": "ecmul", "curve": args.curve, "k": args.k,
                "point": [x1, y1], "result": "infinity" if R is None else list(R),
            })
            if not json_output:
                if R is None:
                    print("Result: point at infinity (identity)")
                else:
                    print(f"Result point: ({R[0]}, {R[1]})")
            if R is not None and b is not None:
                result["result_on_curve"] = is_on_curve(R[0], R[1], a, b, p)
                if not json_output:
                    status = "✓ on curve" if result["result_on_curve"] else "✗ NOT on curve"
                    print(f"  {status}")

        elif cmd == "crt" or args.command == "crt":
            x, m = crt_solve(args.r, args.m)
            result.update({"command": "crt", "r": args.r, "m": args.m, "x": x, "mod": m})
            if not json_output:
                print(f"x ≡ {x} (mod {m})")

        elif cmd == "cf" or args.command == "cf":
            cf = continued_fraction(args.num, args.den, args.terms)
            convs = cf_convergents(cf)
            result.update({
                "command": "cf", "num": args.num, "den": args.den,
                "cf": cf, "convergents": convs,
            })
            if not json_output:
                print(f"[{cf[0]}; {', '.join(str(a) for a in cf[1:])}]")
                for p, q in convs:
                    print(f"  {p}/{q}")

        elif cmd == "egcd" or args.command == "egcd":
            g, x, y = extended_gcd(args.a, args.b)
            result.update({"command": "egcd", "a": args.a, "b": args.b, "gcd": g, "x": x, "y": y})
            if not json_output:
                print(f"gcd({args.a}, {args.b}) = {g} = {args.a}*({x}) + {args.b}*({y})")

        elif cmd == "modinv" or args.command == "modinv":
            inv = mod_inverse(args.a, args.m)
            result.update({"command": "modinv", "a": args.a, "m": args.m, "inverse": inv})
            if not json_output:
                print(f"{args.a}^-1 ≡ {inv} (mod {args.m})")

        elif cmd == "totient" or args.command == "totient":
            phi = euler_phi(args.n)
            lam = carmichael_lambda(args.n)
            cyclic = (phi == lam)
            result.update({
                "command": "totient",
                "n": args.n,
                "phi": phi,
                "lambda": lam,
                "is_cyclic": cyclic,
            })
            if not json_output:
                cyclic_str = " (cyclic)" if cyclic else ""
                print(f"phi({args.n}) = {phi}")
                print(f"lambda({args.n}) = {lam}{cyclic_str}")

        elif cmd == "order" or args.command == "order":
            ord_k = multiplicative_order(args.a, args.m)
            phi = euler_phi(args.m)
            is_prim = (ord_k == phi)
            result.update({
                "command": "order",
                "a": args.a,
                "m": args.m,
                "order": ord_k,
                "is_primitive_root": is_prim,
            })
            if not json_output:
                prim_str = " (primitive root)" if is_prim else ""
                print(f"ord_{args.m}({args.a}) = {ord_k}{prim_str}")

        elif cmd == "dlog" or args.command == "dlog":
            x = discrete_log(args.g, args.h, args.p)
            ord_g = multiplicative_order(args.g, args.p)
            result.update({
                "command": "dlog",
                "g": args.g,
                "h": args.h,
                "p": args.p,
                "x": x,
                "order": ord_g,
            })
            if not json_output:
                print(f"{args.g}^{x} ≡ {args.h} (mod {args.p})")

        elif cmd == "mobius" or args.command == "mobius":
            val = mobius(args.n)
            result.update({"command": "mobius", "n": args.n, "value": val})
            if not json_output:
                print(f"mu({args.n}) = {val}")

        elif cmd == "sigma" or args.command == "sigma":
            val = sigma(args.n, args.k)
            result.update({"command": "sigma", "n": args.n, "k": args.k, "value": val})
            if not json_output:
                prefix = f"sigma_{args.k}" if args.k != 1 else "sigma"
                print(f"{prefix}({args.n}) = {val}")

        elif cmd == "tau" or args.command == "tau":
            val = tau(args.n)
            result.update({"command": "tau", "n": args.n, "value": val})
            if not json_output:
                print(f"tau({args.n}) = {val}")

        elif cmd == "jacobi" or args.command == "jacobi":
            js = jacobi_symbol(args.a, args.n)
            result.update({"command": "jacobi", "a": args.a, "n": args.n, "value": js})
            if not json_output:
                print(f"Jacobi ({args.a}/{args.n}) = {js}")

        elif cmd == "solovay" or args.command == "solovay":
            prime = solovay_strassen(args.n, rounds=args.rounds)
            result.update({"command": "solovay", "n": args.n, "rounds": args.rounds, "is_prime": prime})
            if not json_output:
                if args.n < 2:
                    print(f"{args.n} is not prime")
                elif prime:
                    if args.n in (2, 3):
                        print(f"{args.n} is prime")
                    else:
                        print(f"{args.n} is probably prime (Solovay-Strassen, {args.rounds} rounds)")
                else:
                    print(f"{args.n} is composite")

        elif cmd == "repl" or args.command == "repl":
            cmd_repl()
            return

        elif cmd == "challenge" or args.command == "challenge":
            if args.count < 1:
                raise ValueError(f"--count must be >= 1 (got {args.count})")

            rng = random.Random(args.seed)
            problems = [generate_problem(args.topic, rng) for _ in range(args.count)]

            if args.reveal:
                result.update({"command": "challenge", "topic": args.topic, "problems": problems})
                if not json_output:
                    for i, prob in enumerate(problems, 1):
                        print(f"{i}. {prob['question']}")
                        print(f"   Answer: {prob['answer']}")
            else:
                score = 0
                details = []
                for i, prob in enumerate(problems, 1):
                    if not json_output:
                        print(f"{i}. {prob['question']}")
                    try:
                        raw = input("   Your answer: ")
                    except EOFError:
                        raw = ""
                    try:
                        user_answer = int(raw.strip())
                    except ValueError:
                        user_answer = None
                    correct = user_answer == prob["answer"]
                    score += correct
                    details.append({**prob, "your_answer": user_answer, "correct": correct})
                    if not json_output:
                        if correct:
                            print("   ✓ Correct!")
                        else:
                            print(f"   ✗ Incorrect (answer: {prob['answer']})")
                result.update({
                    "command": "challenge", "topic": args.topic,
                    "score": score, "total": args.count, "details": details,
                })
                if not json_output:
                    print(f"\nScore: {score}/{args.count}")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    if json_output:
        print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
