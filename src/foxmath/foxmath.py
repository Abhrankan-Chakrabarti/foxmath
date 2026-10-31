#!/usr/bin/env python3
"""
FoxMath — Minimalist math & crypto explorer 🦊
One binary. Clean. Educational. Safe.
"""

import argparse
import sys
import json
import random
import time
from collections import Counter
from decimal import Decimal, getcontext
from importlib import metadata
import math
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

def to_sci(n: int, prec: int = 30) -> str:
    """Convert integer n to scientific notation string if it exceeds the given precision.

    Example: to_sci(123456789, prec=6) -> "1.234568e+8"
    """
    if not isinstance(n, int):
        raise TypeError("n must be an integer")
    if prec < 1:
        raise ValueError("prec must be a positive integer >= 1")

    s = str(n)
    if len(s) <= prec:
        return s

    neg = n < 0
    if neg:
        s = str(abs(n))

    exp = len(s) - 1
    sig_len = 1 + prec  # 1 digit before decimal + `prec` digits after

    extended = s if len(s) > sig_len else s + ("0" * (sig_len + 1 - len(s)))
    sig = extended[:sig_len]
    next_digit = extended[sig_len] if len(extended) > sig_len else "0"

    mant_int = int(sig)

    # Round if needed
    if next_digit >= "5":
        mant_int += 1
        if mant_int >= 10 ** sig_len:
            mant_int //= 10
            exp += 1

    mant_str = str(mant_int).rjust(sig_len, "0")
    mantissa = mant_str[0] + "." + mant_str[1:]
    exp_str = f"{exp:+d}"
    sign = "-" if neg else ""
    return f"{sign}{mantissa}e{exp_str}"

def _is_square(d: int) -> bool:
    """Exact integer square test."""
    if d < 0:
        return False
    # Modulo 16 quadratic residues: 0, 1, 4, 9
    if not (0x0213 & (1 << (d & 15))):
        return False
    s = isqrt(d)
    return s * s == d

def fermat_factor(
    n: int,
    quiet: bool = True,
    prec: int = 30,
    print_interval: int = 1000000,
    start_iter: int = 0,
) -> tuple:
    """Factor an integer n using Fermat's difference of squares method (n = x^2 - y^2 = (x - y)(x + y)).

    Returns (p, q, iterations) where p * q == n and p <= q, or None if no factor was found.
    """
    if isinstance(n, str):
        data = n.strip()
        if data.startswith(("0x", "0X")):
            n = int(data, 16)
        else:
            n = int(data)

    if not isinstance(n, int):
        raise TypeError("n must be an integer")
    if n < 1:
        raise ValueError("n must be a positive integer >= 1")
    if start_iter < 0:
        raise ValueError("starting iteration (start_iter) must be non-negative")
    if prec < 1:
        raise ValueError("prec must be a positive integer >= 1")
    if print_interval < 1:
        raise ValueError("print_interval must be a positive integer >= 1")

    if n == 1:
        return (1, 1, 0)

    y0 = isqrt(n)
    if y0 * y0 < n:
        y0 += 1

    y0 += start_iter
    D = y0 * y0 - n
    k = start_iter

    while y0 + (k - start_iter) < n:
        if _is_square(D):
            x = isqrt(D)
            y = y0 + (k - start_iter)
            p = y - x
            q = y + x
            if p * q == n:
                return (p, q, k)

        if not quiet and k % print_interval == 0 and k > 0:
            y = y0 + (k - start_iter)
            approx_x = isqrt(D)
            approx_p = y - approx_x
            approx_q = y + approx_x
            print(
                f"Iteration: {k} p ≈ {to_sci(approx_p, prec)} q ≈ {to_sci(approx_q, prec)}",
                end="\r",
                flush=True,
            )

        k += 1
        D += 2 * y0 + 2 * (k - start_iter) - 1

    return None

def predict_factor(
    n: int,
    quiet: bool = False,
    prec: int = 30,
    print_interval: int = 1000000,
    start_iter: int = 0,
) -> tuple:
    """Alias for fermat_factor matching diffsquare-py naming and defaults."""
    return fermat_factor(
        n,
        quiet=quiet,
        prec=prec,
        print_interval=print_interval,
        start_iter=start_iter,
    )

def pollard_pm1(
    n: int,
    b1: int = 1000,
    b2: int | None = None,
    a: int = 2,
) -> dict:
    """Factor an integer n using Pollard's p - 1 algorithm with two-stage smoothness bounds.

    Returns a dict with factor, cofactor, stage, smoothness information, and steps.
    """
    if isinstance(n, str):
        data = n.strip()
        if data.startswith(("0x", "0X")):
            n = int(data, 16)
        else:
            n = int(data)

    if not isinstance(n, int):
        raise TypeError("n must be an integer")
    if n < 2:
        raise ValueError("n must be an integer >= 2")

    if not isinstance(b1, int) or b1 < 2:
        raise ValueError("b1 must be an integer >= 2")

    if b2 is None:
        b2 = b1 * 10
    if not isinstance(b2, int) or b2 < b1:
        raise ValueError("b2 must be an integer >= b1")

    if not isinstance(a, int) or a <= 1 or a >= n:
        raise ValueError(f"base a must be an integer with 1 < a < {n}")

    start_time = time.perf_counter()

    if n % 2 == 0:
        p = 2
        q = n // 2
        p_factors = factorize(p - 1)
        elapsed_us = int((time.perf_counter() - start_time) * 1_000_000)
        return {
            "n": n,
            "factor": p,
            "cofactor": q,
            "success": True,
            "stage": 0,
            "b1": b1,
            "b2": b2,
            "base": a,
            "p_minus_1": p - 1,
            "p_minus_1_factors": p_factors,
            "p_minus_1_expression": format_factorization(p_factors),
            "steps": 0,
            "message": "n is even; factor 2 found immediately",
            "elapsed_us": elapsed_us,
        }

    if is_prime(n):
        elapsed_us = int((time.perf_counter() - start_time) * 1_000_000)
        return {
            "n": n,
            "factor": None,
            "cofactor": None,
            "success": False,
            "stage": None,
            "b1": b1,
            "b2": b2,
            "base": a,
            "p_minus_1": None,
            "p_minus_1_factors": None,
            "p_minus_1_expression": None,
            "steps": 0,
            "message": f"n = {n} is already prime",
            "elapsed_us": elapsed_us,
        }

    g0 = gcd(a, n)
    if 1 < g0 < n:
        p = min(g0, n // g0)
        q = n // p
        p_factors = factorize(p - 1)
        elapsed_us = int((time.perf_counter() - start_time) * 1_000_000)
        return {
            "n": n,
            "factor": p,
            "cofactor": q,
            "success": True,
            "stage": 0,
            "b1": b1,
            "b2": b2,
            "base": a,
            "p_minus_1": p - 1,
            "p_minus_1_factors": p_factors,
            "p_minus_1_expression": format_factorization(p_factors),
            "steps": 0,
            "message": f"Base {a} shares non-trivial factor {p} with n",
            "elapsed_us": elapsed_us,
        }

    def primes_up_to(limit: int) -> list:
        s = [True] * (limit + 1)
        s[0] = s[1] = False
        for p in range(2, int(limit**0.5) + 1):
            if s[p]:
                for i in range(p * p, limit + 1, p):
                    s[i] = False
        return [i for i, v in enumerate(s) if v]

    all_primes = primes_up_to(b2)
    s1_primes = [p for p in all_primes if p <= b1]
    s2_primes = [p for p in all_primes if p > b1]

    x = a
    steps = 0

    # Stage 1
    for q in s1_primes:
        steps += 1
        pk = q
        k = 1
        while pk * q <= b1:
            pk *= q
            k += 1

        prev_x = x
        x = pow(x, pk, n)
        g = gcd(x - 1, n)

        if 1 < g < n:
            p = min(g, n // g)
            q_factor = n // p
            p_factors = factorize(p - 1)
            elapsed_us = int((time.perf_counter() - start_time) * 1_000_000)
            return {
                "n": n,
                "factor": p,
                "cofactor": q_factor,
                "success": True,
                "stage": 1,
                "b1": b1,
                "b2": b2,
                "base": a,
                "p_minus_1": p - 1,
                "p_minus_1_factors": p_factors,
                "p_minus_1_expression": format_factorization(p_factors),
                "steps": steps,
                "message": f"Found non-trivial factor p = {p} in Stage 1 at prime {q}",
                "elapsed_us": elapsed_us,
            }

        if g == n:
            curr = prev_x
            for _ in range(k):
                curr = pow(curr, q, n)
                g_sub = gcd(curr - 1, n)
                if 1 < g_sub < n:
                    p = min(g_sub, n // g_sub)
                    q_factor = n // p
                    p_factors = factorize(p - 1)
                    elapsed_us = int((time.perf_counter() - start_time) * 1_000_000)
                    return {
                        "n": n,
                        "factor": p,
                        "cofactor": q_factor,
                        "success": True,
                        "stage": 1,
                        "b1": b1,
                        "b2": b2,
                        "base": a,
                        "p_minus_1": p - 1,
                        "p_minus_1_factors": p_factors,
                        "p_minus_1_expression": format_factorization(p_factors),
                        "steps": steps,
                        "message": f"Found non-trivial factor p = {p} in Stage 1 backtrack at prime {q}",
                        "elapsed_us": elapsed_us,
                    }
                if g_sub == n:
                    break
            elapsed_us = int((time.perf_counter() - start_time) * 1_000_000)
            return {
                "n": n,
                "factor": None,
                "cofactor": None,
                "success": False,
                "stage": 1,
                "b1": b1,
                "b2": b2,
                "base": a,
                "p_minus_1": None,
                "p_minus_1_factors": None,
                "p_minus_1_expression": None,
                "steps": steps,
                "message": f"Stage 1 failed: all factors congruent to 1 mod n at prime {q}. Try a smaller B1 or different base.",
                "elapsed_us": elapsed_us,
            }

    # Stage 2
    if s2_primes:
        a0 = x
        q1 = s2_primes[0]
        cur_x = pow(a0, q1, n)
        steps += 1
        g = gcd(cur_x - 1, n)
        if 1 < g < n:
            p = min(g, n // g)
            q_factor = n // p
            p_factors = factorize(p - 1)
            elapsed_us = int((time.perf_counter() - start_time) * 1_000_000)
            return {
                "n": n,
                "factor": p,
                "cofactor": q_factor,
                "success": True,
                "stage": 2,
                "b1": b1,
                "b2": b2,
                "base": a,
                "p_minus_1": p - 1,
                "p_minus_1_factors": p_factors,
                "p_minus_1_expression": format_factorization(p_factors),
                "steps": steps,
                "message": f"Found non-trivial factor p = {p} in Stage 2 at prime {q1}",
                "elapsed_us": elapsed_us,
            }

        a0_d = {}
        for i in range(len(s2_primes) - 1):
            d = s2_primes[i + 1] - s2_primes[i]
            if d not in a0_d:
                a0_d[d] = pow(a0, d, n)
            cur_x = (cur_x * a0_d[d]) % n
            steps += 1
            g = gcd(cur_x - 1, n)
            if 1 < g < n:
                p = min(g, n // g)
                q_factor = n // p
                p_factors = factorize(p - 1)
                elapsed_us = int((time.perf_counter() - start_time) * 1_000_000)
                return {
                    "n": n,
                    "factor": p,
                    "cofactor": q_factor,
                    "success": True,
                    "stage": 2,
                    "b1": b1,
                    "b2": b2,
                    "base": a,
                    "p_minus_1": p - 1,
                    "p_minus_1_factors": p_factors,
                    "p_minus_1_expression": format_factorization(p_factors),
                    "steps": steps,
                    "message": f"Found non-trivial factor p = {p} in Stage 2 at prime {s2_primes[i + 1]}",
                    "elapsed_us": elapsed_us,
                }

    elapsed_us = int((time.perf_counter() - start_time) * 1_000_000)
    return {
        "n": n,
        "factor": None,
        "cofactor": None,
        "success": False,
        "stage": None,
        "b1": b1,
        "b2": b2,
        "base": a,
        "p_minus_1": None,
        "p_minus_1_factors": None,
        "p_minus_1_expression": None,
        "steps": steps,
        "message": f"No non-trivial factor found with B1 = {b1}, B2 = {b2}. Try increasing bounds or changing base.",
        "elapsed_us": elapsed_us,
    }

def lenstra_ecm(
    n: int,
    b1: int = 1000,
    b2: int | None = None,
    max_curves: int = 20,
    seed: int = 42,
) -> dict:
    """Factor an integer n using Lenstra's Elliptic Curve Method (ECM).

    Tests random elliptic curves with two-stage smoothness bounds (B1, B2).
    Returns a dict with factor, cofactor, stage, curve details, and steps.
    """
    if isinstance(n, str):
        data = n.strip()
        if data.startswith(("0x", "0X")):
            n = int(data, 16)
        else:
            n = int(data)

    if not isinstance(n, int):
        raise TypeError("n must be an integer")
    if n < 2:
        raise ValueError("n must be an integer >= 2")

    if not isinstance(b1, int) or b1 < 2:
        raise ValueError("b1 must be an integer >= 2")

    if b2 is None:
        b2 = b1 * 10
    elif not isinstance(b2, int) or b2 < b1:
        raise ValueError(f"b2 must be an integer >= b1 ({b1})")

    if not isinstance(max_curves, int) or max_curves < 1:
        raise ValueError("max_curves must be an integer >= 1")

    start_time = time.perf_counter()

    if n % 2 == 0:
        elapsed_us = int((time.perf_counter() - start_time) * 1_000_000)
        return {
            "n": n,
            "factor": 2,
            "cofactor": n // 2,
            "success": True,
            "curve_index": 0,
            "curves_tested": 0,
            "stage": 0,
            "b1": b1,
            "b2": b2,
            "a": None,
            "b": None,
            "steps": 0,
            "message": "n is even; factor 2 found immediately",
            "elapsed_us": elapsed_us,
        }

    if is_prime(n):
        elapsed_us = int((time.perf_counter() - start_time) * 1_000_000)
        return {
            "n": n,
            "factor": None,
            "cofactor": None,
            "success": False,
            "curve_index": 0,
            "curves_tested": 0,
            "stage": None,
            "b1": b1,
            "b2": b2,
            "a": None,
            "b": None,
            "steps": 0,
            "message": f"n = {n} is already prime",
            "elapsed_us": elapsed_us,
        }

    def primes_up_to(limit: int) -> list:
        s = [True] * (limit + 1)
        s[0] = s[1] = False
        for p in range(2, int(limit**0.5) + 1):
            if s[p]:
                for i in range(p * p, limit + 1, p):
                    s[i] = False
        return [i for i, v in enumerate(s) if v]

    all_primes = primes_up_to(b2)
    s1_primes = [p for p in all_primes if p <= b1]
    s2_primes = [p for p in all_primes if p > b1]

    state = seed & 0xFFFFFFFFFFFFFFFF if seed != 0 else 0x853C49E6748FEA9B

    def next_u64() -> int:
        nonlocal state
        state = (state + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
        z = state
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
        return (z ^ (z >> 31)) & 0xFFFFFFFFFFFFFFFF

    def gen_range(low: int, high: int) -> int:
        if low >= high:
            return low
        span = high - low + 1
        return low + (next_u64() % span)

    def ec_add(P, Q, a_val, mod_n):
        if P is None:
            return 1, Q
        if Q is None:
            return 1, P
        x1, y1 = P
        x2, y2 = Q
        if x1 == x2:
            if (y1 + y2) % mod_n == 0:
                return 1, None
            den = (2 * y1) % mod_n
            g, x_inv, _ = extended_gcd(den, mod_n)
            if 1 < g < mod_n:
                return g, None
            if g == mod_n:
                return 1, None
            lam = ((3 * x1 * x1 + a_val) * (x_inv % mod_n)) % mod_n
        else:
            den = (x2 - x1) % mod_n
            g, x_inv, _ = extended_gcd(den, mod_n)
            if 1 < g < mod_n:
                return g, None
            if g == mod_n:
                return 1, None
            lam = ((y2 - y1) * (x_inv % mod_n)) % mod_n
        x3 = (lam * lam - x1 - x2) % mod_n
        y3 = (lam * (x1 - x3) - y1) % mod_n
        return 1, (x3, y3)

    def ec_mul(k_val, pt, a_val, mod_n):
        R = None
        curr_Q = pt
        while k_val > 0:
            if k_val & 1:
                factor_or_one, R = ec_add(R, curr_Q, a_val, mod_n)
                if 1 < factor_or_one < mod_n:
                    return factor_or_one, None
            k_val >>= 1
            if k_val > 0:
                factor_or_one, curr_Q = ec_add(curr_Q, curr_Q, a_val, mod_n)
                if 1 < factor_or_one < mod_n:
                    return factor_or_one, None
        return 1, R

    total_steps = 0

    for curve_idx in range(1, max_curves + 1):
        x0 = gen_range(1, n - 1)
        y0 = gen_range(1, n - 1)
        a_param = gen_range(0, n - 1)
        b_param = (y0 * y0 - x0 * x0 * x0 - a_param * x0) % n

        disc = (4 * a_param * a_param * a_param + 27 * b_param * b_param) % n
        g = gcd(disc, n)
        if 1 < g < n:
            p = min(g, n // g)
            elapsed_us = int((time.perf_counter() - start_time) * 1_000_000)
            return {
                "n": n,
                "factor": p,
                "cofactor": n // p,
                "success": True,
                "curve_index": curve_idx,
                "curves_tested": curve_idx,
                "stage": 0,
                "b1": b1,
                "b2": b2,
                "a": a_param,
                "b": b_param,
                "steps": total_steps,
                "message": f"Discriminant of curve {curve_idx} shares non-trivial factor {p} with n",
                "elapsed_us": elapsed_us,
            }
        if g == n:
            continue

        pt = (x0, y0)

        # Stage 1
        for q in s1_primes:
            total_steps += 1
            pk = q
            while pk * q <= b1:
                pk *= q

            if pt is not None:
                factor_or_one, pt = ec_mul(pk, pt, a_param, n)
                if 1 < factor_or_one < n:
                    p = min(factor_or_one, n // factor_or_one)
                    elapsed_us = int((time.perf_counter() - start_time) * 1_000_000)
                    return {
                        "n": n,
                        "factor": p,
                        "cofactor": n // p,
                        "success": True,
                        "curve_index": curve_idx,
                        "curves_tested": curve_idx,
                        "stage": 1,
                        "b1": b1,
                        "b2": b2,
                        "a": a_param,
                        "b": b_param,
                        "steps": total_steps,
                        "message": f"Found non-trivial factor p = {p} on curve {curve_idx} in Stage 1 at prime {q}",
                        "elapsed_us": elapsed_us,
                    }
                if pt is None:
                    break
            else:
                break

        # Stage 2
        if pt is not None and s2_primes:
            for q in s2_primes:
                total_steps += 1
                factor_or_one, _ = ec_mul(q, pt, a_param, n)
                if 1 < factor_or_one < n:
                    p = min(factor_or_one, n // factor_or_one)
                    elapsed_us = int((time.perf_counter() - start_time) * 1_000_000)
                    return {
                        "n": n,
                        "factor": p,
                        "cofactor": n // p,
                        "success": True,
                        "curve_index": curve_idx,
                        "curves_tested": curve_idx,
                        "stage": 2,
                        "b1": b1,
                        "b2": b2,
                        "a": a_param,
                        "b": b_param,
                        "steps": total_steps,
                        "message": f"Found non-trivial factor p = {p} on curve {curve_idx} in Stage 2 at prime {q}",
                        "elapsed_us": elapsed_us,
                    }

    elapsed_us = int((time.perf_counter() - start_time) * 1_000_000)
    return {
        "n": n,
        "factor": None,
        "cofactor": None,
        "success": False,
        "curve_index": max_curves,
        "curves_tested": max_curves,
        "stage": None,
        "b1": b1,
        "b2": b2,
        "a": None,
        "b": None,
        "steps": total_steps,
        "message": f"No non-trivial factor found across {max_curves} curves with B1 = {b1}, B2 = {b2}. Try increasing curves or bounds.",
        "elapsed_us": elapsed_us,
    }

def dixon_factor(
    n: int,
    b: int = 50,
    max_steps: int = 50000,
) -> dict:
    """Factor an integer n using Dixon's random squares method with factor base and linear algebra over GF(2)."""
    if isinstance(n, str):
        data = n.strip()
        if data.startswith(("0x", "0X")):
            n = int(data, 16)
        else:
            n = int(data)

    if not isinstance(n, int):
        raise TypeError("n must be an integer")
    if n < 2:
        raise ValueError("n must be an integer >= 2")

    if not isinstance(b, int) or b < 10:
        raise ValueError("b must be an integer >= 10")
    if not isinstance(max_steps, int) or max_steps < 10:
        raise ValueError("max_steps must be an integer >= 10")

    start_time = time.perf_counter()

    if n % 2 == 0:
        elapsed_us = int((time.perf_counter() - start_time) * 1_000_000)
        return {
            "n": n,
            "factor": 2,
            "cofactor": n // 2,
            "success": True,
            "b": b,
            "factor_base_size": 1,
            "relations_needed": 1,
            "relations_found": 1,
            "x": None,
            "y": None,
            "congruence": None,
            "steps": 0,
            "message": "n is even; factor 2 found immediately",
            "elapsed_us": elapsed_us,
        }

    root = int(math.isqrt(n))
    if root * root == n:
        elapsed_us = int((time.perf_counter() - start_time) * 1_000_000)
        return {
            "n": n,
            "factor": root,
            "cofactor": root,
            "success": True,
            "b": b,
            "factor_base_size": 1,
            "relations_needed": 1,
            "relations_found": 1,
            "x": root,
            "y": 0,
            "congruence": f"{root}² ≡ 0² (mod {n})",
            "steps": 0,
            "message": f"n = {n} is a perfect square: {root}²",
            "elapsed_us": elapsed_us,
        }

    if is_prime(n):
        elapsed_us = int((time.perf_counter() - start_time) * 1_000_000)
        return {
            "n": n,
            "factor": None,
            "cofactor": None,
            "success": False,
            "b": b,
            "factor_base_size": 0,
            "relations_needed": 0,
            "relations_found": 0,
            "x": None,
            "y": None,
            "congruence": None,
            "steps": 0,
            "message": f"n = {n} is already prime",
            "elapsed_us": elapsed_us,
        }

    # Factor base: [-1, 2] + odd primes p <= b with (n/p) == 1
    def primes_up_to(limit: int) -> list:
        s = [True] * (limit + 1)
        s[0] = s[1] = False
        for p in range(2, int(limit**0.5) + 1):
            if s[p]:
                for i in range(p * p, limit + 1, p):
                    s[i] = False
        return [i for i, v in enumerate(s) if v]

    fb = [-1, 2]
    for p in primes_up_to(b):
        if p > 2 and legendre_symbol(n % p, p) == 1:
            fb.append(p)

    k = len(fb)
    relations_needed = k + 10
    relations = []
    x = root + 1
    steps = 0

    while len(relations) < relations_needed and steps < max_steps:
        steps += 1
        x_sq = (x * x) % n
        r = x_sq - n if x_sq > n // 2 else x_sq

        temp = r
        exps = [0] * k
        if temp < 0:
            exps[0] = 1
            temp = -temp

        for idx in range(1, k):
            p = fb[idx]
            while temp > 0 and temp % p == 0:
                exps[idx] += 1
                temp //= p

        if temp == 1:
            parity = [e % 2 for e in exps]
            relations.append((x, r, exps, parity))

        x += 1

    relations_found = len(relations)
    if relations_found <= k:
        elapsed_us = int((time.perf_counter() - start_time) * 1_000_000)
        return {
            "n": n,
            "factor": None,
            "cofactor": None,
            "success": False,
            "b": b,
            "factor_base_size": k,
            "relations_needed": relations_needed,
            "relations_found": relations_found,
            "x": None,
            "y": None,
            "congruence": None,
            "steps": steps,
            "message": f"Insufficient smooth relations found (found {relations_found} of {relations_needed} needed with B = {b}). Try increasing B or max steps.",
            "elapsed_us": elapsed_us,
        }

    # Gaussian elimination over GF(2)
    num_rels = len(relations)
    matrix = [list(r[3]) for r in relations]
    tracking = [[1 if i == j else 0 for j in range(num_rels)] for i in range(num_rels)]
    is_pivot = [False] * num_rels

    for col in range(k):
        found = None
        for r in range(num_rels):
            if not is_pivot[r] and matrix[r][col] == 1:
                found = r
                break
        if found is not None:
            is_pivot[found] = True
            pivot_row = list(matrix[found])
            pivot_track = list(tracking[found])
            for r in range(num_rels):
                if r != found and matrix[r][col] == 1:
                    matrix[r] = [a ^ b for a, b in zip(matrix[r], pivot_row)]
                    tracking[r] = [a ^ b for a, b in zip(tracking[r], pivot_track)]

    # Search null-space rows
    for r in range(num_rels):
        if all(v == 0 for v in matrix[r]):
            subset = [i for i, bit in enumerate(tracking[r]) if bit == 1]
            if not subset:
                continue

            x_prod = 1
            for idx in subset:
                x_prod = (x_prod * relations[idx][0]) % n

            total_exp = [0] * k
            for idx in subset:
                for j in range(k):
                    total_exp[j] += relations[idx][2][j]

            if total_exp[0] % 2 != 0:
                continue

            y_prod = 1
            all_even = True
            for j in range(1, k):
                if total_exp[j] % 2 != 0:
                    all_even = False
                    break
                y_prod = (y_prod * pow(fb[j], total_exp[j] // 2, n)) % n

            if not all_even:
                continue

            x_val = x_prod
            y_val = y_prod

            for diff in [abs(x_val - y_val), x_val + y_val]:
                g = math.gcd(diff, n)
                if 1 < g < n:
                    p = min(g, n // g)
                    q = n // p
                    elapsed_us = int((time.perf_counter() - start_time) * 1_000_000)
                    return {
                        "n": n,
                        "factor": p,
                        "cofactor": q,
                        "success": True,
                        "b": b,
                        "factor_base_size": k,
                        "relations_needed": relations_needed,
                        "relations_found": relations_found,
                        "x": x_val,
                        "y": y_val,
                        "congruence": f"{x_val}² ≡ {y_val}² (mod {n})",
                        "steps": steps,
                        "message": f"Found factor p = {p} via congruent squares {x_val}² ≡ {y_val}² (mod {n})",
                        "elapsed_us": elapsed_us,
                    }

    elapsed_us = int((time.perf_counter() - start_time) * 1_000_000)
    return {
        "n": n,
        "factor": None,
        "cofactor": None,
        "success": False,
        "b": b,
        "factor_base_size": k,
        "relations_needed": relations_needed,
        "relations_found": relations_found,
        "x": None,
        "y": None,
        "congruence": None,
        "steps": steps,
        "message": f"Collected {relations_found} smooth relations with B = {b}, but null-space congruences yielded trivial factors. Try increasing B.",
        "elapsed_us": elapsed_us,
    }

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

def is_perfect(n: int) -> bool:
    """Return True if n is a perfect number (sum of proper divisors equals n)."""
    if not isinstance(n, int):
        raise TypeError("n must be an integer")
    if n <= 1:
        return False
    return sigma(n) == 2 * n

def mersenne_number(p: int) -> int:
    """Return the Mersenne number M_p = 2^p - 1 for exponent p >= 1."""
    if not isinstance(p, int):
        raise TypeError("Exponent p must be an integer")
    if p < 1:
        raise ValueError("Exponent p must be an integer >= 1")
    return (1 << p) - 1

def lucas_lehmer(p: int) -> bool:
    """Lucas-Lehmer deterministic primality test for Mersenne number M_p = 2^p - 1.

    Returns True if M_p is prime, False otherwise.
    Raises ValueError if p < 2.
    """
    if not isinstance(p, int):
        raise TypeError("Exponent p must be an integer")
    if p < 2:
        raise ValueError("Exponent p must be an integer >= 2")
    if p == 2:
        return True
    if p % 2 == 0 or not is_prime(p):
        return False

    mp = (1 << p) - 1
    s = 4
    for _ in range(p - 2):
        s = (s * s - 2) % mp
    return s == 0

def mersenne(p: int) -> dict:
    """Inspect Mersenne number M_p = 2^p - 1 and verify its primality."""
    if not isinstance(p, int):
        raise TypeError("Exponent p must be an integer")
    if p < 2:
        raise ValueError("Exponent p must be an integer >= 2")

    mp = mersenne_number(p)
    p_prime = is_prime(p)
    prime = lucas_lehmer(p) if p_prime else False

    info = {
        "p": p,
        "mp": mp,
        "p_is_prime": p_prime,
        "is_prime": prime,
        "bits": p,
        "digits": len(str(mp)),
    }
    if prime:
        info["perfect_number"] = (1 << (p - 1)) * mp
    return info

def mersenne_primes(limit_p: int = 127) -> list:
    """Return a list of prime exponents p <= limit_p such that 2^p - 1 is prime."""
    if not isinstance(limit_p, int):
        raise TypeError("limit_p must be an integer")
    if limit_p < 2:
        return []

    primes = []
    for p in range(2, limit_p + 1):
        if is_prime(p) and lucas_lehmer(p):
            primes.append(p)
    return primes

def perfect_number(p: int) -> int:
    """Generate the even perfect number 2^(p-1) * (2^p - 1) for Mersenne prime exponent p.

    Raises ValueError if 2^p - 1 is not prime.
    """
    if not isinstance(p, int):
        raise TypeError("Exponent p must be an integer")
    if p < 2:
        raise ValueError("Exponent p must be an integer >= 2")
    if not lucas_lehmer(p):
        raise ValueError(f"Mersenne number 2^{p} - 1 is not prime; no perfect number corresponds to p = {p}")
    return (1 << (p - 1)) * ((1 << p) - 1)

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

# ------------------- Cryptosystems (RSA & Diffie-Hellman) -------------------

def _generate_prime(bits: int, rng: random.Random) -> int:
    """Generate a probable prime with the exact bit length specified."""
    if bits < 2:
        raise ValueError("Prime bit length must be at least 2")
    if bits == 2:
        return rng.choice([2, 3])
    while True:
        val = rng.getrandbits(bits)
        val |= (1 << (bits - 1)) | 1
        if is_prime(val):
            return val

def rsa_keygen(bits: int = 32, e: int = 65537, seed: int = None) -> dict:
    """Generate an RSA key pair (p, q, n, e, d) with CRT components.

    bits: Total modulus bit length (n = p * q). Must be >= 8.
    e: Public exponent (default 65537). Must be coprime to lambda(n).
    seed: Optional random seed for deterministic key generation.
    """
    if not isinstance(bits, int) or not isinstance(e, int):
        raise TypeError("bits and e must be integers")
    if bits < 8:
        raise ValueError("RSA modulus bit length must be at least 8")
    if e < 3 or e % 2 == 0:
        raise ValueError("Public exponent e must be an odd integer >= 3")

    rng = random.Random(seed)
    p_bits = bits // 2
    q_bits = bits - p_bits

    for _ in range(5000):
        p = _generate_prime(p_bits, rng)
        q = _generate_prime(q_bits, rng)
        if p == q:
            continue
        if p < q:
            p, q = q, p

        n = p * q
        phi = (p - 1) * (q - 1)
        lam = phi // gcd(p - 1, q - 1)

        if gcd(e, lam) == 1:
            d = mod_inverse(e, lam)
            dp = d % (p - 1)
            dq = d % (q - 1)
            qinv = mod_inverse(q, p)
            return {
                "p": p,
                "q": q,
                "n": n,
                "e": e,
                "d": d,
                "phi": phi,
                "lambda": lam,
                "dp": dp,
                "dq": dq,
                "qinv": qinv,
            }

    raise RuntimeError(f"Failed to generate RSA key with bits={bits} and e={e}")

def rsa_encrypt(m: int, e: int, n: int) -> int:
    """Encrypt message m with public key (e, n): c = m^e (mod n)."""
    if not isinstance(m, int) or not isinstance(e, int) or not isinstance(n, int):
        raise TypeError("m, e, and n must be integers")
    if n <= 1:
        raise ValueError("Modulus n must be greater than 1")
    if e < 1:
        raise ValueError("Exponent e must be >= 1")
    if m < 0 or m >= n:
        raise ValueError(f"Message m must satisfy 0 <= m < n (got m={m}, n={n})")
    return pow(m, e, n)

def rsa_decrypt(c: int, d: int, n: int, p: int = None, q: int = None) -> int:
    """Decrypt ciphertext c with private key (d, n), optionally using CRT acceleration.

    If factors p and q are provided (where p * q == n), evaluates decryption
    via the Chinese Remainder Theorem:
      m1 = c^(d mod p-1) mod p
      m2 = c^(d mod q-1) mod q
      m = m2 + (q^-1 * (m1 - m2) mod p) * q
    """
    if not isinstance(c, int) or not isinstance(d, int) or not isinstance(n, int):
        raise TypeError("c, d, and n must be integers")
    if n <= 1:
        raise ValueError("Modulus n must be greater than 1")
    if d < 1:
        raise ValueError("Exponent d must be >= 1")
    if c < 0 or c >= n:
        raise ValueError(f"Ciphertext c must satisfy 0 <= c < n (got c={c}, n={n})")

    if (p is None) != (q is None):
        raise ValueError("Both p and q must be provided for CRT acceleration, or neither")

    if p is not None and q is not None:
        if not isinstance(p, int) or not isinstance(q, int):
            raise TypeError("p and q must be integers")
        if p * q != n:
            raise ValueError(f"Invalid factors: p * q != n ({p} * {q} != {n})")
        dp = d % (p - 1)
        dq = d % (q - 1)
        qinv = mod_inverse(q, p)
        m1 = pow(c, dp, p)
        m2 = pow(c, dq, q)
        h = (qinv * (m1 - m2)) % p
        return m2 + h * q

    return pow(c, d, n)

def rsa_crack(n: int, e: int, c: int = None, method: str = "auto") -> dict:
    """Recover RSA private key (d) by factoring n, and decrypt c if given.

    method: 'auto' (tries Fermat up to 10k steps, then Pollard's rho),
            'fermat' (uses Fermat's difference of squares), or
            'pollard' (uses Pollard's rho & trial division).
    """
    if not isinstance(n, int) or not isinstance(e, int):
        raise TypeError("n and e must be integers")
    if n <= 1:
        raise ValueError("Modulus n must be greater than 1")
    if e < 1:
        raise ValueError("Exponent e must be >= 1")
    if method not in ("auto", "pollard", "fermat"):
        raise ValueError(f"Unknown factorization method '{method}'; choose 'auto', 'pollard', or 'fermat'")
    if c is not None:
        if not isinstance(c, int):
            raise TypeError("c must be an integer")
        if c < 0 or c >= n:
            raise ValueError(f"Ciphertext c must satisfy 0 <= c < n (got c={c}, n={n})")

    p = q = None
    used_method = method

    if method in ("auto", "fermat"):
        max_f_steps = 10000 if method == "auto" else 1000000
        if n % 2 != 0:
            fermat_res = fermat_factor(n, quiet=True)
            if fermat_res is not None and fermat_res[2] <= max_f_steps:
                fp, fq, _ = fermat_res
                if fp > 1 and fq > 1 and is_prime(fp) and is_prime(fq):
                    p, q = fp, fq
                    used_method = "fermat"

        if method == "fermat" and p is None:
            raise ValueError(f"Fermat factorization failed to factor modulus {n}")

    if p is None:
        factors = factorize(n)
        primes = []
        for prime, power in factors:
            primes.extend([prime] * power)

        if len(primes) != 2:
            raise ValueError(f"Modulus {n} must be a product of two primes (found {len(primes)} factors)")

        p, q = primes[0], primes[1]
        used_method = "pollard"

    if p < q:
        p, q = q, p

    phi = (p - 1) * (q - 1)
    lam = phi // gcd(p - 1, q - 1)

    if gcd(e, lam) != 1:
        raise ValueError(f"Cannot invert e={e} modulo lambda(n)={lam}")

    d = mod_inverse(e, lam)
    res = {
        "n": n,
        "e": e,
        "p": p,
        "q": q,
        "phi": phi,
        "lambda": lam,
        "d": d,
        "method": used_method,
    }
    if c is not None:
        res["c"] = c
        res["m"] = rsa_decrypt(c, d, n, p, q)
    return res

def dh_keygen(p: int, g: int, private_key: int = None, seed: int = None) -> dict:
    """Generate a Diffie-Hellman key pair: private scalar a and public key A = g^a mod p."""
    if not isinstance(p, int) or not isinstance(g, int):
        raise TypeError("p and g must be integers")
    if p <= 2 or not is_prime(p):
        raise ValueError("Diffie-Hellman modulus p must be a prime > 2")
    if g <= 1 or g >= p:
        raise ValueError(f"Generator g must satisfy 1 < g < p (got g={g})")

    if private_key is not None:
        if not isinstance(private_key, int):
            raise TypeError("private_key must be an integer")
        if private_key < 2 or private_key >= p - 1:
            raise ValueError(f"private_key must satisfy 2 <= private_key <= p-2 (got {private_key})")
        a = private_key
    else:
        rng = random.Random(seed)
        a = rng.randint(2, p - 2)

    A = pow(g, a, p)
    return {
        "p": p,
        "g": g,
        "private_key": a,
        "public_key": A,
    }

def dh_derive(private_key: int, peer_public_key: int, p: int) -> int:
    """Derive Diffie-Hellman shared secret: S = peer_public_key^private_key mod p."""
    if not isinstance(private_key, int) or not isinstance(peer_public_key, int) or not isinstance(p, int):
        raise TypeError("private_key, peer_public_key, and p must be integers")
    if p <= 2:
        raise ValueError("Modulus p must be greater than 2")
    if private_key < 1:
        raise ValueError("private_key must be >= 1")
    if peer_public_key < 1 or peer_public_key >= p:
        raise ValueError(f"peer_public_key must satisfy 1 <= key < p (got {peer_public_key})")
    return pow(peer_public_key, private_key, p)

def dh_crack(g: int, public_key: int, p: int) -> int:
    """Recover Diffie-Hellman private key from public key using BSGS discrete logarithm solver."""
    return discrete_log(g, public_key, p)

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

_RSA_PRIMES = [3, 5, 7, 11, 13, 17, 19, 23]

def _generate_rsa_problem(rng: random.Random) -> dict:
    p, q = rng.sample(_RSA_PRIMES, 2)
    n = p * q
    phi = (p - 1) * (q - 1)
    possible_e = [e for e in (3, 5, 7, 11, 13, 17) if gcd(e, phi) == 1]
    e = rng.choice(possible_e)
    m = rng.randint(2, min(n - 1, 15))
    c = pow(m, e, n)
    return {
        "topic": "rsa",
        "question": f"In RSA with e = {e} and n = {n}, find ciphertext c for message m = {m} (c ≡ m^{e} mod {n})",
        "answer": c,
    }

_DH_PARAMS = [
    (11, 2), (13, 2), (17, 3), (19, 2), (23, 5), (29, 2), (31, 3), (41, 6), (43, 3)
]

def _generate_dh_problem(rng: random.Random) -> dict:
    p, g = rng.choice(_DH_PARAMS)
    a = rng.randint(2, 6)
    if rng.choice([True, False]):
        A = pow(g, a, p)
        return {
            "topic": "dh",
            "question": f"In Diffie-Hellman with p = {p} and generator g = {g}, find public key A for private key a = {a} (A ≡ {g}^{a} mod {p})",
            "answer": A,
        }
    else:
        b = rng.randint(2, 6)
        B = pow(g, b, p)
        s = pow(B, a, p)
        return {
            "topic": "dh",
            "question": f"In Diffie-Hellman with p = {p}, Alice's private key is a = {a} and Bob's public key is B = {B}. Find shared secret S (S ≡ {B}^{a} mod {p})",
            "answer": s,
        }

_MERSENNE_CHALLENGE_EXPONENTS = [2, 3, 4, 5, 6, 7, 11, 13, 17, 19, 23]

def _generate_mersenne_problem(rng: random.Random) -> dict:
    p = rng.choice(_MERSENNE_CHALLENGE_EXPONENTS)
    ans = 1 if (p >= 2 and is_prime(p) and lucas_lehmer(p)) else 0
    return {
        "topic": "mersenne",
        "question": f"Is the Mersenne number M_{p} = 2^{p} - 1 prime? (Enter 1 for prime, 0 for composite)",
        "answer": ans,
    }

_PERFECT_CHALLENGE_NUMBERS = [
    (6, 1),
    (12, 0),
    (28, 1),
    (36, 0),
    (496, 1),
    (512, 0),
    (8128, 1),
    (8192, 0),
]

def _generate_perfect_problem(rng: random.Random) -> dict:
    if rng.choice([True, False]):
        n, ans = rng.choice(_PERFECT_CHALLENGE_NUMBERS)
        return {
            "topic": "perfect",
            "question": f"Is n = {n} a perfect number? (Enter 1 for yes, 0 for no)",
            "answer": ans,
        }
    else:
        p = rng.choice([2, 3, 5, 7])
        ans = perfect_number(p)
        return {
            "topic": "perfect",
            "question": f"What is the even perfect number corresponding to Mersenne prime M_{p} = 2^{p} - 1 for p = {p}?",
            "answer": ans,
        }

def _generate_pm1_problem(rng: random.Random) -> dict:
    smooth_primes = [29, 37, 41, 73, 97, 101, 109, 113, 137, 181, 193, 241, 257, 337, 433, 577, 601, 641, 769, 1009]
    other_primes = [103, 107, 127, 131, 139, 149, 151, 157, 163, 167, 173, 179]
    p = rng.choice(smooth_primes)
    q = rng.choice(other_primes)
    n = p * q
    return {
        "topic": "pm1",
        "question": f"Smallest prime factor of composite n = {n} (hint: p - 1 is smooth)",
        "answer": min(p, q),
    }

def _generate_ecm_problem(rng: random.Random) -> dict:
    primes_p = [43, 47, 53, 59, 61, 67, 71, 79, 83, 89, 97, 101, 103, 107, 109, 113]
    primes_q = [127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199]
    p = rng.choice(primes_p)
    q = rng.choice(primes_q)
    n = p * q
    return {
        "topic": "ecm",
        "question": f"Smallest prime factor of composite n = {n} (hint: use Lenstra elliptic curve factorization)",
        "answer": min(p, q),
    }

def _generate_dixon_problem(rng: random.Random) -> dict:
    primes_p = [23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71]
    primes_q = [73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131]
    p = rng.choice(primes_p)
    q = rng.choice(primes_q)
    n = p * q
    return {
        "topic": "dixon",
        "question": f"Smallest prime factor of composite n = {n} (hint: use Dixon's random squares method)",
        "answer": min(p, q),
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
    "rsa": _generate_rsa_problem,
    "dh": _generate_dh_problem,
    "mersenne": _generate_mersenne_problem,
    "perfect": _generate_perfect_problem,
    "pm1": _generate_pm1_problem,
    "ecm": _generate_ecm_problem,
    "dixon": _generate_dixon_problem,
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
        "Functions: lucas_lehmer, mersenne, perfect_number, is_perfect,\n"
        "           rsa_keygen, rsa_encrypt, rsa_decrypt, rsa_crack,\n"
        "           dh_keygen, dh_derive, dh_crack, jacobi_symbol, solovay_strassen,\n"
        "           discrete_log (dlog), mobius, sigma, tau, factorize, is_prime, ..."
    )
    print(banner)
    namespace = {
        "__name__": "__console__",
        "__doc__": None,
        "is_prime": is_prime,
        "factorize": factorize,
        "format_factorization": format_factorization,
        "fermat_factor": fermat_factor,
        "predict_factor": predict_factor,
        "to_sci": to_sci,
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
        "is_perfect": is_perfect,
        "mersenne_number": mersenne_number,
        "lucas_lehmer": lucas_lehmer,
        "mersenne": mersenne,
        "mersenne_primes": mersenne_primes,
        "perfect_number": perfect_number,
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
        "rsa_keygen": rsa_keygen,
        "rsa_encrypt": rsa_encrypt,
        "rsa_decrypt": rsa_decrypt,
        "rsa_crack": rsa_crack,
        "dh_keygen": dh_keygen,
        "dh_derive": dh_derive,
        "dh_crack": dh_crack,
        "pollard_pm1": pollard_pm1,
        "pm1": pollard_pm1,
        "lenstra_ecm": lenstra_ecm,
        "ecm": lenstra_ecm,
        "dixon_factor": dixon_factor,
        "dixon": dixon_factor,
    }
    code.interact(banner="", local=namespace, exitmsg="")

# ------------------- CLI -------------------

def _has_stdin_data() -> bool:
    """Return True if stdin has data waiting to be read without blocking."""
    if sys.stdin.isatty():
        return False
    if sys.platform == "win32":
        try:
            import msvcrt
            import ctypes
            handle = msvcrt.get_osfhandle(sys.stdin.fileno())
            avail = ctypes.c_ulong(0)
            if ctypes.windll.kernel32.PeekNamedPipe(handle, None, 0, None, ctypes.byref(avail), None):
                return avail.value > 0
        except Exception:
            pass
        return False
    try:
        import select
        return bool(select.select([sys.stdin], [], [], 0.0)[0])
    except Exception:
        return False

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

    # Fermat difference of squares
    fermat_cmd = subparsers.add_parser("fermat", help="Fermat's difference of squares factorization")
    fermat_cmd.add_argument("n", type=str, nargs="?", default=None, help="Modulus to factor (decimal or hex with 0x)")
    fermat_cmd.add_argument("-n", "--mod", dest="mod", type=str, default=None, help="Modulus to factor (decimal or hex with 0x)")
    fermat_cmd.add_argument("-q", "--quiet", action="store_true", help="Suppress progress output")
    fermat_cmd.add_argument("-p", "--prec", type=int, default=30, help="Digits to show for approximate factors (default: 30)")
    fermat_cmd.add_argument("-i", "--interval", type=int, default=1000000, help="Print interval for progress (default: 1e6)")
    fermat_cmd.add_argument("-s", "--start", type=int, default=0, help="Starting iteration offset (default: 0)")

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

    # RSA
    rsa_parser = subparsers.add_parser("rsa", help="Toy RSA cryptosystem (keygen, encrypt, decrypt, crack)")
    rsa_sub = rsa_parser.add_subparsers(dest="rsa_action", required=True)

    rsa_keygen_cmd = rsa_sub.add_parser("keygen", help="Generate RSA keypair")
    rsa_keygen_cmd.add_argument("--bits", type=int, default=32, help="Modulus bit length (default: 32)")
    rsa_keygen_cmd.add_argument("--e", type=int, default=65537, help="Public exponent e (default: 65537)")
    rsa_keygen_cmd.add_argument("--seed", type=int, default=None, help="Deterministic random seed")

    rsa_enc_cmd = rsa_sub.add_parser("encrypt", help="Encrypt message: c = m^e mod n")
    rsa_enc_cmd.add_argument("m", type=int, help="Message integer (0 <= m < n)")
    rsa_enc_cmd.add_argument("e", type=int, help="Public exponent e")
    rsa_enc_cmd.add_argument("n", type=int, help="Modulus n")

    rsa_dec_cmd = rsa_sub.add_parser("decrypt", help="Decrypt ciphertext: m = c^d mod n")
    rsa_dec_cmd.add_argument("c", type=int, help="Ciphertext integer (0 <= c < n)")
    rsa_dec_cmd.add_argument("d", type=int, help="Private exponent d")
    rsa_dec_cmd.add_argument("n", type=int, help="Modulus n")
    rsa_dec_cmd.add_argument("--p", type=int, default=None, help="Prime factor p for CRT acceleration")
    rsa_dec_cmd.add_argument("--q", type=int, default=None, help="Prime factor q for CRT acceleration")

    rsa_crack_cmd = rsa_sub.add_parser("crack", help="Factor modulus n, recover d, and decrypt c if given")
    rsa_crack_cmd.add_argument("n", type=int, help="Modulus n")
    rsa_crack_cmd.add_argument("e", type=int, help="Public exponent e")
    rsa_crack_cmd.add_argument("c", type=int, nargs="?", default=None, help="Optional ciphertext c to decrypt")
    rsa_crack_cmd.add_argument("--c", dest="opt_c", type=int, default=None, help="Optional ciphertext c to decrypt")
    rsa_crack_cmd.add_argument("--method", choices=["auto", "pollard", "fermat"], default="auto", help="Factorization attack method (default: auto)")

    # Diffie-Hellman
    dh_parser = subparsers.add_parser("dh", help="Diffie-Hellman key exchange (keygen, derive, crack)")
    dh_sub = dh_parser.add_subparsers(dest="dh_action", required=True)

    dh_keygen_cmd = dh_sub.add_parser("keygen", help="Generate private/public keypair")
    dh_keygen_cmd.add_argument("p", type=int, help="Prime modulus p (> 2)")
    dh_keygen_cmd.add_argument("g", type=int, help="Generator g (1 < g < p)")
    dh_keygen_cmd.add_argument("--priv", "--private-key", dest="priv", type=int, default=None, help="Custom private scalar a (2 <= a <= p-2)")
    dh_keygen_cmd.add_argument("--seed", type=int, default=None, help="Deterministic random seed")

    dh_derive_cmd = dh_sub.add_parser("derive", help="Derive shared secret: S = peer_public^private mod p")
    dh_derive_cmd.add_argument("private_key", type=int, help="Your private key scalar a")
    dh_derive_cmd.add_argument("peer_public_key", type=int, help="Peer's public key B")
    dh_derive_cmd.add_argument("p", type=int, help="Prime modulus p")

    dh_crack_cmd = dh_sub.add_parser("crack", help="Recover private key a from public key A = g^a mod p via BSGS")
    dh_crack_cmd.add_argument("g", type=int, help="Generator g")
    dh_crack_cmd.add_argument("public_key", type=int, help="Public key A")
    dh_crack_cmd.add_argument("p", type=int, help="Prime modulus p")

    # Mersenne numbers & Lucas–Lehmer primality test
    mersenne_cmd = subparsers.add_parser("mersenne", help="Mersenne number M_p = 2^p - 1 and Lucas–Lehmer test")
    mersenne_cmd.add_argument("p", type=int, nargs="?", default=None, help="Mersenne exponent p >= 2")
    mersenne_cmd.add_argument("--list", action="store_true", help="List Mersenne primes up to --limit")
    mersenne_cmd.add_argument("--limit", type=int, default=127, help="Maximum exponent for --list (default: 127)")

    # Perfect numbers
    perfect_cmd = subparsers.add_parser("perfect", help="Perfect number generator and test (Euclid–Euler theorem)")
    perfect_cmd.add_argument("n", type=int, help="Integer n to test, or exponent p if --from-p is set")
    perfect_cmd.add_argument("-p", "--from-p", action="store_true", help="Generate even perfect number from Mersenne exponent p")

    # Pollard p - 1
    pm1_cmd = subparsers.add_parser("pm1", help="Pollard's p - 1 integer factorization")
    pm1_cmd.add_argument("n", type=int, help="Integer to factor (n >= 2)")
    pm1_cmd.add_argument("--b1", type=int, default=1000, help="Stage 1 smoothness bound (default: 1000)")
    pm1_cmd.add_argument("--b2", type=int, default=None, help="Stage 2 smoothness bound (default: b1 * 10)")
    pm1_cmd.add_argument("-a", "--base", type=int, default=2, help="Base a (default: 2)")

    # Lenstra's ECM
    ecm_cmd = subparsers.add_parser("ecm", help="Lenstra's elliptic curve factorization (ECM)")
    ecm_cmd.add_argument("n", type=int, help="Integer to factor (n >= 2)")
    ecm_cmd.add_argument("--b1", type=int, default=1000, help="Stage 1 smoothness bound (default: 1000)")
    ecm_cmd.add_argument("--b2", type=int, default=None, help="Stage 2 smoothness bound (default: b1 * 10)")
    ecm_cmd.add_argument("--curves", type=int, default=20, help="Maximum curves to test (default: 20)")
    ecm_cmd.add_argument("--seed", type=int, default=42, help="Random seed for reproducible curve generation (default: 42)")

    # Dixon's random squares
    dixon_cmd = subparsers.add_parser("dixon", help="Dixon's random squares integer factorization")
    dixon_cmd.add_argument("n", type=int, help="Integer to factor (n >= 2)")
    dixon_cmd.add_argument("--b", type=int, default=50, help="Factor base smoothness bound (default: 50)")
    dixon_cmd.add_argument("--max-steps", type=int, default=50000, help="Maximum search steps (default: 50000)")

    # Interactive REPL
    subparsers.add_parser("repl", help="Interactive FoxMath REPL")

    # Challenge / quiz mode
    challenge = subparsers.add_parser("challenge", help="Practice problems (Legendre, Jacobi, CRT, SqrtMod, Factor, Totient, DLog, RSA, DH, Mersenne, Perfect, PM1, ECM, Dixon)")
    challenge.add_argument("--topic", choices=["legendre", "crt", "ecdlp", "sqrtmod", "factor", "totient", "dlog", "jacobi", "rsa", "dh", "mersenne", "perfect", "pm1", "ecm", "dixon", "mixed"], default="mixed")
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

        elif cmd == "fermat" or args.command == "fermat":
            raw_input = args.n if args.n is not None else args.mod
            if raw_input is None:
                if _has_stdin_data():
                    raw_input = sys.stdin.read().strip()
                else:
                    raise ValueError("Modulus to factor must be provided via argument (-n/positional) or stdin")

            raw_input = raw_input.strip()
            if raw_input == "-":
                raw_input = sys.stdin.read().strip()

            if not raw_input:
                raise ValueError("No modulus provided to factor")

            if raw_input.startswith(("0x", "0X")):
                target_n = int(raw_input, 16)
            else:
                target_n = int(raw_input)

            if args.start < 0:
                raise ValueError("starting iteration (--start) must be non-negative")

            start_time = time.perf_counter()
            res = fermat_factor(
                target_n,
                quiet=args.quiet,
                prec=args.prec,
                print_interval=args.interval,
                start_iter=args.start,
            )
            elapsed = time.perf_counter() - start_time

            if res is None:
                if json_output:
                    result.update({
                        "command": "fermat",
                        "n": target_n,
                        "success": False,
                        "error": f"Failed to factor {target_n}",
                    })
                    print(json.dumps(result, indent=2))
                else:
                    print(f"❌ Failed to factor {target_n}.")
                sys.exit(1)

            p, q, k = res
            result.update({
                "command": "fermat",
                "n": target_n,
                "p": p,
                "q": q,
                "iterations": k,
                "elapsed_seconds": round(elapsed, 6),
                "success": True,
            })
            if not json_output:
                print(f"\n✅ Factors of {target_n}:")
                print(f"p = {p}")
                print(f"q = {q}")
                print(f"Iterations = {k}")
                print(f"⏱️  Execution time: {elapsed:.6f} seconds")

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

        elif cmd == "rsa" or args.command == "rsa":
            if args.rsa_action == "keygen":
                keys = rsa_keygen(bits=args.bits, e=args.e, seed=args.seed)
                result.update({"command": "rsa", "action": "keygen", **keys})
                if not json_output:
                    print(f"RSA {args.bits}-bit Keypair:")
                    print(f"  Modulus   (n) = {keys['n']}")
                    print(f"  Public  e     = {keys['e']}")
                    print(f"  Private d     = {keys['d']}")
                    print(f"  Prime p       = {keys['p']}")
                    print(f"  Prime q       = {keys['q']}")
                    print(f"  Euler phi     = {keys['phi']}")
                    print(f"  Carmichael λ  = {keys['lambda']}")

            elif args.rsa_action == "encrypt":
                c = rsa_encrypt(args.m, args.e, args.n)
                result.update({"command": "rsa", "action": "encrypt", "m": args.m, "e": args.e, "n": args.n, "c": c})
                if not json_output:
                    print(f"Ciphertext c = {c}")

            elif args.rsa_action == "decrypt":
                m = rsa_decrypt(args.c, args.d, args.n, p=args.p, q=args.q)
                crt_used = args.p is not None and args.q is not None
                result.update({"command": "rsa", "action": "decrypt", "c": args.c, "d": args.d, "n": args.n, "m": m, "crt": crt_used})
                if not json_output:
                    suffix = " (CRT accelerated)" if crt_used else ""
                    print(f"Decrypted message m = {m}{suffix}")

            elif args.rsa_action == "crack":
                c_val = args.c if args.c is not None else args.opt_c
                cracked = rsa_crack(args.n, args.e, c=c_val, method=args.method)
                result.update({"command": "rsa", "action": "crack", **cracked})
                if not json_output:
                    print(f"Factored n = {args.n}:")
                    print(f"  p = {cracked['p']}")
                    print(f"  q = {cracked['q']}")
                    print(f"  Method       = {cracked['method']}")
                    print(f"  Euler phi    = {cracked['phi']}")
                    print(f"  Carmichael λ = {cracked['lambda']}")
                    print(f"  Recovered d  = {cracked['d']}")
                    if "m" in cracked:
                        print(f"  Decrypted m  = {cracked['m']}")

        elif cmd == "dh" or args.command == "dh":
            if args.dh_action == "keygen":
                keys = dh_keygen(args.p, args.g, private_key=args.priv, seed=args.seed)
                result.update({"command": "dh", "action": "keygen", **keys})
                if not json_output:
                    print(f"Diffie-Hellman Keypair (p={args.p}, g={args.g}):")
                    print(f"  Private key (a) = {keys['private_key']}")
                    print(f"  Public key  (A) = {keys['public_key']}")

            elif args.dh_action == "derive":
                s = dh_derive(args.private_key, args.peer_public_key, args.p)
                result.update({
                    "command": "dh",
                    "action": "derive",
                    "private_key": args.private_key,
                    "peer_public_key": args.peer_public_key,
                    "p": args.p,
                    "shared_secret": s,
                })
                if not json_output:
                    print(f"Shared secret S = {s}")

            elif args.dh_action == "crack":
                priv = dh_crack(args.g, args.public_key, args.p)
                result.update({
                    "command": "dh",
                    "action": "crack",
                    "g": args.g,
                    "public_key": args.public_key,
                    "p": args.p,
                    "private_key": priv,
                })
                if not json_output:
                    print("Discrete logarithm solved via BSGS:")
                    print(f"  Recovered private key a = {priv}")
                    print(f"  Verification: {args.g}^{priv} ≡ {args.public_key} (mod {args.p})")

        elif cmd == "mersenne" or args.command == "mersenne":
            if args.list:
                m_primes = mersenne_primes(args.limit)
                result.update({
                    "command": "mersenne",
                    "action": "list",
                    "limit": args.limit,
                    "exponents": m_primes,
                    "count": len(m_primes),
                })
                if not json_output:
                    print(f"Mersenne prime exponents p <= {args.limit}:")
                    print(f"  {', '.join(str(p) for p in m_primes)} ({len(m_primes)} primes)")
            else:
                if args.p is None:
                    raise ValueError("Specify exponent p or use --list")
                m_info = mersenne(args.p)
                result.update({"command": "mersenne", **m_info})
                if not json_output:
                    mp_str = str(m_info["mp"])
                    disp = mp_str if len(mp_str) <= 40 else f"{mp_str[:15]}...{mp_str[-15:]} ({m_info['digits']} digits)"
                    print(f"M_{args.p} = 2^{args.p} - 1 = {disp}")
                    if m_info["is_prime"]:
                        print("  ✓ Prime (Lucas-Lehmer test passed)")
                        print(f"  Corresponding perfect number: {m_info['perfect_number']}")
                    else:
                        if not m_info["p_is_prime"]:
                            print(f"  ✗ Composite (exponent {args.p} is composite)")
                        else:
                            print("  ✗ Composite (Lucas-Lehmer test failed)")

        elif cmd == "perfect" or args.command == "perfect":
            if args.from_p:
                pn = perfect_number(args.n)
                result.update({
                    "command": "perfect",
                    "action": "generate",
                    "p": args.n,
                    "perfect_number": pn,
                    "mersenne_prime": (1 << args.n) - 1,
                })
                if not json_output:
                    print(f"Even perfect number for p = {args.n}:")
                    print(f"  2^{args.n - 1} * (2^{args.n} - 1) = {pn}")
            else:
                perf = is_perfect(args.n)
                result.update({
                    "command": "perfect",
                    "action": "test",
                    "n": args.n,
                    "is_perfect": perf,
                })
                if not json_output:
                    if perf:
                        print(f"{args.n} is a perfect number (sum of proper divisors equals {args.n})")
                    else:
                        print(f"{args.n} is not a perfect number")

        elif cmd == "pm1" or args.command == "pm1":
            res = pollard_pm1(args.n, b1=args.b1, b2=args.b2, a=args.base)
            result.update({"command": "pm1", **res})
            if not json_output:
                if res["success"]:
                    print(f"Factor found: p = {res['factor']}")
                    print(f"Cofactor:     q = {res['cofactor']}")
                    print(f"Stage:        {res['stage']}")
                    print(f"Bounds:       B1 = {res['b1']}, B2 = {res['b2']}")
                    print(f"p - 1:        {res['p_minus_1']} = {res['p_minus_1_expression']}")
                    print(f"Steps:        {res['steps']}")
                    print(f"Time:         {res['elapsed_us']} µs")
                else:
                    print(res["message"])

        elif cmd == "ecm" or args.command == "ecm":
            res = lenstra_ecm(args.n, b1=args.b1, b2=args.b2, max_curves=args.curves, seed=args.seed)
            result.update({"command": "ecm", **res})
            if not json_output:
                if res["success"]:
                    print(f"Factor found: p = {res['factor']}")
                    print(f"Cofactor:     q = {res['cofactor']}")
                    print(f"Curve index:  {res['curve_index']} (of {res['curves_tested']} tested)")
                    print(f"Stage:        {res['stage']}")
                    print(f"Bounds:       B1 = {res['b1']}, B2 = {res['b2']}")
                    if res["a"] is not None:
                        print(f"Curve params: a = {res['a']}, b = {res['b']}")
                    print(f"Steps:        {res['steps']}")
                    print(f"Time:         {res['elapsed_us']} µs")
                else:
                    print(res["message"])

        elif cmd == "dixon" or args.command == "dixon":
            res = dixon_factor(args.n, b=args.b, max_steps=args.max_steps)
            result.update({"command": "dixon", **res})
            if not json_output:
                if res["success"]:
                    print(f"Factor found: p = {res['factor']}")
                    print(f"Cofactor:     q = {res['cofactor']}")
                    print(f"Bound B:      {res['b']} (Factor base size: {res['factor_base_size']})")
                    print(f"Relations:    {res['relations_found']} of {res['relations_needed']} needed")
                    if res["congruence"] is not None:
                        print(f"Congruence:   {res['congruence']}")
                    print(f"Steps:        {res['steps']}")
                    print(f"Time:         {res['elapsed_us']} µs")
                else:
                    print(res["message"])

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
