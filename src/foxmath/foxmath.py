#!/usr/bin/env python3
"""
FoxMath — Minimalist math & crypto explorer 🦊
One binary. Clean. Educational. Safe.
"""

import argparse
import sys
import json
from decimal import Decimal, getcontext
from pathlib import Path

# ------------------- Math Functions -------------------

def legendre_symbol(a: int, p: int) -> int:
    """Compute Legendre symbol (a/p)"""
    if p < 2 or p % 2 == 0:
        raise ValueError("p must be an odd prime")
    res = pow(a, (p - 1) // 2, p)
    return 1 if res == 1 else (-1 if res == p - 1 else 0)

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

def _extended_gcd(a: int, b: int):
    """Returns (g, x, y) such that a*x + b*y == g == gcd(a, b)."""
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r != 0:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    return old_r, old_s, old_t

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

# ------------------- CLI -------------------

def get_command():
    """Support symlinks: foxmath-legendre, foxmath-pi, etc."""
    name = Path(sys.argv[0]).stem.lower()
    if name.startswith("foxmath-"):
        return name.split("-", 1)[1]
    return None

def main():
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
    subparsers = parser.add_subparsers(dest="command", required=not cmd)

    # Legendre
    leg = subparsers.add_parser("legendre", help="Legendre symbol (a/p)")
    leg.add_argument("a", type=int, help="Integer a")
    leg.add_argument("p", type=int, help="Odd prime p")

    # Pi
    pi_cmd = subparsers.add_parser("pi", help="Euler/Hermann Machin-like π approximation")
    pi_cmd.add_argument("--terms", type=int, default=100)

    # Elliptic curve add
    ec = subparsers.add_parser("ecadd", help="Elliptic curve point addition")
    ec.add_argument("--x1", type=int, required=True)
    ec.add_argument("--y1", type=int, required=True)
    ec.add_argument("--x2", type=int, required=True)
    ec.add_argument("--y2", type=int, required=True)
    ec.add_argument("--a", type=int, default=-3, help="Curve parameter a")
    ec.add_argument("--p", type=int, required=True, help="Prime modulus")

    # Chinese Remainder Theorem
    crt = subparsers.add_parser("crt", help="Solve a system of congruences x = r (mod m)")
    crt.add_argument("--r", type=int, nargs="+", required=True, help="Remainders, e.g. --r 2 3 2")
    crt.add_argument("--m", type=int, nargs="+", required=True, help="Moduli, e.g. --m 3 5 7")

    # Continued fractions
    cf_cmd = subparsers.add_parser("cf", help="Continued fraction expansion and convergents")
    cf_cmd.add_argument("--num", type=int, required=True, help="Numerator")
    cf_cmd.add_argument("--den", type=int, required=True, help="Denominator")
    cf_cmd.add_argument("--terms", type=int, default=100, help="Max terms to expand")

    args = parser.parse_args(argv)

    result = {"tool": "foxmath"}

    try:
        if cmd == "legendre" or args.command == "legendre":
            ls = legendre_symbol(args.a, args.p)
            result.update({"command": "legendre", "a": args.a, "p": args.p, "value": ls})
            if not json_output:
                print(f"Legendre ({args.a}/{args.p}) = {ls}")

        elif cmd == "pi" or args.command == "pi":
            approx = euler_hermann_pi_approx(args.terms)
            result.update({"command": "pi", "terms": args.terms, "approx": str(approx)})
            if not json_output:
                print(f"π ≈ {approx}  ({args.terms} terms)")

        elif cmd == "ecadd" or args.command == "ecadd":
            x3, y3 = ec_point_add(args.x1, args.y1, args.x2, args.y2, args.a, args.p)
            result.update({"command": "ecadd", "result": (x3, y3)})
            if not json_output:
                print(f"Result point: ({x3}, {y3})")

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

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    if json_output:
        print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
