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

def catalan_pi_approx(terms: int = 100) -> Decimal:
    """
    Fast-converging π approximation via the Euler/Hermann Machin-like
    identity: pi/4 = arctan(1/2) + arctan(1/3).
    """
    getcontext().prec = 50
    s = _arctan_series(Decimal(1) / 2, terms) + _arctan_series(Decimal(1) / 3, terms)
    return 4 * s

def ec_point_add(x1, y1, x2, y2, a, p):
    """Simple elliptic curve point addition over finite field y² = x³ + a x + b (mod p)"""
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
    pi_cmd = subparsers.add_parser("pi", help="Catalan-inspired π approximation")
    pi_cmd.add_argument("--terms", type=int, default=100)

    # Elliptic curve add
    ec = subparsers.add_parser("ecadd", help="Elliptic curve point addition")
    ec.add_argument("--x1", type=int, required=True)
    ec.add_argument("--y1", type=int, required=True)
    ec.add_argument("--x2", type=int, required=True)
    ec.add_argument("--y2", type=int, required=True)
    ec.add_argument("--a", type=int, default=-3, help="Curve parameter a")
    ec.add_argument("--p", type=int, required=True, help="Prime modulus")

    args = parser.parse_args(argv)

    result = {"tool": "foxmath"}

    try:
        if cmd == "legendre" or args.command == "legendre":
            ls = legendre_symbol(args.a, args.p)
            result.update({"command": "legendre", "a": args.a, "p": args.p, "value": ls})
            print(f"Legendre ({args.a}/{args.p}) = {ls}")

        elif cmd == "pi" or args.command == "pi":
            approx = catalan_pi_approx(args.terms)
            result.update({"command": "pi", "terms": args.terms, "approx": str(approx)})
            print(f"π ≈ {approx}  ({args.terms} terms)")

        elif cmd == "ecadd" or args.command == "ecadd":
            x3, y3 = ec_point_add(args.x1, args.y1, args.x2, args.y2, args.a, args.p)
            result.update({"command": "ecadd", "result": (x3, y3)})
            print(f"Result point: ({x3}, {y3})")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    if json_output:
        print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
