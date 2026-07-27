#!/usr/bin/env python3
"""
FoxMath — Minimalist math & crypto explorer 🦊
One binary. Clean. Educational. Safe.
"""

import argparse
import sys
import json
import random
from decimal import Decimal, getcontext
from math import gcd
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
}

def is_on_curve(x: int, y: int, a: int, b: int, p: int) -> bool:
    """Check y^2 = x^3 + a*x + b (mod p)."""
    return (y * y - (x ** 3 + a * x + b)) % p == 0

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

_TOPIC_GENERATORS = {
    "legendre": _generate_legendre_problem,
    "crt": _generate_crt_problem,
}

def generate_problem(topic: str, rng: random.Random) -> dict:
    if topic == "mixed":
        topic = rng.choice(list(_TOPIC_GENERATORS.keys()))
    if topic not in _TOPIC_GENERATORS:
        raise ValueError(f"Unknown challenge topic: {topic}")
    return _TOPIC_GENERATORS[topic](rng)

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
    ec.add_argument("--curve", choices=sorted(CURVES.keys()),
                     help="Named curve (provides a, b, p, and generator G)")
    ec.add_argument("--x1", type=int, help="Defaults to G's x if --curve is given")
    ec.add_argument("--y1", type=int, help="Defaults to G's y if --curve is given")
    ec.add_argument("--x2", type=int, help="Defaults to (x1,y1) if omitted (i.e. doubling)")
    ec.add_argument("--y2", type=int, help="Defaults to (x1,y1) if omitted (i.e. doubling)")
    ec.add_argument("--a", type=int, help="Curve parameter a (overrides --curve if given)")
    ec.add_argument("--p", type=int, help="Prime modulus (overrides --curve if given)")

    # Chinese Remainder Theorem
    crt = subparsers.add_parser("crt", help="Solve a system of congruences x = r (mod m)")
    crt.add_argument("--r", type=int, nargs="+", required=True, help="Remainders, e.g. --r 2 3 2")
    crt.add_argument("--m", type=int, nargs="+", required=True, help="Moduli, e.g. --m 3 5 7")

    # Continued fractions
    cf_cmd = subparsers.add_parser("cf", help="Continued fraction expansion and convergents")
    cf_cmd.add_argument("--num", type=int, required=True, help="Numerator")
    cf_cmd.add_argument("--den", type=int, required=True, help="Denominator")
    cf_cmd.add_argument("--terms", type=int, default=100, help="Max terms to expand")

    # Challenge / quiz mode
    challenge = subparsers.add_parser("challenge", help="Practice problems (Legendre, CRT)")
    challenge.add_argument("--topic", choices=["legendre", "crt", "mixed"], default="mixed")
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
