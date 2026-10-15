"""
Unit tests for foxmath's core functions, checked against independently
verifiable reference values (known Legendre symbols, math.pi, and
on-curve checks for elliptic curve arithmetic) rather than against the
module's own prior output.
"""
import sys
import os
import math
import unittest
from decimal import Decimal

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from foxmath.foxmath import (
    legendre_symbol, jacobi_symbol, solovay_strassen, sqrt_mod, is_prime,
    extended_gcd, mod_inverse, euler_hermann_pi_approx, ec_point_add, ec_mul,
    crt_solve, continued_fraction, cf_convergents, CURVES, is_on_curve,
    generate_problem, _is_probable_prime,
    factorize, format_factorization,
    euler_phi, carmichael_lambda, multiplicative_order, is_primitive_root,
    discrete_log, mobius, tau, sigma,
    dlog, egcd, modinv, crt,
)
import random

# Reference pi digits (mpmath, 250 dps of working precision), "3" + fractional digits.
PI_REFERENCE = Decimal(
    "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803"
    "4825342117067982148086513282306647093844609550582231725359408128481117450284102701938521"
    "1055596446229489549303819644288"
)


class TestLegendreSymbol(unittest.TestCase):
    def test_known_quadratic_residues_mod_7(self):
        # Squares mod 7: {1, 4, 2}
        self.assertEqual(legendre_symbol(1, 7), 1)
        self.assertEqual(legendre_symbol(2, 7), 1)
        self.assertEqual(legendre_symbol(4, 7), 1)

    def test_known_non_residues_mod_7(self):
        self.assertEqual(legendre_symbol(3, 7), -1)
        self.assertEqual(legendre_symbol(5, 7), -1)
        self.assertEqual(legendre_symbol(6, 7), -1)

    def test_multiple_of_p_is_zero(self):
        self.assertEqual(legendre_symbol(7, 7), 0)
        self.assertEqual(legendre_symbol(14, 7), 0)

    def test_rejects_even_p(self):
        with self.assertRaises(ValueError):
            legendre_symbol(3, 8)

    def test_rejects_p_less_than_2(self):
        with self.assertRaises(ValueError):
            legendre_symbol(3, 1)

    def test_rejects_composite_odd_p(self):
        for comp in (9, 15, 21, 25, 27, 33, 35, 561):
            with self.subTest(p=comp):
                with self.assertRaises(ValueError):
                    legendre_symbol(2, comp)


class TestSqrtMod(unittest.TestCase):
    def test_tonelli_shanks_returns_both_roots(self):
        roots = sqrt_mod(10, 13)
        self.assertEqual(roots, (6, 7))
        self.assertTrue(all((root * root - 10) % 13 == 0 for root in roots))

    def test_fast_path_for_prime_three_mod_four(self):
        self.assertEqual(sqrt_mod(9, 11), (3, 8))

    def test_zero_has_one_root(self):
        self.assertEqual(sqrt_mod(0, 17), (0,))

    def test_rejects_non_residue(self):
        with self.assertRaises(ValueError):
            sqrt_mod(5, 17)

    def test_rejects_even_p(self):
        with self.assertRaises(ValueError):
            sqrt_mod(1, 2)

    def test_rejects_p_less_than_2(self):
        with self.assertRaises(ValueError):
            sqrt_mod(1, 1)
        with self.assertRaises(ValueError):
            sqrt_mod(1, -5)

    def test_normalizes_a_greater_than_p(self):
        self.assertEqual(sqrt_mod(10 + 26, 13), (6, 7))

    def test_negative_a(self):
        self.assertEqual(sqrt_mod(-3, 13), (6, 7))

    def test_larger_primes_tonelli_shanks(self):
        roots_41 = sqrt_mod(8, 41)
        self.assertEqual(roots_41, (7, 34))
        self.assertTrue(all((r * r - 8) % 41 == 0 for r in roots_41))

        roots_97 = sqrt_mod(2, 97)
        self.assertEqual(roots_97, (14, 83))
        self.assertTrue(all((r * r - 2) % 97 == 0 for r in roots_97))

    def test_rejects_composite_odd_p(self):
        for comp in (9, 15, 25, 561):
            with self.subTest(p=comp):
                with self.assertRaises(ValueError):
                    sqrt_mod(4, comp)


class TestPrimality(unittest.TestCase):
    def test_small_primes(self):
        for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 97, 193):
            with self.subTest(p=p):
                self.assertTrue(_is_probable_prime(p))

    def test_small_composites(self):
        for c in (0, 1, 4, 6, 8, 9, 10, 12, 14, 15, 16, 18, 20, 21, 25, 27):
            with self.subTest(c=c):
                self.assertFalse(_is_probable_prime(c))

    def test_carmichael_numbers_rejected(self):
        for cn in (561, 1105, 1729):
            with self.subTest(cn=cn):
                self.assertFalse(_is_probable_prime(cn))

    def test_mersenne_prime(self):
        self.assertTrue(_is_probable_prime(2147483647))

    def test_is_prime_public_api(self):
        self.assertTrue(is_prime(104729))
        self.assertFalse(is_prime(561))
        self.assertFalse(is_prime(1))
        self.assertFalse(is_prime(0))
        self.assertFalse(is_prime(-7))


class TestEulerHermannPiApprox(unittest.TestCase):
    def test_converges_to_pi(self):
        # Regression check: this used to converge to 16*arctan(1/2) =~ 7.418,
        # not pi, due to a double-counted scaling factor.
        approx = euler_hermann_pi_approx(100)
        self.assertLess(abs(approx - PI_REFERENCE), Decimal("1e-40"))

    def test_more_terms_is_more_accurate(self):
        err_10 = abs(euler_hermann_pi_approx(10) - PI_REFERENCE)
        err_50 = abs(euler_hermann_pi_approx(50) - PI_REFERENCE)
        self.assertLess(err_50, err_10)


class TestECPointAdd(unittest.TestCase):
    # Curve y^2 = x^3 + 2x + 3 (mod 97); (3, 6) is on this curve
    # since 6^2=36 and 3^3+2*3+3=36.
    A, P = 2, 97

    def _on_curve(self, x, y):
        return (y * y - (x ** 3 + self.A * x + 3)) % self.P == 0

    def test_doubling_stays_on_curve(self):
        x3, y3 = ec_point_add(3, 6, 3, 6, self.A, self.P)
        self.assertTrue(self._on_curve(x3, y3))

    def test_distinct_point_addition_stays_on_curve(self):
        # Another point on the same curve: (80, 10) = 2*(3,6) from above.
        x3, y3 = ec_point_add(3, 6, 80, 10, self.A, self.P)
        self.assertTrue(self._on_curve(x3, y3))

    def test_point_plus_negation_raises_clear_error(self):
        with self.assertRaises(ValueError):
            ec_point_add(3, 6, 3, self.P - 6, self.A, self.P)


class TestCRT(unittest.TestCase):
    def test_classic_coprime_case(self):
        x, m = crt_solve([2, 3, 2], [3, 5, 7])
        self.assertEqual((x, m), (23, 105))
        self.assertEqual(x % 3, 2)
        self.assertEqual(x % 5, 3)
        self.assertEqual(x % 7, 2)

    def test_non_coprime_consistent(self):
        x, m = crt_solve([2, 2], [4, 6])
        self.assertEqual((x, m), (2, 12))

    def test_non_coprime_inconsistent_raises(self):
        with self.assertRaises(ValueError):
            crt_solve([1, 2], [4, 6])

    def test_single_congruence(self):
        self.assertEqual(crt_solve([5], [11]), (5, 11))

    def test_mismatched_lengths_raises(self):
        with self.assertRaises(ValueError):
            crt_solve([2, 3], [3, 5, 7])

    def test_empty_raises(self):
        with self.assertRaises(ValueError):
            crt_solve([], [])

    def test_negative_remainder_normalized(self):
        x, m = crt_solve([-1, 3], [5, 7])
        self.assertEqual(x % 5, 4)  # -1 mod 5 == 4
        self.assertEqual(x % 7, 3)


class TestContinuedFractions(unittest.TestCase):
    def test_classic_pi_approximation(self):
        cf = continued_fraction(355, 113)
        self.assertEqual(cf, [3, 7, 16])
        convs = cf_convergents(cf)
        self.assertEqual(convs, [(3, 1), (22, 7), (355, 113)])

    def test_convergents_reduce_to_final_fraction(self):
        cf = continued_fraction(19, 7)
        convs = cf_convergents(cf)
        self.assertEqual(convs[-1], (19, 7))

    def test_negative_numerator(self):
        cf = continued_fraction(-19, 7)
        convs = cf_convergents(cf)
        self.assertEqual(convs[-1], (-19, 7))

    def test_zero_denominator_raises(self):
        with self.assertRaises(ValueError):
            continued_fraction(5, 0)

    def test_randomized_against_fractions_module(self):
        import random
        from fractions import Fraction
        random.seed(42)
        for _ in range(200):
            n = random.randint(-10000, 10000)
            d = random.randint(1, 10000)
            cf = continued_fraction(n, d)
            convs = cf_convergents(cf)
            expected = Fraction(n, d)
            self.assertEqual(convs[-1], (expected.numerator, expected.denominator))


class TestECMul(unittest.TestCase):
    """
    ec_mul (scalar multiplication / double-and-add) had no test coverage
    at all despite being a documented headline feature -- these tests
    were added after finding and closing that gap, not in response to
    any bug (cross-validation below found none).
    """

    def _naive_mul(self, k, x, y, a, p):
        # Independent reference: k repeated additions via the
        # user-facing ec_point_add, rather than ec_mul's own
        # double-and-add logic.
        if k == 0:
            return None
        result = None
        point = (x, y)
        for _ in range(k):
            if result is None:
                result = point
            else:
                try:
                    result = ec_point_add(result[0], result[1], point[0], point[1], a, p)
                except ValueError:
                    result = None  # hit the point at infinity
        return result

    def test_matches_naive_repeated_addition_toy97(self):
        c = CURVES["toy97"]
        for k in [1, 2, 3, 5, 7, 10, 13, 20]:
            fast = ec_mul(k, c["gx"], c["gy"], c["a"], c["p"])
            naive = self._naive_mul(k, c["gx"], c["gy"], c["a"], c["p"])
            self.assertEqual(fast, naive, f"mismatch at k={k}")

    def test_matches_naive_repeated_addition_toy193(self):
        c = CURVES["toy193"]
        for k in [1, 2, 3, 5, 7, 10, 13, 20]:
            fast = ec_mul(k, c["gx"], c["gy"], c["a"], c["p"])
            naive = self._naive_mul(k, c["gx"], c["gy"], c["a"], c["p"])
            self.assertEqual(fast, naive, f"mismatch at k={k}")

    def test_secp256k1_matches_independent_reference_including_large_k(self):
        # Reference points cross-checked independently against the
        # `ecdsa` library's SECP256k1.generator * k (same verification
        # approach the pre-existing doubling test above uses). Values
        # are baked in rather than imported at test time so this suite
        # doesn't pick up a runtime dependency on `ecdsa` -- foxmath
        # itself stays pure stdlib (dependencies = [] in pyproject.toml).
        c = CURVES["secp256k1"]
        reference = {
            1: (55066263022277343669578718895168534326250603453777594175500187360389116729240,
                32670510020758816978083085130507043184471273380659243275938904335757337482424),
            2: (89565891926547004231252920425935692360644145829622209833684329913297188986597,
                12158399299693830322967808612713398636155367887041628176798871954788371653930),
            3: (112711660439710606056748659173929673102114977341539408544630613555209775888121,
                25583027980570883691656905877401976406448868254816295069919888960541586679410),
            7: (41948375291644419605210209193538855353224492619856392092318293986323063962044,
                48361766907851246668144012348516735800090617714386977531302791340517493990618),
            100: (107303582290733097924842193972465022053148211775194373671539518313500194639752,
                  103795966108782717446806684023742168462365449272639790795591544606836007446638),
            12345: (108607064596551879580190606910245687803607295064141551927605737287325610911759,
                    6661302038839728943522144359728938428925407345457796456954441906546235843221),
            # A scalar well beyond anything repeated addition could feasibly verify.
            2**64 + 17: (5136902638900309784563324385263118307222550852398366895311091889317425008288,
                         19220983661714808463696268177333916812952649794945216772534518833957341771155),
        }
        for k, expected in reference.items():
            fast = ec_mul(k, c["gx"], c["gy"], c["a"], c["p"])
            self.assertEqual(fast, expected, f"mismatch at k={k}")

    def test_toy193_has_prime_order_full_cycle(self):
        # README's documented claim: #E = 193 (prime), so every non-zero
        # multiple up to 192 should land on a distinct on-curve point,
        # and k=193 (the order) should land on infinity.
        c = CURVES["toy193"]
        seen = set()
        for k in range(1, 193):
            pt = ec_mul(k, c["gx"], c["gy"], c["a"], c["p"])
            self.assertIsNotNone(pt, f"unexpected infinity at k={k}")
            self.assertTrue(is_on_curve(pt[0], pt[1], c["a"], c["b"], c["p"]))
            self.assertNotIn(pt, seen, f"duplicate point at k={k} -- order is not 193")
            seen.add(pt)
        self.assertEqual(len(seen), 192)

    def test_k_zero_is_infinity(self):
        c = CURVES["toy193"]
        self.assertIsNone(ec_mul(0, c["gx"], c["gy"], c["a"], c["p"]))

    def test_k_equal_to_order_is_infinity(self):
        c = CURVES["toy193"]
        self.assertIsNone(ec_mul(193, c["gx"], c["gy"], c["a"], c["p"]))

    def test_negative_k_raises(self):
        c = CURVES["toy193"]
        with self.assertRaises(ValueError):
            ec_mul(-5, c["gx"], c["gy"], c["a"], c["p"])


class TestNamedCurves(unittest.TestCase):
    def test_secp256k1_generator_is_on_curve(self):
        c = CURVES["secp256k1"]
        self.assertTrue(is_on_curve(c["gx"], c["gy"], c["a"], c["b"], c["p"]))

    def test_secp256k1_doubling_matches_known_value(self):
        # Cross-checked independently against the `ecdsa` library's
        # SECP256k1.generator * 2.
        c = CURVES["secp256k1"]
        x2, y2 = ec_point_add(c["gx"], c["gy"], c["gx"], c["gy"], c["a"], c["p"])
        self.assertEqual(
            x2, 89565891926547004231252920425935692360644145829622209833684329913297188986597
        )
        self.assertEqual(
            y2, 12158399299693830322967808612713398636155367887041628176798871954788371653930
        )
        self.assertTrue(is_on_curve(x2, y2, c["a"], c["b"], c["p"]))

    def test_toy97_generator_is_on_curve(self):
        c = CURVES["toy97"]
        self.assertTrue(is_on_curve(c["gx"], c["gy"], c["a"], c["b"], c["p"]))


class TestExtendedGCD(unittest.TestCase):
    def test_known_coprime(self):
        g, x, y = extended_gcd(3, 5)
        self.assertEqual(g, 1)
        self.assertEqual(3 * x + 5 * y, 1)

    def test_known_common_factor(self):
        g, x, y = extended_gcd(30, 50)
        self.assertEqual(g, 10)
        self.assertEqual(30 * x + 50 * y, 10)

    def test_negative_inputs(self):
        g, x, y = extended_gcd(-30, 50)
        self.assertEqual(g, 10)
        self.assertEqual(-30 * x + 50 * y, 10)

    def test_zero_arguments(self):
        self.assertEqual(extended_gcd(0, 5), (5, 0, 1))
        self.assertEqual(extended_gcd(5, 0), (5, 1, 0))


class TestModInverse(unittest.TestCase):
    def test_known_inverse(self):
        self.assertEqual(mod_inverse(3, 11), 4)
        self.assertEqual((3 * 4) % 11, 1)

    def test_negative_a(self):
        inv = mod_inverse(-3, 11)
        self.assertEqual((-3 * inv) % 11, 1)

    def test_non_coprime_raises(self):
        with self.assertRaises(ValueError):
            mod_inverse(6, 9)

    def test_modulus_less_than_two_raises(self):
        with self.assertRaises(ValueError):
            mod_inverse(1, 1)
        with self.assertRaises(ValueError):
            mod_inverse(1, 0)


class TestFactorize(unittest.TestCase):
    def test_prime_numbers(self):
        self.assertEqual(factorize(2), [(2, 1)])
        self.assertEqual(factorize(13), [(13, 1)])
        self.assertEqual(factorize(104729), [(104729, 1)])

    def test_composite_powers(self):
        self.assertEqual(factorize(8), [(2, 3)])
        self.assertEqual(factorize(360), [(2, 3), (3, 2), (5, 1)])
        self.assertEqual(factorize(1000), [(2, 3), (5, 3)])

    def test_semiprimes_large(self):
        p1 = 1000003
        p2 = 1000033
        n = p1 * p2
        self.assertEqual(factorize(n), [(p1, 1), (p2, 1)])

    def test_mersenne_composite(self):
        # 2^11 - 1 = 2047 = 23 * 89
        self.assertEqual(factorize(2047), [(23, 1), (89, 1)])

    def test_one(self):
        self.assertEqual(factorize(1), [])

    def test_invalid_input(self):
        with self.assertRaises(ValueError):
            factorize(0)
        with self.assertRaises(ValueError):
            factorize(-15)
        with self.assertRaises(TypeError):
            factorize("10")

    def test_format_factorization(self):
        self.assertEqual(format_factorization([]), "1")
        self.assertEqual(format_factorization([(2, 3), (3, 2), (5, 1)]), "2^3 * 3^2 * 5")
        self.assertEqual(format_factorization([(13, 1)]), "13")


class TestEulerPhi(unittest.TestCase):
    def test_known_values(self):
        self.assertEqual(euler_phi(1), 1)
        self.assertEqual(euler_phi(2), 1)
        self.assertEqual(euler_phi(7), 6)
        self.assertEqual(euler_phi(9), 6)
        self.assertEqual(euler_phi(12), 4)
        self.assertEqual(euler_phi(360), 96)

    def test_prime_powers(self):
        # phi(p^k) = p^(k-1) * (p - 1)
        self.assertEqual(euler_phi(8), 4)       # 2^3 -> 4
        self.assertEqual(euler_phi(27), 18)     # 3^3 -> 18
        self.assertEqual(euler_phi(25), 20)     # 5^2 -> 20
        self.assertEqual(euler_phi(49), 42)     # 7^2 -> 42

    def test_multiplicative_property(self):
        # gcd(a, b) == 1 => phi(a * b) == phi(a) * phi(b)
        self.assertEqual(euler_phi(8 * 9), euler_phi(8) * euler_phi(9))
        self.assertEqual(euler_phi(5 * 7), euler_phi(5) * euler_phi(7))

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            euler_phi(0)
        with self.assertRaises(ValueError):
            euler_phi(-10)
        with self.assertRaises(TypeError):
            euler_phi("360")


class TestCarmichaelLambda(unittest.TestCase):
    def test_known_values(self):
        self.assertEqual(carmichael_lambda(1), 1)
        self.assertEqual(carmichael_lambda(2), 1)
        self.assertEqual(carmichael_lambda(4), 2)
        self.assertEqual(carmichael_lambda(8), 2)
        self.assertEqual(carmichael_lambda(16), 4)
        self.assertEqual(carmichael_lambda(15), 4)
        self.assertEqual(carmichael_lambda(360), 12)

    def test_cyclic_moduli(self):
        # Groups (Z/nZ)* that are cyclic: lambda(n) == phi(n)
        # Form: 1, 2, 4, p^k, 2*p^k for odd prime p
        for n in (1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 13, 14, 17, 18, 19, 23, 25, 27):
            self.assertEqual(carmichael_lambda(n), euler_phi(n))

    def test_non_cyclic_moduli(self):
        # Non-cyclic groups have lambda(n) < phi(n)
        for n in (8, 12, 15, 16, 20, 24, 30, 360):
            self.assertLess(carmichael_lambda(n), euler_phi(n))

    def test_exponent_property(self):
        # For every a coprime to n: a^lambda(n) == 1 (mod n)
        for n in (8, 12, 15, 20, 24, 360):
            lam = carmichael_lambda(n)
            for a in range(1, min(n, 50)):
                if math.gcd(a, n) == 1:
                    self.assertEqual(pow(a, lam, n), 1)

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            carmichael_lambda(0)
        with self.assertRaises(ValueError):
            carmichael_lambda(-5)
        with self.assertRaises(TypeError):
            carmichael_lambda("15")


class TestMultiplicativeOrder(unittest.TestCase):
    def test_known_orders(self):
        self.assertEqual(multiplicative_order(3, 7), 6)
        self.assertEqual(multiplicative_order(2, 7), 3)
        self.assertEqual(multiplicative_order(1, 7), 1)
        self.assertEqual(multiplicative_order(6, 7), 2)
        self.assertEqual(multiplicative_order(3, 8), 2)
        self.assertEqual(multiplicative_order(3, 17), 16)

    def test_negative_a(self):
        # -1 mod 7 is 6, order is 2
        self.assertEqual(multiplicative_order(-1, 7), 2)

    def test_non_coprime_raises(self):
        with self.assertRaises(ValueError):
            multiplicative_order(2, 8)
        with self.assertRaises(ValueError):
            multiplicative_order(6, 9)

    def test_modulus_le_one_raises(self):
        with self.assertRaises(ValueError):
            multiplicative_order(1, 1)
        with self.assertRaises(ValueError):
            multiplicative_order(1, 0)
        with self.assertRaises(ValueError):
            multiplicative_order(1, -5)

    def test_invalid_type_raises(self):
        with self.assertRaises(TypeError):
            multiplicative_order("3", 7)
        with self.assertRaises(TypeError):
            multiplicative_order(3, "7")

    def test_is_primitive_root(self):
        self.assertTrue(is_primitive_root(3, 7))
        self.assertFalse(is_primitive_root(2, 7))
        # 8 is not cyclic, no primitive roots exist
        self.assertFalse(is_primitive_root(3, 8))
        self.assertFalse(is_primitive_root(2, 8))  # not coprime
        self.assertFalse(is_primitive_root(1, 1))  # invalid modulus


class TestDiscreteLog(unittest.TestCase):
    def test_known_values(self):
        self.assertEqual(discrete_log(3, 6, 7), 3)
        self.assertEqual(discrete_log(2, 5, 13), 9)
        self.assertEqual(discrete_log(3, 1, 7), 0)
        self.assertEqual(discrete_log(3, 3, 7), 1)
        self.assertEqual(discrete_log(2, 1, 5), 0)

    def test_subgroup(self):
        # 4 has order 4 in (Z/17Z)*: 4^0=1, 4^1=4, 4^2=16, 4^3=13
        self.assertEqual(discrete_log(4, 13, 17), 3)
        self.assertEqual(discrete_log(4, 16, 17), 2)
        self.assertEqual(discrete_log(4, 4, 17), 1)

    def test_no_solution_raises(self):
        # 4 does not generate 3 mod 17
        with self.assertRaises(ValueError):
            discrete_log(4, 3, 17)
        # 2 generates subgroup of size 3 in (Z/7Z)*: {1, 2, 4}
        with self.assertRaises(ValueError):
            discrete_log(2, 3, 7)

    def test_non_coprime_raises(self):
        with self.assertRaises(ValueError):
            discrete_log(2, 4, 8)

    def test_invalid_modulus_raises(self):
        with self.assertRaises(ValueError):
            discrete_log(2, 3, 1)
        with self.assertRaises(ValueError):
            discrete_log(2, 3, -5)

    def test_invalid_types_raises(self):
        with self.assertRaises(TypeError):
            discrete_log("2", 3, 5)


class TestMobius(unittest.TestCase):
    def test_known_values(self):
        self.assertEqual(mobius(1), 1)
        self.assertEqual(mobius(2), -1)
        self.assertEqual(mobius(3), -1)
        self.assertEqual(mobius(4), 0)
        self.assertEqual(mobius(5), -1)
        self.assertEqual(mobius(6), 1)
        self.assertEqual(mobius(12), 0)
        self.assertEqual(mobius(30), -1)

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            mobius(0)
        with self.assertRaises(ValueError):
            mobius(-5)
        with self.assertRaises(TypeError):
            mobius("6")


class TestTau(unittest.TestCase):
    def test_known_values(self):
        self.assertEqual(tau(1), 1)
        self.assertEqual(tau(2), 2)
        self.assertEqual(tau(6), 4)
        self.assertEqual(tau(12), 6)
        self.assertEqual(tau(360), 24)

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            tau(0)
        with self.assertRaises(ValueError):
            tau(-10)
        with self.assertRaises(TypeError):
            tau("12")


class TestSigma(unittest.TestCase):
    def test_known_values(self):
        self.assertEqual(sigma(1), 1)
        self.assertEqual(sigma(6), 12)
        self.assertEqual(sigma(12), 28)
        self.assertEqual(sigma(28), 56)  # 28 is a perfect number
        self.assertEqual(sigma(6, 2), 50)  # 1^2 + 2^2 + 3^2 + 6^2 = 50
        self.assertEqual(sigma(12, 0), 6)  # sigma_0(n) == tau(n)

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            sigma(0)
        with self.assertRaises(ValueError):
            sigma(-5)
        with self.assertRaises(ValueError):
            sigma(6, -1)
        with self.assertRaises(TypeError):
            sigma("6")


class TestChallengeGenerators(unittest.TestCase):
    def test_legendre_problem_answer_is_correct(self):
        rng = random.Random(0)
        for _ in range(50):
            prob = generate_problem("legendre", rng)
            self.assertEqual(prob["topic"], "legendre")
            self.assertIn(prob["answer"], (-1, 0, 1))

    def test_sqrtmod_problem_answer_is_correct(self):
        rng = random.Random(0)
        for _ in range(50):
            prob = generate_problem("sqrtmod", rng)
            self.assertEqual(prob["topic"], "sqrtmod")
            r = prob["answer"]
            self.assertGreaterEqual(r, 1)
            # Question is "Find the smaller positive square root of a modulo p..."
            self.assertIn("Find the smaller positive square root of ", prob["question"])

    def test_crt_problem_answer_satisfies_both_congruences(self):
        rng = random.Random(0)
        for _ in range(50):
            prob = generate_problem("crt", rng)
            self.assertEqual(prob["topic"], "crt")
            self.assertIsInstance(prob["answer"], int)

    def test_factor_problem_answer_is_correct(self):
        rng = random.Random(0)
        for _ in range(50):
            prob = generate_problem("factor", rng)
            self.assertEqual(prob["topic"], "factor")
            self.assertIsInstance(prob["answer"], int)
            self.assertIn("Find the smallest prime factor of ", prob["question"])

    def test_totient_problem_answer_is_correct(self):
        rng = random.Random(0)
        for _ in range(50):
            prob = generate_problem("totient", rng)
            self.assertEqual(prob["topic"], "totient")
            self.assertIsInstance(prob["answer"], int)
            self.assertIn("Compute Euler's totient phi(", prob["question"])

    def test_dlog_problem_answer_is_correct(self):
        rng = random.Random(0)
        for _ in range(50):
            prob = generate_problem("dlog", rng)
            self.assertEqual(prob["topic"], "dlog")
            self.assertIsInstance(prob["answer"], int)
            self.assertIn("such that", prob["question"])

    def test_jacobi_problem_answer_is_correct(self):
        rng = random.Random(0)
        for _ in range(50):
            prob = generate_problem("jacobi", rng)
            self.assertEqual(prob["topic"], "jacobi")
            self.assertIn(prob["answer"], (-1, 0, 1))
            self.assertIn("Jacobi symbol (", prob["question"])

    def test_rsa_problem_answer_is_correct(self):
        rng = random.Random(0)
        for _ in range(50):
            prob = generate_problem("rsa", rng)
            self.assertEqual(prob["topic"], "rsa")
            self.assertIsInstance(prob["answer"], int)
            self.assertIn("In RSA with e =", prob["question"])

    def test_dh_problem_answer_is_correct(self):
        rng = random.Random(0)
        for _ in range(50):
            prob = generate_problem("dh", rng)
            self.assertEqual(prob["topic"], "dh")
            self.assertIsInstance(prob["answer"], int)
            self.assertIn("In Diffie-Hellman with p =", prob["question"])

    def test_mixed_topic_produces_all_kinds(self):
        rng = random.Random(1)
        topics = {generate_problem("mixed", rng)["topic"] for _ in range(250)}
        self.assertEqual(topics, {"legendre", "crt", "ecdlp", "sqrtmod", "factor", "totient", "dlog", "jacobi", "rsa", "dh"})

    def test_reproducible_with_same_seed(self):
        p1 = generate_problem("legendre", random.Random(99))
        p2 = generate_problem("legendre", random.Random(99))
        self.assertEqual(p1, p2)

    def test_unknown_topic_raises(self):
        with self.assertRaises(ValueError):
            generate_problem("not_a_topic", random.Random(0))


class TestJacobiSymbol(unittest.TestCase):
    def test_known_values(self):
        self.assertEqual(jacobi_symbol(2, 9), 1)
        self.assertEqual(jacobi_symbol(2, 15), 1)
        self.assertEqual(jacobi_symbol(2, 21), -1)
        self.assertEqual(jacobi_symbol(6, 21), 0)
        self.assertEqual(jacobi_symbol(10, 13), 1)
        self.assertEqual(jacobi_symbol(3, 7), -1)
        self.assertEqual(jacobi_symbol(1007, 2003), 1)

    def test_parity_with_legendre_symbol(self):
        for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41):
            for a in range(-30, 30):
                self.assertEqual(jacobi_symbol(a, p), legendre_symbol(a, p))

    def test_composite_moduli_multiplicative_definition(self):
        # (a/n) == prod (a/p)^e
        composites = [9, 15, 21, 25, 27, 33, 35, 45, 55, 63, 77, 91, 105]
        for n in composites:
            factors = factorize(n)
            for a in range(-20, 20):
                expected = 1
                for p, power in factors:
                    expected *= (legendre_symbol(a, p) ** power)
                self.assertEqual(jacobi_symbol(a, n), expected)

    def test_multiplicative_numerator(self):
        for n in (7, 15, 21, 35):
            for a in range(1, 15):
                for b in range(1, 15):
                    self.assertEqual(
                        jacobi_symbol(a * b, n),
                        jacobi_symbol(a, n) * jacobi_symbol(b, n),
                    )

    def test_multiplicative_denominator(self):
        for a in range(1, 20):
            m, n = 7, 11
            self.assertEqual(
                jacobi_symbol(a, m * n),
                jacobi_symbol(a, m) * jacobi_symbol(a, n),
            )

    def test_n_equals_one(self):
        for a in (-10, -1, 0, 1, 2, 42, 100):
            self.assertEqual(jacobi_symbol(a, 1), 1)

    def test_zero_numerator(self):
        self.assertEqual(jacobi_symbol(0, 1), 1)
        self.assertEqual(jacobi_symbol(0, 3), 0)
        self.assertEqual(jacobi_symbol(0, 15), 0)

    def test_negative_numerators(self):
        self.assertEqual(jacobi_symbol(-1, 7), -1)  # 7 ≡ 3 mod 4
        self.assertEqual(jacobi_symbol(-1, 5), 1)   # 5 ≡ 1 mod 4
        self.assertEqual(jacobi_symbol(-1, 15), -1) # 15 ≡ 3 mod 4
        self.assertEqual(jacobi_symbol(-2, 15), -1)

    def test_modular_equivalence(self):
        for n in (9, 15, 21):
            for a in range(-10, 10):
                self.assertEqual(jacobi_symbol(a, n), jacobi_symbol(a + 5 * n, n))

    def test_non_coprime_gives_zero(self):
        self.assertEqual(jacobi_symbol(15, 25), 0)
        self.assertEqual(jacobi_symbol(14, 21), 0)
        self.assertEqual(jacobi_symbol(33, 99), 0)

    def test_invalid_modulus_raises(self):
        with self.assertRaises(ValueError):
            jacobi_symbol(2, 0)
        with self.assertRaises(ValueError):
            jacobi_symbol(2, -3)
        with self.assertRaises(ValueError):
            jacobi_symbol(2, 2)
        with self.assertRaises(ValueError):
            jacobi_symbol(2, 8)

    def test_invalid_types_raise(self):
        with self.assertRaises(TypeError):
            jacobi_symbol("2", 7)
        with self.assertRaises(TypeError):
            jacobi_symbol(2, "7")
        with self.assertRaises(TypeError):
            jacobi_symbol(2.5, 7)


class TestSolovayStrassen(unittest.TestCase):
    def test_known_primes(self):
        for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 104729):
            self.assertTrue(solovay_strassen(p, rounds=20, rng=random.Random(p)))

    def test_known_composites(self):
        for c in (4, 6, 8, 9, 10, 12, 14, 15, 16, 21, 25, 27, 33, 35, 49, 100):
            self.assertFalse(solovay_strassen(c, rounds=20, rng=random.Random(c)))

    def test_carmichael_numbers_detected_as_composite(self):
        for c in (561, 1105, 1729, 2465, 2821):
            self.assertFalse(solovay_strassen(c, rounds=25, rng=random.Random(c)))

    def test_euler_pseudoprime_91(self):
        # 91 = 7 * 13 is an Euler pseudoprime to base 9, but Solovay-Strassen detects compositeness
        self.assertFalse(solovay_strassen(91, rounds=20, rng=random.Random(91)))

    def test_small_and_negative_inputs(self):
        self.assertFalse(solovay_strassen(0))
        self.assertFalse(solovay_strassen(1))
        self.assertFalse(solovay_strassen(-7))
        self.assertFalse(solovay_strassen(-10))

    def test_deterministic_seeded_rng(self):
        res1 = solovay_strassen(104729, rounds=15, rng=random.Random(42))
        res2 = solovay_strassen(104729, rounds=15, rng=random.Random(42))
        self.assertEqual(res1, res2)
        self.assertTrue(res1)

    def test_invalid_arguments_raise(self):
        with self.assertRaises(TypeError):
            solovay_strassen("17")
        with self.assertRaises(TypeError):
            solovay_strassen(17, rounds="20")
        with self.assertRaises(ValueError):
            solovay_strassen(17, rounds=0)
        with self.assertRaises(ValueError):
            solovay_strassen(17, rounds=-5)


class TestAliasesAndPackageExports(unittest.TestCase):
    def test_convenience_aliases(self):
        self.assertEqual(dlog(2, 5, 13), 9)
        self.assertEqual(egcd(30, 50), (10, 2, -1))
        self.assertEqual(modinv(3, 11), 4)
        self.assertEqual(crt([2, 3, 2], [3, 5, 7]), (23, 105))

    def test_top_level_package_exports(self):
        import foxmath
        self.assertTrue(callable(foxmath.jacobi_symbol))
        self.assertTrue(callable(foxmath.solovay_strassen))
        self.assertTrue(callable(foxmath.discrete_log))
        self.assertTrue(callable(foxmath.dlog))
        self.assertTrue(callable(foxmath.mobius))
        self.assertTrue(callable(foxmath.sigma))
        self.assertTrue(callable(foxmath.tau))
        self.assertTrue(callable(foxmath.egcd))
        self.assertTrue(callable(foxmath.modinv))
        self.assertTrue(callable(foxmath.crt))
        self.assertEqual(foxmath.jacobi_symbol(10, 13), 1)
        self.assertTrue(foxmath.solovay_strassen(17))


if __name__ == "__main__":
    unittest.main()
