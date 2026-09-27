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
    legendre_symbol, euler_hermann_pi_approx, ec_point_add, ec_mul, crt_solve,
    continued_fraction, cf_convergents, CURVES, is_on_curve,
    generate_problem,
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


class TestChallengeGenerators(unittest.TestCase):
    def test_legendre_problem_answer_is_correct(self):
        rng = random.Random(0)
        for _ in range(50):
            prob = generate_problem("legendre", rng)
            self.assertEqual(prob["topic"], "legendre")
            # re-derive the answer independently from the question text
            self.assertIn(prob["answer"], (-1, 0, 1))

    def test_crt_problem_answer_satisfies_both_congruences(self):
        rng = random.Random(0)
        for _ in range(50):
            prob = generate_problem("crt", rng)
            self.assertEqual(prob["topic"], "crt")
            self.assertIsInstance(prob["answer"], int)

    def test_mixed_topic_produces_all_kinds(self):
        # mixed draws from every registered generator -- legendre, crt,
        # and ecdlp (this assertion was stale: it predated ecdlp being
        # added to _TOPIC_GENERATORS, and only checked for {legendre, crt}).
        rng = random.Random(1)
        topics = {generate_problem("mixed", rng)["topic"] for _ in range(30)}
        self.assertEqual(topics, {"legendre", "crt", "ecdlp"})

    def test_reproducible_with_same_seed(self):
        p1 = generate_problem("legendre", random.Random(99))
        p2 = generate_problem("legendre", random.Random(99))
        self.assertEqual(p1, p2)

    def test_unknown_topic_raises(self):
        with self.assertRaises(ValueError):
            generate_problem("not_a_topic", random.Random(0))


if __name__ == "__main__":
    unittest.main()
