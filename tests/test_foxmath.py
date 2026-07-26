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

from foxmath.foxmath import legendre_symbol, catalan_pi_approx, ec_point_add

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


class TestCatalanPiApprox(unittest.TestCase):
    def test_converges_to_pi(self):
        # Regression check: this used to converge to 16*arctan(1/2) =~ 7.418,
        # not pi, due to a double-counted scaling factor.
        approx = catalan_pi_approx(100)
        self.assertLess(abs(approx - PI_REFERENCE), Decimal("1e-40"))

    def test_more_terms_is_more_accurate(self):
        err_10 = abs(catalan_pi_approx(10) - PI_REFERENCE)
        err_50 = abs(catalan_pi_approx(50) - PI_REFERENCE)
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


if __name__ == "__main__":
    unittest.main()
