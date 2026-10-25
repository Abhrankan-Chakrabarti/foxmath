import unittest
from foxmath import (
    lucas_lehmer,
    mersenne_number,
    mersenne,
    mersenne_primes,
    perfect_number,
    is_perfect,
    rsa_keygen,
    rsa_encrypt,
    rsa_crack,
)

class TestLucasLehmer(unittest.TestCase):
    def test_known_mersenne_primes(self):
        known_exponents = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127]
        for p in known_exponents:
            with self.subTest(p=p):
                self.assertTrue(lucas_lehmer(p))

    def test_known_mersenne_composites(self):
        # 11: 2047 = 23 * 89
        # 23: 8388607 = 47 * 178481
        # 29: 536870911 = 233 * 1103 * 2089
        composite_exponents = [4, 6, 8, 9, 10, 11, 23, 29, 37, 41, 43]
        for p in composite_exponents:
            with self.subTest(p=p):
                self.assertFalse(lucas_lehmer(p))

    def test_invalid_exponents(self):
        with self.assertRaises(ValueError):
            lucas_lehmer(1)
        with self.assertRaises(ValueError):
            lucas_lehmer(0)
        with self.assertRaises(ValueError):
            lucas_lehmer(-5)
        with self.assertRaises(TypeError):
            lucas_lehmer("7")


class TestMersenneNumber(unittest.TestCase):
    def test_values(self):
        self.assertEqual(mersenne_number(1), 1)
        self.assertEqual(mersenne_number(2), 3)
        self.assertEqual(mersenne_number(3), 7)
        self.assertEqual(mersenne_number(5), 31)
        self.assertEqual(mersenne_number(7), 127)
        self.assertEqual(mersenne_number(13), 8191)

    def test_invalid_input(self):
        with self.assertRaises(ValueError):
            mersenne_number(0)
        with self.assertRaises(ValueError):
            mersenne_number(-1)
        with self.assertRaises(TypeError):
            mersenne_number("5")


class TestMersenneDict(unittest.TestCase):
    def test_prime_mersenne(self):
        info = mersenne(7)
        self.assertEqual(info["p"], 7)
        self.assertEqual(info["mp"], 127)
        self.assertTrue(info["p_is_prime"])
        self.assertTrue(info["is_prime"])
        self.assertEqual(info["bits"], 7)
        self.assertEqual(info["digits"], 3)
        self.assertEqual(info["perfect_number"], 8128)

    def test_composite_mersenne_with_prime_exponent(self):
        info = mersenne(11)
        self.assertEqual(info["p"], 11)
        self.assertEqual(info["mp"], 2047)
        self.assertTrue(info["p_is_prime"])
        self.assertFalse(info["is_prime"])
        self.assertNotIn("perfect_number", info)

    def test_composite_mersenne_with_composite_exponent(self):
        info = mersenne(4)
        self.assertEqual(info["p"], 4)
        self.assertEqual(info["mp"], 15)
        self.assertFalse(info["p_is_prime"])
        self.assertFalse(info["is_prime"])

    def test_invalid_input(self):
        with self.assertRaises(ValueError):
            mersenne(1)
        with self.assertRaises(TypeError):
            mersenne(3.14)


class TestMersennePrimesList(unittest.TestCase):
    def test_list_up_to_31(self):
        self.assertEqual(mersenne_primes(31), [2, 3, 5, 7, 13, 17, 19, 31])

    def test_list_edge_cases(self):
        self.assertEqual(mersenne_primes(1), [])
        self.assertEqual(mersenne_primes(0), [])
        self.assertEqual(mersenne_primes(-10), [])
        with self.assertRaises(TypeError):
            mersenne_primes("31")


class TestPerfectNumbers(unittest.TestCase):
    def test_is_perfect(self):
        self.assertTrue(is_perfect(6))
        self.assertTrue(is_perfect(28))
        self.assertTrue(is_perfect(496))
        self.assertTrue(is_perfect(8128))
        self.assertTrue(is_perfect(33550336))

    def test_is_perfect_non_perfect(self):
        self.assertFalse(is_perfect(1))
        self.assertFalse(is_perfect(2))
        self.assertFalse(is_perfect(12))
        self.assertFalse(is_perfect(100))
        self.assertFalse(is_perfect(8192))
        self.assertFalse(is_perfect(0))
        self.assertFalse(is_perfect(-6))

    def test_is_perfect_type_error(self):
        with self.assertRaises(TypeError):
            is_perfect("28")

    def test_perfect_number_generator(self):
        self.assertEqual(perfect_number(2), 6)
        self.assertEqual(perfect_number(3), 28)
        self.assertEqual(perfect_number(5), 496)
        self.assertEqual(perfect_number(7), 8128)
        self.assertEqual(perfect_number(13), 33550336)
        self.assertEqual(perfect_number(17), 8589869056)
        self.assertEqual(perfect_number(19), 137438691328)

    def test_perfect_number_invalid(self):
        with self.assertRaises(ValueError):
            perfect_number(11)  # 2^11 - 1 = 2047 is composite
        with self.assertRaises(ValueError):
            perfect_number(4)   # exponent is composite
        with self.assertRaises(ValueError):
            perfect_number(1)   # exponent < 2
        with self.assertRaises(TypeError):
            perfect_number("3")


class TestRsaCrackMethods(unittest.TestCase):
    def test_crack_fermat_close_primes(self):
        # p = 61, q = 53 (difference of 8, very close)
        n = 61 * 53
        e = 17
        c = rsa_encrypt(42, e, n)
        cracked = rsa_crack(n, e, c, method="fermat")
        self.assertEqual(cracked["p"], 61)
        self.assertEqual(cracked["q"], 53)
        self.assertEqual(cracked["m"], 42)
        self.assertEqual(cracked["method"], "fermat")

    def test_crack_auto_uses_fermat_for_close_primes(self):
        n = 61 * 53
        e = 17
        cracked = rsa_crack(n, e, method="auto")
        self.assertEqual(cracked["method"], "fermat")

    def test_crack_pollard_explicit(self):
        n = 61 * 53
        e = 17
        cracked = rsa_crack(n, e, method="pollard")
        self.assertEqual(cracked["method"], "pollard")
        self.assertEqual(cracked["p"], 61)
        self.assertEqual(cracked["q"], 53)

    def test_invalid_method(self):
        with self.assertRaises(ValueError):
            rsa_crack(3233, 17, method="unknown_method")


if __name__ == "__main__":
    unittest.main()
