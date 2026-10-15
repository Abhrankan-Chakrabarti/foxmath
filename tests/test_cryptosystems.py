"""
Unit tests for FoxMath's textbook cryptosystems suite (RSA & Diffie-Hellman).
"""
import os
import sys
import unittest
import random

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src")))

from foxmath.foxmath import (
    rsa_keygen, rsa_encrypt, rsa_decrypt, rsa_crack,
    dh_keygen, dh_derive, dh_crack,
    is_prime, gcd,
)
import foxmath


class TestRsaKeygen(unittest.TestCase):
    def test_default_keygen_32_bit(self):
        keys = rsa_keygen(bits=32, seed=100)
        self.assertIn("p", keys)
        self.assertIn("q", keys)
        self.assertIn("n", keys)
        self.assertIn("e", keys)
        self.assertIn("d", keys)
        self.assertEqual(keys["n"], keys["p"] * keys["q"])
        self.assertTrue(is_prime(keys["p"]))
        self.assertTrue(is_prime(keys["q"]))
        self.assertNotEqual(keys["p"], keys["q"])
        self.assertEqual((keys["e"] * keys["d"]) % keys["lambda"], 1)

    def test_configurable_bit_lengths(self):
        for b in (16, 24, 32, 40):
            keys = rsa_keygen(bits=b, seed=b)
            self.assertEqual(keys["n"], keys["p"] * keys["q"])
            self.assertGreaterEqual(keys["n"].bit_length(), b - 2)

    def test_deterministic_seed(self):
        k1 = rsa_keygen(bits=32, seed=42)
        k2 = rsa_keygen(bits=32, seed=42)
        self.assertEqual(k1, k2)
        k3 = rsa_keygen(bits=32, seed=99)
        self.assertNotEqual(k1["n"], k3["n"])

    def test_crt_components(self):
        keys = rsa_keygen(bits=32, seed=123)
        p, q, d = keys["p"], keys["q"], keys["d"]
        self.assertEqual(keys["dp"], d % (p - 1))
        self.assertEqual(keys["dq"], d % (q - 1))
        self.assertEqual((keys["qinv"] * q) % p, 1)

    def test_invalid_arguments_raise(self):
        with self.assertRaises(ValueError):
            rsa_keygen(bits=7)
        with self.assertRaises(ValueError):
            rsa_keygen(bits=32, e=2)  # even e
        with self.assertRaises(ValueError):
            rsa_keygen(bits=32, e=1)  # e < 3
        with self.assertRaises(TypeError):
            rsa_keygen(bits="32")
        with self.assertRaises(TypeError):
            rsa_keygen(bits=32, e="65537")


class TestRsaEncryptDecrypt(unittest.TestCase):
    def setUp(self):
        self.keys = rsa_keygen(bits=32, seed=777)

    def test_encrypt_decrypt_roundtrip(self):
        e, d, n = self.keys["e"], self.keys["d"], self.keys["n"]
        for m in (0, 1, 2, 42, 12345, n - 2, n - 1):
            c = rsa_encrypt(m, e, n)
            dec = rsa_decrypt(c, d, n)
            self.assertEqual(dec, m)

    def test_crt_acceleration_matches_standard(self):
        e, d, n = self.keys["e"], self.keys["d"], self.keys["n"]
        p, q = self.keys["p"], self.keys["q"]
        m = 987654
        c = rsa_encrypt(m, e, n)
        dec_standard = rsa_decrypt(c, d, n)
        dec_crt = rsa_decrypt(c, d, n, p=p, q=q)
        self.assertEqual(dec_standard, m)
        self.assertEqual(dec_crt, m)

    def test_invalid_message_bounds_raise(self):
        e, n = self.keys["e"], self.keys["n"]
        with self.assertRaises(ValueError):
            rsa_encrypt(-1, e, n)
        with self.assertRaises(ValueError):
            rsa_encrypt(n, e, n)
        with self.assertRaises(ValueError):
            rsa_encrypt(n + 10, e, n)

    def test_invalid_ciphertext_bounds_raise(self):
        d, n = self.keys["d"], self.keys["n"]
        with self.assertRaises(ValueError):
            rsa_decrypt(-1, d, n)
        with self.assertRaises(ValueError):
            rsa_decrypt(n, d, n)

    def test_invalid_factors_raise(self):
        c, d, n = 1234, self.keys["d"], self.keys["n"]
        with self.assertRaises(ValueError):
            rsa_decrypt(c, d, n, p=7, q=11)  # 7 * 11 != n

    def test_type_errors(self):
        with self.assertRaises(TypeError):
            rsa_encrypt("42", 65537, 100)
        with self.assertRaises(TypeError):
            rsa_decrypt("42", 123, 100)


class TestRsaCrack(unittest.TestCase):
    def test_crack_known_textbook_example(self):
        # p = 61, q = 53, n = 3233, e = 17, c = 2790 (m = 65)
        cracked = rsa_crack(3233, 17, c=2790)
        self.assertEqual(cracked["p"], 61)
        self.assertEqual(cracked["q"], 53)
        self.assertEqual(cracked["phi"], 3120)
        self.assertEqual(cracked["lambda"], 780)
        self.assertEqual(cracked["d"], 413)
        self.assertEqual(cracked["m"], 65)

    def test_crack_generated_32bit_key(self):
        keys = rsa_keygen(bits=32, seed=888)
        m = 543210
        c = rsa_encrypt(m, keys["e"], keys["n"])
        cracked = rsa_crack(keys["n"], keys["e"], c=c)
        self.assertEqual(cracked["p"], keys["p"])
        self.assertEqual(cracked["q"], keys["q"])
        self.assertEqual(cracked["m"], m)

    def test_crack_without_ciphertext(self):
        keys = rsa_keygen(bits=32, seed=999)
        cracked = rsa_crack(keys["n"], keys["e"])
        self.assertNotIn("c", cracked)
        self.assertNotIn("m", cracked)
        self.assertEqual(cracked["d"], keys["d"])

    def test_crack_invalid_modulus_raises(self):
        with self.assertRaises(ValueError):
            rsa_crack(17, 3)  # prime n
        with self.assertRaises(ValueError):
            rsa_crack(30, 3)  # 2 * 3 * 5 has 3 prime factors


class TestDhKeygen(unittest.TestCase):
    def test_default_keygen(self):
        p, g = 104729, 3
        keypair = dh_keygen(p, g, seed=42)
        self.assertEqual(keypair["p"], p)
        self.assertEqual(keypair["g"], g)
        self.assertGreaterEqual(keypair["private_key"], 2)
        self.assertLess(keypair["private_key"], p - 1)
        self.assertEqual(keypair["public_key"], pow(g, keypair["private_key"], p))

    def test_explicit_private_key(self):
        p, g = 104729, 3
        keypair = dh_keygen(p, g, private_key=1234)
        self.assertEqual(keypair["private_key"], 1234)
        self.assertEqual(keypair["public_key"], pow(3, 1234, 104729))

    def test_deterministic_seed(self):
        k1 = dh_keygen(104729, 3, seed=55)
        k2 = dh_keygen(104729, 3, seed=55)
        self.assertEqual(k1, k2)

    def test_invalid_parameters_raise(self):
        with self.assertRaises(ValueError):
            dh_keygen(561, 3)  # composite p
        with self.assertRaises(ValueError):
            dh_keygen(2, 3)    # p <= 2
        with self.assertRaises(ValueError):
            dh_keygen(17, 1)   # g <= 1
        with self.assertRaises(ValueError):
            dh_keygen(17, 17)  # g >= p
        with self.assertRaises(ValueError):
            dh_keygen(17, 3, private_key=16)  # priv >= p - 1


class TestDhDerive(unittest.TestCase):
    def test_shared_secret_agreement(self):
        p, g = 104729, 3
        alice = dh_keygen(p, g, private_key=1111)
        bob = dh_keygen(p, g, private_key=2222)
        s_alice = dh_derive(alice["private_key"], bob["public_key"], p)
        s_bob = dh_derive(bob["private_key"], alice["public_key"], p)
        self.assertEqual(s_alice, s_bob)
        self.assertEqual(s_alice, pow(g, 1111 * 2222, p))

    def test_invalid_parameters_raise(self):
        with self.assertRaises(ValueError):
            dh_derive(12, 100, 100)  # peer pub >= p
        with self.assertRaises(ValueError):
            dh_derive(12, 0, 100)    # peer pub < 1
        with self.assertRaises(ValueError):
            dh_derive(0, 50, 100)    # priv < 1
        with self.assertRaises(TypeError):
            dh_derive("12", 50, 100)


class TestDhCrack(unittest.TestCase):
    def test_crack_recovers_private_key(self):
        p, g = 104729, 3
        alice = dh_keygen(p, g, private_key=4321)
        recovered_priv = dh_crack(g, alice["public_key"], p)
        self.assertEqual(recovered_priv, 4321)

    def test_eavesdropper_computes_shared_secret(self):
        p, g = 104729, 3
        alice = dh_keygen(p, g, private_key=123)
        bob = dh_keygen(p, g, private_key=456)
        legit_secret = dh_derive(alice["private_key"], bob["public_key"], p)
        # Eavesdropper observes public values g, p, alice["public_key"], bob["public_key"]
        cracked_priv_a = dh_crack(g, alice["public_key"], p)
        stolen_secret = dh_derive(cracked_priv_a, bob["public_key"], p)
        self.assertEqual(stolen_secret, legit_secret)


class TestTopLevelCryptosystemImports(unittest.TestCase):
    def test_all_cryptosystem_functions_exported(self):
        for name in (
            "rsa_keygen", "rsa_encrypt", "rsa_decrypt", "rsa_crack",
            "dh_keygen", "dh_derive", "dh_crack",
        ):
            self.assertTrue(hasattr(foxmath, name), f"foxmath missing {name}")
            self.assertTrue(callable(getattr(foxmath, name)), f"foxmath.{name} not callable")


if __name__ == "__main__":
    unittest.main()
