"""
End-to-end tests for foxmath's CLI, run as a subprocess to match real
usage (including simulated symlink invocation).
"""
import subprocess
import sys
import os
import shutil
import tempfile
import unittest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src", "foxmath"))
SCRIPT = os.path.join(ROOT, "foxmath.py")


def run(*args, script=SCRIPT, stdin_input=None):
    result = subprocess.run(
        [sys.executable, script, *args],
        capture_output=True, text=True, encoding="utf-8", timeout=15, input=stdin_input,
    )
    return result.returncode, result.stdout, result.stderr


class TestJsonFlag(unittest.TestCase):
    def test_json_after_subcommand(self):
        # This is the exact usage shown in the README's own example.
        code, out, err = run("legendre", "2", "7", "--json")
        self.assertEqual(code, 0)
        self.assertIn('"value": 1', out)

    def test_json_before_subcommand(self):
        code, out, err = run("--json", "legendre", "2", "7")
        self.assertEqual(code, 0)
        self.assertIn('"value": 1', out)

    def test_no_json_flag_omits_json_block(self):
        code, out, err = run("legendre", "2", "7")
        self.assertEqual(code, 0)
        self.assertNotIn("{", out)

    def test_sqrtmod(self):
        code, out, err = run("sqrtmod", "10", "13")
        self.assertEqual(code, 0)
        self.assertIn("sqrt(10) = +/-7 (mod 13)", out)

    def test_sqrtmod_json(self):
        code, out, err = run("sqrtmod", "10", "13", "--json")
        self.assertEqual(code, 0)
        self.assertIn('"roots": [', out)

    def test_json_output_is_pure_json_no_preamble(self):
        # Regression check: --json used to still print the human-readable
        # line first, breaking anything piping stdout into a JSON parser.
        code, out, err = run("pi", "--terms", "20", "--json")
        self.assertEqual(code, 0)
        import json
        parsed = json.loads(out)  # raises if there's any leading text
        self.assertEqual(parsed["command"], "pi")


class TestSymlinkInvocation(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _make_symlink_script(self, name):
        path = os.path.join(self.tmpdir, f"{name}.py")
        shutil.copy(SCRIPT, path)
        return path

    def test_legendre_symlink(self):
        script = self._make_symlink_script("foxmath-legendre")
        code, out, err = run("2", "7", script=script)
        self.assertEqual(code, 0)
        self.assertIn("Legendre (2/7) = 1", out)

    def test_sqrtmod_symlink(self):
        script = self._make_symlink_script("foxmath-sqrtmod")
        code, out, err = run("10", "13", script=script)
        self.assertEqual(code, 0)
        self.assertIn("sqrt(10) = +/-7 (mod 13)", out)

    def test_isprime_symlink(self):
        script = self._make_symlink_script("foxmath-isprime")
        code, out, err = run("104729", script=script)
        self.assertEqual(code, 0)
        self.assertIn("104729 is prime", out)

    def test_factor_symlink(self):
        script = self._make_symlink_script("foxmath-factor")
        code, out, err = run("360", script=script)
        self.assertEqual(code, 0)
        self.assertIn("360 = 2^3 * 3^2 * 5", out)

    def test_fermat_symlink(self):
        script = self._make_symlink_script("foxmath-fermat")
        code, out, err = run("5959", script=script)
        self.assertEqual(code, 0)
        self.assertIn("Factors of 5959", out)
        self.assertIn("p = 59", out)
        self.assertIn("q = 101", out)

    def test_totient_symlink(self):
        script = self._make_symlink_script("foxmath-totient")
        code, out, err = run("360", script=script)
        self.assertEqual(code, 0)
        self.assertIn("phi(360) = 96", out)
        self.assertIn("lambda(360) = 12", out)

    def test_order_symlink(self):
        script = self._make_symlink_script("foxmath-order")
        code, out, err = run("3", "7", script=script)
        self.assertEqual(code, 0)
        self.assertIn("ord_7(3) = 6 (primitive root)", out)

    def test_dlog_symlink(self):
        script = self._make_symlink_script("foxmath-dlog")
        code, out, err = run("3", "6", "7", script=script)
        self.assertEqual(code, 0)
        self.assertIn("3^3 ≡ 6 (mod 7)", out)

    def test_mobius_symlink(self):
        script = self._make_symlink_script("foxmath-mobius")
        code, out, err = run("30", script=script)
        self.assertEqual(code, 0)
        self.assertIn("mu(30) = -1", out)

    def test_sigma_symlink(self):
        script = self._make_symlink_script("foxmath-sigma")
        code, out, err = run("12", script=script)
        self.assertEqual(code, 0)
        self.assertIn("sigma(12) = 28", out)

    def test_tau_symlink(self):
        script = self._make_symlink_script("foxmath-tau")
        code, out, err = run("12", script=script)
        self.assertEqual(code, 0)
        self.assertIn("tau(12) = 6", out)

    def test_jacobi_symlink(self):
        script = self._make_symlink_script("foxmath-jacobi")
        code, out, err = run("10", "13", script=script)
        self.assertEqual(code, 0)
        self.assertIn("Jacobi (10/13) = 1", out)

    def test_solovay_symlink(self):
        script = self._make_symlink_script("foxmath-solovay")
        code, out, err = run("104729", "--rounds", "10", script=script)
        self.assertEqual(code, 0)
        self.assertIn("104729 is probably prime (Solovay-Strassen, 10 rounds)", out)

    def test_repl_symlink(self):
        script = self._make_symlink_script("foxmath-repl")
        code, out, err = run(script=script, stdin_input="print(jacobi_symbol(10, 13))\nexit()\n")
        self.assertEqual(code, 0)
        self.assertIn("FoxMath", out)
        self.assertIn("Interactive REPL", out)
        self.assertIn("1", out)

    def test_egcd_symlink(self):
        script = self._make_symlink_script("foxmath-egcd")
        code, out, err = run("30", "50", script=script)
        self.assertEqual(code, 0)
        self.assertIn("gcd(30, 50) = 10", out)

    def test_modinv_symlink(self):
        script = self._make_symlink_script("foxmath-modinv")
        code, out, err = run("3", "11", script=script)
        self.assertEqual(code, 0)
        self.assertIn("3^-1 ≡ 4 (mod 11)", out)

    def test_pi_symlink(self):
        script = self._make_symlink_script("foxmath-pi")
        code, out, err = run("--terms", "10", script=script)
        self.assertEqual(code, 0)
        self.assertIn("π ≈", out)

    def test_ecadd_symlink(self):
        script = self._make_symlink_script("foxmath-ecadd")
        code, out, err = run("--x1", "3", "--y1", "6", "--x2", "3", "--y2", "6",
                              "--a", "2", "--p", "97", script=script)
        self.assertEqual(code, 0)
        self.assertIn("Result point:", out)

    def test_crt_symlink(self):
        script = self._make_symlink_script("foxmath-crt")
        code, out, err = run("--r", "2", "3", "2", "--m", "3", "5", "7", script=script)
        self.assertEqual(code, 0)
        self.assertIn("x ≡ 23 (mod 105)", out)

    def test_cf_symlink(self):
        script = self._make_symlink_script("foxmath-cf")
        code, out, err = run("--num", "355", "--den", "113", script=script)
        self.assertEqual(code, 0)
        self.assertIn("[3; 7, 16]", out)

    def test_challenge_symlink(self):
        script = self._make_symlink_script("foxmath-challenge")
        code, out, err = run("--topic", "legendre", "--count", "2", "--seed", "1",
                              "--reveal", script=script)
        self.assertEqual(code, 0)
        self.assertIn("Answer:", out)

    def test_ecmul_symlink(self):
        script = self._make_symlink_script("foxmath-ecmul")
        code, out, err = run("--curve", "toy193", "--k", "7", script=script)
        self.assertEqual(code, 0)
        self.assertIn("Result point: (66, 24)", out)

    def test_rsa_symlink(self):
        script = self._make_symlink_script("foxmath-rsa")
        code, out, err = run("encrypt", "42", "17", "3233", script=script)
        self.assertEqual(code, 0)
        self.assertIn("Ciphertext c = 2557", out)

    def test_dh_symlink(self):
        script = self._make_symlink_script("foxmath-dh")
        code, out, err = run("derive", "4", "8", "23", script=script)
        self.assertEqual(code, 0)
        self.assertIn("Shared secret S = 2", out)

    def test_mersenne_symlink(self):
        script = self._make_symlink_script("foxmath-mersenne")
        code, out, err = run("7", script=script)
        self.assertEqual(code, 0)
        self.assertIn("M_7 = 2^7 - 1 = 127", out)
        self.assertIn("Prime (Lucas-Lehmer test passed)", out)

    def test_perfect_symlink(self):
        script = self._make_symlink_script("foxmath-perfect")
        code, out, err = run("28", script=script)
        self.assertEqual(code, 0)
        self.assertIn("28 is a perfect number", out)

    def test_symlink_with_json(self):
        script = self._make_symlink_script("foxmath-legendre")
        code, out, err = run("2", "7", "--json", script=script)
        self.assertEqual(code, 0)
        self.assertIn('"value": 1', out)


class TestNormalUsage(unittest.TestCase):
    def test_pi(self):
        code, out, err = run("pi", "--terms", "50")
        self.assertEqual(code, 0)
        self.assertIn("3.14159265358979", out)

    def test_ecadd_point_at_infinity_errors_clearly(self):
        code, out, err = run("ecadd", "--x1", "3", "--y1", "6",
                              "--x2", "3", "--y2", "91", "--a", "2", "--p", "97")
        self.assertEqual(code, 1)
        self.assertIn("point at infinity", err)

    def test_ecadd_curve_secp256k1_defaults_to_doubling_generator(self):
        code, out, err = run("ecadd", "--curve", "secp256k1")
        self.assertEqual(code, 0)
        self.assertIn("89565891926547004231252920425935692360644145829622209833684329913297188986597", out)
        self.assertIn("on curve", out)

    def test_ecadd_curve_manual_a_override(self):
        # --a explicitly provided should override the curve preset's a.
        code, out, err = run("ecadd", "--curve", "toy97", "--a", "2")
        self.assertEqual(code, 0)

    def test_ecadd_no_args_at_all_errors_clearly(self):
        code, out, err = run("ecadd")
        self.assertEqual(code, 1)
        self.assertIn("Error", err)

    def test_ecadd_backward_compatible_manual_call(self):
        # Original documented usage: --a omitted, defaults to -3.
        code, out, err = run("ecadd", "--x1", "1", "--y1", "2",
                              "--x2", "3", "--y2", "4", "--p", "17")
        self.assertEqual(code, 0)
        self.assertIn("Result point: (14, 2)", out)

    def test_ecmul_readme_example_toy193(self):
        # Exact usage/output shown in the README.
        code, out, err = run("ecmul", "--curve", "toy193", "--k", "7")
        self.assertEqual(code, 0)
        self.assertIn("Result point: (66, 24)", out)
        self.assertIn("on curve", out)

    def test_ecmul_readme_example_identity(self):
        code, out, err = run("ecmul", "--curve", "toy193", "--k", "193")
        self.assertEqual(code, 0)
        self.assertIn("point at infinity", out)

    def test_ecmul_secp256k1_matches_known_doubling(self):
        code, out, err = run("ecmul", "--curve", "secp256k1", "--k", "2")
        self.assertEqual(code, 0)
        self.assertIn("89565891926547004231252920425935692360644145829622209833684329913297188986597", out)

    def test_ecmul_json_output(self):
        code, out, err = run("ecmul", "--curve", "toy193", "--k", "7", "--json")
        self.assertEqual(code, 0)
        import json
        parsed = json.loads(out)
        self.assertEqual(parsed["result"], [66, 24])
        self.assertTrue(parsed["result_on_curve"])

    def test_ecmul_json_output_infinity(self):
        code, out, err = run("ecmul", "--curve", "toy193", "--k", "193", "--json")
        self.assertEqual(code, 0)
        import json
        parsed = json.loads(out)
        self.assertEqual(parsed["result"], "infinity")

    def test_ecmul_manual_mode_no_curve(self):
        code, out, err = run("ecmul", "--p", "97", "--x", "3", "--y", "6",
                              "--a", "2", "--k", "3")
        self.assertEqual(code, 0)
        self.assertIn("Result point:", out)

    def test_ecmul_a_override(self):
        code, out, err = run("ecmul", "--curve", "toy97", "--a", "2", "--k", "3")
        self.assertEqual(code, 0)

    def test_ecmul_missing_curve_and_incomplete_manual_args_errors(self):
        code, out, err = run("ecmul", "--k", "5")
        self.assertEqual(code, 1)
        self.assertIn("Error", err)

    def test_ecmul_negative_k_errors(self):
        code, out, err = run("ecmul", "--curve", "toy193", "--k", "-5")
        self.assertEqual(code, 1)
        self.assertIn("Error", err)

    def test_crt(self):
        code, out, err = run("crt", "--r", "2", "3", "2", "--m", "3", "5", "7")
        self.assertEqual(code, 0)
        self.assertIn("x ≡ 23 (mod 105)", out)

    def test_crt_json_is_pure_json(self):
        code, out, err = run("crt", "--r", "2", "3", "2", "--m", "3", "5", "7", "--json")
        self.assertEqual(code, 0)
        import json
        parsed = json.loads(out)
        self.assertEqual(parsed["x"], 23)
        self.assertEqual(parsed["mod"], 105)

    def test_crt_inconsistent_system_errors(self):
        code, out, err = run("crt", "--r", "1", "2", "--m", "4", "6")
        self.assertEqual(code, 1)
        self.assertIn("inconsistent", err)

    def test_cf(self):
        code, out, err = run("cf", "--num", "355", "--den", "113")
        self.assertEqual(code, 0)
        self.assertIn("[3; 7, 16]", out)
        self.assertIn("355/113", out)

    def test_cf_json_is_pure_json(self):
        code, out, err = run("cf", "--num", "355", "--den", "113", "--json")
        self.assertEqual(code, 0)
        import json
        parsed = json.loads(out)
        self.assertEqual(parsed["cf"], [3, 7, 16])
        self.assertEqual(parsed["convergents"][-1], [355, 113])

    def test_cf_zero_denominator_errors(self):
        code, out, err = run("cf", "--num", "5", "--den", "0")
        self.assertEqual(code, 1)
        self.assertIn("Error", err)

    def test_challenge_reveal_is_deterministic_with_seed(self):
        code1, out1, err1 = run("challenge", "--topic", "legendre", "--count", "3",
                                 "--seed", "42", "--reveal")
        code2, out2, err2 = run("challenge", "--topic", "legendre", "--count", "3",
                                 "--seed", "42", "--reveal")
        self.assertEqual(code1, 0)
        self.assertEqual(out1, out2)

    def test_challenge_interactive_scores_correct_answers(self):
        # Known answers for --topic legendre --seed 42 --count 3 (from reveal mode).
        code, out, err = run("challenge", "--topic", "legendre", "--count", "3",
                              "--seed", "42", stdin_input="-1\n-1\n1\n")
        self.assertEqual(code, 0)
        self.assertIn("Score: 3/3", out)

    def test_challenge_interactive_scores_wrong_answers(self):
        code, out, err = run("challenge", "--topic", "legendre", "--count", "3",
                              "--seed", "42", stdin_input="0\n0\n0\n")
        self.assertEqual(code, 0)
        self.assertIn("Score: 0/3", out)

    def test_challenge_interactive_handles_non_numeric_input(self):
        code, out, err = run("challenge", "--topic", "legendre", "--count", "1",
                              "--seed", "42", stdin_input="not-a-number\n")
        self.assertEqual(code, 0)
        self.assertIn("Score: 0/1", out)

    def test_challenge_reveal_json(self):
        code, out, err = run("challenge", "--topic", "crt", "--count", "2",
                              "--seed", "1", "--reveal", "--json")
        self.assertEqual(code, 0)
        import json
        parsed = json.loads(out)
        self.assertEqual(len(parsed["problems"]), 2)

    def test_challenge_zero_count_errors(self):
        code, out, err = run("challenge", "--count", "0")
        self.assertEqual(code, 1)
        self.assertIn("Error", err)

    def test_sqrtmod_non_residue_errors(self):
        code, out, err = run("sqrtmod", "5", "17")
        self.assertEqual(code, 1)
        self.assertIn("no square root", err)

    def test_sqrtmod_zero(self):
        code, out, err = run("sqrtmod", "0", "17")
        self.assertEqual(code, 0)
        self.assertIn("sqrt(0) = 0 (mod 17)", out)

    def test_legendre_composite_p_errors(self):
        code, out, err = run("legendre", "2", "9")
        self.assertEqual(code, 1)
        self.assertIn("p must be an odd prime", err)

    def test_sqrtmod_composite_p_errors(self):
        code, out, err = run("sqrtmod", "2", "9")
        self.assertEqual(code, 1)
        self.assertIn("p must be an odd prime", err)

    def test_isprime_prime(self):
        code, out, err = run("isprime", "104729")
        self.assertEqual(code, 0)
        self.assertIn("104729 is prime", out)

    def test_isprime_composite(self):
        code, out, err = run("isprime", "561")
        self.assertEqual(code, 0)
        self.assertIn("561 is composite", out)

    def test_isprime_less_than_two(self):
        code, out, err = run("isprime", "1")
        self.assertEqual(code, 0)
        self.assertIn("1 is not prime", out)

    def test_isprime_json(self):
        code, out, err = run("isprime", "104729", "--json")
        self.assertEqual(code, 0)
        import json
        data = json.loads(out)
        self.assertEqual(data["command"], "isprime")
        self.assertEqual(data["n"], 104729)
        self.assertTrue(data["is_prime"])

    def test_egcd(self):
        code, out, err = run("egcd", "30", "50")
        self.assertEqual(code, 0)
        self.assertIn("gcd(30, 50) = 10", out)
        self.assertIn("30*(2) + 50*(-1)", out)

    def test_egcd_json(self):
        code, out, err = run("egcd", "30", "50", "--json")
        self.assertEqual(code, 0)
        import json
        data = json.loads(out)
        self.assertEqual(data["gcd"], 10)
        self.assertEqual(data["x"], 2)
        self.assertEqual(data["y"], -1)

    def test_modinv(self):
        code, out, err = run("modinv", "3", "11")
        self.assertEqual(code, 0)
        self.assertIn("3^-1 ≡ 4 (mod 11)", out)

    def test_modinv_json(self):
        code, out, err = run("modinv", "3", "11", "--json")
        self.assertEqual(code, 0)
        import json
        data = json.loads(out)
        self.assertEqual(data["inverse"], 4)

    def test_modinv_non_coprime_errors(self):
        code, out, err = run("modinv", "6", "9")
        self.assertEqual(code, 1)
        self.assertIn("no modular inverse", err)

    def test_challenge_sqrtmod(self):
        code, out, err = run("challenge", "--topic", "sqrtmod", "--count", "2", "--seed", "42", "--reveal")
        self.assertEqual(code, 0)
        self.assertIn("square root", out)
        self.assertIn("Answer:", out)

    def test_factor_composite(self):
        code, out, err = run("factor", "360")
        self.assertEqual(code, 0)
        self.assertIn("360 = 2^3 * 3^2 * 5", out)

    def test_factor_prime(self):
        code, out, err = run("factor", "13")
        self.assertEqual(code, 0)
        self.assertIn("13 is prime", out)

    def test_factor_one(self):
        code, out, err = run("factor", "1")
        self.assertEqual(code, 0)
        self.assertIn("1 has no prime factors", out)

    def test_factor_negative_or_zero_errors(self):
        code, out, err = run("factor", "0")
        self.assertEqual(code, 1)
        self.assertIn("positive integer >= 1", err)

    def test_factor_json(self):
        code, out, err = run("factor", "360", "--json")
        self.assertEqual(code, 0)
        import json
        data = json.loads(out)
        self.assertEqual(data["command"], "factor")
        self.assertEqual(data["n"], 360)
        self.assertEqual(data["expression"], "2^3 * 3^2 * 5")
        self.assertEqual(data["prime_factors"], [2, 3, 5])
        self.assertEqual(data["factors"], [
            {"prime": 2, "power": 3},
            {"prime": 3, "power": 2},
            {"prime": 5, "power": 1},
        ])
        self.assertFalse(data["is_prime"])

    def test_challenge_factor(self):
        code, out, err = run("challenge", "--topic", "factor", "--count", "2", "--seed", "42", "--reveal")
        self.assertEqual(code, 0)
        self.assertIn("smallest prime factor", out)
        self.assertIn("Answer:", out)

    def test_totient(self):
        code, out, err = run("totient", "360")
        self.assertEqual(code, 0)
        self.assertIn("phi(360) = 96", out)
        self.assertIn("lambda(360) = 12", out)

    def test_totient_cyclic(self):
        code, out, err = run("totient", "7")
        self.assertEqual(code, 0)
        self.assertIn("phi(7) = 6", out)
        self.assertIn("lambda(7) = 6 (cyclic)", out)

    def test_totient_json(self):
        code, out, err = run("totient", "360", "--json")
        self.assertEqual(code, 0)
        import json
        data = json.loads(out)
        self.assertEqual(data["command"], "totient")
        self.assertEqual(data["n"], 360)
        self.assertEqual(data["phi"], 96)
        self.assertEqual(data["lambda"], 12)
        self.assertFalse(data["is_cyclic"])

    def test_totient_invalid_errors(self):
        code, out, err = run("totient", "0")
        self.assertEqual(code, 1)
        self.assertIn("positive integer >= 1", err)

    def test_order(self):
        code, out, err = run("order", "3", "7")
        self.assertEqual(code, 0)
        self.assertIn("ord_7(3) = 6 (primitive root)", out)

    def test_order_non_primitive(self):
        code, out, err = run("order", "2", "7")
        self.assertEqual(code, 0)
        self.assertIn("ord_7(2) = 3", out)
        self.assertNotIn("primitive root", out)

    def test_order_json(self):
        code, out, err = run("order", "3", "7", "--json")
        self.assertEqual(code, 0)
        import json
        data = json.loads(out)
        self.assertEqual(data["command"], "order")
        self.assertEqual(data["a"], 3)
        self.assertEqual(data["m"], 7)
        self.assertEqual(data["order"], 6)
        self.assertTrue(data["is_primitive_root"])

    def test_order_non_coprime_errors(self):
        code, out, err = run("order", "2", "8")
        self.assertEqual(code, 1)
        self.assertIn("not coprime", err)

    def test_challenge_totient(self):
        code, out, err = run("challenge", "--topic", "totient", "--count", "2", "--seed", "42", "--reveal")
        self.assertEqual(code, 0)
        self.assertIn("Euler's totient", out)
        self.assertIn("Answer:", out)

    def test_dlog(self):
        code, out, err = run("dlog", "3", "6", "7")
        self.assertEqual(code, 0)
        self.assertIn("3^3 ≡ 6 (mod 7)", out)

    def test_dlog_json(self):
        code, out, err = run("dlog", "3", "6", "7", "--json")
        self.assertEqual(code, 0)
        import json
        data = json.loads(out)
        self.assertEqual(data["command"], "dlog")
        self.assertEqual(data["g"], 3)
        self.assertEqual(data["h"], 6)
        self.assertEqual(data["p"], 7)
        self.assertEqual(data["x"], 3)
        self.assertEqual(data["order"], 6)

    def test_dlog_no_solution_errors(self):
        code, out, err = run("dlog", "2", "3", "7")
        self.assertEqual(code, 1)
        self.assertIn("No discrete logarithm exists", err)

    def test_dlog_non_coprime_errors(self):
        code, out, err = run("dlog", "2", "4", "8")
        self.assertEqual(code, 1)
        self.assertIn("must be coprime", err)

    def test_mobius(self):
        code, out, err = run("mobius", "30")
        self.assertEqual(code, 0)
        self.assertIn("mu(30) = -1", out)

    def test_mobius_json(self):
        code, out, err = run("mobius", "30", "--json")
        self.assertEqual(code, 0)
        import json
        data = json.loads(out)
        self.assertEqual(data["command"], "mobius")
        self.assertEqual(data["n"], 30)
        self.assertEqual(data["value"], -1)

    def test_sigma(self):
        code, out, err = run("sigma", "12")
        self.assertEqual(code, 0)
        self.assertIn("sigma(12) = 28", out)

    def test_sigma_k(self):
        code, out, err = run("sigma", "6", "--k", "2")
        self.assertEqual(code, 0)
        self.assertIn("sigma_2(6) = 50", out)

    def test_sigma_json(self):
        code, out, err = run("sigma", "12", "--json")
        self.assertEqual(code, 0)
        import json
        data = json.loads(out)
        self.assertEqual(data["command"], "sigma")
        self.assertEqual(data["n"], 12)
        self.assertEqual(data["k"], 1)
        self.assertEqual(data["value"], 28)

    def test_tau(self):
        code, out, err = run("tau", "12")
        self.assertEqual(code, 0)
        self.assertIn("tau(12) = 6", out)

    def test_tau_json(self):
        code, out, err = run("tau", "12", "--json")
        self.assertEqual(code, 0)
        import json
        data = json.loads(out)
        self.assertEqual(data["command"], "tau")
        self.assertEqual(data["n"], 12)
        self.assertEqual(data["value"], 6)

    def test_challenge_dlog(self):
        code, out, err = run("challenge", "--topic", "dlog", "--count", "2", "--seed", "42", "--reveal")
        self.assertEqual(code, 0)
        self.assertIn("such that", out)
        self.assertIn("Answer:", out)

    def test_jacobi(self):
        code, out, err = run("jacobi", "10", "13")
        self.assertEqual(code, 0)
        self.assertIn("Jacobi (10/13) = 1", out)

    def test_jacobi_json(self):
        code, out, err = run("jacobi", "10", "13", "--json")
        self.assertEqual(code, 0)
        import json
        data = json.loads(out)
        self.assertEqual(data["command"], "jacobi")
        self.assertEqual(data["a"], 10)
        self.assertEqual(data["n"], 13)
        self.assertEqual(data["value"], 1)

    def test_jacobi_even_modulus_error(self):
        code, out, err = run("jacobi", "2", "8")
        self.assertEqual(code, 1)
        self.assertIn("positive odd integer", err)

    def test_solovay_prime(self):
        code, out, err = run("solovay", "104729", "--rounds", "15")
        self.assertEqual(code, 0)
        self.assertIn("104729 is probably prime (Solovay-Strassen, 15 rounds)", out)

    def test_solovay_prime_small(self):
        code, out, err = run("solovay", "3")
        self.assertEqual(code, 0)
        self.assertIn("3 is prime", out)

    def test_solovay_composite(self):
        code, out, err = run("solovay", "561")
        self.assertEqual(code, 0)
        self.assertIn("561 is composite", out)

    def test_solovay_json(self):
        code, out, err = run("solovay", "17", "--json")
        self.assertEqual(code, 0)
        import json
        data = json.loads(out)
        self.assertEqual(data["command"], "solovay")
        self.assertEqual(data["n"], 17)
        self.assertEqual(data["is_prime"], True)

    def test_repl(self):
        code, out, err = run("repl", stdin_input="print(jacobi_symbol(2, 9))\nexit()\n")
        self.assertEqual(code, 0)
        self.assertIn("FoxMath", out)
        self.assertIn("Interactive REPL", out)
        self.assertIn("1", out)

    def test_challenge_jacobi(self):
        code, out, err = run("challenge", "--topic", "jacobi", "--count", "3", "--seed", "42", "--reveal")
        self.assertEqual(code, 0)
        self.assertIn("Jacobi symbol (", out)
        self.assertIn("Answer:", out)

    def test_rsa_keygen(self):
        code, out, err = run("rsa", "keygen", "--bits", "32", "--seed", "42")
        self.assertEqual(code, 0)
        self.assertIn("RSA 32-bit Keypair:", out)
        self.assertIn("Modulus   (n) =", out)
        self.assertIn("Public  e     =", out)
        self.assertIn("Private d     =", out)

    def test_rsa_keygen_json(self):
        code, out, err = run("rsa", "keygen", "--bits", "32", "--seed", "42", "--json")
        self.assertEqual(code, 0)
        import json
        data = json.loads(out)
        self.assertEqual(data["command"], "rsa")
        self.assertEqual(data["action"], "keygen")
        self.assertEqual(data["e"], 65537)
        self.assertIn("p", data)
        self.assertIn("q", data)
        self.assertIn("d", data)
        self.assertIn("n", data)

    def test_rsa_encrypt_decrypt(self):
        # encrypt
        code, out, err = run("rsa", "encrypt", "42", "17", "3233")
        self.assertEqual(code, 0)
        self.assertIn("Ciphertext c = 2557", out)

        # decrypt
        code, out, err = run("rsa", "decrypt", "2557", "413", "3233")
        self.assertEqual(code, 0)
        self.assertIn("Decrypted message m = 42", out)

        # decrypt with CRT
        code, out, err = run("rsa", "decrypt", "2557", "413", "3233", "--p", "61", "--q", "53")
        self.assertEqual(code, 0)
        self.assertIn("Decrypted message m = 42 (CRT accelerated)", out)

    def test_rsa_encrypt_decrypt_json(self):
        code, out, err = run("rsa", "encrypt", "42", "17", "3233", "--json")
        self.assertEqual(code, 0)
        import json
        data = json.loads(out)
        self.assertEqual(data["c"], 2557)

        code, out, err = run("rsa", "decrypt", "2557", "413", "3233", "--json")
        self.assertEqual(code, 0)
        data = json.loads(out)
        self.assertEqual(data["m"], 42)
        self.assertEqual(data["crt"], False)

    def test_rsa_crack(self):
        code, out, err = run("rsa", "crack", "3233", "17", "2790")
        self.assertEqual(code, 0)
        self.assertIn("Factored n = 3233:", out)
        self.assertIn("p = 61", out)
        self.assertIn("q = 53", out)
        self.assertIn("Recovered d  = 413", out)
        self.assertIn("Decrypted m  = 65", out)

    def test_rsa_crack_json(self):
        code, out, err = run("rsa", "crack", "3233", "17", "--c", "2790", "--json")
        self.assertEqual(code, 0)
        import json
        data = json.loads(out)
        self.assertEqual(data["command"], "rsa")
        self.assertEqual(data["action"], "crack")
        self.assertEqual(data["p"], 61)
        self.assertEqual(data["q"], 53)
        self.assertEqual(data["d"], 413)
        self.assertEqual(data["m"], 65)

    def test_dh_keygen(self):
        code, out, err = run("dh", "keygen", "23", "5", "--priv", "4")
        self.assertEqual(code, 0)
        self.assertIn("Diffie-Hellman Keypair (p=23, g=5):", out)
        self.assertIn("Private key (a) = 4", out)
        self.assertIn("Public key  (A) = 4", out)

    def test_dh_keygen_json(self):
        code, out, err = run("dh", "keygen", "23", "5", "--seed", "42", "--json")
        self.assertEqual(code, 0)
        import json
        data = json.loads(out)
        self.assertEqual(data["command"], "dh")
        self.assertEqual(data["action"], "keygen")
        self.assertEqual(data["p"], 23)
        self.assertEqual(data["g"], 5)
        self.assertIn("private_key", data)
        self.assertIn("public_key", data)

    def test_dh_derive_and_crack(self):
        # derive
        code, out, err = run("dh", "derive", "4", "8", "23")
        self.assertEqual(code, 0)
        self.assertIn("Shared secret S = 2", out)

        # crack
        code, out, err = run("dh", "crack", "5", "4", "23")
        self.assertEqual(code, 0)
        self.assertIn("Discrete logarithm solved via BSGS:", out)
        self.assertIn("Recovered private key a = 4", out)

    def test_dh_derive_and_crack_json(self):
        code, out, err = run("dh", "derive", "4", "8", "23", "--json")
        self.assertEqual(code, 0)
        import json
        data = json.loads(out)
        self.assertEqual(data["command"], "dh")
        self.assertEqual(data["action"], "derive")
        self.assertEqual(data["shared_secret"], 2)

        code, out, err = run("dh", "crack", "5", "4", "23", "--json")
        self.assertEqual(code, 0)
        data = json.loads(out)
        self.assertEqual(data["command"], "dh")
        self.assertEqual(data["action"], "crack")
        self.assertEqual(data["private_key"], 4)

    def test_challenge_rsa(self):
        code, out, err = run("challenge", "--topic", "rsa", "--count", "3", "--seed", "42", "--reveal")
        self.assertEqual(code, 0)
        self.assertIn("In RSA with e =", out)
        self.assertIn("Answer:", out)

    def test_challenge_dh(self):
        code, out, err = run("challenge", "--topic", "dh", "--count", "3", "--seed", "42", "--reveal")
        self.assertEqual(code, 0)
        self.assertIn("In Diffie-Hellman with p =", out)
        self.assertIn("Answer:", out)


class TestFermatCli(unittest.TestCase):
    def test_fermat_positional(self):
        code, out, err = run("fermat", "5959")
        self.assertEqual(code, 0)
        self.assertIn("Factors of 5959", out)
        self.assertIn("p = 59", out)
        self.assertIn("q = 101", out)
        self.assertIn("Iterations = 2", out)

    def test_fermat_flag_mod(self):
        code, out, err = run("fermat", "-n", "5959")
        self.assertEqual(code, 0)
        self.assertIn("p = 59", out)
        self.assertIn("q = 101", out)

    def test_fermat_hex(self):
        code, out, err = run("fermat", "-n", "0x1747")
        self.assertEqual(code, 0)
        self.assertIn("p = 59", out)
        self.assertIn("q = 101", out)

    def test_fermat_stdin(self):
        code, out, err = run("fermat", stdin_input="5959\n")
        self.assertEqual(code, 0)
        self.assertIn("p = 59", out)
        self.assertIn("q = 101", out)

    def test_fermat_quiet(self):
        code, out, err = run("fermat", "-n", "5959", "--quiet")
        self.assertEqual(code, 0)
        self.assertIn("p = 59", out)
        self.assertIn("q = 101", out)

    def test_fermat_start_iter(self):
        code, out, err = run("fermat", "-n", "5959", "-s", "1")
        self.assertEqual(code, 0)
        self.assertIn("p = 59", out)
        self.assertIn("q = 101", out)

    def test_fermat_json(self):
        code, out, err = run("fermat", "-n", "5959", "--json")
        self.assertEqual(code, 0)
        import json
        data = json.loads(out)
        self.assertEqual(data["command"], "fermat")
        self.assertEqual(data["n"], 5959)
        self.assertEqual(data["p"], 59)
        self.assertEqual(data["q"], 101)
        self.assertEqual(data["iterations"], 2)
        self.assertTrue(data["success"])

    def test_fermat_failure(self):
        code, out, err = run("fermat", "-n", "10")
        self.assertEqual(code, 1)
        self.assertIn("Failed to factor 10", out)

    def test_fermat_failure_json(self):
        code, out, err = run("fermat", "-n", "10", "--json")
        self.assertEqual(code, 1)
        import json
        data = json.loads(out)
        self.assertEqual(data["command"], "fermat")
        self.assertFalse(data["success"])

    def test_fermat_negative_start_error(self):
        code, out, err = run("fermat", "-n", "5959", "-s", "-1")
        self.assertEqual(code, 1)
        self.assertIn("must be non-negative", err)

    def test_fermat_no_args_error(self):
        code, out, err = run("fermat")
        self.assertEqual(code, 1)
        self.assertIn("Modulus to factor must be provided", err)

    def test_mersenne_prime(self):
        code, out, err = run("mersenne", "7")
        self.assertEqual(code, 0)
        self.assertIn("M_7 = 2^7 - 1 = 127", out)
        self.assertIn("Prime (Lucas-Lehmer test passed)", out)
        self.assertIn("Corresponding perfect number: 8128", out)

    def test_mersenne_composite(self):
        code, out, err = run("mersenne", "11")
        self.assertEqual(code, 0)
        self.assertIn("M_11 = 2^11 - 1 = 2047", out)
        self.assertIn("Composite (Lucas-Lehmer test failed)", out)

    def test_mersenne_composite_exponent(self):
        code, out, err = run("mersenne", "4")
        self.assertEqual(code, 0)
        self.assertIn("M_4 = 2^4 - 1 = 15", out)
        self.assertIn("Composite (exponent 4 is composite)", out)

    def test_mersenne_list(self):
        code, out, err = run("mersenne", "--list", "--limit", "31")
        self.assertEqual(code, 0)
        self.assertIn("Mersenne prime exponents p <= 31:", out)
        self.assertIn("2, 3, 5, 7, 13, 17, 19, 31", out)

    def test_mersenne_json(self):
        code, out, err = run("mersenne", "5", "--json")
        self.assertEqual(code, 0)
        import json
        data = json.loads(out)
        self.assertEqual(data["command"], "mersenne")
        self.assertEqual(data["p"], 5)
        self.assertEqual(data["mp"], 31)
        self.assertTrue(data["is_prime"])
        self.assertEqual(data["perfect_number"], 496)

    def test_perfect_test_true(self):
        code, out, err = run("perfect", "496")
        self.assertEqual(code, 0)
        self.assertIn("496 is a perfect number", out)

    def test_perfect_test_false(self):
        code, out, err = run("perfect", "12")
        self.assertEqual(code, 0)
        self.assertIn("12 is not a perfect number", out)

    def test_perfect_from_p(self):
        code, out, err = run("perfect", "7", "--from-p")
        self.assertEqual(code, 0)
        self.assertIn("Even perfect number for p = 7:", out)
        self.assertIn("8128", out)

    def test_perfect_from_p_json(self):
        code, out, err = run("perfect", "5", "--from-p", "--json")
        self.assertEqual(code, 0)
        import json
        data = json.loads(out)
        self.assertEqual(data["command"], "perfect")
        self.assertEqual(data["action"], "generate")
        self.assertEqual(data["perfect_number"], 496)

    def test_rsa_crack_fermat_method(self):
        code, out, err = run("rsa", "crack", "3233", "17", "--method", "fermat")
        self.assertEqual(code, 0)
        self.assertIn("Factored n = 3233:", out)
        self.assertIn("Method       = fermat", out)
        self.assertIn("p = 61", out)
        self.assertIn("q = 53", out)

    def test_challenge_mersenne(self):
        code, out, err = run("challenge", "--topic", "mersenne", "--count", "2", "--seed", "42", "--reveal")
        self.assertEqual(code, 0)
        self.assertIn("Is the Mersenne number M_", out)
        self.assertIn("Answer:", out)

    def test_challenge_perfect(self):
        code, out, err = run("challenge", "--topic", "perfect", "--count", "2", "--seed", "42", "--reveal")
        self.assertEqual(code, 0)
        self.assertIn("perfect number", out)
        self.assertIn("Answer:", out)

    def test_pm1_cli(self):
        code, out, err = run("pm1", "1393429", "--b1", "20", "--b2", "100")
        self.assertEqual(code, 0)
        self.assertIn("Factor found: p = 1009", out)
        self.assertIn("Cofactor:     q = 1381", out)
        self.assertIn("Stage:        1", out)

    def test_pm1_json(self):
        code, out, err = run("pm1", "1393429", "--b1", "20", "--b2", "100", "--json")
        self.assertEqual(code, 0)
        import json
        data = json.loads(out)
        self.assertEqual(data["command"], "pm1")
        self.assertEqual(data["factor"], 1009)
        self.assertEqual(data["cofactor"], 1381)
        self.assertEqual(data["stage"], 1)

    def test_challenge_pm1(self):
        code, out, err = run("challenge", "--topic", "pm1", "--count", "2", "--seed", "42", "--reveal")
        self.assertEqual(code, 0)
        self.assertIn("hint: p - 1 is smooth", out)
        self.assertIn("Answer:", out)

    def test_ecm_cli(self):
        code, out, err = run("ecm", "1393429", "--b1", "50", "--curves", "20", "--seed", "42")
        self.assertEqual(code, 0)
        self.assertIn("Factor found:", out)
        self.assertIn("Curve index:", out)

    def test_ecm_json(self):
        code, out, err = run("ecm", "1393429", "--b1", "50", "--curves", "20", "--seed", "42", "--json")
        self.assertEqual(code, 0)
        import json
        data = json.loads(out)
        self.assertEqual(data["command"], "ecm")
        self.assertTrue(data["success"])
        self.assertIn(data["factor"], (1009, 1381))
        self.assertEqual(data["factor"] * data["cofactor"], 1393429)

    def test_challenge_ecm(self):
        code, out, err = run("challenge", "--topic", "ecm", "--count", "2", "--seed", "42", "--reveal")
        self.assertEqual(code, 0)
        self.assertIn("hint: use Lenstra elliptic curve factorization", out)
        self.assertIn("Answer:", out)

    def test_dixon_cli(self):
        code, out, err = run("dixon", "84923", "--b", "30")
        self.assertEqual(code, 0)
        self.assertIn("Factor found:", out)
        self.assertIn("Congruence:", out)

    def test_dixon_json(self):
        code, out, err = run("dixon", "84923", "--b", "30", "--json")
        self.assertEqual(code, 0)
        import json
        data = json.loads(out)
        self.assertEqual(data["command"], "dixon")
        self.assertTrue(data["success"])
        self.assertIn(data["factor"], (163, 521))
        self.assertEqual(data["factor"] * data["cofactor"], 84923)

    def test_challenge_dixon(self):
        code, out, err = run("challenge", "--topic", "dixon", "--count", "2", "--seed", "42", "--reveal")
        self.assertEqual(code, 0)
        self.assertIn("hint: use Dixon's random squares method", out)
        self.assertIn("Answer:", out)



if __name__ == "__main__":
    unittest.main()

