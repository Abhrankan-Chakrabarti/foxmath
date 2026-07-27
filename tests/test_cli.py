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
        capture_output=True, text=True, timeout=15, input=stdin_input,
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


if __name__ == "__main__":
    unittest.main()
