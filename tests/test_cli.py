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


def run(*args, script=SCRIPT):
    result = subprocess.run(
        [sys.executable, script, *args],
        capture_output=True, text=True, timeout=15,
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


if __name__ == "__main__":
    unittest.main()
