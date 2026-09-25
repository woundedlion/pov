"""Behavioral checks for tracked Python suite discovery."""

from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

RUNNER = Path(__file__).resolve().parents[1] / "run_python_tests.py"


class PythonRunner(unittest.TestCase):
    def run_fixture(self, source=None, extra=None):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            subprocess.run(["git", "-C", str(root), "init", "--quiet"], check=True)
            if source is not None:
                path = root / "new_suite" / "test_sample.py"
                path.parent.mkdir()
                path.write_text(source, encoding="utf-8")
                subprocess.run(["git", "-C", str(root), "add", "--",
                                "new_suite/test_sample.py"], check=True)
            if extra is not None:
                empty = root / "new_suite/test_extra.py"
                empty.write_text(extra, encoding="utf-8")
                subprocess.run(["git", "-C", str(root), "add", "--",
                                "new_suite/test_extra.py"], check=True)
            return subprocess.run([sys.executable, str(RUNNER), "--root", str(root)],
                                  capture_output=True, text=True, check=False)

    def test_empty_repository_fails(self):
        self.assertNotEqual(self.run_fixture().returncode, 0)

    def test_empty_suite_fails(self):
        self.assertNotEqual(self.run_fixture("# no cases\n").returncode, 0)

    def test_empty_file_fails_beside_live_tests(self):
        live = ("import unittest\nclass Live(unittest.TestCase):\n"
                "    def test_live(self): self.assertTrue(True)\n")
        self.assertNotEqual(self.run_fixture(live, extra="# empty\n").returncode, 0)

    def test_optional_skipped_file_passes_beside_live_tests(self):
        live = ("import unittest\nclass Live(unittest.TestCase):\n"
                "    def test_live(self): self.assertTrue(True)\n")
        skipped = ("import unittest\n@unittest.skip('fixture')\n"
                   "class Skipped(unittest.TestCase):\n"
                   "    def test_skipped(self): pass\n")
        result = self.run_fixture(live, extra=skipped)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_entirely_skipped_suite_fails(self):
        skipped = ("import unittest\n@unittest.skip('fixture')\n"
                   "class Skipped(unittest.TestCase):\n"
                   "    def test_skipped(self): pass\n")
        result = self.run_fixture(skipped, extra=skipped)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("no unskipped test cases in new_suite", result.stderr)

    def test_new_suite_is_discovered_and_failure_propagates(self):
        source = ("import unittest\nclass Sample(unittest.TestCase):\n"
                  "    def test_value(self):\n        self.assertEqual(1, VALUE)\n")
        result = self.run_fixture(source.replace("VALUE", "1"))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotEqual(self.run_fixture(source.replace("VALUE", "2")).returncode, 0)
