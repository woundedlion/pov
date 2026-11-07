"""Behavioral checks for tracked Python suite discovery."""

from pathlib import Path

import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from git_test_env import git_run, isolated_env

LIVE = ("import unittest\nclass Live(unittest.TestCase):\n"
        "    def test_live(self): self.assertTrue(True)\n")

RUNNER = Path(__file__).resolve().parents[1] / "run_python_tests.py"


class PythonRunner(unittest.TestCase):
    def run_fixture(self, source=None, extra=None, helper=False):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            git_run(root, "init", "--quiet")
            def track(name, text):
                path = root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(text, encoding="utf-8")
                git_run(root, "add", "--", name)

            if source is not None:
                track("new_suite/test_sample.py", source)
            if extra is not None:
                track("new_suite/test_extra.py", extra)
            if helper:
                track("new_suite/tests/data/helper.py", "# fixture helper\n")
            return subprocess.run([sys.executable, str(RUNNER), "--root", str(root)],
                                  capture_output=True, text=True, check=False, env=isolated_env())

    def test_empty_repository_fails(self):
        result = self.run_fixture()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("no tracked Python test suites discovered", result.stderr)
        self.assertNotIn("Traceback", result.stderr)

    def test_empty_suite_fails(self):
        result = self.run_fixture("# no cases\n")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("no test cases discovered in", result.stderr)
        self.assertIn("test_sample.py", result.stderr)
        self.assertNotIn("Traceback", result.stderr)

    def test_empty_file_fails_beside_live_tests(self):
        result = self.run_fixture(LIVE, extra="# empty\n")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("no test cases discovered in", result.stderr)
        self.assertIn("test_extra.py", result.stderr)
        self.assertNotIn("Traceback", result.stderr)

    def test_optional_skipped_file_passes_beside_live_tests(self):
        skipped = ("import unittest\n@unittest.skip('fixture')\n"
                   "class Skipped(unittest.TestCase):\n"
                   "    def test_skipped(self): pass\n")
        result = self.run_fixture(LIVE, extra=skipped)
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

    def test_nested_helpers_do_not_become_suites(self):
        result = self.run_fixture(LIVE, helper=True)
        self.assertEqual(result.returncode, 0, result.stderr)
