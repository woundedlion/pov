#!/usr/bin/env python3
"""End-to-end tests for the glob-discovered test-suite non-empty guard."""

import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
CHECK = REPO / "tools" / "require_test_files.sh"
# The MSYS runtime expands wildcard argv entries before Bash sees them, which
# would split the glob into several arguments on Windows. Ignored elsewhere.
ENV = {**os.environ, "MSYS": "noglob"}


class RequireTestFiles(unittest.TestCase):
    def setUp(self):
        if shutil.which("bash") is None:
            self.skipTest("no Bash shell")
        tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.addCleanup(tmp.cleanup)
        self.tree = Path(tmp.name)
        populated = self.tree / "tools" / "sample_tests"
        populated.mkdir(parents=True)
        (populated / "test_sample.py").write_text("", encoding="utf-8")
        (populated / "test_other.py").write_text("", encoding="utf-8")
        (populated / "helper.py").write_text("", encoding="utf-8")
        (self.tree / "tools" / "empty_tests").mkdir(parents=True)

    def run_check(self, *args: str) -> subprocess.CompletedProcess:
        return subprocess.run(
            ["bash", CHECK.as_posix(), *args], cwd=self.tree, env=ENV,
            capture_output=True, text=True, check=False)

    def test_a_populated_glob_passes_and_lists_its_matches(self):
        done = self.run_check("tools/sample_tests/test*.py")
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertIn("2 test file(s) discovered", done.stdout)
        self.assertIn("tools/sample_tests/test_sample.py", done.stdout)
        self.assertIn("tools/sample_tests/test_other.py", done.stdout)
        # Discovery is the glob, not the directory: non-test files stay out.
        self.assertNotIn("helper.py", done.stdout)

    def test_an_empty_glob_fails(self):
        done = self.run_check("tools/empty_tests/test*.py")
        self.assertEqual(done.returncode, 1, done.stdout + done.stderr)
        self.assertIn("no test files match", done.stdout)

    def test_a_missing_directory_fails(self):
        done = self.run_check("tools/renamed_tests/test*.py")
        self.assertEqual(done.returncode, 1, done.stdout + done.stderr)
        self.assertIn("no test files match", done.stdout)

    def test_wrong_argument_count_is_a_usage_error(self):
        for args in ((), ("a", "b")):
            with self.subTest(args=args):
                done = self.run_check(*args)
                self.assertEqual(done.returncode, 2, done.stdout + done.stderr)
                self.assertIn("usage:", done.stderr)


if __name__ == "__main__":
    unittest.main()
