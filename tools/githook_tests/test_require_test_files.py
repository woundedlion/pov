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
        (populated / "sample.test.mjs").write_text("", encoding="utf-8")
        (populated / "other.test.mjs").write_text("", encoding="utf-8")
        (populated / "helper.mjs").write_text("", encoding="utf-8")
        (self.tree / "tools" / "empty_tests").mkdir(parents=True)

    def run_check(self, *args: str) -> subprocess.CompletedProcess:
        return subprocess.run(
            ["bash", CHECK.as_posix(), *args], cwd=self.tree, env=ENV,
            capture_output=True, text=True, check=False)

    def test_a_populated_glob_passes_and_lists_its_matches(self):
        done = self.run_check("tools/sample_tests/*.test.mjs")
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertIn("2 test file(s) discovered", done.stdout)
        self.assertIn("tools/sample_tests/sample.test.mjs", done.stdout)
        self.assertIn("tools/sample_tests/other.test.mjs", done.stdout)
        # Discovery is the glob, not the directory: non-test files stay out.
        self.assertNotIn("helper.mjs", done.stdout)

    def test_nested_test_outside_the_glob_fails(self):
        nested = self.tree / "tools" / "sample_tests" / "nested"
        nested.mkdir()
        (nested / "hidden.test.mjs").write_text("", encoding="utf-8")
        done = self.run_check("tools/sample_tests/*.test.mjs")
        self.assertEqual(done.returncode, 1, done.stdout + done.stderr)
        self.assertIn("hidden.test.mjs", done.stdout)
        self.assertIn("unreachable", done.stdout)

    def test_javascript_spec_with_unmatched_extension_fails(self):
        scripts = self.tree / "scripts"
        scripts.mkdir()
        (scripts / "active.test.mjs").write_text("", encoding="utf-8")
        (scripts / "hidden.spec.js").write_text("", encoding="utf-8")
        done = self.run_check("scripts/*.test.mjs")
        self.assertEqual(done.returncode, 1, done.stdout + done.stderr)
        self.assertIn("hidden.spec.js", done.stdout)

    def test_git_repository_scan_reports_tests_outside_the_script_glob(self):
        subprocess.run(["git", "-C", str(self.tree), "init", "--quiet"], check=True)
        scripts = self.tree / "scripts"
        scripts.mkdir()
        (scripts / "active.test.mjs").write_text("", encoding="utf-8")
        (self.tree / "tools" / "stray.test.mjs").write_text("", encoding="utf-8")
        done = self.run_check("scripts/*.test.mjs")
        self.assertEqual(done.returncode, 1, done.stdout + done.stderr)
        self.assertIn("tools/stray.test.mjs", done.stdout)
        self.assertIn("unreachable", done.stdout)

    def test_python_glob_with_a_wildcard_directory_passes(self):
        source = self.tree / "tools" / "sample_tests" / "test_sample.py"
        source.write_text("", encoding="utf-8")
        done = self.run_check("tools/*_tests/test_*.py")
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertIn("tools/sample_tests/test_sample.py", done.stdout)
        self.assertIn("1 test file(s) discovered", done.stdout)

    def test_javascript_glob_with_a_wildcard_directory_scans_a_real_parent(self):
        done = self.run_check("tools/*_tests/*.test.mjs")
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        hidden = self.tree / "tools" / "sample_tests" / "nested"
        hidden.mkdir()
        (hidden / "hidden.test.mjs").write_text("", encoding="utf-8")
        done = self.run_check("tools/*_tests/*.test.mjs")
        self.assertEqual(done.returncode, 1, done.stdout + done.stderr)
        self.assertIn("hidden.test.mjs", done.stdout)
        self.assertIn("unreachable", done.stdout)

    def test_an_empty_glob_fails(self):
        done = self.run_check("tools/empty_tests/*.test.mjs")
        self.assertEqual(done.returncode, 1, done.stdout + done.stderr)
        self.assertIn("no test files match", done.stdout)

    def test_a_missing_directory_fails(self):
        done = self.run_check("tools/renamed_tests/*.test.mjs")
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
