#!/usr/bin/env python3
"""End-to-end tests for the glob-discovered test-suite non-empty guard."""

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from git_test_env import isolated_env, require_tool

REPO = Path(__file__).resolve().parents[2]
CHECK = REPO / "tools" / "require_test_files.sh"
# The MSYS runtime expands wildcard argv entries before Bash sees them, which
# would split the glob into several arguments on Windows. Ignored elsewhere.
ENV = {**isolated_env(), "MSYS": "noglob"}


class RequireTestFiles(unittest.TestCase):
    def setUp(self):
        require_tool(self, "bash")
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

    def test_an_empty_glob_fails(self):
        done = self.run_check("tools/empty_tests/*.test.mjs")
        self.assertEqual(done.returncode, 1, done.stdout + done.stderr)
        self.assertIn("no test files match", done.stdout)


if __name__ == "__main__":
    unittest.main()
