#!/usr/bin/env python3
"""End-to-end tests for the CI test-suite directory pin check."""

import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
CHECK = REPO / "tools" / "check_test_dir_pins.sh"


class TestDirPins(unittest.TestCase):
    def setUp(self):
        if shutil.which("bash") is None:
            self.skipTest("no Bash shell")
        tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.addCleanup(tmp.cleanup)
        self.repo = Path(tmp.name)
        (self.repo / ".github" / "workflows").mkdir(parents=True)
        tests = self.repo / "tools" / "sample_tests"
        tests.mkdir(parents=True)
        (tests / "test_sample.py").write_text("", encoding="utf-8")
        subprocess.run(["git", "init", "--quiet"], cwd=self.repo, check=True)
        subprocess.run(["git", "add", "tools/sample_tests/test_sample.py"],
                       cwd=self.repo, check=True)

    def run_check(self, workflow: str) -> subprocess.CompletedProcess[str]:
        path = self.repo / ".github" / "workflows" / "ci.yml"
        path.write_text(workflow, encoding="utf-8")
        return subprocess.run(
            ["bash", CHECK.as_posix()], cwd=self.repo,
            capture_output=True, text=True, check=False)

    def test_discovery_step_is_required(self):
        required = (
            "bash tools/require_test_files.sh 'tools/sample_tests/test*.py'\n"
        )
        done = self.run_check(
            required + "python -m unittest discover -s tools/sample_tests -v\n")
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)

        done = self.run_check(required)
        self.assertEqual(done.returncode, 1, done.stdout + done.stderr)
        self.assertIn("tools/sample_tests", done.stdout)

    def test_a_longer_sibling_directory_does_not_satisfy_the_pin(self):
        done = self.run_check(
            "bash tools/require_test_files.sh 'tools/sample_tests/test*.py'\n"
            "python -m unittest discover -s tools/sample_tests_extra -v\n")
        self.assertEqual(done.returncode, 1, done.stdout + done.stderr)
        self.assertIn("tools/sample_tests", done.stdout)


if __name__ == "__main__":
    unittest.main()
