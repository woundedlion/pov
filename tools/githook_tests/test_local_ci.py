#!/usr/bin/env python3
"""Tests for the local CI job table and the pre-push hook that drives it."""

import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from git_test_env import git_run, isolated_env, require_tool

REPO = Path(__file__).resolve().parents[2]
LOCAL_CI = REPO / "tools" / "local_ci.sh"
PRE_PUSH = REPO / ".githooks" / "pre-push"


class JobTable(unittest.TestCase):
    def setUp(self):
        self.bash = require_tool(self, "bash")

    def local_ci(self, *args: str, **env: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run([self.bash, LOCAL_CI.as_posix(), *args], cwd=REPO,
                              env=dict(isolated_env(), **env),
                              capture_output=True, text=True, check=False)

    def test_the_table_names_every_workflow_job(self):
        done = self.local_ci("--check")
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)

    def test_every_listed_runner_is_a_job_function_or_a_skip_reason(self):
        source = LOCAL_CI.read_text(encoding="utf-8")
        rows = self.local_ci("--list").stdout.splitlines()
        self.assertTrue(rows)
        for job, runner in (row.split("|", 1) for row in rows):
            with self.subTest(job=job):
                if runner.startswith("skip:"):
                    self.assertTrue(runner.removeprefix("skip:").strip())
                else:
                    self.assertIn(f"\n{runner}() {{", source)

    def test_a_workflow_job_missing_from_the_table_fails_the_check(self):
        with tempfile.TemporaryDirectory() as directory:
            workflow = Path(directory) / "ci.yml"
            workflow.write_text(
                (REPO / ".github/workflows/ci.yml").read_text(encoding="utf-8")
                + "  brand-new-job:\n    runs-on: ubuntu-24.04\n",
                encoding="utf-8")
            done = self.local_ci("--check", HS_CI_WORKFLOW=str(workflow))
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("brand-new-job", done.stderr)

    def test_an_unknown_job_is_rejected(self):
        done = self.local_ci("no-such-job")
        self.assertEqual(done.returncode, 2)
        self.assertIn("unknown job no-such-job", done.stderr)


class PrePushHook(unittest.TestCase):
    def setUp(self):
        self.shell = require_tool(self, "sh")
        tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.addCleanup(tmp.cleanup)
        self.repo = Path(tmp.name)
        self.env = isolated_env()
        git_run(self.repo, "init", "--quiet", env=self.env)
        tools = self.repo / "tools"
        tools.mkdir()
        (tools / "local_ci.sh").write_text('echo "LOCAL_CI $*"\n', encoding="utf-8", newline="\n")
        shutil.copyfile(PRE_PUSH, self.repo / "pre-push")

    def push(self, refs: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run([self.shell, "pre-push", "origin", "url"], cwd=self.repo,
                              input=refs, env=self.env, capture_output=True, text=True,
                              check=False)

    def test_each_distinct_pushed_commit_is_run_and_deletions_are_skipped(self):
        a, b, zero = "a" * 40, "b" * 40, "0" * 40
        done = self.push(f"refs/heads/x {a} refs/heads/x {zero}\n"
                         f"refs/heads/y {a} refs/heads/y {b}\n"
                         f"(delete) {zero} refs/heads/z {b}\n"
                         f"refs/heads/w {b} refs/heads/w {a}\n")
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertEqual(done.stdout.strip().splitlines(), [
            f"LOCAL_CI --new-guards {a} {zero}",
            f"LOCAL_CI --new-guards {a} {b}",
            f"LOCAL_CI --new-guards {b} {a}",
            f"LOCAL_CI --commit {a} {b}",
        ])

    def test_a_push_of_only_deletions_runs_nothing(self):
        done = self.push(f"(delete) {'0' * 40} refs/heads/z {'b' * 40}\n")
        self.assertEqual(done.returncode, 0, done.stderr)
        self.assertEqual(done.stdout, "")


if __name__ == "__main__":
    unittest.main()
