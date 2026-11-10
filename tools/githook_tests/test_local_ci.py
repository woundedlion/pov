#!/usr/bin/env python3
"""Tests for the local CI job table and the pre-push hook that drives it."""

import os
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

    def test_only_local_only_table_jobs_may_be_absent_from_the_workflow(self):
        with tempfile.TemporaryDirectory() as directory:
            workflow = Path(directory) / "ci.yml"
            workflow.write_text(
                (REPO / ".github/workflows/ci.yml").read_text(encoding="utf-8")
                .replace("\n  lint:\n", "\n  renamed-lint:\n"),
                encoding="utf-8")
            done = self.local_ci("--check", HS_CI_WORKFLOW=str(workflow))
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("table jobs absent from ci.yml: lint", done.stderr)
        self.assertNotIn("windows-tests", done.stderr)


class CommitJobs(unittest.TestCase):
    def setUp(self):
        self.bash = require_tool(self, "bash")
        tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.repo = self.root / "main"
        self.repo.mkdir()
        self.tree = self.root / "ci"
        self.env = dict(isolated_env(), HS_LOCAL_CI_DIR=self.tree.as_posix())
        git_run(self.repo, "init", "--quiet", env=self.env)
        (self.repo / "tools").mkdir()
        shutil.copyfile(LOCAL_CI, self.repo / "tools/local_ci.sh")
        workflow = self.repo / ".github/workflows"
        workflow.mkdir(parents=True)
        shutil.copyfile(REPO / ".github/workflows/ci.yml", workflow / "ci.yml")
        (self.repo / ".gitignore").write_text(".pio/\nbuild/\n", encoding="utf-8")
        (self.repo / "package-lock.json").write_text("{}\n", encoding="utf-8")
        git_run(self.repo, "add", "tools/local_ci.sh", ".github/workflows/ci.yml",
                ".gitignore", "package-lock.json", env=self.env)
        git_run(self.repo, "commit", "--quiet", "-m", "fixture", env=self.env)
        self.sha = git_run(self.repo, "rev-parse", "HEAD", env=self.env).stdout.strip()
        git_run(self.repo, "worktree", "add", "--quiet", "--detach", str(self.tree),
                self.sha, env=self.env)
        stamp = self.tree / "build/local-ci/package-lock.blob"
        stamp.parent.mkdir(parents=True)
        stamp.write_text(git_run(self.repo, "rev-parse", "HEAD:package-lock.json",
                                 env=self.env).stdout, encoding="utf-8")
        self.passed = self.repo / ".git/hs-local-ci-passed"
        bin_dir = self.root / "bin"
        bin_dir.mkdir()
        just = bin_dir / "just"
        just.write_text('#!/usr/bin/env bash\n[ "$1" = lint ] || exit 1\n'
                        'git rev-parse HEAD > build/ran-job\n', encoding="utf-8", newline="\n")
        just.chmod(0o755)
        self.env["PATH"] = str(bin_dir) + os.pathsep + self.env["PATH"]

    def run_ci(self, *args):
        return subprocess.run([self.bash, "tools/local_ci.sh", "--commit", *args],
                              cwd=self.repo, env=self.env, capture_output=True,
                              text=True, check=False)

    def test_named_jobs_ignore_passed_record_and_fill_missing_libdeps(self):
        self.passed.write_text(self.sha + "\n", encoding="utf-8")
        for tree, env_name, content in ((self.repo, "existing", "main"),
                                       (self.repo, "missing", "copied"),
                                       (self.tree, "existing", "retained")):
            directory = tree / ".pio/libdeps" / env_name
            directory.mkdir(parents=True)
            (directory / "library").write_text(content, encoding="utf-8")
        (self.tree / "dirty").write_text("clean me", encoding="utf-8")
        done = self.run_ci(self.sha, "--", "lint")
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertEqual((self.tree / "build/ran-job").read_text().strip(), self.sha)
        self.assertFalse((self.tree / "dirty").exists())
        self.assertEqual((self.tree / ".pio/libdeps/missing/library").read_text(), "copied")
        self.assertEqual((self.tree / ".pio/libdeps/existing/library").read_text(), "retained")
        self.assertEqual(self.passed.read_text(), self.sha + "\n")

    def test_named_jobs_run_each_commit_without_creating_passed_record(self):
        (self.repo / "marker").write_text("next", encoding="utf-8")
        git_run(self.repo, "add", "marker", env=self.env)
        git_run(self.repo, "commit", "--quiet", "-m", "next", env=self.env)
        tip = git_run(self.repo, "rev-parse", "HEAD", env=self.env).stdout.strip()
        done = self.run_ci(self.sha, tip, "--", "lint")
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertEqual(done.stdout.count("PASS lint"), 2)
        self.assertEqual((self.tree / "build/ran-job").read_text().strip(), tip)
        self.assertFalse(self.passed.exists())

    def test_full_runs_still_skip_passed_commits(self):
        self.passed.write_text(self.sha + "\n", encoding="utf-8")
        done = self.run_ci(self.sha)
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertIn("already passed local CI", done.stdout)
        self.assertFalse((self.tree / "build/ran-job").exists())

    def test_missing_commits_or_jobs_and_unknown_jobs_fail(self):
        for args in ((), ("--", "lint"), (self.sha, "--"), (self.sha, "--", "unknown")):
            with self.subTest(args=args):
                done = self.run_ci(*args)
                self.assertEqual(done.returncode, 2, done.stdout + done.stderr)


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
        (tools / "local_ci.sh").write_text('echo "LOCAL_CI $*"\n[ "$1" != "${FAIL_MODE:-}" ]\n',
                                           encoding="utf-8", newline="\n")
        shutil.copyfile(PRE_PUSH, self.repo / "pre-push")

    def push(self, refs: str, **extra: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run([self.shell, "pre-push", "origin", "url"], cwd=self.repo,
                              input=refs, env=dict(self.env, **extra), capture_output=True, text=True,
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

    def test_a_failing_guard_or_commit_check_fails_the_push(self):
        refs = f"refs/heads/x {'a' * 40} refs/heads/x {'0' * 40}\n"
        for mode in ("--new-guards", "--commit"):
            with self.subTest(mode=mode):
                self.assertNotEqual(self.push(refs, FAIL_MODE=mode).returncode, 0)


if __name__ == "__main__":
    unittest.main()
