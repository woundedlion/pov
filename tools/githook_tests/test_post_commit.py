"""End-to-end controls for the per-worktree, non-failing size-trail hook."""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from git_test_env import isolated_env  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
HOOK = REPO / ".githooks/post-commit"
TOOL = """import json, os, pathlib, sys
pending = pathlib.Path(sys.argv[sys.argv.index('--pending') + 1])
trail = pathlib.Path(sys.argv[sys.argv.index('--trail') + 1])
if os.environ.get('FAIL_SIZE_TRAIL'):
    raise SystemExit(7)
trail.write_text(json.dumps({'pending': str(pending), 'value': pending.read_text()}))
pending.unlink()
"""


class PostCommitHook(unittest.TestCase):
    def setUp(self):
        if shutil.which("sh") is None:
            if os.environ.get("CI"):
                self.fail("POSIX shell required in CI")
            self.skipTest("no POSIX shell")
        tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.repo = self.root / "main"
        self.repo.mkdir()
        self.worktree = self.root / "linked"
        neutral = str(self.root / "empty-config")
        Path(neutral).write_text("", encoding="utf-8")
        self.env = dict(isolated_env(), GIT_CONFIG_GLOBAL=neutral,
                        GIT_CONFIG_SYSTEM=neutral,
                        GIT_AUTHOR_NAME="hook test", GIT_AUTHOR_EMAIL="hook@test.invalid",
                        GIT_COMMITTER_NAME="hook test", GIT_COMMITTER_EMAIL="hook@test.invalid",
                        HS_PYTHON=Path(sys.executable).as_posix())
        self.git("init", "--quiet", "-b", "master")
        (self.repo / "tools").mkdir()
        (self.repo / "tools/teensy_size_trail.py").write_text(TOOL, encoding="utf-8")
        self.git("add", "--", "tools/teensy_size_trail.py")
        self.git("commit", "--quiet", "-m", "fixture")
        self.git("worktree", "add", "--quiet", "-b", "linked", str(self.worktree))
        gitdir = self.git("rev-parse", "--absolute-git-dir", cwd=self.worktree).stdout.strip()
        self.pending = Path(gitdir) / "teensy-size-trail.pending.json"
        self.main_pending = self.repo / ".git/teensy-size-trail.pending.json"
        self.trail = self.repo / ".git/teensy-size-trail.tsv"
        self.pending.write_text("linked capture", encoding="utf-8")
        self.main_pending.write_text("main capture", encoding="utf-8")

    def git(self, *args, cwd=None):
        return subprocess.run(["git", "-C", str(cwd or self.repo), *args],
                              env=self.env, capture_output=True, text=True, check=True)

    def run_hook(self, **extra):
        return subprocess.run(["sh", HOOK.as_posix()], cwd=self.worktree,
                              env=dict(self.env, **extra), capture_output=True,
                              text=True, check=False)

    def test_linked_capture_uses_private_pending_and_common_trail(self):
        done = self.run_hook()
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertFalse(self.pending.exists())
        self.assertEqual(self.main_pending.read_text(encoding="utf-8"), "main capture")
        stamp = json.loads(self.trail.read_text(encoding="utf-8"))
        self.assertEqual(Path(stamp["pending"]), self.pending)
        self.assertEqual(stamp["value"], "linked capture")
        self.assertFalse((self.pending.parent / self.trail.name).exists())

    def test_failed_capture_keeps_pending_and_exits_zero(self):
        done = self.run_hook(FAIL_SIZE_TRAIL="1")
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertIn("size trail not updated", done.stdout + done.stderr)
        self.assertEqual(self.pending.read_text(encoding="utf-8"), "linked capture")
        self.assertFalse(self.trail.exists())
        self.assertEqual(self.run_hook().returncode, 0)
        self.assertTrue(self.trail.exists())

    def test_absent_tool_removes_only_private_pending(self):
        tool = self.worktree / "tools/teensy_size_trail.py"
        tool.unlink()
        done = self.run_hook()
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertFalse(self.pending.exists())
        self.assertTrue(self.main_pending.exists())
        self.assertFalse(self.trail.exists())
        tool.write_text(TOOL, encoding="utf-8")
        self.pending.write_text("fresh capture", encoding="utf-8")
        self.assertEqual(self.run_hook().returncode, 0)
        self.assertEqual(json.loads(self.trail.read_text(encoding="utf-8"))["value"],
                         "fresh capture")

    def test_unsupported_python_preserves_pending_and_exits_zero(self):
        python = self.root / "old-python"
        python.write_text('#!/bin/sh\n[ "$1" = --version ] && exit 0\nexit 1\n',
                          encoding="utf-8", newline="\n")
        python.chmod(0o755)
        done = self.run_hook(HS_PYTHON=python.as_posix())
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertIn("Python 3.11 or newer is required", done.stdout + done.stderr)
        self.assertTrue(self.pending.exists())
        self.assertFalse(self.trail.exists())
        self.assertEqual(self.run_hook().returncode, 0)
        self.assertTrue(self.trail.exists())

    def test_no_pending_is_quiet_and_does_not_stamp_another_worktree(self):
        self.pending.unlink()
        done = self.run_hook()
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertEqual(done.stdout + done.stderr, "")
        self.assertTrue(self.main_pending.exists())
        self.assertFalse(self.trail.exists())
        self.pending.write_text("new capture", encoding="utf-8")
        self.assertEqual(self.run_hook().returncode, 0)
        self.assertTrue(self.trail.exists())
