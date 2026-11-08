#!/usr/bin/env python3
"""Host tests for .githooks/reference-transaction fast-forward enforcement.

The hook refuses every non-fast-forward move of refs/heads/master unless a
one-shot token names the target commit (full or abbreviated SHA), the current
commit for deletion, or ANY.

Run:  python -m unittest discover -s tools/githook_tests
"""

import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from git_test_env import git_run, isolated_env

REPO = Path(__file__).resolve().parents[2]
HOOK = REPO / ".githooks" / "reference-transaction"
ZERO = "0" * 40
MASTER = "refs/heads/master"


class ReferenceTransactionHook(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.addCleanup(tmp.cleanup)
        self.repo = Path(tmp.name)
        self.git("init", "--quiet", "-b", "master")

        # base <- one <- two (branch ahead); side branches from base.
        self.git("commit", "--allow-empty", "--quiet", "-m", "base")
        self.base = self.rev("HEAD")
        self.git("commit", "--allow-empty", "--quiet", "-m", "one")
        self.one = self.rev("HEAD")
        self.git("checkout", "--quiet", "-b", "ahead")
        self.git("commit", "--allow-empty", "--quiet", "-m", "two")
        self.two = self.rev("HEAD")
        self.git("checkout", "--quiet", "-b", "other", self.base)
        self.git("commit", "--allow-empty", "--quiet", "-m", "side")
        self.side = self.rev("HEAD")
        self.git("checkout", "--quiet", "master")

        # Installed last: the fixture commits above must not be gated.
        hooks = self.repo / ".git" / "hooks"
        hooks.mkdir(exist_ok=True)
        self.installed = hooks / "reference-transaction"
        shutil.copy(HOOK, self.installed)
        os.chmod(self.installed, 0o755)
        self.git("config", "core.hooksPath", str(hooks))
        self.token = self.repo / ".git" / "hs-allow-nonff"
        self.log = self.repo / ".git" / "hs-nonff.log"

    def git(self, *args, check=True):
        return git_run(self.repo, *args, check=check)

    def rev(self, ref: str) -> str:
        return self.git("rev-parse", "--verify", ref).stdout.strip()

    def hook(self, *lines: str, phase: str = "prepared"):
        """Feed transaction lines to the hook the way git does.

        stdin is written as bytes: text mode rewrites each newline to CRLF on
        Windows, which leaves a stray CR on the ref name the hook matches.
        """
        done = subprocess.run(
            ["sh", HOOK.as_posix(), phase], cwd=self.repo,
            env=isolated_env(),
            input="".join(f"{line}\n" for line in lines).encode("utf-8"),
            capture_output=True)
        return subprocess.CompletedProcess(
            done.args, done.returncode,
            done.stdout.decode("utf-8", "replace"),
            done.stderr.decode("utf-8", "replace"))

    def log_text(self) -> str:
        return self.log.read_text(encoding="utf-8") if self.log.exists() else ""

    # ── End to end, through git's own transaction ────────────────────────────

    def test_fast_forward_of_master_is_allowed(self):
        self.assertEqual(self.git("update-ref", MASTER, self.two,
                                  check=False).returncode, 0)
        self.assertEqual(self.rev(MASTER), self.two)
        self.assertEqual(self.log_text(), "")

    def test_deletion_of_packed_master_is_refused(self):
        self.git("pack-refs", "--all")
        done = self.git("update-ref", "-d", MASTER, check=False)
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("REFUSED deletion", done.stderr)
        self.assertEqual(self.rev(MASTER), self.one)

    def test_rewind_of_master_is_refused(self):
        done = self.git("update-ref", MASTER, self.base, check=False)
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("non-fast-forward", done.stderr)
        self.assertEqual(self.rev(MASTER), self.one)
        self.assertIn(f"REFUSED master {self.one} -> {self.base}",
                      self.log_text())

    def test_an_exact_token_authorizes_one_rewind(self):
        self.token.write_text(f"{self.base}\n", encoding="utf-8")
        self.assertEqual(self.git("update-ref", MASTER, self.base,
                                  check=False).returncode, 0)
        self.assertEqual(self.rev(MASTER), self.base)
        self.assertFalse(self.token.exists(), "the token must be consumed")
        self.assertIn(f"OVERRIDE master {self.one} -> {self.base}",
                      self.log_text())
        # Consumed: the next non-fast-forward move is refused again.
        self.assertNotEqual(self.git("update-ref", "-d", MASTER,
                                     check=False).returncode, 0)

    # ── Directly, for the lines git will not produce on demand ───────────────

    def test_a_token_for_another_commit_does_not_authorize(self):
        self.token.write_text(self.side, encoding="utf-8")
        done = self.hook(f"{ZERO} {self.base} {MASTER}")
        self.assertEqual(done.returncode, 1)
        self.assertTrue(self.token.exists(), "an unused token must survive")
        self.assertIn("REFUSED master", self.log_text())

    def test_the_stdin_old_field_is_not_trusted(self):
        # git reports a zero old OID to this hook, and a caller could report any
        # value; the fast-forward check reads the ref itself.
        self.assertEqual(
            self.hook(f"{self.base} {self.base} {MASTER}").returncode, 1)


if __name__ == "__main__":
    unittest.main()
