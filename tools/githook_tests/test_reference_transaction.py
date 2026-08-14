#!/usr/bin/env python3
"""Host tests for .githooks/reference-transaction fast-forward enforcement.

The hook is what makes the landing model safe: it refuses every
non-fast-forward move of refs/heads/master, and an intentional rewind needs a
one-shot token naming the exact new commit. A regression there silently permits
a clobber of landed work on the one branch every session lands on. The hook is
driven here against a scratch repository -- end to end through `git update-ref`,
and directly against crafted transaction lines for the cases git will not
produce on demand.

Run:  python -m unittest discover -s tools/githook_tests
"""

import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
HOOK = REPO / ".githooks" / "reference-transaction"
ZERO = "0" * 40
MASTER = "refs/heads/master"


def git_env(repo: Path) -> dict[str, str]:
    """Environment that keeps the developer's git config out of the fixture."""
    neutral = str(repo / "empty-config")
    return dict(
        os.environ,
        GIT_CONFIG_GLOBAL=neutral,
        GIT_CONFIG_SYSTEM=neutral,
        GIT_AUTHOR_NAME="hook test",
        GIT_AUTHOR_EMAIL="hook@test.invalid",
        GIT_COMMITTER_NAME="hook test",
        GIT_COMMITTER_EMAIL="hook@test.invalid",
    )


class ReferenceTransactionHook(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.addCleanup(tmp.cleanup)
        self.repo = Path(tmp.name)
        (self.repo / "empty-config").write_text("", encoding="utf-8")
        self.git("init", "--quiet", "-b", "master")

        # base <- one is the fast-forward line; side diverges from base, so it
        # is neither an ancestor nor a descendant of one.
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
        return subprocess.run(
            ["git", *args], cwd=self.repo, env=git_env(self.repo),
            capture_output=True, text=True, check=check)

    def rev(self, ref: str) -> str:
        return self.git("rev-parse", "--verify", ref).stdout.strip()

    def hook(self, *lines: str, phase: str = "prepared"):
        """Feed transaction lines to the hook the way git does.

        stdin is written as bytes: text mode rewrites each newline to CRLF on
        Windows, which leaves a stray CR on the ref name the hook matches.
        """
        done = subprocess.run(
            ["sh", HOOK.as_posix(), phase], cwd=self.repo,
            env=git_env(self.repo),
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

    def test_rewind_of_master_is_refused(self):
        done = self.git("update-ref", MASTER, self.base, check=False)
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("non-fast-forward", done.stderr)
        self.assertEqual(self.rev(MASTER), self.one)
        self.assertIn(f"REFUSED master {self.one} -> {self.base}",
                      self.log_text())

    def test_sideways_move_of_master_is_refused(self):
        done = self.git("update-ref", MASTER, self.side, check=False)
        self.assertNotEqual(done.returncode, 0)
        self.assertEqual(self.rev(MASTER), self.one)

    def test_a_commit_on_master_is_allowed(self):
        self.assertEqual(
            self.git("commit", "--allow-empty", "--quiet", "-m", "next",
                     check=False).returncode, 0)
        self.assertEqual(self.rev(MASTER), self.rev("HEAD"))

    def test_another_branch_may_move_backwards(self):
        self.assertEqual(self.git("update-ref", "refs/heads/ahead", self.base,
                                  check=False).returncode, 0)
        self.assertEqual(self.rev("refs/heads/ahead"), self.base)
        self.assertEqual(self.log_text(), "")

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

    def test_a_zero_token_authorizes_deletion_and_recreation_is_free(self):
        self.token.write_text(f"{ZERO}\n", encoding="utf-8")
        self.assertEqual(self.git("update-ref", "-d", MASTER,
                                  check=False).returncode, 0)
        self.assertEqual(self.git("rev-parse", "--verify", "-q", MASTER,
                                  check=False).returncode, 1)
        # Creating a ref has nothing to fast-forward from, so no token is spent.
        self.assertEqual(self.git("update-ref", MASTER, self.side,
                                  check=False).returncode, 0)
        self.assertEqual(self.rev(MASTER), self.side)

    def test_deletion_of_master_is_refused_without_a_token(self):
        done = self.git("update-ref", "-d", MASTER, check=False)
        self.assertNotEqual(done.returncode, 0)
        self.assertEqual(self.rev(MASTER), self.one)

    # ── Directly, for the lines git will not produce on demand ───────────────

    def test_an_abbreviated_token_authorizes_its_own_commit(self):
        self.token.write_text(self.base[:10], encoding="utf-8")
        self.assertEqual(self.hook(f"{ZERO} {self.base} {MASTER}").returncode, 0)
        self.assertFalse(self.token.exists())

    def test_a_token_for_another_commit_does_not_authorize(self):
        self.token.write_text(self.side, encoding="utf-8")
        done = self.hook(f"{ZERO} {self.base} {MASTER}")
        self.assertEqual(done.returncode, 1)
        self.assertTrue(self.token.exists(), "an unused token must survive")
        self.assertIn("REFUSED master", self.log_text())

    def test_a_token_too_short_to_name_a_commit_does_not_authorize(self):
        self.token.write_text("abc", encoding="utf-8")
        self.assertEqual(self.hook(f"{ZERO} {self.base} {MASTER}").returncode, 1)

    def test_an_empty_token_does_not_authorize(self):
        self.token.write_text("\n", encoding="utf-8")
        self.assertEqual(self.hook(f"{ZERO} {self.base} {MASTER}").returncode, 1)

    def test_a_deletion_token_does_not_authorize_a_rewind(self):
        self.token.write_text(ZERO, encoding="utf-8")
        self.assertEqual(self.hook(f"{ZERO} {self.base} {MASTER}").returncode, 1)
        self.assertTrue(self.token.exists())

    def test_the_any_token_authorizes_any_rewind(self):
        self.token.write_text("ANY\n", encoding="utf-8")
        self.assertEqual(self.hook(f"{ZERO} {self.base} {MASTER}").returncode, 0)
        self.assertFalse(self.token.exists())
        self.assertIn("OVERRIDE master", self.log_text())

    def test_the_stdin_old_field_is_not_trusted(self):
        # git reports a zero old OID to this hook, and a caller could report any
        # value; the fast-forward check reads the ref itself.
        self.assertEqual(
            self.hook(f"{self.base} {self.base} {MASTER}").returncode, 1)

    def test_a_malformed_transaction_line_is_refused(self):
        done = self.hook(f"{ZERO} {self.one}")
        self.assertEqual(done.returncode, 1)
        self.assertIn("malformed", done.stderr)
        self.assertIn("REFUSED malformed", self.log_text())

    def test_a_refused_master_line_fails_the_whole_transaction(self):
        done = self.hook(f"{ZERO} {self.base} refs/heads/other",
                         f"{ZERO} {self.base} {MASTER}",
                         f"{ZERO} {self.base} refs/heads/ahead")
        self.assertEqual(done.returncode, 1)

    def test_no_phase_but_prepared_is_inspected(self):
        for phase in ("committed", "aborted"):
            with self.subTest(phase=phase):
                done = self.hook(f"{ZERO} {self.base} {MASTER}", phase=phase)
                self.assertEqual(done.returncode, 0)
                self.assertEqual(self.log_text(), "")


if __name__ == "__main__":
    unittest.main()
