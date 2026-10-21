#!/usr/bin/env python3
"""End-to-end tests for commit_checks.sh over a commit range."""

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from git_test_env import git_run, isolated_env, require_tool  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
CHECKS = REPO / "tools" / "commit_checks.sh"


class CommitRangeChecks(unittest.TestCase):
    def setUp(self):
        require_tool(self, "sh")
        tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.addCleanup(tmp.cleanup)
        self.repo = Path(tmp.name) / "repo"
        self.repo.mkdir()
        self.env = isolated_env()
        self.git("init", "--quiet", "-b", "main")
        tools = self.repo / "tools"
        tools.mkdir()
        (self.repo / "README.md").write_bytes(b"valid\n")
        (tools / "docs_check.py").write_text(
            "import pathlib, sys\n"
            "root = pathlib.Path(sys.argv[sys.argv.index('--root') + 1])\n"
            "raise SystemExit('BROKEN' in (root / 'README.md').read_text())\n",
            encoding="utf-8")
        for name in ["docs_images.py", "build_pins.py", "license_check.py", "include_paths.py"]:
            (tools / name).write_text("raise SystemExit(0)\n", encoding="utf-8")
        (self.repo / ".clang-format").write_bytes(b"BasedOnStyle: LLVM\n")
        self.formatter = Path(tmp.name) / "clang-format"
        self.formatter.write_text(
            "#!/bin/sh\ngrep -q UNFORMATTED && exit 1\nexit 0\n",
            encoding="utf-8", newline="\n")
        self.formatter.chmod(0o755)
        self.commit("base", "README.md", "tools", ".clang-format")

    def git(self, *args: str) -> subprocess.CompletedProcess[str]:
        return git_run(self.repo, *args, env=self.env)

    def commit(self, message: str, *paths: str) -> str:
        self.git("add", "--", *paths)
        self.git("commit", "--quiet", "-m", message)
        return self.git("rev-parse", "HEAD").stdout.strip()

    def write(self, name: str, text: str) -> str:
        (self.repo / name).write_text(text, encoding="utf-8", newline="\n")
        return name

    def check(self, *args: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            ["sh", CHECKS.as_posix(), *args], cwd=self.repo,
            env=dict(self.env, CLANG_FORMAT=self.formatter.as_posix()),
            capture_output=True, text=True, check=False)

    def test_an_unformatted_file_in_the_range_fails(self):
        self.commit("good", self.write("good.cpp", "int good;\n"))
        self.commit("bad", self.write("bad.cpp", "int UNFORMATTED;\n"))
        done = self.check("main~2..main")
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("clang-format failed: bad.cpp", done.stdout + done.stderr)
        self.assertNotIn("good.cpp", done.stdout + done.stderr)

    def test_a_formatted_range_passes(self):
        self.commit("good", self.write("good.cpp", "int good;\n"))
        done = self.check("main~1..main")
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertIn("range main~1..main formatting, lint, and documentation checks passed",
                      done.stdout)

    def test_files_are_checked_as_of_the_tip(self):
        self.commit("bad", self.write("sample.cpp", "int UNFORMATTED;\n"))
        self.commit("fixed", self.write("sample.cpp", "int fixed;\n"))
        done = self.check("main~2..main")
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)

    def test_the_index_and_working_tree_are_ignored(self):
        self.commit("good", self.write("sample.cpp", "int good;\n"))
        self.write("sample.cpp", "int UNFORMATTED;\n")
        self.git("add", "sample.cpp")
        self.write("README.md", "BROKEN\n")
        done = self.check("main~1..main")
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertEqual(self.git("diff", "--cached", "--name-only").stdout, "sample.cpp\n")

    def test_tree_wide_checks_run_on_the_tip_snapshot(self):
        self.commit("docs", self.write("README.md", "BROKEN\n"))
        self.assertNotEqual(self.check("main~1..main").returncode, 0)

    def test_a_diverged_base_checks_only_the_tip_side(self):
        fork = self.git("rev-parse", "HEAD").stdout.strip()
        self.commit("peer", self.write("peer.cpp", "int UNFORMATTED;\n"))
        self.git("checkout", "--quiet", "-b", "topic", fork)
        self.commit("topic", self.write("topic.cpp", "int good;\n"))
        done = self.check("main..topic")
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.commit("topic bad", self.write("topic_bad.cpp", "int UNFORMATTED;\n"))
        done = self.check("main..topic")
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("clang-format failed: topic_bad.cpp", done.stdout + done.stderr)

    def test_whitespace_errors_in_the_range_fail(self):
        self.commit("trailing", self.write("notes.md", "trailing  \n"))
        done = self.check("main~1..main")
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("range main~1..main whitespace errors", done.stdout + done.stderr)

    def test_an_empty_range_passes(self):
        done = self.check("main..main")
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertIn("changes no files", done.stdout)

    def test_malformed_ranges_are_rejected(self):
        for spec in ("main", "main..", "..main", "main...main", "main..absent"):
            with self.subTest(spec=spec):
                done = self.check(spec)
                self.assertEqual(done.returncode, 2, done.stdout + done.stderr)


if __name__ == "__main__":
    unittest.main()
