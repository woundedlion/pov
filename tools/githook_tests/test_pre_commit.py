#!/usr/bin/env python3
"""End-to-end tests for the staged pre-commit gates."""

import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
HOOK = REPO / ".githooks" / "pre-commit"
LOCAL_GIT_ENV = {
    "GIT_ALTERNATE_OBJECT_DIRECTORIES",
    "GIT_COMMON_DIR",
    "GIT_CONFIG",
    "GIT_CONFIG_COUNT",
    "GIT_CONFIG_PARAMETERS",
    "GIT_DIR",
    "GIT_GRAFT_FILE",
    "GIT_IMPLICIT_WORK_TREE",
    "GIT_INDEX_FILE",
    "GIT_NO_REPLACE_OBJECTS",
    "GIT_OBJECT_DIRECTORY",
    "GIT_PREFIX",
    "GIT_REPLACE_REF_BASE",
    "GIT_SHALLOW_FILE",
    "GIT_WORK_TREE",
}


def isolated_env() -> dict[str, str]:
    env = {key: value for key, value in os.environ.items()
           if key not in LOCAL_GIT_ENV}
    env.update({
        "GIT_CONFIG_GLOBAL": os.devnull,
        "GIT_CONFIG_SYSTEM": os.devnull,
        "GIT_AUTHOR_NAME": "fixture",
        "GIT_AUTHOR_EMAIL": "fixture@example.invalid",
        "GIT_COMMITTER_NAME": "fixture",
        "GIT_COMMITTER_EMAIL": "fixture@example.invalid",
        "PYTHON": sys.executable,
    })
    return env


class PreCommitHook(unittest.TestCase):
    def setUp(self):
        if shutil.which("sh") is None:
            self.skipTest("no POSIX shell")
        tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.addCleanup(tmp.cleanup)
        self.repo = Path(tmp.name)
        self.env = isolated_env()
        self.git("init", "--quiet")

        tools = self.repo / "tools"
        tools.mkdir()
        (self.repo / "README.md").write_bytes(b"valid\n")
        (tools / "docs_check.py").write_text(
            "import pathlib, sys\n"
            "root = pathlib.Path(sys.argv[sys.argv.index('--root') + 1])\n"
            "raise SystemExit('BROKEN' in (root / 'README.md').read_text())\n",
            encoding="utf-8")
        for name in ["docs_images.py", "build_pins.py", "license_check.py"]:
            (tools / name).write_text("raise SystemExit(0)\n", encoding="utf-8")
        self.git("add", "README.md", "tools")
        self.git("commit", "--quiet", "-m", "base")

    def git(self, *args: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            ["git", *args], cwd=self.repo, env=self.env,
            capture_output=True, text=True, check=True)

    def run_hook(self, **extra: str) -> subprocess.CompletedProcess[str]:
        env = dict(self.env, GIT_DIR=str(self.repo / ".git"),
                   GIT_WORK_TREE=str(self.repo), **extra)
        return subprocess.run(
            ["sh", HOOK.as_posix()], cwd=self.repo, env=env,
            capture_output=True, text=True, check=False)

    def test_documentation_reads_the_index(self):
        readme = self.repo / "README.md"
        readme.write_bytes(b"BROKEN\n")
        self.git("add", "README.md")
        readme.write_bytes(b"valid working tree\n")
        self.assertNotEqual(self.run_hook().returncode, 0)

        self.git("add", "README.md")
        readme.write_bytes(b"BROKEN working tree\n")
        done = self.run_hook()
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)

    def test_whitespace_errors_fail_before_documentation(self):
        (self.repo / "README.md").write_bytes(b"trailing  \n")
        self.git("add", "README.md")
        done = self.run_hook()
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("staged whitespace errors", done.stdout + done.stderr)

    def test_clang_format_reads_the_index(self):
        formatter = self.repo / "clang-format"
        formatter.write_text(
            "#!/bin/sh\n"
            "case \"${1:-}\" in --version) echo 'clang-format version 22.1.8'; exit 0;; esac\n"
            "grep -q UNFORMATTED && exit 1\n"
            "exit 0\n",
            encoding="utf-8")
        formatter.chmod(0o755)
        source = self.repo / "sample.cpp"
        source.write_bytes(b"int good;\n")
        self.git("add", "sample.cpp")
        self.git("commit", "--quiet", "-m", "source")

        source.write_bytes(b"int UNFORMATTED;\n")
        self.git("add", "sample.cpp")
        source.write_bytes(b"int good_working_tree;\n")
        done = self.run_hook(CLANG_FORMAT=formatter.as_posix())
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("clang-format failed", done.stdout + done.stderr)

    def test_an_unreadable_staged_blob_fails_the_lint(self):
        bin_dir = self.repo / "fakebin"
        bin_dir.mkdir()
        ruff = bin_dir / "ruff"
        ruff.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        ruff.chmod(0o755)

        source = self.repo / "sample.py"
        source.write_bytes(b"x = 1\n")
        self.git("add", "sample.py")
        blob = self.git("rev-parse", ":sample.py").stdout.strip()
        loose = self.repo / ".git" / "objects" / blob[:2] / blob[2:]
        if not loose.exists():
            self.skipTest("staged blob is not a loose object")
        loose.chmod(0o644)
        loose.unlink()

        done = self.run_hook(
            PATH=os.pathsep.join([str(bin_dir), self.env["PATH"]]))
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("cannot read staged sample.py",
                      done.stdout + done.stderr)


if __name__ == "__main__":
    unittest.main()
