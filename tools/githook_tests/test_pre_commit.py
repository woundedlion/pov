#!/usr/bin/env python3
"""End-to-end tests for the staged pre-commit gates and the range-only whole-tree gates."""

import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from git_test_env import git_run, isolated_env, require_tool

REPO = Path(__file__).resolve().parents[2]
HOOK = REPO / ".githooks" / "pre-commit"


class PreCommitHook(unittest.TestCase):
    def setUp(self):
        self.shell = require_tool(self, "sh")
        git = require_tool(self, "git")
        tmp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.addCleanup(tmp.cleanup)
        self.repo = Path(tmp.name)
        self.env = isolated_env()
        self.git("init", "--quiet")

        tools = self.repo / "tools"
        tools.mkdir()
        (self.repo / "README.md").write_bytes(b"valid\n")
        (tools / "docs_check.py").write_text(
            "import pathlib, subprocess, sys\n"
            "root = pathlib.Path(sys.argv[sys.argv.index('--root') + 1])\n"
            f"tracked = subprocess.run([{git!r}, '-C', str(root), 'ls-files'],\n"
            "                         check=True, capture_output=True, text=True)\n"
            "print('TRACKED', tracked.stdout.split())\n"
            "raise SystemExit('BROKEN' in (root / 'README.md').read_text())\n",
            encoding="utf-8")
        (self.repo / "tests").mkdir()
        gates = [(tools, name) for name in ["docs_images.py", "build_pins.py", "license_check.py",
                                            "include_paths.py", "encoding_check.py"]]
        for directory, name in gates + [(self.repo / "tests", "check_death_pins.py")]:
            (directory / name).write_text(
                "import os, sys\n"
                "staged = sys.stdin.read().split() if '--staged' in sys.argv else None\n"
                f"print('GATE {name}', staged)\n"
                f"raise SystemExit(int(os.environ.get('FAIL_{name.removesuffix('.py').upper()}', '0')))\n",
                encoding="utf-8")
        requirements = self.repo / "requirements"
        requirements.mkdir()
        shutil.copyfile(REPO / "requirements" / "ruff.txt", requirements / "ruff.txt")
        shutil.copyfile(REPO / "ruff.toml", self.repo / "ruff.toml")
        shutil.copyfile(REPO / ".clang-format", self.repo / ".clang-format")
        self.git("add", "README.md", "tools", "tests", "requirements", "ruff.toml", ".clang-format")
        self.git("commit", "--quiet", "-m", "base")

    def git(self, *args: str) -> subprocess.CompletedProcess[str]:
        return git_run(self.repo, *args, env=self.env)

    def run_hook(self, **extra: str) -> subprocess.CompletedProcess[str]:
        env = dict(self.env, GIT_DIR=str(self.repo / ".git"),
                   GIT_WORK_TREE=str(self.repo), **extra)
        return subprocess.run(
            [self.shell, HOOK.as_posix()], cwd=self.repo, env=env,
            capture_output=True, text=True, check=False)

    def run_range(self, spec="HEAD~1..HEAD", **extra: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [self.shell, (REPO / "tools" / "commit_checks.sh").as_posix(), spec],
            cwd=self.repo, env=dict(self.env, **extra),
            capture_output=True, text=True, check=False)

    def commit_all(self, message: str):
        self.git("commit", "--quiet", "-m", message)

    def formatter(self):
        formatter = self.repo / "clang-format"
        formatter.write_text(
            "#!/bin/sh\n"
            "grep -q UNFORMATTED && exit 1\n"
            "exit 0\n",
            encoding="utf-8")
        formatter.chmod(0o755)
        return formatter

    def test_staged_gates_get_the_staged_paths_and_failures_propagate(self):
        (self.repo / "README.md").write_bytes(b"valid updated\n")
        (self.repo / "unstaged.md").write_bytes(b"unstaged\n")
        self.git("add", "README.md")
        for gate in ("license_check", "include_paths", "encoding_check"):
            with self.subTest(gate=gate):
                clean = self.run_hook()
                self.assertEqual(clean.returncode, 0, clean.stdout + clean.stderr)
                self.assertIn(f"GATE {gate}.py ['README.md']", clean.stdout)
                failed = self.run_hook(**{f"FAIL_{gate.upper()}": "7"})
                self.assertNotEqual(failed.returncode, 0)

    def test_whole_tree_gates_run_over_a_range_not_per_commit(self):
        (self.repo / "README.md").write_bytes(b"BROKEN\n")
        self.git("add", "README.md")
        staged = self.run_hook()
        self.assertEqual(staged.returncode, 0, staged.stdout + staged.stderr)
        for gate in ("TRACKED", "GATE docs_images.py", "GATE build_pins.py"):
            self.assertNotIn(gate, staged.stdout)
        self.commit_all("broken docs")
        self.assertNotEqual(self.run_range().returncode, 0)

        (self.repo / "README.md").write_bytes(b"valid\n")
        self.git("add", "README.md")
        self.commit_all("fixed docs")
        for gate in ("docs_images", "build_pins", "license_check", "include_paths", "encoding_check",
                     "check_death_pins"):
            with self.subTest(gate=gate):
                clean = self.run_range()
                self.assertEqual(clean.returncode, 0, clean.stdout + clean.stderr)
                self.assertIn(f"GATE {gate}.py None", clean.stdout)
                self.assertNotEqual(self.run_range(**{f"FAIL_{gate.upper()}": "7"}).returncode, 0)

    def test_python_uses_staged_lint_configuration(self):
        require_tool(self, "ruff")
        config = self.repo / "ruff.toml"
        config.write_text('lint.select = ["F401"]\n', encoding="utf-8", newline="\n")
        self.git("add", "ruff.toml")
        self.git("commit", "--quiet", "-m", "strict config")
        (self.repo / "sample.py").write_text("import os\n", encoding="utf-8", newline="\n")
        self.git("add", "sample.py")
        config.write_text('lint.ignore = ["F401"]\n', encoding="utf-8", newline="\n")
        done = self.run_hook()
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("F401", done.stdout + done.stderr)

    def test_python_lints_temporary_commit_indexes(self):
        require_tool(self, "ruff")
        source = self.repo / "sample.py"
        source.write_bytes(b"value = 1\n")
        self.git("add", "sample.py")
        self.git("commit", "--quiet", "-m", "source")
        self.git("config", "core.hooksPath", HOOK.parent.as_posix())
        source.write_bytes(b"import os\n")
        for args in [("-a",), ("--", "sample.py")]:
            with self.subTest(args=args):
                done = subprocess.run(
                    ["git", "-C", str(self.repo), "commit", "-m", "invalid", *args],
                    env=self.env, capture_output=True, text=True, check=False)
                self.assertNotEqual(done.returncode, 0)
                self.assertIn("F401", done.stdout + done.stderr)

    def profile_fixture(self):
        require_tool(self, "node")
        scripts = self.repo / "scripts"
        scripts.mkdir()
        for name in ["check_profiles.mjs", "effect_roster.mjs", "exit.mjs"]:
            shutil.copyfile(REPO / "scripts" / name, scripts / name)
        targets = self.repo / "targets"
        (targets / "Phantasm").mkdir(parents=True)
        (targets / "effects.h").write_text(
            "#define HS_EFFECT_LIST(X) X(Example) X(Galaxies)\n",
            encoding="utf-8", newline="\n")
        (targets / "Phantasm" / "phantasm_playlist.def").write_text(
            "show(Example, 1)\nshow(Galaxies, 1)\n",
            encoding="utf-8", newline="\n")
        profiles = self.repo / "docs" / "profiles"
        reports = {name: f"profile_{name.lower()}_teensy_2026-09-28.md"
                   for name in ("Example", "Galaxies")}
        for timing_set in ["shipping", "O3"]:
            directory = profiles / timing_set
            directory.mkdir(parents=True)
            for name, filename in reports.items():
                (directory / filename).write_text(
                    f"# {name} on-device profile (2026-09-28)\n\n"
                    "## Setup\n\nHardware.\n\n## Frame cadence\n\nCadence.\n\n"
                    "## Summary ranking\n\nRanking.\n\n## Harness\n\nHarness.\n",
                    encoding="utf-8", newline="\n")
            (directory / "README.md").write_text(
                "Reports covering the 2 effects.\n\n" + "".join(
                    f"[{name}]({filename})\n" for name, filename in reports.items()),
                encoding="utf-8", newline="\n")
        (profiles / "README.md").write_text(
            "On-device timing for the **2 effects in the Phantasm image**\n\n" + "".join(
                f"[{timing_set} {name}]({timing_set}/{filename})\n"
                for timing_set in ("shipping", "O3")
                for name, filename in reports.items()),
            encoding="utf-8", newline="\n")
        self.git("add", "scripts", "targets", "docs/profiles")
        self.git("commit", "--quiet", "-m", "profile fixture")
        return profiles

    def test_profile_archive_is_validated_as_of_the_range_tip(self):
        profiles = self.profile_fixture()
        index = profiles / "README.md"
        valid = index.read_bytes()
        broken = valid.replace(b"shipping/profile_", b"shipping/missing_")
        index.write_bytes(broken)
        self.git("add", "docs/profiles/README.md")
        staged = self.run_hook()
        self.assertEqual(staged.returncode, 0, staged.stdout + staged.stderr)
        self.assertNotIn("Profile archive structure", staged.stdout)
        self.commit_all("broken profile index")
        index.write_bytes(valid)
        done = self.run_range()
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("main shipping index is missing", done.stdout + done.stderr)

        index.write_bytes(valid + b"\nCapture notes.\n")
        self.git("add", "docs/profiles/README.md")
        self.commit_all("profile notes")
        index.write_bytes(broken)
        done = self.run_range()
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertIn("Profile archive structure checks all 2", done.stdout)

    def test_whitespace_errors_fail_before_the_staged_gates(self):
        (self.repo / "README.md").write_bytes(b"trailing  \n")
        self.git("add", "README.md")
        done = self.run_hook()
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("staged whitespace errors", done.stdout + done.stderr)
        self.assertNotIn("GATE", done.stdout)

    def test_clang_format_reads_the_index(self):
        formatter = self.formatter()
        source = self.repo / "sample.cpp"
        source.write_bytes(b"int good;\n")
        self.git("add", "sample.cpp")
        done = self.run_hook(CLANG_FORMAT=formatter.as_posix())
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.git("commit", "--quiet", "-m", "source")

        source.write_bytes(b"int UNFORMATTED;\n")
        self.git("add", "sample.cpp")
        source.write_bytes(b"int good_working_tree;\n")
        done = self.run_hook(CLANG_FORMAT=formatter.as_posix())
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("clang-format failed", done.stdout + done.stderr)

    def check_missing_linter(self, tool, name, content):
        directory = self.repo / "linter-bin"
        directory.mkdir()
        for command in ("sh", "git", "dirname", "mktemp", "rm", "mkdir",
                        "cp", "mv", "grep", "sed"):
            resolved = subprocess.run(
                [self.shell, "-c", 'command -v "$1"', "fixture", command],
                env=self.env, capture_output=True, text=True, check=True).stdout.strip()
            proxy = directory / command
            proxy.write_text(f'#!/bin/sh\nexec {shlex.quote(resolved)} "$@"\n',
                             encoding="utf-8", newline="\n")
            proxy.chmod(0o755)
        (self.repo / name).write_text(content, encoding="utf-8", newline="\n")
        self.git("add", name)
        result = self.run_hook(PATH=str(directory))
        self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn(f"{tool} not found", result.stdout + result.stderr)

        linter = directory / tool
        linter.write_text('#!/bin/sh\necho FIXTURE_LINTER_RAN\n'
                          'exit "${FIXTURE_LINTER_EXIT:-0}"\n',
                          encoding="utf-8", newline="\n")
        linter.chmod(0o755)
        present = self.run_hook(PATH=str(directory))
        self.assertEqual(present.returncode, 0, present.stdout + present.stderr)
        self.assertIn("FIXTURE_LINTER_RAN", present.stdout)
        failed = self.run_hook(PATH=str(directory), FIXTURE_LINTER_EXIT="7")
        self.assertEqual(failed.returncode, 7, failed.stdout + failed.stderr)
        self.assertIn("FIXTURE_LINTER_RAN", failed.stdout)

    def test_missing_ruff_fails_the_commit(self):
        self.check_missing_linter("ruff", "sample.py", "value = 1\n")

    def eslint_stub(self):
        path = self.repo / "node_modules" / ".bin" / "eslint"
        path.parent.mkdir(parents=True)
        path.write_text('#!/bin/sh\nprintf "%s\\n" "$*" >> "$ESLINT_LOG"\n'
                        'last=; for arg do last=$arg; done\n'
                        '[ "$last" != . ] || exit "${ESLINT_TREE_FAIL:-0}"\n'
                        'if grep -q ESLINT_REJECT; then exit 1; fi\nexit 0\n',
                        encoding="utf-8", newline="\n")
        path.chmod(0o755)
        return self.repo / "eslint.log"

    def test_javascript_lint_reads_the_staged_blob(self):
        log = self.eslint_stub()
        source = self.repo / "sample.mjs"
        source.write_text("ESLINT_REJECT\n", encoding="utf-8", newline="\n")
        self.git("add", "sample.mjs")
        source.write_text("export const value = 1;\n", encoding="utf-8", newline="\n")
        self.assertNotEqual(self.run_hook(ESLINT_LOG=log.as_posix()).returncode, 0)
        self.git("add", "sample.mjs")
        source.write_text("ESLINT_REJECT\n", encoding="utf-8", newline="\n")
        done = self.run_hook(ESLINT_LOG=log.as_posix())
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)

    def test_eslint_configuration_change_lints_the_whole_tree_over_a_range(self):
        log = self.eslint_stub()
        (self.repo / "eslint.config.mjs").write_text("export default [];\n", encoding="utf-8", newline="\n")
        self.git("add", "eslint.config.mjs")
        staged = self.run_hook(ESLINT_LOG=log.as_posix(), ESLINT_TREE_FAIL="7")
        self.assertEqual(staged.returncode, 0, staged.stdout + staged.stderr)
        self.commit_all("eslint config")
        failed = self.run_range(ESLINT_LOG=log.as_posix(), ESLINT_TREE_FAIL="7")
        self.assertNotEqual(failed.returncode, 0)
        self.assertTrue(any(line == "." for line in log.read_text().splitlines()))
        passed = self.run_range(ESLINT_LOG=log.as_posix())
        self.assertEqual(passed.returncode, 0, passed.stdout + passed.stderr)

    def test_shell_lint_reads_the_index(self):
        require_tool(self, "shellcheck")
        script = self.repo / "sample.sh"
        script.write_text('#!/bin/sh\nprintf "%s" $1\n', encoding="utf-8", newline="\n")
        self.git("add", "sample.sh")
        script.write_text('#!/bin/sh\nprintf "%s" "$1"\n', encoding="utf-8", newline="\n")
        done = self.run_hook()
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("SC2086", done.stdout + done.stderr)
        self.git("add", "sample.sh")
        script.write_text('#!/bin/sh\nprintf "%s" $1\n', encoding="utf-8", newline="\n")
        done = self.run_hook()
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)

    def test_an_unreadable_staged_blob_fails_the_snapshot(self):
        source = self.repo / "sample.py"
        source.write_bytes(b"x = 1\n")
        self.git("add", "sample.py")
        blob = self.git("rev-parse", ":sample.py").stdout.strip()
        loose = self.repo / ".git" / "objects" / blob[:2] / blob[2:]
        if not loose.exists():
            self.skipTest("staged blob is not a loose object")
        loose.chmod(0o644)
        loose.unlink()

        done = self.run_hook()
        self.assertNotEqual(done.returncode, 0)
        self.assertIn("cannot materialize staged snapshot",
                      done.stdout + done.stderr)


if __name__ == "__main__":
    unittest.main()
