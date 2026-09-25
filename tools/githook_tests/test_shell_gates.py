"""Behavioral fixtures for shell selection, build, and profiling gates."""

import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[2]


class ShellGateTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(ignore_cleanup_errors=True)
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        (self.root / "stubs").mkdir()
        self.env = dict(os.environ, GIT_CONFIG_GLOBAL=str(self.root / "no-config"),
                        GIT_CONFIG_SYSTEM=str(self.root / "no-config"))
        self.git("init", "--quiet")
        self.git("config", "core.autocrlf", "false")

    def git(self, *args):
        return subprocess.run(["git", "-C", str(self.root), *args], env=self.env,
                              check=True, capture_output=True)

    def stub(self, name, body):
        path = self.root / "stubs" / name
        path.write_text("#!/usr/bin/env bash\n" + body + "\n", encoding="utf-8", newline="\n")
        path.chmod(0o755)

    def gate(self, name, *args, script=None):
        return subprocess.run(
            ["bash", "-c", 'export PATH="$PWD/stubs:$PATH"; exec bash "$@"',
             "fixture", str(script or REPO / "tools" / name), *args],
            cwd=self.root, env=self.env, capture_output=True, text=True, timeout=30)

    def test_eol_checks_worktree_bytes_and_refuses_an_empty_selection(self):
        empty = self.gate("eol_gate.sh")
        self.assertNotEqual(empty.returncode, 0)
        self.assertIn("no tracked file", empty.stderr)
        (self.root / ".gitattributes").write_text("*.txt text eol=lf\n", encoding="utf-8")
        payload = self.root / "payload.txt"
        payload.write_bytes(b"line\n")
        self.git("add", "--", ".gitattributes", "payload.txt")
        good = self.gate("eol_gate.sh")
        self.assertEqual(good.returncode, 0, good.stdout + good.stderr)
        payload.write_bytes(b"line\r\n")
        bad = self.gate("eol_gate.sh")
        self.assertNotEqual(bad.returncode, 0)
        self.assertIn("worktree line endings are crlf", bad.stdout)

    def test_crlf_attribute_requires_lf_index_and_crlf_checkout(self):
        (self.root / ".gitattributes").write_text("*.txt text eol=crlf\n", encoding="utf-8")
        payload = self.root / "payload.txt"
        payload.write_bytes(b"line\r\n")
        self.git("add", "--", ".gitattributes", "payload.txt")
        self.assertEqual(self.git("show", ":payload.txt").stdout, b"line\n")
        good = self.gate("eol_gate.sh")
        self.assertEqual(good.returncode, 0, good.stdout + good.stderr)
        blob = self.git("hash-object", "-w", "--no-filters", "payload.txt").stdout.decode().strip()
        self.git("update-index", "--cacheinfo", f"100644,{blob},payload.txt")
        bad = self.gate("eol_gate.sh")
        self.assertNotEqual(bad.returncode, 0)
        self.assertIn("index line endings are crlf, expected lf", bad.stdout)

    def test_eol_repair_preserves_unstaged_edits_and_index(self):
        (self.root / ".gitattributes").write_text("*.txt text eol=lf\n", encoding="utf-8")
        payload = self.root / "payload.txt"
        payload.write_bytes(b"staged\n")
        self.git("add", "--", ".gitattributes", "payload.txt")
        payload.write_bytes(b"staged\r\nuser edit\n")
        repaired = self.gate("eol_gate.sh", "--fix-worktree")
        self.assertEqual(repaired.returncode, 0, repaired.stdout + repaired.stderr)
        self.assertEqual(payload.read_bytes(), b"staged\nuser edit\n")
        self.assertEqual(self.git("show", ":payload.txt").stdout, b"staged\n")

    def test_ruff_selection_includes_failed_lints_but_not_empty_reports(self):
        (self.root / "source.py").touch()
        self.git("add", "--", "source.py")
        self.stub("ruff", "exit 1")
        self.assertNotEqual(self.gate("ruff_selection_guard.sh").returncode, 0)
        self.stub("ruff", "echo source.py; exit 1")
        selected = self.gate("ruff_selection_guard.sh")
        self.assertEqual(selected.returncode, 0, selected.stdout + selected.stderr)
        self.stub("ruff", "echo unrelated.py; exit 1")
        self.assertNotEqual(self.gate("ruff_selection_guard.sh").returncode, 0)
        self.stub("ruff", "echo source.py; exit 1")
        (self.root / "omitted.py").touch()
        self.git("add", "--", "omitted.py")
        self.assertNotEqual(self.gate("ruff_selection_guard.sh").returncode, 0)

    def test_eslint_selection_requires_a_file_path(self):
        (self.root / "source.js").touch()
        self.git("add", "--", "source.js")
        self.stub("npx", "echo '[]'")
        self.assertNotEqual(self.gate("eslint_selection_guard.sh").returncode, 0)

        self.stub("npx", "echo '[{\"filePath\":\"source.js\"}]'; exit 1")
        selected = self.gate("eslint_selection_guard.sh")
        self.assertEqual(selected.returncode, 0, selected.stdout + selected.stderr)
        self.stub("npx", "echo '[{\"filePath\":\"unrelated.js\"}]'; exit 1")
        self.assertNotEqual(self.gate("eslint_selection_guard.sh").returncode, 0)
        self.stub("npx", "echo '[{\"filePath\":\"source.js\"}]'; exit 1")
        (self.root / "omitted.js").touch()
        self.git("add", "--", "omitted.js")
        self.assertNotEqual(self.gate("eslint_selection_guard.sh").returncode, 0)

    def test_eslint_selection_reports_a_failed_lint_run(self):
        for payload in ("", "not-json", "{}"):
            self.stub("npx", f"echo '{payload}'; exit 2")
            result = self.gate("eslint_selection_guard.sh")
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("the lint run itself failed", result.stderr)

    def test_clang_format_requires_version_sources_and_success(self):
        self.stub("clang-format", "echo 'clang-format version 21.0.0'")
        wrong = self.gate("clang_format_gate.sh")
        self.assertNotEqual(wrong.returncode, 0)
        self.assertIn("clang-format 22 is required", wrong.stderr)
        self.stub("clang-format", "if [ \"$1\" = --version ]; then echo 'clang-format version 22.0.0'; else printf '%s\\n' \"$@\"; fi")
        empty = self.gate("clang_format_gate.sh")
        self.assertNotEqual(empty.returncode, 0)
        self.assertIn("no files selected", empty.stdout)
        generated = self.root / "core" / "color" / "gamut_lut.h"
        generated.parent.mkdir(parents=True)
        generated.touch()
        (self.root / "selected file.cpp").touch()
        self.git("add", "--", ".")
        selected = self.gate("clang_format_gate.sh")
        self.assertEqual(selected.returncode, 0, selected.stderr)
        self.assertIn("selected file.cpp", selected.stdout)
        self.assertNotIn("gamut_lut.h", selected.stdout)
        self.assertIn("--dry-run\n--Werror\n--style=file", selected.stdout)
        self.stub("clang-format", "if [ \"$1\" = --version ]; then echo 'clang-format version 22.0.0'; else exit 9; fi")
        self.assertNotEqual(self.gate("clang_format_gate.sh").returncode, 0)

    def test_shellcheck_refuses_empty_selection_and_propagates_lint_failure(self):
        empty = self.gate("shellcheck_gate.sh")
        self.assertNotEqual(empty.returncode, 0)
        self.assertIn("no shell files selected", empty.stdout)
        (self.root / "selected.sh").write_text("#!/bin/sh\nexit 0\n", encoding="utf-8", newline="\n")
        self.git("add", "--", "selected.sh")
        self.stub("shellcheck", "printf '%s\\n' \"$*\"; exit 9")
        bad = self.gate("shellcheck_gate.sh")
        self.assertNotEqual(bad.returncode, 0)
        self.assertIn("selected.sh", bad.stdout)

    def test_composite_steps_are_isolated_and_unsupported_scalars_fail(self):
        if shutil.which("shellcheck") is None:
            self.skipTest("shellcheck is not installed")
        (self.root / "selected.sh").write_text("#!/bin/sh\nexit 0\n", encoding="utf-8", newline="\n")
        action = self.root / ".github/actions/fixture/action.yml"
        action.parent.mkdir(parents=True)
        self.git("add", "--", "selected.sh")
        for run in ('|\n        echo "$value"', 'echo ok', '>\n        echo ok'):
            action.write_text('runs:\n  using: composite\n  steps:\n'
                              '    - shell: bash\n      run: |\n'
                              '        value=ok; echo "$value"\n'
                              '    - shell: sh\n      run: ' + run + '\n',
                              encoding="utf-8", newline="\n")
            self.git("add", "--", ".github/actions/fixture/action.yml")
            result = self.gate("shellcheck_gate.sh")
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
            if run.startswith("|"):
                self.assertIn("SC2154", result.stdout + result.stderr)

    def test_cold_build_removes_only_fixture_caches_and_preserves_pio_failure(self):
        script = self.root / "tools" / "teensy_cold_build.sh"
        script.parent.mkdir()
        shutil.copyfile(REPO / "tools" / script.name, script)
        for name in ("build", "build_cache"):
            directory = self.root / ".pio" / name
            directory.mkdir(parents=True)
            (directory / "stale").write_text("cached", encoding="utf-8")
        sentinel = self.root / ".pio" / "keep"
        sentinel.write_text("keep", encoding="utf-8")
        self.stub("pio", "echo fixture-build-failure; exit 7")
        failed = self.gate(script.name, "capture.log", script=script)
        self.assertEqual(failed.returncode, 7, failed.stdout + failed.stderr)
        self.assertFalse((self.root / ".pio" / "build").exists())
        self.assertFalse((self.root / ".pio" / "build_cache").exists())
        self.assertEqual(sentinel.read_text(encoding="utf-8"), "keep")
        self.assertIn("fixture-build-failure", (self.root / "capture.log").read_text())

    def test_profile_sweep_detects_a_missing_playlist_member_without_hardware(self):
        tools = self.root / "tools"
        tools.mkdir()
        script = tools / "profile_sweep.sh"
        shutil.copyfile(REPO / "tools" / script.name, script)
        playlist = self.root / "targets" / "Phantasm" / "phantasm_playlist.h"
        playlist.parent.mkdir(parents=True)
        shutil.copyfile(REPO / "targets" / "Phantasm" / playlist.name, playlist)
        valid = self.gate(script.name, "check", script=script)
        self.assertEqual(valid.returncode, 0, valid.stdout + valid.stderr)
        text = playlist.read_text(encoding="utf-8")
        text = text.replace("#define HS_PHANTASM_EFFECT_LIST(X)",
                            "#define HS_PHANTASM_EFFECT_LIST(X) \\\n  X(FixtureMissingEffect, 1) ", 1)
        playlist.write_text(text, encoding="utf-8", newline="\n")
        bad = self.gate(script.name, "check", script=script)
        self.assertNotEqual(bad.returncode, 0)
        self.assertIn("FixtureMissingEffect", bad.stderr)


if __name__ == "__main__":
    unittest.main()
