#!/usr/bin/env python3
"""Host tests for the build-pin gate (tools/build_pins.py).

Run:  python -m unittest discover -s tools/build_pins_tests
"""

import contextlib
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
import unittest.mock
from pathlib import Path
from types import SimpleNamespace

TOOLS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TOOLS))

from git_test_env import git_run, isolated_env  # noqa: E402
import build_pins as bp  # noqa: E402


class AuthorityPresence(unittest.TestCase):
    def test_missing_spelling_in_one_authority_is_rejected(self):
        for check, relative, original, replacement in (
                (bp.check_shared_literals, 'tools/commit_checks.sh',
                 "grep -vE '", "grep -v -E '"),
                (bp.check_shared_literals, 'tools/commit_checks.sh',
                 'core.whitespace=', 'core.whitespace '),
                (bp.check_inline_pins, 'requirements/ruff.in',
                 'ruff==', 'ruff >= ')):
            path = bp.ROOT / relative
            read = bp.read_scanned
            def changed(candidate, errors):
                text = read(candidate, errors)
                return text.replace(original, replacement) if candidate == path else text
            with self.subTest(path=relative), unittest.mock.patch.object(
                    bp, 'read_scanned', side_effect=changed):
                self.assertTrue(any('missing' in error and relative in error
                                    for error in check()))


class RequirementPins(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        (self.root / "tools").mkdir()
        (self.root / "requirements").mkdir()
        for name in ("build_pins.py", "teensy_gate.py"):
            shutil.copyfile(TOOLS / name, self.root / "tools" / name)
        for path in (bp.ROOT / "requirements").glob("*.in"):
            shutil.copyfile(path, self.root / "requirements" / path.name)
        self.source = self.root / "requirements/ruff.in"
        self.lock = self.root / "requirements/ruff.txt"

    def _run(self, *args, **kwargs):
        return subprocess.run(
            [sys.executable, str(self.root / "tools/build_pins.py"), *args],
            capture_output=True, text=True, **kwargs)

    def _check_pair(self):
        uses = tuple(use for use in bp.INLINE_USES if use[1] == "ruff")
        with unittest.mock.patch.object(bp, "ROOT", self.root), \
                unittest.mock.patch.object(bp, "INLINE_SCAN",
                                           (self.source, self.lock)), \
                unittest.mock.patch.object(bp, "INLINE_USES", uses), \
                unittest.mock.patch.dict(bp.PINS,
                                         ruff=self._run("ruff").stdout.strip()):
            return bp.check_inline_pins()

    def test_a_coordinated_dependency_bump_needs_no_script_edit(self):
        self.source.write_text("ruff==99.0.0\n", encoding="utf-8")
        self.lock.write_text("ruff==99.0.0 \\\n    --hash=sha256:abc\n",
                             encoding="utf-8")
        result = self._run("ruff")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), "99.0.0")
        self.assertEqual(self._check_pair(), [])
        output = self.root / "github-output"
        result = self._run("--github-output",
                           env=dict(os.environ, GITHUB_OUTPUT=str(output)))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("ruff=99.0.0\n", output.read_text(encoding="utf-8"))

    def test_a_stale_lock_still_fails(self):
        self.source.write_text("ruff==99.0.0\n", encoding="utf-8")
        self.lock.write_text("ruff==98.0.0\n", encoding="utf-8")
        errors = self._check_pair()
        self.assertEqual(len(errors), 1)
        self.assertIn("requirements/ruff.txt", errors[0].replace("\\", "/"))
        self.assertIn("99.0.0", errors[0])
        self.assertIn("98.0.0", errors[0])


class DuplicatesPin(unittest.TestCase):
    """An injected pin re-spelled in a build file is a second source of truth."""

    def test_value_on_the_naming_line_is_a_duplicate(self):
        self.assertTrue(bp.duplicates_pin(
            "        node-version: 24.13.0\n", "node", "24.13.0"))

    def test_a_commented_out_value_is_not_a_duplicate(self):
        self.assertFalse(bp.duplicates_pin(
            "        node: latest   # was 24.13.0\n", "node", "24.13.0"))


class EngineRanges(unittest.TestCase):
    """package.json declares the floor; the pin is what CI installs."""

    def _check(self, engines, pin):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "package.json"
            manifest.write_text(json.dumps({"engines": engines}),
                                encoding="utf-8")
            with unittest.mock.patch.object(bp, "ROOT", root), \
                    unittest.mock.patch.object(bp, "ENGINE_RANGES", ((manifest, ("engines", "node"), "node"),)), \
                    unittest.mock.patch.object(bp, "PINS", dict(bp.PINS, node=pin)):
                return bp.check_engine_ranges()

    def test_pin_above_the_floor_passes(self):
        self.assertEqual(self._check({"node": ">=22"}, "24.13.0"), [])

    def test_pin_below_the_floor_is_reported(self):
        errors = self._check({"node": ">=26"}, "24.13.0")
        self.assertEqual(len(errors), 1)
        self.assertIn("24.13.0", errors[0])

    def test_a_range_the_gate_cannot_compare_is_reported(self):
        # A caret/tilde/OR range silently stops declaring a floor; the gate must
        # say so rather than skip the manifest.
        for spec in ("^22", "~24.13.0", ">=22 <25", "*"):
            errors = self._check({"node": spec}, "24.13.0")
            self.assertEqual(len(errors), 1, spec)
            self.assertIn("expected '>=X'", errors[0])


class FastMathFlags(unittest.TestCase):
    def test_fast_math_axis_marks_its_test_contract(self):
        definitions = [match.groups() for flag in bp.FAST_MATH_TEST_FLAGS
                       if (match := re.fullmatch(r"-D(\w+)=(\d+)", flag))]
        self.assertTrue(definitions)
        harness = (bp.ROOT / "tests/test_harness.h").read_text(encoding="utf-8")
        for name, value in definitions:
            self.assertNotEqual(int(value), 0)
            self.assertRegex(harness, rf"defined\s*\(\s*{re.escape(name)}\s*\)")


class InlinePins(unittest.TestCase):
    def test_repeated_consistent_pins_need_no_count_update(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "ci.yml"
            pattern = next(row for row in bp.INLINE_USES if row[1] == "python")
            with unittest.mock.patch.object(bp, "ROOT", root), \
                    unittest.mock.patch.object(bp, "INLINE_SCAN", (path,)), \
                    unittest.mock.patch.object(bp, "INLINE_USES", (pattern,)), \
                    unittest.mock.patch.dict(bp.INLINE_PINS, python="3.12"):
                for count in (1, 3):
                    path.write_text("python-version: '3.12'\n" * count,
                                    encoding="utf-8")
                    self.assertEqual(bp.check_inline_pins(), [])
                path.write_text("python-version: '3.13'\n", encoding="utf-8")
                self.assertIn("written '3.13'", bp.check_inline_pins()[0])
                path.write_text("", encoding="utf-8")
                self.assertIn("was not found", bp.check_inline_pins()[0])


class CheckTool(unittest.TestCase):
    """--check-tool holds PATH to the pin, so it must be able to reach it."""

    def _check(self, name, stdout):

        def run(command, **kwargs):
            if stdout is None:
                raise OSError("not found")
            return subprocess.CompletedProcess(command, 0, stdout, "")

        with unittest.mock.patch.object(bp.subprocess, "run", run):
            with contextlib.redirect_stdout(io.StringIO()) as out:
                status = bp.check_tool(name)
        return status, out.getvalue()

    def test_clang_format_honors_the_binary_override(self):

        with unittest.mock.patch.dict(bp.os.environ, CLANG_FORMAT="/custom/clang-format"), \
                unittest.mock.patch.object(bp.subprocess, "run") as run:
            run.return_value = SimpleNamespace(stdout="clang-format version " + bp.INLINE_PINS["clang-format"])
            self.assertEqual(bp.check_tool("clang-format"), 0)
        self.assertEqual(run.call_args.args[0], ["/custom/clang-format", "--version"])

    def test_a_dependency_bump_updates_the_runtime_version_check(self):
        with unittest.mock.patch.dict(bp.PINS, ruff="99.0.0"):
            self.assertEqual(self._check("ruff", "ruff 99.0.0")[0], 0)
            status, message = self._check("ruff", "ruff 98.0.0")
        self.assertEqual(status, 1)
        self.assertIn("pip install --require-hashes -r requirements/ruff.txt",
                      message)

    def test_a_missing_tool_reports_how_to_install_that_tool(self):
        for name, want in (("just",
                            "pip install --require-hashes -r requirements/just.txt"),
                           ("shellcheck",
                            "pip install --require-hashes -r "
                            "requirements/shellcheck.txt"),
                           ("node", f"install Node {bp.PINS['node']}"),
                           ("doxygen", f"install Doxygen {bp.INLINE_PINS['doxygen']}")):
            status, message = self._check(name, None)
            self.assertEqual(status, 1, name)
            self.assertIn("nothing runnable", message)
            self.assertIn(want, message)


class InstallSet(unittest.TestCase):
    """Repository sources copied into the sibling Daydream checkout."""

    def test_an_unpinned_file_is_reported(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            git_run(root, 'init', '--quiet')
            attributes = root / ".gitattributes"
            attributes.write_text(
                "* text=auto eol=lf\nunpinned/file.txt !eol !binary\n",
                encoding="utf-8", newline="\n")
            with unittest.mock.patch.object(bp, "ROOT", root), \
                    unittest.mock.patch.dict(os.environ, isolated_env(), clear=True):
                self.assertEqual(bp.check_install_eol(["pinned/file.txt"]), [])
                self.assertEqual(len(bp.check_install_eol(["unpinned/file.txt"])), 1)
                attributes.write_text("* text=auto eol=lf\n", encoding="utf-8", newline="\n")
                self.assertEqual(bp.check_install_eol(["unpinned/file.txt"]), [])

    def test_an_empty_set_fails_instead_of_passing_vacuously(self):
        self.assertTrue(bp.check_install_eol([]))


class FlexRamGeometry(unittest.TestCase):
    """The budgets, gate and linker script spell one FlexRAM bank geometry."""

    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        (self.root / "tools").mkdir()
        for name in ("teensy_gate.py", "teensy_budgets.json", "phantasm.ld"):
            shutil.copyfile(bp.ROOT / "tools" / name, self.root / "tools" / name)
        self.linker = self.root / "tools/phantasm.ld"
        self.gate = self.root / "tools/teensy_gate.py"

    def check(self):
        with unittest.mock.patch.object(bp, "ROOT", self.root):
            return bp.check_flexram_geometry()

    def rewrite(self, path, text):
        path.write_text(text, encoding="utf-8", newline="")

    def test_the_committed_tree_agrees(self):
        self.assertEqual(self.check(), [])

    def test_a_lowercased_hex_literal_is_the_same_geometry(self):
        text = self.linker.read_text(encoding="utf-8")
        self.rewrite(self.linker, re.sub(
            r"0[xX][0-9A-Fa-f]+", lambda found: found.group(0).lower(), text))
        self.assertEqual(self.check(), [])

    def test_a_reflowed_expression_is_the_same_geometry(self):
        text = self.linker.read_text(encoding="utf-8")
        self.rewrite(self.linker, text.replace(") >> ", ")\n\t>> "))
        self.assertEqual(self.check(), [])

    def test_a_linker_script_without_the_geometry_is_reported(self):
        self.rewrite(self.linker, "MEMORY { }\n")
        errors = self.check()
        self.assertEqual(len(errors), 2)
        self.assertTrue(all("phantasm.ld" in error for error in errors), errors)

    def test_a_gate_bank_size_that_disagrees_is_reported(self):
        text = self.gate.read_text(encoding="utf-8")
        self.rewrite(self.gate, re.sub(
            r"^FLEXRAM_BANK_BYTES = .*$", "FLEXRAM_BANK_BYTES = 1", text,
            flags=re.MULTILINE))
        errors = self.check()
        self.assertTrue(any("FlexRAM bank size differs" in error
                            for error in errors), errors)


if __name__ == "__main__":
    unittest.main()
