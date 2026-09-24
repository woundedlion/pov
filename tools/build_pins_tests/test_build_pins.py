#!/usr/bin/env python3
"""Host tests for the build-pin gate (tools/build_pins.py).

The gate is the single source for every externally-installed version, and each
of its three checks reads a foreign file shape: duplicates_pin scans workflow
YAML two lines at a time, check_engine_ranges parses package.json's `>=X`
string, and _version_tuple compares versions of unequal width. A drift in any
of those shapes makes the check detect nothing while still printing PASS.

The install-set check reads a fourth shape, CMakeLists.txt's install() rules; a
rule it stops recognising silently exempts those files from their line-ending
pin. The FlexRAM check reads a fifth, tools/phantasm.ld's derived symbols, and
is the only tie between the budgets, the size gate and the linker script.

Run:  python -m unittest discover -s tools/build_pins_tests
"""

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

TOOLS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TOOLS))

import build_pins as bp  # noqa: E402


class AuthorityPresence(unittest.TestCase):
    def test_missing_spelling_in_one_authority_is_rejected(self):
        for check, relative, original, replacement in (
                (bp.check_shared_literals, '.githooks/pre-commit',
                 "grep -vE '", "grep -v -E '"),
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

    def test_comments_and_blank_lines_are_allowed(self):
        self.source.write_text("# Linter\n\nruff==99.0.0 # exact pin\n",
                               encoding="utf-8")
        result = self._run("ruff")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), "99.0.0")

    def test_malformed_sources_fail_without_a_traceback(self):
        for content in ("", "ruff>=99.0.0", "ruff==99.*", "ruff==",
                        "other==99.0.0", "ruff==99.0.0\nruff==98.0.0"):
            with self.subTest(content=content):
                self.source.write_text(content, encoding="utf-8")
                result = self._run("ruff")
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("expected one exact ruff==VERSION pin",
                              result.stderr)
                self.assertNotIn("Traceback", result.stderr)

    def test_a_missing_source_fails_without_a_traceback(self):
        self.source.unlink()
        result = self._run("ruff")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("ruff.in: cannot be read", result.stderr)
        self.assertNotIn("Traceback", result.stderr)


class DuplicatesPin(unittest.TestCase):
    """An injected pin re-spelled in a build file is a second source of truth."""

    def test_value_on_the_naming_line_is_a_duplicate(self):
        self.assertTrue(bp.duplicates_pin(
            "        node-version: 24.13.0\n", "node", "24.13.0"))

    def test_name_two_lines_above_is_still_in_context(self):
        self.assertTrue(bp.duplicates_pin(
            "      - name: Set up Node\n"
            "        with:\n"
            "          version: 24.13.0\n", "node", "24.13.0"))

    def test_name_three_lines_above_is_out_of_context(self):
        # The window is the two preceding lines plus the line itself; a value
        # further from any mention of the tool is not a pin duplicate.
        self.assertFalse(bp.duplicates_pin(
            "      - name: Set up Node\n"
            "        with:\n"
            "          cache: npm\n"
            "          version: 24.13.0\n", "node", "24.13.0"))

    def test_underscore_alias_matches(self):
        self.assertTrue(bp.duplicates_pin(
            "        clang_format: 22.1.8\n", "clang-format", "22.1.8"))

    def test_a_commented_out_value_is_not_a_duplicate(self):
        self.assertFalse(bp.duplicates_pin(
            "        node: latest   # was 24.13.0\n", "node", "24.13.0"))

    def test_a_longer_version_is_not_the_pin(self):
        for line in ("        node: 24.13.0.1\n", "        node: 124.13.0\n"):
            self.assertFalse(bp.duplicates_pin(line, "node", "24.13.0"), line)

    def test_the_tracked_workflows_duplicate_no_pin(self):
        # The live gate, so a real duplicate fails here as well as in CI.
        for path in sorted((bp.ROOT / ".github/workflows").glob("*.yml")):
            text = path.read_text(encoding="utf-8")
            for name, value in bp.PINS.items():
                self.assertFalse(bp.duplicates_pin(text, name, value),
                                 f"{path.name} duplicates {name} {value}")


class VersionTuple(unittest.TestCase):
    def test_parses_each_dotted_component(self):
        self.assertEqual(bp._version_tuple("24.13.0"), (24, 13, 0))
        self.assertEqual(bp._version_tuple("22"), (22,))

    def test_components_compare_numerically_not_lexically(self):
        self.assertLess(bp._version_tuple("24.9.0"), bp._version_tuple("24.13.0"))
        self.assertGreater(bp._version_tuple("1.10.0"), bp._version_tuple("1.9.0"))


class EngineRanges(unittest.TestCase):
    """package.json declares the floor; the pin is what CI installs."""

    def _check(self, engines, pin):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "package.json"
            manifest.write_text(json.dumps({"engines": engines}),
                                encoding="utf-8")
            saved = (bp.ROOT, bp.ENGINE_RANGES, bp.PINS)
            bp.ROOT = root
            bp.ENGINE_RANGES = ((manifest, ("engines", "node"), "node"),)
            bp.PINS = dict(bp.PINS, node=pin)
            try:
                return bp.check_engine_ranges()
            finally:
                bp.ROOT, bp.ENGINE_RANGES, bp.PINS = saved

    def test_pin_above_the_floor_passes(self):
        self.assertEqual(self._check({"node": ">=22"}, "24.13.0"), [])

    def test_pin_equal_to_the_floor_passes(self):
        self.assertEqual(self._check({"node": ">=24.13.0"}, "24.13.0"), [])

    def test_versions_of_unequal_width_compare_by_component(self):
        # Raw tuples order a prefix below its extension, so an unpadded (24,)
        # would read as below a (24, 0, 0) floor it exactly meets.
        self.assertEqual(self._check({"node": ">=24.0.0"}, "24"), [])
        self.assertEqual(self._check({"node": ">=24"}, "24.13.0"), [])
        # A real shortfall is still caught at every width.
        self.assertEqual(len(self._check({"node": ">=24.14"}, "24.13.0")), 1)
        self.assertEqual(len(self._check({"node": ">=24.0.1"}, "24")), 1)

    def test_pin_below_the_floor_is_reported(self):
        errors = self._check({"node": ">=26"}, "24.13.0")
        self.assertEqual(len(errors), 1)
        self.assertIn("24.13.0", errors[0])

    def test_whitespace_in_the_range_is_tolerated(self):
        self.assertEqual(self._check({"node": ">= 22"}, "24.13.0"), [])

    def test_a_range_the_gate_cannot_compare_is_reported(self):
        # A caret/tilde/OR range silently stops declaring a floor; the gate must
        # say so rather than skip the manifest.
        for spec in ("^22", "~24.13.0", ">=22 <25", "*"):
            errors = self._check({"node": spec}, "24.13.0")
            self.assertEqual(len(errors), 1, spec)
            self.assertIn("expected '>=X'", errors[0])

    def test_a_missing_range_is_reported(self):
        errors = self._check({"npm": ">=10"}, "24.13.0")
        self.assertEqual(len(errors), 1)
        self.assertIn("no engines.node range", errors[0])

    def test_the_tracked_manifest_satisfies_its_pin(self):
        self.assertEqual(bp.check_engine_ranges(), [])


class SharedLiterals(unittest.TestCase):
    """Strings ci.yml, the justfile and the pre-commit hook must spell alike."""

    def test_the_tracked_copies_agree(self):
        self.assertEqual(bp.check_shared_literals(), [])

    def test_every_scanned_file_the_check_counts_is_tracked(self):
        for path in bp.INLINE_SCAN:
            self.assertTrue(path.is_file(), path)

    def test_the_format_exclude_copies_are_where_the_count_expects(self):
        want = bp.SHARED_LITERALS["format-exclude"]
        carriers = [path.name for path in bp.INLINE_SCAN
                    if want in path.read_text(encoding="utf-8")]
        self.assertEqual(sorted(carriers),
                         ["clang_format_gate.sh", "pre-commit"])

    def test_fast_math_axis_marks_its_test_contract(self):
        self.assertEqual(
            bp.FAST_MATH_TEST_FLAGS,
            (*bp.FLOAT_FLAGS, "-DHS_TEST_FAST_MATH=1"),
        )
        self.assertEqual(
            bp.SHARED_LITERALS["float-test-flags"],
            "-ffast-math -fno-finite-math-only -DHS_TEST_FAST_MATH=1",
        )


class InlinePins(unittest.TestCase):
    def test_the_tracked_spellings_agree(self):
        self.assertEqual(bp.check_inline_pins(), [])

    def test_repeated_consistent_pins_need_no_count_update(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = root / "ci.yml"
            pattern = next(row for row in bp.INLINE_USES if row[1] == "python")
            with unittest.mock.patch.object(bp, "ROOT", root), \
                    unittest.mock.patch.object(bp, "INLINE_SCAN", (path,)), \
                    unittest.mock.patch.object(bp, "INLINE_USES", (pattern,)):
                for count in (1, 3):
                    path.write_text("python-version: '3.12'\n" * count,
                                    encoding="utf-8")
                    self.assertEqual(bp.check_inline_pins(), [])
                path.write_text("python-version: '3.13'\n", encoding="utf-8")
                self.assertIn("written '3.13'", bp.check_inline_pins()[0])
                path.write_text("", encoding="utf-8")
                self.assertIn("was not found", bp.check_inline_pins()[0])

    def test_the_python_version_input_is_read_in_either_quoting(self):
        pattern = next(p for p, name, _ in bp.INLINE_USES if name == "python")
        for line in ("          python-version: '3.13'",
                     '          python-version: "3.13"',
                     "          python-version: 3.13"):
            with self.subTest(line=line):
                self.assertEqual(re.findall(pattern, line), ["3.13"])

    def test_the_hook_compares_against_its_pinned_format_major(self):
        lines = (bp.ROOT / ".githooks/pre-commit").read_text(
            encoding="utf-8").splitlines()
        major = bp.INLINE_PINS["clang-format"].split(".")[0]
        self.assertIn(f"HS_CLANG_FORMAT_MAJOR={major}", lines)
        self.assertIn(
            '  if [ "$major" != "$HS_CLANG_FORMAT_MAJOR" ]; then', lines)


class UnreadableScannedFile(unittest.TestCase):
    """A renamed scanned path is an error line, not a traceback out of the hook.

    Every scanned path is named in a table here, so a rename leaves the table
    pointing at nothing; the pre-commit hook runs these checks, and a
    FileNotFoundError there reports no finding at all.
    """

    def test_a_missing_inline_scan_entry_is_reported(self):
        missing = bp.ROOT / "no-such-build-file.yml"
        with unittest.mock.patch.object(bp, "INLINE_SCAN", (missing,)):
            inline = bp.check_inline_pins()
            shared = bp.check_shared_literals()
        self.assertTrue(any("no-such-build-file.yml" in e for e in inline))
        self.assertTrue(any("no-such-build-file.yml" in e for e in shared))

    def test_a_missing_consumer_is_reported(self):
        import contextlib
        import io
        missing = bp.ROOT / "no-such-consumer.yml"
        out = io.StringIO()
        with unittest.mock.patch.object(bp, "CONSUMERS", {missing: ("x",)}), \
                contextlib.redirect_stdout(out):
            status = bp.check_consumers()
        self.assertEqual(status, 1)
        self.assertIn("no-such-consumer.yml", out.getvalue())

    def test_a_manifest_that_is_not_json_is_reported(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "package.json"
            manifest.write_text("{not json", encoding="utf-8")
            with unittest.mock.patch.object(
                    bp, "ENGINE_RANGES",
                    ((manifest, ("engines", "node"), "node"),)), \
                    unittest.mock.patch.object(bp, "ROOT", Path(tmp)):
                errors = bp.check_engine_ranges()
        self.assertEqual(len(errors), 1)
        self.assertIn("not valid JSON", errors[0])

    def test_a_missing_manifest_is_reported(self):
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "package.json"
            with unittest.mock.patch.object(
                    bp, "ENGINE_RANGES",
                    ((manifest, ("engines", "node"), "node"),)), \
                    unittest.mock.patch.object(bp, "ROOT", Path(tmp)):
                errors = bp.check_engine_ranges()
        self.assertEqual(len(errors), 1)
        self.assertIn("cannot be read", errors[0])


class ConsumerCallSites(unittest.TestCase):
    def test_every_build_pin_gate_entry_point_is_required(self):
        paths = (bp.ROOT / ".github/workflows/ci.yml",
                 bp.ROOT / "justfile",
                 bp.ROOT / ".githooks/pre-commit")
        for path in paths:
            calls = [line.strip()
                     for line in path.read_text(encoding="utf-8").splitlines()
                     if not line.lstrip().startswith("#")
                     and "tools/build_pins.py" in line
                     and line.strip().endswith("--check")]
            self.assertEqual(len(calls), 1, path)
            self.assertIn(calls[0], bp.CONSUMERS.get(path, ()), path)


class CheckTool(unittest.TestCase):
    """--check-tool holds PATH to the pin, so it must be able to reach it.

    A pin naming a git ref, a file digest or an SDK has no `--version` to
    compare, and the install command differs per pin: the PyPI distribution of
    `just` is rust-just, of `shellcheck` is shellcheck-py, and clang, Node,
    Doxygen and Python do not come from pip at all.
    """

    def _check(self, name, stdout):
        import contextlib
        import io
        import subprocess

        def run(command, **kwargs):
            if stdout is None:
                raise OSError("not found")
            return subprocess.CompletedProcess(command, 0, stdout, "")

        with unittest.mock.patch.object(bp.subprocess, "run", run):
            with contextlib.redirect_stdout(io.StringIO()) as out:
                status = bp.check_tool(name)
        return status, out.getvalue()

    def test_every_target_names_a_pin(self):
        self.assertEqual(set(bp.CHECK_TOOLS) - set(bp.PINS | bp.INLINE_PINS),
                         set())
        self.assertEqual(bp.check_consumers(), 0)

    def test_justfile_tool_checks_are_supported(self):
        source = (bp.ROOT / "justfile").read_text(encoding="utf-8")
        names = set(re.findall(r"--check-tool ([a-z-]+)", source))
        self.assertTrue(names)
        self.assertEqual(names - bp.CHECK_TOOLS.keys(), set())

    def test_actionlint_checks_the_binary_release_and_command(self):
        pin = bp.PINS["actionlint"]
        self.assertEqual(bp.CHECK_TOOLS["actionlint"][0], ["actionlint", "-version"])
        self.assertEqual(self._check("actionlint", pin.rsplit(".", 1)[0] + "\nbuilt with go")[0], 0)
        self.assertEqual(self._check("actionlint", "0.0.0")[0], 1)

    def test_pins_with_no_version_to_report_are_not_targets(self):
        for name in ("daydream", "doxygen-awesome", "doxygen-sha256",
                     "llvm-key-sha256", "emsdk", "kicad"):
            self.assertIn(name, bp.PINS | bp.INLINE_PINS)
            self.assertNotIn(name, bp.CHECK_TOOLS)

    def test_a_major_only_pin_is_met_by_a_release_of_that_major(self):
        # clang's pin is a major; the binary reports the full version, which
        # an equality test could never satisfy.
        self.assertEqual(
            self._check("clang", "Ubuntu clang version 22.1.8 (tags/x)")[0], 0)
        self.assertEqual(self._check("python", "Python 3.12.9")[0], 0)

    def test_a_different_major_still_fails(self):
        status, message = self._check("clang", "clang version 21.1.0")
        self.assertEqual(status, 1)
        self.assertIn("apt install clang-22", message)

    def test_a_packaging_suffix_is_not_expected_from_the_binary(self):
        # shellcheck-py's version is the release plus a suffix; shellcheck
        # reports the release, so the pin was unsatisfiable by equality.
        pin = bp.PINS["shellcheck"]
        status, message = self._check(
            "shellcheck", "ShellCheck - shell script analysis tool\n"
                          f"version: {pin.rsplit('.', 1)[0]}\n")
        self.assertEqual(status, 0)
        self.assertIn(pin, message)

    def test_a_dependency_bump_updates_the_runtime_version_check(self):
        with unittest.mock.patch.dict(bp.PINS, ruff="99.0.0"):
            self.assertEqual(self._check("ruff", "ruff 99.0.0")[0], 0)
            status, message = self._check("ruff", "ruff 98.0.0")
        self.assertEqual(status, 1)
        self.assertIn("pip install ruff==99.0.0", message)

    def test_a_missing_tool_reports_how_to_install_that_tool(self):
        for name, want in (("just", f"pip install rust-just=={bp.PINS['just']}"),
                           ("shellcheck",
                            f"pip install shellcheck-py=={bp.PINS['shellcheck']}"),
                           ("node", "install Node 24.13.0"),
                           ("doxygen", "install Doxygen 1.17.0")):
            status, message = self._check(name, None)
            self.assertEqual(status, 1, name)
            self.assertIn("nothing runnable", message)
            self.assertIn(want, message)

    def test_the_remediation_is_not_a_blanket_pip_install(self):
        for name in bp.CHECK_TOOLS:
            _, message = self._check(name, "0.0.0")
            pin = (bp.PINS | bp.INLINE_PINS)[name]
            if name in ("clang", "doxygen", "just", "node", "python",
                        "shellcheck"):
                self.assertNotIn(f"pip install {name}=={pin}", message, name)


class InstallSet(unittest.TestCase):
    """Repository sources copied into the sibling Daydream checkout."""

    def test_the_runtime_install_has_each_source_class(self):
        installed = bp.installed_sources()
        for path in ("hardware/pov_segment_map.json",
                     "scripts/shader_workbench.mjs", "scripts/sha256.mjs"):
            self.assertIn(path, installed)
        self.assertTrue(
            any(path.startswith("patterns/") for path in installed))
        self.assertIn("README.md", installed)
        self.assertTrue(
            any(path.startswith("docs/screenshots/") for path in installed))

    def test_documentation_is_part_of_the_simulator_install(self):
        installed = bp.installed_sources()
        self.assertIn("README.md", installed)
        self.assertTrue(
            all(path.endswith(".png") for path in installed
                if path.startswith("docs/screenshots/")))

    def test_a_directory_rule_selects_only_its_patterns(self):
        # patterns/ also holds a README the FILES_MATCHING patterns exclude.
        self.assertNotIn("patterns/README.md", bp.installed_sources())

    def test_a_generated_artifact_is_not_a_repository_source(self):
        # The module, glue, and exported catalog are generated during install.
        self.assertNotIn("scripts/engine_catalog.json", bp.installed_sources())
        for path in bp.installed_sources():
            self.assertTrue((bp.ROOT / path).is_file(), path)

    def test_the_live_install_set_is_pinned(self):
        self.assertEqual(bp.check_install_eol(bp.installed_sources()), [])

    def test_an_unpinned_file_is_reported(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            subprocess.run(["git", "init", "--quiet", str(root)], check=True)
            attributes = root / ".gitattributes"
            attributes.write_text(
                "* text=auto eol=lf\nunpinned/file.txt !eol !binary\n",
                encoding="utf-8", newline="\n")
            with unittest.mock.patch.object(bp, "ROOT", root):
                self.assertEqual(bp.check_install_eol(["pinned/file.txt"]), [])
                self.assertEqual(len(bp.check_install_eol(["unpinned/file.txt"])), 1)
                attributes.write_text("* text=auto eol=lf\n", encoding="utf-8", newline="\n")
                self.assertEqual(bp.check_install_eol(["unpinned/file.txt"]), [])

    def test_an_empty_set_fails_instead_of_passing_vacuously(self):
        self.assertTrue(bp.check_install_eol([]))


class FlexRamGeometry(unittest.TestCase):
    """tools/teensy_budgets.json, tools/teensy_gate.py and tools/phantasm.ld all
    spell one FlexRAM bank geometry, and this check is what ties them."""

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
