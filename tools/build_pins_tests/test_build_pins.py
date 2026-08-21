#!/usr/bin/env python3
"""Host tests for the build-pin gate (tools/build_pins.py).

The gate is the single source for every externally-installed version, and each
of its three checks reads a foreign file shape: duplicates_pin scans workflow
YAML two lines at a time, check_engine_ranges parses package.json's `>=X`
string, and _version_tuple compares versions of unequal width. A drift in any
of those shapes makes the check detect nothing while still printing PASS.

The install-set check reads a fourth shape, CMakeLists.txt's install() rules; a
rule it stops recognising silently exempts those files from their line-ending
pin.

Run:  python -m unittest discover -s tools/build_pins_tests
"""

import json
import sys
import tempfile
import unittest
import unittest.mock
from pathlib import Path

TOOLS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TOOLS))

import build_pins as bp  # noqa: E402


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
                         ["ci.yml", "justfile", "pre-commit"])


class InlinePins(unittest.TestCase):
    def test_the_tracked_spellings_agree(self):
        self.assertEqual(bp.check_inline_pins(), [])


class TestFilePins(unittest.TestCase):
    """check_test_files.sh counts, mirrored between ci.yml and the justfile."""

    def _check(self, sources, files):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            scan = []
            for name, text in sources.items():
                path = root / name
                path.write_text(text, encoding="utf-8")
                scan.append(path)
            for name in files:
                path = root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("", encoding="utf-8")
            saved = (bp.ROOT, bp.TEST_FILE_PIN_SCAN)
            bp.ROOT = root
            bp.TEST_FILE_PIN_SCAN = tuple(scan)
            try:
                return bp.check_test_file_pins()
            finally:
                bp.ROOT, bp.TEST_FILE_PIN_SCAN = saved

    def test_copies_that_agree_with_the_suite_pass(self):
        errors = self._check(
            {"justfile": 'bash tools/check_test_files.sh 2 "suite/test*.py"',
             "ci.yml": "bash tools/check_test_files.sh 2 'suite/test*.py'"},
            ["suite/test_a.py", "suite/test_b.py"])
        self.assertEqual(errors, [])

    def test_a_copy_left_behind_is_reported_at_both_sites(self):
        # The finding this gate exists for: CI bumped 1 -> 2 with the new file,
        # the justfile mirror kept 1, and nothing executed the mirror.
        errors = self._check(
            {"justfile": 'bash tools/check_test_files.sh 1 "suite/test*.py"',
             "ci.yml": "bash tools/check_test_files.sh 2 'suite/test*.py'"},
            ["suite/test_a.py", "suite/test_b.py"])
        self.assertEqual(len(errors), 2)
        self.assertIn("justfile:1", errors[0])
        self.assertIn("ci.yml:1", errors[1])

    def test_copies_that_agree_with_each_other_but_not_the_suite_are_reported(self):
        errors = self._check(
            {"justfile": 'bash tools/check_test_files.sh 1 "suite/test*.py"',
             "ci.yml": "bash tools/check_test_files.sh 1 'suite/test*.py'"},
            ["suite/test_a.py", "suite/test_b.py"])
        self.assertEqual(len(errors), 2)
        self.assertIn("matches 2 file(s)", errors[0])

    def test_a_glob_matching_nothing_is_reported(self):
        errors = self._check(
            {"ci.yml": "bash tools/check_test_files.sh 1 'suite/test*.py'"}, [])
        self.assertEqual(len(errors), 1)
        self.assertIn("matches 0 file(s)", errors[0])

    def test_both_quotings_and_a_single_sided_pin_are_read(self):
        errors = self._check(
            {"justfile": 'bash tools/check_test_files.sh 1 "suite/test*.py"'},
            ["suite/test_a.py"])
        self.assertEqual(errors, [])

    def test_the_tracked_copies_agree(self):
        self.assertEqual(bp.check_test_file_pins(), [])

    def test_the_mirrored_globs_are_pinned_in_both_files(self):
        # A recipe that mirrors a CI step but drops its pin is the same hole
        # from the other side.
        pins = bp.collect_test_file_pins()
        self.assertTrue(pins)
        mirrored = [glob for glob, uses in pins.items() if len(uses) > 1]
        self.assertGreater(len(mirrored), 1)


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

    def test_pins_with_no_version_to_report_are_not_targets(self):
        for name in ("doxygen-awesome", "doxygen-sha256", "llvm-key-sha256",
                     "emsdk", "kicad"):
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
        status, message = self._check(
            "shellcheck", "ShellCheck - shell script analysis tool\n"
                          "version: 0.11.0\n")
        self.assertEqual(status, 0)
        self.assertIn("0.11.0.1", message)

    def test_a_missing_tool_reports_how_to_install_that_tool(self):
        for name, want in (("just", "pip install rust-just==1.52.0"),
                           ("shellcheck",
                            "pip install shellcheck-py==0.11.0.1"),
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
    """The files the daydream install mirrors, and their line-ending pins."""

    def test_every_cross_repo_rule_contributes_its_sources(self):
        installed = bp.installed_sources()
        for path in ("README.md", "hardware/pov_segment_map.json",
                     "scripts/shader_workbench.mjs", "scripts/sha256.mjs",
                     "scripts/engine_catalog.json",
                     "patterns/example.shader.json",
                     "patterns/shaderball_migration.json"):
            self.assertIn(path, installed)
        self.assertTrue(
            any(path.startswith("docs/screenshots/") for path in installed))

    def test_a_directory_rule_selects_only_its_patterns(self):
        # patterns/ also holds a README the FILES_MATCHING patterns exclude.
        self.assertNotIn("patterns/README.md", bp.installed_sources())

    def test_a_generated_artifact_is_not_a_repository_source(self):
        # The module and its glue are installed from the build directory.
        for path in bp.installed_sources():
            self.assertTrue((bp.ROOT / path).is_file(), path)

    def test_the_live_install_set_is_pinned(self):
        self.assertEqual(bp.check_install_eol(bp.installed_sources()), [])

    def test_an_unpinned_file_is_reported(self):
        self.assertEqual(len(bp.check_install_eol(["unpinned/file.txt"])), 1)

    def test_an_empty_set_fails_instead_of_passing_vacuously(self):
        self.assertTrue(bp.check_install_eol([]))


if __name__ == "__main__":
    unittest.main()
