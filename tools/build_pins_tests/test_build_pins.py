#!/usr/bin/env python3
"""Host tests for the build-pin gate (tools/build_pins.py).

The gate is the single source for every externally-installed version, and each
of its three checks reads a foreign file shape: duplicates_pin scans workflow
YAML two lines at a time, check_engine_ranges parses package.json's `>=X`
string, and _version_tuple compares versions of unequal width. A drift in any
of those shapes makes the check detect nothing while still printing PASS.

Run:  python -m unittest discover -s tools/build_pins_tests
"""

import json
import sys
import tempfile
import unittest
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


if __name__ == "__main__":
    unittest.main()
