#!/usr/bin/env python3
"""Host tests for the test-pin updater (tools/update_test_pins.py).

The updater rewrites the ratchet the whole native suite is gated by, so the
rules that keep a partial measurement from writing a bad floor are what these
drive: the minimum has to come from at least two applying sources, an entry no
source observed has to survive untouched, and a red or whole-suite-skipped run
is not a measurement at all. The parsing side is driven too, because the
roster's shape (a backslash-continued X-macro whose effects rows select
between tier constants) is what tells the updater which number to move. The
drift allowance is driven from both sides, since a floor no module can fall
under is as useless as one it clears many times over.

Run:  python -m unittest discover -s tools/update_test_pins_tests
"""

import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TOOLS))

import update_test_pins as up  # noqa: E402


ROSTER = '''constexpr int EFFECTS_QUICK_MIN_ASSERTIONS = 100;
constexpr int EFFECTS_FULL_MIN_ASSERTIONS = 200;

#define HS_TEST_MODULE_LIST(X)                                                 \\
  X("alpha", 7, hs_test::alpha_tests::run_alpha_tests, 30014)                  \\
  X("beta", 3, hs_test::beta_tests::run_beta_tests,                            \\
    77)                                                                        \\
  X("effects", 5, hs_test::effects_tests::run_effects_tests,                   \\
    hs_test::effects_tests::effects_full_suite()                               \\
        ? EFFECTS_FULL_MIN_ASSERTIONS                                          \\
        : EFFECTS_QUICK_MIN_ASSERTIONS)

#define HS_TEST_MODULE_ENTRY(name, case_sites, fn, min_assertions)             \\
'''


def source(label, tier, counts, skips=()):
    """A Source standing in for one configuration's run."""
    return up.Source(label, tier, label, dict(counts), list(skips))


def floor(name, value, module, tier=None):
    """A Floor over a standalone pin, detached from any file."""
    return up.Floor(up.Pin(Path("run_tests.cpp"), name, 0, 0, value), module, tier)


class RosterParsing(unittest.TestCase):
    """The roster's shape decides which literal each pin names."""

    def setUp(self):
        self.original = up.scan_headers
        up.scan_headers = lambda: {
            "test_alpha.h": ("run_alpha_tests", 7),
            "test_beta.h": ("run_beta_tests", 4),
            "test_effects.h": ("run_effects_tests", 5),
        }
        self.bound = up.bind_off_roster
        up.bind_off_roster = lambda headers, covered, cases: []

    def tearDown(self):
        up.scan_headers = self.original
        up.bind_off_roster = self.bound

    def test_case_pins_carry_the_count_scanned_from_their_header(self):
        cases, _, errors = up.parse_roster(ROSTER)
        self.assertEqual(errors, [])
        self.assertEqual(
            {pin.value: sites for pin, sites in cases}, {7: 7, 3: 4, 5: 5})

    def test_wrapped_row_floor_is_found_across_the_continuation(self):
        _, floors, _ = up.parse_roster(ROSTER)
        beta = [f for f in floors if f.module == "beta"]
        self.assertEqual([f.pin.value for f in beta], [77])
        self.assertEqual(ROSTER[beta[0].pin.start:beta[0].pin.end], "77")

    def test_tiered_row_yields_one_floor_per_tier_constant(self):
        _, floors, _ = up.parse_roster(ROSTER)
        tiered = {f.tier: f.pin.value for f in floors if f.module == "effects"}
        self.assertEqual(tiered, {"full": 200, "quick": 100})

    def test_quick_tier_floor_names_the_constant_the_ternary_selects_last(self):
        _, floors, _ = up.parse_roster(ROSTER)
        quick = next(f for f in floors if f.tier == "quick")
        self.assertEqual(quick.pin.name, "EFFECTS_QUICK_MIN_ASSERTIONS")


class MinimumOverSources(unittest.TestCase):
    """A floor has to clear every configuration, so it is the minimum."""

    def test_minimum_wins_over_the_larger_configuration(self):
        sources = [source("a", "quick", {"m": 900}), source("b", "quick", {"m": 700})]
        self.assertEqual(up.measure_floor(floor("m", 500, "m"), sources)[0], 700)

    def test_one_source_cannot_move_a_floor(self):
        value, note = up.measure_floor(
            floor("m", 500, "m"),
            [source("a", "quick", {"m": 900}), source("b", "quick", {})])
        self.assertIsNone(value)
        self.assertIn("observed only by a", note)

    def test_tier_floor_ignores_the_other_tier(self):
        sources = [
            source("a", "quick", {"e": 100}), source("b", "quick", {"e": 120}),
            source("c", "full", {"e": 9}),
        ]
        self.assertEqual(
            up.measure_floor(floor("e", 50, "e", "quick"), sources)[0], 100)

    def test_tier_floor_with_one_matching_source_is_held(self):
        sources = [source("a", "quick", {"e": 100}), source("b", "full", {"e": 900})]
        value, note = up.measure_floor(floor("e", 50, "e", "full"), sources)
        self.assertIsNone(value)
        self.assertIn("1 full-tier source(s)", note)

    def test_unfloored_module_stays_unfloored(self):
        sources = [source("a", "quick", {"m": 900}), source("b", "quick", {"m": 700})]
        value, note = up.measure_floor(floor("m", 0, "m"), sources)
        self.assertIsNone(value)
        self.assertIn("unfloored", note)


class DriftAllowance(unittest.TestCase):
    """A floor several times under its module has stopped bounding it."""

    def test_a_floor_the_module_clears_many_times_over_has_drifted(self):
        self.assertTrue(up.drifted(189868, 1186754))

    def test_headroom_inside_the_fraction_is_not_drift(self):
        self.assertFalse(up.drifted(1000000, 1200000))

    def test_a_small_floor_is_held_to_the_absolute_bound_not_the_fraction(self):
        self.assertFalse(up.drifted(12, 24))
        self.assertTrue(up.drifted(12, 300))

    def test_a_wide_fraction_under_the_absolute_bound_is_not_drift(self):
        self.assertFalse(up.drifted(100, 250))

    def test_a_floor_the_measurement_matches_has_no_drift(self):
        self.assertFalse(up.drifted(5000, 5000))


class WriteMargin(unittest.TestCase):
    """--margin writes the floor under the measurement, never over it."""

    def test_no_margin_writes_the_measurement(self):
        self.assertEqual(up.with_margin(1186754, 0.0), 1186754)

    def test_margin_is_a_percentage_of_the_measurement(self):
        self.assertEqual(up.with_margin(1000000, 1.0), 990000)

    def test_margin_cannot_round_a_small_measurement_down(self):
        self.assertEqual(up.with_margin(12, 1.0), 12)

    def test_a_margined_floor_does_not_read_as_drifted(self):
        measured = 1186754
        self.assertFalse(up.drifted(up.with_margin(measured, 1.0), measured))

    def test_a_margin_outside_a_percentage_is_refused(self):
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(up.main(["--margin", "100"]), 1)
            self.assertEqual(up.main(["--margin", "-1"]), 1)


class UnobservedEntriesSurvive(unittest.TestCase):
    """A narrowed run must never drop what it did not run."""

    def test_module_no_source_ran_is_held_and_named(self):
        sources = [source("a", "quick", {"m": 900}), source("b", "quick", {"m": 700})]
        value, note = up.measure_floor(floor("other", 12345, "other"), sources)
        self.assertIsNone(value)
        self.assertEqual(note, "unobserved by a, b")

    def test_a_partially_observed_floor_is_held_not_raised(self):
        sources = [
            source("a", "quick", {"m": 900}), source("b", "quick", {}),
            source("c", "quick", {}),
        ]
        self.assertIsNone(up.measure_floor(floor("m", 500, "m"), sources)[0])


class RunOutputParsing(unittest.TestCase):
    """Only a green run that actually asserted is a measurement."""

    def test_footer_totals_passed_plus_failed(self):
        parsed = up.parse_run_output(
            "a", "quick", "log", "=== alpha: 40 passed, 0 failed ===\n")
        self.assertEqual(parsed.counts, {"alpha": 40})

    def test_skipped_variant_of_the_footer_still_counts(self):
        parsed = up.parse_run_output(
            "a", "quick", "log", "=== alpha: 40 passed, 0 failed, 2 SKIPPED ===\n")
        self.assertEqual(parsed.counts, {"alpha": 40})

    def test_whole_suite_skip_is_recorded_not_measured(self):
        parsed = up.parse_run_output(
            "a", "quick", "log",
            "=== alpha: 4 passed, 0 failed ===\n"
            "=== death: 0 passed, 0 failed, 1 SKIPPED ===\n")
        self.assertEqual(parsed.counts, {"alpha": 4})
        self.assertEqual(parsed.suite_skips, ["death"])

    def test_a_red_run_is_rejected(self):
        with self.assertRaises(SystemExit):
            up.parse_run_output(
                "a", "quick", "log", "=== alpha: 40 passed, 1 failed ===\n")

    def test_output_with_no_footer_is_rejected(self):
        with self.assertRaises(SystemExit):
            up.parse_run_output("a", "quick", "log", "=== alpha ===\n")

    def test_repeated_module_takes_the_lower_count(self):
        parsed = up.parse_run_output(
            "a", "quick", "log",
            "=== alpha: 40 passed, 0 failed ===\n"
            "=== alpha: 31 passed, 0 failed ===\n")
        self.assertEqual(parsed.counts, {"alpha": 31})


class CaseScanning(unittest.TestCase):
    """The case count mirrors tests/check_case_calls.cmake's scan."""

    def test_a_case_named_only_in_a_comment_is_not_a_definition(self):
        text = up.strip_code(
            "// void test_ghost()\n/* void check_ghost() */\nvoid test_real() {}\n")
        self.assertEqual(up.CASE_DEF.findall(text), ["test_real"])

    def test_a_case_named_in_a_string_is_not_a_definition(self):
        text = up.strip_code('const char *s = "void test_ghost(";\n')
        self.assertEqual(up.CASE_DEF.findall(text), [])

    def test_a_wrapped_signature_still_counts(self):
        self.assertEqual(
            up.CASE_DEF.findall("\ninline static void\ntest_wrapped(int a) {}\n"),
            ["test_wrapped"])


class DriftReporting(unittest.TestCase):
    """--check-drift only fails when it is asked for, and only on drift."""

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.log = Path(self.directory.name) / "drift.log"
        self.log.write_text(
            "=== scan: 1186754 passed, 0 failed ===\n"
            "=== mesh: 69923 passed, 0 failed ===\n", encoding="utf-8")
        self.roster = up.parse_roster
        up.parse_roster = lambda text: (
            [],
            [floor("scan floor", 189868, "scan"), floor("mesh floor", 69923, "mesh")],
            [])

    def tearDown(self):
        up.parse_roster = self.roster
        self.directory.cleanup()

    def run_main(self, *flags):
        """main() over one drifted floor, with its report captured."""
        argv = ["--log", f"quick={self.log}", "--log", f"full={self.log}", *flags]
        printed = io.StringIO()
        with contextlib.redirect_stdout(printed):
            status = up.main(argv)
        return status, printed.getvalue()

    def test_drift_is_ignored_unless_the_flag_is_passed(self):
        status, printed = self.run_main()
        self.assertEqual(status, 0)
        self.assertIn("0 drifted", printed)

    def test_drift_fails_the_run_and_names_the_floor(self):
        status, printed = self.run_main("--check-drift")
        self.assertEqual(status, 1)
        self.assertIn("scan floor: pinned 189868, measured 1186754", printed)

    def test_remeasure_offers_the_measurement_instead_of_failing_on_drift(self):
        status, printed = self.run_main("--check-drift", "--remeasure")
        self.assertEqual(status, 1, "the pin still has to be written")
        self.assertIn("189868 -> 1186754 RAISED", printed)
        self.assertNotIn("have drifted so far", printed)

    def test_margin_writes_under_the_measurement(self):
        _, printed = self.run_main("--check-drift", "--remeasure", "--margin", "1")
        self.assertIn("189868 -> 1174887 RAISED", printed)

    def test_repairing_drift_leaves_a_floor_on_its_measurement_alone(self):
        _, printed = self.run_main("--check-drift", "--remeasure", "--margin", "1")
        self.assertNotIn("mesh floor", printed)

    def test_remeasure_without_the_drift_flag_margins_every_floor_down(self):
        _, printed = self.run_main("--remeasure", "--margin", "1")
        self.assertIn("mesh floor: 69923 -> 69224 LOWERED", printed)


class RealTreeAgreement(unittest.TestCase):
    """The committed pins are what this tool reads back out of the tree."""

    def test_every_case_pin_reproduces(self):
        cases, floors, errors = up.parse_roster(
            up.RUN_TESTS.read_text(encoding="utf-8"))
        self.assertEqual(errors, [])
        self.assertTrue(floors)
        drift = [(pin.name, pin.value, sites)
                 for pin, sites in cases if pin.value != sites]
        self.assertEqual(drift, [])


if __name__ == "__main__":
    unittest.main()
