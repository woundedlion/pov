import contextlib
import io
import os
import sys
import tempfile
import unittest
import unittest.mock
from pathlib import Path

from tools import check_domain_ratchets


class ParseFloors(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def write(self, name: str, text: str) -> Path:
        path = Path(self.tmp.name) / name
        path.write_text(text, encoding="utf-8")
        return path

    def test_rejects_harness_without_floor(self):
        path = self.write("harness.cpp", "")
        with self.assertRaisesRegex(SystemExit, "MIN_RELAX_BAKES_VERIFIED"):
            check_domain_ratchets.floors(path)

    def test_ignores_a_commented_floor(self):
        path = self.write(
            "harness.cpp",
            "constexpr int MIN_RELAX_BAKES_VERIFIED = 21;\n"
            "// constexpr int MIN_RELAX_BAKES_VERIFIED = 5;\n"
            "/* constexpr int MIN_RELAX_BAKES_VERIFIED = 5; */\n",
        )
        self.assertEqual(
            check_domain_ratchets.floors(path),
            {"MIN_RELAX_BAKES_VERIFIED": 21},
        )

    def test_rejects_a_floor_declared_twice(self):
        path = self.write(
            "harness.cpp",
            "constexpr int MIN_RELAX_BAKES_VERIFIED = 21;\n"
            "constexpr int MIN_RELAX_BAKES_VERIFIED = 5;\n",
        )
        with self.assertRaisesRegex(SystemExit, "declared twice"):
            check_domain_ratchets.floors(path)


class ParseGaps(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def write(self, text: str) -> Path:
        path = Path(self.tmp.name) / "death.h"
        path.write_text(text, encoding="utf-8")
        return path

    def test_ignores_a_commented_row(self):
        path = self.write(GAP_TABLE.format(
            rows='    {"sdf.h", 3},\n    // {"sdf.h", 99},\n'))
        self.assertEqual(
            check_domain_ratchets.death_pins(path), {"guard_gap.sdf.h": 3})

    def test_rejects_a_row_declared_twice(self):
        path = self.write(GAP_TABLE.format(
            rows='    {"sdf.h", 3},\n    {"sdf.h", 99},\n'))
        with self.assertRaisesRegex(SystemExit, "declared twice"):
            check_domain_ratchets.death_pins(path)


HARNESS = "constexpr int MIN_RELAX_BAKES_VERIFIED = {floor};\n"
GAP_TABLE = "constexpr GuardGap GUARD_GAP_ALLOW[] = {{\n{rows}}};\n"


def death_source(gaps):
    """A death harness, carrying no GUARD_GAP_ALLOW table at all when gaps is None."""
    if gaps is None:
        return ""
    rows = "".join('    {{"{0}", {1}}},\n'.format(name, gap)
                   for name, gap in gaps.items())
    return GAP_TABLE.format(rows=rows)


class Ratchets(unittest.TestCase):
    """main()'s comparison, over temp harness and death files.

    A state is (bakes verified, guard-gap rows).
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def run_check(self, previous, current, allow=""):
        paths = []
        for tag, (bakes, gaps) in (("prev", previous), ("cur", current)):
            harness = Path(self.tmp.name) / (tag + "_harness.cpp")
            harness.write_text(HARNESS.format(floor=bakes), encoding="utf-8")
            deaths = Path(self.tmp.name) / (tag + "_death.h")
            deaths.write_text(death_source(gaps), encoding="utf-8")
            paths.append((harness, deaths))
        argv = [
            "check_domain_ratchets.py",
            str(paths[0][0]), str(paths[1][0]),
            str(paths[0][1]), str(paths[1][1]),
            "--previous-ref", "origin/master",
        ]
        out = io.StringIO()
        with unittest.mock.patch.object(sys, "argv", argv), \
                unittest.mock.patch.dict(
                    os.environ, {"DOMAIN_RATCHET_ALLOW_WEAKEN": allow}), \
                contextlib.redirect_stdout(out):
            status = check_domain_ratchets.main()
        return status, out.getvalue()

    def test_compares_floors_and_guard_gaps(self):
        base = (4, {"sdf.h": 2, "mesh.h": 1})
        cases = [
            ("unchanged", base, 0, "0 approved weakening(s)"),
            ("raised", (5, {"sdf.h": 1}), 0, "domain ratchets: 1 floors"),
            ("fewer bakes verified", (3, {"sdf.h": 2, "mesh.h": 1}), 1,
             "MIN_RELAX_BAKES_VERIFIED weakened (4 -> 3)"),
            ("widened gap", (4, {"sdf.h": 5, "mesh.h": 1}), 1,
             "guard_gap.sdf.h weakened (2 -> 5)"),
            ("new gap row", (4, {"sdf.h": 2, "mesh.h": 1, "new.h": 1}), 1,
             "guard_gap.new.h weakened (0 -> 1)"),
        ]
        for name, current, expected, message in cases:
            with self.subTest(name):
                status, output = self.run_check(base, current)
                self.assertEqual(status, expected, output)
                self.assertIn(message, output)

    def test_listed_weakening_is_approved(self):
        status, output = self.run_check(
            (4, {"sdf.h": 2}), (3, {"sdf.h": 5}),
            allow="MIN_RELAX_BAKES_VERIFIED\nguard_gap.sdf.h\n",
        )
        self.assertEqual(status, 0, output)
        self.assertNotIn("::error", output)
        self.assertIn("2 approved weakening(s)", output)

    def test_reports_an_allowance_that_approves_nothing(self):
        status, output = self.run_check(
            (4, {"sdf.h": 2}), (4, {"sdf.h": 2}),
            allow="guard_gap.sdf.h\n",
        )
        self.assertEqual(status, 0, output)
        self.assertIn('"guard_gap.sdf.h" approves no weakening', output)

    def test_skips_the_gap_check_when_the_previous_side_has_no_table(self):
        status, output = self.run_check((4, None), (4, {"sdf.h": 9}))
        self.assertEqual(status, 0, output)
        self.assertIn("gap-widening check skipped", output)

    def test_rejects_a_death_harness_without_gap_rows(self):
        with self.assertRaisesRegex(SystemExit, "no GUARD_GAP_ALLOW rows parsed"):
            self.run_check((4, {"sdf.h": 2}), (4, {}))


if __name__ == "__main__":
    unittest.main()
