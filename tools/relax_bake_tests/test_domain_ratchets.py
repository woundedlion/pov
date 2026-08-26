import contextlib
import io
import os
import subprocess
import sys
import tempfile
import unittest
import unittest.mock
from pathlib import Path

TOOLS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TOOLS))

import check_domain_ratchets  # noqa: E402


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
            allow=(
                "MIN_RELAX_BAKES_VERIFIED=4->3\n"
                "guard_gap.sdf.h=2->5\n"
            ),
        )
        self.assertEqual(status, 0, output)
        self.assertNotIn("::error", output)
        self.assertIn("2 approved weakening(s)", output)

    def test_reports_an_allowance_that_approves_nothing(self):
        status, output = self.run_check(
            (4, {"sdf.h": 2}), (4, {"sdf.h": 2}),
            allow="guard_gap.sdf.h=2->5\n",
        )
        self.assertEqual(status, 0, output)
        self.assertIn("transition was not exercised", output)

    def test_an_allowance_only_approves_its_exact_transition(self):
        status, output = self.run_check(
            (4, {"sdf.h": 2}), (4, {"sdf.h": 6}),
            allow="guard_gap.sdf.h=2->5\n",
        )
        self.assertEqual(status, 1, output)
        self.assertIn("guard_gap.sdf.h weakened (2 -> 6)", output)

    def test_rejects_a_bare_allowance_name(self):
        with self.assertRaisesRegex(SystemExit, "expected name=before->after"):
            self.run_check(
                (4, {"sdf.h": 2}), (4, {"sdf.h": 5}),
                allow="guard_gap.sdf.h\n",
            )

    def test_skips_the_gap_check_when_the_previous_side_has_no_table(self):
        status, output = self.run_check((4, None), (4, {"sdf.h": 9}))
        self.assertEqual(status, 0, output)
        self.assertIn("gap-widening check skipped", output)

    def test_rejects_a_death_harness_without_gap_rows(self):
        with self.assertRaisesRegex(SystemExit, "no GUARD_GAP_ALLOW rows parsed"):
            self.run_check((4, {"sdf.h": 2}), (4, {}))


class GitRange(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.repo = Path(self.tmp.name)
        self.git("init", "--initial-branch=master")
        self.git("config", "user.name", "Domain Ratchet Test")
        self.git("config", "user.email", "ratchet@example.invalid")
        (self.repo / "tools").mkdir()
        (self.repo / "tests").mkdir()
        (self.repo / ".github" / "workflows").mkdir(parents=True)

    def git(self, *args: str) -> str:
        return subprocess.run(
            ["git", "-C", str(self.repo), *args], check=True,
            capture_output=True, text=True
        ).stdout.strip()

    def commit(self, gap: int, workflow: str, message: str) -> str:
        (self.repo / "tools" / "relax_bake_harness.cpp").write_text(
            HARNESS.format(floor=4), encoding="utf-8"
        )
        (self.repo / "tests" / "test_death.h").write_text(
            death_source({"sdf.h": gap}), encoding="utf-8"
        )
        (self.repo / ".github" / "workflows" / "ci.yml").write_text(
            workflow, encoding="utf-8"
        )
        self.git("add", "tools/relax_bake_harness.cpp", "tests/test_death.h",
                 ".github/workflows/ci.yml")
        self.git("commit", "-m", message)
        return self.git("rev-parse", "HEAD")

    def test_checks_each_gated_commit_edge_in_a_batch(self):
        base = self.commit(2, "jobs:\n  lint:\n", "base")
        self.commit(4, "jobs:\n  lint:\n", "ungated weakening")
        head = self.commit(
            5, "jobs:\n  domain-ratchets:\n", "gated weakening"
        )
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            status = check_domain_ratchets.check_git_range(
                self.repo, base, head, {("guard_gap.sdf.h", "4", "5")}
            )
        self.assertEqual(status, 0, output.getvalue())
        self.assertIn("domain ratchet range: 1 commit edge(s) checked",
                      output.getvalue())


if __name__ == "__main__":
    unittest.main()
