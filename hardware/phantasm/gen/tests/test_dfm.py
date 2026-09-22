"""Mask export regressions and placement constraints for regenerated projects."""
import contextlib
import io
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

GEN = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN))

import board  # noqa: E402
import fab  # noqa: E402
import pcb  # noqa: E402
import sexp  # noqa: E402
from kicad_common import F, kicad_cli  # noqa: E402
from test_pcb_generation import GENERATES, GENERATES_REASON, generate, read  # noqa: E402


MASK_BOARD = """(kicad_pcb
    (setup
        (pad_to_mask_clearance 0)
        (solder_mask_min_width 0.1)
        (allow_soldermask_bridges_in_footprints no)
        (tenting (front yes) (back yes)))
    (via (at 1 2) (size 0.45) (drill 0.2) (layers "F.Cu" "B.Cu")))"""


class SolderMaskTests(unittest.TestCase):
    def check(self, text):
        fab.validate_solder_mask("fixture.kicad_pcb", board=sexp.parse_one(text))

    def test_corrected_board_and_inherited_via_tenting_pass(self):
        self.check(MASK_BOARD)
        fab.validate_solder_mask(fab.PCB)

    def test_missing_or_small_mask_web_fails(self):
        for replacement in ("", "(solder_mask_min_width 0.09)",
                            "(solder_mask_min_width nan)"):
            with self.subTest(replacement=replacement):
                with self.assertRaisesRegex(fab.SolderMaskError, "web"):
                    self.check(MASK_BOARD.replace("(solder_mask_min_width 0.1)", replacement))

    def test_missing_or_enabled_bridge_exception_fails(self):
        for replacement in ("", "(allow_soldermask_bridges_in_footprints yes)"):
            with self.subTest(replacement=replacement):
                with self.assertRaisesRegex(fab.SolderMaskError, "bridges"):
                    self.check(MASK_BOARD.replace(
                        "(allow_soldermask_bridges_in_footprints no)", replacement))

    def test_missing_or_open_via_tenting_fails(self):
        for replacement in ("", "(tenting (front no) (back yes))",
                            "(tenting (front yes) (back no))"):
            with self.subTest(replacement=replacement):
                with self.assertRaisesRegex(fab.SolderMaskError, "tenting"):
                    self.check(MASK_BOARD.replace(
                        "(tenting (front yes) (back yes))", replacement))

    def test_a_via_cannot_override_board_tenting(self):
        root = sexp.parse_one(MASK_BOARD)
        F(root, "via")[0].append(sexp.parse_one("(tenting (back no))"))
        with self.assertRaisesRegex(fab.SolderMaskError, "via: back"):
            fab.validate_solder_mask("fixture.kicad_pcb", board=root)

    def test_failed_mask_gate_stops_before_exports(self):
        root = sexp.parse_one(MASK_BOARD.replace("(solder_mask_min_width 0.1)", ""))
        with mock.patch.object(fab, "kicad_cli", return_value="fixture-cli"), \
                mock.patch.object(fab, "read_board", return_value=root), \
                mock.patch.object(fab, "run_export") as export, \
                contextlib.redirect_stdout(io.StringIO()):
            with self.assertRaisesRegex(SystemExit, "solder mask web"):
                fab.main()
        export.assert_not_called()


class NewProjectMarginsTests(unittest.TestCase):
    def test_seed_sets_clearances_and_rejects_silk_on_pads(self):
        settings = json.loads(board.project_seed("root"))["board"]["design_settings"]
        self.assertEqual(settings["rules"]["min_hole_clearance"], 0.1016)
        self.assertEqual(settings["rules"]["min_silk_clearance"], 0.15)
        self.assertEqual(settings["rules"]["solder_mask_to_copper_clearance"], 0.1)
        self.assertEqual(settings["rule_severities"]["silk_over_copper"], "error")

    def test_regeneration_updates_floors_and_preserves_other_settings(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "phantasm.kicad_pro"
            project = json.loads(board.project_seed("root"))
            settings = project["board"]["design_settings"]
            settings["rules"]["min_silk_clearance"] = 0
            settings["rules"]["min_track_width"] = 0.4
            settings["rule_severities"]["silk_over_copper"] = "ignore"
            project["text_variables"] = {"CUSTOM": "preserved"}
            path.write_text(json.dumps(project), encoding="utf-8")
            board.write_project(path, "new-root")
            updated = json.loads(path.read_text(encoding="utf-8"))
            rules = updated["board"]["design_settings"]["rules"]
            self.assertEqual(rules["min_silk_clearance"], 0.15)
            self.assertEqual(rules["min_track_width"], 0.4)
            self.assertEqual(updated["text_variables"], {"CUSTOM": "preserved"})
            self.assertEqual(updated["sheets"], [["new-root", "Root"]])

    def test_changed_footprints_do_not_inherit_fixed_placements(self):
        comps = {ref: (ref, libid, "", False)
                 for ref, libid in pcb.QUILTER_FIXED_FOOTPRINTS.items()}
        self.assertEqual(set(pcb.fixed_placements(comps)), set(comps))
        for ref in comps:
            changed = dict(comps, **{ref: (ref, "Different:Footprint", "", False)})
            self.assertNotIn(ref, pcb.fixed_placements(changed))


@unittest.skipUnless(GENERATES, GENERATES_REASON)
class GeneratedMaskTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.directory = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.directory.cleanup)
        cls.paths = {}
        for unplaced in (False, True):
            target = Path(cls.directory.name) / str(unplaced)
            target.mkdir()
            cls.paths[unplaced] = Path(generate(str(target), unplaced=unplaced))

    def test_placed_and_unplaced_boards_preserve_mask_policy(self):
        for unplaced, path in self.paths.items():
            with self.subTest(unplaced=unplaced):
                root = read(path)
                fab.validate_solder_mask("generated.kicad_pcb", board=root)
                setup = F(root, "setup")[0]
                self.assertEqual(float(sexp.val(setup, "pad_to_mask_clearance")[0]), 0)
                self.assertEqual(float(sexp.val(setup, "solder_mask_min_width")[0]), 0.1)

    def test_both_projects_carry_new_layout_rules(self):
        for unplaced, path in self.paths.items():
            with self.subTest(unplaced=unplaced):
                project = json.loads(path.with_suffix(".kicad_pro").read_text(encoding="utf-8"))
                self.assertEqual(project["meta"]["filename"], path.with_suffix(".kicad_pro").name)
                self.assertTrue(project["sheets"][0][0])
                rules = project["board"]["design_settings"]["rules"]
                self.assertEqual(rules["min_silk_clearance"], 0.15)
                self.assertEqual(rules["solder_mask_to_copper_clearance"], 0.1)
                self.assertEqual(rules["min_hole_clearance"], 0.25 if unplaced else 0.1016)

    def test_generated_placements_have_no_physical_drc_violations(self):
        for unplaced, path in self.paths.items():
            with self.subTest(unplaced=unplaced):
                report = path.with_suffix(".drc.json")
                result = subprocess.run(
                    [kicad_cli(), "pcb", "drc", "--format", "json",
                     "--severity-error", "--severity-warning", "-o", str(report), str(path)],
                    capture_output=True, text=True, timeout=120)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                findings = json.loads(report.read_text(encoding="utf-8"))
                self.assertIn("violations", findings)
                physical = [item for item in findings["violations"]
                            if item["type"] not in {"lib_footprint_issues", "lib_footprint_mismatch"}]
                self.assertEqual(physical, [])


if __name__ == "__main__":
    unittest.main()
