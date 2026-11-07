"""Revision selection, electrical partitions and current revision generator provenance."""
import contextlib
import copy
import math
import io
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

GEN = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN))

import analyze_candidates  # noqa: E402
import board  # noqa: E402
import check  # noqa: E402
import fab  # noqa: E402
import pcb  # noqa: E402
import sexp  # noqa: E402
import shorts  # noqa: E402
from kicad_common import F, export_netlist, kicad_cli, rotate_point  # noqa: E402
from courtyard_bounds import courtyard_box  # noqa: E402
from pin_landing import dangling_pins  # noqa: E402
from kicad_env import committed_board_nets, generate, footprints_by_ref  # noqa: E402
from kicad_env import GENERATES, GENERATES_REASON  # noqa: E402

PROTOTYPE = GEN.parent / "1.3"
GENERATED_PROJECT_FILES = (
    "fp-lib-table",
    "phantasm.kicad_pcb",
    "phantasm.kicad_pro",
    "phantasm.kicad_sch",
    "phantasm.kicad_sym",
    "phantasm.pretty/Teensy4.0.kicad_mod",
    "phantasm.pretty/Teensy4.0.wrl",
    "phantasm.pretty/TerminalBlock_GCT_TBC05-02-1-G-G.kicad_mod",
    "phantasm.pretty/TerminalBlock_GCT_TBC05-03-1-G-G.kicad_mod",
    "sym-lib-table",
)




class PrototypeContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.schematic = sexp.parse_one((PROTOTYPE / "phantasm.kicad_sch").read_text(encoding="utf-8"))
        cls.board = sexp.parse_one((PROTOTYPE / "phantasm.kicad_pcb").read_text(encoding="utf-8"))

    def test_default_generation_target_is_1_2(self):
        self.assertEqual(board.parse_args([]).revision, "1.2")
        self.assertEqual(pcb.parse_args([]).revision, "1.2")

    def test_positional_schematic_revision_mismatch_is_rejected(self):
        root = sexp.parse_one('(export (design (sheet (name "/") (title_block (rev "1.2")))))')
        with mock.patch.object(check, "export_netlist", return_value=root), \
                mock.patch.object(check, "kicad_cli", return_value="unused"):
            with self.assertRaisesRegex(SystemExit, "does not match requested 1.3"):
                check.main(["--revision", "1.3", "rev_12.kicad_sch"])
        with contextlib.redirect_stderr(io.StringIO()):
            result = shorts.main(["--revision", "1.3", str(GEN.parent / "1.2" / "phantasm.kicad_sch")])
        self.assertEqual(result, 2)

    def test_committed_board_matches_differential_electrical_partition(self):
        self.assertTrue(check.check(committed_board_nets("1.3"), "1.3"))

    def test_led_resistors_follow_their_output_channels(self):
        footprints = footprints_by_ref(self.board)
        data_x = float(sexp.val(footprints["R_D1"], "at")[0])
        clock_x = float(sexp.val(footprints["R_D2"], "at")[0])
        self.assertGreater(data_x, clock_x)
        nets = committed_board_nets("1.3")
        for ref, pin, net in (("R_D1", "6", "DATA_SRC"), ("R_D2", "3", "CLK_SRC")):
            self.assertIn(f"U1.{pin}", nets[net])
            source_pad = next(p for p in F(footprints[ref], "pad") if str(p[1]) == "2")
            self.assertEqual(sexp.val(source_pad, "net")[-1], f"/{net}")
            self.assertLess(math.dist(self.pad_position(footprints[ref], "2"),
                                      self.pad_position(footprints["U1"], pin)), 3)

    @staticmethod
    def pad_position(footprint, number):
        pad = next(p for p in F(footprint, "pad") if str(p[1]) == number)
        x, y, angle = map(float, sexp.val(footprint, "at"))
        dx, dy = rotate_point(tuple(map(float, sexp.val(pad, "at")[:2])), angle)
        return x + dx, y + dy

    def test_front_layout_preserves_full_terminal_approaches(self):
        footprints = footprints_by_ref(self.board)
        self.assertEqual(len(footprints), 36)
        self.assertTrue(all(sexp.val(f, "layer") == ["F.Cu"] for f in footprints.values()))
        for ref, count in (("J1", 2), ("J2", 3), ("J3A", 4), ("J3B", 4)):
            footprint = footprints[ref]
            self.assertEqual(float(sexp.val(footprint, "at")[2]), 270)
            zones = {str(sexp.val(z, "name", [""])[0]): z for z in F(footprint, "zone")}
            for pin in range(1, count + 1):
                zone = zones[f"Terminal wire access pin {pin}"]
                points = F(F(F(zone, "polygon")[0], "pts")[0], "xy")
                xs, ys = zip(*[(float(p[1]), float(p[2])) for p in points])
                self.assertAlmostEqual(max(xs) - min(xs), 2.4)
                self.assertAlmostEqual(max(ys) - min(ys), 10)
                px, py = self.pad_position(footprint, str(pin))
                self.assertAlmostEqual((min(xs) + max(xs)) / 2, px)
                self.assertAlmostEqual(max(ys), py - 3.25)
        pcb.validate_front_access(list(footprints.values()), 58.28)

    def test_ground_and_shield_approaches_reject_off_center_obstruction(self):
        for pin in ("3", "4"):
            footprints = footprints_by_ref(copy.deepcopy(self.board))
            x, _ = self.pad_position(footprints["J3B"], pin)
            probe = sexp.parse_one(
                f'(footprint "probe" (layer "F.Cu") (at {x + 1.15} -5 0) '
                '(property "Reference" "R_PROBE") '
                '(fp_rect (start -0.025 -0.025) (end 0.025 0.025) '
                '(stroke (width 0.01) (type solid)) (fill none) (layer "F.CrtYd")))')
            with self.assertRaisesRegex(ValueError, f"J3B Terminal wire access pin {pin} obstructed"):
                pcb.validate_front_access([*footprints.values(), probe], 58.28)

    def test_usb_mating_volume_rejects_front_obstruction(self):
        footprints = footprints_by_ref(copy.deepcopy(self.board))
        probe = sexp.parse_one(
            '(footprint "probe" (layer "F.Cu") (at 0 20 0) '
            '(property "Reference" "R_PROBE") '
            '(fp_rect (start -0.1 -0.1) (end 0.1 0.1) '
            '(stroke (width 0.01) (type solid)) (fill none) (layer "F.CrtYd")))')
        with self.assertRaisesRegex(ValueError, "U_MCU USB mating plug obstructed by R_PROBE"):
            pcb.validate_front_access([*footprints.values(), probe], 58.28)

    def test_bypass_pads_remain_within_three_mm_of_supply_pins(self):
        footprints = footprints_by_ref(self.board)
        for cap, ic, pin in (("C_DEC1", "U_MCU", "VIN"), ("C_DEC2", "U1", "14"),
                             ("C_DEC3", "U_SYNC", "8")):
            self.assertLess(math.dist(self.pad_position(footprints[cap], "1"),
                                      self.pad_position(footprints[ic], pin)), 3)

    def test_termination_legend_is_revision_specific(self):
        self.assertIn("TERM=CABLE ENDS", [text[1] for text in F(self.board, "gr_text")])
        legacy = sexp.parse_one((GEN.parent / "1.2" / "phantasm.kicad_pcb").read_text(encoding="utf-8"))
        self.assertNotIn("TERM=CABLE ENDS", [text[1] for text in F(legacy, "gr_text")])

    def test_sync_parts_have_orderable_manufacturer_fields(self):
        expected = {
            "U_SYNC": "THVD2410DR",
            "D_SYNC": "CDSOT23-SM712",
            "R_A": "CRCW060310R0FKEAHP",
            "R_B": "CRCW060310R0FKEAHP",
            "R_TERM": "CRCW0805120RFKEAHP",
            "C_DEC3": "CC0603KRX7R9BB104",
            "C_BULK3": "CC0603KRX7R8BB105",
            "C_IN": "UPW1H101MPD",
        }
        symbols = {props["Reference"]: props for symbol in F(self.schematic, "symbol")
                   if (props := {p[1]: p[2] for p in F(symbol, "property")})}
        footprints = {ref: {p[1]: p[2] for p in F(footprint, "property")}
                      for ref, footprint in footprints_by_ref(self.board).items()}
        for ref, mpn in expected.items():
            self.assertEqual(symbols[ref]["MPN"], mpn, ref)
            self.assertTrue(symbols[ref]["Manufacturer"], ref)
            self.assertTrue(symbols[ref]["Datasheet"].startswith("https://"), ref)
            for field in ("MPN", "Manufacturer", "Datasheet"):
                self.assertEqual(footprints[ref][field], symbols[ref][field], (ref, field))

    def test_led_reset_pulldown_open_is_rejected(self):
        nets = committed_board_nets("1.3")
        for ref, net in (("R_DATA_PD", "DATA_IN"), ("R_CLK_PD", "CLK_IN")):
            nets[net].remove(ref)
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertFalse(check.check(nets, "1.3"))
            nets[net].add(ref)

    def test_schematic_has_no_shorts_or_unconnected_pins(self):
        self.assertEqual(shorts.analyze(self.schematic)[0], [])
        self.assertEqual(dangling_pins(self.schematic), [])

    def test_rev_12_board_fails_differential_gate(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            self.assertFalse(check.check(committed_board_nets("1.2"), "1.3"))
        self.assertIn("FAIL MASTER_EN", output.getvalue())
        self.assertIn("FAIL SYNC_A", output.getvalue())

    def test_connector_pad_swap_is_rejected(self):
        nets = committed_board_nets("1.3")
        nets["SYNC_A"].remove("J3B.1")
        nets["SYNC_A"].add("J3B.2")
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertFalse(check.check(nets, "1.3"))

    def test_locked_courtyards_clear_each_other_and_mounting_reservations(self):
        positions, bounds, courtyards = {}, {}, {}
        for ref, footprint in footprints_by_ref(self.board).items():
            if sexp.val(footprint, "locked") != ["yes"]:
                continue
            if ref.startswith("H"):
                continue
            positions[ref] = tuple(map(float, sexp.val(footprint, "at")))
            bounds[ref] = pcb.fp_bbox(footprint, graphic_layers=pcb.COURTYARD_LAYERS)
            courtyards[ref] = courtyard_box(footprint)
            self.assertIsNotNone(courtyards[ref], ref)
        self.assertEqual(pcb.keepout_clashes(positions, bounds, pcb.QUILTER_LENGTH), [])
        refs = list(positions)
        for index, ref in enumerate(refs):
            for other in refs[index + 1:]:
                self.assertFalse(pcb._boxes_overlap(courtyards[ref], courtyards[other]),
                                 f"{ref}/{other}")

    def test_schematic_revision_mismatch_is_rejected_before_placement(self):
        root = sexp.parse_one('(export (design (sheet (name "/") (title_block (rev "1.2")))))')
        with tempfile.TemporaryDirectory() as out, \
                mock.patch.object(pcb, "export_netlist", return_value=root), \
                mock.patch.object(pcb, "kicad_cli", return_value="unused"):
            with self.assertRaisesRegex(SystemExit, "schematic revision 1.2.*revision 1.3"):
                pcb.main(unplaced=True, revision="1.3", output_dir=out)
        self.assertEqual(pcb._GENERATION.get(), ("1.2", None))

    def test_prototype_cannot_be_exported_with_rev_11_bom(self):
        result = subprocess.run([sys.executable, str(GEN / "fab.py"), "--revision", "1.3"],
                                capture_output=True, text=True, check=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("revision-specific BOM are not validated", result.stderr)

    def test_prototype_cannot_receive_single_ended_si_score(self):
        with self.assertRaisesRegex(ValueError, "no validated differential-bus scoring model"):
            analyze_candidates.analyze(str(PROTOTYPE / "phantasm.kicad_pcb"))


@unittest.skipUnless(GENERATES, GENERATES_REASON)
class RevisionGenerationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.directory = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.directory.cleanup)
        cls.rev_12 = Path(cls.directory.name) / "1.2"
        cls.prototype = Path(cls.directory.name) / "1.3"
        generate(cls.rev_12, unplaced=True, revision="1.2")
        cls.before = {path: (cls.rev_12 / path).read_bytes() for path in GENERATED_PROJECT_FILES}
        generate(cls.prototype, unplaced=True, revision="1.3")
        generate(cls.rev_12, unplaced=True, revision="1.2")

    def test_current_schematics_pass_real_erc(self):
        for revision in ("1.2", "1.3"):
            with self.subTest(revision=revision), tempfile.TemporaryDirectory() as directory:
                project = Path(directory) / revision
                shutil.copytree(GEN.parent / revision, project,
                                ignore=shutil.ignore_patterns("routed", "candidates"))
                report = str(Path(directory) / "erc.json")
                with mock.patch.object(fab, "SCH", str(project / "phantasm.kicad_sch")):
                    self.assertEqual(fab.run_erc(report), 0)

    def test_inner_planes_pass_between_teensy_header_pads(self):
        for project in (self.rev_12, self.prototype):
            root = sexp.parse_one((project / "phantasm.kicad_pcb").read_text(encoding="utf-8"))
            teensy = footprints_by_ref(root)["U_MCU"]
            rows = {}
            for pad in F(teensy, "pad"):
                x, y = map(float, sexp.val(pad, "at")[:2])
                rows.setdefault(y, []).append((x, float(sexp.val(pad, "size")[0])))
            for zone in F(root, "zone"):
                if sexp.val(zone, "layer") not in (["In1.Cu"], ["In2.Cu"]):
                    continue
                clearance = float(sexp.val(F(zone, "connect_pads")[0], "clearance")[0])
                thickness = float(sexp.val(zone, "min_thickness")[0])
                for pads in rows.values():
                    pads.sort()
                    for (x0, width0), (x1, width1) in zip(pads, pads[1:]):
                        with self.subTest(revision=project.name, layer=sexp.val(zone, "layer")):
                            self.assertGreaterEqual(x1 - x0 - (width0 + width1) / 2
                                                    - 2 * clearance, thickness)

    def test_rev_12_output_is_byte_identical_after_revision_switch(self):
        for path in GENERATED_PROJECT_FILES:
            actual = (self.rev_12 / path).read_bytes()
            self.assertEqual(actual, self.before[path], path)
            self.assertEqual(actual, (GEN.parent / "1.2" / path).read_bytes(), path)

    def test_prototype_schematic_exports_expected_pin_partition(self):
        root = export_netlist(kicad_cli(), str(self.prototype / "phantasm.kicad_sch"))
        self.assertEqual(check.netlist_revision(root), "1.3")
        self.assertTrue(check.check(check.netlist_nets(root), "1.3"))

    def test_rev_13_generation_reproduces_committed_project(self):
        for path in GENERATED_PROJECT_FILES:
            self.assertEqual((self.prototype / path).read_bytes(),
                             (PROTOTYPE / path).read_bytes(), path)


if __name__ == "__main__":
    unittest.main()
