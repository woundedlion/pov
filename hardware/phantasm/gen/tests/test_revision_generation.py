"""Revision selection, electrical partitions and current revision generator provenance."""
import contextlib
import io
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
import pcb  # noqa: E402
import sexp  # noqa: E402
import shorts  # noqa: E402
from kicad_common import F, export_netlist, kicad_cli  # noqa: E402
from courtyard_bounds import courtyard_box  # noqa: E402
from pin_landing import dangling_pins  # noqa: E402
from test_check import committed_board_nets  # noqa: E402
from test_pcb_generation import GENERATES, GENERATES_REASON  # noqa: E402

PROTOTYPE = GEN.parent / "1.3"
REV_12_FILES = (
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


def generate(out, revision):
    with contextlib.redirect_stdout(io.StringIO()):
        board.main(force=True, revision=revision, output_dir=str(out))
        pcb.main(unplaced=True, force=True, force_teensy_library=True,
                 revision=revision, output_dir=str(out))


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
        footprints = {props["Reference"]: props for footprint in F(self.board, "footprint")
                      if (props := {p[1]: p[2] for p in F(footprint, "property")})}
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
        for footprint in F(self.board, "footprint"):
            if sexp.val(footprint, "locked") != ["yes"]:
                continue
            ref = next(str(p[2]) for p in F(footprint, "property") if p[1] == "Reference")
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
        generate(cls.rev_12, "1.2")
        cls.before = {path: (cls.rev_12 / path).read_bytes() for path in REV_12_FILES}
        generate(cls.prototype, "1.3")
        generate(cls.rev_12, "1.2")

    def test_rev_12_output_is_byte_identical_after_revision_switch(self):
        for path in REV_12_FILES:
            actual = (self.rev_12 / path).read_bytes()
            self.assertEqual(actual, self.before[path], path)
            self.assertEqual(actual, (GEN.parent / "1.2" / path).read_bytes(), path)

    def test_prototype_schematic_exports_expected_pin_partition(self):
        root = export_netlist(kicad_cli(), str(self.prototype / "phantasm.kicad_sch"))
        self.assertEqual(check.netlist_revision(root), "1.3")
        self.assertTrue(check.check(check.netlist_nets(root), "1.3"))

    def test_rev_13_generation_reproduces_committed_project(self):
        for path in REV_12_FILES:
            self.assertEqual((self.prototype / path).read_bytes(),
                             (PROTOTYPE / path).read_bytes(), path)


if __name__ == "__main__":
    unittest.main()
