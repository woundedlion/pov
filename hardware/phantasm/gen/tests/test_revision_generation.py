"""Revision selection, electrical partitions and legacy generator provenance."""
import contextlib
import hashlib
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
from test_board import dangling_pins  # noqa: E402
from test_check import committed_board_nets  # noqa: E402
from test_pcb_generation import GENERATES, GENERATES_REASON  # noqa: E402

PROTOTYPE = GEN.parent / "1.3"
LEGACY_HASHES = {
    "fp-lib-table": "fadeef1477fd022f6d1f2cef7a2f6a4e0133f193d770ca68fb6ed8eb2c46f7c9",
    "phantasm.kicad_pcb": "e94f57b9acb25da9535abd8dfae65cd8b94008185308cb501c974b1c4c1e4c28",
    "phantasm.kicad_pro": "5f7017ac920f088bd626b39c705ccf037a992e7208d65dc951187affccfa5534",
    "phantasm.kicad_sch": "c1ee258c1767aa53e900880f429d7f30eb9b7222b695aa783df32e7430fa601c",
    "phantasm.kicad_sym": "225312bff8e7626bba6c2ba271aac410cf397eacbca540304a6323146590135f",
    "phantasm.pretty/Teensy4.0.kicad_mod": "5f9f953b323356b7f498b73df4a8bd371c3759b7ebb4c2a65685175f9c3c5c95",
    "phantasm.pretty/Teensy4.0.wrl": "f120ec83682c6d9ab152f7671a7b713dc5044dbcb37e8b23f3393c5dc9b477a5",
    "phantasm.pretty/TerminalBlock_GCT_TBC05-02-1-G-G.kicad_mod": "de30f2ca5859eed9402b2cc86f64d889fefe2c808cdafae25f68cc173ed11c74",
    "phantasm.pretty/TerminalBlock_GCT_TBC05-03-1-G-G.kicad_mod": "526b838ef8340373cfc319b7f11f9fe393202d974492d1ce78636609b2550360",
    "sym-lib-table": "d3f2a4c416f8e75b72735784e2679123d55312be551e8b1d92884e98d7b83baf"
}


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
                check.main(["--revision", "1.3", "legacy.kicad_sch"])
        with contextlib.redirect_stderr(io.StringIO()):
            result = shorts.main(["--revision", "1.3", str(GEN.parent / "1.2" / "phantasm.kicad_sch")])
        self.assertEqual(result, 2)

    def test_committed_board_matches_differential_electrical_partition(self):
        self.assertTrue(check.check(committed_board_nets("1.3"), "1.3"))

    def test_schematic_has_no_shorts_or_unconnected_pins(self):
        self.assertEqual(shorts.analyze(self.schematic)[0], [])
        self.assertEqual(dangling_pins(self.schematic), [])

    def test_legacy_board_fails_differential_gate(self):
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
        positions, bounds = {}, {}
        for footprint in F(self.board, "footprint"):
            if sexp.val(footprint, "locked") != ["yes"]:
                continue
            ref = next(str(p[2]) for p in F(footprint, "property") if p[1] == "Reference")
            if ref.startswith("H"):
                continue
            positions[ref] = tuple(map(float, sexp.val(footprint, "at")))
            bounds[ref] = pcb.fp_bbox(footprint, graphic_layers=pcb.COURTYARD_LAYERS)
        self.assertEqual(pcb.keepout_clashes(positions, bounds, pcb.QUILTER_LENGTH), [])
        refs = list(positions)
        for index, ref in enumerate(refs):
            x, y, angle = positions[ref]
            x0, y0, x1, y1 = pcb._rot_bb(bounds[ref], angle)
            box = (x + x0, y + y0, x + x1, y + y1)
            for other in refs[index + 1:]:
                ox, oy, rotation = positions[other]
                a, b, c, d = pcb._rot_bb(bounds[other], rotation)
                self.assertFalse(pcb._boxes_overlap(box, (ox + a, oy + b, ox + c, oy + d)),
                                 f"{ref}/{other}")

    def test_schematic_revision_mismatch_is_rejected_before_placement(self):
        root = sexp.parse_one('(export (design (sheet (name "/") (title_block (rev "1.2")))))')
        with tempfile.TemporaryDirectory() as out, \
                mock.patch.object(pcb, "export_netlist", return_value=root), \
                mock.patch.object(pcb, "kicad_cli", return_value="unused"):
            with self.assertRaisesRegex(SystemExit, "schematic revision 1.2.*revision 1.3"):
                pcb.main(unplaced=True, revision="1.3", output_dir=out)
        self.assertEqual(pcb._GENERATION.get(), ("1.2", None))

    def test_prototype_cannot_be_exported_with_legacy_bom(self):
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
        cls.legacy = Path(cls.directory.name) / "1.2"
        cls.prototype = Path(cls.directory.name) / "1.3"
        generate(cls.legacy, "1.2")
        cls.before = {path: (cls.legacy / path).read_bytes() for path in LEGACY_HASHES}
        generate(cls.prototype, "1.3")
        generate(cls.legacy, "1.2")

    def test_legacy_output_is_byte_identical_after_revision_switch(self):
        for path, expected in LEGACY_HASHES.items():
            actual = (self.legacy / path).read_bytes()
            self.assertEqual(actual, self.before[path], path)
            self.assertEqual(hashlib.sha256(actual).hexdigest(), expected, path)

    def test_prototype_schematic_exports_expected_pin_partition(self):
        root = export_netlist(kicad_cli(), str(self.prototype / "phantasm.kicad_sch"))
        self.assertEqual(check.netlist_revision(root), "1.3")
        self.assertTrue(check.check(check.netlist_nets(root), "1.3"))

    def test_rev_13_generation_reproduces_committed_board(self):
        self.assertEqual((self.prototype / "phantasm.kicad_pcb").read_bytes(),
                         (PROTOTYPE / "phantasm.kicad_pcb").read_bytes())


if __name__ == "__main__":
    unittest.main()
