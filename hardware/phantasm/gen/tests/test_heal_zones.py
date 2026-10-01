import contextlib
import io
from pathlib import Path
import sys
import tempfile
import unittest

GEN = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN))

import fab  # noqa: E402
import heal_zones  # noqa: E402
import sexp  # noqa: E402
from kicad_common import F, is_copper_pour  # noqa: E402


ZONE = '''(zone (net "GND") (layer "In1.Cu")
    (connect_pads (clearance 0.2)) (min_thickness 0.0254)
    (fill yes (thermal_gap 0.0254) (thermal_bridge_width 0.2))
    (polygon (pts (xy 0 0) (xy 10 0) (xy 10 10) (xy 0 10)))
    (filled_polygon (layer "In1.Cu") (pts (xy 0 0) (xy 10 0) (xy 10 10))))'''
KEEPOUT = '''(zone (net "GND") (layer "F.Cu") (min_thickness 0.01)
    (keepout (tracks allowed) (copperpour not_allowed)))'''
ROUTING = '''(footprint "Test" (at 1 2 90) (property "Reference" "J1")
    (pad "1" thru_hole circle (at 0 0) (size 2 2) (drill 1) (net "GND")))
  (segment (start 0 0) (end 10 0) (width 0.3) (layer "F.Cu") (net "GND"))
  (via (at 5 0) (size 0.6) (drill 0.3) (layers "F.Cu" "B.Cu") (net "GND"))'''
BOARD = '(kicad_pcb\n  ' + ROUTING + '\n  ' + ZONE + '\n  ' + \
    ZONE.replace('In1.Cu', 'In2.Cu') + '\n  ' + KEEPOUT + '\n)\n'


class HealZoneTests(unittest.TestCase):
    def test_repairs_both_planes_and_preserves_other_geometry(self):
        repaired, count = heal_zones.heal_zones(BOARD)
        self.assertEqual(count, 2)
        self.assertIn(ROUTING, repaired)
        self.assertIn(KEEPOUT, repaired)
        before, after = sexp.parse_one(BOARD), sexp.parse_one(repaired)
        for key in ("segment", "via", "footprint"):
            self.assertEqual(F(before, key), F(after, key))
        for old, new in zip(F(before, "zone"), F(after, "zone")):
            self.assertEqual(F(old, "polygon"), F(new, "polygon"))
        for zone in F(after, "zone"):
            if not is_copper_pour(zone):
                continue
            self.assertFalse(F(zone, "filled_polygon"))
            self.assertEqual(sexp.val(zone, "min_thickness"), ["0.25"])
            self.assertEqual(sexp.val(F(zone, "fill")[0], "thermal_gap"), ["0.5"])
            self.assertEqual(sexp.val(F(zone, "fill")[0], "thermal_bridge_width"), ["0.5"])
        self.assertEqual(fab.validate_zone_geometry("fixture", board=after), 2)
        self.assertEqual(heal_zones.heal_zones(repaired), (repaired, 0))

    def test_preserves_larger_features_and_invalidates_all_fill_caches(self):
        safe = ZONE.replace("0.0254", "0.8").replace(
            "(thermal_bridge_width 0.2)", "(thermal_bridge_width 0.8)")
        source = '(kicad_pcb ' + ZONE + safe + ')'
        repaired, count = heal_zones.heal_zones(source)
        self.assertEqual(count, 1)
        self.assertNotIn("filled_polygon", repaired)
        self.assertIn("(thermal_gap 0.8)", repaired)

    def test_safe_board_remains_byte_identical_including_fills(self):
        source = BOARD.replace("0.0254", "0.5").replace(
            "(thermal_bridge_width 0.2)", "(thermal_bridge_width 0.5)")
        self.assertEqual(heal_zones.heal_zones(source), (source, 0))

    def test_preserves_crlf_line_endings(self):
        repaired, _ = heal_zones.heal_zones(BOARD.replace("\n", "\r\n"))
        self.assertNotIn("\n", repaired.replace("\r\n", ""))

    def test_rejects_ambiguous_or_invalid_fields(self):
        for replacement in ("", "(min_thickness 0.2) (min_thickness 0.3)",
                            "(min_thickness nan)", "(min_thickness inf)",
                            "(min_thickness 0)", "(min_thickness -1)",
                            "(min_thickness foo)", "(min_thickness 0.2 0.3)"):
            with self.subTest(replacement=replacement), self.assertRaises(ValueError):
                heal_zones.heal_zones(BOARD.replace("(min_thickness 0.0254)", replacement))

    def test_rejects_disabled_or_solid_pours_and_missing_plane(self):
        for source in (BOARD.replace("(fill yes", "(fill no"),
                       BOARD.replace("(connect_pads ", "(connect_pads yes "),
                       '(kicad_pcb ' + ZONE + ')', '(kicad_sch)',
                       BOARD + BOARD):
            with self.subTest(source=source), self.assertRaises(ValueError):
                heal_zones.heal_zones(source)


class HealZoneCliTests(unittest.TestCase):
    def setUp(self):
        directory = self.enterContext(tempfile.TemporaryDirectory())
        self.source = Path(directory) / "source.kicad_pcb"
        self.output = Path(directory) / "healed.kicad_pcb"
        self.source.write_text(BOARD, encoding="utf-8")
        self.stdout = self.enterContext(contextlib.redirect_stdout(io.StringIO()))
        self.enterContext(contextlib.redirect_stderr(io.StringIO()))

    def run_cli(self, *args):
        return heal_zones.main([str(self.source), *map(str, args)])

    def test_explicit_copy_requires_refill_and_preserves_source(self):
        original = self.source.read_bytes()
        self.assertEqual(self.run_cli("--output", self.output), 0)
        self.assertEqual(self.source.read_bytes(), original)
        self.assertIn("Refill all zones", self.stdout.getvalue())
        self.assertEqual(heal_zones.main([str(self.output), "--check"]), 0)
        self.assertEqual(self.run_cli("--check"), 1)

    def test_overwriting_an_output_requires_force(self):
        self.output.write_text("previous output", encoding="utf-8")
        self.assertEqual(self.run_cli("-o", self.output), 1)
        self.assertEqual(self.output.read_text(encoding="utf-8"), "previous output")
        self.assertEqual(self.run_cli("-o", self.output, "--force"), 0)

    def test_never_overwrites_source_even_with_force(self):
        original = self.source.read_bytes()
        self.assertEqual(self.run_cli("-o", self.source, "--force"), 1)
        self.assertEqual(self.source.read_bytes(), original)

    def test_never_writes_into_manifested_directory(self):
        (self.output.parent / "SHA256SUMS.txt").touch()
        self.assertEqual(self.run_cli("-o", self.output, "--force"), 1)
        self.assertFalse(self.output.exists())

    def test_bad_source_does_not_replace_existing_output(self):
        self.output.write_text("previous output", encoding="utf-8")
        self.source.write_text("(kicad_pcb)", encoding="utf-8")
        self.assertEqual(self.run_cli("-o", self.output, "--force"), 1)
        self.assertEqual(self.output.read_text(encoding="utf-8"), "previous output")


if __name__ == "__main__":
    unittest.main()
