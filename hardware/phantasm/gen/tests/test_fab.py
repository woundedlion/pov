import sys
import tempfile
import unittest
from pathlib import Path

GEN = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN))

import fab  # noqa: E402


class ViaGeometryTests(unittest.TestCase):
    def validate_source(self, source):
        with tempfile.TemporaryDirectory() as directory:
            pcb = Path(directory) / "test.kicad_pcb"
            pcb.write_text(source, encoding="utf-8")
            return fab.validate_via_geometry(pcb)

    def validate(self, size, drill):
        return self.validate_source(
            "(kicad_pcb "
            f"(via (at 1 2) (size {size}) (drill {drill}) "
            '(layers "F.Cu" "B.Cu")))'
        )

    def test_accepts_standard_cost_via(self):
        self.assertEqual(self.validate("0.45", "0.20"), 1)

    def test_rejects_small_via_diameter(self):
        with self.assertRaisesRegex(
                fab.ViaGeometryError, "0.25 mm diameter is below 0.45 mm"):
            self.validate("0.25", "0.20")

    def test_rejects_small_via_drill(self):
        with self.assertRaisesRegex(
                fab.ViaGeometryError, "0.15 mm drill is below 0.2 mm"):
            self.validate("0.45", "0.15")

    def test_rejects_close_via_copper(self):
        source = (
            "(kicad_pcb "
            '(via (at 1 2) (size 0.45) (drill 0.20) '
            '(layers "F.Cu" "B.Cu")) '
            '(via (at 1.48 2) (size 0.45) (drill 0.20) '
            '(layers "F.Cu" "B.Cu")))'
        )
        with self.assertRaisesRegex(
                fab.ViaGeometryError,
                "0.03 mm copper spacing is below 0.15 mm"):
            self.validate_source(source)


class ZoneGeometryTests(unittest.TestCase):
    def validate(self, min_thickness, thermal_gap, bridge_width="0.1016"):
        with tempfile.TemporaryDirectory() as directory:
            pcb = Path(directory) / "test.kicad_pcb"
            pcb.write_text(
                '(kicad_pcb (zone (net "GND") (name "GND_IN1") '
                f'(min_thickness {min_thickness}) '
                f'(fill yes (thermal_gap {thermal_gap}) '
                f'(thermal_bridge_width {bridge_width}))))',
                encoding="utf-8")
            return fab.validate_zone_geometry(pcb)

    def test_accepts_process_floor_fill(self):
        self.assertEqual(self.validate("0.1016", "0.1016"), 1)

    def test_rejects_thin_pour_sliver(self):
        with self.assertRaisesRegex(
                fab.ZoneGeometryError,
                "GND_IN1: min_thickness 0.0254 mm is below 0.1016 mm"):
            self.validate("0.0254", "0.1016")

    def test_rejects_unresolvable_thermal_gap(self):
        with self.assertRaisesRegex(
                fab.ZoneGeometryError,
                "GND_IN1: thermal_gap 0.0254 mm is below 0.1016 mm"):
            self.validate("0.1016", "0.0254")

    def test_rejects_thin_thermal_spoke(self):
        with self.assertRaisesRegex(
                fab.ZoneGeometryError,
                "GND_IN1: thermal_bridge_width 0.05 mm is below 0.1016 mm"):
            self.validate("0.1016", "0.1016", "0.05")

    def test_ignores_keepout_zones(self):
        with tempfile.TemporaryDirectory() as directory:
            pcb = Path(directory) / "test.kicad_pcb"
            pcb.write_text(
                '(kicad_pcb (zone (name "H1 screw head") '
                '(min_thickness 0.0254) (keepout (tracks not_allowed))))',
                encoding="utf-8")
            self.assertEqual(fab.validate_zone_geometry(pcb), 0)


class AssemblyMetadataTests(unittest.TestCase):
    def setUp(self):
        self.ref = "X1"
        self.comps = {self.ref: {"lcsc": " C123 "}}
        self.posrows = {
            self.ref: {
                "PosX": " 1.25 ",
                "PosY": "-2",
                "Rot": "370",
                "Side": " TOP ",
            }
        }

    def validate(self):
        return fab.validate_assembly_metadata(
            self.comps, self.posrows, [self.ref])

    def assert_diagnostic(self, expected):
        with self.assertRaises(fab.AssemblyMetadataError) as caught:
            self.validate()
        self.assertIn(expected, caught.exception.diagnostics)

    def test_accepts_complete_metadata_and_normalizes_text(self):
        metadata = self.validate()

        self.assertEqual(
            metadata[self.ref],
            {
                "lcsc": "C123",
                "pos_x": "1.25",
                "pos_y": "-2",
                "rotation": 370.0,
                "side": "top",
            },
        )

    def test_rejects_missing_supplier_part_number(self):
        del self.comps[self.ref]["lcsc"]

        self.assert_diagnostic(
            f"{self.ref}: supplier part number (LCSC) is missing")

    def test_rejects_blank_supplier_part_number(self):
        self.comps[self.ref]["lcsc"] = "  "

        self.assert_diagnostic(
            f"{self.ref}: supplier part number (LCSC) is blank")

    def test_rejects_missing_centroid_row(self):
        self.posrows = {}

        self.assert_diagnostic(f"{self.ref}: centroid row is missing")

    def test_rejects_missing_centroid_fields(self):
        for field, label in (("PosX", "X"), ("PosY", "Y"),
                             ("Rot", "rotation"), ("Side", "side")):
            with self.subTest(field=field):
                original = self.posrows[self.ref].pop(field)
                self.assert_diagnostic(
                    f"{self.ref}: centroid {label} is missing")
                self.posrows[self.ref][field] = original

    def test_rejects_blank_centroid_fields(self):
        for field, label in (("PosX", "X"), ("PosY", "Y"),
                             ("Rot", "rotation"), ("Side", "side")):
            with self.subTest(field=field):
                original = self.posrows[self.ref][field]
                self.posrows[self.ref][field] = "  "
                self.assert_diagnostic(
                    f"{self.ref}: centroid {label} is blank")
                self.posrows[self.ref][field] = original

    def test_rejects_nonnumeric_centroid_coordinates_and_rotation(self):
        for field, label in (("PosX", "X"), ("PosY", "Y"),
                             ("Rot", "rotation")):
            with self.subTest(field=field):
                original = self.posrows[self.ref][field]
                self.posrows[self.ref][field] = "invalid"
                self.assert_diagnostic(
                    f"{self.ref}: centroid {label} is not numeric: 'invalid'")
                self.posrows[self.ref][field] = original

    def test_rejects_nonfinite_centroid_coordinates_and_rotation(self):
        for field, label, value in (("PosX", "X", "nan"),
                                    ("PosY", "Y", "inf"),
                                    ("Rot", "rotation", "-inf")):
            with self.subTest(field=field):
                original = self.posrows[self.ref][field]
                self.posrows[self.ref][field] = value
                self.assert_diagnostic(
                    f"{self.ref}: centroid {label} is not finite: {value!r}")
                self.posrows[self.ref][field] = original

    def test_rejects_wrong_centroid_side(self):
        self.posrows[self.ref]["Side"] = "bottom"

        self.assert_diagnostic(
            f"{self.ref}: centroid side is 'bottom'; expected 'top'")

    def test_aggregates_diagnostics(self):
        self.comps[self.ref]["lcsc"] = ""
        self.posrows = {}

        with self.assertRaises(fab.AssemblyMetadataError) as caught:
            self.validate()
        self.assertEqual(
            caught.exception.diagnostics,
            (
                f"{self.ref}: supplier part number (LCSC) is blank",
                f"{self.ref}: centroid row is missing",
            ),
        )


class PartCatalogTests(unittest.TestCase):
    def test_accepts_every_assigned_part(self):
        metadata = {
            ref: {"lcsc": lcsc}
            for ref, lcsc in fab.LCSC_BY_REF.items()
        }

        fab.validate_part_catalog(metadata)

    def test_rejects_unmapped_supplier_part(self):
        with self.assertRaisesRegex(
                fab.PartCatalogError,
                "X1: LCSC C999999999 has no catalog entry"):
            fab.validate_part_catalog({"X1": {"lcsc": "C999999999"}})


if __name__ == "__main__":
    unittest.main()
