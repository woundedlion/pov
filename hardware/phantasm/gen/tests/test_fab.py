import json
import re
import subprocess
import sys
import tempfile
import unittest
import unittest.mock
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


class PlotOriginTests(unittest.TestCase):
    def validate_source(self, source):
        with tempfile.TemporaryDirectory() as directory:
            pcb = Path(directory) / "test.kicad_pcb"
            pcb.write_text(source, encoding="utf-8")
            return fab.validate_plot_origin(pcb)

    def test_accepts_board_without_plot_origin(self):
        source = "(kicad_pcb (setup (pad_to_mask_clearance 0)))"

        self.assertEqual(self.validate_source(source), (0.0, 0.0))

    def test_accepts_plot_origin_at_absolute_zero(self):
        self.assertEqual(
            self.validate_source("(kicad_pcb (setup (aux_axis_origin 0 0)))"),
            (0.0, 0.0))

    def test_rejects_offset_plot_origin(self):
        with self.assertRaisesRegex(
                fab.PlotOriginError,
                r"drill/place origin is 10,-5 mm"):
            self.validate_source("(kicad_pcb (setup (aux_axis_origin 10 -5)))")

    def test_rejects_invalid_plot_origin(self):
        with self.assertRaisesRegex(
                fab.PlotOriginError, "aux_axis_origin is invalid: 10"):
            self.validate_source("(kicad_pcb (setup (aux_axis_origin 10)))")

    def test_rejects_unreadable_board(self):
        with self.assertRaisesRegex(
                fab.PlotOriginError, "cannot read PCB setup"):
            fab.validate_plot_origin("does-not-exist.kicad_pcb")

    def test_committed_board_exports_in_absolute_coordinates(self):
        self.assertEqual(fab.validate_plot_origin(fab.PCB), (0.0, 0.0))


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


class AssemblyPolicyTests(unittest.TestCase):
    """JLC reflows top-side SMD only; every hand-soldered part stays out."""

    def component(self, footprint, value="", dnp=False):
        return {"footprint": footprint, "value": value, "dnp": dnp}

    def test_excludes_hand_soldered_packages(self):
        for footprint in (
                "Connector_PinHeader_2.54mm:PinHeader_1x03_P2.54mm_Vertical",
                "Connector_JST:JST_XA_B02B-XASK-1-A_1x02_P2.50mm_Vertical",
                "Jumper:SolderJumper-2_P1.3mm_Open_TrianglePad1.0x1.5mm",
                "Capacitor_THT:CP_Radial_D8.0mm_P3.50mm"):
            with self.subTest(footprint=footprint):
                self.assertFalse(fab.is_assembled(self.component(footprint)))

    def test_excludes_teensy_module(self):
        self.assertFalse(fab.is_assembled(self.component("", "Teensy4.0")))

    def test_includes_top_side_smd(self):
        for footprint in ("Package_TO_SOT_SMD:SOT-23",
                          "Resistor_SMD:R_0603_1608Metric",
                          "Package_SO:SOIC-14_3.9x8.7mm_P1.27mm"):
            with self.subTest(footprint=footprint):
                self.assertTrue(fab.is_assembled(self.component(footprint)))

    def test_excludes_dnp(self):
        self.assertFalse(fab.is_assembled(
            self.component("Resistor_SMD:R_0603_1608Metric", dnp=True)))

    def test_accepts_exact_assigned_part_set(self):
        fab.validate_assembled_refs(fab.LCSC_BY_REF)

    def test_rejects_assigned_part_missing_from_assembly(self):
        assembled = set(fab.LCSC_BY_REF) - {"U1"}

        with self.assertRaisesRegex(
                fab.AssemblyMetadataError,
                "assigned parts missing from assembly: U1"):
            fab.validate_assembled_refs(assembled)

    def test_rejects_rotation_correction_missing_from_assembly(self):
        assembled = set(fab.LCSC_BY_REF) - {"U1"}

        with self.assertRaisesRegex(
                fab.AssemblyMetadataError,
                "rotation corrections missing from assembly: U1"):
            fab.validate_rotation_refs(assembled)

    def test_excludes_the_j1_footprint_the_schematic_generator_emits(self):
        source = (GEN / "board.py").read_text(encoding="utf-8")
        match = re.search(r'"J1".*?fp="([^"]+)"', source, re.DOTALL)
        self.assertIsNotNone(match, "board.py assigns J1 no footprint")
        self.assertFalse(fab.is_assembled(self.component(match.group(1))))


PARITY_DESCRIPTIONS = {
    "extra_footprint": "Extra footprint",
    "footprint_symbol_mismatch": "Footprint attributes don't match symbol",
}


class SchematicParityTests(unittest.TestCase):
    KNOWN = [
        {"type": kind, "description": PARITY_DESCRIPTIONS[kind],
         "items": [{"description": f"Footprint {ref}"}]}
        for kind, ref in sorted(fab.KNOWN_PARITY_ITEMS)
    ]

    def require_report(self, report):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "parity.json"
            path.write_text(json.dumps(report), encoding="utf-8")
            return fab.require_schematic_parity(path)

    def require(self, entries):
        return self.require_report({"schematic_parity": entries})

    def test_accepts_known_differences(self):
        self.assertEqual(self.require(self.KNOWN), len(self.KNOWN))

    def test_rejects_report_without_a_parity_section(self):
        with self.assertRaisesRegex(
                fab.SchematicParityError,
                "no schematic_parity section"):
            self.require_report({"coordinate_units": "mm"})

    def test_rejects_report_missing_known_differences(self):
        with self.assertRaisesRegex(
                fab.SchematicParityError,
                r"lists 0 items, fewer than the \d+ KiCad reports"):
            self.require([])

    def test_rejects_net_drift(self):
        entries = self.KNOWN + [{
            "type": "net_conflict",
            "description": "Pad net (/DATA_SRC) doesn't match net given by "
                           "schematic (/DATA)",
            "items": [{"description":
                       "Pad 1 [/DATA_SRC] of R_D1 on Top Layer"}],
        }]

        with self.assertRaisesRegex(
                fab.SchematicParityError,
                r"net_conflict: Pad net .* doesn't match"):
            self.require(entries)

    def test_rejects_unknown_extra_footprint(self):
        entries = self.KNOWN + [{"type": "extra_footprint",
                                 "description": "Extra footprint",
                                 "items": [{"description": "Footprint R_PDX"}]}]

        with self.assertRaisesRegex(
                fab.SchematicParityError, "re-route it in Quilter"):
            self.require(entries)

    def test_rejects_missing_footprint(self):
        entries = self.KNOWN + [{"type": "missing_footprint",
                                 "description": "Missing footprint R_PD (10k)"}]

        with self.assertRaisesRegex(
                fab.SchematicParityError,
                r"missing_footprint: Missing footprint R_PD"):
            self.require(entries)

    def test_rejects_unreadable_report(self):
        with self.assertRaisesRegex(
                fab.SchematicParityError, "cannot read parity report"):
            fab.require_schematic_parity("does-not-exist.json")


class DesignRuleTests(unittest.TestCase):
    def require_report(self, report):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "drc.json"
            path.write_text(json.dumps(report), encoding="utf-8")
            return fab.require_clean_drc(path)

    def test_accepts_clean_report(self):
        self.assertEqual(
            self.require_report({"violations": [], "unconnected_items": []}),
            (0, 0))

    def test_rejects_violations(self):
        with self.assertRaisesRegex(
                fab.DesignRuleError,
                "1 error-severity violations, 0 unconnected items"):
            self.require_report({
                "violations": [{"type": "clearance",
                                "description": "Clearance violation"}],
                "unconnected_items": [],
            })

    def test_rejects_unconnected_items(self):
        with self.assertRaisesRegex(
                fab.DesignRuleError,
                "0 error-severity violations, 2 unconnected items"):
            self.require_report({
                "violations": [],
                "unconnected_items": [{"type": "unconnected_items"},
                                      {"type": "unconnected_items"}],
            })

    def test_rejects_report_without_violation_sections(self):
        with self.assertRaisesRegex(
                fab.DesignRuleError,
                "no violations/unconnected_items sections"):
            self.require_report({"schematic_parity": []})

    def test_rejects_unreadable_report(self):
        with self.assertRaisesRegex(
                fab.DesignRuleError, "cannot read DRC report"):
            fab.require_clean_drc("does-not-exist.json")

    def test_rejects_nonzero_exit_on_a_clean_report(self):
        def fake_run(args, check=True, **kw):
            Path(args[args.index("-o") + 1]).write_text(
                json.dumps({"violations": [], "unconnected_items": []}),
                encoding="utf-8")
            return subprocess.CompletedProcess(args, 5, "", "")

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "drc.json"
            with unittest.mock.patch.object(fab, "run", fake_run):
                with self.assertRaisesRegex(
                        fab.DesignRuleError, "kicad-cli drc exited 5"):
                    fab.run_drc(path)


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
