import contextlib
import io
import json
import re
from decimal import Decimal
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path

GEN = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN))

import fab  # noqa: E402
import heal_clearance  # noqa: E402
import sexp  # noqa: E402
from constraints import NEW_LAYOUT_RULES  # noqa: E402
from kicad_common import F, property_value  # noqa: E402
from kicad_env import expected_nodes, footprints_by_ref  # noqa: E402


def checked_rule_count(project):
    document = json.loads(Path(project).read_text(encoding="utf-8"))
    return len(fab.rule_checks(document, fab.RULE_MINIMUMS, fab.DEFAULT_CLASS_MINIMUMS))


class ProjectRulesTests(unittest.TestCase):
    def write_project(self, directory, rules=None, default=None):
        document = {
            "board": {"design_settings": {
                "rules": {**fab.RULE_MINIMUMS, **NEW_LAYOUT_RULES},
                "rule_severities": {"silk_over_copper": "error"}}},
            "net_settings": {
                "classes": [dict(fab.DEFAULT_CLASS_MINIMUMS, name="Default")]
            },
        }
        document["board"]["design_settings"]["rules"].update(rules or {})
        document["net_settings"]["classes"][0].update(default or {})
        project = Path(directory) / "phantasm.kicad_pro"
        project.write_text(json.dumps(document), encoding="utf-8")
        return project

    def test_accepts_project_at_the_floor(self):
        with tempfile.TemporaryDirectory() as directory:
            project = self.write_project(directory)
            self.assertEqual(fab.validate_project_rules(project),
                             checked_rule_count(project))

    def test_rejects_gui_rezeroed_clearance(self):
        with tempfile.TemporaryDirectory() as directory:
            project = self.write_project(directory, rules={"min_clearance": 0})
            with self.assertRaises(fab.ProjectRulesError) as caught:
                fab.validate_project_rules(project)
        self.assertIn("min_clearance", str(caught.exception))
        self.assertIn("heal_clearance.py", str(caught.exception))
        self.assertIn("on a working copy of the project, then retry", str(caught.exception))

    def test_rejects_relaxed_default_net_class(self):
        with tempfile.TemporaryDirectory() as directory:
            project = self.write_project(directory, default={"via_drill": 0.1})
            with self.assertRaises(fab.ProjectRulesError) as caught:
                fab.validate_project_rules(project)
        self.assertIn("Default.via_drill", str(caught.exception))

    def test_rejects_project_without_default_net_class(self):
        with tempfile.TemporaryDirectory() as directory:
            project = Path(directory) / "phantasm.kicad_pro"
            project.write_text(
                json.dumps({"net_settings": {"classes": []}}), encoding="utf-8")
            with self.assertRaises(fab.ProjectRulesError) as caught:
                fab.validate_project_rules(project)
        self.assertIn("Default net class", str(caught.exception))

    def test_rejects_unreadable_project(self):
        with tempfile.TemporaryDirectory() as directory:
            project = Path(directory) / "phantasm.kicad_pro"
            project.write_text("{", encoding="utf-8")
            with self.assertRaises(fab.ProjectRulesError):
                fab.validate_project_rules(project)

    def test_committed_project_meets_the_fabrication_floors(self):
        self.assertEqual(fab.validate_project_rules(fab.PRO), checked_rule_count(fab.PRO))

    def test_rejects_disabled_layout_constraints(self):
        for field in NEW_LAYOUT_RULES:
            with self.subTest(field=field), tempfile.TemporaryDirectory() as directory:
                project = self.write_project(directory, rules={field: 0})
                with self.assertRaisesRegex(fab.ProjectRulesError, field):
                    fab.validate_project_rules(project)

    def test_heal_restores_silk_severity(self):
        with tempfile.TemporaryDirectory() as directory:
            project = self.write_project(directory)
            document = json.loads(project.read_text(encoding="utf-8"))
            document["board"]["design_settings"]["rule_severities"][
                "silk_over_copper"] = "warning"
            project.write_text(json.dumps(document), encoding="utf-8")
            with self.assertRaisesRegex(fab.ProjectRulesError, "silk_over_copper"):
                fab.validate_project_rules(project)
            self.assertTrue(heal_clearance.heal_project(project))
            fab.validate_project_rules(project)
            self.assertFalse(heal_clearance.heal_project(project))


class ViaGeometryTests(unittest.TestCase):
    def validate_source(self, source, min_vias=1):
        with tempfile.TemporaryDirectory() as directory:
            pcb = Path(directory) / "test.kicad_pcb"
            pcb.write_text(source, encoding="utf-8")
            return fab.validate_via_geometry(pcb, min_vias=min_vias)

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
                "0.03 mm via-to-via copper spacing is below 0.15 mm"):
            self.validate_source(source)

    def test_accepts_close_same_net_stitching_vias(self):
        source = (
            "(kicad_pcb "
            '(via (at 1 2) (size 0.45) (drill 0.20) (net 3) '
            '(layers "F.Cu" "B.Cu")) '
            '(via (at 1.48 2) (size 0.45) (drill 0.20) (net 3) '
            '(layers "F.Cu" "B.Cu")))'
        )
        self.assertEqual(self.validate_source(source), 2)

    def test_accepts_root_sheet_aliases_of_the_same_via_net(self):
        root = sexp.parse_one(
            '(kicad_pcb '
            '(via (at 1 2) (size 0.45) (drill 0.20) (net "/GND")) '
            '(via (at 1.48 2) (size 0.45) (drill 0.20) (net "GND")))')
        self.assertEqual(fab.validate_via_geometry("fixture", min_vias=0, board=root), 2)

    def test_rejects_close_unnamed_vias(self):
        root = sexp.parse_one(
            '(kicad_pcb '
            '(via (at 1 2) (size 0.45) (drill 0.20) (net 0)) '
            '(via (at 1.48 2) (size 0.45) (drill 0.20) (net 0)))')
        with self.assertRaises(fab.ViaGeometryError):
            fab.validate_via_geometry("fixture", min_vias=0, board=root)

    def test_rejects_close_via_copper_on_different_nets(self):
        source = (
            "(kicad_pcb "
            '(via (at 1 2) (size 0.45) (drill 0.20) (net 3) '
            '(layers "F.Cu" "B.Cu")) '
            '(via (at 1.48 2) (size 0.45) (drill 0.20) (net 4) '
            '(layers "F.Cu" "B.Cu")))'
        )
        with self.assertRaisesRegex(
                fab.ViaGeometryError,
                "0.03 mm via-to-via copper spacing is below 0.15 mm"):
            self.validate_source(source)

    def test_rejects_board_without_vias(self):
        with self.assertRaisesRegex(
                fab.ViaGeometryError,
                r"board lists 0 vias, fewer than the 99"):
            self.validate_source("(kicad_pcb)",
                                 min_vias=fab.MIN_BOARD_VIAS)

    def test_committed_board_meets_the_via_floor(self):
        self.assertGreaterEqual(fab.validate_via_geometry(fab.PCB),
                                fab.MIN_BOARD_VIAS)


class ZoneGeometryTests(unittest.TestCase):
    def validate(self, min_thickness, thermal_gap, bridge_width="0.4",
                 connect_pads="", footprint="", fill='(fill yes (thermal_gap {gap}) '
                                      '(thermal_bridge_width {bridge}))'):
        with tempfile.TemporaryDirectory() as directory:
            pcb = Path(directory) / "test.kicad_pcb"
            pcb.write_text(
                '(kicad_pcb (zone (net "GND") (name "GND_IN1") '
                f'{connect_pads}'
                f'(min_thickness {min_thickness}) '
                + fill.format(gap=thermal_gap, bridge=bridge_width) + ")" + footprint + ")",
                encoding="utf-8")
            return fab.validate_zone_geometry(pcb, min_pours=1, plane_layers=())

    def test_accepts_zone_feature_minimums(self):
        self.assertEqual(self.validate("0.13", "0.3"), 1)

    def test_rejects_solid_through_hole_overrides(self):
        for footprint_mode, pad_mode in (("", "2"), ("2", "")):
            with self.subTest(footprint=footprint_mode, pad=pad_mode), \
                    self.assertRaisesRegex(fab.ZoneGeometryError, "zone_connect 2"):
                self.validate("0.13", "0.3", footprint=
                              f'(footprint "Terminal" (property "Reference" "J1") '
                              f'{"(zone_connect " + footprint_mode + ")" if footprint_mode else ""}'
                              f'(pad "1" thru_hole circle (net "GND") '
                              f'{"(zone_connect " + pad_mode + ")" if pad_mode else ""}))')

    def test_accepts_relief_override_and_unpoured_or_smd_pads(self):
        for pad_type, net, pad_mode in (("thru_hole", "GND", "0"),
                                        ("thru_hole", "GND", "1"),
                                        ("thru_hole", "OTHER", "2"),
                                        ("smd", "GND", "2")):
            with self.subTest(pad_type=pad_type, net=net, pad_mode=pad_mode):
                self.assertEqual(self.validate("0.13", "0.3", footprint=
                    f'(footprint "Terminal" (zone_connect 2) '
                    f'(pad "1" {pad_type} circle (net "{net}") '
                    f'(zone_connect {pad_mode})))'), 1)

    def test_rejects_solid_through_hole_overrides_on_numeric_nets(self):
        for zone_name in ('(net_name "GND")', ''):
            for pad_net in ('(net 23 "GND")', '(net 23)'):
                for footprint_mode, pad_mode in (("", "2"), ("2", "")):
                    with self.subTest(zone_name=zone_name, pad_net=pad_net,
                                      footprint=footprint_mode, pad=pad_mode):
                        root = sexp.parse_one(
                            '(kicad_pcb (net 23 "GND") '
                            f'(zone (net 23) {zone_name} (min_thickness 0.13) '
                            '(fill yes (thermal_gap 0.3) (thermal_bridge_width 0.4))))')
                        root.append(sexp.parse_one(
                            '(footprint "Terminal" (property "Reference" "J1") '
                            f'{"(zone_connect " + footprint_mode + ")" if footprint_mode else ""}'
                            f'(pad "1" thru_hole circle {pad_net} '
                            f'{"(zone_connect " + pad_mode + ")" if pad_mode else ""}))'))
                        with self.assertRaisesRegex(fab.ZoneGeometryError, "zone_connect 2"):
                            fab.validate_zone_geometry("fixture", min_pours=1, board=root,
                                                       plane_layers=())

    def test_rejects_thin_pour_sliver(self):
        with self.assertRaisesRegex(
                fab.ZoneGeometryError,
                "GND_IN1: min_thickness 0.1016 mm is below 0.13 mm"):
            self.validate("0.1016", "0.3")

    def test_rejects_unresolvable_thermal_gap(self):
        with self.assertRaisesRegex(
                fab.ZoneGeometryError,
                "GND_IN1: thermal_gap 0.0254 mm is below 0.3 mm"):
            self.validate("0.13", "0.0254")

    def test_rejects_a_thermal_gap_at_the_bare_process_floor(self):
        with self.assertRaisesRegex(fab.ZoneGeometryError, "below 0.3 mm"):
            self.validate("0.13", "0.1016")

    def test_rejects_thin_thermal_spoke(self):
        with self.assertRaisesRegex(
                fab.ZoneGeometryError,
                "GND_IN1: thermal_bridge_width 0.1016 mm is below 0.4 mm"):
            self.validate("0.13", "0.3", "0.1016")

    def test_accepts_the_thermal_relief_default(self):
        self.assertEqual(
            self.validate("0.13", "0.3",
                          connect_pads="(connect_pads (clearance 0.1016)) "), 1)

    def test_accepts_reliefs_on_through_hole_pads_only(self):
        self.assertEqual(
            self.validate("0.13", "0.3",
                          connect_pads="(connect_pads thru_hole_only "
                                       "(clearance 0.1016)) "), 1)

    def test_rejects_solidly_connected_pads(self):
        with self.assertRaisesRegex(
                fab.ZoneGeometryError,
                "GND_IN1: connect_pads yes solders every pad"):
            self.validate("0.13", "0.3",
                          connect_pads="(connect_pads yes (clearance 0.1016)) ")

    def test_rejects_an_unfilled_pour(self):
        with self.assertRaisesRegex(
                fab.ZoneGeometryError, "GND_IN1: fill no pours no copper"):
            self.validate("0.13", "0.3",
                          fill="(fill no (thermal_gap {gap}) "
                               "(thermal_bridge_width {bridge}))")

    def test_rejects_a_pour_with_no_fill_enable_token(self):
        with self.assertRaisesRegex(
                fab.ZoneGeometryError, "GND_IN1: fill no pours no copper"):
            self.validate("0.13", "0.3",
                          fill="(fill (thermal_gap {gap}) "
                               "(thermal_bridge_width {bridge}))")

    def test_rejects_a_pour_with_no_fill_node(self):
        with self.assertRaisesRegex(
                fab.ZoneGeometryError,
                "GND_IN1: no .fill .... node"):
            self.validate("0.13", "0.3", fill="")

    def keepout_count(self, net_node):
        with tempfile.TemporaryDirectory() as directory:
            pcb = Path(directory) / "test.kicad_pcb"
            pcb.write_text(
                f'(kicad_pcb (zone {net_node}(name "H1 screw head") '
                '(min_thickness 0.0254) (keepout (tracks not_allowed))))',
                encoding="utf-8")
            return fab.validate_zone_geometry(pcb, min_pours=0, plane_layers=())

    def test_ignores_keepout_zones(self):
        self.assertEqual(self.keepout_count(""), 0)

    def test_ignores_keepout_zones_carrying_a_net(self):
        self.assertEqual(self.keepout_count('(net 0) (net_name "") '), 0)

    def test_rejects_board_without_copper_pours(self):
        with tempfile.TemporaryDirectory() as directory:
            pcb = Path(directory) / "test.kicad_pcb"
            pcb.write_text("(kicad_pcb)", encoding="utf-8")
            with self.assertRaisesRegex(
                    fab.ZoneGeometryError,
                    r"board lists 0 copper pours, fewer than the 2"):
                fab.validate_zone_geometry(pcb)

    def test_committed_board_pours_the_reference_planes(self):
        self.assertEqual(fab.validate_zone_geometry(fab.PCB),
                         fab.MIN_COPPER_POURS)

    def test_duplicate_and_outer_pours_do_not_replace_inner_planes(self):
        stack = '(layers (0 "F.Cu" signal) (2 "In1.Cu" power) (4 "In2.Cu" power) (31 "B.Cu" signal))'
        zone = ('(zone (net "GND") (layer "{layer}") (min_thickness 0.25) '
                '(fill yes (thermal_gap 0.5) (thermal_bridge_width 0.5)))')
        for first, second in (("In1.Cu", "In1.Cu"), ("F.Cu", "B.Cu")):
            with self.subTest(first=first, second=second):
                board = sexp.parse_one(
                    f'(kicad_pcb {stack} {zone.format(layer=first)} {zone.format(layer=second)})')
                with self.assertRaisesRegex(fab.ZoneGeometryError,
                                            r"required plane In2\.Cu has no copper pour"):
                    fab.validate_zone_geometry("fixture", board=board)


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


class SharedBoardParseTests(unittest.TestCase):
    """One parse feeds every geometry gate; no gate re-reads the file."""

    def test_gates_accept_a_preparsed_board(self):
        board = fab.read_board(fab.PCB)
        missing = "does-not-exist.kicad_pcb"

        self.assertEqual(fab.validate_plot_origin(missing, board=board),
                         (0.0, 0.0))
        self.assertGreaterEqual(
            fab.validate_via_geometry(missing, board=board),
            fab.MIN_BOARD_VIAS)
        self.assertEqual(fab.validate_zone_geometry(missing, board=board),
                         fab.MIN_COPPER_POURS)

    def test_rejects_unreadable_board(self):
        with self.assertRaisesRegex(fab.BoardReadError, "cannot read board"):
            fab.read_board("does-not-exist.kicad_pcb")


class AssemblyMetadataTests(unittest.TestCase):
    def setUp(self):
        self.ref = "U1"
        self.posrows = {
            self.ref: {
                "PosX": " 1.25 ",
                "PosY": "-2",
                "Rot": "370",
                "Side": " TOP ",
            }
        }

    def validate(self):
        return fab.validate_assembly_metadata(self.posrows, [self.ref])

    def assert_diagnostic(self, expected):
        with self.assertRaises(fab.AssemblyMetadataError) as caught:
            self.validate()
        self.assertIn(expected, caught.exception.diagnostics)

    def test_accepts_complete_metadata_and_normalizes_text(self):
        metadata = self.validate()

        self.assertEqual(
            metadata[self.ref],
            {
                "lcsc": fab.LCSC_BY_REF[self.ref],
                "pos_x": "1.25",
                "pos_y": "-2",
                "rotation": 370.0,
                "side": "top",
            },
        )

    def test_rejects_a_reference_without_a_part_assignment(self):
        self.ref = "X1"
        self.posrows = {self.ref: dict(self.posrows["U1"])}

        self.assert_diagnostic(
            f"{self.ref}: supplier part number (LCSC) is missing")

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
        self.ref = "X1"
        self.posrows = {}

        with self.assertRaises(fab.AssemblyMetadataError) as caught:
            self.validate()
        self.assertEqual(
            caught.exception.diagnostics,
            (
                f"{self.ref}: supplier part number (LCSC) is missing",
                f"{self.ref}: centroid row is missing",
            ),
        )


class AssemblyExclusionTests(unittest.TestCase):
    """fab.py's assembly set and the board's own attributes must agree."""

    BOARD = sexp.parse_one(
        "(kicad_pcb"
        ' (footprint "L:smd" (attr smd)'
        ' (property "Reference" "R1"))'
        ' (footprint "L:header" (attr through_hole exclude_from_pos_files'
        ' exclude_from_bom) (property "Reference" "J1")))')
    COMPS = {"R1": {}, "J1": {}}

    def test_accepts_an_assembly_set_matching_the_board(self):
        self.assertEqual(
            fab.validate_assembly_exclusions(self.COMPS, ["R1"], self.BOARD), 1)

    def test_rejects_assembling_a_part_the_board_excludes(self):
        with self.assertRaisesRegex(
                fab.AssemblyMetadataError,
                "assembled parts the board excludes from the BOM or "
                "centroid: J1"):
            fab.validate_assembly_exclusions(self.COMPS, ["R1", "J1"],
                                             self.BOARD)

    def test_rejects_holding_back_a_part_the_board_assembles(self):
        with self.assertRaisesRegex(
                fab.AssemblyMetadataError,
                "parts held back from assembly that the board does not "
                "exclude: R1"):
            fab.validate_assembly_exclusions(self.COMPS, [], self.BOARD)


class AssemblyPolicyTests(unittest.TestCase):
    """JLC reflows top-side SMD only; every hand-soldered part stays out."""

    NATIVE_EXCLUDED_REFS = {"U_MCU", "J1", "J2", "J3A", "J3B", "J4", "C_IN"}
    NATIVE_EXCLUSION_FLAGS = {"exclude_from_pos_files", "exclude_from_bom"}

    def component(self, footprint, value="", dnp=False):
        return {"footprint": footprint, "value": value, "dnp": dnp}

    def test_excludes_hand_soldered_packages(self):
        for footprint in (
                "Connector_PinHeader_2.54mm:PinHeader_1x03_P2.54mm_Vertical",
                "phantasm:TerminalBlock_GCT_TBC05-03-1-G-G",
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

    def test_netlist_parser_preserves_dnp_properties(self):
        with tempfile.TemporaryDirectory() as directory:
            netlist = Path(directory) / "parts.net"
            netlist.write_text('''(export (components
                (comp (ref "R1") (value "10k")
                    (footprint "Resistor_SMD:R_0603_1608Metric")
                    (property (name "dnp")))
                (comp (ref "R2") (value "10k")
                    (footprint "Resistor_SMD:R_0603_1608Metric"))))''', encoding="utf-8")
            components = fab.parse_components(netlist)
        self.assertEqual(components, {
            "R1": self.component("Resistor_SMD:R_0603_1608Metric", "10k", dnp=True),
            "R2": self.component("Resistor_SMD:R_0603_1608Metric", "10k"),
        })
        self.assertFalse(fab.is_assembled(components["R1"]))
        self.assertTrue(fab.is_assembled(components["R2"]))

    def test_board_marks_hand_soldered_parts_as_excluded(self):
        root = fab.read_board(fab.PCB)
        attrs_by_ref = {}
        for ref, footprint in footprints_by_ref(root).items():
            if ref in self.NATIVE_EXCLUDED_REFS:
                attrs_by_ref[ref] = {
                    str(value) for value in sexp.val(footprint, "attr", [])
                }

        self.assertEqual(set(attrs_by_ref), self.NATIVE_EXCLUDED_REFS)
        for ref, attrs in attrs_by_ref.items():
            with self.subTest(ref=ref):
                self.assertLessEqual(self.NATIVE_EXCLUSION_FLAGS, attrs)

    def test_board_excludes_no_assembled_part(self):
        excluded = fab.board_assembly_exclusions(fab.read_board(fab.PCB))
        self.assertLessEqual(self.NATIVE_EXCLUDED_REFS, excluded)
        self.assertEqual(excluded & set(fab.LCSC_BY_REF), set())

    def test_accepts_exact_assigned_part_set(self):
        fab.validate_assembled_refs(fab.LCSC_BY_REF, "1.2")

    def test_revision_1_1_assembly_has_no_transmit_pulldown(self):
        assembled = set(fab.LCSC_BY_REF) - {"R_TX"}
        fab.validate_assembled_refs(assembled, "1.1")
        with self.assertRaisesRegex(fab.AssemblyMetadataError,
                                    "assigned parts missing from assembly: R_TX"):
            fab.validate_assembled_refs(assembled, "1.2")

    def test_rejects_assigned_part_missing_from_assembly(self):
        assembled = set(fab.LCSC_BY_REF) - {"U1"}

        with self.assertRaisesRegex(
                fab.AssemblyMetadataError,
                "assigned parts missing from assembly: U1"):
            fab.validate_assembled_refs(assembled, "1.2")

    def test_rejects_rotation_correction_missing_from_assembly(self):
        assembled = set(fab.LCSC_BY_REF) - {"U1"}

        with self.assertRaisesRegex(
                fab.AssemblyMetadataError,
                "rotation corrections missing from assembly: U1"):
            fab.validate_rotation_refs(assembled)

    def test_d_bus_exports_committed_rotation_and_cathode_net(self):
        footprint = footprints_by_ref(fab.read_board(fab.PCB))["D_BUS"]
        rotation = float(sexp.val(footprint, "at")[2])
        self.assertEqual(rotation, 90)
        self.assertEqual(fab.cpl_rotation("D_BUS", rotation), 90)
        self.assertEqual(
            {str(pad[1]): sexp.val(pad, "net")[-1]
             for pad in F(footprint, "pad")},
            {"1": "/SYNC_BUS", "2": "GND"})

    def test_nonpolarized_resistor_rotation_is_unchanged(self):
        self.assertEqual(fab.cpl_rotation("R_D1", 180), 180)

class SchematicParityTests(unittest.TestCase):

    KNOWN = [
        {"type": "extra_footprint", "description": "Extra footprint",
         "items": [{"description": f"Footprint H{index}"}]}
        for index in range(1, 5)
    ]


    KNOWN_WARNINGS = [
        {"type": kind, "description": "Known warning",
         "items": [{"description": f"Footprint {ref}"}]}
        for (kind, ref), count in sorted(fab.KNOWN_PARITY_WARNING_COUNTS.items())
        for _ in range(count)
    ]


    def require_report(self, report):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "parity.json"
            path.write_text(json.dumps(report), encoding="utf-8")
            return fab.require_schematic_parity(path)


    def require(self, entries, violations=None):
        return self.require_report({
            "schematic_parity": entries,
            "violations": (self.KNOWN_WARNINGS if violations is None
                           else violations),
        })

    def test_nonzero_exit_fails_even_with_clean_report(self):
        with mock.patch.object(fab, "kicad_cli", return_value="kicad-cli"), \
                mock.patch.object(fab, "run", return_value=mock.Mock(returncode=2)), \
                mock.patch.object(fab, "require_schematic_parity", return_value=0):
            with self.assertRaisesRegex(fab.SchematicParityError, "status 2"):
                fab.run_parity("missing-parity-fixture.json")

    def test_committed_artifacts_pass_pinned_kicad_parity(self):
        try:
            fab.kicad_cli()
        except SystemExit as error:
            self.skipTest(str(error))
        with tempfile.TemporaryDirectory() as directory:
            report = str(Path(directory) / "parity.json")
            self.assertEqual(fab.run_parity(report), len(fab.KNOWN_PARITY_ITEMS))


    def test_accepts_known_differences(self):
        self.assertEqual(self.require(self.KNOWN), len(self.KNOWN))

    def test_rejects_report_without_a_parity_section(self):
        with self.assertRaisesRegex(
                fab.SchematicParityError,
                "no schematic_parity section"):
            self.require_report({"coordinate_units": "mm"})

    def test_rejects_report_without_a_violations_section(self):
        with self.assertRaisesRegex(
                fab.SchematicParityError,
                "no violations section"):
            self.require_report({"schematic_parity": self.KNOWN})

    def test_rejects_unexpected_warning_type(self):
        violations = self.KNOWN_WARNINGS + [{
            "type": "via_dangling",
            "description": "Through via has no connection",
        }]

        with self.assertRaisesRegex(
                fab.SchematicParityError,
                r"via_dangling: .* reported 1 times .* expected 0") as caught:
            self.require(self.KNOWN, violations)
        self.assertIn("resolve in Pcbnew", str(caught.exception))
        self.assertNotIn("re-route", str(caught.exception))

    def test_rejects_known_warning_count_drift(self):
        with self.assertRaisesRegex(
                fab.SchematicParityError,
                r"lib_footprint_mismatch: .* reported 0 times .* expected 1") as caught:
            self.require(self.KNOWN, self.KNOWN_WARNINGS[:-1])
        self.assertIn("Land edits", str(caught.exception))
        self.assertNotIn("re-route", str(caught.exception))

    def test_optional_warning_accepts_up_to_its_maximum(self):
        extra = {"type": "lib_footprint_mismatch", "description": "Optional",
                 "items": [{"description": "Footprint J1"}]}
        with mock.patch.object(fab, "OPTIONAL_PARITY_WARNING_MAX",
                               {("lib_footprint_mismatch", "J1"): 1}):
            for count in (0, 1):
                self.assertEqual(
                    self.require(self.KNOWN, self.KNOWN_WARNINGS + [extra] * count),
                    len(fab.KNOWN_PARITY_ITEMS))
            with self.assertRaisesRegex(fab.SchematicParityError,
                                        r"J1 reported 2 times .* expected at most 1"):
                self.require(self.KNOWN, self.KNOWN_WARNINGS + [extra] * 2)

    def test_rejects_replacement_warning_with_same_total(self):
        violations = self.KNOWN_WARNINGS[:-1] + [{
            "type": "lib_footprint_mismatch", "description": "New warning",
            "items": [{"description": "Footprint C_DEC1"}]}]
        with self.assertRaisesRegex(fab.SchematicParityError,
                                    r"C_DEC1 reported 1 times .* expected 0"):
            self.require(self.KNOWN, violations)

    def test_rejects_report_missing_known_differences(self):
        with self.assertRaisesRegex(
                fab.SchematicParityError,
                r"extra_footprint: H1 reported 0 times"):
            self.require([])

    def test_rejects_a_duplicate_standing_in_for_a_vanished_item(self):
        entries = [e for e in self.KNOWN
                   if e["items"][0]["description"] != "Footprint H2"]
        entries.append({"type": "extra_footprint",
                        "description":
                            fab.KNOWN_PARITY_ITEMS[("extra_footprint", "H1")],
                        "items": [{"description": "Footprint H1"}]})

        with self.assertRaisesRegex(
                fab.SchematicParityError,
                r"extra_footprint: H2 reported 0 times"):
            self.require(entries)

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
                fab.SchematicParityError, "Update PCB from Schematic"):
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
            with mock.patch.object(
                    fab, "kicad_cli", return_value="kicad-cli"), \
                    mock.patch.object(fab, "run", fake_run):
                with self.assertRaisesRegex(
                        fab.DesignRuleError, "kicad-cli drc exited 5"):
                    fab.run_drc(path)


class PartCatalogTests(unittest.TestCase):
    def test_passive_catalog_values_match_committed_schematics(self):
        def normalized(value):
            token = value.split()[0].removesuffix("ohm").removesuffix("F")
            match = re.fullmatch(r"(\d+(?:\.\d+)?)([RkMpun]?)(\d*)", token)
            self.assertIsNotNone(match, value)
            whole, prefix, fraction = match.groups()
            number = Decimal(whole + "." + fraction if fraction else whole)
            return number * {"": Decimal(1), "R": Decimal(1), "k": Decimal("1e3"),
                             "M": Decimal("1e6"), "p": Decimal("1e-12"),
                             "n": Decimal("1e-9"), "u": Decimal("1e-6")}[prefix]

        for value, expected in (("1R5", "1.5"), ("10k", "10000"),
                                ("0.1uF", "1e-7"), ("100nF", "1e-7")):
            self.assertEqual(normalized(value), Decimal(expected))
        expected_refs = {ref for ref in fab.LCSC_BY_REF if ref.startswith(("R", "C"))}
        for revision in ("1.1", "1.2"):
            schematic = GEN.parent / revision / "phantasm.kicad_sch"
            root = sexp.parse_one(schematic.read_text(encoding="utf-8"))
            checked = set()
            for symbol in F(root, "symbol"):
                ref = property_value(symbol, "Reference", "")
                if ref not in expected_refs:
                    continue
                checked.add(ref)
                value = property_value(symbol, "Value")
                part = fab.PART_BY_LCSC[fab.LCSC_BY_REF[ref]]
                with self.subTest(revision=revision, ref=ref, value=value, lcsc=fab.LCSC_BY_REF[ref]):
                    self.assertEqual(normalized(value), normalized(part["description"]))
            self.assertEqual(expected_refs - checked, {"R_TX"} if revision == "1.1" else set(), revision)

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

    def test_rejects_a_catalog_entry_no_reference_names(self):
        catalog = dict(fab.PART_BY_LCSC)
        catalog["C999999999"] = {
            "manufacturer": "ACME",
            "mpn": "NOPART",
            "description": "unassigned",
        }
        metadata = {
            ref: {"lcsc": lcsc}
            for ref, lcsc in fab.LCSC_BY_REF.items()
        }
        with mock.patch.object(fab, "PART_BY_LCSC", catalog):
            with self.assertRaisesRegex(
                    fab.PartCatalogError,
                    "LCSC C999999999: catalog entry no reference is assigned"):
                fab.validate_part_catalog(metadata)


class CommandRunTests(unittest.TestCase):
    def test_netlist_export_checks_annotation(self):
        result = subprocess.CompletedProcess(["kicad-cli"], 0, "", "")
        with mock.patch.object(fab, "run", return_value=result), \
                mock.patch.object(fab, "require_annotated_export") as annotated:
            self.assertIs(fab.run_export("netlist", ["kicad-cli", "source.kicad_sch"]), result)
        annotated.assert_called_once_with(result, "source.kicad_sch")

    def test_nonraising_failure_surfaces_stderr(self):
        result = subprocess.CompletedProcess(["kicad-cli"], 1, "", "cannot load board\n")
        stderr = io.StringIO()
        with mock.patch.object(subprocess, "run", return_value=result), \
                contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(stderr):
            self.assertIs(fab.run(["kicad-cli"], check=False), result)
        self.assertEqual(stderr.getvalue(), "cannot load board\n")

    def invoke(self, error):
        def fake_run(*args, **kw):
            raise error

        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout), \
                mock.patch.object(subprocess, "run", fake_run):
            with self.assertRaises(SystemExit) as caught:
                fab.run(["kicad-cli", "version"])
        return str(caught.exception.code)

    def test_missing_binary_reports_a_one_line_diagnostic(self):
        message = self.invoke(
            FileNotFoundError(2, "No such file or directory"))
        self.assertEqual(message.splitlines(), [message])
        self.assertIn("kicad-cli", message)
        self.assertIn("KICAD_CLI", message)


class CommandLineTests(unittest.TestCase):
    def test_verify_reports_the_matching_digest_count(self):
        output = io.StringIO()
        with mock.patch.object(fab, "verify_package", return_value=7) as verify, \
                contextlib.redirect_stdout(output):
            fab.verify_main()
        verify.assert_called_once_with(fab.JLC, fab.SHIPPED_SUMS)
        self.assertIn(f"fab package verified: 7 digests match {fab.SHIPPED_SUMS}", output.getvalue())

    def test_verify_exits_with_the_package_error(self):
        with mock.patch.object(fab, "verify_package",
                               side_effect=fab.PackageVerificationError("digest mismatch")):
            with self.assertRaises(SystemExit) as caught:
                fab.verify_main()
        self.assertEqual(caught.exception.code, "digest mismatch")

    def test_fabrication_uses_the_requested_revision(self):
        with tempfile.TemporaryDirectory() as directory, \
                mock.patch.object(fab, "OUT", directory), \
                mock.patch.object(fab, "fabricate") as fabricate:
            fab.main("1.2")
        fabricate.assert_called_once_with("1.2")

    def parse(self, argv):
        stdout, stderr = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(stdout), \
                contextlib.redirect_stderr(stderr):
            with self.assertRaises(SystemExit) as caught:
                fab.parse_args(argv)
        return caught.exception.code, stdout.getvalue() + stderr.getvalue()

    def test_accepts_no_arguments(self):
        self.assertEqual(vars(fab.parse_args([])), {"verify": False, "revision": "1.1"})

    def test_explicit_revision_1_1_preserves_the_export_target(self):
        self.assertEqual(fab.parse_args(["--revision", "1.1"]).revision, "1.1")

    def test_other_revisions_are_rejected_by_the_parser(self):
        for revision in ("1.2", "1.3"):
            with self.subTest(revision=revision):
                code, text = self.parse(["--revision", revision])
                self.assertEqual(code, 2)
                self.assertIn("revision-specific BOM are not validated", text)
                self.assertIn("rev 1.1 only", text)

    def test_verify_selects_the_digest_check(self):
        self.assertTrue(fab.parse_args(["--verify"]).verify)

    def test_help_exits_without_running_the_workflow(self):
        code, text = self.parse(["--help"])
        self.assertEqual(code, 0)
        self.assertIn("usage:", text)

    def test_rejects_unknown_option(self):
        code, text = self.parse(["--refresh"])
        self.assertEqual(code, 2)
        self.assertIn("unrecognized arguments: --refresh", text)

    def test_rejects_positional_argument(self):
        self.assertEqual(self.parse(["board.kicad_pcb"])[0], 2)


class NetlistSpecTests(unittest.TestCase):
    """The fab run holds its own netlist export to check.py's named-net table."""

    def validate(self, nets):
        blocks = "".join(
            '(net (name "{}") {})'.format(
                name, " ".join(f'(node (ref "{ref}") (pin "{pin}"))'
                               for ref, pin in nodes))
            for name, nodes in nets.items())
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "test.net"
            path.write_text(
                '(export (design (sheet (name "/") (title_block (rev "1.2"))))'
                f'(nets {blocks}))', encoding="utf-8")
            with contextlib.redirect_stdout(io.StringIO()):
                return fab.validate_netlist_spec(path)

    def test_accepts_the_specified_partition(self):
        self.assertEqual(self.validate(expected_nodes()),
                         len(fab.netlist_spec.EXPECT_REV_1_2))

    def test_rejects_a_net_missing_a_member(self):
        nets = expected_nodes()
        nets["SYNC_BUS"].pop()
        with self.assertRaisesRegex(fab.NetlistSpecError,
                                    "does not match the electrical"):
            self.validate(nets)

    def test_rejects_a_permuted_pinout(self):
        nets = expected_nodes()
        nets["DATA_IN"] = [("U_MCU", "11"), ("U1", "5")]
        with self.assertRaisesRegex(fab.NetlistSpecError,
                                    "does not match the electrical"):
            self.validate(nets)


class ElectricalRuleTests(unittest.TestCase):
    def require_report(self, document):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "erc.json"
            path.write_text(json.dumps(document), encoding="utf-8")
            return fab.require_clean_erc(path)

    def test_accepts_clean_sheets(self):
        self.assertEqual(self.require_report({"sheets": [{"violations": []}]}), 0)

    def test_rejects_violation_on_later_sheet(self):
        with self.assertRaisesRegex(fab.ElectricalRuleError, "1 error-severity"):
            self.require_report({"sheets": [{"violations": []},
                                           {"violations": [{"type": "pin_not_connected"}]}]})

    def test_rejects_malformed_sections(self):
        for document in [[], {}, {"sheets": []}, {"sheets": [{}]},
                         {"sheets": [{"violations": None}]}, {"sheets": [None]}]:
            with self.subTest(document=document), self.assertRaises(fab.ElectricalRuleError):
                self.require_report(document)

    def test_missing_report_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(fab.ElectricalRuleError, "cannot read ERC"):
                fab.require_clean_erc(Path(directory) / "missing.json")

    def test_runner_rejects_stale_reports_and_nonzero_clean_exit(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "erc.json"
            path.write_text('{"sheets": [{"violations": []}]}', encoding="utf-8")
            def fake_run(args, check=True):
                self.assertFalse(path.exists())
                self.assertFalse(check)
                self.assertEqual(args[1:8], ["sch", "erc", "--severity-error", "--format",
                                             "json", "--exit-code-violations", "-o"])
                path.write_text('{"sheets": [{"violations": []}]}', encoding="utf-8")
                return subprocess.CompletedProcess(args, 5)
            with mock.patch.object(fab, "kicad_cli", return_value="kicad-cli"), \
                    mock.patch.object(fab, "run", side_effect=fake_run):
                with self.assertRaisesRegex(fab.ElectricalRuleError, "erc exited 5"):
                    fab.run_erc(path)


if __name__ == "__main__":
    unittest.main()
