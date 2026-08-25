import contextlib
import hashlib
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
import unittest.mock
import zipfile
from pathlib import Path

GEN = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN))

import fab  # noqa: E402


class ProjectRulesTests(unittest.TestCase):
    def write_project(self, directory, rules=None, default=None):
        document = {
            "board": {"design_settings": {"rules": dict(fab.RULE_MINIMUMS)}},
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
            self.assertEqual(
                fab.validate_project_rules(project),
                len(fab.RULE_MINIMUMS) + len(fab.DEFAULT_CLASS_MINIMUMS))

    def test_rejects_gui_rezeroed_clearance(self):
        with tempfile.TemporaryDirectory() as directory:
            project = self.write_project(directory, rules={"min_clearance": 0})
            with self.assertRaises(fab.ProjectRulesError) as caught:
                fab.validate_project_rules(project)
        self.assertIn("min_clearance", str(caught.exception))
        self.assertIn("heal_clearance.py", str(caught.exception))

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
        self.assertEqual(
            fab.validate_project_rules(),
            len(fab.RULE_MINIMUMS) + len(fab.DEFAULT_CLASS_MINIMUMS))


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

    def test_spacing_threshold_is_named_for_via_pairs(self):
        self.assertEqual(fab.MIN_VIA_TO_VIA_COPPER_SPACING_MM, 0.15)
        self.assertFalse(hasattr(fab, "MIN_VIA_COPPER_SPACING_MM"))

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
    def validate(self, min_thickness, thermal_gap, bridge_width="0.13",
                 connect_pads="", fill='(fill yes (thermal_gap {gap}) '
                                      '(thermal_bridge_width {bridge}))'):
        with tempfile.TemporaryDirectory() as directory:
            pcb = Path(directory) / "test.kicad_pcb"
            pcb.write_text(
                '(kicad_pcb (zone (net "GND") (name "GND_IN1") '
                f'{connect_pads}'
                f'(min_thickness {min_thickness}) '
                + fill.format(gap=thermal_gap, bridge=bridge_width) + "))",
                encoding="utf-8")
            return fab.validate_zone_geometry(pcb, min_pours=1)

    def test_accepts_zone_feature_minimums(self):
        self.assertEqual(self.validate("0.13", "0.1016"), 1)

    def test_rejects_thin_pour_sliver(self):
        with self.assertRaisesRegex(
                fab.ZoneGeometryError,
                "GND_IN1: min_thickness 0.1016 mm is below 0.13 mm"):
            self.validate("0.1016", "0.1016")

    def test_rejects_unresolvable_thermal_gap(self):
        with self.assertRaisesRegex(
                fab.ZoneGeometryError,
                "GND_IN1: thermal_gap 0.0254 mm is below 0.1016 mm"):
            self.validate("0.13", "0.0254")

    def test_rejects_thin_thermal_spoke(self):
        with self.assertRaisesRegex(
                fab.ZoneGeometryError,
                "GND_IN1: thermal_bridge_width 0.1016 mm is below 0.13 mm"):
            self.validate("0.13", "0.1016", "0.1016")

    def test_accepts_the_thermal_relief_default(self):
        self.assertEqual(
            self.validate("0.13", "0.1016",
                          connect_pads="(connect_pads (clearance 0.1016)) "), 1)

    def test_accepts_reliefs_on_through_hole_pads_only(self):
        self.assertEqual(
            self.validate("0.13", "0.1016",
                          connect_pads="(connect_pads thru_hole_only "
                                       "(clearance 0.1016)) "), 1)

    def test_rejects_solidly_connected_pads(self):
        with self.assertRaisesRegex(
                fab.ZoneGeometryError,
                "GND_IN1: connect_pads yes solders every pad"):
            self.validate("0.13", "0.1016",
                          connect_pads="(connect_pads yes (clearance 0.1016)) ")

    def test_rejects_an_unfilled_pour(self):
        with self.assertRaisesRegex(
                fab.ZoneGeometryError, "GND_IN1: fill no pours no copper"):
            self.validate("0.13", "0.1016",
                          fill="(fill no (thermal_gap {gap}) "
                               "(thermal_bridge_width {bridge}))")

    def test_rejects_a_pour_with_no_fill_enable_token(self):
        with self.assertRaisesRegex(
                fab.ZoneGeometryError, "GND_IN1: fill no pours no copper"):
            self.validate("0.13", "0.1016",
                          fill="(fill (thermal_gap {gap}) "
                               "(thermal_bridge_width {bridge}))")

    def test_rejects_a_pour_with_no_fill_node(self):
        with self.assertRaisesRegex(
                fab.ZoneGeometryError,
                "GND_IN1: no .fill .... node"):
            self.validate("0.13", "0.1016", fill="")

    def keepout_count(self, net_node):
        with tempfile.TemporaryDirectory() as directory:
            pcb = Path(directory) / "test.kicad_pcb"
            pcb.write_text(
                f'(kicad_pcb (zone {net_node}(name "H1 screw head") '
                '(min_thickness 0.0254) (keepout (tracks not_allowed))))',
                encoding="utf-8")
            return fab.validate_zone_geometry(pcb, min_pours=0)

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


class AssemblyPolicyTests(unittest.TestCase):
    """JLC reflows top-side SMD only; every hand-soldered part stays out."""

    NATIVE_EXCLUDED_REFS = {"U_MCU", "J1", "J2", "J3A", "J3B", "J4", "C_IN"}
    NATIVE_EXCLUSION_FLAGS = {"exclude_from_pos_files", "exclude_from_bom"}

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

    def test_board_marks_hand_soldered_parts_as_excluded(self):
        root = fab.read_board(fab.PCB)
        attrs_by_ref = {}
        for footprint in fab.F(root, "footprint"):
            ref = next(
                (prop[2] for prop in fab.F(footprint, "property")
                 if prop[1] == "Reference"),
                None,
            )
            if ref in self.NATIVE_EXCLUDED_REFS:
                attrs_by_ref[ref] = {
                    str(value) for value in fab.sexp.val(footprint, "attr", [])
                }

        self.assertEqual(set(attrs_by_ref), self.NATIVE_EXCLUDED_REFS)
        for ref, attrs in attrs_by_ref.items():
            with self.subTest(ref=ref):
                self.assertLessEqual(self.NATIVE_EXCLUSION_FLAGS, attrs)

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

    def test_d_bus_cathode_rotation(self):
        self.assertEqual(fab.cpl_rotation("D_BUS", 90), 270)

    def test_nonpolarized_resistor_rotation_is_unchanged(self):
        self.assertEqual(fab.cpl_rotation("R_D1", 180), 180)

class SchematicParityTests(unittest.TestCase):
    KNOWN = [
        {"type": "extra_footprint", "description": "Extra footprint",
         "items": [{"description": f"Footprint H{index}"}]}
        for index in range(1, 5)
    ] + [
        {"type": "footprint_symbol_mismatch",
         "description": "Exclude from bill of materials",
         "items": [{"description": f"Footprint {ref}"}]}
        for ref in ("JP_ID0", "JP_ID1", "JP_ID2", "JP_SHLD")
    ]
    KNOWN_WARNINGS = [
        {"type": kind, "description": "Known warning"}
        for kind, count in sorted(fab.KNOWN_PARITY_WARNING_COUNTS.items())
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
                r"via_dangling: reported 1 times .* expected 0"):
            self.require(self.KNOWN, violations)

    def test_rejects_known_warning_count_drift(self):
        with self.assertRaisesRegex(
                fab.SchematicParityError,
                r"lib_footprint_mismatch: reported 10 times .* expected 11"):
            self.require(self.KNOWN, self.KNOWN_WARNINGS[:-1])

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
        with unittest.mock.patch.object(fab, "PART_BY_LCSC", catalog):
            with self.assertRaisesRegex(
                    fab.PartCatalogError,
                    "LCSC C999999999: catalog entry no reference is assigned"):
                fab.validate_part_catalog(metadata)


class CommandRunTests(unittest.TestCase):
    def invoke(self, error):
        def fake_run(*args, **kw):
            raise error

        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout), \
                unittest.mock.patch.object(subprocess, "run", fake_run):
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
    def parse(self, argv):
        stdout, stderr = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(stdout), \
                contextlib.redirect_stderr(stderr):
            with self.assertRaises(SystemExit) as caught:
                fab.parse_args(argv)
        return caught.exception.code, stdout.getvalue() + stderr.getvalue()

    def test_accepts_no_arguments(self):
        self.assertEqual(vars(fab.parse_args([])), {"verify": False})

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


class ZipMemberTests(unittest.TestCase):
    def build(self, mtime):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "phantasm-F_Cu.gtl"
            source.write_text("G04 gerber*\n", encoding="utf-8")
            os.utime(source, (mtime, mtime))
            buffer = io.BytesIO()
            with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
                archive.writestr(fab.zip_member(source.name), source.read_bytes())
            return buffer.getvalue()

    def test_archive_is_independent_of_source_mtime(self):
        self.assertEqual(self.build(315532800), self.build(1700000000))

    def test_member_carries_no_host_metadata(self):
        info = fab.zip_member("phantasm-F_Cu.gtl")
        self.assertEqual(info.date_time, (1980, 1, 1, 0, 0, 0))
        self.assertEqual(info.create_system, 3)
        self.assertEqual(info.compress_type, zipfile.ZIP_DEFLATED)


class TimestampNormalizationTests(unittest.TestCase):
    """Every artifact records the export wall clock, so without normalization
    two runs over an unchanged board differ in every file."""

    # Gerbers are CRLF, Excellon and the job file LF.
    GERBER = ("%TF.GenerationSoftware,KiCad,Pcbnew,10.0.4*%\r\n"
              "%TF.CreationDate,{stamp}-07:00*%\r\n"
              "G04 Created by KiCad (PCBNEW 10.0.4) date {date}*\r\n"
              "D10*\r\nX1000Y2000D02*\r\nM02*\r\n")
    DRILL = ("M48\n; DRILL file KiCad 10.0.4 date {stamp}\n"
             "; #@! TF.CreationDate,{stamp}-07:00\nFMAT,2\nM30\n")
    JOB = ('{{\n  "GeneralSpecs": {{\n'
           '    "CreationDate": "{stamp}-07:00"\n  }}\n}}\n')

    def export(self, stamp, date):
        directory = Path(self.enterContext(tempfile.TemporaryDirectory()))
        for name, template in (("phantasm-Top Layer.gtl", self.GERBER),
                               ("phantasm-PTH.drl", self.DRILL),
                               ("phantasm-job.gbrjob", self.JOB)):
            (directory / name).write_bytes(
                template.format(stamp=stamp, date=date).encode("utf-8"))
        rewritten = fab.normalize_fab_timestamps(directory)
        return directory, rewritten

    def test_two_exports_of_one_board_become_identical(self):
        first, rewritten = self.export("2026-08-12T17:25:19", "2026-08-12 17:25:19")
        second, _ = self.export("2019-03-04T01:02:03", "2019-03-04 01:02:03")
        self.assertEqual(len(rewritten), 3)
        for path in sorted(first.iterdir()):
            self.assertEqual(path.read_bytes(),
                             (second / path.name).read_bytes(), path.name)

    def test_line_endings_and_tool_version_survive(self):
        directory, _ = self.export("2026-08-12T17:25:19", "2026-08-12 17:25:19")
        gerber = (directory / "phantasm-Top Layer.gtl").read_bytes()
        self.assertNotIn(b"\n", gerber.replace(b"\r\n", b""))
        self.assertIn(f"%TF.CreationDate,{fab.FAB_TIMESTAMP}*%\r\n".encode(), gerber)
        self.assertIn(b"KiCad (PCBNEW 10.0.4) date", gerber)
        drill = (directory / "phantasm-PTH.drl").read_bytes()
        self.assertNotIn(b"\r", drill)
        self.assertIn(f"date {fab.FAB_TIMESTAMP}\n".encode(), drill)

    def test_an_undecodable_artifact_is_loud(self):
        directory = Path(self.enterContext(tempfile.TemporaryDirectory()))
        (directory / "phantasm-F_Cu.gtl").write_bytes(
            b"%TF.CreationDate,\xff\xfe*%\r\n")
        with self.assertRaisesRegex(fab.TimestampNormalizationError,
                                    "cannot read exported artifact"):
            fab.normalize_fab_timestamps(directory)

    def test_an_artifact_no_pattern_reaches_is_loud(self):
        directory = Path(self.enterContext(tempfile.TemporaryDirectory()))
        (directory / "phantasm-F_Cu.gtl").write_bytes(
            b"%TF.Creation-Date,2026-08-12T17:25:19-07:00*%\r\n")
        with self.assertRaisesRegex(fab.TimestampNormalizationError,
                                    "phantasm-F_Cu.gtl"):
            fab.normalize_fab_timestamps(directory)

    def test_a_respelled_banner_cannot_report_a_clean_run(self):
        directory, _ = self.export("2026-08-12T17:25:19", "2026-08-12 17:25:19")
        with unittest.mock.patch.object(fab, "TIMESTAMP_SUBSTITUTIONS", ()):
            with self.assertRaises(fab.TimestampNormalizationError) as caught:
                fab.normalize_fab_timestamps(directory)
        for name in sorted(path.name for path in directory.iterdir()):
            self.assertIn(name, str(caught.exception))

    def test_a_gerber_carrying_only_one_of_its_two_stamps_is_loud(self):
        directory = Path(self.enterContext(tempfile.TemporaryDirectory()))
        gerber = self.GERBER.format(stamp="2026-08-12T17:25:19",
                                    date="2026-08-12 17:25:19")
        (directory / "phantasm-F_Cu.gtl").write_bytes(
            gerber.replace("G04 Created by KiCad", "G04 Made by KiCad").encode())
        with self.assertRaisesRegex(fab.TimestampNormalizationError,
                                    "gerber banner"):
            fab.normalize_fab_timestamps(directory)

    def test_an_unknown_artifact_type_is_loud(self):
        directory = Path(self.enterContext(tempfile.TemporaryDirectory()))
        (directory / "phantasm-notes.txt").write_bytes(b"G04 gerber*\r\n")
        with self.assertRaisesRegex(fab.TimestampNormalizationError,
                                    "phantasm-notes.txt"):
            fab.normalize_fab_timestamps(directory)

    def test_every_zipped_artifact_type_states_its_stamps(self):
        for name in sorted(fab.ZIP_MEMBERS):
            with self.subTest(name=name):
                self.assertTrue(
                    fab.REQUIRED_STAMPS[os.path.splitext(name)[1]])

    def test_normalizing_an_already_stamped_export_stays_clean(self):
        directory, first = self.export("2026-08-12T17:25:19",
                                       "2026-08-12 17:25:19")
        self.assertEqual(fab.normalize_fab_timestamps(directory), first)


class ZipMembershipTests(unittest.TestCase):
    """A full export: GERBER_LAYERS in Protel naming, both drills, the job."""

    EXPORTED = ("phantasm-F_Cu.gtl", "phantasm-In1_Cu.g1", "phantasm-In2_Cu.g2",
                "phantasm-B_Cu.gbl", "phantasm-F_Silkscreen.gto",
                "phantasm-B_Silkscreen.gbo", "phantasm-F_Mask.gts",
                "phantasm-B_Mask.gbs", "phantasm-F_Paste.gtp",
                "phantasm-B_Paste.gbp", "phantasm-Edge_Cuts.gm1",
                "phantasm-PTH.drl", "phantasm-NPTH.drl", "phantasm-job.gbrjob")

    def test_the_fixture_covers_every_required_member(self):
        self.assertEqual(set(self.EXPORTED), fab.ZIP_MEMBERS)

    def test_every_exported_artifact_is_zipped(self):
        self.assertEqual(fab.zip_members(self.EXPORTED + tuple(fab.ZIP_EXCLUDED)),
                         sorted(self.EXPORTED))

    def test_an_extension_outside_the_allowlist_is_rejected(self):
        with self.assertRaisesRegex(
                fab.UploadPackageError, "phantasm-User_2.gm2"):
            fab.zip_members(self.EXPORTED + ("phantasm-User_2.gm2",))

    def test_an_unexported_layer_is_rejected(self):
        for missing in self.EXPORTED:
            extension = os.path.splitext(missing)[1]
            partial = tuple(n for n in self.EXPORTED
                            if os.path.splitext(n)[1] != extension)
            with self.subTest(extension=extension):
                with self.assertRaisesRegex(fab.UploadPackageError, extension):
                    fab.zip_members(partial)

    def test_missing_npth_drill_is_rejected(self):
        partial = tuple(n for n in self.EXPORTED
                        if n != "phantasm-NPTH.drl")
        with self.assertRaisesRegex(fab.UploadPackageError,
                                    "phantasm-NPTH.drl"):
            fab.zip_members(partial)

    def test_the_upload_zip_itself_is_not_a_member(self):
        self.assertNotIn("phantasm-jlc-gerbers.zip",
                         fab.zip_members(self.EXPORTED
                                         + ("phantasm-jlc-gerbers.zip",)))


class PackageManifestTests(unittest.TestCase):
    """The fab run builds a byte-reproducible package; the manifest is what
    lets a rebuild be checked against what was ordered."""

    ARCHIVE = "phantasm-jlc-gerbers.zip"

    def test_manifest_covers_every_member_and_the_archive(self):
        directory = Path(self.enterContext(tempfile.TemporaryDirectory()))
        payloads = {"phantasm-F_Cu.gtl": b"G04 gerber*",
                    "phantasm-PTH.drl": b"M48 M30",
                    self.ARCHIVE: b"zip bytes"}
        for name, data in payloads.items():
            (directory / name).write_bytes(data)
        members = ["phantasm-F_Cu.gtl", "phantasm-PTH.drl"]
        lines = fab.package_manifest(str(directory), members,
                                     self.ARCHIVE).splitlines()
        self.assertEqual(
            lines,
            [f"{hashlib.sha256(payloads[name]).hexdigest()}  {name}"
             for name in members + [self.ARCHIVE]])

    def test_the_manifest_is_not_an_upload_member(self):
        self.assertIn(fab.SUMS_FILE, fab.ZIP_EXCLUDED)
        self.assertNotIn(fab.SUMS_FILE, fab.ZIP_MEMBERS)


class FabContentTests(unittest.TestCase):
    """Zip membership is a filename test; this reads the exported bytes."""

    # An aperture and a flash that uses it, then the same layer plotting nothing.
    GERBER = ("%FSLAX46Y46*%\r\n%MOMM*%\r\n%ADD10C,0.500000*%\r\n"
              "D10*\r\nX1000000Y2000000D03*\r\nM02*\r\n")
    APERTURELESS = "%FSLAX46Y46*%\r\n%MOMM*%\r\nM02*\r\n"
    UNDRAWN = "%FSLAX46Y46*%\r\n%MOMM*%\r\n%ADD10C,0.500000*%\r\nM02*\r\n"
    JOB = '{"GeneralSpecs": {}}\n'
    # One via and one plated pad, one unplated mounting hole, one SMD pad.
    BOARD = ('(kicad_pcb (via (at 1 2) (size 0.45) (drill 0.2)) '
             '(footprint "R" (pad "1" thru_hole circle (at 0 0) (drill 1)) '
             '(pad "2" smd roundrect (at 1 0)) '
             '(pad "" np_thru_hole circle (at 2 0) (drill 2.7))))')

    def drill(self, holes):
        body = "".join(f"X{index}.0Y1.0\n" for index in range(1, holes + 1))
        return "M48\nFMAT,2\nMETRIC\nT1C0.200\n%\nG90\nT1\n" + body + "M30\n"

    def export(self, overrides=None):
        """A staged export of the fixture board, with `overrides` swapped in."""
        directory = Path(self.enterContext(tempfile.TemporaryDirectory()))
        pcb = directory.parent / "content.kicad_pcb"
        pcb.write_text(self.BOARD, encoding="utf-8")
        board = fab.read_board(pcb)
        holes = fab.board_hole_counts(board)
        payloads = {}
        for name in ZipMembershipTests.EXPORTED:
            extension = os.path.splitext(name)[1]
            if extension == ".gbrjob":
                payloads[name] = self.JOB
            elif extension == ".drl":
                payloads[name] = self.drill(holes[fab.DRILL_MEMBERS[name]])
            elif name in fab.APERTURELESS_MEMBERS:
                payloads[name] = self.APERTURELESS
            else:
                payloads[name] = self.GERBER
        payloads.update(overrides or {})
        for name, text in payloads.items():
            (directory / name).write_bytes(text.encode())
        return directory, board

    def validate(self, overrides=None):
        directory, board = self.export(overrides)
        return fab.validate_fab_content(directory, board)

    def test_accepts_an_export_that_carries_the_board(self):
        self.assertEqual(self.validate(), {"plated": 2, "unplated": 1})

    def test_rejects_a_layer_that_plots_nothing(self):
        with self.assertRaisesRegex(fab.FabContentError,
                                    "phantasm-F_Cu.gtl: defines no apertures"):
            self.validate({"phantasm-F_Cu.gtl": self.APERTURELESS})

    def test_rejects_a_layer_that_draws_with_none_of_its_apertures(self):
        with self.assertRaisesRegex(fab.FabContentError, "draws with none"):
            self.validate({"phantasm-In1_Cu.g1": self.UNDRAWN})

    def test_allows_only_the_bottom_stencil_to_be_apertureless(self):
        self.assertEqual(fab.APERTURELESS_MEMBERS, {"phantasm-B_Paste.gbp"})
        self.assertLess(fab.APERTURELESS_MEMBERS, fab.ZIP_MEMBERS)

    def test_rejects_a_drill_file_short_a_hole(self):
        with self.assertRaisesRegex(
                fab.FabContentError,
                r"phantasm-PTH.drl: drills 1 holes, the board carries 2 plated"):
            self.validate({"phantasm-PTH.drl": self.drill(1)})

    def test_rejects_a_drill_file_the_board_does_not_account_for(self):
        with self.assertRaisesRegex(fab.FabContentError,
                                    "no board hole count covers"):
            self.validate({"phantasm-Extra.drl": self.drill(1)})

    def test_counts_a_slot_as_one_hole(self):
        self.assertEqual(
            len(fab.EXCELLON_HOLE.findall("X1.0Y1.0G85X2.0Y1.0\nX3.0Y1.0\n")),
            2)

    def test_committed_board_holes(self):
        # 99 vias and 45 plated pads; the four mounting holes are unplated.
        # Re-measure after promoting a re-route, as with the via floor.
        self.assertEqual(fab.board_hole_counts(fab.read_board(fab.PCB)),
                         {"plated": 144, "unplated": 4})


class PackageVerificationTests(unittest.TestCase):
    """The manifest is only worth writing if something reads it back."""

    def package(self):
        directory = Path(self.enterContext(tempfile.TemporaryDirectory()))
        for name in ZipMembershipTests.EXPORTED + (fab.ARCHIVE,
                                                  "phantasm-BOM.csv"):
            (directory / name).write_bytes(name.encode())
        members = fab.zip_members(os.listdir(directory))
        (directory / fab.SUMS_FILE).write_text(
            fab.package_manifest(str(directory), members, fab.ARCHIVE),
            encoding="utf-8", newline="\n")
        ordered = Path(self.enterContext(tempfile.TemporaryDirectory()))
        baseline = ordered / "ordered-SHA256SUMS.txt"
        baseline.write_bytes((directory / fab.SUMS_FILE).read_bytes())
        return directory, baseline

    def test_accepts_a_package_matching_manifest_and_baseline(self):
        directory, baseline = self.package()
        self.assertEqual(fab.verify_package(str(directory), str(baseline)),
                         len(fab.ZIP_MEMBERS) + 1)

    def test_rejects_an_edited_artifact(self):
        directory, baseline = self.package()
        (directory / "phantasm-F_Cu.gtl").write_bytes(b"edited")
        with self.assertRaisesRegex(fab.PackageVerificationError,
                                    "phantasm-F_Cu.gtl"):
            fab.verify_package(str(directory), str(baseline))

    def test_rejects_a_manifest_that_skips_an_artifact(self):
        directory, baseline = self.package()
        manifest = directory / fab.SUMS_FILE
        kept = [line for line in manifest.read_text(encoding="utf-8").splitlines()
                if not line.endswith("phantasm-PTH.drl")]
        manifest.write_text("\n".join(kept) + "\n", encoding="utf-8",
                            newline="\n")
        with self.assertRaisesRegex(fab.PackageVerificationError,
                                    "records no digest for: phantasm-PTH.drl"):
            fab.verify_package(str(directory), str(baseline))

    def test_rejects_a_package_that_is_not_what_was_ordered(self):
        directory, baseline = self.package()
        lines = baseline.read_text(encoding="utf-8").splitlines()
        flipped = ("1" if lines[0][0] == "0" else "0") + lines[0][1:]
        baseline.write_text("\n".join([flipped] + lines[1:]) + "\n",
                            encoding="utf-8", newline="\n")
        with self.assertRaisesRegex(fab.PackageVerificationError, "ordered"):
            fab.verify_package(str(directory), str(baseline))

    def test_reports_an_absent_baseline(self):
        directory, baseline = self.package()
        baseline.unlink()
        with self.assertRaisesRegex(fab.PackageVerificationError,
                                    "no committed digest baseline"):
            fab.verify_package(str(directory), str(baseline))

    def test_checks_the_package_alone_without_a_baseline(self):
        directory, _ = self.package()
        self.assertEqual(fab.verify_package(str(directory), None),
                         len(fab.ZIP_MEMBERS) + 1)

    def test_rejects_a_malformed_manifest_line(self):
        directory, baseline = self.package()
        (directory / fab.SUMS_FILE).write_text("not a digest line\n",
                                               encoding="utf-8", newline="\n")
        with self.assertRaisesRegex(fab.PackageVerificationError,
                                    "malformed manifest line"):
            fab.verify_package(str(directory), str(baseline))

    def test_reports_an_ungenerated_package(self):
        directory = Path(self.enterContext(tempfile.TemporaryDirectory()))
        with self.assertRaisesRegex(fab.PackageVerificationError,
                                    "no fab package to verify"):
            fab.verify_package(str(directory / "jlc"), None)

    def test_the_committed_baseline_records_the_shipped_package(self):
        if not os.path.exists(fab.SHIPPED_SUMS):
            self.skipTest(f"no committed baseline at {fab.SHIPPED_SUMS}")
        recorded = fab.read_manifest(fab.SHIPPED_SUMS)
        self.assertEqual(set(recorded), fab.ZIP_MEMBERS | {fab.ARCHIVE})


class NetlistSpecTests(unittest.TestCase):
    """The fab run holds its own netlist export to check.py's named-net table."""

    @staticmethod
    def spec_nodes():
        """EXPECT -> {net: [(ref, pin)]}; a symmetric key names the ref alone."""
        return {name: [tuple(key.split(".")) if "." in key else (key, "1")
                       for key in sorted(keys)]
                for name, keys in fab.netlist_spec.EXPECT.items()}

    def validate(self, nets):
        blocks = "".join(
            '(net (name "{}") {})'.format(
                name, " ".join(f'(node (ref "{ref}") (pin "{pin}"))'
                               for ref, pin in nodes))
            for name, nodes in nets.items())
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "test.net"
            path.write_text(f"(export (nets {blocks}))", encoding="utf-8")
            with contextlib.redirect_stdout(io.StringIO()):
                return fab.validate_netlist_spec(path)

    def test_accepts_the_specified_partition(self):
        self.assertEqual(self.validate(self.spec_nodes()),
                         len(fab.netlist_spec.EXPECT))

    def test_rejects_a_net_missing_a_member(self):
        nets = self.spec_nodes()
        nets["SYNC_BUS"].pop()
        with self.assertRaisesRegex(fab.NetlistSpecError,
                                    "does not match the electrical"):
            self.validate(nets)

    def test_rejects_a_permuted_pinout(self):
        nets = self.spec_nodes()
        nets["DATA_IN"] = [("U_MCU", "11"), ("U1", "5")]
        with self.assertRaisesRegex(fab.NetlistSpecError,
                                    "does not match the electrical"):
            self.validate(nets)


if __name__ == "__main__":
    unittest.main()
