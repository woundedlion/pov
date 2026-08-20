import sys
from decimal import Decimal
from pathlib import Path
import unittest


GEN_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = GEN_DIR.parents[2]
sys.path.insert(0, str(GEN_DIR))

import board_metadata   # noqa: E402
import fab   # noqa: E402


class BoardMetadataTests(unittest.TestCase):
    def test_extracts_committed_routed_board(self):
        metadata = board_metadata.load_board(
            REPO_ROOT / "hardware" / "phantasm" / "phantasm.kicad_pcb"
        )

        self.assertEqual(
            (metadata.width_mm, metadata.height_mm),
            (Decimal("58.28"), Decimal("32")),
        )
        self.assertEqual(metadata.footprint_sides, (("F.Cu", 32), ("B.Cu", 0)))
        self.assertEqual(metadata.track_segments, 339)
        self.assertEqual(metadata.vias, 100)
        self.assertEqual(metadata.copper_pours, 2)
        self.assertEqual(metadata.pour_layers, (("In1.Cu", 1), ("In2.Cu", 1)))
        self.assertEqual(metadata.rule_areas, 4)
        self.assertEqual(
            metadata.rule_area_layers,
            (("F.Cu", 4), ("In1.Cu", 4), ("In2.Cu", 4), ("B.Cu", 4)),
        )
        self.assertEqual(metadata.copper_layers, ("F.Cu", "In1.Cu", "In2.Cu", "B.Cu"))
        self.assertEqual(
            metadata.copper_stackup,
            (
                ("F.Cu", Decimal("0.035001")),
                ("In1.Cu", Decimal("0.015189")),
                ("In2.Cu", Decimal("0.015189")),
                ("B.Cu", Decimal("0.035001")),
            ),
        )
        self.assertEqual(metadata.copper_finish, "Lead-Free")

    def test_pour_count_matches_the_fab_gate(self):
        """The facts block and gen/fab.py must count copper the same way."""
        board = REPO_ROOT / "hardware" / "phantasm" / "phantasm.kicad_pcb"
        metadata = board_metadata.load_board(board)

        self.assertEqual(metadata.copper_pours,
                         fab.validate_zone_geometry(board))

    def test_counts_a_keepout_rule_area_separately_from_a_pour(self):
        metadata = board_metadata.parse_board("""
            (kicad_pcb
              (gr_rect (start 0 0) (end 10 10) (layer "Edge.Cuts"))
              (general (thickness 1.6))
              (layers (0 "F.Cu" signal) (2 "B.Cu" signal))
              (setup (stackup
                (layer "F.Cu" (type "copper") (thickness 0.035))
                (layer "B.Cu" (type "copper") (thickness 0.035))
                (copper_finish "ENIG")))
              (zone (net 1) (layer "F.Cu") (min_thickness 0.25))
              (zone (net 0) (layers "F.Cu" "B.Cu")
                (keepout (tracks not_allowed) (copperpour not_allowed))))
        """)

        self.assertEqual(metadata.copper_pours, 1)
        self.assertEqual(metadata.pour_layers, (("F.Cu", 1),))
        self.assertEqual(metadata.rule_areas, 1)
        self.assertEqual(metadata.rule_area_layers, (("F.Cu", 1), ("B.Cu", 1)))

    def test_rejects_malformed_board(self):
        with self.assertRaisesRegex(board_metadata.MetadataError, "invalid KiCad S-expression"):
            board_metadata.parse_board("(kicad_pcb (general")

    def test_curved_outline_bounds_include_arc_extrema(self):
        root = board_metadata.sexp.parse("""
            (kicad_pcb
              (gr_arc (start 0 0) (mid 5 -5) (end 10 0) (layer "Edge.Cuts"))
              (gr_line (start 10 0) (end 10 5) (layer "Edge.Cuts"))
              (gr_line (start 10 5) (end 0 5) (layer "Edge.Cuts"))
              (gr_line (start 0 5) (end 0 0) (layer "Edge.Cuts")))
        """)[0]
        self.assertEqual(
            board_metadata._outline_bounds(root),
            (Decimal("10"), Decimal("10")),
        )

    def test_circle_outline_uses_full_radius(self):
        root = board_metadata.sexp.parse("""
            (kicad_pcb
              (gr_circle (center 1 2) (end 4 6) (layer "Edge.Cuts")))
        """)[0]
        self.assertEqual(
            board_metadata._outline_bounds(root),
            (Decimal("10"), Decimal("10")),
        )

    def test_rejects_footprint_level_edge_cuts_geometry(self):
        """Board edges drawn as footprint graphics would under-report bounds."""
        root = board_metadata.sexp.parse("""
            (kicad_pcb
              (gr_rect (start 0 0) (end 4 4) (layer "Edge.Cuts"))
              (footprint "Board:Outline"
                (layer "F.Cu")
                (fp_line (start 0 0) (end 60 0) (layer "Edge.Cuts"))))
        """)[0]
        with self.assertRaisesRegex(board_metadata.MetadataError,
                                    "Edge.Cuts graphics inside footprint: Board:Outline"):
            board_metadata._outline_bounds(root)

    def test_footprint_graphics_off_edge_cuts_do_not_reject(self):
        root = board_metadata.sexp.parse("""
            (kicad_pcb
              (gr_rect (start 0 0) (end 10 10) (layer "Edge.Cuts"))
              (footprint "Lib:R"
                (layer "F.Cu")
                (fp_line (start 0 0) (end 1 0) (layer "F.SilkS"))
                (fp_poly (pts (xy 0 0) (xy 1 1)) (layer "F.Fab"))))
        """)[0]
        self.assertEqual(
            board_metadata._outline_bounds(root),
            (Decimal("10"), Decimal("10")),
        )

    def test_rejects_readme_drift(self):
        metadata = board_metadata.load_board(
            REPO_ROOT / "hardware" / "phantasm" / "phantasm.kicad_pcb"
        )
        facts = board_metadata.render_facts(metadata)
        readme = (REPO_ROOT / "hardware" / "phantasm" / "README.md").read_text(encoding="utf-8")
        current = f"| Track segments | {metadata.track_segments} |"
        stale = readme.replace(
            current,
            f"| Track segments | {metadata.track_segments + 1} |",
            1,
        )

        with self.assertRaisesRegex(board_metadata.MetadataError, "facts are stale"):
            board_metadata.check_facts(stale, facts)


if __name__ == "__main__":
    unittest.main()
