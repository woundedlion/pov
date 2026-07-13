import sys
from decimal import Decimal
from pathlib import Path
import unittest


GEN_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = GEN_DIR.parents[2]
sys.path.insert(0, str(GEN_DIR))

import board_metadata


class BoardMetadataTests(unittest.TestCase):
    def test_extracts_committed_routed_board(self):
        metadata = board_metadata.load_board(
            REPO_ROOT / "hardware" / "phantasm" / "phantasm.kicad_pcb"
        )

        self.assertEqual(
            (metadata.width_mm, metadata.height_mm),
            (Decimal("58.28"), Decimal("32")),
        )
        self.assertEqual(metadata.footprint_sides, (("F.Cu", 29), ("B.Cu", 0)))
        self.assertEqual(metadata.track_segments, 384)
        self.assertEqual(metadata.vias, 83)
        self.assertEqual(metadata.zones, 2)
        self.assertEqual(metadata.zone_layers, (("In1.Cu", 1), ("In2.Cu", 1)))
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

    def test_rejects_malformed_board(self):
        with self.assertRaisesRegex(board_metadata.MetadataError, "invalid KiCad S-expression"):
            board_metadata.parse_board("(kicad_pcb (general")

    def test_rejects_readme_drift(self):
        metadata = board_metadata.load_board(
            REPO_ROOT / "hardware" / "phantasm" / "phantasm.kicad_pcb"
        )
        facts = board_metadata.render_facts(metadata)
        readme = (REPO_ROOT / "hardware" / "phantasm" / "README.md").read_text(encoding="utf-8")
        stale = readme.replace("| Track segments | 384 |", "| Track segments | 385 |", 1)

        with self.assertRaisesRegex(board_metadata.MetadataError, "facts are stale"):
            board_metadata.check_facts(stale, facts)


if __name__ == "__main__":
    unittest.main()
