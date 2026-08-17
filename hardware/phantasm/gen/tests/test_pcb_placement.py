import sys
import unittest
from pathlib import Path

GEN_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN_DIR))

import pcb  # noqa: E402
import sexp  # noqa: E402
from kicad_common import F  # noqa: E402

ROUTED = GEN_DIR.parent / pcb.PCB_FILE
UNPLACED = GEN_DIR.parent / pcb.UNPLACED_FILE


def placements(path):
    """ref -> (x, y, rotation) for every referenced footprint on a board."""
    root = sexp.parse(path.read_text(encoding="utf-8"))[0]
    out = {}
    for node in F(root, "footprint"):
        at = sexp.val(node, "at")
        for prop in F(node, "property"):
            if prop[1] == "Reference":
                out[prop[2]] = (float(at[0]), float(at[1]),
                                float(at[2]) if len(at) > 2 else 0.0)
    return out


class LockedPlacementTests(unittest.TestCase):
    """QUILTER_FIXED is what `pcb.py --unplaced` locks before an autoplacer run,
    so a constant that disagrees with the committed boards hands Quilter a
    different mechanical placement than the routing was produced under."""

    def test_every_locked_placement_matches_the_committed_boards(self):
        routed = placements(ROUTED)
        unplaced = placements(UNPLACED)
        self.assertTrue(pcb.QUILTER_FIXED)
        for ref, (x, y, rot) in sorted(pcb.QUILTER_FIXED.items()):
            with self.subTest(ref=ref):
                self.assertEqual(routed.get(ref), (float(x), float(y), float(rot)))
                self.assertEqual(unplaced.get(ref), (float(x), float(y), float(rot)))


class OutlineBoundsTests(unittest.TestCase):
    """Both board paths run this gate: DRC never reports copper hanging past
    Edge.Cuts, and the fab routs it away."""

    PADS = {"R1": (-2.0, -0.5, 2.0, 0.5)}
    LENGTH = 20.0

    def overflows(self, x, y, rot=0, refs=("R1",)):
        return pcb.outline_overflows({"R1": (x, y, rot)}, self.PADS,
                                     self.LENGTH, refs)

    def test_accepts_a_part_inside_the_outline(self):
        self.assertEqual(self.overflows(10.0, 10.0), [])

    def test_reports_a_part_past_the_far_end(self):
        self.assertEqual(self.overflows(self.LENGTH, 10.0),
                         ["R1 [18.00,9.50]-[22.00,10.50]"])

    def test_reports_a_part_past_the_hub_end(self):
        self.assertEqual(len(self.overflows(1.0, 10.0)), 1)

    def test_reports_a_part_past_the_board_width(self):
        self.assertEqual(len(self.overflows(10.0, pcb.PCB_W)), 1)

    def test_measures_a_rotated_part_rotated(self):
        self.assertEqual(self.overflows(10.0, 1.0), [])
        self.assertEqual(len(self.overflows(10.0, 1.0, rot=90)), 1)

    def test_ignores_refs_outside_the_checked_set(self):
        self.assertEqual(self.overflows(self.LENGTH, 10.0, refs=()), [])


class OverwriteAuthorizationTests(unittest.TestCase):
    """The board overwrite and the Teensy footprint library overwrite are
    unrelated destructive writes and take separate opt-ins."""

    def test_board_force_does_not_authorize_the_footprint_library(self):
        args = pcb.parse_args(["--unplaced", "--force"])
        self.assertTrue(args.force)
        self.assertFalse(args.force_teensy_library)

    def test_footprint_library_force_does_not_authorize_the_board(self):
        args = pcb.parse_args(["--force-teensy-library"])
        self.assertFalse(args.force)
        self.assertTrue(args.force_teensy_library)


if __name__ == "__main__":
    unittest.main()
