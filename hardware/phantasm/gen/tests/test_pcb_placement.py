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


if __name__ == "__main__":
    unittest.main()
