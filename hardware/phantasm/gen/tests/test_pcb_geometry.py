import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pcb
import sexp


class FootprintBoundsTests(unittest.TestCase):
    def test_circle_uses_full_radius(self):
        footprint = sexp.parse("(footprint (fp_circle (center 2 3) (end 5 7)))")[0]
        self.assertEqual(pcb.fp_bbox(footprint), (-3.0, -2.0, 7.0, 8.0))

    def test_arc_includes_cardinal_extrema(self):
        footprint = sexp.parse(
            "(footprint (fp_arc (start 1 0) (mid 0 1) (end -1 0)))")[0]
        bounds = pcb.fp_bbox(footprint)
        for actual, expected in zip(bounds, (-1.0, 0.0, 1.0, 1.0)):
            self.assertAlmostEqual(actual, expected)


if __name__ == "__main__":
    unittest.main()
