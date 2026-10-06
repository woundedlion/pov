"""Courtyard-bounds helper tests."""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import sexp  # noqa: E402
from courtyard_bounds import courtyard_box  # noqa: E402


class CourtyardBoundsTests(unittest.TestCase):
    def test_arc_includes_cardinal_extrema(self):
        footprint = sexp.parse(
            '(footprint (at 0 0) (fp_arc (start 0.8 0.6) (mid -0.6 0.8) '
            '(end -0.8 -0.6) (layer "F.CrtYd")))')[0]
        for actual, expected in zip(courtyard_box(footprint), (-1.0, -0.6, 0.8, 1.0)):
            self.assertAlmostEqual(actual, expected)
