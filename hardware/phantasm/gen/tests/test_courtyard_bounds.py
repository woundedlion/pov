"""Courtyard-bounds helper tests."""

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import sexp
from courtyard_bounds import courtyard_box


class CourtyardBoundsTests(unittest.TestCase):
    def test_rectangle_corners_are_rotated_and_translated(self):
        footprint = sexp.parse_one(
            '(footprint (at 10 20 90) (fp_rect (start -1 -2) (end 3 4) '
            '(layer "F.CrtYd")))')
        self.assertEqual(courtyard_box(footprint), (8.0, 17.0, 14.0, 21.0))

    def test_circle_uses_its_full_radius(self):
        footprint = sexp.parse_one(
            '(footprint (at 10 -2) (fp_circle (center 2 3) (end 5 7) '
            '(layer "B.CrtYd")))')
        self.assertEqual(courtyard_box(footprint), (7.0, -4.0, 17.0, 6.0))

    def test_polygon_includes_all_vertices(self):
        footprint = sexp.parse_one(
            '(footprint (at 5 7) (fp_poly (pts (xy -1 0) (xy 2 1) (xy 0 4)) '
            '(layer "F.CrtYd")))')
        self.assertEqual(courtyard_box(footprint), (4.0, 7.0, 7.0, 11.0))

    def test_other_layers_do_not_expand_the_courtyard(self):
        footprint = sexp.parse_one(
            '(footprint (at 0 0) (fp_rect (start 0 0) (end 1 2) '
            '(layer "F.CrtYd")) (fp_rect (start -100 -100) (end 100 100) '
            '(layer "F.Fab")))')
        self.assertEqual(courtyard_box(footprint), (0.0, 0.0, 1.0, 2.0))

    def test_footprint_without_a_courtyard_has_no_bounds(self):
        footprint = sexp.parse_one(
            '(footprint (at 0 0) (fp_rect (start -1 -1) (end 1 1) '
            '(layer "F.Fab")))')
        self.assertIsNone(courtyard_box(footprint))

    def test_arc_includes_cardinal_extrema(self):
        footprint = sexp.parse(
            '(footprint (at 0 0) (fp_arc (start 0.8 0.6) (mid -0.6 0.8) '
            '(end -0.8 -0.6) (layer "F.CrtYd")))')[0]
        for actual, expected in zip(courtyard_box(footprint), (-1.0, -0.6, 0.8, 1.0)):
            self.assertAlmostEqual(actual, expected)
