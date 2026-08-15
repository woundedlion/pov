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


class RotatedBoundsTests(unittest.TestCase):
    BOX = (-1.0, -3.0, 2.0, 4.0)

    def test_right_angles_swap_the_extent(self):
        self.assertEqual(pcb._rot_bb(self.BOX, 0), self.BOX)
        self.assertEqual(pcb._rot_bb(self.BOX, 90), (-3.0, -2.0, 4.0, 1.0))
        self.assertEqual(pcb._rot_bb(self.BOX, 180), (-2.0, -4.0, 1.0, 3.0))
        self.assertEqual(pcb._rot_bb(self.BOX, 270), (-4.0, -1.0, 3.0, 2.0))

    def test_rejects_an_angle_it_cannot_rotate(self):
        for rot in (45, "90", sexp.Sym("90")):
            with self.subTest(rot=rot):
                with self.assertRaisesRegex(ValueError, "0/90/180/270"):
                    pcb._rot_bb(self.BOX, rot)


class MountingKeepoutTests(unittest.TestCase):
    """Nothing may be packed onto a mounting hole: the screw head sits there."""

    BOXES = {
        "J1": (-1.8, -3.0, 1.8, 3.0),
        "J4": (-1.8, -2.0, 1.8, 2.0),
        "J2": (-1.8, -4.3, 1.8, 4.3),
        "J3A": (-1.8, -4.3, 1.8, 4.3),
        "J3B": (-1.8, -4.3, 1.8, 4.3),
        "U_MCU": (-24.2, -9.5, 24.2, 9.5),
        "C_IN": (-1.5, -0.8, 1.5, 0.8),
        "R1": (-1.5, -0.8, 1.5, 0.8),
    }

    def test_pack_places_every_part_clear_of_every_hole(self):
        place, length = pcb.pack(dict(self.BOXES), pcb.PCB_W)
        self.assertEqual(sorted(place), sorted(self.BOXES))
        self.assertEqual(pcb.keepout_clashes(place, self.BOXES, length), [])

    def test_packed_length_reserves_a_tail_for_the_far_holes(self):
        place, length = pcb.pack(dict(self.BOXES), pcb.PCB_W)
        extent = max(x + pcb._rot_bb(self.BOXES[ref], rot)[2]
                     for ref, (x, _, rot) in place.items())
        self.assertGreaterEqual(length - extent,
                                pcb.MOUNTING_HOLE_INSET + pcb.MOUNTING_KEEPOUT_RADIUS)

    def test_pack_keeps_every_part_inside_the_outline(self):
        place, length = pcb.pack(dict(self.BOXES), pcb.PCB_W)
        self.assertEqual(
            pcb.outline_overflows(place, self.BOXES, length, list(self.BOXES)), [])

    def test_a_part_too_tall_to_pack_overflows_the_outline(self):
        boxes = dict(self.BOXES, TALL=(-1.5, -pcb.PCB_W, 1.5, pcb.PCB_W))
        place, length = pcb.pack(boxes, pcb.PCB_W)
        self.assertEqual(
            [entry.split()[0]
             for entry in pcb.outline_overflows(place, boxes, length, list(boxes))],
            ["TALL"])

    def test_clash_report_names_the_hole(self):
        place = {"J1": (pcb.MOUNTING_HOLE_INSET, pcb.MOUNTING_HOLE_INSET, 0)}
        self.assertEqual(pcb.keepout_clashes(place, self.BOXES, 40.0), ["J1/H1"])

    def test_a_part_flush_with_the_keepout_edge_is_not_a_clash(self):
        edge = pcb.MOUNTING_HOLE_INSET + pcb.MOUNTING_KEEPOUT_RADIUS
        place = {"C_IN": (edge + 1.5, pcb.MOUNTING_HOLE_INSET, 0)}
        self.assertEqual(pcb.keepout_clashes(place, self.BOXES, 40.0), [])


if __name__ == "__main__":
    unittest.main()
