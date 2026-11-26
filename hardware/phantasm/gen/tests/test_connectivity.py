import contextlib
import io
import math
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path

GEN = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN))

import connectivity  # noqa: E402
import pcb  # noqa: E402
import sexp  # noqa: E402

from committed_boards import committed_board  # noqa: E402

ROUTED = committed_board()

# One F.Cu track joining the two SMD pads of one net.
BOARD = """(kicad_pcb
\t(layers
\t\t(0 "F.Cu" signal)
\t\t(4 "In1.Cu" signal)
\t\t(2 "B.Cu" signal)
\t\t(25 "Edge.Cuts" user)
\t)
\t(footprint "R"
\t\t(layer "F.Cu")
\t\t(at 0 0)
\t\t(property "Reference" "R1")
\t\t(pad "1" smd rect (at 0 0) (size 0.5 0.5) (layers "F.Cu") (net "/DATA"))
\t\t(pad "2" smd rect (at 10 0) (size 0.5 0.5) (layers "F.Cu") (net "/DATA"))
\t)
\t(segment (start 0 0) (end 10 0) (width 0.2) (layer "F.Cu") (net "/DATA"))
)"""

NATIVE_BOARD = (BOARD
                .replace("\t(layers", '\t(net 1 "/DATA")\n\t(layers')
                .replace('(net "/DATA"))', '(net 1 "/DATA"))', 2)
                .replace('(net "/DATA"))', '(net 1))'))

# The same net routed as a curved track, which Quilter emits for corners.
ARC_BOARD = BOARD.replace(
    '(segment (start 0 0) (end 10 0) (width 0.2) (layer "F.Cu") (net "/DATA"))',
    '(arc (start 0 0) (mid 5 3) (end 10 0) (width 0.2) (layer "F.Cu")'
    ' (net "/DATA"))')

ROUND_PAD_BOARD = """(kicad_pcb
\t(layers (0 "F.Cu" signal))
\t(footprint "J"
\t\t(layer "F.Cu")
\t\t(at 0 0)
\t\t(property "Reference" "J1")
\t\t(pad "1" smd circle (at 0 0) (size 2.7 2.7) (layers "F.Cu") (net "/ROUND"))
\t\t(pad "2" smd circle (at 3.818 0) (size 2.7 2.7) (layers "F.Cu") (net "/ROUND"))
\t)
)"""

# A chip resistor placed at 90 degrees: KiCad folds the footprint rotation into
# the pad angle, so the land's long axis runs across the (size) x axis.
ROTATED_LAND_BOARD = """(kicad_pcb
	(layers (0 "F.Cu" signal))
	(footprint "R"
		(layer "F.Cu")
		(at 29.7 22.8 90)
		(property "Reference" "R2")
		(pad "1" smd roundrect (at -0.825 0 90) (size 1.2 0.95) (layers "F.Cu") (net "/R"))
	)
)"""

# Two SOIC-14 lands on one net, with the track stopping 0.535 mm short of the
# second: bare laminate, yet inside the disc that circumscribes a 1.95 x 0.6
# land.
SOIC_BOARD = """(kicad_pcb
	(layers (0 "F.Cu" signal))
	(footprint "SOIC"
		(layer "F.Cu")
		(at 0 0)
		(property "Reference" "U1")
		(pad "1" smd roundrect (at 0 0) (size 1.95 0.6) (layers "F.Cu") (net "/SIG"))
		(pad "11" smd roundrect (at 0 -3) (size 1.95 0.6) (layers "F.Cu") (net "/SIG"))
	)
	(segment (start 0 0) (end 0 -2.065) (width 0.2) (layer "F.Cu") (net "/SIG"))
)"""

# The solder-jumper land: a custom pad whose primitives reach past (size),
# with KiCad's own nesting under (primitives ...).
CUSTOM_PAD_BOARD = """(kicad_pcb
	(layers (0 "F.Cu" signal))
	(footprint "SJ"
		(layer "F.Cu")
		(at 0 0)
		(property "Reference" "SJ1")
		(pad "1" smd custom
			(at 0 0 90)
			(size 1 0.5)
			(layers "F.Cu" "F.Mask")
			(net "/SJ")
			(primitives
				(gr_circle (center 0 0.25) (end 0.5 0.25) (width 0) (fill yes))
				(gr_circle (center 0 -0.25) (end 0.5 -0.25) (width 0) (fill yes))
				(gr_poly
					(pts (xy 0.5 0.75) (xy 0 0.75) (xy 0 -0.75) (xy 0.5 -0.75))
					(width 0)
					(fill yes)
				)
			)
		)
	)
)"""

# A through-hole pad on GND whose only copper is the In1 pour it sits in, with
# a solid pour spanning the through-hole pads.
POUR_BOARD = """(kicad_pcb
\t(layers
\t\t(0 "F.Cu" signal)
\t\t(4 "In1.Cu" signal)
\t\t(2 "B.Cu" signal)
\t)
\t(footprint "J"
\t\t(layer "F.Cu")
\t\t(at 5 5)
\t\t(property "Reference" "J1")
\t\t(pad "1" thru_hole circle (at 0 0) (size 1.6 1.6) (layers "*.Cu") (net "/GND"))
\t\t(pad "2" thru_hole circle (at 20 0) (size 1.6 1.6) (layers "*.Cu") (net "/GND"))
\t)
\t(zone
\t\t(net "/GND")
\t\t(layer "In1.Cu")
\t\t(filled_polygon
\t\t\t(layer "In1.Cu")
\t\t\t(pts (xy 0 0) (xy 30 0) (xy 30 10) (xy 0 10))
\t\t)
\t)
)"""

# A pour spanning two layers: KiCad spells that `(layers ...)`, and the fills of
# an unfilled-then-refilled zone can carry no layer of their own.
SPAN_POUR_BOARD = POUR_BOARD.replace(
    '\t(layer "In1.Cu")\n',
    '\t(layers "In1.Cu" "B.Cu")\n', 1).replace(
    '\t\t\t(layer "In1.Cu")\n', "", 1)

STAR_POUR_BOARD = POUR_BOARD.replace(
    '\t(layer "In1.Cu")\n', '\t(layers "*.Cu")\n', 1).replace(
    '\t\t\t(layer "In1.Cu")\n', "", 1)

EMPTY_BOARD = """(kicad_pcb
	(layers
		(0 "F.Cu" signal)
		(2 "B.Cu" signal)
	)
)"""


def parse(text):
    return sexp.parse(text)[0]


def drop(root, kind, net):
    return [node for node in root
            if not (isinstance(node, list) and node and node[0] == kind
                    and connectivity.net_name(node) == net)]


def probe(point, layer="F.Cu"):
    """A point-sized copper contact on one layer."""
    return connectivity.Capsule(point, point, 0, {layer})


class ExactLandTests(unittest.TestCase):
    def test_oval_uses_its_minor_radius_and_rotated_major_axis(self):
        for size in ("4 1", "1 4"):
            pad = sexp.parse_one(
                f'(pad "1" smd oval (at 2 1 90) (size {size}) (layers "F.Cu"))')
            land = connectivity.pad_copper(pad, (10, 20), 90, ["F.Cu"])
            inside = (11, 16) if size == "4 1" else (13, 18)
            outside = (11.45, 16.1) if size == "4 1" else (12.9, 18.45)
            with self.subTest(size=size):
                self.assertTrue(land.touches(probe(inside)))
                self.assertFalse(land.touches(probe(outside)))

    def test_roundrect_corner_contacts_use_the_minkowski_radius(self):
        for angle, radius_field in ((0, ""), (0, "(roundrect_rratio 0.25)"),
                                    (37, "(roundrect_rratio 0.25)"),
                                    (90, "(roundrect_rratio 0.25)")):
            pad = sexp.parse_one(
                f'(pad "1" smd roundrect (at 0 0 {angle}) (size 4 2) '
                f'{radius_field} (layers "F.Cu"))')
            land = connectivity.pad_copper(pad, (0, 0), 0, ["F.Cu"])
            for point, touches in (((1.8, 0.9), True), ((1.95, 0.95), False)):
                x, y = point
                radians = math.radians(angle)
                position = (x * math.cos(radians) + y * math.sin(radians),
                            y * math.cos(radians) - x * math.sin(radians))
                for contact in (probe(position), connectivity.Polygon([position] * 4, {"F.Cu"})):
                    with self.subTest(angle=angle, point=point, shape=type(contact).__name__):
                        self.assertEqual(land.touches(contact), touches)
                        self.assertEqual(contact.touches(land), touches)

    def test_unfilled_custom_circle_preserves_its_hole(self):
        for fill, inside in (("none", False), ("yes", True)):
            pad = sexp.parse_one(
                '(pad "1" smd custom (at 0 0) (size 0.1 0.1) (layers "F.Cu") '
                '(primitives (gr_circle (center 3 0) (end 4 0) '
                f'(width 0.2) (fill {fill}))))')
            land = connectivity.pad_copper(pad, (0, 0), 0, ["F.Cu"])
            self.assertEqual(land.touches(probe((3, 0))), inside)
            self.assertTrue(land.touches(probe((4.09, 0))))
            self.assertFalse(land.touches(probe((4.12, 0))))

    def test_custom_polygon_retains_concavity_and_stroke(self):
        pad = sexp.parse_one(
            '(pad "1" smd custom (at 0 0) (size 0.1 0.1) (layers "F.Cu") '
            '(primitives (gr_poly (pts (xy 0 0) (xy 3 0) (xy 3 1) '
            '(xy 1 1) (xy 1 3) (xy 0 3)) (width 0.2) (fill yes))))')
        land = connectivity.pad_copper(pad, (0, 0), 0, ["F.Cu"])
        self.assertTrue(land.touches(probe((0.5, 2))))
        self.assertTrue(land.touches(probe((1.09, 2))))
        self.assertFalse(land.touches(probe((1.12, 2))))
        self.assertFalse(land.touches(probe((2, 2))))

    def test_custom_anchor_shape_controls_corner_copper(self):
        for anchor, expected in (("circle", False), ("rect", True)):
            pad = sexp.parse_one(
                '(pad "1" smd custom (at 0 0) (size 2 2) (layers "F.Cu") '
                f'(options (anchor {anchor})))')
            land = connectivity.pad_copper(pad, (0, 0), 0, ["F.Cu"])
            self.assertEqual(land.touches(probe((0.9, 0.9))), expected)

    def test_unsupported_land_geometry_is_rejected(self):
        for shape, extra in (("unsupported", ""),
                             ("custom", "(primitives (gr_curve))")):
            pad = sexp.parse_one(
                f'(pad "1" smd {shape} (at 0 0) (size 1 1) (layers "F.Cu") {extra})')
            with self.subTest(shape=shape), self.assertRaisesRegex(ValueError, "unsupported"):
                connectivity.pad_copper(pad, (0, 0), 0, ["F.Cu"])

    def test_arc_contacts_follow_major_and_reversed_sweeps(self):
        on_arc = probe((math.sqrt(50), math.sqrt(50)))
        on_chord = probe((5, 5))
        for start, mid, end in (((10, 0), (0, 10), (-10, 0)),
                                ((-10, 0), (0, 10), (10, 0)),
                                ((10, 0), (-10, 0), (0, -10))):
            arc = sexp.parse_one(
                '(arc ' + ' '.join(f'({key} {p[0]} {p[1]})'
                                  for key, p in zip(("start", "mid", "end"), (start, mid, end)))
                + ' (width 0.02) (layer "F.Cu"))')
            copper = connectivity.CopperUnion(connectivity.arc_capsules(arc), {"F.Cu"})
            with self.subTest(start=start, end=end):
                self.assertTrue(copper.touches(on_arc))
                self.assertFalse(copper.touches(on_chord))


class CustomPrimitiveBoundsTests(unittest.TestCase):
    def test_custom_strokes_preserve_empty_space(self):
        for primitive, on_stroke, empty in (
                ("(gr_line (start -2 -1) (end 3 4) (width 0.2))", (-1, 0), (-1, 3)),
                ("(gr_rect (start -2 -1) (end 3 4) (width 0.2))", (-2, 2), (1, 1)),
                ("(gr_arc (start 2 0) (mid 0 2) (end -2 0) (width 0.2))",
                 (math.sqrt(2), math.sqrt(2)), (1, 0.5))):
            with self.subTest(primitive=primitive):
                pad = sexp.parse_one(
                    '(pad "1" smd custom (at 0 0) (size 0.1 0.1) '
                    '(layers "F.Cu") (primitives ' + primitive + '))')
                copper = connectivity.pad_copper(pad, (0, 0), 0, ["F.Cu"])
                self.assertTrue(copper.touches(probe(on_stroke)))
                self.assertFalse(copper.touches(probe(empty)))


class SyntheticBoardTests(unittest.TestCase):
    def test_circle_pad_uses_its_declared_radius(self):
        root = parse(ROUND_PAD_BOARD)
        pads = list(connectivity.F(
            connectivity.F(root, "footprint")[0], "pad"))
        first = connectivity.pad_copper(pads[0], (0, 0), 0, ["F.Cu"])
        self.assertAlmostEqual(first.radius, 1.35)

    def test_a_rectangular_land_is_the_rotated_rectangle(self):
        root = parse(ROTATED_LAND_BOARD.replace("roundrect", "rect"))
        pad = connectivity.F(
            connectivity.F(root, "footprint")[0], "pad")[0]
        land = connectivity.pad_copper(pad, (29.7, 22.8), 90, ["F.Cu"])
        xs = sorted(round(x, 3) for x, _ in land.polygon)
        ys = sorted(round(y, 3) for _, y in land.polygon)
        self.assertEqual((xs[0], xs[-1]), (29.225, 30.175))
        self.assertEqual((ys[0], ys[-1]), (23.025, 24.225))

    def test_trapezoid_taper_excludes_the_bounding_box_corners(self):
        for delta, inside, outside in (
                ("0.8 0", (0.9, 0.5), (0.9, 1.2)),
                ("-0.8 0", (-0.9, 0.5), (-0.9, 1.2)),
                ("0 0.8", (0.5, 0.9), (1.2, -0.9)),
                ("0 -0.8", (-0.5, 0.9), (-1.2, 0.9))):
            with self.subTest(delta=delta):
                pad = sexp.parse_one(
                    '(pad "1" smd trapezoid (at 0 0) (size 2 2) '
                    f'(rect_delta {delta}) (layers "F.Cu"))')
                land = connectivity.pad_copper(pad, (0, 0), 0, ["F.Cu"])
                self.assertTrue(land.touches(probe(inside)))
                self.assertFalse(land.touches(probe(outside)))

    def test_custom_pad_rotates_its_primitive_union(self):
        pad = connectivity.F(
            connectivity.F(parse(CUSTOM_PAD_BOARD), "footprint")[0], "pad")[0]
        land = connectivity.pad_copper(pad, (0, 0), 0, ["F.Cu"])
        self.assertFalse(land.touches(probe((0.70, 0.45))))
        self.assertTrue(land.touches(probe((0.65, -0.4))))
        short_track = connectivity.Capsule((0, 2), (0, 0.775), 0.1, {"F.Cu"})
        touching_track = connectivity.Capsule((0, 2), (0, 0.6), 0.1, {"F.Cu"})
        self.assertFalse(land.touches(short_track))
        self.assertTrue(land.touches(touching_track))

    def test_rotated_square_pad_reaches_the_outline(self):
        pad = sexp.parse_one('(pad "1" smd rect (at 0 0 45) (size 2 2) (layers "F.Cu"))')
        bounds = {"R1": pcb.fp_bbox(["footprint", pad], pads_only=True)}
        self.assertEqual(pcb.outline_overflows({"R1": (2, 2, 0)}, bounds, 10, ["R1"]), [])
        self.assertEqual(len(pcb.outline_overflows({"R1": (1.1, 2, 0)}, bounds, 10, ["R1"])), 1)

    def test_custom_pad_primitive_bounds_reach_the_outline_gate(self):
        footprint = connectivity.F(parse(CUSTOM_PAD_BOARD), "footprint")[0]
        bounds = pcb.fp_bbox(footprint, pads_only=True)
        for actual, expected in zip(bounds, (-0.75, -0.5, 0.75, 0.5)):
            self.assertAlmostEqual(actual, expected)
        self.assertEqual(pcb.outline_overflows(
            {"SJ1": (9.25, 1, 0)}, {"SJ1": bounds}, 10, ["SJ1"]), [])
        self.assertEqual(len(pcb.outline_overflows(
            {"SJ1": (9.4, 1, 0)}, {"SJ1": bounds}, 10, ["SJ1"])), 1)

    def test_fill_overlap_is_symmetric(self):
        big = connectivity.Polygon(
            [(0, 0), (10, 0), (10, 10), (0, 10)], {"In1.Cu"})
        small = connectivity.Polygon(
            [(4, 4), (6, 4), (6, 6), (4, 6)], {"In1.Cu"})
        self.assertTrue(big.touches(small))
        self.assertTrue(small.touches(big))

    def test_fills_crossing_edge_on_overlap(self):
        across = connectivity.Polygon(
            [(0, 4), (10, 4), (10, 6), (0, 6)], {"In1.Cu"})
        down = connectivity.Polygon(
            [(4, 0), (6, 0), (6, 10), (4, 10)], {"In1.Cu"})
        self.assertTrue(across.touches(down))
        self.assertTrue(down.touches(across))

    def test_disjoint_fills_do_not_overlap(self):
        left = connectivity.Polygon(
            [(0, 0), (1, 0), (1, 1), (0, 1)], {"In1.Cu"})
        right = connectivity.Polygon(
            [(5, 0), (6, 0), (6, 1), (5, 1)], {"In1.Cu"})
        self.assertFalse(left.touches(right))
        self.assertFalse(right.touches(left))

    def test_touch_tolerance_is_one_micron(self):
        left = connectivity.Capsule((0, 0), (0, 0), 0.5, {"F.Cu"})
        within = connectivity.Capsule((1.0005, 0), (1.0005, 0), 0.5,
                                      {"F.Cu"})
        outside = connectivity.Capsule((1.0015, 0), (1.0015, 0), 0.5,
                                       {"F.Cu"})
        self.assertTrue(left.touches(within))
        self.assertFalse(left.touches(outside))

    def test_a_track_between_two_pads_is_one_island(self):
        self.assertEqual(connectivity.opens(parse(BOARD)), {})

    def test_deleting_the_track_splits_the_net(self):
        broken = connectivity.opens(drop(parse(BOARD), "segment", "DATA"))
        self.assertEqual(broken, {"DATA": [[("R1", "1")], [("R1", "2")]]})

    def test_an_arc_track_between_two_pads_is_one_island(self):
        self.assertEqual(connectivity.opens(parse(ARC_BOARD)), {})

    def test_deleting_the_arc_splits_the_net(self):
        broken = connectivity.opens(drop(parse(ARC_BOARD), "arc", "DATA"))
        self.assertEqual(broken, {"DATA": [[("R1", "1")], [("R1", "2")]]})

    def test_arc_chords_stay_within_one_micron(self):
        copper, _, _ = connectivity.board_copper(parse(ARC_BOARD))
        chords = [item for item in copper["DATA"] if isinstance(item, connectivity.Capsule)
                  and item.start != item.end]
        self.assertEqual(chords[0].start, (0, 0))
        self.assertEqual(chords[-1].end, (10, 0))
        self.assertIn((5, 3), [item.end for item in chords])
        for first, second in zip(chords, chords[1:]):
            self.assertEqual(first.end, second.start)
        for chord in chords:
            midpoint = tuple((a + b) / 2 for a, b in zip(chord.start, chord.end))
            sagitta = 17 / 3 - math.dist(midpoint, (5, -8 / 3))
            self.assertLessEqual(abs(sagitta), connectivity.TOUCH_TOLERANCE)

    def test_native_numeric_net_references_join_named_pads(self):
        self.assertEqual(connectivity.opens(parse(NATIVE_BOARD)), {})

    def test_a_numeric_board_reports_the_open_by_declared_name(self):
        broken = connectivity.opens(drop(parse(NATIVE_BOARD), "segment", "1"))
        self.assertEqual(broken, {"DATA": [[("R1", "1")], [("R1", "2")]]})

    def test_round_pads_do_not_bridge_bare_laminate(self):
        broken = connectivity.opens(parse(ROUND_PAD_BOARD))
        self.assertEqual(broken,
                         {"ROUND": [[("J1", "1")], [("J1", "2")]]})

    def test_rectangular_lands_do_not_bridge_bare_laminate(self):
        broken = connectivity.opens(parse(SOIC_BOARD))
        self.assertEqual(broken, {"SIG": [[("U1", "1")], [("U1", "11")]]})

    def test_a_track_landing_on_a_rectangular_land_connects(self):
        text = SOIC_BOARD.replace("(end 0 -2.065)", "(end 0 -2.7)")
        self.assertEqual(connectivity.opens(parse(text)), {})

    def test_a_solid_pour_connects_through_hole_pads(self):
        self.assertEqual(connectivity.opens(parse(POUR_BOARD)), {})

    def test_emptying_the_pour_splits_its_pads(self):
        root = [node if not (isinstance(node, list) and node and node[0] == "zone")
                else [child for child in node
                      if not (isinstance(child, list) and child
                              and child[0] == "filled_polygon")]
                for node in parse(POUR_BOARD)]
        self.assertEqual(sorted(connectivity.opens(root)), ["GND"])

    def test_a_multi_layer_pour_carries_its_pads(self):
        # KiCad writes a zone spanning several layers as `(layers ...)`; reading
        # only the singular spelling breaks the walk instead of scoring it.
        self.assertEqual(connectivity.opens(parse(SPAN_POUR_BOARD)), {})
        self.assertEqual(connectivity.opens(parse(STAR_POUR_BOARD)), {})

    def test_a_multi_layer_pour_is_placed_on_every_layer_it_spans(self):
        copper, _, _ = connectivity.board_copper(parse(SPAN_POUR_BOARD))
        fills = [item for item in copper["GND"]
                 if isinstance(item, connectivity.Polygon)]
        self.assertEqual(sorted(layer for fill in fills for layer in fill.layers),
                         ["B.Cu", "In1.Cu"])

    def test_a_named_fill_only_connects_pads_on_its_own_layer(self):
        text = SPAN_POUR_BOARD.replace(
            '(filled_polygon', '(filled_polygon (layer "B.Cu")', 1)
        self.assertEqual(connectivity.opens(parse(text)), {})
        text = text.replace(
            '(pad "1" thru_hole circle (at 0 0) (size 1.6 1.6) (layers "*.Cu")',
            '(pad "1" smd circle (at 0 0) (size 1.6 1.6) (layers "In1.Cu")')
        self.assertEqual(connectivity.opens(parse(text)),
                         {"GND": [[("J1", "1")], [("J1", "2")]]})

    def test_unnetted_pads_do_not_form_an_open(self):
        root = parse(NATIVE_BOARD)
        self.assertEqual(connectivity.opens(root), {})
        root.append(sexp.parse_one('(net 0 "")'))
        root.append(sexp.parse_one('''(footprint "R" (layer "F.Cu") (at 0 0)
            (property "Reference" "R0")
            (pad "1" smd rect (at 0 10) (size 0.5 0.5) (layers "F.Cu") (net 0 ""))
            (pad "2" smd rect (at 10 10) (size 0.5 0.5) (layers "F.Cu") (net 0 "")))'''))
        self.assertEqual(connectivity.opens(root), {})

    def test_a_zone_naming_no_copper_layer_is_refused(self):
        text = SPAN_POUR_BOARD.replace('(layers "In1.Cu" "B.Cu")', "(layers)")
        with self.assertRaisesRegex(ValueError, "outside the copper stack"):
            connectivity.opens(parse(text))

    def test_a_zone_declaring_both_spellings_is_refused(self):
        text = SPAN_POUR_BOARD.replace(
            '(layers "In1.Cu" "B.Cu")',
            '(layer "In1.Cu") (layers "In1.Cu" "B.Cu")')
        with self.assertRaisesRegex(ValueError, "2 times, not once"):
            connectivity.opens(parse(text))

    def test_a_back_side_footprint_is_refused(self):
        text = BOARD.replace("(layer \"F.Cu\")", "(layer \"B.Cu\")", 1)
        with self.assertRaisesRegex(ValueError, "R1 layer is B.Cu"):
            connectivity.opens(parse(text))

    def test_a_footprint_without_a_layer_is_refused(self):
        text = BOARD.replace("\t\t(layer \"F.Cu\")\n", "", 1)
        with self.assertRaisesRegex(ValueError, "R1 layer is missing"):
            connectivity.opens(parse(text))

    def test_a_via_bridges_the_layers_it_spans(self):
        text = BOARD.replace(
            '(pad "2" smd rect (at 10 0) (size 0.5 0.5) (layers "F.Cu")',
            '(pad "2" smd rect (at 10 0) (size 0.5 0.5) (layers "B.Cu")')
        self.assertEqual(sorted(connectivity.opens(parse(text))), ["DATA"])
        bridged = text.replace(
            "(segment (start 0 0) (end 10 0) (width 0.2) (layer \"F.Cu\") (net \"/DATA\"))",
            "(segment (start 0 0) (end 10 0) (width 0.2) (layer \"F.Cu\") (net \"/DATA\"))\n"
            "\t(via (at 10 0) (size 0.45) (drill 0.2) (layers \"F.Cu\" \"B.Cu\")"
            " (net \"/DATA\"))")
        self.assertEqual(connectivity.opens(parse(bridged)), {})


class PadPresenceTests(unittest.TestCase):
    def test_unused_inner_lands_do_not_connect_through_cached_pours(self):
        for remove, expected in (("no", {}),
                                 ("yes", {"GND": [[("J1", "1")], [("J1", "2")]]}),
                                 ("", {"GND": [[("J1", "1")], [("J1", "2")]]})):
            text = POUR_BOARD.replace('(layers "*.Cu")',
                                      f'(layers "*.Cu") (remove_unused_layers {remove})')
            with self.subTest(remove=remove):
                self.assertEqual(connectivity.opens(parse(text)), expected)

    def test_keep_end_layers_retains_only_external_copper(self):
        text = POUR_BOARD.replace('(layers "*.Cu")',
                                  '(layers "*.Cu") (remove_unused_layers yes) (keep_end_layers yes)')
        _, pads, _ = connectivity.board_copper(parse(text))
        for _, item in pads['GND']:
            self.assertEqual(item.layers, {"F.Cu", "B.Cu"})
        self.assertIn("GND", connectivity.opens(parse(text)))
        self.assertEqual(connectivity.opens(parse(text.replace('"In1.Cu")', '"B.Cu")'))), {})

    def test_saved_zone_connections_retain_the_named_inner_lands(self):
        text = POUR_BOARD.replace('(layers "*.Cu")',
                                  '(layers "*.Cu") (remove_unused_layers yes) '
                                  '(zone_layer_connections "In1.Cu")')
        _, pads, _ = connectivity.board_copper(parse(text))
        for _, item in pads['GND']:
            self.assertEqual(item.layers, {"In1.Cu"})
        self.assertEqual(connectivity.opens(parse(text)), {})

    def test_removed_layers_follow_direct_nonzone_contacts(self):
        for contact in (
                '(segment (start 5 5) (end 25 5) (width 0.2) (layer "In1.Cu") (net "/GND"))',
                '(arc (start 5 5) (mid 15 8) (end 25 5) (width 0.2) (layer "In1.Cu") (net "/GND"))',
                '(via (at 5 5) (size 0.6) (layers "In1.Cu" "B.Cu") (net "/GND"))'
                '(via (at 25 5) (size 0.6) (layers "In1.Cu" "B.Cu") (net "/GND"))',
                '(footprint "contact" (layer "F.Cu") (at 0 0) '
                '(property "Reference" "J2") '
                '(pad "1" smd circle (at 5 5) (size 1 1) (layers "In1.Cu") (net "/GND"))'
                '(pad "2" smd circle (at 25 5) (size 1 1) (layers "In1.Cu") (net "/GND")))'):
            text = POUR_BOARD.replace('(layers "*.Cu")',
                                      '(layers "*.Cu") (remove_unused_layers yes)')
            text = text.rstrip()[:-1] + contact + ')'
            with self.subTest(kind=contact.split()[0]):
                _, pads, _ = connectivity.board_copper(parse(text))
                for (ref, _), item in pads['GND']:
                    if ref == 'J1':
                        self.assertIn('In1.Cu', item.layers)
                        self.assertNotIn('F.Cu', item.layers)
                self.assertEqual(connectivity.opens(parse(text)), {})

    def test_disconnected_custom_primitives_do_not_bridge_a_gap(self):
        for bridge, expected_open in (("", True),
                                      ('(gr_line (start 0 0) (end 10 0) (width 0.2))', False)):
            text = BOARD.replace(
                'rect (at 0 0) (size 0.5 0.5)',
                'custom (at 0 0) (size 0.1 0.1) (primitives '
                '(gr_rect (start 9.8 -0.2) (end 10.2 0.2) (width 0) (fill yes))' + bridge + ')')
            root = drop(parse(text), 'segment', 'DATA')
            with self.subTest(bridge=bool(bridge)):
                self.assertEqual('DATA' in connectivity.opens(root), expected_open)


class EmptyScanTests(unittest.TestCase):
    def test_empty_board_fails(self):
        with self.assertRaisesRegex(ValueError, "nothing to analyze"):
            connectivity.opens(parse(EMPTY_BOARD))

    def test_main_exits_nonzero_on_an_empty_board(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "empty.kicad_pcb"
            path.write_text(EMPTY_BOARD, encoding="utf-8")
            with contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(connectivity.main([str(path)]), 2)

    def test_main_reports_missing_or_unreadable_board(self):
        for error in (FileNotFoundError("missing board"), PermissionError("unreadable board")):
            output = io.StringIO()
            with self.subTest(error=error), mock.patch("builtins.open", side_effect=error), \
                    contextlib.redirect_stderr(output):
                self.assertEqual(connectivity.main(["test.kicad_pcb"]), 2)
            self.assertIn(f"test.kicad_pcb: {error}", output.getvalue())
            self.assertNotIn("Traceback", output.getvalue())

    def test_main_exits_two_on_a_schematic_or_layerless_board(self):
        documents = {
            "schematic": '(kicad_sch (version 20231120) (generator "eeschema"))',
            "layerless": BOARD[:BOARD.index("\t(layers")]
            + BOARD[BOARD.index("\t(footprint"):],
        }
        for label, text in documents.items():
            output = io.StringIO()
            with self.subTest(document=label), tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp) / "input.kicad_pcb"
                path.write_text(text, encoding="utf-8")
                with contextlib.redirect_stderr(output):
                    self.assertEqual(connectivity.main([str(path)]), 2)
            self.assertNotIn("Traceback", output.getvalue())

    def test_main_exits_two_on_a_file_holding_no_document(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "blank.kicad_pcb"
            path.write_text("   ", encoding="utf-8")
            with contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(connectivity.main([str(path)]), 2)


class RoutedBoardTests(unittest.TestCase):
    """The shipped board's copper."""

    @classmethod
    def setUpClass(cls):
        cls.root = parse(ROUTED.read_text(encoding="utf-8"))

    def test_every_net_on_the_routed_board_is_one_island(self):
        broken = connectivity.opens(self.root)
        self.assertEqual(broken, {}, connectivity.report(broken))

    def test_deleting_the_strip_data_track_is_caught(self):
        broken = connectivity.opens(drop(self.root, "segment", "DATA"))
        self.assertIn("DATA", broken)


if __name__ == "__main__":
    unittest.main()
