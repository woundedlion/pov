import contextlib
import io
import math
import sys
import tempfile
import unittest
from pathlib import Path

GEN = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN))

import connectivity  # noqa: E402
import sexp  # noqa: E402

ROUTED = GEN.parent / "phantasm.kicad_pcb"

# One F.Cu track between two SMD pads, plus a second pad pair the track misses.
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
# the pour voided around the drill the way KiCad thermally relieves it.
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


class SyntheticBoardTests(unittest.TestCase):
    def test_pad_geometry_uses_the_circumscribed_radius(self):
        root = parse(ROUND_PAD_BOARD)
        pads = list(connectivity.F(
            connectivity.F(root, "footprint")[0], "pad"))
        first = connectivity.pad_copper(pads[0], (0, 0), 0, ["F.Cu"])
        self.assertAlmostEqual(first.radius, 1.35)

    def test_a_rectangular_land_is_the_rotated_rectangle(self):
        root = parse(ROTATED_LAND_BOARD)
        pad = connectivity.F(
            connectivity.F(root, "footprint")[0], "pad")[0]
        land = connectivity.pad_copper(pad, (29.7, 22.8), 90, ["F.Cu"])
        xs = sorted(round(x, 3) for x, _ in land.polygon)
        ys = sorted(round(y, 3) for _, y in land.polygon)
        self.assertEqual((xs[0], xs[-1]), (29.225, 30.175))
        self.assertEqual((ys[0], ys[-1]), (23.025, 24.225))

    def test_a_trapezoid_land_boxes_in_its_rect_delta(self):
        text = ROTATED_LAND_BOARD.replace(
            'roundrect (at -0.825 0 90) (size 1.2 0.95)',
            'trapezoid (at -0.825 0) (size 1 2) (rect_delta 0.4 0)')
        pad = connectivity.F(
            connectivity.F(parse(text), "footprint")[0], "pad")[0]
        land = connectivity.pad_copper(pad, (0, 0), 0, ["F.Cu"])
        xs = sorted(round(x, 3) for x, _ in land.polygon)
        ys = sorted(round(y, 3) for _, y in land.polygon)
        self.assertEqual((xs[0], xs[-1]), (-1.325, -0.325))
        self.assertEqual((ys[0], ys[-1]), (-1.2, 1.2))

    def test_a_custom_pad_reaches_its_furthest_primitive(self):
        root = parse(CUSTOM_PAD_BOARD)
        pad = connectivity.F(
            connectivity.F(root, "footprint")[0], "pad")[0]
        capsule = connectivity.pad_copper(pad, (0, 0), 0, ["F.Cu"])
        # (size 1 0.5) circumscribes to 0.559; the poly corner is 0.901 out.
        self.assertAlmostEqual(capsule.radius, math.hypot(0.5, 0.75))

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

    def test_a_pour_carries_through_hole_pads_it_voids_around(self):
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


class RoutedBoardTests(unittest.TestCase):
    """The shipped board's copper, which no net-label gate reads."""

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
