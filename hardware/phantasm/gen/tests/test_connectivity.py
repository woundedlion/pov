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
        first = connectivity.pad_capsule(pads[0], (0, 0), 0, ["F.Cu"])
        self.assertAlmostEqual(first.radius, 1.35)

    def test_a_custom_pad_reaches_its_furthest_primitive(self):
        root = parse(CUSTOM_PAD_BOARD)
        pad = connectivity.F(
            connectivity.F(root, "footprint")[0], "pad")[0]
        capsule = connectivity.pad_capsule(pad, (0, 0), 0, ["F.Cu"])
        # (size 1 0.5) circumscribes to 0.559; the poly corner is 0.901 out.
        self.assertAlmostEqual(capsule.radius, math.hypot(0.5, 0.75))

    def test_fill_overlap_is_symmetric(self):
        big = connectivity.Fill([(0, 0), (10, 0), (10, 10), (0, 10)], "In1.Cu")
        small = connectivity.Fill([(4, 4), (6, 4), (6, 6), (4, 6)], "In1.Cu")
        self.assertTrue(big.touches(small))
        self.assertTrue(small.touches(big))

    def test_fills_crossing_edge_on_overlap(self):
        across = connectivity.Fill([(0, 4), (10, 4), (10, 6), (0, 6)], "In1.Cu")
        down = connectivity.Fill([(4, 0), (6, 0), (6, 10), (4, 10)], "In1.Cu")
        self.assertTrue(across.touches(down))
        self.assertTrue(down.touches(across))

    def test_disjoint_fills_do_not_overlap(self):
        left = connectivity.Fill([(0, 0), (1, 0), (1, 1), (0, 1)], "In1.Cu")
        right = connectivity.Fill([(5, 0), (6, 0), (6, 1), (5, 1)], "In1.Cu")
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

    def test_round_pads_do_not_bridge_bare_laminate(self):
        broken = connectivity.opens(parse(ROUND_PAD_BOARD))
        self.assertEqual(broken,
                         {"ROUND": [[("J1", "1")], [("J1", "2")]]})

    def test_a_pour_carries_through_hole_pads_it_voids_around(self):
        self.assertEqual(connectivity.opens(parse(POUR_BOARD)), {})

    def test_emptying_the_pour_splits_its_pads(self):
        root = [node if not (isinstance(node, list) and node and node[0] == "zone")
                else [child for child in node
                      if not (isinstance(child, list) and child
                              and child[0] == "filled_polygon")]
                for node in parse(POUR_BOARD)]
        self.assertEqual(sorted(connectivity.opens(root)), ["GND"])

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
