import sys
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
\t\t(at 0 0)
\t\t(property "Reference" "R1")
\t\t(pad "1" smd rect (at 0 0) (size 0.5 0.5) (layers "F.Cu") (net "/DATA"))
\t\t(pad "2" smd rect (at 10 0) (size 0.5 0.5) (layers "F.Cu") (net "/DATA"))
\t)
\t(segment (start 0 0) (end 10 0) (width 0.2) (layer "F.Cu") (net "/DATA"))
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


def parse(text):
    return sexp.parse(text)[0]


def drop(root, kind, net):
    return [node for node in root
            if not (isinstance(node, list) and node and node[0] == kind
                    and connectivity.net_name(node) == net)]


class SyntheticBoardTests(unittest.TestCase):
    def test_a_track_between_two_pads_is_one_island(self):
        self.assertEqual(connectivity.opens(parse(BOARD)), {})

    def test_deleting_the_track_splits_the_net(self):
        broken = connectivity.opens(drop(parse(BOARD), "segment", "DATA"))
        self.assertEqual(broken, {"DATA": [[("R1", "1")], [("R1", "2")]]})

    def test_a_pour_carries_through_hole_pads_it_voids_around(self):
        self.assertEqual(connectivity.opens(parse(POUR_BOARD)), {})

    def test_emptying_the_pour_splits_its_pads(self):
        root = [node if not (isinstance(node, list) and node and node[0] == "zone")
                else [child for child in node
                      if not (isinstance(child, list) and child
                              and child[0] == "filled_polygon")]
                for node in parse(POUR_BOARD)]
        self.assertEqual(sorted(connectivity.opens(root)), ["GND"])

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
