import math
import sys
import unittest
from pathlib import Path
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pcb
import sexp


REPO_ROOT = Path(__file__).resolve().parents[4]
ROUTED = REPO_ROOT / "hardware" / "phantasm" / "phantasm.kicad_pcb"
UNPLACED = (REPO_ROOT / "hardware" / "phantasm" / "unplaced" /
            "phantasm_unplaced.kicad_pcb")


def _child(node, key):
    return next(
        child[1:]
        for child in node
        if isinstance(child, list) and child and child[0] == key
    )


def _footprint(board, ref):
    return next(
        node
        for node in board
        if isinstance(node, list) and node and node[0] == "footprint"
        and any(isinstance(child, list) and child
                and child[0] == "property"
                and child[1:3] == ["Reference", ref]
                for child in node)
    )


def _without_uuids(node):
    if not isinstance(node, list):
        return node
    return tuple(
        _without_uuids(child)
        for child in node
        if not (isinstance(child, list) and child and child[0] == "uuid")
    )


def _silk_outline(footprint):
    return [
        node
        for node in footprint
        if isinstance(node, list) and node
        and node[0] in ("fp_line", "fp_arc")
        and _child(node, "layer")[0] == "F.SilkS"
    ]


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


class RoutedTraceTests(unittest.TestCase):
    def test_master_enable_bottom_route_uses_an_obtuse_corner(self):
        board_path = REPO_ROOT / "hardware" / "phantasm" / "phantasm.kicad_pcb"
        board = sexp.parse(board_path.read_text(encoding="utf-8"))[0]
        route_uuids = {
            "3dca6b57-636f-4af4-a67b-e4e8047ab42c",
            "8660a7d1-70bc-4768-a0fb-1074dad66d08",
        }
        segments = [
            node
            for node in board
            if isinstance(node, list)
            and node
            and node[0] == "segment"
            and str(_child(node, "uuid")[0]) in route_uuids
        ]
        self.assertEqual(len(segments), 2)

        endpoints = [
            {
                tuple(map(float, _child(segment, "start")[:2])),
                tuple(map(float, _child(segment, "end")[:2])),
            }
            for segment in segments
        ]
        corner, = endpoints[0] & endpoints[1]
        outer = [next(iter(points - {corner})) for points in endpoints]
        vectors = [(x - corner[0], y - corner[1]) for x, y in outer]
        cosine = sum(a * b for a, b in zip(*vectors)) / math.prod(
            math.hypot(*vector) for vector in vectors
        )
        angle = math.degrees(math.acos(max(-1.0, min(1.0, cosine))))
        self.assertGreater(angle, 90.0)


class CathodeMarkTests(unittest.TestCase):
    """D_BUS is a unidirectional TVS: reversed, it clamps the bus below the
    AHCT HIGH threshold, and the silk bar beside pad 1 is the only way to read
    its orientation on an assembled board."""

    def test_d_bus_carries_a_silk_bar_beside_its_cathode(self):
        board_path = REPO_ROOT / "hardware" / "phantasm" / "phantasm.kicad_pcb"
        board = sexp.parse(board_path.read_text(encoding="utf-8"))[0]
        footprint, = [
            node
            for node in board
            if isinstance(node, list) and node and node[0] == "footprint"
            and any(isinstance(child, list) and child
                    and child[0] == "property" and child[1:3] == ["Reference",
                                                                  "D_BUS"]
                    for child in node)
        ]
        cathode_x, = [
            float(_child(pad, "at")[0])
            for pad in footprint
            if isinstance(pad, list) and pad and pad[0] == "pad"
            and str(pad[1]) == "1"
        ]
        bars = [
            node
            for node in footprint
            if isinstance(node, list) and node and node[0] == "fp_line"
            and _child(node, "layer")[0] == "F.SilkS"
            and max(float(_child(node, "start")[0]),
                    float(_child(node, "end")[0])) <= cathode_x
        ]
        self.assertEqual(len(bars), 1)


class IdStrapSilkscreenTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.boards = {
            path: sexp.parse(path.read_text(encoding="utf-8"))[0]
            for path in (ROUTED, UNPLACED)
        }

    def test_every_id_strap_has_the_same_outline(self):
        for path, board in self.boards.items():
            expected = [_without_uuids(node)
                        for node in _silk_outline(_footprint(board, "JP_ID0"))]
            self.assertEqual(
                [node[0] for node in expected].count("fp_line"), 4)
            self.assertEqual(
                [node[0] for node in expected].count("fp_arc"), 4)
            for ref in ("JP_ID1", "JP_ID2", "JP_SHLD"):
                with self.subTest(board=path.name, ref=ref):
                    actual = [_without_uuids(node) for node in
                              _silk_outline(_footprint(board, ref))]
                    self.assertEqual(actual, expected)

    def test_generator_keeps_jp_id2_outline_on_silkscreen(self):
        source = _footprint(self.boards[ROUTED], "JP_ID0")
        libid = "Test:Jumper"
        with mock.patch.dict(pcb._MOD_CACHE, {libid: source}):
            generated = pcb.embed(libid, "JP_ID2", "Jumper", 0, 0, 0,
                                  {}, {})
        self.assertEqual(len(_silk_outline(generated)), 8)

    def test_boards_describe_the_id2_half_of_the_truth_table(self):
        expected = "N8 ID2 OPEN=0-3 GND=4-7; M=OPEN; SHLD=M"
        for path, board in self.boards.items():
            texts = [str(node[1]) for node in board
                     if isinstance(node, list) and node
                     and node[0] == "gr_text"
                     and _child(node, "layer")[0] == "B.SilkS"]
            with self.subTest(board=path.name):
                self.assertIn(expected, texts)


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
