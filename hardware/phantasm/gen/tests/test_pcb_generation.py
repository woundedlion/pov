"""Self-tests for the board generator.

pcb.py is the only writer of phantasm.kicad_pcb and refuses to overwrite the
committed file, so nothing else in the repo executes it: every other gate reads
the committed board. These tests run the generator into a temporary directory
and push what it wrote back through the repo's KiCad-free readers.

Generating needs KiCad's stock symbol and footprint libraries plus a kicad-cli
on the pin for the netlist export; without them the whole class is skipped.
"""
import contextlib
import io
import math
import os
import sys
import tempfile
import unittest
import unittest.mock as mock
from decimal import Decimal
from pathlib import Path

GEN = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN))

import board_metadata   # noqa: E402
import connectivity     # noqa: E402
import pcb              # noqa: E402
import sexp             # noqa: E402
from kicad_common import F, is_copper_pour, kicad_cli  # noqa: E402

STOCK_SYMBOLS = os.path.isdir(sexp.KICAD_SHARE)
STOCK_FOOTPRINTS = os.path.isdir(pcb.FP_DIR)


def pinned_kicad_cli():
    """True when a kicad-cli on the pin resolves; kicad_cli() exits when not."""
    try:
        return bool(kicad_cli())
    except SystemExit:
        return False


GENERATES = STOCK_SYMBOLS and STOCK_FOOTPRINTS and pinned_kicad_cli()
GENERATES_REASON = (
    f"KiCad {sexp.KICAD_MAJOR} stock libraries or kicad-cli not found")


def generate(out, unplaced=False):
    """Run the generator into `out`; return the board path."""
    with mock.patch.object(pcb, "OUT", out), \
            contextlib.redirect_stdout(io.StringIO()):
        pcb.main(unplaced=unplaced, force=True, force_teensy_library=True)
    return os.path.join(out, pcb.UNPLACED_FILE if unplaced else pcb.PCB_FILE)


def read(path):
    return sexp.parse(Path(path).read_text(encoding="utf-8"))[0]


def reference(footprint):
    for child in F(footprint, "property"):
        if len(child) > 2 and child[1] == "Reference":
            return str(child[2])
    return "?"


def zone_polygon(zone):
    points = F(F(F(zone, "polygon")[0], "pts")[0], "xy")
    return [(float(point[1]), float(point[2])) for point in points]


def _graphic_points(node):
    """Local-frame corner points of one footprint graphic."""
    if str(node[0]) == "fp_circle":
        centre = sexp.val(node, "center")
        rim = sexp.val(node, "end")
        x, y = float(centre[0]), float(centre[1])
        radius = math.hypot(float(rim[0]) - x, float(rim[1]) - y)
        return [(x - radius, y - radius), (x + radius, y - radius),
                (x + radius, y + radius), (x - radius, y + radius)]
    points = [(float(value[0]), float(value[1]))
              for value in (sexp.val(node, key)
                            for key in ("start", "mid", "end", "center"))
              if value]
    for vertex in (F(F(node, "pts")[0], "xy") if F(node, "pts") else []):
        points.append((float(vertex[1]), float(vertex[2])))
    if str(node[0]) == "fp_rect" and len(points) == 2:
        (x0, y0), (x1, y1) = points
        points = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
    return points


def courtyard_box(footprint):
    """Placed bounding box of a footprint's courtyard, or None if it draws none.

    The box circumscribes the courtyard, so it over-reports a clash between two
    interlocking outlines; nothing on this board is placed that tightly.
    """
    placement = sexp.val(footprint, "at")
    origin = (float(placement[0]), float(placement[1]))
    rotation = float(placement[2]) if len(placement) > 2 else 0.0
    xs, ys = [], []
    for child in footprint:
        if not (isinstance(child, list) and child):
            continue
        layer = sexp.val(child, "layer")
        if not layer or str(layer[0]) not in pcb.COURTYARD_LAYERS:
            continue
        for point in _graphic_points(child):
            x, y = connectivity._rotate(point, rotation)
            xs.append(origin[0] + x)
            ys.append(origin[1] + y)
    return (min(xs), min(ys), max(xs), max(ys)) if xs else None


@unittest.skipUnless(GENERATES, GENERATES_REASON)
class GeneratedBoardTests(unittest.TestCase):
    """The placed draft `pcb.py --force` emits, read back without KiCad."""

    @classmethod
    def setUpClass(cls):
        cls.out = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.out.cleanup)
        cls.path = generate(cls.out.name)
        cls.root = read(cls.path)
        cls.metadata = board_metadata.load_board(Path(cls.path))
        cls.length = float(cls.metadata.width_mm)

    def test_writes_the_board_and_the_teensy_library(self):
        for name in (pcb.PCB_FILE, "fp-lib-table",
                     os.path.join("phantasm.pretty", "Teensy4.0.kicad_mod")):
            with self.subTest(artifact=name):
                self.assertTrue(
                    os.path.exists(os.path.join(self.out.name, name)), name)

    def test_the_teensy_library_parses_as_a_footprint(self):
        module = read(os.path.join(self.out.name, "phantasm.pretty",
                                   "Teensy4.0.kicad_mod"))
        self.assertEqual(str(module[0]), "footprint")
        self.assertEqual(str(module[1]), "Teensy4.0")

    def test_emits_the_declared_stackup(self):
        self.assertEqual(self.metadata.height_mm, Decimal(str(pcb.PCB_W)))
        self.assertEqual(self.metadata.thickness_mm, Decimal("1.6"))
        self.assertEqual(self.metadata.copper_layers, pcb.copper_layer_names())
        self.assertEqual(self.metadata.copper_finish, "ENIG")

    def test_the_draft_carries_no_routing(self):
        self.assertEqual(self.metadata.track_segments, 0)
        self.assertEqual(self.metadata.vias, 0)

    def test_places_every_footprint_on_the_front(self):
        self.assertEqual(self.metadata.footprint_sides[1], ("B.Cu", 0))
        self.assertGreater(self.metadata.footprint_sides[0][1], 0)

    def test_every_padded_net_resolves_to_a_declaration(self):
        copper, pads, names = connectivity.board_copper(self.root)
        self.assertTrue(pads)
        self.assertEqual(sorted(set(pads) - set(names)), [])

    def test_net_ids_are_unique_and_dense(self):
        ids = [int(node[1]) for node in F(self.root, "net")]
        self.assertEqual(sorted(ids), list(range(len(ids))))

    def test_pours_both_inner_reference_planes(self):
        self.assertEqual(self.metadata.copper_pours,
                         len(pcb.GROUND_PLANE_LAYERS))
        self.assertEqual(
            self.metadata.pour_layers,
            tuple((layer, 1) for layer in pcb.GROUND_PLANE_LAYERS))
        outline = [(0.0, 0.0), (self.length, 0.0),
                   (self.length, pcb.PCB_W), (0.0, pcb.PCB_W)]
        for zone in F(self.root, "zone"):
            if is_copper_pour(zone):
                with self.subTest(zone=str(sexp.val(zone, "name")[0])):
                    self.assertEqual(str(sexp.val(zone, "net_name")[0]),
                                     pcb.GROUND_NET)
                    self.assertEqual(zone_polygon(zone), outline)

    def test_reserves_a_rule_area_at_every_mounting_hole(self):
        expected = pcb.keepout_rects(self.length)
        self.assertEqual(self.metadata.rule_areas, len(expected))
        self.assertEqual(
            self.metadata.rule_area_layers,
            tuple((layer, len(expected)) for layer in pcb.copper_layer_names()))
        found = {}
        for zone in F(self.root, "zone"):
            if is_copper_pour(zone):
                continue
            xs = [x for x, _ in zone_polygon(zone)]
            ys = [y for _, y in zone_polygon(zone)]
            found[str(sexp.val(zone, "name")[0]).split()[0]] = (
                min(xs), min(ys), max(xs), max(ys))
        for ref, rect in expected.items():
            with self.subTest(hole=ref):
                self.assertEqual(
                    tuple(round(value, 3) for value in found[ref]),
                    tuple(round(value, 3) for value in rect))

    def test_stamps_the_revision_on_the_back_silkscreen(self):
        texts = [str(node[1]) for node in F(self.root, "gr_text")
                 if str(sexp.val(node, "layer")[0]) == "B.SilkS"]
        self.assertIn(pcb.SILK_REVISION, texts)

    def test_no_two_courtyards_overlap(self):
        boxes = {}
        for footprint in F(self.root, "footprint"):
            box = courtyard_box(footprint)
            self.assertIsNotNone(box, reference(footprint))
            boxes[reference(footprint)] = box
        refs = sorted(boxes)
        overlaps = [f"{a}/{b}"
                    for index, a in enumerate(refs) for b in refs[index + 1:]
                    if pcb._boxes_overlap(boxes[a], boxes[b])]
        self.assertEqual(overlaps, [])


@unittest.skipUnless(GENERATES, GENERATES_REASON)
class UnplacedBoardTests(unittest.TestCase):
    """The autoplacer upload `pcb.py --unplaced` emits."""

    @classmethod
    def setUpClass(cls):
        cls.out = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.out.cleanup)
        cls.root = read(generate(cls.out.name, unplaced=True))

    def test_locks_the_mechanical_placements(self):
        placed = {}
        for footprint in F(self.root, "footprint"):
            at = sexp.val(footprint, "at")
            placed[reference(footprint)] = (
                float(at[0]), float(at[1]),
                float(at[2]) if len(at) > 2 else 0.0)
        for ref, (x, y, rot) in pcb.QUILTER_FIXED.items():
            with self.subTest(ref=ref):
                self.assertEqual(placed[ref], (float(x), float(y), float(rot)))

    def test_labels_the_id_straps_on_the_front_silkscreen(self):
        texts = [str(node[1]) for node in F(self.root, "gr_text")
                 if str(sexp.val(node, "layer")[0]) == "F.SilkS"]
        self.assertTrue({"ID0", "ID1", "ID2", "SHLD", "SYNC IN", "SYNC OUT"}
                        <= set(texts), sorted(texts))


@unittest.skipUnless(GENERATES, GENERATES_REASON)
class OrphanPadTests(unittest.TestCase):
    """A netlist pin with no pad of that name would drop its net silently."""

    def test_a_netlist_pin_with_no_pad_is_refused(self):
        build_nets = pcb.build_nets

        def with_a_ghost_pin(nlroot):
            pad_net, netid = build_nets(nlroot)
            return pad_net | {("U_MCU", "999"): pcb.GROUND_NET}, netid

        out = self.enterContext(tempfile.TemporaryDirectory())
        with mock.patch.object(pcb, "build_nets", with_a_ghost_pin), \
                self.assertRaises(SystemExit) as caught:
            generate(out)
        self.assertIn("U_MCU.999", str(caught.exception))


if __name__ == "__main__":
    unittest.main()
