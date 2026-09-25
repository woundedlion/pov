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
import board            # noqa: E402
import builder          # noqa: E402
import check            # noqa: E402
import fab              # noqa: E402
import connectivity     # noqa: E402
import pcb              # noqa: E402
import sexp             # noqa: E402
from kicad_common import F, export_netlist, is_copper_pour, kicad_cli  # noqa: E402

COMMITTED_PCB = GEN.parent / pcb.PCB_FILE

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
    schematic = os.path.join(out, "phantasm.kicad_sch")
    with mock.patch.object(board, "OUT", out), \
            mock.patch.object(board, "SCH", schematic), \
            contextlib.redirect_stdout(io.StringIO()):
        board.main(force=True)
    with mock.patch.object(pcb, "OUT", out), \
            mock.patch.object(pcb, "SCH", schematic), \
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


def assembly_exclusions(root):
    """{ref: sorted attr flags} for every footprint kept off the assembly."""
    excluded = {}
    for footprint in F(root, "footprint"):
        flags = [str(flag) for attr in F(footprint, "attr") for flag in attr[1:]]
        if "exclude_from_bom" in flags:
            excluded[reference(footprint)] = sorted(flags)
    return excluded


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


class OverwriteProtectionTests(unittest.TestCase):
    def test_existing_board_is_preserved_without_kicad(self):
        for unplaced in (False, True):
            with self.subTest(unplaced=unplaced), tempfile.TemporaryDirectory() as directory:
                target = Path(directory) / (pcb.UNPLACED_FILE if unplaced else pcb.PCB_FILE)
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(b"existing routed board\n")
                with mock.patch.object(pcb, "OUT", directory), \
                        mock.patch.object(pcb, "kicad_cli", side_effect=AssertionError("KiCad called")):
                    with self.assertRaisesRegex(SystemExit, "refusing to overwrite"):
                        pcb.main(unplaced=unplaced)
                self.assertEqual(target.read_bytes(), b"existing routed board\n")


@unittest.skipUnless(GENERATES, GENERATES_REASON)
class GeneratedBoardTests(unittest.TestCase):
    """The placed draft `pcb.py --force` emits, read back without KiCad."""

    def test_generated_pair_passes_schematic_parity(self):
        with tempfile.TemporaryDirectory() as directory:
            board_path = generate(directory)
            warnings = {("lib_footprint_mismatch", ref): 1 for ref in ("D_BUS", "U_MCU")}
            with mock.patch.object(fab, "PCB", board_path), \
                    mock.patch.object(fab, "SCH", str(Path(directory) / "phantasm.kicad_sch")), \
                    mock.patch.object(fab, "KNOWN_PARITY_WARNING_COUNTS", warnings):
                self.assertEqual(fab.run_parity(str(Path(directory) / "parity.json")),
                                 len(fab.KNOWN_PARITY_ITEMS))

    def test_refused_library_overwrite_preserves_existing_board(self):
        with tempfile.TemporaryDirectory() as directory:
            board = Path(generate(directory))
            board.write_bytes(b"routed board source of truth\r\n")
            library = Path(directory) / "phantasm.pretty" / "Teensy4.0.kicad_mod"
            library.write_bytes(b"hand-maintained footprint\n")
            before = board.read_bytes()
            with mock.patch.object(pcb, "OUT", directory), \
                    mock.patch.object(pcb, "SCH", str(Path(directory) / "phantasm.kicad_sch")), \
                    contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaises(SystemExit):
                    pcb.main(force=True)
            self.assertEqual(board.read_bytes(), before)
            self.assertEqual(library.read_bytes(), b"hand-maintained footprint\n")

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

    def test_sync_transmit_is_isolated_and_pulled_low_on_every_board(self):
        nets = {}
        for footprint in F(self.root, "footprint"):
            ref = reference(footprint)
            for pad in F(footprint, "pad"):
                name = sexp.val(pad, "net")
                if name:
                    nets.setdefault(str(name[-1]).lstrip("/"), set()).add(
                        check.node_key(ref, str(pad[1])))
        self.assertEqual(nets["FRAME_SYNC"], {"U_MCU.3", "R1", "R2", "C_SYNC"})
        self.assertEqual(nets["SYNC_TX"], {"U_MCU.4", "U1.9", "R_TX"})
        self.assertIn("R_TX", nets["GND"])
        self.assertEqual(nets["MASTER_EN"], check.EXPECT["MASTER_EN"])
        resistor = next(fp for fp in F(self.root, "footprint")
                        if reference(fp) == "R_TX")
        self.assertIn(["property", "Value", "10k"],
                      [p[:3] for p in F(resistor, "property")])
        self.assertNotIn("R_TX", assembly_exclusions(self.root))

    def test_generated_assembly_has_all_revision_parts(self):
        netlist = Path(self.out.name) / "assembly.net"
        root = export_netlist(kicad_cli(), str(Path(self.out.name) / "phantasm.kicad_sch"))
        netlist.write_text(sexp.dumps(root), encoding="utf-8")
        self.assertEqual(fab.validate_netlist_spec(netlist), len(check.EXPECT))
        components = fab.parse_components(netlist)
        fab.validate_assembled_refs(
            [ref for ref, component in components.items() if fab.is_assembled(component)],
            builder.REVISION)

    def test_old_schematic_cannot_be_stamped_with_the_new_revision(self):
        with tempfile.TemporaryDirectory() as directory, \
                mock.patch.object(pcb, "OUT", directory):
            with self.assertRaisesRegex(SystemExit, "regenerate the schematic first"):
                pcb.main()

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

    def test_board_id_has_large_legible_digits(self):
        labels = [node for node in F(self.root, "gr_text")
                  if str(node[1]) == "BOARD ID: ____"]
        self.assertEqual(len(labels), 1)
        label = labels[0]
        self.assertEqual(str(sexp.val(label, "layer")[0]), "B.SilkS")
        font = F(F(label, "effects")[0], "font")[0]
        self.assertEqual([float(v) for v in sexp.val(font, "size")], [2.0, 2.0])
        self.assertEqual(float(sexp.val(label, "at")[1]), 23.5)

    def test_reproduces_the_routed_board_assembly_exclusions(self):
        self.assertEqual(assembly_exclusions(self.root),
                         assembly_exclusions(read(COMMITTED_PCB)))

    def test_stamps_the_revision_in_the_title_block(self):
        blocks = F(self.root, "title_block")
        self.assertEqual(len(blocks), 1)
        self.assertEqual(str(sexp.val(blocks[0], "rev")[0]), builder.REVISION)

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
        comps = {}
        locked = set()
        for footprint in F(self.root, "footprint"):
            ref = reference(footprint)
            comps[ref] = (ref, str(footprint[1]), "", False)
            if sexp.val(footprint, "locked", []) == ["yes"]:
                locked.add(ref)
            at = sexp.val(footprint, "at")
            placed[ref] = (
                float(at[0]), float(at[1]),
                float(at[2]) if len(at) > 2 else 0.0)
        fixed = pcb.fixed_placements(comps)
        self.assertEqual(locked - {"H1", "H2", "H3", "H4"}, set(fixed))
        for ref, (x, y, rot) in fixed.items():
            with self.subTest(ref=ref):
                self.assertEqual(placed[ref], (float(x), float(y), float(rot)))
        for ref in pcb.QUILTER_FIXED.keys() - fixed.keys():
            with self.subTest(staged=ref):
                self.assertGreater(placed[ref][1], pcb.PCB_W)

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
