"""Self-tests for the board generator.

Generating needs KiCad's stock symbol and footprint libraries plus a kicad-cli
on the pin for the netlist export; without them classes gated on GENERATES are
skipped.
"""
import contextlib
import io
import json
import math
import os
import subprocess
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
from courtyard_bounds import courtyard_box  # noqa: E402
from kicad_common import F, export_netlist, is_copper_pour, kicad_cli  # noqa: E402

COMMITTED_PCB = GEN.parent / "1.1" / pcb.PCB_FILE

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
    return os.path.join(out, pcb.UNPLACED_FILE if unplaced else pcb.DRAFT_FILE)


def read(path):
    return sexp.parse(Path(path).read_text(encoding="utf-8"))[0]


def assembly_exclusions(root):
    """{ref: sorted attr flags} for every footprint kept off the assembly."""
    excluded = {}
    for footprint in F(root, "footprint"):
        flags = [str(flag) for attr in F(footprint, "attr") for flag in attr[1:]]
        if "exclude_from_bom" in flags:
            excluded[connectivity.footprint_reference(footprint)] = sorted(flags)
    return excluded


def zone_polygon(zone):
    points = F(F(F(zone, "polygon")[0], "pts")[0], "xy")
    return [(float(point[1]), float(point[2])) for point in points]


class TerminalBodyChecks:
    def test_terminal_bodies_have_component_reservations(self):
        footprints = {connectivity.footprint_reference(fp): fp for fp in F(self.root, "footprint")}
        for ref, count in (("J1", 2), ("J2", 3), ("J3A", 3), ("J3B", 3)):
            with self.subTest(ref=ref):
                fp = footprints[ref]
                name = f"TerminalBlock_GCT_TBC05-0{count}-1-G-G"
                self.assertEqual(str(fp[1]), f"phantasm:{name}")
                pads = F(fp, "pad")
                self.assertEqual([str(p[1]) for p in pads],
                                 [str(i + 1) for i in range(count)])
                for i, pad in enumerate(pads):
                    self.assertEqual([float(v) for v in sexp.val(pad, "at")[:2]],
                                     [0.0, i * 2.54])
                    self.assertEqual(float(sexp.val(pad, "drill")[0]), 1.3)
                reservation = (-3.75, -1.97, 3.75, (count - 1) * 2.54 + 1.97)
                self.assertEqual(pcb.fp_bbox(fp, graphic_layers=("F.CrtYd",)),
                                 reservation)
                zones = F(fp, "zone")
                self.assertEqual(len(zones), 2 if ref == "J1" else 1)
                zone = zones[0]
                self.assertEqual(sexp.val(zone, "layer"), ["F.Cu"])
                x0, y0, x1, y1 = reservation
                local_points = {(x0, y0), (x1, y0), (x1, y1), (x0, y1)}
                module = pcb.load_mod(f"phantasm:{name}")
                self.assertEqual(set(zone_polygon(F(module, "zone")[0])), local_points)
                x, y, angle = map(float, sexp.val(fp, "at"))
                angle = math.radians(angle)
                placed_points = {
                    (round(x + px * math.cos(angle) + py * math.sin(angle), 6),
                     round(y - px * math.sin(angle) + py * math.cos(angle), 6))
                    for px, py in local_points}
                self.assertEqual(set(zone_polygon(zone)), placed_points)
                keepout = F(zone, "keepout")[0]
                self.assertEqual(sexp.val(keepout, "footprints"), ["not_allowed"])
                for item in ("tracks", "vias", "pads", "copperpour"):
                    self.assertEqual(sexp.val(keepout, item), ["allowed"])
                self.assertIn(ref, assembly_exclusions(self.root))

    def test_terminal_keepout_rotates_and_moves_with_footprint(self):
        local_points = [(-3.75, -1.97), (3.75, -1.97), (3.75, 4.51), (-3.75, 4.51)]
        transforms = {0: lambda x, y: (x, y), 90: lambda x, y: (y, -x),
                      180: lambda x, y: (-x, -y), 270: lambda x, y: (-y, x)}
        for angle, transform in transforms.items():
            with self.subTest(angle=angle):
                fp = pcb.embed(pcb.TERMINAL_LIBIDS[0], "J1", "power", 20, 30,
                               angle, {}, {"": 0})
                expected = {(round(20 + x, 6), round(30 + y, 6))
                            for x, y in map(lambda point: transform(*point), local_points)}
                self.assertEqual(set(zone_polygon(F(fp, "zone")[0])), expected)

    def test_terminal_libraries_are_available_in_generated_project(self):
        for count in (2, 3):
            name = f"TerminalBlock_GCT_TBC05-0{count}-1-G-G"
            path = Path(self.out.name) / "phantasm.pretty" / f"{name}.kicad_mod"
            self.assertEqual(str(read(path)[1]), name)


class OverwriteProtectionTests(unittest.TestCase):
    def test_existing_board_is_preserved_without_kicad(self):
        for unplaced in (False, True):
            with self.subTest(unplaced=unplaced), tempfile.TemporaryDirectory() as directory:
                target = Path(directory) / (pcb.UNPLACED_FILE if unplaced else pcb.DRAFT_FILE)
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(b"existing routed board\n")
                with mock.patch.object(pcb, "OUT", directory), \
                        mock.patch.object(pcb, "kicad_cli", side_effect=AssertionError("KiCad called")):
                    with self.assertRaisesRegex(SystemExit, "refusing to overwrite"):
                        pcb.main(unplaced=unplaced)
                self.assertEqual(target.read_bytes(), b"existing routed board\n")


class TerminalFootprintTests(unittest.TestCase):
    def test_stock_library_discovery_reports_the_override_and_keeps_local_loading(self):
        with tempfile.TemporaryDirectory() as directory:
            stock = Path(directory) / "stock"
            library = stock / "Fixture.pretty"
            library.mkdir(parents=True)
            (library / "sample.kicad_mod").write_text('(footprint "sample")', encoding="utf-8")
            local = Path(directory) / "local"
            local.mkdir()
            (local / "custom.kicad_mod").write_text('(footprint "custom")', encoding="utf-8")
            with mock.patch.object(pcb, "_MOD_CACHE", {}), \
                    mock.patch.object(pcb, "FP_DIR", str(stock)), \
                    mock.patch.object(pcb, "LOCAL_FOOTPRINT_DIR", str(local)):
                self.assertEqual(str(pcb.load_mod("Fixture:sample")[1]), "sample")
                with mock.patch.object(pcb, "FP_DIR", str(stock / "missing")):
                    with self.assertRaisesRegex(RuntimeError, "KICAD_FOOTPRINT_DIR") as caught:
                        pcb.load_mod("Fixture:sample")
                    for pattern in sexp.kicad_data_dir_patterns("footprints"):
                        self.assertIn(pattern, str(caught.exception))
                    self.assertEqual(str(pcb.load_mod("phantasm:custom")[1]), "custom")

    def test_revision_13_rejects_legacy_connector_footprints(self):
        token = pcb._GENERATION.set(("1.3", "test"))
        try:
            comps = {ref: (ref, pcb.REVISION_LAYOUTS["1.3"]["terminal_by_ref"][ref], "", False)
                     for ref in pcb.TERMINAL_EDGE_PLACEMENTS_1_3}
            self.assertEqual(pcb.fixed_placements(comps), pcb.TERMINAL_EDGE_PLACEMENTS_1_3)
            comps = {ref: (ref, pcb.QUILTER_FIXED_FOOTPRINTS[ref], "", False)
                     for ref in pcb.TERMINAL_EDGE_PLACEMENTS_1_3}
            self.assertFalse(pcb.TERMINAL_EDGE_PLACEMENTS_1_3.keys() &
                             pcb.fixed_placements(comps).keys())
        finally:
            pcb._GENERATION.reset(token)


@unittest.skipUnless(GENERATES, GENERATES_REASON)
class GeneratedBoardTests(TerminalBodyChecks, unittest.TestCase):
    """The placed draft `pcb.py --force` emits, read back without KiCad."""

    def test_placed_draft_preserves_the_quilter_board_and_rules(self):
        with tempfile.TemporaryDirectory() as directory:
            upload = Path(generate(directory, unplaced=True))
            project = upload.with_suffix(".kicad_pro")
            original = {path: path.read_bytes() for path in (upload, project)}
            draft = Path(generate(directory))
            self.assertNotEqual(draft, upload)
            self.assertEqual(draft.name, "phantasm-draft.kicad_pcb")
            for path, data in original.items():
                self.assertEqual(path.read_bytes(), data)
            self.assertEqual(json.loads(project.read_text())["text_variables"]["PHANTASM_LAYOUT"], "unplaced")
            draft_project = json.loads(draft.with_suffix(".kicad_pro").read_text())
            self.assertEqual(draft_project["text_variables"]["PHANTASM_LAYOUT"], "placed")

    def test_drc_rejects_a_component_under_the_terminal_body(self):
        root = read(self.path)
        footprints = {connectivity.footprint_reference(fp): fp for fp in F(root, "footprint")}
        terminal = footprints["J1"]
        resistor = footprints["R_MEN"]
        x, y = map(float, sexp.val(terminal, "at")[:2])
        at = F(resistor, "at")[0]
        at[1:] = [x + 2.0, y + 1.27, 0]
        board_path = Path(self.out.name) / "blocked.kicad_pcb"
        report_path = Path(self.out.name) / "blocked-drc.json"
        board_path.write_text(sexp.dumps(root), encoding="utf-8")
        result = subprocess.run(
            [kicad_cli(), "pcb", "drc", "--format", "json", "-o",
             str(report_path), str(board_path)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        report = json.loads(report_path.read_text(encoding="utf-8"))
        resistor_uuid = str(sexp.val(resistor, "uuid")[0])
        self.assertTrue(any(
            violation["type"] == "items_not_allowed" and
            any(item["uuid"] == resistor_uuid for item in violation["items"])
            for violation in report["violations"]), report["violations"])

    def test_generated_pair_passes_schematic_parity(self):
        with tempfile.TemporaryDirectory() as directory:
            board_path = generate(directory, unplaced=True)
            subprocess.run(
                [kicad_cli(), "pcb", "drc", "--refill-zones", "--save-board",
                 "-o", str(Path(directory) / "refill.rpt"), board_path],
                capture_output=True, text=True, check=True)
            warnings = {("lib_footprint_mismatch", ref): 1 for ref in ("D_BUS", "J1")}
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
        for name in (pcb.DRAFT_FILE, "fp-lib-table",
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
            ref = connectivity.footprint_reference(footprint)
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
                        if connectivity.footprint_reference(fp) == "R_TX")
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
                mock.patch.object(pcb, "OUT", directory), \
                mock.patch.object(pcb, "SCH", str(COMMITTED_PCB.with_suffix(".kicad_sch"))):
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
        expected = assembly_exclusions(read(COMMITTED_PCB))
        del expected["J4"]
        self.assertEqual(assembly_exclusions(self.root), expected)

    def test_stamps_the_revision_in_the_title_block(self):
        blocks = F(self.root, "title_block")
        self.assertEqual(len(blocks), 1)
        self.assertEqual(str(sexp.val(blocks[0], "rev")[0]), builder.REVISION)

    def test_no_two_courtyards_overlap(self):
        boxes = {}
        for footprint in F(self.root, "footprint"):
            box = courtyard_box(footprint)
            self.assertIsNotNone(box, connectivity.footprint_reference(footprint))
            boxes[connectivity.footprint_reference(footprint)] = box
        refs = sorted(boxes)
        overlaps = [f"{a}/{b}"
                    for index, a in enumerate(refs) for b in refs[index + 1:]
                    if pcb._boxes_overlap(boxes[a], boxes[b])]
        self.assertEqual(overlaps, [])


class TerminalEdgePlacementChecks:
    def locked_front_track_length(self, footprints, ref, pad_number):
        net_id = sexp.val(next(pad for pad in F(footprints[ref], "pad")
                              if pad[1] == pad_number), "net")[0]
        tracks = [track for track in F(self.root, "segment")
                  if sexp.val(track, "net") == [net_id]]
        self.assertTrue(tracks)
        length = 0.0
        for track in tracks:
            self.assertEqual(sexp.val(track, "locked"), ["yes"])
            self.assertEqual(sexp.val(track, "layer"), ["F.Cu"])
            a = tuple(map(float, sexp.val(track, "start")))
            b = tuple(map(float, sexp.val(track, "end")))
            length += math.dist(a, b)
        return length

    def test_vin_bypass_has_a_direct_locked_connection(self):
        footprints = {connectivity.footprint_reference(fp): fp for fp in F(self.root, "footprint")}
        self.assertLess(self.locked_front_track_length(footprints, "C_DEC1", "1"), 3.0)
        groups = connectivity.opens(self.root)["+5V_LOGIC"]
        self.assertTrue(any({("C_DEC1", "1"), ("U_MCU", "VIN")} <= set(group)
                            for group in groups), groups)

    def test_receive_filter_is_locked_and_prerouted_at_d3(self):
        footprints = {connectivity.footprint_reference(fp): fp for fp in F(self.root, "footprint")}
        for ref in ("R1", "R2", "C_SYNC"):
            self.assertEqual(sexp.val(footprints[ref], "locked"), ["yes"], ref)
        self.assertNotIn("FRAME_SYNC", connectivity.opens(self.root))
        self.assertLess(self.locked_front_track_length(footprints, "C_SYNC", "1"), 10.0)

    def test_locked_decouplers_are_within_three_mm_of_supply_pins(self):
        footprints = {connectivity.footprint_reference(fp): fp for fp in F(self.root, "footprint")}

        def pad_center(fp, number):
            pad = next(p for p in F(fp, "pad") if p[1] == number)
            x, y, angle = map(float, sexp.val(fp, "at"))
            dx, dy = connectivity._rotate(tuple(map(float, sexp.val(pad, "at")[:2])), angle)
            return x + dx, y + dy

        for cap, parent, pin in (("C_DEC1", "U_MCU", "VIN"), ("C_DEC2", "U1", "14")):
            with self.subTest(cap=cap):
                self.assertEqual(sexp.val(footprints[cap], "locked"), ["yes"])
                self.assertEqual(sexp.val(footprints[cap], "layer"),
                                 sexp.val(footprints[parent], "layer"))
                self.assertLess(math.dist(pad_center(footprints[cap], "1"),
                                          pad_center(footprints[parent], pin)), 3.0)

    def test_debug_header_and_serial_connection_are_absent(self):
        footprints = {connectivity.footprint_reference(fp): fp for fp in F(self.root, "footprint")}
        self.assertNotIn("J4", footprints)
        nets = {sexp.val(pad, "net")[-1]
                for fp in footprints.values() for pad in F(fp, "pad")
                if sexp.val(pad, "net")}
        self.assertNotIn("/SERIAL1_TX", nets)
        tx = next(pad for pad in F(footprints["U_MCU"], "pad") if pad[1] == "1")
        tx_net = sexp.val(tx, "net")[-1]
        self.assertTrue(tx_net.startswith("unconnected-"), tx_net)
        self.assertEqual(sum(sexp.val(pad, "net", [""])[-1] == tx_net
                             for fp in footprints.values() for pad in F(fp, "pad")), 1)

    def test_connectors_are_locked_inside_the_outline(self):
        footprints = {connectivity.footprint_reference(fp): fp for fp in F(self.root, "footprint")}
        edge = F(self.root, "gr_rect")[0]
        length, width = map(float, sexp.val(edge, "end"))
        self.assertEqual((length, width), (58.28, 32.0))
        row = []
        for ref in ("J1", "J2", "J3A", "J3B"):
            fp = footprints[ref]
            self.assertEqual(sexp.val(fp, "locked"), ["yes"], ref)
            x, y, angle = map(float, sexp.val(fp, "at"))
            self.assertEqual(angle, 0, ref)
            x0, y0, x1, y1 = courtyard_box(fp)
            self.assertGreaterEqual(y0, 0, ref)
            self.assertLessEqual(x1, length, ref)
            self.assertLessEqual(y1, width, ref)
            if ref == "J1":
                self.assertLess(x, 5)
                body = pcb.fp_bbox(fp, graphic_layers=("F.Fab",))
                self.assertGreaterEqual(x + body[0], 0)
                outlines = [line for line in F(fp, "fp_line")
                            if sexp.val(line, "layer") == ["F.SilkS"]]
                self.assertEqual(len(outlines), 3)
                for line in outlines:
                    stroke_width = float(sexp.val(F(line, "stroke")[0], "width")[0])
                    for endpoint in ("start", "end"):
                        self.assertGreaterEqual(
                            x + float(sexp.val(line, endpoint)[0]) - stroke_width / 2,
                            pcb.NEW_LAYOUT_RULES["min_silk_clearance"] - 1e-6)
                module = pcb.load_mod(pcb.REVISION_LAYOUTS["1.2"]["terminal_by_ref"][ref])
                outline = next(rect for rect in F(module, "fp_rect")
                               if sexp.val(rect, "layer") == ["F.SilkS"])
                self.assertEqual(float(sexp.val(outline, "start")[0]), -3.37)
            else:
                self.assertGreaterEqual(x0, 0, ref)
                self.assertLess(length - x, 11)
                row.append((x, y0, y1))
        self.assertEqual(len({entry[0] for entry in row}), 1)
        for first, second in zip(row, row[1:]):
            self.assertGreaterEqual(second[1] - first[2], 0.4)

    def test_mounting_centers_match_the_routed_revision(self):
        routed = {connectivity.footprint_reference(fp): sexp.val(fp, "at")[:2]
                  for fp in F(read(COMMITTED_PCB), "footprint")}
        current = {connectivity.footprint_reference(fp): sexp.val(fp, "at")[:2]
                   for fp in F(self.root, "footprint")}
        for ref in ("H1", "H2", "H3", "H4"):
            self.assertEqual(list(map(float, current[ref])),
                             list(map(float, routed[ref])), ref)

    def test_terminal_body_and_wire_access_keepouts_are_clear(self):
        footprints = {connectivity.footprint_reference(fp): fp for fp in F(self.root, "footprint")}
        boxes = {ref: courtyard_box(fp) for ref, fp in footprints.items()}
        for ref in ("J1", "J2", "J3A", "J3B"):
            zones = {sexp.val(zone, "name")[0]: zone
                     for zone in F(footprints[ref], "zone")}
            self.assertIn("Terminal body and assembly clearance", zones)
            if ref == "J1":
                self.assertIn("Terminal wire access", zones)
            for zone in zones.values():
                self.assertEqual(sexp.val(F(zone, "keepout")[0], "footprints"),
                                 ["not_allowed"])
                points = zone_polygon(zone)
                xs, ys = zip(*points)
                area = (min(xs), min(ys), max(xs), max(ys))
                for other, box in boxes.items():
                    if other != ref and box is not None:
                        self.assertFalse(pcb._boxes_overlap(area, box),
                                         f"{ref} access / {other}")


class CommittedTerminalEdgeTests(TerminalEdgePlacementChecks, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = read(GEN.parent / "1.2" / "phantasm.kicad_pcb")

    def test_front_silk_anchors_clear_mounting_reservations(self):
        for revision in ('1.2', '1.3'):
            board_path = GEN.parent / revision / pcb.PCB_FILE
            root = sexp.parse(board_path.read_text(encoding='utf-8'))[0]
            for node in F(root, 'gr_text'):
                if str(sexp.val(node, 'layer')[0]) != 'F.SilkS':
                    continue
                x, y = (float(v) for v in sexp.val(node, 'at')[:2])
                for ref, (x0, y0, x1, y1) in pcb.mounting_reserve_rects(pcb.QUILTER_LENGTH).items():
                    with self.subTest(revision=revision, label=str(node[1]), hole=ref):
                        self.assertFalse(x0 <= x <= x1 and y0 <= y <= y1)


@unittest.skipUnless(GENERATES, GENERATES_REASON)
class UnplacedBoardTests(TerminalBodyChecks, TerminalEdgePlacementChecks, unittest.TestCase):
    """The autoplacer upload `pcb.py --unplaced` emits."""

    @classmethod
    def setUpClass(cls):
        cls.out = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.out.cleanup)
        cls.root = read(generate(cls.out.name, unplaced=True))

    def test_committed_upload_matches_generator(self):
        for filename in (
            "phantasm.kicad_sch", "phantasm.kicad_pcb", "phantasm.kicad_sym",
            "phantasm.kicad_pro", "phantasm.pretty/Teensy4.0.kicad_mod",
        ):
            with self.subTest(filename=filename):
                self.assertEqual((Path(self.out.name) / filename).read_bytes(),
                                 (GEN.parent / "1.2" / filename).read_bytes())

    def test_locks_the_mechanical_placements(self):
        placed = {}
        comps = {}
        locked = set()
        for footprint in F(self.root, "footprint"):
            ref = connectivity.footprint_reference(footprint)
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
        for ref in (pcb.QUILTER_FIXED.keys() & placed.keys()) - fixed.keys():
            with self.subTest(staged=ref):
                self.assertGreater(placed[ref][1], pcb.PCB_W)

    def test_labels_the_id_straps_on_the_front_silkscreen(self):
        texts = [str(node[1]) for node in F(self.root, "gr_text")
                 if str(sexp.val(node, "layer")[0]) == "F.SilkS"]
        self.assertTrue({"ID0", "ID1", "ID2", "SHLD"}
                        <= set(texts), sorted(texts))
        self.assertTrue({"SYNC IN", "SYNC OUT", "LED OUT"}.isdisjoint(texts))

    def test_labels_terminal_pin_functions(self):
        marks = [(str(node[1]), *(float(v) for v in sexp.val(node, "at")[:2]))
                 for node in F(self.root, "gr_text")
                 if str(sexp.val(node, "layer")[0]) == "F.SilkS"]
        reserves = pcb.mounting_reserve_rects(pcb.QUILTER_LENGTH)
        for ref, labels in (("J2", "DGC"), ("J3A", "SGH"), ("J3B", "SGH")):
            x, y, _ = pcb.TERMINAL_EDGE_PLACEMENTS[ref]
            for pin, label in enumerate(labels):
                row = y + pin * 2.54
                with self.subTest(ref=ref, pin=pin, label=label):
                    hits = [(mx, my) for text, mx, my in marks if text == label]
                    self.assertTrue(hits, sorted(marks))
                    label_x, label_y = min(
                        hits, key=lambda at: math.hypot(at[0] - x, at[1] - row))
                    self.assertAlmostEqual(label_y, row, delta=1e-6)
                    self.assertTrue(x < label_x <= x + 4.5, (label_x, x))
                    for hole, (x0, y0, x1, y1) in reserves.items():
                        self.assertFalse(x0 <= label_x <= x1 and y0 <= label_y <= y1,
                                         f"{ref} {label} inside {hole} reservation")


@unittest.skipUnless(GENERATES, GENERATES_REASON)
class ConnectorEdgePlacementGateTests(unittest.TestCase):
    """Connector placement and footprint gates."""

    def test_missing_connector_placement_is_rejected(self):
        fixed_placements = pcb.fixed_placements

        def without_led(comps):
            fixed = fixed_placements(comps)
            fixed.pop("J2")
            return fixed

        out = self.enterContext(tempfile.TemporaryDirectory())
        with mock.patch.object(pcb, "fixed_placements", without_led), \
                self.assertRaisesRegex(SystemExit, "connectors require verified edge placements: J2"):
            generate(out, unplaced=True)

    def test_legacy_connector_coordinates_are_rejected(self):
        fixed_placements = pcb.fixed_placements

        def legacy_led(comps):
            fixed = fixed_placements(comps)
            fixed["J2"] = pcb.QUILTER_FIXED["J2"]
            return fixed

        out = self.enterContext(tempfile.TemporaryDirectory())
        with mock.patch.object(pcb, "fixed_placements", legacy_led), \
                self.assertRaisesRegex(SystemExit, "connectors require verified edge placements: J2"):
            generate(out, unplaced=True)

    def test_legacy_connector_footprints_are_rejected(self):
        components = pcb.schematic_components

        def legacy_led(schematic):
            return [(ref, pcb.QUILTER_FIXED_FOOTPRINTS[ref] if ref == "J2" else fp, value, dnp)
                    for ref, fp, value, dnp in components(schematic)]

        out = self.enterContext(tempfile.TemporaryDirectory())
        with mock.patch.object(pcb, "schematic_components", legacy_led), \
                self.assertRaisesRegex(SystemExit, "connectors require verified edge placements: J2"):
            generate(out, unplaced=True)

@unittest.skipUnless(GENERATES, GENERATES_REASON)
class OrphanPadTests(unittest.TestCase):
    """A netlist pin must have a matching pad."""

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
