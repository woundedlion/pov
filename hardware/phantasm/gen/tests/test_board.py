"""Self-tests for the schematic generator.

board.py is the only writer of phantasm.kicad_sch and refuses to overwrite the
committed file. check.py and shorts.py read the committed schematic; board,
PCB-generation and revision-generation tests run the generator into a
temporary directory and assert on what it wrote.

Generating needs KiCad's stock symbol libraries (sexp.KICAD_SHARE); the checks
that do not are kept outside that guard and also run against the committed
rev 1.2 schematic, which is generator output. The rev 1.1 schematic carries
KiCad-authored content and is checked by test_shorts and test_builder.
"""
import contextlib
import io
import importlib
import json
import os
import sys
import tempfile
import unittest
import unittest.mock
from pathlib import Path

GEN = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN))

import board    # noqa: E402
import builder  # noqa: E402
import sexp     # noqa: E402
import shorts   # noqa: E402
from constraints import (DEFAULT_CLASS_MINIMUMS, NEW_LAYOUT_RULES,  # noqa: E402
                         RULE_MINIMUMS)
from kicad_common import F  # noqa: E402

STOCK_SYMBOLS = os.path.isdir(sexp.KICAD_SHARE)

PROJECT_FILES = ("phantasm.kicad_sch", "phantasm.kicad_sym", "sym-lib-table",
                 "phantasm.kicad_pro")

# One resistor across a wire span and a net label, so both of its pins land.
LANDED = """(kicad_sch
\t(lib_symbols
\t\t(symbol "Device:R"
\t\t\t(symbol "R_1_1"
\t\t\t\t(pin passive line (at 0 3.81 270) (length 1.27)
\t\t\t\t\t(name "~") (number "1"))
\t\t\t\t(pin passive line (at 0 -3.81 90) (length 1.27)
\t\t\t\t\t(name "~") (number "2")))))
\t(symbol (lib_id "Device:R") (at 100 100 0) (unit 1)
\t\t(property "Reference" "R1") (property "Value" "10k"))
\t(wire (pts (xy 100 96.19) (xy 110 96.19)))
\t(label "NET_A" (at 110 96.19 0))
\t(label "NET_B" (at 100 103.81 0))
)"""

# Same board with pin 2's label dropped: the pin now connects to nothing.
DANGLING = LANDED.replace('\t(label "NET_B" (at 100 103.81 0))\n', "")
WIRE_OVER_PIN = DANGLING.replace(
    '\t(label "NET_A" (at 110 96.19 0))\n',
    '\t(label "NET_A" (at 110 96.19 0))\n'
    '\t(wire (pts (xy 95 103.81) (xy 105 103.81)))\n')


def generate(out):
    """Run the generator into `out`; return the schematic path."""
    sch = os.path.join(out, "phantasm.kicad_sch")
    with unittest.mock.patch.object(board, "OUT", out), \
            unittest.mock.patch.object(board, "SCH", sch), \
            contextlib.redirect_stdout(io.StringIO()):
        board.main(force=True)
    return sch


def dangling_pins(root):
    """[(ref, pin number, point)] for pins with no wiring or no-connect marker.

    Pin coordinates come from the schematic's own lib_symbols, so a stock
    symbol whose pin moved is measured where the generated file placed it.
    """
    libs = {node[1]: builder._index_unit_pins(node)
            for node in sexp.val(root, "lib_symbols", [])
            if isinstance(node, list) and node and node[0] == "symbol"}
    named, wires, junctions = shorts.geometry(root)
    anchors = set(named) | set(junctions)
    anchors.update(shorts.R(tuple(map(float, sexp.val(node, "at"))))
                   for node in F(root, "no_connect"))
    for a, b in wires:
        anchors.add(a)
        anchors.add(b)
    loose = []
    for inst in F(root, "symbol"):
        at = sexp.val(inst, "at")
        mirror = sexp.val(inst, "mirror")
        units = libs[sexp.val(inst, "lib_id")[0]]
        pins = dict(units.get(0, {}))
        pins.update(units.get(int(sexp.val(inst, "unit", [1])[0]), {}))
        for number, pin in pins.items():
            point = shorts.R(builder.transform(
                float(at[0]), float(at[1]),
                float(at[2]) if len(at) > 2 else 0.0,
                mirror[0] if mirror else None, pin["x"], pin["y"]))
            if point in anchors:
                continue
            ref = next((p[2] for p in F(inst, "property")
                        if p[1] == "Reference"), None)
            loose.append((ref, number, point))
    return loose


class BoardEntryPointTests(unittest.TestCase):
    def test_import_ignores_host_arguments(self):
        with unittest.mock.patch.object(sys, "argv", ["host", "--unknown"]):
            importlib.reload(board)

        self.assertTrue(callable(board.main))

    def test_refuses_to_overwrite_an_existing_schematic(self):
        out = self.enterContext(tempfile.TemporaryDirectory())
        sch = Path(out) / "phantasm.kicad_sch"
        sch.write_text("(kicad_sch)", encoding="utf-8")
        with unittest.mock.patch.object(board, "OUT", out), \
                unittest.mock.patch.object(board, "SCH", str(sch)):
            with self.assertRaises(SystemExit) as caught:
                board.main()
        self.assertIn(str(sch), str(caught.exception))
        self.assertEqual(sch.read_text(encoding="utf-8"), "(kicad_sch)")


class ProjectSeedTests(unittest.TestCase):
    def test_uses_all_fabrication_rule_minimums(self):
        project = json.loads(board.project_seed("root-uuid"))
        self.assertEqual(project["board"]["design_settings"]["rules"],
                         {**RULE_MINIMUMS, **NEW_LAYOUT_RULES})
        default = project["net_settings"]["classes"][0]
        for field, minimum in DEFAULT_CLASS_MINIMUMS.items():
            with self.subTest(field=field):
                self.assertEqual(default[field], minimum)

    def test_includes_new_fabrication_rule_minimums(self):
        with unittest.mock.patch.dict(RULE_MINIMUMS, {"min_test_clearance": 0.4}):
            project = json.loads(board.project_seed("root-uuid"))
            self.assertEqual(project["board"]["design_settings"]["rules"],
                             {**RULE_MINIMUMS, **NEW_LAYOUT_RULES})

    def test_links_the_root_sheet(self):
        project = json.loads(board.project_seed("root-uuid"))
        self.assertEqual(project["sheets"], [["root-uuid", "Root"]])


class DanglingPinTests(unittest.TestCase):
    """The landing check must fire on a broken board and stay quiet on a
    connected one, or it proves nothing about the generated schematic."""

    def test_a_connected_pin_is_not_reported(self):
        self.assertEqual(dangling_pins(sexp.parse(LANDED)[0]), [])

    def test_an_unconnected_pin_is_reported(self):
        self.assertEqual(dangling_pins(sexp.parse(DANGLING)[0]),
                         [("R1", "2", (100.0, 103.81))])

    def test_an_explicit_no_connect_is_not_reported(self):
        root = sexp.parse(DANGLING)[0]
        root.append(["no_connect", ["at", "100", "103.81"]])
        self.assertEqual(dangling_pins(root), [])

    def test_a_wire_crossing_a_pin_does_not_connect_it(self):
        self.assertEqual(dangling_pins(sexp.parse(WIRE_OVER_PIN)[0]),
                         [("R1", "2", (100.0, 103.81))])


class BypassConnectionChecks:
    def test_bypass_caps_are_wired_directly_to_their_parent_pins(self):
        libs = {node[1]: builder._index_unit_pins(node)
                for node in sexp.val(self.root, "lib_symbols")}
        positions = {}
        for inst in F(self.root, "symbol"):
            ref = next(p[2] for p in F(inst, "property") if p[1] == "Reference")
            units = libs[sexp.val(inst, "lib_id")[0]]
            pins = {**units.get(0, {}), **units.get(int(sexp.val(inst, "unit")[0]), {})}
            x, y, angle = map(float, sexp.val(inst, "at"))
            mirror = sexp.val(inst, "mirror", [None])[0]
            for number, pin in pins.items():
                positions[ref, number] = shorts.R(builder.transform(
                    x, y, angle, mirror, pin["x"], pin["y"]))
        _, wires, _ = shorts.geometry(self.root)
        connections = {frozenset((a, b)) for a, b in wires}
        for cap, parent, pin in (("C_DEC1", "U_MCU", "VIN"), ("C_DEC2", "U1", "14")):
            self.assertIn(frozenset((positions[cap, "1"], positions[parent, pin])), connections)


class CommittedSchematicTests(BypassConnectionChecks, unittest.TestCase):
    """Both checks read the file's own lib_symbols, so no stock library and no
    kicad-cli is needed: they gate the shipped schematic on every push."""

    @classmethod
    def setUpClass(cls):
        cls.root = sexp.parse(Path(board.SCH).read_text(encoding="utf-8"))[0]

    def test_no_two_named_nets_share_a_group(self):
        conflicts, _ = shorts.analyze(self.root)
        self.assertEqual(conflicts, [])

    def test_every_placed_pin_lands_on_the_wiring(self):
        self.assertEqual(dangling_pins(self.root), [])


@unittest.skipUnless(STOCK_SYMBOLS,
                     f"KiCad stock symbol libraries not found ({sexp.KICAD_SHARE})")
class GeneratedSchematicTests(BypassConnectionChecks, unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.out.cleanup)
        cls.sch = generate(cls.out.name)
        cls.root = sexp.parse(Path(cls.sch).read_text(encoding="utf-8"))[0]

    def test_writes_every_project_file(self):
        for name in PROJECT_FILES:
            self.assertTrue(os.path.exists(os.path.join(self.out.name, name)), name)

    def test_custom_power_descriptions_name_the_generated_net(self):
        symbols = sexp.val(self.root, "lib_symbols")
        for net in ("+5V_RAW", "+5V_LOGIC"):
            with self.subTest(net=net):
                symbol = next(node for node in symbols if node[:2] == ["symbol", f"phantasm:{net}"])
                description = next(prop[2] for prop in F(symbol, "property") if prop[1] == "Description")
                self.assertEqual(description, f'Power symbol creates a global label with name "{net}"')

    def test_teensy_pins_match_footprint_and_unused_pins_are_no_connect(self):
        footprint = GEN.parent / "1.2" / "phantasm.pretty" / "Teensy4.0.kicad_mod"
        pads = F(sexp.parse(footprint.read_text(encoding="utf-8"))[0], "pad")
        pad_numbers = {pad[1] for pad in pads if pad[1]}
        lib = next(node for node in sexp.val(self.root, "lib_symbols")
                   if node[0:2] == ["symbol", "phantasm:Teensy4.0"])
        pins = builder._index_unit_pins(lib)[1]
        self.assertEqual(set(pins), pad_numbers)
        inst = next(node for node in F(self.root, "symbol")
                    if sexp.val(node, "lib_id") == ["phantasm:Teensy4.0"])
        self.assertEqual({node[1] for node in F(inst, "pin")}, pad_numbers)
        x, y, rot = map(float, sexp.val(inst, "at"))
        positions = {number: builder.transform(x, y, rot, None, pin["x"], pin["y"])
                     for number, pin in pins.items()}
        self.assertEqual(len(set(positions.values())), len(pins))
        markers = {tuple(map(float, sexp.val(node, "at")))
                   for node in F(self.root, "no_connect")}
        connected = {"VIN", "3V3", "GND", "3", "4", "5", "11", "13",
                     "21", "22", "23"}
        self.assertEqual(markers, {positions[number] for number in pad_numbers - connected})
        named, wires, junctions = shorts.geometry(self.root)
        self.assertTrue(markers.isdisjoint(set(named) | set(junctions)))
        for start, end in wires:
            for point in markers:
                self.assertFalse(min(start[0], end[0]) <= point[0] <= max(start[0], end[0])
                                 and min(start[1], end[1]) <= point[1] <= max(start[1], end[1]))

    def test_project_uses_fabrication_rule_minimums(self):
        project = Path(self.out.name, "phantasm.kicad_pro")
        rules = json.loads(project.read_text(encoding="utf-8"))[
            "board"]["design_settings"]["rules"]
        for field, minimum in RULE_MINIMUMS.items():
            with self.subTest(field=field):
                self.assertEqual(rules[field], minimum)

    def test_no_two_named_nets_share_a_group(self):
        conflicts, _ = shorts.analyze(self.root)
        self.assertEqual(conflicts, [])

    def test_every_placed_pin_lands_on_the_wiring(self):
        self.assertEqual(dangling_pins(self.root), [])

    def test_places_the_whole_schematic(self):
        refs = {p[2] for inst in F(self.root, "symbol")
                for p in F(inst, "property") if p[1] == "Reference"}
        self.assertTrue({"U_MCU", "U1", "J1", "J2", "J3A", "J3B",
                         "D_BUS", "Q_REV", "F1", "FB"} <= refs, sorted(refs))
        self.assertNotIn("J4", refs)

    def test_d_bus_uses_a_polarized_symbol(self):
        instances = [
            inst for inst in F(self.root, "symbol")
            if any(p[1:3] == ["Reference", "D_BUS"]
                   for p in F(inst, "property"))
        ]
        self.assertEqual(len(instances), 1)
        self.assertEqual(sexp.val(instances[0], "lib_id"), ["Device:D_Zener"])


if __name__ == "__main__":
    unittest.main()
