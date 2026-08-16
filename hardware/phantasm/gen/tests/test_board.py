"""Self-tests for the schematic generator.

board.py is the only writer of phantasm.kicad_sch and refuses to overwrite the
committed file, so nothing else in the repo executes it: check.py and shorts.py
both read the committed schematic. These tests run the generator into a
temporary directory and assert on what it wrote.

Generating needs KiCad's stock symbol libraries (sexp.KICAD_SHARE); the checks
that do not are kept outside that guard.
"""
import contextlib
import io
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


def generate(out):
    """Run the generator into `out`; return the schematic path."""
    sch = os.path.join(out, "phantasm.kicad_sch")
    with unittest.mock.patch.object(board, "OUT", out), \
            unittest.mock.patch.object(board, "SCH", sch), \
            contextlib.redirect_stdout(io.StringIO()):
        board.main(force=True)
    return sch


def dangling_pins(root):
    """[(ref, pin number, point)] for pins touching no wire, junction or label.

    Pin coordinates come from the schematic's own lib_symbols, so a stock
    symbol whose pin moved is measured where the generated file placed it.
    """
    libs = {node[1]: builder._index_unit_pins(node)
            for node in sexp.val(root, "lib_symbols", [])
            if isinstance(node, list) and node and node[0] == "symbol"}
    named, wires, junctions = shorts.geometry(root)
    anchors = set(named) | set(junctions)
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
            if point in anchors or any(shorts.on_seg(point, a, b)
                                       for a, b in wires):
                continue
            ref = next((p[2] for p in F(inst, "property")
                        if p[1] == "Reference"), None)
            loose.append((ref, number, point))
    return loose


class BoardEntryPointTests(unittest.TestCase):
    def test_import_ignores_host_arguments(self):
        with unittest.mock.patch.object(sys, "argv", ["host", "--unknown"]):
            import board

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


class DanglingPinTests(unittest.TestCase):
    """The landing check must fire on a broken board and stay quiet on a
    connected one, or it proves nothing about the generated schematic."""

    def test_a_connected_pin_is_not_reported(self):
        self.assertEqual(dangling_pins(sexp.parse(LANDED)[0]), [])

    def test_an_unconnected_pin_is_reported(self):
        self.assertEqual(dangling_pins(sexp.parse(DANGLING)[0]),
                         [("R1", "2", (100.0, 103.81))])


@unittest.skipUnless(STOCK_SYMBOLS,
                     f"KiCad stock symbol libraries not found ({sexp.KICAD_SHARE})")
class GeneratedSchematicTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.out = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.out.cleanup)
        cls.sch = generate(cls.out.name)
        cls.root = sexp.parse(Path(cls.sch).read_text(encoding="utf-8"))[0]

    def test_writes_every_project_file(self):
        for name in PROJECT_FILES:
            self.assertTrue(os.path.exists(os.path.join(self.out.name, name)), name)

    def test_no_two_named_nets_share_a_group(self):
        conflicts, _ = shorts.analyze(self.root)
        self.assertEqual(conflicts, [])

    def test_every_placed_pin_lands_on_the_wiring(self):
        self.assertEqual(dangling_pins(self.root), [])

    def test_places_the_whole_schematic(self):
        refs = {p[2] for inst in F(self.root, "symbol")
                for p in F(inst, "property") if p[1] == "Reference"}
        self.assertTrue({"U_MCU", "U1", "J1", "J2", "J3A", "J3B", "J4",
                         "D_BUS", "Q_REV", "F1", "FB"} <= refs, sorted(refs))

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
