import sys
import tempfile
import unittest
import unittest.mock
from pathlib import Path

GEN_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN_DIR))

import fab  # noqa: E402
import pcb  # noqa: E402
import sexp  # noqa: E402
from kicad_common import F  # noqa: E402

UNPLACED = GEN_DIR.parent / "unplaced" / "phantasm_unplaced.kicad_pcb"
ROUTED = GEN_DIR.parent / "phantasm.kicad_pcb"
TEENSY_LIBRARY = GEN_DIR.parent / "phantasm.pretty" / "Teensy4.0.kicad_mod"

# Every reference the generator embeds straight from its footprint library.
LIBRARY_LAND_REFS = ("R1", "R2", "R_PD", "R_S", "R_MEN", "R_D1", "R_D2",
                     "R_LF", "C_SYNC", "C_DEC1", "C_DEC2", "C_LF")


def chip(width, pitch, height=0.95):
    return (("1", -pitch, 0.0, width, height), ("2", pitch, 0.0, width, height))


# Pad geometry the routed board ships, per reference. R1/R2/R_PD/R_S carry the
# stock footprint id with the pads widened in place, so KiCad reports no parity
# difference and "Update Footprints from Library" would silently restore the
# library land.
SHIPPED_CHIP_LANDS = {
    "R1": ("Resistor_SMD:R_0603_1608Metric", chip(1.2, 0.825)),
    "R2": ("Resistor_SMD:R_0603_1608Metric", chip(1.2, 0.825)),
    "R_PD": ("Resistor_SMD:R_0603_1608Metric", chip(1.2, 0.825)),
    "R_S": ("Resistor_SMD:R_0805_2012Metric", chip(1.4, 0.9125, 1.4)),
    "R_MEN": ("Resistor_SMD:R_0603_1608Metric", chip(0.8, 0.825)),
    "R_D1": ("Resistor_SMD:R_0805_2012Metric", chip(1.025, 0.9125, 1.4)),
    "R_D2": ("Resistor_SMD:R_0805_2012Metric", chip(1.025, 0.9125, 1.4)),
    "R_LF": ("Resistor_SMD:R_0805_2012Metric", chip(1.025, 0.9125, 1.4)),
    "C_SYNC": ("Capacitor_SMD:C_0603_1608Metric", chip(0.9, 0.775)),
    "C_DEC1": ("Capacitor_SMD:C_0603_1608Metric", chip(0.9, 0.775)),
    "C_DEC2": ("Capacitor_SMD:C_0603_1608Metric", chip(0.9, 0.775)),
}

CHIP_LIBID = "Test:R_chip"
CHIP_MOD = """(footprint "R_chip"
\t(layer "F.Cu")
\t(pad "1" smd roundrect (at -0.825 0) (size 0.8 0.95)
\t\t(layers "F.Cu" "F.Mask" "F.Paste") (uuid "a"))
\t(pad "2" smd roundrect (at 0.825 0) (size 0.8 0.95)
\t\t(layers "F.Cu" "F.Mask" "F.Paste") (uuid "b"))
)"""

SCH_TEENSY = """(kicad_sch
\t(symbol (lib_id "phantasm:Teensy4.0")
\t\t(property "Reference" "U_MCU")
\t\t(property "Value" "Teensy4.0")
\t\t(property "Footprint" "")
\t)
)"""

SCH_BLANK_PASSIVE = """(kicad_sch
\t(symbol (lib_id "Device:R")
\t\t(property "Reference" "R1")
\t\t(property "Value" "10k")
\t\t(property "Footprint" "")
\t)
)"""


def pad_lands(node):
    """Sorted (pad number, at x, at y, size x, size y) tuples — pad numbers repeat
    on multi-pad nets. Pad angle is placement state (the footprint rotation folded
    in), not part of the land pattern."""
    lands = []
    for pad in F(node, "pad"):
        at = sexp.val(pad, "at")
        size = sexp.val(pad, "size")
        lands.append((str(pad[1]), float(at[0]), float(at[1]),
                      float(size[0]), float(size[1])))
    return tuple(sorted(lands))


def reference(node):
    for child in node:
        if (isinstance(child, list) and child and child[0] == "property"
                and child[1] == "Reference"):
            return child[2]
    return None


class EmbedLandTests(unittest.TestCase):
    """A per-reference pad override makes two boards carry the same part on
    different lands, so the footprint library is the only land source."""

    def setUp(self):
        pcb._MOD_CACHE[CHIP_LIBID] = sexp.parse(CHIP_MOD)[0]
        self.addCleanup(pcb._MOD_CACHE.pop, CHIP_LIBID, None)

    def test_embed_keeps_the_library_land(self):
        expected = pad_lands(sexp.parse(CHIP_MOD)[0])
        for ref in LIBRARY_LAND_REFS:
            for rot in (0, 90, 180):
                node = pcb.embed(CHIP_LIBID, ref, "10k", 10.0, 20.0, rot,
                                 {}, {"": 0})
                self.assertEqual(pad_lands(node), expected, f"{ref} at {rot}")


class UnplacedBoardLandTests(unittest.TestCase):
    """The unplaced board is the placement/routing input, so it is where a land
    deviation would reach fabrication."""

    def test_one_land_per_footprint_library_id(self):
        root = sexp.parse(UNPLACED.read_text(encoding="utf-8"))[0]
        families = {}
        for node in F(root, "footprint"):
            ref = reference(node)
            self.assertIsNotNone(ref)
            families.setdefault(str(node[1]), {})[ref] = pad_lands(node)
        self.assertTrue(families)
        for libid, lands in sorted(families.items()):
            self.assertEqual(len(set(lands.values())), 1, f"{libid}: {lands}")


class RoutedBoardLandTests(unittest.TestCase):
    """The routed board is the fabrication source of truth, and its chip lands are
    not all the library ones: the four bench-tuned sync resistors carry widened
    pads under a stock footprint id, which no parity or DRC gate can see."""

    def test_routed_board_ships_the_pinned_lands(self):
        root = sexp.parse(ROUTED.read_text(encoding="utf-8"))[0]
        shipped = {}
        for node in F(root, "footprint"):
            ref = reference(node)
            if ref in SHIPPED_CHIP_LANDS:
                shipped[ref] = (str(node[1]), pad_lands(node))
        self.assertEqual(shipped, SHIPPED_CHIP_LANDS)


class PowerInletTests(unittest.TestCase):
    """J1 ships unkeyed on both artifacts while gen/board.py names a keyed JST XA,
    and fab.EXCLUDE_FP_SUBSTR holds both spellings, so the assembly gate cannot
    see the substitution."""

    SHIPPED = "Connector_PinHeader_2.54mm:PinHeader_1x02_P2.54mm_Vertical"
    KEYED = "Connector_JST:JST_XA_B02B-XASK-1-A_1x02_P2.50mm_Vertical"

    def test_schematic_ships_the_unkeyed_header(self):
        footprints = {ref: fp for ref, fp, _, _ in pcb.schematic_components()}
        self.assertEqual(footprints["J1"], self.SHIPPED)

    def test_routed_board_ships_the_unkeyed_header(self):
        root = sexp.parse(ROUTED.read_text(encoding="utf-8"))[0]
        footprints = {reference(node): str(node[1])
                      for node in F(root, "footprint")}
        self.assertEqual(footprints["J1"], self.SHIPPED)

    def test_the_assembly_gate_excludes_both_spellings(self):
        for footprint in (self.SHIPPED, self.KEYED):
            with self.subTest(footprint=footprint):
                self.assertFalse(fab.is_assembled(
                    {"footprint": footprint, "value": "", "dnp": False}))


class TeensyLibraryTests(unittest.TestCase):
    """The committed, routed board resolves its Teensy pads against the
    committed library, so generator drift there invalidates the routing."""

    def test_generator_matches_the_committed_library(self):
        library = sexp.parse(TEENSY_LIBRARY.read_text(encoding="utf-8"))[0]

        self.assertEqual(pad_lands(pcb.teensy_footprint()), pad_lands(library))


class SchematicFootprintTests(unittest.TestCase):
    """A blank Footprint field is the Teensy's marker, so any other symbol that
    loses one would be embedded as the 37x19 mm through-hole Teensy land."""

    def components(self, text):
        path = Path(self.enterContext(tempfile.TemporaryDirectory())) / "s.kicad_sch"
        path.write_text(text, encoding="utf-8")
        self.enterContext(unittest.mock.patch.object(pcb, "SCH", str(path)))
        return pcb.schematic_components()

    def test_blank_footprint_resolves_to_the_teensy(self):
        self.assertEqual(self.components(SCH_TEENSY),
                         [("U_MCU", pcb.TEENSY_LIBID, "Teensy4.0", False)])

    def test_blank_footprint_on_another_symbol_is_rejected(self):
        with self.assertRaises(SystemExit) as caught:
            self.components(SCH_BLANK_PASSIVE)
        self.assertIn("R1", str(caught.exception))
        self.assertIn("Device:R", str(caught.exception))

    def test_committed_schematic_names_every_footprint(self):
        for ref, fp, _, _ in pcb.schematic_components():
            self.assertTrue(fp, ref)


if __name__ == "__main__":
    unittest.main()
