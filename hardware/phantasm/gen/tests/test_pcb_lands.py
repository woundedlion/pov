import sys
import unittest
from pathlib import Path

GEN_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN_DIR))

import pcb  # noqa: E402
import sexp  # noqa: E402
from kicad_common import F  # noqa: E402

UNPLACED = GEN_DIR.parent / "unplaced" / "phantasm_unplaced.kicad_pcb"
TEENSY_LIBRARY = GEN_DIR.parent / "phantasm.pretty" / "Teensy4.0.kicad_mod"

# Every reference whose land must come from its footprint library unchanged.
LIBRARY_LAND_REFS = ("R1", "R2", "R_PD", "R_S", "R_MEN", "R_D1", "R_D2",
                     "R_LF", "C_SYNC", "C_DEC1", "C_DEC2", "C_LF")

CHIP_LIBID = "Test:R_chip"
CHIP_MOD = """(footprint "R_chip"
\t(layer "F.Cu")
\t(pad "1" smd roundrect (at -0.825 0) (size 0.8 0.95)
\t\t(layers "F.Cu" "F.Mask" "F.Paste") (uuid "a"))
\t(pad "2" smd roundrect (at 0.825 0) (size 0.8 0.95)
\t\t(layers "F.Cu" "F.Mask" "F.Paste") (uuid "b"))
)"""


def pad_lands(node):
    """Sorted (pad number, at x, at y, size x, size y) tuples — pad numbers repeat
    on multi-pad nets. Pad angle is placement state (the footprint rotation folded
    in), not part of the land pattern."""
    lands = []
    for pad in F(node, "pad"):
        at = sexp._val(pad, "at")
        size = sexp._val(pad, "size")
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


class TeensyLibraryTests(unittest.TestCase):
    """The committed, routed board resolves its Teensy pads against the
    committed library, so generator drift there invalidates the routing."""

    def test_generator_matches_the_committed_library(self):
        library = sexp.parse(TEENSY_LIBRARY.read_text(encoding="utf-8"))[0]

        self.assertEqual(pad_lands(pcb.teensy_footprint()), pad_lands(library))


if __name__ == "__main__":
    unittest.main()
