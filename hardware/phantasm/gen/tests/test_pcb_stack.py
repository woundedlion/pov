import sys
import unittest
from pathlib import Path

GEN = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN))

import fab  # noqa: E402
import pcb  # noqa: E402
import sexp  # noqa: E402
from kicad_common import F  # noqa: E402


BOARDS = {
    "routed": GEN.parent / "phantasm.kicad_pcb",
    "unplaced": GEN.parent / "unplaced" / "phantasm_unplaced.kicad_pcb",
}


def declared_layer_names(path):
    """Every name a board's `(layers ...)` block declares -- canonical and alias."""
    root = sexp.parse(path.read_text(encoding="utf-8"))[0]
    names = []
    for entry in F(root, "layers")[0][1:]:
        names.extend(str(value) for value in entry[1:]
                     if str(value) not in ("signal", "power", "user", "mixed",
                                           "jumper"))
    return names


def declared_layer_types(path):
    """Canonical layer name -> KiCad layer type from a committed board."""
    root = sexp.parse(path.read_text(encoding="utf-8"))[0]
    return {str(entry[1]): str(entry[2])
            for entry in F(root, "layers")[0][1:]}


def stackup_copper():
    """(name, thickness) of every copper layer the emitted stackup declares."""
    root = sexp.parse("(stackup " + " ".join(pcb.STACKUP) + ")")[0]
    return tuple(
        (str(layer[1]), str(sexp.val(layer, "thickness")[0]))
        for layer in F(root, "layer")
        if sexp.val(layer, "type")[0] == "copper"
    )


class CopperStackTests(unittest.TestCase):
    """A layer the board does not declare plots as an empty Gerber, so the fab
    builds a board whose inner planes and keepouts silently vanish."""

    def test_declares_the_four_layer_signal_ground_stack(self):
        self.assertEqual(pcb.copper_layer_names(),
                         ("F.Cu", "In1.Cu", "In2.Cu", "B.Cu"))
        types = {name: kind for _, name, kind in pcb.COPPER_LAYERS}
        self.assertEqual(types["F.Cu"], "signal")
        self.assertEqual(types["B.Cu"], "signal")
        for name in pcb.GROUND_PLANE_LAYERS:
            self.assertEqual(types[name], "power")

    def test_layer_numbers_and_names_are_unique(self):
        layers = pcb.COPPER_LAYERS + pcb.TECHNICAL_LAYERS
        self.assertEqual(len({number for number, _, _ in layers}), len(layers))
        self.assertEqual(len({name for _, name, _ in layers}), len(layers))

    def test_stackup_copper_matches_the_declared_layers(self):
        self.assertEqual(tuple(name for name, _ in stackup_copper()),
                         pcb.copper_layer_names())

    def test_outer_copper_is_thicker_than_inner(self):
        thickness = dict(stackup_copper())
        for outer in ("F.Cu", "B.Cu"):
            for inner in pcb.GROUND_PLANE_LAYERS:
                self.assertGreater(float(thickness[outer]),
                                   float(thickness[inner]))

    def test_ground_planes_sit_on_declared_copper(self):
        for name in pcb.GROUND_PLANE_LAYERS:
            self.assertIn(name, pcb.copper_layer_names())

    def test_committed_boards_declare_inner_planes_as_power(self):
        for board, path in BOARDS.items():
            types = declared_layer_types(path)
            with self.subTest(board=board):
                self.assertEqual(types["F.Cu"], "signal")
                self.assertEqual(types["B.Cu"], "signal")
                for name in pcb.GROUND_PLANE_LAYERS:
                    self.assertEqual(types[name], "power")

    def test_every_plotted_copper_layer_is_declared(self):
        plotted = [name for name in fab.GERBER_LAYERS.split(",")
                   if name.endswith(".Cu")]

        self.assertEqual(tuple(plotted), pcb.copper_layer_names())


class LayerNameTests(unittest.TestCase):
    """KiCad builds each Gerber's filename from the layer's name, so a name
    carrying a space ships an upload zip the fab has to guess at."""

    def test_no_declared_layer_name_carries_whitespace(self):
        for board, path in BOARDS.items():
            for name in declared_layer_names(path):
                with self.subTest(board=board, name=name):
                    self.assertEqual(name.split(), [name])


if __name__ == "__main__":
    unittest.main()
