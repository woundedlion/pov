import sys
import unittest
from pathlib import Path

GEN = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN))

import sexp  # noqa: E402
from kicad_common import F  # noqa: E402


class RevisionProjectTests(unittest.TestCase):
    def test_each_project_carries_its_directory_revision(self):
        for revision in ("1.1", "1.2"):
            for suffix in ("kicad_sch", "kicad_pcb"):
                with self.subTest(revision=revision, suffix=suffix):
                    path = GEN.parent / revision / f"phantasm.{suffix}"
                    root = sexp.parse_one(path.read_text(encoding="utf-8"))
                    self.assertEqual(sexp.val(F(root, "title_block")[0], "rev"),
                                     [revision])

    def test_upload_component_pin_numbers_match(self):
        project = GEN.parent / "1.2"
        schematic = sexp.parse_one((project / "phantasm.kicad_sch").read_text(encoding="utf-8"))
        board = sexp.parse_one((project / "phantasm.kicad_pcb").read_text(encoding="utf-8"))
        pins = {}
        for symbol in F(schematic, "symbol"):
            properties = {prop[1]: prop[2] for prop in F(symbol, "property")}
            if properties.get("Footprint"):
                pins.setdefault(properties["Reference"], set()).update(
                    str(pin[1]) for pin in F(symbol, "pin"))
        pads = {}
        for footprint in F(board, "footprint"):
            numbers = {str(pad[1]) for pad in F(footprint, "pad") if pad[1]}
            if not numbers:
                continue
            properties = {prop[1]: prop[2] for prop in F(footprint, "property")}
            reference = properties.get("Reference") or next(
                text[2] for text in F(footprint, "fp_text") if text[1] == "reference")
            pads[reference] = numbers
        self.assertEqual(len(pins), 28)
        self.assertNotIn("J4", pins)
        self.assertEqual(pins, pads)


if __name__ == "__main__":
    unittest.main()
