import sys
from committed_boards import committed_project
import unittest
from pathlib import Path

GEN = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN))

import sexp  # noqa: E402
from kicad_common import F  # noqa: E402
from kicad_env import footprints_by_ref  # noqa: E402


class FootprintLookupTests(unittest.TestCase):
    def test_lookup_uses_reference_properties_and_retains_nodes(self):
        resistor = sexp.parse_one('(footprint "R_0603" (property "Reference" "R1"))')
        connector = sexp.parse_one('(footprint "Terminal" (property "Reference" "J1"))')
        root = ['kicad_pcb', ['gr_text', 'ignored', ['property', 'Reference', 'TEXT']],
                resistor, connector]
        found = footprints_by_ref(root)
        self.assertEqual(set(found), {"R1", "J1"})
        self.assertIs(found["R1"], resistor)
        self.assertIs(found["J1"], connector)

    def test_missing_and_duplicate_references_fail_loudly(self):
        for nodes, message in ((['(footprint "R_0603")'], 'footprint has no reference'),
                               (['(footprint "R_0603" (property "Reference" "R1"))',
                                 '(footprint "R_0805" (property "Reference" "R1"))'],
                                'duplicate footprint reference: R1')):
            with self.subTest(message=message), self.assertRaisesRegex(AssertionError, message):
                footprints_by_ref(['kicad_pcb', *(sexp.parse_one(node) for node in nodes)])


class RevisionProjectTests(unittest.TestCase):
    def test_upload_component_pin_numbers_match(self):
        for revision, expected_count in (("1.2", 28), ("1.3", 32)):
            with self.subTest(revision=revision):
                project = committed_project(revision)
                schematic = sexp.parse_one((project / "phantasm.kicad_sch").read_text(encoding="utf-8"))
                board = sexp.parse_one((project / "phantasm.kicad_pcb").read_text(encoding="utf-8"))
                pins = {}
                for symbol in F(schematic, "symbol"):
                    properties = {prop[1]: prop[2] for prop in F(symbol, "property")}
                    if properties.get("Footprint"):
                        pins.setdefault(properties["Reference"], set()).update(
                            str(pin[1]) for pin in F(symbol, "pin"))
                pads = {}
                for reference, footprint in footprints_by_ref(board).items():
                    numbers = {str(pad[1]) for pad in F(footprint, "pad") if pad[1]}
                    if not numbers:
                        continue
                    pads[reference] = numbers
                self.assertEqual(len(pins), expected_count)
                self.assertNotIn("J4", pins)
                self.assertEqual(pins, pads)


if __name__ == "__main__":
    unittest.main()
