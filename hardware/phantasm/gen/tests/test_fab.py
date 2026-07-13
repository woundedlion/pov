import sys
import unittest
from pathlib import Path

GEN = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN))

import fab  # noqa: E402


class AssemblyMetadataTests(unittest.TestCase):
    def setUp(self):
        self.ref = "X1"
        self.comps = {self.ref: {"lcsc": " C123 "}}
        self.posrows = {
            self.ref: {
                "PosX": " 1.25 ",
                "PosY": "-2",
                "Rot": "370",
                "Side": " TOP ",
            }
        }

    def validate(self):
        return fab.validate_assembly_metadata(
            self.comps, self.posrows, [self.ref])

    def assert_diagnostic(self, expected):
        with self.assertRaises(fab.AssemblyMetadataError) as caught:
            self.validate()
        self.assertIn(expected, caught.exception.diagnostics)

    def test_accepts_complete_metadata_and_normalizes_text(self):
        metadata = self.validate()

        self.assertEqual(
            metadata[self.ref],
            {
                "lcsc": "C123",
                "pos_x": "1.25",
                "pos_y": "-2",
                "rotation": 370.0,
                "side": "top",
            },
        )

    def test_rejects_missing_supplier_part_number(self):
        del self.comps[self.ref]["lcsc"]

        self.assert_diagnostic(
            f"{self.ref}: supplier part number (LCSC) is missing")

    def test_rejects_blank_supplier_part_number(self):
        self.comps[self.ref]["lcsc"] = "  "

        self.assert_diagnostic(
            f"{self.ref}: supplier part number (LCSC) is blank")

    def test_rejects_missing_centroid_row(self):
        self.posrows = {}

        self.assert_diagnostic(f"{self.ref}: centroid row is missing")

    def test_rejects_missing_centroid_fields(self):
        for field, label in (("PosX", "X"), ("PosY", "Y"),
                             ("Rot", "rotation"), ("Side", "side")):
            with self.subTest(field=field):
                original = self.posrows[self.ref].pop(field)
                self.assert_diagnostic(
                    f"{self.ref}: centroid {label} is missing")
                self.posrows[self.ref][field] = original

    def test_rejects_blank_centroid_fields(self):
        for field, label in (("PosX", "X"), ("PosY", "Y"),
                             ("Rot", "rotation"), ("Side", "side")):
            with self.subTest(field=field):
                original = self.posrows[self.ref][field]
                self.posrows[self.ref][field] = "  "
                self.assert_diagnostic(
                    f"{self.ref}: centroid {label} is blank")
                self.posrows[self.ref][field] = original

    def test_rejects_nonnumeric_centroid_coordinates_and_rotation(self):
        for field, label in (("PosX", "X"), ("PosY", "Y"),
                             ("Rot", "rotation")):
            with self.subTest(field=field):
                original = self.posrows[self.ref][field]
                self.posrows[self.ref][field] = "invalid"
                self.assert_diagnostic(
                    f"{self.ref}: centroid {label} is not numeric: 'invalid'")
                self.posrows[self.ref][field] = original

    def test_rejects_nonfinite_centroid_coordinates_and_rotation(self):
        for field, label, value in (("PosX", "X", "nan"),
                                    ("PosY", "Y", "inf"),
                                    ("Rot", "rotation", "-inf")):
            with self.subTest(field=field):
                original = self.posrows[self.ref][field]
                self.posrows[self.ref][field] = value
                self.assert_diagnostic(
                    f"{self.ref}: centroid {label} is not finite: {value!r}")
                self.posrows[self.ref][field] = original

    def test_rejects_wrong_centroid_side(self):
        self.posrows[self.ref]["Side"] = "bottom"

        self.assert_diagnostic(
            f"{self.ref}: centroid side is 'bottom'; expected 'top'")

    def test_aggregates_diagnostics(self):
        self.comps[self.ref]["lcsc"] = ""
        self.posrows = {}

        with self.assertRaises(fab.AssemblyMetadataError) as caught:
            self.validate()
        self.assertEqual(
            caught.exception.diagnostics,
            (
                f"{self.ref}: supplier part number (LCSC) is blank",
                f"{self.ref}: centroid row is missing",
            ),
        )


if __name__ == "__main__":
    unittest.main()
