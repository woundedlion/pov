# Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
# Licensed under the PolyForm Noncommercial License 1.0.0
"""Typed capture fixture generation contracts."""

import copy
import json
from pathlib import Path
import struct
import unittest

import gen_chain_capture_fixtures as generator


class ChainCaptureFixtureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = (
            Path(__file__).resolve().parents[1]
            / "tests/data/chain_capture_fixtures.jsonl"
        )
        cls.records = [
            json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()
        ]

    def test_authored_corpus_is_complete_and_emits_typed_state(self):
        self.assertEqual(len(self.records), 488)
        for width, height in ((96, 20), (288, 144)):
            selected = [
                r for r in self.records if (r["width"], r["height"]) == (width, height)
            ]
            self.assertEqual(sum(r["kind"] == "frame" for r in selected), 227)
            self.assertEqual(sum(r["kind"] == "oracle" for r in selected), 17)
        output = generator.generate(self.records)
        self.assertIn("SpatialWalkSnapshot", output)
        self.assertIn("GeneratedPaletteBank::Snapshot", output)
        self.assertIn('"sample.grid.v3"', output)
        self.assertNotIn("reinterpret_cast", output)

    def test_duplicate_case_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "duplicate"):
            generator.generate([self.records[0], self.records[0]])

    def test_wire_and_metadata_refusals(self):
        for field, value in (
            ("operation", 65536),
            ("preset", -1),
            ("width", 0),
            ("kind", "other"),
        ):
            record = copy.deepcopy(self.records[0])
            record[field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                generator.generate([record])
        for version in (0, 1, 3):
            record = copy.deepcopy(self.records[0])
            record["snapshot"]["schemaVersion"] = version
            with self.subTest(version=version), self.assertRaisesRegex(
                ValueError, "version-two"
            ):
                generator.generate([record])

    def test_float_payload_is_binary32_exact(self):
        for bits in (0, 0x80000000, 0x3F49BA5E, 1, 0x7F7FFFFF):
            value = struct.unpack("<f", struct.pack("<I", bits))[0]
            self.assertEqual(generator.float_cpp(value), f"f(0x{bits:08x}u)")
        for value in (float("nan"), float("inf")):
            with self.assertRaisesRegex(ValueError, "finite"):
                generator.float_cpp(value)


if __name__ == "__main__":
    unittest.main()
