# Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
# Licensed under the PolyForm Noncommercial License 1.0.0
"""Typed capture fixture generation contracts."""

import copy
import json
from pathlib import Path
import struct
import unittest

import gen_chain_capture_fixtures as generator
from generate_pullback_manifest_header import load_and_validate, protocol_definition
from pullback_capture import operation_specs, oracle_operation_specs


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
        directory = Path(__file__).resolve().parents[1] / "tests/data/pullback"
        programs, oracles, _ = load_and_validate(directory)
        codes = protocol_definition()[1]
        expected = set()
        for width, height in programs["corpus"]["resolutions"]:
            for spec in operation_specs(programs):
                expected.add(("frame", spec["name"], width, height,
                              spec["preset"], codes[spec["mapping"]],
                              struct.pack("<f", 0.0)))
            for spec in oracle_operation_specs(oracles):
                expected.add(("oracle", spec["oracle"], width, height,
                              spec["preset"], codes[spec["mapping"]],
                              struct.pack("<f", spec["hue_noise_phase"])))
        actual = [
            (record["kind"], record["name"], record["width"], record["height"],
             record["preset"], record["operation"],
             struct.pack("<f", record.get("hueNoisePhase", 0.0)))
            for record in self.records
        ]
        self.assertEqual(set(actual), expected)
        self.assertEqual(len(actual), len(expected))
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
