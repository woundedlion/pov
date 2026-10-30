# Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
# Licensed under the PolyForm Noncommercial License 1.0.0
"""Typed capture fixture generation contracts."""

import copy
import json
from pathlib import Path
import struct
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import gen_chain_capture_fixtures as generator  # noqa: E402
from generate_pullback_manifest_header import load_and_validate, protocol_definition  # noqa: E402
from pullback_capture import operation_specs, oracle_operation_specs  # noqa: E402


class ChainCaptureFixtureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = (
            Path(__file__).resolve().parents[2]
            / "tests/data/chain_capture_fixtures.jsonl"
        )
        cls.records = [
            json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()
        ]

    def test_authored_corpus_is_complete_and_emits_typed_state(self):
        directory = Path(__file__).resolve().parents[2] / "tests/data/pullback"
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
        self.assertIn("ColorModulation::Config", output)
        self.assertIn("Pullback::Modulation::State", output)
        self.assertIn("Pullback::Interp::ModulationClockSnapshot", output)
        self.assertNotIn("reinterpret_cast", output)

    def test_corpus_uses_canonical_recipes_and_typed_clocks(self):
        for record in self.records:
            snapshot = record["snapshot"]
            with self.subTest(name=record["name"]):
                self.assertEqual(snapshot["schemaVersion"], 3)
                recipe = snapshot["modulation"][0]["recipe"]
                clock = next(entry for entry in snapshot["runtime"]
                             if entry["kind"] == "color-clock-v2")
                self.assertEqual(recipe["seed"], 6047)
                self.assertEqual(clock["state"]["modulation"]["source"], recipe["source"])
                self.assertEqual(clock["state"]["modulation"]["version"], recipe["version"])
                self.assertNotIn("hueNoiseSeed", clock["state"])
                self.assertNotIn("hueNoisePhase", clock["state"])
                fields = {entry["name"]: entry["value"] for entry in snapshot["parameters"]}
                self.assertNotIn("colorize.hue-noise-scale", fields)
                self.assertNotIn("colorize.hue-noise-speed", fields)
                path = fields["colorize.hue-shift-mode"] == 1
                self.assertEqual(recipe["enabled"], not path)
                if not path:
                    self.assertEqual(fields["colorize.hue-shift-amount"], 0)

    def test_duplicate_case_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "duplicate"):
            generator.generate([self.records[0], self.records[0]])

    def test_frozen_noise_controls_keep_their_authored_values(self):
        snapshot = self.records[0]["snapshot"]
        recipe = snapshot["modulation"][0]["recipe"]
        self.assertEqual(recipe["noise"]["scale"], 0.8300312757492065)
        self.assertEqual(recipe["noise"]["speed"], 0.00021200001356191933)
        self.assertEqual(recipe["mechanic"]["amount"], 0.39800000190734863)
        clock = next(entry for entry in snapshot["runtime"]
                     if entry["kind"] == "color-clock-v2")
        self.assertEqual(clock["state"]["modulation"]["phaseTurns"], 0.78100002)

    def test_legacy_snapshots_split_recipes_and_clocks_without_mutating_the_corpus(self):
        original = self.records[0]["snapshot"]
        before = copy.deepcopy(original)
        snapshot = generator.migrate_legacy_snapshot(original)
        self.assertEqual(original, before)
        self.assertEqual(snapshot["schemaVersion"], 4)
        self.assertEqual(snapshot["chain"][-3:], [
            {"instance": "colorize", "operator": "colorize.generated-palette.v4"},
            {"instance": "displacement", "operator": "color.modulate.v1"},
            {"instance": "modulation", "operator": "color.modulate.v1"},
        ])
        path, spatial = snapshot["modulation"]
        self.assertEqual(path["recipe"]["source"], 4)
        self.assertEqual(path["recipe"]["evaluation"], 0)
        self.assertFalse(path["recipe"]["enabled"])
        self.assertEqual(spatial["recipe"], original["modulation"][0]["recipe"])
        clocks = {entry["instance"]: entry for entry in snapshot["runtime"]}
        self.assertEqual(clocks["modulation"]["state"],
                         next(entry for entry in original["runtime"]
                              if entry["instance"] == "colorize")["state"]["modulation"])
        self.assertEqual(clocks["displacement"]["state"]["phaseTurns"], 0)
        self.assertNotIn("modulation", clocks["colorize"]["state"])
        self.assertTrue(all("hue-shift" not in entry["name"] for entry in snapshot["parameters"]))
        self.assertIs(generator.migrate_legacy_snapshot(snapshot), snapshot)
        self.assertIn("color.modulate.v1", generator.generate([{**self.records[0], "snapshot": snapshot}]))

    def test_unknown_pattern_and_retired_runtime_are_rejected(self):
        self.assertIn("ColorModulation::Config", generator.generate([self.records[0]]))
        record = copy.deepcopy(self.records[0])
        record["snapshot"]["modulation"][0]["recipe"]["source"] = 5
        with self.assertRaisesRegex(ValueError, "invalid Pattern::SourceKind"):
            generator.generate([record])
        record = copy.deepcopy(self.records[0])
        clock = next(entry for entry in record["snapshot"]["runtime"]
                     if entry["kind"] == "color-clock-v2")
        clock["kind"] = "color-clock-v1"
        with self.assertRaisesRegex(ValueError, "unknown capture runtime"):
            generator.generate([record])

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
        for version in (0, 1, 2, 5):
            record = copy.deepcopy(self.records[0])
            record["snapshot"]["schemaVersion"] = version
            with self.subTest(version=version), self.assertRaisesRegex(
                ValueError, "version-four"
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
