"""Self-tests for pullback manifest generation and capture backends."""

from __future__ import annotations

import copy
import contextlib
import io
import json
import re
import shutil
import struct
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


TOOLS = Path(__file__).resolve().parents[1]
ROOT = TOOLS.parent
sys.path.insert(0, str(TOOLS))

import generate_pullback_manifest_header as generator  # noqa: E402
import pullback_capture as capture  # noqa: E402


MANIFEST_DIR = ROOT / "tests/data/pullback"
_, ORACLES, _ = generator.load_and_validate(MANIFEST_DIR)


def _backend_stream(programs):
    specs = capture.operation_specs(programs)
    header = b"HSPB" + struct.pack("<HHHIH", 3, 16, 16, len(specs), len(ORACLES))
    records = []
    for spec in specs:
        name = spec["name"].encode()
        mapping = spec["mapping"]
        selected = 256 if (mapping.startswith("CASE_") or mapping.startswith("THROUGH_CLEAR_")
                           or mapping == "FULL_FRAME") else sum(
            capture._selected_pixel(mapping, x, y, 16, 16)
            for y in range(16) for x in range(16))
        records.append(
            struct.pack("<HHH", spec["preset"],
                        generator.protocol_definition()[1][spec["mapping"]], len(name))
            + name + struct.pack("<HHHHI", 0, 1, 0, 60, selected)
            + struct.pack("<HHH", spec["preset"], 3, 4) * 256
        )
    metrics = []
    for oracle in ORACLES:
        name = oracle["oracle_id"].encode()
        metrics.append(struct.pack("<H", len(name)) + name + struct.pack("<HI", 0, 3))
    return header, records, metrics


def _test_manifest():
    programs, _, _ = generator.load_and_validate(MANIFEST_DIR)
    programs = copy.deepcopy(programs)
    programs["corpus"]["resolutions"] = [[16, 16], [32, 24]]
    return programs



@contextlib.contextmanager
def _staged_manifest():
    with tempfile.TemporaryDirectory() as directory:
        staged = Path(directory)
        for source in MANIFEST_DIR.glob("*.json"):
            shutil.copy(source, staged / source.name)
        yield staged


def _replace_framebuffer_maximum(manifest_dir):
    for path in manifest_dir.glob("*.json"):
        if path.name in {"schema.json", "programs.json"}:
            continue
        oracle = generator._load(path)
        for metric in oracle["metrics"]:
            if (metric["domain"], metric["aggregation"]) == (
                "FRAMEBUFFER", "MAXIMUM"
            ):
                metric["aggregation"] = "MEAN"
        path.write_text(json.dumps(oracle), encoding="utf-8")


class DuplicateJsonKeys(unittest.TestCase):
    def test_duplicate_keys_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "duplicate.json"
            for content in ('{"key": 1, "key": 2}', '{"nested": {"key": 1, "key": 2}}'):
                with self.subTest(content=content):
                    path.write_text(content, encoding="utf-8")
                    with self.assertRaisesRegex(generator.ManifestError, "duplicate JSON key"):
                        generator._load(path)


class ManifestValidation(unittest.TestCase):
    def test_invalid_protocol_definition_uses_the_cli_diagnostic(self):
        for content in (None, 'invalid protocol'):
            with self.subTest(content=content), tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / 'operations.def'
                if content is not None:
                    path.write_text(content, encoding='utf-8')
                generator.protocol_definition.cache_clear()
                stderr = io.StringIO()
                try:
                    with mock.patch.object(generator, 'PROTOCOL_DEFINITION', path), \
                            contextlib.redirect_stderr(stderr), \
                            self.assertRaises(SystemExit) as raised:
                        generator.main(['--manifest-dir', str(MANIFEST_DIR), '--validate-only'])
                    self.assertEqual(raised.exception.code, 2)
                    self.assertIn(str(path), stderr.getvalue())
                    self.assertNotIn('Traceback', stderr.getvalue())
                finally:
                    generator.protocol_definition.cache_clear()

    def test_generator_main_returns_success_for_validation(self):
        self.assertEqual(generator.main([
            "--manifest-dir", str(MANIFEST_DIR), "--validate-only"
        ]), 0)

    def test_header_rows_encode_each_programs_presets(self):
        programs, oracles, _ = generator.load_and_validate(MANIFEST_DIR)
        header = generator.generate_header(programs, oracles)
        rows = {
            program_id: {bit for bit in range(int(mask, 16).bit_length())
                         if int(mask, 16) >> bit & 1}
            for program_id, mask in re.findall(
                r'^    \{"([A-Z0-9_]+)", 0x([0-9a-f]+), \{\{', header, re.MULTILINE)
        }
        expected = {program["id"]: set(program["presets"])
                    for program in programs["programs"]}
        self.assertEqual(rows, expected)
        self.assertTrue(any(len(presets) > 1 for presets in expected.values()))
        preset_count, operations = generator.protocol_definition()
        self.assertEqual(
            re.findall(r"PRESET_COUNT = (\d+);", header), [str(preset_count)])
        self.assertEqual(set().union(*expected.values()), set(range(preset_count)))
        self.assertIn("FULL_FRAME", operations)

    def test_validate_only_runs_header_generation(self):
        programs, oracles, schema = generator.load_and_validate(MANIFEST_DIR)
        broken = copy.deepcopy(oracles)
        broken[0]["metrics"][0].pop("aggregation")
        stderr = io.StringIO()
        with (
            mock.patch.object(
                generator,
                "load_and_validate",
                return_value=(programs, broken, schema),
            ),
            mock.patch.object(
                sys,
                "argv",
                ["generate_pullback_manifest_header.py", "--validate-only"],
            ),
            contextlib.redirect_stderr(stderr),
            self.assertRaises(SystemExit) as caught,
        ):
            generator.main()
        self.assertEqual(caught.exception.code, 2)
        self.assertIn("header generation failed", stderr.getvalue())
        self.assertIn("aggregation", stderr.getvalue())

    def test_program_completeness_is_checked_through_loader(self):
        programs, _, _ = generator.load_and_validate(MANIFEST_DIR)
        broken = copy.deepcopy(programs)
        broken["programs"] = broken["programs"][:1]
        broken["programs"][0]["presets"] = [0]
        with _staged_manifest() as manifest_dir:
            path = manifest_dir / "programs.json"
            path.write_text(json.dumps(broken), encoding="utf-8")
            with self.assertRaisesRegex(
                generator.ManifestError, "presets must be covered exactly once"
            ):
                generator.load_and_validate(manifest_dir)

    def test_new_oracle_is_discovered_and_validated(self):
        with _staged_manifest() as manifest_dir:
            path = manifest_dir / "new_oracle.json"
            path.write_text(json.dumps({"kind": "pullback-oracle"}),
                            encoding="utf-8")
            with self.assertRaisesRegex(generator.ManifestError,
                                        "new_oracle.json"):
                generator.load_and_validate(manifest_dir)

    def test_oracle_set_must_be_nonempty(self):
        with _staged_manifest() as manifest_dir:
            for path in manifest_dir.glob("*.json"):
                if path.name not in {"schema.json", "programs.json"}:
                    path.unlink()
            with self.assertRaisesRegex(generator.ManifestError,
                                        "at least one oracle manifest"):
                generator.load_and_validate(manifest_dir)

    def test_nested_schema_rejects_unknown_fields(self):
        programs, _, _ = generator.load_and_validate(MANIFEST_DIR)
        schema = generator._load(MANIFEST_DIR / "schema.json")
        broken = copy.deepcopy(programs)
        broken["programs"][0]["tolerances"]["unlicensed"] = 1
        with self.assertRaisesRegex(
            generator.ManifestError, "unexpected field unlicensed"
        ):
            generator._validate_schema(broken, schema, schema, "programs")

    def test_schema_shape_rejects_missing_definitions(self):
        schema_path = MANIFEST_DIR / "schema.json"
        schema = generator._load(schema_path)
        schema["$defs"].pop("program")
        with self.assertRaisesRegex(generator.ManifestError,
                                    "missing schema definitions"):
            generator._validate_schema_shape(schema, schema_path)

    def test_schema_shape_rejects_unsupported_keywords(self):
        schema_path = MANIFEST_DIR / "schema.json"
        schema = generator._load(schema_path)
        schema["$defs"]["program"]["properties"]["id"]["format"] = "uuid"
        with self.assertRaisesRegex(generator.ManifestError,
                                    "unsupported schema keywords.*format"):
            generator._validate_schema_shape(schema, schema_path)

    def test_schema_metric_enums_follow_the_approximation_contract(self):
        schema = generator._load(MANIFEST_DIR / "schema.json")
        contract = (ROOT / "core/render/pullback/contract.h").read_text(encoding="utf-8")
        metric = schema["$defs"]["metric"]["properties"]
        for field, enum in (("domain", "ApproximationDomain"),
                            ("aggregation", "ApproximationAggregation")):
            with self.subTest(enum=enum):
                body = re.search(rf"enum class {enum} : uint8_t \{{([^}}]*)\}};", contract)
                self.assertIsNotNone(body, f"{enum} not found in contract.h")
                members = re.sub(r"//[^\n]*|/\*.*?\*/", "", body[1], flags=re.S)
                names = [name.strip() for name in members.split(",") if name.strip()]
                self.assertEqual(metric[field]["enum"], names)

    def test_schema_operations_follow_the_protocol_table(self):
        schema_path = MANIFEST_DIR / "schema.json"
        schema = generator._load(schema_path)
        generator._validate_schema_shape(schema, schema_path)
        preset_count, operations = generator.protocol_definition()
        for name in ("COLUMN_ZERO", "FULL_FRAME"):
            changed = {key: value for key, value in operations.items() if key != name}
            with self.subTest(removed=name), mock.patch.object(
                    generator, "protocol_definition", return_value=(preset_count, changed)):
                with self.assertRaisesRegex(generator.ManifestError,
                                            "must enumerate protocol operations"):
                    generator._validate_schema_shape(schema, schema_path)

    def test_schema_const_rejects_boolean_versions(self):
        programs, _, _ = generator.load_and_validate(MANIFEST_DIR)
        schema = generator._load(MANIFEST_DIR / "schema.json")
        for field_path in (("schema_version",), ("corpus", "version")):
            broken = copy.deepcopy(programs)
            target = broken
            for field in field_path[:-1]:
                target = target[field]
            target[field_path[-1]] = True
            with self.subTest(field=".".join(field_path)):
                with self.assertRaises(generator.ManifestError):
                    generator._validate_schema(broken, schema, schema, "programs")
        for keyword in ("const", "enum"):
            rule = {keyword: 1 if keyword == "const" else [1]}
            with self.subTest(keyword=keyword):
                with self.assertRaises(generator.ManifestError):
                    generator._validate_schema(True, rule, {"$defs": {}}, "numeric")

    def test_non_json_numeric_constants_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "invalid.json"
            for token in ("Infinity", "NaN"):
                with self.subTest(token=token):
                    path.write_text(f'{{"value": {token}}}', encoding="utf-8")
                    with self.assertRaisesRegex(
                        generator.ManifestError, f"non-JSON numeric constant {token}"
                    ):
                        generator._load(path)

    def test_schema_shape_is_asserted(self):
        path = MANIFEST_DIR / "schema.json"
        schema = generator._load(path)
        generator._validate_schema_shape(schema, path)
        gutted = {"$schema": schema["$schema"], "$defs": {}}
        with self.assertRaisesRegex(
            generator.ManifestError, "missing schema definitions"
        ):
            generator._validate_schema_shape(gutted, path)
        mutations = {
            "root selector": lambda doc: doc.pop("oneOf"),
            "closed objects": lambda doc: doc["$defs"]["program"].pop(
                "additionalProperties"
            ),
            "required fields": lambda doc: doc["$defs"]["programManifest"].pop(
                "required"
            ),
            "properties": lambda doc: doc["$defs"]["metric"].pop("properties"),
            "topology fields": lambda doc: doc["$defs"]["topologyKey"][
                "required"
            ].pop(),
        }
        for name, mutate in mutations.items():
            broken = copy.deepcopy(schema)
            mutate(broken)
            with self.subTest(dropped=name):
                with self.assertRaises(generator.ManifestError):
                    generator._validate_schema_shape(broken, path)

    def test_gutted_schema_directory_is_refused(self):
        with _staged_manifest() as staged:
            (staged / "schema.json").write_text(
                json.dumps({"$schema": generator.SCHEMA_DRAFT, "$defs": {}}),
                encoding="utf-8",
            )
            with self.assertRaises(generator.ManifestError):
                generator.load_and_validate(staged)

    def test_oracle_aggregation_is_checked_through_loader(self):
        _, oracles, _ = generator.load_and_validate(MANIFEST_DIR)
        broken = copy.deepcopy(next(
            oracle for oracle in oracles
            if oracle["oracle_id"] == "HUE_ROTATION_AND_NOISE_LUTS"
        ))
        framebuffer = next(
            metric for metric in broken["metrics"] if metric["domain"] == "FRAMEBUFFER"
        )
        framebuffer["configuration_baselines"]["native-debug"] += 2
        with _staged_manifest() as manifest_dir:
            path = manifest_dir / "hue_rotation_noise_luts.json"
            path.write_text(json.dumps(broken), encoding="utf-8")
            with self.assertRaisesRegex(generator.ManifestError, "aggregate"):
                generator.load_and_validate(manifest_dir)

    def test_framebuffer_maximum_is_required(self):
        with _staged_manifest() as manifest_dir:
            _replace_framebuffer_maximum(manifest_dir)
            with self.assertRaisesRegex(generator.ManifestError,
                                        "FRAMEBUFFER/MAXIMUM"):
                generator.load_and_validate(manifest_dir)

    def test_capture_reports_missing_framebuffer_maximum(self):
        with _staged_manifest() as manifest_dir:
            _replace_framebuffer_maximum(manifest_dir)
            stderr = io.StringIO()
            with (
                mock.patch.object(sys, "argv", [
                    "pullback_capture.py",
                    "--manifest-dir", str(manifest_dir),
                    "--operations-output", str(manifest_dir / "operations.txt"),
                ]),
                contextlib.redirect_stderr(stderr),
            ):
                code = capture.main()
            # rc 1, not argparse's rc 2, and no usage dump over the diagnostic:
            # the manifest is wrong, the command line is not.
            self.assertEqual(code, 1)
            self.assertIn("FRAMEBUFFER/MAXIMUM", stderr.getvalue())
            self.assertNotIn("usage:", stderr.getvalue())

    def test_capture_still_reports_a_bad_command_line_as_a_usage_error(self):
        with _staged_manifest() as manifest_dir:
            stderr = io.StringIO()
            with (
                mock.patch.object(sys, "argv", [
                    "pullback_capture.py",
                    "--manifest-dir", str(manifest_dir),
                    "--backend-audit", str(manifest_dir / "audit.bin"),
                ]),
                contextlib.redirect_stderr(stderr),
                self.assertRaises(SystemExit) as caught,
            ):
                capture.main()
            self.assertEqual(caught.exception.code, 2)
            self.assertIn("requires --configuration", stderr.getvalue())



class CaptureBackend(unittest.TestCase):
    def test_backend_oracle_resolution_drift_within_limit_is_accepted(self):
        metrics = {}
        for oracle in ORACLES:
            framebuffer = next(
                metric
                for metric in oracle["metrics"]
                if metric["domain"] == "FRAMEBUFFER"
            )
            metrics[oracle["oracle_id"]] = {
                "value": framebuffer["resolution_baselines"]["native-debug"][
                    "96x20"
                ],
            }
        metrics["HUE_ROTATION_AND_NOISE_LUTS"]["value"] += 1
        capture.validate_backend_oracles(
            metrics, ORACLES, "native-debug", (96, 20)
        )

    def test_backend_oracle_limit_violation_is_refused(self):
        metrics = {
            oracle["oracle_id"]: {
                "value": next(
                    metric
                    for metric in oracle["metrics"]
                    if metric["domain"] == "FRAMEBUFFER"
                )["accepted_limit"]
                + 1,
            }
            for oracle in ORACLES
        }
        with self.assertRaisesRegex(capture.CaptureError, "accepted limit"):
            capture.validate_backend_oracles(
                metrics, ORACLES, "native-debug", (96, 20)
            )


    def test_backend_rejects_malformed_streams(self):
        programs = _test_manifest()
        specs = capture.operation_specs(programs)
        header, records, metrics = _backend_stream(programs)
        body = b"".join(records)
        tail = b"".join(metrics)
        valid = header + body + tail
        unknown_oracle = b"unknown"
        cases = {
            "magic": b"FAIL" + valid[4:],
            "truncated_header": valid[:5],
            "version": valid[:4] + struct.pack("<H", 2) + valid[6:],
            "record_count": valid[:10] + struct.pack("<I", len(specs) - 1) + valid[14:],
            "truncated_name": header + records[0][:7],
            "invalid_name": header + records[0][:6] + b"\xff" + records[0][7:] + b"".join(records[1:]) + tail,
            "duplicate_record": header + records[0] + records[0] + b"".join(records[2:]) + tail,
            "unknown_record": header + struct.pack("<H", 65535) + body[2:] + tail,
            "operation_code": header + body[:2] + struct.pack("<H", 65535) + body[4:] + tail,
            "truncated_pixels": header + records[0][:-1],
            "missing_record": header + b"".join(records[:-1]) + tail,
            "duplicate_oracle": header + body + metrics[0] + metrics[0] + b"".join(metrics[2:]),
            "unknown_oracle": header + body + struct.pack("<H", len(unknown_oracle)) + unknown_oracle + struct.pack("<HI", 0, 3) + b"".join(metrics[1:]),
            "zero_samples": header + body + metrics[0][:-4] + struct.pack("<I", 0) + b"".join(metrics[1:]),
            "truncated_oracle_name": header + body + metrics[0][:3],
            "invalid_oracle_name": header + body + metrics[0][:2] + b"\xff" + metrics[0][3:] + b"".join(metrics[1:]),
            "missing_oracle": header[:-2] + struct.pack("<H", len(metrics) - 1) + body + b"".join(metrics[:-1]),
            "trailing_bytes": valid + b"trailing",
        }
        expected_errors = {
            "magic": "capture backend magic mismatch",
            "truncated_header": "capture backend stream is truncated",
            "version": "capture backend header mismatch",
            "record_count": "capture backend header mismatch",
            "truncated_name": "capture backend operation name is truncated",
            "invalid_name": "capture backend operation name is invalid",
            "duplicate_record": "invalid capture backend record",
            "unknown_record": "invalid capture backend record",
            "operation_code": "capture backend operation mismatch",
            "truncated_pixels": "capture backend stream is truncated",
            "missing_record": "capture backend operation name is truncated",
            "duplicate_oracle": "invalid capture backend oracle metric",
            "unknown_oracle": "invalid capture backend oracle metric",
            "zero_samples": "invalid capture backend oracle metric",
            "truncated_oracle_name": "capture backend oracle name is truncated",
            "invalid_oracle_name": "capture backend oracle name is invalid",
            "missing_oracle": "capture backend stream is incomplete or has trailing bytes",
            "trailing_bytes": "capture backend stream is incomplete or has trailing bytes",
        }
        self.assertEqual(cases.keys(), expected_errors.keys())
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "capture.bin"
            path.write_bytes(valid)
            resolution, loaded_metrics = capture.load_backend(
                path, programs, ORACLES
            )
            self.assertEqual(resolution, (16, 16))
            self.assertEqual(set(loaded_metrics), {oracle["oracle_id"] for oracle in ORACLES})
            for name, data in cases.items():
                with self.subTest(name=name):
                    path.write_bytes(data)
                    with self.assertRaisesRegex(
                        capture.CaptureError, expected_errors[name]
                    ):
                        capture.load_backend(path, programs, ORACLES)

    def test_backend_rejects_drifted_spatial_selection(self):
        programs = _test_manifest()
        specs = capture.operation_specs(programs)
        header, records, metrics = _backend_stream(programs)
        index = next(i for i, spec in enumerate(specs)
                     if spec["mapping"] == "OCTANT_COLUMNS")
        count_at = 6 + len(specs[index]["name"].encode()) + 8
        record = records[index]
        records[index] = record[:count_at] + struct.pack("<I", 0) + record[count_at + 4:]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "capture.bin"
            path.write_bytes(header + b"".join(records) + b"".join(metrics))
            with self.assertRaisesRegex(capture.CaptureError,
                                        "capture backend selected-pixel mismatch"):
                capture.load_backend(path, programs, ORACLES)

if __name__ == "__main__":
    unittest.main()
