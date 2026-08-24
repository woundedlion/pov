"""Self-tests for pullback manifest generation and capture comparison."""

from __future__ import annotations

import copy
import contextlib
import hashlib
import io
import json
import shutil
import struct
import sys
import tempfile
import unittest
import weakref
from pathlib import Path
from unittest import mock


TOOLS = Path(__file__).resolve().parents[1]
ROOT = TOOLS.parent
sys.path.insert(0, str(TOOLS))

import generate_pullback_manifest_header as generator  # noqa: E402
import pullback_crosscheck as crosscheck  # noqa: E402
import pullback_capture as capture  # noqa: E402


MANIFEST_DIR = ROOT / "tests/data/pullback"
_, ORACLES, SCHEMA = generator.load_and_validate(MANIFEST_DIR)


def _frame(programs, program, preset, case, resolution, probe, value=1):
    total = resolution[0] * resolution[1]
    preset_count = sum(len(entry["presets"])
                       for entry in programs["programs"])
    pixels = [[value, value, value, 65535] for _ in range(total)]
    if probe == "steady":
        operation = {"kind": "parameter-case", "selected_pixels": total}
    else:
        mapping = programs["corpus"]["probe_operations"][probe]
        if mapping in ("THROUGH_CLEAR_FROM", "THROUGH_CLEAR_TO"):
            endpoint = "from" if mapping.endswith("FROM") else "to"
            operation = {
                "kind": "through-clear",
                "endpoint": endpoint,
                "source_preset": preset if endpoint == "from"
                else (preset + preset_count - 1) % preset_count,
                "destination_preset": (preset + 1) % preset_count
                if endpoint == "from"
                else preset,
                "elapsed": 0 if endpoint == "from" else 60,
                "duration": 60,
                "selected_pixels": total,
            }
        else:
            selected = 0
            for index, pixel in enumerate(pixels):
                x = index % resolution[0]
                y = index // resolution[0]
                if crosscheck._selected_pixel(mapping, x, y, *resolution):
                    selected += 1
                else:
                    pixel[:3] = [0, 0, 0]
            operation = {
                "kind": "spatial-extract",
                "mapping": mapping,
                "selected_pixels": selected,
            }
    frame = {
        "program": program,
        "preset": preset,
        "case": case,
        "resolution": list(resolution),
        "probe": probe,
        "operation": operation,
        "pixels": pixels,
    }
    frame["sha256"] = hashlib.sha256(
        crosscheck._canonical_frame_bytes(frame)
    ).hexdigest()
    return frame


def _capture(
    programs, digest, configuration="native-debug", value=1, checkout_sha="a" * 40
):
    frames = [
        _frame(programs, *key, value=value)
        for key in sorted(crosscheck._expected_frame_keys(programs))
    ]
    return {
        "schema_version": 1,
        "checkout_sha": checkout_sha,
        "manifest_sha256": digest,
        "toolchain": crosscheck._expected_toolchain(programs, configuration),
        "configuration": configuration,
        "corpus": programs["corpus"],
        "oracle_metrics": [
            {
                "oracle_id": oracle["oracle_id"],
                "domain": "FRAMEBUFFER",
                "aggregation": "MAXIMUM",
                "value": next(
                    metric["configuration_baselines"][configuration]
                    for metric in oracle["metrics"]
                    if metric["domain"] == "FRAMEBUFFER"
                ),
                "unit": next(
                    metric["unit"]
                    for metric in oracle["metrics"]
                    if metric["domain"] == "FRAMEBUFFER"
                ),
                "sample_count": 3,
                "resolution_values": next(
                    metric["resolution_baselines"][configuration]
                    for metric in oracle["metrics"]
                    if metric["domain"] == "FRAMEBUFFER"
                ),
            }
            for oracle in ORACLES
        ],
        "frames": frames,
    }


def _test_manifest():
    programs, oracles, schema = generator.load_and_validate(MANIFEST_DIR)
    programs = copy.deepcopy(programs)
    programs["base_sha"] = "a" * 40
    programs["corpus"]["resolutions"] = [[16, 16], [32, 24]]
    for program in programs["programs"]:
        program["tolerances"]["release_wasm_framebuffer"][
            "max_differing_pixel_fraction"
        ] = 1.0
    return programs, generator.manifest_sha256(programs, oracles, schema)


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


class ManifestValidation(unittest.TestCase):
    def test_generator_main_returns_success_for_validation(self):
        self.assertEqual(generator.main([
            "--manifest-dir", str(MANIFEST_DIR), "--validate-only"
        ]), 0)

    def test_generation_is_deterministic(self):
        programs, oracles, schema = generator.load_and_validate(MANIFEST_DIR)
        first = generator.generate_header(programs, oracles, schema)
        second = generator.generate_header(programs, oracles, schema)
        self.assertEqual(first, second)
        self.assertIn("0x40000", first)
        self.assertEqual(generator.PRESET_COUNT, 24)
        self.assertEqual(generator.OPERATION_CODES["FULL_FRAME"], 18)

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

    def test_common_constraints_are_checked_through_loader(self):
        with _staged_manifest() as manifest_dir:
            schema_path = manifest_dir / "schema.json"
            schema = generator._load(schema_path)
            schema["$defs"]["sha"]["pattern"] = ".*"
            schema_path.write_text(json.dumps(schema), encoding="utf-8")
            for path in manifest_dir.glob("*.json"):
                if path.name == "schema.json":
                    continue
                document = generator._load(path)
                document["base_sha"] = "invalid"
                path.write_text(json.dumps(document), encoding="utf-8")
            with self.assertRaisesRegex(
                generator.ManifestError, "base_sha must be a full lowercase Git SHA"
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
        broken["programs"][0]["tolerances"]["release_wasm_framebuffer"][
            "unlicensed"
        ] = 1
        with self.assertRaisesRegex(
            generator.ManifestError, "unexpected field unlicensed"
        ):
            generator._validate_schema(broken, schema, schema, "programs")

    def test_schema_shape_rejects_missing_definitions(self):
        schema_path = MANIFEST_DIR / "schema.json"
        schema = generator._load(schema_path)
        schema["$defs"].pop("sha")
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
            "sha pattern": lambda doc: doc["$defs"]["sha"].pop("pattern"),
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

    def test_digest_covers_the_schema(self):
        programs, oracles, schema = generator.load_and_validate(MANIFEST_DIR)
        relaxed = copy.deepcopy(schema)
        relaxed["$defs"]["program"]["properties"]["id"].pop("pattern")
        self.assertNotEqual(
            generator.manifest_sha256(programs, oracles, schema),
            generator.manifest_sha256(programs, oracles, relaxed),
        )
        self.assertNotEqual(
            generator.generate_header(programs, oracles, schema),
            generator.generate_header(programs, oracles, relaxed),
        )

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
                self.assertRaises(SystemExit) as caught,
            ):
                capture.main()
            self.assertEqual(caught.exception.code, 2)
            self.assertIn("FRAMEBUFFER/MAXIMUM", stderr.getvalue())

    def test_crosscheck_reports_missing_framebuffer_maximum(self):
        with _staged_manifest() as manifest_dir:
            _replace_framebuffer_maximum(manifest_dir)
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                status = crosscheck.main([
                    "compare",
                    "--manifest-dir", str(manifest_dir),
                    "--base", str(manifest_dir / "base.json"),
                    "--candidate", str(manifest_dir / "candidate.json"),
                    "--base-sha", "a" * 40,
                    "--candidate-sha", "b" * 40,
                ])
            self.assertEqual(status, 1)
            self.assertIn("FRAMEBUFFER/MAXIMUM", stderr.getvalue())


class CaptureComparison(unittest.TestCase):
    def test_capture_path_comparison_releases_and_reloads_documents(self):
        class TrackedCapture(dict):
            pass

        loaded = []
        references = []
        comparisons = []

        def load(path):
            configuration = "strict" if "strict" in path.name else "release"
            document = TrackedCapture(
                name=path.name,
                configuration=configuration,
                toolchain={"configuration": configuration},
            )
            loaded.append(path.name)
            references.append(weakref.ref(document))
            return document

        def compare(base, candidate, *args, **kwargs):
            self.assertEqual(kwargs, {"oracles": []})
            documents = (base, candidate, *args[-2:])
            comparisons.append(
                tuple(
                    document["name"] if document is not None else None
                    for document in documents
                )
            )
            if len(comparisons) == 1:
                raise crosscheck.StrictFpRequired("strict-FP required")

        paths = [
            Path("base-release.json"),
            Path("candidate-release.json"),
            Path("base-strict.json"),
            Path("candidate-strict.json"),
        ]

        def compare_paths(selected):
            return crosscheck._compare_capture_paths(
                selected[0],
                selected[1],
                {},
                "digest",
                "a" * 40,
                "b" * 40,
                *selected[2:],
                oracles=[],
            )

        with (
            mock.patch.object(crosscheck, "_load_capture", new=load),
            mock.patch.object(crosscheck, "compare_captures", new=compare),
        ):
            strict_required, toolchains = compare_paths(paths[:2])
            self.assertTrue(strict_required)
            self.assertEqual(toolchains, {"release": {"configuration": "release"}})
            self.assertTrue(all(reference() is None for reference in references))
            strict_required, toolchains = compare_paths(paths)
            self.assertFalse(strict_required)
            self.assertEqual(
                toolchains,
                {
                    "release": {"configuration": "release"},
                    "strict": {"configuration": "strict"},
                },
            )
            self.assertTrue(all(reference() is None for reference in references))

        self.assertEqual(
            loaded,
            [path.name for path in paths[:2] + paths],
        )
        self.assertEqual(
            comparisons,
            [
                (*[path.name for path in paths[:2]], None, None),
                tuple(path.name for path in paths),
            ],
        )

    def test_crosscheck_main_returns_failure_for_vacuous_run(self):
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            status = crosscheck.main([
                "run", "--repository", ".", "--controller", ".",
                "--base-sha", "a" * 40, "--candidate-sha", "a" * 40,
                "--output", "output",
            ])
        self.assertEqual(status, 1)
        self.assertIn("must differ", stderr.getvalue())

    def test_orchestration_rejects_identical_revisions(self):
        with self.assertRaisesRegex(crosscheck.CrosscheckError, "must differ"):
            crosscheck.orchestrate(
                Path("repository"), Path("controller"), "a" * 40, "a" * 40,
                Path("output"))

    def test_toolchain_mismatch_is_refused(self):
        programs, digest = _test_manifest()
        base = _capture(programs, digest)
        candidate = copy.deepcopy(base)
        candidate["checkout_sha"] = "b" * 40
        candidate["toolchain"]["compiler"] = "different"
        with self.assertRaisesRegex(crosscheck.CrosscheckError, "toolchain"):
            crosscheck.compare_captures(
                base, candidate, programs, digest, "a" * 40, "b" * 40,
                oracles=ORACLES,
            )

    def test_oracle_metric_drift_is_refused(self):
        programs, digest = _test_manifest()
        base = _capture(programs, digest)
        candidate = copy.deepcopy(base)
        candidate["checkout_sha"] = "b" * 40
        candidate["oracle_metrics"][0]["value"] += 1
        with self.assertRaisesRegex(crosscheck.CrosscheckError, "oracle"):
            crosscheck.compare_captures(
                base,
                candidate,
                programs,
                digest,
                "a" * 40,
                "b" * 40,
                oracles=ORACLES,
            )

    def test_oracle_metrics_cannot_be_left_uncompared(self):
        programs, digest = _test_manifest()
        base = _capture(programs, digest)
        base["oracle_metrics"] = []
        candidate = copy.deepcopy(base)
        candidate["checkout_sha"] = "b" * 40
        with self.assertRaises(TypeError):
            crosscheck.compare_captures(
                base, candidate, programs, digest, "a" * 40, "b" * 40
            )
        with self.assertRaisesRegex(crosscheck.CrosscheckError, "coverage"):
            crosscheck.compare_captures(
                base, candidate, programs, digest, "a" * 40, "b" * 40,
                oracles=ORACLES,
            )

    def test_observed_toolchain_mismatch_is_refused(self):
        programs, _ = _test_manifest()
        cache = {
            "CMAKE_CXX_COMPILER": "clang++",
            "CMAKE_COMMAND": "cmake",
            "CMAKE_BUILD_TYPE": "Debug",
            "HS_PULLBACK_CAPTURE_PRESET": "tests",
        }
        with (
            mock.patch.object(capture, "_cache_values", return_value=cache),
            mock.patch.object(
                capture,
                "_compiler_metadata",
                return_value=(Path("clang++"), "22.0.0"),
            ),
            mock.patch.object(
                capture,
                "_first_line",
                side_effect=["clang version 22.0.0", "cmake version 4.2.3"],
            ),
            self.assertRaisesRegex(capture.CaptureError, "manifest pin"),
        ):
            capture.attest_toolchain(Path("build"), "native-debug", programs)

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
                "sample_count": 3,
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
                "sample_count": 3,
            }
            for oracle in ORACLES
        }
        with self.assertRaisesRegex(capture.CaptureError, "accepted limit"):
            capture.validate_backend_oracles(
                metrics, ORACLES, "native-debug", (96, 20)
            )

    def test_native_deterministic_replay_passes(self):
        programs, digest = _test_manifest()
        base = _capture(programs, digest)
        candidate = copy.deepcopy(base)
        candidate["checkout_sha"] = "b" * 40
        crosscheck.compare_captures(
            base, candidate, programs, digest, "a" * 40, "b" * 40,
            oracles=ORACLES,
        )

    def test_release_mismatch_requires_strict_fp_identity(self):
        programs, digest = _test_manifest()
        base = _capture(programs, digest, "wasm-release", 1)
        candidate = copy.deepcopy(base)
        candidate["checkout_sha"] = "b" * 40
        candidate["frames"][0] = _frame(
            programs,
            candidate["frames"][0]["program"],
            candidate["frames"][0]["preset"],
            candidate["frames"][0]["case"],
            candidate["frames"][0]["resolution"],
            candidate["frames"][0]["probe"],
            2,
        )
        with self.assertRaisesRegex(crosscheck.StrictFpRequired, "strict-FP"):
            crosscheck.compare_captures(
                base, candidate, programs, digest, "a" * 40, "b" * 40,
                oracles=ORACLES,
            )
        strict_base = _capture(programs, digest, "wasm-strict-fp", 3)
        strict_candidate = copy.deepcopy(strict_base)
        strict_candidate["checkout_sha"] = "b" * 40
        crosscheck.compare_captures(
            base,
            candidate,
            programs,
            digest,
            "a" * 40,
            "b" * 40,
            strict_base,
            strict_candidate,
            oracles=ORACLES,
        )

    def test_release_pixel_fraction_is_rounded_to_whole_pixels(self):
        self.assertTrue(crosscheck._within_differing_pixel_limit(2, 1920, 0.001))
        self.assertFalse(
            crosscheck._within_differing_pixel_limit(3, 1920, 0.001)
        )
        self.assertFalse(crosscheck._within_differing_pixel_limit(1, 1920, 0.0))

    def test_raw_hash_pixel_shape_and_provenance_are_checked(self):
        programs, digest = _test_manifest()
        base = _capture(programs, digest)
        candidate = copy.deepcopy(base)
        candidate["checkout_sha"] = "b" * 40
        candidate["frames"][0]["sha256"] = "0" * 64
        with self.assertRaisesRegex(crosscheck.CrosscheckError, "raw hash"):
            crosscheck.compare_captures(
                base, candidate, programs, digest, "a" * 40, "b" * 40,
                oracles=ORACLES,
            )
        candidate = copy.deepcopy(base)
        candidate["checkout_sha"] = "b" * 40
        candidate["frames"][0]["pixels"].pop()
        with self.assertRaisesRegex(crosscheck.CrosscheckError, "pixels"):
            crosscheck.compare_captures(
                base, candidate, programs, digest, "a" * 40, "b" * 40,
                oracles=ORACLES,
            )
        candidate = copy.deepcopy(base)
        candidate["checkout_sha"] = "b" * 40
        candidate["manifest_sha256"] = "0" * 64
        with self.assertRaisesRegex(crosscheck.CrosscheckError, "manifest hash"):
            crosscheck.compare_captures(
                base, candidate, programs, digest, "a" * 40, "b" * 40,
                oracles=ORACLES,
            )

    def test_probe_operations_are_targeted_and_not_aliases(self):
        programs, digest = _test_manifest()
        capture_doc = _capture(programs, digest)
        spatial = [
            frame
            for frame in capture_doc["frames"]
            if frame["operation"]["kind"] == "spatial-extract"
        ]
        self.assertGreater(len({frame["sha256"] for frame in spatial}), 1)
        self.assertTrue(
            all(
                frame["operation"]["selected_pixels"]
                < frame["resolution"][0] * frame["resolution"][1]
                for frame in spatial
            )
        )
        candidate = copy.deepcopy(capture_doc)
        candidate["checkout_sha"] = "b" * 40
        target = next(
            frame
            for frame in candidate["frames"]
            if frame["operation"]["kind"] == "spatial-extract"
        )
        mapping = target["operation"]["mapping"]
        width, height = target["resolution"]
        outside = next(
            index
            for index in range(width * height)
            if not crosscheck._selected_pixel(
                mapping, index % width, index // width, width, height
            )
        )
        target["pixels"][outside][0] = 1
        target["sha256"] = hashlib.sha256(
            crosscheck._canonical_frame_bytes(target)
        ).hexdigest()
        with self.assertRaisesRegex(crosscheck.CrosscheckError, "outside mapping"):
            crosscheck.compare_captures(
                capture_doc, candidate, programs, digest, "a" * 40, "b" * 40,
                oracles=ORACLES,
            )

    def test_strict_capture_provenance_and_identity_are_checked(self):
        programs, digest = _test_manifest()
        base = _capture(programs, digest, "wasm-release")
        candidate = copy.deepcopy(base)
        candidate["checkout_sha"] = "b" * 40
        strict_base = _capture(programs, digest, "wasm-strict-fp")
        strict_candidate = copy.deepcopy(strict_base)
        strict_candidate["checkout_sha"] = "b" * 40
        strict_candidate["corpus"] = {"wrong": True}
        with self.assertRaisesRegex(crosscheck.CrosscheckError, "corpus"):
            crosscheck.compare_captures(
                base,
                candidate,
                programs,
                digest,
                "a" * 40,
                "b" * 40,
                strict_base,
                strict_candidate,
                oracles=ORACLES,
            )
        strict_candidate = copy.deepcopy(strict_base)
        strict_candidate["checkout_sha"] = "b" * 40
        strict_candidate["frames"][0] = _frame(
            programs,
            strict_candidate["frames"][0]["program"],
            strict_candidate["frames"][0]["preset"],
            strict_candidate["frames"][0]["case"],
            strict_candidate["frames"][0]["resolution"],
            strict_candidate["frames"][0]["probe"],
            2,
        )
        with self.assertRaisesRegex(
            crosscheck.CrosscheckError, "strict-FP frame differs"
        ):
            crosscheck.compare_captures(
                base,
                candidate,
                programs,
                digest,
                "a" * 40,
                "b" * 40,
                strict_base,
                strict_candidate,
                oracles=ORACLES,
            )

    def test_backend_stream_is_complete(self):
        programs, _ = _test_manifest()
        specs = capture.operation_specs(programs)
        data = bytearray(b"HSPB")
        data.extend(struct.pack("<HHHIH", 3, 1, 1, len(specs), len(ORACLES)))
        for spec in specs:
            name = spec["name"].encode()
            mapping = spec["mapping"]
            data.extend(
                struct.pack(
                    "<HHH",
                    spec["preset"],
                    capture.OPERATION_CODES[mapping],
                    len(name),
                )
            )
            data.extend(name)
            data.extend(struct.pack("<HHHHI", 0, 1, 0, 60, 1))
            data.extend(struct.pack("<HHH", spec["preset"], 3, 4))
        for oracle in ORACLES:
            name = oracle["oracle_id"].encode()
            data.extend(struct.pack("<H", len(name)))
            data.extend(name)
            data.extend(struct.pack("<HI", 0, 3))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "capture.bin"
            path.write_bytes(data)
            resolution, records, metrics = capture.load_backend(
                path, programs, ORACLES
            )
        self.assertEqual(resolution, (1, 1))
        self.assertEqual(len(records), len(specs))
        self.assertEqual(len(metrics), len(ORACLES))
        pixels = records[(11, "case:interior")]["pixels"]
        self.assertEqual(pixels, struct.pack("<HHHH", 11, 3, 4, 65535))
        # The record's byte string is what the frame hash is taken over, so it
        # must equal what the crosschecker rebuilds from the published lists.
        self.assertEqual(
            pixels,
            crosscheck._canonical_frame_bytes(
                {"resolution": [1, 1], "pixels": [[11, 3, 4, 65535]]}
            ),
        )

    def test_streamed_capture_matches_a_materialized_document(self):
        # The corpus is hashed and diffed run to run, so its encoding is part
        # of the contract: streaming the frames must not move a single byte.
        frames = [
            {
                "program": "P",
                "preset": 3,
                "case": "default",
                "resolution": [2, 1],
                "probe": "steady",
                "operation": {"kind": "parameter-case", "selected_pixels": 2},
                "pixels": struct.pack("<HHHHHHHH", 1, 2, 3, 65535,
                                      65535, 0, 40000, 65535),
                "sha256": "a" * 64,
            }
        ]
        head = {"schema_version": 1, "configuration": "native-debug"}
        expanded = dict(head, frames=[
            dict(frames[0], pixels=[[1, 2, 3, 65535], [65535, 0, 40000, 65535]])
        ])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "nested" / "capture.json"
            capture.write_capture(head, frames, path)
            written = path.read_text(encoding="utf-8")
        self.assertEqual(
            written, json.dumps(expanded, separators=(",", ":")) + "\n"
        )

    def test_build_directories_are_isolated(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sources = {"base": root / "base", "candidate": root / "candidate"}
            builds = {
                name: {
                    configuration: root / f"build-{name}-{configuration}"
                    for configuration in ("native", "release", "strict")
                }
                for name in sources
            }
            crosscheck.validate_build_isolation(sources, builds)
            builds["candidate"]["strict"] = builds["base"]["strict"]
            with self.assertRaisesRegex(crosscheck.CrosscheckError, "not isolated"):
                crosscheck.validate_build_isolation(sources, builds)


if __name__ == "__main__":
    unittest.main()
