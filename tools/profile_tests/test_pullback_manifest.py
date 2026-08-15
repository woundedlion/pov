"""Self-tests for pullback manifest generation and capture comparison."""

from __future__ import annotations

import copy
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path


TOOLS = Path(__file__).resolve().parents[1]
ROOT = TOOLS.parent
sys.path.insert(0, str(TOOLS))

import generate_pullback_manifest_header as generator  # noqa: E402
import pullback_crosscheck as crosscheck  # noqa: E402


MANIFEST_DIR = ROOT / "tests/data/pullback"


def _frame(program, case, resolution, probe, value=1):
    pixels = [[value, value, value, 65535]
              for _ in range(resolution[0] * resolution[1])]
    frame = {
        "program": program,
        "case": case,
        "resolution": list(resolution),
        "probe": probe,
        "pixels": pixels,
    }
    frame["sha256"] = hashlib.sha256(
        crosscheck._canonical_frame_bytes(frame)).hexdigest()
    return frame


def _capture(programs, digest, configuration="native-debug", value=1,
             checkout_sha="a" * 40):
    frames = [_frame(*key, value=value)
              for key in sorted(crosscheck._expected_frame_keys(programs))]
    return {
        "schema_version": 1,
        "checkout_sha": checkout_sha,
        "manifest_sha256": digest,
        "toolchains": programs["toolchains"],
        "configuration": configuration,
        "corpus": programs["corpus"],
        "frames": frames,
    }


def _test_manifest():
    programs, oracles = generator.load_and_validate(MANIFEST_DIR)
    programs = copy.deepcopy(programs)
    programs["base_sha"] = "a" * 40
    programs["corpus"]["resolutions"] = [[1, 1], [2, 1]]
    for program in programs["programs"]:
        program["tolerances"]["release_wasm_framebuffer"][
            "max_differing_pixel_fraction"] = 1.0
    return programs, generator.manifest_sha256(programs, oracles)


class ManifestValidation(unittest.TestCase):
    def test_generation_is_deterministic(self):
        programs, oracles = generator.load_and_validate(MANIFEST_DIR)
        first = generator.generate_header(programs, oracles)
        second = generator.generate_header(programs, oracles)
        self.assertEqual(first, second)
        self.assertIn("0x300", first)
        runtime = json.loads(generator.generate_runtime_manifest(
            programs, oracles, "a" * 40))
        self.assertEqual(runtime["capture_sha"], "a" * 40)
        self.assertEqual(len(runtime["manifest_sha256"]), 64)

    def test_manifest_completeness_is_checked(self):
        programs, _ = generator.load_and_validate(MANIFEST_DIR)
        broken = copy.deepcopy(programs)
        broken["programs"][0]["parameter_cases"].pop()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "programs.json"
            path.write_text(json.dumps(broken), encoding="utf-8")
            with self.assertRaises(generator.ManifestError):
                generator._validate_programs(broken, path)

    def test_nested_schema_rejects_unknown_fields(self):
        programs, _ = generator.load_and_validate(MANIFEST_DIR)
        schema = generator._load(MANIFEST_DIR / "schema.json")
        broken = copy.deepcopy(programs)
        broken["programs"][0]["tolerances"]["release_wasm_framebuffer"][
            "unlicensed"] = 1
        with self.assertRaisesRegex(generator.ManifestError,
                                    "unexpected field unlicensed"):
            generator._validate_schema(broken, schema, schema, "programs")

    def test_schema_const_rejects_boolean_versions(self):
        programs, _ = generator.load_and_validate(MANIFEST_DIR)
        schema = generator._load(MANIFEST_DIR / "schema.json")
        for field_path in (("schema_version",), ("corpus", "version")):
            broken = copy.deepcopy(programs)
            target = broken
            for field in field_path[:-1]:
                target = target[field]
            target[field_path[-1]] = True
            with self.subTest(field=".".join(field_path)):
                with self.assertRaises(generator.ManifestError):
                    generator._validate_schema(
                        broken, schema, schema, "programs")
                with self.assertRaises(generator.ManifestError):
                    generator._validate_programs(
                        broken, MANIFEST_DIR / "programs.json")
        for keyword in ("const", "enum"):
            rule = {keyword: 1 if keyword == "const" else [1]}
            with self.subTest(keyword=keyword):
                with self.assertRaises(generator.ManifestError):
                    generator._validate_schema(
                        True, rule, {"$defs": {}}, "numeric")

    def test_non_json_numeric_constants_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "invalid.json"
            for token in ("Infinity", "NaN"):
                with self.subTest(token=token):
                    path.write_text(f'{{"value": {token}}}', encoding="utf-8")
                    with self.assertRaisesRegex(
                            generator.ManifestError,
                            f"non-JSON numeric constant {token}"):
                        generator._load(path)


class CaptureComparison(unittest.TestCase):
    def test_toolchain_mismatch_is_refused(self):
        programs, digest = _test_manifest()
        base = _capture(programs, digest)
        candidate = copy.deepcopy(base)
        candidate["checkout_sha"] = "b" * 40
        candidate["toolchains"]["native"]["compiler"] = "different"
        with self.assertRaisesRegex(crosscheck.CrosscheckError, "toolchains"):
            crosscheck.compare_captures(
                base, candidate, programs, digest, "a" * 40, "b" * 40)

    def test_native_deterministic_replay_passes(self):
        programs, digest = _test_manifest()
        base = _capture(programs, digest)
        candidate = copy.deepcopy(base)
        candidate["checkout_sha"] = "b" * 40
        crosscheck.compare_captures(
            base, candidate, programs, digest, "a" * 40, "b" * 40)

    def test_release_mismatch_requires_strict_fp_identity(self):
        programs, digest = _test_manifest()
        base = _capture(programs, digest, "wasm-release", 1)
        candidate = copy.deepcopy(base)
        candidate["checkout_sha"] = "b" * 40
        candidate["frames"][0] = _frame(
            candidate["frames"][0]["program"],
            candidate["frames"][0]["case"],
            candidate["frames"][0]["resolution"],
            candidate["frames"][0]["probe"], 2)
        with self.assertRaisesRegex(crosscheck.CrosscheckError, "strict-FP"):
            crosscheck.compare_captures(
                base, candidate, programs, digest, "a" * 40, "b" * 40)
        strict_base = _capture(programs, digest, "wasm-strict-fp", 3)
        strict_candidate = copy.deepcopy(strict_base)
        strict_candidate["checkout_sha"] = "b" * 40
        crosscheck.compare_captures(
            base, candidate, programs, digest, "a" * 40, "b" * 40,
            strict_base, strict_candidate)

    def test_raw_hash_pixel_shape_and_provenance_are_checked(self):
        programs, digest = _test_manifest()
        base = _capture(programs, digest)
        candidate = copy.deepcopy(base)
        candidate["checkout_sha"] = "b" * 40
        candidate["frames"][0]["sha256"] = "0" * 64
        with self.assertRaisesRegex(crosscheck.CrosscheckError, "raw hash"):
            crosscheck.compare_captures(
                base, candidate, programs, digest, "a" * 40, "b" * 40)
        candidate = copy.deepcopy(base)
        candidate["checkout_sha"] = "b" * 40
        candidate["frames"][0]["pixels"].pop()
        with self.assertRaisesRegex(crosscheck.CrosscheckError, "pixels"):
            crosscheck.compare_captures(
                base, candidate, programs, digest, "a" * 40, "b" * 40)
        candidate = copy.deepcopy(base)
        candidate["checkout_sha"] = "b" * 40
        candidate["manifest_sha256"] = "0" * 64
        with self.assertRaisesRegex(crosscheck.CrosscheckError, "manifest hash"):
            crosscheck.compare_captures(
                base, candidate, programs, digest, "a" * 40, "b" * 40)

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
                base, candidate, programs, digest, "a" * 40, "b" * 40,
                strict_base, strict_candidate)
        strict_candidate = copy.deepcopy(strict_base)
        strict_candidate["checkout_sha"] = "b" * 40
        strict_candidate["frames"][0] = _frame(
            strict_candidate["frames"][0]["program"],
            strict_candidate["frames"][0]["case"],
            strict_candidate["frames"][0]["resolution"],
            strict_candidate["frames"][0]["probe"], 2)
        with self.assertRaisesRegex(crosscheck.CrosscheckError,
                                    "strict-FP frame differs"):
            crosscheck.compare_captures(
                base, candidate, programs, digest, "a" * 40, "b" * 40,
                strict_base, strict_candidate)

if __name__ == "__main__":
    unittest.main()
