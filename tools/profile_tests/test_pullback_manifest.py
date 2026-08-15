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
    pixels = [[value, value, value, 65535]]
    payload = (f"{program}:{case}:{resolution}:{probe}:{value}".encode() +
               bytes(channel & 0xff for pixel in pixels for channel in pixel))
    return {
        "program": program,
        "case": case,
        "resolution": list(resolution),
        "probe": probe,
        "sha256": hashlib.sha256(payload).hexdigest(),
        "pixels": pixels,
    }


def _capture(configuration="native-debug", value=1):
    programs, _ = generator.load_and_validate(MANIFEST_DIR)
    frames = [_frame(*key, value=value)
              for key in sorted(crosscheck._expected_frame_keys(programs))]
    return {
        "toolchains": programs["toolchains"],
        "configuration": configuration,
        "corpus_version": 1,
        "frames": frames,
    }


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


class CaptureComparison(unittest.TestCase):
    def test_toolchain_mismatch_is_refused(self):
        programs, _ = generator.load_and_validate(MANIFEST_DIR)
        base = _capture()
        candidate = copy.deepcopy(base)
        candidate["toolchains"]["native"]["compiler"] = "different"
        with self.assertRaisesRegex(crosscheck.CrosscheckError, "toolchains"):
            crosscheck.compare_captures(base, candidate, programs)

    def test_native_deterministic_replay_passes(self):
        programs, _ = generator.load_and_validate(MANIFEST_DIR)
        capture = _capture()
        crosscheck.compare_captures(capture, copy.deepcopy(capture), programs)

    def test_release_mismatch_requires_strict_fp_identity(self):
        programs, _ = generator.load_and_validate(MANIFEST_DIR)
        base = _capture("wasm-release", 1)
        candidate = copy.deepcopy(base)
        candidate["frames"][0]["sha256"] = "different"
        with self.assertRaisesRegex(crosscheck.CrosscheckError, "strict-FP"):
            crosscheck.compare_captures(base, candidate, programs)
        strict_base = _capture("wasm-strict-fp", 3)
        strict_candidate = copy.deepcopy(strict_base)
        crosscheck.compare_captures(base, candidate, programs, strict_base,
                                    strict_candidate)

    def test_build_directories_are_checkout_local(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            base = root / "base"
            candidate = root / "candidate"
            crosscheck.validate_isolation(base, candidate, base / "build",
                                          candidate / "build")
            with self.assertRaisesRegex(crosscheck.CrosscheckError, "escapes"):
                crosscheck.validate_isolation(base, candidate, root / "shared",
                                              candidate / "build")


if __name__ == "__main__":
    unittest.main()
