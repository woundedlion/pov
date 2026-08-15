#!/usr/bin/env python3
"""Validate and compare canonical pullback capture files.

Capture production is intentionally outside this comparator. The framebuffer
gate remains pending until a repository-owned native/WASM producer exists.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

from generate_pullback_manifest_header import (
    ManifestError,
    SHA_RE,
    load_and_validate,
    manifest_sha256,
)


class CrosscheckError(RuntimeError):
    """Cross-checkout capture provenance or output is invalid."""


def _load_capture(path: Path) -> dict:
    try:
        capture = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise CrosscheckError(f"{path}: {error}") from error
    if not isinstance(capture, dict):
        raise CrosscheckError(f"{path}: capture root must be an object")
    return capture


def _canonical_frame_bytes(frame: dict) -> bytes:
    resolution = frame.get("resolution")
    if not isinstance(resolution, list) or len(resolution) != 2 or not all(
            isinstance(value, int) and value > 0 for value in resolution):
        raise CrosscheckError("frame resolution must contain two positive integers")
    pixels = frame.get("pixels")
    expected_count = resolution[0] * resolution[1]
    if not isinstance(pixels, list) or len(pixels) != expected_count:
        raise CrosscheckError(
            f"frame has {len(pixels) if isinstance(pixels, list) else 'invalid'} "
            f"pixels; expected {expected_count}")
    raw = bytearray()
    for pixel in pixels:
        if not isinstance(pixel, list) or len(pixel) != 4:
            raise CrosscheckError("pixels must contain four 16-bit channels")
        for channel in pixel:
            if not isinstance(channel, int) or isinstance(channel, bool) or not (
                    0 <= channel <= 65535):
                raise CrosscheckError("pixel channels must be 16-bit integers")
            raw.extend(channel.to_bytes(2, "little"))
    return bytes(raw)


def _frame_map(capture: dict) -> dict[tuple, tuple[dict, str]]:
    frames = capture.get("frames")
    if not isinstance(frames, list):
        raise CrosscheckError("capture.frames must be an array")
    mapped = {}
    for frame in frames:
        if not isinstance(frame, dict):
            raise CrosscheckError("capture frame must be an object")
        if set(frame) != {
                "program", "case", "resolution", "probe", "sha256", "pixels"}:
            raise CrosscheckError("capture frame fields are incomplete or unknown")
        resolution = tuple(frame.get("resolution", ()))
        key = (frame.get("program"), frame.get("case"), resolution,
               frame.get("probe"))
        if key in mapped:
            raise CrosscheckError(f"duplicate capture frame {key}")
        digest = hashlib.sha256(_canonical_frame_bytes(frame)).hexdigest()
        reported = frame["sha256"]
        if reported != digest:
            raise CrosscheckError(f"frame raw hash is invalid: {key}")
        mapped[key] = (frame, digest)
    return mapped


def _expected_frame_keys(programs: dict) -> set[tuple]:
    expected = set()
    resolutions = programs["corpus"]["resolutions"]
    for program in programs["programs"]:
        program_id = program["id"]
        for resolution in resolutions:
            dimensions = tuple(resolution)
            for case in program["parameter_cases"]:
                expected.add((program_id, case["id"], dimensions, "steady"))
            for category, probes in program["probes"].items():
                for probe in probes:
                    expected.add((program_id, "default", dimensions,
                                  f"{category}:{probe}"))
    return expected


def _compare_pixels(base: list, candidate: list) -> tuple[int, int, int]:
    if len(base) != len(candidate):
        raise CrosscheckError("frame pixel counts differ")
    maximum = 0
    differing = 0
    for base_pixel, candidate_pixel in zip(base, candidate, strict=True):
        if len(base_pixel) != 4 or len(candidate_pixel) != 4:
            raise CrosscheckError("pixels must contain four 16-bit channels")
        deltas = [abs(a - b) for a, b in zip(base_pixel, candidate_pixel,
                                             strict=True)]
        maximum = max(maximum, *deltas)
        differing += any(deltas)
    return maximum, differing, len(base)


def _validate_capture(capture: dict, programs: dict, digest: str,
                      configuration: str, checkout_sha: str) -> dict:
    required = {
        "schema_version", "checkout_sha", "manifest_sha256", "toolchains",
        "configuration", "corpus", "frames",
    }
    if set(capture) != required:
        raise CrosscheckError("capture provenance fields are incomplete or unknown")
    if capture["schema_version"] != 1:
        raise CrosscheckError("capture schema version mismatch")
    if capture["checkout_sha"] != checkout_sha:
        raise CrosscheckError("capture checkout SHA mismatch")
    if capture["manifest_sha256"] != digest:
        raise CrosscheckError("capture manifest hash mismatch")
    if capture["toolchains"] != programs["toolchains"]:
        raise CrosscheckError("capture toolchains differ from manifest pins")
    if capture["configuration"] != configuration:
        raise CrosscheckError("capture configuration mismatch")
    if capture["corpus"] != programs["corpus"]:
        raise CrosscheckError("capture corpus differs from manifest")
    frames = _frame_map(capture)
    if frames.keys() != _expected_frame_keys(programs):
        raise CrosscheckError("capture does not cover every manifest case and probe")
    return frames


def compare_captures(base: dict, candidate: dict, programs: dict,
                     digest: str, base_sha: str, candidate_sha: str,
                     strict_base: dict | None = None,
                     strict_candidate: dict | None = None) -> None:
    if base_sha != programs["base_sha"]:
        raise CrosscheckError("base checkout SHA differs from manifest pin")
    if SHA_RE.fullmatch(candidate_sha) is None:
        raise CrosscheckError("candidate checkout SHA must be a full lowercase SHA")
    configuration = base.get("configuration")
    if configuration not in ("native-debug", "wasm-release"):
        raise CrosscheckError(f"unsupported configuration {configuration}")
    base_frames = _validate_capture(base, programs, digest, configuration,
                                    base_sha)
    candidate_frames = _validate_capture(
        candidate, programs, digest, configuration, candidate_sha)
    strict_frames = None
    if strict_base is not None or strict_candidate is not None:
        if strict_base is None or strict_candidate is None:
            raise CrosscheckError("strict-FP captures must be supplied as a pair")
        strict_frames = (
            _validate_capture(strict_base, programs, digest, "wasm-strict-fp",
                              base_sha),
            _validate_capture(strict_candidate, programs, digest,
                              "wasm-strict-fp", candidate_sha),
        )
        for key in strict_frames[0]:
            if strict_frames[0][key][1] != strict_frames[1][key][1]:
                raise CrosscheckError(f"strict-FP frame differs: {key}")
    program_tolerances = {
        entry["id"]: entry["tolerances"] for entry in programs["programs"]
    }
    release_differences = []
    for key, (base_frame, base_digest) in base_frames.items():
        candidate_frame, candidate_digest = candidate_frames[key]
        if configuration == "native-debug":
            if base_digest != candidate_digest:
                raise CrosscheckError(f"native frame differs: {key}")
            continue
        if base_digest == candidate_digest:
            continue
        maximum, differing, count = _compare_pixels(
            base_frame["pixels"], candidate_frame["pixels"])
        limits = program_tolerances[key[0]]["release_wasm_framebuffer"]
        release_differences.append((key, maximum, differing, count, limits))
    if release_differences:
        if strict_frames is None:
            raise CrosscheckError("release mismatch requires strict-FP captures")
        for key, maximum, differing, count, limits in release_differences:
            if maximum > limits["max_channel_delta_u16"]:
                raise CrosscheckError(
                    f"release channel delta exceeds limit: {key}")
            if differing / count > limits["max_differing_pixel_fraction"]:
                raise CrosscheckError(
                    f"release differing-pixel share exceeds limit: {key}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest-dir", type=Path, required=True)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--base-sha", required=True)
    parser.add_argument("--candidate-sha", required=True)
    parser.add_argument("--strict-base", type=Path)
    parser.add_argument("--strict-candidate", type=Path)
    args = parser.parse_args()
    try:
        programs, oracles = load_and_validate(args.manifest_dir)
        compare_captures(
            _load_capture(args.base), _load_capture(args.candidate), programs,
            manifest_sha256(programs, oracles), args.base_sha,
            args.candidate_sha,
            _load_capture(args.strict_base) if args.strict_base else None,
            _load_capture(args.strict_candidate)
            if args.strict_candidate else None)
    except (CrosscheckError, ManifestError) as error:
        print(f"pullback_crosscheck: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
