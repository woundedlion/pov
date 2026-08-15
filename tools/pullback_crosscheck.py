#!/usr/bin/env python3
"""Compare isolated base/candidate pullback framebuffer captures."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shlex
import subprocess
import sys
import tempfile
from pathlib import Path

from generate_pullback_manifest_header import ManifestError, load_and_validate


CONFIGURATIONS = ("native-debug", "wasm-release")


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


def _frame_map(capture: dict) -> dict[tuple, dict]:
    frames = capture.get("frames")
    if not isinstance(frames, list):
        raise CrosscheckError("capture.frames must be an array")
    mapped = {}
    for frame in frames:
        if not isinstance(frame, dict):
            raise CrosscheckError("capture frame must be an object")
        resolution = tuple(frame.get("resolution", ()))
        key = (frame.get("program"), frame.get("case"), resolution,
               frame.get("probe"))
        if key in mapped:
            raise CrosscheckError(f"duplicate capture frame {key}")
        mapped[key] = frame
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


def compare_captures(base: dict, candidate: dict, programs: dict,
                     strict_base: dict | None = None,
                     strict_candidate: dict | None = None) -> None:
    for field in ("toolchains", "configuration", "corpus_version"):
        if base.get(field) != candidate.get(field):
            raise CrosscheckError(f"capture {field} mismatch")
    configuration = base.get("configuration")
    base_frames = _frame_map(base)
    candidate_frames = _frame_map(candidate)
    if base.get("toolchains") != programs["toolchains"]:
        raise CrosscheckError("capture toolchains differ from manifest pins")
    if base_frames.keys() != candidate_frames.keys():
        raise CrosscheckError("capture manifests are incomplete or differ")
    expected_keys = _expected_frame_keys(programs)
    if base_frames.keys() != expected_keys:
        raise CrosscheckError("capture does not cover every manifest case and probe")
    program_tolerances = {
        entry["id"]: entry["tolerances"] for entry in programs["programs"]
    }
    release_differences = []
    for key, base_frame in base_frames.items():
        candidate_frame = candidate_frames[key]
        if configuration == "native-debug":
            if base_frame.get("sha256") != candidate_frame.get("sha256"):
                raise CrosscheckError(f"native frame differs: {key}")
            continue
        if configuration != "wasm-release":
            raise CrosscheckError(f"unsupported configuration {configuration}")
        if base_frame.get("sha256") == candidate_frame.get("sha256"):
            continue
        pixels = base_frame.get("pixels"), candidate_frame.get("pixels")
        if not all(isinstance(value, list) for value in pixels):
            raise CrosscheckError(f"release mismatch lacks pixels: {key}")
        maximum, differing, count = _compare_pixels(*pixels)
        limits = program_tolerances[key[0]]["release_wasm_framebuffer"]
        release_differences.append((key, maximum, differing, count, limits))
    if release_differences:
        if strict_base is None or strict_candidate is None:
            raise CrosscheckError("release mismatch requires strict-FP captures")
        if strict_base.get("configuration") != "wasm-strict-fp" or \
                strict_candidate.get("configuration") != "wasm-strict-fp":
            raise CrosscheckError("strict-FP captures use the wrong configuration")
        strict_frames = _frame_map(strict_base), _frame_map(strict_candidate)
        if strict_frames[0].keys() != base_frames.keys() or \
                strict_frames[1].keys() != base_frames.keys():
            raise CrosscheckError("strict-FP capture manifests differ")
        for key in base_frames:
            if strict_frames[0][key].get("sha256") != \
                    strict_frames[1][key].get("sha256"):
                raise CrosscheckError(f"strict-FP frame differs: {key}")
        for key, maximum, differing, count, limits in release_differences:
            if maximum > limits["max_channel_delta_u16"]:
                raise CrosscheckError(
                    f"release channel delta exceeds limit: {key}")
            if differing / count > limits["max_differing_pixel_fraction"]:
                raise CrosscheckError(
                    f"release differing-pixel share exceeds limit: {key}")


def validate_isolation(base_checkout: Path, candidate_checkout: Path,
                       base_build: Path, candidate_build: Path) -> None:
    paths = [path.resolve() for path in (base_checkout, candidate_checkout,
                                        base_build, candidate_build)]
    base_checkout, candidate_checkout, base_build, candidate_build = paths
    if base_checkout == candidate_checkout:
        raise CrosscheckError("base and candidate checkouts are identical")
    if base_checkout in candidate_checkout.parents or \
            candidate_checkout in base_checkout.parents:
        raise CrosscheckError("one checkout is nested inside the other")
    if base_checkout not in base_build.parents:
        raise CrosscheckError("base build directory escapes its checkout")
    if candidate_checkout not in candidate_build.parents:
        raise CrosscheckError("candidate build directory escapes its checkout")


def _run_capture(command: str, checkout: Path, build: Path, manifest: Path,
                 configuration: str, output: Path) -> None:
    rendered = command.format(checkout=checkout, build=build, manifest=manifest,
                              configuration=configuration, output=output)
    environment = os.environ.copy()
    environment["CCACHE_DISABLE"] = "1"
    environment["CMAKE_CXX_COMPILER_LAUNCHER"] = ""
    subprocess.run(shlex.split(rendered), cwd=checkout, env=environment,
                   check=True)


def orchestrate(repository: Path, base_sha: str, candidate_sha: str,
                manifest: Path, capture_command: str, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="pullback-crosscheck-") as directory:
        root = Path(directory)
        base_checkout = root / "base"
        candidate_checkout = root / "candidate"
        subprocess.run(["git", "-C", str(repository), "worktree", "add",
                        "--detach", str(base_checkout), base_sha], check=True)
        try:
            subprocess.run(["git", "-C", str(repository), "worktree", "add",
                            "--detach", str(candidate_checkout), candidate_sha],
                           check=True)
            try:
                base_build = base_checkout / "build" / "pullback-crosscheck"
                candidate_build = (candidate_checkout / "build" /
                                   "pullback-crosscheck")
                validate_isolation(base_checkout, candidate_checkout, base_build,
                                   candidate_build)
                for configuration in CONFIGURATIONS:
                    for name, checkout, build in (
                            ("base", base_checkout, base_build),
                            ("candidate", candidate_checkout, candidate_build)):
                        capture = output / f"{name}-{configuration}.json"
                        _run_capture(capture_command, checkout, build, manifest,
                                     configuration, capture)
                programs, _ = load_and_validate(manifest.parent)
                base_native = _load_capture(output / "base-native-debug.json")
                candidate_native = _load_capture(
                    output / "candidate-native-debug.json")
                compare_captures(base_native, candidate_native, programs)
                base_release = _load_capture(output / "base-wasm-release.json")
                candidate_release = _load_capture(
                    output / "candidate-wasm-release.json")
                try:
                    compare_captures(base_release, candidate_release, programs)
                except CrosscheckError as error:
                    if "requires strict-FP" not in str(error):
                        raise
                    for name, checkout, build in (
                            ("base", base_checkout, base_build),
                            ("candidate", candidate_checkout, candidate_build)):
                        capture = output / f"{name}-wasm-strict-fp.json"
                        _run_capture(capture_command, checkout, build, manifest,
                                     "wasm-strict-fp", capture)
                    compare_captures(
                        base_release, candidate_release, programs,
                        _load_capture(output / "base-wasm-strict-fp.json"),
                        _load_capture(output /
                                      "candidate-wasm-strict-fp.json"))
                captures = {}
                for capture in sorted(output.glob("*.json")):
                    captures[capture.name] = hashlib.sha256(
                        capture.read_bytes()).hexdigest()
                summary = {
                    "schema_version": 1,
                    "base_sha": base_sha,
                    "head_sha": candidate_sha,
                    "toolchains": programs["toolchains"],
                    "manifest": str(manifest),
                    "capture_sha256": captures,
                }
                (output / "summary.json").write_text(
                    json.dumps(summary, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8", newline="\n")
            finally:
                subprocess.run(["git", "-C", str(repository), "worktree",
                                "remove", str(candidate_checkout)], check=True)
        finally:
            subprocess.run(["git", "-C", str(repository), "worktree", "remove",
                            str(base_checkout)], check=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    compare = subparsers.add_parser("compare")
    compare.add_argument("--manifest-dir", type=Path, required=True)
    compare.add_argument("--base", type=Path, required=True)
    compare.add_argument("--candidate", type=Path, required=True)
    compare.add_argument("--strict-base", type=Path)
    compare.add_argument("--strict-candidate", type=Path)
    run = subparsers.add_parser("run")
    run.add_argument("--repository", type=Path, required=True)
    run.add_argument("--base-sha", required=True)
    run.add_argument("--candidate-sha", required=True)
    run.add_argument("--manifest", type=Path, required=True)
    run.add_argument("--capture-command", required=True)
    run.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command == "run":
            orchestrate(args.repository, args.base_sha, args.candidate_sha,
                        args.manifest, args.capture_command, args.output)
        else:
            programs, _ = load_and_validate(args.manifest_dir)
            compare_captures(
                _load_capture(args.base), _load_capture(args.candidate), programs,
                _load_capture(args.strict_base) if args.strict_base else None,
                _load_capture(args.strict_candidate)
                if args.strict_candidate else None)
    except (CrosscheckError, ManifestError, subprocess.CalledProcessError) as error:
        print(f"pullback_crosscheck: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
