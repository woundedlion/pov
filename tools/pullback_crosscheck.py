#!/usr/bin/env python3
"""Produce, validate, and compare canonical pullback capture files."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import tempfile
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
    if (
        not isinstance(resolution, list)
        or len(resolution) != 2
        or not all(isinstance(value, int) and value > 0 for value in resolution)
    ):
        raise CrosscheckError("frame resolution must contain two positive integers")
    pixels = frame.get("pixels")
    expected_count = resolution[0] * resolution[1]
    if not isinstance(pixels, list) or len(pixels) != expected_count:
        raise CrosscheckError(
            f"frame has {len(pixels) if isinstance(pixels, list) else 'invalid'} "
            f"pixels; expected {expected_count}"
        )
    raw = bytearray()
    for pixel in pixels:
        if not isinstance(pixel, list) or len(pixel) != 4:
            raise CrosscheckError("pixels must contain four 16-bit channels")
        for channel in pixel:
            if (
                not isinstance(channel, int)
                or isinstance(channel, bool)
                or not (0 <= channel <= 65535)
            ):
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
            "program",
            "preset",
            "case",
            "resolution",
            "probe",
            "sha256",
            "pixels",
        }:
            raise CrosscheckError("capture frame fields are incomplete or unknown")
        resolution = tuple(frame.get("resolution", ()))
        key = (
            frame.get("program"),
            frame.get("preset"),
            frame.get("case"),
            resolution,
            frame.get("probe"),
        )
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
        for preset in program["presets"]:
            for resolution in resolutions:
                dimensions = tuple(resolution)
                for case in program["parameter_cases"]:
                    expected.add((program_id, preset, case["id"], dimensions, "steady"))
                for category, probes in program["probes"].items():
                    for probe in probes:
                        expected.add(
                            (
                                program_id,
                                preset,
                                "default",
                                dimensions,
                                f"{category}:{probe}",
                            )
                        )
    return expected


def _compare_pixels(base: list, candidate: list) -> tuple[int, int, int]:
    if len(base) != len(candidate):
        raise CrosscheckError("frame pixel counts differ")
    maximum = 0
    differing = 0
    for base_pixel, candidate_pixel in zip(base, candidate, strict=True):
        if len(base_pixel) != 4 or len(candidate_pixel) != 4:
            raise CrosscheckError("pixels must contain four 16-bit channels")
        deltas = [abs(a - b) for a, b in zip(base_pixel, candidate_pixel, strict=True)]
        maximum = max(maximum, *deltas)
        differing += any(deltas)
    return maximum, differing, len(base)


def _validate_capture(
    capture: dict, programs: dict, digest: str, configuration: str, checkout_sha: str
) -> dict:
    required = {
        "schema_version",
        "checkout_sha",
        "manifest_sha256",
        "toolchains",
        "configuration",
        "corpus",
        "frames",
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


def compare_captures(
    base: dict,
    candidate: dict,
    programs: dict,
    digest: str,
    base_sha: str,
    candidate_sha: str,
    strict_base: dict | None = None,
    strict_candidate: dict | None = None,
) -> None:
    if base_sha != programs["base_sha"]:
        raise CrosscheckError("base checkout SHA differs from manifest pin")
    if SHA_RE.fullmatch(candidate_sha) is None:
        raise CrosscheckError("candidate checkout SHA must be a full lowercase SHA")
    configuration = base.get("configuration")
    if configuration not in ("native-debug", "wasm-release"):
        raise CrosscheckError(f"unsupported configuration {configuration}")
    base_frames = _validate_capture(base, programs, digest, configuration, base_sha)
    candidate_frames = _validate_capture(
        candidate, programs, digest, configuration, candidate_sha
    )
    strict_frames = None
    if strict_base is not None or strict_candidate is not None:
        if strict_base is None or strict_candidate is None:
            raise CrosscheckError("strict-FP captures must be supplied as a pair")
        strict_frames = (
            _validate_capture(
                strict_base, programs, digest, "wasm-strict-fp", base_sha
            ),
            _validate_capture(
                strict_candidate, programs, digest, "wasm-strict-fp", candidate_sha
            ),
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
            base_frame["pixels"], candidate_frame["pixels"]
        )
        limits = program_tolerances[key[0]]["release_wasm_framebuffer"]
        release_differences.append((key, maximum, differing, count, limits))
    if release_differences:
        if strict_frames is None:
            raise CrosscheckError("release mismatch requires strict-FP captures")
        for key, maximum, differing, count, limits in release_differences:
            if maximum > limits["max_channel_delta_u16"]:
                raise CrosscheckError(f"release channel delta exceeds limit: {key}")
            if differing / count > limits["max_differing_pixel_fraction"]:
                raise CrosscheckError(
                    f"release differing-pixel share exceeds limit: {key}"
                )


def _run(command: list[str], cwd: Path, environment: dict) -> None:
    subprocess.run(command, cwd=cwd, env=environment, check=True)


def validate_build_isolation(
    sources: dict[str, Path], builds: dict[str, dict[str, Path]]
) -> None:
    source_paths = {name: path.resolve() for name, path in sources.items()}
    if len(set(source_paths.values())) != len(source_paths):
        raise CrosscheckError("capture source checkouts are not isolated")
    build_paths = [
        path.resolve()
        for configurations in builds.values()
        for path in configurations.values()
    ]
    if len(set(build_paths)) != len(build_paths):
        raise CrosscheckError("capture build directories are not isolated")
    for build in build_paths:
        for source in source_paths.values():
            if build == source or build in source.parents or source in build.parents:
                raise CrosscheckError("capture build and source checkout are nested")


def _build_capture_backend(
    controller: Path, source: Path, build: Path, configuration: str, environment: dict
) -> None:
    build_type = {
        "native-debug": "Debug",
        "wasm-release": "Release",
        "wasm-strict-fp": "PullbackStrictFP",
    }[configuration]
    if configuration == "native-debug":
        toolchain = controller / "cmake/toolchain-native-clang.cmake"
    else:
        emsdk = environment.get("EMSDK")
        if not emsdk:
            raise CrosscheckError("EMSDK is required for WASM capture builds")
        toolchain = (
            Path(emsdk) / "upstream/emscripten/cmake/Modules/Platform/Emscripten.cmake"
        )
    command = [
        "cmake",
        "-S",
        str(controller),
        "-B",
        str(build),
        "-G",
        "Ninja",
        f"-DCMAKE_BUILD_TYPE={build_type}",
        f"-DCMAKE_TOOLCHAIN_FILE={toolchain}",
        f"-DHS_PULLBACK_CAPTURE_SOURCE_ROOT={source}",
        "-DHS_INSTALL_GIT_HOOKS=OFF",
        "-DCMAKE_CXX_COMPILER_LAUNCHER=",
    ]
    _run(command, controller, environment)
    target = (
        "pullback_capture_native"
        if configuration == "native-debug"
        else "holosphere_wasm"
    )
    _run(["cmake", "--build", str(build), "--target", target], controller, environment)


def _produce_capture(
    controller: Path,
    source: Path,
    build: Path,
    manifest_dir: Path,
    configuration: str,
    output: Path,
    environment: dict,
) -> None:
    _run(
        [
            sys.executable,
            str(controller / "tools/pullback_capture.py"),
            "--configuration",
            configuration,
            "--checkout",
            str(source),
            "--build-dir",
            str(build),
            "--producer-root",
            str(controller),
            "--manifest-dir",
            str(manifest_dir),
            "--output",
            str(output),
        ],
        controller,
        environment,
    )


def orchestrate(
    repository: Path, controller: Path, base_sha: str, candidate_sha: str, output: Path
) -> None:
    programs, oracles = load_and_validate(controller / "tests/data/pullback")
    digest = manifest_sha256(programs, oracles)
    environment = os.environ.copy()
    environment["CCACHE_DISABLE"] = "1"
    environment["CMAKE_CXX_COMPILER_LAUNCHER"] = ""
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="pullback-crosscheck-") as directory:
        root = Path(directory)
        sources = {"base": root / "base", "candidate": root / "candidate"}
        builds = {
            name: {
                configuration: root / f"build-{name}-{configuration}"
                for configuration in ("native-debug", "wasm-release", "wasm-strict-fp")
            }
            for name in sources
        }
        validate_build_isolation(sources, builds)
        _run(
            [
                "git",
                "-C",
                str(repository),
                "worktree",
                "add",
                "--detach",
                str(sources["base"]),
                base_sha,
            ],
            repository,
            environment,
        )
        try:
            _run(
                [
                    "git",
                    "-C",
                    str(repository),
                    "worktree",
                    "add",
                    "--detach",
                    str(sources["candidate"]),
                    candidate_sha,
                ],
                repository,
                environment,
            )
            try:
                captures = {}
                capture_paths = {}
                for configuration in ("native-debug", "wasm-release"):
                    for name in ("base", "candidate"):
                        _build_capture_backend(
                            controller,
                            sources[name],
                            builds[name][configuration],
                            configuration,
                            environment,
                        )
                        path = output / f"{name}-{configuration}.json"
                        _produce_capture(
                            controller,
                            sources[name],
                            builds[name][configuration],
                            controller / "tests/data/pullback",
                            configuration,
                            path,
                            environment,
                        )
                        captures[(name, configuration)] = _load_capture(path)
                        capture_paths[(name, configuration)] = path
                compare_captures(
                    captures[("base", "native-debug")],
                    captures[("candidate", "native-debug")],
                    programs,
                    digest,
                    base_sha,
                    candidate_sha,
                )
                try:
                    compare_captures(
                        captures[("base", "wasm-release")],
                        captures[("candidate", "wasm-release")],
                        programs,
                        digest,
                        base_sha,
                        candidate_sha,
                    )
                except CrosscheckError as error:
                    if "requires strict-FP" not in str(error):
                        raise
                    for name in ("base", "candidate"):
                        configuration = "wasm-strict-fp"
                        _build_capture_backend(
                            controller,
                            sources[name],
                            builds[name][configuration],
                            configuration,
                            environment,
                        )
                        path = output / f"{name}-{configuration}.json"
                        _produce_capture(
                            controller,
                            sources[name],
                            builds[name][configuration],
                            controller / "tests/data/pullback",
                            configuration,
                            path,
                            environment,
                        )
                        captures[(name, configuration)] = _load_capture(path)
                        capture_paths[(name, configuration)] = path
                    compare_captures(
                        captures[("base", "wasm-release")],
                        captures[("candidate", "wasm-release")],
                        programs,
                        digest,
                        base_sha,
                        candidate_sha,
                        captures[("base", "wasm-strict-fp")],
                        captures[("candidate", "wasm-strict-fp")],
                    )
                summary = {
                    "schema_version": 1,
                    "base_sha": base_sha,
                    "candidate_sha": candidate_sha,
                    "manifest_sha256": digest,
                    "toolchains": programs["toolchains"],
                    "configurations": sorted(
                        {configuration for _, configuration in capture_paths}
                    ),
                    "capture_sha256": {
                        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                        for path in capture_paths.values()
                    },
                }
                (output / "summary.json").write_text(
                    json.dumps(summary, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8",
                )
            finally:
                _run(
                    [
                        "git",
                        "-C",
                        str(repository),
                        "worktree",
                        "remove",
                        str(sources["candidate"]),
                    ],
                    repository,
                    environment,
                )
        finally:
            _run(
                [
                    "git",
                    "-C",
                    str(repository),
                    "worktree",
                    "remove",
                    str(sources["base"]),
                ],
                repository,
                environment,
            )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    compare = commands.add_parser("compare")
    compare.add_argument("--manifest-dir", type=Path, required=True)
    compare.add_argument("--base", type=Path, required=True)
    compare.add_argument("--candidate", type=Path, required=True)
    compare.add_argument("--base-sha", required=True)
    compare.add_argument("--candidate-sha", required=True)
    compare.add_argument("--strict-base", type=Path)
    compare.add_argument("--strict-candidate", type=Path)
    run = commands.add_parser("run")
    run.add_argument("--repository", type=Path, required=True)
    run.add_argument(
        "--controller", type=Path, default=Path(__file__).resolve().parents[1]
    )
    run.add_argument("--base-sha", required=True)
    run.add_argument("--candidate-sha", required=True)
    run.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command == "run":
            orchestrate(
                args.repository.resolve(),
                args.controller.resolve(),
                args.base_sha,
                args.candidate_sha,
                args.output.resolve(),
            )
        else:
            programs, oracles = load_and_validate(args.manifest_dir)
            compare_captures(
                _load_capture(args.base),
                _load_capture(args.candidate),
                programs,
                manifest_sha256(programs, oracles),
                args.base_sha,
                args.candidate_sha,
                _load_capture(args.strict_base) if args.strict_base else None,
                _load_capture(args.strict_candidate) if args.strict_candidate else None,
            )
    except (
        CrosscheckError,
        ManifestError,
        OSError,
        subprocess.CalledProcessError,
    ) as error:
        print(f"pullback_crosscheck: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
