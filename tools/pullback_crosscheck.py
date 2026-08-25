#!/usr/bin/env python3
"""Produce, validate, and compare canonical pullback capture files."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
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


class StrictFpRequired(CrosscheckError):
    """Release captures differ and require strict-FP attribution."""


def _load_capture(path: Path) -> dict:
    try:
        capture = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise CrosscheckError(f"{path}: {error}") from error
    if not isinstance(capture, dict):
        raise CrosscheckError(f"{path}: capture root must be an object")
    return capture


def _file_sha256(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


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
            "operation",
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


def _selected_pixel(mapping: str, x: int, y: int, width: int, height: int) -> bool:
    if mapping == "COLUMN_ZERO":
        return x == 0
    if mapping == "WRAP_COLUMNS":
        return x in (0, width - 1)
    if mapping == "NORTH_ROW":
        return y == 0
    if mapping == "SOUTH_ROW":
        return y == height - 1
    if mapping == "POLE_BANDS":
        return y < 2 or y >= height - 2
    if mapping == "HORIZON_COLUMNS":
        return x in (width // 4, 3 * width // 4)
    if mapping == "OCTANT_COLUMNS":
        return x in (width // 8, 3 * width // 8, 5 * width // 8, 7 * width // 8)
    if mapping == "MIRROR_GRID":
        return x % (width // 8) == 0 or y in (
            height // 4,
            height // 2,
            3 * height // 4,
        )
    if mapping == "CARDINAL_POINTS":
        return y == height // 2 and x % (width // 4) == 0
    if mapping == "FRAME_PERIMETER":
        return x in (0, width - 1) or y in (0, height - 1)
    if mapping == "EQUATOR_ROW":
        return y == height // 2
    if mapping == "FRONT_AXIS_POINT":
        return x == 0 and y == height // 2
    raise CrosscheckError(f"unknown spatial probe mapping {mapping}")


def _validate_frame_operations(frames: dict, programs: dict) -> None:
    mappings = programs["corpus"]["probe_operations"]
    preset_count = sum(len(program["presets"])
                       for program in programs["programs"])
    spatial_hashes = set()
    for key, (frame, digest) in frames.items():
        _, preset, _, resolution, probe = key
        operation = frame["operation"]
        if not isinstance(operation, dict):
            raise CrosscheckError(f"frame operation must be an object: {key}")
        pixels = frame["pixels"]
        total = resolution[0] * resolution[1]
        if probe == "steady":
            if operation != {"kind": "parameter-case", "selected_pixels": total}:
                raise CrosscheckError(f"parameter-case operation mismatch: {key}")
            continue
        mapping = mappings[probe]
        if mapping in ("THROUGH_CLEAR_FROM", "THROUGH_CLEAR_TO"):
            endpoint = "from" if mapping.endswith("FROM") else "to"
            expected = {
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
            if operation != expected:
                raise CrosscheckError(f"transition endpoint operation mismatch: {key}")
            steady_key = (key[0], preset, "default", resolution, "steady")
            if digest != frames[steady_key][1]:
                raise CrosscheckError(
                    f"transition endpoint is not the steady endpoint: {key}"
                )
            continue
        selected = 0
        for index, pixel in enumerate(pixels):
            x = index % resolution[0]
            y = index // resolution[0]
            if _selected_pixel(mapping, x, y, *resolution):
                selected += 1
            elif any(pixel[:3]):
                raise CrosscheckError(
                    f"spatial probe contains pixels outside mapping: {key}"
                )
        expected = {
            "kind": "spatial-extract",
            "mapping": mapping,
            "selected_pixels": selected,
        }
        if operation != expected or not 0 < selected < total:
            raise CrosscheckError(f"spatial probe operation mismatch: {key}")
        spatial_hashes.add(digest)
    if len(spatial_hashes) < 2:
        raise CrosscheckError("spatial probe frames are aliases")


def expected_toolchain(programs: dict, configuration: str) -> dict:
    if configuration == "native-debug":
        pin = programs["toolchains"]["native"]
        keys = ("compiler", "cmake", "cmake_preset", "configuration")
        return {key: pin[key] for key in keys}
    pin = programs["toolchains"]["wasm"]
    if configuration == "wasm-strict-fp":
        return {
            "compiler": pin["compiler"],
            "cmake": pin["cmake"],
            "cmake_preset": pin["strict_cmake_preset"],
            "configuration": pin["strict_configuration"],
        }
    keys = ("compiler", "cmake", "cmake_preset", "configuration")
    return {key: pin[key] for key in keys}


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


def _within_differing_pixel_limit(
    differing: int, count: int, fraction: float
) -> bool:
    return differing <= math.ceil(count * fraction)


def _oracle_metric_map(capture: dict, oracles: list[dict]) -> dict:
    metrics = capture.get("oracle_metrics")
    if not isinstance(metrics, list):
        raise CrosscheckError("capture.oracle_metrics must be an array")
    mapped = {}
    for metric in metrics:
        required = {
            "oracle_id",
            "domain",
            "aggregation",
            "value",
            "unit",
            "sample_count",
            "resolution_values",
        }
        if not isinstance(metric, dict) or set(metric) != required:
            raise CrosscheckError("capture oracle metric fields are invalid")
        key = (metric["oracle_id"], metric["domain"], metric["aggregation"])
        if (
            key in mapped
            or key[1:] != ("FRAMEBUFFER", "MAXIMUM")
            or type(metric["value"]) is not int
            or metric["value"] < 0
            or type(metric["sample_count"]) is not int
            or metric["sample_count"] <= 0
            or not isinstance(metric["unit"], str)
            or not metric["unit"]
            or not isinstance(metric["resolution_values"], dict)
            or set(metric["resolution_values"]) != {"96x20", "288x144"}
            or any(
                type(value) is not int or value < 0
                for value in metric["resolution_values"].values()
            )
            or max(metric["resolution_values"].values()) != metric["value"]
        ):
            raise CrosscheckError("capture oracle metric is invalid")
        mapped[key] = metric
    expected = {}
    for oracle in oracles:
        manifest_metric = next(
            metric
            for metric in oracle["metrics"]
            if metric["domain"] == "FRAMEBUFFER"
            and metric["aggregation"] == "MAXIMUM"
        )
        key = (oracle["oracle_id"], "FRAMEBUFFER", "MAXIMUM")
        expected[key] = manifest_metric
    if mapped.keys() != expected.keys():
        raise CrosscheckError("capture oracle metric coverage differs from manifests")
    for key, metric in mapped.items():
        manifest_metric = expected[key]
        if (
            metric["value"]
            != manifest_metric["configuration_baselines"][capture["configuration"]]
            or metric["unit"] != manifest_metric["unit"]
            or metric["resolution_values"]
            != manifest_metric["resolution_baselines"][capture["configuration"]]
        ):
            raise CrosscheckError("capture oracle metric differs from manifest")
    return mapped


def _validate_capture(
    capture: dict,
    programs: dict,
    digest: str,
    configuration: str,
    checkout_sha: str,
    oracles: list[dict],
) -> tuple[dict, dict]:
    required = {
        "schema_version",
        "checkout_sha",
        "manifest_sha256",
        "toolchain",
        "configuration",
        "corpus",
        "oracle_metrics",
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
    if capture["toolchain"] != expected_toolchain(programs, configuration):
        raise CrosscheckError("observed capture toolchain differs from manifest pin")
    if capture["configuration"] != configuration:
        raise CrosscheckError("capture configuration mismatch")
    if capture["corpus"] != programs["corpus"]:
        raise CrosscheckError("capture corpus differs from manifest")
    frames = _frame_map(capture)
    if frames.keys() != _expected_frame_keys(programs):
        raise CrosscheckError("capture does not cover every manifest case and probe")
    _validate_frame_operations(frames, programs)
    return frames, _oracle_metric_map(capture, oracles)


def compare_captures(
    base: dict,
    candidate: dict,
    programs: dict,
    digest: str,
    base_sha: str,
    candidate_sha: str,
    strict_base: dict | None = None,
    strict_candidate: dict | None = None,
    *,
    oracles: list[dict],
) -> None:
    if base_sha != programs["base_sha"]:
        raise CrosscheckError("base checkout SHA differs from manifest pin")
    if SHA_RE.fullmatch(candidate_sha) is None:
        raise CrosscheckError("candidate checkout SHA must be a full lowercase SHA")
    configuration = base.get("configuration")
    if configuration not in ("native-debug", "wasm-release"):
        raise CrosscheckError(f"unsupported configuration {configuration}")
    base_frames, base_oracles = _validate_capture(
        base, programs, digest, configuration, base_sha, oracles
    )
    candidate_frames, candidate_oracles = _validate_capture(
        candidate, programs, digest, configuration, candidate_sha, oracles
    )
    if base_oracles != candidate_oracles:
        raise CrosscheckError("capture oracle metrics differ")
    strict_frames = None
    if strict_base is not None or strict_candidate is not None:
        if strict_base is None or strict_candidate is None:
            raise CrosscheckError("strict-FP captures must be supplied as a pair")
        strict_captures = (
            _validate_capture(
                strict_base, programs, digest, "wasm-strict-fp", base_sha, oracles
            ),
            _validate_capture(
                strict_candidate,
                programs,
                digest,
                "wasm-strict-fp",
                candidate_sha,
                oracles,
            ),
        )
        strict_frames = (strict_captures[0][0], strict_captures[1][0])
        if strict_captures[0][1] != strict_captures[1][1]:
            raise CrosscheckError("strict-FP oracle metrics differ")
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
            raise StrictFpRequired("release mismatch requires strict-FP captures")
        for key, maximum, differing, count, limits in release_differences:
            if maximum > limits["max_channel_delta_u16"]:
                raise CrosscheckError(f"release channel delta exceeds limit: {key}")
            if not _within_differing_pixel_limit(
                differing, count, limits["max_differing_pixel_fraction"]
            ):
                raise CrosscheckError(
                    f"release differing-pixel share exceeds limit: {key}"
                )


def _compare_capture_paths(
    base_path: Path,
    candidate_path: Path,
    programs: dict,
    digest: str,
    base_sha: str,
    candidate_sha: str,
    strict_base_path: Path | None = None,
    strict_candidate_path: Path | None = None,
    *,
    oracles: list[dict],
) -> tuple[bool, dict[str, dict]]:
    """Compare captures loaded for this call; return strict request and toolchains."""
    base = _load_capture(base_path)
    candidate = _load_capture(candidate_path)
    strict_base = _load_capture(strict_base_path) if strict_base_path else None
    strict_candidate = (
        _load_capture(strict_candidate_path) if strict_candidate_path else None
    )
    try:
        compare_captures(
            base,
            candidate,
            programs,
            digest,
            base_sha,
            candidate_sha,
            strict_base,
            strict_candidate,
            oracles=oracles,
        )
    except StrictFpRequired:
        return True, {base["configuration"]: dict(base["toolchain"])}
    observed_toolchains = {base["configuration"]: dict(base["toolchain"])}
    if strict_base is not None:
        observed_toolchains[strict_base["configuration"]] = dict(
            strict_base["toolchain"]
        )
    return False, observed_toolchains


def _run(command: list[str], cwd: Path, environment: dict) -> None:
    subprocess.run(
        command, cwd=cwd, env=environment, check=True, timeout=1800)


def _remove_worktree(
    repository: Path,
    worktree: Path,
    environment: dict,
    active_error: BaseException | None,
) -> None:
    try:
        _run(
            ["git", "-C", str(repository), "worktree", "remove", str(worktree)],
            repository,
            environment,
        )
    except (OSError, subprocess.SubprocessError) as cleanup_error:
        if active_error is None:
            raise
        active_error.add_note(f"worktree cleanup also failed: {cleanup_error}")


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
    preset = {
        "native-debug": "tests",
        "wasm-release": "wasm-release",
        "wasm-strict-fp": "wasm-strict-fp",
    }[configuration]
    if configuration != "native-debug" and not environment.get("EMSDK"):
        raise CrosscheckError("EMSDK is required for WASM capture builds")
    command = [
        "cmake",
        "--preset",
        preset,
        "-B",
        str(build),
        f"-DHS_PULLBACK_CAPTURE_SOURCE_ROOT={source}",
        "-DHS_INSTALL_GIT_HOOKS=OFF",
        "-DCMAKE_CXX_COMPILER_LAUNCHER=",
    ]
    _run(command, controller, environment)
    target = (
        "pullback_capture_native"
        if configuration == "native-debug"
        else "pullback_capture_wasm"
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
    if base_sha == candidate_sha:
        raise CrosscheckError("base and candidate revisions must differ")
    programs, oracles, schema = load_and_validate(controller / "tests/data/pullback")
    digest = manifest_sha256(programs, oracles, schema)
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
                capture_paths = {}
                observed_toolchains = {}
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
                        capture_paths[(name, configuration)] = path
                    strict_required, compared_toolchains = _compare_capture_paths(
                        capture_paths[("base", configuration)],
                        capture_paths[("candidate", configuration)],
                        programs,
                        digest,
                        base_sha,
                        candidate_sha,
                        oracles=oracles,
                    )
                    observed_toolchains.update(compared_toolchains)
                    if configuration == "native-debug" and strict_required:
                        raise CrosscheckError(
                            "native comparison requested strict-FP captures"
                        )
                if strict_required:
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
                        capture_paths[(name, configuration)] = path
                    strict_required, compared_toolchains = _compare_capture_paths(
                        capture_paths[("base", "wasm-release")],
                        capture_paths[("candidate", "wasm-release")],
                        programs,
                        digest,
                        base_sha,
                        candidate_sha,
                        capture_paths[("base", "wasm-strict-fp")],
                        capture_paths[("candidate", "wasm-strict-fp")],
                        oracles=oracles,
                    )
                    observed_toolchains.update(compared_toolchains)
                    if strict_required:
                        raise CrosscheckError(
                            "strict-FP comparison requested strict-FP captures"
                        )
                summary = {
                    "schema_version": 1,
                    "base_sha": base_sha,
                    "candidate_sha": candidate_sha,
                    "manifest_sha256": digest,
                    "toolchains": programs["toolchains"],
                    "observed_toolchains": observed_toolchains,
                    "configurations": sorted(
                        {configuration for _, configuration in capture_paths}
                    ),
                    "capture_sha256": {
                        path.name: _file_sha256(path)
                        for path in capture_paths.values()
                    },
                }
                (output / "summary.json").write_text(
                    json.dumps(summary, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8",
                )
            finally:
                _remove_worktree(
                    repository,
                    sources["candidate"],
                    environment,
                    sys.exception(),
                )
        finally:
            _remove_worktree(
                repository,
                sources["base"],
                environment,
                sys.exception(),
            )


def main(argv=None) -> int:
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
    args = parser.parse_args(argv)
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
            programs, oracles, schema = load_and_validate(args.manifest_dir)
            compare_captures(
                _load_capture(args.base),
                _load_capture(args.candidate),
                programs,
                manifest_sha256(programs, oracles, schema),
                args.base_sha,
                args.candidate_sha,
                _load_capture(args.strict_base) if args.strict_base else None,
                _load_capture(args.strict_candidate) if args.strict_candidate else None,
                oracles=oracles,
            )
    except (
        CrosscheckError,
        ManifestError,
        OSError,
        subprocess.SubprocessError,
    ) as error:
        print(f"pullback_crosscheck: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
