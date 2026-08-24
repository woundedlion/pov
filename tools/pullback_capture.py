#!/usr/bin/env python3
"""Produce canonical pullback captures from native or WASM ShaderWorkbench builds."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import struct
import subprocess
import tempfile
from pathlib import Path

from generate_pullback_manifest_header import (
    ManifestError,
    OPERATION_CODES,
    load_and_validate,
    manifest_sha256,
)
from pullback_crosscheck import _expected_toolchain


CASE_OPERATIONS = {
    "default": "CASE_DEFAULT",
    "endpoint_min": "CASE_ENDPOINT_MIN",
    "endpoint_max": "CASE_ENDPOINT_MAX",
    "interior": "CASE_INTERIOR",
}


class CaptureError(RuntimeError):
    """A capture backend, corpus, or toolchain attestation is invalid."""


def _read_u16(data: bytes, offset: int) -> tuple[int, int]:
    if offset + 2 > len(data):
        raise CaptureError("capture backend stream is truncated")
    return struct.unpack_from("<H", data, offset)[0], offset + 2


def _read_u32(data: bytes, offset: int) -> tuple[int, int]:
    low, offset = _read_u16(data, offset)
    high, offset = _read_u16(data, offset)
    return low | high << 16, offset


def operation_specs(programs: dict) -> list[dict]:
    mappings = programs["corpus"]["probe_operations"]
    specs = []
    for program in programs["programs"]:
        for preset in program["presets"]:
            for case in program["parameter_cases"]:
                case_id = case["id"]
                specs.append(
                    {
                        "preset": preset,
                        "name": f"case:{case_id}",
                        "mapping": CASE_OPERATIONS[case_id],
                    }
                )
            for category, probes in program["probes"].items():
                for probe in probes:
                    name = f"{category}:{probe}"
                    specs.append(
                        {"preset": preset, "name": name, "mapping": mappings[name]}
                    )
    keys = [(spec["preset"], spec["name"]) for spec in specs]
    if len(keys) != len(set(keys)):
        raise CaptureError("capture corpus operations are not unique")
    return specs


def oracle_operation_specs(oracles: list[dict]) -> list[dict]:
    specs = []
    for oracle in oracles:
        for probe in oracle["corpus"]["framebuffer"]["probes"]:
            for operation in probe["operations"]:
                specs.append(
                    {
                        "oracle": oracle["oracle_id"],
                        "preset": probe["preset"],
                        "hue_noise_phase": probe["hue_noise_phase"],
                        "mapping": operation,
                    }
                )
    keys = [
        (spec["oracle"], spec["preset"], spec["mapping"]) for spec in specs
    ]
    if len(keys) != len(set(keys)):
        raise CaptureError("oracle capture operations are not unique")
    return specs


def write_operations(programs: dict, oracles: list[dict], path: Path) -> None:
    specs = operation_specs(programs)
    oracle_specs = oracle_operation_specs(oracles)
    data = bytearray(b"HSPO")
    data.extend(struct.pack("<HHH", 2, len(specs), len(oracle_specs)))
    for spec in specs:
        name = spec["name"].encode("utf-8")
        mapping = spec["mapping"]
        if mapping not in OPERATION_CODES:
            raise CaptureError(f"unknown capture operation mapping {mapping}")
        data.extend(
            struct.pack("<HHH", spec["preset"], OPERATION_CODES[mapping], len(name))
        )
        data.extend(name)
    for spec in oracle_specs:
        name = spec["oracle"].encode("utf-8")
        mapping = spec["mapping"]
        if mapping not in OPERATION_CODES:
            raise CaptureError(f"unknown oracle operation mapping {mapping}")
        data.extend(struct.pack("<H", len(name)))
        data.extend(name)
        data.extend(
            struct.pack(
                "<HHf",
                spec["preset"],
                OPERATION_CODES[mapping],
                spec["hue_noise_phase"],
            )
        )
    path.write_bytes(data)


def load_backend(
    path: Path, programs: dict, oracles: list[dict]
) -> tuple[tuple[int, int], dict, dict]:
    data = path.read_bytes()
    if data[:4] != b"HSPB":
        raise CaptureError("capture backend magic mismatch")
    offset = 4
    version, offset = _read_u16(data, offset)
    width, offset = _read_u16(data, offset)
    height, offset = _read_u16(data, offset)
    record_count, offset = _read_u32(data, offset)
    oracle_count, offset = _read_u16(data, offset)
    specs = {(spec["preset"], spec["name"]): spec for spec in operation_specs(programs)}
    if version != 3 or record_count != len(specs):
        raise CaptureError("capture backend header mismatch")
    records = {}
    for _ in range(record_count):
        preset, offset = _read_u16(data, offset)
        operation_code, offset = _read_u16(data, offset)
        name_length, offset = _read_u16(data, offset)
        if offset + name_length > len(data):
            raise CaptureError("capture backend operation name is truncated")
        try:
            name = data[offset : offset + name_length].decode("utf-8")
        except UnicodeDecodeError as error:
            raise CaptureError("capture backend operation name is invalid") from error
        offset += name_length
        source, offset = _read_u16(data, offset)
        destination, offset = _read_u16(data, offset)
        elapsed, offset = _read_u16(data, offset)
        duration, offset = _read_u16(data, offset)
        selected_pixels, offset = _read_u32(data, offset)
        key = (preset, name)
        if key in records or key not in specs:
            raise CaptureError(f"invalid capture backend record {key}")
        expected_code = OPERATION_CODES[specs[key]["mapping"]]
        if operation_code != expected_code:
            raise CaptureError(f"capture backend operation mismatch for {key}")
        end = offset + width * height * 6
        if end > len(data):
            raise CaptureError("capture backend stream is truncated")
        # Corpus pixels are the backend's three little-endian 16-bit channels
        # plus an opaque alpha, held as that byte string rather than a list per
        # pixel: it is what the frame hash is taken over, and a full-resolution
        # roster of four-element lists costs gigabytes.
        channels = data[offset:end]
        offset = end
        pixels = bytearray(width * height * 8)
        for byte in range(6):
            pixels[byte::8] = channels[byte::6]
        pixels[6::8] = b"\xff" * (width * height)
        pixels[7::8] = b"\xff" * (width * height)
        records[key] = {
            "pixels": bytes(pixels),
            "mapping": specs[key]["mapping"],
            "selected_pixels": selected_pixels,
            "source": source,
            "destination": destination,
            "elapsed": elapsed,
            "duration": duration,
        }
    metrics = {}
    expected_oracles = {oracle["oracle_id"] for oracle in oracles}
    for _ in range(oracle_count):
        name_length, offset = _read_u16(data, offset)
        if offset + name_length > len(data):
            raise CaptureError("capture backend oracle name is truncated")
        try:
            oracle = data[offset : offset + name_length].decode("utf-8")
        except UnicodeDecodeError as error:
            raise CaptureError("capture backend oracle name is invalid") from error
        offset += name_length
        maximum, offset = _read_u16(data, offset)
        samples, offset = _read_u32(data, offset)
        if oracle in metrics or oracle not in expected_oracles or samples == 0:
            raise CaptureError(f"invalid capture backend oracle metric {oracle}")
        metrics[oracle] = {"value": maximum, "sample_count": samples}
    if (
        offset != len(data)
        or records.keys() != specs.keys()
        or metrics.keys() != expected_oracles
    ):
        raise CaptureError("capture backend stream is incomplete or has trailing bytes")
    return (width, height), records, metrics


def validate_backend_oracles(
    metrics: dict,
    oracles: list[dict],
    configuration: str,
    resolution: tuple[int, int] | None = None,
) -> None:
    for oracle in oracles:
        oracle_id = oracle["oracle_id"]
        manifest_metric = next(
            metric
            for metric in oracle["metrics"]
            if metric["domain"] == "FRAMEBUFFER"
            and metric["aggregation"] == "MAXIMUM"
        )
        observed = metrics[oracle_id]["value"]
        limit = manifest_metric["accepted_limit"]
        if observed > limit:
            raise CaptureError(
                f"{oracle_id} framebuffer error exceeds its accepted limit: "
                f"{observed} > {limit}"
            )
        dimensions = "aggregate" if resolution is None else f"{resolution[0]}x{resolution[1]}"
        print(
            f"pullback metric: {configuration} {dimensions} {oracle_id} "
            f"observed={observed} limit={limit}"
        )


def _cache_values(build_dir: Path) -> dict[str, str]:
    cache = build_dir / "CMakeCache.txt"
    try:
        lines = cache.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        raise CaptureError(f"cannot read CMake cache: {cache}") from error
    values = {}
    for line in lines:
        if not line or line.startswith(("#", "//")) or "=" not in line:
            continue
        name_type, value = line.split("=", 1)
        name = name_type.split(":", 1)[0]
        values[name] = value
    return values


def _compiler_metadata(build_dir: Path) -> tuple[Path, str]:
    files = list((build_dir / "CMakeFiles").glob("*/CMakeCXXCompiler.cmake"))
    if len(files) != 1:
        raise CaptureError("cannot identify CMake compiler metadata")
    contents = files[0].read_text(encoding="utf-8")
    path_match = re.search(r'set\(CMAKE_CXX_COMPILER "([^"]+)"\)', contents)
    version_match = re.search(
        r'set\(CMAKE_CXX_COMPILER_VERSION "([0-9]+\.[0-9]+\.[0-9]+)"\)',
        contents,
    )
    if path_match is None or version_match is None:
        raise CaptureError("CMake compiler metadata is incomplete")
    return Path(path_match.group(1)), version_match.group(1)


def _first_line(command: list[str]) -> str:
    try:
        result = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=30,
        )
    except (OSError, subprocess.SubprocessError) as error:
        raise CaptureError(f"toolchain query failed: {command[0]}") from error
    lines = result.stdout.splitlines()
    if not lines:
        raise CaptureError(f"toolchain query produced no version: {command[0]}")
    return lines[0]


def attest_toolchain(build_dir: Path, configuration: str, programs: dict) -> dict:
    cache = _cache_values(build_dir)
    compiler, clang_version = _compiler_metadata(build_dir)
    cmake = cache.get("CMAKE_COMMAND", "cmake")
    compiler_line = _first_line([str(compiler), "--version"])
    cmake_line = _first_line([cmake, "--version"])
    cmake_match = re.fullmatch(r"cmake version ([0-9]+\.[0-9]+\.[0-9]+)", cmake_line)
    if cmake_match is None:
        raise CaptureError("CMake version output is unrecognized")
    if configuration == "native-debug":
        if (
            re.search(
                rf"clang version {re.escape(clang_version)}(?:git)?\b", compiler_line
            )
            is None
        ):
            raise CaptureError("native compiler version disagrees with CMake metadata")
        compiler_version = f"clang {clang_version}"
    else:
        emscripten_match = re.search(r"\) ([0-9]+\.[0-9]+\.[0-9]+) \(", compiler_line)
        if emscripten_match is None:
            raise CaptureError("Emscripten compiler version output is unrecognized")
        compiler_version = (
            f"emscripten {emscripten_match.group(1)} clang {clang_version}"
        )
    observed = {
        "compiler": compiler_version,
        "cmake": cmake_match.group(1),
        "cmake_preset": cache.get("HS_PULLBACK_CAPTURE_PRESET", ""),
        "configuration": cache.get("CMAKE_BUILD_TYPE", ""),
    }
    expected = _expected_toolchain(programs, configuration)
    if observed != expected:
        raise CaptureError(
            f"observed toolchain differs from manifest pin: {observed} != {expected}"
        )
    return observed


def _backend_command(
    configuration: str,
    build_dir: Path,
    resolution: tuple[int, int],
    operations: Path,
    output: Path,
) -> list[str]:
    dimensions = f"{resolution[0]}x{resolution[1]}"
    arguments = [
        "--resolution",
        dimensions,
        "--operations",
        str(operations),
        str(output),
    ]
    if configuration == "native-debug":
        executable = build_dir / "tests/pullback_capture_native"
        if not executable.exists():
            executable = executable.with_suffix(".exe")
        return [str(executable), *arguments]
    return ["node", str(build_dir / "pullback_capture_wasm.cjs"), *arguments]


def _frame_operation(record: dict) -> dict:
    mapping = record["mapping"]
    if mapping.startswith("CASE_"):
        return {"kind": "parameter-case", "selected_pixels": record["selected_pixels"]}
    if mapping == "THROUGH_CLEAR_FROM" or mapping == "THROUGH_CLEAR_TO":
        return {
            "kind": "through-clear",
            "endpoint": "from" if mapping.endswith("FROM") else "to",
            "source_preset": record["source"],
            "destination_preset": record["destination"],
            "elapsed": record["elapsed"],
            "duration": record["duration"],
            "selected_pixels": record["selected_pixels"],
        }
    return {
        "kind": "spatial-extract",
        "mapping": mapping,
        "selected_pixels": record["selected_pixels"],
    }


def _pixels_json(pixels: bytes) -> str:
    """The corpus byte string as its `[[r,g,b,a],...]` JSON array."""
    return "[" + ",".join(
        f"[{r},{g},{b},{a}]" for r, g, b, a in struct.iter_unpack("<HHHH", pixels)
    ) + "]"


def write_capture(capture: dict, frames: list[dict], path: Path) -> None:
    """Write the capture JSON, expanding each frame's pixels as it streams.

    Byte-identical to `json.dumps(capture | {"frames": frames}, separators=
    (",", ":"))` with each frame's pixel bytes written as four-channel lists.
    The corpus is a quarter of a gigabyte of text over ten million pixels, so
    neither those lists nor the document are ever materialized.
    """
    compact = {"separators": (",", ":")}
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as stream:
        head = json.dumps(capture, **compact)
        stream.write(head[:-1] + ("," if capture else ""))
        stream.write('"frames":[')
        for index, frame in enumerate(frames):
            if index:
                stream.write(",")
            stream.write("{" + ",".join(
                json.dumps(key, **compact) + ":"
                + (_pixels_json(value) if key == "pixels"
                   else json.dumps(value, **compact))
                for key, value in frame.items()) + "}")
        stream.write("]}\n")


def produce(
    configuration: str,
    checkout: Path,
    build_dir: Path,
    manifest_dir: Path,
    output: Path,
) -> None:
    programs, oracles, schema = load_and_validate(manifest_dir)
    toolchain = attest_toolchain(build_dir, configuration, programs)
    checkout_sha = subprocess.run(
        ["git", "-C", str(checkout), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=30,
    ).stdout.strip()
    backend_records = {}
    oracle_totals = {
        oracle["oracle_id"]: {
            "value": 0,
            "sample_count": 0,
            "resolution_values": {},
        }
        for oracle in oracles
    }
    with tempfile.TemporaryDirectory(prefix="pullback-capture-") as directory:
        operations = Path(directory) / "operations.bin"
        write_operations(programs, oracles, operations)
        for resolution in programs["corpus"]["resolutions"]:
            dimensions = tuple(resolution)
            backend_output = Path(directory) / f"{dimensions[0]}x{dimensions[1]}.bin"
            subprocess.run(
                _backend_command(
                    configuration, build_dir, dimensions, operations, backend_output
                ),
                cwd=checkout,
                check=True,
                timeout=900,
            )
            actual_resolution, records, oracle_metrics = load_backend(
                backend_output, programs, oracles
            )
            if actual_resolution != dimensions:
                raise CaptureError("capture backend resolution mismatch")
            validate_backend_oracles(
                oracle_metrics, oracles, configuration, dimensions
            )
            for key, record in records.items():
                backend_records[(dimensions, *key)] = record
            for oracle, metric in oracle_metrics.items():
                oracle_totals[oracle]["value"] = max(
                    oracle_totals[oracle]["value"], metric["value"]
                )
                oracle_totals[oracle]["sample_count"] += metric["sample_count"]
                oracle_totals[oracle]["resolution_values"][
                    f"{dimensions[0]}x{dimensions[1]}"
                ] = metric["value"]
    frames = []
    for program in programs["programs"]:
        for preset in program["presets"]:
            for resolution in programs["corpus"]["resolutions"]:
                dimensions = tuple(resolution)
                for case in program["parameter_cases"]:
                    case_id = case["id"]
                    record = backend_records[(dimensions, preset, f"case:{case_id}")]
                    frame = {
                        "program": program["id"],
                        "preset": preset,
                        "case": case_id,
                        "resolution": resolution,
                        "probe": "steady",
                        "operation": _frame_operation(record),
                        "pixels": record["pixels"],
                    }
                    frame["sha256"] = hashlib.sha256(record["pixels"]).hexdigest()
                    frames.append(frame)
                for category, probes in program["probes"].items():
                    for probe in probes:
                        name = f"{category}:{probe}"
                        record = backend_records[(dimensions, preset, name)]
                        frame = {
                            "program": program["id"],
                            "preset": preset,
                            "case": "default",
                            "resolution": resolution,
                            "probe": name,
                            "operation": _frame_operation(record),
                            "pixels": record["pixels"],
                        }
                        frame["sha256"] = hashlib.sha256(
                            record["pixels"]).hexdigest()
                        frames.append(frame)
    oracle_metrics = []
    for oracle in oracles:
        manifest_metric = next(
            metric
            for metric in oracle["metrics"]
            if metric["domain"] == "FRAMEBUFFER"
            and metric["aggregation"] == "MAXIMUM"
        )
        observed = oracle_totals[oracle["oracle_id"]]
        accepted_limit = manifest_metric["accepted_limit"]
        if observed["value"] > accepted_limit:
            raise CaptureError(
                f'{oracle["oracle_id"]} framebuffer error exceeds its accepted '
                f'limit: {observed["value"]} > {accepted_limit}'
            )
        print(
            f'pullback metric: {configuration} aggregate {oracle["oracle_id"]} '
            f'observed={observed["value"]} limit={accepted_limit}'
        )
        oracle_metrics.append(
            {
                "oracle_id": oracle["oracle_id"],
                "domain": "FRAMEBUFFER",
                "aggregation": "MAXIMUM",
                "value": observed["value"],
                "unit": manifest_metric["unit"],
                "sample_count": observed["sample_count"],
                "resolution_values": observed["resolution_values"],
            }
        )
    capture = {
        "schema_version": 1,
        "checkout_sha": checkout_sha,
        "manifest_sha256": manifest_sha256(programs, oracles, schema),
        "toolchain": toolchain,
        "configuration": configuration,
        "corpus": programs["corpus"],
        "oracle_metrics": oracle_metrics,
    }
    write_capture(capture, frames, output)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--configuration",
        choices=("native-debug", "wasm-release", "wasm-strict-fp"),
    )
    parser.add_argument("--checkout", type=Path)
    parser.add_argument("--build-dir", type=Path)
    parser.add_argument("--manifest-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--operations-output", type=Path)
    parser.add_argument("--backend-audit", type=Path, nargs="+")
    args = parser.parse_args()
    try:
        if args.operations_output is not None:
            programs, oracles, _ = load_and_validate(args.manifest_dir.resolve())
            write_operations(programs, oracles, args.operations_output.resolve())
            return 0
        if args.backend_audit is not None:
            if args.configuration is None:
                parser.error("backend audit requires --configuration")
            programs, oracles, _ = load_and_validate(args.manifest_dir.resolve())
            metrics = {
                oracle["oracle_id"]: {"value": 0, "sample_count": 0}
                for oracle in oracles
            }
            resolutions = set()
            for path in args.backend_audit:
                resolution, _, observed = load_backend(
                    path.resolve(), programs, oracles
                )
                resolutions.add(resolution)
                validate_backend_oracles(
                    observed, oracles, args.configuration, resolution
                )
                for oracle, metric in observed.items():
                    metrics[oracle]["value"] = max(
                        metrics[oracle]["value"], metric["value"]
                    )
                    metrics[oracle]["sample_count"] += metric["sample_count"]
            expected_resolutions = {
                tuple(resolution) for resolution in programs["corpus"]["resolutions"]
            }
            if resolutions != expected_resolutions:
                raise CaptureError("backend audit does not cover the resolution roster")
            validate_backend_oracles(metrics, oracles, args.configuration)
            return 0
        if None in (args.configuration, args.checkout, args.build_dir, args.output):
            parser.error(
                "capture requires --configuration, --checkout, --build-dir, and --output"
            )
        produce(
            args.configuration,
            args.checkout.resolve(),
            args.build_dir.resolve(),
            args.manifest_dir.resolve(),
            args.output.resolve(),
        )
    except (
        CaptureError,
        ManifestError,
        OSError,
        subprocess.SubprocessError,
    ) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
