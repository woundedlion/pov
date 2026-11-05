#!/usr/bin/env python3
"""Generate pullback operations and audit binary backend metrics."""

from __future__ import annotations

import argparse
import struct
import sys
from pathlib import Path

from generate_pullback_manifest_header import (
    ManifestError,
    protocol_definition,
    load_and_validate,
)


class CaptureError(RuntimeError):
    """A capture backend stream or corpus is invalid."""


CASE_OPERATIONS = {
    "default": "CASE_DEFAULT",
    "endpoint_min": "CASE_ENDPOINT_MIN",
    "endpoint_max": "CASE_ENDPOINT_MAX",
    "interior": "CASE_INTERIOR",
}


def corpus_operations(programs: dict):
    """Yield the canonical preset, parameter-case and probe operations."""
    mappings = programs["corpus"]["probe_operations"]
    for program in programs["programs"]:
        for preset in program["presets"]:
            for case in program["parameter_cases"]:
                case_id = case["id"]
                yield {"program": program["id"], "preset": preset, "case": case_id,
                       "probe": "steady", "name": f"case:{case_id}",
                       "mapping": CASE_OPERATIONS[case_id]}
            for category, probes in program["probes"].items():
                for probe in probes:
                    name = f"{category}:{probe}"
                    yield {"program": program["id"], "preset": preset, "case": "default",
                           "probe": name, "name": name, "mapping": mappings[name]}


def framebuffer_metric(oracle: dict) -> dict:
    """Return the maximum framebuffer metric from an oracle manifest."""
    return next(metric for metric in oracle["metrics"]
                if metric["domain"] == "FRAMEBUFFER" and metric["aggregation"] == "MAXIMUM")


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
    raise CaptureError(f"unknown spatial probe mapping {mapping}")


def _read_u16(data: bytes, offset: int) -> tuple[int, int]:
    if offset + 2 > len(data):
        raise CaptureError("capture backend stream is truncated")
    return struct.unpack_from("<H", data, offset)[0], offset + 2


def _read_u32(data: bytes, offset: int) -> tuple[int, int]:
    low, offset = _read_u16(data, offset)
    high, offset = _read_u16(data, offset)
    return low | high << 16, offset


def operation_specs(programs: dict) -> list[dict]:
    specs = [{name: op[name] for name in ("preset", "name", "mapping")}
             for op in corpus_operations(programs)]
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
        if mapping not in protocol_definition()[1]:
            raise CaptureError(f"unknown capture operation mapping {mapping}")
        data.extend(
            struct.pack("<HHH", spec["preset"], protocol_definition()[1][mapping], len(name))
        )
        data.extend(name)
    for spec in oracle_specs:
        name = spec["oracle"].encode("utf-8")
        mapping = spec["mapping"]
        if mapping not in protocol_definition()[1]:
            raise CaptureError(f"unknown oracle operation mapping {mapping}")
        data.extend(struct.pack("<H", len(name)))
        data.extend(name)
        data.extend(
            struct.pack(
                "<HHf",
                spec["preset"],
                protocol_definition()[1][mapping],
                spec["hue_noise_phase"],
            )
        )
    path.write_bytes(data)


def load_backend(
    path: Path, programs: dict, oracles: list[dict]
) -> tuple[tuple[int, int], dict]:
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
    records = set()
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
        for _ in range(4):
            _, offset = _read_u16(data, offset)
        selected_pixels, offset = _read_u32(data, offset)
        key = (preset, name)
        if key in records or key not in specs:
            raise CaptureError(f"invalid capture backend record {key}")
        expected_code = protocol_definition()[1][specs[key]["mapping"]]
        if operation_code != expected_code:
            raise CaptureError(f"capture backend operation mismatch for {key}")
        mapping = specs[key]["mapping"]
        if not (mapping.startswith("CASE_") or mapping.startswith("THROUGH_CLEAR_")
                or mapping == "FULL_FRAME"):
            expected_pixels = sum(_selected_pixel(mapping, x, y, width, height)
                                  for y in range(height) for x in range(width))
            if selected_pixels != expected_pixels:
                raise CaptureError(f"capture backend selected-pixel mismatch for {key}")
        end = offset + width * height * 6
        if end > len(data):
            raise CaptureError("capture backend stream is truncated")
        offset = end
        records.add(key)
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
        or metrics.keys() != expected_oracles
    ):
        raise CaptureError("capture backend stream is incomplete or has trailing bytes")
    return (width, height), metrics


def validate_backend_oracles(
    metrics: dict,
    oracles: list[dict],
    configuration: str,
    resolution: tuple[int, int] | None = None,
) -> None:
    for oracle in oracles:
        oracle_id = oracle["oracle_id"]
        manifest_metric = framebuffer_metric(oracle)
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--configuration",
        choices=("native-debug", "wasm-release", "wasm-strict-fp"),
    )
    parser.add_argument("--manifest-dir", type=Path, required=True)
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
                resolution, observed = load_backend(
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
        parser.error("requires --operations-output or --backend-audit")
    except (
        CaptureError,
        ManifestError,
        OSError,
    ) as error:
        print(f"pullback_capture: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
