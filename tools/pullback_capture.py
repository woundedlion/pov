#!/usr/bin/env python3
"""Produce canonical pullback captures from native or WASM ShaderBall builds."""

from __future__ import annotations

import argparse
import hashlib
import json
import struct
import subprocess
import tempfile
from pathlib import Path

from generate_pullback_manifest_header import load_and_validate, manifest_sha256
from pullback_crosscheck import _canonical_frame_bytes


CASE_IDS = ("default", "endpoint_min", "endpoint_max", "interior")


class CaptureError(RuntimeError):
    """A capture backend or its binary stream is invalid."""


def _read_u16(data: bytes, offset: int) -> tuple[int, int]:
    if offset + 2 > len(data):
        raise CaptureError("capture backend stream is truncated")
    return struct.unpack_from("<H", data, offset)[0], offset + 2


def load_backend(path: Path) -> tuple[tuple[int, int], dict[tuple, list]]:
    data = path.read_bytes()
    if data[:4] != b"HSPB":
        raise CaptureError("capture backend magic mismatch")
    offset = 4
    header = []
    for _ in range(6):
        value, offset = _read_u16(data, offset)
        header.append(value)
    version, width, height, preset_count, case_count, records_low = header
    records_high, offset = _read_u16(data, offset)
    record_count = records_low | records_high << 16
    if (version, preset_count, case_count, record_count) != (1, 12, 4, 48):
        raise CaptureError("capture backend header mismatch")
    pixels = width * height
    records = {}
    for _ in range(record_count):
        preset, offset = _read_u16(data, offset)
        case_index, offset = _read_u16(data, offset)
        if preset >= preset_count or case_index >= case_count:
            raise CaptureError("capture backend record key is invalid")
        key = (preset, CASE_IDS[case_index])
        if key in records:
            raise CaptureError(f"duplicate capture backend record {key}")
        frame = []
        for _pixel in range(pixels):
            channels = []
            for _channel in range(3):
                channel, offset = _read_u16(data, offset)
                channels.append(channel)
            frame.append([*channels, 65535])
        records[key] = frame
    if offset != len(data):
        raise CaptureError("capture backend stream has trailing bytes")
    return (width, height), records


def _backend_command(
    configuration: str,
    producer_root: Path,
    build_dir: Path,
    resolution: tuple[int, int],
    output: Path,
) -> list[str]:
    dimensions = f"{resolution[0]}x{resolution[1]}"
    if configuration == "native-debug":
        executable = build_dir / "tests/pullback_capture_native"
        if not executable.exists():
            executable = executable.with_suffix(".exe")
        return [str(executable), "--resolution", dimensions, str(output)]
    return [
        "node",
        str(producer_root / "scripts/pullback_capture_wasm.mjs"),
        str(build_dir / "holosphere_wasm.js"),
        dimensions,
        str(output),
    ]


def produce(
    configuration: str,
    checkout: Path,
    build_dir: Path,
    producer_root: Path,
    manifest_dir: Path,
    output: Path,
) -> None:
    programs, oracles = load_and_validate(manifest_dir)
    checkout_sha = subprocess.run(
        ["git", "-C", str(checkout), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    raw = {}
    with tempfile.TemporaryDirectory(prefix="pullback-capture-") as directory:
        for resolution in programs["corpus"]["resolutions"]:
            dimensions = tuple(resolution)
            backend_output = Path(directory) / f"{dimensions[0]}x{dimensions[1]}.bin"
            subprocess.run(
                _backend_command(
                    configuration, producer_root, build_dir, dimensions, backend_output
                ),
                cwd=checkout,
                check=True,
            )
            actual_resolution, records = load_backend(backend_output)
            if actual_resolution != dimensions:
                raise CaptureError("capture backend resolution mismatch")
            for key, pixels in records.items():
                raw[(dimensions, *key)] = pixels
    frames = []
    for program in programs["programs"]:
        for preset in program["presets"]:
            for resolution in programs["corpus"]["resolutions"]:
                dimensions = tuple(resolution)
                for case in program["parameter_cases"]:
                    frame = {
                        "program": program["id"],
                        "preset": preset,
                        "case": case["id"],
                        "resolution": resolution,
                        "probe": "steady",
                        "pixels": raw[(dimensions, preset, case["id"])],
                    }
                    frame["sha256"] = hashlib.sha256(
                        _canonical_frame_bytes(frame)
                    ).hexdigest()
                    frames.append(frame)
                for category, probes in program["probes"].items():
                    for probe in probes:
                        frame = {
                            "program": program["id"],
                            "preset": preset,
                            "case": "default",
                            "resolution": resolution,
                            "probe": f"{category}:{probe}",
                            "pixels": raw[(dimensions, preset, "default")],
                        }
                        frame["sha256"] = hashlib.sha256(
                            _canonical_frame_bytes(frame)
                        ).hexdigest()
                        frames.append(frame)
    capture = {
        "schema_version": 1,
        "checkout_sha": checkout_sha,
        "manifest_sha256": manifest_sha256(programs, oracles),
        "toolchains": programs["toolchains"],
        "configuration": configuration,
        "corpus": programs["corpus"],
        "frames": frames,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(capture, separators=(",", ":")) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--configuration",
        required=True,
        choices=("native-debug", "wasm-release", "wasm-strict-fp"),
    )
    parser.add_argument("--checkout", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument(
        "--producer-root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--manifest-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        produce(
            args.configuration,
            args.checkout.resolve(),
            args.build_dir.resolve(),
            args.producer_root.resolve(),
            args.manifest_dir.resolve(),
            args.output.resolve(),
        )
    except (CaptureError, OSError, subprocess.CalledProcessError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
