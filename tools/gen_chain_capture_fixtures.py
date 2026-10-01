# Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
# Licensed under the PolyForm Noncommercial License 1.0.0
"""Emit typed C++ snapshots for the authored pullback capture cases."""

import argparse
import json
import math
from pathlib import Path
import struct


def float_cpp(value):
    if not math.isfinite(value):
        raise ValueError("capture fixture requires finite floats")
    bits = struct.unpack("<I", struct.pack("<f", value))[0]
    return f"f(0x{bits:08x}u)"


def string_cpp(value):
    return json.dumps(value, ensure_ascii=True)


def walk_cpp(state):
    def vector(values):
        return "math::Vector(" + ",".join(map(float_cpp, values)) + ")"

    def quaternion(values):
        return "math::Quaternion(" + ",".join(map(float_cpp, values)) + ")"

    return (
        "Pullback::Interp::SpatialWalkSnapshot{"
        + ",".join(
            [
                str(state["noiseSeed"]),
                str(state["walkTime"]),
                vector(state["position"]),
                vector(state["direction"]),
                quaternion(state["wander"]),
                float_cpp(state["angularVelocity"]),
                float_cpp(state["spinPhase"]),
            ]
        )
        + "}"
    )


def runtime_cpp(entry):
    state = entry["state"]
    kind = entry["kind"]
    if kind == "spatial-walk-v2":
        return walk_cpp(state)
    names = {
        "source-clock-v1": ("SourceClockSnapshot", ["primary", "secondary", "angle"]),
        "noise-clock-v1": ("NoiseClockSnapshot", ["phase", "noiseSeed"]),
        "phase-clock-v1": ("PhaseClockSnapshot", ["phase"]),
        "ripple-clock-v1": ("RippleClockSnapshot", ["phase"]),
        "affine-clock-v1": ("AffineClockSnapshot", ["phase", "rotation"]),
        "color-clock-v1": (
            "ColorClockSnapshot",
            ["oscillationPhase", "hueNoisePhase", "hueNoiseSeed"],
        ),
    }
    if kind == "spherical-rings-v2":
        return (
            "Pullback::Interp::SphericalRingsSnapshot{"
            + walk_cpp(state["walk"])
            + ","
            + float_cpp(state["phase"])
            + "}"
        )
    typename, fields = names[kind]
    values = [
        str(state[field]) if field.endswith("Seed") else float_cpp(state[field])
        for field in fields
    ]
    return "Pullback::Interp::" + typename + "{" + ",".join(values) + "}"


def snapshot_cpp(snapshot):
    chain = ",".join(
        "{" + string_cpp(entry["instance"]) + "," + string_cpp(entry["operator"]) + "}"
        for entry in snapshot["chain"]
    )
    parameters = ",".join(
        "{" + string_cpp(entry["name"]) + "," + float_cpp(entry["value"]) + "}"
        for entry in snapshot["parameters"]
    )
    runtime = ",".join(
        "{" + string_cpp(entry["instance"]) + "," + runtime_cpp(entry) + "}"
        for entry in snapshot["runtime"]
    )
    bank = snapshot["paletteBank"]
    cycles = ",".join(
        "{"
        + ",".join(
            [
                str(clock["frame"]),
                str(clock["nextSequence"]),
                str(clock["fadeActive"]).lower(),
                str(clock["displayDirty"]).lower(),
            ]
        )
        + "}"
        for clock in bank["cycles"]
    )
    palette = (
        "GeneratedPaletteBank::Snapshot{"
        + float_cpp(bank["chroma"])
        + ",{ "
        + ",".join(map(str, bank["hues"]))
        + "},{{"
        + cycles
        + "}}}"
    )
    return (
        "ChainSnapshot{ChainSnapshot::SCHEMA_VERSION,{"
        + chain
        + "},{"
        + parameters
        + "},std::vector<ChainSnapshot::Runtime>{"
        + runtime
        + "},"
        + palette
        + ","
        + str(snapshot["animationsPaused"]).lower()
        + "}"
    )


def generate(records):
    snapshots = []
    snapshot_indices = {}
    cases = []
    keys = set()
    for record in records:
        if record["kind"] not in ("frame", "oracle"):
            raise ValueError("unknown capture fixture kind")
        for field in (
            "preset",
            "operation",
            "source",
            "destination",
            "elapsed",
            "duration",
        ):
            value = record[field]
            if type(value) is not int or not 0 <= value <= 65535:
                raise ValueError(f"capture fixture {field} exceeds uint16")
        for field in ("width", "height"):
            value = record[field]
            if type(value) is not int or not 0 < value <= 2147483647:
                raise ValueError(f"capture fixture {field} is not a positive int32")
        snapshot = record["snapshot"]
        if snapshot["schemaVersion"] != 2 or not snapshot["chain"]:
            raise ValueError("capture fixture requires a version-two chain")
        canonical = json.dumps(snapshot, sort_keys=True, separators=(",", ":"))
        if canonical not in snapshot_indices:
            snapshot_indices[canonical] = len(snapshots)
            snapshots.append(snapshot)
        phase = record.get("hueNoisePhase", 0.0)
        key = (
            record["kind"],
            record["name"],
            record["width"],
            record["height"],
            record["preset"],
            record["operation"],
            struct.pack("<f", phase),
        )
        if key in keys:
            raise ValueError(f"duplicate capture fixture {key}")
        keys.add(key)
        cases.append(
            "{"
            + ",".join(
                [
                    string_cpp(record["kind"]),
                    string_cpp(record["name"]),
                    str(record["width"]),
                    str(record["height"]),
                    str(record["preset"]),
                    str(record["operation"]),
                    float_cpp(phase),
                    str(record["source"]),
                    str(record["destination"]),
                    str(record["elapsed"]),
                    str(record["duration"]),
                    str(snapshot_indices[canonical]),
                ]
            )
            + "}"
        )
    lines = [
        "/* Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.",
        " * Licensed under the PolyForm Noncommercial License 1.0.0 */",
        "#pragma once",
        "#include <bit>",
        "#include <string_view>",
        '#include "workbench/shader/chain_snapshot.h"',
        "namespace ChainCaptureFixtures {",
        "constexpr float f(uint32_t bits) { return std::bit_cast<float>(bits); }",
        "struct Case { const char *kind,*name; int width,height; uint16_t preset,operation; float phase; uint16_t source,destination,elapsed,duration; size_t snapshot_index; };",
        "inline constexpr Case CASES[] = {",
        ",\n".join(cases),
        "};",
        "inline ChainSnapshot snapshot(size_t index) { switch(index) {",
    ]
    lines.extend(
        f"case {index}: return {snapshot_cpp(snapshot)};"
        for index, snapshot in enumerate(snapshots)
    )
    lines.extend(
        [
            'default: HS_CHECK(false,"invalid capture fixture index"); return {}; } }',
            "}",
            "",
        ]
    )
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    records = [
        json.loads(line)
        for line in args.input.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(generate(records), encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
