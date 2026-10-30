# Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
# Licensed under the PolyForm Noncommercial License 1.0.0
"""Emit typed C++ snapshots for the frozen, hand-edited pullback capture corpus.

The unit_chain_capture_fixtures test validates the authored corpus.
"""

import argparse
import copy
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


def vector_cpp(values):
    return "math::Vector(" + ",".join(map(float_cpp, values)) + ")"


def sine_cpp(params):
    return "Pattern::SineParams{" + ",".join([
        vector_cpp(params["axis"]), float_cpp(params["frequency"]),
        float_cpp(params["offset"]), float_cpp(params["speed"]),
    ]) + "}"


def interference_cpp(params):
    return "Pattern::InterferenceParams{" + ",".join([
        sine_cpp(params), vector_cpp(params["axisB"]),
        float_cpp(params["frequencyB"]), float_cpp(params["offsetB"]),
    ]) + "}"


def enum_cpp(typename, value, count):
    if type(value) is not int or not 0 <= value < count:
        raise ValueError(f"capture fixture has an invalid {typename}")
    return f"static_cast<{typename}>({value})"


def recipe_cpp(recipe):
    if recipe["version"] != 1 or type(recipe["enabled"]) is not bool:
        raise ValueError("capture fixture requires a version-one pattern recipe")
    noise, mechanic = recipe["noise"], recipe["mechanic"]
    return "ColorModulation::Config{" + ",".join([
        str(recipe["enabled"]).lower(),
        enum_cpp("Pattern::SourceKind", recipe["source"], 5), "1",
        enum_cpp("Pattern::Evaluation", recipe["evaluation"], 2),
        "Pattern::NoiseParams{" + ",".join([
            float_cpp(noise["scale"]), float_cpp(noise["speed"]),
            enum_cpp("Pattern::NoiseBasis", noise["basis"], 3),
            float_cpp(noise["contrast"]),
            enum_cpp("Pattern::NoiseVariant", noise["variant"], 2),
        ]) + "}", sine_cpp(recipe["sine"]),
        interference_cpp(recipe["interference"]),
        "Pattern::GridParams{" + interference_cpp(recipe["grid"]) + "}",
        str(recipe["seed"]),
        "ColorModulation::Params{" + enum_cpp(
            "ColorModulation::Mechanic", mechanic["kind"], 3
        ) + "," + float_cpp(mechanic["amount"]) + "}",
    ]) + "}"


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
    if kind == "color-clock-v3":
        return "Pullback::Interp::ColorClockSnapshot{" + float_cpp(state["oscillationPhase"]) + "}"
    if kind == "modulation-clock-v1":
        if state["version"] != 1:
            raise ValueError("capture fixture requires a version-one pattern clock")
        return "Pullback::Interp::ModulationClockSnapshot{" + ",".join([
            "Pullback::Modulation::State{" + ",".join([
                enum_cpp("Pattern::SourceKind", state["source"], 5), "1",
                "Pattern::PhaseState{" + float_cpp(state["phaseTurns"]) + "}",
                str(state["initialized"]).lower(),
            ]) + "}",
        ]) + "}"
    names = {
        "source-clock-v1": ("SourceClockSnapshot", ["primary", "secondary", "angle"]),
        "noise-clock-v1": ("NoiseClockSnapshot", ["phase", "noiseSeed"]),
        "phase-clock-v1": ("PhaseClockSnapshot", ["phase"]),
        "ripple-clock-v1": ("RippleClockSnapshot", ["phase"]),
        "affine-clock-v1": ("AffineClockSnapshot", ["phase", "rotation"]),
    }
    if kind == "spherical-rings-v2":
        return (
            "Pullback::Interp::SphericalRingsSnapshot{"
            + walk_cpp(state["walk"])
            + ","
            + float_cpp(state["phase"])
            + "}"
        )
    if kind not in names:
        raise ValueError("unknown capture runtime kind")
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
    modulation = ",".join(
        "{" + string_cpp(entry["instance"]) + "," + recipe_cpp(entry["recipe"]) + "}"
        for entry in snapshot["modulation"]
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
        + "},{"
        + modulation
        + "},std::vector<ChainSnapshot::Runtime>{"
        + runtime
        + "},"
        + palette
        + ","
        + str(snapshot["animationsPaused"]).lower()
        + "}"
    )


def migrate_legacy_snapshot(snapshot):
    if snapshot["schemaVersion"] != 3:
        return snapshot
    result = copy.deepcopy(snapshot)
    result["schemaVersion"] = 4
    used = {entry["instance"] for entry in result["chain"]}

    def label_for(base):
        label = base
        suffix = 2
        while label in used:
            label = f"{base}-{suffix}"
            suffix += 1
        used.add(label)
        return label

    chain = []
    recipes = []
    for entry in result["chain"]:
        chain.append(entry)
        if entry["operator"] != "colorize.generated-palette.v3":
            continue
        owner = entry["instance"]
        displacement, spatial = label_for("displacement"), label_for("modulation")
        entry["operator"] = "colorize.generated-palette.v4"
        chain.extend([
            {"instance": displacement, "operator": "color.modulate.v1"},
            {"instance": spatial, "operator": "color.modulate.v1"},
        ])
        values = {parameter["name"]: parameter["value"] for parameter in result["parameters"]}
        result["parameters"] = [parameter for parameter in result["parameters"]
                                if parameter["name"] not in
                                (owner + ".hue-shift-mode", owner + ".hue-shift-amount")]
        original = next(recipe["recipe"] for recipe in result["modulation"]
                        if recipe["instance"] == owner)
        wave = {"axis": [0, 1, 0], "frequency": 1, "offset": 0, "speed": 0}
        interference = {**wave, "axisB": [1, 0, 0], "frequencyB": 1, "offsetB": 0}
        path = {
            "enabled": values.get(owner + ".hue-shift-mode", 0) == 1,
            "source": 4, "version": 1, "evaluation": 0, "seed": 6047,
            "noise": {"scale": 1, "speed": 0, "contrast": 0, "basis": 0, "variant": 0},
            "sine": wave, "interference": interference, "grid": interference,
            "mechanic": {"kind": 0, "amount": values.get(owner + ".hue-shift-amount", 0)},
        }
        recipes.extend([{"instance": displacement, "recipe": path},
                        {"instance": spatial, "recipe": original}])
        clock = next(value for value in result["runtime"] if value["instance"] == owner)
        if clock["kind"] != "color-clock-v2":
            raise ValueError("unknown capture runtime kind")
        old_clock = clock["state"].pop("modulation")
        clock["kind"] = "color-clock-v3"
        result["runtime"].extend([
            {"instance": displacement, "kind": "modulation-clock-v1", "state": {
                "source": 4, "version": 1, "phaseTurns": 0, "initialized": True}},
            {"instance": spatial, "kind": "modulation-clock-v1", "state": old_clock},
        ])
    result["chain"] = chain
    result["modulation"] = recipes
    return result


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
        snapshot = migrate_legacy_snapshot(record["snapshot"])
        if snapshot["schemaVersion"] != 4 or not snapshot["chain"]:
            raise ValueError("capture fixture requires a version-four chain")
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
        "#include <cstddef>",
        "#include <cstdint>",
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
