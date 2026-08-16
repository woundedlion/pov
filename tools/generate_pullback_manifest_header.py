#!/usr/bin/env python3
"""Validate pullback manifests and generate their native-test header."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import re
from pathlib import Path


TOPOLOGY_FIELDS = (
    "function", "projection", "projection_frame", "surface_lens",
    "signal_weight", "value_transfer", "coverage", "peirce_layout",
    "airocean_layout", "bonne_hemisphere", "gnomonic_hemisphere",
    "surface_noise", "surface_noise_placement", "surface_noise_basis",
    "surface_curl_integrator", "source_noise_basis", "outer_warp",
    "outer_warp_basis", "outer_warp_envelope", "outer_polar_mode",
    "outer_curl_integrator", "outer_polar_harmonic", "inner_warp",
    "inner_warp_basis", "inner_warp_envelope", "inner_polar_mode",
    "inner_curl_integrator", "inner_polar_harmonic",
)
ORACLE_FILES = ("peirce_fast_square.json", "hue_rotation_noise_luts.json")
SHA_RE = re.compile(r"^[0-9a-f]{40}$")
CASE_IDS = {"default", "endpoint_min", "endpoint_max", "interior"}
PRESET_COUNT = 19


class ManifestError(ValueError):
    """A manifest violates the checked schema or cross-file contract."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ManifestError(message)


def _reject_json_constant(token: str):
    raise ValueError(f"non-JSON numeric constant {token}")


def _load(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"),
                           parse_constant=_reject_json_constant)
    except (OSError, ValueError) as error:
        raise ManifestError(f"{path}: {error}") from error
    _require(isinstance(value, dict), f"{path}: root must be an object")
    return value


def _schema_type_matches(value, expected: str) -> bool:
    if expected == "object":
        return isinstance(value, dict)
    if expected == "array":
        return isinstance(value, list)
    if expected == "string":
        return isinstance(value, str)
    if expected == "integer":
        return isinstance(value, int) and not isinstance(value, bool)
    if expected == "number":
        return (isinstance(value, (int, float)) and
                not isinstance(value, bool) and math.isfinite(value))
    if expected == "boolean":
        return isinstance(value, bool)
    if expected == "null":
        return value is None
    raise ManifestError(f"schema uses unsupported type {expected}")


def _json_equal(first, second) -> bool:
    return type(first) is type(second) and first == second


def _is_number(value) -> bool:
    return (isinstance(value, (int, float)) and not isinstance(value, bool) and
            math.isfinite(value))


def _resolve_ref(root: dict, reference: str) -> dict:
    _require(reference.startswith("#/$defs/"),
             f"unsupported schema reference {reference}")
    name = reference.removeprefix("#/$defs/")
    resolved = root.get("$defs", {}).get(name)
    _require(isinstance(resolved, dict), f"missing schema definition {name}")
    return resolved


def _validate_schema(value, schema: dict, root: dict, location: str) -> None:
    if "$ref" in schema:
        _validate_schema(value, _resolve_ref(root, schema["$ref"]), root,
                         location)
        return
    if "oneOf" in schema:
        matches = 0
        errors = []
        for alternative in schema["oneOf"]:
            try:
                _validate_schema(value, alternative, root, location)
                matches += 1
            except ManifestError as error:
                errors.append(str(error))
        _require(matches == 1,
                 f"{location}: expected exactly one schema alternative; "
                 f"matched {matches}: {'; '.join(errors)}")
        return
    if "const" in schema:
        _require(_json_equal(value, schema["const"]),
                 f"{location}: expected {schema['const']!r}")
    if "enum" in schema:
        _require(any(_json_equal(value, allowed) for allowed in schema["enum"]),
                 f"{location}: value is outside the allowed enum")
    expected_types = schema.get("type")
    if expected_types is not None:
        if isinstance(expected_types, str):
            expected_types = [expected_types]
        _require(any(_schema_type_matches(value, expected)
                     for expected in expected_types),
                 f"{location}: expected type {' or '.join(expected_types)}")
    if isinstance(value, str):
        _require(len(value) >= schema.get("minLength", 0),
                 f"{location}: string is too short")
        if "pattern" in schema:
            _require(re.fullmatch(schema["pattern"], value) is not None,
                     f"{location}: string does not match {schema['pattern']}")
    if _is_number(value):
        _require(value >= schema.get("minimum", value),
                 f"{location}: value is below minimum")
        _require(value <= schema.get("maximum", value),
                 f"{location}: value is above maximum")
    if isinstance(value, list):
        _require(len(value) >= schema.get("minItems", 0),
                 f"{location}: array is too short")
        _require(len(value) <= schema.get("maxItems", len(value)),
                 f"{location}: array is too long")
        if schema.get("uniqueItems"):
            canonical = [json.dumps(item, sort_keys=True) for item in value]
            _require(len(canonical) == len(set(canonical)),
                     f"{location}: array entries must be unique")
        if "items" in schema:
            for index, item in enumerate(value):
                _validate_schema(item, schema["items"], root,
                                 f"{location}[{index}]")
    if isinstance(value, dict):
        required = schema.get("required", [])
        missing = sorted(set(required) - value.keys())
        _require(not missing, f"{location}: missing fields {missing}")
        properties = schema.get("properties", {})
        for name, item in value.items():
            if name in properties:
                _validate_schema(item, properties[name], root,
                                 f"{location}.{name}")
                continue
            additional = schema.get("additionalProperties", True)
            _require(additional is not False,
                     f"{location}: unexpected field {name}")
            if isinstance(additional, dict):
                _validate_schema(item, additional, root,
                                 f"{location}.{name}")
        _require(len(value) >= schema.get("minProperties", 0),
                 f"{location}: object has too few fields")


def _validate_common(document: dict, path: Path) -> None:
    _require(type(document.get("schema_version")) is int and
             document["schema_version"] == 1,
             f"{path}: schema_version must be 1")
    _require(SHA_RE.fullmatch(document.get("base_sha", "")) is not None,
             f"{path}: base_sha must be a full lowercase Git SHA")
    toolchains = document.get("toolchains")
    _require(isinstance(toolchains, dict), f"{path}: toolchains must be an object")
    for name in ("native", "wasm", "teensy"):
        _require(isinstance(toolchains.get(name), dict),
                 f"{path}: toolchains.{name} must be an object")
        _require(bool(toolchains[name].get("compiler")),
                 f"{path}: toolchains.{name}.compiler is required")


def _validate_programs(document: dict, path: Path) -> None:
    _validate_common(document, path)
    _require(document.get("kind") == "pullback-programs",
             f"{path}: kind must be pullback-programs")
    corpus = document.get("corpus")
    _require(isinstance(corpus, dict), f"{path}: corpus must be an object")
    _require(corpus.get("generator") == "shaderball-pullback-corpus",
             f"{path}: unknown corpus generator")
    _require(type(corpus.get("version")) is int and corpus["version"] == 1,
             f"{path}: corpus.version must be 1")
    _require(corpus.get("resolutions") == [[96, 20], [288, 144]],
             f"{path}: both roster resolutions are required")
    probe_operations = corpus.get("probe_operations")
    _require(isinstance(probe_operations, dict) and probe_operations,
             f"{path}: corpus.probe_operations must be a non-empty object")
    programs = document.get("programs")
    _require(isinstance(programs, list), f"{path}: programs must be an array")
    ids = [program.get("id") for program in programs
           if isinstance(program, dict)]
    _require(len(ids) == len(programs) and len(ids) == len(set(ids)),
             f"{path}: program IDs must be unique")
    covered_presets: list[int] = []
    referenced_probes: set[str] = set()
    for program in programs:
        program_id = program["id"]
        presets = program.get("presets")
        _require(isinstance(presets, list) and presets,
                 f"{path}: {program_id}.presets must be non-empty")
        _require(all(type(index) is int and 0 <= index < PRESET_COUNT
                     for index in presets),
                 f"{path}: {program_id}.presets contains an invalid index")
        covered_presets.extend(presets)
        topology = program.get("topology_key")
        _require(set(topology) == set(TOPOLOGY_FIELDS),
                 f"{path}: {program_id}.topology_key fields are incomplete")
        cases = program.get("parameter_cases")
        _require(isinstance(cases, list),
                 f"{path}: {program_id}.parameter_cases must be an array")
        _require({case.get("id") for case in cases if isinstance(case, dict)} ==
                 CASE_IDS,
                 f"{path}: {program_id} must define default, endpoints, and interior cases")
        probes = program.get("probes")
        _require(isinstance(probes, dict),
                 f"{path}: {program_id}.probes must be an object")
        for category in ("seam", "singularity", "boundary", "transition"):
            _require(isinstance(probes.get(category), list) and probes[category],
                     f"{path}: {program_id}.probes.{category} must be non-empty")
            referenced_probes.update(
                f"{category}:{probe}" for probe in probes[category])
        tolerances = program.get("tolerances")
        _require(isinstance(tolerances, dict),
                 f"{path}: {program_id}.tolerances must be an object")
        _require(tolerances.get("native_bit_identity") is True,
                 f"{path}: {program_id} native output must be bit-identical")
        _require(isinstance(tolerances.get("carriers"), dict),
                 f"{path}: {program_id}.tolerances.carriers must be an object")
        framebuffer = tolerances.get("release_wasm_framebuffer")
        _require(isinstance(framebuffer, dict),
                 f"{path}: {program_id} release WASM tolerance is required")
        delta = framebuffer.get("max_channel_delta_u16")
        fraction = framebuffer.get("max_differing_pixel_fraction")
        _require(type(delta) is int and 0 <= delta <= 4,
                 f"{path}: {program_id} exceeds the channel-delta license")
        _require(_is_number(fraction) and 0 <= fraction <= 0.001,
                 f"{path}: {program_id} exceeds the differing-pixel license")
    _require(sorted(covered_presets) == list(range(PRESET_COUNT)),
             f"{path}: presets must be covered exactly once")
    _require(set(probe_operations) == referenced_probes,
             f"{path}: probe operation mappings must exactly cover referenced probes")


def _validate_oracle(document: dict, path: Path) -> None:
    _validate_common(document, path)
    _require(document.get("kind") == "pullback-oracle",
             f"{path}: kind must be pullback-oracle")
    _require(bool(document.get("oracle_id")), f"{path}: oracle_id is required")
    _require(bool(document.get("exact_callable")),
             f"{path}: exact_callable is required")
    _require(document.get("non_floating_fields_exact") is True,
             f"{path}: non-floating fields must be exact")
    corpus = document.get("corpus")
    _require(isinstance(corpus, dict) and
             type(corpus.get("version")) is int and corpus["version"] == 1,
             f"{path}: a versioned deterministic corpus is required")
    _require(corpus.get("deterministic") is True,
             f"{path}: oracle corpus must be deterministic")
    framebuffer = corpus.get("framebuffer")
    _require(isinstance(framebuffer, dict),
             f"{path}: framebuffer corpus is required")
    probes = framebuffer.get("probes")
    _require(framebuffer.get("sampling") == "pixel centers" and
             framebuffer.get("aggregation") ==
             "maximum absolute premultiplied linear-u16 channel delta" and
             bool(framebuffer.get("exact_callable")) and
             bool(framebuffer.get("optimized_callable")) and
             isinstance(probes, list) and probes,
             f"{path}: framebuffer pixel-center probes are required")
    corpus_presets = corpus.get("presets")
    _require(isinstance(corpus_presets, list) and corpus_presets,
             f"{path}: oracle presets are required")
    probe_presets = []
    for probe in probes:
        _require(isinstance(probe, dict),
                 f"{path}: framebuffer probe must be an object")
        preset = probe.get("preset")
        operations = probe.get("operations")
        _require(type(preset) is int and preset in corpus_presets,
                 f"{path}: framebuffer probe preset is outside the oracle corpus")
        _require(_is_number(probe.get("hue_noise_phase")),
                 f"{path}: framebuffer hue-noise phase must be numeric")
        _require(isinstance(operations, list) and "FULL_FRAME" in operations,
                 f"{path}: each framebuffer probe must render a full frame")
        probe_presets.append(preset)
    _require(sorted(probe_presets) == sorted(corpus_presets),
             f"{path}: framebuffer probes must cover each oracle preset once")
    metrics = document.get("metrics")
    _require(isinstance(metrics, list) and metrics,
             f"{path}: metrics must be non-empty")
    domains = set()
    for metric in metrics:
        _require(isinstance(metric, dict), f"{path}: metric must be an object")
        domains.add(metric.get("domain"))
        measured = metric.get("measured_baseline")
        accepted = metric.get("accepted_limit")
        _require(_is_number(measured) and measured >= 0,
                 f"{path}: measured_baseline must be nonnegative")
        _require(_is_number(accepted) and accepted >= measured,
                 f"{path}: accepted_limit must cover measured_baseline")
        _require(bool(metric.get("unit")) and
                 bool(metric.get("limit_provenance")) and
                 bool(metric.get("measurement_provenance")),
                 f"{path}: metric unit and provenance are required")
        configurations = metric.get("configuration_baselines")
        resolutions = metric.get("resolution_baselines")
        if metric.get("domain") == "FRAMEBUFFER":
            _require(isinstance(configurations, dict) and
                     set(configurations) == {"native-debug", "wasm-release",
                                             "wasm-strict-fp"} and
                     all(_is_number(value) and value >= 0
                         for value in configurations.values()),
                     f"{path}: framebuffer configuration baselines are required")
            _require(max(configurations.values()) == measured,
                     f"{path}: framebuffer baseline must aggregate configurations")
            _require(isinstance(resolutions, dict) and
                     set(resolutions) == set(configurations) and
                     all(isinstance(values, dict) and
                         set(values) == {"96x20", "288x144"} and
                         all(_is_number(value) and value >= 0
                             for value in values.values())
                         for values in resolutions.values()),
                     f"{path}: framebuffer resolution baselines are required")
            _require(all(max(resolutions[configuration].values()) == baseline
                         for configuration, baseline in configurations.items()),
                     f"{path}: configuration baselines must aggregate resolutions")
        else:
            _require(configurations is None and resolutions is None,
                     f"{path}: only framebuffer metrics have configuration baselines")
    _require("FRAMEBUFFER" in domains,
             f"{path}: a final-framebuffer metric is required")


def load_and_validate(directory: Path) -> tuple[dict, list[dict]]:
    schema_path = directory / "schema.json"
    schema = _load(schema_path)
    _require(schema.get("$schema") ==
             "https://json-schema.org/draft/2020-12/schema",
             f"{schema_path}: unsupported JSON Schema draft")
    _require(isinstance(schema.get("$defs"), dict),
             f"{schema_path}: $defs must be an object")
    programs_path = directory / "programs.json"
    programs = _load(programs_path)
    _validate_schema(programs, schema, schema, str(programs_path))
    _validate_programs(programs, programs_path)
    oracles = []
    for filename in ORACLE_FILES:
        path = directory / filename
        oracle = _load(path)
        _validate_schema(oracle, schema, schema, str(path))
        _validate_oracle(oracle, path)
        _require(oracle["corpus"]["framebuffer"]["resolutions"] ==
                 programs["corpus"]["resolutions"],
                 f"{path}: framebuffer resolutions must match the capture roster")
        oracles.append(oracle)
    for oracle in oracles:
        _require(oracle["base_sha"] == programs["base_sha"],
                 f"{directory}: all manifests must use the same base_sha")
        _require(oracle["toolchains"] == programs["toolchains"],
                 f"{directory}: all manifests must use the same toolchain pins")
    return programs, oracles


def _canonical_bytes(documents: list[dict]) -> bytes:
    return b"\n".join(json.dumps(document, sort_keys=True,
                                 separators=(",", ":")).encode("utf-8")
                      for document in documents)


def manifest_sha256(programs: dict, oracles: list[dict]) -> str:
    return hashlib.sha256(_canonical_bytes([programs, *oracles])).hexdigest()


def generate_header(programs: dict, oracles: list[dict]) -> str:
    digest = manifest_sha256(programs, oracles)
    entries = []
    for program in programs["programs"]:
        mask = sum(1 << preset for preset in program["presets"])
        topology = ", ".join(str(program["topology_key"][field])
                             for field in TOPOLOGY_FIELDS)
        entries.append(
            f'    {{"{program["id"]}", 0x{mask:03x}, {{{{{topology}}}}}}},')
    metric_entries = []
    for oracle in oracles:
        for metric in oracle["metrics"]:
            measured = metric["measured_baseline"]
            metric_entries.append(
                "    {"
                f'"{oracle["oracle_id"]}", "{metric["domain"]}", '
                f'"{metric["aggregation"]}", '
                f'{str(measured is not None).lower()}, '
                f'{measured if measured is not None else 0.0}, '
                f'{metric["accepted_limit"]}'
                "},")
    return "\n".join((
        "#pragma once",
        "",
        "#include <array>",
        "#include <cstdint>",
        "#include <string_view>",
        "",
        "namespace PullbackManifest {",
        "struct ProgramEntry {",
        "  std::string_view id;",
        "  uint32_t preset_mask;",
        f"  std::array<uint8_t, {len(TOPOLOGY_FIELDS)}> topology_key;",
        "};",
        "struct OracleMetric {",
        "  std::string_view oracle_id;",
        "  std::string_view domain;",
        "  std::string_view aggregation;",
        "  bool measured;",
        "  float measured_baseline;",
        "  float accepted_limit;",
        "};",
        "",
        f'inline constexpr std::string_view BASE_SHA = "{programs["base_sha"]}";',
        f'inline constexpr std::string_view MANIFEST_SHA256 = "{digest}";',
        f"inline constexpr std::array<ProgramEntry, {len(entries)}> PROGRAMS{{{{",
        *entries,
        "}};",
        f"inline constexpr std::array<OracleMetric, {len(metric_entries)}> ORACLE_METRICS{{{{",
        *metric_entries,
        "}};",
        "}",
        "",
    ))


def generate_runtime_manifest(programs: dict, oracles: list[dict],
                              capture_sha: str) -> str:
    _require(SHA_RE.fullmatch(capture_sha) is not None,
             "capture_sha must be a full lowercase Git SHA")
    runtime = copy.deepcopy(programs)
    runtime["capture_sha"] = capture_sha
    runtime["manifest_sha256"] = manifest_sha256(programs, oracles)
    return json.dumps(runtime, indent=2, sort_keys=True) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest-dir", type=Path,
                        default=Path(__file__).resolve().parents[1] /
                        "tests/data/pullback")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--runtime-output", type=Path)
    parser.add_argument("--capture-sha")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    try:
        programs, oracles = load_and_validate(args.manifest_dir)
    except ManifestError as error:
        parser.error(str(error))
    if args.validate_only:
        return 0
    if args.output is None and args.runtime_output is None:
        parser.error("--output or --runtime-output is required")
    if args.output is not None:
        output = generate_header(programs, oracles)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output, encoding="utf-8", newline="\n")
    if args.runtime_output is not None:
        if args.capture_sha is None:
            parser.error("--runtime-output requires --capture-sha")
        try:
            output = generate_runtime_manifest(programs, oracles,
                                               args.capture_sha)
        except ManifestError as error:
            parser.error(str(error))
        args.runtime_output.parent.mkdir(parents=True, exist_ok=True)
        args.runtime_output.write_text(output, encoding="utf-8", newline="\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
