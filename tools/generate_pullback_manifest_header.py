#!/usr/bin/env python3
"""Validate pullback manifests and generate their native-test header."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
from pathlib import Path


PROGRAM_IDS = (
    "BONNE_KALEIDOSCOPE_LATTICE_MIRROR",
    "GLITCH_NOISE_GRID_WAVE_SHEAR",
    "KALEIDOSCOPE_TWIN_WAVE_INNER_MIRROR",
    "GNOMONIC_KALEIDOSCOPE_GRID_MIRROR",
    "GNOMONIC_GLITCH_GRID_MIRROR",
    "PEIRCE_DODECAHEDRAL_GRID",
    "GNOMONIC_DODECAHEDRAL_GRID_WAVE_MIRROR",
    "GNOMONIC_AFFINE_LATTICE_CONTOUR",
    "SINUSOIDAL_CURL_LATTICE",
    "STEREOGRAPHIC_PRISM_POLAR_WAVE_LATTICE",
    "GNOMONIC_DODECAHEDRAL_GRID_VECTOR_MIRROR",
)
ORACLE_FILES = ("peirce_fast_square.json", "hue_rotation_noise_luts.json")
SHA_RE = re.compile(r"^[0-9a-f]{40}$")
CASE_IDS = {"default", "endpoint_min", "endpoint_max", "interior"}


class ManifestError(ValueError):
    """A manifest violates the checked schema or cross-file contract."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ManifestError(message)


def _load(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ManifestError(f"{path}: {error}") from error
    _require(isinstance(value, dict), f"{path}: root must be an object")
    return value


def _validate_common(document: dict, path: Path) -> None:
    _require(document.get("schema_version") == 1,
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
    _require(corpus.get("version") == 1, f"{path}: corpus.version must be 1")
    _require(corpus.get("resolutions") == [[96, 20], [288, 144]],
             f"{path}: both roster resolutions are required")
    programs = document.get("programs")
    _require(isinstance(programs, list), f"{path}: programs must be an array")
    _require(len(programs) == len(PROGRAM_IDS),
             f"{path}: expected {len(PROGRAM_IDS)} programs")
    ids = [program.get("id") for program in programs
           if isinstance(program, dict)]
    _require(tuple(ids) == PROGRAM_IDS,
             f"{path}: program IDs or stable order differ from the compiled manifest")
    covered_presets: list[int] = []
    for program in programs:
        program_id = program["id"]
        presets = program.get("presets")
        _require(isinstance(presets, list) and presets,
                 f"{path}: {program_id}.presets must be non-empty")
        _require(all(isinstance(index, int) and 0 <= index < 12
                     for index in presets),
                 f"{path}: {program_id}.presets contains an invalid index")
        covered_presets.extend(presets)
        topology = program.get("topology_key")
        _require(isinstance(topology, dict) and topology,
                 f"{path}: {program_id}.topology_key must be non-empty")
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
        _require(0 <= framebuffer.get("max_channel_delta_u16", -1) <= 4,
                 f"{path}: {program_id} exceeds the channel-delta license")
        _require(0 <= framebuffer.get("max_differing_pixel_fraction", -1) <= 0.001,
                 f"{path}: {program_id} exceeds the differing-pixel license")
    _require(sorted(covered_presets) == list(range(12)),
             f"{path}: presets must be covered exactly once")


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
    _require(isinstance(corpus, dict) and corpus.get("version") == 1,
             f"{path}: a versioned deterministic corpus is required")
    _require(corpus.get("deterministic") is True,
             f"{path}: oracle corpus must be deterministic")
    metrics = document.get("metrics")
    _require(isinstance(metrics, list) and metrics,
             f"{path}: metrics must be non-empty")
    domains = set()
    for metric in metrics:
        _require(isinstance(metric, dict), f"{path}: metric must be an object")
        domains.add(metric.get("domain"))
        measured = metric.get("measured_baseline")
        accepted = metric.get("accepted_limit")
        if measured is None:
            _require(metric.get("domain") == "FRAMEBUFFER" and
                     metric.get("measurement_status") == "requires capture",
                     f"{path}: only uncaptured framebuffer metrics may be pending")
            _require(isinstance(accepted, (int, float)) and accepted >= 0,
                     f"{path}: pending metric accepted_limit must be nonnegative")
        else:
            _require(isinstance(measured, (int, float)) and measured >= 0,
                     f"{path}: measured_baseline must be nonnegative")
            _require(isinstance(accepted, (int, float)) and accepted >= measured,
                     f"{path}: accepted_limit must cover measured_baseline")
        _require(bool(metric.get("unit")) and bool(metric.get("limit_provenance")),
                 f"{path}: metric unit and limit provenance are required")
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
    _validate_programs(programs, programs_path)
    oracles = []
    for filename in ORACLE_FILES:
        path = directory / filename
        oracle = _load(path)
        _validate_oracle(oracle, path)
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


def generate_header(programs: dict, oracles: list[dict]) -> str:
    digest = hashlib.sha256(_canonical_bytes([programs, *oracles])).hexdigest()
    entries = []
    for program in programs["programs"]:
        mask = sum(1 << preset for preset in program["presets"])
        entries.append(f'    {{"{program["id"]}", 0x{mask:03x}}},')
    oracle_entries = [f'    "{oracle["oracle_id"]}",' for oracle in oracles]
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
        "  uint16_t preset_mask;",
        "};",
        "",
        f'inline constexpr std::string_view BASE_SHA = "{programs["base_sha"]}";',
        f'inline constexpr std::string_view MANIFEST_SHA256 = "{digest}";',
        f"inline constexpr std::array<ProgramEntry, {len(entries)}> PROGRAMS{{{{",
        *entries,
        "}};",
        f"inline constexpr std::array<std::string_view, {len(oracle_entries)}> ORACLES{{{{",
        *oracle_entries,
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
    runtime["manifest_sha256"] = hashlib.sha256(
        _canonical_bytes([programs, *oracles])).hexdigest()
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
