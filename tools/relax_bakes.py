#!/usr/bin/env python3
"""Validate or regenerate authoritative Teensy relax-bake assets."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
ASSET = ROOT / "core/mesh/relax_bakes_generated.h"
MANIFEST = ROOT / "core/mesh/relax_bakes_manifest.json"
REBAKE_COMMAND = "python tools/relax_bakes.py rebake --port COM4"
COMPILER_IDENTITY = "Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86); 15.2.1 20251203"
ARITHMETIC_FLAGS = (
    "-O3 -ffast-math -fno-finite-math-only -fno-unswitch-loops "
    "-std=gnu++20"
)
INPUT_PATHS = (
    "core/mesh/conway.h",
    "core/mesh/hankin.h",
    "core/mesh/mesh.h",
    "core/mesh/solids.h",
    "core/math/3dmath.h",
    "core/math/geometry.h",
    "core/engine/platform.h",
    "platformio.ini",
    "tools/relax_bake_check.py",
    "tools/relax_bakes.py",
)

DATA_RE = re.compile(r"^RELAX_BAKE_DATA(?:4|2) (.+)$")

# Authoritative source identities from the pinned Teensy extraction. The
# generator rejects additions or changes until a fresh device run updates this
# inventory deliberately.
SOURCES = {
    0xCDE6A3EF: ("truncated_icosahedron_ambo_converged", 90, 92, 360, 13,
                 0x0A20BCE1, 0xC91B1EDF),
    0x12ACE029: ("dodecahedron_ambo_bevel33_converged", 240, 122, 720, 50,
                 0x0E6068E9, 0x9A817BFD),
    0x06F43533: ("dodecahedron_bevel20_converged", 120, 62, 360, 17,
                 0x3B53F3C9, 0xE8B39AB3),
    0xA8661B09: ("truncated_icosidodecahedron_converged", 120, 62, 360, 17,
                 0x3B53F3C9, 0x95C0C28B),
    0x06CFFC5D: ("truncated_icosidodecahedron_bevel50_converged", 360, 362,
                 1440, 1064, 0x1F875855, 0x88D55C89),
    0x7F1CC361: ("dodecahedron_hankin_ambo_hankin_ambo_converged", 720, 362,
                 2160, 73, 0x236DE54D, 0x2F1A35E1),
    0x04C71362: ("snub_dodecahedron_converged", 60, 92, 300, 21,
                 0x004BC741, 0x2783988C),
    0x855BA2BF: ("icosahedron_snub_converged", 60, 92, 300, 21,
                 0x75B29ECD, 0x322E178F),
}

EXPECTED_CALLS = {
    (100, 0xCDE6A3EF),
    (217, 0xCDE6A3EF),
    (8, 0xCDE6A3EF),
    (100, 0x12ACE029),
    (100, 0x06F43533),
    (50, 0xA8661B09),
    (100, 0x06CFFC5D),
    (100, 0x7F1CC361),
    (50, 0x04C71362),
    (8, 0x855BA2BF),
}


def sha256_text(path: Path) -> str:
    """Hash text with Git-independent newline normalization."""
    contents = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(contents).hexdigest()


def vertex_hash(words: tuple[int, ...]) -> int:
    value = 2166136261
    for word in words:
        value = ((value ^ word) * 16777619) & 0xFFFFFFFF
    return value


def parse_capture(path: Path) -> tuple[dict[tuple[int, int], tuple], int]:
    complete: dict[tuple[int, int], tuple] = {}
    incomplete = 0
    meta = None
    words: list[int] = []

    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split()
        if parts[:1] == ["RELAX_BAKE_BEGIN"] and len(parts) in (7, 8, 9):
            if meta is not None:
                incomplete += 1
            if len(parts) == 7:
                requested, vertices, faces, indices = (
                    int(value) for value in parts[1:5]
                )
                actual = requested
                hashes = [int(value, 16) for value in parts[5:]]
            else:
                requested, actual, vertices, faces, indices = (
                    int(value) for value in parts[1:6]
                )
                hashes = [int(value, 16) for value in parts[6:]]
            if len(hashes) == 2:
                source_hash, output_hash = hashes
                topology_hash = None
            else:
                source_hash, topology_hash, output_hash = hashes
            meta = (requested, actual, vertices, faces, indices, source_hash,
                    topology_hash, output_hash)
            words = []
            continue
        data = DATA_RE.match(line)
        if data and meta is not None:
            words.extend(int(word, 16) for word in data.group(1).split())
            continue
        if line != "RELAX_BAKE_END" or meta is None:
            continue

        (requested, actual, vertices, faces, indices, source_hash,
         topology_hash, output_hash) = meta
        payload = tuple(words)
        meta = None
        if len(payload) != 3 * vertices:
            incomplete += 1
            continue
        if vertex_hash(payload) != output_hash:
            raise ValueError(f"{path}: output hash mismatch for {source_hash:08x}")
        key = (requested, source_hash)
        record = (actual, vertices, faces, indices, topology_hash, output_hash,
                  payload)
        previous = complete.get(key)
        if previous is not None and previous != record:
            raise ValueError(f"{path}: nondeterministic duplicate {key}")
        complete[key] = record

    if meta is not None:
        incomplete += 1
    missing = EXPECTED_CALLS - complete.keys()
    extra = complete.keys() - EXPECTED_CALLS
    if missing or extra:
        raise ValueError(f"{path}: capture inventory missing={missing}, extra={extra}")
    return complete, incomplete


def validate_captures(shipping: Path, o3: Path) -> dict[int, tuple[int, ...]]:
    ship, ship_incomplete = parse_capture(shipping)
    opt, opt_incomplete = parse_capture(o3)
    if ship != opt:
        differing = sorted(key for key in EXPECTED_CALLS if ship[key] != opt[key])
        raise ValueError(f"shipping/O3 payload mismatch: {differing}")

    payloads: dict[int, tuple[int, ...]] = {}
    for (requested, source_hash), record in ship.items():
        (actual, vertices, faces, indices, topology_hash, output_hash,
         payload) = record
        (name, expected_v, expected_f, expected_i, expected_actual,
         expected_topology, expected_out) = SOURCES[source_hash]
        observed = (vertices, faces, indices, actual, output_hash)
        expected = (expected_v, expected_f, expected_i, expected_actual,
                    expected_out)
        if observed != expected:
            raise ValueError(
                f"{name}: converged metadata changed: observed={observed}, "
                f"expected={expected}"
            )
        if topology_hash is not None and topology_hash != expected_topology:
            raise ValueError(
                f"{name}: topology hash changed: {topology_hash:08x} != "
                f"{expected_topology:08x}"
            )
        previous = payloads.get(source_hash)
        if previous is not None and previous != payload:
            raise ValueError(f"{name}: former caps did not converge identically")
        payloads[source_hash] = payload

    print(
        f"validated {len(EXPECTED_CALLS)} call signatures, "
        f"{len(payloads)} unique payloads; ignored redundant incomplete blocks "
        f"shipping={ship_incomplete}, o3={opt_incomplete}"
    )
    return payloads


def emit_header(payloads: dict[int, tuple[int, ...]]) -> str:
    lines = [
        "/* Generated by tools/relax_bakes.py from pinned Teensy 4.0 captures. */",
        "#pragma once",
        "",
        "// clang-format off",
        "namespace RelaxBakes {",
        "",
    ]
    for source_hash, values in SOURCES.items():
        name, vertices, faces, indices, actual, topology_hash, output_hash = values
        words = payloads[source_hash]
        lines.append(f"inline const uint32_t {name}_bits[] PROGMEM = {{")
        for start in range(0, len(words), 8):
            chunk = ", ".join(f"0x{word:08x}u" for word in words[start:start + 8])
            lines.append(f"    {chunk},")
        lines.extend(
            [
                "};",
                f"inline constexpr MeshOps::RelaxBake {name} = {{",
                f"    {name}_bits, {vertices}, {faces}, {indices}, {actual},",
                f"    0x{topology_hash:08x}u, 0x{source_hash:08x}u,",
                f"    0x{output_hash:08x}u}};",
                "",
            ]
        )
    lines.extend(["} // namespace RelaxBakes", "// clang-format on", ""])
    return "\n".join(lines)


def write_manifest(shipping: Path, o3: Path) -> None:
    manifest = {
        "schema": 1,
        "rebake_command": REBAKE_COMMAND,
        "compiler_identity": COMPILER_IDENTITY,
        "arithmetic_flags": ARITHMETIC_FLAGS,
        "asset_sha256": sha256_text(ASSET),
        "shipping_capture_sha256": sha256_text(shipping),
        "o3_capture_sha256": sha256_text(o3),
        "inputs": {path: sha256_text(ROOT / path) for path in INPUT_PATHS},
        "datasets": {
            f"{source_hash:08x}": {
                "name": values[0],
                "vertices": values[1],
                "faces": values[2],
                "indices": values[3],
                "convergence_iterations": values[4],
                "topology_hash": f"{values[5]:08x}",
                "output_hash": f"{values[6]:08x}",
            }
            for source_hash, values in SOURCES.items()
        },
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def emit(args: argparse.Namespace) -> None:
    shipping = Path(args.shipping_log).resolve()
    o3 = Path(args.o3_log).resolve()
    payloads = validate_captures(shipping, o3)
    ASSET.write_text(emit_header(payloads), encoding="utf-8", newline="\n")
    write_manifest(shipping, o3)
    raw_bytes = sum(len(payload) * 4 for payload in payloads.values())
    print(f"wrote {ASSET.relative_to(ROOT)}: {raw_bytes} raw payload bytes")
    print(f"wrote {MANIFEST.relative_to(ROOT)}")


def stale(message: str) -> None:
    raise SystemExit(
        f"relax bake assets are stale: {message}\n"
        f"Rebuild on the pinned Teensy 4.0 toolchain:\n  {REBAKE_COMMAND}"
    )


def check(_: argparse.Namespace) -> None:
    if not ASSET.exists() or not MANIFEST.exists():
        stale("generated asset or manifest is missing")
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if manifest.get("compiler_identity") != COMPILER_IDENTITY:
        stale("compiler identity changed")
    if manifest.get("arithmetic_flags") != ARITHMETIC_FLAGS:
        stale("relax arithmetic flags changed")
    if manifest.get("asset_sha256") != sha256_text(ASSET):
        stale("generated asset hash differs")
    for path in INPUT_PATHS:
        expected = manifest.get("inputs", {}).get(path)
        if expected != sha256_text(ROOT / path):
            stale(f"input changed: {path}")
    print("relax bake freshness: OK")


def verify_compiler(compiler: str) -> None:
    version = subprocess.check_output([compiler, "--version"], text=True)
    if ("Arm GNU Toolchain 15.2.Rel1 (Build arm-15.86)" not in version or
            "15.2.1 20251203" not in version):
        raise SystemExit(f"unexpected PlatformIO ARM compiler {compiler}:\n{version}")


def compiler(args: argparse.Namespace) -> None:
    verify_compiler(args.path)
    print(f"relax bake compiler: OK ({args.path})")


def rebake(args: argparse.Namespace) -> None:
    env = os.environ.copy()
    env["HS_TEENSY_PORT"] = args.port
    env["HS_PROFILE_TREE"] = str(ROOT).replace("\\", "/")
    env["HS_DEVICE_WAIT"] = "900"
    flags = [
        "-D HS_PROFILE_TRANS_SPEED=8",
        "-D HS_PROFILE_EPOCH_REVS=2560",
        "-D HS_RELAX_BAKE_EXTRACT",
    ]
    captures = {}
    for profile_env, tag in (("profile", "ship"), ("profile_o3", "o3")):
        command = [
            "bash", "tools/profile_one.sh", "IslamicStars", profile_env,
            "80", "16", *flags,
        ]
        subprocess.run(command, cwd=ROOT, env=env, check=True)
        source = ROOT / f"build/prof/islamicstars_{tag}.log"
        destination = ROOT / f"build/prof/islamicstars_{tag}_relax_converged_extract.log"
        shutil.copyfile(source, destination)
        captures[tag] = destination
    emit(argparse.Namespace(shipping_log=captures["ship"], o3_log=captures["o3"]))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("check").set_defaults(func=check)
    compiler_parser = subparsers.add_parser("compiler")
    compiler_parser.add_argument("--path", required=True)
    compiler_parser.set_defaults(func=compiler)
    emit_parser = subparsers.add_parser("emit")
    emit_parser.add_argument("--shipping-log", required=True)
    emit_parser.add_argument("--o3-log", required=True)
    emit_parser.set_defaults(func=emit)
    rebake_parser = subparsers.add_parser("rebake")
    rebake_parser.add_argument("--port", required=True)
    rebake_parser.set_defaults(func=rebake)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
