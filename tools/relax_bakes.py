#!/usr/bin/env python3
"""Regenerate the host relax-bake asset (core/mesh/relax_bakes_generated.h).

The bake payloads are produced on the host: the `relax_bake_gen` harness
(tools/relax_bake_harness.cpp, built with HS_RELAX_BAKE_EXTRACT) runs every
bake-bearing SolidBuilder recipe and logs a RELAX_BAKE block per payload. This
script parses that stream and writes the generated header. Because host relax
is deterministic, the emitted bits load unchanged on host and device; the
`relax_bake_verify` ctest (HS_RELAX_BAKE_VERIFY) re-derives them and asserts
bit-exact equality, so a stale asset fails the suite.

That gate covers the payload values only. `check` covers the file's form —
banner, chunking, declaration layout — by re-emitting the header from a fresh
dump and diffing the full text, so a legitimate regeneration cannot arrive
buried in a reformat the emitter drifted into meanwhile.

Usage:
    <build>/relax_bake_gen | python tools/relax_bakes.py emit --stdin
    python tools/relax_bakes.py emit --dump captured_dump.txt
    <build>/relax_bake_gen | python tools/relax_bakes.py check --stdin
"""

from __future__ import annotations

import argparse
import difflib
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
ASSET = ROOT / "core/mesh/relax_bakes_generated.h"


def unterminated_block(meta: dict, words: list[int]) -> ValueError:
    """Error for a bake the harness died inside; dropping it emits a short asset."""
    return ValueError(
        f"{meta['name']}: unterminated block ({len(words)} of "
        f"{3 * meta['vertices']} words, no RELAX_BAKE_END)"
    )


def parse_dump(text: str) -> list[dict]:
    """Parse RELAX_BAKE blocks into per-name payloads (first-seen order)."""
    bakes: dict[str, dict] = {}
    order: list[str] = []
    meta: dict | None = None
    words: list[int] = []
    for line in text.splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "RELAX_BAKE_BEGIN":
            if meta is not None:
                raise unterminated_block(meta, words)
            name, iterations, v, f, i, topo, out = parts[1:8]
            meta = {
                "name": name,
                "iterations": int(iterations),
                "vertices": int(v),
                "faces": int(f),
                "indices": int(i),
                "topology_hash": int(topo, 16),
                "output_hash": int(out, 16),
            }
            words = []
        elif parts[0] in ("RELAX_BAKE_DATA", "RELAX_BAKE_END") and meta is None:
            raise ValueError(f"{parts[0]} outside a block: dump starts mid-bake")
        elif parts[0] == "RELAX_BAKE_DATA":
            words.extend(int(w, 16) for w in parts[1:4])
        elif parts[0] == "RELAX_BAKE_END":
            if len(words) != 3 * meta["vertices"]:
                raise ValueError(
                    f"{meta['name']}: {len(words)} words for "
                    f"{meta['vertices']} vertices"
                )
            if vertex_hash(words) != meta["output_hash"]:
                raise ValueError(f"{meta['name']}: output hash mismatch")
            meta["bits"] = words
            existing = bakes.get(meta["name"])
            if existing is None:
                order.append(meta["name"])
                bakes[meta["name"]] = meta
            elif existing["bits"] != meta["bits"]:
                raise ValueError(f"{meta['name']}: nondeterministic duplicate")
            meta = None
    if meta is not None:
        raise unterminated_block(meta, words)
    if not order:
        raise ValueError("no RELAX_BAKE blocks found in dump")
    return [bakes[name] for name in order]


def vertex_hash(words: list[int]) -> int:
    value = 2166136261
    for word in words:
        value = ((value ^ word) * 16777619) & 0xFFFFFFFF
    return value


def emit_header(bakes: list[dict]) -> str:
    lines = [
        "/*",
        " * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.",
        " * Licensed under the PolyForm Noncommercial License 1.0.0",
        " * GENERATED FILE - DO NOT EDIT. Regenerate from host relax bakes with:",
        " *   <build>/relax_bake_gen | python tools/relax_bakes.py emit --stdin",
        " */",
        "#pragma once",
        "",
        "#include <cstdint>",
        "#include <iterator>",
        '#include "mesh/conway.h"',
        "",
        "// clang-format off",
        "namespace Solids {",
        "namespace RelaxBakes {",
        "",
    ]
    for bake in bakes:
        name = bake["name"]
        words = bake["bits"]
        lines.append(f"inline const uint32_t {name}_bits[] HS_PROGMEM_UNIQUE({name}_bits) = {{")
        for start in range(0, len(words), 8):
            chunk = ", ".join(f"0x{w:08x}u" for w in words[start:start + 8])
            lines.append(f"    {chunk},")
        lines.extend([
            "};",
            f"inline constexpr MeshOps::RelaxBake {name} = {{",
            f'    .name = "{name}", .vertex_bits = {name}_bits,',
            f'    .vertex_count = {bake["vertices"]}, '
            f'.face_count = {bake["faces"]}, '
            f'.index_count = {bake["indices"]}, '
            f'.iterations = {bake["iterations"]},',
            f'    .topology_hash = 0x{bake["topology_hash"]:08x}u,',
            f'    .output_hash = 0x{bake["output_hash"]:08x}u}};',
            f"static_assert(std::size({name}_bits) == 3u * {name}.vertex_count);",
            "",
        ])
    lines.extend([
        "} // namespace RelaxBakes",
        "} // namespace Solids",
        "// clang-format on",
        "",
    ])
    return "\n".join(lines)


def read_dump(args: argparse.Namespace) -> str:
    """Read the harness dump from stdin or the named file."""
    if args.stdin:
        return sys.stdin.read()
    return Path(args.dump).read_text(encoding="utf-8", errors="replace")


def emit(args: argparse.Namespace) -> None:
    bakes = parse_dump(read_dump(args))
    ASSET.write_text(emit_header(bakes), encoding="utf-8", newline="\n")
    raw = sum(len(b["bits"]) * 4 for b in bakes)
    print(f"wrote {ASSET.relative_to(ROOT)}: {len(bakes)} bakes, {raw} bytes")


def check(args: argparse.Namespace) -> None:
    """Re-emit the header from a dump and diff the full text against the asset."""
    bakes = parse_dump(read_dump(args))
    generated = emit_header(bakes)
    # read_text decodes universal newlines, so a CRLF checkout compares equal.
    committed = ASSET.read_text(encoding="utf-8")
    name = ASSET.relative_to(ROOT).as_posix()
    if generated != committed:
        diff = difflib.unified_diff(committed.splitlines(),
                                    generated.splitlines(),
                                    fromfile=name, tofile="regenerated",
                                    lineterm="")
        for line in list(diff)[:40]:
            print(line, file=sys.stderr)
        sys.exit(
            f"{name} is out of sync with tools/relax_bakes.py; regenerate with: "
            "relax_bake_gen | python tools/relax_bakes.py emit --stdin")
    print(f"{name} matches the regenerated header in full: "
          f"{len(bakes)} bake payload(s).")


def add_dump_source(parser: argparse.ArgumentParser) -> None:
    """Add the mutually exclusive --stdin / --dump input selection."""
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--stdin", action="store_true",
                       help="read the harness dump from stdin")
    group.add_argument("--dump", help="read the harness dump from a file")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    emit_parser = sub.add_parser("emit", help="write the header from a dump")
    add_dump_source(emit_parser)
    emit_parser.set_defaults(func=emit)
    check_parser = sub.add_parser(
        "check", help="diff the header re-emitted from a dump against the asset")
    add_dump_source(check_parser)
    check_parser.set_defaults(func=check)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
