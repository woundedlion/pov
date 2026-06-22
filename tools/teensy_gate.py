#!/usr/bin/env python3
"""Teensy 4 firmware size + memory-layout gate (pure parser / classifier).

This module is the testable heart of docs/teensy_ci_gate_spec.md. It holds NO
PlatformIO dependency so it runs as an ordinary host Python module (stdlib only,
matching the repo's scripts/ convention) and is exercised by the golden /
deliberately-broken fixtures under tools/teensy_gate_tests/. The thin PlatformIO
post-build wrapper that feeds it real toolchain output is tools/teensy_gate_extra.py.

What it does (spec §7.3, §7.4, §8):
  * region totals  — parse `teensy_size` (or `arm-none-eabi-size -A`) and compare
    each region's used bytes against a per-target ceiling, and the DTCM "free for
    local variables" stack headroom against a floor.
  * layout         — classify the framebuffer / arena / reaction-graph symbols by
    their LOAD ADDRESS against the Teensy 4 memory map (NOT an `nm` type letter:
    DTCM .bss and OCRAM .dmabuffers are both NOBITS and `nm` cannot tell them
    apart), and assert each lands in the region it must, with the arena's
    MAGNITUDE pinned near 335 KB so a leaked HS_TEST_BUILD (8 MB) arena fails.
  * fail-loud       — a configured layout symbol that is NOT FOUND in the ELF is a
    violation, never a silent skip: a name that never matches would make the
    invariant never fire (false-green), the exact trap spec §7.4 warns about.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

# ---------------------------------------------------------------------------
# Teensy 4 (i.MX RT1062) memory map — the load-bearing address buckets (§7.4).
# Classify a symbol by where it LIVES (its VMA), not by section name or nm letter.
#   ITCM  0x0000_0000  fast code              (part of RAM1)
#   DTCM  0x2000_0000  fast data + stack      (part of RAM1) — the 335 KB arena
#   OCRAM 0x2020_0000  DMAMEM + heap          (RAM2)         — the framebuffers
#   FLASH 0x6000_0000  code + const/PROGMEM   (2 MiB on T4.0)— the reaction graph
# Half-open [lo, hi) ranges. RAM1 in the budget == ITCM ∪ DTCM.
# ---------------------------------------------------------------------------
MEMORY_MAP: tuple[tuple[str, int, int], ...] = (
    ("ITCM", 0x00000000, 0x00080000),   # 512 KiB window; ITCM is carved from RAM1
    ("DTCM", 0x20000000, 0x20080000),   # 512 KiB window
    ("OCRAM", 0x20200000, 0x20280000),  # 512 KiB RAM2
    ("FLASH", 0x60000000, 0x60200000),  # 2 MiB T4.0 flash
)


def region_for_address(addr: int) -> str:
    """Bucket a load address into a Teensy 4 memory region, or 'OTHER'."""
    for name, lo, hi in MEMORY_MAP:
        if lo <= addr < hi:
            return name
    return "OTHER"


# ---------------------------------------------------------------------------
# Parsed-data containers
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Symbol:
    """One `readelf -s` symbol-table row."""

    num: int
    value: int          # VMA (load address)
    size: int
    type: str
    bind: str
    ndx: str            # section index as printed (number, or UND/ABS/COM)
    name: str           # possibly-mangled linkage name

    @property
    def region(self) -> str:
        return region_for_address(self.value)


@dataclass
class Violation:
    """One gate failure, rendered as a GitHub `::error::` annotation."""

    code: str
    message: str


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------
# teensy_size summary lines (§7.3), e.g.
#   teensy_size: FLASH: code:62788, data:13684, headers:8460   free for files: 1979136
#   teensy_size: RAM1: variables:343040, code:62240, padding:30496   free for local variables: 88512
#   teensy_size: RAM2: variables:497920   free for malloc/new: 26368
# The component blob and the "free for ..." figure are separated by run(s) of
# whitespace whose width is not contractual — teensy_size has used 2+ spaces but
# a single-space variant is a valid format the parser must not silently miss
# (yielding a "region-missing" gate failure). `.*?` is lazy and "free for" is a
# unique literal on the line, so `\s+` cannot over-consume the blob's own single
# spaces (e.g. ", data:").
_TS_REGION_RE = re.compile(
    r"\bteensy_size:\s*(FLASH|RAM1|RAM2):\s*(.*?)\s+free for [^:]+:\s*(\d+)",
    re.IGNORECASE,
)
_TS_PAIR_RE = re.compile(r"(\w+)\s*:\s*(\d+)")


def parse_teensy_size(text: str) -> dict[str, dict[str, int]]:
    """Parse `teensy_size` stdout into per-region {used, free} byte totals.

    `used` is the sum of the region's components (code/data/headers/...); `free`
    is the region's "free for ..." figure. RAM1's free is the DTCM stack
    headroom (§7.4 #4); RAM2's free is the OCRAM heap room (§8).
    """
    out: dict[str, dict[str, int]] = {}
    for m in _TS_REGION_RE.finditer(text):
        region = m.group(1).lower()
        components = {k.lower(): int(v) for k, v in _TS_PAIR_RE.findall(m.group(2))}
        out[region] = {"used": sum(components.values()), "free": int(m.group(3)),
                       "components": components}
    return out


# `arm-none-eabi-size -A -x` section row, e.g.
#   .text.progmem   0xb948   0x60000200
_SIZE_A_RE = re.compile(r"^(\.\S+)\s+(0x[0-9a-fA-F]+|\d+)\s+(0x[0-9a-fA-F]+|\d+)\s*$")


def parse_size_a(text: str) -> list[tuple[str, int, int]]:
    """Parse `arm-none-eabi-size -A` into [(section, size, addr)] rows."""
    rows: list[tuple[str, int, int]] = []
    for line in text.splitlines():
        m = _SIZE_A_RE.match(line.strip())
        if m:
            rows.append((m.group(1), int(m.group(2), 0), int(m.group(3), 0)))
    return rows


def region_totals_from_size_a(text: str) -> dict[str, int]:
    """Sum `size -A` sections into ITCM/DTCM/OCRAM/FLASH totals by VMA.

    A documented fallback / cross-check for parse_teensy_size (§7.3, §9). NOTE:
    it buckets by VMA, so an ITCM/.fast code section is counted under RAM1 here
    even though its initialized copy also consumes flash — i.e. this UNDERCOUNTS
    flash relative to teensy_size, which does the correct LMA accounting. Use it
    for RAM1/RAM2 placement cross-checks, not as the authoritative flash ceiling.
    """
    totals: dict[str, int] = {}
    for _name, size, addr in parse_size_a(text):
        if addr == 0 and size == 0:
            continue
        totals[region_for_address(addr)] = totals.get(region_for_address(addr), 0) + size
    return totals


# `readelf -s`/`-sW` symbol row, e.g.
#    124: 20200000 248832 OBJECT  GLOBAL DEFAULT    8 _ZN6Effect8buffer_aE
# The Value is hex (no prefix); the Size is decimal for small objects but HEX
# (0x-prefixed) for large ones in arm-none-eabi readelf 11.3.1 — accept either.
_READELF_SYM_RE = re.compile(
    r"^\s*(\d+):\s+([0-9a-fA-F]+)\s+(0x[0-9a-fA-F]+|\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$"
)


def parse_readelf_symbols(text: str) -> list[Symbol]:
    """Parse `arm-none-eabi-readelf -s` (or -sW) symbol-table output."""
    syms: list[Symbol] = []
    for line in text.splitlines():
        m = _READELF_SYM_RE.match(line)
        if not m:
            continue
        syms.append(Symbol(
            num=int(m.group(1)),
            value=int(m.group(2), 16),
            size=int(m.group(3), 0),  # decimal or 0x-hex
            type=m.group(4),
            bind=m.group(5),
            ndx=m.group(7),
            name=m.group(8),
        ))
    return syms


# `readelf -S`/`-SW` section-header row, e.g.
#   [ 2] .text.progmem     PROGBITS        60000200 000200 00b948 00  AX  0   0 16
_READELF_SEC_RE = re.compile(
    r"^\s*\[\s*(\d+)\]\s+(\S+)\s+(\S+)\s+([0-9a-fA-F]+)\b"
)


def parse_readelf_sections(text: str) -> dict[str, tuple[str, int]]:
    """Parse `arm-none-eabi-readelf -S` into {ndx: (name, addr)}."""
    secs: dict[str, tuple[str, int]] = {}
    for line in text.splitlines():
        m = _READELF_SEC_RE.match(line)
        if m:
            secs[m.group(1)] = (m.group(2), int(m.group(4), 16))
    return secs


def section_name(symbol: Symbol, sections: dict[str, tuple[str, int]]) -> str:
    """Human-readable section for a symbol's Ndx (for the report message)."""
    entry = sections.get(symbol.ndx)
    return entry[0] if entry else symbol.ndx


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------
@dataclass
class GateResult:
    env: str
    violations: list[Violation] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def passed(self) -> bool:
        return not self.violations


def evaluate(
    env: str,
    budget: dict,
    sizes: dict[str, dict[str, int]],
    symbols: list[Symbol],
    sections: dict[str, tuple[str, int]] | None = None,
) -> GateResult:
    """Apply a target's budget to its measured sizes + symbols (§7.3, §7.4, §8)."""
    sections = sections or {}
    result = GateResult(env=env)
    v = result.violations

    # --- Region ceilings + DTCM stack-headroom floor (§7.3, §7.4 #4, §8) ---
    for region, spec in budget.get("regions", {}).items():
        measured = sizes.get(region)
        if measured is None:
            v.append(Violation(
                "region-missing",
                f"{env}: no measured size for region '{region}' - "
                f"teensy_size output did not report it (parser or build broke)."))
            continue
        cap = spec.get("max_bytes")
        if cap is not None and measured["used"] > cap:
            v.append(Violation(
                "region-over-budget",
                f"{env}: {region.upper()} uses {measured['used']:,} B, "
                f"over the {cap:,} B ceiling (by {measured['used'] - cap:,} B)."))
        floor = spec.get("free_min_bytes")
        if floor is not None and measured.get("free", 0) < floor:
            v.append(Violation(
                "headroom-below-floor",
                f"{env}: {region.upper()} free-for-local-variables "
                f"{measured.get('free', 0):,} B is below the {floor:,} B floor "
                f"(stack headroom squeezed)."))

    # --- Layout invariants: symbol -> region (+ magnitude) (§7.4 #1-#3) ---
    for key, spec in budget.get("symbols", {}).items():
        name = spec["name"]
        # Consider only DEFINED symbol-table rows (real section + non-zero size).
        # readelf -s can also carry a same-named UND reference row (size 0, ndx
        # UND, null value); including it would spuriously trip "symbol-too-small"
        # on its 0 size and mis-derive a region from address 0. The definition is
        # the row that carries the object, so a name present ONLY as UND means the
        # definition was removed/renamed -> still a hard "symbol-not-found" below.
        matches = [s for s in symbols
                   if s.name == name and s.ndx != "UND" and s.size > 0]
        if not matches:
            # Fail loud: a never-matching name silently disables its invariant.
            v.append(Violation(
                "symbol-not-found",
                f"{env}: layout symbol '{key}' ({name}) not found in the ELF - "
                f"renamed/removed? A missing symbol is a hard failure, never a "
                f"silent pass (spec 7.4)."))
            continue
        want_region = spec.get("region")
        for sym in matches:
            got_region = sym.region
            where = section_name(sym, sections)
            if want_region is not None and got_region != want_region:
                v.append(Violation(
                    "symbol-wrong-region",
                    f"{env}: '{key}' ({name}) is in {got_region} (section "
                    f"{where}, addr 0x{sym.value:08x}) but must be in "
                    f"{want_region}."))
            lo, hi = spec.get("min_bytes"), spec.get("max_bytes")
            if lo is not None and sym.size < lo:
                v.append(Violation(
                    "symbol-too-small",
                    f"{env}: '{key}' ({name}) is {sym.size:,} B, below its "
                    f"{lo:,} B floor - expected magnitude regression."))
            if hi is not None and sym.size > hi:
                v.append(Violation(
                    "symbol-too-large",
                    f"{env}: '{key}' ({name}) is {sym.size:,} B, above its "
                    f"{hi:,} B cap (e.g. an HS_TEST_BUILD 8 MB arena leak, spec 7.4 #1)."))

    return result


# ---------------------------------------------------------------------------
# CLI (standalone / local debugging on captured toolchain output)
# ---------------------------------------------------------------------------
def _strip_jsonc_comments(text: str) -> str:
    """Remove // line and /* */ block comments from JSONC text, leaving any such
    sequences that occur INSIDE a JSON string value untouched.

    A character scanner is used rather than a regex on purpose: the naive
    ``re.sub`` this replaced would corrupt a ``//`` or ``/*`` that appears inside
    a string (a path, URL, or note), silently mangling a future budgets entry.
    """
    out = []
    i, n = 0, len(text)
    in_string = False
    while i < n:
        c = text[i]
        if in_string:
            out.append(c)
            if c == "\\" and i + 1 < n:
                # Escaped char: copy verbatim; it can't terminate the string.
                out.append(text[i + 1])
                i += 2
                continue
            if c == '"':
                in_string = False
            i += 1
        elif c == '"':
            in_string = True
            out.append(c)
            i += 1
        elif c == "/" and i + 1 < n and text[i + 1] == "/":
            i += 2  # // line comment: drop through end of line (keep the \n)
            while i < n and text[i] != "\n":
                i += 1
        elif c == "/" and i + 1 < n and text[i + 1] == "*":
            i += 2  # /* */ block comment: drop through the closing */
            while i + 1 < n and not (text[i] == "*" and text[i + 1] == "/"):
                i += 1
            i += 2
        else:
            out.append(c)
            i += 1
    return "".join(out)


def load_budgets(path: str | Path) -> dict:
    """Load tools/teensy_budgets.json, tolerating // and /* */ comments."""
    raw = Path(path).read_text(encoding="utf-8")
    return json.loads(_strip_jsonc_comments(raw))


def render_report(result: GateResult, *, github: bool = False) -> str:
    """Render a pass/fail report; GitHub mode prefixes `::error::` annotations."""
    lines = []
    if result.passed:
        lines.append(f"[teensy-gate] {result.env}: PASS - all region budgets and "
                     f"layout invariants satisfied.")
    else:
        lines.append(f"[teensy-gate] {result.env}: FAIL "
                     f"({len(result.violations)} violation(s)):")
        for viol in result.violations:
            prefix = "::error::" if github else "  - "
            lines.append(f"{prefix}{viol.message}")
    for note in result.notes:
        lines.append(f"  note: {note}")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Teensy 4 size/layout gate (parser).")
    p.add_argument("--env", required=True, help="budget key, e.g. holosphere")
    p.add_argument("--budgets", default="tools/teensy_budgets.json")
    p.add_argument("--teensy-size", help="file with captured teensy_size stdout")
    p.add_argument("--size-a", help="file with captured `size -A` output (fallback)")
    p.add_argument("--readelf-syms", required=True, help="file with `readelf -s` output")
    p.add_argument("--readelf-secs", help="file with `readelf -S` output")
    p.add_argument("--github", action="store_true", help="emit ::error:: annotations")
    args = p.parse_args(argv)

    budgets = load_budgets(args.budgets)
    if args.env not in budgets:
        print(f"::error::no budget for env '{args.env}' in {args.budgets}",
              file=sys.stderr)
        return 2

    if args.teensy_size:
        sizes = parse_teensy_size(Path(args.teensy_size).read_text(encoding="utf-8"))
    elif args.size_a:
        totals = region_totals_from_size_a(
            Path(args.size_a).read_text(encoding="utf-8"))
        ram1 = totals.get("ITCM", 0) + totals.get("DTCM", 0)
        sizes = {
            "flash": {"used": totals.get("FLASH", 0), "free": 0},
            "ram1": {"used": ram1, "free": 0x80000 - ram1},
            "ram2": {"used": totals.get("OCRAM", 0), "free": 0x80000 - totals.get("OCRAM", 0)},
        }
    else:
        p.error("one of --teensy-size or --size-a is required")

    symbols = parse_readelf_symbols(Path(args.readelf_syms).read_text(encoding="utf-8"))
    sections = (parse_readelf_sections(Path(args.readelf_secs).read_text(encoding="utf-8"))
                if args.readelf_secs else {})

    result = evaluate(args.env, budgets[args.env], sizes, symbols, sections)
    print(render_report(result, github=args.github))
    return 0 if result.passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
