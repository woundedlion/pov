#!/usr/bin/env python3
"""Teensy 4 firmware size and memory-layout gate; stdlib only.

tools/teensy_budgets.json defines region ceilings, DTCM stack-headroom floors
and layout requirements. Layout checks use symbol VMAs, not nm type letters.
Missing configured symbols, malformed region lines, unknown budget keys and
missing required entries all fail.
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import NotRequired, TypedDict

from teensy_memory import (FLASH_PROGRAM_BYTES as FLASH_PROGRAM_BYTES,
                           FLEXRAM_BANK_BYTES as FLEXRAM_BANK_BYTES,
                           MEMORY_MAP as MEMORY_MAP,
                           region_for_address as region_for_address)
from teensy_budget_schema import (
    BudgetSchemaError as BudgetSchemaError,
    _strip_jsonc_comments as _strip_jsonc_comments,
    declares_components as declares_components,
    load_budgets as load_budgets,
    validate_budgets as validate_budgets,
)


def sibling_tool(compiler: str, suffix: str, search_path: str | None = None) -> str:
    """Resolve a GCC sibling tool, or an ARM binutils name on PATH."""
    resolved = shutil.which(compiler, path=search_path) or compiler
    stem, matched = re.subn(r"-(gcc|g\+\+|c\+\+)(\.exe)?$", "", resolved,
                            flags=re.IGNORECASE)
    if not matched:
        return "arm-none-eabi-" + suffix
    extension = ".exe" if resolved.lower().endswith(".exe") else ""
    return stem + "-" + suffix + extension


# ---------------------------------------------------------------------------
# Parsed-data containers
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Symbol:
    """One `readelf -s` symbol-table row."""

    num: int
    value: int          # VMA (runtime address)
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
# teensy_size summary lines, e.g.
#   teensy_size: FLASH: code:62788, data:13684, headers:8460   free for files: 1979136
#   teensy_size: RAM1: variables:343040, code:62240, padding:30496   free for local variables: 88512
#   teensy_size: RAM2: variables:497920   free for malloc/new: 26368
# The lazy `.*?` stops at "free for", a unique literal on the line, so `\s+`
# cannot over-consume the blob's own spaces.
_TS_REGION_RE = re.compile(
    r"\bteensy_size:\s*(FLASH|RAM1|RAM2):\s*(.*?)\s+free for [^:]+:\s*(-?\d+)",
    re.IGNORECASE,
)
_TS_PAIR_RE = re.compile(r"(\w+)\s*:\s*(\d+)")


class TeensySizeFormatError(ValueError):
    pass


class RegionSizes(TypedDict):
    """One measured memory region.

    `components` is present only when the figures came from `teensy_size`, which
    breaks each region down; the `size -A` fallback can only total sections, so
    it omits the key entirely.
    """

    used: int
    free: int
    components: NotRequired[dict[str, int]]


def parse_teensy_size(text: str) -> dict[str, RegionSizes]:
    """Parse `teensy_size` stdout into per-region {used, free} byte totals.

    `used` is the sum of the region's components (code/data/headers/...); `free`
    is the region's "free for ..." figure. RAM1's free is the DTCM stack
    headroom; RAM2's free is the OCRAM heap room.

    A region whose component blob yields no `name:bytes` pair raises: summing
    nothing gives `used` 0, which clears every ceiling and reports PASS.
    """
    out: dict[str, RegionSizes] = {}
    for m in _TS_REGION_RE.finditer(text):
        region = m.group(1).lower()
        if region in out:
            raise TeensySizeFormatError(
                f"region {region.upper()} appears more than once")
        components = {k.lower(): int(v) for k, v in _TS_PAIR_RE.findall(m.group(2))}
        if not components:
            raise TeensySizeFormatError(
                f"region {region.upper()} reports no component:bytes pairs in "
                f"{m.group(2)!r}")
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


# Non-allocated metadata sections report a VMA of 0 in `size -A`. Filter them by
# name, not by VMA == 0: real ITCM sections (.text.itcm) also run at 0x0, so a
# VMA test would drop live code while these consume no target memory.
_NON_ALLOC_SECTIONS = (
    ".ARM.attributes", ".comment", ".debug", ".note",
    ".symtab", ".strtab", ".shstrtab", ".stab",
)

# `size -A`'s `<file> :` header often starts with '.' (`.pio/build/...`); no
# section row ends in ':'.
_SIZE_A_FILE_HEADER_RE = re.compile(r"^\S+\s*:$")


class SizeAFormatError(ValueError):
    pass


def fallback_sizes_from_size_a(text: str) -> dict[str, RegionSizes]:
    """Validate `size -A` output and synthesize Teensy budget regions."""
    malformed = [line.strip() for line in text.splitlines()
                 if line.strip().startswith(".")
                 and not _SIZE_A_RE.match(line.strip())
                 and not _SIZE_A_FILE_HEADER_RE.match(line.strip())]
    if malformed:
        raise SizeAFormatError(
            f"malformed section row {malformed[0]!r}")

    rows = parse_size_a(text)
    allocated = [
        (name, size, addr) for name, size, addr in rows
        if size > 0 and not name.startswith(_NON_ALLOC_SECTIONS)
        and region_for_address(addr) != "OTHER"
    ]
    if not allocated:
        raise SizeAFormatError(
            "no recognized positive allocated section rows")

    totals: dict[str, int] = {}
    for _, size, addr in allocated:
        region = region_for_address(addr)
        totals[region] = totals.get(region, 0) + size

    missing = [region for region in ("FLASH", "ITCM", "DTCM", "OCRAM")
               if totals.get(region, 0) <= 0]
    if missing:
        raise SizeAFormatError(
            "missing positive Teensy memory bucket(s): " + ", ".join(missing))

    # FlexRAM hands ITCM whole banks (tools/phantasm.ld _itcm_block_count), so
    # the raw section total under-counts RAM1.
    itcm = -(-totals["ITCM"] // FLEXRAM_BANK_BYTES) * FLEXRAM_BANK_BYTES
    ram1 = itcm + totals["DTCM"]
    initialized_data = sum(size for name, size, addr in allocated
                           if name == ".data" and region_for_address(addr) == "DTCM")
    flash = totals["FLASH"] + totals["ITCM"] + initialized_data
    return {
        "flash": {"used": flash, "free": FLASH_PROGRAM_BYTES - flash},
        "ram1": {"used": ram1, "free": 0x80000 - ram1},
        "ram2": {"used": totals["OCRAM"],
                 "free": 0x80000 - totals["OCRAM"]},
    }


# `readelf -s`/`-sW` symbol row, e.g.
#    124: 20200000 248832 OBJECT  GLOBAL DEFAULT    8 _ZN6Effect8buffer_aE
# readelf prints Value in hex without a prefix. Size is decimal for small
# objects and 0x-prefixed hex for large ones.
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
        if m and m.group(1) != "0":
            # ndx 0 is the NULL section: its empty Name column shifts every \S+
            # group one field left (Type read as Name, Addr as Type). Skip it.
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


def _check_derived_component_ceiling(
    result: GateResult,
    env: str,
    region: str,
    region_spec: dict,
    cname: str,
    cmeasured: int,
    components: dict[str, int],
    derived: dict,
) -> None:
    """Enforce a stack-floor-derived component ceiling.

    FlexRAM splits RAM1's `total_banks` banks of `bank_bytes` between ITCM and
    DTCM. DTCM keeps ceil((variables + free_min_bytes) / bank_bytes) banks; the
    remaining banks less `min_headroom_bytes` cap the component. A missing input
    is a hard failure.
    """
    v = result.violations
    bank = derived["bank_bytes"]
    total_banks = derived["total_banks"]
    floor = region_spec.get("free_min_bytes")
    if floor is None:
        v.append(Violation(
            "component-floor-missing",
            f"{env}: {region.upper()} component '{cname}' has a "
            f"stack-floor-derived ceiling but the region sets no "
            f"free_min_bytes - the derivation needs the stack floor. A missing "
            f"floor is a hard failure, never a silent pass."))
        return
    variables = components.get("variables")
    if variables is None:
        v.append(Violation(
            "component-missing",
            f"{env}: {region.upper()} component 'variables' not reported by "
            f"the size output but required to derive the '{cname}' ceiling - "
            f"renamed/removed field? A missing component is a hard failure, "
            f"never a silent pass."))
        return
    dtcm_banks = -(-(variables + floor) // bank)  # ceil
    itcm_banks = total_banks - dtcm_banks
    bank_ceiling = itcm_banks * bank
    min_headroom = derived["min_headroom_bytes"]
    ceiling = bank_ceiling - min_headroom
    if cmeasured > ceiling:
        v.append(Violation(
            "component-over-derived-ceiling",
            f"{env}: {region.upper()} component '{cname}' uses {cmeasured:,} B, "
            f"over the derived {ceiling:,} B ceiling (by "
            f"{cmeasured - ceiling:,} B). The {bank_ceiling:,} B bank "
            f"allocation reserves {min_headroom:,} B of boundary headroom. "
            f"The binding constraint is the DTCM "
            f"stack floor: {variables:,} B of variables + the {floor:,} B floor "
            f"need {dtcm_banks} of {total_banks} FlexRAM banks, leaving "
            f"{itcm_banks} bank(s) x {bank:,} B for '{cname}'."))
    to_boundary = (bank - cmeasured % bank) % bank
    result.notes.append(
        f"{env}: {region.upper()} '{cname}' derived ceiling: measured "
        f"{cmeasured:,} B of {ceiling:,} B ({itcm_banks} x {bank:,} B banks; "
        f"{dtcm_banks} DTCM banks cover {variables:,} B variables + the "
        f"{floor:,} B stack floor); remaining {ceiling - cmeasured:,} B; "
        f"boundary headroom ratchet {min_headroom:,} B; "
        f"{to_boundary:,} B to the next bank boundary.")


def _evaluate_region_budgets(env, budget, sizes, result):
    v = result.violations
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
            free_label = {
                "ram1": "free for local variables (stack)",
                "ram2": "free for malloc/new (heap)",
                "flash": "free for files",
            }[region]
            v.append(Violation(
                "headroom-below-floor",
                f"{env}: {region.upper()} {free_label} "
                f"{measured.get('free', 0):,} B is below the {floor:,} B floor."))
        # Per-component ceilings: static max_bytes, stack-floor-derived, or
        # both. A configured component absent from the output is a hard failure.
        for cname, cspec in spec.get("components", {}).items():
            cmeasured = measured.get("components", {}).get(cname)
            if cmeasured is None:
                v.append(Violation(
                    "component-missing",
                    f"{env}: {region.upper()} component '{cname}' not reported "
                    f"by the size output - renamed/removed field? A missing "
                    f"component is a hard failure, never a silent pass."))
                continue
            ccap = cspec.get("max_bytes")
            if ccap is not None and cmeasured > ccap:
                v.append(Violation(
                    "component-over-budget",
                    f"{env}: {region.upper()} component '{cname}' uses "
                    f"{cmeasured:,} B, over the {ccap:,} B ceiling "
                    f"(by {cmeasured - ccap:,} B)."))
            derived = cspec.get("max_banks_from_stack_floor")
            if derived is not None:
                _check_derived_component_ceiling(
                    result, env, region, spec, cname, cmeasured,
                    measured.get("components", {}), derived)


def _evaluate_layout_symbols(env, budget, symbols, sections, result):
    v = result.violations
    for key, spec in budget.get("symbols", {}).items():
        name = spec["name"]
        # Only DEFINED rows count; readelf -s can also list a same-named UND row.
        matches = [s for s in symbols
                   if s.name == name and s.ndx != "UND" and s.size > 0]
        if not matches:
            # Fail loud: a never-matching name silently disables its invariant.
            v.append(Violation(
                "symbol-not-found",
                f"{env}: layout symbol '{key}' ({name}) not found in the ELF - "
                f"renamed/removed? A missing symbol is a hard failure, never a "
                f"silent pass."))
            continue
        # Vague-linkage copies share a VMA and size and collapse to one.
        distinct = {(s.value, s.size) for s in matches}
        if len(distinct) > 1:
            v.append(Violation(
                "symbol-duplicate",
                f"{env}: layout symbol '{key}' ({name}) has {len(distinct)} "
                f"distinct definitions in the ELF - expected exactly one."))
            continue
        sym = matches[0]
        want_region = spec.get("region")
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
                f"{hi:,} B cap (e.g. an 8 MiB host arena leak)."))


def evaluate(
    env: str,
    budget: dict,
    sizes: dict[str, RegionSizes],
    symbols: list[Symbol],
    sections: dict[str, tuple[str, int]] | None = None,
) -> GateResult:
    """Apply a target's budget to its measured sizes + symbols."""
    sections = sections or {}
    result = GateResult(env=env)
    v = result.violations

    # A budget carrying neither key would satisfy every loop vacuously; that
    # is a hard failure.
    if not budget.get("regions") and not budget.get("symbols"):
        v.append(Violation(
            "budget-empty",
            f"{env}: budget declares no 'regions' and no 'symbols' - every "
            f"check would pass vacuously (key typo or truncated entry?). An "
            f"empty budget is a hard failure, never a silent pass."))

    _evaluate_region_budgets(env, budget, sizes, result)
    _evaluate_layout_symbols(env, budget, symbols, sections, result)
    return result


# ---------------------------------------------------------------------------
# CLI (standalone / local debugging on captured toolchain output)
# ---------------------------------------------------------------------------
def read_capture(path: str | Path) -> str:
    """Read captured toolchain output, replacing undecodable bytes.

    A Windows `tee` capture interleaves cp1252 bytes; the parsers match ASCII.
    """
    return Path(path).read_text(encoding="utf-8", errors="replace")


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


#: Exit code for a PASS produced from the uncalibrated `size -A` fallback.
#: Distinct from 0 (calibrated PASS), 1 (violation) and 2 (the gate could not
#: run: unknown env, unparseable `size -A`).
EXIT_UNCALIBRATED_PASS = 3

#: Note prefixed to any verdict computed from the `size -A` fallback.
UNCALIBRATED_NOTE = (
    "ADVISORY: `size -A` fallback in use (teensy_size unavailable). Flash "
    "includes ITCM and initialized data load images, but section alignment "
    "padding may differ, so this region verdict is NOT calibrated "
    "- run with teensy_size for an authoritative ceiling decision.")


def size_format_annotation(exc: ValueError, github: bool = True) -> str:
    """Describe an unparseable size capture in the requested diagnostic mode."""
    prefix = "::error::teensy-gate: " if github else "teensy-gate: error: "
    source = ("teensy_size output" if isinstance(exc, TeensySizeFormatError)
              else "`size -A` output")
    return (f"{prefix}invalid {source} ({exc}). This is a "
            f"tooling/format error, not a size-budget violation.")


def verdict(env: str, budget: dict, sizes: dict[str, RegionSizes],
            symbols: list[Symbol], sections: dict[str, tuple[str, int]], *,
            uncalibrated: bool, github: bool) -> tuple[str, int]:
    """Evaluate one env and render it; returns (report, exit code).

    An uncalibrated PASS exits with EXIT_UNCALIBRATED_PASS.
    """
    result = evaluate(env, budget, sizes, symbols, sections)
    if uncalibrated:
        result.notes.insert(0, UNCALIBRATED_NOTE)
    if not result.passed:
        code = 1
    else:
        code = EXIT_UNCALIBRATED_PASS if uncalibrated else 0
    return render_report(result, github=github), code


def run(env: str, budgets_path: str, size_text: str, *,
        from_teensy_size: bool, syms_text: str, secs_text: str,
        github: bool) -> tuple[list[str], int]:
    """Evaluate captured tool output; return report lines and gate exit status."""
    error = "::error::teensy-gate: " if github else "teensy-gate: error: "
    warning = "::warning::" if github else "teensy-gate: warning: "
    try:
        budgets = load_budgets(budgets_path)
    except BudgetSchemaError as exc:
        return [f"{error}invalid budgets schema in {budgets_path} "
                f"({exc}). This is a budgets-file error, not a size-budget "
                "violation."], 2
    except (OSError, ValueError) as exc:
        return [f"{error}cannot load budgets from {budgets_path} "
                f"({exc}). This is a budgets-file error, not a size-budget "
                "violation."], 2
    if env not in budgets:
        return [f"{error}no budget for env '{env}' in {budgets_path}"], 2

    try:
        sizes = (parse_teensy_size(size_text) if from_teensy_size
                 else fallback_sizes_from_size_a(size_text))
    except (TeensySizeFormatError, SizeAFormatError) as exc:
        return [size_format_annotation(exc, github)], 2
    if not any(region in sizes for region in ("flash", "ram1", "ram2")):
        return [f"{error}parsed no FLASH/RAM1/RAM2 regions from "
                "the size output. This is a toolchain/format break, not a "
                "size-budget violation."], 2
    if not from_teensy_size and declares_components(budgets[env]):
        return [f"{error}env '{env}' declares per-component "
                "ceilings, which the `size -A` fallback cannot measure. "
                "No component figure was read: this is a tooling break, "
                "not a size-budget violation; re-run with teensy_size."], 2

    symbols = parse_readelf_symbols(syms_text)
    sections = parse_readelf_sections(secs_text) if secs_text else {}
    report, code = verdict(env, budgets[env], sizes, symbols, sections,
                           uncalibrated=not from_teensy_size, github=github)
    lines = []
    if not from_teensy_size:
        lines.append(f"{warning}teensy_size not found; using `size -A` fallback "
                     "(section alignment padding may differ; calibrate against teensy_size).")
    lines.append(report)
    if code == EXIT_UNCALIBRATED_PASS:
        lines.append(f"{error}PASS is UNCALIBRATED; install the "
                     "Teensy platform tools and re-run with teensy_size.")
    return lines, code


def main(argv: list[str] | None = None) -> int:
    """Run the gate. Exit: 0 PASS, 1 violation, 2 cannot-run, 3 advisory PASS."""
    p = argparse.ArgumentParser(description="Teensy 4 size/layout gate (parser).")
    p.add_argument("--env", required=True, help="budget key, e.g. holosphere")
    p.add_argument("--budgets", default="tools/teensy_budgets.json")
    # Exclusive: both measure the same regions.
    sizes_from = p.add_mutually_exclusive_group(required=True)
    sizes_from.add_argument("--teensy-size",
                            help="file with captured teensy_size stdout")
    sizes_from.add_argument("--size-a",
                            help="file with captured `size -A` output (fallback)")
    p.add_argument("--readelf-syms", required=True, help="file with `readelf -s` output")
    p.add_argument("--readelf-secs", help="file with `readelf -S` output")
    p.add_argument("--github", action="store_true", help="emit ::error:: annotations")
    args = p.parse_args(argv)
    error = "::error::teensy-gate: " if args.github else "teensy-gate: error: "

    try:
        size_text = read_capture(args.teensy_size or args.size_a)
    except OSError as exc:
        print(f"{error}cannot read the captured size output "
              f"({exc}). This is a tooling error, not a size-budget violation.",
              file=sys.stderr)
        return 2
    try:
        syms_text = read_capture(args.readelf_syms)
        secs_text = read_capture(args.readelf_secs) if args.readelf_secs else ""
    except OSError as exc:
        print(f"{error}cannot read the captured readelf output "
              f"({exc}). This is a tooling error, not a size-budget violation.",
              file=sys.stderr)
        return 2
    lines, code = run(args.env, args.budgets, size_text,
                      from_teensy_size=bool(args.teensy_size),
                      syms_text=syms_text, secs_text=secs_text, github=args.github)
    print("\n".join(lines), file=sys.stderr if code == 2 else sys.stdout)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
