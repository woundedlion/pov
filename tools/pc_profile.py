#!/usr/bin/env python3
"""Symbolize the PC-sampler dump of an HS_PROFILE_PC capture.

Reads the one `pc profile` block an HS_PROFILE_PC image logs
(tools/pc_sampler.h), resolves every sampled address bin against the exact
firmware ELF of that run with the ARM toolchain's nm, and prints each
function's share of the render samples, its memory region (ITCM, flash or
other) and size, a per-region summary, and the hottest address ranges inside
the hottest functions. Ticks serviced late (behind a higher-priority interrupt
or a masked section) are not binned and report as one deferred row.

A bin that straddles a function boundary is split between the functions by
the bytes each covers. Folded functions sharing one address are reported once,
with their alias count.

  python tools/pc_profile.py build/prof/shapeshifter_ship_pc.log

The ELF defaults to the capture's `artifact_profile_elf` provenance line.
"""

import argparse
import bisect
import os
import re
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

ITCM_END = 0x0008_0000
FLASH_START = 0x6000_0000
FLASH_END = 0x8000_0000
TEXT_TYPES = frozenset("tTwW")

HEADER_RE = re.compile(r"^pc profile: version=(\d+) (.*)$")
COUNTERS_RE = re.compile(r"^pc profile: (ticks=.*)$")
LATENCY_RE = re.compile(r"^pc profile latency:((?: \d+)+)$")
BIN_RE = re.compile(r"^pc ([0-9a-f]{8}) (\d+)$")
END_LINE = "pc profile end"
ELF_RE = re.compile(r"^profile provenance: artifact_profile_elf=(.+)$")
NM_RE = re.compile(
    r"^([0-9a-fA-F]{8})\s+(?:([0-9a-fA-F]{8})\s+)?([A-Za-z])\s+(.+)$")
HEADER_KEYS = ("hz", "bin_bytes", "capacity", "defer_counts", "counts_per_us",
               "elapsed_ms")
COUNTER_KEYS = ("ticks", "idle", "deferred", "missed", "handler", "recorded",
                "dropped", "bins", "service_avg", "service_max")


class DumpError(ValueError):
    """A capture holds no usable PC dump."""


@dataclass
class Dump:
    """One parsed `pc profile` block."""
    header: dict
    counters: dict
    latency: list
    bins: dict = field(default_factory=dict)

    @property
    def render(self):
        """Ticks taken with the render gate open."""
        c = self.counters
        return c["recorded"] + c["dropped"] + c["deferred"] + c["missed"]


@dataclass
class Function:
    """A text symbol, with any folded aliases at the same address."""
    address: int
    size: int
    name: str
    aliases: int = 0

    @property
    def end(self):
        return self.address + self.size


def region_of(address):
    """Memory region an address executes from."""
    if address < ITCM_END:
        return "ITCM"
    if FLASH_START <= address < FLASH_END:
        return "flash"
    return "other"


def parse_fields(text, keys):
    """Parse `k=v` integer fields, requiring exactly @p keys."""
    fields = dict(item.split("=", 1) for item in text.split())
    if set(fields) != set(keys):
        raise DumpError(f"pc profile fields {sorted(fields)} != {sorted(keys)}")
    return {k: int(v) for k, v in fields.items()}


def parse_dump(lines):
    """Parse the single PC dump in a capture's lines."""
    dump = None
    ended = False
    for raw in lines:
        line = raw.rstrip("\r\n")
        m = HEADER_RE.match(line)
        if m:
            if dump is not None:
                raise DumpError("more than one pc profile dump")
            if m.group(1) != "1":
                raise DumpError(f"unsupported pc profile version {m.group(1)}")
            dump = Dump(parse_fields(m.group(2), HEADER_KEYS), {}, [])
            continue
        if dump is None or ended:
            continue
        m = COUNTERS_RE.match(line)
        if m:
            dump.counters = parse_fields(m.group(1), COUNTER_KEYS)
            continue
        m = LATENCY_RE.match(line)
        if m:
            dump.latency = [int(v) for v in m.group(1).split()]
            continue
        m = BIN_RE.match(line)
        if m:
            address = int(m.group(1), 16)
            if address in dump.bins:
                raise DumpError(f"duplicate pc bin {address:08x}")
            dump.bins[address] = int(m.group(2))
            continue
        if line == END_LINE:
            ended = True
    if dump is None:
        raise DumpError("no pc profile dump (built without -D HS_PROFILE_PC?)")
    if not ended:
        raise DumpError("pc profile dump has no end line (capture cut short?)")
    if not dump.counters:
        raise DumpError("pc profile dump has no counter line")
    recorded = sum(dump.bins.values())
    if recorded != dump.counters["recorded"]:
        raise DumpError(f"pc bins sum to {recorded}, counters say "
                        f"{dump.counters['recorded']}")
    if len(dump.bins) != dump.counters["bins"]:
        raise DumpError(f"{len(dump.bins)} pc bins listed, counters say "
                        f"{dump.counters['bins']}")
    return dump


def elf_from_log(lines):
    """The profile ELF a capture's provenance names, or None."""
    for line in lines:
        m = ELF_RE.match(line.rstrip("\r\n"))
        if m:
            return m.group(1)
    return None


def parse_nm(text):
    """Text symbols from `nm -S -n -C --defined-only` output.

    Returns (address, size or None, type, name) tuples with the Thumb bit
    cleared.
    """
    symbols = []
    for line in text.splitlines():
        m = NM_RE.match(line)
        if not m or m.group(3) not in TEXT_TYPES:
            continue
        size = int(m.group(2), 16) if m.group(2) else None
        symbols.append((int(m.group(1), 16) & ~1, size, m.group(3),
                        m.group(4)))
    return symbols


def build_functions(symbols):
    """Fold symbols into address-sorted functions.

    Symbols sharing an address become one function named after its first
    global symbol (else its first symbol), sized by the largest alias. A
    symbol without a size extends to the next function's start.
    """
    by_address = {}
    for address, size, kind, name in symbols:
        by_address.setdefault(address, []).append((size, kind, name))
    addresses = sorted(by_address)
    functions = []
    for i, address in enumerate(addresses):
        group = by_address[address]
        named = [g for g in group if g[1].isupper()] or group
        sizes = [g[0] for g in group if g[0]]
        if sizes:
            size = max(sizes)
        elif i + 1 < len(addresses):
            size = addresses[i + 1] - address
        else:
            size = 0
        functions.append(Function(address, size, named[0][2], len(group) - 1))
    return functions


class Symbolizer:
    """Maps address ranges onto functions."""

    def __init__(self, functions):
        self.functions = functions
        self.starts = [f.address for f in functions]

    def overlaps(self, start, end):
        """Yield (function or None, bytes) covering [start, end).

        None stands for bytes no function covers.
        """
        i = bisect.bisect_right(self.starts, start) - 1
        cursor = start
        i = max(i, 0)
        while cursor < end and i < len(self.functions):
            f = self.functions[i]
            if f.end <= cursor:
                i += 1
                continue
            if f.address >= end:
                break
            if f.address > cursor:
                yield None, f.address - cursor
                cursor = f.address
            covered = min(f.end, end) - cursor
            yield f, covered
            cursor += covered
            i += 1
        if cursor < end:
            yield None, end - cursor


@dataclass
class Row:
    """One function's (or one unknown region's) samples."""
    name: str
    region: str
    size: int
    aliases: int
    address: int
    samples: float = 0.0
    bins: dict = field(default_factory=dict)


def attribute(dump, symbolizer):
    """Split every bin's samples across the functions it overlaps."""
    width = dump.header["bin_bytes"]
    rows = {}
    for start, count in dump.bins.items():
        for f, covered in symbolizer.overlaps(start, start + width):
            share = count * covered / width
            if f is None:
                region = region_of(start)
                key = ("?", region)
                row = rows.get(key) or rows.setdefault(
                    key, Row(f"<no symbol: {region}>", region, 0, 0, 0))
            else:
                key = f.address
                row = rows.get(key) or rows.setdefault(
                    key, Row(f.name, region_of(f.address), f.size, f.aliases,
                             f.address))
            row.samples += share
            row.bins[start] = row.bins.get(start, 0.0) + share
    return sorted(rows.values(), key=lambda r: (-r.samples, r.name))


def region_totals(rows):
    """Samples per region, in a fixed region order."""
    totals = {"ITCM": 0.0, "flash": 0.0, "other": 0.0}
    for row in rows:
        totals[row.region] += row.samples
    return totals


def pct(part, whole):
    return 100.0 * part / whole if whole else 0.0


def short(name, width):
    if width <= 0 or len(name) <= width:
        return name
    return name[:max(width - 3, 1)] + "..."


def hot_ranges(row, bin_bytes, limit):
    """The row's hottest bins as (start, end, samples), contiguous bins merged."""
    merged = []
    for start in sorted(row.bins):
        samples = row.bins[start]
        if merged and merged[-1][1] == start:
            s, _, total = merged[-1]
            merged[-1] = (s, start + bin_bytes, total + samples)
        else:
            merged.append((start, start + bin_bytes, samples))
    merged.sort(key=lambda r: -r[2])
    return merged[:limit]


def report(dump, rows, *, top, ranges, range_rows, name_width, lines=None):
    """Render the text report as a list of lines.

    Percentages are of the render samples: the binned rows plus one row for
    deferred and missed ticks (time in higher-priority handlers or masked
    sections) and one for dropped ticks. @p lines maps an address to a
    `file:line` string for the hot ranges.
    """
    c = dump.counters
    h = dump.header
    out = []
    recorded = c["recorded"]
    render = dump.render
    measured_hz = c["ticks"] * 1000.0 / h["elapsed_ms"] if h["elapsed_ms"] else 0
    out.append(f"ticks {c['ticks']} over {h['elapsed_ms']} ms "
               f"(nominal {h['hz']} Hz, measured {measured_hz:.0f} Hz); "
               f"service {c['service_avg']} cyc avg, {c['service_max']} max")
    out.append(f"render samples {render} ({pct(render, render + c['idle']):.1f}%"
               f" of all): binned {recorded}, deferred {c['deferred']}, "
               f"missed {c['missed']}, dropped {c['dropped']}, "
               f"handler-mode {c['handler']}")
    if c["dropped"]:
        out.append(f"WARNING: {c['dropped']} samples found no histogram slot; "
                   "raise HS_PROFILE_PC_BINS")
    totals = region_totals(rows)
    late = c["deferred"] + c["missed"]
    out.append("regions (% of render): " + "  ".join(
        f"{region} {pct(samples, render):.1f}%"
        for region, samples in totals.items())
        + f"  deferred {pct(late, render):.1f}%"
        + f"  dropped {pct(c['dropped'], render):.1f}%")
    out.append("regions (% of binned): " + "  ".join(
        f"{region} {pct(samples, recorded):.1f}%"
        for region, samples in totals.items()))
    table = list(rows)
    if late:
        table.append(Row("<deferred: higher-priority IRQ or IRQs masked>",
                         "-", 0, 0, 0, late))
    if c["dropped"]:
        table.append(Row("<dropped: histogram full>", "-", 0, 0, 0,
                         c["dropped"]))
    table.sort(key=lambda r: (-r.samples, r.name))
    out.append("")
    out.append(f"{'pct':>6} {'samples':>8} {'region':<6} {'size':>6}  function")
    for row in table[:top]:
        alias = f" [+{row.aliases} folded]" if row.aliases else ""
        out.append(f"{pct(row.samples, render):5.1f}% {row.samples:8.0f} "
                   f"{row.region:<6} {row.size:6d}  "
                   f"{short(row.name, name_width)}{alias}")
    if ranges:
        out.append("")
        out.append("hot ranges (function offset, samples, % of render):")
        for row in rows[:ranges]:
            if row.name.startswith("<no symbol"):
                continue
            out.append(f"  {short(row.name, name_width)}")
            for start, end, samples in hot_ranges(row, h["bin_bytes"],
                                                  range_rows):
                where = (lines or {}).get(start, "")
                out.append(f"    +0x{start - row.address:04x}..+0x"
                           f"{end - row.address - 1:04x} {samples:8.0f} "
                           f"{pct(samples, render):5.1f}%"
                           + (f"  {where}" if where else ""))
    return out


def toolchain_tool(name, override):
    """Path of an ARM binutils tool, honoring an explicit override."""
    if override:
        return override
    core = os.environ.get("PLATFORMIO_CORE_DIR",
                          str(Path.home() / ".platformio"))
    base = Path(core) / "packages" / "toolchain-gccarmnoneeabi-teensy" / "bin"
    exe = base / f"arm-none-eabi-{name}.exe"
    return str(exe if exe.exists() else base / f"arm-none-eabi-{name}")


def run_nm(nm, elf):
    return subprocess.run([nm, "-S", "-n", "-C", "--defined-only", elf],
                          check=True, capture_output=True, text=True).stdout


def run_addr2line(addr2line, elf, addresses):
    """Map addresses to `file:line`, omitting unresolved ones."""
    if not addresses:
        return {}
    try:
        result = subprocess.run(
            [addr2line, "-e", elf, *(f"0x{a:x}" for a in addresses)],
            check=True, capture_output=True, text=True).stdout.splitlines()
    except (OSError, subprocess.CalledProcessError):
        return {}
    return {a: loc for a, loc in zip(addresses, result)
            if loc and not loc.startswith("??")}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("log", help="capture log holding the pc profile dump")
    ap.add_argument("--elf", help="firmware ELF (default: the log's "
                    "artifact_profile_elf provenance line)")
    ap.add_argument("--nm-file", help="read `nm -S -n -C --defined-only` "
                    "output from this file instead of running nm")
    ap.add_argument("--nm", help="nm executable (default: PlatformIO's ARM "
                    "toolchain, or $HS_ARM_NM)",
                    default=os.environ.get("HS_ARM_NM"))
    ap.add_argument("--addr2line", default=os.environ.get("HS_ARM_ADDR2LINE"),
                    help="addr2line executable (default: PlatformIO's)")
    ap.add_argument("--top", type=int, default=40, help="functions listed")
    ap.add_argument("--ranges", type=int, default=8,
                    help="hottest functions given address ranges")
    ap.add_argument("--range-rows", type=int, default=4,
                    help="ranges per function")
    ap.add_argument("--name-width", type=int, default=150,
                    help="truncate names past this width (0: never)")
    args = ap.parse_args(argv)

    lines = Path(args.log).read_text(encoding="utf-8",
                                     errors="replace").splitlines()
    try:
        dump = parse_dump(lines)
    except DumpError as exc:
        print(f"pc_profile: {args.log}: {exc}", file=sys.stderr)
        return 1
    elf = args.elf or elf_from_log(lines)
    if args.nm_file:
        nm_text = Path(args.nm_file).read_text(encoding="utf-8")
    else:
        if not elf or not Path(elf).is_file():
            print(f"pc_profile: no firmware ELF ({elf or 'none in the log'}); "
                  "pass --elf", file=sys.stderr)
            return 1
        nm_text = run_nm(toolchain_tool("nm", args.nm), elf)
    functions = build_functions(parse_nm(nm_text))
    rows = attribute(dump, Symbolizer(functions))
    locations = {}
    if elf and Path(elf).is_file() and args.ranges:
        starts = sorted({start for row in rows[:args.ranges]
                         for start, _, _ in hot_ranges(
                             row, dump.header["bin_bytes"], args.range_rows)})
        locations = run_addr2line(toolchain_tool("addr2line", args.addr2line),
                                  elf, starts)
    print(f"pc profile: {args.log}")
    print(f"elf: {elf}")
    print("\n".join(report(dump, rows, top=args.top, ranges=args.ranges,
                           range_rows=args.range_rows,
                           name_width=args.name_width, lines=locations)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
