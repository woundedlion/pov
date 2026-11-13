#!/usr/bin/env python3
"""Instruction-cache layout of Phantasm's flash-resident hot code; stdlib only.

The Teensy 4.0 instruction cache is 32 KiB, 2-way, with 32-byte lines: 512 sets
of two lines, one 16 KiB way. Code outside ITCM runs from flash through it, so
an effect slows when its hot flash code maps more lines to one set than the
cache has ways.

`plan` reads a probe link (ELF, linker map and object section tables), derives
each effect's hot set and writes the linker-script fragment tools/phantasm.ld
INCLUDEs. `check` measures the per-set occupancy of every effect's hot set in a
linked image and fails when any effect exceeds `MAX_LINES_PER_SET`.

An effect's roots are the functions whose mangled names name its class or its
`<Class>Spec`, plus the POVSegmented driver. Its hot set is the transitive
closure over direct branches and code-address literals, through ITCM, startup,
hot, cold and toolchain-library code and linker veneers, and through the flash
objects a function references that are at least half function pointers,
restricted to flash `.text.hot*` sections, toolchain-library sections and
cold sections that hot, ITCM or startup code calls. Explicit linker startup
ranges retain roots and traversal without placing the startup code. Members
rank by call depth from the roots, then by depth-first call order; the gated
set takes them in rank order, skipping any that would exceed `GATED_BYTES`.

Run:  python tools/hot_layout.py check --elf firmware.elf --map firmware.map
"""

from __future__ import annotations

import argparse
import bisect
import heapq
import os
import re
import shutil
import struct
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

from cpp_roster import (PHANTASM_PLAYLIST_SOURCE, macro_body, named_rows,
                        playlist_shows)
from teensy_gate import parse_readelf_symbols, sibling_tool

#: Instruction-cache line size in bytes.
CACHE_LINE_BYTES = 32
#: Instruction-cache sets.
CACHE_SETS = 512
#: Instruction-cache ways.
CACHE_WAYS = 2
#: Bytes one cache way spans.
CACHE_WAY_BYTES = CACHE_LINE_BYTES * CACHE_SETS
#: Most lines of one effect's gated hot set that may map to one set.
MAX_LINES_PER_SET = CACHE_WAYS
#: Hot bytes per effect the gate covers, highest rank first. It stays well below
#: the cache size so the planner has slack to fit shared blocks between effects.
GATED_BYTES = 20 * 1024

#: Most planning passes before the plan with fewest lines over the bound wins.
PLAN_ROUNDS = 32

#: File name of the generated fragment that tools/phantasm.ld INCLUDEs.
FRAGMENT_NAME = "phantasm_hot_layout.ld"
#: Linker symbol tools/phantasm.ld defines where the fragment starts.
START_SYMBOL = "hs_hot_layout_start"
#: Linker ranges whose startup code retains its graph role after leaving ITCM.
STARTUP_RANGE_SYMBOLS = (
    ("hs_hot_layout_startup_ctors_start", "hs_hot_layout_startup_ctors_end"),
    ("hs_hot_layout_startup_raster_start", "hs_hot_layout_startup_raster_end"),
)
#: Toolchain archives whose code the layout may place: arithmetic, libm, and
#: the libc memory primitives (`libc_a-mem*.o`).
LIBRARY_ARCHIVES = ("libm.a", "libgcc.a")
LIBRARY_MEMBERS = re.compile(r"^(?:libc_a-mem\w*\.o)$")
LIBC_ARCHIVE = "libc_nano.a"
#: Mangled class names of the code every effect runs.
DRIVER_CLASSES = ("POVSegmented",)

ITCM_END = 0x20000000
FLASH_START = 0x60000000
FLASH_OUTPUT = ".text.code"
ITCM_OUTPUT = ".text.itcm"

# Function kinds.
ITCM = "itcm"
HOT = "hot"
COLD = "cold"
LIBRARY = "library"
VENEER = "veneer"
SETUP = "setup"
STARTUP = "startup"
#: Kinds the closure walks through.
TRAVERSED = frozenset({ITCM, HOT, COLD, LIBRARY, VENEER, STARTUP})
#: Kinds the layout places.
PLACED = frozenset({HOT, LIBRARY})
#: Kinds that root an effect when its name names the effect.
ROOTED = frozenset({ITCM, HOT, COLD, STARTUP})


class LayoutError(Exception):
    """Malformed input or a toolchain failure."""


@dataclass(frozen=True, order=True)
class Origin:
    """Input file of a section: object path, or archive path plus member."""

    path: str
    member: str = ""

    @property
    def archive(self) -> str:
        """Archive file name, or empty for a plain object."""
        return _basename(self.path) if self.member else ""


@dataclass(frozen=True)
class Section:
    """One input section of the linked image."""

    output: str
    name: str
    origin: Origin
    address: int
    size: int

    @property
    def key(self) -> tuple[Origin, str]:
        """Identity of the section across links."""
        return (self.origin, self.name)


@dataclass(frozen=True)
class Function:
    """One FUNC symbol."""

    name: str
    address: int
    size: int


@dataclass
class EffectHotSet:
    """Placed members of one effect's closure, ranked."""

    effect: str
    ranked: list[tuple[int, tuple[Origin, str]]] = field(default_factory=list)
    gated: list[tuple[Origin, str]] = field(default_factory=list)
    overflow: list[tuple[Origin, str]] = field(default_factory=list)
    order: dict[tuple[Origin, str], int] = field(default_factory=dict)


def _basename(path: str) -> str:
    return re.split(r"[\\/]", path)[-1]


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------
_OUTPUT_HEADER = re.compile(r"^(\.\S+)(?:\s+0x([0-9a-f]+)\s+0x([0-9a-f]+))?")
_INPUT_HEADER = re.compile(
    r"^ (\.\S+)(?:\s+0x([0-9a-f]+)\s+0x([0-9a-f]+)\s+(.+))?$")
_INPUT_TAIL = re.compile(r"^\s+0x([0-9a-f]+)\s+0x([0-9a-f]+)\s+(.+)$")
_SYMBOL_ASSIGN = re.compile(r"^\s+0x([0-9a-f]+)\s+(\w+) = ")
_ARCHIVE_MEMBER = re.compile(r"^(.*\.a)\((.+)\)$")


def parse_origin(text: str) -> Origin:
    """Origin from a linker-map file column such as `lib.a(member.o)`."""
    text = text.strip()
    member = _ARCHIVE_MEMBER.match(text)
    if member:
        return Origin(member.group(1), member.group(2))
    return Origin(text)


def parse_map(text: str) -> tuple[list[Section], dict[str, int]]:
    """Non-empty code input sections and assigned symbols of a GNU ld map."""
    marker = text.find("Linker script and memory map")
    if marker < 0:
        raise LayoutError("linker map has no 'Linker script and memory map'")
    sections: list[Section] = []
    symbols: dict[str, int] = {}
    output = ""
    pending = ""
    for line in text[marker:].splitlines():
        if pending:
            tail = _INPUT_TAIL.match(line)
            if tail:
                _add_section(sections, output, pending, tail.group(1),
                             tail.group(2), tail.group(3))
            pending = ""
            continue
        if line.startswith(" ."):
            head = _INPUT_HEADER.match(line)
            if head and head.group(2) is None:
                pending = head.group(1)
            elif head:
                _add_section(sections, output, head.group(1), head.group(2),
                             head.group(3), head.group(4))
            continue
        if line.startswith("."):
            head = _OUTPUT_HEADER.match(line)
            if head:
                output = head.group(1)
            continue
        assign = _SYMBOL_ASSIGN.match(line)
        if assign:
            symbols[assign.group(2)] = int(assign.group(1), 16)
    return sections, symbols


def _add_section(sections, output, name, address, size, origin):
    size = int(size, 16)
    if size and output in (FLASH_OUTPUT, ITCM_OUTPUT):
        sections.append(Section(output, name, parse_origin(origin),
                                int(address, 16), size))


_OBJDUMP_MEMBER = re.compile(r"^(\S+):\s+file format ")
_OBJDUMP_SECTION = re.compile(
    r"^\s*\d+\s+(\S+)\s+[0-9a-f]+\s+[0-9a-f]+\s+[0-9a-f]+\s+[0-9a-f]+\s+2\*\*(\d+)\s")


def parse_section_alignments(text: str, path: str) -> dict[tuple[Origin, str], int]:
    """Section alignments from `objdump -h -w` of one object or archive.

    readelf truncates long section names, so the alignments come from objdump.
    """
    alignments: dict[tuple[Origin, str], int] = {}
    archive = text.lstrip().startswith("In archive ")
    origin = Origin(path)
    for line in text.splitlines():
        header = _OBJDUMP_MEMBER.match(line)
        if header:
            origin = Origin(path, header.group(1)) if archive else Origin(path)
            continue
        row = _OBJDUMP_SECTION.match(line)
        if row:
            alignments[(origin, row.group(1))] = 1 << int(row.group(2))
    return alignments


def parse_functions(text: str) -> list[Function]:
    """FUNC symbols from `readelf -sW`, one per address (largest size)."""
    by_address: dict[int, Function] = {}
    for symbol in parse_readelf_symbols(text):
        if symbol.type != "FUNC":
            continue
        address = symbol.value & ~1
        known = by_address.get(address)
        if known is None or symbol.size > known.size or (
                symbol.size == known.size and symbol.name < known.name):
            by_address[address] = Function(symbol.name, address, symbol.size)
    return sorted(by_address.values(), key=lambda f: f.address)


_INSN = re.compile(r"^\s*([0-9a-f]+):\s+(\S+)\s*(.*)$")
_BRANCH_TARGET = re.compile(r"^(?:[^,]*,\s*)?([0-9a-f]+) <")
_WORD = re.compile(r"^0x([0-9a-f]+)")
_MOVE = re.compile(r"^(\w+), #-?\d+\s+@ 0x([0-9a-f]+)")


def parse_references(text: str, function_starts: set[int] | None = None
                     ) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    """(site, target) code and data references in `objdump -d` output.

    Branch targets and odd literal-pool words or movw/movt pairs are code
    references; even ones are data references. Targets are raw values,
    resolved against functions and data objects by the caller.
    """
    references: list[tuple[int, int]] = []
    data: list[tuple[int, int]] = []
    low_halves: dict[str, int] = {}
    for line in text.splitlines():
        insn = _INSN.match(line)
        if not insn:
            if line.endswith(">:"):
                low_halves.clear()
            continue
        site = int(insn.group(1), 16)
        mnemonic, operands = insn.group(2), insn.group(3)
        if mnemonic == ".word":
            word = _WORD.match(operands)
            if word:
                value = int(word.group(1), 16)
                if value & 1:
                    if function_starts is None or (value & ~1) in function_starts:
                        references.append((site, value & ~1))
                else:
                    data.append((site, value))
        elif mnemonic.startswith(("b", "cb")):
            target = _BRANCH_TARGET.match(operands)
            if target:
                references.append((site, int(target.group(1), 16)))
        elif mnemonic in ("movw", "movt"):
            move = _MOVE.match(operands)
            if move:
                register, value = move.group(1), int(move.group(2), 16)
                if mnemonic == "movw":
                    low_halves[register] = value
                elif register in low_halves:
                    full = (value << 16) | low_halves.pop(register)
                    if full & 1:
                        if function_starts is None or (full & ~1) in function_starts:
                            references.append((site, full & ~1))
                    else:
                        data.append((site, full))
    return references, data


def parse_objects(text: str) -> list[tuple[int, int]]:
    """(address, size) of every sized flash OBJECT symbol in `readelf -sW`."""
    return sorted({(symbol.value, symbol.size) for symbol in parse_readelf_symbols(text)
                   if symbol.type == "OBJECT" and symbol.size >= 4
                   and symbol.value >= FLASH_START})


def read_words(elf: bytes, ranges: list[tuple[int, int]]) -> dict[int, list[int]]:
    """Aligned 32-bit words of each (address, size) range of a little-endian ELF32."""
    shoff, = struct.unpack_from("<I", elf, 0x20)
    shentsize, shnum = struct.unpack_from("<HH", elf, 0x2E)
    loaded = []
    for index in range(shnum):
        _, kind, _, address, offset, size = struct.unpack_from(
            "<IIIIII", elf, shoff + index * shentsize)
        if kind == 1 and address:
            loaded.append((address, offset, size))
    words: dict[int, list[int]] = {}
    for start, size in ranges:
        for address, offset, length in loaded:
            if address <= start and start + size <= address + length:
                first = (start + 3) & ~3
                count = (start + size - first) // 4
                base = offset + first - address
                words[start] = list(struct.unpack_from(f"<{count}I", elf, base))
                break
    return words


_IDENTIFIER = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")


def mangled_identifiers(name: str) -> set[str]:
    """Length-prefixed source identifiers embedded in an Itanium-mangled name."""
    found: set[str] = set()
    for number in re.finditer(r"(?<![0-9])([1-9][0-9]*)", name):
        length = int(number.group(1))
        word = name[number.end():number.end() + length]
        if len(word) == length and _IDENTIFIER.fullmatch(word):
            found.add(word)
    return found


def parse_effect_classes(text: str) -> list[str]:
    """Class names of the X(Class, ...) rows of an effect list macro."""
    return named_rows(text)


def parse_effect_list(text: str, macro: str) -> list[str]:
    """Class names in the body of `#define <macro>(X)`."""
    body = macro_body(text, macro)
    if body is None:
        raise LayoutError(f"no #define {macro}(X)")
    return parse_effect_classes(body)


# ---------------------------------------------------------------------------
# Program model
# ---------------------------------------------------------------------------
class Program:
    """Functions, sections and the code-reference graph of one linked image."""

    def __init__(self, functions: list[Function], sections: list[Section],
                 references: list[tuple[int, int]],
                 data: list[tuple[int, int]] = (),
                 objects: dict[int, tuple[int, list[int]]] | None = None,
                 startup_ranges: tuple[tuple[int, int], ...] = ()):
        self.functions = sorted(functions, key=lambda f: f.address)
        self.starts = [f.address for f in self.functions]
        self.by_address = {f.address: f for f in functions}
        self.sections = sorted(sections, key=lambda s: s.address)
        self.section_starts = [s.address for s in self.sections]
        self.by_key = {s.key: s for s in self.sections}
        self.section_of = {f.address: self._containing_section(f.address)
                           for f in self.functions}
        self.startup_ranges = startup_ranges
        self.kind = {f.address: self._classify(f) for f in self.functions}
        objects = objects or {}
        object_starts = sorted(objects)
        targets = list(references)
        for site, address in data:
            index = bisect.bisect_right(object_starts, address) - 1
            if index < 0:
                continue
            start = object_starts[index]
            size, words = objects[start]
            pointers = [word & ~1 for word in words
                        if word & 1 and (word & ~1) in self.by_address]
            if address < start + size and 2 * len(pointers) >= len(words):
                targets.extend((site, pointer) for pointer in pointers)
        self.edges: dict[int, set[int]] = {}
        self.calls: dict[int, list[int]] = {}
        for site, target in sorted(targets, key=lambda reference: reference[0]):
            caller = self.function_at(site)
            callee = self.function_at(target)
            if caller is not None and callee is not None and caller != callee:
                callees = self.edges.setdefault(caller, set())
                if callee not in callees:
                    callees.add(callee)
                    self.calls.setdefault(caller, []).append(callee)
        self.hot_called = self._hot_called()
        self.identifiers = {f.address: mangled_identifiers(f.name)
                            for f in self.functions}

    def function_at(self, address: int) -> int | None:
        """Start address of the function containing `address`."""
        index = bisect.bisect_right(self.starts, address) - 1
        if index < 0:
            return None
        start = self.starts[index]
        if address < start + max(self.by_address[start].size, 2):
            return start
        return None

    def _containing_section(self, address: int) -> Section | None:
        index = bisect.bisect_right(self.section_starts, address) - 1
        if index < 0:
            return None
        section = self.sections[index]
        return section if address < section.address + section.size else None

    def _classify(self, function: Function) -> str:
        if function.name.endswith("_veneer"):
            return VENEER
        if function.address < ITCM_END:
            return ITCM
        section = self.section_of[function.address]
        if section is None or section.output != FLASH_OUTPUT:
            return SETUP
        if any(start <= function.address < end for start, end in self.startup_ranges):
            return STARTUP
        if section.name.startswith(".text.hot"):
            return HOT
        origin = section.origin
        if origin.archive in LIBRARY_ARCHIVES or (
                origin.archive == LIBC_ARCHIVE and LIBRARY_MEMBERS.match(origin.member)):
            return LIBRARY
        if section.name.startswith(".text.unlikely"):
            return COLD
        return SETUP

    def _hot_called(self) -> set[int]:
        """Cold functions called by hot, ITCM or startup code, including through veneers."""
        found: set[int] = set()
        for caller, callees in self.edges.items():
            if self.kind[caller] not in (HOT, ITCM, STARTUP):
                continue
            for callee in callees:
                if self.kind[callee] == VENEER:
                    found.update(c for c in self.edges.get(callee, ()) if self.kind[c] == COLD)
                elif self.kind[callee] == COLD:
                    found.add(callee)
        return found

    def placeable(self, function: int) -> bool:
        """True for hot and library code and for cold code that hot code calls."""
        return self.kind[function] in PLACED or function in self.hot_called

    def roots(self, classes: tuple[str, ...]) -> list[int]:
        """Functions whose mangled names name any of `classes`."""
        names = set(classes)
        return [f.address for f in self.functions
                if self.kind[f.address] in ROOTED
                and self.identifiers[f.address] & names]

    def depths(self, roots: list[int]) -> dict[int, int]:
        """Call depth of every function reachable from `roots`; veneers are free."""
        depth: dict[int, int] = {}
        queue = [(0, root) for root in roots]
        heapq.heapify(queue)
        while queue:
            level, function = heapq.heappop(queue)
            if function in depth:
                continue
            depth[function] = level
            for callee in self.edges.get(function, ()):
                if callee not in depth and self.kind[callee] in TRAVERSED:
                    step = 0 if self.kind[callee] == VENEER else 1
                    heapq.heappush(queue, (level + step, callee))
        return depth

    def call_order(self, roots: list[int]) -> dict[int, int]:
        """Depth-first preorder from `roots`, callees in call-site order.

        Roots that no reached function calls go first, by name.
        """
        reached = set(self.depths(roots))
        called = {callee for caller in reached for callee in self.edges.get(caller, ())}
        index: dict[int, int] = {}
        for root in sorted(roots, key=lambda r: (r in called, self.by_address[r].name)):
            if root in index:
                continue
            index[root] = len(index)
            stack = [iter(self.calls.get(root, ()))]
            while stack:
                callee = next(stack[-1], None)
                if callee is None:
                    stack.pop()
                elif callee not in index and self.kind[callee] in TRAVERSED:
                    index[callee] = len(index)
                    stack.append(iter(self.calls.get(callee, ())))
        return index

    def hot_set(self, effect: str, driver: list[int]) -> EffectHotSet:
        """Placed members reachable from `effect`'s roots and `driver`."""
        result = EffectHotSet(effect)
        roots = self.roots((effect, effect + "Spec"))
        if not roots:
            return result
        depth = self.depths(roots + driver)
        visit = self.call_order(roots + driver)
        rank: dict[tuple[Origin, str], int] = {}
        for function, level in depth.items():
            section = self.section_of[function]
            if not self.placeable(function) or section is None:
                continue
            key = section.key
            rank[key] = min(level, rank.get(key, level))
            result.order[key] = min(visit[function], result.order.get(key, visit[function]))
        result.ranked = sorted(((level, key) for key, level in rank.items()),
                               key=lambda item: (item[0], result.order[item[1]]))
        total = 0
        for _, key in result.ranked:
            size = self.by_key[key].size
            if total + size <= GATED_BYTES:
                result.gated.append(key)
                total += size
            else:
                result.overflow.append(key)
        return result

    def hot_sets(self, effects: list[str]) -> dict[str, EffectHotSet]:
        """Hot sets of every effect in `effects` with roots in the image."""
        driver = self.roots(DRIVER_CLASSES)
        sets = {}
        for effect in effects:
            hot = self.hot_set(effect, driver)
            if hot.ranked:
                sets[effect] = hot
        return sets


# ---------------------------------------------------------------------------
# Occupancy and the gate
# ---------------------------------------------------------------------------
def lines_of(address: int, size: int) -> range:
    """Cache-line indices a byte range touches."""
    return range(address // CACHE_LINE_BYTES,
                 (address + max(size, 1) - 1) // CACHE_LINE_BYTES + 1)


def occupancy(sections: list[Section]) -> dict[int, list[int]]:
    """Distinct lines per cache set, keyed by set index."""
    lines: set[int] = set()
    for section in sections:
        lines.update(lines_of(section.address, section.size))
    per_set: dict[int, list[int]] = {}
    for line in sorted(lines):
        per_set.setdefault(line % CACHE_SETS, []).append(line)
    return per_set


@dataclass
class EffectReport:
    """Gate result for one effect."""

    effect: str
    gated_bytes: int
    overflow_bytes: int
    worst: int
    full_worst: int
    over_sets: int
    offenders: list[str]

    @property
    def passed(self) -> bool:
        """True when no set holds more than `MAX_LINES_PER_SET` gated lines."""
        return self.worst <= MAX_LINES_PER_SET


def evaluate(program: Program, hot_sets: dict[str, EffectHotSet]) -> list[EffectReport]:
    """Per-effect occupancy of the gated and full hot sets in `program`."""
    reports = []
    for effect, hot in hot_sets.items():
        gated = [program.by_key[key] for key in hot.gated]
        full = gated + [program.by_key[key] for key in hot.overflow]
        per_set = occupancy(gated)
        worst = max((len(lines) for lines in per_set.values()), default=0)
        over = [s for s, lines in per_set.items() if len(lines) > MAX_LINES_PER_SET]
        offenders: list[str] = []
        if over:
            worst_set = min(s for s, lines in per_set.items() if len(lines) == worst)
            touched = set(per_set[worst_set])
            for section in gated:
                if touched & set(lines_of(section.address, section.size)):
                    offenders.append(_function_label(section))
        reports.append(EffectReport(
            effect, sum(s.size for s in gated),
            sum(program.by_key[k].size for k in hot.overflow), worst,
            max((len(v) for v in occupancy(full).values()), default=0),
            len(over), offenders))
    return reports


def _function_label(section: Section) -> str:
    name = section.name
    for prefix in (".text.hot.", ".text."):
        if name.startswith(prefix):
            name = name[len(prefix):]
            break
    if section.origin.member:
        return f"{section.origin.archive}:{section.origin.member}({name})"
    return name


def format_reports(env: str, reports: list[EffectReport], missing: list[str]) -> list[str]:
    """Gate lines: one per effect, offenders under a failing one."""
    failed = [r for r in reports if not r.passed]
    verdict = "FAIL" if failed or missing else "PASS"
    lines = [f"[hot-layout] {env}: {verdict} - {len(reports)} effects, at most "
             f"{MAX_LINES_PER_SET} lines per instruction-cache set over each "
             f"effect's first {GATED_BYTES // 1024} KiB of hot code"]
    for effect in missing:
        lines.append(f"  FAIL {effect}: no code in the image names this effect")
    for report in sorted(reports, key=lambda r: (r.passed, -r.worst, r.effect)):
        status = "ok  " if report.passed else "FAIL"
        lines.append(
            f"  {status} {report.effect:<24} worst {report.worst} lines/set "
            f"({report.over_sets} sets over)  gated {report.gated_bytes:>6,} B  "
            f"overflow {report.overflow_bytes:>6,} B  full-set worst {report.full_worst}")
        for offender in report.offenders:
            lines.append(f"         {offender}")
    return lines


# ---------------------------------------------------------------------------
# Planner
# ---------------------------------------------------------------------------
def _align(address: int, alignment: int) -> int:
    return (address + alignment - 1) // alignment * alignment


class _Planner:
    """One placement pass: hot units at cache phases where their gated code fits."""

    def __init__(self, program: Program, hot_sets: dict[str, EffectHotSet],
                 alignments: dict[tuple[Origin, str], int], start: int,
                 promoted: dict[tuple[Origin, str], int]):
        self.program = program
        self.alignments = alignments
        self.hot_sets = hot_sets
        self.placed: dict[tuple[Origin, str], tuple[int, int]] = {}
        self.address = start
        self.order: list[tuple[Origin, str]] = []
        members: dict[tuple[Origin, str], set[str]] = {}
        self.gates: dict[tuple[Origin, str], frozenset[str]] = {}
        for effect, hot in hot_sets.items():
            for _, key in hot.ranked:
                members.setdefault(key, set()).add(effect)
            for key in hot.gated:
                self.gates[key] = self.gates.get(key, frozenset()) | {effect}
        grouped: dict[frozenset[str], list[tuple[Origin, str]]] = {}
        for key, effects in members.items():
            grouped.setdefault(frozenset(effects), []).append(key)
        self.units: list[tuple[frozenset[str], int, list[tuple[Origin, str]]]] = []
        for effects, keys in grouped.items():
            lead = hot_sets[min(effects)].order
            keys.sort(key=lambda k: (lead[k], k[1], k[0]))
            chunk: list[tuple[Origin, str]] = []
            size = 0
            for key in keys:
                section_size = program.by_key[key].size
                if chunk and size + section_size > CACHE_WAY_BYTES:
                    self.units.append((effects, len(self.units), chunk))
                    chunk, size = [], 0
                chunk.append(key)
                size += section_size
            self.units.append((effects, len(self.units), chunk))
        self.units.sort(key=lambda unit: (-promoted.get(unit[2][0], 0), -len(unit[0]),
                                          -self._gated_size(unit[2]), sorted(unit[0]),
                                          unit[1]))
        self.unit_of = {key: unit[2][0] for unit in self.units for key in unit[2]}
        self.filler = sorted(
            (s.key for s in program.sections
             if s.output == FLASH_OUTPUT and s.key not in members
             and s.name.startswith((".text.hot", ".text.unlikely"))),
            key=lambda k: (-program.by_key[k].size, k[1], k[0]))
        self.lines = {effect: [0] * CACHE_SETS for effect in hot_sets}
        self.last = {effect: -1 for effect in hot_sets}

    def _span(self, key: tuple[Origin, str], at: int) -> tuple[int, int]:
        begin = _align(at, self.alignments[key])
        return begin, begin + self.program.by_key[key].size

    def _gated_size(self, keys: list[tuple[Origin, str]]) -> int:
        return sum(self.program.by_key[key].size for key in keys if key in self.gates)

    def _gated_lines(self, keys, at: int) -> tuple[list[tuple[int, frozenset[str]]], int]:
        """(line, gating effects) of every gated section of `keys` placed from `at`."""
        lines: list[tuple[int, frozenset[str]]] = []
        for key in keys:
            begin, at = self._span(key, at)
            effects = self.gates.get(key)
            if effects:
                lines.extend((line, effects) for line in lines_of(begin, at - begin))
        return lines, at

    def _fits(self, keys, at: int) -> bool:
        lines, _ = self._gated_lines(keys, at)
        added: dict[tuple[str, int], int] = {}
        counted: dict[str, int] = {}
        for line, effects in lines:
            for effect in effects:
                if line == self.last[effect] or counted.get(effect) == line:
                    continue
                counted[effect] = line
                index = line % CACHE_SETS
                total = added.get((effect, index), 0) + 1
                added[(effect, index)] = total
                if self.lines[effect][index] + total > MAX_LINES_PER_SET:
                    return False
        return True

    def _shift(self, keys) -> int | None:
        """Fewest whole lines to advance before the gated code of `keys` fits."""
        lines, _ = self._gated_lines(keys, self.address)
        if not lines:
            return 0
        if self._fits(keys, self.address):
            return 0
        room = [MAX_LINES_PER_SET] * CACHE_SETS
        for effect in {e for _, effects in lines for e in effects}:
            counts = self.lines[effect]
            for index in range(CACHE_SETS):
                room[index] = min(room[index], MAX_LINES_PER_SET - counts[index])
        first = lines[0][0]
        need: dict[int, int] = {}
        for line in sorted({line for line, _ in lines}):
            offset = (line - first) % CACHE_SETS
            need[offset] = need.get(offset, 0) + 1
        start = first % CACHE_SETS
        for shift in range(1, CACHE_SETS):
            base = start + shift
            if all(count <= room[(base + offset) % CACHE_SETS]
                   for offset, count in need.items()):
                if self._fits(keys, self.address + shift * CACHE_LINE_BYTES):
                    return shift
        return None

    def _place(self, keys) -> None:
        at = self.address
        for key in keys:
            self.placed[key] = self._span(key, at)
            at = self.placed[key][1]
        lines, end = self._gated_lines(keys, self.address)
        for line, effects in lines:
            for effect in effects:
                if line != self.last[effect]:
                    self.lines[effect][line % CACHE_SETS] += 1
                    self.last[effect] = line
        self.order.extend(keys)
        self.address = end

    def violations(self) -> tuple[int, set[tuple[Origin, str]]]:
        """Lines over the bound across effects, and the units holding them."""
        over = 0
        units: set[tuple[Origin, str]] = set()
        for effect, counts in self.lines.items():
            full = {index for index, count in enumerate(counts)
                    if count > MAX_LINES_PER_SET}
            over += sum(counts[index] - MAX_LINES_PER_SET for index in full)
            for key in self.hot_sets[effect].gated if full else ():
                begin, end = self.placed[key]
                if any(line % CACHE_SETS in full for line in lines_of(begin, end - begin)):
                    units.add(self.unit_of[key])
        return over, units

    def _gap_unit(self, users: frozenset[str], gap: int) -> int | None:
        """Largest other unit that fits here within `gap` bytes, sharers first."""
        best = None
        for index in range(1, len(self.units)):
            others, _, keys = self.units[index]
            end = self.address
            for key in keys:
                end = self._span(key, end)[1]
            if end - self.address > gap or not self._fits(keys, self.address):
                continue
            rank = (-len(others & users), -(end - self.address), index)
            if best is None or rank < best[0]:
                best = (rank, index)
        return None if best is None else best[1]

    def _gap_filler(self, gap: int) -> tuple[Origin, str] | None:
        for key in self.filler:
            if self._span(key, self.address)[1] - self.address <= gap:
                return key
        return None

    def run(self) -> list[tuple[Origin, str]]:
        while self.units:
            users, serial, keys = self.units[0]
            shift = self._shift(keys)
            if shift is None and len(keys) > 1:
                half = len(keys) // 2
                self.units[0:1] = [(users, serial, keys[:half]), (users, serial, keys[half:])]
                continue
            if not shift:
                self._place(keys)
                self.units.pop(0)
                continue
            gap = shift * CACHE_LINE_BYTES
            unit = self._gap_unit(users, gap)
            if unit is not None:
                self._place(self.units.pop(unit)[2])
                continue
            filler = self._gap_filler(gap)
            if filler is None:
                self._place(keys)
                self.units.pop(0)
                continue
            self.filler.remove(filler)
            self._place([filler])
        return self.order


def plan(program: Program, hot_sets: dict[str, EffectHotSet],
         alignments: dict[tuple[Origin, str], int], start: int) -> list[tuple[Origin, str]]:
    """Order of the gated sections, and any filler between them, from `start`.

    Every hot-set member groups into units by the set of effects whose hot sets
    hold it, in the first such effect's depth-first call order, split at one
    cache way, so each call subtree stays contiguous. Units go most-shared first,
    then largest gated code first: a unit whose gated code fits at the current
    address is placed there; otherwise the gap to the nearest address where it
    fits is filled with another unit that fits, sharers first, else with a hot
    or cold flash section outside every hot set. A unit with no fitting address
    splits in two; a single section with none is placed where it is.
    While lines exceed the bound, the plan repeats with the units holding them
    moved ahead, up to `PLAN_ROUNDS` times; the plan with fewest such lines wins.
    """
    promoted: dict[tuple[Origin, str], int] = {}
    best: tuple[int, list[tuple[Origin, str]]] | None = None
    for _ in range(PLAN_ROUNDS):
        planner = _Planner(program, hot_sets, alignments, start, promoted)
        order = planner.run()
        over, units = planner.violations()
        if best is None or over < best[0]:
            best = (over, order)
        if not over:
            break
        for unit in units:
            promoted[unit] = promoted.get(unit, 0) + 1
    return best[1]


def render_fragment(program: Program, order: list[tuple[Origin, str]], end: int) -> str:
    """Linker-script statements placing `order`, then an end-address check."""
    out = ["/* Generated by tools/hot_layout.py from the probe link. */"]
    for key in order:
        section = program.by_key[key]
        origin = section.origin
        if origin.member:
            out.append(f"*{origin.archive}:{origin.member}({section.name})")
        else:
            out.append(f"*{_basename(origin.path)}({section.name})")
    out.append(f"ASSERT(ABSOLUTE(.) == 0x{end:08x}, "
               f"\"hot code layout differs from its plan\");")
    return "\n".join(out) + "\n"


def planned_end(program: Program, order: list[tuple[Origin, str]],
                alignments: dict[tuple[Origin, str], int], start: int) -> int:
    """Address after the last section of `order` placed from `start`."""
    at = start
    for key in order:
        at = _align(at, alignments[key]) + program.by_key[key].size
    return at


# ---------------------------------------------------------------------------
# Toolchain driving
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Toolchain:
    """Paths of the binutils this module runs."""

    objdump: str
    readelf: str

    @classmethod
    def beside(cls, compiler: str, search_path: str | None = None) -> Toolchain:
        """Binutils next to an `arm-none-eabi-gcc`/`g++` named by path or PATH."""
        tools = cls(sibling_tool(compiler, "objdump", search_path),
                    sibling_tool(compiler, "readelf", search_path))
        for tool in (tools.objdump, tools.readelf):
            if not Path(tool).is_file() and not shutil.which(tool, path=search_path):
                raise LayoutError(f"no binutils tool at {tool}")
        return tools

    def run(self, args: list[str]) -> str:
        """Stdout of one tool invocation."""
        result = subprocess.run(args, capture_output=True, text=True,
                                encoding="utf-8", errors="replace", check=False)
        if result.returncode:
            raise LayoutError(f"{' '.join(args[:2])} failed: {result.stderr.strip()}")
        return result.stdout


def startup_ranges(symbols: dict[str, int]) -> tuple[tuple[int, int], ...]:
    """Explicit startup graph ranges, or none for links without the metadata."""
    ranges = []
    for begin, end in STARTUP_RANGE_SYMBOLS:
        if begin not in symbols and end not in symbols:
            continue
        if begin not in symbols or end not in symbols:
            raise LayoutError(f"incomplete startup range: {begin}, {end}")
        if not FLASH_START <= symbols[begin] <= symbols[end]:
            raise LayoutError(f"invalid startup range: {begin}, {end}")
        ranges.append((symbols[begin], symbols[end]))
    return tuple(ranges)


def load_program(tools: Toolchain, elf: str, map_path: str) -> tuple[Program, dict[str, int]]:
    """Program model and map symbols of a linked image."""
    sections, symbols = parse_map(Path(map_path).read_text(encoding="utf-8", errors="replace"))
    symbol_table = tools.run([tools.readelf, "-sW", elf])
    functions = parse_functions(symbol_table)
    references, data = parse_references(tools.run(
        [tools.objdump, "-d", "--no-show-raw-insn", "-j", FLASH_OUTPUT,
         "-j", ITCM_OUTPUT, elf]), {f.address for f in functions})
    objects = parse_objects(symbol_table)
    words = read_words(Path(elf).read_bytes(), objects)
    contents = {start: (size, words[start]) for start, size in objects if start in words}
    return Program(functions, sections, references, data, contents, startup_ranges(symbols)), symbols


def load_alignments(tools: Toolchain, keys: set[tuple[Origin, str]],
                    root: str) -> dict[tuple[Origin, str], int]:
    """Alignment of every section in `keys`, read from its input file."""
    alignments: dict[tuple[Origin, str], int] = {}
    for path in sorted({key[0].path for key in keys}):
        resolved = path if os.path.isabs(path) else os.path.join(root, path)
        alignments.update(parse_section_alignments(
            tools.run([tools.objdump, "-h", "-w", resolved]), path))
    missing = sorted(key[1] for key in keys if key not in alignments)
    if missing:
        raise LayoutError(f"no section alignment for {missing[0]} "
                          f"(+{len(missing) - 1} more)")
    return alignments


def write_plan(tools: Toolchain, elf: str, map_path: str, effects: list[str],
               fragment: str, root: str) -> None:
    """Plan from a probe link and write the fragment."""
    program, symbols = load_program(tools, elf, map_path)
    if START_SYMBOL not in symbols:
        raise LayoutError(f"probe map defines no {START_SYMBOL}")
    hot_sets = program.hot_sets(effects)
    start = symbols[START_SYMBOL]
    placeable = {key for hot in hot_sets.values() for _, key in hot.ranked}
    placeable |= {s.key for s in program.sections if s.output == FLASH_OUTPUT
                  and s.name.startswith((".text.hot", ".text.unlikely"))}
    alignments = load_alignments(tools, placeable, root)
    order = plan(program, hot_sets, alignments, start)
    end = planned_end(program, order, alignments, start)
    Path(fragment).write_text(render_fragment(program, order, end), encoding="utf-8")


def check(tools: Toolchain, elf: str, map_path: str, effects: list[str],
          expected: list[str], env: str) -> tuple[list[str], bool]:
    """Gate lines for a linked image and whether it passes."""
    program, _ = load_program(tools, elf, map_path)
    hot_sets = program.hot_sets(effects)
    missing = [e for e in expected if e not in hot_sets]
    reports = evaluate(program, hot_sets)
    return (format_reports(env, reports, missing),
            not missing and all(r.passed for r in reports))


def roster(root: str) -> list[str]:
    """Every effect class in core/platform/effect_identity.h."""
    text = Path(root, "core/platform/effect_identity.h").read_text(encoding="utf-8")
    return parse_effect_list(text, "HS_EFFECT_IDENTITY_LIST")


def phantasm_playlist(root: str) -> list[str]:
    """Effect classes of the Phantasm device playlist."""
    text = Path(root, PHANTASM_PLAYLIST_SOURCE).read_text(encoding="utf-8")
    names = playlist_shows(text)
    if not names:
        raise LayoutError(f"no show rows in {PHANTASM_PLAYLIST_SOURCE}")
    return names


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("command", choices=("check", "plan"))
    parser.add_argument("--elf", required=True)
    parser.add_argument("--map", required=True)
    parser.add_argument("--compiler", default="arm-none-eabi-gcc",
                        help="gcc whose sibling binutils to run")
    parser.add_argument("--root", default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument("--expect", action="append", default=[],
                        help="effect that must have code in the image; "
                             "'playlist' expands to the Phantasm playlist")
    parser.add_argument("--fragment", help="plan: fragment path to write")
    args = parser.parse_args(argv)
    try:
        tools = Toolchain.beside(args.compiler)
        expected = []
        for name in args.expect:
            expected.extend(phantasm_playlist(args.root) if name == "playlist" else [name])
        if args.command == "plan":
            if not args.fragment:
                parser.error("plan needs --fragment")
            write_plan(tools, args.elf, args.map, roster(args.root), args.fragment, args.root)
            return 0
        lines, passed = check(tools, args.elf, args.map, roster(args.root), expected,
                              Path(args.elf).stem)
    except (LayoutError, OSError) as exc:
        print(f"hot_layout: {exc}", file=sys.stderr)
        return 2
    print("\n".join(lines))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
