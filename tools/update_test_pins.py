#!/usr/bin/env python3
"""Re-measure the test-suite pins in tests/run_tests.cpp from real runs.

Two kinds of number live in the roster, and only one of them can be read out
of the tree.

CASE PINS -- the second column of every HS_TEST_MODULE_LIST row, plus every
entry in tests/off_roster_headers.cmake -- are the EXACT count of
`void test_*(` / `void check_*(` / `void case_*(` definitions in a header.
This script counts them the way tests/check_case_calls.cmake counts them, so
no run is needed.

ASSERTION FLOORS -- the fourth column, or one of the tier constants the
effects rows select between -- are LOWER BOUNDS that have to hold in every
configuration that moves a module's count: the HS_SMOKE_FRAMES window, Debug
-O0 against RelWithDebInfo -O2, the effects depth tier, and the NDEBUG-gated
cases in test_memory.h. They are therefore MEASURED, never derived. Each
--run or --log is one configuration's run_tests output, and a floor's
measurement is the MINIMUM over the sources that apply to it.

A floor being a bound and not an expectation, the default is to CHECK it:
measured-below-pinned is reported as a violation and everything else is left
alone, headroom included. --remeasure writes the measured minimum instead, for
a module whose floor was never measured or whose cases were deliberately cut,
less --margin percent so a benign shift does not have to red CI.

Headroom is expected, but a floor a module clears several times over has
stopped bounding it: most of the module can be deleted and the gate stays
green. --check-drift turns that into a failure, on the two-sided allowance
below. --remeasure --check-drift --margin <pct> --write is its repair, and
re-pins the drifted floors ONLY -- plain --remeasure rewrites every measured
floor, which would walk one sitting exactly on its measurement back down by the
margin. Drift is measured the same way the floor is, so it is only as sound as
the sources: a floor tracking a configuration none of them ran reads as
drifted.

Two rules keep a partial measurement from writing a floor no other
configuration clears. A floor moves only when at least two applying sources
observed it, so one configuration cannot lift it over what another one
produces. An entry no source observed -- every module outside a narrowed run's
list -- is carried through unchanged and named in the report, never dropped or
reset.

The script reports and exits unless --write is passed, and --write refuses
under CI: a ratchet a workflow can regenerate gates nothing.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
TESTS = ROOT / "tests"
RUN_TESTS = TESTS / "run_tests.cpp"
OFF_ROSTER = TESTS / "off_roster_headers.cmake"

TIERS = ("quick", "full")

ROSTER_BEGIN = "#define HS_TEST_MODULE_LIST(X)"
ROSTER_END = "#define HS_TEST_MODULE_ENTRY"

# Mirrors of tests/check_case_calls.cmake's scan, applied in the same order it
# applies them: string bodies are blanked before comment spans, so a `//` or
# `/*` inside a message cannot open a span.
STRING_BODY = re.compile(r'"([^"\\\n]|\\.)*"')
BLOCK_COMMENT = re.compile(r"/\*[^*]*\*+([^/*][^*]*\*+)*/")
LINE_COMMENT = re.compile(r"//[^\n]*")
CASE_DEF = re.compile(
    r"\n[ \t]*(?:template[ \t]*<[^\n]*>[ \t]*)?(?:(?:inline|static)[ \t]+)*void[ \t\r\n]+"
    r"((?:test|check|case)_[A-Za-z0-9_]+)\("
)
ENTRY_FN = re.compile(r"\n[ \t]*(?:inline )?(?:static )?int (run_[A-Za-z0-9_]+_tests)\(")

ROSTER_ROW = re.compile(r'X\("([A-Za-z0-9_]+)"')
TIER_CONSTANT = re.compile(r"^constexpr int ([A-Za-z_][A-Za-z0-9_]*) = (\d+);", re.M)
TIER_SELECT = re.compile(
    r"\?[\s\\]*([A-Za-z_][A-Za-z0-9_]*)[\s\\]*:[\s\\]*([A-Za-z_][A-Za-z0-9_]*)"
)
INTEGER_FIELD = re.compile(r"[\s\\]*(\d+)[\s\\]*")
OFF_ROSTER_ENTRY = re.compile(r'"([A-Za-z0-9_.]+)=([0-9]+)"')

# end_module()'s per-module footer. The SKIPPED variant is printed whenever any
# case skipped, so it does not by itself mark a whole-suite skip; a module that
# reported zero assertions does, and its count is no measurement.
MODULE_FOOTER = re.compile(
    r"^=== ([A-Za-z0-9_]+): (\d+) passed, (\d+) failed(?:, (\d+) SKIPPED)? ===\s*$",
    re.M,
)
NO_ASSERTIONS = "NO ASSERTIONS RAN"

# Drift allowance. Both bounds have to be cleared before a floor counts as
# drifted: the fraction is what makes a floor stop bounding its module, and the
# absolute keeps a small module (a floor of 12, a measurement of 24) out of it.
DRIFT_FRACTION = 0.25
DRIFT_MINIMUM = 200


@dataclass
class Pin:
    """One editable integer literal in a pin source file."""

    path: Path
    name: str
    start: int
    end: int
    value: int


@dataclass
class Floor:
    """A roster assertion floor and the runs that can measure it."""

    pin: Pin
    module: str
    tier: str | None = None


@dataclass
class Source:
    """One configuration's run_tests output."""

    label: str
    tier: str
    origin: str
    counts: dict[str, int] = field(default_factory=dict)
    suite_skips: list[str] = field(default_factory=list)


def strip_code(text: str) -> str:
    """Blank string bodies, then drop block and line comment spans."""
    text = STRING_BODY.sub('""', text)
    text = BLOCK_COMMENT.sub("\n", text)
    return LINE_COMMENT.sub("", text)


def scan_headers() -> dict[str, tuple[str | None, int]]:
    """Return header name -> (its run_*_tests entry point or None, case sites)."""
    headers: dict[str, tuple[str | None, int]] = {}
    for path in sorted(TESTS.rglob("*.h")) + sorted(TESTS.rglob("*.hpp")):
        text = strip_code(path.read_text(encoding="utf-8"))
        seen: list[str] = []
        for name in CASE_DEF.findall(text):
            if name not in seen:
                seen.append(name)
        entry = ENTRY_FN.search(text)
        headers[path.name] = (entry.group(1) if entry else None, len(seen))
    return headers


def field_spans(row: str) -> list[tuple[int, int]]:
    """Return the spans of a roster row's top-level comma-separated fields."""
    spans: list[tuple[int, int]] = []
    depth = 0
    start = 0
    in_string = False
    for index, char in enumerate(row):
        if in_string:
            if char == '"':
                in_string = False
        elif char == '"':
            in_string = True
        elif char == "(":
            depth += 1
            if depth == 1:
                start = index + 1
        elif char == ")":
            depth -= 1
            if depth == 0:
                spans.append((start, index))
                break
        elif char == "," and depth == 1:
            spans.append((start, index))
            start = index + 1
    return spans


def integer_pin(name: str, text: str, offset: int, span: tuple[int, int]) -> Pin | None:
    """Return the pin for a roster field holding a bare integer, else None."""
    found = INTEGER_FIELD.fullmatch(text, offset + span[0], offset + span[1])
    if found is None:
        return None
    return Pin(RUN_TESTS, name, found.start(1), found.end(1), int(found.group(1)))


def parse_roster(text: str) -> tuple[list[tuple[Pin, int]], list[Floor], list[str]]:
    """Return (case pins paired with the count scanned out of their header,
    assertion floors, errors) read from run_tests.cpp."""
    begin = text.find(ROSTER_BEGIN)
    end = text.find(ROSTER_END)
    if begin < 0 or end < begin:
        return [], [], [f"cannot find the {ROSTER_BEGIN} block in {RUN_TESTS.name}"]

    constants = {
        match.group(1): Pin(
            RUN_TESTS, match.group(1), match.start(2), match.end(2), int(match.group(2))
        )
        for match in TIER_CONSTANT.finditer(text[:begin])
    }
    headers = scan_headers()
    by_entry_fn = {
        entry_fn: (name, sites)
        for name, (entry_fn, sites) in headers.items()
        if entry_fn
    }

    cases: list[tuple[Pin, int]] = []
    floors: list[Floor] = []
    errors: list[str] = []
    covered: set[str] = set()
    for head in ROSTER_ROW.finditer(text, begin, end):
        module = head.group(1)
        spans = field_spans(text[head.start() : end])
        if len(spans) != 4:
            errors.append(
                f"roster row for '{module}' has {len(spans)} field(s); every row "
                f'must read X("<name>", <case sites>, <entry point>, <floor>)'
            )
            continue
        offset = head.start()

        entry_fn = text[offset + spans[2][0] : offset + spans[2][1]]
        entry_fn = entry_fn.strip(" \t\r\n\\").split("::")[-1]
        case_pin = integer_pin(f"{module} cases", text, offset, spans[1])
        if case_pin is None:
            errors.append(f"roster row for '{module}' has no integer case-site pin")
        elif entry_fn not in by_entry_fn:
            errors.append(
                f"roster module '{module}' names {entry_fn}, which no header "
                f"under {TESTS.name}/ defines"
            )
        else:
            header, sites = by_entry_fn[entry_fn]
            covered.add(header)
            case_pin.name = f"{module} cases ({header})"
            cases.append((case_pin, sites))

        floor_pin = integer_pin(f"{module} floor", text, offset, spans[3])
        if floor_pin is not None:
            floors.append(Floor(floor_pin, module))
            continue
        floor_text = text[offset + spans[3][0] : offset + spans[3][1]]
        tiered = TIER_SELECT.search(floor_text)
        if tiered is None:
            errors.append(
                f"roster row for '{module}' has a floor this script cannot read: "
                f"{' '.join(floor_text.split())}"
            )
            continue
        for name, tier in ((tiered.group(1), "full"), (tiered.group(2), "quick")):
            if name in constants:
                floors.append(Floor(constants[name], module, tier))
            else:
                errors.append(f"roster module '{module}' selects undefined {name}")

    errors += bind_off_roster(headers, covered, cases)
    return cases, floors, errors


def bind_off_roster(
    headers: dict[str, tuple[str | None, int]],
    covered: set[str],
    cases: list[tuple[Pin, int]],
) -> list[str]:
    """Add the off_roster_headers.cmake pins for every header the roster misses."""
    errors: list[str] = []
    text = OFF_ROSTER.read_text(encoding="utf-8")
    pinned = {match.group(1): match for match in OFF_ROSTER_ENTRY.finditer(text)}
    for header, (_, sites) in sorted(headers.items()):
        if header in covered:
            continue
        match = pinned.get(header)
        if match is None:
            errors.append(
                f"{OFF_ROSTER.name} pins no case count for {header}; add "
                f'"{header}={sites}" to HS_OFF_ROSTER_HEADERS'
            )
            continue
        pin = Pin(
            OFF_ROSTER,
            f"{header} cases",
            match.start(2),
            match.end(2),
            int(match.group(2)),
        )
        cases.append((pin, sites))
    return errors


def parse_run_output(label: str, tier: str, origin: str, text: str) -> Source:
    """Turn one run_tests output into a Source, or exit if it is no measurement."""
    source = Source(label, tier, origin)
    failed_modules: set[str] = set()
    for name, passed, failed, _ in MODULE_FOOTER.findall(text):
        total = int(passed) + int(failed)
        if int(failed):
            failed_modules.add(name)
        if total == 0:
            source.suite_skips.append(name)
            continue
        source.counts[name] = min(total, source.counts.get(name, total))
    if failed_modules:
        raise SystemExit(
            f"[test-pins] {origin}: {', '.join(sorted(failed_modules))} reported "
            f"failures; a red run is not a measurement"
        )
    if NO_ASSERTIONS in text:
        raise SystemExit(
            f"[test-pins] {origin}: a module ran no assertions at all; fix the run "
            f"before measuring from it"
        )
    if not source.counts:
        raise SystemExit(
            f"[test-pins] {origin}: no module footer found; is this run_tests output?"
        )
    return source


def read_source(spec: str, index: int, modules: list[str], execute: bool) -> Source:
    """Build a Source from a TIER=PATH spec, running the binary when asked."""
    tier, separator, target = spec.partition("=")
    if not separator or tier not in TIERS:
        raise SystemExit(
            f"[test-pins] '{spec}' is not <tier>=<path> with tier one of "
            f"{'/'.join(TIERS)}"
        )
    path = Path(target)
    label = f"{tier}#{index}"
    if not execute:
        return parse_run_output(
            label, tier, str(path), path.read_text(encoding="utf-8", errors="replace")
        )
    environment = os.environ.copy()
    environment["HS_EFFECTS_FULL"] = "1" if tier == "full" else "0"
    origin = " ".join([str(path), *modules])
    print(f"[test-pins] {label}: running {origin} with HS_EFFECTS_FULL="
          f"{environment['HS_EFFECTS_FULL']}")
    result = subprocess.run(
        [str(path), *modules],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        env=environment,
    )
    if result.returncode != 0:
        sys.stdout.write(result.stdout)
        raise SystemExit(f"[test-pins] {origin} exited {result.returncode}")
    return parse_run_output(label, tier, origin, result.stdout)


def measure_floor(floor: Floor, sources: list[Source]) -> tuple[int | None, str]:
    """Return (measured minimum, note); None leaves the committed value alone."""
    if floor.pin.value == 0:
        return None, "deliberately unfloored (0)"
    scope = "" if floor.tier is None else f" {floor.tier}-tier"
    applying = [s for s in sources if floor.tier is None or s.tier == floor.tier]
    if len(applying) < 2:
        return None, f"{len(applying)}{scope} source(s); two are required"
    observed = [s for s in applying if floor.module in s.counts]
    missing = [s.label for s in applying if floor.module not in s.counts]
    if not observed:
        return None, f"unobserved by {', '.join(missing)}"
    if len(observed) < 2:
        return None, (f"observed only by {observed[0].label} of "
                      f"{len(applying)}{scope} source(s); two are required")
    note = f"min over {len(observed)} source(s)"
    if missing:
        note += f"; unobserved by {', '.join(missing)}"
    return min(s.counts[floor.module] for s in observed), note


def drifted(pinned: int, measured: int) -> bool:
    """True when a measurement clears its floor by more than the allowance."""
    slack = measured - pinned
    return slack >= DRIFT_MINIMUM and slack > pinned * DRIFT_FRACTION


def with_margin(measured: int, margin: float) -> int:
    """The value --remeasure writes: the measurement less margin percent."""
    return measured - int(measured * margin / 100.0)


def report(title: str, changes: list[str], held: list[str]) -> None:
    """Print one section of the report."""
    print(f"[test-pins] {title}: {len(changes)} change(s), "
          f"{len(held)} carried through unchanged")
    for line in changes:
        print(f"  {line}")
    for line in held:
        print(f"  = {line}")


def apply(edits: list[tuple[Pin, int]]) -> None:
    """Splice new values into the pin sources, last offset first."""
    by_path: dict[Path, list[tuple[Pin, int]]] = {}
    for pin, value in edits:
        by_path.setdefault(pin.path, []).append((pin, value))
    for path, group in by_path.items():
        text = path.read_text(encoding="utf-8")
        for pin, value in sorted(group, key=lambda e: e[0].start, reverse=True):
            text = text[: pin.start] + str(value) + text[pin.end :]
        path.write_text(text, encoding="utf-8", newline="\n")
        print(f"[test-pins] wrote {len(group)} pin(s) to {path.relative_to(ROOT)}")


def reformat() -> None:
    """Re-run clang-format over run_tests.cpp; a widened pin rewraps its row."""
    want = subprocess.run(
        [sys.executable, str(TOOLS / "build_pins.py"), "clang-format"],
        capture_output=True, text=True, check=True,
    ).stdout.strip()
    try:
        reported = subprocess.run(
            ["clang-format", "--version"], capture_output=True, text=True, check=True
        ).stdout
    except (OSError, subprocess.SubprocessError):
        reported = ""
    found = re.search(r"\d[\w.]*", reported)
    if found is None or found.group() != want:
        print(f"[test-pins] clang-format {want} is not on PATH, so a pin whose "
              f"width changed leaves its roster row unwrapped; run "
              f"`just clang-format` before committing")
        return
    subprocess.run(["clang-format", "-i", "--style=file", str(RUN_TESTS)], check=True)
    print(f"[test-pins] reformatted {RUN_TESTS.relative_to(ROOT)}")


def build_sources(runs: list[str], logs: list[str], modules: list[str]) -> list[Source]:
    """Read every measurement source named on the command line."""
    specs = [(spec, True) for spec in runs] + [(spec, False) for spec in logs]
    return [
        read_source(spec, index, modules, execute)
        for index, (spec, execute) in enumerate(specs, 1)
    ]


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Re-measure the tests/run_tests.cpp case pins and assertion "
                    "floors from two or more configurations."
    )
    parser.add_argument(
        "--run", action="append", default=[], metavar="TIER=PATH",
        help="run_tests binary to execute with HS_EFFECTS_FULL set from TIER "
             "(quick or full); repeatable")
    parser.add_argument(
        "--log", action="append", default=[], metavar="TIER=PATH",
        help="captured run_tests output for one configuration; repeatable")
    parser.add_argument(
        "--module", action="append", default=[], metavar="NAME",
        help="restrict every --run to these roster modules; repeatable")
    parser.add_argument(
        "--remeasure", action="store_true",
        help="set every measured floor to its minimum instead of checking that "
             "the committed one holds")
    parser.add_argument(
        "--check-drift", action="store_true",
        help=f"also fail on a floor its module clears by both {DRIFT_MINIMUM} "
             f"assertions and {DRIFT_FRACTION * 100:.0f} percent, which no longer "
             f"bounds it; narrows --remeasure to those floors")
    parser.add_argument(
        "--margin", type=float, default=0.0, metavar="PCT",
        help="write each --remeasure floor this far under the measured "
             "minimum (default 0)")
    parser.add_argument(
        "--write", action="store_true",
        help="rewrite the pins in place (default: report the diff only)")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not 0.0 <= args.margin < 100.0:
        print(f"[test-pins] --margin {args.margin} is not a percentage in [0, 100)")
        return 1
    if args.write and os.environ.get("CI"):
        print("[test-pins] --write is refused under CI: the pins are a ratchet, "
              "and one a workflow can regenerate gates nothing. Re-measure on a "
              "developer machine and commit the result.")
        return 1

    cases, floors, errors = parse_roster(RUN_TESTS.read_text(encoding="utf-8"))
    if errors:
        for line in errors:
            print(f"[test-pins] {line}")
        return 1

    sources = build_sources(args.run, args.log, args.module)
    if len(sources) == 1:
        print("[test-pins] one measurement source: a floor is the minimum over "
              "the configurations that move a module's count, so at least two "
              "are required. Pass another --run/--log.")
        return 1
    for source in sources:
        print(f"[test-pins] {source.label}: {len(source.counts)} module(s) "
              f"measured from {source.origin}")
        if source.suite_skips:
            print(f"[test-pins]   whole-suite skips, not measured: "
                  f"{', '.join(source.suite_skips)}")

    edits: list[tuple[Pin, int]] = []
    drift = [(pin, sites) for pin, sites in cases if pin.value != sites]
    changes = [f"{pin.name}: {pin.value} -> {sites}" for pin, sites in drift]
    edits += drift
    report(f"{len(cases)} case pin(s)", changes, [])

    changes: list[str] = []
    held: list[str] = []
    violations: list[str] = []
    drifts: list[str] = []
    headroom: list[int] = []
    for floor in floors if sources else []:
        measured, note = measure_floor(floor, sources)
        if measured is None:
            held.append(f"{floor.pin.name} {floor.pin.value}: {note}")
            continue
        is_drift = False
        if measured < floor.pin.value:
            violations.append(
                f"{floor.pin.name}: pinned {floor.pin.value}, measured "
                f"{measured} ({note})")
        elif args.check_drift and drifted(floor.pin.value, measured):
            is_drift = True
            drifts.append(
                f"{floor.pin.name}: pinned {floor.pin.value}, measured "
                f"{measured}, {measured - floor.pin.value} of slack ({note})")
        target = with_margin(measured, args.margin)
        repairing = is_drift or not args.check_drift
        if args.remeasure and repairing and target != floor.pin.value:
            direction = "RAISED" if target > floor.pin.value else "LOWERED"
            edits.append((floor.pin, target))
            changes.append(
                f"{floor.pin.name}: {floor.pin.value} -> {target} {direction} "
                f"by {abs(target - floor.pin.value)} ({note})")
        elif measured >= floor.pin.value:
            headroom.append(measured - floor.pin.value)
    if not sources:
        print(f"[test-pins] {len(floors)} assertion floor(s): unmeasured, no "
              f"--run/--log given")
    else:
        span = f" (headroom {min(headroom)}..{max(headroom)})" if headroom else ""
        print(f"[test-pins] {len(floors)} assertion floor(s): {len(headroom)} "
              f"hold{span}, {len(violations)} violated, {len(drifts)} drifted, "
              f"{len(changes)} re-measured, {len(held)} carried through unchanged")
        for line in violations:
            print(f"  ! {line}")
        for line in drifts:
            print(f"  ~ {line}")
        for line in changes:
            print(f"  {line}")
        for line in held:
            print(f"  = {line}")

    if violations and (args.check_drift or not args.remeasure):
        print(f"[test-pins] {len(violations)} floor(s) no longer hold. Restore the "
              f"cases that went missing, or lower them deliberately with "
              f"--remeasure --write.")
        return 1
    if drifts and not args.remeasure:
        print(f"[test-pins] {len(drifts)} floor(s) have drifted so far under "
              f"their module that most of it could be deleted green. Re-pin them "
              f"with --remeasure --check-drift --margin <pct> --write.")
        return 1
    if not edits:
        print("[test-pins] every pin already matches; nothing to write")
        return 0
    if not args.write:
        print(f"[test-pins] {len(edits)} pin(s) would move; re-run with --write to "
              f"apply. A LOWERED floor belongs only alongside the case removal "
              f"that caused it.")
        return 1
    apply(edits)
    reformat()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
