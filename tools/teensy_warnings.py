#!/usr/bin/env python3
"""Teensy firmware warning-hygiene ratchet (docs/teensy_ci_gate_spec.md §7.2).

Policy is a BASELINE RATCHET, not -Werror: capture the current first-party
warning set once, then fail only on warnings NOT in that baseline. New warnings
are blocked; fixing a baseline warning is fine.

Two properties the spec calls load-bearing:
  * SET-based, not line-ordered. PlatformIO builds in parallel (-j) so warning
    emission order is nondeterministic; an ordered diff would flap green/red on
    identical inputs. We compare normalized SETS.
  * Normalized to drop volatile bits (absolute path prefix, line:col) so the
    fingerprint is stable across unrelated edits, while the file + message + flag
    that identify the warning are kept.

First-party only: warnings outside core/ effects/ hardware/ targets/ (i.e.
FastLED and the Teensy core) are dropped — the independent backstop to the
`-isystem` plan, so the ratchet is robust even if the -isystem demotion proves
awkward under PlatformIO (§7.2).

This ratchet must run on a COLD build: a cached TU emits no warnings, so a warm
build hides header-introduced warnings (§7.2, §15). `PLATFORMIO_BUILD_CACHE_DIR=`
does not achieve that — PlatformIO ignores an empty sysenvvar and falls back to
platformio.ini's `build_cache_dir` — so the capture must DELETE `.pio/build_cache`
and `.pio/build` first. The capture audit below enforces coldness from the log
itself rather than trusting the caller to have done so.
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path

FIRST_PARTY = ("core/", "effects/", "hardware/", "targets/")

# Library/toolchain roots: a path through any of these is third-party even when a
# nested dir reuses a first-party name (e.g. .platformio/lib/Foo/effects/x.h).
THIRD_PARTY = ("lib/", ".platformio/", "packages/")

# gcc: "<path>:<line>[:<col>]: warning: <message> [-Wflag]"
_WARNING_RE = re.compile(r"^(.*?):(\d+):(?:\d+:)?\s*warning:\s*(.*)$")

# PlatformIO's non-verbose step line: "Compiling <object>". `pio run -v` prints
# the raw compiler command instead, and never these.
_PIO_STEP_RE = re.compile(r"^\s*Compiling\s+(\S+)")

# `-c` as a standalone flag on a compiler command line. In the `-v` echo it is
# NOT adjacent to the source: the flags follow it and the source comes last.
_DASH_C_RE = re.compile(r"(?:^|\s)-c(?:\s|$)")

# A positional source argument anywhere on that command line. Tokens starting
# with `-` are flags (`-Icore`, `-x c++`), never the translation unit.
_SOURCE_RE = re.compile(r"(?:^|\s)(?!-)(\S+\.(?:cpp|cc|cxx|c|S))(?=\s|$)")

# PlatformIO's per-environment banner. It opens every env's section AND prints
# that env's resolved options, `build_src_filter` among them — which is what
# makes the expected translation-unit set derivable from the log alone:
#   Processing phantasm (board: teensy40; build_src_filter: -<*>, +<core/...>; ...)
_ENV_HEADER_RE = re.compile(r"^Processing\s+(\S+)\s+\((.*)\)\s*$")

# `key: value` pairs in that banner are `; `-separated, but a value may itself
# contain `, ` and spaces, so the value runs to the next `; <key>: `.
_SRC_FILTER_RE = re.compile(
    r"(?:^|;\s*)build_src_filter:\s*(.*?)(?=;\s*[\w.]+:\s|$)")

# One src-filter term: `+<path>` includes, `-<path>` excludes.
_SRC_FILTER_TERM_RE = re.compile(r"([+-])<([^>]*)>")

# SCons object-cache hit:  Retrieved `.pio/build/x/src/core/engine/memory.cpp.o' from cache
_CACHE_HIT_RE = re.compile(r"^\s*Retrieved\s+[`'\"](.+?)['\"]\s+from cache\s*$")

_GLOB_CHARS = "*?["


_FIRST_PARTY_DIRS = frozenset(fp.rstrip("/") for fp in FIRST_PARTY)
_THIRD_PARTY_DIRS = frozenset(tp.rstrip("/") for tp in THIRD_PARTY)


def _relativize(path: str) -> str | None:
    """Return the repo-root-relative path if first-party, else None.

    Anchored at the FIRST segment naming a first-party top-level dir, so the full
    nested path is kept (targets/Phantasm/effects/Foo.h != top-level effects/Foo.h)
    and not collapsed to the first matching segment.
    """
    segs = path.replace("\\", "/").split("/")
    if _THIRD_PARTY_DIRS & set(segs):
        return None
    for i, seg in enumerate(segs):
        if seg in _FIRST_PARTY_DIRS:
            return "/".join(segs[i:])
    return None


def normalize(line: str) -> str | None:
    """Normalize one compiler line to a stable first-party warning key, or None."""
    m = _WARNING_RE.match(line.strip())
    if not m:
        return None
    rel = _relativize(m.group(1))
    if rel is None:
        return None
    return f"{rel}: warning: {m.group(3)}".rstrip()


def extract_warnings(build_log: str) -> set[str]:
    """The deduplicated, normalized, first-party warning set from a build log."""
    out: set[str] = set()
    for line in build_log.splitlines():
        key = normalize(line)
        if key is not None:
            out.add(key)
    return out


def compiled_paths(line: str) -> list[str]:
    """The path(s) one build-log line shows being compiled, or an empty list.

    Two log shapes: PlatformIO's non-verbose step line names the OBJECT, while
    the `pio run -v` echo of the raw compiler command names the SOURCE.
    """
    m = _PIO_STEP_RE.match(line)
    if m:
        return [m.group(1)]
    if not _DASH_C_RE.search(line):
        return []
    return _SOURCE_RE.findall(line)


def count_first_party_compiles(build_log: str) -> int:
    """Compiler invocations on first-party sources visible in the log.

    The ratchet's green and a broken capture both yield an empty warning set, so
    the comparison is only meaningful once the log is known to hold a build that
    could have emitted first-party warnings at all. A third-party-only build
    (FastLED, the Teensy core) is NOT evidence: those TUs cannot emit a warning
    the ratchet would ever look at.
    """
    n = 0
    for line in build_log.splitlines():
        if any(_relativize(p) is not None for p in compiled_paths(line)):
            n += 1
    return n


class CaptureError(Exception):
    """The log cannot be audited for coldness, so its warning set proves nothing."""


@dataclass(frozen=True)
class EnvSection:
    """One `Processing <env> (...)` banner and the build lines that follow it."""

    name: str
    header: str
    lines: tuple[str, ...]


@dataclass(frozen=True)
class EnvAudit:
    """Expected vs actual first-party translation units for one environment."""

    name: str
    declared: frozenset[str]
    compiled: frozenset[str]
    from_cache: frozenset[str]

    @property
    def missing(self) -> frozenset[str]:
        return self.declared - self.compiled


@dataclass(frozen=True)
class CaptureAudit:
    """Whole-log coldness evidence."""

    envs: tuple[EnvAudit, ...]
    compiles: int
    cache_hits: int

    @property
    def declared(self) -> int:
        return sum(len(e.declared) for e in self.envs)

    @property
    def missing(self) -> tuple[EnvAudit, ...]:
        return tuple(e for e in self.envs if e.missing)

    @property
    def missing_count(self) -> int:
        return sum(len(e.missing) for e in self.envs)

    @property
    def first_party_cache_hits(self) -> int:
        return sum(len(e.from_cache) for e in self.envs)


def parse_env_sections(build_log: str) -> list[EnvSection]:
    """Split a `pio run` log into its per-environment sections.

    PlatformIO builds environments sequentially, so a banner closes the previous
    section. Lines before the first banner belong to no environment.
    """
    sections: list[EnvSection] = []
    name: str | None = None
    header = ""
    body: list[str] = []
    for line in build_log.splitlines():
        m = _ENV_HEADER_RE.match(line)
        if m:
            if name is not None:
                sections.append(EnvSection(name, header, tuple(body)))
            name, header, body = m.group(1), m.group(2), []
        elif name is not None:
            body.append(line)
    if name is not None:
        sections.append(EnvSection(name, header, tuple(body)))
    return sections


def declared_first_party_sources(section: EnvSection) -> set[str]:
    """The first-party translation units this env's `build_src_filter` selects.

    Derived from PlatformIO's own banner, so adding a TU or an environment moves
    the expectation with no second place to edit. A glob that could select a
    first-party source is not countable from the log and raises rather than
    silently lowering the bar.
    """
    m = _SRC_FILTER_RE.search(section.header)
    if m is None:
        raise CaptureError(
            f"environment '{section.name}' banner has no build_src_filter, so its "
            f"expected first-party translation units cannot be derived")
    sources: set[str] = set()
    for sign, pattern in _SRC_FILTER_TERM_RE.findall(m.group(1)):
        rel = _relativize(pattern)
        if any(ch in pattern for ch in _GLOB_CHARS):
            if rel is not None:
                raise CaptureError(
                    f"environment '{section.name}' build_src_filter term "
                    f"'{sign}<{pattern}>' is a first-party glob; the expected "
                    f"translation-unit set cannot be counted from the log")
            continue
        if rel is None:
            continue
        if sign == "+":
            sources.add(rel)
        else:
            sources.discard(rel)
    return sources


def _source_key(path: str) -> str | None:
    """A first-party path as a translation-unit key, or None.

    Both log shapes and the cache-hit line reduce to the same key: the object
    suffix is dropped so `<tu>.cpp.o` and `<tu>.cpp` compare equal against the
    sources build_src_filter declares.
    """
    rel = _relativize(path)
    if rel is None:
        return None
    return rel[:-2] if rel.endswith(".o") else rel


def compiled_first_party_sources(lines: tuple[str, ...] | list[str]) -> set[str]:
    """First-party translation units a compiler was actually invoked on."""
    out: set[str] = set()
    for line in lines:
        for path in compiled_paths(line):
            key = _source_key(path)
            if key is not None:
                out.add(key)
    return out


def cached_first_party_sources(lines: tuple[str, ...] | list[str]) -> set[str]:
    """First-party translation units whose object SCons served from the cache."""
    out: set[str] = set()
    for line in lines:
        m = _CACHE_HIT_RE.match(line)
        if m:
            key = _source_key(m.group(1))
            if key is not None:
                out.add(key)
    return out


def count_cache_hits(build_log: str) -> int:
    """Every SCons cache retrieval in the log, first-party or not."""
    return sum(1 for line in build_log.splitlines() if _CACHE_HIT_RE.match(line))


def audit_capture(build_log: str) -> CaptureAudit:
    """Compare each environment's declared first-party TUs against the compiles."""
    envs = []
    for section in parse_env_sections(build_log):
        envs.append(EnvAudit(
            name=section.name,
            declared=frozenset(declared_first_party_sources(section)),
            compiled=frozenset(compiled_first_party_sources(section.lines)),
            from_cache=frozenset(cached_first_party_sources(section.lines)),
        ))
    return CaptureAudit(
        envs=tuple(envs),
        compiles=count_first_party_compiles(build_log),
        cache_hits=count_cache_hits(build_log),
    )


def load_baseline(path: str | Path) -> set[str]:
    """Read the committed baseline set (ignoring blank and #-comment lines)."""
    text = Path(path).read_text(encoding="utf-8") if Path(path).exists() else ""
    return {ln.strip() for ln in text.splitlines()
            if ln.strip() and not ln.lstrip().startswith("#")}


def render_baseline(warnings: set[str]) -> str:
    """Render a baseline file: header + sorted, deduplicated warning set."""
    header = [
        "# Teensy firmware first-party warning baseline (tools/teensy_warnings.py).",
        "# Sorted, deduplicated, normalized (path:line:col stripped). The ratchet",
        "# fails only on warnings NOT listed here. Regenerate with --update-baseline",
        "# in the SAME PR as the change that legitimately alters the set (§7.2).",
        "",
    ]
    return "\n".join(header + sorted(warnings)) + "\n"


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Teensy warning-baseline ratchet.")
    p.add_argument("--build-log", required=True, help="compiler output of a COLD build")
    p.add_argument("--baseline", default="tools/teensy_warning_baseline.txt")
    p.add_argument("--update-baseline", action="store_true",
                   help="rewrite the baseline from this build's warning set")
    p.add_argument("--github", action="store_true", help="emit ::error:: annotations")
    args = p.parse_args(argv)

    build_log = Path(args.build_log).read_text(encoding="utf-8")
    prefix = "::error::" if args.github else ""
    compiles = count_first_party_compiles(build_log)
    if compiles == 0:
        print(f"{prefix}[teensy-warnings] FAIL - no first-party compiler "
              f"invocation in {args.build_log}: the capture broke (or the build "
              f"was fully cached), so its empty warning set proves nothing. "
              f"Drive this from a cold `pio run -v 2>&1 | tee` log.")
        return 1

    try:
        audit = audit_capture(build_log)
    except CaptureError as exc:
        print(f"{prefix}[teensy-warnings] FAIL - {exc} ({args.build_log}).")
        return 1
    if not audit.envs:
        print(f"{prefix}[teensy-warnings] FAIL - no `Processing <env> (...)` banner "
              f"in {args.build_log}: this is not a `pio run` capture, so the "
              f"expected first-party translation-unit set cannot be derived and a "
              f"partially cached build would read as green.")
        return 1
    if audit.missing:
        print(f"{prefix}[teensy-warnings] FAIL - the capture is not cold: "
              f"{audit.missing_count} of {audit.declared} first-party "
              f"translation unit(s) across {len(audit.envs)} environment(s) never "
              f"reached the compiler, so any warning they emit is absent from "
              f"{args.build_log}.")
        for env in audit.missing:
            print(f"  - {env.name}: {len(env.missing)} of {len(env.declared)} not "
                  f"compiled: {', '.join(sorted(env.missing))}")
        if audit.cache_hits:
            print(f"Cause: PlatformIO's object cache served {audit.cache_hits} "
                  f"object(s) ({audit.first_party_cache_hits} first-party) from "
                  f"platformio.ini's build_cache_dir.")
        else:
            print("Cause: an incremental build - the objects were already up to "
                  "date in .pio/build.")
        print("Delete .pio/build_cache and .pio/build before the capture. Note "
              "PLATFORMIO_BUILD_CACHE_DIR= does NOT disable the cache: PlatformIO "
              "ignores an empty sysenvvar, so build_cache_dir still wins.")
        return 1
    current = extract_warnings(build_log)

    if args.update_baseline:
        Path(args.baseline).write_text(
            render_baseline(current), encoding="utf-8", newline="\n"
        )
        print(f"[teensy-warnings] wrote {len(current)} warning(s) to {args.baseline}")
        return 0

    baseline = load_baseline(args.baseline)
    new = sorted(current - baseline)
    if not new:
        print(f"[teensy-warnings] PASS - {len(current)} warning(s), none new "
              f"(baseline has {len(baseline)}; cold capture: all {audit.declared} "
              f"first-party translation unit(s) across {len(audit.envs)} "
              f"environment(s) compiled, {compiles} invocation(s)).")
        return 0

    item_prefix = "::error::" if args.github else "  - "
    print(f"[teensy-warnings] FAIL - {len(new)} new first-party warning(s) not in "
          f"the baseline:")
    for w in new:
        print(f"{item_prefix}{w}")
    print("If intentional, regenerate the baseline in this PR: "
          "python tools/teensy_warnings.py --build-log <log> --update-baseline")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
