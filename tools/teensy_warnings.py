#!/usr/bin/env python3
"""Gate cold Teensy captures against the zero first-party warning policy."""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path

FIRST_PARTY = ("core/", "effects/", "workbench/", "hardware/", "targets/")

# Library/toolchain roots: a path through any of these is third-party even when a
# nested dir reuses a first-party name (e.g. .platformio/lib/Foo/effects/x.h).
THIRD_PARTY = ("lib/", "libdeps/", ".platformio/", "packages/")

# gcc: "<path>:<line>[:<col>]: warning: <message> [-Wflag]"
_WARNING_RE = re.compile(r"^(.*?):(\d+):(?:\d+:)?\s*warning:\s*(.*)$")
_FILELESS_WARNING_RE = re.compile(
    r"^(<command-line>|(?:\S*[/\\])?(?:cc1plus|ld)(?:\.exe)?):\s*warning:\s*(.*)$")

# PlatformIO's non-verbose step line: "Compiling <object>" (absent under -v).
_PIO_STEP_RE = re.compile(r"^\s*Compiling\s+(\S+)")

# `-c` as a standalone flag; in the `-v` echo the source comes last, not after it.
_DASH_C_RE = re.compile(r"(?:^|\s)-c(?:\s|$)")

# A positional (not `-`-prefixed) source argument anywhere on that command line.
_SOURCE_RE = re.compile(r"(?:^|\s)(?!-)(\S+\.(?:cpp|cc|cxx|c|S))(?=\s|$)")

# PlatformIO's per-environment banner, which also prints the env's resolved
# options, `build_src_filter` among them:
#   Processing phantasm (board: teensy40; build_src_filter: -<*>, +<core/...>; ...)
_ENV_HEADER_RE = re.compile(r"^Processing\s+(\S+)\s+\((.*)\)\s*$")

# `key: value` pairs in that banner are `; `-separated, but a value may itself
# contain `, ` and spaces, so the value runs to the next `; <key>: `.
_SRC_FILTER_RE = re.compile(
    r"(?:^|;\s*)build_src_filter:\s*(.*?)(?=;\s*[\w.]+:\s|$)")

# One src-filter term: `+<path>` includes, `-<path>` excludes.
_SRC_FILTER_TERM_RE = re.compile(r"([+-])<([^>]*)>")

# An `[env:<name>]` declaration in platformio.ini; the shared `[env]` is not one.
_INI_ENV_RE = re.compile(r"^\s*\[env:([^\]\s]+)\]\s*$", re.MULTILINE)

# SCons object-cache hit:  Retrieved `.pio/build/x/src/core/memory.cpp.o' from cache
_CACHE_HIT_RE = re.compile(r"^\s*Retrieved\s+[`'\"](.+?)['\"]\s+from cache\s*$")

_GLOB_CHARS = "*?["


_FIRST_PARTY_DIRS = frozenset(fp.rstrip("/") for fp in FIRST_PARTY)
_THIRD_PARTY_DIRS = frozenset(tp.rstrip("/") for tp in THIRD_PARTY)


def _relativize(path: str) -> str | None:
    """Return the repo-root-relative path if first-party, else None.

    Anchored at the FIRST first-party segment, so the full nested path is kept
    (targets/Phantasm/effects/Foo.h, not effects/Foo.h).
    """
    segs = path.replace("\\", "/").split("/")
    if ".platformio" in segs or "libdeps" in segs:
        return None
    for i, seg in enumerate(segs):
        if seg in _FIRST_PARTY_DIRS:
            remainder = segs[i:]
            if _THIRD_PARTY_DIRS & set(remainder):
                return None
            return "/".join(remainder)
    return None


def normalize(line: str) -> str | None:
    """Normalize one compiler line to a stable first-party warning key, or None."""
    stripped = line.strip()
    m = _WARNING_RE.match(stripped)
    if not m:
        m = _FILELESS_WARNING_RE.match(stripped)
        if m:
            tool = m.group(1).replace("\\", "/").rsplit("/", 1)[-1].removesuffix(".exe")
            return f"{tool}: warning: {m.group(2)}".rstrip()
        return None
    rel = _relativize(m.group(1))
    if rel is None:
        return None
    return f"{rel}: warning: {m.group(3)}".rstrip()


def extract_warnings(build_log: str) -> set[str]:
    """The deduplicated, normalized, first-party warning set from a build log."""
    out: set[str] = set()
    origin = None
    for line in build_log.splitlines():
        if line.startswith("In file included from") or _PIO_STEP_RE.match(line):
            origin = None
        if "inlined from" in line:
            location = re.search(r" at (.*?):\d+(?::\d+)?[, :]?$", line)
            if location and origin is None:
                origin = _relativize(location.group(1))
        key = normalize(line)
        warning = _WARNING_RE.match(line.strip())
        if key is None and warning and origin:
            path = warning.group(1).replace("\\", "/")
            if re.search(r"(?:^|/)packages/toolchain-[^/]+/", path):
                key = f"{origin}: warning: {warning.group(3)}".rstrip()
        if key is not None:
            out.add(key)
        if warning:
            origin = None
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

    Zero means an empty warning set proves nothing.
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

    Each banner closes the previous section; lines before the first belong to
    no environment.
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

    A first-party glob is not countable from the log and raises.
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

    The object suffix is dropped, so `<tu>.cpp.o` and `<tu>.cpp` compare equal.
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


def declared_environments(ini_path: str | Path) -> tuple[str, ...]:
    """Every environment platformio.ini defines, in file order."""
    path = Path(ini_path)
    if not path.exists():
        raise CaptureError(
            f"'{ini_path}' does not exist, so the set of environments the build "
            f"was expected to cover is unknown; pass --env explicitly")
    envs = tuple(_INI_ENV_RE.findall(path.read_text(encoding="utf-8")))
    if not envs:
        raise CaptureError(
            f"'{ini_path}' declares no [env:<name>] section, so the set of "
            f"environments the build was expected to cover is unknown")
    return envs


def read_build_log(path: str | Path) -> str:
    """Read a captured build log, replacing undecodable bytes.

    A Windows `tee` capture interleaves cp1252 bytes; the matched warning
    fingerprints are ASCII, so substitution cannot alter the set.
    """
    return Path(path).read_text(encoding="utf-8", errors="replace")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Teensy first-party warning gate.")
    p.add_argument("--build-log", required=True, help="compiler output of a COLD build")
    p.add_argument("--env", action="append", metavar="NAME",
                   help="environment the capture was expected to cover; "
                        "repeatable. Defaults to every [env:<name>] in "
                        "--platformio-ini")
    p.add_argument("--platformio-ini", default="platformio.ini")
    p.add_argument("--github", action="store_true", help="emit ::error:: annotations")
    args = p.parse_args(argv)

    prefix = "::error::" if args.github else ""
    try:
        build_log = read_build_log(args.build_log)
    except OSError as exc:
        print(f"{prefix}[teensy-warnings] FAIL - cannot read {args.build_log} "
              f"({exc}): the warning gate has nothing to check.")
        return 1
    compiles = count_first_party_compiles(build_log)
    if compiles == 0:
        print(f"{prefix}[teensy-warnings] FAIL - no first-party compiler "
              f"invocation in {args.build_log}: the capture broke (or the build "
              f"was fully cached), so its empty warning set proves nothing. "
              f"Drive this from a cold `pio run -v 2>&1 | tee` log.")
        return 1

    try:
        audit = audit_capture(build_log)
        expected = (tuple(args.env) if args.env
                    else declared_environments(args.platformio_ini))
    except CaptureError as exc:
        print(f"{prefix}[teensy-warnings] FAIL - {exc} ({args.build_log}).")
        return 1
    if not audit.envs:
        print(f"{prefix}[teensy-warnings] FAIL - no `Processing <env> (...)` banner "
              f"in {args.build_log}: this is not a `pio run` capture, so the "
              f"expected first-party translation-unit set cannot be derived and a "
              f"partially cached build would read as green.")
        return 1
    absent = [e for e in expected if e not in {a.name for a in audit.envs}]
    if absent:
        print(f"{prefix}[teensy-warnings] FAIL - {len(absent)} of {len(expected)} "
              f"expected environment(s) have no section in {args.build_log}: "
              f"{', '.join(absent)}. The build stopped early (or never started "
              f"them), so their warnings are absent and the remaining "
              f"environment(s) cannot vouch for them.")
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

    if not current:
        print(f"[teensy-warnings] PASS - no first-party warnings "
              f"(cold capture: all {audit.declared} first-party translation "
              f"unit(s) across {len(audit.envs)} environment(s) compiled, "
              f"{compiles} invocation(s)).")
        return 0

    item_prefix = "::error::" if args.github else "  - "
    print(f"[teensy-warnings] FAIL - {len(current)} first-party warning(s):")
    for warning in sorted(current):
        print(f"{item_prefix}{warning}")
    print("The firmware policy is zero first-party warnings: fix them at the source.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
