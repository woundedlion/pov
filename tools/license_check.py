#!/usr/bin/env python3
"""Check that every tracked C/C++ source carries the license header LICENSE grants it.

LICENSE grants PolyForm Noncommercial to everything outside `effects/` and
reserves all rights over `effects/`, with a named exception list.

Scope is C/C++ sources; other files may carry the header but are ungated.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path, PurePosixPath

from git_files import git_nul_list, staged_blobs

NOTICE = "Required Notice: Copyright 2025 Gabriel Levy. All rights reserved."
POLYFORM = "Licensed under the PolyForm Noncommercial License 1.0.0"
RESERVED = "LICENSE: ALL RIGHTS RESERVED."
# Grant text of the third-party licenses LICENSE names. A path marked with one
# of these carries its upstream copyright line instead of NOTICE.
MIT_GRANT = "Permission is hereby granted"
MIT_TITLE = "MIT License"
THIRD_PARTY = frozenset({MIT_GRANT, MIT_TITLE})

SOURCE_SUFFIXES = (".h", ".hpp", ".c", ".cc", ".cpp", ".def", ".inl",
                   ".ino")

# Bytes of a file the header must appear within.
HEAD_BYTES = 600

# Paths LICENSE names as carrying their own terms, each with the marker its
# header must hold.
EXCEPTIONS = {
    "core/engine/effects_legacy.h": RESERVED,
    "core/math/projections.h": MIT_GRANT,
    "core/vendor/FastNoiseLite.h": MIT_TITLE,
}

UPSTREAM_HOLDERS = {
    "core/math/projections.h": ("Frank Warmerdam", "PROJ contributors", "Gerald I. Evenden",
                                 "Kristian Evers", "Toby C Wilkinson"),
    "core/vendor/FastNoiseLite.h": ("Jordan Peck", "Contributors"),
}

# Markers that cannot stand beside each other. The two MIT markers spell one
# grant, so they do not contradict each other.
CONTRADICTIONS = {
    POLYFORM: (RESERVED, MIT_GRANT, MIT_TITLE),
    RESERVED: (POLYFORM, MIT_GRANT, MIT_TITLE),
    MIT_GRANT: (POLYFORM, RESERVED),
    MIT_TITLE: (POLYFORM, RESERVED),
}

LICENSE_EXCEPTIONS_HEADING = (
    "Exceptions outside `effects/`, each carrying its own terms in its own header:"
)


def license_exception_paths(text: str) -> set[str]:
    """Paths named in LICENSE's exceptions block."""
    _, separator, remainder = text.partition(LICENSE_EXCEPTIONS_HEADING)
    if not separator:
        return set()

    paths = set()
    for line in remainder.splitlines():
        match = re.match(r"^- `([^`]+)`(?: |$)", line)
        if match:
            paths.add(match.group(1))
        elif paths and line and not line.startswith("  "):
            break
    return paths


def license_exception_issues(text: str) -> list[str]:
    """Differences between LICENSE and the executable exception map."""
    licensed = license_exception_paths(text)
    configured = set(EXCEPTIONS)
    issues = []
    missing = licensed - configured
    if missing:
        issues.append("LICENSE names exceptions absent from EXCEPTIONS: "
                      + ", ".join(sorted(missing)))
    extra = configured - licensed
    if extra:
        issues.append("EXCEPTIONS names paths absent from LICENSE's exceptions "
                      "block: " + ", ".join(sorted(extra)))
    return issues


class ToolingError(Exception):
    """A condition under which the checker can certify nothing."""


def tracked_sources(root: Path) -> list[str]:
    """Repo-relative paths of every tracked C/C++ source under root."""
    try:
        listed = git_nul_list(root, "ls-files", "-z")
    except (OSError, subprocess.SubprocessError, UnicodeError) as error:
        raise ToolingError(f"cannot run git: {error}") from error
    return sorted(p for p in listed
                  if p and PurePosixPath(p).suffix in SOURCE_SUFFIXES)


def expected_marker(path: str) -> str:
    """The header marker a path must carry."""
    return EXCEPTIONS.get(path) or (RESERVED if path.startswith("effects/") else POLYFORM)


def header_issue(path: str, head: str) -> str | None:
    """The complaint a file's header earns, or None where it satisfies LICENSE."""
    marker = expected_marker(path)
    # Case-folded: the license name is spelled both PolyForm and Polyform.
    folded = head.lower()
    if NOTICE.lower() not in folded and marker not in THIRD_PARTY:
        return f"no copyright notice in the first {HEAD_BYTES} bytes"
    if marker.lower() not in folded:
        return f"header does not carry {marker!r}, which LICENSE grants this path"
    for contradiction in CONTRADICTIONS[marker]:
        if contradiction.lower() in folded:
            return (f"header carries {contradiction!r} alongside {marker!r}; "
                    f"LICENSE gives this path only the latter")
    if marker in THIRD_PARTY:
        holders = UPSTREAM_HOLDERS[path]
        for holder in holders:
            if not any("copyright" in line.lower() and holder.lower() in line.lower()
                       for line in head.splitlines()):
                return f"missing upstream copyright notice for {holder}"
    return None


def check_staged(root: Path) -> int:
    """Check the staged headers of the paths read from stdin; return the exit code."""
    paths = [path for path in sys.stdin.read().splitlines()
             if PurePosixPath(path).suffix in SOURCE_SUFFIXES]
    try:
        blobs = staged_blobs(root, paths)
    except (OSError, subprocess.SubprocessError) as error:
        print(f"[license-check] tooling error: cannot run git: {error}", file=sys.stderr)
        return 2
    issues = [f"{path}:1: {issue}" for path, data in sorted(blobs.items())
              if (issue := header_issue(
                  path, data[:HEAD_BYTES].decode("utf-8", errors="replace")))]
    for issue in issues:
        print(issue)
    if issues:
        print(f"[license-check] FAIL - {len(issues)} issue(s)", file=sys.stderr)
        return 1
    print(f"[license-check] PASS - {len(blobs)} staged source(s)")
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."),
                        help="repository root to check")
    parser.add_argument("--staged", action="store_true",
                        help="check only the staged headers of the paths read from stdin")
    args = parser.parse_args(argv)

    if args.staged:
        return check_staged(args.root)

    try:
        sources = tracked_sources(args.root)
    except ToolingError as error:
        print(f"[license-check] tooling error: {error}", file=sys.stderr)
        return 2

    # Zero tracked sources would certify nothing.
    if not sources:
        print(f"[license-check] tooling error: no tracked C/C++ sources under "
              f"{args.root.resolve()}", file=sys.stderr)
        return 2

    try:
        license_text = (args.root / "LICENSE").read_text(encoding="utf-8")
    except OSError as error:
        print(f"[license-check] tooling error: LICENSE is unreadable: {error}",
              file=sys.stderr)
        return 2

    issues = license_exception_issues(license_text)
    stale = set(EXCEPTIONS) - set(sources)
    for path in sources:
        # A tracked source can be absent from a sparse or interrupted checkout.
        try:
            head = (args.root / path).read_bytes()[:HEAD_BYTES]
        except OSError as error:
            issues.append(f"{path}:1: tracked source is unreadable: {error}")
            continue
        issue = header_issue(path, head.decode("utf-8", errors="replace"))
        if issue:
            issues.append(f"{path}:1: {issue}")
    if stale:
        issues.append("EXCEPTIONS names paths with no tracked C/C++ source: "
                      + ", ".join(sorted(stale)))
    if issues:
        for issue in issues:
            print(issue)
        print(f"[license-check] FAIL - {len(issues)} issue(s)", file=sys.stderr)
        return 1
    print(f"[license-check] PASS - {len(sources)} tracked source(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
