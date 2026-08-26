#!/usr/bin/env python3
"""Check that every tracked C/C++ source carries the license header LICENSE grants it.

LICENSE grants PolyForm Noncommercial to everything outside `effects/` and
reserves all rights over `effects/`, with a named exception list. A header that
disagrees with LICENSE is what a licensee acts on, so the two are gated against
each other here.

Scope is C/C++ sources: the Python and JavaScript tooling carries no license
header today, and this checker does not invent one for it.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path, PurePosixPath

NOTICE = "Required Notice: Copyright 2025 Gabriel Levy. All rights reserved."
POLYFORM = "Licensed under the PolyForm Noncommercial License 1.0.0"
RESERVED = "LICENSE: ALL RIGHTS RESERVED."
# Grant text of the third-party licenses LICENSE names. A path marked with one
# of these carries its upstream copyright line instead of NOTICE.
MIT_GRANT = "Permission is hereby granted"
MIT_TITLE = "MIT License"
THIRD_PARTY = frozenset({MIT_GRANT, MIT_TITLE})

SOURCE_SUFFIXES = (".h", ".hpp", ".c", ".cc", ".cpp", ".inl", ".ino")

# Bytes of a file the header must appear within. Generated files name their
# generator right below the notice, so this is a few lines rather than one.
HEAD_BYTES = 600

# Paths LICENSE names as carrying their own terms, each with the marker its
# header must hold. A directory prefix covers everything beneath it. Every entry
# names a marker: an exempt path asserts nothing, so stripping its banner would
# pass. core/vendor/FastNoiseLite_config.h is first-party and is not here.
EXCEPTIONS = {
    "core/engine/effects_legacy.h": RESERVED,
    "workbench/": RESERVED,
    "core/math/projections.h": MIT_GRANT,
    "core/vendor/FastNoiseLite.h": MIT_TITLE,
}

# Markers that cannot stand beside each other: a header granting PolyForm while
# reserving all rights tells a licensee two incompatible things, and so does a
# third-party header carrying either. The two MIT markers spell one grant, so
# they do not contradict each other.
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
        listed = subprocess.run(["git", "-C", str(root), "ls-files", "-z"],
                                capture_output=True, text=True)
    except OSError as error:
        raise ToolingError(f"cannot run git: {error}") from error
    if listed.returncode != 0:
        detail = listed.stderr.strip() or f"git exited {listed.returncode}"
        raise ToolingError(
            f"cannot list tracked sources under {root}: {detail}")
    out = listed.stdout
    return sorted(p for p in out.split("\0")
                  if p and PurePosixPath(p).suffix in SOURCE_SUFFIXES)


def expected_marker(path: str) -> str:
    """The header marker a path must carry."""
    for prefix, marker in EXCEPTIONS.items():
        if path == prefix or (prefix.endswith("/") and path.startswith(prefix)):
            return marker
    return RESERVED if path.startswith("effects/") else POLYFORM


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
    return None


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."),
                        help="repository root to check")
    args = parser.parse_args(argv)

    try:
        sources = tracked_sources(args.root)
    except ToolingError as error:
        print(f"[license-check] tooling error: {error}", file=sys.stderr)
        return 2

    try:
        license_text = (args.root / "LICENSE").read_text(encoding="utf-8")
    except OSError as error:
        print(f"[license-check] tooling error: LICENSE is unreadable: {error}",
              file=sys.stderr)
        return 2

    issues = license_exception_issues(license_text)
    stale = []
    for path in sources:
        # A tracked source can be absent from the working tree (an interrupted
        # checkout, a sparse one). That is an issue to report, not a traceback.
        try:
            head = (args.root / path).read_bytes()[:HEAD_BYTES]
        except OSError as error:
            issues.append(f"{path}:1: tracked source is unreadable: {error}")
            continue
        issue = header_issue(path, head.decode("utf-8", errors="replace"))
        if issue:
            issues.append(f"{path}:1: {issue}")
    for prefix in EXCEPTIONS:
        if not any(p == prefix or (prefix.endswith("/") and p.startswith(prefix))
                   for p in sources):
            stale.append(prefix)

    if stale:
        issues.append("EXCEPTIONS names paths with no tracked C/C++ source: "
                      + ", ".join(sorted(stale)))
    if issues:
        for issue in issues:
            print(issue)
        print(f"[license-check] FAIL - {len(issues)} issue(s)", file=sys.stderr)
        return 1
    # No tracked sources means the checker was pointed somewhere it cannot see
    # the repository; passing would certify nothing.
    if not sources:
        print(f"[license-check] tooling error: no tracked C/C++ sources under "
              f"{args.root.resolve()}", file=sys.stderr)
        return 2
    print(f"[license-check] PASS - {len(sources)} tracked source(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
