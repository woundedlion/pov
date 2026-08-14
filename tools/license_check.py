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
import subprocess
import sys
from pathlib import Path, PurePosixPath

NOTICE = "Required Notice: Copyright 2025 Gabriel Levy. All rights reserved."
POLYFORM = "Licensed under the PolyForm Noncommercial License 1.0.0"
RESERVED = "LICENSE: ALL RIGHTS RESERVED."

SOURCE_SUFFIXES = (".h", ".hpp", ".c", ".cpp", ".ino")

# Bytes of a file the header must appear within. Generated files name their
# generator right below the notice, so this is a few lines rather than one.
HEAD_BYTES = 600

# Paths LICENSE names as carrying their own terms, each with the marker its
# header must hold. A directory prefix covers everything beneath it.
EXCEPTIONS = {
    "core/engine/effects_legacy.h": RESERVED,
    "core/math/projections.h": "Permission is hereby granted",
    "core/vendor/": None,
}

# Markers that cannot stand beside each other: a header granting PolyForm while
# reserving all rights tells a licensee two incompatible things.
CONTRADICTIONS = {POLYFORM: RESERVED, RESERVED: POLYFORM}


def tracked_sources(root: Path) -> list[str]:
    """Repo-relative paths of every tracked C/C++ source under root."""
    out = subprocess.run(["git", "-C", str(root), "ls-files", "-z"],
                         capture_output=True, text=True, check=True).stdout
    return sorted(p for p in out.split("\0")
                  if p and PurePosixPath(p).suffix in SOURCE_SUFFIXES)


def expected_marker(path: str) -> str | None:
    """The header marker a path must carry, or None where it is exempt."""
    for prefix, marker in EXCEPTIONS.items():
        if path == prefix or (prefix.endswith("/") and path.startswith(prefix)):
            return marker
    return RESERVED if path.startswith("effects/") else POLYFORM


def header_issue(path: str, head: str) -> str | None:
    """The complaint a file's header earns, or None where it satisfies LICENSE."""
    marker = expected_marker(path)
    if marker is None:
        return None
    # Case-folded: the license name is spelled both PolyForm and Polyform.
    folded = head.lower()
    if NOTICE.lower() not in folded and marker != "Permission is hereby granted":
        return f"no copyright notice in the first {HEAD_BYTES} bytes"
    if marker.lower() not in folded:
        return f"header does not carry {marker!r}, which LICENSE grants this path"
    contradiction = CONTRADICTIONS.get(marker)
    if contradiction is not None and contradiction.lower() in folded:
        return (f"header carries {contradiction!r} alongside {marker!r}; "
                f"LICENSE gives this path only the latter")
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."),
                        help="repository root to check")
    args = parser.parse_args()

    sources = tracked_sources(args.root)
    issues = []
    stale = []
    for path in sources:
        head = (args.root / path).read_bytes()[:HEAD_BYTES].decode(
            "utf-8", errors="replace")
        issue = header_issue(path, head)
        if issue:
            issues.append(f"{path}:1: {issue}")
    for prefix in EXCEPTIONS:
        if not any(p == prefix or (prefix.endswith("/") and p.startswith(prefix))
                   for p in sources):
            stale.append(prefix)

    if stale:
        print("::warning::EXCEPTIONS in tools/license_check.py names paths that "
              f"no longer exist: {', '.join(sorted(stale))}")
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
