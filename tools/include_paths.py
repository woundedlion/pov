#!/usr/bin/env python3
"""Check repository-root spelling of first-party C++ includes."""

from __future__ import annotations

import argparse
from pathlib import Path
import re

from git_files import git_nul_list

SOURCE_SUFFIXES = {".h", ".hpp", ".cpp", ".cc", ".inl", ".ino"}
CONFIG_HEADER = "core/vendor/FastNoiseLite_config.h"
INCLUDE = re.compile(r'^[ \t]*#[ \t]*include[ \t]*"([^"\n]+)"', re.MULTILINE)
COMMENTS = re.compile(r"/\*.*?\*/|//[^\n]*", re.DOTALL)


def check_includes(root: Path, files: list[str]) -> list[str]:
    """Return diagnostics for noncanonical includes of tracked headers."""
    root = root.resolve()
    known = set(files)
    diagnostics = []
    for name in files:
        if Path(name).suffix not in SOURCE_SUFFIXES:
            continue
        if name.startswith("core/vendor/") and name != CONFIG_HEADER:
            continue
        path = root / name
        try:
            text = path.read_text(encoding="utf-8")
        except OSError as exc:
            diagnostics.append(f"{name}:1: cannot read source ({exc})")
            continue
        text = COMMENTS.sub(lambda match: "\n" * match.group().count("\n"), text)
        for match in INCLUDE.finditer(text):
            include = match.group(1)
            candidates = [path.parent / include, root / include,
                          root / "core" / include, root / "effects" / include,
                          root / "hardware" / include, root / "tests" / include]
            target = next((candidate.resolve() for candidate in candidates
                           if candidate.is_file()), None)
            if target is None or not target.is_relative_to(root):
                continue
            canonical = target.relative_to(root).as_posix()
            if canonical in known and include != canonical:
                line = text.count("\n", 0, match.start()) + 1
                diagnostics.append(f'{name}:{line}: include "{canonical}" '
                                   f'instead of "{include}"')
    return diagnostics


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path,
                        default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    files = git_nul_list(args.root, "ls-files", "-z", timeout=None)
    diagnostics = check_includes(args.root, files)
    if diagnostics:
        print("\n".join(diagnostics))
        return 1
    print("[include-paths] PASS - first-party includes use repository-root paths")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
