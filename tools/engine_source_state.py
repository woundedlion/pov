#!/usr/bin/env python3
"""Report tracked working-tree edits against HEAD and untracked non-ignored files."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from git_files import git_nul_list


def changed_sources(root: Path) -> list[str]:
    return [*git_nul_list(root, "diff", "HEAD", "--name-only", "-z"),
            *git_nul_list(root, "ls-files", "--others", "--exclude-standard", "-z")]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    try:
        changed = changed_sources(args.root.resolve())
    except (OSError, UnicodeError, ValueError, subprocess.SubprocessError) as error:
        print(f"Cannot determine engine source state: {error}", file=sys.stderr)
        return 1
    if changed:
        print("\n".join(changed))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
