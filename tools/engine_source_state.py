#!/usr/bin/env python3
"""Report working-tree, index and untracked source edits against HEAD."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


def changed_sources(root: Path) -> list[str]:
    changed = subprocess.check_output(
        ["git", "-C", str(root), "diff", "HEAD", "--name-only", "-z"],
        timeout=30).decode("utf-8").split("\0")
    changed = [name for name in changed if name]
    untracked = subprocess.check_output(
        ["git", "-C", str(root), "ls-files", "--others", "--exclude-standard", "-z"],
        timeout=30).decode("utf-8").split("\0")
    changed.extend(name for name in untracked if name)
    return changed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    try:
        changed = changed_sources(args.root.resolve())
    except (OSError, UnicodeError, ValueError, subprocess.SubprocessError) as error:
        print(f"Cannot determine engine source state: {error}")
        return 1
    if changed:
        print("\n".join(changed))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
