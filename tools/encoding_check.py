#!/usr/bin/env python3
"""Reject invalid UTF-8 and reversible cp1252 mojibake in tracked text."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import subprocess
import sys

from git_files import git_nul_list, staged_blobs


def mojibake_lines(text: str) -> list[int]:
    """Return line numbers containing a reversible non-ASCII run."""
    lines = []
    for number, line in enumerate(text.splitlines(), 1):
        if line.isascii():
            continue
        projected = bytearray()
        for char in line:
            if ord(char) in {0x81, 0x8D, 0x8F, 0x90, 0x9D}:
                projected.append(ord(char))
            else:
                try:
                    projected.extend(char.encode("cp1252"))
                except UnicodeError:
                    projected.append(0)
        for match in re.finditer(
                rb"[\xc2-\xdf][\x80-\xbf]|[\xe0-\xef][\x80-\xbf]{2}|[\xf0-\xf4][\x80-\xbf]{3}",
                projected):
            try:
                match[0].decode("utf-8")
            except UnicodeError:
                continue
            lines.append(number)
            break
    return lines


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--staged", action="store_true",
                        help="check the staged content of the paths read from stdin")
    args = parser.parse_args()
    root = args.root.resolve()
    try:
        if args.staged:
            contents = staged_blobs(root, sys.stdin.read().splitlines())
        else:
            files = git_nul_list(root, "ls-files", "-z")
            if not files:
                print("[encoding-check] tooling error: no tracked files", file=sys.stderr)
                return 2
            contents = {name: (root / name).read_bytes()
                        for name in files if (root / name).is_file()}
    except (OSError, subprocess.SubprocessError, UnicodeError) as error:
        print(f"[encoding-check] tooling error: cannot read tracked files ({error})",
              file=sys.stderr)
        return 2
    failed = False
    for name, data in contents.items():
        if b"\0" in data:
            continue
        try:
            text = data.decode("utf-8")
        except UnicodeDecodeError:
            print(f"{name}:1: not valid UTF-8", file=sys.stderr)
            failed = True
            continue
        for line in mojibake_lines(text):
            print(f"{name}:{line}: reversible cp1252 mojibake", file=sys.stderr)
            failed = True
    if not failed:
        print("[encoding-check] PASS - tracked text uses valid UTF-8")
    return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())
