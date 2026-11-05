#!/usr/bin/env python3
"""Reject reversible cp1252 mojibake in tracked UTF-8 text."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import sys

from git_files import git_nul_list


def mojibake_lines(text: str) -> list[int]:
    """Return line numbers containing a reversible non-ASCII run."""
    lines = []
    for number, line in enumerate(text.splitlines(), 1):
        projected = bytearray()
        for char in line:
            if ord(char) in {0x81, 0x8D, 0x8F, 0x90, 0x9D}:
                projected.append(ord(char))
            else:
                try:
                    projected.extend(char.encode("cp1252"))
                except UnicodeError:
                    projected.append(0)
        for match in re.finditer(rb"[\xc2-\xf4][\x80-\xbf]+", projected):
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
    args = parser.parse_args()
    root = args.root.resolve()
    failed = False
    for name in git_nul_list(root, "ls-files", "-z"):
        path = root / name
        if not path.is_file():
            continue
        data = path.read_bytes()
        if b"\0" in data:
            continue
        try:
            text = data.decode("utf-8")
        except UnicodeDecodeError:
            continue
        for line in mojibake_lines(text):
            print(f"{name}:{line}: reversible cp1252 mojibake", file=sys.stderr)
            failed = True
    return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())
