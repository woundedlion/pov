#!/usr/bin/env python3
"""Single source for external build-tool versions shared by CI and just."""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path


PINS = {
    "doxygen-awesome": "568f56cde6ac78b6dfcc14acd380b2e745c301ea",
    "emsdk": "5.0.0",
    "node": "24.13.0",
    "platformio": "6.1.19",
}

ROOT = Path(__file__).resolve().parents[1]
CONSUMERS = {
    ROOT / ".github/workflows/ci.yml": (
        "build_pins.py --github-output",
        "build_pins.py platformio",
    ),
    ROOT / ".github/workflows/docs.yml": ("build_pins.py doxygen-awesome",),
    ROOT / "justfile": ("build_pins.py doxygen-awesome",),
}


def duplicates_pin(text: str, name: str, value: str) -> bool:
    """Return whether a dependency context contains its literal pinned value."""
    aliases = (name.lower(), name.lower().replace("-", "_"))
    value_pattern = re.compile(
        rf"(?<![0-9A-Za-z.]){re.escape(value)}(?![0-9A-Za-z.])"
    )
    lines = text.splitlines()
    for index, line in enumerate(lines):
        code = line.split("#", 1)[0]
        if not value_pattern.search(code):
            continue
        context = "\n".join(lines[max(0, index - 2) : index + 1]).lower()
        if any(alias in context for alias in aliases):
            return True
    return False


def check_consumers() -> int:
    errors: list[str] = []
    for path, references in CONSUMERS.items():
        text = path.read_text(encoding="utf-8")
        for reference in references:
            if reference not in text:
                errors.append(f"{path.relative_to(ROOT)}: missing {reference!r}")
    workflow_paths = set((ROOT / ".github/workflows").glob("*.yml"))
    workflow_paths.update((ROOT / ".github/workflows").glob("*.yaml"))
    duplicate_paths = workflow_paths | set(CONSUMERS)
    for path in sorted(duplicate_paths):
        text = path.read_text(encoding="utf-8")
        for name, value in PINS.items():
            if duplicates_pin(text, name, value):
                errors.append(
                    f"{path.relative_to(ROOT)}: duplicates {name} pin {value}"
                )
    if errors:
        for error in errors:
            print(error)
        return 1
    print("build pins are single-sourced")
    return 0


def write_github_output() -> None:
    output = os.environ.get("GITHUB_OUTPUT")
    if not output:
        raise SystemExit("GITHUB_OUTPUT is not set")
    with open(output, "a", encoding="utf-8", newline="\n") as stream:
        for name, value in PINS.items():
            stream.write(f"{name.replace('-', '_')}={value}\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("name", nargs="?", choices=sorted(PINS))
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--github-output", action="store_true")
    args = parser.parse_args()
    if args.check:
        return check_consumers()
    if args.github_output:
        write_github_output()
        return 0
    if args.name:
        print(PINS[args.name])
        return 0
    parser.error("specify a pin name, --check, or --github-output")


if __name__ == "__main__":
    raise SystemExit(main())
