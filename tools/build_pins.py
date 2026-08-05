#!/usr/bin/env python3
"""Single source for external build-tool versions shared by CI and just.

PINS are injected (--github-output / a named lookup) and must appear in no
build file literally. INLINE_PINS cannot be injected -- they name an action
input, an apt package or a pip requirement resolved before this script could
run -- so --check pins them by consistency instead: every spelling of them in
INLINE_SCAN must derive from the value here.
"""

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

# Versions the build files must spell out literally: a setup-action input, an
# apt package name and a pip requirement are all resolved before this script
# could inject anything, and the interpreter pin bootstraps the script itself.
# They are single-sourced by CONSISTENCY instead of substitution -- --check
# asserts every occurrence equals the value here, so a partial bump fails.
INLINE_PINS = {
    "python": "3.12",
    "numpy": "2.4.3",
    "clang": "18",
    "clang-format": "18.1.8",
    "doxygen": "1.17.0",
    "doxygen-sha256":
        "75419ef4f446fc1c24ef12514b574e66e898ee6f527c6ae2ad84f91a905823c2",
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

# Files scanned for INLINE_PINS occurrences.
INLINE_SCAN = (
    ROOT / ".github/workflows/ci.yml",
    ROOT / ".github/workflows/docs.yml",
    ROOT / "justfile",
    ROOT / "scripts/generate_luts.py",
    ROOT / ".githooks/pre-commit",
)

# (pattern, pin name, expected form of the pin value). The pattern's single
# capture group is the version as that spelling writes it; `form` maps the pin
# value to that spelling, so every occurrence is derived from one string.
INLINE_USES = (
    (r"python-version:\s*'([^']*)'", "python", lambda v: v),
    (r"\bnumpy==([\w.]+)", "numpy", lambda v: v),
    (r"\b(?:clang\+\+|clang|llvm)-(\d+)\b", "clang", lambda v: v),
    (r"\bllvm-\w+-(\d+)\b", "clang", lambda v: v),
    (r"\bclang-format==([\w.]+)", "clang-format", lambda v: v),
    (r"\bclang-format-(\d+)\b", "clang-format", lambda v: v.split(".")[0]),
    (r"EXPECTED_CLANG_FORMAT_MAJOR = (\d+)", "clang-format",
     lambda v: v.split(".")[0]),
    (r"HS_CLANG_FORMAT_MAJOR=(\d+)", "clang-format", lambda v: v.split(".")[0]),
    (r"Install Doxygen ([\w.]+) ", "doxygen", lambda v: v),
    (r"/Release_(\w+)/", "doxygen", lambda v: v.replace(".", "_")),
    (r"\bdoxygen-([\w.]+)\.linux", "doxygen", lambda v: v),
    (r"\bdoxygen-([\w.]+)/bin", "doxygen", lambda v: v),
    (r"([0-9a-f]{64})  doxygen\.tar\.gz", "doxygen-sha256", lambda v: v),
)


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


def check_inline_pins() -> list[str]:
    """Return one error per occurrence of an INLINE_PINS value that disagrees
    with the pin, plus one per pin no scanned file spells at all."""
    errors: list[str] = []
    seen: dict[str, int] = {name: 0 for name in INLINE_PINS}
    for path in INLINE_SCAN:
        text = path.read_text(encoding="utf-8")
        for index, line in enumerate(text.splitlines(), 1):
            for pattern, name, form in INLINE_USES:
                want = form(INLINE_PINS[name])
                for found in re.findall(pattern, line):
                    seen[name] += 1
                    if found != want:
                        errors.append(
                            f"{path.relative_to(ROOT)}:{index}: {name} pinned to "
                            f"{want!r} in build_pins.py but written {found!r}"
                        )
    for name, count in seen.items():
        if count == 0:
            errors.append(f"{name} pin is spelled in no scanned file")
    return errors


def check_consumers() -> int:
    errors: list[str] = check_inline_pins()
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
    print(f"build pins are single-sourced ({len(PINS)} injected, "
          f"{len(INLINE_PINS)} inline)")
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
    parser.add_argument("name", nargs="?", choices=sorted(PINS | INLINE_PINS))
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--github-output", action="store_true")
    args = parser.parse_args()
    if args.check:
        return check_consumers()
    if args.github_output:
        write_github_output()
        return 0
    if args.name:
        print((PINS | INLINE_PINS)[args.name])
        return 0
    parser.error("specify a pin name, --check, or --github-output")


if __name__ == "__main__":
    raise SystemExit(main())
