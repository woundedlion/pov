#!/usr/bin/env python3
# Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
# Licensed under the PolyForm Noncommercial License 1.0.0
"""Check death-pin uniqueness and cross-directory breadcrumb aliases."""

import ast
from pathlib import Path
import re
import sys

STRING = r'"(?:[^"\\]|\\.)*"'
TOKEN = re.compile(STRING + r"|//[^\n]*|/\*.*?\*/|[(),]|[^(),\"/]+|/", re.S)


def literals(text):
    """Decode adjacent C++ string literals."""
    return "".join(ast.literal_eval(s) for s in re.findall(STRING, text))


def pins(text):
    """Read the literal fields of the death-case table."""
    pattern = (r'\{\s*"([^"\n]+)"\s*,\s*case_\w+\s*,\s*'
               r'"([^"\n]+)"\s*,\s*((?:' + STRING + r'\s*)+)\}')
    return [(m[1], m[2], literals(m[3])) for m in re.finditer(pattern, text)]


def guard_texts(text):
    """Read macro condition/message pairs with balanced parentheses."""
    text = re.sub(STRING + r"|//[^\n]*|/\*.*?\*/",
                  lambda m: m[0] if m[0].startswith('"') else " ", text,
                  flags=re.S)
    for match in re.finditer(r"HS_(?:AUDIT_)?CHECK\s*\(", text):
        depth = 1
        condition = []
        message = []
        field = condition
        for token in TOKEN.finditer(text, match.end()):
            value = token[0]
            if value == "(":
                depth += 1
            elif value == ")":
                depth -= 1
                if depth == 0:
                    break
            elif value == "," and depth == 1:
                if field is message:
                    break
                field = message
                continue
            field.append(value)
        yield "(" + "".join(condition).strip() + ") " + literals("".join(message))


def compact(text):
    return re.sub(r"\s+", "", text)


def main():
    root = Path(sys.argv[1])
    case_pins = pins((root / "tests/test_death.h").read_text(encoding="utf-8"))
    pinned_names = {Path(file).name for _, file, _ in case_pins}
    sources = {}
    for directory in ("core", "effects", "workbench", "hardware", "targets", "tools"):
        for path in (root / directory).rglob("*"):
            if path.is_file() and path.name in pinned_names:
                sources[path.relative_to(root).as_posix()] = list(
                    guard_texts(path.read_text(encoding="utf-8")))
    for name, file, text in case_pins:
        for other, candidates in sources.items():
            if other == file or Path(other).name != Path(file).name:
                continue
            if any(compact(guard).startswith(compact(text)) for guard in candidates):
                raise SystemExit(f"{name}: breadcrumb for {file} also matches {other}")
    guards = list(guard_texts((root / "core/animation/opleg.h").read_text(encoding="utf-8")))
    checked = 0
    for name, _file, text in case_pins:
        if not name.startswith("opleg_"):
            continue
        matches = sum(compact(guard).startswith(compact(text)) for guard in guards)
        if matches != 1:
            raise SystemExit(f"{name}: pin identifies {matches} OpLeg guard sites, expected 1")
        checked += 1
    if checked == 0:
        raise SystemExit("no OpLeg death pins found")
    print(f"death pins: {checked} OpLeg cases each identify one guard site")


if __name__ == "__main__":
    main()
