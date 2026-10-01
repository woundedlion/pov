#!/usr/bin/env python3
"""Omit the host process environment from a PlatformIO SCons dump."""

import ast
import re
import sys


def sanitize_envdump(text: str) -> str:
    matches = list(re.finditer(r"(?m)^  ['\"]ENV['\"]: *", text))
    if len(matches) != 1:
        raise ValueError("expected one SCons ENV field")
    start = matches[0].end()
    cursor = start
    wrapped = text.startswith("environ(", cursor)
    if wrapped:
        cursor += len("environ(")
    if text[cursor:cursor + 1] != "{":
        raise ValueError("expected an ENV dictionary")
    dictionary_start = cursor
    depth = 0
    quote = ""
    escaped = False
    for cursor in range(cursor, len(text)):
        char = text[cursor]
        if quote:
            if escaped:
                escaped = False
            elif char == "\\":
                escaped = True
            elif char == quote:
                quote = ""
        elif char in "'\"":
            quote = char
        elif char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                break
    else:
        raise ValueError("unterminated ENV dictionary")
    try:
        environment = ast.literal_eval(text[dictionary_start:cursor + 1])
    except (SyntaxError, ValueError) as exc:
        raise ValueError("invalid ENV dictionary") from exc
    if not isinstance(environment, dict) or any(
        not isinstance(key, str) or not isinstance(value, str)
        for key, value in environment.items()
    ):
        raise ValueError("expected string ENV entries")
    end = cursor + 1
    if wrapped:
        if text[end:end + 1] != ")":
            raise ValueError("unterminated environ wrapper")
        end += 1
    if text[end:].lstrip()[:1] not in (",", "}"):
        raise ValueError("invalid ENV field boundary")
    return text[:start] + "{}" + text[end:]


def main() -> int:
    try:
        text = sys.stdin.buffer.read().decode("utf-8", errors="surrogateescape")
        sanitized = sanitize_envdump(text)
    except ValueError as exc:
        print(f"profile envdump: {exc}", file=sys.stderr)
        return 1
    sys.stdout.buffer.write(sanitized.encode("utf-8", errors="surrogateescape"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
