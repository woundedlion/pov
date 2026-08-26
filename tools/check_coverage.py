#!/usr/bin/env python3
"""Enforce line-coverage floors from llvm-cov JSON: aggregate and per-directory."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def line_coverage(document: object) -> float:
    """Return the aggregate line percentage from an llvm-cov export."""
    if not isinstance(document, dict):
        raise ValueError("coverage document is not an object")
    data = document.get("data")
    if not isinstance(data, list) or len(data) != 1:
        raise ValueError("coverage document must contain one data entry")
    try:
        percent = data[0]["totals"]["lines"]["percent"]
    except (KeyError, TypeError) as error:
        raise ValueError("coverage document has no aggregate line percentage") from error
    if not isinstance(percent, (int, float)) or not math.isfinite(percent):
        raise ValueError("aggregate line percentage is not finite")
    if not 0.0 <= percent <= 100.0:
        raise ValueError("aggregate line percentage is outside [0, 100]")
    return float(percent)


def directory_coverage(document: object, directory: str) -> float:
    """Return the line percentage over one repository-relative directory.

    A directory the report names no file under -- or names only files with no
    instrumented line -- is fatal, not vacuously covered: either would report
    100% and satisfy any floor, and the aggregate floor already survives a
    subsystem dropping out.
    """
    if not isinstance(document, dict):
        raise ValueError("coverage document is not an object")
    data = document.get("data")
    if not isinstance(data, list) or len(data) != 1:
        raise ValueError("coverage document must contain one data entry")
    files = data[0].get("files")
    if not isinstance(files, list):
        raise ValueError("coverage document has no per-file summaries")
    needle = f"/{directory.strip('/')}/"
    total = covered = 0
    matched = 0
    for entry in files:
        try:
            name = entry["filename"]
            lines = entry["summary"]["lines"]
            count, hit = lines["count"], lines["covered"]
        except (KeyError, TypeError) as error:
            raise ValueError("per-file summary is malformed") from error
        if needle not in "/" + str(name).replace("\\", "/").lstrip("/"):
            continue
        matched += 1
        total += int(count)
        covered += int(hit)
    if matched == 0:
        raise ValueError(f"no file under {directory} appears in the report")
    if total == 0:
        raise ValueError(
            f"the {matched} file(s) under {directory} report no instrumented line")
    return 100.0 * covered / total


def directory_floor(text: str) -> tuple[str, float]:
    """Parse a `<directory>=<percent>` floor."""
    directory, separator, percent = text.partition("=")
    if not separator or not directory.strip("/"):
        raise argparse.ArgumentTypeError(f"{text!r} is not <directory>=<percent>")
    try:
        floor = float(percent)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{percent!r} is not a percentage") from None
    if not 0.0 <= floor <= 100.0:
        raise argparse.ArgumentTypeError("a floor must be between 0 and 100")
    return directory.strip("/"), floor


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("report", type=Path)
    parser.add_argument("--min-lines", type=float, required=True)
    parser.add_argument("--min-directory", type=directory_floor, action="append",
                        default=[], metavar="DIR=PERCENT",
                        help="line floor for one directory; repeatable")
    args = parser.parse_args()
    if not 0.0 <= args.min_lines <= 100.0:
        parser.error("--min-lines must be between 0 and 100")

    try:
        document = json.loads(args.report.read_text(encoding="utf-8"))
        percent = line_coverage(document)
        measured = [(directory, directory_coverage(document, directory), floor)
                    for directory, floor in args.min_directory]
    except (OSError, json.JSONDecodeError, ValueError) as error:
        parser.error(str(error))

    print(f"line coverage: {percent:.2f}% (minimum {args.min_lines:.2f}%)")
    failed = percent < args.min_lines
    for directory, found, floor in measured:
        print(f"  {directory}: {found:.2f}% (minimum {floor:.2f}%)")
        failed = failed or found < floor
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
