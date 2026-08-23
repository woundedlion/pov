#!/usr/bin/env python3
"""Enforce a catastrophic line-coverage floor from llvm-cov JSON."""

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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("report", type=Path)
    parser.add_argument("--min-lines", type=float, required=True)
    args = parser.parse_args()
    if not 0.0 <= args.min_lines <= 100.0:
        parser.error("--min-lines must be between 0 and 100")

    try:
        document = json.loads(args.report.read_text(encoding="utf-8"))
        percent = line_coverage(document)
    except (OSError, json.JSONDecodeError, ValueError) as error:
        parser.error(str(error))

    print(f"line coverage: {percent:.2f}% (minimum {args.min_lines:.2f}%)")
    return 0 if percent >= args.min_lines else 1


if __name__ == "__main__":
    raise SystemExit(main())
