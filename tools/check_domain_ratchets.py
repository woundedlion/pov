#!/usr/bin/env python3
"""Reject reductions in concrete relax-bake and death-harness coverage."""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path


HARNESS_FLOOR = "MIN_RELAX_BAKES_VERIFIED"
FLOOR_RES = {
    name: re.compile(rf"constexpr int ({name})\s*=\s*(\d+)\s*;")
    for name in (HARNESS_FLOOR,)
}
GAP_TABLE_RE = re.compile(r"GUARD_GAP_ALLOW\[\]\s*=\s*\{(.*?)\};", re.S)
GAP_ROW_RE = re.compile(r'\{\s*"([^"]+)"\s*,\s*(\d+)\s*\}')
COMMENT_RE = re.compile(r"/\*.*?\*/|//[^\n]*", re.S)


def source(path: Path) -> str:
    """The file's text with comments stripped, so no comment can be matched."""
    return COMMENT_RE.sub("", path.read_text(encoding="utf-8"))


def counts(matches, path: Path, what: str) -> dict[str, int]:
    """(name, count) matches as a mapping; a repeated name is fatal, not last-wins."""
    result: dict[str, int] = {}
    for key, value in matches:
        if key in result:
            raise SystemExit(f"{path}: {what} {key} is declared twice")
        result[key] = int(value)
    return result


def floors(path: Path, name: str) -> dict[str, int]:
    """The named coverage floor, keyed by its name; a rename is fatal."""
    result = counts(
        FLOOR_RES[name].findall(source(path)), path, "coverage floor"
    )
    if not result:
        raise SystemExit(f"{path}: coverage floor {name} is missing or renamed")
    return result


def death_pins(path: Path) -> dict[str, int]:
    """Each census file's approved count of guard sites no case pins.

    Every file is gated exactly and in both directions, and a file with no row
    must be fully pinned, so these rows subsume a whole-suite coverage total.
    """
    table = GAP_TABLE_RE.search(source(path))
    rows = counts(
        GAP_ROW_RE.findall(table.group(1) if table else ""),
        path,
        "guard-gap row",
    )
    return {f"guard_gap.{name}": gap for name, gap in rows.items()}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("previous_harness", type=Path)
    parser.add_argument("current_harness", type=Path)
    parser.add_argument("previous_death", type=Path)
    parser.add_argument("current_death", type=Path)
    parser.add_argument("--previous-ref", required=True)
    args = parser.parse_args()

    current = floors(args.current_harness, HARNESS_FLOOR)
    gaps_current = death_pins(args.current_death)
    if not gaps_current:
        raise SystemExit(f"no GUARD_GAP_ALLOW rows parsed from {args.current_death}")

    previous = floors(args.previous_harness, HARNESS_FLOOR)
    gaps_previous = death_pins(args.previous_death)
    allow = {
        line.strip()
        for line in os.environ.get("DOMAIN_RATCHET_ALLOW_WEAKEN", "").splitlines()
        if line.strip()
    }
    weakened = {
        key: f"{value} -> {current.get(key, 'gone')}"
        for key, value in previous.items()
        if key not in current or current[key] < value
    }
    if gaps_previous:
        weakened.update(
            {
                key: f"{gaps_previous.get(key, 0)} -> {value}"
                for key, value in gaps_current.items()
                if value > gaps_previous.get(key, 0)
            }
        )
    else:
        print(
            f"::warning::no GUARD_GAP_ALLOW table in {args.previous_ref}:"
            f"{args.current_death} - gap-widening check skipped"
        )

    unapproved = sorted(set(weakened) - allow)
    for key in unapproved:
        source = (
            args.current_death
            if key.startswith("guard_gap.")
            else args.current_harness
        )
        print(
            f"::error file={source}::{key} weakened ({weakened[key]}) - "
            f"restore its coverage or list {key} in "
            "DOMAIN_RATCHET_ALLOW_WEAKEN in this commit"
        )
    if unapproved:
        return 1

    for key in sorted(allow - set(weakened)):
        print(
            f'::warning::DOMAIN_RATCHET_ALLOW_WEAKEN is stale - "{key}" '
            "approves no weakening in this commit"
        )
    print(
        f"domain ratchets: {len(current)} floors and {len(gaps_current)} "
        f"guard-gap rows, {len(weakened)} approved weakening(s) against "
        f"{args.previous_ref}."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
