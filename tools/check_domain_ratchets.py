#!/usr/bin/env python3
"""Reject reductions in concrete relax-bake and death-harness coverage."""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path


HARNESS_FLOOR_RE = re.compile(
    r"constexpr int (MIN_RELAX_BAKES_VERIFIED)\s*=\s*(\d+)\s*;"
)
DEATH_FLOOR_RE = re.compile(
    r"constexpr int (MIN_COVERED_GUARD_SITES)\s*=\s*(\d+)\s*;"
)
GAP_TABLE_RE = re.compile(r"GUARD_GAP_ALLOW\[\]\s*=\s*\{(.*?)\};", re.S)
GAP_ROW_RE = re.compile(r'\{\s*"([^"]+)"\s*,\s*(\d+)\s*\}')


def floors(path: Path, pattern: re.Pattern[str]) -> dict[str, int]:
    text = path.read_text(encoding="utf-8")
    return {name: int(value) for name, value in pattern.findall(text)}


def death_pins(path: Path) -> tuple[dict[str, int], dict[str, int]]:
    text = path.read_text(encoding="utf-8")
    result = {name: int(value) for name, value in DEATH_FLOOR_RE.findall(text)}
    table = GAP_TABLE_RE.search(text)
    gaps = {
        f"guard_gap.{name}": int(gap)
        for name, gap in GAP_ROW_RE.findall(table.group(1) if table else "")
    }
    return result, gaps


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("previous_harness", type=Path)
    parser.add_argument("current_harness", type=Path)
    parser.add_argument("previous_death", type=Path)
    parser.add_argument("current_death", type=Path)
    parser.add_argument("--previous-ref", required=True)
    args = parser.parse_args()

    current = floors(args.current_harness, HARNESS_FLOOR_RE)
    death_current, gaps_current = death_pins(args.current_death)
    current.update(death_current)
    required = {"MIN_RELAX_BAKES_VERIFIED", "MIN_COVERED_GUARD_SITES"}
    if set(current) != required:
        raise SystemExit("domain coverage floor missing or renamed")
    if not gaps_current:
        raise SystemExit(f"no GUARD_GAP_ALLOW rows parsed from {args.current_death}")

    previous = floors(args.previous_harness, HARNESS_FLOOR_RE)
    death_previous, gaps_previous = death_pins(args.previous_death)
    previous.update(death_previous)
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
            if key.startswith(("guard_gap.", "MIN_COVERED"))
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
