#!/usr/bin/env python3
"""Reject reductions in concrete relax-bake and death-harness coverage."""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import tempfile
from pathlib import Path


HARNESS_FLOOR = "MIN_RELAX_BAKES_VERIFIED"
FLOOR_RE = re.compile(rf"constexpr int ({HARNESS_FLOOR})\s*=\s*(\d+)\s*;")
GAP_TABLE_RE = re.compile(r"GUARD_GAP_ALLOW\[\]\s*=\s*\{(.*?)\};", re.S)
GAP_ROW_RE = re.compile(r'\{\s*"([^"]+)"\s*,\s*(\d+)\s*\}')
COMMENT_RE = re.compile(r"/\*.*?\*/|//[^\n]*", re.S)
ALLOW_RE = re.compile(r"^([A-Za-z0-9_.-]+)=(\d+|gone)->(\d+|gone)$")
WORKFLOW = ".github/workflows/ci.yml"
HARNESS_PATH = "tools/relax_bake_harness.cpp"
DEATH_PATH = "tests/test_death.h"
# Weakenings already on master, each exact and each expiring at ALLOWANCE_EXPIRY;
# DOMAIN_RATCHET_ALLOW_WEAKEN adds to them.
HISTORICAL_ALLOWANCES = """
guard_gap.canvas.h=29->30
guard_gap.carousel.h=3->4
guard_gap.raster.h=9->10
guard_gap.volume.h=7->8
guard_gap.memory.h=1->2
guard_gap.Profile.ino=3->4
guard_gap.chain_host.h=0->1
guard_gap.shader_host.h=0->13
guard_gap.kernels.h=0->3
guard_gap.param_host.h=0->26
guard_gap.preset_host.h=0->2
guard_gap.conway.h=30->33
guard_gap.engine_bindings.h=5->7
guard_gap.Holosphere.ino=0->1
"""
ALLOWANCE_EXPIRY = "713c53e423133e5cc07594c8b2f33362b706f7a5"
Transition = tuple[str, str, str]


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


def floors(path: Path) -> dict[str, int]:
    """The harness coverage floor, keyed by its name; a rename is fatal."""
    result = counts(FLOOR_RE.findall(source(path)), path, "coverage floor")
    if not result:
        raise SystemExit(
            f"{path}: coverage floor {HARNESS_FLOOR} is missing or renamed")
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


def weakening_allowances(value: str) -> set[Transition]:
    """Exact before/after transitions the allowance text approves."""
    result: set[Transition] = set()
    for raw in value.splitlines():
        line = raw.strip()
        if not line:
            continue
        match = ALLOW_RE.fullmatch(line)
        if match is None:
            raise SystemExit(
                f"invalid weakening allowance entry: {line!r}; "
                "expected name=before->after"
            )
        entry = match.groups()
        if entry in result:
            raise SystemExit(f"weakening allowance is repeated: {line}")
        result.add(entry)
    return result


def compare_files(previous_harness: Path, current_harness: Path,
                  previous_death: Path, current_death: Path,
                  previous_ref: str,
                  allow: set[Transition]) -> tuple[int, set[Transition]]:
    """Compare one commit edge and return its status and used allowances."""
    current = floors(current_harness)
    gaps_current = death_pins(current_death)
    if not gaps_current:
        raise SystemExit(f"no GUARD_GAP_ALLOW rows parsed from {current_death}")

    previous = floors(previous_harness)
    gaps_previous = death_pins(previous_death)
    weakened: dict[str, tuple[str, str]] = {
        key: (str(value), str(current.get(key, "gone")))
        for key, value in previous.items()
        if key not in current or current[key] < value
    }
    if gaps_previous:
        weakened.update(
            {
                key: (str(gaps_previous.get(key, 0)), str(value))
                for key, value in gaps_current.items()
                if value > gaps_previous.get(key, 0)
            }
        )
    else:
        raise SystemExit(
            f"no GUARD_GAP_ALLOW rows parsed from {previous_death} "
            f"at {previous_ref}")

    transitions = {
        (key, before, after) for key, (before, after) in weakened.items()
    }
    approved = transitions & allow
    unapproved = sorted(transitions - allow)
    for key, before, after in unapproved:
        source_path = (
            current_death if key.startswith("guard_gap.") else current_harness
        )
        print(
            f"::error file={source_path}::{key} weakened ({before} -> {after}) - "
            f"restore its coverage or list {key}={before}->{after} in "
            "HISTORICAL_ALLOWANCES"
        )
    print(
        f"domain ratchets: {len(current)} floors and {len(gaps_current)} "
        f"guard-gap rows, {len(approved)} approved weakening(s) against "
        f"{previous_ref}."
    )
    return (1 if unapproved else 0), approved


def git(repo: Path, *args: str, required: bool = True) -> str | None:
    """Run a read-only git query in `repo`."""
    result = subprocess.run(
        ["git", "-C", str(repo), *args], capture_output=True, text=True
    )
    if result.returncode == 0:
        return result.stdout
    if not required:
        return None
    raise SystemExit(result.stderr.strip() or f"git {' '.join(args)} failed")


def gated_sources(repo: Path, ref: str, root: Path,
                  tag: str) -> tuple[Path, Path] | None:
    """The gated harness and death files at @p ref, copied under @p root.

    None once a gated path is absent there, which is a renamed or deleted file
    rather than a weakening, and is annotated as such.
    """
    written = []
    for gated, name in ((HARNESS_PATH, "harness.cpp"), (DEATH_PATH, "death.h")):
        text = git(repo, "show", f"{ref}:{gated}", required=False)
        if text is None:
            print(f"::error file={gated}::{gated} does not exist at {ref} - a "
                  "ratchet-gated file was renamed or deleted without updating "
                  "the gate")
            return None
        local = root / f"{tag}-{name}"
        local.write_text(text, encoding="utf-8")
        written.append(local)
    return (written[0], written[1])


def check_git_range(repo: Path, previous_ref: str, current_ref: str,
                    allow: set[Transition], allow_through: str | None = None) -> int:
    """Check each first-parent edge whose commit contains the ratchet job."""
    if allow:
        if allow_through is None:
            raise SystemExit("range allowances require --allow-through COMMIT")
        allow_through = git(repo, "rev-parse", "--verify",
                            f"{allow_through}^{{commit}}").strip()
    used = set()
    in_window = 0
    commits = git(
        repo, "rev-list", "--first-parent", "--reverse",
        f"{previous_ref}..{current_ref}", "--", HARNESS_PATH, DEATH_PATH
    ).splitlines()
    checked = 0
    failed = False
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        for commit in commits:
            workflow = git(repo, "show", f"{commit}:{WORKFLOW}", required=False)
            if workflow is None or not re.search(
                    r"^  domain-ratchets:\s*$", workflow, re.M):
                continue
            parent = git(repo, "rev-parse", f"{commit}^1").strip()
            paths = [gated_sources(repo, ref, root, tag)
                     for tag, ref in (("previous", parent),
                                      ("current", commit))]
            if any(pair is None for pair in paths):
                return 1
            print(f"checking domain ratchets for {commit} against {parent}")
            within = bool(allow) and git(
                repo, "merge-base", "--is-ancestor", commit, allow_through,
                required=False) is not None
            in_window += within
            edge_allow = allow if within else set()
            status, edge_used = compare_files(
                paths[0][0], paths[1][0], paths[0][1], paths[1][1],
                parent, edge_allow
            )
            used.update(edge_used)
            checked += 1
            failed |= status != 0
    if checked == 0:
        print("::warning::no commits with the domain-ratchets job in range")
    # Only a range reaching into the allowance window can under-use one.
    if in_window:
        for key, before, after in sorted(allow - used):
            print("::warning::allowed transition was not exercised: "
                  f"{key}={before}->{after}")
    print(f"domain ratchet range: {checked} commit edge(s) checked")
    return 1 if failed else 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("previous_harness", type=Path, nargs="?")
    parser.add_argument("current_harness", type=Path, nargs="?")
    parser.add_argument("previous_death", type=Path, nargs="?")
    parser.add_argument("current_death", type=Path, nargs="?")
    parser.add_argument("--previous-ref")
    parser.add_argument("--git-range", nargs=2, metavar=("PREVIOUS", "CURRENT"))
    parser.add_argument("--allow-through")
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    args = parser.parse_args(argv)
    allow = weakening_allowances(
        os.environ.get("DOMAIN_RATCHET_ALLOW_WEAKEN", "")
    )
    if args.git_range is not None:
        if any((args.previous_harness, args.current_harness,
                args.previous_death, args.current_death, args.previous_ref)):
            parser.error("--git-range cannot be combined with file arguments")
        return check_git_range(
            args.repo, *args.git_range,
            allow | weakening_allowances(HISTORICAL_ALLOWANCES),
            args.allow_through or ALLOWANCE_EXPIRY,
        )
    paths = (
        args.previous_harness, args.current_harness,
        args.previous_death, args.current_death
    )
    if any(path is None for path in paths) or args.previous_ref is None:
        parser.error("four files and --previous-ref are required")
    status, used = compare_files(*paths, args.previous_ref, allow)
    for key, before, after in sorted(allow - used):
        print(
            "::warning::DOMAIN_RATCHET_ALLOW_WEAKEN transition was not "
            f"exercised: {key}={before}->{after}"
        )
    return status


if __name__ == "__main__":
    raise SystemExit(main())
