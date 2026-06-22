#!/usr/bin/env python3
"""Teensy firmware warning-hygiene ratchet (docs/teensy_ci_gate_spec.md §7.2).

Policy is a BASELINE RATCHET, not -Werror: capture the current first-party
warning set once, then fail only on warnings NOT in that baseline. New warnings
are blocked; fixing a baseline warning is fine.

Two properties the spec calls load-bearing:
  * SET-based, not line-ordered. PlatformIO builds in parallel (-j) so warning
    emission order is nondeterministic; an ordered diff would flap green/red on
    identical inputs. We compare normalized SETS.
  * Normalized to drop volatile bits (absolute path prefix, line:col) so the
    fingerprint is stable across unrelated edits, while the file + message + flag
    that identify the warning are kept.

First-party only: warnings outside core/ effects/ hardware/ targets/ (i.e.
FastLED and the Teensy core) are dropped — the independent backstop to the
`-isystem` plan, so the ratchet is robust even if the -isystem demotion proves
awkward under PlatformIO (§7.2).

This ratchet must run on a CACHE-DISABLED build: a cached TU emits no warnings,
so a warm build hides header-introduced warnings (§7.2, §15). Drive it from the
build log of a clean build, not the cached size build.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

FIRST_PARTY = ("core/", "effects/", "hardware/", "targets/")

# gcc: "<path>:<line>[:<col>]: warning: <message> [-Wflag]"
_WARNING_RE = re.compile(r"^(.*?):(\d+):(?:\d+:)?\s*warning:\s*(.*)$")


def _relativize(path: str) -> str | None:
    """Return a first-party repo-relative path, or None if not first-party."""
    p = path.replace("\\", "/")
    for fp in FIRST_PARTY:
        if p.startswith(fp):
            return p
        idx = p.find("/" + fp)
        if idx != -1:
            return p[idx + 1:]
    return None


def normalize(line: str) -> str | None:
    """Normalize one compiler line to a stable first-party warning key, or None."""
    m = _WARNING_RE.match(line.strip())
    if not m:
        return None
    rel = _relativize(m.group(1))
    if rel is None:
        return None
    return f"{rel}: warning: {m.group(3)}".rstrip()


def extract_warnings(build_log: str) -> set[str]:
    """The deduplicated, normalized, first-party warning set from a build log."""
    out: set[str] = set()
    for line in build_log.splitlines():
        key = normalize(line)
        if key is not None:
            out.add(key)
    return out


def load_baseline(path: str | Path) -> set[str]:
    """Read the committed baseline set (ignoring blank and #-comment lines)."""
    text = Path(path).read_text(encoding="utf-8") if Path(path).exists() else ""
    return {ln.strip() for ln in text.splitlines()
            if ln.strip() and not ln.lstrip().startswith("#")}


def render_baseline(warnings: set[str]) -> str:
    """Render a baseline file: header + sorted, deduplicated warning set."""
    header = [
        "# Teensy firmware first-party warning baseline (tools/teensy_warnings.py).",
        "# Sorted, deduplicated, normalized (path:line:col stripped). The ratchet",
        "# fails only on warnings NOT listed here. Regenerate with --update-baseline",
        "# in the SAME PR as the change that legitimately alters the set (§7.2).",
        "",
    ]
    return "\n".join(header + sorted(warnings)) + "\n"


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description="Teensy warning-baseline ratchet.")
    p.add_argument("--build-log", required=True, help="compiler output of a CLEAN build")
    p.add_argument("--baseline", default="tools/teensy_warning_baseline.txt")
    p.add_argument("--update-baseline", action="store_true",
                   help="rewrite the baseline from this build's warning set")
    p.add_argument("--github", action="store_true", help="emit ::error:: annotations")
    args = p.parse_args(argv)

    current = extract_warnings(Path(args.build_log).read_text(encoding="utf-8"))

    if args.update_baseline:
        Path(args.baseline).write_text(render_baseline(current), encoding="utf-8")
        print(f"[teensy-warnings] wrote {len(current)} warning(s) to {args.baseline}")
        return 0

    baseline = load_baseline(args.baseline)
    new = sorted(current - baseline)
    if not new:
        print(f"[teensy-warnings] PASS - {len(current)} warning(s), none new "
              f"(baseline has {len(baseline)}).")
        return 0

    prefix = "::error::" if args.github else "  - "
    print(f"[teensy-warnings] FAIL - {len(new)} new first-party warning(s) not in "
          f"the baseline:")
    for w in new:
        print(f"{prefix}{w}")
    print("If intentional, regenerate the baseline in this PR: "
          "python tools/teensy_warnings.py --build-log <log> --update-baseline")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
