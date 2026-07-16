#!/usr/bin/env python3
"""Run `pio run` for the given Teensy envs and append a combined memory table.

`just teensy-size` drives this wrapper instead of calling pio directly. It
streams the pio output unchanged, then re-prints the per-env teensy_size
FLASH/RAM1/RAM2 details — which otherwise scroll past one env at a time — as a
single side-by-side table after PlatformIO's own summary. The pio exit code is
propagated, so a size-gate or compile failure still fails the recipe.

Stdlib only; the teensy_size line parsing is teensy_gate.parse_teensy_size, so
the table and the gate can never disagree about what a line means.

Run:  python tools/teensy_size_table.py holosphere phantasm ...
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import teensy_gate  # noqa: E402

# PlatformIO's per-env banner, e.g. "Processing phantasm (board: teensy40; ...)".
_PROCESSING_RE = re.compile(r"^Processing (\S+) \(")

# Table rows: (region, key, label). Keys name either a teensy_size component,
# the region's "free for ..." figure, or "used" — the component sum the size
# gate budgets against (tools/teensy_budgets.json). RAM2 has a single component,
# so its sum row would be a duplicate and is omitted.
_ROWS = (
    ("flash", "code", "FLASH code"),
    ("flash", "data", "FLASH data"),
    ("flash", "headers", "FLASH headers"),
    ("flash", "used", "FLASH used (gate total)"),
    ("flash", "free", "FLASH free (files)"),
    ("ram1", "code", "RAM1 code (ITCM)"),
    ("ram1", "variables", "RAM1 variables (DTCM)"),
    ("ram1", "padding", "RAM1 padding"),
    ("ram1", "used", "RAM1 used (gate total)"),
    ("ram1", "free", "RAM1 free (stack)"),
    ("ram2", "variables", "RAM2 variables (OCRAM)"),
    ("ram2", "free", "RAM2 free (malloc/new)"),
)


def split_by_env(lines: list[str]) -> tuple[list[str], dict[str, str]]:
    """Split a full `pio run` log into per-env chunks, in build order."""
    order: list[str] = []
    chunks: dict[str, list[str]] = {}
    cur: str | None = None
    for line in lines:
        m = _PROCESSING_RE.match(line)
        if m:
            cur = m.group(1)
            if cur not in chunks:
                order.append(cur)
                chunks[cur] = []
        if cur is not None:
            chunks[cur].append(line)
    return order, {env: "\n".join(ls) for env, ls in chunks.items()}


def collect_sizes(lines: list[str]) -> tuple[list[str], dict[str, dict]]:
    """Per-env {region: {used, free, components}} parsed from a pio run log."""
    order, chunks = split_by_env(lines)
    return order, {env: teensy_gate.parse_teensy_size(text)
                   for env, text in chunks.items()}


def _cell(sizes: dict, region: str, key: str) -> str:
    measured = sizes.get(region)
    if measured is None:
        return "-"
    value = measured.get(key) if key in ("used", "free") \
        else measured.get("components", {}).get(key)
    return f"{value:,}" if value is not None else "-"


def render_table(order: list[str], sizes_by_env: dict[str, dict]) -> str:
    """Side-by-side memory table: one row per teensy_size figure, one column
    per env. An env that produced no teensy_size output (compile/link failure)
    renders as '-' cells rather than being dropped, so the column set always
    matches what was requested of pio."""
    rows = [[label] + [_cell(sizes_by_env.get(env, {}), region, key)
                       for env in order]
            for region, key, label in _ROWS]
    header = ["Memory (bytes)"] + order
    widths = [max(len(r[i]) for r in [header] + rows) for i in range(len(header))]

    def fmt(cells: list[str]) -> str:
        # Left-align the label column, right-align the numeric columns.
        first = cells[0].ljust(widths[0])
        rest = (c.rjust(w) for c, w in zip(cells[1:], widths[1:]))
        return "  ".join([first, *rest]).rstrip()

    rule = "  ".join("-" * w for w in widths)
    return "\n".join([fmt(header), rule] + [fmt(r) for r in rows])


def main(argv: list[str] | None = None) -> int:
    envs = sys.argv[1:] if argv is None else argv
    if not envs:
        print("usage: teensy_size_table.py <env> [<env> ...]", file=sys.stderr)
        return 2
    pio = shutil.which("pio") or shutil.which("platformio")
    if pio is None:
        print("error: pio not found on PATH (pip install platformio)",
              file=sys.stderr)
        return 2

    cmd = [pio, "run"]
    for env in envs:
        cmd += ["-e", env]
    # Merge stderr into stdout: teensy_size prints to stderr, and interleaving
    # the streams in real order is what lets split_by_env attribute its lines.
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True,
                            encoding="utf-8", errors="replace")
    lines: list[str] = []
    assert proc.stdout is not None
    for line in proc.stdout:
        sys.stdout.write(line)
        lines.append(line.rstrip("\n"))
    rc = proc.wait()

    order, sizes_by_env = collect_sizes(lines)
    if order:
        print()
        print(render_table(order, sizes_by_env))
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
