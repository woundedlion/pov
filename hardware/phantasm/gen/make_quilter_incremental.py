"""Snapshot the corrected placed board for an incremental Quilter routing job."""

from __future__ import annotations

import hashlib
import shutil
from pathlib import Path


PROJECT = Path(__file__).resolve().parent.parent
OUTPUT = PROJECT / "quilter_incremental"
FILES = (
    "phantasm.kicad_pcb",
    "phantasm.kicad_sch",
    "phantasm.kicad_pro",
    "phantasm.kicad_sym",
    "fp-lib-table",
    "sym-lib-table",
)


def require_incremental_input(board: str, schematic: str) -> None:
    expected_board = (
        '(net "/SYNC_PULLDOWN")',
        '(net "/MASTER_EN")',
        '(name "GND_IN1")',
        '(name "GND_IN2")',
    )
    expected_schematic = (
        '(label "SYNC_PULLDOWN"',
        '(label "MASTER_EN"',
    )
    missing = [item for item in expected_board if item not in board]
    missing += [item for item in expected_schematic if item not in schematic]
    if missing:
        raise RuntimeError("incremental Quilter input is missing: " + ", ".join(missing))

    if board.count('(net "/SYNC_PULLDOWN")') != 2:
        raise RuntimeError(
            "SYNC_PULLDOWN must have exactly two unrouted PCB pads before snapshotting"
        )


def main() -> None:
    board = (PROJECT / "phantasm.kicad_pcb").read_text(encoding="utf-8")
    schematic = (PROJECT / "phantasm.kicad_sch").read_text(encoding="utf-8")
    require_incremental_input(board, schematic)

    OUTPUT.mkdir(parents=True, exist_ok=True)
    for name in FILES:
        shutil.copy2(PROJECT / name, OUTPUT / name)
    shutil.copytree(
        PROJECT / "phantasm.pretty",
        OUTPUT / "phantasm.pretty",
        dirs_exist_ok=True,
    )

    hashes = []
    for name in FILES:
        digest = hashlib.sha256((OUTPUT / name).read_bytes()).hexdigest()
        hashes.append(f"{digest}  {name}")
    (OUTPUT / "SHA256SUMS.txt").write_text("\n".join(hashes) + "\n", encoding="utf-8")
    print(f"wrote protected incremental Quilter project: {OUTPUT}")


if __name__ == "__main__":
    main()
