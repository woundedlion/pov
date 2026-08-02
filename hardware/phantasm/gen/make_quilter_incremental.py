"""Snapshot the corrected placed board for an incremental Quilter routing job."""

from __future__ import annotations

import argparse
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


def snapshot_digest(path: Path) -> str:
    data = path.read_bytes().replace(b"\r\n", b"\n")
    return hashlib.sha256(data).hexdigest()


def snapshot_paths(root: Path) -> list[Path]:
    paths = [Path(name) for name in FILES]
    paths.extend(
        path.relative_to(root)
        for path in sorted((root / "phantasm.pretty").rglob("*"))
        if path.is_file()
    )
    return paths


def verify_snapshot(root: Path = OUTPUT) -> None:
    manifest = root / "SHA256SUMS.txt"
    expected = {}
    for line in manifest.read_text(encoding="utf-8").splitlines():
        digest, name = line.split("  ", 1)
        expected[Path(name)] = digest

    actual_paths = snapshot_paths(root)
    if set(expected) != set(actual_paths):
        missing = sorted(str(path) for path in set(actual_paths) - set(expected))
        extra = sorted(str(path) for path in set(expected) - set(actual_paths))
        raise RuntimeError(f"snapshot manifest mismatch: missing={missing}, extra={extra}")

    for path in actual_paths:
        digest = snapshot_digest(root / path)
        if digest != expected[path]:
            raise RuntimeError(f"snapshot hash mismatch: {path}")


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


def make_snapshot() -> None:
    board = (PROJECT / "phantasm.kicad_pcb").read_text(encoding="utf-8")
    schematic = (PROJECT / "phantasm.kicad_sch").read_text(encoding="utf-8")
    require_incremental_input(board, schematic)

    OUTPUT.mkdir(parents=True, exist_ok=True)
    for name in FILES:
        shutil.copy2(PROJECT / name, OUTPUT / name)
    footprint_output = OUTPUT / "phantasm.pretty"
    if footprint_output.exists():
        shutil.rmtree(footprint_output)
    shutil.copytree(PROJECT / "phantasm.pretty", footprint_output)

    hashes = []
    for path in snapshot_paths(OUTPUT):
        digest = snapshot_digest(OUTPUT / path)
        hashes.append(f"{digest}  {path.as_posix()}")
    (OUTPUT / "SHA256SUMS.txt").write_text("\n".join(hashes) + "\n", encoding="utf-8")
    verify_snapshot()
    print(f"wrote protected incremental Quilter project: {OUTPUT}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.verify:
        verify_snapshot()
        print(f"verified protected incremental Quilter project: {OUTPUT}")
    else:
        make_snapshot()


if __name__ == "__main__":
    main()
