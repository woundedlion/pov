"""Snapshot the corrected placed board for an incremental Quilter routing job."""

from __future__ import annotations

import argparse
import hashlib
import shutil
from fnmatch import fnmatch
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

# Paths under the snapshot that carry no routing payload and so are never
# hashed: the manifest cannot cover itself, the README is prose, and the rest
# are KiCad's per-user lock/autosave/backup droppings plus editor scratch,
# all of which .gitignore already keeps out of the tree. Every other file the
# snapshot directory holds must appear in the manifest.
UNMANAGED = (
    "SHA256SUMS.txt",
    "README.md",
    "*.lck",
    "*.bak",
    "*-bak",
    "*.kicad_prl",
    "fp-info-cache",
    ".history",
)


def snapshot_digest(path: Path) -> str:
    data = path.read_bytes().replace(b"\r\n", b"\n")
    return hashlib.sha256(data).hexdigest()


def is_unmanaged(path: Path) -> bool:
    return any(fnmatch(part, pattern) for part in path.parts for pattern in UNMANAGED)


def snapshot_paths(root: Path) -> list[Path]:
    found = {
        path.relative_to(root)
        for path in root.rglob("*")
        if path.is_file() and not is_unmanaged(path.relative_to(root))
    }
    # Manifest order: the project files first, then everything else sorted.
    paths = [Path(name) for name in FILES if Path(name) in found]
    paths.extend(sorted(found.difference(paths)))
    return paths


def require_no_strays(root: Path) -> None:
    """Reject managed files a snapshot holds that this run does not rewrite.

    Hashing them would adopt them into the manifest, and verify_snapshot() would
    then demand them forever.
    """
    written = {Path(name) for name in FILES}
    strays = sorted(
        str(path) for path in snapshot_paths(root)
        if path not in written and path.parts[0] != "phantasm.pretty"
    )
    if strays:
        raise RuntimeError(
            "snapshot holds files this run does not write: " + ", ".join(strays)
        )


def verify_snapshot(root: Path = OUTPUT) -> None:
    manifest = root / "SHA256SUMS.txt"
    expected = {}
    for line in manifest.read_text(encoding="utf-8").splitlines():
        digest, _, name = line.partition("  ")
        if not name:
            raise RuntimeError(f"snapshot manifest line is malformed: {line!r}")
        expected[Path(name)] = digest

    actual_paths = snapshot_paths(root)
    if set(expected) != set(actual_paths):
        missing = sorted(str(path) for path in set(expected) - set(actual_paths))
        extra = sorted(str(path) for path in set(actual_paths) - set(expected))
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
    require_no_strays(OUTPUT)
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
    try:
        if args.verify:
            verify_snapshot()
            print(f"verified protected incremental Quilter project: {OUTPUT}")
        else:
            make_snapshot()
    except RuntimeError as exc:
        raise SystemExit(str(exc)) from None


if __name__ == "__main__":
    main()
