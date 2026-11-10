"""Shared KiCad prerequisites, generators and board fixtures."""

import contextlib
import io
import os
import sys
from committed_boards import ROUTED_REVISION, committed_board
from pathlib import Path

GEN = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN))

import board  # noqa: E402
import pcb  # noqa: E402
import sexp  # noqa: E402
import check  # noqa: E402
from connectivity import footprint_reference  # noqa: E402
import kicad_common  # noqa: E402
from kicad_common import F, kicad_cli, rotate_point  # noqa: E402

STOCK_SYMBOLS = os.path.isdir(sexp.KICAD_SHARE)
STOCK_FOOTPRINTS = os.path.isdir(pcb.FP_DIR)


def pinned_kicad_cli_error():
    """None when a kicad-cli on the pin resolves, else why kicad_cli() exited."""
    try:
        kicad_cli()
    except SystemExit as error:
        return str(error)
    return None


def generation_blockers():
    """Every missing prerequisite for generating, as human-readable causes."""
    blockers = []
    if not STOCK_SYMBOLS:
        blockers.append(f"KiCad {sexp.KICAD_MAJOR} stock symbols not found; "
                        "set KICAD_SYMBOL_DIR")
    if not STOCK_FOOTPRINTS:
        blockers.append(f"KiCad {sexp.KICAD_MAJOR} stock footprints not found; "
                        "set KICAD_FOOTPRINT_DIR")
    cli_error = pinned_kicad_cli_error()
    if cli_error:
        blockers.append(cli_error)
    return blockers


GENERATION_BLOCKERS = generation_blockers()
GENERATES = not GENERATION_BLOCKERS
GENERATES_REASON = "; ".join(GENERATION_BLOCKERS)


def generate_schematic(out, revision=board.B.REVISION):
    """Generate a schematic into `out`; return its path."""
    with contextlib.redirect_stdout(io.StringIO()):
        board.main(force=True, revision=revision, output_dir=str(out))
    return os.path.join(out, "phantasm.kicad_sch")


def generate(out, unplaced=False, revision=board.B.REVISION):
    """Generate a schematic and board into `out`; return the board path."""
    generate_schematic(out, revision=revision)
    with contextlib.redirect_stdout(io.StringIO()):
        pcb.main(unplaced=unplaced, force=True, force_teensy_library=True,
                 revision=revision, output_dir=str(out))
    return os.path.join(out, pcb.UNPLACED_FILE if unplaced else pcb.DRAFT_FILE)


def read(path):
    return sexp.parse(Path(path).read_text(encoding="utf-8"))[0]


def footprints_by_ref(root):
    """Referenced board footprints, rejecting missing or duplicate identities."""
    footprints = {}
    for footprint in F(root, "footprint"):
        reference = footprint_reference(footprint)
        if not reference or reference == "?":
            raise AssertionError("footprint has no reference")
        if reference in footprints:
            raise AssertionError(f"duplicate footprint reference: {reference}")
        footprints[reference] = footprint
    return footprints


def footprint_pose(footprint):
    """Board position and rotation, defaulting an omitted rotation to zero."""
    at = sexp.val(footprint, "at")
    return float(at[0]), float(at[1]), float(at[2]) if len(at) > 2 else 0.0


def pad_position(footprint, number):
    """Board coordinates of a numbered footprint pad."""
    pad = next(p for p in F(footprint, "pad") if str(p[1]) == number)
    x, y, angle = footprint_pose(footprint)
    dx, dy = rotate_point(tuple(map(float, sexp.val(pad, "at")[:2])), angle)
    return x + dx, y + dy


def committed_board_nets(revision=ROUTED_REVISION):
    root = read(committed_board(revision))
    nets = {}
    for reference, footprint in footprints_by_ref(root).items():
        for pad in F(footprint, "pad"):
            name = kicad_common.net_name(pad)
            if name:
                nets.setdefault(name, set()).add(check.node_key(reference, str(pad[1])))
    return nets
