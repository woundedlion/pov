"""Shared KiCad prerequisites, generated boards and committed net fixtures."""

import contextlib
import io
import os
import sys
from pathlib import Path

GEN = Path(__file__).resolve().parent.parent
REPO_ROOT = GEN.parents[2]
sys.path.insert(0, str(GEN))

import board  # noqa: E402
import pcb  # noqa: E402
import sexp  # noqa: E402
import check  # noqa: E402
import kicad_common  # noqa: E402
from kicad_common import kicad_cli  # noqa: E402

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


def committed_board_nets(revision="1.1"):
    path = REPO_ROOT / "hardware" / "phantasm" / revision / "phantasm.kicad_pcb"
    root = sexp.parse(path.read_text(encoding="utf-8"))[0]
    nets = {}
    for footprint in (
            node for node in root if isinstance(node, list) and node
            and node[0] == "footprint"):
        reference = next(
            (str(node[2]) for node in footprint
             if isinstance(node, list) and len(node) > 2
             and node[0] == "property" and node[1] == "Reference"),
            None,
        )
        for pad in (
                node for node in footprint if isinstance(node, list) and node
                and node[0] == "pad"):
            name = kicad_common.net_name(pad)
            if reference and name:
                nets.setdefault(name, set()).add(check.node_key(reference, str(pad[1])))
    return nets
