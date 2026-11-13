"""Paths of committed Phantasm projects used by the test fixtures."""

from pathlib import Path
import sys

GEN = Path(__file__).resolve().parents[1]
if str(GEN) not in sys.path:
    sys.path.insert(0, str(GEN))

from constraints import ROUTED_REVISION  # noqa: E402

PROJECTS = Path(__file__).resolve().parents[2]


def committed_project(revision=ROUTED_REVISION):
    return PROJECTS / revision


def committed_board(revision=ROUTED_REVISION):
    return committed_project(revision) / "phantasm.kicad_pcb"
