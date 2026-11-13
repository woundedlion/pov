"""Paths of committed Phantasm projects used by the test fixtures."""

from pathlib import Path


from constraints import ROUTED_REVISION
PROJECTS = Path(__file__).resolve().parents[2]


def committed_project(revision=ROUTED_REVISION):
    return PROJECTS / revision


def committed_board(revision=ROUTED_REVISION):
    return committed_project(revision) / "phantasm.kicad_pcb"
