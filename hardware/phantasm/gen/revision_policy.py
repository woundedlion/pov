"""Revision paths and typed schematic/placement policy data."""

from dataclasses import dataclass, field
import os

from parts_rev13 import PARTS as REV13_PARTS


GENERATION_REVISIONS = ("1.2", "1.3")
DEFAULT_REVISION = "1.2"
DEFAULT_OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                           DEFAULT_REVISION)
DEFAULT_SCH = os.path.join(DEFAULT_OUT, "phantasm.kicad_sch")


def revision_paths(revision, output_dir=None, *, default_out=DEFAULT_OUT,
                   default_sch=DEFAULT_SCH):
    """Resolve a generation revision with an optional current-project override."""
    if revision not in GENERATION_REVISIONS:
        raise ValueError(f"unsupported board revision: {revision!r}")
    out = output_dir or (default_out if revision == DEFAULT_REVISION else os.path.join(
        os.path.dirname(default_out), revision))
    sch = default_sch if output_dir is None and revision == DEFAULT_REVISION else os.path.join(
        out, "phantasm.kicad_sch")
    return out, sch


@dataclass(frozen=True)
class RevisionLayout:
    """Fixed component locations, local wiring and placement annotations."""

    edge_placements: dict[str, tuple[float, float, float]]
    fixed: dict[str, tuple[float, float, float]]
    local_routes: tuple
    ground_routes: tuple
    jumper_marks: tuple
    back_silk: list[tuple[str, float, float]]
    back_silk_x: dict[str, float | str]
    sync_marks: str
    terminal_libids: tuple[str, ...]
    terminal_by_ref: dict[str, str]
    part_properties: dict[str, dict[str, str]]


@dataclass(frozen=True)
class CircuitPolicy:
    """Electrical and drawing choices for one generated revision."""

    differential: bool
    paper: str
    logic_budget: str
    enable_net: str
    enable_ref: str
    enable_power: str
    tx_idle_power: str
    data_y: float
    clock_y: float
    channels: tuple[str, str]
    part_properties: dict[str, dict[str, str]] = field(default_factory=dict)
    property_offsets: dict[str, dict[str, tuple[float, float]]] = field(default_factory=dict)
    unit_property_offsets: dict[tuple[str, int], dict[str, tuple[float, float]]] = field(default_factory=dict)


CIRCUIT_POLICIES = {
    "1.2": CircuitPolicy(False, "A3", "~0.15 A", "MASTER_EN", "R_MEN", "+3V3", "GND",
                          137.16, 167.64, ("DATA", "CLK")),
    "1.3": CircuitPolicy(True, "A2", "0.25 A", "SYNC_DE", "R_DE", "GND", "+3V3",
                          167.64, 137.16, ("CLK", "DATA"), REV13_PARTS,
                          {"U_MCU": {"Reference": (0, -45.72), "Value": (0, 45.72)},
                           "U_SYNC": {"Reference": (-5.08, -15.24), "Value": (8.89, -15.24)},
                           "D_SYNC": {"Reference": (0, -7.62), "Value": (15.24, 5.08)},
                           "U1": {"Value": (-12.7, 5.08)},
                           "Q_REV": {"Reference": (0, -10.16)}},
                          {("U1", 5): {"Reference": (-12.7, -2.54)},
                           ("U1", 3): {"Value": (-25.4, 5.08)},
                           ("U1", 4): {"Value": (-25.4, 5.08)}}),
}


def circuit_policy(revision):
    """Return the supported revision's schematic policy."""
    try:
        return CIRCUIT_POLICIES[revision]
    except KeyError as error:
        raise ValueError(f"unsupported board revision: {revision!r}") from error
