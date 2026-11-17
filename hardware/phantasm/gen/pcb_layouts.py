"""Revision placements, labels and locked local routes."""
import connectivity
import sexp
from kicad_common import (uid, fmt, F, rotate_point)
from parts_rev13 import PARTS as REV13_PARTS
from pcb_model import GROUND_NET
from pcb_model import TERMINAL_LIBIDS
from pcb_model import GenerationContext
from revision_policy import DEFAULT_REVISION, RevisionLayout, revision_paths
from pcb_pack import mounting_reserve_rects
from pcb_model import silk_revision


TERMINAL_EDGE_PLACEMENTS = {
    "J1": (3.25, 18.9, 0),
    "J2": (48.0, 3.96, 0),
    "J3A": (48.0, 13.46, 0),
    "J3B": (48.0, 22.96, 0),
}


TERMINAL_EDGE_PLACEMENTS_1_3 = {
    "J1": (13.11, 3.43, 270),
    "J2": (22.23, 3.43, 270),
    "J3A": (33.89, 3.43, 270),
    "J3B": (45.55, 3.43, 270),
}


TERMINAL_FIXED_1_3 = {
    "U_MCU": (27.15, 20.2, 0),
    "C_IN": (2.55, 10.75, 0),
    "C_DEC1": (11.5, 9.8, 0),
    "U1": (50.8, 19.0, 90),
    "C_DEC2": (47.19, 14.5, 0),
    "R_D1": (52.27, 23.7, 0),
    "R_D2": (48.3, 23.7, 0),
    "JP_ID0": (24.55, 8.93, 0),
    "JP_ID1": (28.35, 8.93, 0),
    "JP_ID2": (32.15, 8.93, 0),
    "JP_SHLD": (35.95, 8.93, 0),
    "JP_TERM": (39.35, 8.93, 90),
    "F1": (15.76, 8.93, 0),
    "Q_REV": (20.47, 8.93, 0),
    "D_SYNC": (43.03, 8.93, 0),
    "U_SYNC": (50.2, 10.15, 0),
    "C_DEC3": (55.5, 8.5, 0),
    "C_BULK3": (55.5, 10.25, 0),
    "R_A": (55.6, 12.8, 0),
    "R_B": (55.6, 14.5, 0),
    "R_TERM": (50.75, 14.25, 0),
    "FB": (56.82, 22.5, 90),
    "C_LF": (48.7, 30.0, 0),
    "R_LF": (44, 30.75, 0),
    "R_TX": (13, 30.95, 0),
    "R_DE": (20, 30.95, 0),
    "R_DATA_PD": (48.5, 26.5, 270),
    "R_CLK_PD": (46.5, 26.5, 270),
}


SYNC_FILTER_PLACEMENTS = {
    "R1": (15.7, 22.0, 0),
    "R2": (22.5, 22.0, 0),
    "C_SYNC": (19.15, 22.0, 0),
}


TERMINAL_FIXED = {
    "JP_ID0": (56.0, 8.5, 90),
    "JP_ID1": (56.0, 11.9, 90),
    "JP_ID2": (56.0, 15.3, 90),
    "JP_SHLD": (56.0, 19.0, 90),
    "U_MCU": (25.5, 11.7, 0),
    "C_IN": (20.5, 27.0, 0),
    "C_DEC1": (9.75, 1.35, 0),
    "U1": (36.6, 26.5, 180),
    "C_DEC2": (31.3, 30.0, 180),
    "R_D1": (42.2, 27.77, 180),
    "R_D2": (42.2, 23.96, 180),
    **SYNC_FILTER_PLACEMENTS,
}


REVISION_LAYOUTS = {
    "1.2": RevisionLayout(
        front_fixed_layout=False,
        edge_placements=TERMINAL_EDGE_PLACEMENTS,
        fixed=TERMINAL_FIXED,
        local_routes=(
            (("C_DEC1", "1"), ("U_MCU", "VIN"), ()),
            (("U_MCU", "3"), ("C_SYNC", "1"), ((19.15, 21.225),)),
            (("R1", "2"), ("C_SYNC", "1"), ()),
            (("R2", "1"), ("C_SYNC", "1"),
             ((20.7875, 21.2), (19.175, 21.2))),
        ),
        ground_routes=((("C_DEC1", "2"), (11.55, 1.35)),
                          (("C_SYNC", "2"), (19.925, 23.0)),
                          (("R2", "2"), (24.4125, 22.0))),
        jumper_marks=(("ID0", 54.3, 8.9, 90, "F.SilkS"),
                         ("ID1", 54.3, 11.9, 90, "F.SilkS"),
                         ("ID2", 54.3, 15.3, 90, "F.SilkS"),
                         ("SHLD", 54.3, 19.0, 90, "F.SilkS")),
        back_silk=[
            ("ID  ID0   ID1   ROLE(N4)", 7.0, 0.9),
            ("0   OPEN  OPEN  MASTER", 9.0, 0.9),
            ("1   GND   OPEN  A-SOUTH", 11.0, 0.9),
            ("2   OPEN  GND   B-NORTH", 13.0, 0.9),
            ("3   GND   GND   B-SOUTH", 15.0, 0.9),
            ("N8 ID2 OPEN=0-3 GND=4-7; M=OPEN; SHLD=M", 17.0, 0.8),
            ("CUT TEENSY VIN-VUSB", 21.3, 0.8),
            ("BOARD ID: ____", 23.5, 2.0),
            (silk_revision("1.2"), 29.5, 1.0),
        ],
        back_silk_x={"BOARD ID: ____": 25.5},
        sync_marks="SGH",
        terminal_libids=TERMINAL_LIBIDS[:2],
        terminal_by_ref={
            "J1": TERMINAL_LIBIDS[0],
            **dict.fromkeys(("J2", "J3A", "J3B"), TERMINAL_LIBIDS[1]),
        },
        part_properties={},
    ),
    "1.3": RevisionLayout(
        front_fixed_layout=True,
        edge_placements=TERMINAL_EDGE_PLACEMENTS_1_3,
        fixed=TERMINAL_FIXED_1_3,
        local_routes=((("C_DEC1", "1"), ("U_MCU", "VIN"), ()),),
        ground_routes=((("C_DEC1", "2"), (12.275, 10.95)),),
        jumper_marks=(("ID0", "JP_ID0", 8.93, 0, "B.SilkS"),
                         ("ID1", "JP_ID1", 8.93, 0, "B.SilkS"),
                         ("ID2", "JP_ID2", 8.93, 0, "B.SilkS"),
                         ("SHLD", "JP_SHLD", 8.93, 0, "B.SilkS")),
        back_silk=[
            ("TERM=CABLE ENDS", 7.0, 0.8),
            ("ID  ID0   ID1   ROLE(N4)", 17.0, 0.9),
            ("0   OPEN  OPEN  MASTER", 19.0, 0.9),
            ("1   GND   OPEN  A-SOUTH", 21.0, 0.9),
            ("2   OPEN  GND   B-NORTH", 23.0, 0.9),
            ("3   GND   GND   B-SOUTH", 25.0, 0.9),
            ("N8 ID2 OPEN=0-3 GND=4-7; M=OPEN; SHLD=M", 29.5, 0.8),
            ("CUT TEENSY VIN-VUSB", 10.6, 0.8),
            ("BOARD ID: ____", 14.8, 2.0),
            (silk_revision("1.3"), 31.0, 0.8),
            ("TERM", 8.93, 0.8),
        ],
        back_silk_x={"BOARD ID: ____": 20.0, "TERM": "JP_TERM"},
        sync_marks="ABGH",
        terminal_libids=TERMINAL_LIBIDS,
        terminal_by_ref={
            "J1": TERMINAL_LIBIDS[0],
            "J2": TERMINAL_LIBIDS[1],
            **dict.fromkeys(("J3A", "J3B"), TERMINAL_LIBIDS[2]),
        },
        part_properties=REV13_PARTS,
    ),
}


def fixed_placements(comps, ctx=None):
    """Verified fixed origins for the context's component records."""
    layout = (ctx or generation_context()).layout
    fixed = {ref: placement for ref, placement in layout.fixed.items()
             if ref in comps}
    fixed.update({ref: placement for ref, placement in layout.edge_placements.items()
                  if ref in comps and comps[ref][1] == layout.terminal_by_ref[ref]})
    return fixed


def terminal_edge_distances(layout):
    """Resolve one board-edge distance per terminal library footprint."""
    distances = {}
    for ref, (x, y, rot) in layout.edge_placements.items():
        libid = layout.terminal_by_ref[ref]
        distance = y if rot == 270 else x
        if libid in distances and distances[libid] != distance:
            raise ValueError(f"terminal library {libid} has conflicting board-edge distances")
        distances[libid] = distance
    return distances


def pin_silk_marks(ctx, fixed, length):
    """Resolve jumper and terminal-pin silkscreen labels."""
    layout = ctx.layout
    silk_marks = [(label, fixed[x][0] if isinstance(x, str) else x, y, angle, layer)
                  for label, x, y, angle, layer in layout.jumper_marks]
    for ref, pin_marks in (("J1", "+-"), ("J2", "DGC"),
                           ("J3A", layout.sync_marks), ("J3B", layout.sync_marks)):
        x, y, rot = fixed[ref]
        back = ref == "J1" or ctx.front_fixed_layout or any(
            x0 <= x + 4.5 <= x1 and y0 <= y + pin * 2.54 <= y1
            for pin in range(len(pin_marks))
            for x0, y0, x1, y1 in mounting_reserve_rects(length).values())
        offset = 2.2 if ctx.front_fixed_layout else (2 if ref == "J1" else (2.5 if back else 4.5))
        for pin, mark in enumerate(pin_marks):
            dx, dy = rotate_point((offset, pin * 2.54), rot)
            silk_marks.append((mark, x + dx, y + dy, 0, "B.SilkS" if back else "F.SilkS"))
    return silk_marks


def local_routes(ctx, footprints):
    """Pre-route the revision's locked local tracks and ground-via stubs on front copper."""
    by_ref = {connectivity.footprint_reference(footprint): footprint
              for footprint in footprints}

    def terminal(ref, number):
        footprint = by_ref[ref]
        if sexp.val(footprint, "locked") != ["yes"]:
            raise ValueError(f"local route requires locked {ref}")
        pad = next(pad for pad in F(footprint, "pad") if str(pad[1]) == number)
        x, y, angle = map(float, sexp.val(footprint, "at"))
        dx, dy = rotate_point(tuple(map(float, sexp.val(pad, "at")[:2])), angle)
        return (x + dx, y + dy), sexp.val(pad, "net")

    layout = ctx.layout
    lines = []
    for source, destination, bends in layout.local_routes:
        start, net = terminal(*source)
        end, other_net = terminal(*destination)
        if not net or net != other_net:
            raise ValueError(f"local route net mismatch: {source} / {destination}")
        points = (start, *bends, end)
        for a, b in zip(points, points[1:]):
            lines.append(f'\t(segment (start {fmt(a[0])} {fmt(a[1])}) '
                         f'(end {fmt(b[0])} {fmt(b[1])}) (width 0.2) '
                         f'(layer "F.Cu") (net {net[0]}) '
                         f'(uuid "{uid(site="pcb.py:local_routes")}") (locked yes))')
    for source, end in layout.ground_routes:
        start, net = terminal(*source)
        if str(net[-1]).lstrip("/") != GROUND_NET:
            raise ValueError(f"local ground via net mismatch: {source}")
        lines.append(f'\t(segment (start {fmt(start[0])} {fmt(start[1])}) '
                     f'(end {fmt(end[0])} {fmt(end[1])}) (width 0.3) '
                     f'(layer "F.Cu") (net {net[0]}) '
                     f'(uuid "{uid(site="pcb.py:local_routes")}") (locked yes))')
        lines.append(f'\t(via (at {fmt(end[0])} {fmt(end[1])}) '
                     '(size 0.6) (drill 0.3) (layers "F.Cu" "B.Cu") '
                     f'(net {net[0]}) (uuid "{uid(site="pcb.py:local_routes")}") (locked yes))')
    return lines


def generation_context(revision=DEFAULT_REVISION, unplaced=False, output_dir=None):
    """Resolve the revision layout and output destination."""
    if revision not in REVISION_LAYOUTS:
        raise ValueError(f"unsupported board revision: {revision!r}")
    out, _ = revision_paths(revision, output_dir)
    return GenerationContext(revision, unplaced, REVISION_LAYOUTS[revision], out)
