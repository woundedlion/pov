"""Placed courtyard bounds for board-generation tests."""

import math

import connectivity
import pcb
import sexp
from kicad_common import F, arc_extrema


def _graphic_points(node):
    """Local-frame corner points of one footprint graphic."""
    if str(node[0]) == "fp_circle":
        centre = sexp.val(node, "center")
        rim = sexp.val(node, "end")
        x, y = float(centre[0]), float(centre[1])
        radius = math.hypot(float(rim[0]) - x, float(rim[1]) - y)
        return [(x - radius, y - radius), (x + radius, y - radius),
                (x + radius, y + radius), (x - radius, y + radius)]
    points = [(float(value[0]), float(value[1]))
              for value in (sexp.val(node, key)
                            for key in ("start", "mid", "end", "center"))
              if value]
    if str(node[0]) == "fp_arc":
        points.extend(arc_extrema(*points[:3]))
    for vertex in (F(F(node, "pts")[0], "xy") if F(node, "pts") else []):
        points.append((float(vertex[1]), float(vertex[2])))
    if str(node[0]) == "fp_rect" and len(points) == 2:
        (x0, y0), (x1, y1) = points
        points = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
    return points


def courtyard_box(footprint):
    """Placed bounding box of a footprint's courtyard, or None if it draws none.

    The box circumscribes the courtyard, so it over-reports a clash between two
    interlocking outlines; nothing on this board is placed that tightly.
    """
    placement = sexp.val(footprint, "at")
    origin = (float(placement[0]), float(placement[1]))
    rotation = float(placement[2]) if len(placement) > 2 else 0.0
    xs, ys = [], []
    for child in footprint:
        if not (isinstance(child, list) and child):
            continue
        layer = sexp.val(child, "layer")
        if not layer or str(layer[0]) not in pcb.COURTYARD_LAYERS:
            continue
        for point in _graphic_points(child):
            x, y = connectivity._rotate(point, rotation)
            xs.append(origin[0] + x)
            ys.append(origin[1] + y)
    return (min(xs), min(ys), max(xs), max(ys)) if xs else None
