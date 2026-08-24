"""Gate: walk the routed board's copper and verify every net is one island.

The named-net gate in check.py reads pad net ATTRIBUTES -- what a pad is meant
to sit on. Those survive deleting the tracks, vias or pour that realize them, so
a board can pass every net-label gate and still ship an open circuit. This one
unions the copper geometry instead: tracks, vias, pads and zone fills that touch
on a shared layer, per net. A net whose pads land in more than one island is an
open.

Pads are modelled by their circumscribed disc, so a marginal touch reads as
connected. The error runs towards passing: this gate answers "is the copper
still there", not "does it meet clearance" -- kicad-cli DRC owns that.
"""
import argparse
import math
import os
import sys

import sexp
from kicad_common import F, is_copper_pour

# Copper that lands this close counts as touching. KiCad snaps track ends to
# pad and via anchors, so the slack absorbs export rounding, not a real gap.
TOUCH_TOLERANCE = 0.001


def copper_layers(root):
    """Copper layer names in stack order, which is what a via span enumerates."""
    return [str(entry[1]) for entry in F(root, "layers")[0][1:]
            if str(entry[1]).endswith(".Cu")]


def net_name(node):
    net = sexp.val(node, "net")
    return str(net[-1]).lstrip("/") if net else None


def _xy(values):
    return float(values[0]), float(values[1])


def _rotate(point, degrees):
    if not degrees:
        return point
    angle = math.radians(degrees)
    cos, sin = math.cos(angle), math.sin(angle)
    x, y = point
    # KiCad footprint rotation is counter-clockwise on a y-down canvas.
    return x * cos + y * sin, -x * sin + y * cos


def _point_segment_distance(point, start, end):
    (px, py), (ax, ay), (bx, by) = point, start, end
    dx, dy = bx - ax, by - ay
    span = dx * dx + dy * dy
    t = 0.0 if span == 0 else max(
        0.0, min(1.0, ((px - ax) * dx + (py - ay) * dy) / span))
    return math.hypot(px - (ax + dx * t), py - (ay + dy * t))


def _segment_distance(a1, a2, b1, b2):
    """Distance between two closed segments; zero when they cross."""
    def side(p, q, r):
        return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])

    d1, d2 = side(b1, b2, a1), side(b1, b2, a2)
    d3, d4 = side(a1, a2, b1), side(a1, a2, b2)
    if ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0)):
        return 0.0
    return min(_point_segment_distance(a1, b1, b2),
               _point_segment_distance(a2, b1, b2),
               _point_segment_distance(b1, a1, a2),
               _point_segment_distance(b2, a1, a2))


def _point_in_polygon(point, polygon):
    x, y = point
    inside = False
    previous = polygon[-1]
    for current in polygon:
        (xi, yi), (xj, yj) = current, previous
        if (yi > y) != (yj > y) and x < (xj - xi) * (y - yi) / (yj - yi) + xi:
            inside = not inside
        previous = current
    return inside


class Capsule:
    """A track segment, a via or a pad: a thick line, degenerate to a disc."""

    def __init__(self, start, end, radius, layers):
        self.start, self.end, self.radius = start, end, radius
        self.layers = layers

    def points(self):
        mid = ((self.start[0] + self.end[0]) / 2,
               (self.start[1] + self.end[1]) / 2)
        return (self.start, mid, self.end)

    def touches(self, other):
        if self.layers.isdisjoint(other.layers):
            return False
        if isinstance(other, Fill):
            return other.touches(self)
        gap = _segment_distance(self.start, self.end, other.start, other.end)
        return gap <= self.radius + other.radius + TOUCH_TOLERANCE


class Fill:
    """One filled polygon of a copper pour, on one layer.

    A pour voids around a through-hole pad and reaches it through thermal
    spokes, so the pad centre sits outside the fill while the annulus overlaps
    it. Overlap is therefore "a point inside" OR "the fill boundary runs through
    the item".
    """

    def __init__(self, polygon, layer):
        self.polygon, self.layers = polygon, {layer}

    def edges(self):
        return zip(self.polygon, self.polygon[1:] + self.polygon[:1])

    def touches(self, other):
        if self.layers.isdisjoint(other.layers):
            return False
        if isinstance(other, Fill):
            return any(_point_in_polygon(p, self.polygon) for p in other.polygon)
        if any(_point_in_polygon(p, self.polygon) for p in other.points()):
            return True
        reach = other.radius + TOUCH_TOLERANCE
        return any(_segment_distance(other.start, other.end, a, b) <= reach
                   for a, b in self.edges())


def pad_capsule(pad, origin, rotation, stack):
    """Pad copper as a disc around its centre, sized to circumscribe the land.

    Custom pads carry primitives outside the base size, so their reach is the
    furthest primitive vertex or circle rim.
    """
    offset = _rotate(_xy(sexp.val(pad, "at")), rotation)
    centre = (origin[0] + offset[0], origin[1] + offset[1])
    width, height = _xy(sexp.val(pad, "size"))
    radius = math.hypot(width, height) / 2
    if str(pad[3]) == "custom":
        for primitive in pad[1:]:
            if not isinstance(primitive, list) or not primitive:
                continue
            if str(primitive[0]) == "gr_circle":
                rim = _xy(sexp.val(primitive, "center"))
                edge = _xy(sexp.val(primitive, "end"))
                radius = max(radius, math.hypot(*rim) + math.dist(rim, edge))
            elif str(primitive[0]) == "gr_poly":
                for vertex in F(F(primitive, "pts")[0], "xy"):
                    radius = max(radius, math.hypot(*_xy(vertex[1:])))
    names = [str(value) for value in sexp.val(pad, "layers")]
    layers = (set(stack) if "*.Cu" in names
              else {name for name in names if name.endswith(".Cu")})
    return Capsule(centre, centre, radius, layers)


def via_capsule(via, stack):
    centre = _xy(sexp.val(via, "at"))
    span = [str(value) for value in sexp.val(via, "layers")]
    first, last = stack.index(span[0]), stack.index(span[-1])
    layers = set(stack[min(first, last):max(first, last) + 1])
    return Capsule(centre, centre, float(sexp.val(via, "size")[0]) / 2, layers)


def segment_capsule(segment):
    return Capsule(_xy(sexp.val(segment, "start")),
                   _xy(sexp.val(segment, "end")),
                   float(sexp.val(segment, "width")[0]) / 2,
                   {str(sexp.val(segment, "layer")[0])})


def footprint_reference(footprint):
    for child in footprint:
        if (isinstance(child, list) and len(child) > 2
                and child[0] == "property" and child[1] == "Reference"):
            return str(child[2])
    return "?"


def board_copper(root):
    """({net: [copper item]}, {net: [((ref, pad), item)]}) over the board."""
    stack = copper_layers(root)
    copper, pads = {}, {}
    for footprint in F(root, "footprint"):
        reference = footprint_reference(footprint)
        placement = sexp.val(footprint, "at")
        origin = _xy(placement)
        rotation = float(placement[2]) if len(placement) > 2 else 0.0
        for pad in F(footprint, "pad"):
            net = net_name(pad)
            if net is None:
                continue
            item = pad_capsule(pad, origin, rotation, stack)
            copper.setdefault(net, []).append(item)
            pads.setdefault(net, []).append(((reference, str(pad[1])), item))
    for segment in F(root, "segment"):
        net = net_name(segment)
        if net is not None:
            copper.setdefault(net, []).append(segment_capsule(segment))
    for via in F(root, "via"):
        net = net_name(via)
        if net is not None:
            copper.setdefault(net, []).append(via_capsule(via, stack))
    for zone in F(root, "zone"):
        net = net_name(zone)
        if net is None or not is_copper_pour(zone):
            continue
        layer = str(sexp.val(zone, "layer")[0])
        for filled in F(zone, "filled_polygon"):
            polygon = [_xy(vertex[1:]) for vertex in F(F(filled, "pts")[0], "xy")]
            if polygon:
                copper.setdefault(net, []).append(Fill(polygon, layer))
    return copper, pads


def islands(items):
    """Union-find over one net's copper; returns each item's island root."""
    parent = list(range(len(items)))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, first in enumerate(items):
        for j in range(i + 1, len(items)):
            if root(i) == root(j):
                continue
            if first.touches(items[j]):
                parent[root(i)] = root(j)
    return [root(i) for i in range(len(items))]


def opens(root):
    """{net: [[(ref, pad), ...], ...]} per net whose pads span several islands.

    Raises ValueError when the board contains no netted pads.
    """
    copper, pads = board_copper(root)
    if not pads:
        raise ValueError("nothing to analyze: no netted pads")
    broken = {}
    for net, net_pads in pads.items():
        if len(net_pads) < 2:
            continue
        items = copper[net]
        position = {id(item): index for index, item in enumerate(items)}
        roots = islands(items)
        grouped = {}
        for node, item in net_pads:
            grouped.setdefault(roots[position[id(item)]], []).append(node)
        if len(grouped) > 1:
            broken[net] = sorted(sorted(group) for group in grouped.values())
    return broken


def report(broken):
    lines = []
    for net, groups in sorted(broken.items()):
        lines.append(f"OPEN {net}")
        for group in groups:
            lines.append("   island "
                         + ", ".join(f"{ref}.{pad}" for ref, pad in group))
    lines.append("COPPER OPEN" if broken else "COPPER CONNECTED")
    return "\n".join(lines)


def main(argv=None):
    default = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "phantasm.kicad_pcb")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("board", nargs="?", default=default)
    args = parser.parse_args(argv)
    try:
        with open(args.board, encoding="utf-8") as handle:
            broken = opens(sexp.parse(handle.read())[0])
    except ValueError as error:
        print(f"{args.board}: {error}", file=sys.stderr)
        return 2
    print(report(broken))
    return 1 if broken else 0


if __name__ == "__main__":
    sys.exit(main())
