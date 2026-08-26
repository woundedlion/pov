"""Gate: walk the routed board's copper and verify every net is one island.

The named-net gate in check.py reads pad net ATTRIBUTES -- what a pad is meant
to sit on. Those survive deleting the tracks, vias or pour that realize them, so
a board can pass every net-label gate and still ship an open circuit. This one
unions the copper geometry instead: tracks (straight and arc), vias, pads and
zone fills that touch on a shared layer, per net. A net whose pads land in more than one island is an
open.

A rectangular land is modelled by its rotated rectangle; every other pad by its
circumscribed disc, so a marginal touch reads as connected. The error runs
towards passing: this gate answers "is the copper still there", not "does it
meet clearance" -- kicad-cli DRC owns that.
"""
import argparse
import math
import os
import sys

import sexp
from kicad_common import F, is_copper_pour, net_name

# Copper that lands this close counts as touching. KiCad snaps track ends to
# pad and via anchors, so the slack absorbs export rounding, not a real gap.
TOUCH_TOLERANCE = 0.001


def copper_layers(root):
    """Copper layer names in stack order, which is what a via span enumerates."""
    return [str(entry[1]) for entry in F(root, "layers")[0][1:]
            if str(entry[1]).endswith(".Cu")]


def net_id(node):
    net = sexp.val(node, "net")
    if not net or str(net[0]) == "0":
        return None
    return str(net[0]).lstrip("/")


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
    """A track segment, a via or a round pad: a thick line, degenerate to a
    disc."""

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
        if isinstance(other, Polygon):
            return other.touches(self)
        gap = _segment_distance(self.start, self.end, other.start, other.end)
        return gap <= self.radius + other.radius + TOUCH_TOLERANCE


class Polygon:
    """A copper outline: one layer of a pour fill, or a rectangular pad land.

    A pour voids around a through-hole pad and reaches it through thermal
    spokes, so the pad centre sits outside the fill while the annulus overlaps
    it. Overlap is therefore "a point inside" OR "the fill boundary runs through
    the item".
    """

    def __init__(self, polygon, layers):
        self.polygon, self.layers = polygon, layers

    def edges(self):
        return zip(self.polygon, self.polygon[1:] + self.polygon[:1])

    def touches(self, other):
        if self.layers.isdisjoint(other.layers):
            return False
        if isinstance(other, Polygon):
            if any(_point_in_polygon(p, self.polygon) for p in other.polygon):
                return True
            if any(_point_in_polygon(p, other.polygon) for p in self.polygon):
                return True
            # Two pours can overlap edge-on with no vertex of either inside the
            # other (a cross), which neither containment pass sees.
            return any(_segment_distance(a, b, c, d) <= TOUCH_TOLERANCE
                       for a, b in self.edges() for c, d in other.edges())
        if any(_point_in_polygon(p, self.polygon) for p in other.points()):
            return True
        reach = other.radius + TOUCH_TOLERANCE
        return any(_segment_distance(other.start, other.end, a, b) <= reach
                   for a, b in self.edges())


def _rectangle(centre, width, height, degrees):
    """Corners of a rectangle centred on `centre`, rotated into the board."""
    dx, dy = width / 2, height / 2
    return [(centre[0] + x, centre[1] + y)
            for x, y in (_rotate(corner, degrees) for corner in
                         ((-dx, -dy), (dx, -dy), (dx, dy), (-dx, dy)))]


def pad_copper(pad, origin, rotation, stack):
    """Pad copper: a rectangular land as its rotated rectangle, anything else as
    a disc around its centre, sized to circumscribe the land.

    Custom pads carry primitives outside the base size, so their reach is the
    furthest primitive vertex or circle rim.
    """
    placement = sexp.val(pad, "at")
    offset = _rotate(_xy(placement), rotation)
    centre = (origin[0] + offset[0], origin[1] + offset[1])
    width, height = _xy(sexp.val(pad, "size"))
    shape = str(pad[3])
    names = [str(value) for value in sexp.val(pad, "layers")]
    layers = (set(stack) if "*.Cu" in names
              else {name for name in names if name.endswith(".Cu")})
    if shape in {"rect", "roundrect", "trapezoid"}:
        # A pad angle is absolute: KiCad folds the footprint rotation into it.
        angle = float(placement[2]) if len(placement) > 2 else 0.0
        # A trapezoid's rect_delta pushes corners outside (size); box them in.
        skew = sexp.val(pad, "rect_delta")
        if skew:
            width, height = (width + abs(float(skew[1])),
                             height + abs(float(skew[0])))
        return Polygon(_rectangle(centre, width, height, angle), layers)
    radius = (max(width, height) / 2 if shape in {"circle", "oval"}
              else math.hypot(width, height) / 2)
    if shape == "custom":
        # KiCad nests the primitives one level down, under (primitives ...).
        blocks = F(pad, "primitives")
        for primitive in (blocks[0][1:] if blocks else []):
            if not isinstance(primitive, list) or not primitive:
                continue
            if str(primitive[0]) == "gr_circle":
                rim = _xy(sexp.val(primitive, "center"))
                edge = _xy(sexp.val(primitive, "end"))
                radius = max(radius, math.hypot(*rim) + math.dist(rim, edge))
            elif str(primitive[0]) == "gr_poly":
                for vertex in F(F(primitive, "pts")[0], "xy"):
                    radius = max(radius, math.hypot(*_xy(vertex[1:])))
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


def arc_capsules(arc):
    """A curved track as its two chords, start->mid and mid->end.

    The chords cut inside the arc by the sagitta, so a third item landing
    mid-span reads as further away than it is. Both endpoints are exact, which
    is where KiCad joins an arc to the rest of the net.
    """
    radius = float(sexp.val(arc, "width")[0]) / 2
    layers = {str(sexp.val(arc, "layer")[0])}
    start, mid, end = (_xy(sexp.val(arc, key))
                       for key in ("start", "mid", "end"))
    return [Capsule(start, mid, radius, layers),
            Capsule(mid, end, radius, layers)]


def footprint_reference(footprint):
    for child in footprint:
        if (isinstance(child, list) and len(child) > 2
                and child[0] == "property" and child[1] == "Reference"):
            return str(child[2])
    return "?"


def zone_layers(zone, stack):
    """Copper layers a zone occupies, from either KiCad spelling.

    A one-layer zone declares `(layer "In1.Cu")`; a multi-layer one declares
    `(layers "F.Cu" "B.Cu")`, or `(layers "*.Cu")` for the whole stack.
    """
    nodes = F(zone, "layer") + F(zone, "layers")
    if len(nodes) != 1:
        raise ValueError(
            f"zone on net {net_name(zone) or '?'} declares layer or layers "
            f"{len(nodes)} times, not once")
    names = [str(value) for value in nodes[0][1:]]
    layers = list(dict.fromkeys(
        entry for name in names
        for entry in (stack if name == "*.Cu" else (name,))))
    if not layers or any(name not in stack for name in layers):
        raise ValueError(
            f"zone on net {net_name(zone) or '?'} names a layer outside the "
            f"copper stack: {', '.join(names) or '(none)'}")
    return layers


def board_copper(root):
    """Returns copper, netted pads, and native net-id-to-name mappings."""
    stack = copper_layers(root)
    copper, pads = {}, {}
    # A `(net id "name")` declaration is itself the node net_id and net_name
    # look up as a child, so read it through a one-element parent.
    names = {}
    for declaration in F(root, "net"):
        key = net_id([declaration])
        if key is not None:
            names[key] = net_name([declaration])
    for footprint in F(root, "footprint"):
        reference = footprint_reference(footprint)
        side = sexp.val(footprint, "layer", [None])[0]
        if str(side) != "F.Cu":
            raise ValueError(
                f"footprint {reference} layer is {side or 'missing'}, not F.Cu; "
                "pads are placed by rotation alone, never mirrored")
        placement = sexp.val(footprint, "at")
        origin = _xy(placement)
        rotation = float(placement[2]) if len(placement) > 2 else 0.0
        for pad in F(footprint, "pad"):
            net = net_id(pad)
            if net is None:
                continue
            item = pad_copper(pad, origin, rotation, stack)
            copper.setdefault(net, []).append(item)
            pads.setdefault(net, []).append(((reference, str(pad[1])), item))
    for segment in F(root, "segment"):
        net = net_id(segment)
        if net is not None:
            copper.setdefault(net, []).append(segment_capsule(segment))
    for arc in F(root, "arc"):
        net = net_id(arc)
        if net is not None:
            copper.setdefault(net, []).extend(arc_capsules(arc))
    for via in F(root, "via"):
        net = net_id(via)
        if net is not None:
            copper.setdefault(net, []).append(via_capsule(via, stack))
    for zone in F(root, "zone"):
        net = net_id(zone)
        if net is None or not is_copper_pour(zone):
            continue
        declared = zone_layers(zone, stack)
        for filled in F(zone, "filled_polygon"):
            polygon = [_xy(vertex[1:]) for vertex in F(F(filled, "pts")[0], "xy")]
            if not polygon:
                continue
            # A multi-layer pour splits into one filled_polygon per layer, each
            # naming its own; a fill that names none covers the whole zone.
            own = sexp.val(filled, "layer")
            for layer in ([str(own[0])] if own else declared):
                copper.setdefault(net, []).append(Polygon(polygon, {layer}))
    return copper, pads, names


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
    copper, pads, names = board_copper(root)
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
            broken[names.get(net, net)] = sorted(
                sorted(group) for group in grouped.values())
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
            broken = opens(sexp.parse_one(handle.read()))
    except ValueError as error:
        print(f"{args.board}: {error}", file=sys.stderr)
        return 2
    print(report(broken))
    return 1 if broken else 0


if __name__ == "__main__":
    sys.exit(main())
