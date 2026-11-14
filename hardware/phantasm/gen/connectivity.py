"""Gate: walk the routed board's copper and verify every net is one island.

Unions the copper geometry -- tracks (straight and arc), vias, pads and zone
fills that touch on a shared layer, per net -- rather than pad net attributes. A
net whose pads land in more than one island is an open.

Lands retain their copper outlines; curved track chords deviate by at most
TOUCH_TOLERANCE from the arc centreline.
This gate checks copper continuity; kicad-cli DRC owns clearance.
"""
import argparse
import math
import os
import sys

import sexp
from constraints import ROUTED_REVISION
from kicad_common import F, is_copper_pour, net_name, property_value, rotate_point

# Copper that lands this close counts as touching. KiCad snaps track ends to
# pad and via anchors, so the slack absorbs export rounding, not a real gap.
TOUCH_TOLERANCE = 0.001


def copper_layers(root):
    """Copper layer names in stack order, which is what a via span enumerates.

    Raises ValueError when the board declares no copper layers.
    """
    blocks = F(root, "layers")
    layers = [str(entry[1]) for entry in (blocks[0][1:] if blocks else [])
              if isinstance(entry, list) and len(entry) > 1
              and str(entry[1]).endswith(".Cu")]
    if not layers:
        raise ValueError("board declares no copper layers")
    return layers


def net_id(node):
    net = sexp.val(node, "net")
    if not net or str(net[0]) == "0":
        return None
    return str(net[0]) if len(net) > 1 else net_name(node) or None


def _xy(values):
    return float(values[0]), float(values[1])


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
        if isinstance(other, (Polygon, CopperUnion, Ring)):
            return other.touches(self)
        gap = _segment_distance(self.start, self.end, other.start, other.end)
        return gap <= self.radius + other.radius + TOUCH_TOLERANCE


class Polygon:
    """A filled polygon Minkowski-grown by a circular radius."""

    def __init__(self, polygon, layers, radius=0):
        self.polygon, self.layers, self.radius = polygon, layers, radius

    def edges(self):
        return zip(self.polygon, self.polygon[1:] + self.polygon[:1])

    def touches(self, other):
        if self.layers.isdisjoint(other.layers):
            return False
        if isinstance(other, (CopperUnion, Ring)):
            return other.touches(self)
        if isinstance(other, Polygon):
            if any(_point_in_polygon(p, self.polygon) for p in other.polygon):
                return True
            if any(_point_in_polygon(p, other.polygon) for p in self.polygon):
                return True
            # Two pours can overlap edge-on with no vertex of either inside the
            # other (a cross), which neither containment pass sees.
            return any(_segment_distance(a, b, c, d) <= self.radius + other.radius + TOUCH_TOLERANCE
                       for a, b in self.edges() for c, d in other.edges())
        if any(_point_in_polygon(p, self.polygon) for p in other.points()):
            return True
        reach = self.radius + other.radius + TOUCH_TOLERANCE
        return any(_segment_distance(other.start, other.end, a, b) <= reach
                   for a, b in self.edges())


class CopperUnion:
    """Copper occupied by the union of a custom pad's anchor and primitives."""

    def __init__(self, shapes, layers):
        """Store member shapes sharing the pad's layer set."""
        self.shapes, self.layers = shapes, layers

    def touches(self, other):
        """Whether any member touches the other copper on a shared layer."""
        if self.layers.isdisjoint(other.layers):
            return False
        return any(shape.touches(other) for shape in self.shapes)


class Ring:
    """An unfilled circular primitive with a finite stroke width."""

    def __init__(self, centre, radius, half_stroke, layers):
        """Store the circle centre, centreline radius and half stroke width."""
        self.centre, self.radius = centre, radius
        self.half_stroke, self.layers = half_stroke, layers

    def touches(self, other):
        """Intersect radial intervals of connected copper with the ring."""
        if self.layers.isdisjoint(other.layers):
            return False
        if isinstance(other, CopperUnion):
            return other.touches(self)
        if isinstance(other, Ring):
            distance = math.dist(self.centre, other.centre)
            gap = max(0, distance - self.radius - other.radius,
                      abs(self.radius - other.radius) - distance)
            return gap <= self.half_stroke + other.half_stroke + TOUCH_TOLERANCE
        if isinstance(other, Polygon):
            near = (0 if _point_in_polygon(self.centre, other.polygon) else
                    min(_point_segment_distance(self.centre, a, b) for a, b in other.edges()))
            far = max(math.dist(self.centre, p) for p in other.polygon)
        else:
            near = _point_segment_distance(self.centre, other.start, other.end)
            far = max(math.dist(self.centre, p) for p in (other.start, other.end))
        return (max(0, near - other.radius) <= self.radius + self.half_stroke + TOUCH_TOLERANCE
                and far + other.radius >= self.radius - self.half_stroke - TOUCH_TOLERANCE)


def arc_polyline(start, mid, end):
    """Preserve declared arc points and bound each chord's radial error."""
    bx, by = mid[0] - start[0], mid[1] - start[1]
    cx, cy = end[0] - start[0], end[1] - start[1]
    determinant = 2 * (bx * cy - by * cx)
    if abs(determinant) < 1e-12:
        return [start, mid, end]
    b_sq, c_sq = bx * bx + by * by, cx * cx + cy * cy
    ux = (cy * b_sq - by * c_sq) / determinant
    uy = (bx * c_sq - cx * b_sq) / determinant
    centre = (start[0] + ux, start[1] + uy)
    radius = math.hypot(ux, uy)
    angles = [math.atan2(y - centre[1], x - centre[0]) for x, y in (start, mid, end)]
    first, middle, last = angles
    direction = 1 if (middle - first) % math.tau <= (last - first) % math.tau else -1
    max_step = min(math.pi, 4 * math.asin(math.sqrt(min(1, TOUCH_TOLERANCE / (2 * radius)))))
    points = [start]
    for a, b, finish in ((first, middle, mid), (middle, last, end)):
        sweep = ((b - a) * direction) % math.tau
        count = max(1, math.ceil(sweep / max_step))
        for i in range(1, count):
            angle = a + direction * sweep * i / count
            points.append((centre[0] + radius * math.cos(angle),
                           centre[1] + radius * math.sin(angle)))
        points.append(finish)
    return points


def primitive_copper(primitive, layers, transform):
    """Copper occupied by a custom graphic, preserving fill and stroke geometry."""
    kind = str(primitive[0])
    half_stroke = float(sexp.val(primitive, "width", [0])[0]) / 2
    filled = str(sexp.val(primitive, "fill", ["none"])[0]) in {"yes", "solid"}
    if kind == "gr_circle":
        centre = transform(_xy(sexp.val(primitive, "center")))
        radius = math.dist(centre, transform(_xy(sexp.val(primitive, "end"))))
        return ([Capsule(centre, centre, radius + half_stroke, layers)] if filled
                else [Ring(centre, radius, half_stroke, layers)])
    if kind == "gr_arc":
        points = arc_polyline(*(_xy(sexp.val(primitive, key)) for key in ("start", "mid", "end")))
    elif kind == "gr_line":
        points = [_xy(sexp.val(primitive, key)) for key in ("start", "end")]
    elif kind == "gr_rect":
        (x0, y0), (x1, y1) = (_xy(sexp.val(primitive, key)) for key in ("start", "end"))
        points = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
    elif kind == "gr_poly":
        points = [_xy(vertex[1:]) for vertex in F(F(primitive, "pts")[0], "xy")]
    else:
        raise ValueError(f"unsupported custom pad primitive: {kind}")
    points = [transform(p) for p in points]
    if kind in {"gr_rect", "gr_poly"}:
        if filled:
            return [Polygon(points, layers, half_stroke)]
        points.append(points[0])
    return [Capsule(a, b, half_stroke, layers) for a, b in zip(points, points[1:])]


def _rectangle(centre, width, height, degrees):
    """Corners of a rectangle centred on `centre`, rotated into the board."""
    dx, dy = width / 2, height / 2
    return [(centre[0] + x, centre[1] + y)
            for x, y in (rotate_point(corner, degrees) for corner in
                         ((-dx, -dy), (dx, -dy), (dx, dy), (-dx, dy)))]


def copper_bounds(item):
    """Axis-aligned bounds of a land or a custom pad's primitive union."""
    if isinstance(item, CopperUnion):
        bounds = [copper_bounds(shape) for shape in item.shapes]
        return (min(b[0] for b in bounds), min(b[1] for b in bounds),
                max(b[2] for b in bounds), max(b[3] for b in bounds))
    if isinstance(item, Ring):
        points, radius = [item.centre], item.radius + item.half_stroke
    elif isinstance(item, Polygon):
        points, radius = item.polygon, item.radius
    else:
        points, radius = [item.start, item.end], item.radius
    xs, ys = zip(*points)
    return min(xs) - radius, min(ys) - radius, max(xs) + radius, max(ys) + radius


def pad_copper(pad, origin, rotation, stack):
    """Copper land in board coordinates, with absolute KiCad pad orientation."""
    placement = sexp.val(pad, "at")
    offset = rotate_point(_xy(placement), rotation)
    centre = (origin[0] + offset[0], origin[1] + offset[1])
    width, height = _xy(sexp.val(pad, "size"))
    angle = float(placement[2]) if len(placement) > 2 else 0.0
    shape = str(pad[3])
    names = [str(value) for value in sexp.val(pad, "layers")]
    layers = (set(stack) if "*.Cu" in names
              else {name for name in names if name.endswith(".Cu")})

    def transform(point):
        """Map pad-local geometry into board coordinates."""
        x, y = rotate_point(point, angle)
        return centre[0] + x, centre[1] + y

    if shape == "circle":
        return Capsule(centre, centre, width / 2, layers)
    if shape == "oval":
        radius = min(width, height) / 2
        delta = (width / 2 - radius, height / 2 - radius)
        return Capsule(transform((-delta[0], -delta[1])), transform(delta), radius, layers)
    if shape in {"rect", "roundrect"}:
        ratio = float(sexp.val(pad, "roundrect_rratio", [0.25])[0]) if shape == "roundrect" else 0
        radius = min(width, height) * ratio
        if not 0 <= ratio <= 0.5:
            raise ValueError(f"invalid roundrect radius ratio: {ratio}")
        return Polygon(_rectangle(centre, width - 2 * radius, height - 2 * radius, angle),
                       layers, radius)
    if shape == "trapezoid":
        dx, dy = width / 2, height / 2
        sx, sy = (float(v) / 2 for v in sexp.val(pad, "rect_delta", [0, 0]))
        corners = ((-dx + sy, -dy - sx), (dx - sy, -dy + sx),
                   (dx + sy, dy - sx), (-dx - sy, dy + sx))
        return Polygon([transform(point) for point in corners], layers)
    if shape == "custom":
        options = F(pad, "options")
        anchor = str(sexp.val(options[0], "anchor", ["circle"])[0]) if options else "circle"
        if anchor == "rect":
            shapes = [Polygon(_rectangle(centre, width, height, angle), layers)]
        elif anchor == "circle":
            shapes = [Capsule(centre, centre, width / 2, layers)]
        else:
            raise ValueError(f"unsupported custom pad anchor: {anchor}")
        blocks = F(pad, "primitives")
        for primitive in (blocks[0][1:] if blocks else []):
            shapes.extend(primitive_copper(primitive, layers, transform))
        return CopperUnion(shapes, layers)
    raise ValueError(f"unsupported pad shape: {shape}")


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
    """Subdivide a track arc with centreline sagitta at most TOUCH_TOLERANCE."""
    radius = float(sexp.val(arc, "width")[0]) / 2
    layers = {str(sexp.val(arc, "layer")[0])}
    points = arc_polyline(*(_xy(sexp.val(arc, key)) for key in ("start", "mid", "end")))
    return [Capsule(a, b, radius, layers) for a, b in zip(points, points[1:])]


def footprint_reference(footprint):
    return property_value(footprint, "Reference", "?")


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


def land_components(land):
    """Keep disconnected custom primitives in separate electrical islands."""
    if not isinstance(land, CopperUnion):
        return [land]
    roots = islands(land.shapes)
    groups = {}
    for index, shape in zip(roots, land.shapes):
        groups.setdefault(index, []).append(shape)
    result = []
    for shapes in groups.values():
        layers = set(land.layers)
        for shape in shapes:
            shape.layers = layers
        result.append(CopperUnion(shapes, layers))
    return result


def flag_enabled(node, name):
    """Read a KiCad boolean, including the legacy empty-node spelling."""
    fields = F(node, name)
    return bool(fields) and (len(fields[0]) == 1 or str(fields[0][1]) == "yes")


def pad_presence(pad, land, stack, nonzone_copper):
    """Retain end layers, saved zone flashes and direct non-zone connections."""
    declared = set(land.layers)
    if str(pad[2]) != "thru_hole" or not flag_enabled(pad, "remove_unused_layers"):
        return declared
    present = {str(layer) for layer in sexp.val(pad, "zone_layer_connections", [])}
    if flag_enabled(pad, "keep_end_layers"):
        present.update((stack[0], stack[-1]))
    for other in nonzone_copper:
        if other is not land and land.touches(other):
            present.update(declared.intersection(other.layers))
    return declared.intersection(present)


def board_copper(root):
    """Returns copper, netted pads, and native net-id-to-name mappings."""
    if not isinstance(root, list) or not root or str(root[0]) != "kicad_pcb":
        raise ValueError("document is not a kicad_pcb")
    stack = copper_layers(root)
    copper, pads = {}, {}
    pad_policies = []
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
            land = pad_copper(pad, origin, rotation, stack)
            components = land_components(land)
            for item in components:
                copper.setdefault(net, []).append(item)
                pads.setdefault(net, []).append(((reference, str(pad[1])), item))
                pad_policies.append((net, pad, item))
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
    presence = [(item, pad_presence(pad, item, stack, copper[net]))
                for net, pad, item in pad_policies]
    for item, layers in presence:
        item.layers.intersection_update(layers)
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
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))), ROUTED_REVISION,
        "phantasm.kicad_pcb")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("board", nargs="?", default=default)
    args = parser.parse_args(argv)
    try:
        with open(args.board, encoding="utf-8") as handle:
            broken = opens(sexp.parse_one(handle.read()))
    except (OSError, ValueError) as error:
        print(f"{args.board}: {error}", file=sys.stderr)
        return 2
    print(report(broken))
    return 1 if broken else 0


if __name__ == "__main__":
    sys.exit(main())
