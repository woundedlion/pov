"""Placed schematic pins and their wiring anchors for generator tests."""

import builder
import sexp
import shorts
from kicad_common import F


def dangling_pins(root):
    """[(ref, pin number, point)] for pins with no wiring or no-connect marker.

    Pin coordinates come from the schematic's own lib_symbols, so a stock
    symbol whose pin moved is measured where the generated file placed it.
    """
    libs = {node[1]: builder._index_unit_pins(node)
            for node in sexp.val(root, "lib_symbols", [])
            if isinstance(node, list) and node and node[0] == "symbol"}
    _, wires, junctions = shorts.geometry(root)
    anchors = set(junctions)
    anchors.update(shorts.R(sexp.val(node, "at"))
                   for kind in ("label", "global_label", "hierarchical_label")
                   for node in F(root, kind))
    anchors.update(shorts.R(tuple(map(float, sexp.val(node, "at"))))
                   for node in F(root, "no_connect"))
    for a, b in wires:
        anchors.add(a)
        anchors.add(b)
    placed = []
    touching = {}
    for inst in F(root, "symbol"):
        at = sexp.val(inst, "at")
        mirror = sexp.val(inst, "mirror")
        units = libs[sexp.val(inst, "lib_id")[0]]
        pins = dict(units.get(0, {}))
        pins.update(units.get(int(sexp.val(inst, "unit", [1])[0]), {}))
        for number, pin in pins.items():
            point = shorts.R(builder.transform(
                float(at[0]), float(at[1]),
                float(at[2]) if len(at) > 2 else 0.0,
                mirror[0] if mirror else None, pin["x"], pin["y"]))
            ref = next((p[2] for p in F(inst, "property")
                        if p[1] == "Reference"), None)
            placed.append((ref, number, point))
            touching[point] = touching.get(point, 0) + 1
    return [(ref, number, point) for ref, number, point in placed
            if point not in anchors and touching[point] < 2]


def pin_positions(root):
    """(reference, pin number) -> schematic point, from the file's lib_symbols."""
    libs = {node[1]: builder._index_unit_pins(node)
            for node in sexp.val(root, "lib_symbols")}
    positions = {}
    for inst in F(root, "symbol"):
        ref = next(p[2] for p in F(inst, "property") if p[1] == "Reference")
        units = libs[sexp.val(inst, "lib_id")[0]]
        pins = {**units.get(0, {}), **units.get(int(sexp.val(inst, "unit")[0]), {})}
        x, y, angle = map(float, sexp.val(inst, "at"))
        mirror = sexp.val(inst, "mirror", [None])[0]
        for number, pin in pins.items():
            positions[ref, number] = shorts.R(builder.transform(
                x, y, angle, mirror, pin["x"], pin["y"]))
    return positions
