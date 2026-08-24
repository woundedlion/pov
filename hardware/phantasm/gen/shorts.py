"""Union-find connectivity check over raw schematic geometry.
Reports electrical groups that merge two DIFFERENT named nets (labels / power)."""
import argparse
import math
import os
import sys
import sexp
from kicad_common import F

# Power flags mark a rail as driven for ERC; they adopt whatever net they sit on
# and never form one, so they are removed before counting nets in a group.
FLAG_NET = "PWR_FLAG"

DEFAULT_SCH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "phantasm.kicad_sch")

TOL = 0.02   # mm, perpendicular distance from the wire


def R(v):
    return (round(float(v[0]), 3), round(float(v[1]), 3))


def prop(sym, name):
    for d in sym:
        if isinstance(d, list) and d and d[0] == "property" and d[1] == name:
            return d[2]
    return None


def on_seg(p, a, b):
    (x, y), (x1, y1), (x2, y2) = p, a, b
    if not (min(x1, x2) - TOL <= x <= max(x1, x2) + TOL and
            min(y1, y2) - TOL <= y <= max(y1, y2) + TOL):
        return False
    cross = abs((x2 - x1) * (y - y1) - (y2 - y1) * (x - x1))
    length = math.hypot(x2 - x1, y2 - y1)
    return cross < TOL * length if length else True


def named_points(root):
    """point -> set of net names, from labels and power symbols.

    Names at one coordinate accumulate; a point carrying two different ones is a
    short on its own, with no wire between them.
    """
    named = {}
    for kind in ("label", "global_label", "hierarchical_label"):
        for c in F(root, kind):
            at = sexp.val(c, "at", [])
            if len(c) < 2 or len(at) < 2:
                continue
            named.setdefault(R(at), set()).add(c[1])
    for c in F(root, "symbol"):
        lib_id = sexp.val(c, "lib_id", [])
        at = sexp.val(c, "at", [])
        if not lib_id or len(at) < 2:
            continue
        lib = lib_id[0]
        if lib.startswith("power:") or "phantasm:+" in lib:
            value = prop(c, "Value")
            if value is not None:
                named.setdefault(R(at), set()).add(value)
    return named


def geometry(root):
    """Scan the top level; return (named points, wire spans, junction points).

    Malformed elements are skipped. Raises ValueError on a hierarchical sheet or
    a bus: the scan reads plain top-level wires only, so a sheet's contents and
    every bus member would be invisible.
    """
    unscanned = [kind for kind in ("sheet", "bus", "bus_entry") if F(root, kind)]
    if unscanned:
        raise ValueError(f"{', '.join(unscanned)} present: this checker scans "
                         "top-level wires only and would miss the connectivity "
                         "carried there")
    named = named_points(root)
    wires = []
    for c in F(root, "wire"):
        xy = [R((p[1], p[2])) for p in sexp.val(c, "pts", [])
              if isinstance(p, list) and len(p) >= 3 and p[0] == "xy"]
        if len(xy) >= 2:
            wires.append((xy[0], xy[1]))
    # junction dots explicitly connect crossing/T wires that only touch mid-span.
    junctions = [R(at) for c in F(root, "junction")
                 if len(at := sexp.val(c, "at", [])) >= 2]
    return named, wires, junctions


def analyze(root):
    """Union-find over the geometry; return (conflicts, bridges).

    conflicts: [(sorted real net names, [(net, point), ...])] per merged group.
    bridges:   human-readable description of every point landing mid-span.

    Raises ValueError if the schematic yields no named points or no wires -- an
    empty scan is a broken input, not a clean result.
    """
    named, wires, junctions = geometry(root)
    if not named or not wires:
        raise ValueError(f"nothing to analyze: {len(named)} named points, "
                         f"{len(wires)} wires")

    parent = {}

    def find(p):
        parent.setdefault(p, p)
        while parent[p] != p:
            parent[p] = parent[parent[p]]
            p = parent[p]
        return p

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    allpts = set(named) | set(junctions)
    for a, b in wires:
        allpts.add(a); allpts.add(b)
    for p in allpts:
        find(p)

    bridges = []
    for a, b in wires:
        union(a, b)
        # Any point on this span is on the wire: a named pin, another wire's
        # endpoint landing mid-span (a T-junction), or an explicit junction dot.
        for p in allpts:
            if p != a and p != b and on_seg(p, a, b):
                tag = (f"label/power {'+'.join(sorted(named[p]))}"
                       if p in named else "junction")
                union(p, a)
                bridges.append(f"{tag} on wire {a}-{b}")

    groups = {}
    for p, nets in named.items():
        groups.setdefault(find(p), set()).update(nets)
    conflicts = []
    for g, nets in groups.items():
        real = sorted(nets - {FLAG_NET})
        if len(real) > 1:
            members = sorted((net, p) for p, ns in named.items() if find(p) == g
                             for net in ns)
            conflicts.append((real, members))
    return conflicts, bridges


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("schematic", nargs="?", default=DEFAULT_SCH)
    path = parser.parse_args(argv).schematic
    try:
        with open(path, encoding="utf-8") as fh:
            conflicts, bridges = analyze(sexp.parse(fh.read())[0])
    except ValueError as e:
        print(f"{path}: {e}", file=sys.stderr)
        return 2
    print("=== suspect bridging edges (touching 2 named pts) ===")
    for reason in bridges:
        print("  ", reason)
    print(f"=== conflict groups (excl {FLAG_NET}) ===")
    for nets, members in conflicts:
        print("  MERGED:", nets)
        for net, p in members:
            print("     ", net, "@", p)
    print("conflict groups:", len(conflicts))
    return 1 if conflicts else 0


if __name__ == "__main__":
    sys.exit(main())
