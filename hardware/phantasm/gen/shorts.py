"""Union-find connectivity check over raw schematic geometry.
Reports electrical groups that merge two DIFFERENT named nets (labels / power)."""
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
    """point -> net name, from local labels and power symbols."""
    named = {}
    for c in F(root, "label"):
        named[R(sexp._val(c, "at"))] = c[1]
    for c in F(root, "symbol"):
        lib = sexp._val(c, "lib_id")[0]
        if lib.startswith("power:") or "phantasm:+" in lib:
            named[R(sexp._val(c, "at"))] = prop(c, "Value")
    return named


def analyze(root):
    """Union-find over the geometry; return (conflicts, bridges).

    conflicts: [(sorted real net names, [(net, point), ...])] per merged group.
    bridges:   human-readable description of every point landing mid-span.
    """
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

    named = named_points(root)
    wires = []
    for c in F(root, "wire"):
        xy = [R((p[1], p[2])) for p in sexp._val(c, "pts")
              if isinstance(p, list) and p[0] == "xy"]
        wires.append((xy[0], xy[1]))
    # junction dots explicitly connect crossing/T wires that only touch mid-span.
    junctions = [R(sexp._val(c, "at")) for c in F(root, "junction")
                 if sexp._val(c, "at")]

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
                tag = f"label/power {named[p]}" if p in named else "junction"
                union(p, a)
                bridges.append(f"{tag} on wire {a}-{b}")

    groups = {}
    for p, net in named.items():
        groups.setdefault(find(p), set()).add(net)
    conflicts = []
    for g, nets in groups.items():
        real = sorted(nets - {FLAG_NET})
        if len(real) > 1:
            members = sorted((net, p) for p, net in named.items() if find(p) == g)
            conflicts.append((real, members))
    return conflicts, bridges


def main(argv):
    path = argv[1] if len(argv) > 1 else DEFAULT_SCH
    conflicts, bridges = analyze(sexp.parse(open(path, encoding="utf-8").read())[0])
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
    sys.exit(main(sys.argv))
