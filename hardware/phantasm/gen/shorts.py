"""Union-find connectivity check over raw schematic geometry.
Reports edges that merge two DIFFERENT named nets (labels / power)."""
import os
import sys
import sexp
from kicad_common import F

path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "phantasm.kicad_sch")
root = sexp.parse(open(path, encoding="utf-8").read())[0]


def R(v):
    return (round(float(v[0]), 3), round(float(v[1]), 3))


parent = {}


def find(p):
    parent.setdefault(p, p)
    while parent[p] != p:
        parent[p] = parent[parent[p]]
        p = parent[p]
    return p


named = {}   # point -> net name (label/power)
causes = []  # (a,b,reason)


def union(a, b, reason):
    ra, rb = find(a), find(b)
    if ra != rb:
        parent[ra] = rb
    causes.append((a, b, reason))


# collect labels & power points
for c in root:
    if not isinstance(c, list) or not c:
        continue
    if c[0] == "label":
        at = sexp._val(c, "at"); p = R(at); named[p] = c[1]; find(p)
    if c[0] == "symbol":
        lib = sexp._val(c, "lib_id")[0]
        if lib.startswith("power:") or "phantasm:+" in lib:
            at = sexp._val(c, "at"); p = R(at)
            net = None
            for d in c:
                if isinstance(d, list) and d and d[0] == "property" and d[1] == "Value":
                    net = d[2]
            named[p] = net; find(p)

# wires: union endpoints; also any named point colinear on the segment
wires = []
for c in F(root, "wire"):
    pts = [R((p[1], p[2])) for p in sexp._val(c, "pts") if isinstance(p, list) and p[0] == "xy"]
    wires.append((pts[0], pts[1]))


def on_seg(p, a, b):
    (x, y), (x1, y1), (x2, y2) = p, a, b
    if not (min(x1, x2) - 0.02 <= x <= max(x1, x2) + 0.02 and
            min(y1, y2) - 0.02 <= y <= max(y1, y2) + 0.02):
        return False
    return abs((x2 - x1) * (y - y1) - (y2 - y1) * (x - x1)) < 0.02


# junction dots explicitly connect crossing/T wires that only touch mid-span.
junctions = []
for c in F(root, "junction"):
    at = sexp._val(c, "at")
    if at:
        junctions.append(R(at))

allpts = set(named)  # plus wire endpoints and junction dots
for a, b in wires:
    allpts.add(a); allpts.add(b)
allpts.update(junctions)

for a, b in wires:
    union(a, b, f"wire {a}-{b}")
    # Any point colinear on this span is on the wire: a named pin, another wire's
    # endpoint landing mid-span (a T-junction), or an explicit junction dot.
    for p in allpts:
        if p != a and p != b and on_seg(p, a, b):
            tag = f"label/power {named[p]}" if p in named else "junction"
            union(p, a, f"{tag} on wire {a}-{b}")

# same-coordinate merges already same key

# now report: groups with >1 distinct named net
groups = {}
for p, net in named.items():
    groups.setdefault(find(p), set()).add(net)
print("=== nets per electrical group (only conflicts) ===")
conf = 0
for g, nets in groups.items():
    if len(nets) > 1:
        conf += 1
        print("  MERGED:", sorted(nets))
print("conflict groups:", conf)
print("=== suspect bridging edges (touching 2 named pts) ===")
for a, b, reason in causes:
    if "on wire" in reason:
        print("  ", reason)

print("\n=== conflict group members (excl PWR_FLAG) ===")
real = 0
for g, nets in groups.items():
    if len(nets) > 1 and nets != {'GND','PWR_FLAG'} and 'PWR_FLAG' not in nets:
        real += 1
        print("GROUP", sorted(nets))
        for p, net in sorted(named.items()):
            if find(p) == g:
                print("   ", net, "@", p)

sys.exit(1 if real else 0)
