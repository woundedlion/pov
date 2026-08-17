"""Generate phantasm.kicad_pcb from the schematic + its netlist.

Embeds each component's footprint (from KiCad stock libs, or a generated Teensy
footprint), assigns pad nets by name from the exported netlist, places everything
linearly inside a <=35 mm-wide board outline (R-MECH-6), and declares the nets.
Emits a 4-layer SIG/GND/GND/SIG board: the physical stackup and the inner GND
planes are encoded in the file, so an autoplacer/fab reads them on upload.
Placement is a rough starting arrangement; route/refine interactively in Pcbnew.
"""
import argparse
import math
import os
import sys
import sexp
from kicad_common import (uid, reset_uid_sequence, fmt, F, arc_extrema,
                          export_netlist, find_kicad_cli, require_writable)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
SCH = os.path.join(OUT, "phantasm.kicad_sch")
PCB_FILE = "phantasm.kicad_pcb"
UNPLACED_FILE = os.path.join("unplaced", "phantasm_unplaced.kicad_pcb")
UNPLACED_REASON = (
    "The committed unplaced board is the Quilter upload input and carries\n"
    "  KiCad GUI edits these generators do not reproduce.")
FP_DIR = sexp.find_kicad_data_dir("footprints", "KICAD_FOOTPRINT_DIR")
KCLI = find_kicad_cli()
PCB_W = 32.0  # board width (mm); within the R-MECH-6 cap (<=35), trimmed to part extent
QUILTER_LENGTH = 58.28
TEENSY_LIBRARY_REASON = (
    "The committed, routed phantasm.kicad_pcb resolves its Teensy pads\n"
    "  against this library file.")
TEENSY_LIBID = "phantasm:Teensy4.0"
MOUNTING_HOLE_FOOTPRINT = "MountingHole:MountingHole_2.7mm_M2.5"
MOUNTING_KEEPOUT_RADIUS = 2.7
MOUNTING_HOLE_INSET = 3.5

# R-SI-1: SIG/GND/GND/SIG, so a trace on either outer layer has an adjacent
# reference plane. The inner planes are poured GND (GROUND_PLANE_LAYERS).
COPPER_LAYERS = ((0, "F.Cu", "signal"), (4, "In1.Cu", "power"),
                 (6, "In2.Cu", "power"), (2, "B.Cu", "signal"))
TECHNICAL_LAYERS = ((9, "F.Adhes", "user"), (11, "B.Adhes", "user"),
                    (13, "F.Paste", "user"), (15, "B.Paste", "user"),
                    (5, "F.SilkS", "user"), (7, "B.SilkS", "user"),
                    (1, "F.Mask", "user"), (3, "B.Mask", "user"),
                    (17, "Dwgs.User", "user"), (19, "Cmts.User", "user"),
                    (21, "Eco1.User", "user"), (23, "Eco2.User", "user"),
                    (25, "Edge.Cuts", "user"), (27, "Margin", "user"),
                    (31, "F.CrtYd", "user"), (29, "B.CrtYd", "user"),
                    (35, "F.Fab", "user"), (33, "B.Fab", "user"))
GROUND_NET = "GND"
GROUND_PLANE_LAYERS = ("In1.Cu", "In2.Cu")

# Standard 1.6 mm JLCPCB 4-layer build (JLC04161H-7628): 1 oz outer copper,
# 0.5 oz inner, thin outer prepreg over a thick core, ENIG. Encoded in the file
# so an autoplacer/fab reads the intended stackup without hand-entered values.
STACKUP = (
    '(layer "F.SilkS" (type "Top Silk Screen"))',
    '(layer "F.Paste" (type "Top Solder Paste"))',
    '(layer "F.Mask" (type "Top Solder Mask") (thickness 0.01))',
    '(layer "F.Cu" (type "copper") (thickness 0.035))',
    '(layer "dielectric 1" (type "prepreg") (thickness 0.2104) (material "FR4")'
    ' (epsilon_r 4.5) (loss_tangent 0.02))',
    '(layer "In1.Cu" (type "copper") (thickness 0.0152))',
    '(layer "dielectric 2" (type "core") (thickness 1.065) (material "FR4")'
    ' (epsilon_r 4.5) (loss_tangent 0.02))',
    '(layer "In2.Cu" (type "copper") (thickness 0.0152))',
    '(layer "dielectric 3" (type "prepreg") (thickness 0.2104) (material "FR4")'
    ' (epsilon_r 4.5) (loss_tangent 0.02))',
    '(layer "B.Cu" (type "copper") (thickness 0.035))',
    '(layer "B.Mask" (type "Bottom Solder Mask") (thickness 0.01))',
    '(layer "B.Paste" (type "Bottom Solder Paste"))',
    '(layer "B.SilkS" (type "Bottom Silk Screen"))',
    '(copper_finish "ENIG")',
    '(dielectric_constraints no)',
)


def copper_layer_names():
    return tuple(name for _, name, _ in COPPER_LAYERS)


def mounting_holes(L):
    """R-MECH-2: four M2.5 clearance holes, centres 3.5 mm in from each corner."""
    d = MOUNTING_HOLE_INSET
    return {"H1": (d, d), "H2": (d, PCB_W - d),
            "H3": (L - d, d), "H4": (L - d, PCB_W - d)}


def keepout_rects(L):
    """All-copper keepout square (screw head) around each mounting hole, as
    ref -> (min x, min y, max x, max y)."""
    kr = MOUNTING_KEEPOUT_RADIUS
    return {ref: (x - kr, y - kr, x + kr, y + kr)
            for ref, (x, y) in mounting_holes(L).items()}


def _boxes_overlap(a, b, tol=1e-6):
    return (a[0] < b[2] - tol and a[2] > b[0] + tol
            and a[1] < b[3] - tol and a[3] > b[1] + tol)


def keepout_clashes(place, pad_bxs, L):
    """Refs whose copper lands inside a mounting-hole keepout, as "REF/HOLE"."""
    keepouts = keepout_rects(L)
    out = []
    for ref, (x, y, rot) in sorted(place.items()):
        mnx, mny, mxx, mxy = _rot_bb(pad_bxs[ref], rot)
        box = (x + mnx, y + mny, x + mxx, y + mxy)
        out += [f"{ref}/{hole}" for hole in sorted(keepouts)
                if _boxes_overlap(box, keepouts[hole])]
    return out


def outline_overflows(place, pad_bxs, L, refs):
    """Refs of `refs` whose copper reaches past the L x PCB_W outline."""
    out = []
    for ref in sorted(refs):
        x, y, rot = place[ref]
        mnx, mny, mxx, mxy = _rot_bb(pad_bxs[ref], rot)
        if x + mnx < 0 or y + mny < 0 or x + mxx > L or y + mxy > PCB_W:
            out.append(f"{ref} [{x+mnx:.2f},{y+mny:.2f}]-"
                       f"[{x+mxx:.2f},{y+mxy:.2f}]")
    return out


def build_nets(nlroot):
    pad_net = {}            # (ref, pad) -> netname
    names = []
    for nb in F(nlroot, "nets"):
        for net in F(nb, "net"):
            name = sexp.val(net, "name")[0]
            names.append(name)
            for nd in F(net, "node"):
                pad_net[(sexp.val(nd, "ref")[0], sexp.val(nd, "pin")[0])] = name
    # stable net ids: 0 = "", then unique names in first-seen order
    netid = {"": 0}
    for n in names:
        if n not in netid:
            netid[n] = len(netid)
    return pad_net, netid


def build_paths(nlroot):
    """ref -> footprint (path) linking to the schematic symbol, from netlist tstamps."""
    out = {}
    for cb in F(nlroot, "components"):
        for comp in F(cb, "comp"):
            ref = sexp.val(comp, "ref")[0]
            ts = sexp.val(comp, "tstamps")
            if ts:
                out[ref] = "/" + ts[0].strip("/")
    return out


# ---------------------------------------------------------------- components
def schematic_components():
    """Return ordered unique component records, skipping power and flag symbols.

    The Teensy carries no Footprint field; every other symbol must, or its part
    would be embedded as the generated Teensy land.
    """
    with open(SCH, encoding="utf-8") as f:
        root = sexp.parse(f.read())[0]
    seen = {}
    order = []
    for c in root:
        if not (isinstance(c, list) and c and c[0] == "symbol"):
            continue
        ref = val = fp = None
        for d in c:
            if isinstance(d, list) and d and d[0] == "property":
                if d[1] == "Reference":
                    ref = d[2]
                elif d[1] == "Value":
                    val = d[2]
                elif d[1] == "Footprint":
                    fp = d[2]
        if not ref or ref.startswith("#"):
            continue
        if ref not in seen:
            dnp = sexp.val(c, "dnp", [sexp.Sym("no")])[0] == "yes"
            if not fp:
                libid = str(sexp.val(c, "lib_id", [""])[0])
                if libid != TEENSY_LIBID:
                    sys.exit(f"ERROR {ref} ({libid or 'no lib_id'}) has an empty "
                             "Footprint property; set one in Eeschema")
                fp = TEENSY_LIBID
            seen[ref] = (ref, fp, val or "", dnp)
            order.append(ref)
    return [seen[r] for r in order]


# ---------------------------------------------------------------- footprints
_MOD_CACHE = {}


def load_mod(libid):
    if libid not in _MOD_CACHE:
        lib, name = libid.split(":", 1)
        path = os.path.join(FP_DIR, lib + ".pretty", name + ".kicad_mod")
        with open(path, encoding="utf-8") as f:
            _MOD_CACHE[libid] = sexp.parse(f.read())[0]
    return _MOD_CACHE[libid]


def set_d_bus_land_pattern(node):
    geometry = {"1": (-1.1, 0), "2": (1.1, 0)}
    for pad in F(node, "pad"):
        number = str(pad[1])
        if number not in geometry:
            continue
        for field in pad:
            if not isinstance(field, list) or not field:
                continue
            if field[0] == "at":
                field[1:3] = geometry[number]
            elif field[0] == "size":
                field[1:3] = (0.8, 0.5)


def move_silk_graphics_to_fab(node):
    for item in node:
        if not isinstance(item, list) or not item or item[0] not in (
                "fp_line", "fp_arc", "fp_rect", "fp_poly"):
            continue
        for field in item:
            if (isinstance(field, list) and field and field[0] == "layer"
                    and field[1] == "F.SilkS"):
                field[1] = "F.Fab"


def set_pad_orientations(node, footprint_rotation):
    for pad in F(node, "pad"):
        for field in pad:
            if not isinstance(field, list) or not field or field[0] != "at":
                continue
            pad_rotation = float(field[3]) if len(field) > 3 else 0.0
            absolute_rotation = (pad_rotation + footprint_rotation) % 360
            if absolute_rotation:
                if len(field) > 3:
                    field[3] = sexp.Sym(fmt(absolute_rotation))
                else:
                    field.append(sexp.Sym(fmt(absolute_rotation)))
            elif len(field) > 3:
                del field[3:]


def refresh_uuids(node):
    for child in node:
        if not isinstance(child, list) or not child:
            continue
        if child[0] == "uuid":
            child[1] = uid()
        else:
            refresh_uuids(child)


def teensy_footprint(model_path="${KIPRJMOD}/phantasm.pretty/Teensy4.0.wrl"):
    """2x14 THT module, pads named by Teensy pin label; long axis along X.
    Pad map is the top view (component side up) with the USB end at -X:
    top row (-Y) VIN,GND,3V3,23..13 / bottom row (+Y) GND,0..12."""
    TOP = ["VIN", "GND", "3V3", "23", "22", "21", "20", "19", "18", "17", "16", "15", "14", "13"]
    BOTTOM = ["GND", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"]
    pads = []
    x0 = -((14 - 1) / 2) * 2.54
    for i, nm in enumerate(TOP):
        pads.append(f'\t(pad "{nm}" thru_hole oval (at {fmt(x0 + i*2.54)} -7.62) '
                    f'(size 1.7 1.7) (drill 1.0) (layers "*.Cu" "*.Mask"))')
    for i, nm in enumerate(BOTTOM):
        pads.append(f'\t(pad "{nm}" thru_hole oval (at {fmt(x0 + i*2.54)} 7.62) '
                    f'(size 1.7 1.7) (drill 1.0) (layers "*.Cu" "*.Mask"))')
    body = ('\t(fp_rect (start -18.5 -9.5) (end 18.5 9.5) (stroke (width 0.12) (type solid))'
            ' (fill none) (layer "F.SilkS"))\n'
            '\t(fp_rect (start -18.5 -9.5) (end 18.5 9.5) (stroke (width 0.05) (type solid))'
            ' (fill none) (layer "F.CrtYd"))\n'
            '\t(fp_rect (start -30 -5) (end -18.5 5) (stroke (width 0.25) (type dash))'
            ' (fill none) (layer "Dwgs.User"))')
    keepout = (
        '\t(zone\n'
        '\t\t(layer "F.Cu")\n'
        f'\t\t(uuid "{uid()}")\n'
        '\t\t(name "USB mating plug")\n'
        '\t\t(hatch full 0.5)\n'
        '\t\t(connect_pads (clearance 0))\n'
        '\t\t(min_thickness 0.25)\n'
        '\t\t(keepout (tracks allowed) (vias allowed) (pads allowed)'
        ' (copperpour allowed) (footprints not_allowed))\n'
        '\t\t(placement (enabled no) (sheetname ""))\n'
        '\t\t(fill (thermal_gap 0.3) (thermal_bridge_width 0.3))\n'
        '\t\t(polygon (pts (xy -30 -5) (xy -18.5 -5)'
        ' (xy -18.5 5) (xy -30 5)))\n'
        '\t)\n')
    model = (
        f'\t(model "{model_path}"\n'
        '\t\t(offset (xyz 0 0 0))\n'
        '\t\t(scale (xyz 1 1 1))\n'
        '\t\t(rotate (xyz 0 0 0))\n'
        '\t)\n')
    text = ('(footprint "phantasm:Teensy4.0"\n'
            '\t(layer "F.Cu")\n'
            '\t(descr "Teensy 4.0 module, 2x14 0.1in headers")\n'
            '\t(property "Reference" "U" (at 0 -11 0) (layer "F.SilkS")\n'
            '\t\t(effects (font (size 1 1) (thickness 0.15))))\n'
            '\t(property "Value" "Teensy4.0" (at 0 11 0) (layer "F.Fab")\n'
            '\t\t(effects (font (size 1 1) (thickness 0.15))))\n'
            '\t(attr through_hole)\n'
            + body + "\n" + keepout + "\n".join(pads) + "\n" + model + ")\n")
    return sexp.parse(text)[0]


def embedded_footprint(ref, libid,
                       teensy_model_path="${KIPRJMOD}/phantasm.pretty/Teensy4.0.wrl"):
    """A fresh footprint node carrying the land `embed` emits for `ref`,
    including any per-reference land override."""
    node = teensy_footprint(teensy_model_path) if libid == TEENSY_LIBID else \
        sexp.parse(sexp.dumps(load_mod(libid)))[0]  # deep copy via round-trip
    if ref == "D_BUS":
        set_d_bus_land_pattern(node)
    return node


def embed(libid, ref, value, x, y, rot, pad_net, netid, path=None, locked=False,
          dnp=False, hide_reference=False,
          teensy_model_path="${KIPRJMOD}/phantasm.pretty/Teensy4.0.wrl",
          consumed=None):
    node = embedded_footprint(ref, libid, teensy_model_path)
    refresh_uuids(node)
    if ref in ("D_BUS", "JP_ID2"):
        move_silk_graphics_to_fab(node)
    set_pad_orientations(node, rot)
    node[1] = libid
    # strip lib-file-only headers
    node[:] = [c for c in node if not (isinstance(c, list) and c and
               c[0] in ("version", "generator", "generator_version", "tedit"))]
    # ensure layer / at / uuid present and correct
    def setkv(key, value_list):
        for c in node:
            if isinstance(c, list) and c and c[0] == key:
                c[1:] = value_list
                return
        node.insert(2, [sexp.Sym(key)] + value_list)
    setkv("layer", ["F.Cu"])
    if locked:
        node.insert(2, [sexp.Sym("locked"), sexp.Sym("yes")])
    if dnp:
        # `(attr ...)` is optional in the footprint format, so a stock library
        # part may carry none; an empty one is the same "no flags" default.
        attr = next((c for c in node if isinstance(c, list) and c and c[0] == "attr"),
                    None)
        if attr is None:
            attr = [sexp.Sym("attr")]
            node.insert(2, attr)
        attr.append(sexp.Sym("dnp"))
    # insert at + uuid after layer
    node.insert(3, [sexp.Sym("at"), sexp.Sym(fmt(x)), sexp.Sym(fmt(y)), sexp.Sym(fmt(rot))])
    node.insert(4, [sexp.Sym("uuid"), uid()])
    # link to the schematic symbol so the board matches the schematic
    if path:
        node.insert(5, [sexp.Sym("path"), path])
    # set Reference / Value property text
    for c in node:
        if isinstance(c, list) and c and c[0] == "property":
            if c[1] == "Reference":
                c[2] = ref
                if hide_reference and not any(
                        isinstance(d, list) and d and d[0] == "hide" for d in c):
                    c.insert(-1, [sexp.Sym("hide"), sexp.Sym("yes")])
            elif c[1] == "Value":
                c[2] = value
    # assign pad nets by name; `consumed` collects the (ref, pad) keys that matched
    for c in node:
        if isinstance(c, list) and c and c[0] == "pad":
            padname = c[1]
            nn = pad_net.get((ref, padname))
            if nn is not None:
                if consumed is not None:
                    consumed.add((ref, padname))
                # remove any existing net, then add
                c[:] = [d for d in c if not (isinstance(d, list) and d and d[0] == "net")]
                c.append([sexp.Sym("net"), netid[nn], nn])
    return node


# ---------------------------------------------------------------- layout
# 2-D shelf strip-packer: minimise board LENGTH within the fixed WIDTH cap by
# stacking parts across the width in shelves (first-fit-decreasing). The Teensy
# (~37x19) sets the floor; small SMD parts pack into the leftover width beside it.
# Draft placement — refine orientation / push connectors to the edges in Pcbnew.
def _arc_points(start, mid, end):
    return [start, mid, end, *arc_extrema(start, mid, end)]


def fp_bbox(node, pads_only=False):
    """Footprint bounding box (minx,miny,maxx,maxy) in its local (origin) frame,
    over pads plus (unless `pads_only`) graphic outlines. Pad rotation is folded
    into a max-dim radius."""
    xs = []; ys = []
    for c in node:
        if not (isinstance(c, list) and c):
            continue
        if c[0] == "pad":
            at = sexp.val(c, "at"); sz = sexp.val(c, "size")
            x = float(at[0]); y = float(at[1])
            r = max(float(sz[0]), float(sz[1])) / 2 if sz else 0.5
            xs += [x - r, x + r]; ys += [y - r, y + r]
        elif not pads_only and c[0] == "fp_circle":
            center = sexp.val(c, "center")
            end = sexp.val(c, "end")
            if center and end:
                x = float(center[0]); y = float(center[1])
                radius = math.hypot(float(end[0]) - x, float(end[1]) - y)
                xs += [x - radius, x + radius]
                ys += [y - radius, y + radius]
        elif not pads_only and c[0] == "fp_arc":
            values = [sexp.val(c, key) for key in ("start", "mid", "end")]
            if all(values):
                points = [(float(value[0]), float(value[1])) for value in values]
                for x, y in _arc_points(*points):
                    xs.append(x); ys.append(y)
        elif not pads_only and c[0] in ("fp_rect", "fp_line", "fp_poly"):
            for k in ("start", "end", "center"):
                v = sexp.val(c, k)
                if v:
                    xs.append(float(v[0])); ys.append(float(v[1]))
            pts = sexp.val(c, "pts")
            if pts:
                for p in pts:
                    if isinstance(p, list) and p and p[0] == "xy":
                        xs.append(float(p[1])); ys.append(float(p[2]))
    if not xs:
        return (-1.0, -1.0, 1.0, 1.0)
    return (min(xs), min(ys), max(xs), max(ys))


def _rot_bb(bb, rot):
    """Bounding box of `bb` rotated by a right angle.

    The packer and the mounting-hole keepout gate measure parts through this,
    so an angle it cannot rotate is an error, not the unrotated box.
    """
    mnx, mny, mxx, mxy = bb
    pts = [(mnx, mny), (mxx, mny), (mxx, mxy), (mnx, mxy)]
    if rot == 90:
        pts = [(y, -x) for x, y in pts]
    elif rot == 180:
        pts = [(-x, -y) for x, y in pts]
    elif rot == 270:
        pts = [(-y, x) for x, y in pts]
    elif rot != 0:
        raise ValueError(f"placement rotation must be 0/90/180/270: {rot!r}")
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    return (min(xs), min(ys), max(xs), max(ys))


def _rotatable(ref):
    # Placement aligns the rotated bbox corner to the cell, so any origin packs
    # correctly. Rotate chip passives (R/C) and pin-header connectors (J*). Keep
    # solder jumpers (JP*, custom pads + clearance-outline trip DRC at 90), diodes,
    # fuse, bead and ICs at rot 0.
    if ref.startswith("JP"):
        return False
    return ref[0] in "RC" or ref[0] == "J"


# Through-hole connectors are pinned to the board ends (not skyline-packed) so the
# mating cables are accessible at the edge, per the spec's signal flow: power/debug
# at the hub end, strip + sync daisy at the far end (R-CON-4). Each 0.1in header is
# stood with its pin-row across the width (rot 0 for these 1xN vertical headers), so it
# hugs the end edge and adds almost no length. (v1 opposite-edge split — best SI; the
# one-edge variant crowds the narrow edge and lengthens the fast nets.)
HUB_CONNS = ("J1", "J4")            # logic power in, debug — hub end (left)
FAR_CONNS = ("J2", "J3A", "J3B")    # strip signal, sync daisy in/out — far end (right)

QUILTER_FIXED = {
    # J1's JST XA body runs along the length; it sits in the hub pocket between
    # the USB approach corridor and J4, clear of the Teensy courtyard.
    "J1": (8.5, 18.8, 0),
    "J4": (8.5, 28.5, 90),
    "J2": (49.8, 2.77, 0),
    "J3A": (49.8, 11.7, 0),
    "J3B": (49.8, 20.6, 0),
    "JP_ID0": (56.0, 8.5, 90),
    "JP_ID1": (56.0, 11.9, 90),
    "JP_ID2": (56.0, 15.3, 90),
    "JP_SHLD": (56.0, 19.0, 90),
    "U_MCU": (29.0, 11.7, 0),
    "C_IN": (20.5, 27.0, 0),
    "C_DEC1": (12.5, 1.35, 0),
    "U1": (36.6, 26.5, 180),
    "C_DEC2": (31.3, 30.0, 180),
    "R_D1": (42.2, 27.77, 180),
    "R_D2": (42.2, 23.96, 180),
    "R_S": (31.7, 24.0, 90),
    "R1": (29.7, 26.2, 90),
    "R2": (29.7, 22.8, 90),
    "C_SYNC": (27.8, 22.8, 90),
    "D_BUS": (46.0, 26.5, 90),
    "R_PD": (46.0, 23.0, 90),
}


def _column_height(refs, bxs, gap):
    """Y extent of a rot-0 column of `refs` separated by `gap`."""
    heights = [_rot_bb(bxs[r], 0)[3] - _rot_bb(bxs[r], 0)[1] for r in refs]
    return sum(heights) + gap * max(0, len(heights) - 1)


def _stack(refs, bxs, x0, edge, gap, place, y0=0.0):
    """Place `refs` as a column stacked along Y from y0 at left edge x0 (rot 0).
    Returns (column right-edge x, total stacked height)."""
    y = y0
    w = 0.0
    for ref in refs:
        rb = _rot_bb(bxs[ref], 0)
        place[ref] = (round(edge + x0 - rb[0], 3), round(edge + y - rb[1], 3), 0)
        y += (rb[3] - rb[1]) + gap
        w = max(w, rb[2] - rb[0])
    return x0 + w, y - gap - y0


def pack(bxs, width, edge=1.0, gap=1.2):
    """Connectors pinned to the ends + skyline (bottom-left-fill) interior pack.
    Board WIDTH is fixed; minimise LENGTH. Returns ({ref:(x,y,rot)}, length mm).

    The hub-end mounting-hole keepouts sit at a fixed x, so they seed the skyline;
    the far-end pair follows the length, so it is reserved as tail length past the
    rightmost part instead. Either way no part lands under a screw head."""
    usable = width - 2 * edge
    hub = [r for r in HUB_CONNS if r in bxs]
    far = [r for r in FAR_CONNS if r in bxs]
    pinned = set(hub) | set(far)
    place = {}
    sky = [(0.0, usable, 0.0)]    # skyline over the width: (y0, y1, x_right)

    def free_x(yb, h):            # leftmost x a [yb,yb+h] slot can start at
        x = 0.0
        for s0, s1, sx in sky:
            if s1 > yb + 1e-9 and s0 < yb + h - 1e-9:
                x = max(x, sx)
        return x

    def best_pos(h):             # (x, yb) bottom-left slot for a band of height h
        cands = {0.0}
        for s0, s1, _ in sky:
            cands.add(s0); cands.add(s1 - h)
        best = None
        for yb in sorted(cands):
            if yb < -1e-6 or yb + h > usable + 1e-6:
                continue
            r = (free_x(yb, h), yb)
            if best is None or r < best:
                best = r
        return best

    def reserve(yb, w, h, x):    # raise the skyline over [yb,yb+h] to x+w
        out = []
        for s0, s1, sx in sky:
            if s1 <= yb or s0 >= yb + h:
                out.append((s0, s1, sx)); continue
            if s0 < yb:
                out.append((s0, yb, sx))
            out.append((max(s0, yb), min(s1, yb + h), x + w))
            if s1 > yb + h:
                out.append((yb + h, s1, sx))
        out.sort()
        merged = []
        for seg in out:
            if merged and abs(merged[-1][2] - seg[2]) < 1e-6 \
                    and abs(merged[-1][1] - seg[0]) < 1e-6:
                merged[-1] = (merged[-1][0], seg[1], seg[2])
            else:
                merged.append(seg)
        return merged

    # 0) hub-end screw heads: block their squares before anything is placed
    kr = MOUNTING_KEEPOUT_RADIUS
    for cy in (MOUNTING_HOLE_INSET, width - MOUNTING_HOLE_INSET):
        yb = max(0.0, cy - kr - edge)
        yt = min(usable, cy + kr - edge)
        if yt > yb:
            sky = reserve(yb, 2 * kr, yt - yb, MOUNTING_HOLE_INSET - kr - edge)

    # 1) hub connectors: column at the left edge; reserve it so the interior packs right
    if hub:
        hubh = _column_height(hub, bxs, gap) + gap
        x0, y0 = best_pos(hubh) or (free_x(0.0, hubh), 0.0)
        hubw, _ = _stack(hub, bxs, x0, edge, gap, place, y0)
        sky = reserve(y0, hubw + gap, hubh, 0.0)

    # 2) interior parts: Teensy + SMD + jumpers, bottom-left skyline pack
    interior = sorted((r for r in bxs if r not in pinned),
                      key=lambda r: -max(bxs[r][2] - bxs[r][0], bxs[r][3] - bxs[r][1]))
    for ref in interior:
        bb = bxs[ref]
        rots = (0, 90) if _rotatable(ref) else (0,)
        choice = None
        for rot in rots:
            rb = _rot_bb(bb, rot)
            wg, hg = rb[2] - rb[0] + gap, rb[3] - rb[1] + gap
            if hg > usable + 1e-6:
                continue
            bp = best_pos(hg)
            if bp is None:
                continue
            x, yb = bp
            if choice is None or x < choice[0] - 1e-9:
                choice = (x, yb, rot, rb, wg, hg)
        if choice is None:        # too tall even rotated — force rot 0 at the end
            rb = _rot_bb(bb, 0)
            wg, hg = rb[2] - rb[0] + gap, rb[3] - rb[1] + gap
            x = max(s[2] for s in sky); yb = 0.0
            choice = (x, yb, 0, rb, wg, hg)
        x, yb, rot, rb, wg, hg = choice
        place[ref] = (round(edge + x - rb[0], 3), round(edge + yb - rb[1], 3), rot)
        sky = reserve(yb, wg, hg, x)

    # 3) far connectors: single column on the far edge. All I/O (power + strip +
    #    sync) is grouped here, so the column can be tall — fit the inter-connector
    #    gap to the usable width so the last pad keeps board-edge clearance.
    right = max(s[2] for s in sky)
    if far:
        sumh = _column_height(far, bxs, 0.0)
        fgap = gap
        if len(far) > 1:
            fgap = max(0.6, min(gap, (usable - sumh) / (len(far) - 1)))
        right, _ = _stack(far, bxs, right + gap, edge, fgap, place)

    # 4) length: enough tail for the far-end screw heads to clear every part
    extent = max((x + _rot_bb(bxs[ref], rot)[2] for ref, (x, _, rot) in place.items()),
                 default=edge + right)
    tail = max(edge, MOUNTING_HOLE_INSET + MOUNTING_KEEPOUT_RADIUS)
    return place, math.ceil((extent + tail) * 100 - 1e-9) / 100


def unplaced_layout(bxs, L, width, margin=2.0, gap=2.0):
    """Stage every footprint in a grid BELOW the board outline (rot 0), so an
    autoplacer (Quilter) sees them as unplaced. Returns {ref:(x,y,rot)}."""
    place = {}
    x = margin
    y = width + 6.0           # start clear of the outline (which is 0..width in Y)
    rowh = 0.0
    wrap = max(L, 40.0)
    for ref in sorted(bxs):
        mnx, mny, mxx, mxy = bxs[ref]
        w, h = mxx - mnx, mxy - mny
        if x + w > wrap - margin:
            x = margin; y += rowh + gap; rowh = 0.0
        place[ref] = (round(x - mnx, 3), round(y - mny, 3), 0)
        x += w + gap; rowh = max(rowh, h)
    return place


def main(unplaced=False, force=False, force_teensy_library=False):
    reset_uid_sequence()
    if unplaced:
        require_writable(os.path.join(OUT, UNPLACED_FILE), force, UNPLACED_REASON)
    else:
        require_writable(os.path.join(OUT, PCB_FILE), force)
    nlroot = export_netlist(KCLI, SCH)
    pad_net, netid = build_nets(nlroot)
    paths = build_paths(nlroot)                   # ref -> schematic-symbol path
    comps = {r: (r, fp, v, dnp) for r, fp, v, dnp in schematic_components()}
    # footprint bounding boxes -> 2-D shelf-pack to minimise length
    bxs = {}
    pad_bxs = {}
    for ref, (_, fp, _, _) in comps.items():
        node = embedded_footprint(ref, fp)
        bxs[ref] = fp_bbox(node)
        pad_bxs[ref] = fp_bbox(node, pads_only=True)
    if unplaced:
        L = QUILTER_LENGTH
        staged = unplaced_layout(bxs, L, PCB_W)
        PLACE = {r: QUILTER_FIXED.get(r, staged[r]) for r in bxs}
        # the staged grid sits below the outline by design; only locked parts are on it
        bounded = [r for r in bxs if r in QUILTER_FIXED]
        OUTFILE = UNPLACED_FILE
        NOTE = (f'PHANTASM segment board UNPLACED  -  {fmt(L)}x{fmt(PCB_W)}mm outline '
                '(width <=35mm); mechanical and signal-integrity placements locked')
    else:
        PLACE, L = pack(bxs, PCB_W)
        bounded = list(bxs)
        OUTFILE = PCB_FILE
        NOTE = (f'PHANTASM segment board  -  {fmt(L)}x{fmt(PCB_W)}mm (width <=35mm, R-MECH-6); '
                'shelf-packed draft, route in Pcbnew')

    outside = outline_overflows(PLACE, pad_bxs, L, bounded)
    if outside:
        sys.exit(f"ERROR placements outside the {fmt(L)}x{fmt(PCB_W)}mm outline: "
                 + ", ".join(outside))

    clashes = keepout_clashes(PLACE, pad_bxs, L)
    if clashes:
        sys.exit("ERROR pads inside a mounting-hole keepout: " + ", ".join(clashes))

    HOLES = mounting_holes(L)
    teensy_model_path = "${KIPRJMOD}/../phantasm.pretty/Teensy4.0.wrl" if unplaced else \
        "${KIPRJMOD}/phantasm.pretty/Teensy4.0.wrl"
    foot_nodes = []
    consumed = set()
    for ref, (x, y, rot) in PLACE.items():
        _, fp, val, dnp = comps[ref]
        lock = unplaced and ref in QUILTER_FIXED
        foot_nodes.append(embed(fp, ref, val, x, y, rot, pad_net, netid,
                                path=paths.get(ref), locked=lock, dnp=dnp,
                                hide_reference=lock,
                                teensy_model_path=teensy_model_path,
                                consumed=consumed))
    for ref, (x, y) in HOLES.items():
        foot_nodes.append(embed(MOUNTING_HOLE_FOOTPRINT, ref, "M2.5",
                                x, y, 0, pad_net, netid, locked=unplaced,
                                hide_reference=unplaced, consumed=consumed))
    # A netlist pin with no pad of that name would drop its connection silently;
    # refs that were never embedded drop their nets by design.
    embedded = set(PLACE) | set(HOLES)
    orphan = sorted(k for k in pad_net if k[0] in embedded and k not in consumed)
    if orphan:
        sys.exit("ERROR unmatched netlist pins (no footprint pad of that name): "
                 + ", ".join(f"{r}.{p}" for r, p in orphan))

    lines = []
    lines.append("(kicad_pcb")
    lines.append("\t(version 20250513)")
    lines.append('\t(generator "phantasm-gen")')
    lines.append('\t(generator_version "10.0")')
    lines.append("\t(general (thickness 1.6) (legacy_teardrops no))")
    lines.append('\t(paper "A2")')
    lines.append("\t(layers")
    for n, nm, ty in COPPER_LAYERS + TECHNICAL_LAYERS:
        lines.append(f'\t\t({n} "{nm}" {ty})')
    lines.append("\t)")
    lines.append("\t(setup")
    lines.append("\t\t(stackup")
    lines += [f"\t\t\t{layer}" for layer in STACKUP]
    lines.append("\t\t)")
    lines.append("\t\t(pad_to_mask_clearance 0)")
    lines.append("\t)")
    # nets
    for nm, i in sorted(netid.items(), key=lambda kv: kv[1]):
        lines.append(f'\t(net {i} {sexp.quote(nm)})')
    # board outline (Edge.Cuts) — <=35 mm wide strip, length minimised by packer
    lines.append(f'\t(gr_rect (start 0 0) (end {fmt(L)} {fmt(PCB_W)}) '
                 '(stroke (width 0.15) (type solid)) (fill none) (layer "Edge.Cuts") '
                 f'(uuid "{uid()}"))')
    # inner reference planes over the whole outline (R-SI-1), left unfilled for
    # the router to pour
    if GROUND_NET not in netid:
        sys.exit(f"ERROR no {GROUND_NET} net to pour the inner planes on")
    for layer in GROUND_PLANE_LAYERS:
        lines.append('\t(zone')
        lines.append(f'\t\t(net {netid[GROUND_NET]})')
        lines.append(f'\t\t(net_name "{GROUND_NET}")')
        lines.append(f'\t\t(layer "{layer}")')
        lines.append(f'\t\t(uuid "{uid()}")')
        lines.append(f'\t\t(name "{GROUND_NET} plane {layer}")')
        lines.append('\t\t(hatch edge 0.5)')
        lines.append('\t\t(connect_pads (clearance 0.5))')
        lines.append('\t\t(min_thickness 0.25)')
        lines.append('\t\t(fill yes (thermal_gap 0.5) (thermal_bridge_width 0.5))')
        lines.append('\t\t(polygon (pts '
                     f'(xy 0 0) (xy {fmt(L)} 0) '
                     f'(xy {fmt(L)} {fmt(PCB_W)}) (xy 0 {fmt(PCB_W)})))')
        lines.append('\t)')
    for ref, (kx0, ky0, kx1, ky1) in keepout_rects(L).items():
        lines.append('\t(zone')
        lines.append('\t\t(net 0)')
        lines.append('\t\t(net_name "")')
        lines.append('\t\t(layers '
                     + " ".join(f'"{nm}"' for nm in copper_layer_names()) + ')')
        lines.append(f'\t\t(uuid "{uid()}")')
        lines.append(f'\t\t(name "{ref} screw head")')
        lines.append('\t\t(hatch full 0.5)')
        lines.append('\t\t(connect_pads (clearance 0))')
        lines.append('\t\t(min_thickness 0.25)')
        lines.append('\t\t(keepout (tracks not_allowed) (vias not_allowed) (pads allowed)'
                     ' (copperpour not_allowed) (footprints allowed))')
        lines.append('\t\t(placement (enabled no) (sheetname ""))')
        lines.append('\t\t(fill (thermal_gap 0.3) (thermal_bridge_width 0.3))')
        lines.append('\t\t(polygon (pts '
                     f'(xy {fmt(kx0)} {fmt(ky0)}) (xy {fmt(kx1)} {fmt(ky0)}) '
                     f'(xy {fmt(kx1)} {fmt(ky1)}) (xy {fmt(kx0)} {fmt(ky1)})))')
        lines.append('\t)')
    lines.append(f'\t(gr_text {sexp.quote(NOTE)}'
                 f' (at 4 -4 0) (layer "Cmts.User") (uuid "{uid()}") '
                 '(effects (font (size 2 2) (thickness 0.3)) (justify left bottom)))')
    # Coordinates below annotate the QUILTER_FIXED placement; a fresh pack moves
    # those parts, so they are emitted only with the locked layout.
    if unplaced:
        lines.append(f'\t(gr_text "USB PLUG KEEP-OUT" (at 5.3 11.5 90)'
                     f' (layer "Dwgs.User") (uuid "{uid()}") '
                     '(effects (font (size 0.8 0.8) (thickness 0.15))))')
        front_silk = [
            ("S", 48.0, 11.7, 0),
            ("G", 48.0, 14.24, 0),
            ("H", 48.0, 16.78, 0),
            ("S", 48.0, 20.6, 0),
            ("G", 48.0, 23.14, 0),
            ("H", 48.0, 25.68, 0),
            ("SYNC IN", 52.0, 14.24, 90),
            ("SYNC OUT", 52.0, 23.14, 90),
            ("ID0", 53.9, 8.9, 90),
            ("ID1", 54.3, 11.9, 90),
            ("ID2", 54.1, 15.3, 90),
            ("SHLD", 54.3, 19.0, 90),
        ]
        for text, x, y, angle in front_silk:
            lines.append(f'\t(gr_text "{text}" (at {fmt(x)} {fmt(y)} {angle})'
                         f' (layer "F.SilkS") (uuid "{uid()}") '
                         '(effects (font (size 0.8 0.8) (thickness 0.15))))')
    back_silk = [
        ("ID  ID0   ID1   ROLE", 7.0, 0.9),
        ("0   OPEN  OPEN  MASTER", 9.0, 0.9),
        ("1   GND   OPEN  A-SOUTH", 11.0, 0.9),
        ("2   OPEN  GND   B-NORTH", 13.0, 0.9),
        ("3   GND   GND   B-SOUTH", 15.0, 0.9),
        ("MASTER = ALL ID OPEN   SHLD = MASTER ONLY", 17.0, 0.8),
        ("BOARD ID: ____", 23.5, 2.0),
        ("Phantasm Rev 1.1", 29.5, 1.0),
    ]
    for text, y, size in back_silk:
        lines.append(f'\t(gr_text "{text}" (at {fmt(L/2)} {fmt(y)} 0)'
                     f' (layer "B.SilkS") (uuid "{uid()}") '
                     f'(effects (font (size {fmt(size)} {fmt(size)})'
                     f' (thickness {fmt(max(0.15, size*0.15))})) (justify mirror)))')
    for fn in foot_nodes:
        lines.append(sexp.dumps(fn, indent=1))
    lines.append(")")
    outpath = os.path.join(OUT, OUTFILE)
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    with open(outpath, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")

    # --- custom footprint library (Teensy) + fp-lib-table ---
    pretty = os.path.join(OUT, "phantasm.pretty")
    os.makedirs(pretty, exist_ok=True)
    mod = teensy_footprint()
    mod[1] = "Teensy4.0"
    mod.insert(2, [sexp.Sym("version"), sexp.Sym("20240108")])
    mod.insert(3, [sexp.Sym("generator"), "phantasm-gen"])
    mod.insert(4, [sexp.Sym("generator_version"), "10.0"])
    mod_path = os.path.join(pretty, "Teensy4.0.kicad_mod")
    mod_text = sexp.dumps(mod) + "\n"
    if os.path.exists(mod_path):
        with open(mod_path, encoding="utf-8") as f:
            existing_mod_text = f.read()
    else:
        existing_mod_text = None
    if existing_mod_text != mod_text:
        require_writable(mod_path, force_teensy_library, TEENSY_LIBRARY_REASON,
                         flag="--force-teensy-library")
        with open(mod_path, "w", encoding="utf-8", newline="\n") as f:
            f.write(mod_text)
    fplt = os.path.join(OUT, "fp-lib-table")
    if not os.path.exists(fplt):
        with open(fplt, "w", encoding="utf-8", newline="\n") as f:
            f.write(
                '(fp_lib_table\n\t(version 7)\n'
                '\t(lib (name "phantasm")(type "KiCad")(uri "${KIPRJMOD}/phantasm.pretty")'
                '(options "")(descr "PHANTASM custom footprints"))\n)\n')
    if unplaced:
        with open(os.path.join(OUT, "unplaced", "fp-lib-table"), "w", encoding="utf-8",
                  newline="\n") as f:
            f.write(
                '(fp_lib_table\n\t(version 7)\n'
                '\t(lib (name "phantasm")(type "KiCad")(uri "${KIPRJMOD}/../phantasm.pretty")'
                '(options "")(descr "PHANTASM custom footprints"))\n)\n')
    print(f"wrote {OUTFILE}  footprints:{len(foot_nodes)} nets:{len(netid)} length:{L:.0f}mm")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--unplaced", action="store_true",
                        help="write unplaced/phantasm_unplaced.kicad_pcb for "
                             "the autoplacer instead")
    parser.add_argument("--force", action="store_true",
                        help=f"overwrite the committed {PCB_FILE} / "
                             f"{UNPLACED_FILE}")
    parser.add_argument("--force-teensy-library", action="store_true",
                        help="overwrite phantasm.pretty/Teensy4.0.kicad_mod "
                             "when the generated footprint differs; the routed "
                             "board resolves its Teensy pads against it")
    return parser.parse_args(argv)


if __name__ == "__main__":
    ARGS = parse_args()
    main(unplaced=ARGS.unplaced, force=ARGS.force,
         force_teensy_library=ARGS.force_teensy_library)
