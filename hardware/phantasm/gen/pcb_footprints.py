"""Footprint loading, land overrides and embedding."""
import copy
import math
import os
import sys
from pcb_emit import zone_sexpr
from constraints import KEEPOUT_FLOORS
import connectivity
import sexp
from constraints import (NEW_LAYOUT_RULES)
from kicad_common import (uid, fmt, F, arc_points,
                          rotate_point)
from pcb_model import COURTYARD_LAYERS
from pcb_model import HERE
from pcb_model import OUT
from pcb_model import PCB_W
from pcb_model import TEENSY_LIBID
from pcb_model import TEENSY_MODEL_PATH
from pcb_model import TERMINAL_LIBIDS
from pcb_layouts import generation_context
from pcb_pack import _boxes_overlap
from pcb_pack import _rot_bb
from pcb_pack import mounting_reserve_rects
from pcb_pack import mounting_holes
from pcb_model import MOUNTING_HOLE_FOOTPRINT
from constraints import is_hand_assembled


LOCAL_FOOTPRINT_DIR = os.path.join(OUT, "phantasm.pretty")


FP_DIR = sexp.find_kicad_data_dir("footprints", "KICAD_FOOTPRINT_DIR")


_MOD_CACHE = {}


def terminal_library_dir(libid):
    if libid == TERMINAL_LIBIDS[2]:
        return os.path.join(HERE, "footprints")
    return LOCAL_FOOTPRINT_DIR


def load_mod(libid):
    lib, name = libid.split(":", 1)
    if lib != "phantasm" and not os.path.isdir(FP_DIR):
        raise RuntimeError(
            "KiCad stock footprint libraries not found; set KICAD_FOOTPRINT_DIR "
            "to a share/kicad/footprints directory, or install KiCad in one "
            "of: " + ", ".join(sexp.kicad_data_dir_patterns("footprints")))
    directory = terminal_library_dir(libid) if lib == "phantasm" else \
        os.path.join(FP_DIR, lib + ".pretty")
    key = (directory, libid)
    if key not in _MOD_CACHE:
        path = os.path.join(directory, name + ".kicad_mod")
        with open(path, encoding="utf-8") as f:
            _MOD_CACHE[key] = sexp.parse(f.read())[0]
    return _MOD_CACHE[key]


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
            child[1] = uid(site="pcb.py:refresh_uuids")
        else:
            refresh_uuids(child)


def component_keepout(name, x0, y0, x1, y1, uuid):
    """Serialize a front-side component reservation."""
    return zone_sexpr(
        ("F.Cu",), name, uuid, ((x0, y0), (x1, y0), (x1, y1), (x0, y1)),
        KEEPOUT_FLOORS, keepout="(tracks allowed) (vias allowed) (pads allowed)"
        " (copperpour allowed) (footprints not_allowed)")


def teensy_footprint(ctx=None):
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
    outline_layer = 'F.Fab' if (ctx or generation_context()).is_rev13 else 'F.SilkS'
    body = ('\t(fp_rect (start -18.5 -9.5) (end 18.5 9.5) (stroke (width 0.12) (type solid))'
            f' (fill none) (layer "{outline_layer}"))\n'
            '\t(fp_rect (start -18.5 -9.5) (end 18.5 9.5) (stroke (width 0.05) (type solid))'
            ' (fill none) (layer "F.CrtYd"))\n'
            '\t(fp_rect (start -30 -5) (end -18.5 5) (stroke (width 0.25) (type dash))'
            ' (fill none) (layer "Dwgs.User"))'
            '\n\t(fp_text user "USB" (at -17.3 0 90) (layer "F.SilkS")'
            ' (effects (font (size 1 1) (thickness 0.15))))')
    keepout = component_keepout("USB mating plug", -30, -5, -18.5, 5, uid(site="pcb.py:teensy_footprint"))
    model = (
        f'\t(model {sexp.quote(TEENSY_MODEL_PATH)}\n'
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


def embedded_footprint(ref, libid, ctx=None):
    """A fresh footprint node carrying the land `embed` emits for `ref`,
    including any per-reference land override."""
    node = teensy_footprint(ctx) if libid == TEENSY_LIBID else \
        copy.deepcopy(load_mod(libid))
    if ref == "D_BUS":
        set_d_bus_land_pattern(node)
    if libid not in TERMINAL_LIBIDS:
        position_reference(node, ref)
    return node


WIRE_APPROACH_MM = 10.0


WIRE_ENVELOPE_MM = 2.4


def reserve_terminal_wire_access(node):
    """Reserve each terminal's nominal insulated-wire insertion envelope."""
    node[:] = [n for n in node if not (isinstance(n, list) and n and n[0] == "zone"
               and str(sexp.val(n, "name", [""])[0]).startswith("Terminal wire access"))]
    node[:] = [n for n in node if not (
        isinstance(n, list) and n and n[0] == "fp_rect"
        and sexp.val(n, "layer") == ["Dwgs.User"]
        and max(float(sexp.val(n, "start")[0]),
                float(sexp.val(n, "end")[0])) <= -3.25)]
    for pad in F(node, "pad"):
        y = float(sexp.val(pad, "at")[1])
        x0, x1 = -3.25 - WIRE_APPROACH_MM, -3.25
        y0, y1 = y - WIRE_ENVELOPE_MM / 2, y + WIRE_ENVELOPE_MM / 2
        node.append(sexp.parse_one(component_keepout(
            f"Terminal wire access pin {pad[1]}", x0, y0, x1, y1, uid(site="pcb.py:reserve_terminal_wire_access"))))


def validate_front_access(footprints, length):
    """Reject front components inside wire-entry and USB mating reservations."""
    boxes = {}
    for node in footprints:
        ref = connectivity.footprint_reference(node)
        if sexp.val(node, "layer") != ["F.Cu"]:
            raise ValueError(f"front-only layout requires F.Cu: {ref}")
        if ref.startswith("H"):
            continue
        x, y, angle = map(float, sexp.val(node, "at"))
        x0, y0, x1, y1 = _rot_bb(fp_bbox(node, graphic_layers=COURTYARD_LAYERS), angle)
        boxes[ref] = (x + x0, y + y0, x + x1, y + y1)
    for node in footprints:
        ref = connectivity.footprint_reference(node)
        for zone in F(node, "zone"):
            name = str(sexp.val(zone, "name", [""])[0])
            if not (name.startswith("Terminal wire access") or name == "USB mating plug"):
                continue
            points = F(F(F(zone, "polygon")[0], "pts")[0], "xy")
            xs, ys = zip(*[(float(p[1]), float(p[2])) for p in points])
            box = min(xs), min(ys), max(xs), max(ys)
            for other, bounds in {**boxes, **mounting_reserve_rects(length)}.items():
                if other != ref and _boxes_overlap(box, bounds):
                    raise ValueError(f"{ref} {name} obstructed by {other}")


def position_reference(node, ref):
    bounds = fp_bbox(node, graphic_layers=("F.SilkS", "F.CrtYd"))
    for prop in F(node, "property"):
        if prop[1] != "Reference":
            continue
        prop[2] = ref
        for child in prop:
            if isinstance(child, list) and child and child[0] == "at":
                child[1:] = [sexp.Sym(fmt((bounds[0] + bounds[2]) / 2)),
                             sexp.Sym(fmt(bounds[1] - 0.8)), 0]
            elif isinstance(child, list) and child and child[0] == "effects":
                child[1:] = sexp.parse_one(
                    "(effects (font (size 0.8 0.8) (thickness 0.15)))")[1:]


def reference_bbox(node, ref):
    bounds = fp_bbox(node)
    for prop in F(node, "property"):
        if prop[1] == "Reference":
            x, y = map(float, sexp.val(prop, "at")[:2])
            half_width = len(ref) * 0.8 * 0.6 + 0.075
            half_height = 0.475
            at = sexp.val(prop, "at")
            if len(at) > 2 and float(at[2]) % 180 == 90:
                half_width, half_height = half_height, half_width
            bounds = (min(bounds[0], x - half_width), min(bounds[1], y - half_height),
                      max(bounds[2], x + half_width), max(bounds[3], y + half_height))
    return bounds


def embed(libid, ref, value, x, y, rot, pad_net, netid, path=None, locked=False,
          dnp=False, hand_assembled=False, hide_reference=False,
          consumed=None, ctx=None):
    ctx = ctx or generation_context()
    node = embedded_footprint(ref, libid, ctx)
    if ctx.is_rev13 and libid in TERMINAL_LIBIDS:
        reserve_terminal_wire_access(node)
    refresh_uuids(node)
    set_pad_orientations(node, rot)
    # Embedded zone vertices use board coordinates; library vertices are local.
    for zone in F(node, "zone"):
        for polygon in F(zone, "polygon"):
            for points in F(polygon, "pts"):
                for point in F(points, "xy"):
                    px, py = rotate_point((float(point[1]), float(point[2])), rot)
                    point[1:3] = [sexp.Sym(fmt(x + px)), sexp.Sym(fmt(y + py))]
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
    flags = ["dnp"] if dnp else []
    if hand_assembled:
        flags += ["exclude_from_pos_files", "exclude_from_bom"]
    if flags:
        # `(attr ...)` is optional in the footprint format, so a stock library
        # part may carry none; an empty one is the same "no flags" default.
        attr = next((c for c in node if isinstance(c, list) and c and c[0] == "attr"),
                    None)
        if attr is None:
            attr = [sexp.Sym("attr")]
            node.insert(2, attr)
        carried = {str(flag) for flag in attr[1:]}
        attr += [sexp.Sym(flag) for flag in flags if flag not in carried]
    node.insert(3, [sexp.Sym("at"), sexp.Sym(fmt(x)), sexp.Sym(fmt(y)), sexp.Sym(fmt(rot))])
    node.insert(4, [sexp.Sym("uuid"), uid(site="pcb.py:embed")])
    # link to the schematic symbol so the board matches the schematic
    if path:
        node.insert(5, [sexp.Sym("path"), path])
    # set Reference / Value property text
    for c in node:
        if isinstance(c, list) and c and c[0] == "property":
            if c[1] == "Reference":
                c[2] = ref
                for at in F(c, "at"):
                    local_rotation = float(at[3]) if len(at) > 3 else 0
                    at[3:] = [sexp.Sym(fmt((local_rotation + rot) % 360))]
                if hide_reference and not any(
                        isinstance(d, list) and d and d[0] == "hide" for d in c):
                    c.insert(-1, [sexp.Sym("hide"), sexp.Sym("yes")])
            elif c[1] == "Value":
                c[2] = value
    properties = ctx.layout["part_properties"]
    for name, value in properties.get(ref, {}).items():
        existing = next((p for p in F(node, "property") if p[1] == name), None)
        if existing is not None:
            existing[2] = value
        else:
            node.append(sexp.parse_one(
                f'(property {sexp.quote(name)} {sexp.quote(value)} (at 0 0 0) '
                '(layer "F.Fab") (hide yes) (effects (font (size 1 1) (thickness 0.15))))'))
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


def trim_terminal_silkscreen(node, edge_distance):
    """Open the local -X outline at the adjacent board edge."""
    x = edge_distance
    for rectangle in list(F(node, "fp_rect")):
        if sexp.val(rectangle, "layer") != ["F.SilkS"]:
            continue
        x0, y0 = map(float, sexp.val(rectangle, "start"))
        x1, y1 = map(float, sexp.val(rectangle, "end"))
        stroke = F(rectangle, "stroke")[0]
        width = float(sexp.val(stroke, "width")[0])
        limit = NEW_LAYOUT_RULES["min_silk_clearance"] + width / 2 - x
        if x0 >= limit:
            continue
        node.remove(rectangle)
        if x1 <= limit:
            continue
        for start, end in (((limit, y0), (x1, y0)),
                           ((x1, y0), (x1, y1)),
                           ((x1, y1), (limit, y1))):
            node.append([sexp.Sym("fp_line"),
                         [sexp.Sym("start"), *(sexp.Sym(fmt(v)) for v in start)],
                         [sexp.Sym("end"), *(sexp.Sym(fmt(v)) for v in end)],
                         copy.deepcopy(stroke), [sexp.Sym("layer"), "F.SilkS"]])


def move_clipped_silkscreen_to_fab(node, length):
    """Keep component geometry that crosses the silk margin on the fab layer."""
    x, y, angle = map(float, sexp.val(node, "at"))
    for kind in ("fp_line", "fp_arc", "fp_circle", "fp_rect", "fp_poly"):
        for graphic in F(node, kind):
            if sexp.val(graphic, "layer") != ["F.SilkS"]:
                continue
            bounds = fp_bbox([sexp.Sym("footprint"), graphic], graphic_layers=("F.SilkS",))
            x0, y0, x1, y1 = _rot_bb(bounds, angle)
            stroke = float(sexp.val(F(graphic, "stroke")[0], "width")[0])
            margin = NEW_LAYOUT_RULES["min_silk_clearance"] + stroke / 2
            if (x + x0 < margin - 1e-6 or y + y0 < margin - 1e-6
                    or x + x1 > length - margin + 1e-6 or y + y1 > PCB_W - margin + 1e-6):
                F(graphic, "layer")[0][1] = "F.Fab"


def _on_layers(node, layers):
    if layers is None:
        return True
    names = sexp.val(node, "layer") or sexp.val(node, "layers") or []
    return any(str(name) in layers for name in names)


def fp_bbox(node, pads_only=False, graphic_layers=None):
    """Footprint bounding box (minx,miny,maxx,maxy) in its local (origin) frame,
    over pads plus (unless `pads_only`) graphic outlines, restricted to
    `graphic_layers` when given. Custom pads use primitive copper bounds;
    other pads use a max-dimension radius."""
    xs = []; ys = []
    for c in node:
        if not (isinstance(c, list) and c):
            continue
        if c[0] == "pad":
            if str(c[3]) == "custom":
                copper = connectivity.pad_copper(c, (0, 0), 0, ["F.Cu"])
                x0, y0, x1, y1 = connectivity.copper_bounds(copper)
                xs += [x0, x1]; ys += [y0, y1]
                continue
            at = sexp.val(c, "at"); sz = sexp.val(c, "size")
            x = float(at[0]); y = float(at[1])
            r = max(float(sz[0]), float(sz[1])) / 2 if sz else 0.5
            xs += [x - r, x + r]; ys += [y - r, y + r]
            continue
        if pads_only or not _on_layers(c, graphic_layers):
            continue
        if c[0] == "fp_circle":
            center = sexp.val(c, "center")
            end = sexp.val(c, "end")
            if center and end:
                x = float(center[0]); y = float(center[1])
                radius = math.hypot(float(end[0]) - x, float(end[1]) - y)
                xs += [x - radius, x + radius]
                ys += [y - radius, y + radius]
        elif c[0] == "fp_arc":
            values = [sexp.val(c, key) for key in ("start", "mid", "end")]
            if all(values):
                points = [(float(value[0]), float(value[1])) for value in values]
                for x, y in arc_points(*points):
                    xs.append(x); ys.append(y)
        elif c[0] in ("fp_rect", "fp_line", "fp_poly"):
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


def embed_components(ctx, netlist, placement):
    """Embed placed parts and mounting holes, rejecting unmatched netlist pins."""
    placements, length, fixed = placement.positions, placement.length, placement.fixed
    comps, pad_net, netid, paths = netlist.components, netlist.pad_net, netlist.netid, netlist.paths
    unplaced, layout = ctx.unplaced, ctx.layout
    holes = mounting_holes(length)
    foot_nodes = []
    consumed = set()
    for ref, (x, y, rot) in placements.items():
        _, fp, val, dnp = comps[ref]
        lock = unplaced and ref in fixed
        foot_nodes.append(embed(fp, ref, val, x, y, rot, pad_net, netid,
                                path=paths.get(ref), locked=lock, dnp=dnp,
                                hand_assembled=is_hand_assembled(fp, val),
                                hide_reference=lock,
                                consumed=consumed, ctx=ctx))
        if unplaced and ref in layout["edge_placements"]:
            if ctx.is_rev13 or ref == "J1":
                trim_terminal_silkscreen(foot_nodes[-1], y if rot == 270 else x)
        if ctx.is_rev13:
            move_clipped_silkscreen_to_fab(foot_nodes[-1], length)
    for ref, (x, y) in holes.items():
        foot_nodes.append(embed(MOUNTING_HOLE_FOOTPRINT, ref, "M2.5",
                                x, y, 0, pad_net, netid, locked=unplaced,
                                hide_reference=True, consumed=consumed, ctx=ctx))
    if unplaced and ctx.is_rev13:
        validate_front_access(foot_nodes, length)
    # A netlist pin with no pad of that name would drop its connection silently;
    # refs that were never embedded drop their nets by design.
    embedded = set(placements) | set(holes)
    orphan = sorted(k for k in pad_net if k[0] in embedded and k not in consumed)
    if orphan:
        sys.exit("ERROR unmatched netlist pins (no footprint pad of that name): "
                 + ", ".join(f"{r}.{p}" for r, p in orphan))
    return foot_nodes
