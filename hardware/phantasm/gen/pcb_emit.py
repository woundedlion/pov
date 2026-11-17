"""KiCad board serialization in deterministic generation order."""
import sys
import connectivity
import sexp
from constraints import (GROUND_PLANE_LAYERS, ZONE_DEFAULTS, UNPLACED_RULES,
                         MIN_SOLDER_MASK_WEB_MM, KEEPOUT_FLOORS, ZoneFloors)
from kicad_common import (uid, fmt, F)
from pcb_model import BOARD_THICKNESS_MM as BOARD_THICKNESS_MM
from pcb_model import COPPER_LAYERS as COPPER_LAYERS
from pcb_model import GROUND_NET as GROUND_NET
from pcb_model import LAYER_NAMES as LAYER_NAMES
from pcb_model import PCB_W as PCB_W
from pcb_model import STACKUP as STACKUP
from pcb_model import TECHNICAL_LAYERS as TECHNICAL_LAYERS
from pcb_model import copper_layer_names as copper_layer_names
from pcb_pack import keepout_rects as keepout_rects
from pcb_layouts import local_routes as local_routes
from pcb_layouts import pin_silk_marks as pin_silk_marks


def header(ctx):
    """Serialize the board's header."""
    selected = ctx.revision
    lines = []
    lines.append("(kicad_pcb")
    lines.append(f"\t(version {sexp.PCB_FORMAT})")
    lines.append('\t(generator "phantasm-gen")')
    lines.append(f'\t(generator_version "{sexp.GENERATOR_VERSION}")')
    # KiCad writes this rev into the Gerber X2 ProjectId attribute, which reads
    # `rev?` when the board carries no title block.
    lines.append("\t(title_block")
    lines.append(f'\t\t(rev {sexp.quote(selected)})')
    lines.append("\t)")
    lines.append(f"\t(general (thickness {BOARD_THICKNESS_MM}) (legacy_teardrops no))")
    lines.append('\t(paper "A2")')
    return lines


def layers_setup():
    """Serialize the board's layers setup."""
    lines = []
    lines.append("\t(layers")
    for n, nm, ty in COPPER_LAYERS + TECHNICAL_LAYERS:
        alias = f' {sexp.quote(LAYER_NAMES[nm])}' if nm in LAYER_NAMES else ""
        lines.append(f'\t\t({n} "{nm}" {ty}{alias})')
    lines.append("\t)")
    lines.append("\t(setup")
    lines.append("\t\t(stackup")
    lines += [f"\t\t\t{layer}" for layer in STACKUP]
    lines.append("\t\t)")
    lines.append("\t\t(pad_to_mask_clearance 0)")
    lines.append(f"\t\t(solder_mask_min_width {fmt(MIN_SOLDER_MASK_WEB_MM)})")
    lines.append("\t\t(allow_soldermask_bridges_in_footprints no)")
    lines.append("\t\t(tenting (front yes) (back yes))")
    lines.append("\t)")
    return lines


def nets_outline(netlist, placement):
    """Serialize the board's nets outline."""
    netid, length = netlist.netid, placement.length
    lines = []
    # nets
    for nm, i in sorted(netid.items(), key=lambda kv: kv[1]):
        lines.append(f'\t(net {i} {sexp.quote(nm)})')
    # board outline (Edge.Cuts)
    lines.append(f'\t(gr_rect (start 0 0) (end {fmt(length)} {fmt(PCB_W)}) '
                 '(stroke (width 0.15) (type solid)) (fill none) (layer "Edge.Cuts") '
                 f'(uuid "{uid(site="pcb.py:main")}"))')
    return lines


def ground_planes(netlist, placement):
    """Serialize unfilled inner reference planes over the outline."""
    if GROUND_NET not in netlist.netid:
        sys.exit(f"ERROR no {GROUND_NET} net to pour the inner planes on")
    points = ((0, 0), (placement.length, 0),
              (placement.length, PCB_W), (0, PCB_W))
    return [line for layer in GROUND_PLANE_LAYERS
            for line in zone_sexpr(
                (layer,), f"{GROUND_NET} plane {layer}", uid(site="pcb.py:main"), points,
                ZONE_DEFAULTS, net=(netlist.netid[GROUND_NET], GROUND_NET),
                clearance=UNPLACED_RULES["min_clearance"]).splitlines()]


def keepout_zones(placement):
    """Serialize the mounting-hole screw-head copper reservations."""
    lines = []
    for ref, (x0, y0, x1, y1) in keepout_rects(placement.length).items():
        lines.extend(zone_sexpr(
            copper_layer_names(), f"{ref} screw head", uid(site="pcb.py:main"),
            ((x0, y0), (x1, y0), (x1, y1), (x0, y1)), KEEPOUT_FLOORS,
            net=(0, ""), keepout="(tracks not_allowed) (vias not_allowed) (pads allowed)"
            " (copperpour not_allowed) (footprints allowed)").splitlines())
    return lines


def annotations(ctx, placement, foot_nodes):
    """Serialize the board's annotations."""
    unplaced, layout = ctx.unplaced, ctx.layout
    length, note, fixed, placements = placement.length, placement.note, placement.fixed, placement.positions
    lines = []
    lines.append(f'\t(gr_text {sexp.quote(note)}'
                 f' (at 4 -4 0) (layer "Cmts.User") (uuid "{uid(site="pcb.py:main")}") '
                 '(effects (font (size 2 2) (thickness 0.3)) (justify left bottom)))')
    if unplaced:
        usb_label_y = fmt(fixed["U_MCU"][1] - 0.2)
        lines.append(f'\t(gr_text {sexp.quote("USB PLUG KEEP-OUT")} (at 5.3 {usb_label_y} 90)'
                     f' (layer "Dwgs.User") (uuid "{uid(site="pcb.py:main")}") '
                     '(effects (font (size 0.8 0.8) (thickness 0.15))))')
        silk_marks = pin_silk_marks(ctx, fixed, length)
        for text, x, y, angle, layer in silk_marks:
            justify = " (justify mirror)" if layer == "B.SilkS" else ""
            lines.append(f'\t(gr_text {sexp.quote(text)} (at {fmt(x)} {fmt(y)} {angle})'
                         f' (layer "{layer}") (uuid "{uid(site="pcb.py:main")}") '
                         f'(effects (font (size 0.8 0.8) (thickness 0.15)){justify}))')
        for ref, label in (("J1", "PWR"), ("J2", "LED OUT"),
                           ("J3A", "SYNC IN"), ("J3B", "SYNC OUT")):
            x, y, rot = fixed[ref]
            if ctx.front_fixed_layout:
                pins = len(F(next(fn for fn in foot_nodes
                                  if connectivity.footprint_reference(fn) == ref), "pad"))
                label_x, label_y, angle = x - (pins - 1) * 1.27, 1.2, 0
            else:
                label_x = x + 3.2 if ref == "J1" else x - 2
                label_y = y + (1.27 if ref == "J1" else 2.54)
                angle = 90
            lines.append(f'\t(gr_text {sexp.quote(label)} '
                         f'(at {fmt(label_x)} {fmt(label_y)} {angle}) '
                         f'(layer "B.SilkS") (uuid "{uid(site="pcb.py:main")}") '
                         '(effects (font (size 0.8 0.8) (thickness 0.15)) (justify mirror)))')
    legend_x = placements["U_MCU"][0]
    for text, y, size in layout.back_silk:
        text_x = layout.back_silk_x.get(text, legend_x)
        if isinstance(text_x, str):
            text_x = fixed[text_x][0]
        lines.append(f'\t(gr_text {sexp.quote(text)} (at {fmt(text_x)} {fmt(y)} 0)'
                     f' (layer "B.SilkS") (uuid "{uid(site="pcb.py:main")}") '
                     f'(effects (font (size {fmt(size)} {fmt(size)})'
                     f' (thickness {fmt(max(0.15, size*0.15))})) (justify mirror)))')
    return lines


def footprints(ctx, foot_nodes):
    """Serialize the board's footprints."""
    unplaced = ctx.unplaced
    lines = []
    for fn in foot_nodes:
        lines.append(sexp.dumps(fn, indent=1))
    if unplaced:
        lines.extend(local_routes(ctx, foot_nodes))
    return lines


def board_text(ctx, netlist, placement, foot_nodes):
    """Serialize one complete board with its placed footprints."""
    lines = header(ctx)
    lines.extend(layers_setup())
    lines.extend(nets_outline(netlist, placement))
    lines.extend(ground_planes(netlist, placement))
    lines.extend(keepout_zones(placement))
    lines.extend(annotations(ctx, placement, foot_nodes))
    lines.extend(footprints(ctx, foot_nodes))
    lines.append(")")
    return "\n".join(lines) + "\n"


def zone_sexpr(layers, name, uuid, points, floors: ZoneFloors, *, net=None,
               clearance=0, keepout=None):
    """Serialize a rectangular pour or reservation with explicit feature floors."""
    lines = ["\t(zone"]
    if net is not None:
        lines.extend((f"\t\t(net {net[0]})", f"\t\t(net_name {sexp.quote(net[1])})"))
    layer_key = "layer" if len(layers) == 1 else "layers"
    lines.extend((f"\t\t({layer_key} " + " ".join(sexp.quote(layer) for layer in layers) + ")",
                  f'\t\t(uuid "{uuid}")', f"\t\t(name {sexp.quote(name)})",
                  f"\t\t(hatch {'full' if keepout else 'edge'} 0.5)",
                  f"\t\t(connect_pads (clearance {clearance:g}))",
                  f"\t\t(min_thickness {floors.min_thickness:g})"))
    if keepout:
        lines.extend((f"\t\t(keepout {keepout})",
                      '\t\t(placement (enabled no) (sheetname ""))'))
    lines.append(f"\t\t(fill{'' if keepout else ' yes'} (thermal_gap {floors.thermal_gap:g}) "
                 f"(thermal_bridge_width {floors.thermal_bridge_width:g}))")
    lines.append("\t\t(polygon (pts " + " ".join(f"(xy {fmt(x)} {fmt(y)})" for x, y in points)
                 + "))")
    lines.append("\t)")
    return "\n".join(lines) + "\n"
