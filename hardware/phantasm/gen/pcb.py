"""Generate a placed draft or Quilter board from the schematic + its netlist.

Embeds footprints from KiCad stock libraries, hand-maintained terminal blocks
(TBC05-02/03 in the default revision's phantasm.pretty, TBC05-04 in gen/footprints),
or the generated Teensy footprint. Assigns pad nets by name
from the exported netlist and declares the nets. Draft mode skyline-packs
footprints into a PCB_W-wide outline (R-MECH-6). --unplaced locks fixed
parts on the QUILTER_LENGTH outline; rev 1.3 fixes the complete front placement,
while rev 1.2 stages its remaining parts below the outline for the autoplacer.
Emits a 4-layer SIG/GND/GND/SIG board: the physical stackup and the inner GND
planes are encoded in the file, so an autoplacer/fab reads them on upload.
Draft placement is a rough starting arrangement; route/refine it in Pcbnew.
"""
import argparse
import os
import shutil
import sys
import pcb_model
import pcb_layouts
import board as schematic_generator
import builder
import check
import sexp
from constraints import GROUND_PLANE_LAYERS as GROUND_PLANE_LAYERS
from constraints import NEW_LAYOUT_RULES as NEW_LAYOUT_RULES
from kicad_common import atomic_write_text
from kicad_common import (reset_uid_sequence, fmt, export_netlist, kicad_cli, require_writable)
from pcb_model import BOARD_THICKNESS_MM as BOARD_THICKNESS_MM
from pcb_model import COPPER_LAYERS as COPPER_LAYERS
from pcb_model import COURTYARD_LAYERS as COURTYARD_LAYERS
from pcb_model import DRAFT_FILE as DRAFT_FILE
from pcb_model import FAR_CONNS as FAR_CONNS
from pcb_footprints import FP_DIR as FP_DIR
from pcb_model import GENERATED_LAND_DEVIATIONS as GENERATED_LAND_DEVIATIONS
from pcb_model import GENERATED_OPTIONAL_DEVIATIONS as GENERATED_OPTIONAL_DEVIATIONS
from pcb_model import GROUND_NET as GROUND_NET
from pcb_model import HUB_CONNS as HUB_CONNS
from pcb_footprints import LOCAL_FOOTPRINT_DIR as LOCAL_FOOTPRINT_DIR
from pcb_model import MOUNTING_COURTYARD_RADIUS as MOUNTING_COURTYARD_RADIUS
from pcb_model import MOUNTING_HOLE_INSET as MOUNTING_HOLE_INSET
from pcb_model import MOUNTING_KEEPOUT_HALF_SIDE as MOUNTING_KEEPOUT_HALF_SIDE
from pcb_model import MOUNTING_RESERVE_HALF_SIDE as MOUNTING_RESERVE_HALF_SIDE
from pcb_model import PCB_W as PCB_W
from pcb_model import PCB_W_MAX as PCB_W_MAX
from pcb_model import QUILTER_LENGTH as QUILTER_LENGTH
from pcb_layouts import REVISION_LAYOUTS as REVISION_LAYOUTS
from pcb_model import STACKUP as STACKUP
from pcb_model import TECHNICAL_LAYERS as TECHNICAL_LAYERS
from pcb_model import TEENSY_LIBRARY_REASON as TEENSY_LIBRARY_REASON
from pcb_layouts import TERMINAL_EDGE_PLACEMENTS_1_3 as TERMINAL_EDGE_PLACEMENTS_1_3
from pcb_model import TERMINAL_LIBIDS as TERMINAL_LIBIDS
from pcb_model import UNPLACED_FILE as UNPLACED_FILE
from pcb_model import UNPLACED_REASON as UNPLACED_REASON
from pcb_model import Netlist, Placement, Boxes
from pcb_layouts import generation_context as generation_context
from pcb_footprints import embed_components
from pcb_emit import board_text
from pcb_pack import _boxes_overlap as _boxes_overlap
from pcb_pack import _rot_bb as _rot_bb
from pcb_pack import _rotatable as _rotatable
from pcb_model import build_nets as build_nets
from pcb_model import build_paths as build_paths
from pcb_model import copper_layer_names as copper_layer_names
from pcb_footprints import embed as embed
from pcb_footprints import embedded_footprint as embedded_footprint
from pcb_layouts import fixed_placements as fixed_placements
from pcb_footprints import fp_bbox as fp_bbox
from pcb_pack import keepout_clashes as keepout_clashes
from pcb_pack import keepout_rects as keepout_rects
from pcb_footprints import load_mod as load_mod
from pcb_pack import mounting_reserve_rects as mounting_reserve_rects
from pcb_pack import outline_overflows as outline_overflows
from pcb_pack import pack as pack
from pcb_footprints import reference_bbox as reference_bbox
from pcb_footprints import reserve_terminal_wire_access as reserve_terminal_wire_access
from pcb_model import schematic_components as schematic_components
from pcb_model import silk_revision as silk_revision
from pcb_footprints import teensy_footprint as teensy_footprint
from pcb_layouts import terminal_edge_distances as terminal_edge_distances
from pcb_footprints import terminal_library_dir as terminal_library_dir
from pcb_footprints import trim_terminal_silkscreen as trim_terminal_silkscreen
from pcb_pack import unplaced_layout as unplaced_layout
from pcb_footprints import validate_front_access as validate_front_access


def place_components(ctx, netlist):
    """Pack or stage the netlist and measure its copper and courtyard bounds."""
    comps, layout, unplaced = netlist.components, ctx.layout, ctx.unplaced
    bxs = {}
    pad_bxs = {}
    crt_bxs = {}
    for ref, (_, fp, _, _) in comps.items():
        node = embedded_footprint(ref, fp, ctx)
        bxs[ref] = reference_bbox(node, ref)
        pad_bxs[ref] = fp_bbox(node, pads_only=True)
        crt_bxs[ref] = fp_bbox(node, graphic_layers=COURTYARD_LAYERS)
    if unplaced:
        length = QUILTER_LENGTH
        fixed = pcb_layouts.fixed_placements(comps, ctx)
        edge_map = layout.edge_placements
        missing_edges = sorted(ref for ref, placement in edge_map.items()
                               if fixed.get(ref) != placement or ref not in comps or
                               comps[ref][1] != layout.terminal_by_ref[ref])
        if missing_edges:
            sys.exit("ERROR connectors require verified edge placements: "
                     + ", ".join(missing_edges))
        staged = unplaced_layout(bxs, length, PCB_W)
        placements = {r: fixed.get(r, staged[r]) for r in bxs}
        # the staged grid sits below the outline by design; only locked parts are on it
        bounded = list(fixed)
        outfile = UNPLACED_FILE
        status = 'UNROUTED' if ctx.front_fixed_layout else 'UNPLACED'
        note = (f'PHANTASM segment board {status}  -  {fmt(length)}x{fmt(PCB_W)}mm outline '
                f'(width <={fmt(PCB_W_MAX)}mm); mechanical and signal-integrity '
                'placements locked')
    else:
        placements, length = pack(bxs, PCB_W)
        bounded = list(bxs)
        outfile = DRAFT_FILE
        note = (f'PHANTASM segment board  -  {fmt(length)}x{fmt(PCB_W)}mm '
                f'(width <={fmt(PCB_W_MAX)}mm, R-MECH-6); '
                'skyline-packed draft, route in Pcbnew')
    return (Placement(placements, length, bounded, outfile, note, fixed if unplaced else {}),
            Boxes(pad_bxs, crt_bxs))


def validate_placements(placement, boxes):
    """Reject copper outside the outline and courtyards inside reservations."""
    placements, length, bounded = placement.positions, placement.length, placement.bounded
    pad_bxs, crt_bxs = boxes.pads, boxes.courtyards
    outside = outline_overflows(placements, pad_bxs, length, bounded)
    if outside:
        sys.exit(f"ERROR placements outside the {fmt(length)}x{fmt(PCB_W)}mm outline: "
                 + ", ".join(outside))

    clashes = keepout_clashes(placements, crt_bxs, length)
    if clashes:
        sys.exit("ERROR footprints inside a mounting-hole reservation: "
                 + ", ".join(clashes))


def prepare_teensy_library(ctx, force_teensy_library):
    """Prepare the module library and enforce its overwrite policy."""
    out = ctx.out_dir
    pretty = os.path.join(out, "phantasm.pretty")
    mod = teensy_footprint(ctx)
    mod[1] = "Teensy4.0"
    mod.insert(2, [sexp.Sym("version"), sexp.Sym(sexp.FOOTPRINT_FORMAT)])
    mod.insert(3, [sexp.Sym("generator"), "phantasm-gen"])
    mod.insert(4, [sexp.Sym("generator_version"), sexp.GENERATOR_VERSION])
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
    return pretty, mod_path, mod_text, existing_mod_text


def write_footprint_library(ctx, pretty, mod_path, mod_text, existing_mod_text):
    """Write the generated project's local footprint libraries."""
    out, layout = ctx.out_dir, ctx.layout
    os.makedirs(pretty, exist_ok=True)
    if existing_mod_text != mod_text:
        atomic_write_text(mod_path, mod_text)
    edge_distances = terminal_edge_distances(layout) if ctx.front_fixed_layout else {}
    for libid in layout.terminal_libids:
        name = libid.split(":", 1)[1] + ".kicad_mod"
        source_dir = terminal_library_dir(libid)
        source = os.path.join(source_dir, name)
        destination = os.path.join(pretty, name)
        if ctx.front_fixed_layout:
            with open(source, encoding="utf-8") as f:
                terminal = sexp.parse_one(f.read())
            reserve_terminal_wire_access(terminal)
            trim_terminal_silkscreen(terminal, edge_distances[libid])
            atomic_write_text(destination, sexp.dumps(terminal) + "\n")
        elif os.path.abspath(source) != os.path.abspath(destination):
            with open(source, encoding="utf-8") as f:
                atomic_write_text(destination, f.read())
    model_source = os.path.join(LOCAL_FOOTPRINT_DIR, "Teensy4.0.wrl")
    model_target = os.path.join(pretty, "Teensy4.0.wrl")
    if os.path.abspath(model_source) != os.path.abspath(model_target):
        shutil.copyfile(model_source, model_target)
    fplt = os.path.join(out, "fp-lib-table")
    if not os.path.exists(fplt):
        atomic_write_text(fplt, '(fp_lib_table\n\t(version 7)\n'
                '\t(lib (name "phantasm")(type "KiCad")(uri "${KIPRJMOD}/phantasm.pretty")'
                '(options "")(descr "PHANTASM custom footprints"))\n)\n')


def main(unplaced=False, force=False, force_teensy_library=False,
         revision=builder.REVISION, output_dir=None):
    """Generate the selected revision's draft or fixed-placement project."""
    out, sch = schematic_generator.revision_paths(revision, output_dir)
    ctx = generation_context(revision, unplaced, out)
    selected = ctx.revision
    if ctx.front_fixed_layout and not unplaced:
        sys.exit("ERROR revision 1.3 has a fixed placement and requires --unplaced")
    out = ctx.out_dir
    reset_uid_sequence()
    if unplaced:
        require_writable(os.path.join(out, UNPLACED_FILE), force, UNPLACED_REASON)
    else:
        require_writable(os.path.join(out, DRAFT_FILE), force)
    nlroot = export_netlist(kicad_cli(), sch)
    revision = check.netlist_revision(nlroot)
    if revision != selected:
        sys.exit(f"ERROR schematic revision {revision} does not match generator "
                 f"revision {selected}; regenerate the schematic first")
    if not check.check(check.netlist_nets(nlroot), revision):
        sys.exit("ERROR schematic does not match the electrical specification")
    pad_net, netid = pcb_model.build_nets(nlroot)
    paths = build_paths(nlroot)                   # ref -> schematic-symbol path
    comps = {r: (r, fp, v, dnp) for r, fp, v, dnp in pcb_model.schematic_components(sch)}
    netlist = Netlist(comps, pad_net, netid, paths)
    placement, boxes = place_components(ctx, netlist)
    validate_placements(placement, boxes)
    foot_nodes = embed_components(ctx, netlist, placement)
    text = board_text(ctx, netlist, placement, foot_nodes)
    pretty, mod_path, mod_text, existing_mod_text = prepare_teensy_library(
        ctx, force_teensy_library)
    outpath = os.path.join(out, placement.outfile)
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    atomic_write_text(outpath, text)
    project_path = os.path.splitext(outpath)[0] + ".kicad_pro"
    with open(sch, encoding="utf-8") as f:
        root_uuid = sexp.val(sexp.parse_one(f.read()), "uuid", [""])[0]
    schematic_generator.write_project(project_path, root_uuid, unplaced=unplaced)

    write_footprint_library(ctx, pretty, mod_path, mod_text, existing_mod_text)
    print(f"wrote {placement.outfile}  footprints:{len(foot_nodes)} nets:{len(netid)} length:{fmt(placement.length)}mm")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--revision", choices=schematic_generator.GENERATION_REVISIONS, default=builder.REVISION)
    parser.add_argument("--unplaced", action="store_true",
                        help="write the revisioned unplaced project for "
                             "the autoplacer instead")
    parser.add_argument("--force", action="store_true",
                        help=f"overwrite {DRAFT_FILE} / "
                             f"{UNPLACED_FILE}")
    parser.add_argument("--force-teensy-library", action="store_true",
                        help="overwrite phantasm.pretty/Teensy4.0.kicad_mod "
                             "when the generated footprint differs; this revision's "
                             "committed phantasm.kicad_pcb resolves its pads against it")
    return parser.parse_args(argv)


if __name__ == "__main__":
    ARGS = parse_args()
    main(unplaced=ARGS.unplaced, force=ARGS.force,
         force_teensy_library=ARGS.force_teensy_library, revision=ARGS.revision)
