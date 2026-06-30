"""Generate phantasm.kicad_pcb from the schematic + its netlist.

Embeds each component's footprint (from KiCad stock libs, or a generated Teensy
footprint), assigns pad nets by name from the exported netlist, places everything
linearly inside a <=30 mm-wide board outline (R-MECH-6), and declares the nets.
Placement is a rough starting arrangement; route/refine interactively in Pcbnew.
"""
import os
import subprocess
import sys
import uuid as _uuid
import sexp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
SCH = os.path.join(OUT, "phantasm.kicad_sch")
FP_DIR = os.environ.get("KICAD_FOOTPRINT_DIR",
                        r"C:\Program Files\KiCad\10.0\share\kicad\footprints")
KCLI = r"C:\Program Files\KiCad\10.0\bin\kicad-cli.exe"
PCB_W = 32.0  # board width (mm); within the R-MECH-6 cap (<=35), trimmed to part extent


def uid():
    return str(_uuid.uuid4())


def fmt(v):
    return str(int(v)) if float(v) == int(v) else f"{v:.4f}".rstrip("0").rstrip(".")


# ---------------------------------------------------------------- netlist
def export_netlist():
    net = os.path.join(OUT, "_pcb.net")
    subprocess.run([KCLI, "sch", "export", "netlist", "--format", "kicadsexpr",
                    "-o", net, SCH], check=True, capture_output=True)
    root = sexp.parse(open(net, encoding="utf-8").read())[0]
    os.remove(net)
    return root


def build_nets(nlroot):
    def F(n, k):
        return [c for c in n if isinstance(c, list) and c and c[0] == k]
    pad_net = {}            # (ref, pad) -> netname
    names = []
    for nb in F(nlroot, "nets"):
        for net in F(nb, "net"):
            name = sexp._val(net, "name")[0].lstrip("/")
            names.append(name)
            for nd in F(net, "node"):
                pad_net[(sexp._val(nd, "ref")[0], sexp._val(nd, "pin")[0])] = name
    # stable net ids: 0 = "", then unique names in first-seen order
    netid = {"": 0}
    for n in names:
        if n not in netid:
            netid[n] = len(netid)
    return pad_net, netid


def build_paths(nlroot):
    """ref -> footprint (path) linking to the schematic symbol, from netlist tstamps."""
    def F(n, k):
        return [c for c in n if isinstance(c, list) and c and c[0] == k]
    out = {}
    for cb in F(nlroot, "components"):
        for comp in F(cb, "comp"):
            ref = sexp._val(comp, "ref")[0]
            ts = sexp._val(comp, "tstamps")
            if ts:
                out[ref] = "/" + ts[0].strip("/")
    return out


# ---------------------------------------------------------------- components
def schematic_components():
    """Return ordered unique [(ref, footprint_libid, value)] (skip power/flag)."""
    root = sexp.parse(open(SCH, encoding="utf-8").read())[0]
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
            seen[ref] = (ref, fp or "", val or "")
            order.append(ref)
    return [seen[r] for r in order]


# ---------------------------------------------------------------- footprints
_MOD_CACHE = {}


def load_mod(libid):
    if libid not in _MOD_CACHE:
        lib, name = libid.split(":", 1)
        path = os.path.join(FP_DIR, lib + ".pretty", name + ".kicad_mod")
        _MOD_CACHE[libid] = sexp.parse(open(path, encoding="utf-8").read())[0]
    return _MOD_CACHE[libid]


def teensy_footprint():
    """2x14 THT module, pads named by Teensy pin label; long axis along X."""
    LEFT = ["GND", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"]
    RIGHT = ["VIN", "GND", "3V3", "23", "22", "21", "20", "19", "18", "17", "16", "15", "14", "13"]
    pads = []
    x0 = -((14 - 1) / 2) * 2.54
    for i, nm in enumerate(LEFT):
        pads.append(f'\t(pad "{nm}" thru_hole oval (at {fmt(x0 + i*2.54)} -7.62) '
                    f'(size 1.7 1.7) (drill 1.0) (layers "*.Cu" "*.Mask"))')
    for i, nm in enumerate(RIGHT):
        pads.append(f'\t(pad "{nm}" thru_hole oval (at {fmt(x0 + i*2.54)} 7.62) '
                    f'(size 1.7 1.7) (drill 1.0) (layers "*.Cu" "*.Mask"))')
    body = ('\t(fp_rect (start -18.5 -9.5) (end 18.5 9.5) (stroke (width 0.12) (type solid))'
            ' (fill none) (layer "F.SilkS"))\n'
            '\t(fp_rect (start -18.5 -9.5) (end 18.5 9.5) (stroke (width 0.05) (type solid))'
            ' (fill none) (layer "F.CrtYd"))')
    text = ('(footprint "phantasm:Teensy4.0"\n'
            '\t(layer "F.Cu")\n'
            '\t(descr "Teensy 4.0 module, 2x14 0.1in headers")\n'
            '\t(property "Reference" "U" (at 0 -11 0) (layer "F.SilkS")\n'
            '\t\t(effects (font (size 1 1) (thickness 0.15))))\n'
            '\t(property "Value" "Teensy4.0" (at 0 11 0) (layer "F.Fab")\n'
            '\t\t(effects (font (size 1 1) (thickness 0.15))))\n'
            '\t(attr through_hole)\n'
            + body + "\n" + "\n".join(pads) + "\n)\n")
    return sexp.parse(text)[0]


def embed(libid, ref, value, x, y, rot, pad_net, netid, path=None, locked=False):
    node = teensy_footprint() if libid in ("", "phantasm:Teensy4.0") else \
        sexp.parse(sexp.dumps(load_mod(libid)))[0]  # deep copy via round-trip
    fid = libid if libid else "phantasm:Teensy4.0"
    node[1] = fid
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
            elif c[1] == "Value":
                c[2] = value
    # assign pad nets by name
    for c in node:
        if isinstance(c, list) and c and c[0] == "pad":
            padname = c[1]
            nn = pad_net.get((ref, padname))
            if nn is not None:
                # remove any existing net, then add
                c[:] = [d for d in c if not (isinstance(d, list) and d and d[0] == "net")]
                c.append([sexp.Sym("net"), netid[nn], nn])
    return node


# ---------------------------------------------------------------- layout
# 2-D shelf strip-packer: minimise board LENGTH within the fixed WIDTH cap by
# stacking parts across the width in shelves (first-fit-decreasing). The Teensy
# (~37x19) sets the floor; small SMD parts pack into the leftover width beside it.
# Draft placement — refine orientation / push connectors to the edges in Pcbnew.
def fp_bbox(node):
    """Footprint bounding box (minx,miny,maxx,maxy) in its local (origin) frame,
    over pads + graphic outlines. Pad rotation is folded into a max-dim radius."""
    xs = []; ys = []
    for c in node:
        if not (isinstance(c, list) and c):
            continue
        if c[0] == "pad":
            at = sexp._val(c, "at"); sz = sexp._val(c, "size")
            x = float(at[0]); y = float(at[1])
            r = max(float(sz[0]), float(sz[1])) / 2 if sz else 0.5
            xs += [x - r, x + r]; ys += [y - r, y + r]
        elif c[0] in ("fp_rect", "fp_line", "fp_poly", "fp_circle"):
            for k in ("start", "end", "center"):
                v = sexp._val(c, k)
                if v:
                    xs.append(float(v[0])); ys.append(float(v[1]))
    if not xs:
        return (-1.0, -1.0, 1.0, 1.0)
    return (min(xs), min(ys), max(xs), max(ys))


def _rot_bb(bb, rot):
    mnx, mny, mxx, mxy = bb
    pts = [(mnx, mny), (mxx, mny), (mxx, mxy), (mnx, mxy)]
    if rot == 90:
        pts = [(-y, x) for x, y in pts]
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    return (min(xs), min(ys), max(xs), max(ys))


def _rotatable(ref, bb):
    # Placement aligns the rotated bbox corner to the cell, so any origin packs
    # correctly. Rotate chip passives (R/C) and pin-header connectors (J*). Keep
    # solder jumpers (JP*, custom pads + clearance-outline trip DRC at 90), diodes,
    # fuse, bead and ICs at rot 0.
    if ref.startswith("JP"):
        return False
    return ref[0] in "RC" or ref[0] == "J"


# Through-hole connectors are pinned to the board ends (not skyline-packed) so the
# mating cables are accessible at the edge, per the spec's signal flow: power/debug
# at the hub end, strip + sync daisy at the far end (R-CON-4). Each is stood with its
# pin-row across the width (rot 0 for these 1xN vertical headers), so it hugs the end
# edge and adds almost no length. (v1 opposite-edge split — best SI; the one-edge
# variant crowds the narrow edge and lengthens the fast nets.)
HUB_CONNS = ("J1", "J4")            # logic power in, debug — hub end (left)
FAR_CONNS = ("J2", "J3A", "J3B")    # strip signal, sync daisy in/out — far end (right)


def _stack(refs, bxs, x0, edge, gap, place):
    """Place `refs` as a column stacked along Y at left edge x0 (rot 0). Returns
    (column right-edge x, total stacked height)."""
    y = 0.0
    w = 0.0
    for ref in refs:
        rb = _rot_bb(bxs[ref], 0)
        place[ref] = (round(edge + x0 - rb[0], 3), round(edge + y - rb[1], 3), 0)
        y += (rb[3] - rb[1]) + gap
        w = max(w, rb[2] - rb[0])
    return x0 + w, y - gap


def pack(bxs, width, edge=1.0, gap=1.2):
    """Connectors pinned to the ends + skyline (bottom-left-fill) interior pack.
    Board WIDTH is fixed; minimise LENGTH. Returns ({ref:(x,y,rot)}, length mm)."""
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

    # 1) hub connectors: column at the left edge; reserve it so the interior packs right
    if hub:
        hubw, hubh = _stack(hub, bxs, 0.0, edge, gap, place)
        sky = reserve(0.0, hubw + gap, hubh + gap, 0.0)

    # 2) interior parts: Teensy + SMD + jumpers, bottom-left skyline pack
    interior = sorted((r for r in bxs if r not in pinned),
                      key=lambda r: -max(bxs[r][2] - bxs[r][0], bxs[r][3] - bxs[r][1]))
    for ref in interior:
        bb = bxs[ref]
        rots = (0, 90) if _rotatable(ref, bb) else (0,)
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
        sumh = sum(_rot_bb(bxs[r], 0)[3] - _rot_bb(bxs[r], 0)[1] for r in far)
        fgap = gap
        if len(far) > 1:
            fgap = max(0.6, min(gap, (usable - sumh) / (len(far) - 1)))
        right, _ = _stack(far, bxs, right + gap, edge, fgap, place)
    return place, round(edge + right + edge, 2)


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


def main(unplaced=False):
    nlroot = export_netlist()
    pad_net, netid = build_nets(nlroot)
    paths = build_paths(nlroot)                   # ref -> schematic-symbol path
    comps = {r: (r, fp, v) for r, fp, v in schematic_components()}
    # footprint bounding boxes -> 2-D shelf-pack to minimise length
    bxs = {}
    for ref, (_, fp, _) in comps.items():
        node = teensy_footprint() if fp in ("", "phantasm:Teensy4.0") else load_mod(fp)
        bxs[ref] = fp_bbox(node)
    PLACE, L = pack(bxs, PCB_W)
    conns = set(HUB_CONNS) | set(FAR_CONNS)
    # Ring-mount clearance (R-MECH): all through-hole parts (Teensy, the 0.1" headers)
    # and the tall electrolytic C_IN must stay on the TOP side so the board's back is
    # clear of bodies/the big cap. They're pre-placed + locked on F.Cu so Quilter can't
    # flip them to the bottom; only the low-profile SMD is staged below for Quilter.
    TOP_LOCKED = conns | {"U_MCU", "C_IN"}
    if unplaced:
        # top-locked parts pre-placed (locked) on F.Cu; low SMD staged below for Quilter
        staged = unplaced_layout(bxs, L, PCB_W)
        PLACE = {r: (PLACE[r] if r in TOP_LOCKED else staged[r]) for r in bxs}
        OUTFILE = os.path.join("unplaced", "phantasm_unplaced.kicad_pcb")
        NOTE = (f'PHANTASM segment board UNPLACED  -  {fmt(L)}x{fmt(PCB_W)}mm outline '
                '(width <=35mm); THT + tall C_IN locked on TOP, SMD staged below for Quilter')
    else:
        OUTFILE = "phantasm.kicad_pcb"
        NOTE = (f'PHANTASM segment board  -  {fmt(L)}x{fmt(PCB_W)}mm (width <=35mm, R-MECH-6); '
                'shelf-packed draft, route in Pcbnew')

    foot_nodes = []
    for ref, (x, y, rot) in PLACE.items():
        _, fp, val = comps[ref]
        lock = unplaced and ref in TOP_LOCKED
        foot_nodes.append(embed(fp, ref, val, x, y, rot, pad_net, netid,
                                path=paths.get(ref), locked=lock))

    lines = []
    lines.append("(kicad_pcb")
    lines.append("\t(version 20250513)")
    lines.append('\t(generator "phantasm-gen")')
    lines.append('\t(generator_version "10.0")')
    lines.append("\t(general (thickness 1.6) (legacy_teardrops no))")
    lines.append('\t(paper "A2")')
    lines.append("\t(layers")
    for n, nm, ty in [(0, "F.Cu", "signal"), (2, "B.Cu", "signal"),
                      (9, "F.Adhes", "user"), (11, "B.Adhes", "user"),
                      (13, "F.Paste", "user"), (15, "B.Paste", "user"),
                      (5, "F.SilkS", "user"), (7, "B.SilkS", "user"),
                      (1, "F.Mask", "user"), (3, "B.Mask", "user"),
                      (17, "Dwgs.User", "user"), (19, "Cmts.User", "user"),
                      (21, "Eco1.User", "user"), (23, "Eco2.User", "user"),
                      (25, "Edge.Cuts", "user"), (27, "Margin", "user"),
                      (31, "F.CrtYd", "user"), (29, "B.CrtYd", "user"),
                      (35, "F.Fab", "user"), (33, "B.Fab", "user")]:
        lines.append(f'\t\t({n} "{nm}" {ty})')
    lines.append("\t)")
    lines.append("\t(setup (pad_to_mask_clearance 0))")
    # nets
    for nm, i in sorted(netid.items(), key=lambda kv: kv[1]):
        esc = nm.replace("\\", "\\\\").replace('"', '\\"')
        lines.append(f'\t(net {i} "{esc}")')
    # board outline (Edge.Cuts) — <=35 mm wide strip, length minimised by packer
    lines.append(f'\t(gr_rect (start 0 0) (end {fmt(L)} {fmt(PCB_W)}) '
                 '(stroke (width 0.15) (type solid)) (fill none) (layer "Edge.Cuts") '
                 f'(uuid "{uid()}"))')
    esc_note = NOTE.replace("\\", "\\\\").replace('"', '\\"')
    lines.append(f'\t(gr_text "{esc_note}"'
                 f' (at 4 -4 0) (layer "Cmts.User") (uuid "{uid()}") '
                 '(effects (font (size 2 2) (thickness 0.3)) (justify left bottom)))')
    for fn in foot_nodes:
        lines.append(sexp.dumps(fn, indent=1))
    lines.append(")")
    outpath = os.path.join(OUT, OUTFILE)
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    open(outpath, "w", encoding="utf-8").write("\n".join(lines) + "\n")

    # --- custom footprint library (Teensy) + fp-lib-table ---
    pretty = os.path.join(OUT, "phantasm.pretty")
    os.makedirs(pretty, exist_ok=True)
    mod = teensy_footprint()
    mod[1] = "Teensy4.0"
    mod.insert(2, [sexp.Sym("version"), sexp.Sym("20240108")])
    mod.insert(3, [sexp.Sym("generator"), "phantasm-gen"])
    mod.insert(4, [sexp.Sym("generator_version"), "10.0"])
    open(os.path.join(pretty, "Teensy4.0.kicad_mod"), "w", encoding="utf-8").write(
        sexp.dumps(mod) + "\n")
    fplt = os.path.join(OUT, "fp-lib-table")
    if not os.path.exists(fplt):
        open(fplt, "w", encoding="utf-8").write(
            '(fp_lib_table\n\t(version 7)\n'
            '\t(lib (name "phantasm")(type "KiCad")(uri "${KIPRJMOD}/phantasm.pretty")'
            '(options "")(descr "PHANTASM custom footprints"))\n)\n')
    print(f"wrote {OUTFILE}  footprints:{len(foot_nodes)} nets:{len(netid)} length:{L:.0f}mm")


if __name__ == "__main__":
    main(unplaced="--unplaced" in sys.argv)
