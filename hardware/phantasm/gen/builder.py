"""KiCad 10 schematic builder: place stock symbols, draw wires/labels, emit .kicad_sch."""
import copy
import math
import uuid as _uuid
import sexp

VERSION = "20250114"


def uid():
    return str(_uuid.uuid4())


# ---------- transform: library local coords -> schematic screen coords ----------
def transform(sx, sy, rot, mirror, lx, ly):
    """Map a symbol-local pin coord to its global schematic coord.

    Verified empirically against kicad-cli netlist export (see calib.py).
    """
    # Library Y is up; schematic Y is down -> flip Y first.
    x, y = lx, -ly
    # mirror is applied in the symbol (screen) frame, before rotation.
    if mirror == "x":      # mirror across X axis (flip vertically)
        y = -y
    elif mirror == "y":    # mirror across Y axis (flip horizontally)
        x = -x
    # rotation: KiCad rotates counter-clockwise on screen (y-down frame).
    a = math.radians(rot)
    c, s = math.cos(a), math.sin(a)
    rx = x * c + y * s
    ry = -x * s + y * c
    return (round(sx + rx, 4), round(sy + ry, 4))


class Symbol:
    def __init__(self, lib_id, ref, value, x, y, rot=0, mirror=None,
                 unit=1, footprint="", dnp=False, fields=None, hide_value=False):
        self.lib_id = lib_id
        self.ref = ref
        self.value = value
        # snap placement to the 1.27 mm (50 mil) connection grid so pins land on grid
        self.x = round(x / 1.27) * 1.27
        self.y = round(y / 1.27) * 1.27
        self.rot = rot
        self.mirror = mirror
        self.unit = unit
        self.footprint = footprint
        self.dnp = dnp
        self.fields = fields or {}
        self.hide_value = hide_value
        self.uuid = uid()
        self._pins = None  # filled by builder

    def pin(self, number):
        """Global (x,y) of pin `number` for this instance's unit."""
        p = self._pins[str(number)]
        return transform(self.x, self.y, self.rot, self.mirror, p["x"], p["y"])

    def pin_dir(self, number):
        """Outward unit axis (dx,dy) in screen coords for a stub.

        Derived from the pin's own orientation (angle points tip->body, so
        outward = angle+180), mapped through the verified placement transform.
        Robust for tall symbols where distance-to-center would mislead.
        """
        p = self._pins[str(number)]
        th = math.radians(p["angle"] + 180.0)
        ox, oy = math.cos(th), math.sin(th)
        sx, sy = transform(0, 0, self.rot, self.mirror, ox, oy)
        if abs(sx) >= abs(sy):
            return (1 if sx > 0 else -1, 0)
        return (0, 1 if sy > 0 else -1)


class Builder:
    def __init__(self, title, paper="A3"):
        self.title = title
        self.paper = paper
        self.uuid = uid()
        self.symbols = []
        self.wires = []
        self.labels = []
        self.junctions = []
        self.no_connects = []
        self.texts = []          # (point, string, size)
        self.lib_defs = {}       # lib_id -> symbol node (for lib_symbols)
        self._unit_pins = {}     # lib_id -> {unit -> {num -> pin}}

    # ---- library symbol management ----
    def _resolve(self, lib, name):
        """Return a fully-resolved symbol node (extends flattened)."""
        node = copy.deepcopy(sexp.get_symbol(lib, name))
        ext = sexp._val(node, "extends")
        if ext:
            base = self._resolve(lib, ext[0])
            base = copy.deepcopy(base)
            # rename base sub-symbols from base-root to this name
            _rename_subsymbols(base, ext[0], name)
            base[1] = name
            # overlay this node's properties onto base
            _overlay_props(base, node)
            node = base
        return node

    def ensure_lib(self, lib, name, lib_id=None, value=None, footprint=None,
                   reference=None):
        lib_id = lib_id or f"{lib}:{name}"
        if lib_id not in self.lib_defs:
            node = self._resolve(lib, name)
            node[1] = lib_id  # top name becomes the full lib_id
            self.lib_defs[lib_id] = node
            # index per-unit pins
            self._unit_pins[lib_id] = _index_unit_pins(node)
        return lib_id

    def register_custom(self, node, lib_id):
        """Register an already-built custom symbol node under lib_id."""
        node = copy.deepcopy(node)
        node[1] = lib_id
        self.lib_defs[lib_id] = node
        self._unit_pins[lib_id] = _index_unit_pins(node)
        return lib_id

    # ---- placement / drawing ----
    def place(self, sym: Symbol):
        sym._pins = self._unit_pins[sym.lib_id].get(sym.unit, {})
        # include common (unit 0) pins
        common = self._unit_pins[sym.lib_id].get(0, {})
        merged = dict(common)
        merged.update(sym._pins)
        sym._pins = merged
        self.symbols.append(sym)
        return sym

    def wire(self, p1, p2):
        # orthogonal: if not aligned, route as an L through a corner
        x1, y1 = p1
        x2, y2 = p2
        if x1 == x2 or y1 == y2:
            self._seg(p1, p2)
        else:
            corner = (x2, y1)
            self._seg(p1, corner)
            self._seg(corner, p2)
            self.junctions.append((x2, y1))  # may be pruned later
        return p2

    def _seg(self, p1, p2):
        if p1 != p2:
            self.wires.append((p1, p2))

    def stub(self, sym, number, length=2.54):
        """Draw a stub from a pin outward; return the stub end point."""
        tip = sym.pin(number)
        dx, dy = sym.pin_dir(number)
        end = (round(tip[0] + dx * length, 4), round(tip[1] + dy * length, 4))
        self._seg(tip, end)
        return end

    def label(self, p, text, angle=0):
        self.labels.append((p, text, angle))

    def text(self, p, s, size=2.0):
        self.texts.append((p, s, size))

    def junction(self, p):
        self.junctions.append(p)

    def nc(self, p):
        self.no_connects.append(p)

    # ---- emit ----
    def dumps(self):
        out = []
        out.append(f'(kicad_sch')
        out.append(f'\t(version {VERSION})')
        out.append(f'\t(generator "eeschema")')
        out.append(f'\t(generator_version "10.0")')
        out.append(f'\t(uuid "{self.uuid}")')
        out.append(f'\t(paper "{self.paper}")')
        out.append('\t(title_block')
        out.append(f'\t\t(title "{self.title}")')
        out.append('\t\t(rev "A")')
        out.append('\t)')
        # lib_symbols
        out.append('\t(lib_symbols')
        for lib_id in sorted(self.lib_defs):
            out.append(sexp.dumps(self.lib_defs[lib_id], indent=2))
        out.append('\t)')
        # junctions (dedup)
        for p in _dedup(self.junctions):
            out.append('\t(junction')
            out.append(f'\t\t(at {fmt(p[0])} {fmt(p[1])})')
            out.append('\t\t(diameter 0)')
            out.append('\t\t(color 0 0 0 0)')
            out.append(f'\t\t(uuid "{uid()}")')
            out.append('\t)')
        # no connects
        for p in self.no_connects:
            out.append('\t(no_connect')
            out.append(f'\t\t(at {fmt(p[0])} {fmt(p[1])})')
            out.append(f'\t\t(uuid "{uid()}")')
            out.append('\t)')
        # wires
        for (p1, p2) in self.wires:
            out.append('\t(wire')
            out.append(f'\t\t(pts (xy {fmt(p1[0])} {fmt(p1[1])}) (xy {fmt(p2[0])} {fmt(p2[1])}))')
            out.append('\t\t(stroke (width 0) (type default))')
            out.append(f'\t\t(uuid "{uid()}")')
            out.append('\t)')
        # labels
        for (p, text, angle) in self.labels:
            just = "left" if angle in (0, 90) else "right"
            out.append(f'\t(label "{text}"')
            out.append(f'\t\t(at {fmt(p[0])} {fmt(p[1])} {angle})')
            out.append(f'\t\t(effects (font (size 1.27 1.27)) (justify {just} bottom))')
            out.append(f'\t\t(uuid "{uid()}")')
            out.append('\t)')
        # text notes (block headers)
        for (p, s, size) in self.texts:
            esc = s.replace("\\", "\\\\").replace('"', '\\"')
            out.append(f'\t(text "{esc}"')
            out.append(f'\t\t(at {fmt(p[0])} {fmt(p[1])} 0)')
            out.append(f'\t\t(effects (font (size {size} {size}) (bold yes)) (justify left bottom))')
            out.append(f'\t\t(uuid "{uid()}")')
            out.append('\t)')
        # symbols
        for s in self.symbols:
            out.append(self._dump_symbol(s))
        out.append('\t(sheet_instances')
        out.append('\t\t(path "/" (page "1"))')
        out.append('\t)')
        out.append('\t(embedded_fonts no)')
        out.append(')')
        return "\n".join(out) + "\n"

    def _dump_symbol(self, s):
        L = []
        L.append('\t(symbol')
        L.append(f'\t\t(lib_id "{s.lib_id}")')
        rot = s.rot
        L.append(f'\t\t(at {fmt(s.x)} {fmt(s.y)} {rot})')
        if s.mirror:
            L.append(f'\t\t(mirror {s.mirror})')
        L.append(f'\t\t(unit {s.unit})')
        L.append('\t\t(exclude_from_sim no)')
        L.append('\t\t(in_bom yes)')
        L.append('\t\t(on_board yes)')
        L.append(f'\t\t(dnp {"yes" if s.dnp else "no"})')
        L.append(f'\t\t(uuid "{s.uuid}")')
        # properties
        refy = s.y - 5.08
        valy = s.y + 5.08
        L += _prop("Reference", s.ref, s.x, refy)
        L += _prop("Value", s.value, s.x, valy, hide=s.hide_value)
        L += _prop("Footprint", s.footprint, s.x, s.y + 7.62, hide=True)
        L += _prop("Datasheet", "", s.x, s.y, hide=True)
        for k, v in s.fields.items():
            L += _prop(k, v, s.x, s.y, hide=True)
        # pin uuids
        for num in self._unit_pins[s.lib_id].get(s.unit, {}):
            L.append(f'\t\t(pin "{num}" (uuid "{uid()}"))')
        for num in self._unit_pins[s.lib_id].get(0, {}):
            L.append(f'\t\t(pin "{num}" (uuid "{uid()}"))')
        # instances
        L.append('\t\t(instances')
        L.append('\t\t\t(project "phantasm"')
        L.append(f'\t\t\t\t(path "/{self.uuid}"')
        L.append(f'\t\t\t\t\t(reference "{s.ref}") (unit {s.unit})')
        L.append('\t\t\t\t)')
        L.append('\t\t\t)')
        L.append('\t\t)')
        L.append('\t)')
        return "\n".join(L)


def _prop(name, value, x, y, hide=False, angle=0):
    v = value.replace('\\', '\\\\').replace('"', '\\"')
    out = [f'\t\t(property "{name}" "{v}"',
           f'\t\t\t(at {fmt(x)} {fmt(y)} {angle})',
           '\t\t\t(effects (font (size 1.27 1.27))' + (' (hide yes)' if hide else '') + ')',
           '\t\t)']
    return out


def fmt(v):
    if v == int(v):
        return str(int(v))
    return f"{v:.4f}".rstrip("0").rstrip(".")


def _dedup(pts):
    seen = set(); out = []
    for p in pts:
        k = (round(p[0], 3), round(p[1], 3))
        if k not in seen:
            seen.add(k); out.append(p)
    return out


# ---------- lib symbol helpers ----------
def _rename_subsymbols(node, old_root, new_root):
    for c in node:
        if isinstance(c, list) and c and c[0] == "symbol" and isinstance(c[1], str):
            if c[1].startswith(old_root + "_"):
                c[1] = new_root + c[1][len(old_root):]


def _overlay_props(base, node):
    # replace Value/Footprint/Reference/Datasheet/Description from node
    wanted = {}
    for c in node:
        if isinstance(c, list) and c and c[0] == "property":
            wanted[c[1]] = c
    newchildren = []
    seen = set()
    for c in base:
        if isinstance(c, list) and c and c[0] == "property" and c[1] in wanted:
            newchildren.append(copy.deepcopy(wanted[c[1]]))
            seen.add(c[1])
        else:
            newchildren.append(c)
    base[:] = newchildren


def _index_unit_pins(node):
    """Return {unit -> {number -> {x,y,angle,length,name}}}."""
    res = {}

    def walk(n):
        for c in n:
            if isinstance(c, list) and c:
                if c[0] == "symbol" and isinstance(c[1], str) and "_" in c[1]:
                    parts = c[1].rsplit("_", 2)
                    try:
                        unit = int(parts[-2])
                    except ValueError:
                        unit = 1
                    for d in c:
                        if isinstance(d, list) and d and d[0] == "pin":
                            at = sexp._val(d, "at")
                            num = None; name = None
                            for e in d:
                                if isinstance(e, list) and e and e[0] == "number":
                                    num = e[1]
                                if isinstance(e, list) and e and e[0] == "name":
                                    name = e[1]
                            res.setdefault(unit, {})[str(num)] = dict(
                                x=float(at[0]), y=float(at[1]),
                                angle=float(at[2]) if len(at) > 2 else 0.0,
                                name=str(name))
                else:
                    walk(c)
    walk(node)
    return res
