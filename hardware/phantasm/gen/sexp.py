"""Minimal KiCad S-expression parser + symbol-library helpers.

Just enough to: load a .kicad_sym, pull a top-level (symbol "Name" ...) block,
serialize a node back to text, and enumerate a symbol's pins with their
*local* connection coordinates (the pin tip you wire to).
"""
import os

# Stock KiCad symbol libraries. Override with env KICAD_SYMBOL_DIR if installed
# elsewhere or on a newer/older KiCad version.
KICAD_SHARE = os.environ.get(
    "KICAD_SYMBOL_DIR", r"C:\Program Files\KiCad\10.0\share\kicad\symbols")


# ---------- tokenizer / parser ----------
class Sym(str):
    """An unquoted atom (token), distinct from a quoted string."""


def tokenize(s):
    toks, i, n = [], 0, len(s)
    while i < n:
        c = s[i]
        if c in " \t\r\n":
            i += 1
        elif c == "(" or c == ")":
            toks.append(c); i += 1
        elif c == '"':
            j = i + 1; buf = []
            while j < n:
                if s[j] == "\\":
                    buf.append(s[j + 1]); j += 2
                elif s[j] == '"':
                    break
                else:
                    buf.append(s[j]); j += 1
            toks.append(("STR", "".join(buf))); i = j + 1
        else:
            j = i
            while j < n and s[j] not in ' \t\r\n()"':
                j += 1
            toks.append(("ATOM", s[i:j])); i = j
    return toks


def parse(s):
    toks = tokenize(s)
    pos = [0]

    def rd():
        t = toks[pos[0]]; pos[0] += 1
        if t == "(":
            lst = []
            while toks[pos[0]] != ")":
                lst.append(rd())
            pos[0] += 1
            return lst
        if t == ")":
            raise ValueError("unexpected )")
        if isinstance(t, tuple):
            return t[1] if t[0] == "STR" else Sym(t[1])
        raise ValueError(t)

    out = []
    while pos[0] < len(toks):
        out.append(rd())
    return out


def dumps(node, indent=0):
    pad = "\t" * indent
    if isinstance(node, list):
        if not node:
            return pad + "()"
        head = node[0]
        # keep simple short lists on one line
        if all(not isinstance(x, list) for x in node):
            return pad + "(" + " ".join(_atom(x) for x in node) + ")"
        parts = [pad + "(" + _atom(head)]
        for child in node[1:]:
            if isinstance(child, list):
                parts.append(dumps(child, indent + 1))
            else:
                parts[-1] += " " + _atom(child)
        parts.append(pad + ")")
        return "\n".join(parts)
    return pad + _atom(node)


def _atom(x):
    if isinstance(x, Sym):
        return str(x)
    if isinstance(x, str):
        esc = x.replace("\\", "\\\\").replace('"', '\\"')
        return '"' + esc + '"'
    return str(x)


# ---------- library helpers ----------
_LIB_CACHE = {}


def load_lib(libname):
    if libname not in _LIB_CACHE:
        path = os.path.join(KICAD_SHARE, libname + ".kicad_sym")
        with open(path, "r", encoding="utf-8") as f:
            _LIB_CACHE[libname] = parse(f.read())[0]
    return _LIB_CACHE[libname]


def get_symbol(libname, name):
    root = load_lib(libname)
    for node in root:
        if isinstance(node, list) and node and node[0] == "symbol" and node[1] == name:
            return node
    raise KeyError(f"{libname}:{name} not found")


def _find(node, key):
    return [c for c in node if isinstance(c, list) and c and c[0] == key]


def _val(node, key, default=None):
    for c in node:
        if isinstance(c, list) and c and c[0] == key:
            return c[1:]
    return default


def pins(sym):
    """Return list of dicts: number, name, x, y, angle, length (local coords)."""
    out = []

    def walk(node):
        for c in node:
            if isinstance(c, list) and c:
                if c[0] == "pin":
                    at = _val(c, "at")
                    length = _val(c, "length", [0])[0]
                    num = None; name = None
                    for d in c:
                        if isinstance(d, list) and d and d[0] == "number":
                            num = d[1]
                        if isinstance(d, list) and d and d[0] == "name":
                            name = d[1]
                    out.append(dict(number=str(num), name=str(name),
                                    x=float(at[0]), y=float(at[1]),
                                    angle=float(at[2]) if len(at) > 2 else 0.0,
                                    length=float(length)))
                else:
                    walk(c)
    walk(sym)
    return out


if __name__ == "__main__":
    import sys
    lib, name = sys.argv[1], sys.argv[2]
    sym = get_symbol(lib, name)
    print(f"# {lib}:{name}")
    for p in pins(sym):
        print(f"  pin {p['number']:>3} {p['name']:<8} at ({p['x']:>7.3f},{p['y']:>7.3f}) "
              f"ang={p['angle']:>5.0f} len={p['length']}")
