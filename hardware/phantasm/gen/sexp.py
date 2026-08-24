"""Minimal KiCad S-expression parser + symbol-library helpers.

Just enough to: load a .kicad_sym, pull a top-level (symbol "Name" ...) block,
serialize a node back to text, and enumerate a symbol's pins with their
*local* connection coordinates (the pin tip you wire to).
"""
import glob
import os
import re
from dataclasses import dataclass

def kicad_version_key(path):
    """(major, minor) of a KiCad install path; (0, 0) if unversioned."""
    m = re.search(r"KiCad[\\/](\d+)\.(\d+)", path)
    return (int(m.group(1)), int(m.group(2))) if m else (0, 0)


# The KiCad major every committed board, schematic and fab output was generated
# with. Pinned in tools/build_pins.py, which asserts this spelling.
KICAD_MAJOR = 10

# The (generator_version) string stamped into every file the generators emit.
GENERATOR_VERSION = f"{KICAD_MAJOR}.0"

# (version) epochs, one per file format. KiCad revises each format on its own
# schedule, so these differ from each other and cannot be derived from
# KICAD_MAJOR; each is the epoch KICAD_MAJOR reads and writes for that format.
SCHEMATIC_FORMAT = "20251024"
PCB_FORMAT = "20250513"
FOOTPRINT_FORMAT = "20240108"


def kicad_data_dir_patterns(kind):
    """Stock install globs searched for a KiCad data directory, in order."""
    return [
        fr"C:\Program Files\KiCad\*\share\kicad\{kind}",
        fr"C:\Program Files (x86)\KiCad\*\share\kicad\{kind}",
        f"/Applications/KiCad/KiCad.app/Contents/SharedSupport/{kind}",
        f"/usr/share/kicad/{kind}",
        f"/usr/local/share/kicad/{kind}",
    ]


# Stock KiCad data directories (symbols, footprints, 3dmodels). Override with
# the matching env var if installed elsewhere or on a newer/older KiCad version.
def find_kicad_data_dir(kind, env_name):
    """Locate a stock KICAD_MAJOR data directory, falling back to the bare `kind`.

    An install whose path names another major is skipped: its land patterns and
    symbol graphics differ from what the pinned kicad-cli validates, so drawing
    from one would embed foreign geometry into a board that CLI then passes.
    Paths carrying no version (the Unix prefixes) name no major and are taken.
    The env override is an explicit choice and is not version-checked.

    The fallback is deliberately not a directory: callers probe the result with
    isdir() to skip the KiCad-dependent work rather than fail at import.
    """
    env = os.environ.get(env_name)
    if env and os.path.isdir(env):
        return env
    for pattern in kicad_data_dir_patterns(kind):
        hits = [hit for hit in glob.glob(pattern)
                if kicad_version_key(hit)[0] in (0, KICAD_MAJOR)]
        if hits:
            return max(hits, key=kicad_version_key)   # newest minor
    return kind


KICAD_SHARE = find_kicad_data_dir("symbols", "KICAD_SYMBOL_DIR")


# ---------- tokenizer / parser ----------
class Sym(str):
    """An unquoted atom (token), distinct from a quoted string."""


class SList(list):
    """List carrying the source layout used to parse it."""

    def __init__(self, values=(), layout=None):
        super().__init__(values)
        self.layout = layout


@dataclass(frozen=True)
class _Token:
    kind: str
    value: str
    start: int
    end: int


@dataclass(frozen=True)
class _Layout:
    prefix: str
    separators: tuple[str, ...]
    suffix: str
    atoms: tuple[tuple[type, str, str] | None, ...]
    indent: str


ESCAPE_DECODE = {
    "\\": "\\",
    '"': '"',
    "a": "\a",
    "b": "\b",
    "f": "\f",
    "n": "\n",
    "r": "\r",
    "t": "\t",
    "v": "\v",
}
ESCAPE_ENCODE = {value: key for key, value in ESCAPE_DECODE.items()}


def tokenize(s):
    toks, i, n = [], 0, len(s)
    while i < n:
        c = s[i]
        if c in " \t\r\n":
            i += 1
        elif c == "(" or c == ")":
            toks.append(_Token(c, c, i, i + 1)); i += 1
        elif c == '"':
            j = i + 1; buf = bytearray()
            while j < n:
                if s[j] == "\\":
                    if j + 1 >= n:
                        raise ValueError("trailing backslash in string")
                    escaped = s[j + 1]
                    if escaped in ESCAPE_DECODE:
                        buf.extend(ESCAPE_DECODE[escaped].encode()); j += 2
                    elif escaped == "x":
                        digits = re.match(r"[0-9a-fA-F]{1,2}", s[j + 2:])
                        if digits:
                            value = digits.group(0)
                            buf.append(int(value, 16))
                            j += 2 + len(value)
                        else:
                            buf.append(ord("x")); j += 2
                    elif escaped in "01234567":
                        value = re.match(r"[0-7]{1,3}", s[j + 1:]).group(0)
                        buf.append(int(value, 8) & 0xff)
                        j += 1 + len(value)
                    else:
                        buf.append(ord("\\")); j += 1
                elif s[j] == '"':
                    break
                else:
                    buf.extend(s[j].encode()); j += 1
            if j == n:
                raise ValueError("unterminated string")
            toks.append(_Token("STR", buf.decode(), i, j + 1)); i = j + 1
        else:
            j = i
            while j < n and s[j] not in ' \t\r\n()"':
                j += 1
            toks.append(_Token("ATOM", s[i:j], i, j)); i = j
    return toks


def parse(s):
    toks = tokenize(s)
    pos = [0]

    def rd():
        if pos[0] >= len(toks):
            raise ValueError(f"unexpected end of input at token {pos[0]}")
        token = toks[pos[0]]; pos[0] += 1
        if token.kind == "(":
            values = []
            spans = []
            while True:
                if pos[0] >= len(toks):
                    raise ValueError(f"unexpected end of input at token {pos[0]}")
                if toks[pos[0]].kind == ")":
                    break
                child, start, end = rd()
                values.append(child)
                spans.append((start, end))
            close = toks[pos[0]]
            pos[0] += 1
            if spans:
                prefix = s[token.end:spans[0][0]]
                separators = tuple(s[left[1]:right[0]]
                                   for left, right in zip(spans, spans[1:]))
                suffix = s[spans[-1][1]:close.start]
            else:
                prefix = s[token.end:close.start]
                separators = ()
                suffix = ""
            line_start = s.rfind("\n", 0, token.start) + 1
            line_prefix = s[line_start:token.start]
            indent = re.match(r"[ \t]*", line_prefix).group(0)
            atoms = tuple(
                None if isinstance(value, list)
                else (type(value), value, s[start:end])
                for value, (start, end) in zip(values, spans)
            )
            value = SList(values, _Layout(prefix, separators, suffix,
                                          atoms, indent))
            return value, token.start, close.end
        if token.kind == ")":
            raise ValueError("unexpected )")
        value = token.value if token.kind == "STR" else Sym(token.value)
        return value, token.start, token.end

    out = []
    while pos[0] < len(toks):
        out.append(rd()[0])
    return out


def dumps(node, indent=0):
    if isinstance(node, SList) and node.layout is not None:
        rendered = _preserved_list(node)
        return _reindent(rendered, node.layout.indent, "\t" * indent)
    return _canonical_dumps(node, indent)


def _canonical_dumps(node, indent=0):
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
                parts.append(_canonical_dumps(child, indent + 1))
            else:
                parts[-1] += " " + _atom(child)
        parts.append(pad + ")")
        return "\n".join(parts)
    return pad + _atom(node)


def _preserved_list(node):
    layout = node.layout
    if len(node) != len(layout.atoms):
        canonical_indent = (layout.indent
                            if set(layout.indent) <= {"\t"} else "")
        level = len(canonical_indent)
        rendered = _canonical_dumps(node, level)
        return rendered[len(canonical_indent):]
    parts = ["(", layout.prefix]
    for index, child in enumerate(node):
        if index:
            parts.append(layout.separators[index - 1])
        original = layout.atoms[index]
        if isinstance(child, SList) and child.layout is not None:
            parts.append(_preserved_list(child))
        elif isinstance(child, list):
            parts.append(_canonical_dumps(child))
        elif (original is not None and type(child) is original[0]
              and child == original[1]):
            parts.append(original[2])
        else:
            parts.append(_atom(child))
    parts.extend((layout.suffix, ")"))
    return "".join(parts)


def _reindent(text, original, desired):
    if original == desired:
        return desired + text
    lines = text.split("\n")
    for index in range(1, len(lines)):
        if lines[index].startswith(original):
            lines[index] = lines[index][len(original):]
        lines[index] = desired + lines[index]
    return desired + "\n".join(lines)


def quote(s):
    """Render `s` as a quoted atom with KiCad-compatible control escapes."""
    esc = "".join(_escaped_char(c) for c in s)
    return '"' + esc + '"'


def _escaped_char(c):
    if c in ESCAPE_ENCODE:
        return "\\" + ESCAPE_ENCODE[c]
    if ord(c) < 0x20 or ord(c) == 0x7f:
        return f"\\x{ord(c):02x}"
    return c


def _atom(x):
    if isinstance(x, Sym):
        return str(x)
    if isinstance(x, str):
        return quote(x)
    return str(x)


# ---------- library helpers ----------
_LIB_CACHE = {}


def load_lib(libname):
    if libname not in _LIB_CACHE:
        if not os.path.isdir(KICAD_SHARE):
            raise RuntimeError(
                "KiCad stock symbol libraries not found; set KICAD_SYMBOL_DIR "
                "to a share/kicad/symbols directory, or install KiCad in one "
                "of: " + ", ".join(kicad_data_dir_patterns("symbols")))
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


def val(node, key, default=None):
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
                    at = val(c, "at")
                    length = val(c, "length", [0])[0]
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
