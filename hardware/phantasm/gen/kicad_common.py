"""Shared KiCad-gen helpers: kicad-cli discovery + pure s-expr/format/geometry utilities.

Single source of truth for finding kicad-cli and for the netlist-export call,
plus the small helpers reused across pcb.py / check.py / shorts.py / builder.py
/ board_metadata.py, so a KiCad-flag or schema change touches one place.
"""
import glob
import math
import os
import re
import subprocess
import sys
import tempfile
import uuid as _uuid
import sexp


# The KiCad major every committed board, schematic and fab output was generated
# with. Pinned in tools/build_pins.py, which asserts this spelling.
KICAD_MAJOR = 10

KICAD_CLI_PATTERNS = (
    r"C:\Program Files\KiCad\*\bin\kicad-cli.exe",
    r"C:\Program Files (x86)\KiCad\*\bin\kicad-cli.exe",
    "/Applications/KiCad/KiCad.app/Contents/MacOS/kicad-cli",
    "/usr/bin/kicad-cli",
    "/usr/local/bin/kicad-cli",
)


def kicad_cli_major(path):
    """Major version of the kicad-cli at `path`, or None if unreadable.

    Stock Windows installs carry the version in the path; the Unix ones do not,
    so those are asked directly.
    """
    major = sexp.kicad_version_key(path)[0]
    if major:
        return major
    try:
        reported = subprocess.run([path, "--version"], capture_output=True,
                                  text=True, check=True).stdout
    except (OSError, subprocess.SubprocessError):
        return None
    found = re.search(r"(\d+)\.\d+", reported)
    return int(found.group(1)) if found else None


def find_kicad_cli():
    """Path to kicad-cli: $KICAD_CLI, else a KICAD_MAJOR install, else the PATH name.

    Exits when installs are present but none is KICAD_MAJOR: a newer KiCad
    upgrades the board on open and formats the fab outputs differently.
    """
    env = os.environ.get("KICAD_CLI")
    if env and os.path.exists(env):
        return env
    hits = [hit for pattern in KICAD_CLI_PATTERNS for hit in glob.glob(pattern)]
    pinned = [hit for hit in hits if kicad_cli_major(hit) == KICAD_MAJOR]
    if pinned:
        return max(pinned, key=sexp.kicad_version_key)
    if hits:
        sys.exit(f"the fab gates are pinned to KiCad {KICAD_MAJOR}, which is not "
                 "installed\n"
                 f"  found: {', '.join(hits)}\n"
                 f"  Install KiCad {KICAD_MAJOR} or set KICAD_CLI to its "
                 "kicad-cli.")
    return "kicad-cli"                 # assume on PATH


UID_NAMESPACE = _uuid.UUID("73e63115-f1ba-4688-a7e5-5a689c53d6fd")
_uid_seq = {}


def reset_uid_sequence():
    _uid_seq.clear()


def uid():
    """UUID for a generated element, derived from the call site and occurrence.

    Regenerating unchanged inputs yields a byte-identical file. The source file
    is part of the key, so the schematic and board generators draw from
    disjoint id spaces.
    """
    code = sys._getframe(1).f_code
    site = f"{os.path.basename(code.co_filename)}:{code.co_name}"
    n = _uid_seq[site] = _uid_seq.get(site, 0) + 1
    return str(_uuid.uuid5(UID_NAMESPACE, f"{site}#{n}"))


def fmt(v):
    return str(int(v)) if float(v) == int(v) else f"{v:.4f}".rstrip("0").rstrip(".")


def F(n, k):
    return [c for c in n if isinstance(c, list) and c and c[0] == k]


def is_copper_pour(zone):
    """True for a zone that pours copper.

    A rule area carries a keepout node and pours nothing. Its net is not the
    discriminator: gen/pcb.py writes keepouts with (net 0).
    """
    return not F(zone, "keepout")


def arc_extrema(start, mid, end, collinear_error=None):
    """Axis-aligned extreme points of the arc through three (x, y) floats.

    Solves the circumcircle, then keeps whichever axis crossings (0, 90, 180,
    270 degrees) lie on the start->mid->end sweep. Collinear inputs have no
    circumcircle: raise `collinear_error` if given, else return no extrema.
    """
    x1, y1 = start
    x2, y2 = mid
    x3, y3 = end
    denominator = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
    if abs(denominator) < 1e-12:
        if collinear_error is not None:
            raise collinear_error
        return []

    center_x = (
        (x1 * x1 + y1 * y1) * (y2 - y3)
        + (x2 * x2 + y2 * y2) * (y3 - y1)
        + (x3 * x3 + y3 * y3) * (y1 - y2)
    ) / denominator
    center_y = (
        (x1 * x1 + y1 * y1) * (x3 - x2)
        + (x2 * x2 + y2 * y2) * (x1 - x3)
        + (x3 * x3 + y3 * y3) * (x2 - x1)
    ) / denominator
    radius = math.hypot(x1 - center_x, y1 - center_y)
    angles = tuple(math.atan2(y - center_y, x - center_x) for x, y in (start, mid, end))

    def ccw_between(angle, first, last):
        return (angle - first) % math.tau <= (last - first) % math.tau + 1e-12

    first, middle, last = angles
    counterclockwise = ccw_between(middle, first, last)
    points = []
    for angle in (0, math.pi / 2, math.pi, 3 * math.pi / 2):
        on_arc = ccw_between(angle, first, last) if counterclockwise else \
            ccw_between(angle, last, first)
        if on_arc:
            points.append((center_x + radius * math.cos(angle),
                           center_y + radius * math.sin(angle)))
    return points


DEFAULT_OVERWRITE_REASON = (
    "It holds KiCad-authored work these generators do not reproduce\n"
    "  (routing, vias, silk, hand edits).")


def require_writable(path, force, reason=DEFAULT_OVERWRITE_REASON, flag="--force"):
    """Exit unless `path` is absent or `force` is set.

    `reason` states what the overwrite would cost, defaulting to the
    KiCad-authored routing and hand edits a committed board or schematic holds.
    `flag` names the option that authorizes this particular overwrite.
    """
    if force or not os.path.exists(path):
        return
    sys.exit(f"refusing to overwrite {path}\n"
             f"  {reason}\n"
             f"  Re-run with {flag} to regenerate it from scratch anyway.")


def export_netlist(kcli, sch):
    """Export `sch` to a kicadsexpr netlist via kicad-cli; return its parsed root.

    Exits with a diagnostic when kicad-cli is absent or the export fails: every
    caller is a command-line gate, for which a traceback says less.
    """
    fd, net = tempfile.mkstemp(suffix=".net")
    os.close(fd)
    try:
        subprocess.run([kcli, "sch", "export", "netlist", "--format", "kicadsexpr",
                        "-o", net, sch], check=True, capture_output=True, text=True)
        with open(net, encoding="utf-8") as fh:
            return sexp.parse(fh.read())[0]
    except FileNotFoundError:
        sys.exit(f"kicad-cli not found: {kcli}; install KiCad, put kicad-cli on "
                 "PATH, or set KICAD_CLI to its full path")
    except subprocess.CalledProcessError as e:
        sys.stderr.write(e.stderr or "")
        sys.exit(f"netlist export failed: kicad-cli exited {e.returncode} on {sch}")
    finally:
        if os.path.exists(net):
            os.remove(net)
