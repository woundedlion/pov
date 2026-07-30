"""Shared KiCad-gen helpers: netlist export + pure s-expr/format utilities.

Single source of truth for the kicad-cli netlist-export call and the small
helpers reused across pcb.py / check.py / shorts.py / builder.py, so a
KiCad-flag or schema change touches one place.
"""
import os
import subprocess
import sys
import tempfile
import uuid as _uuid
import sexp


UID_NAMESPACE = _uuid.UUID("73e63115-f1ba-4688-a7e5-5a689c53d6fd")
_uid_seq = {}


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


def require_writable(path, force):
    """Exit unless `path` is absent or `force` is set.

    The committed board and schematic carry KiCad-authored routing, silk and
    hand edits that these generators do not reproduce.
    """
    if force or not os.path.exists(path):
        return
    sys.exit(f"refusing to overwrite {path}\n"
             "  It holds KiCad-authored work these generators do not reproduce\n"
             "  (routing, vias, silk, hand edits). Re-run with --force to\n"
             "  regenerate it from scratch and discard that work.")


def export_netlist(kcli, sch):
    """Export `sch` to a kicadsexpr netlist via kicad-cli; return its parsed root."""
    fd, net = tempfile.mkstemp(suffix=".net")
    os.close(fd)
    try:
        subprocess.run([kcli, "sch", "export", "netlist", "--format", "kicadsexpr",
                        "-o", net, sch], check=True, capture_output=True, text=True)
        return sexp.parse(open(net, encoding="utf-8").read())[0]
    except subprocess.CalledProcessError as e:
        sys.stderr.write(e.stderr or "")
        raise
    finally:
        if os.path.exists(net):
            os.remove(net)
