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


def uid():
    return str(_uuid.uuid4())


def fmt(v):
    return str(int(v)) if float(v) == int(v) else f"{v:.4f}".rstrip("0").rstrip(".")


def F(n, k):
    return [c for c in n if isinstance(c, list) and c and c[0] == k]


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
