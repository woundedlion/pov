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


DEFAULT_OVERWRITE_REASON = (
    "It holds KiCad-authored work these generators do not reproduce\n"
    "  (routing, vias, silk, hand edits).")


def require_writable(path, force, reason=DEFAULT_OVERWRITE_REASON):
    """Exit unless `path` is absent or `force` is set.

    `reason` states what the overwrite would cost, defaulting to the
    KiCad-authored routing and hand edits a committed board or schematic holds.
    """
    if force or not os.path.exists(path):
        return
    sys.exit(f"refusing to overwrite {path}\n"
             f"  {reason}\n"
             "  Re-run with --force to regenerate it from scratch anyway.")


def export_netlist(kcli, sch):
    """Export `sch` to a kicadsexpr netlist via kicad-cli; return its parsed root."""
    fd, net = tempfile.mkstemp(suffix=".net")
    os.close(fd)
    try:
        subprocess.run([kcli, "sch", "export", "netlist", "--format", "kicadsexpr",
                        "-o", net, sch], check=True, capture_output=True, text=True)
        with open(net, encoding="utf-8") as fh:
            return sexp.parse(fh.read())[0]
    except subprocess.CalledProcessError as e:
        sys.stderr.write(e.stderr or "")
        raise
    finally:
        if os.path.exists(net):
            os.remove(net)
