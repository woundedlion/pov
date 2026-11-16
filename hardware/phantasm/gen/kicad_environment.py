"""KiCad generation prerequisites shared by CLI setup and tests."""

import os

import sexp
from kicad_common import kicad_cli


def pinned_kicad_cli_error():
    """Return the CLI pin failure, or None when the pinned executable resolves."""
    try:
        kicad_cli()
    except SystemExit as error:
        return str(error)
    return None


def generation_blockers():
    """Return missing stock libraries and CLI prerequisites."""
    blockers = []
    if not os.path.isdir(sexp.KICAD_SHARE):
        blockers.append(f"KiCad {sexp.KICAD_MAJOR} stock symbols not found; "
                        "set KICAD_SYMBOL_DIR")
    footprints = sexp.find_kicad_data_dir("footprints", "KICAD_FOOTPRINT_DIR")
    if not os.path.isdir(footprints):
        blockers.append(f"KiCad {sexp.KICAD_MAJOR} stock footprints not found; "
                        "set KICAD_FOOTPRINT_DIR")
    cli_error = pinned_kicad_cli_error()
    if cli_error:
        blockers.append(cli_error)
    return blockers
