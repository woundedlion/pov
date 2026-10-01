"""Restore routed copper-zone feature sizes to the generator's defaults.

Writes a separate board with cached fills removed. Refill all zones in KiCad,
save, and rerun DRC and the fabrication gates before using the output.
"""
import argparse
import math
from pathlib import Path
import sys

import fab
import sexp
from kicad_common import F, atomic_write_text, is_copper_pour


ZONE_DEFAULTS = {
    "min_thickness": 0.25,
    "thermal_gap": 0.5,
    "thermal_bridge_width": 0.5,
}
FILL_CACHE = {"filled_polygon", "fill_segments"}


def single_field(node, name):
    fields = F(node, name)
    if len(fields) != 1:
        raise ValueError(f"expected exactly one {name}, found {len(fields)}")
    return fields[0]


def heal_zones(source):
    """Return repaired board text and changed-zone count; preserve routed copper."""
    root = sexp.parse_one(source)
    if not isinstance(root, list) or not root or root[0] != "kicad_pcb":
        raise ValueError("expected a kicad_pcb document")
    zones = [zone for zone in F(root, "zone") if is_copper_pour(zone)]
    changed = 0
    for zone in zones:
        fill = single_field(zone, "fill")
        zone_changed = False
        for name, minimum in ZONE_DEFAULTS.items():
            field = single_field(zone if name == "min_thickness" else fill, name)
            try:
                value = float(field[1]) if len(field) == 2 else float("nan")
            except (TypeError, ValueError):
                value = float("nan")
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be one finite positive number")
            if value < minimum:
                field[1] = sexp.Sym(f"{minimum:g}")
                zone_changed = True
        changed += zone_changed
    fab.validate_zone_geometry("<repaired board>", board=root)
    if not changed:
        return source, 0
    for zone in zones:
        zone[:] = [child for child in zone
                   if not isinstance(child, list) or not child
                   or child[0] not in FILL_CACHE]
    prefix = source[:len(source) - len(source.lstrip())]
    suffix = source[len(source.rstrip()):]
    repaired = prefix + sexp.dumps(root) + suffix
    if "\r\n" in source and "\n" not in source.replace("\r\n", ""):
        repaired = repaired.replace("\r\n", "\n").replace("\n", "\r\n")
    return repaired, changed


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("board", type=Path, help="routed .kicad_pcb to inspect")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true",
                      help="check fabrication zone geometry without writing")
    mode.add_argument("-o", "--output", type=Path, help="separate repaired board")
    parser.add_argument("--force", action="store_true",
                        help="replace an existing output, never the source board")
    args = parser.parse_args(argv)
    try:
        if args.check:
            fab.validate_zone_geometry(args.board)
            print("Zone settings pass; refill and DRC are separate required checks.")
            return 0
        source_path = args.board.resolve()
        output_path = args.output.resolve()
        if source_path == output_path or (
                output_path.exists() and source_path.samefile(output_path)):
            raise ValueError("output must be separate from the source board")
        if (output_path.parent / "SHA256SUMS.txt").exists():
            raise ValueError("output directory is hash-manifested")
        if output_path.exists() and not args.force:
            raise ValueError("output exists; use --force to replace it")
        with source_path.open(encoding="utf-8", newline="") as handle:
            source = handle.read()
        repaired, changed = heal_zones(source)
        atomic_write_text(output_path, repaired, newline="")
        print(f"Wrote {output_path}: repaired {changed} copper zone(s).")
        if changed:
            print("Cached copper fills removed. Refill all zones in KiCad, save, "
                  "then rerun DRC and fabrication gates.")
        return 0
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
