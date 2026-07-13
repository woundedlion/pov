"""Extract routed-board manufacturing facts from a KiCad PCB."""

import argparse
from collections import Counter
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path
import sys

import sexp


FACTS_START = "<!-- BEGIN ROUTED PCB FACTS -->"
FACTS_END = "<!-- END ROUTED PCB FACTS -->"


class MetadataError(ValueError):
    pass


@dataclass(frozen=True)
class BoardMetadata:
    width_mm: Decimal
    height_mm: Decimal
    thickness_mm: Decimal
    footprint_sides: tuple[tuple[str, int], ...]
    track_segments: int
    vias: int
    zones: int
    zone_layers: tuple[tuple[str, int], ...]
    copper_layers: tuple[str, ...]
    copper_stackup: tuple[tuple[str, Decimal], ...]
    copper_finish: str


def _children(node, key):
    return [
        child for child in node
        if isinstance(child, list) and child and str(child[0]) == key
    ]


def _one_child(node, key):
    children = _children(node, key)
    if len(children) != 1:
        raise MetadataError(f"expected one {key}, found {len(children)}")
    return children[0]


def _one_value(node, key):
    child = _one_child(node, key)
    if len(child) != 2 or isinstance(child[1], list):
        raise MetadataError(f"expected one value for {key}")
    return str(child[1])


def _decimal(value, field):
    try:
        return Decimal(str(value))
    except InvalidOperation as error:
        raise MetadataError(f"invalid decimal for {field}: {value}") from error


def _point(node, key):
    child = _one_child(node, key)
    if len(child) < 3:
        raise MetadataError(f"expected x/y coordinates for {key}")
    return _decimal(child[1], key), _decimal(child[2], key)


def _outline_bounds(root):
    points = []
    supported = {"gr_line", "gr_rect", "gr_poly"}
    for node in root[1:]:
        if not isinstance(node, list) or not node:
            continue
        key = str(node[0])
        if not key.startswith("gr_") or not _children(node, "layer"):
            continue
        if _one_value(node, "layer") != "Edge.Cuts":
            continue
        if key not in supported:
            raise MetadataError(f"unsupported Edge.Cuts primitive: {key}")
        if key in {"gr_line", "gr_rect"}:
            points.extend((_point(node, "start"), _point(node, "end")))
            continue
        pts = _one_child(node, "pts")
        polygon_points = [_point([child], "xy") for child in _children(pts, "xy")]
        if not polygon_points:
            raise MetadataError("Edge.Cuts polygon has no points")
        points.extend(polygon_points)
    if not points:
        raise MetadataError("board has no supported Edge.Cuts geometry")
    xs = [point[0] for point in points]
    ys = [point[1] for point in points]
    width = max(xs) - min(xs)
    height = max(ys) - min(ys)
    if width <= 0 or height <= 0:
        raise MetadataError("board outline has non-positive dimensions")
    return width, height


def parse_board(text):
    try:
        documents = sexp.parse(text)
    except (IndexError, ValueError) as error:
        raise MetadataError(f"invalid KiCad S-expression: {error}") from error
    if len(documents) != 1:
        raise MetadataError(f"expected one document, found {len(documents)}")
    root = documents[0]
    if not isinstance(root, list) or not root or str(root[0]) != "kicad_pcb":
        raise MetadataError("document is not a kicad_pcb")

    width, height = _outline_bounds(root)
    general = _one_child(root, "general")
    thickness = _decimal(_one_value(general, "thickness"), "thickness")

    layers_node = _one_child(root, "layers")
    copper_layers = tuple(
        str(layer[1]) for layer in layers_node[1:]
        if isinstance(layer, list) and len(layer) >= 3 and str(layer[1]).endswith(".Cu")
    )
    if not copper_layers:
        raise MetadataError("board has no copper layers")

    setup = _one_child(root, "setup")
    stackup = _one_child(setup, "stackup")
    copper_stackup = []
    for layer in _children(stackup, "layer"):
        if _one_value(layer, "type") == "copper":
            if len(layer) < 2 or isinstance(layer[1], list):
                raise MetadataError("stackup copper layer has no name")
            copper_stackup.append((
                str(layer[1]),
                _decimal(_one_value(layer, "thickness"), "copper thickness"),
            ))
    if tuple(name for name, _ in copper_stackup) != copper_layers:
        raise MetadataError("declared copper layers do not match the stackup")

    footprint_counts = Counter({"F.Cu": 0, "B.Cu": 0})
    footprints = _children(root, "footprint")
    for footprint in footprints:
        side = _one_value(footprint, "layer")
        if side not in footprint_counts:
            raise MetadataError(f"footprint is on unsupported layer: {side}")
        footprint_counts[side] += 1

    zone_counts = Counter()
    zones = _children(root, "zone")
    for zone in zones:
        layer_nodes = _children(zone, "layer")
        layers_nodes = _children(zone, "layers")
        if len(layer_nodes) + len(layers_nodes) != 1:
            raise MetadataError("zone must declare layer or layers exactly once")
        if layer_nodes:
            zone_layers = [str(value) for value in layer_nodes[0][1:]]
        else:
            zone_layers = [str(value) for value in layers_nodes[0][1:]]
        if not zone_layers or any(layer not in copper_layers for layer in zone_layers):
            raise MetadataError("zone has an invalid copper layer")
        zone_counts.update(zone_layers)

    return BoardMetadata(
        width_mm=width,
        height_mm=height,
        thickness_mm=thickness,
        footprint_sides=tuple((side, footprint_counts[side]) for side in ("F.Cu", "B.Cu")),
        track_segments=len(_children(root, "segment")),
        vias=len(_children(root, "via")),
        zones=len(zones),
        zone_layers=tuple((layer, zone_counts[layer]) for layer in copper_layers if zone_counts[layer]),
        copper_layers=copper_layers,
        copper_stackup=tuple(copper_stackup),
        copper_finish=_one_value(stackup, "copper_finish"),
    )


def load_board(path):
    try:
        return parse_board(path.read_text(encoding="utf-8"))
    except OSError as error:
        raise MetadataError(f"cannot read {path}: {error}") from error


def _format_decimal(value):
    formatted = format(value, "f").rstrip("0").rstrip(".")
    return formatted or "0"


def render_facts(metadata):
    dimensions = f"{_format_decimal(metadata.width_mm)} × {_format_decimal(metadata.height_mm)} mm"
    footprint_sides = ", ".join(f"{side}: {count}" for side, count in metadata.footprint_sides)
    zone_layers = ", ".join(f"{layer}: {count}" for layer, count in metadata.zone_layers)
    copper_stackup = "; ".join(
        f"{layer}: {_format_decimal(thickness)} mm"
        for layer, thickness in metadata.copper_stackup
    )
    rows = [
        ("Board dimensions", dimensions),
        ("Board thickness", f"{_format_decimal(metadata.thickness_mm)} mm"),
        ("Footprints by side", f"{sum(count for _, count in metadata.footprint_sides)} ({footprint_sides})"),
        ("Track segments", str(metadata.track_segments)),
        ("Vias", str(metadata.vias)),
        ("Copper zones", f"{metadata.zones} ({zone_layers})"),
        ("Copper layers", f"{len(metadata.copper_layers)} ({', '.join(metadata.copper_layers)})"),
        ("Copper thicknesses", copper_stackup),
        ("Copper finish", metadata.copper_finish),
    ]
    lines = [
        FACTS_START,
        "<!-- Generated by `python gen/board_metadata.py --write-readme`; do not edit. -->",
        "| Routed-board fact | Extracted value |",
        "|---|---|",
    ]
    lines.extend(f"| {name} | {value} |" for name, value in rows)
    lines.append(FACTS_END)
    return "\n".join(lines)


def replace_facts(readme, facts):
    if readme.count(FACTS_START) != 1 or readme.count(FACTS_END) != 1:
        raise MetadataError("README must contain exactly one routed PCB facts block")
    start = readme.index(FACTS_START)
    end = readme.index(FACTS_END, start) + len(FACTS_END)
    return readme[:start] + facts + readme[end:]


def check_facts(readme, facts):
    if replace_facts(readme, facts) != readme:
        raise MetadataError(
            "README routed PCB facts are stale; run "
            "python hardware/phantasm/gen/board_metadata.py --write-readme"
        )


def main(argv=None):
    base = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--board", type=Path, default=base / "phantasm.kicad_pcb")
    parser.add_argument("--readme", type=Path, default=base / "README.md")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write-readme", action="store_true")
    args = parser.parse_args(argv)

    try:
        facts = render_facts(load_board(args.board))
        if args.check or args.write_readme:
            readme = args.readme.read_text(encoding="utf-8")
            if args.check:
                check_facts(readme, facts)
            else:
                args.readme.write_text(replace_facts(readme, facts), encoding="utf-8")
        else:
            print(facts)
    except (MetadataError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
