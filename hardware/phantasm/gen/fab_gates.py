"""KiCad board, project and DRC report validation without fabrication exports."""

import json
import math

import sexp
from connectivity import copper_layers, footprint_reference, zone_layers
from constraints import (DEFAULT_CLASS_MINIMUMS, GROUND_PLANE_LAYERS,
                         MIN_SOLDER_MASK_WEB_MM, MIN_THERMAL_GAP_MM,
                         MIN_THERMAL_SPOKE_MM, MIN_VIA_TO_VIA_COPPER_SPACING_MM,
                         RULE_MINIMUMS, rule_checks)
from kicad_common import F, is_copper_pour, net_name, net_name_resolver

MIN_STANDARD_VIA_DIAMETER_MM = RULE_MINIMUMS["min_via_diameter"]


MIN_STANDARD_VIA_DRILL_MM = DEFAULT_CLASS_MINIMUMS["via_drill"]


MIN_ZONE_GAP_MM = MIN_THERMAL_GAP_MM


MIN_ZONE_WIDTH_MM = RULE_MINIMUMS["min_track_width"]


ZONE_FILL_FEATURES = ("thermal_gap", "thermal_bridge_width")


ZONE_FEATURE_MINIMUMS = {
    "min_thickness": MIN_ZONE_WIDTH_MM,
    "thermal_gap": MIN_ZONE_GAP_MM,
    "thermal_bridge_width": MIN_THERMAL_SPOKE_MM,
}


MIN_BOARD_VIAS = 99


MIN_COPPER_POURS = len(GROUND_PLANE_LAYERS)


class DesignRuleError(ValueError):
    pass


def read_drc_report(report_path):
    """Return (violations, unconnected_items) from a kicad-cli JSON DRC report.

    Raises unless both sections are present as lists, so a report whose shape
    changed reads as an error rather than as zero violations.
    """
    try:
        with open(report_path, encoding="utf-8") as fh:
            report = json.load(fh)
    except (OSError, UnicodeError, ValueError) as exc:
        raise DesignRuleError(
            f"cannot read DRC report: {report_path}") from exc

    if not isinstance(report, dict) or any(
            not isinstance(report.get(key), list)
            for key in ("violations", "unconnected_items")):
        raise DesignRuleError(
            "DRC report has no violations/unconnected_items sections: "
            f"{report_path}")
    return report["violations"], report["unconnected_items"]


class ProjectRulesError(ValueError):
    pass


def validate_project_rules(project_path):
    """Require the project's DRC floors to meet the fabrication minimums.

    kicad-cli reads its constraints from the .kicad_pro, and KiCad re-zeroes
    min_clearance every time the project is opened in the GUI, so a DRC run
    can report clean against a disabled floor.
    """
    try:
        with open(project_path, encoding="utf-8") as project_file:
            document = json.load(project_file)
    except (OSError, UnicodeError, ValueError) as exc:
        raise ProjectRulesError(
            f"cannot read project file: {project_path}") from exc

    try:
        checks = rule_checks(document, RULE_MINIMUMS, DEFAULT_CLASS_MINIMUMS)
    except ValueError as exc:
        raise ProjectRulesError(f"{project_path}: {exc}") from exc

    shortfalls = sorted((field, current, required)
                        for field, (current, required, met) in checks.items() if not met)
    if shortfalls:
        summary = "; ".join(
            f"{field} is {current!r}, requires {required!r}"
            for field, current, required in shortfalls)
        raise ProjectRulesError(
            f"{project_path} would compute DRC under relaxed constraints: "
            f"{summary}. Run gen/heal_clearance.py on a working copy of the project, then retry.")
    return len(checks)


class ViaGeometryError(ValueError):
    pass


class ZoneGeometryError(ValueError):
    pass


class PlotOriginError(ValueError):
    pass


class BoardReadError(ValueError):
    pass


def read_board(pcb_path, error=BoardReadError, what="board"):
    """Parse a KiCad board; raise `error` naming `what` when it cannot be read."""
    try:
        with open(pcb_path, encoding="utf-8") as fh:
            return sexp.parse_one(fh.read())
    except (OSError, UnicodeError, ValueError) as exc:
        raise error(f"cannot read {what}: {pcb_path}") from exc


class SolderMaskError(ValueError):
    pass


def validate_solder_mask(pcb_path, board=None):
    """Require the corrected board's mask web and tented via defaults."""
    root = board if board is not None else read_board(
        pcb_path, SolderMaskError, "PCB solder mask")
    setups = F(root, "setup")
    if len(setups) != 1:
        raise SolderMaskError(f"{pcb_path}: expected one PCB setup")
    setup = setups[0]
    diagnostics = []
    try:
        width = float(sexp.val(setup, "solder_mask_min_width", [0])[0])
    except (ValueError, TypeError, IndexError):
        width = float("nan")
    if not math.isfinite(width) or width < MIN_SOLDER_MASK_WEB_MM:
        diagnostics.append(f"solder mask web must be at least {MIN_SOLDER_MASK_WEB_MM:g} mm")
    if sexp.val(setup, "allow_soldermask_bridges_in_footprints", []) != ["no"]:
        diagnostics.append("solder mask bridges within footprints must be disabled")
    for node in [setup] + F(root, "via"):
        tenting = F(node, "tenting")
        if not tenting and node is not setup:
            continue
        for side in ("front", "back"):
            value = sexp.val(tenting[0], side, []) if tenting else []
            if value == ["no"] or (node is setup and value != ["yes"]):
                diagnostics.append(f"{node[0]}: {side} via tenting must be enabled")
    if diagnostics:
        raise SolderMaskError(f"{pcb_path}: " + "; ".join(diagnostics))


def validate_plot_origin(pcb_path, board=None):
    """Return the board plot origin, or raise unless it is absolute (0, 0).

    Gerbers, drill, and centroid all export in absolute board coordinates. A
    non-zero `aux_axis_origin` shifts only the origin-relative exports, so the
    frames would diverge and JLC would place every part off-board.
    """
    root = board if board is not None else read_board(
        pcb_path, PlotOriginError, "PCB setup")

    origin = (0.0, 0.0)
    for setup in F(root, "setup"):
        raw = sexp.val(setup, "aux_axis_origin", [])
        if not raw:
            continue
        try:
            origin = (float(raw[0]), float(raw[1]))
        except (IndexError, TypeError, ValueError):
            raise PlotOriginError(
                f"{pcb_path}: aux_axis_origin is invalid: "
                + " ".join(str(value) for value in raw)) from None
        if origin != (0.0, 0.0):
            raise PlotOriginError(
                f"{pcb_path}: drill/place origin is {origin[0]:g},"
                f"{origin[1]:g} mm; the fab exports are pinned to absolute "
                "board coordinates. Reset the drill/place origin in Pcbnew "
                "before regenerating the fab package.")
    return origin


def validate_via_geometry(pcb_path, min_vias=MIN_BOARD_VIAS, board=None):
    root = board if board is not None else read_board(
        pcb_path, ViaGeometryError, "PCB vias")

    diagnostics = []
    vias = F(root, "via")
    if len(vias) < min_vias:
        diagnostics.append(
            f"board lists {len(vias)} vias, fewer than the {min_vias} the "
            "routed board carries; re-measure with gen/board_metadata.py "
            "after promoting a re-route")
    valid_vias = []
    for index, via in enumerate(vias, 1):
        at = sexp.val(via, "at", [])
        size = sexp.val(via, "size", [])
        drill = sexp.val(via, "drill", [])
        location = ",".join(str(value) for value in at[:2]) or f"index {index}"
        try:
            diameter_mm = float(size[0])
            drill_mm = float(drill[0])
        except (IndexError, TypeError, ValueError):
            diagnostics.append(f"via at {location}: size or drill is invalid")
            continue
        try:
            x_mm = float(at[0])
            y_mm = float(at[1])
        except (IndexError, TypeError, ValueError):
            diagnostics.append(f"via at {location}: position is invalid")
            continue
        net = sexp.val(via, "net", [])
        name = net_name(via) if net and str(net[0]) != "0" else None
        valid_vias.append((x_mm, y_mm, diameter_mm, name or None))
        if diameter_mm < MIN_STANDARD_VIA_DIAMETER_MM:
            diagnostics.append(
                f"via at {location}: {diameter_mm:g} mm diameter is below "
                f"{MIN_STANDARD_VIA_DIAMETER_MM:g} mm")
        if drill_mm < MIN_STANDARD_VIA_DRILL_MM:
            diagnostics.append(
                f"via at {location}: {drill_mm:g} mm drill is below "
                f"{MIN_STANDARD_VIA_DRILL_MM:g} mm")

    # Copper spacing is a different-net rule: same-net stitching vias may sit
    # as close as the drills allow. An unnetted via belongs to no net, so it
    # pairs against everything.
    for index, (x1, y1, diameter1, net1) in enumerate(valid_vias):
        for x2, y2, diameter2, net2 in valid_vias[index + 1:]:
            if net1 is not None and net1 == net2:
                continue
            spacing_mm = (
                math.hypot(x2 - x1, y2 - y1)
                - (diameter1 + diameter2) / 2
            )
            if spacing_mm + 1e-9 < MIN_VIA_TO_VIA_COPPER_SPACING_MM:
                diagnostics.append(
                    f"vias at {x1:g},{y1:g} and {x2:g},{y2:g}: "
                    f"{spacing_mm:g} mm via-to-via copper spacing is below "
                    f"{MIN_VIA_TO_VIA_COPPER_SPACING_MM:g} mm")

    if diagnostics:
        raise ViaGeometryError(
            "via geometry validation failed:\n  " + "\n  ".join(diagnostics))
    return len(vias)


def validate_zone_geometry(pcb_path, min_pours=MIN_COPPER_POURS, board=None,
                           *, plane_layers=GROUND_PLANE_LAYERS):
    """Gate required plane coverage, copper-pour pad connections and fill features.

    KiCad DRC never flags these: thermal reliefs are same-net geometry. A
    sub-process gap fabricates as a solid pour, defeating the reliefs R-ASM-6
    requires; `connect_pads yes` is a solid connection, and a zone with no
    `(fill ...)` node states no relief features. A zone whose fill is not
    enabled pours nothing.
    """
    root = board if board is not None else read_board(
        pcb_path, ZoneGeometryError, "PCB zones")

    diagnostics = []
    zones = [zone for zone in F(root, "zone") if is_copper_pour(zone)]
    if len(zones) < min_pours:
        diagnostics.append(
            f"board lists {len(zones)} copper pours, fewer than the "
            f"{min_pours} reference planes the design pours")
    if plane_layers:
        try:
            stack = copper_layers(root) if zones else ()
            poured = {layer for zone in zones for layer in zone_layers(zone, stack)}
        except ValueError as exc:
            diagnostics.append(str(exc))
        else:
            for layer in plane_layers:
                if layer not in poured:
                    diagnostics.append(f"required plane {layer} has no copper pour")
    for index, zone in enumerate(zones, 1):
        name = sexp.val(zone, "name", [])
        label = name[0] if name else f"index {index}"
        # (connect_pads [yes|no|thru_hole_only] (clearance X)): an absent mode
        # token is KiCad's thermal-relief default. `yes` is a solid connection.
        connect = sexp.val(zone, "connect_pads", [])
        if connect and not isinstance(connect[0], list) and connect[0] == "yes":
            diagnostics.append(
                f"zone {label}: connect_pads yes solders every pad straight "
                f"into the pour - through-hole GND joints need thermal "
                f"reliefs (R-ASM-6)")
        features = [("min_thickness", zone)]
        fills = F(zone, "fill")
        if not fills:
            diagnostics.append(
                f"zone {label}: no (fill ...) node, so "
                f"{' and '.join(ZONE_FILL_FEATURES)} are unstated")
        for fill in fills:
            # (fill [yes|no] (thermal_gap X) ...): an absent enable token leaves
            # the zone unfilled, so the plane plots as empty copper.
            enable = fill[1] if len(fill) > 1 and not isinstance(fill[1], list) \
                else "no"
            if enable != "yes":
                diagnostics.append(
                    f"zone {label}: fill {enable} pours no copper")
            features += [(key, fill) for key in ZONE_FILL_FEATURES]
        for key, node in features:
            value = sexp.val(node, key, [])
            try:
                feature_mm = float(value[0])
            except (IndexError, TypeError, ValueError):
                diagnostics.append(f"zone {label}: {key} is missing or invalid")
                continue
            minimum_mm = ZONE_FEATURE_MINIMUMS[key]
            if feature_mm + 1e-9 < minimum_mm:
                diagnostics.append(
                    f"zone {label}: {key} {feature_mm:g} mm is below "
                    f"{minimum_mm:g} mm")

    resolved_net_name = net_name_resolver(root)

    poured_nets = {resolved_net_name(zone) for zone in zones}
    for footprint in F(root, "footprint"):
        for pad in F(footprint, "pad"):
            if len(pad) < 3 or pad[2] != "thru_hole" or resolved_net_name(pad) not in poured_nets:
                continue
            connection = sexp.val(pad, "zone_connect", [])
            if not connection:
                connection = sexp.val(footprint, "zone_connect", [])
            if connection and str(connection[0]) == "2":
                diagnostics.append(
                    f"{footprint_reference(footprint)} pad {pad[1]}: zone_connect 2 "
                    "solders the through-hole joint straight into the pour (R-ASM-6)")

    if diagnostics:
        raise ZoneGeometryError(
            "zone geometry validation failed:\n  " + "\n  ".join(diagnostics))
    return len(zones)
