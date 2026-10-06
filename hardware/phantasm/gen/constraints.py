"""Shared PHANTASM fabrication constraints (millimeters) and assembly policy."""

ZONE_DEFAULTS = {
    "min_thickness": 0.25,
    "thermal_gap": 0.5,
    "thermal_bridge_width": 0.5,
}

RULE_MINIMUMS = {
    "min_clearance": 0.1016,
    "min_copper_edge_clearance": 0.3,
    "min_hole_clearance": 0.1016,
    "min_hole_to_hole": 0.5,
    "min_resolved_spokes": 2,
    "min_through_hole_diameter": 0.2,
    "min_track_width": 0.13,
    "min_via_annular_width": 0.125,
    "min_via_diameter": 0.45,
}

NEW_LAYOUT_RULES = {
    "min_silk_clearance": 0.15,
    "solder_mask_to_copper_clearance": 0.1,
}

MIN_SOLDER_MASK_WEB_MM = 0.1
MIN_THERMAL_GAP_MM = 0.3
MIN_THERMAL_SPOKE_MM = 0.4
MAX_BOARD_WIDTH_MM = 35.0
MIN_VIA_TO_VIA_COPPER_SPACING_MM = 0.15

DEFAULT_CLASS_MINIMUMS = {
    "via_diameter": RULE_MINIMUMS["min_via_diameter"],
    "via_drill": 0.2,
}

# Quilter placement starts with wider constraints than the accepted routed board.
UNPLACED_RULES = {
    "min_clearance": 0.2,
    "min_hole_clearance": 0.25,
    "min_hole_to_hole": 0.5,
    "min_resolved_spokes": 2,
    "min_through_hole_diameter": 0.3,
    "min_via_annular_width": 0.125,
    "min_via_diameter": 0.6,
    "min_track_width": 0.2,
    "min_copper_edge_clearance": 0.5,
}

UNPLACED_DEFAULT_CLASS = {
    "clearance": 0.2,
    "track_width": 0.3,
    "via_diameter": 0.6,
    "via_drill": 0.3,
}

# Assembly policy: JLC reflows only top-side SMD. Exclude hand-soldered
# through-hole (connectors, electrolytic, Teensy) and solder jumpers.
EXCLUDE_FP_SUBSTR = ("TerminalBlock", "PinHeader", "JST_", "Molex_KK-254", "SolderJumper", "CP_Radial")
EXCLUDE_VAL_SUBSTR = ("Teensy",)


def rule_checks(d, rule_minimums, class_minimums):
    """Every fabrication floor and severity checked on project document d.

    Maps the field name -- Default net class fields prefixed "Default." -- to
    (current, required, met). A field KiCad has dropped reads as 0, the same
    as one it re-zeroed.
    """
    rules = d.get("board", {}).get("design_settings", {}).get("rules", {})
    checks = {}
    for field, minimum in {**rule_minimums, **NEW_LAYOUT_RULES}.items():
        current = rules.get(field, 0) or 0
        checks[field] = (current, minimum, current >= minimum)

    severity = d.get("board", {}).get("design_settings", {}).get(
        "rule_severities", {}).get("silk_over_copper")
    checks["rule_severities.silk_over_copper"] = (severity, "error", severity == "error")

    classes = d.get("net_settings", {}).get("classes", [])
    default = next((item for item in classes if item.get("name") == "Default"), None)
    if default is None:
        raise ValueError("missing Default net class")
    for field, minimum in class_minimums.items():
        current = default.get(field, 0) or 0
        checks[f"Default.{field}"] = (current, minimum, current >= minimum)
    return checks


def rule_shortfalls(d, rule_minimums, class_minimums):
    """The `rule_checks` entries of d that are not met, as field -> (current, required)."""
    return {field: (current, required)
            for field, (current, required, met)
            in rule_checks(d, rule_minimums, class_minimums).items() if not met}


def apply_project_floors(document, rule_minimums, class_minimums):
    """Apply shared fabrication floors to a project, creating its default class."""
    settings = document.setdefault("board", {}).setdefault("design_settings", {})
    rules = settings.setdefault("rules", {})
    classes = document.setdefault("net_settings", {}).setdefault("classes", [])
    default = next((item for item in classes if item.get("name") == "Default"), None)
    if default is None:
        default = {"name": "Default"}
        classes.append(default)
    for field, (_, minimum) in rule_shortfalls(document, rule_minimums, class_minimums).items():
        if field.startswith("Default."):
            default[field[len("Default."):]] = minimum
        elif field.startswith("rule_severities."):
            settings.setdefault("rule_severities", {})[field[len("rule_severities."):]] = minimum
        else:
            rules[field] = minimum
