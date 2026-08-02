"""Shared PHANTASM fabrication constraints in millimeters."""

RULE_MINIMUMS = {
    "min_clearance": 0.1016,
    "min_through_hole_diameter": 0.2,
    "min_via_annular_width": 0.125,
    "min_via_diameter": 0.45,
}

DEFAULT_CLASS_MINIMUMS = {
    "via_diameter": RULE_MINIMUMS["min_via_diameter"],
    "via_drill": 0.2,
}
