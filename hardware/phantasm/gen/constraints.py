"""Shared PHANTASM fabrication constraints (millimeters) and assembly policy."""

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
    "min_via_diameter": 0.5,
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
# gen/pcb.py stamps the matching board attributes; gen/fab.py keeps the same
# parts out of the assembly BOM and centroid.
EXCLUDE_FP_SUBSTR = ("PinHeader", "JST_", "Molex_KK-254", "SolderJumper", "CP_Radial")
EXCLUDE_VAL_SUBSTR = ("Teensy",)
