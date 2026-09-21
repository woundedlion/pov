"""Shared PHANTASM fabrication constraints (millimeters) and assembly policy."""

RULE_MINIMUMS = {
    "min_clearance": 0.1016,
    "min_copper_edge_clearance": 0.3,
    "min_through_hole_diameter": 0.2,
    "min_track_width": 0.13,
    "min_via_annular_width": 0.125,
    "min_via_diameter": 0.45,
}

DEFAULT_CLASS_MINIMUMS = {
    "via_diameter": RULE_MINIMUMS["min_via_diameter"],
    "via_drill": 0.2,
}

# unplaced/phantasm_unplaced.kicad_pro is a captured artifact: no generator
# writes it, and `pcb.py --unplaced` produces only the board and fp-lib-table.
# Its constraints are wider than the routed board's; these are the values to
# restore.
UNPLACED_RULES = {
    "min_clearance": 0.2,
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
EXCLUDE_FP_SUBSTR = ("PinHeader", "JST_", "SolderJumper", "CP_Radial")
EXCLUDE_VAL_SUBSTR = ("Teensy",)
