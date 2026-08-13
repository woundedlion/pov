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

# unplaced/phantasm_unplaced.kicad_pro is a captured artifact: no generator
# writes it, and `pcb.py --unplaced` produces only the board and fp-lib-table.
# Its constraints are wider than the routed board's, so heal_clearance.py --
# which only raises a value that is BELOW a minimum -- would restore the routed
# floors into a recreated file and hand Quilter a different constraint set than
# the candidate boards were produced under. These are the values to restore.
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
