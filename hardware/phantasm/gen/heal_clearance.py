"""Restore Quilter and standard-cost via constraints in phantasm project files.

Quilter rejects an upload whose project has min_clearance == 0 ("min clearance must
be greater than zero"). KiCad re-zeroes that field every time the project is opened
in the GUI, so this heal must run as the LAST step before any Quilter upload --
for either the placed board (phantasm.kicad_pro) or the unplaced board
(unplaced/phantasm_unplaced.kicad_pro). Idempotent; safe to run anytime.

Hash-manifested snapshot directories -- those carrying SHA256SUMS.txt -- are skipped;
rewriting one breaks its manifest.

    python gen/heal_clearance.py
"""
import glob
import json
import os
import sys

from constraints import DEFAULT_CLASS_MINIMUMS, RULE_MINIMUMS

OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def is_manifested(path):
    """True when path sits in a directory whose contents are hash-manifested."""
    return os.path.exists(os.path.join(os.path.dirname(path), "SHA256SUMS.txt"))


def project_files():
    candidates = glob.glob(os.path.join(OUT, "phantasm*.kicad_pro")) \
        + glob.glob(os.path.join(OUT, "unplaced", "phantasm*.kicad_pro")) \
        + glob.glob(os.path.join(OUT, "quilter_incremental", "phantasm*.kicad_pro"))
    return [p for p in candidates if not is_manifested(p)]


def heal_project(p):
    with open(p, encoding="utf-8") as project_file:
        d = json.load(project_file)
    rules = d.setdefault("board", {}).setdefault("design_settings", {}).setdefault("rules", {})
    changes = {}
    for field, minimum in RULE_MINIMUMS.items():
        current = rules.get(field, 0) or 0
        if current < minimum:
            rules[field] = minimum
            changes[field] = (current, minimum)

    classes = d.setdefault("net_settings", {}).setdefault("classes", [])
    default = next((item for item in classes if item.get("name") == "Default"), None)
    if default is not None:
        diameter = default.get("via_diameter", 0) or 0
        drill = default.get("via_drill", 0) or 0
        minimum_diameter = DEFAULT_CLASS_MINIMUMS["via_diameter"]
        minimum_drill = DEFAULT_CLASS_MINIMUMS["via_drill"]
        if diameter < minimum_diameter:
            default["via_diameter"] = minimum_diameter
            changes["Default.via_diameter"] = (diameter, minimum_diameter)
        if drill < minimum_drill:
            default["via_drill"] = minimum_drill
            changes["Default.via_drill"] = (drill, minimum_drill)

    if changes:
        with open(p, "w", encoding="utf-8") as project_file:
            json.dump(d, project_file, indent=2)
        summary = ", ".join(
            f"{field} {old} -> {new}"
            for field, (old, new) in changes.items()
        )
        print(f"healed {os.path.relpath(p, OUT)}: {summary}")
    else:
        print(f"ok     {os.path.relpath(p, OUT)}")
    return bool(changes)


def main():
    pros = project_files()
    if not pros:
        print(f"error: no uploadable project files found under {OUT}", file=sys.stderr)
        return 1

    healed = 0
    for p in pros:
        try:
            healed += heal_project(p)
        except (OSError, json.JSONDecodeError) as error:
            print(f"error: cannot process {os.path.relpath(p, OUT)}: {error}",
                  file=sys.stderr)
            return 1

    print(f"\n{healed} file(s) healed. Upload-ready for Quilter.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
