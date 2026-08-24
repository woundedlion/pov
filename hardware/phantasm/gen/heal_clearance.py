"""Restore Quilter and standard-cost via constraints in phantasm project files.

Quilter rejects an upload whose project has min_clearance == 0 ("min clearance must
be greater than zero"). KiCad re-zeroes that field every time the project is opened
in the GUI, so this heal must run as the LAST step before any Quilter upload --
for either the placed board (phantasm.kicad_pro) or the unplaced board
(unplaced/phantasm_unplaced.kicad_pro). Idempotent; safe to run anytime.

The unplaced board is restored to the wider constraints its candidate boards were
produced under (constraints.UNPLACED_RULES / UNPLACED_DEFAULT_CLASS); every other
project gets the routed board's fabrication floors.

Hash-manifested snapshot directories -- those carrying SHA256SUMS.txt -- are skipped;
rewriting one breaks its manifest.

    python gen/heal_clearance.py
"""
import glob
import json
import os
import sys

from constraints import (DEFAULT_CLASS_MINIMUMS, RULE_MINIMUMS,
                         UNPLACED_DEFAULT_CLASS, UNPLACED_RULES)

OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def is_manifested(path):
    """True when path sits in a directory whose contents are hash-manifested."""
    return os.path.exists(os.path.join(os.path.dirname(path), "SHA256SUMS.txt"))


def project_files():
    candidates = glob.glob(os.path.join(OUT, "phantasm*.kicad_pro")) \
        + glob.glob(os.path.join(OUT, "unplaced", "phantasm*.kicad_pro")) \
        + glob.glob(os.path.join(OUT, "quilter_incremental", "phantasm*.kicad_pro"))
    return [p for p in candidates if not is_manifested(p)]


def minimums_for(p):
    """The (rule, Default net class) floors project p must be restored to."""
    if os.path.basename(os.path.dirname(p)) == "unplaced":
        return UNPLACED_RULES, UNPLACED_DEFAULT_CLASS
    return RULE_MINIMUMS, DEFAULT_CLASS_MINIMUMS


def heal_project(p):
    rule_minimums, class_minimums = minimums_for(p)
    with open(p, encoding="utf-8") as project_file:
        d = json.load(project_file)
        # Rewrite in the file's own convention; a mixed file gets the repo's.
        seen = project_file.newlines
        newline = seen if isinstance(seen, str) else "\n"
    rules = d.setdefault("board", {}).setdefault("design_settings", {}).setdefault("rules", {})
    changes = {}
    for field, minimum in rule_minimums.items():
        current = rules.get(field, 0) or 0
        if current < minimum:
            rules[field] = minimum
            changes[field] = (current, minimum)

    classes = d.setdefault("net_settings", {}).setdefault("classes", [])
    default = next((item for item in classes if item.get("name") == "Default"), None)
    if default is None:
        raise ValueError("missing Default net class")
    for field, minimum in class_minimums.items():
        current = default.get(field, 0) or 0
        if current < minimum:
            default[field] = minimum
            changes[f"Default.{field}"] = (current, minimum)

    if changes:
        with open(p, "w", encoding="utf-8", newline=newline) as project_file:
            json.dump(d, project_file, indent=2)
            project_file.write("\n")
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
        except (OSError, ValueError) as error:
            print(f"error: cannot process {os.path.relpath(p, OUT)}: {error}",
                  file=sys.stderr)
            return 1

    print(f"\n{healed} file(s) healed. Upload-ready for Quilter.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
