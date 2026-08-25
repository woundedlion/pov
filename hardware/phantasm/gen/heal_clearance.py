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
import argparse
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


def project_files(paths=()):
    if paths:
        return [os.path.abspath(path) for path in paths if not is_manifested(path)]
    candidates = glob.glob(os.path.join(OUT, "phantasm*.kicad_pro")) \
        + glob.glob(os.path.join(OUT, "unplaced", "phantasm*.kicad_pro")) \
        + glob.glob(os.path.join(OUT, "quilter_incremental", "phantasm*.kicad_pro"))
    return [p for p in candidates if not is_manifested(p)]


def minimums_for(p):
    """The (rule, Default net class) floors project p must be restored to."""
    if os.path.basename(os.path.dirname(p)) == "unplaced":
        return UNPLACED_RULES, UNPLACED_DEFAULT_CLASS
    return RULE_MINIMUMS, DEFAULT_CLASS_MINIMUMS


def rule_shortfalls(d, rule_minimums, class_minimums):
    """Fields of project document d sitting below their fabrication floor.

    Maps the field name -- Default net class fields prefixed "Default." -- to
    (current, minimum). A field KiCad has dropped reads as 0, the same as one
    it re-zeroed.
    """
    rules = d.get("board", {}).get("design_settings", {}).get("rules", {})
    shortfalls = {}
    for field, minimum in rule_minimums.items():
        current = rules.get(field, 0) or 0
        if current < minimum:
            shortfalls[field] = (current, minimum)

    classes = d.get("net_settings", {}).get("classes", [])
    default = next((item for item in classes if item.get("name") == "Default"), None)
    if default is None:
        raise ValueError("missing Default net class")
    for field, minimum in class_minimums.items():
        current = default.get(field, 0) or 0
        if current < minimum:
            shortfalls[f"Default.{field}"] = (current, minimum)
    return shortfalls


def heal_project(p, dry_run=False):
    rule_minimums, class_minimums = minimums_for(p)
    with open(p, encoding="utf-8") as project_file:
        d = json.load(project_file)
        # Rewrite in the file's own convention; a mixed file gets the repo's.
        seen = project_file.newlines
        newline = seen if isinstance(seen, str) else "\n"
    changes = rule_shortfalls(d, rule_minimums, class_minimums)

    if changes:
        rules = d.setdefault("board", {}).setdefault(
            "design_settings", {}).setdefault("rules", {})
        default = next(item for item in d["net_settings"]["classes"]
                       if item.get("name") == "Default")
        for field, (_, minimum) in changes.items():
            if field.startswith("Default."):
                default[field[len("Default."):]] = minimum
            else:
                rules[field] = minimum
        if not dry_run:
            with open(p, "w", encoding="utf-8", newline=newline) as project_file:
                json.dump(d, project_file, indent=2)
                project_file.write("\n")
        summary = ", ".join(
            f"{field} {old} -> {new}"
            for field, (old, new) in changes.items()
        )
        action = "would heal" if dry_run else "healed"
        print(f"{action} {os.path.relpath(p, OUT)}: {summary}")
    else:
        print(f"ok     {os.path.relpath(p, OUT)}")
    return bool(changes)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("projects", nargs="*", help="project files to heal")
    parser.add_argument("--dry-run", action="store_true",
                        help="report changes without rewriting files")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    pros = project_files(args.projects)
    protected = [os.path.abspath(path) for path in args.projects
                 if is_manifested(path)]
    for path in protected:
        print(f"skip   {os.path.relpath(path, OUT)}: hash-manifested by "
              f"SHA256SUMS.txt; healing it would break the manifest",
              file=sys.stderr)
    if not pros:
        if protected:
            print("error: every named project file is hash-manifested",
                  file=sys.stderr)
        else:
            print(f"error: no uploadable project files found under {OUT}",
                  file=sys.stderr)
        return 1

    healed = 0
    for p in pros:
        try:
            healed += heal_project(p, args.dry_run)
        except (OSError, ValueError) as error:
            print(f"error: cannot process {os.path.relpath(p, OUT)}: {error}",
                  file=sys.stderr)
            return 1

    if args.dry_run:
        print(f"\n{healed} file(s) would be healed.")
    else:
        print(f"\n{healed} file(s) healed. Upload-ready for Quilter.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
