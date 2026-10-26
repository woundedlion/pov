"""Restore Quilter and standard-cost via constraints in phantasm project files.

Quilter rejects an upload whose project has min_clearance == 0 ("min clearance must
be greater than zero"). KiCad re-zeroes that field every time the project is opened
in the GUI, so this heal must run as the LAST step before any Quilter upload --
for any revision's phantasm.kicad_pro. Idempotent; safe to run anytime.

The unplaced board is restored to the wider constraints its candidate boards were
produced under (constraints.UNPLACED_RULES / UNPLACED_DEFAULT_CLASS); every other
project gets the routed board's fabrication floors.

Directories carrying candidate.json or SHA256SUMS.txt, and release.json trees,
are hash-manifested and skipped.

    python hardware/phantasm/gen/heal_clearance.py
"""
import argparse
import glob
import json
import os
from pathlib import Path
from kicad_common import atomic_write_text, is_manifested
import sys

from constraints import (DEFAULT_CLASS_MINIMUMS, RULE_MINIMUMS,
                         UNPLACED_DEFAULT_CLASS, UNPLACED_RULES,
                         apply_project_floors, rule_shortfalls)

OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def project_files(paths=()):
    if paths:
        return [os.path.abspath(path) for path in paths if not is_manifested(path)]
    candidates = glob.glob(os.path.join(OUT, "[0-9]*", "phantasm*.kicad_pro"))
    return sorted(p for p in candidates if not is_manifested(p))


def minimums_for(project=None):
    """The (rule, Default net class) floors the parsed project document must be
    restored to; unplaced projects get the wider Quilter floors."""
    layout = (project or {}).get("text_variables", {}).get("PHANTASM_LAYOUT")
    if layout == "unplaced":
        return UNPLACED_RULES, UNPLACED_DEFAULT_CLASS
    return RULE_MINIMUMS, DEFAULT_CLASS_MINIMUMS



def heal_project(p, dry_run=False):
    with open(p, encoding="utf-8") as project_file:
        d = json.load(project_file)
        # Rewrite in the file's own convention; a mixed file gets the repo's.
        seen = project_file.newlines
        newline = seen if isinstance(seen, str) else "\n"
    rule_minimums, class_minimums = minimums_for(d)
    changes = rule_shortfalls(d, rule_minimums, class_minimums)

    if changes:
        apply_project_floors(d, rule_minimums, class_minimums)
        if not dry_run:
            atomic_write_text(p, json.dumps(d, indent=2) + "\n", newline=newline)
        summary = ", ".join(
            f"{field} {old} -> {new}"
            for field, (old, new) in changes.items()
        )
        action = "would heal" if dry_run else "healed"
        print(f"{action} {display_path(p)}: {summary}")
    else:
        print(f"ok     {display_path(p)}")
    return bool(changes)


def display_path(path):
    resolved = Path(path).resolve()
    try:
        return str(resolved.relative_to(Path(OUT).resolve()))
    except ValueError:
        return str(resolved)


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
        print(f"skip   {display_path(path)}: hash-manifested by "
              f"SHA256SUMS.txt, candidate.json, or release.json; "
              f"healing it would break the manifest",
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
            print(f"error: cannot process {display_path(p)}: {error}",
                  file=sys.stderr)
            return 1

    if args.dry_run:
        print(f"\n{healed} file(s) would be healed.")
    else:
        print(f"\n{healed} file(s) healed. Upload-ready for Quilter.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
