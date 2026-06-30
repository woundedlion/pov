"""Restore design_settings.rules.min_clearance > 0 in every phantasm project file.

Quilter rejects an upload whose project has min_clearance == 0 ("min clearance must
be greater than zero"). KiCad re-zeroes that field every time the project is opened
in the GUI, so this heal must run as the LAST step before any Quilter upload --
for either the placed board (phantasm.kicad_pro) or the unplaced board
(unplaced/phantasm_unplaced.kicad_pro). Idempotent; safe to run anytime.

    python gen/heal_clearance.py
"""
import glob
import json
import os

OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MIN_CLEARANCE = 0.2

# the project's own uploadable pros: top-level + the Quilter-prep subdirs (unplaced/,
# divider_rework/, ...). Quilter-result candidate dirs (v1/, v2/, ...) are not touched.
pros = glob.glob(os.path.join(OUT, "phantasm*.kicad_pro")) \
    + glob.glob(os.path.join(OUT, "unplaced", "phantasm*.kicad_pro")) \
    + glob.glob(os.path.join(OUT, "divider_rework", "phantasm*.kicad_pro"))

healed = 0
for p in pros:
    d = json.load(open(p, encoding="utf-8"))
    rules = d.setdefault("board", {}).setdefault("design_settings", {}).setdefault("rules", {})
    cur = rules.get("min_clearance")
    if not cur or cur <= 0:
        rules["min_clearance"] = MIN_CLEARANCE
        json.dump(d, open(p, "w", encoding="utf-8"), indent=2)
        print(f"healed {os.path.relpath(p, OUT)}: min_clearance {cur} -> {MIN_CLEARANCE}")
        healed += 1
    else:
        print(f"ok     {os.path.relpath(p, OUT)}: min_clearance = {cur}")

print(f"\n{healed} file(s) healed. Upload-ready for Quilter.")
