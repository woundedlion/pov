"""Upgrade phantasm_unplaced.kicad_pcb to a 4-layer SIG/GND/GND/SIG stackup with
inner ground planes, so an autoplacer/fab reads the intended stackup from the file.

Run with KiCad's python:
    "C:\\Program Files\\KiCad\\10.0\\bin\\python.exe" stackup.py
Re-run after `python pcb.py --unplaced` regenerates the board.
"""
import json
import os
import pcbnew

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "unplaced")
F = os.path.join(OUT, "phantasm_unplaced.kicad_pcb")
PRO = os.path.join(OUT, "phantasm_unplaced.kicad_pro")
MIN_CLEARANCE = 0.2   # Quilter rejects the KiCad default of 0 ("min clearance must be > 0")

b = pcbnew.LoadBoard(F)

# --- 4 copper layers: F.Cu / In1.Cu / In2.Cu / B.Cu ---
b.SetCopperLayerCount(4)
en = b.GetEnabledLayers()
for L in (pcbnew.In1_Cu, pcbnew.In2_Cu):
    en.AddLayer(L)
b.SetEnabledLayers(en)
b.SetVisibleLayers(en)

# mark the inner copper as plane (power-type) layers
for L in (pcbnew.In1_Cu, pcbnew.In2_Cu):
    try:
        b.SetLayerType(L, pcbnew.LT_POWER)
    except Exception:
        pass

# Physical dielectric stackup uses KiCad's default 1.6 mm 4-layer (1 oz copper,
# thin outer prepreg, thick core) — applied when the board is opened/saved in KiCad.

# --- inner ground planes: GND zones on In1.Cu and In2.Cu over the board outline ---
gnd = b.GetNetcodeFromNetname("GND")
bb = b.GetBoardEdgesBoundingBox()
x0, y0, x1, y1 = bb.GetLeft(), bb.GetTop(), bb.GetRight(), bb.GetBottom()
corners = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
for layer in (pcbnew.In1_Cu, pcbnew.In2_Cu):
    z = pcbnew.ZONE(b)
    z.SetLayer(layer)
    z.SetNetCode(gnd)
    z.SetAssignedPriority(0)
    poly = z.Outline()
    poly.NewOutline()
    for x, y in corners:
        poly.Append(int(x), int(y))
    z.SetIsFilled(False)
    b.Add(z)

b.Save(F)

# --- inject an explicit physical stackup (standard 1.6 mm 4-layer, JLC ...7628):
#     1 oz copper, ~0.21 mm outer prepreg (tight F.Cu->In1 GND coupling), 1.065 mm core.
STACKUP = """\t\t(stackup
\t\t\t(layer "F.SilkS" (type "Top Silk Screen"))
\t\t\t(layer "F.Paste" (type "Top Solder Paste"))
\t\t\t(layer "F.Mask" (type "Top Solder Mask") (thickness 0.01))
\t\t\t(layer "F.Cu" (type "copper") (thickness 0.035))
\t\t\t(layer "dielectric 1" (type "prepreg") (thickness 0.2104) (material "FR4") (epsilon_r 4.5) (loss_tangent 0.02))
\t\t\t(layer "In1.Cu" (type "copper") (thickness 0.035))
\t\t\t(layer "dielectric 2" (type "core") (thickness 1.065) (material "FR4") (epsilon_r 4.5) (loss_tangent 0.02))
\t\t\t(layer "In2.Cu" (type "copper") (thickness 0.035))
\t\t\t(layer "dielectric 3" (type "prepreg") (thickness 0.2104) (material "FR4") (epsilon_r 4.5) (loss_tangent 0.02))
\t\t\t(layer "B.Cu" (type "copper") (thickness 0.035))
\t\t\t(layer "B.Mask" (type "Bottom Solder Mask") (thickness 0.01))
\t\t\t(layer "B.Paste" (type "Bottom Solder Paste"))
\t\t\t(layer "B.SilkS" (type "Bottom Silk Screen"))
\t\t\t(copper_finish "ENIG")
\t\t\t(dielectric_constraints no)
\t\t)
"""
txt = open(F, encoding="utf-8").read()
if "(stackup" not in txt:
    txt = txt.replace("\t(setup\n", "\t(setup\n" + STACKUP, 1)
    open(F, "w", encoding="utf-8").write(txt)

# --- ensure the project DRC keeps min_clearance > 0 (KiCad resets it to 0 on every
#     open/save; Quilter rejects the upload otherwise). Self-heal the Quilter-prep. ---
if os.path.exists(PRO):
    pro = json.load(open(PRO, encoding="utf-8"))
    rules = pro.setdefault("board", {}).setdefault("design_settings", {}).setdefault("rules", {})
    if not rules.get("min_clearance"):
        rules["min_clearance"] = MIN_CLEARANCE
        json.dump(pro, open(PRO, "w", encoding="utf-8"), indent=2)
        print(f"patched {os.path.basename(PRO)}: min_clearance -> {MIN_CLEARANCE} mm (Quilter)")

print("4-layer SIG/GND/GND/SIG; copper layers =", b.GetCopperLayerCount(),
      "; inner GND planes + explicit 1.6 mm stackup encoded")
