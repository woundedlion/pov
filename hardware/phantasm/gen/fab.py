"""Regenerate the PCB fabrication / assembly outputs from the COMMITTED board.

Produces, into ../gen/out/:
  * jlc/        — Gerbers (Protel ext), Excellon drill, and a JLCPCB upload zip
  * jlc/phantasm-BOM.csv / phantasm-CPL.csv — JLCPCB assembly BOM + centroid
  * phantasm-drc.rpt — DRC report (informational; prints the violation count)

It NEVER runs board.py / pcb.py: those rewrite phantasm.kicad_{sch,pcb} and
discard the routing + silk. This script only reads the committed board and
emits derived artifacts (all of gen/out/ is gitignored).

kicad-cli is found via $KICAD_CLI, else common install paths, else PATH.
"""
import csv
import glob
import os
import re
import subprocess
import sys
import zipfile

GEN = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(GEN)                       # hardware/phantasm
PCB = os.path.join(PROJ, "phantasm.kicad_pcb")
SCH = os.path.join(PROJ, "phantasm.kicad_sch")
OUT = os.path.join(GEN, "out")
JLC = os.path.join(OUT, "jlc")

# Layers JLCPCB needs (this board names silk "F.SilkS" / "B.SilkS").
GERBER_LAYERS = ("F.Cu,In1.Cu,In2.Cu,B.Cu,F.SilkS,B.SilkS,"
                 "F.Mask,B.Mask,F.Paste,B.Paste,Edge.Cuts")
# Fab-layer + drill extensions that belong in the JLC upload zip.
ZIP_EXT = {".gtl", ".g1", ".g2", ".gbl", ".gto", ".gbo", ".gts", ".gbs",
           ".gtp", ".gbp", ".gm1", ".drl", ".gbrjob"}

# Assembly policy: JLC reflows only top-side SMD. Exclude hand-soldered
# through-hole (connectors, electrolytic, Teensy), solder jumpers, and DNP.
EXCLUDE_FP_SUBSTR = ("PinHeader", "SolderJumper", "CP_Radial")
EXCLUDE_VAL_SUBSTR = ("Teensy",)

# JLCPCB part assignments (LCSC #) keyed by reference. Kept here rather than in
# the schematic so the JLC assembly output owns the supplier mapping. R_D1/R_D2
# use an 0805 33R (C17634) to match their 0805 land pattern, not a 0603 part.
LCSC_BY_REF = {
    "C_DEC1": "C14663", "C_DEC2": "C14663",
    "C_LF": "C12891", "C_SYNC": "C1603",
    "F1": "C261952", "FB": "C73732", "Q_REV": "C7420330",
    "R1": "C25804", "R_MEN": "C25804", "R_PD": "C25804",
    "R2": "C22809", "R_D1": "C17634", "R_D2": "C17634",
    "R_LF": "C48928179", "R_S": "C17408", "U1": "C155176",
}
# CPL rotation correction (degrees added to KiCad's angle) for parts whose
# JLCPCB library zero-reference differs from KiCad's. U1 (SOIC-14): the raw
# KiCad angle of 180 reads wrong in JLC's viewer; a +270 correction (180 -> 90)
# lands pin 1 on the silk mark. Verify each against the assembly preview.
ROT_CORRECTION = {
    "U1": 270,    # SOIC-14: KiCad 180 -> 90, pin 1 on silk mark
    "Q_REV": 180, # SOD-123: cathode (band/"-") must land on +5V_PROT, not RAW
}
# JLCPCB catalog description keyed by LCSC #, for the BOM Description column.
DESC_BY_LCSC = {
    "C14663": "100nF 50V X7R ±10% 0603 Multilayer Ceramic Capacitors MLCC ROHS",
    "C12891": "22uF 25V X5R ±10% 1206 Multilayer Ceramic Capacitors MLCC ROHS",
    "C1603": "220pF 50V X7R ±10% 0603 Multilayer Ceramic Capacitors MLCC ROHS",
    "C261952": "13.2V 500mA hold 1A trip 750mΩ 1206 Resettable Fuse PPTC ROHS",
    "C73732": "600Ω@100MHz ±25% 2A 200mΩ 1206 Ferrite Bead ROHS",
    "C7420330": "40V 1A 600mV@1A SOD-123 Schottky Diode ROHS",
    "C25804": "10kΩ ±1% 100mW 0603 Thick Film Resistor ROHS",
    "C22809": "15kΩ ±1% 100mW 0603 Thick Film Resistor ROHS",
    "C17634": "33Ω ±1% 125mW 0805 Thick Film Resistor ROHS",
    "C48928179": "1.5Ω ±5% 125mW 0805 Thick Film Resistor ROHS",
    "C17408": "100Ω ±1% 125mW 0805 Thick Film Resistor ROHS",
    "C155176": "74AHCT125 4.5~5.5V quad bus buffer 3-state SOIC-14 ROHS",
}


def _kicad_version_key(path):
    """(major, minor) of a KiCad install path; (0, 0) if unversioned."""
    m = re.search(r"KiCad[\\/](\d+)\.(\d+)", path)
    return (int(m.group(1)), int(m.group(2))) if m else (0, 0)


def find_kicad_cli():
    env = os.environ.get("KICAD_CLI")
    if env and os.path.exists(env):
        return env
    pats = [
        r"C:\Program Files\KiCad\*\bin\kicad-cli.exe",
        r"C:\Program Files (x86)\KiCad\*\bin\kicad-cli.exe",
        "/Applications/KiCad/KiCad.app/Contents/MacOS/kicad-cli",
        "/usr/bin/kicad-cli",
        "/usr/local/bin/kicad-cli",
    ]
    for p in pats:
        hits = glob.glob(p)
        if hits:
            return max(hits, key=_kicad_version_key)   # newest version
    return "kicad-cli"                 # assume on PATH


KCLI = find_kicad_cli()


def run(args, **kw):
    print("  $", os.path.basename(args[0]), " ".join(args[1:]))
    return subprocess.run(args, check=True, capture_output=True, text=True, **kw)


def parse_components(net_path):
    """ref -> {value, footprint, dnp, lcsc} from a kicadsexpr netlist."""
    src = open(net_path, encoding="utf-8").read()
    comps = {}
    for m in re.finditer(r'\(comp\s', src):
        i = m.start()
        d = 0
        j = i
        while j < len(src):
            c = src[j]
            if c == '(':
                d += 1
            elif c == ')':
                d -= 1
                if d == 0:
                    break
            j += 1
        blk = src[i:j + 1]
        ref = re.search(r'\(ref "([^"]+)"\)', blk)
        if not ref:
            continue
        ref = ref.group(1)
        val = re.search(r'\(value "([^"]*)"\)', blk)
        fp = re.search(r'\(footprint "([^"]*)"\)', blk)
        lcsc = re.search(r'\(property\s*\(name "LCSC"\)\s*\(value "([^"]*)"\)', blk)
        comps[ref] = {
            "value": val.group(1) if val else "",
            "footprint": fp.group(1) if fp else "",
            "dnp": '(name "dnp")' in blk or '(name "DNP")' in blk,
            "lcsc": lcsc.group(1) if lcsc else "",
        }
    return comps


def is_assembled(c):
    if c["dnp"]:
        return False
    if any(s in c["footprint"] for s in EXCLUDE_FP_SUBSTR):
        return False
    if any(s in c["value"] for s in EXCLUDE_VAL_SUBSTR):
        return False
    return True


def main():
    print(f"kicad-cli: {KCLI}")
    if not os.path.exists(PCB):
        sys.exit(f"board not found: {PCB}")
    # clean fab dir, keep the rest of out/
    if os.path.isdir(JLC):
        for f in os.listdir(JLC):
            os.remove(os.path.join(JLC, f))
    os.makedirs(JLC, exist_ok=True)

    print("[1/6] Gerbers")
    run([KCLI, "pcb", "export", "gerbers", "--layers", GERBER_LAYERS,
         "-o", JLC + os.sep, PCB])
    print("[2/6] Drill")
    run([KCLI, "pcb", "export", "drill", "--format", "excellon",
         "--drill-origin", "absolute", "--excellon-units", "mm",
         "--excellon-separate-th", "-o", JLC + os.sep, PCB])

    print("[3/6] DRC report")
    rpt = os.path.join(OUT, "phantasm-drc.rpt")
    drc = subprocess.run([KCLI, "pcb", "drc", "--severity-error",
                          "--severity-warning", "-o", rpt, PCB],
                         capture_output=True, text=True)
    nviol = re.search(r"Found (\d+) violation", drc.stdout or "")
    print(f"  DRC: {nviol.group(1) if nviol else '?'} violations -> {rpt}")

    print("[4/6] Netlist + centroid")
    net = os.path.join(OUT, "_fab.net")
    run([KCLI, "sch", "export", "netlist", "--format", "kicadsexpr",
         "-o", net, SCH])
    pos = os.path.join(OUT, "_fab_pos.csv")
    run([KCLI, "pcb", "export", "pos", "--format", "csv", "--units", "mm",
         "--use-drill-file-origin", "-o", pos, PCB])
    comps = parse_components(net)
    posrows = {}
    with open(pos, newline='', encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            posrows[r["Ref"]] = r
    assembled = sorted(r for r, c in comps.items() if is_assembled(c))

    print("[5/6] BOM + CPL")
    # BOM grouped by (value, footprint)
    groups = {}
    for r in assembled:
        lcsc = LCSC_BY_REF.get(r, comps[r]["lcsc"])
        key = (comps[r]["value"], comps[r]["footprint"].split(":")[-1], lcsc)
        groups.setdefault(key, []).append(r)
    missing = sorted(r for r in assembled
                     if not LCSC_BY_REF.get(r, comps[r]["lcsc"]))
    if missing:
        print(f"  WARNING: no LCSC # for {', '.join(missing)}")
    with open(os.path.join(JLC, "phantasm-BOM.csv"), "w", newline='',
              encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["Comment", "Designator", "Footprint", "LCSC Part #",
                    "Description"])
        for (v, f, lcsc), refs in sorted(groups.items()):
            w.writerow([v, ",".join(sorted(refs)), f, lcsc,
                        DESC_BY_LCSC.get(lcsc, "")])
    with open(os.path.join(JLC, "phantasm-CPL.csv"), "w", newline='',
              encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["Designator", "Mid X", "Mid Y", "Layer", "Rotation"])
        for r in assembled:
            p = posrows.get(r, {})
            rot = p.get("Rot", "")
            try:
                rot = f"{(float(rot) + ROT_CORRECTION.get(r, 0)) % 360:.6f}"
            except ValueError:
                pass  # non-numeric rotation: pass through untouched
            w.writerow([r, p.get("PosX", ""), p.get("PosY", ""),
                        p.get("Side", "top"), rot])

    print("[6/6] JLC upload zip")
    zpath = os.path.join(JLC, "phantasm-jlc-gerbers.zip")
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
        for f in sorted(os.listdir(JLC)):
            if os.path.splitext(f)[1].lower() in ZIP_EXT:
                z.write(os.path.join(JLC, f), f)

    print(f"\nDone. {len(assembled)} assembled SMD parts; "
          f"{len(groups)} BOM lines (fill LCSC where blank).")
    print(f"  fab package: {zpath}")


if __name__ == "__main__":
    main()
