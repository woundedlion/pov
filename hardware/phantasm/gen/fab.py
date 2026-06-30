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
        hits = sorted(glob.glob(p))
        if hits:
            return hits[-1]            # newest version
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
         "--side", "front", "--use-drill-file-origin", "-o", pos, PCB])
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
        key = (comps[r]["value"], comps[r]["footprint"].split(":")[-1],
               comps[r]["lcsc"])
        groups.setdefault(key, []).append(r)
    with open(os.path.join(JLC, "phantasm-BOM.csv"), "w", newline='',
              encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["Comment", "Designator", "Footprint", "LCSC Part #"])
        for (v, f, lcsc), refs in sorted(groups.items()):
            w.writerow([v, ",".join(sorted(refs)), f, lcsc])
    with open(os.path.join(JLC, "phantasm-CPL.csv"), "w", newline='',
              encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["Designator", "Mid X", "Mid Y", "Layer", "Rotation"])
        for r in assembled:
            p = posrows.get(r, {})
            w.writerow([r, p.get("PosX", ""), p.get("PosY", ""),
                        p.get("Side", "top"), p.get("Rot", "")])

    print("[6/6] JLC upload zip")
    zpath = os.path.join(JLC, "phantasm-jlc-gerbers.zip")
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
        for f in sorted(os.listdir(JLC)):
            if os.path.splitext(f)[1].lower() in ZIP_EXT:
                z.write(os.path.join(JLC, f), f)

    print(f"\nDone. {len(assembled)} assembled SMD parts; "
          f"{len(groups)} BOM lines (fill LCSC where blank).")
    print(f"  fab package: {zpath}")
    missing = [r for r in assembled if not comps[r]["lcsc"]]
    if missing:
        print(f"  NOTE: {len(missing)} parts have no LCSC number yet "
              f"(add an 'LCSC' field in the schematic).")


if __name__ == "__main__":
    main()
