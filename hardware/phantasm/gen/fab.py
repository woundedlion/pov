"""Regenerate the PCB fabrication / assembly outputs from the COMMITTED board.

Produces, into ../gen/out/:
  * jlc/        — Gerbers (Protel ext), Excellon drill, and a JLCPCB upload zip
  * jlc/phantasm-BOM.csv / phantasm-CPL.csv — JLCPCB assembly BOM + centroid
  * phantasm-drc.rpt — gating error-severity DRC report
  * phantasm-parity.json — gating board/schematic parity report

It NEVER runs board.py / pcb.py: those rewrite phantasm.kicad_{sch,pcb} and
discard the routing + silk. This script only reads the committed board and
emits derived artifacts (all of gen/out/ is gitignored).

kicad-cli is found via $KICAD_CLI, else common install paths, else PATH.
"""
import csv
import glob
import json
import math
import os
import re
import subprocess
import sys
import tempfile
import zipfile

import sexp
from kicad_common import F

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
EXCLUDE_FP_SUBSTR = ("PinHeader", "JST_", "SolderJumper", "CP_Radial")
EXCLUDE_VAL_SUBSTR = ("Teensy",)
ASSEMBLY_SIDE = "top"
MIN_STANDARD_VIA_DIAMETER_MM = 0.45
MIN_STANDARD_VIA_DRILL_MM = 0.20
MIN_VIA_COPPER_SPACING_MM = 0.15
# Smallest copper feature the fab resolves (4 mil), matching the board's
# min_clearance. Applies to pour fill features as well as tracks.
MIN_ZONE_FEATURE_MM = 0.1016
ZONE_FILL_FEATURES = ("thermal_gap", "thermal_bridge_width")

# Parity items KiCad reports on a board that IS in sync with the schematic:
# the mounting holes have no symbol, the ID/shield jumpers are excluded from
# the BOM on the board only, and the Teensy symbol carries no footprint field.
# Anything else means the routed copper predates the current schematic.
KNOWN_PARITY_ITEMS = {
    ("extra_footprint", "H1"),
    ("extra_footprint", "H2"),
    ("extra_footprint", "H3"),
    ("extra_footprint", "H4"),
    ("footprint_symbol_mismatch", "JP_ID0"),
    ("footprint_symbol_mismatch", "JP_ID1"),
    ("footprint_symbol_mismatch", "JP_ID2"),
    ("footprint_symbol_mismatch", "JP_SHLD"),
    ("footprint_symbol_mismatch", "U_MCU"),
}

# JLCPCB part assignments (LCSC #) keyed by reference. Kept here rather than in
# the schematic so the JLC assembly output owns the supplier mapping. R_D1/R_D2
# use an 0805 33R (C17634) to match their 0805 land pattern, not a 0603 part.
LCSC_BY_REF = {
    "C_DEC1": "C14663", "C_DEC2": "C14663",
    "C_LF": "C12891", "C_SYNC": "C1603",
    "D_BUS": "C1975255", "F1": "C261952", "FB": "C73732",
    "Q_REV": "C15127",
    "R1": "C25804", "R_MEN": "C25804", "R_PD": "C25804",
    "R2": "C22809", "R_D1": "C17634", "R_D2": "C17634",
    "R_LF": "C48928179", "R_S": "C17408", "U1": "C155176",
}
# CPL rotation correction (degrees added to KiCad's angle) for parts whose
# JLCPCB library zero-reference differs from KiCad's. U1 (SOIC-14): the raw
# KiCad angle of 180 reads wrong in JLC's viewer; a +270 correction (180 -> 90)
# lands pin 1 on the silk mark. Verify each against the assembly preview.
ROT_CORRECTION = {
    "Q_REV": 180,  # SOT-23: JLC single-lead side aligned to pad 3
    "U1": 270,    # SOIC-14: KiCad 180 -> 90, pin 1 on silk mark
}
# Exact JLCPCB catalog identity keyed by LCSC #. The LCSC number drives
# automatic matching; manufacturer and MPN make every match independently
# auditable in the uploaded BOM.
PART_BY_LCSC = {
    "C14663": {
        "manufacturer": "YAGEO",
        "mpn": "CC0603KRX7R9BB104",
        "description": "100nF 50V X7R +/-10% 0603 MLCC",
    },
    "C12891": {
        "manufacturer": "Samsung Electro-Mechanics",
        "mpn": "CL31A226KAHNNNE",
        "description": "22uF 25V X5R +/-10% 1206 MLCC",
    },
    "C1603": {
        "manufacturer": "Samsung Electro-Mechanics",
        "mpn": "CL10B221KB8NNNC",
        "description": "220pF 50V X7R +/-10% 0603 MLCC",
    },
    "C1975255": {
        "manufacturer": "Bourns",
        "mpn": "CDSOD323-T05L",
        "description": "5V unidirectional 1pF 9.8V-clamp 350W SOD-323 TVS diode",
    },
    "C261952": {
        "manufacturer": "TLC",
        "mpn": "TLC-NSMD050",
        "description": "13.2V 500mA hold 1A trip 1206 resettable PPTC fuse",
    },
    "C73732": {
        "manufacturer": "FH (Guangdong Fenghua Advanced Tech)",
        "mpn": "CBW321609U601T",
        "description": "600ohm @ 100MHz +/-25% 2A 1206 ferrite bead",
    },
    "C15127": {
        "manufacturer": "Alpha & Omega Semiconductor",
        "mpn": "AO3401A",
        "description": "-30V P-channel MOSFET 60mohm @ -4.5V SOT-23",
    },
    "C25804": {
        "manufacturer": "UNI-ROYAL",
        "mpn": "0603WAF1002T5E",
        "description": "10kohm +/-1% 100mW 0603 thick-film resistor",
    },
    "C22809": {
        "manufacturer": "UNI-ROYAL",
        "mpn": "0603WAF1502T5E",
        "description": "15kohm +/-1% 100mW 0603 thick-film resistor",
    },
    "C17634": {
        "manufacturer": "UNI-ROYAL",
        "mpn": "0805W8F330JT5E",
        "description": "33ohm +/-1% 125mW 0805 thick-film resistor",
    },
    "C48928179": {
        "manufacturer": "HKR (Hong Kong Resistors)",
        "mpn": "RCA051R5JLF",
        "description": "1.5ohm +/-5% 125mW 0805 thick-film resistor",
    },
    "C17408": {
        "manufacturer": "UNI-ROYAL",
        "mpn": "0805W8F1000T5E",
        "description": "100ohm +/-1% 125mW 0805 thick-film resistor",
    },
    "C155176": {
        "manufacturer": "Texas Instruments",
        "mpn": "SN74AHCT125DR",
        "description": "4.5V to 5.5V quad 3-state buffer SOIC-14",
    },
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


def run(args, check=True, **kw):
    print("  $", os.path.basename(args[0]), " ".join(args[1:]))
    return subprocess.run(args, check=check, capture_output=True, text=True, **kw)


def require_clean_drc(report_path):
    """Return DRC counts or raise unless both exact summaries are zero."""
    try:
        with open(report_path, encoding="utf-8") as fh:
            src = fh.read()
    except (OSError, UnicodeError) as exc:
        raise RuntimeError(f"cannot read DRC report: {report_path}") from exc

    violations = re.findall(
        r"^\*\* Found (\d+) DRC violations \*\*$", src, re.MULTILINE)
    unconnected = re.findall(
        r"^\*\* Found (\d+) unconnected pads \*\*$", src, re.MULTILINE)
    if len(violations) != 1 or len(unconnected) != 1:
        raise RuntimeError(f"cannot parse DRC report summaries: {report_path}")

    num_violations = int(violations[0])
    num_unconnected = int(unconnected[0])
    if num_violations or num_unconnected:
        raise RuntimeError(
            f"DRC failed: {num_violations} error-severity violations, "
            f"{num_unconnected} unconnected items")
    return num_violations, num_unconnected


def run_drc(report_path):
    """Generate and require a clean error-severity DRC report."""
    if os.path.exists(report_path):
        os.remove(report_path)
    # --exit-code-violations exits nonzero on a dirty board; the report is the diagnostic
    run([KCLI, "pcb", "drc", "--severity-error", "--exit-code-violations",
         "-o", report_path, PCB], check=False)
    return require_clean_drc(report_path)


class SchematicParityError(ValueError):
    pass


def parity_refs(entry):
    """Footprint references named by a parity entry's items."""
    refs = []
    for item in entry.get("items", []):
        match = re.match(r"Footprint (\S+)$", str(item.get("description", "")))
        if match:
            refs.append(match.group(1))
    return refs


def require_schematic_parity(report_path):
    """Return the parity item count or raise on any unexpected difference."""
    try:
        with open(report_path, encoding="utf-8") as fh:
            report = json.load(fh)
    except (OSError, UnicodeError, ValueError) as exc:
        raise SchematicParityError(
            f"cannot read parity report: {report_path}") from exc

    entries = report.get("schematic_parity", [])
    diagnostics = []
    for entry in entries:
        kind = str(entry.get("type", ""))
        refs = parity_refs(entry)
        if len(refs) == 1 and (kind, refs[0]) in KNOWN_PARITY_ITEMS:
            continue
        detail = "; ".join(str(item.get("description", ""))
                           for item in entry.get("items", []))
        diagnostics.append(f"{kind}: {entry.get('description', '')}"
                           + (f" [{detail}]" if detail else ""))

    if diagnostics:
        raise SchematicParityError(
            f"{PCB} no longer matches {SCH}; regenerate the board with "
            "gen/pcb.py --unplaced, re-route it in Quilter and promote the "
            "result before shipping these gerbers:\n  " +
            "\n  ".join(diagnostics))
    return len(entries)


def run_parity(report_path):
    """Generate and require a board that matches the committed schematic."""
    if os.path.exists(report_path):
        os.remove(report_path)
    # KiCad reports parity items at warning severity, so they are invisible to
    # the error-severity DRC gate and need their own run.
    run([KCLI, "pcb", "drc", "--schematic-parity", "--format", "json",
         "--severity-error", "--severity-warning", "-o", report_path, PCB],
        check=False)
    return require_schematic_parity(report_path)


def parse_components(net_path):
    """ref -> {value, footprint, dnp, lcsc} from a kicadsexpr netlist."""
    root = sexp.parse(open(net_path, encoding="utf-8").read())[0]
    comps = {}
    for block in F(root, "components"):
        for comp in F(block, "comp"):
            ref = sexp._val(comp, "ref")
            if not ref:
                continue
            props = {}
            for p in F(comp, "property"):
                name = sexp._val(p, "name")
                if name:
                    value = sexp._val(p, "value")
                    props[name[0]] = value[0] if value else ""
            val = sexp._val(comp, "value")
            fp = sexp._val(comp, "footprint")
            comps[ref[0]] = {
                "value": val[0] if val else "",
                "footprint": fp[0] if fp else "",
                "dnp": "dnp" in props or "DNP" in props,
                "lcsc": props.get("LCSC", ""),
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


class AssemblyMetadataError(ValueError):
    def __init__(self, diagnostics):
        self.diagnostics = tuple(diagnostics)
        message = "assembly metadata validation failed:\n  " + "\n  ".join(
            self.diagnostics)
        super().__init__(message)


class PartCatalogError(ValueError):
    pass


class ViaGeometryError(ValueError):
    pass


class ZoneGeometryError(ValueError):
    pass


class PlotOriginError(ValueError):
    pass


def validate_plot_origin(pcb_path):
    """Return the board plot origin, or raise unless it is absolute (0, 0).

    Gerbers, drill, and centroid all export in absolute board coordinates. A
    non-zero `aux_axis_origin` shifts only the origin-relative exports, so the
    frames would diverge and JLC would place every part off-board.
    """
    try:
        with open(pcb_path, encoding="utf-8") as fh:
            root = sexp.parse(fh.read())[0]
    except (OSError, UnicodeError, ValueError) as exc:
        raise PlotOriginError(f"cannot read PCB setup: {pcb_path}") from exc

    origin = (0.0, 0.0)
    for setup in F(root, "setup"):
        raw = sexp._val(setup, "aux_axis_origin", [])
        if not raw:
            continue
        try:
            origin = (float(raw[0]), float(raw[1]))
        except (IndexError, TypeError, ValueError):
            raise PlotOriginError(
                f"{pcb_path}: aux_axis_origin is invalid: "
                + " ".join(str(value) for value in raw)) from None
        if origin != (0.0, 0.0):
            raise PlotOriginError(
                f"{pcb_path}: drill/place origin is {origin[0]:g},"
                f"{origin[1]:g} mm; the fab exports are pinned to absolute "
                "board coordinates. Reset the drill/place origin in Pcbnew "
                "before regenerating the fab package.")
    return origin


def validate_via_geometry(pcb_path):
    try:
        with open(pcb_path, encoding="utf-8") as fh:
            root = sexp.parse(fh.read())[0]
    except (OSError, UnicodeError, ValueError) as exc:
        raise ViaGeometryError(f"cannot read PCB vias: {pcb_path}") from exc

    diagnostics = []
    vias = F(root, "via")
    valid_vias = []
    for index, via in enumerate(vias, 1):
        at = sexp._val(via, "at", [])
        size = sexp._val(via, "size", [])
        drill = sexp._val(via, "drill", [])
        location = ",".join(str(value) for value in at[:2]) or f"index {index}"
        try:
            diameter_mm = float(size[0])
            drill_mm = float(drill[0])
        except (IndexError, TypeError, ValueError):
            diagnostics.append(f"via at {location}: size or drill is invalid")
            continue
        try:
            x_mm = float(at[0])
            y_mm = float(at[1])
        except (IndexError, TypeError, ValueError):
            diagnostics.append(f"via at {location}: position is invalid")
            continue
        valid_vias.append((x_mm, y_mm, diameter_mm))
        if diameter_mm < MIN_STANDARD_VIA_DIAMETER_MM:
            diagnostics.append(
                f"via at {location}: {diameter_mm:g} mm diameter is below "
                f"{MIN_STANDARD_VIA_DIAMETER_MM:g} mm")
        if drill_mm < MIN_STANDARD_VIA_DRILL_MM:
            diagnostics.append(
                f"via at {location}: {drill_mm:g} mm drill is below "
                f"{MIN_STANDARD_VIA_DRILL_MM:g} mm")

    for index, (x1, y1, diameter1) in enumerate(valid_vias):
        for x2, y2, diameter2 in valid_vias[index + 1:]:
            spacing_mm = (
                math.hypot(x2 - x1, y2 - y1)
                - (diameter1 + diameter2) / 2
            )
            if spacing_mm + 1e-9 < MIN_VIA_COPPER_SPACING_MM:
                diagnostics.append(
                    f"vias at {x1:g},{y1:g} and {x2:g},{y2:g}: "
                    f"{spacing_mm:g} mm copper spacing is below "
                    f"{MIN_VIA_COPPER_SPACING_MM:g} mm")

    if diagnostics:
        raise ViaGeometryError(
            "via geometry validation failed:\n  " + "\n  ".join(diagnostics))
    return len(vias)


def validate_zone_geometry(pcb_path):
    """Gate copper-pour fill features at the fab's minimum feature size.

    KiCad DRC never flags these: thermal reliefs are same-net geometry, so a
    sub-process gap exports clean Gerbers that the fab resolves as a solid
    pour, tying every through-hole GND pad to the full plane and starving the
    hand-soldered joints R-ASM-6 protects.
    """
    try:
        with open(pcb_path, encoding="utf-8") as fh:
            root = sexp.parse(fh.read())[0]
    except (OSError, UnicodeError, ValueError) as exc:
        raise ZoneGeometryError(f"cannot read PCB zones: {pcb_path}") from exc

    diagnostics = []
    # Only net-bearing zones pour copper; the rest are keepouts.
    zones = [zone for zone in F(root, "zone") if F(zone, "net")]
    for index, zone in enumerate(zones, 1):
        name = sexp._val(zone, "name", [])
        label = name[0] if name else f"index {index}"
        features = [("min_thickness", zone)]
        for fill in F(zone, "fill"):
            features += [(key, fill) for key in ZONE_FILL_FEATURES]
        for key, node in features:
            value = sexp._val(node, key, [])
            try:
                feature_mm = float(value[0])
            except (IndexError, TypeError, ValueError):
                diagnostics.append(f"zone {label}: {key} is missing or invalid")
                continue
            if feature_mm + 1e-9 < MIN_ZONE_FEATURE_MM:
                diagnostics.append(
                    f"zone {label}: {key} {feature_mm:g} mm is below "
                    f"{MIN_ZONE_FEATURE_MM:g} mm")

    if diagnostics:
        raise ZoneGeometryError(
            "zone geometry validation failed:\n  " + "\n  ".join(diagnostics))
    return len(zones)


def validate_part_catalog(assembly_metadata):
    diagnostics = []
    required = ("manufacturer", "mpn", "description")
    for ref, metadata in assembly_metadata.items():
        lcsc = metadata["lcsc"]
        part = PART_BY_LCSC.get(lcsc)
        if part is None:
            diagnostics.append(f"{ref}: LCSC {lcsc} has no catalog entry")
            continue
        for field in required:
            if not str(part.get(field, "")).strip():
                diagnostics.append(
                    f"{ref}: LCSC {lcsc} catalog {field} is blank")
    if diagnostics:
        raise PartCatalogError(
            "supplier catalog validation failed:\n  " +
            "\n  ".join(diagnostics))


def validate_assembly_metadata(comps, posrows, assembled):
    diagnostics = []
    metadata = {}
    for ref in assembled:
        assigned_lcsc = LCSC_BY_REF.get(ref)
        raw_lcsc = assigned_lcsc if assigned_lcsc is not None else comps[ref].get(
            "lcsc")
        if raw_lcsc is None:
            diagnostics.append(f"{ref}: supplier part number (LCSC) is missing")
            lcsc = ""
        else:
            lcsc = str(raw_lcsc).strip()
            if not lcsc:
                diagnostics.append(f"{ref}: supplier part number (LCSC) is blank")

        pos = posrows.get(ref)
        if pos is None:
            diagnostics.append(f"{ref}: centroid row is missing")
            continue

        numeric_values = {}
        for field, label in (("PosX", "X"), ("PosY", "Y"),
                             ("Rot", "rotation")):
            raw_value = pos.get(field)
            if raw_value is None:
                diagnostics.append(f"{ref}: centroid {label} is missing")
                continue
            value_text = str(raw_value).strip()
            if not value_text:
                diagnostics.append(f"{ref}: centroid {label} is blank")
                continue
            try:
                value = float(value_text)
            except ValueError:
                diagnostics.append(
                    f"{ref}: centroid {label} is not numeric: {raw_value!r}")
                continue
            if not math.isfinite(value):
                diagnostics.append(
                    f"{ref}: centroid {label} is not finite: {raw_value!r}")
                continue
            numeric_values[field] = (value_text, value)

        raw_side = pos.get("Side")
        if raw_side is None:
            diagnostics.append(f"{ref}: centroid side is missing")
            side = ""
        else:
            side = str(raw_side).strip().lower()
            if not side:
                diagnostics.append(f"{ref}: centroid side is blank")
            elif side != ASSEMBLY_SIDE:
                diagnostics.append(
                    f"{ref}: centroid side is {raw_side!r}; expected "
                    f"{ASSEMBLY_SIDE!r}")

        if (lcsc and len(numeric_values) == 3 and side == ASSEMBLY_SIDE):
            metadata[ref] = {
                "lcsc": lcsc,
                "pos_x": numeric_values["PosX"][0],
                "pos_y": numeric_values["PosY"][0],
                "rotation": numeric_values["Rot"][1],
                "side": side,
            }

    if diagnostics:
        raise AssemblyMetadataError(diagnostics)
    return metadata


def main():
    print(f"kicad-cli: {KCLI}")
    if not os.path.exists(PCB):
        sys.exit(f"board not found: {PCB}")
    try:
        validate_plot_origin(PCB)
    except PlotOriginError as exc:
        sys.exit(str(exc))
    print("plot origin: absolute board coordinates for gerbers, drill, "
          "and centroid")
    try:
        num_vias = validate_via_geometry(PCB)
    except ViaGeometryError as exc:
        sys.exit(str(exc))
    print(
        f"via geometry: {num_vias} vias meet "
        f"{MIN_STANDARD_VIA_DIAMETER_MM:g}/"
        f"{MIN_STANDARD_VIA_DRILL_MM:g} mm minimum and "
        f"{MIN_VIA_COPPER_SPACING_MM:g} mm copper spacing")
    try:
        num_zones = validate_zone_geometry(PCB)
    except ZoneGeometryError as exc:
        sys.exit(str(exc))
    print(
        f"zone geometry: {num_zones} copper pours meet the "
        f"{MIN_ZONE_FEATURE_MM:g} mm minimum fill feature")
    os.makedirs(OUT, exist_ok=True)
    print("[1/6] DRC report + schematic parity")
    rpt = os.path.join(OUT, "phantasm-drc.rpt")
    num_violations, num_unconnected = run_drc(rpt)
    print(f"  DRC: {num_violations} error-severity violations, "
          f"{num_unconnected} unconnected items -> {rpt}")
    parity_rpt = os.path.join(OUT, "phantasm-parity.json")
    try:
        num_parity = run_parity(parity_rpt)
    except SchematicParityError as exc:
        sys.exit(str(exc))
    print(f"  parity: board matches schematic, {num_parity} known "
          f"differences -> {parity_rpt}")

    with tempfile.TemporaryDirectory(prefix="phantasm-jlc-", dir=OUT) as staged:
        print("[2/6] Gerbers")
        run([KCLI, "pcb", "export", "gerbers", "--layers", GERBER_LAYERS,
             "-o", staged + os.sep, PCB])
        print("[3/6] Drill")
        # Absolute origin, matching the gerbers and the centroid below.
        run([KCLI, "pcb", "export", "drill", "--format", "excellon",
             "--drill-origin", "absolute", "--excellon-units", "mm",
             "--excellon-separate-th", "-o", staged + os.sep, PCB])

        if os.path.isdir(JLC):
            for name in os.listdir(JLC):
                os.remove(os.path.join(JLC, name))
        os.makedirs(JLC, exist_ok=True)
        for name in os.listdir(staged):
            os.replace(os.path.join(staged, name), os.path.join(JLC, name))

    print("[4/6] Netlist + centroid")
    net = os.path.join(OUT, "_fab.net")
    run([KCLI, "sch", "export", "netlist", "--format", "kicadsexpr",
         "-o", net, SCH])
    pos = os.path.join(OUT, "_fab_pos.csv")
    # No --use-drill-file-origin: the CPL stays in the same absolute frame as
    # the gerbers and the drill file.
    run([KCLI, "pcb", "export", "pos", "--format", "csv", "--units", "mm",
         "-o", pos, PCB])
    comps = parse_components(net)
    posrows = {}
    with open(pos, newline='', encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            posrows[r["Ref"]] = r
    assembled = sorted(r for r, c in comps.items() if is_assembled(c))
    try:
        assembly_metadata = validate_assembly_metadata(comps, posrows, assembled)
    except AssemblyMetadataError as exc:
        sys.exit(str(exc))
    try:
        validate_part_catalog(assembly_metadata)
    except PartCatalogError as exc:
        sys.exit(str(exc))

    print("[5/6] BOM + CPL")
    # BOM grouped by (value, footprint)
    groups = {}
    for r in assembled:
        lcsc = assembly_metadata[r]["lcsc"]
        key = (comps[r]["value"], comps[r]["footprint"].split(":")[-1], lcsc)
        groups.setdefault(key, []).append(r)
    with open(os.path.join(JLC, "phantasm-BOM.csv"), "w", newline='',
              encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["Comment", "Designator", "Footprint", "LCSC Part #",
                    "Manufacturer", "Manufacturer Part Number", "Description"])
        for (v, f, lcsc), refs in sorted(groups.items()):
            part = PART_BY_LCSC[lcsc]
            w.writerow([
                v,
                ",".join(sorted(refs)),
                f,
                lcsc,
                part["manufacturer"],
                part["mpn"],
                part["description"],
            ])
    with open(os.path.join(JLC, "phantasm-CPL.csv"), "w", newline='',
              encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["Designator", "Mid X", "Mid Y", "Layer", "Rotation"])
        for r in assembled:
            p = assembly_metadata[r]
            rot = f"{(p['rotation'] + ROT_CORRECTION.get(r, 0)) % 360:.6f}"
            w.writerow([r, p["pos_x"], p["pos_y"], p["side"], rot])

    print("[6/6] JLC upload zip")
    zpath = os.path.join(JLC, "phantasm-jlc-gerbers.zip")
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
        for f in sorted(os.listdir(JLC)):
            if os.path.splitext(f)[1].lower() in ZIP_EXT:
                z.write(os.path.join(JLC, f), f)

    print(f"\nDone. {len(assembled)} assembled SMD parts; "
          f"{len(groups)} BOM lines.")
    print(f"  fab package: {zpath}")


if __name__ == "__main__":
    main()
