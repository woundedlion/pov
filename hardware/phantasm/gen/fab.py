"""Regenerate the PCB fabrication / assembly outputs from the COMMITTED board.

Produces, into ../gen/out/:
  * jlc/        — Gerbers (Protel ext), Excellon drill, and a JLCPCB upload zip
  * jlc/phantasm-BOM.csv / phantasm-CPL.csv — JLCPCB assembly BOM + centroid
  * phantasm-drc.json — gating error-severity DRC report
  * phantasm-parity.json — gating board/schematic parity report

It NEVER runs board.py / pcb.py: those rewrite phantasm.kicad_{sch,pcb} and
discard the routing + silk. This script only reads the committed board and
emits derived artifacts (all of gen/out/ is gitignored).

kicad-cli is found via $KICAD_CLI, else common install paths, else PATH.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile

import check as netlist_spec
import sexp
from constraints import DEFAULT_CLASS_MINIMUMS, RULE_MINIMUMS
from kicad_common import F, is_copper_pour, kicad_cli

GEN = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(GEN)                       # hardware/phantasm
PCB = os.path.join(PROJ, "phantasm.kicad_pcb")
SCH = os.path.join(PROJ, "phantasm.kicad_sch")
OUT = os.path.join(GEN, "out")
JLC = os.path.join(OUT, "jlc")
# Prefix of the per-run staging directory the exports are written into, under OUT.
STAGE_PREFIX = "phantasm-jlc-"

# Layers JLCPCB needs (this board names silk "F.SilkS" / "B.SilkS").
GERBER_LAYERS = ("F.Cu,In1.Cu,In2.Cu,B.Cu,F.SilkS,B.SilkS,"
                 "F.Mask,B.Mask,F.Paste,B.Paste,Edge.Cuts")
# Fab-layer and drill files that belong in the JLC upload zip.
ZIP_MEMBERS = {
    "phantasm-F_Cu.gtl", "phantasm-In1_Cu.g1", "phantasm-In2_Cu.g2",
    "phantasm-B_Cu.gbl", "phantasm-F_Silkscreen.gto",
    "phantasm-B_Silkscreen.gbo", "phantasm-F_Mask.gts",
    "phantasm-B_Mask.gbs", "phantasm-F_Paste.gtp", "phantasm-B_Paste.gbp",
    "phantasm-Edge_Cuts.gm1", "phantasm-PTH.drl", "phantasm-NPTH.drl",
    "phantasm-job.gbrjob",
}
#: Digest manifest written beside the upload zip, never inside it.
SUMS_FILE = "SHA256SUMS.txt"

# Everything else the run writes into jlc/: assembly data, the zip itself,
# and the zip's digest manifest.
ZIP_EXCLUDED = {"phantasm-BOM.csv", "phantasm-CPL.csv",
                "phantasm-jlc-gerbers.zip", SUMS_FILE}


#: Creation stamp written over KiCad's wall clock in every exported artifact,
#: matching the zip member date below.
FAB_TIMESTAMP = "1980-01-01T00:00:00+00:00"
FAB_DATE = "1980-01-01 00:00:00"

#: KiCad's four export-time stamp forms: the Gerber X2 file attribute and its
#: plain banner, the two Excellon header comments, and the job file's JSON
#: field. Each substitution keeps the tool-version text so a toolchain upgrade
#: still shows in the diff, and stops before the line ending: the Gerbers are
#: CRLF and the Excellon files LF.
TIMESTAMP_SUBSTITUTIONS = (
    (re.compile(r"^%TF\.CreationDate,[^\r\n]*\*%(?=\r?$)", re.M),
     f"%TF.CreationDate,{FAB_TIMESTAMP}*%"),
    (re.compile(r"^(G04 Created by KiCad \([^\r\n]*\) date )[^\r\n]*\*(?=\r?$)",
                re.M),
     r"\g<1>" + FAB_DATE + "*"),
    (re.compile(r"^(; DRILL file KiCad [^\r\n]* date )[^\r\n]*(?=\r?$)", re.M),
     r"\g<1>" + FAB_TIMESTAMP),
    (re.compile(r"^; #@! TF\.CreationDate,[^\r\n]*(?=\r?$)", re.M),
     f"; #@! TF.CreationDate,{FAB_TIMESTAMP}"),
    (re.compile(r'^([ \t]*"CreationDate": ")[^"\r\n]*(")', re.M),
     r"\g<1>" + FAB_TIMESTAMP + r"\g<2>"),
)


class TimestampNormalizationError(ValueError):
    """An exported fab artifact would ship KiCad's export stamp."""


def normalize_fab_timestamps(directory):
    """Stamp FAB_TIMESTAMP over KiCad's export time; return the names stamped.

    Every Gerber, drill file and job file records the wall clock, so without
    this two exports of an unchanged board differ in every artifact and a
    re-run tells the reader nothing. The stamps are comments and metadata
    attributes; the copper is untouched.

    Every artifact in the directory must carry a stamp one of
    TIMESTAMP_SUBSTITUTIONS matches. One that cannot be decoded, or that a
    KiCad banner respelling puts out of reach of every pattern, raises rather
    than being skipped: matching nothing leaves the export times in place and
    only shrinks the reported count.
    """
    stamped_names = []
    unstamped = []
    for name in sorted(os.listdir(directory)):
        path = os.path.join(directory, name)
        try:
            with open(path, encoding="utf-8", newline="") as fh:
                text = fh.read()
        except (OSError, UnicodeError) as exc:
            raise TimestampNormalizationError(
                f"cannot read exported artifact {path}: {exc}\n"
                "  Its export stamp would ship unnormalized and every re-run "
                "of an unchanged board would differ.") from exc
        stamped = text
        matches = 0
        for pattern, replacement in TIMESTAMP_SUBSTITUTIONS:
            stamped, hits = pattern.subn(replacement, stamped)
            matches += hits
        if not matches:
            unstamped.append(name)
            continue
        if stamped != text:
            with open(path, "w", encoding="utf-8", newline="") as fh:
                fh.write(stamped)
        stamped_names.append(name)
    if unstamped:
        raise TimestampNormalizationError(
            "exported fab artifacts carry no recognized creation stamp: "
            + ", ".join(unstamped)
            + "\n  They would ship KiCad's export time and every re-run of an "
            "unchanged board would differ. Update TIMESTAMP_SUBSTITUTIONS to "
            "the current kicad-cli banner spellings.")
    return stamped_names


def zip_member(name):
    """Zip entry with fixed metadata, so an unchanged board rezips byte-identically.

    Source mtimes, host permissions and the packing platform all vary between
    runs; none of them belong in the fab package. The members themselves are
    stamp-normalized by normalize_fab_timestamps().
    """
    info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
    info.compress_type = zipfile.ZIP_DEFLATED
    info.create_system = 3
    info.external_attr = 0o644 << 16
    return info


def sha256_file(path):
    """SHA-256 of a file, over the bytes the upload zip stores."""
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def package_manifest(directory, members, archive):
    """`sha256sum -c` lines for every zipped artifact and the archive itself.

    zip_member() and normalize_fab_timestamps() make the package
    byte-reproducible; the manifest is what a rebuild is checked against.
    """
    return "".join(
        f"{sha256_file(os.path.join(directory, name))}  {name}\n"
        for name in list(members) + [archive])


class UploadPackageError(ValueError):
    """An exported fab artifact would not reach the JLC upload zip."""


def zip_members(names):
    """Sorted upload-zip members; every ZIP_MEMBERS file present, no others.

    The names bind every GERBER_LAYERS export plus both plated and unplated drill
    files. An omitted layer or drill would otherwise zip clean and fabricate an
    incomplete board.
    """
    members = sorted(n for n in names if n in ZIP_MEMBERS)
    dropped = sorted(set(names) - set(members) - ZIP_EXCLUDED)
    if dropped:
        raise UploadPackageError(
            "exported fab artifacts the JLC upload zip would drop: "
            + ", ".join(dropped)
            + " - add the name to ZIP_MEMBERS or ZIP_EXCLUDED.")
    absent = sorted(ZIP_MEMBERS - set(members))
    if absent:
        raise UploadPackageError(
            "the export produced no fab artifact for: " + ", ".join(absent)
            + " - check the kicad-cli gerber and drill exports.")
    return members

# Assembly policy: JLC reflows only top-side SMD. Exclude hand-soldered
# through-hole (connectors, electrolytic, Teensy), solder jumpers, and DNP.
EXCLUDE_FP_SUBSTR = ("PinHeader", "JST_", "SolderJumper", "CP_Radial")
EXCLUDE_VAL_SUBSTR = ("Teensy",)
ASSEMBLY_SIDE = "top"
MIN_STANDARD_VIA_DIAMETER_MM = RULE_MINIMUMS["min_via_diameter"]
MIN_STANDARD_VIA_DRILL_MM = DEFAULT_CLASS_MINIMUMS["via_drill"]
MIN_VIA_TO_VIA_COPPER_SPACING_MM = 0.15
MIN_ZONE_GAP_MM = RULE_MINIMUMS["min_clearance"]
MIN_ZONE_WIDTH_MM = RULE_MINIMUMS["min_track_width"]
ZONE_FILL_FEATURES = ("thermal_gap", "thermal_bridge_width")
ZONE_FEATURE_MINIMUMS = {
    "min_thickness": MIN_ZONE_WIDTH_MM,
    "thermal_gap": MIN_ZONE_GAP_MM,
    "thermal_bridge_width": MIN_ZONE_WIDTH_MM,
}
# Floors on what a routed board holds: the committed board carries 99 vias
# (gen/board_metadata.py) and pours the In1/In2 reference planes.
MIN_BOARD_VIAS = 99
MIN_COPPER_POURS = 2

# Parity items KiCad reports on a board that IS in sync with the schematic:
# the mounting holes have no symbol, the ID/shield jumpers are excluded from
# the BOM on the board only.
# Anything else means the routed copper predates the current schematic.
KNOWN_PARITY_ITEMS = {
    ("extra_footprint", "H1"): "Extra footprint",
    ("extra_footprint", "H2"): "Extra footprint",
    ("extra_footprint", "H3"): "Extra footprint",
    ("extra_footprint", "H4"): "Extra footprint",
    ("footprint_symbol_mismatch", "JP_ID0"): "Exclude from bill of materials",
    ("footprint_symbol_mismatch", "JP_ID1"): "Exclude from bill of materials",
    ("footprint_symbol_mismatch", "JP_ID2"): "Exclude from bill of materials",
    ("footprint_symbol_mismatch", "JP_SHLD"): "Exclude from bill of materials",
}
KNOWN_PARITY_WARNING_COUNTS = {
    "lib_footprint_mismatch": 11,
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
    "D_BUS": 180,  # SOD-323: JLC cathode band aligned to pad 1
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


def run(args, check=True, **kw):
    print("  $", os.path.basename(args[0]), " ".join(args[1:]))
    try:
        return subprocess.run(args, check=check, capture_output=True, text=True,
                              **kw)
    except subprocess.CalledProcessError as exc:
        sys.stderr.write(exc.stderr or "")
        raise
    except FileNotFoundError as exc:
        raise SystemExit(
            f"kicad-cli not found: {args[0]}; install KiCad, put kicad-cli on "
            "PATH, or set KICAD_CLI to its full path") from exc


def run_export(stage, args):
    """Run an export step, exiting with a diagnostic like the other stages."""
    try:
        return run(args)
    except subprocess.CalledProcessError as exc:
        sys.exit(f"{stage} export failed: kicad-cli exited {exc.returncode}")


class DesignRuleError(ValueError):
    pass


def read_drc_report(report_path):
    """Return (violations, unconnected_items) from a kicad-cli JSON DRC report.

    Raises unless both sections are present as lists, so a report whose shape
    changed reads as an error rather than as zero violations.
    """
    try:
        with open(report_path, encoding="utf-8") as fh:
            report = json.load(fh)
    except (OSError, UnicodeError, ValueError) as exc:
        raise DesignRuleError(
            f"cannot read DRC report: {report_path}") from exc

    if not isinstance(report, dict) or any(
            not isinstance(report.get(key), list)
            for key in ("violations", "unconnected_items")):
        raise DesignRuleError(
            "DRC report has no violations/unconnected_items sections: "
            f"{report_path}")
    return report["violations"], report["unconnected_items"]


def require_clean_drc(report_path):
    """Return DRC counts or raise unless both report sections are empty."""
    violations, unconnected_items = read_drc_report(report_path)
    num_violations = len(violations)
    num_unconnected = len(unconnected_items)
    if num_violations or num_unconnected:
        raise DesignRuleError(
            f"DRC failed: {num_violations} error-severity violations, "
            f"{num_unconnected} unconnected items -> {report_path}")
    return num_violations, num_unconnected


def run_drc(report_path):
    """Generate and require a clean error-severity DRC report."""
    if os.path.exists(report_path):
        os.remove(report_path)
    result = run([kicad_cli(), "pcb", "drc", "--severity-error", "--format", "json",
                  "--exit-code-violations", "-o", report_path, PCB],
                 check=False)
    counts = require_clean_drc(report_path)
    if result.returncode:
        raise DesignRuleError(
            f"kicad-cli drc exited {result.returncode} on a report that lists "
            f"no error-severity items: {report_path}")
    return counts


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
    """Return the parity item count or raise on any unexpected difference.

    Every KNOWN_PARITY_ITEMS entry must appear exactly once with its expected
    description: a count plus a membership test lets a duplicated item stand in
    for a vanished one.
    """
    try:
        with open(report_path, encoding="utf-8") as fh:
            report = json.load(fh)
    except (OSError, UnicodeError, ValueError) as exc:
        raise SchematicParityError(
            f"cannot read parity report: {report_path}") from exc

    if (not isinstance(report, dict) or
            not isinstance(report.get("schematic_parity"), list)):
        raise SchematicParityError(
            f"parity report has no schematic_parity section: {report_path}")
    if not isinstance(report.get("violations"), list):
        raise SchematicParityError(
            f"parity report has no violations section: {report_path}")
    entries = report["schematic_parity"]
    violations = report["violations"]

    diagnostics = []
    warning_counts = {}
    for violation in violations:
        kind = str(violation.get("type", ""))
        warning_counts[kind] = warning_counts.get(kind, 0) + 1
    for kind in sorted(set(warning_counts) |
                       set(KNOWN_PARITY_WARNING_COUNTS)):
        actual = warning_counts.get(kind, 0)
        expected = KNOWN_PARITY_WARNING_COUNTS.get(kind, 0)
        if actual != expected:
            diagnostics.append(
                f"{kind}: reported {actual} times in {report_path}, "
                f"expected {expected}")

    seen = {}
    for entry in entries:
        kind = str(entry.get("type", ""))
        refs = parity_refs(entry)
        description = str(entry.get("description", ""))
        key = (kind, refs[0]) if len(refs) == 1 else None
        expected_description = KNOWN_PARITY_ITEMS.get(key)
        if expected_description is not None and expected_description in description:
            seen[key] = seen.get(key, 0) + 1
            continue
        detail = "; ".join(str(item.get("description", ""))
                           for item in entry.get("items", []))
        diagnostics.append(f"{kind}: {description}"
                           + (f" [{detail}]" if detail else ""))
    for kind, ref in sorted(KNOWN_PARITY_ITEMS):
        times = seen.get((kind, ref), 0)
        if times != 1:
            diagnostics.append(
                f"{kind}: {ref} reported {times} times in {report_path}, "
                "expected exactly once")

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
    run([kicad_cli(), "pcb", "drc", "--schematic-parity", "--format", "json",
         "--severity-error", "--severity-warning", "-o", report_path, PCB],
        check=False)
    return require_schematic_parity(report_path)


class NetlistSpecError(ValueError):
    pass


def validate_netlist_spec(net_path):
    """Hold an exported netlist to check.py's named-net table; return its size.

    Mismatches are printed by the table's own reporter before this raises.
    """
    with open(net_path, encoding="utf-8") as fh:
        root = sexp.parse(fh.read())[0]
    if not netlist_spec.check(netlist_spec.netlist_nets(root)):
        raise NetlistSpecError(
            "netlist does not match the electrical specification (see the "
            "FAIL lines above)")
    return len(netlist_spec.EXPECT)


def parse_components(net_path):
    """ref -> {value, footprint, dnp} from a kicadsexpr netlist."""
    with open(net_path, encoding="utf-8") as fh:
        root = sexp.parse(fh.read())[0]
    comps = {}
    for block in F(root, "components"):
        for comp in F(block, "comp"):
            ref = sexp.val(comp, "ref")
            if not ref:
                continue
            props = {}
            for p in F(comp, "property"):
                name = sexp.val(p, "name")
                if name:
                    value = sexp.val(p, "value")
                    props[name[0]] = value[0] if value else ""
            val = sexp.val(comp, "value")
            fp = sexp.val(comp, "footprint")
            comps[ref[0]] = {
                "value": val[0] if val else "",
                "footprint": fp[0] if fp else "",
                "dnp": "dnp" in props or "DNP" in props,
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


def validate_assembled_refs(assembled):
    actual = set(assembled)
    expected = set(LCSC_BY_REF)
    diagnostics = []
    if missing := sorted(expected - actual):
        diagnostics.append(
            "assigned parts missing from assembly: " + ", ".join(missing))
    if unexpected := sorted(actual - expected):
        diagnostics.append(
            "assembled parts without assignments: " + ", ".join(unexpected))
    if diagnostics:
        raise AssemblyMetadataError(diagnostics)


def validate_rotation_refs(assembled):
    unknown = sorted(set(ROT_CORRECTION) - set(assembled))
    if unknown:
        raise AssemblyMetadataError([
            "rotation corrections missing from assembly: " +
            ", ".join(unknown)
        ])


def cpl_rotation(ref, rotation):
    return (rotation + ROT_CORRECTION.get(ref, 0)) % 360


class PartCatalogError(ValueError):
    pass


class ViaGeometryError(ValueError):
    pass


class ZoneGeometryError(ValueError):
    pass


class PlotOriginError(ValueError):
    pass


class BoardReadError(ValueError):
    pass


def read_board(pcb_path, error=BoardReadError, what="board"):
    """Parse a KiCad board; raise `error` naming `what` when it cannot be read.

    The board is 400 KB of s-expressions and the geometry gates each need the
    whole tree, so main() reads it once and hands the result to all three.
    """
    try:
        with open(pcb_path, encoding="utf-8") as fh:
            return sexp.parse(fh.read())[0]
    except (OSError, UnicodeError, ValueError) as exc:
        raise error(f"cannot read {what}: {pcb_path}") from exc


def validate_plot_origin(pcb_path, board=None):
    """Return the board plot origin, or raise unless it is absolute (0, 0).

    Gerbers, drill, and centroid all export in absolute board coordinates. A
    non-zero `aux_axis_origin` shifts only the origin-relative exports, so the
    frames would diverge and JLC would place every part off-board.
    """
    root = board if board is not None else read_board(
        pcb_path, PlotOriginError, "PCB setup")

    origin = (0.0, 0.0)
    for setup in F(root, "setup"):
        raw = sexp.val(setup, "aux_axis_origin", [])
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


def validate_via_geometry(pcb_path, min_vias=MIN_BOARD_VIAS, board=None):
    root = board if board is not None else read_board(
        pcb_path, ViaGeometryError, "PCB vias")

    diagnostics = []
    vias = F(root, "via")
    if len(vias) < min_vias:
        diagnostics.append(
            f"board lists {len(vias)} vias, fewer than the {min_vias} the "
            "routed board carries; re-measure with gen/board_metadata.py "
            "after promoting a re-route")
    valid_vias = []
    for index, via in enumerate(vias, 1):
        at = sexp.val(via, "at", [])
        size = sexp.val(via, "size", [])
        drill = sexp.val(via, "drill", [])
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
            if spacing_mm + 1e-9 < MIN_VIA_TO_VIA_COPPER_SPACING_MM:
                diagnostics.append(
                    f"vias at {x1:g},{y1:g} and {x2:g},{y2:g}: "
                    f"{spacing_mm:g} mm via-to-via copper spacing is below "
                    f"{MIN_VIA_TO_VIA_COPPER_SPACING_MM:g} mm")

    if diagnostics:
        raise ViaGeometryError(
            "via geometry validation failed:\n  " + "\n  ".join(diagnostics))
    return len(vias)


def validate_zone_geometry(pcb_path, min_pours=MIN_COPPER_POURS, board=None):
    """Gate copper-pour pad connection and fill features.

    KiCad DRC never flags these: thermal reliefs are same-net geometry, so a
    sub-process gap exports clean Gerbers that the fab resolves as a solid
    pour, tying every through-hole GND pad to the full plane and starving the
    hand-soldered joints R-ASM-6 protects. `connect_pads yes` reaches the same
    end directly, and a zone carrying no `(fill ...)` node states no relief
    features at all. A zone whose fill is not enabled pours nothing.
    """
    root = board if board is not None else read_board(
        pcb_path, ZoneGeometryError, "PCB zones")

    diagnostics = []
    zones = [zone for zone in F(root, "zone") if is_copper_pour(zone)]
    if len(zones) < min_pours:
        diagnostics.append(
            f"board lists {len(zones)} copper pours, fewer than the "
            f"{min_pours} reference planes the design pours")
    for index, zone in enumerate(zones, 1):
        name = sexp.val(zone, "name", [])
        label = name[0] if name else f"index {index}"
        # (connect_pads [yes|no|thru_hole_only] (clearance X)): an absent mode
        # token is KiCad's thermal-relief default. `yes` is a solid connection.
        connect = sexp.val(zone, "connect_pads", [])
        if connect and not isinstance(connect[0], list) and connect[0] == "yes":
            diagnostics.append(
                f"zone {label}: connect_pads yes solders every pad straight "
                f"into the pour - through-hole GND joints need thermal "
                f"reliefs (R-ASM-6)")
        features = [("min_thickness", zone)]
        fills = F(zone, "fill")
        if not fills:
            diagnostics.append(
                f"zone {label}: no (fill ...) node, so "
                f"{' and '.join(ZONE_FILL_FEATURES)} are unstated")
        for fill in fills:
            # (fill [yes|no] (thermal_gap X) ...): an absent enable token leaves
            # the zone unfilled, so the plane plots as empty copper.
            enable = fill[1] if len(fill) > 1 and not isinstance(fill[1], list) \
                else "no"
            if enable != "yes":
                diagnostics.append(
                    f"zone {label}: fill {enable} pours no copper")
            features += [(key, fill) for key in ZONE_FILL_FEATURES]
        for key, node in features:
            value = sexp.val(node, key, [])
            try:
                feature_mm = float(value[0])
            except (IndexError, TypeError, ValueError):
                diagnostics.append(f"zone {label}: {key} is missing or invalid")
                continue
            minimum_mm = ZONE_FEATURE_MINIMUMS[key]
            if feature_mm + 1e-9 < minimum_mm:
                diagnostics.append(
                    f"zone {label}: {key} {feature_mm:g} mm is below "
                    f"{minimum_mm:g} mm")

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


def validate_assembly_metadata(posrows, assembled):
    """Return per-ref assembly metadata, or raise with every diagnostic found.

    LCSC_BY_REF owns the supplier mapping; the schematic's LCSC field is never
    consulted, so an unassigned reference is a missing part number.
    """
    diagnostics = []
    metadata = {}
    for ref in assembled:
        lcsc = LCSC_BY_REF.get(ref, "")
        if not lcsc:
            diagnostics.append(f"{ref}: supplier part number (LCSC) is missing")

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


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    return parser.parse_args(argv)


def main():
    print(f"kicad-cli: {kicad_cli()}")
    if not os.path.exists(PCB):
        sys.exit(f"board not found: {PCB}")
    try:
        board = read_board(PCB)
    except BoardReadError as exc:
        sys.exit(str(exc))
    print("[1/9] Plot origin")
    try:
        validate_plot_origin(PCB, board=board)
    except PlotOriginError as exc:
        sys.exit(str(exc))
    print("  plot origin: absolute board coordinates for gerbers, drill, "
          "and centroid")
    print("[2/9] Via geometry")
    try:
        num_vias = validate_via_geometry(PCB, board=board)
    except ViaGeometryError as exc:
        sys.exit(str(exc))
    print(
        f"  via geometry: {num_vias} vias meet "
        f"{MIN_STANDARD_VIA_DIAMETER_MM:g}/"
        f"{MIN_STANDARD_VIA_DRILL_MM:g} mm diameter/drill minimum and "
        f"{MIN_VIA_TO_VIA_COPPER_SPACING_MM:g} mm minimum via-to-via "
        "copper spacing")
    print("[3/9] Zone geometry")
    try:
        num_zones = validate_zone_geometry(PCB, board=board)
    except ZoneGeometryError as exc:
        sys.exit(str(exc))
    print(
        f"  zone geometry: {num_zones} copper pours relieve their pads and "
        f"meet the {MIN_ZONE_GAP_MM:g} mm gap and "
        f"{MIN_ZONE_WIDTH_MM:g} mm width minimums")
    os.makedirs(OUT, exist_ok=True)
    print("[4/9] DRC report + schematic parity")
    rpt = os.path.join(OUT, "phantasm-drc.json")
    try:
        num_violations, num_unconnected = run_drc(rpt)
    except DesignRuleError as exc:
        sys.exit(str(exc))
    print(f"  DRC: {num_violations} error-severity violations, "
          f"{num_unconnected} unconnected items -> {rpt}")
    parity_rpt = os.path.join(OUT, "phantasm-parity.json")
    try:
        num_parity = run_parity(parity_rpt)
    except SchematicParityError as exc:
        sys.exit(str(exc))
    print(f"  parity: board matches schematic, {num_parity} known "
          f"differences -> {parity_rpt}")

    print("[5/9] Netlist + centroid")
    net = os.path.join(OUT, "_fab.net")
    run_export("netlist", [kicad_cli(), "sch", "export", "netlist", "--format",
                           "kicadsexpr", "-o", net, SCH])
    try:
        num_spec_nets = validate_netlist_spec(net)
    except NetlistSpecError as exc:
        sys.exit(str(exc))
    print(f"  netlist: {num_spec_nets} specified nets match member-for-member")
    pos = os.path.join(OUT, "_fab_pos.csv")
    # No --use-drill-file-origin: the CPL stays in the same absolute frame as
    # the gerbers and the drill file.
    run_export("centroid", [kicad_cli(), "pcb", "export", "pos", "--format", "csv",
                            "--units", "mm", "-o", pos, PCB])
    comps = parse_components(net)
    posrows = {}
    with open(pos, newline='', encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            posrows[r["Ref"]] = r
    assembled = sorted(r for r, c in comps.items() if is_assembled(c))
    # Gate before the jlc/ package is touched: a failure here must leave the
    # previous good package intact.
    try:
        validate_assembled_refs(assembled)
        validate_rotation_refs(assembled)
        assembly_metadata = validate_assembly_metadata(posrows, assembled)
    except AssemblyMetadataError as exc:
        sys.exit(str(exc))
    try:
        validate_part_catalog(assembly_metadata)
    except PartCatalogError as exc:
        sys.exit(str(exc))

    # TemporaryDirectory only removes its own: an aborted run strands its
    # staging directory under OUT, where nothing else clears it.
    for name in os.listdir(OUT):
        stranded = os.path.join(OUT, name)
        if name.startswith(STAGE_PREFIX) and os.path.isdir(stranded) \
                and not os.path.islink(stranded):
            shutil.rmtree(stranded)
    with tempfile.TemporaryDirectory(prefix=STAGE_PREFIX, dir=OUT) as staged:
        print("[6/9] Gerbers")
        run_export("gerber", [kicad_cli(), "pcb", "export", "gerbers",
                              "--check-zones", "--layers", GERBER_LAYERS,
                              "-o", staged + os.sep, PCB])
        print("[7/9] Drill")
        # Absolute origin, matching the gerbers and the centroid above.
        run_export("drill", [kicad_cli(), "pcb", "export", "drill", "--format",
                             "excellon", "--drill-origin", "absolute",
                             "--excellon-units", "mm", "--excellon-separate-th",
                             "-o", staged + os.sep, PCB])

        # Before the package is assembled: every staged artifact must carry a
        # creation stamp, and a failure has to leave jlc/ untouched.
        try:
            stamped = normalize_fab_timestamps(staged)
        except TimestampNormalizationError as exc:
            sys.exit(str(exc))

        if os.path.isdir(JLC):
            for name in os.listdir(JLC):
                stale = os.path.join(JLC, name)
                if os.path.isdir(stale) and not os.path.islink(stale):
                    shutil.rmtree(stale)
                else:
                    os.remove(stale)
        os.makedirs(JLC, exist_ok=True)
        for name in os.listdir(staged):
            os.replace(os.path.join(staged, name), os.path.join(JLC, name))
    print(f"  creation stamps: {len(stamped)} artifact(s) normalized to "
          f"{FAB_TIMESTAMP}")

    print("[8/9] BOM + CPL")
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
            rot = f"{cpl_rotation(r, p['rotation']):.6f}"
            w.writerow([r, p["pos_x"], p["pos_y"], p["side"], rot])

    print("[9/9] JLC upload zip")
    zpath = os.path.join(JLC, "phantasm-jlc-gerbers.zip")
    try:
        members = zip_members(os.listdir(JLC))
    except UploadPackageError as exc:
        sys.exit(str(exc))
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
        for f in members:
            with open(os.path.join(JLC, f), "rb") as fh:
                z.writestr(zip_member(f), fh.read())

    manifest_path = os.path.join(JLC, SUMS_FILE)
    with open(manifest_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(package_manifest(JLC, members, os.path.basename(zpath)))

    print(f"\nDone. {len(assembled)} assembled SMD parts; "
          f"{len(groups)} BOM lines.")
    print(f"  fab package: {zpath}")
    print(f"  package sha256: {sha256_file(zpath)}")
    print(f"  digest manifest: {manifest_path}")


if __name__ == "__main__":
    parse_args()
    main()
