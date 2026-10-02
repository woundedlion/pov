"""Generate phantasm.kicad_sch — PHANTASM per-segment carrier board.

Run from this directory:  python board.py
Writes the project files into hardware/phantasm/<revision>/.
Requires KiCad's stock symbol libs (see sexp.KICAD_SHARE / env KICAD_SYMBOL_DIR).

Layout: left-to-right signal flow in labelled blocks. Power distribution uses
visible horizontal rail wires with vertical component drops + junctions; signal
buses between the Teensy, level shifter and connectors use net labels (ports).
"""
import argparse
import copy
import json
import os
import builder as B
import sexp
from constraints import (DEFAULT_CLASS_MINIMUMS, NEW_LAYOUT_RULES, RULE_MINIMUMS,
                         UNPLACED_DEFAULT_CLASS, UNPLACED_RULES, apply_project_floors)
from kicad_common import atomic_write_text
from kicad_common import require_writable, reset_uid_sequence

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   B.REVISION)
SCH = os.path.join(OUT, "phantasm.kicad_sch")

SCH_REASON = (
    "Regeneration replaces this revision's schematic, including any KiCad\n"
    "  edits. Regenerate its board afterward to keep their symbol links aligned.")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--revision", choices=("1.2", "1.3"), default=B.REVISION)
    parser.add_argument("--force", action="store_true",
                        help="overwrite the committed phantasm.kicad_sch")
    return parser.parse_args(argv)


def project_seed(root_uuid):
    return json.dumps({
        "board": {"design_settings": {
            "rules": {**RULE_MINIMUMS, **NEW_LAYOUT_RULES},
            "rule_severities": {"silk_over_copper": "error"},
        }},
        "boards": [],
        "cvpcb": {"equivalence_files": []},
        "libraries": {"pinned_footprint_libs": [], "pinned_symbol_libs": []},
        "meta": {"filename": "phantasm.kicad_pro", "version": 3},
        "net_settings": {"classes": [{
            "name": "Default", "clearance": 0.2, "track_width": 0.3,
            **DEFAULT_CLASS_MINIMUMS,
        }]},
        "pcbnew": {"page_layout_descr_file": ""},
        "schematic": {
            "annotate_start_num": 0,
            "drawing": {"default_line_thickness": 6.0, "label_size_ratio": 0.375},
        },
        "sheets": [[root_uuid, "Root"]],
        "text_variables": {},
    }, indent=2) + "\n"


def write_project(path, root_uuid="", unplaced=None):
    if os.path.exists(path):
        with open(path, encoding="utf-8") as file:
            project = json.load(file)
    else:
        project = json.loads(project_seed(root_uuid))
    project.setdefault("meta", {})["filename"] = os.path.basename(path)
    if root_uuid:
        sheets = [entry for entry in project.get("sheets", []) if entry[1] != "Root"]
        project["sheets"] = [[root_uuid, "Root"], *sheets]
    if unplaced is None:
        unplaced = project.get("text_variables", {}).get("PHANTASM_LAYOUT") == "unplaced"
    project.setdefault("text_variables", {})["PHANTASM_LAYOUT"] = (
        "unplaced" if unplaced else "placed")
    apply_project_floors(project, UNPLACED_RULES if unplaced else RULE_MINIMUMS,
                         UNPLACED_DEFAULT_CLASS if unplaced else DEFAULT_CLASS_MINIMUMS)
    atomic_write_text(path, json.dumps(project, indent=2) + "\n")


def main(force=False, revision=B.REVISION, output_dir=None):
    if revision not in ("1.2", "1.3"):
        raise ValueError(f"unsupported board revision: {revision!r}")
    out = output_dir or (OUT if revision == B.REVISION else os.path.join(
        os.path.dirname(OUT), revision))
    sch = SCH if output_dir is None and revision == B.REVISION else os.path.join(
        out, "phantasm.kicad_sch")
    # uid() keys on call site + occurrence, so the sequence starts empty or a
    # second call in one process would renumber every generated uuid.
    reset_uid_sequence()
    require_writable(sch, force, SCH_REASON)

    b = B.Builder("PHANTASM Segment Board  -  per-segment carrier (x4, strap-selected role)",
                  paper="A3", revision=revision)
    GND = "GND"; V3 = "+3V3"


    # ---------------------------------------------------------------- custom symbols
    def make_power(net):
        node = copy.deepcopy(sexp.get_symbol("power", "+5V"))
        B._rename_subsymbols(node, "+5V", net)
        for c in node:
            if isinstance(c, list) and c and c[0] == "property":
                if c[1] == "Value":
                    c[2] = net
                elif c[1] == "Description":
                    c[2] = f'Power symbol creates a global label with name "{net}"'
        return b.register_custom(node, f"phantasm:{net}")


    def make_teensy():
        """Teensy 4.0 symbol with every numbered pad in the carrier footprint."""
        # (display-name, pad number), grouped by side.
        LEFT = [("D11/MOSI", "11"), ("D13/SCK", "13"), ("D3", "3"),
                ("D4", "4"), ("D5", "5"),
                ("D1/TX1", "1"), ("D21", "21"), ("D22", "22"), ("D23", "23")]
        RIGHT = [("VIN", "VIN"), ("3V3", "3V3"), ("GND", "GND")]
        UNUSED = ["0", "2", "6", "7", "8", "9", "10", "12",
                  "14", "15", "16", "17", "18", "19", "20"]
        bodyx = 13.97
        pitch = 5.08
        name2num = {}
        pins = []

        def addpin(name, number, x, y, ang, etype):
            name2num[name] = number
            pins.append(
                f'\t\t(pin {etype} line (at {x} {y} {ang}) (length 5.08)\n'
                f'\t\t\t(name "{name}" (effects (font (size 1.27 1.27))))\n'
                f'\t\t\t(number "{number}" (effects (font (size 1.27 1.27)))))')

        topL = (len(LEFT) - 1) / 2 * pitch
        for i, (name, num) in enumerate(LEFT):
            addpin(name, num, -(bodyx + 5.08), topL - i * pitch, 0, "passive")
        topR = (len(RIGHT) - 1) / 2 * pitch
        for i, (name, num) in enumerate(RIGHT):
            et = {"VIN": "power_in", "GND": "power_in", "3V3": "power_out"}[name]
            addpin(name, num, (bodyx + 5.08), topR - i * pitch, 180, et)

        unused_y = [30.48 - i * 2.54 for i in range(7)]
        unused_y += [-17.78 - i * 2.54 for i in range(8)]
        for num, y in zip(UNUSED, unused_y):
            addpin(f"D{num}", num, bodyx + 5.08, round(y, 2), 180, "passive")

        half = 40.64
        text = (
            '(symbol "phantasm:Teensy4.0"\n'
            '\t(pin_names (offset 1.016))\n'
            '\t(exclude_from_sim no) (in_bom yes) (on_board yes)\n'
            f'\t(property "Reference" "U" (at {-bodyx} {half + 2.54} 0)\n'
            '\t\t(effects (font (size 1.27 1.27)) (justify left)))\n'
            f'\t(property "Value" "Teensy4.0" (at {-bodyx} {-half - 2.54} 0)\n'
            '\t\t(effects (font (size 1.27 1.27)) (justify left)))\n'
            '\t(property "Footprint" "phantasm:Teensy4.0" (at 0 0 0)\n'
            '\t\t(effects (font (size 1.27 1.27)) (hide yes)))\n'
            '\t(property "Datasheet" "" (at 0 0 0)\n'
            '\t\t(effects (font (size 1.27 1.27)) (hide yes)))\n'
            '\t(symbol "Teensy4.0_0_1"\n'
            f'\t\t(rectangle (start {-bodyx} {half}) (end {bodyx} {-half})\n'
            '\t\t\t(stroke (width 0.254) (type default)) (fill (type background))))\n'
            '\t(symbol "Teensy4.0_1_1"\n'
            + "\n".join(pins) + ")\n)\n")
        node = sexp.parse(text)[0]
        return b.register_custom(node, "phantasm:Teensy4.0"), name2num, UNUSED


    make_power("+5V_RAW"); make_power("+5V_LOGIC")
    TEENSY, TPN, UNUSED_TEENSY_PINS = make_teensy()

    for lib, name in [("power", "GND"), ("power", "+3V3"), ("power", "PWR_FLAG"),
                      ("Device", "R"), ("Device", "C"), ("Device", "C_Polarized"),
                      ("Device", "FerriteBead"), ("Device", "Fuse"), ("Device", "D_Zener"),
                      ("Transistor_FET", "Q_PMOS_GSD"),
                      ("Connector_Generic", "Conn_01x02"),
                      ("Connector_Generic", "Conn_01x03"),
                      ("Jumper", "SolderJumper_2_Open"), ("74xx", "74AHCT125")]:
        b.ensure_lib(lib, name)

    if revision == "1.3":
        b.ensure_lib("Connector_Generic", "Conn_01x04")
        node = b._resolve("Interface_UART", "THVD1450D")
        B._rename_subsymbols(node, "THVD1450D", "THVD1410DR")
        for prop in [item for item in node if isinstance(item, list) and item and item[0] == "property"]:
            if prop[1] == "Value":
                prop[2] = "THVD1410DR"
            elif prop[1] == "Description":
                prop[2] = "500-kbps 3.3-V to 5-V RS-485 transceiver, SOIC-8"
        b.register_custom(node, "phantasm:THVD1410DR")
        protection = '''(symbol "phantasm:CDSOT23-SM712"
            (pin_names (offset 0.5)) (in_bom yes) (on_board yes)
            (property "Reference" "D" (at 0 7.62 0)
                (effects (font (size 1.27 1.27))))
            (property "Value" "CDSOT23-SM712" (at 0 -7.62 0)
                (effects (font (size 1.27 1.27))))
            (property "Datasheet" "https://www.bourns.com/docs/product-datasheets/cdsot23-sm712.pdf"
                (at 0 0 0) (effects (font (size 1.27 1.27)) (hide yes)))
            (symbol "CDSOT23-SM712_0_1"
                (rectangle (start -5.08 5.08) (end 5.08 -5.08)
                    (stroke (width 0.254) (type default)) (fill (type background))))
            (symbol "CDSOT23-SM712_1_1"
                (pin passive line (at -7.62 2.54 0) (length 2.54)
                    (name "LINE1" (effects (font (size 1.27 1.27))))
                    (number "1" (effects (font (size 1.27 1.27)))))
                (pin passive line (at -7.62 -2.54 0) (length 2.54)
                    (name "LINE2" (effects (font (size 1.27 1.27))))
                    (number "2" (effects (font (size 1.27 1.27)))))
                (pin passive line (at 0 -7.62 90) (length 2.54)
                    (name "GND" (effects (font (size 1.27 1.27))))
                    (number "3" (effects (font (size 1.27 1.27)))))))'''
        b.register_custom(sexp.parse_one(protection), "phantasm:CDSOT23-SM712")


    # ---------------------------------------------------------------- helpers
    def place(lib, ref, val, x, y, rot=0, unit=1, fp="", dnp=False, in_bom=True):
        return b.place(B.Symbol(lib, ref, val, x, y, rot=rot, unit=unit,
                                footprint=fp, dnp=dnp, in_bom=in_bom))


    def hw(x1, x2, y):
        b.wire((x1, y), (x2, y))


    def vw(x, y1, y2):
        b.wire((x, y1), (x, y2))


    power_index = 0


    def pwr_sym(kind, x, y):
        nonlocal power_index
        power_index += 1
        ref = f"#PWR{power_index:02d}"
        lib = {"GND": "power:GND", "+3V3": "power:+3V3"}.get(kind, f"phantasm:{kind}")
        b.place(B.Symbol(lib, ref, kind, x, y))


    def to_power(s, num, kind, length=2.54):
        end = b.stub(s, num, length)
        pwr_sym(kind, end[0], end[1])


    def to_label(s, num, name, length=3.81):
        end = b.stub(s, num, length)
        b.label(end, name)


    def bypass(parent, pin, ref, net="+5V_LOGIC", dx=15.24):
        x, y = parent.pin(pin)
        cap = place("Device:C", ref, "0.1uF", x + dx, y + 3.81, fp=C06)
        b.wire(parent.pin(pin), cap.pin("1"))
        to_power(cap, "1", net)
        to_power(cap, "2", GND)


    def rail_drop(s, num, yrail, span):
        """Drop a pin vertically onto the rail at yrail drawn over x range span.

        A junction past a rail end connects nothing: the pin becomes its own net.
        """
        t = s.pin(num)
        x0, x1 = span
        if not x0 - 1e-6 <= t[0] <= x1 + 1e-6:
            raise ValueError(
                f"{s.ref}.{num} drops at x={t[0]:g}, outside the y={yrail:g} "
                f"rail span {x0:g}..{x1:g}")
        vw(t[0], t[1], yrail)
        b.junction((t[0], yrail))


    def to_gnd_down(s, num, yrail, span):
        """Drop a (bottom) pin straight down onto a GND rail at yrail."""
        rail_drop(s, num, yrail, span)


    def to_rail_up(s, num, yrail, span):
        """Drop a (top) pin straight up onto a rail at yrail."""
        rail_drop(s, num, yrail, span)


    def series_wire(src_sym, src_pin, r_sym, far_label, src_label=None):
        """Wire src pin to the NEARER terminal of r_sym; label the far terminal.
        If src_label is given, name the source-side stub too (else KiCad auto-names
        it Net-(R-PadN), which shows up cryptic in autoplacers like Quilter)."""
        tip = src_sym.pin(src_pin)
        t1, t2 = r_sym.pin("1"), r_sym.pin("2")
        d1 = abs(t1[0] - tip[0]) + abs(t1[1] - tip[1])
        d2 = abs(t2[0] - tip[0]) + abs(t2[1] - tip[1])
        near, far = ("1", "2") if d1 <= d2 else ("2", "1")
        b.wire(tip, r_sym.pin(near))
        to_label(r_sym, far, far_label)
        if src_label:
            b.label(tip, src_label)


    SMD08 = "Resistor_SMD:R_0805_2012Metric"
    SMD06 = "Resistor_SMD:R_0603_1608Metric"
    # Rev 1.1 ships stock ids with pads widened in place (../1.1/README.md).
    # Rev 1.2 carries the toe-extended lands required by spec 11.1.
    SMD08_HAND = "Resistor_SMD:R_0805_2012Metric_Pad1.20x1.40mm_HandSolder"
    SMD06_HAND = "Resistor_SMD:R_0603_1608Metric_Pad0.98x0.95mm_HandSolder"
    C06 = "Capacitor_SMD:C_0603_1608Metric"


    # ============================================================ BLOCK 1: POWER
    b.text((25, 25), "POWER ENTRY / PROTECTION / RAIL FILTER  (logic ~0.15 A; LED 4.3 A off-board)", 2.2)
    Y_LOG, Y_GND = 60.96, 96.52

    # --- light logic feed only; LED 4.3 A power is delivered off-board (spec 2.3) ---
    # Series chain on the rail line: J1 -> F1 -> Q_REV -> FB -> +5V_LOGIC.
    J1 = place("Connector_Generic:Conn_01x02", "J1", "TBC05-02-1-G-G +5V IN ~1A", 25.4, 60.96, in_bom=False,
               rot=180,
               fp="phantasm:TerminalBlock_GCT_TBC05-02-1-G-G")
    # Small logic-only fuse/PTC (R-PWR-8) — the 4.3 A strip current never flows here.
    F1 = place("Device:Fuse", "F1", "0.5A hold", 40.64, 60.96, rot=90,
               fp="Fuse:Fuse_1206_3216Metric")
    # AO3401A pinout: 1=G, 2=S, 3=D. Drain faces the input so its body diode
    # initially raises the source; the grounded gate then enhances the channel.
    QREV = place("Transistor_FET:Q_PMOS_GSD", "Q_REV", "AO3401A", 53.34, 63.5,
                 rot=90, fp="Package_TO_SOT_SMD:SOT-23")
    FB = place("Device:FerriteBead", "FB", "600R@100MHz", 68.58, 60.96, rot=90,
               fp="Inductor_SMD:L_1206_3216Metric")

    # J1.1(+5V) -> F1.1 ; J1.2 -> GND ; name the pre-fuse entry node +5V_IN
    b.wire(J1.pin("1"), F1.pin("1"))
    b.label(F1.pin("1"), "+5V_IN")
    to_power(J1, "2", GND)
    # F1.2 -> Q_REV drain ; label the protected-input node +5V_RAW
    b.wire(F1.pin("2"), QREV.pin("3"))
    b.label(F1.pin("2"), "+5V_RAW")
    # Q_REV source -> FB.1 ; name the reverse-protected node +5V_PROT
    b.wire(QREV.pin("2"), FB.pin("1"))
    b.label(QREV.pin("2"), "+5V_PROT")
    to_power(QREV, "1", GND)

    # --- +5V_LOGIC rail (post-bead) and its drops: C_IN, R_LF/C_LF damper ---
    # The GND rail is drawn below, after the drops that fix its span.
    LOG_L, LOG_R = 76.2, 119.38
    GND_L, GND_R = 88.9, 119.38
    LOG_SPAN, GND_SPAN = (LOG_L, LOG_R), (GND_L, GND_R)
    b.wire(FB.pin("2"), (LOG_L, Y_LOG))
    hw(LOG_L, LOG_R, Y_LOG)
    pwr_sym("+5V_LOGIC", LOG_R, Y_LOG - 5.08); vw(LOG_R, Y_LOG - 5.08, Y_LOG)
    # C_IN: the card's only electrolytic, on the logic feed (R-PWR-3/6, spec 10).
    CIN = place("Device:C_Polarized", "C_IN", "100uF", 88.9, 78.74, in_bom=False,
                fp="Capacitor_THT:CP_Radial_D8.0mm_P3.50mm")
    to_rail_up(CIN, "1", Y_LOG, LOG_SPAN)
    to_gnd_down(CIN, "2", Y_GND, GND_SPAN)
    # R_LF + C_LF bead-LC damper (R-PWR-5): LOG -> R_LF -> CLF_NODE -> C_LF -> GND
    RLF = place("Device:R", "R_LF", "1R5", 109.22, 71.12, fp=SMD08)
    CLF = place("Device:C", "C_LF", "22uF", 109.22, 86.36, fp="Capacitor_SMD:C_1206_3216Metric")
    to_rail_up(RLF, "1", Y_LOG, LOG_SPAN)
    b.wire(RLF.pin("2"), CLF.pin("1"))   # CLF_NODE
    b.label(CLF.pin("1"), "LF_DAMP")
    to_gnd_down(CLF, "2", Y_GND, GND_SPAN)

    # --- GND rail spanning the cap drops + its GND symbol at the right end ---
    hw(GND_L, GND_R, Y_GND)
    pwr_sym(GND, GND_R, Y_GND)

    # ============================================================ BLOCK 2: LOGIC
    b.text((25, 116), "TEENSY 4.0  +  74AHCT125 LEVEL SHIFTER  ->  LED STRIP", 2.2)
    # --- Teensy ---
    U = place(TEENSY, "U_MCU", "Teensy4.0", 60.96, 165.1,
              fp="phantasm:Teensy4.0", in_bom=False)
    tn = lambda d: TPN[d]
    bypass(U, tn("VIN"), "C_DEC1")
    to_power(U, tn("3V3"), V3)
    to_power(U, tn("GND"), GND)
    to_label(U, tn("D11/MOSI"), "DATA_IN")
    to_label(U, tn("D13/SCK"), "CLK_IN")
    to_label(U, tn("D3"), "FRAME_SYNC")
    to_label(U, tn("D4"), "SYNC_TX")
    to_label(U, tn("D5"), "MASTER_EN")
    to_label(U, tn("D21"), "ID0")
    to_label(U, tn("D22"), "ID1")
    to_label(U, tn("D23"), "ID2")   # read by the N=8 firmware profile
    b.no_connect(U.pin(tn("D1/TX1")))
    for num in UNUSED_TEENSY_PINS:
        b.no_connect(U.pin(num))

    # --- U1 buffer units (A/B/C/D) + power unit (E) ---
    ux = 152.4
    U1A = place("74xx:74AHCT125", "U1", "74AHCT125", ux, 137.16, unit=1,
                fp="Package_SO:SOIC-14_3.9x8.7mm_P1.27mm")
    U1B = place("74xx:74AHCT125", "U1", "74AHCT125", ux, 167.64, unit=2,
                fp="Package_SO:SOIC-14_3.9x8.7mm_P1.27mm")
    U1C = place("74xx:74AHCT125", "U1", "74AHCT125", ux, 198.12, unit=3,
                fp="Package_SO:SOIC-14_3.9x8.7mm_P1.27mm")
    U1D = place("74xx:74AHCT125", "U1", "74AHCT125", ux, 228.6, unit=4,
                fp="Package_SO:SOIC-14_3.9x8.7mm_P1.27mm")
    U1E = place("74xx:74AHCT125", "U1", "74AHCT125", 215.9, 137.16, unit=5,
                fp="Package_SO:SOIC-14_3.9x8.7mm_P1.27mm")
    RD1 = place("Device:R", "R_D1", "33R", 177.8, 137.16, rot=270, fp=SMD08)
    RD2 = place("Device:R", "R_D2", "33R", 177.8, 167.64, rot=270, fp=SMD08)
    if revision == "1.2":
        RS = place("Device:R", "R_S", "100R", 177.8, 198.12, rot=270, fp=SMD08_HAND)
    # ch A (DATA) — source stub (U1->R_D1) named DATA_SRC; post-term net is DATA
    to_label(U1A, "2", "DATA_IN"); to_power(U1A, "1", GND); series_wire(U1A, "3", RD1, "DATA", "DATA_SRC")
    # ch B (CLK)
    to_label(U1B, "5", "CLK_IN"); to_power(U1B, "4", GND); series_wire(U1B, "6", RD2, "CLK", "CLK_SRC")
    # ch C (SYNC)
    if revision == "1.2":
        to_label(U1C, "9", "SYNC_TX"); to_label(U1C, "10", "MASTER_EN"); series_wire(U1C, "8", RS, "SYNC_BUS", "SYNC_SRC")
    else:
        to_power(U1C, "9", GND); to_power(U1C, "10", "+5V_LOGIC")
        b.no_connect(U1C.pin("8"))
    RTX = place("Device:R", "R_TX", "10k", 111.76, 213.36, fp=SMD06)
    to_label(RTX, "1", "SYNC_TX"); to_power(RTX, "2", GND if revision == "1.2" else V3)
    # ch D switches the single bus idle pull-down on only when this board is master.
    # MASTER_EN is LOW on the master and HIGH on slaves, matching the active-low OE.
    to_power(U1D, "12", GND)
    if revision == "1.2":
        to_label(U1D, "13", "MASTER_EN")
        to_label(U1D, "11", "SYNC_PULLDOWN")
    else:
        to_power(U1D, "13", "+5V_LOGIC")
        b.no_connect(U1D.pin("11"))
    # power unit
    bypass(U1E, "14", "C_DEC2")
    to_power(U1E, "7", GND)

    # --- J2 strip SIGNAL out (3-pin, no power): DI / SIG_GND / CI (R-CON-1) ---
    # Strip 5 V/GND are injected off-board (spec 2.3); SIG_GND is the card's logic GND,
    # landed on the strip GND pin at the load end (the off-board ground star).
    J2 = place("Connector_Generic:Conn_01x03", "J2", "TBC05-03-1-G-G LED DI/SIG_GND/CI", 281.94, 165.1, in_bom=False,
               fp="phantasm:TerminalBlock_GCT_TBC05-03-1-G-G")
    to_label(J2, "1", "DATA"); to_power(J2, "2", GND); to_label(J2, "3", "CLK")

    # ============================================================ BLOCK 3: SYNC
    if revision == "1.2":
        b.text((220, 180), "SYNC BUS  -  RX DIVIDER + DAISY (Belden 8451)", 2.2)
        # divider: SYNC_BUS -> R1 -> node(FRAME_SYNC) -> R2 -> GND ; C_SYNC at node
        R1 = place("Device:R", "R1", "10k", 246.38, 205.74, fp=SMD06_HAND)
        R2 = place("Device:R", "R2", "15k", 246.38, 228.6, fp=SMD06_HAND)
        CSY = place("Device:C", "C_SYNC", "220pF", 269.24, 217.17, rot=90, fp=C06)
        to_label(R1, "1", "SYNC_BUS")
        nd = R1.pin("2")
        b.wire(nd, R2.pin("1"))             # divider node
        b.label(nd, "FRAME_SYNC")
        b.wire(R2.pin("1"), CSY.pin("1"))   # C_SYNC onto node
        b.junction(R2.pin("1"))
        to_power(R2, "2", GND); to_power(CSY, "2", GND)
        # Master-only bus idle pulldown + bus TVS. U1 ch D drives
        # SYNC_PULLDOWN low on the master and is high-impedance on every slave.
        RPD = place("Device:R", "R_PD", "10k", 292.1, 205.74, fp=SMD06_HAND)
        to_label(RPD, "1", "SYNC_BUS"); to_label(RPD, "2", "SYNC_PULLDOWN")
        DBUS = place("Device:D_Zener", "D_BUS", "CDSOD323-T08L", 292.1, 231.14,
                     fp="Diode_SMD:D_SOD-323")
        to_label(DBUS, "1", "SYNC_BUS"); to_power(DBUS, "2", GND)
        # daisy connectors
        J3A = place("Connector_Generic:Conn_01x03", "J3A", "TBC05-03-1-G-G SYNC in", 330.2, 205.74, in_bom=False,
                    fp="phantasm:TerminalBlock_GCT_TBC05-03-1-G-G")
        J3B = place("Connector_Generic:Conn_01x03", "J3B", "TBC05-03-1-G-G SYNC out", 330.2, 233.68, in_bom=False,
                    fp="phantasm:TerminalBlock_GCT_TBC05-03-1-G-G")
        for J in (J3A, J3B):
            to_label(J, "1", "SYNC_BUS"); to_power(J, "2", GND); to_label(J, "3", "SHIELD")
    else:
        b.text((220, 180), "RS-485 SYNC / 120R TWISTED-PAIR TRUNK", 2.2)
        urs = place("phantasm:THVD1410DR", "U_SYNC", "THVD1410DR", 246.38, 210.82,
                    fp="Package_SO:SOIC-8_3.9x4.9mm_P1.27mm")
        to_label(urs, "1", "FRAME_SYNC"); to_power(urs, "2", GND)
        to_label(urs, "3", "MASTER_EN"); to_label(urs, "4", "SYNC_TX")
        to_power(urs, "5", GND); to_power(urs, "8", V3)
        bypass(urs, "8", "C_DEC3", net=V3, dx=-15.24)
        cbulk = place("Device:C", "C_BULK3", "1uF", 204.47, 187.96, fp=C06)
        to_power(cbulk, "1", V3); to_power(cbulk, "2", GND)
        for ref, pin, net, y in (("R_A", "6", "SYNC_A", 200.66),
                                 ("R_B", "7", "SYNC_B", 217.17)):
            resistor = place("Device:R", ref, "10R CRCW0603010RJNEAHP", 281.94, y,
                             rot=270, fp=SMD06_HAND)
            series_wire(urs, pin, resistor, net, net + "_IC")
        dbus = place("phantasm:CDSOT23-SM712", "D_SYNC", "CDSOT23-SM712", 309.88, 248.92,
                     fp="Package_TO_SOT_SMD:SOT-23")
        to_label(dbus, "1", "SYNC_A"); to_label(dbus, "2", "SYNC_B")
        to_power(dbus, "3", GND)
        for ref, y in (("J3A", 202.0), ("J3B", 229.87)):
            connector = place("Connector_Generic:Conn_01x04", ref,
                              "TBC05-04-1-G-G SYNC " + ("in" if ref == "J3A" else "out"),
                              355.6, y, in_bom=False,
                              fp="phantasm:TerminalBlock_GCT_TBC05-04-1-G-G")
            to_label(connector, "1", "SYNC_A"); to_label(connector, "2", "SYNC_B")
            to_power(connector, "3", GND); to_label(connector, "4", "SHIELD")
        rterm = place("Device:R", "R_TERM", "120R 1% 0.25W", 279.4, 251.46, fp=SMD08_HAND)
        jpterm = place("Jumper:SolderJumper_2_Open", "JP_TERM", "TERM endpoints only", 274.32, 271.78,
                       fp="Jumper:SolderJumper-2_P1.3mm_Open_RoundedPad1.0x1.5mm", in_bom=False)
        to_label(rterm, "1", "SYNC_A"); to_label(rterm, "2", "TERM_LINK")
        to_label(jpterm, "1", "TERM_LINK"); to_label(jpterm, "2", "SYNC_B")
    JPS = place("Jumper:SolderJumper_2_Open", "JP_SHLD", "shield gnd (master only)", 363.22, 248.92,
                fp="Jumper:SolderJumper-2_P1.3mm_Open_RoundedPad1.0x1.5mm", in_bom=False)
    to_label(JPS, "1", "SHIELD"); to_power(JPS, "2", GND)

    # ============================================================ BLOCK 4: STRAPS
    b.text((25, 245), "ID STRAPS / MASTER_EN " + ("PULL-UP" if revision == "1.2" else "PULL-DOWN"), 2.2)
    RMEN = place("Device:R", "R_MEN", "10k", 76.2, 261.62, fp=SMD06)
    to_power(RMEN, "1", V3 if revision == "1.2" else GND); to_label(RMEN, "2", "MASTER_EN")
    JID0 = place("Jumper:SolderJumper_2_Open", "JP_ID0", "ID0->GND", 127.0, 274.32,
                 fp="Jumper:SolderJumper-2_P1.3mm_Open_RoundedPad1.0x1.5mm", in_bom=False)
    to_label(JID0, "1", "ID0"); to_power(JID0, "2", GND)
    JID1 = place("Jumper:SolderJumper_2_Open", "JP_ID1", "ID1->GND", 152.4, 274.32,
                 fp="Jumper:SolderJumper-2_P1.3mm_Open_RoundedPad1.0x1.5mm", in_bom=False)
    to_label(JID1, "1", "ID1"); to_power(JID1, "2", GND)
    # ID2 strap (pin 23) — read only by the N=8 firmware profile (R-ID-1)
    JID2 = place("Jumper:SolderJumper_2_Open", "JP_ID2", "ID2->GND (N=8)", 177.8, 274.32,
                 fp="Jumper:SolderJumper-2_P1.3mm_Open_RoundedPad1.0x1.5mm", in_bom=False)
    to_label(JID2, "1", "ID2"); to_power(JID2, "2", GND)

    # ============================================================ POWER FLAGS (ERC)
    b.text((220, 262), "POWER FLAGS (ERC)", 2.2)
    flag_index = 0


    def flag(net, x, y):
        nonlocal flag_index
        flag_index += 1
        pwr_sym(net, x, y)
        ref = f"#FLG{flag_index:02d}"
        fl = B.Symbol("power:PWR_FLAG", ref, "PWR_FLAG", x, y - 5.08)
        b.place(fl)
        b.wire((x, y), fl.pin("1"))

    fy = 274.32
    flag("+5V_RAW", 246.38, fy)
    flag("+5V_LOGIC", 284.48, fy)
    flag(GND, 317.5, fy)

    # ---------------------------------------------------------------- write files
    os.makedirs(out, exist_ok=True)
    atomic_write_text(sch, b.dumps())

    lib_lines = ['(kicad_symbol_lib', f'\t(version {sexp.SYMBOL_LIB_FORMAT})',
                 '\t(generator "phantasm-gen")',
                 f'\t(generator_version "{sexp.GENERATOR_VERSION}")']
    for lib_id in sorted(b.lib_defs):
        if lib_id.startswith("phantasm:"):
            node = copy.deepcopy(b.lib_defs[lib_id]); node[1] = lib_id.split(":", 1)[1]
            lib_lines.append(sexp.dumps(node, indent=1))
    lib_lines.append(')')
    atomic_write_text(os.path.join(out, "phantasm.kicad_sym"), "\n".join(lib_lines) + "\n")

    atomic_write_text(os.path.join(out, "sym-lib-table"), '(sym_lib_table\n\t(version 7)\n'
            '\t(lib (name "phantasm")(type "KiCad")(uri "${KIPRJMOD}/phantasm.kicad_sym")'
            '(options "")(descr "PHANTASM custom symbols"))\n)\n')

    PRO = os.path.join(out, "phantasm.kicad_pro")
    write_project(PRO, b.uuid, unplaced=None)

    print("wrote files  symbols:", len(b.symbols), "wires:", len(b.wires),
          "labels:", len(b.labels), "texts:", len(b.texts))


if __name__ == "__main__":
    args = parse_args()
    main(force=args.force, revision=args.revision)
