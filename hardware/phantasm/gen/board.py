"""Generate phantasm.kicad_sch — PHANTASM per-segment carrier board.

Run from this directory:  python board.py
Writes the project files into the parent (hardware/phantasm/).
Requires KiCad's stock symbol libs (see sexp.KICAD_SHARE / env KICAD_SYMBOL_DIR).

Layout: left-to-right signal flow in labelled blocks. Power distribution uses
visible horizontal rail wires with vertical component drops + junctions; signal
buses between the Teensy, level shifter and connectors use net labels (ports).
"""
import copy
import os
import builder as B
import sexp

b = B.Builder("PHANTASM Segment Board  -  per-segment carrier (x4, strap-selected role)",
              paper="A3")
GND = "GND"; V3 = "+3V3"


# ---------------------------------------------------------------- custom symbols
def make_power(net):
    node = copy.deepcopy(sexp.get_symbol("power", "+5V"))
    B._rename_subsymbols(node, "+5V", net)
    for c in node:
        if isinstance(c, list) and c and c[0] == "property" and c[1] == "Value":
            c[2] = net
    return b.register_custom(node, f"phantasm:{net}")


def make_teensy():
    """Compact Teensy 4.0 symbol showing only the pins this board uses (the other
    ~17 pads are unconnected on this design and omitted for readability). Pin
    NUMBER = the Teensy pad label; assign a Teensy footprint and verify the pad
    map before PCB layout."""
    # (display-name, number, electrical-type) — left side then right side
    LEFT = [("D11/MOSI", "11"), ("D13/SCK", "13"), ("D3", "3"), ("D5", "5"),
            ("D1/TX1", "1"), ("D21", "21"), ("D22", "22"), ("D23", "23")]
    RIGHT = [("VIN", "VIN"), ("3V3", "3V3"), ("GND", "GND")]
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

    half = topL + pitch
    text = (
        '(symbol "phantasm:Teensy4.0"\n'
        '\t(pin_names (offset 1.016))\n'
        '\t(exclude_from_sim no) (in_bom yes) (on_board yes)\n'
        f'\t(property "Reference" "U" (at {-bodyx} {half + 2.54} 0)\n'
        '\t\t(effects (font (size 1.27 1.27)) (justify left)))\n'
        f'\t(property "Value" "Teensy4.0" (at {-bodyx} {-half - 2.54} 0)\n'
        '\t\t(effects (font (size 1.27 1.27)) (justify left)))\n'
        '\t(property "Footprint" "" (at 0 0 0)\n'
        '\t\t(effects (font (size 1.27 1.27)) (hide yes)))\n'
        '\t(property "Datasheet" "" (at 0 0 0)\n'
        '\t\t(effects (font (size 1.27 1.27)) (hide yes)))\n'
        '\t(symbol "Teensy4.0_0_1"\n'
        f'\t\t(rectangle (start {-bodyx} {half}) (end {bodyx} {-half})\n'
        '\t\t\t(stroke (width 0.254) (type default)) (fill (type background))))\n'
        '\t(symbol "Teensy4.0_1_1"\n'
        + "\n".join(pins) + ")\n)\n")
    node = sexp.parse(text)[0]
    return b.register_custom(node, "phantasm:Teensy4.0"), name2num


make_power("+5V_RAW"); make_power("+5V_LOGIC")
TEENSY, TPN = make_teensy()

for lib, name in [("power", "GND"), ("power", "+3V3"), ("power", "PWR_FLAG"),
                  ("Device", "R"), ("Device", "C"), ("Device", "C_Polarized"),
                  ("Device", "FerriteBead"), ("Device", "Fuse"), ("Device", "D_TVS"),
                  ("Device", "D_Schottky"),
                  ("Connector_Generic", "Conn_01x02"),
                  ("Connector_Generic", "Conn_01x03"),
                  ("Connector_Generic", "Conn_01x04"),
                  ("Jumper", "SolderJumper_2_Open"), ("74xx", "74AHCT125")]:
    b.ensure_lib(lib, name)


# ---------------------------------------------------------------- helpers
def place(lib, ref, val, x, y, rot=0, unit=1, fp="", dnp=False):
    return b.place(B.Symbol(lib, ref, val, x, y, rot=rot, unit=unit,
                            footprint=fp, dnp=dnp))


def hw(x1, x2, y):
    b.wire((x1, y), (x2, y))


def vw(x, y1, y2):
    b.wire((x, y1), (x, y2))


def pwr_sym(kind, x, y):
    ref = "#PWR"
    lib = {"GND": "power:GND", "+3V3": "power:+3V3"}.get(kind, f"phantasm:{kind}")
    b.place(B.Symbol(lib, ref, kind, x, y))


def to_power(s, num, kind, length=2.54):
    end = b.stub(s, num, length)
    pwr_sym(kind, end[0], end[1])


def to_label(s, num, name, length=3.81):
    end = b.stub(s, num, length)
    b.label(end, name)


def to_gnd_down(s, num, yrail):
    """Drop a (bottom) pin straight down onto a GND rail at yrail."""
    t = s.pin(num)
    vw(t[0], t[1], yrail)
    b.junction((t[0], yrail))


def to_rail_up(s, num, yrail):
    """Drop a (top) pin straight up onto a rail at yrail."""
    t = s.pin(num)
    vw(t[0], t[1], yrail)
    b.junction((t[0], yrail))


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
C06 = "Capacitor_SMD:C_0603_1608Metric"


# ============================================================ BLOCK 1: POWER
b.text((25, 25), "POWER ENTRY / PROTECTION / RAIL FILTER  (logic ~0.15 A; LED 4.3 A off-board)", 2.2)
Y_LOG, Y_GND = 60.96, 96.52

# --- light logic feed only; LED 4.3 A power is delivered off-board (spec 2.3) ---
# Series chain on the rail line: J1 -> F1 -> D1 (reverse block) -> FB -> +5V_LOGIC.
J1 = place("Connector_Generic:Conn_01x02", "J1", "+5V IN ~1A", 25.4, 60.96, rot=180,
           fp="Connector_PinHeader_2.54mm:PinHeader_1x02_P2.54mm_Vertical")
# Small logic-only fuse/PTC (R-PWR-8) — the 4.3 A strip current never flows here.
F1 = place("Device:Fuse", "F1", "0.75A", 40.64, 60.96, rot=90,
           fp="Fuse:Fuse_1206_3216Metric")
# Series Schottky reverse-protect (R-PWR-7): at ~0.15 A the ideal-diode P-FET is
# unnecessary, so a small SOD-123 Schottky replaces the TO-220 FET + gate resistor.
# Ref kept as Q_REV per spec 9/10; pin2 = anode (source side), pin1 = cathode (load).
D1 = place("Device:D_Schottky", "Q_REV", "B5819W", 53.34, 60.96, rot=180,
           fp="Diode_SMD:D_SOD-123")
FB = place("Device:FerriteBead", "FB", "600R@100MHz", 68.58, 60.96, rot=90,
           fp="Inductor_SMD:L_1206_3216Metric")

# J1.1(+5V) -> F1.1 ; J1.2 -> GND ; name the pre-fuse entry node +5V_IN
b.wire(J1.pin("1"), F1.pin("1"))
b.label(F1.pin("1"), "+5V_IN")
to_power(J1, "2", GND)
# F1.2 -> D1 anode (pin2) ; label the protected-input node +5V_RAW
b.wire(F1.pin("2"), D1.pin("2"))
b.label(F1.pin("2"), "+5V_RAW")
# D1 cathode (pin1) -> FB.1 (reverse-blocked) ; name the post-diode node +5V_PROT
b.wire(D1.pin("1"), FB.pin("1"))
b.label(D1.pin("1"), "+5V_PROT")

# --- +5V_LOGIC rail (post-bead) and its drops: C_IN, R_LF/C_LF damper, C_DEC1/2 ---
LOG_L, LOG_R = 76.2, 219.71
b.wire(FB.pin("2"), (LOG_L, Y_LOG))
hw(LOG_L, LOG_R, Y_LOG)
pwr_sym("+5V_LOGIC", LOG_R, Y_LOG - 5.08); vw(LOG_R, Y_LOG - 5.08, Y_LOG)
# C_IN: the card's only electrolytic, now on the logic feed (R-PWR-3/6, spec 10).
CIN = place("Device:C_Polarized", "C_IN", "100uF", 88.9, 78.74,
            fp="Capacitor_THT:CP_Radial_D8.0mm_P3.50mm")
to_rail_up(CIN, "1", Y_LOG)
to_gnd_down(CIN, "2", Y_GND)
# R_LF + C_LF bead-LC damper (R-PWR-5): LOG -> R_LF -> CLF_NODE -> C_LF -> GND
RLF = place("Device:R", "R_LF", "1R5", 109.22, 71.12, fp=SMD08)
CLF = place("Device:C", "C_LF", "22uF", 109.22, 86.36, fp="Capacitor_SMD:C_1206_3216Metric")
to_rail_up(RLF, "1", Y_LOG)
b.wire(RLF.pin("2"), CLF.pin("1"))   # CLF_NODE
b.label(CLF.pin("1"), "LF_DAMP")
to_gnd_down(CLF, "2", Y_GND)
# decoupling caps (R-PWR-4): one at Teensy VIN, one at U1 Vcc
CD1 = place("Device:C", "C_DEC1", "0.1uF", 140.97, 78.74, fp=C06)
CD2 = place("Device:C", "C_DEC2", "0.1uF", 166.37, 78.74, fp=C06)
for c in (CD1, CD2):
    to_rail_up(c, "1", Y_LOG)
    to_gnd_down(c, "2", Y_GND)

# --- GND rail spanning the cap drops + its GND symbol at the right end ---
GND_L, GND_R = 88.9, 166.37
hw(GND_L, GND_R, Y_GND)
pwr_sym(GND, GND_R, Y_GND)

# ============================================================ BLOCK 2: LOGIC
b.text((25, 116), "TEENSY 4.0  +  74AHCT125 LEVEL SHIFTER  ->  LED STRIP", 2.2)
# --- Teensy (compact symbol, only used pins) ---
U = place(TEENSY, "U_MCU", "Teensy4.0", 60.96, 165.1)
tn = lambda d: TPN[d]
to_power(U, tn("VIN"), "+5V_LOGIC")
to_power(U, tn("3V3"), V3)
to_power(U, tn("GND"), GND)
to_label(U, tn("D11/MOSI"), "DATA_IN")
to_label(U, tn("D13/SCK"), "CLK_IN")
to_label(U, tn("D3"), "FRAME_SYNC")
to_label(U, tn("D5"), "MASTER_EN")
to_label(U, tn("D21"), "ID0")
to_label(U, tn("D22"), "ID1")
to_label(U, tn("D23"), "ID2")   # read by the N=8 firmware profile
to_label(U, tn("D1/TX1"), "SERIAL1_TX")

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
RS = place("Device:R", "R_S", "100R", 177.8, 198.12, rot=270, fp=SMD08)
# ch A (DATA) — source stub (U1->R_D1) named DATA_SRC; post-term net is DATA
to_label(U1A, "2", "DATA_IN"); to_power(U1A, "1", GND); series_wire(U1A, "3", RD1, "DATA", "DATA_SRC")
# ch B (CLK)
to_label(U1B, "5", "CLK_IN"); to_power(U1B, "4", GND); series_wire(U1B, "6", RD2, "CLK", "CLK_SRC")
# ch C (SYNC)
to_label(U1C, "9", "FRAME_SYNC"); to_label(U1C, "10", "MASTER_EN"); series_wire(U1C, "8", RS, "SYNC_BUS", "SYNC_SRC")
# ch D (unused)
to_power(U1D, "12", GND); to_power(U1D, "13", "+5V_LOGIC"); b.nc(U1D.pin("11"))
# power unit
to_power(U1E, "14", "+5V_LOGIC"); to_power(U1E, "7", GND)

# --- J2 strip SIGNAL out (3-pin, no power): DI / CI / SIG_GND (R-CON-1) ---
# Strip 5 V/GND are injected off-board (spec 2.3); SIG_GND is the card's logic GND,
# landed on the strip GND pin at the load end (the off-board ground star).
J2 = place("Connector_Generic:Conn_01x03", "J2", "LED sig DI/CI/SIG_GND", 281.94, 165.1,
           fp="Connector_PinHeader_2.54mm:PinHeader_1x03_P2.54mm_Vertical")
to_label(J2, "1", "DATA"); to_label(J2, "2", "CLK"); to_power(J2, "3", GND)

# ============================================================ BLOCK 3: SYNC
b.text((220, 180), "SYNC BUS  -  RX DIVIDER + DAISY (Belden 8451)", 2.2)
# divider: SYNC_BUS -> R1 -> node(FRAME_SYNC) -> R2 -> GND ; C_SYNC at node
R1 = place("Device:R", "R1", "10k", 246.38, 205.74, fp=SMD06)
R2 = place("Device:R", "R2", "15k", 246.38, 228.6, fp=SMD06)
CSY = place("Device:C", "C_SYNC", "220pF", 269.24, 217.17, rot=90, fp=C06)
to_label(R1, "1", "SYNC_BUS")
nd = R1.pin("2")
b.wire(nd, R2.pin("1"))             # divider node
b.label(nd, "FRAME_SYNC")
b.wire(R2.pin("1"), CSY.pin("1"))   # C_SYNC onto node
b.junction(R2.pin("1"))
to_power(R2, "2", GND); to_power(CSY, "2", GND)
# bus idle pulldown + optional TVS
RPD = place("Device:R", "R_PD", "10k", 292.1, 205.74, fp=SMD06)
to_label(RPD, "1", "SYNC_BUS"); to_power(RPD, "2", GND)
DBUS = place("Device:D_TVS", "D_BUS", "TVS", 292.1, 231.14, dnp=True,
             fp="Diode_SMD:D_SOD-323")
to_label(DBUS, "1", "SYNC_BUS"); to_power(DBUS, "2", GND)
# daisy connectors
J3A = place("Connector_Generic:Conn_01x03", "J3A", "SYNC in", 330.2, 205.74,
            fp="Connector_PinHeader_2.54mm:PinHeader_1x03_P2.54mm_Vertical")
J3B = place("Connector_Generic:Conn_01x03", "J3B", "SYNC out", 330.2, 233.68,
            fp="Connector_PinHeader_2.54mm:PinHeader_1x03_P2.54mm_Vertical")
for J in (J3A, J3B):
    to_label(J, "1", "SYNC_BUS"); to_power(J, "2", GND); to_label(J, "3", "SHIELD")
JPS = place("Jumper:SolderJumper_2_Open", "JP_SHLD", "shield gnd (master only)", 363.22, 248.92,
            fp="Jumper:SolderJumper-2_P1.3mm_Open_RoundedPad1.0x1.5mm")
to_label(JPS, "1", "SHIELD"); to_power(JPS, "2", GND)

# ============================================================ BLOCK 4: STRAPS/DBG
b.text((25, 245), "ID STRAPS / MASTER_EN PULL-UP / DEBUG", 2.2)
RMEN = place("Device:R", "R_MEN", "10k", 76.2, 261.62, fp=SMD06)
to_power(RMEN, "1", V3); to_label(RMEN, "2", "MASTER_EN")
RID0 = place("Device:R", "R_ID0", "10k DNP", 101.6, 261.62, dnp=True, fp=SMD06)
to_power(RID0, "1", V3); to_label(RID0, "2", "ID0")
JID0 = place("Jumper:SolderJumper_2_Open", "JP_ID0", "ID0->GND", 127.0, 274.32,
             fp="Jumper:SolderJumper-2_P1.3mm_Open_RoundedPad1.0x1.5mm")
to_label(JID0, "1", "ID0"); to_power(JID0, "2", GND)
JID1 = place("Jumper:SolderJumper_2_Open", "JP_ID1", "ID1->GND", 152.4, 274.32,
             fp="Jumper:SolderJumper-2_P1.3mm_Open_RoundedPad1.0x1.5mm")
to_label(JID1, "1", "ID1"); to_power(JID1, "2", GND)
# ID2 strap (pin 23) — read only by the N=8 firmware profile (R-ID-1)
JID2 = place("Jumper:SolderJumper_2_Open", "JP_ID2", "ID2->GND (N=8)", 177.8, 274.32,
             fp="Jumper:SolderJumper-2_P1.3mm_Open_RoundedPad1.0x1.5mm")
to_label(JID2, "1", "ID2"); to_power(JID2, "2", GND)
J4 = place("Connector_Generic:Conn_01x04", "J4", "debug", 38.1, 266.7,
           fp="Connector_PinHeader_2.54mm:PinHeader_1x04_P2.54mm_Vertical")
to_power(J4, "1", V3); to_power(J4, "2", GND)
to_label(J4, "3", "MASTER_EN"); to_label(J4, "4", "SERIAL1_TX")

# ============================================================ POWER FLAGS (ERC)
b.text((220, 262), "POWER FLAGS (ERC)", 2.2)
def flag(net, x, y):
    pwr_sym(net, x, y)
    fl = B.Symbol("power:PWR_FLAG", "#FLG", "PWR_FLAG", x, y - 5.08)
    b.place(fl)
    b.wire((x, y), fl.pin("1"))

fy = 274.32
flag("+5V_RAW", 246.38, fy)
flag("+5V_LOGIC", 284.48, fy)
flag(GND, 317.5, fy)

# ---------------------------------------------------------------- write files
OUT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
open(os.path.join(OUT, "phantasm.kicad_sch"), "w", encoding="utf-8").write(b.dumps())

lib_lines = ['(kicad_symbol_lib', '\t(version 20251024)',
             '\t(generator "phantasm-gen")', '\t(generator_version "10.0")']
for lib_id in sorted(b.lib_defs):
    if lib_id.startswith("phantasm:"):
        node = copy.deepcopy(b.lib_defs[lib_id]); node[1] = lib_id.split(":", 1)[1]
        lib_lines.append(sexp.dumps(node, indent=1))
lib_lines.append(')')
open(os.path.join(OUT, "phantasm.kicad_sym"), "w", encoding="utf-8").write("\n".join(lib_lines) + "\n")

open(os.path.join(OUT, "sym-lib-table"), "w", encoding="utf-8").write(
    '(sym_lib_table\n\t(version 7)\n'
    '\t(lib (name "phantasm")(type "KiCad")(uri "${KIPRJMOD}/phantasm.kicad_sym")'
    '(options "")(descr "PHANTASM custom symbols"))\n)\n')

open(os.path.join(OUT, "phantasm.kicad_pro"), "w", encoding="utf-8").write(
    # design_settings.rules.min_clearance > 0 so Quilter accepts the project on upload
    # (it rejects 0). KiCad re-zeroes it whenever the project is opened in the GUI, so
    # run gen/heal_clearance.py before any Quilter upload to restore it.
    '{\n  "board": { "design_settings": { "rules": { "min_clearance": 0.2 } } },\n'
    '  "boards": [],\n  "cvpcb": { "equivalence_files": [] },\n'
    '  "libraries": { "pinned_footprint_libs": [], "pinned_symbol_libs": [] },\n'
    '  "meta": { "filename": "phantasm.kicad_pro", "version": 3 },\n'
    '  "net_settings": { "classes": [ { "name": "Default", "clearance": 0.2,\n'
    '    "track_width": 0.3, "via_diameter": 0.6, "via_drill": 0.3 } ] },\n'
    '  "pcbnew": { "page_layout_descr_file": "" },\n'
    '  "schematic": { "annotate_start_num": 0,\n'
    '    "drawing": { "default_line_thickness": 6.0, "label_size_ratio": 0.375 } },\n'
    '  "sheets": [ [ "' + b.uuid + '", "Root" ] ],\n  "text_variables": {}\n}\n')

print("wrote files  symbols:", len(b.symbols), "wires:", len(b.wires),
      "labels:", len(b.labels), "texts:", len(b.texts))
