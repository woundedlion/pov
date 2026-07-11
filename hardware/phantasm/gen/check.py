"""Gate: export netlist via kicad-cli and verify named nets match spec section 10
by component-ref membership (pin-number agnostic, so it catches any short/break)."""
import subprocess
import sys
import os
import tempfile
sys.path.insert(0, os.path.dirname(__file__))
import sexp
import fab

SCH = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "phantasm.kicad_sch")
KCLI = fab.find_kicad_cli()

# expected: named net -> set of component refs (per spec section 10)
EXPECT = {
    "+5V_RAW":    {"F1", "Q_REV"},
    "+5V_LOGIC":  {"FB", "C_IN", "R_LF", "C_DEC1", "C_DEC2", "U1", "U_MCU"},
    "+3V3":       {"U_MCU", "R_MEN", "R_ID0", "J4"},
    "DATA_IN":    {"U_MCU", "U1"},
    "CLK_IN":     {"U_MCU", "U1"},
    "FRAME_SYNC": {"U_MCU", "U1", "R1", "R2", "C_SYNC"},
    "MASTER_EN":  {"U_MCU", "U1", "R_MEN", "J4"},
    "DATA":       {"R_D1", "J2"},
    "CLK":        {"R_D2", "J2"},
    "SYNC_BUS":   {"R_S", "R1", "R_PD", "D_BUS", "J3A", "J3B"},
    "ID0":        {"U_MCU", "R_ID0", "JP_ID0"},
    "ID1":        {"U_MCU", "JP_ID1"},
    "ID2":        {"U_MCU", "JP_ID2"},
    "SHIELD":     {"J3A", "J3B", "JP_SHLD"},
    "SERIAL1_TX": {"U_MCU", "J4"},
    "GND":        {"C_IN", "C_LF", "C_DEC1", "C_DEC2", "C_SYNC", "D_BUS",
                   "J1", "J2", "J3A", "J3B", "J4", "JP_ID0", "JP_ID1", "JP_ID2",
                   "JP_SHLD", "R2", "R_PD", "U1", "U_MCU"},
}

_netfd, NET = tempfile.mkstemp(suffix=".net")
os.close(_netfd)
try:
    subprocess.run([KCLI, "sch", "export", "netlist", "--format", "kicadsexpr",
                    "-o", NET, SCH], check=True, capture_output=True, text=True)
    root = sexp.parse(open(NET, encoding="utf-8").read())[0]
except subprocess.CalledProcessError as e:
    sys.stderr.write(e.stderr or "")
    raise
finally:
    if os.path.exists(NET):
        os.remove(NET)


def F(n, k):
    return [c for c in n if isinstance(c, list) and c and c[0] == k]


got = {}
for nb in F(root, "nets"):
    for net in F(nb, "net"):
        name = sexp._val(net, "name")[0].lstrip("/")
        got[name] = set(sexp._val(nd, "ref")[0] for nd in F(net, "node"))

ok = True
for name, refs in EXPECT.items():
    g = got.get(name, set())
    if g != refs:
        ok = False
        print(f"FAIL {name}\n   expected {sorted(refs)}\n   got      {sorted(g)}")
# also flag any named (non-auto) net that merged two expected nets
for name, g in got.items():
    if name.startswith(("unconnected", "Net-", "PWR")):
        continue
    # intentional intermediate-node labels (named so autoplacers show readable nets)
    INTERMEDIATE = {"DATA_SRC", "CLK_SRC", "SYNC_SRC", "+5V_IN", "+5V_PROT", "LF_DAMP"}
    if name not in EXPECT and name not in INTERMEDIATE:
        print(f"NOTE extra named net {name}: {sorted(g)}")
print("NETLIST OK" if ok else "NETLIST MISMATCH")
sys.exit(0 if ok else 1)
