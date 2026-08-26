"""Gate: export the schematic's netlist and verify its named-net partition.

Nodes are keyed on (ref, pin), so a connector or IC pinout permutation fails the
gate as loudly as a short or a break does. gen/tests/test_check.py applies the
same table to the committed board's pad nets; that is the leg CI runs.
"""
import argparse
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))
import sexp
from kicad_common import F, export_netlist, kicad_cli

# Refs whose terminals are electrically interchangeable: plain resistors,
# non-polarized caps, the fuse, the bead, solder jumpers. Their pin numbers
# follow schematic geometry rather than the spec, so these nodes are keyed on
# ref alone. Each such ref is expected on two nets, so a bridged or open
# terminal still fails. D_BUS (CDSOD323-T05L) is a unidirectional TVS and is
# pin-keyed: pin 1 cathode on SYNC_BUS, pin 2 anode on GND.
SYMMETRIC = {"R1", "R2", "R_D1", "R_D2", "R_S", "R_PD", "R_MEN", "R_LF",
             "C_LF", "C_DEC1", "C_DEC2", "C_SYNC", "F1", "FB",
             "JP_ID0", "JP_ID1", "JP_ID2", "JP_SHLD"}


def node_key(ref, pin):
    return ref if ref in SYMMETRIC else f"{ref}.{pin}"


# Expected named net -> set of node keys across the electrical requirements.
EXPECT = {
    "+5V_IN":     {"J1.1", "F1"},
    "+5V_RAW":    {"F1", "Q_REV.3"},
    "+5V_PROT":   {"Q_REV.2", "FB"},
    "+5V_LOGIC":  {"FB", "R_LF", "C_IN.1", "C_DEC1", "C_DEC2",
                   "U1.14", "U_MCU.VIN"},
    "LF_DAMP":    {"R_LF", "C_LF"},
    "+3V3":       {"U_MCU.3V3", "R_MEN", "J4.1"},
    "DATA_IN":    {"U_MCU.11", "U1.2"},
    "CLK_IN":     {"U_MCU.13", "U1.5"},
    "DATA_SRC":   {"U1.3", "R_D1"},
    "CLK_SRC":    {"U1.6", "R_D2"},
    "DATA":       {"R_D1", "J2.1"},
    "CLK":        {"R_D2", "J2.2"},
    "FRAME_SYNC": {"U_MCU.3", "U1.9", "R1", "R2", "C_SYNC"},
    "MASTER_EN":  {"U_MCU.5", "U1.10", "U1.13", "R_MEN", "J4.3"},
    "SYNC_SRC":   {"U1.8", "R_S"},
    "SYNC_BUS":   {"R_S", "R1", "R_PD", "D_BUS.1", "J3A.1", "J3B.1"},
    "SYNC_PULLDOWN": {"U1.11", "R_PD"},
    "ID0":        {"U_MCU.21", "JP_ID0"},
    "ID1":        {"U_MCU.22", "JP_ID1"},
    "ID2":        {"U_MCU.23", "JP_ID2"},
    "SHIELD":     {"J3A.3", "J3B.3", "JP_SHLD"},
    "SERIAL1_TX": {"U_MCU.1", "J4.4"},
    "GND":        {"C_IN.2", "C_LF", "C_DEC1", "C_DEC2", "C_SYNC", "D_BUS.2",
                   "J1.2", "J2.3", "J3A.2", "J3B.2", "J4.2",
                   "JP_ID0", "JP_ID1", "JP_ID2", "JP_SHLD", "Q_REV.1", "R2",
                   "U1.1", "U1.4", "U1.7", "U1.12", "U_MCU.GND"},
}


def netlist_nets(root):
    """Parsed netlist -> {net name: set of node keys}."""
    got = {}
    for nb in F(root, "nets"):
        for net in F(nb, "net"):
            name = sexp.val(net, "name")
            if not name:
                continue
            keys = set()
            for nd in F(net, "node"):
                ref = sexp.val(nd, "ref")
                pin = sexp.val(nd, "pin")
                if not ref or not pin:
                    continue
                keys.add(node_key(ref[0], pin[0]))
            got[name[0].lstrip("/")] = keys
    return got


def check(got):
    """Report every net that differs from EXPECT; return True when all match.

    A named net absent from EXPECT is printed as an advisory NOTE and does not
    move the verdict: the gate partitions the nets it knows, it does not close
    the set.
    """
    ok = True
    for name, keys in sorted(EXPECT.items()):
        g = got.get(name, set())
        if g != keys:
            ok = False
            print(f"FAIL {name}\n   missing {sorted(keys - g)}"
                  f"\n   unexpected {sorted(g - keys)}")
    # note any named net outside the spec table; the skipped prefixes are
    # KiCad's own auto-generated names
    for name, g in got.items():
        if name.startswith(("unconnected", "Net-")):
            continue
        if name not in EXPECT:
            print(f"NOTE extra named net {name}: {sorted(g)}")
    return ok


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("schematic", nargs="?", default=os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "phantasm.kicad_sch"))
    args = parser.parse_args(argv)
    ok = check(netlist_nets(export_netlist(kicad_cli(), args.schematic)))
    print("NETLIST OK" if ok else "NETLIST MISMATCH")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
