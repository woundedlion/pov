"""Board generation data and shared stackup constants."""
from dataclasses import dataclass
import os
import sys
from revision_policy import DEFAULT_OUT, RevisionLayout
import sexp
from constraints import (MAX_BOARD_WIDTH_MM)
from kicad_common import property_value
from kicad_common import (fmt, F)


HERE = os.path.dirname(os.path.abspath(__file__))


OUT = DEFAULT_OUT


TERMINAL_LIBIDS = tuple(
    f"phantasm:TerminalBlock_GCT_TBC05-0{pins}-1-G-G" for pins in (2, 3, 4))


SCH = os.path.join(OUT, "phantasm.kicad_sch")


PCB_FILE = "phantasm.kicad_pcb"


DRAFT_FILE = "phantasm-draft.kicad_pcb"


UNPLACED_FILE = PCB_FILE


GENERATED_LAND_DEVIATIONS = {("lib_footprint_mismatch", "D_BUS"): 1}


GENERATED_OPTIONAL_DEVIATIONS = {("lib_footprint_mismatch", "J1"): 1}


def silk_revision(revision):
    """Bottom-silkscreen stamp for the selected board revision."""
    return f"Phantasm Rev {revision}"


UNPLACED_REASON = (
    "Regeneration replaces this revision's board, including any placement\n"
    "  and routing saved in KiCad.")


PCB_W_MAX = MAX_BOARD_WIDTH_MM


PCB_W = 32.0  # board width (mm), fixed


if PCB_W > PCB_W_MAX:
    raise ValueError(f"PCB_W {fmt(PCB_W)} mm exceeds the R-MECH-6 cap of "
                     f"{fmt(PCB_W_MAX)} mm")


QUILTER_LENGTH = 58.28


TEENSY_LIBRARY_REASON = (
    "This revision's committed phantasm.kicad_pcb resolves its Teensy pads\n"
    "  against this library file.")


TEENSY_MODEL_PATH = "${KIPRJMOD}/phantasm.pretty/Teensy4.0.wrl"


TEENSY_LIBID = "phantasm:Teensy4.0"


MOUNTING_HOLE_FOOTPRINT = "MountingHole:MountingHole_2.7mm_M2.5"


MOUNTING_KEEPOUT_HALF_SIDE = 2.7


MOUNTING_HOLE_INSET = 3.5


MOUNTING_COURTYARD_RADIUS = 2.95


COURTYARD_MARGIN = 0.05


MOUNTING_RESERVE_HALF_SIDE = max(MOUNTING_KEEPOUT_HALF_SIDE,
                                 MOUNTING_COURTYARD_RADIUS + COURTYARD_MARGIN)


COURTYARD_LAYERS = ("F.CrtYd", "B.CrtYd")


COPPER_LAYERS = ((0, "F.Cu", "signal"), (4, "In1.Cu", "power"),
                 (6, "In2.Cu", "power"), (2, "B.Cu", "signal"))


TECHNICAL_LAYERS = ((9, "F.Adhes", "user"), (11, "B.Adhes", "user"),
                    (13, "F.Paste", "user"), (15, "B.Paste", "user"),
                    (5, "F.SilkS", "user"), (7, "B.SilkS", "user"),
                    (1, "F.Mask", "user"), (3, "B.Mask", "user"),
                    (17, "Dwgs.User", "user"), (19, "Cmts.User", "user"),
                    (21, "Eco1.User", "user"), (23, "Eco2.User", "user"),
                    (25, "Edge.Cuts", "user"), (27, "Margin", "user"),
                    (31, "F.CrtYd", "user"), (29, "B.CrtYd", "user"),
                    (35, "F.Fab", "user"), (33, "B.Fab", "user"))


GROUND_NET = "GND"


LAYER_NAMES = {"In1.Cu": "GND", "In2.Cu": "Ground"}


BOARD_THICKNESS_MM = 1.6


STACKUP = (
    '(layer "F.SilkS" (type "Top Silk Screen"))',
    '(layer "F.Paste" (type "Top Solder Paste"))',
    '(layer "F.Mask" (type "Top Solder Mask") (thickness 0.01))',
    '(layer "F.Cu" (type "copper") (thickness 0.035))',
    '(layer "dielectric 1" (type "prepreg") (thickness 0.2104) (material "FR4")'
    ' (epsilon_r 4.5) (loss_tangent 0.02))',
    '(layer "In1.Cu" (type "copper") (thickness 0.0152))',
    '(layer "dielectric 2" (type "core") (thickness 1.065) (material "FR4")'
    ' (epsilon_r 4.5) (loss_tangent 0.02))',
    '(layer "In2.Cu" (type "copper") (thickness 0.0152))',
    '(layer "dielectric 3" (type "prepreg") (thickness 0.2104) (material "FR4")'
    ' (epsilon_r 4.5) (loss_tangent 0.02))',
    '(layer "B.Cu" (type "copper") (thickness 0.035))',
    '(layer "B.Mask" (type "Bottom Solder Mask") (thickness 0.01))',
    '(layer "B.Paste" (type "Bottom Solder Paste"))',
    '(layer "B.SilkS" (type "Bottom Silk Screen"))',
    '(copper_finish "ENIG")',
    '(dielectric_constraints no)',
)


def copper_layer_names():
    return tuple(name for _, name, _ in COPPER_LAYERS)


def build_nets(nlroot):
    pad_net = {}            # (ref, pad) -> netname
    names = []
    for nb in F(nlroot, "nets"):
        for net in F(nb, "net"):
            name = sexp.val(net, "name")
            if not name:
                continue
            names.append(name[0])
            for nd in F(net, "node"):
                ref = sexp.val(nd, "ref")
                pin = sexp.val(nd, "pin")
                if not ref or not pin:
                    continue
                pad_net[(ref[0], pin[0])] = name[0]
    # stable net ids: 0 = "", then unique names in first-seen order
    netid = {"": 0}
    for n in names:
        if n not in netid:
            netid[n] = len(netid)
    return pad_net, netid


def build_paths(nlroot):
    """ref -> footprint (path) linking to the schematic symbol, from netlist tstamps."""
    out = {}
    for cb in F(nlroot, "components"):
        for comp in F(cb, "comp"):
            ref = sexp.val(comp, "ref")
            ts = sexp.val(comp, "tstamps")
            if ref and ts:
                out[ref[0]] = "/" + ts[0].strip("/")
    return out


def schematic_components(path=None):
    """Return ordered unique component records, skipping power and flag symbols."""
    with open(path or SCH, encoding="utf-8") as f:
        root = sexp.parse(f.read())[0]
    seen = {}
    order = []
    for c in root:
        if not (isinstance(c, list) and c and c[0] == "symbol"):
            continue
        ref = property_value(c, "Reference")
        val = property_value(c, "Value")
        fp = property_value(c, "Footprint")
        if not ref or ref.startswith("#"):
            continue
        if ref not in seen:
            dnp = sexp.val(c, "dnp", [sexp.Sym("no")])[0] == "yes"
            if not fp:
                libid = str(sexp.val(c, "lib_id", [""])[0])
                sys.exit(f"ERROR {ref} ({libid or 'no lib_id'}) has an empty "
                         "Footprint property; set one in Eeschema")
            seen[ref] = (ref, fp, val or "", dnp)
            order.append(ref)
    return [seen[r] for r in order]


@dataclass(frozen=True)
class GenerationContext:
    """Revision, placement mode, layout policy and output destination."""

    revision: str
    unplaced: bool
    layout: RevisionLayout
    out_dir: str

    @property
    def is_rev13(self):
        """Whether the fixed differential-sync placement is selected."""
        return self.revision == "1.3"


@dataclass(frozen=True)
class Netlist:
    """Component records, pad nets, native net ids and schematic paths."""

    components: dict
    pad_net: dict
    netid: dict
    paths: dict


@dataclass(frozen=True)
class Placement:
    """Placed origins, outline length and placement annotations."""

    positions: dict
    length: float
    bounded: list
    outfile: str
    note: str
    fixed: dict


@dataclass(frozen=True)
class Boxes:
    """Copper and courtyard bounds keyed by reference."""

    pads: dict
    courtyards: dict


HUB_CONNS = ("J1",)                # logic power in — hub end (left)


FAR_CONNS = ("J2", "J3A", "J3B")    # strip signal, sync daisy in/out — far end (right)
