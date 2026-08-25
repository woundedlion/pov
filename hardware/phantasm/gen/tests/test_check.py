import contextlib
import io
import sys
import unittest
from pathlib import Path

GEN = Path(__file__).resolve().parent.parent
REPO_ROOT = GEN.parents[2]
sys.path.insert(0, str(GEN))

import check  # noqa: E402
import kicad_common  # noqa: E402
import sexp  # noqa: E402


def expected_nodes():
    """EXPECT -> {net: [(ref, pin), ...]}, numbering symmetric terminals 1, 2."""
    seq, nets = {}, {}
    for name, keys in sorted(check.EXPECT.items()):
        nodes = []
        for key in sorted(keys):
            if "." in key:
                ref, pin = key.split(".", 1)
            else:
                ref = key
                seq[ref] = seq.get(ref, 0) + 1
                pin = str(seq[ref])
            nodes.append((ref, pin))
        nets[name] = nodes
    return nets


def netlist(nets):
    """Render {net: [(ref, pin), ...]} as a kicadsexpr netlist root."""
    out = ['(export (version "E") (nets']
    for code, (name, nodes) in enumerate(sorted(nets.items()), 1):
        out.append(f'(net (code "{code}") (name "/{name}")')
        out += [f'(node (ref "{ref}") (pin "{pin}"))' for ref, pin in nodes]
        out.append(")")
    out.append("))")
    return sexp.parse("\n".join(out))[0]


def committed_board_nets():
    path = REPO_ROOT / "hardware" / "phantasm" / "phantasm.kicad_pcb"
    root = sexp.parse(path.read_text(encoding="utf-8"))[0]
    nets = {}
    for footprint in (
            node for node in root if isinstance(node, list) and node
            and node[0] == "footprint"):
        reference = next(
            (str(node[2]) for node in footprint
             if isinstance(node, list) and len(node) > 2
             and node[0] == "property" and node[1] == "Reference"),
            None,
        )
        for pad in (
                node for node in footprint if isinstance(node, list) and node
                and node[0] == "pad"):
            name = kicad_common.net_name(pad)
            if reference and name:
                nets.setdefault(name, set()).add(check.node_key(reference, str(pad[1])))
    return nets


def run(nets):
    """Run the gate over a synthetic netlist; return (ok, printed output)."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        ok = check.check(check.netlist_nets(netlist(nets)))
    return ok, buf.getvalue()


def swap_pins(nets, ref, a, b):
    return {
        name: [(r, b if p == a else a if p == b else p) if r == ref else (r, p)
               for r, p in nodes]
        for name, nodes in nets.items()
    }


class NodeKeyTests(unittest.TestCase):
    def test_keys_fixed_pinouts_on_ref_and_pin(self):
        self.assertEqual(check.node_key("J2", "1"), "J2.1")
        self.assertEqual(check.node_key("U_MCU", "VIN"), "U_MCU.VIN")

    def test_keys_interchangeable_terminals_on_ref_alone(self):
        self.assertEqual(check.node_key("R_D1", "2"), "R_D1")


class ExpectTableTests(unittest.TestCase):
    def test_every_fixed_pinout_node_carries_a_pin(self):
        for name, keys in check.EXPECT.items():
            for key in keys:
                if key.split(".")[0] not in check.SYMMETRIC:
                    self.assertIn(".", key, f"{name}: {key} has no pin")

    def test_every_symmetric_ref_spans_two_nets(self):
        for ref in check.SYMMETRIC:
            hits = [n for n, keys in check.EXPECT.items() if ref in keys]
            self.assertEqual(len(hits), 2, f"{ref} on nets {sorted(hits)}")


class MalformedNetlistTests(unittest.TestCase):
    MALFORMED = sexp.parse(
        '(export (nets'
        ' (net (code "1"))'
        ' (net (code "2") (name "/DATA") (node (ref "J2")))'
        ' (net (code "3") (name "/CLK") (node (ref "J2") (pin "2")))'
        "))")[0]

    def test_skips_entries_a_malformed_export_left_incomplete(self):
        self.assertEqual(check.netlist_nets(self.MALFORMED),
                         {"DATA": set(), "CLK": {"J2.2"}})


class GateTests(unittest.TestCase):
    def test_accepts_committed_board(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            ok = check.check(committed_board_nets())
        self.assertTrue(ok, output.getvalue())
        self.assertEqual(output.getvalue(), "")

    def test_accepts_spec_netlist(self):
        ok, out = run(expected_nodes())
        self.assertTrue(ok, out)
        self.assertEqual(out, "")

    def test_rejects_swapped_strip_data_and_clock(self):
        ok, out = run(swap_pins(expected_nodes(), "J2", "1", "2"))
        self.assertFalse(ok)
        self.assertIn("FAIL CLK\n   missing ['J2.2']", out)
        self.assertIn("FAIL DATA\n   missing ['J2.1']", out)

    def test_rejects_swapped_sync_and_ground(self):
        ok, out = run(swap_pins(expected_nodes(), "J3A", "1", "2"))
        self.assertFalse(ok)
        self.assertIn("FAIL GND\n   missing ['J3A.2']", out)
        self.assertIn("FAIL SYNC_BUS\n   missing ['J3A.1']", out)

    def test_rejects_swapped_debug_header_rails(self):
        ok, out = run(swap_pins(expected_nodes(), "J4", "1", "2"))
        self.assertFalse(ok)
        self.assertIn("FAIL +3V3", out)
        self.assertIn("FAIL GND", out)

    def test_rejects_bridged_series_resistor(self):
        nets = expected_nodes()
        nets["DATA"] = nets["DATA"] + nets["DATA_SRC"]
        nets["DATA_SRC"] = [n for n in nets["DATA_SRC"] if n[0] != "R_D1"]
        ok, out = run(nets)
        self.assertFalse(ok)
        self.assertIn("FAIL DATA_SRC\n   missing ['R_D1']", out)

    def test_notes_named_net_outside_the_spec_table(self):
        nets = expected_nodes()
        nets["STRAY"] = [("J2", "9")]
        ok, out = run(nets)
        self.assertTrue(ok, out)
        self.assertIn("NOTE extra named net STRAY", out)

    def test_notes_a_power_named_net_outside_the_spec_table(self):
        nets = expected_nodes()
        nets["PWR_EN"] = [("J4", "9")]
        ok, out = run(nets)
        self.assertTrue(ok, out)
        self.assertIn("NOTE extra named net PWR_EN", out)

    def test_skips_kicad_auto_generated_net_names(self):
        nets = expected_nodes()
        nets["unconnected-(U1-Pad2)"] = [("U1", "2")]
        nets["Net-(R_S-Pad1)"] = [("R_S", "1")]
        ok, out = run(nets)
        self.assertTrue(ok, out)
        self.assertNotIn("NOTE", out)


if __name__ == "__main__":
    unittest.main()
