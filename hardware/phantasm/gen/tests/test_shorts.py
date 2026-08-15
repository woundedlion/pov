"""Self-tests for the schematic short-detection gate.

Run:  python -m unittest discover -s hardware/phantasm/gen/tests
Every case is a synthetic .kicad_sch fragment, plus one regression case over the
committed schematic, so the gate is proven to fire AND proven not to cry wolf.
"""
import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path

GEN_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = GEN_DIR.parents[2]
sys.path.insert(0, str(GEN_DIR))

import sexp        # noqa: E402
import shorts      # noqa: E402


def _pin(number, x, y, angle):
    return (f'(pin passive line (at {x} {y} {angle}) (length 1.27) '
            f'(name "~" (effects (font (size 1.27 1.27)))) '
            f'(number "{number}" (effects (font (size 1.27 1.27)))))')


# Pin geometry copied from the stock KiCad libraries these tests stand in for.
LIB = {
    "power:GND": '(symbol "power:GND" (power global)'
                 f'  (symbol "GND_1_1" {_pin(1, 0, 0, 90)}))',
    "power:+5V": '(symbol "power:+5V" (power global)'
                 f'  (symbol "+5V_1_1" {_pin(1, 0, 0, 270)}))',
    # The flag's pin is a unit-0 (all-units) pin, as in the stock library.
    "power:PWR_FLAG": '(symbol "power:PWR_FLAG" (power global)'
                      f'  (symbol "PWR_FLAG_0_0" {_pin(1, 0, 0, 90)}))',
    "Device:R": '(symbol "Device:R"'
                f'  (symbol "R_1_1" {_pin(1, 0, 3.81, 270)} {_pin(2, 0, -3.81, 90)}))',
}


def sym(lib_id, ref, x, y, rot=0, mirror=None, unit=1, value=None):
    return dict(lib_id=lib_id, ref=ref, x=x, y=y, rot=rot, mirror=mirror,
                unit=unit, value=value)


def power(lib_id, x, y, rot=0):
    """A power/flag symbol: its Value is the net it names."""
    return sym(lib_id, "#PWR", x, y, rot=rot, value=lib_id.split(":")[1])


def build(symbols=(), wires=(), labels=(), global_labels=(),
          hierarchical_labels=(), junctions=(), sheets=(), buses=(),
          bus_entries=()):
    """Parse a minimal schematic holding exactly the given geometry."""
    out = ["(kicad_sch (version 20250114) (paper \"A3\")", "(lib_symbols"]
    for lib_id in sorted({s["lib_id"] for s in symbols}):
        out.append(LIB[lib_id])
    out.append(")")
    for (a, b) in wires:
        out.append(f"(wire (pts (xy {a[0]} {a[1]}) (xy {b[0]} {b[1]})))")
    for (p, text) in labels:
        out.append(f'(label "{text}" (at {p[0]} {p[1]} 0))')
    for (p, text) in global_labels:
        out.append(f'(global_label "{text}" (at {p[0]} {p[1]} 0))')
    for (p, text) in hierarchical_labels:
        out.append(f'(hierarchical_label "{text}" (at {p[0]} {p[1]} 0))')
    for p in junctions:
        out.append(f"(junction (at {p[0]} {p[1]}))")
    for name in sheets:
        out.append(f'(sheet (at 0 0) (property "Sheetname" "{name}" (at 0 0 0)))')
    for (a, b) in buses:
        out.append(f"(bus (pts (xy {a[0]} {a[1]}) (xy {b[0]} {b[1]})))")
    for p in bus_entries:
        out.append(f"(bus_entry (at {p[0]} {p[1]} 0) (size 2.54 2.54))")
    for s in symbols:
        out.append(f'(symbol (lib_id "{s["lib_id"]}") '
                   f'(at {s["x"]} {s["y"]} {s["rot"]})')
        if s["mirror"]:
            out.append(f'(mirror {s["mirror"]})')
        out.append(f'(unit {s["unit"]})')
        out.append(f'(property "Reference" "{s["ref"]}" (at 0 0 0))')
        out.append(f'(property "Value" "{s["value"] or s["ref"]}" (at 0 0 0))')
        out.append(")")
    out.append(")")
    return sexp.parse(" ".join(out))[0]


def conflicts(**kwargs):
    return shorts.analyze(build(**kwargs))[0]


class PowerFlagTests(unittest.TestCase):
    def test_flag_on_a_single_rail_passes(self):
        """A PWR_FLAG shares its rail's net by design -- not a short."""
        self.assertEqual(conflicts(
            symbols=[power("power:GND", 100, 100),
                     power("power:PWR_FLAG", 110, 100)],
            wires=[((100, 100), (110, 100))]), [])

    def test_flag_bridging_two_rails_fails(self):
        """A wire down the flag row merges +5V into GND: rail-to-rail short."""
        found = conflicts(
            symbols=[power("power:GND", 100, 100),
                     power("power:PWR_FLAG", 110, 100),
                     power("power:+5V", 120, 100)],
            wires=[((100, 100), (120, 100))])
        self.assertEqual([nets for nets, _ in found], [["+5V", "GND"]])

    def test_two_rails_without_a_flag_fails(self):
        found = conflicts(
            symbols=[power("power:GND", 100, 100), power("power:+5V", 120, 100)],
            wires=[((100, 100), (120, 100))])
        self.assertEqual([nets for nets, _ in found], [["+5V", "GND"]])

    def test_separate_rails_pass(self):
        self.assertEqual(conflicts(
            symbols=[power("power:GND", 100, 100), power("power:+5V", 100, 120)],
            wires=[((100, 100), (110, 100)), ((100, 120), (110, 120))]), [])

    def test_global_label_bridging_local_label_fails(self):
        found = conflicts(
            labels=[((100, 100), "GND")],
            global_labels=[((120, 100), "+5V")],
            wires=[((100, 100), (120, 100))])

        self.assertEqual([nets for nets, _ in found], [["+5V", "GND"]])

    def test_hierarchical_label_bridging_power_fails(self):
        found = conflicts(
            symbols=[power("power:GND", 100, 100)],
            hierarchical_labels=[((120, 100), "+5V")],
            wires=[((100, 100), (120, 100))])

        self.assertEqual([nets for nets, _ in found], [["+5V", "GND"]])

    def test_incomplete_geometry_is_ignored(self):
        root = sexp.parse(
            '(kicad_sch (label "GND") (global_label) '
            '(symbol (at 1 2)) (wire (pts (xy 1 2))) (junction (at 1)))'
        )[0]

        self.assertEqual(shorts.geometry(root), ({}, [], []))


class StackedNameTests(unittest.TestCase):
    """Two names at one coordinate are a short with no wire between them."""

    def test_stacked_labels_fail(self):
        found = conflicts(
            labels=[((100, 100), "GND"), ((100, 100), "+5V")],
            wires=[((100, 100), (110, 100))])

        self.assertEqual([nets for nets, _ in found], [["+5V", "GND"]])
        self.assertEqual(found[0][1],
                         [("+5V", (100.0, 100.0)), ("GND", (100.0, 100.0))])

    def test_label_stacked_on_a_power_symbol_fails(self):
        found = conflicts(
            symbols=[power("power:GND", 100, 100)],
            labels=[((100, 100), "+5V")],
            wires=[((100, 100), (110, 100))])

        self.assertEqual([nets for nets, _ in found], [["+5V", "GND"]])

    def test_stacked_label_pulls_its_whole_wire_group_in(self):
        """The second name at the point joins the group the point already had."""
        found = conflicts(
            symbols=[power("power:GND", 100, 100)],
            labels=[((120, 100), "+5V")],
            wires=[((100, 100), (120, 100))],
            global_labels=[((100, 100), "VBUS")])

        self.assertEqual([nets for nets, _ in found], [["+5V", "GND", "VBUS"]])

    def test_stacked_duplicate_of_the_same_net_passes(self):
        self.assertEqual(conflicts(
            labels=[((100, 100), "GND"), ((100, 100), "GND")],
            wires=[((100, 100), (110, 100))]), [])

    def test_flag_stacked_on_a_rail_passes(self):
        self.assertEqual(conflicts(
            symbols=[power("power:GND", 100, 100),
                     power("power:PWR_FLAG", 100, 100)],
            wires=[((100, 100), (110, 100))]), [])

    def test_named_points_keeps_every_name(self):
        named = shorts.named_points(build(
            labels=[((100, 100), "GND"), ((100, 100), "+5V")],
            wires=[((100, 100), (110, 100))]))

        self.assertEqual(named, {(100.0, 100.0): {"GND", "+5V"}})


class EmptyScanTests(unittest.TestCase):
    """An empty scan is a broken input, never a clean bill of health."""

    def test_empty_schematic_fails(self):
        with self.assertRaises(ValueError):
            shorts.analyze(sexp.parse("(kicad_sch)")[0])

    def test_geometry_without_names_fails(self):
        with self.assertRaises(ValueError):
            shorts.analyze(build(wires=[((100, 100), (110, 100))]))

    def test_names_without_wires_fail(self):
        with self.assertRaises(ValueError):
            shorts.analyze(build(labels=[((100, 100), "GND")]))

    def test_main_exits_nonzero_on_an_empty_schematic(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "empty.kicad_sch"
            path.write_text("(kicad_sch)", encoding="utf-8")
            with contextlib.redirect_stderr(io.StringIO()):
                self.assertNotEqual(shorts.main(["shorts.py", str(path)]), 0)


class UnscannedConnectivityTests(unittest.TestCase):
    """The scan reads top-level wires only, so connectivity it cannot follow
    must hard-fail rather than read clean."""

    def clean_top_level(self, **extra):
        return build(symbols=[power("power:GND", 100, 100)],
                     wires=[((100, 100), (110, 100))], **extra)

    def test_sheet_fails_even_when_the_top_level_is_clean(self):
        with self.assertRaisesRegex(ValueError, "sheet present"):
            shorts.analyze(self.clean_top_level(sheets=["power"]))

    def test_bus_fails_even_when_the_top_level_is_clean(self):
        with self.assertRaisesRegex(ValueError, "bus present"):
            shorts.analyze(self.clean_top_level(buses=[((120, 100), (140, 100))]))

    def test_bus_entry_fails_even_when_the_top_level_is_clean(self):
        with self.assertRaisesRegex(ValueError, "bus_entry present"):
            shorts.analyze(self.clean_top_level(bus_entries=[(120, 100)]))


class CommittedSchematicTests(unittest.TestCase):
    def setUp(self):
        path = REPO_ROOT / "hardware" / "phantasm" / "phantasm.kicad_sch"
        self.root = sexp.parse(path.read_text(encoding="utf-8"))[0]

    def test_committed_schematic_scan_is_non_trivial(self):
        """Pins the geometry a passing run must have seen, so a scan that finds
        nothing cannot masquerade as a clean schematic."""
        named, wires, junctions = shorts.geometry(self.root)

        self.assertEqual(len(named), 73)
        self.assertEqual(len(wires), 81)
        self.assertEqual(len(junctions), 9)

    def test_committed_schematic_has_no_conflicts(self):
        conflicts, bridges = shorts.analyze(self.root)

        self.assertEqual(conflicts, [])
        self.assertEqual(len(bridges), 6)


if __name__ == "__main__":
    unittest.main()
