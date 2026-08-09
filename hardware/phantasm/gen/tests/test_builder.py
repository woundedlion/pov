import sys
import unittest
from pathlib import Path

GEN = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN))

import builder  # noqa: E402
import sexp  # noqa: E402

SCHEMATIC = GEN.parent / "phantasm.kicad_sch"

ORIGIN = (10.0, 20.0)
LOCAL = (2.0, 1.0)  # library frame, y up

# Expected globals derived by hand: library (2, 1) flips to screen (2, -1);
# mirror acts in the screen frame before rotation ("x" vertical, "y"
# horizontal); rotation is counter-clockwise as displayed, so 90 deg maps
# screen (x, y) -> (y, -x).
TRANSFORM_CASES = [
    (0, None, (12.0, 19.0)),
    (90, None, (9.0, 18.0)),
    (180, None, (8.0, 21.0)),
    (270, None, (11.0, 22.0)),
    (0, "x", (12.0, 21.0)),
    (90, "x", (11.0, 18.0)),
    (180, "x", (8.0, 19.0)),
    (270, "x", (9.0, 22.0)),
    (0, "y", (8.0, 19.0)),
    (90, "y", (9.0, 22.0)),
    (180, "y", (12.0, 21.0)),
    (270, "y", (11.0, 18.0)),
]

# Screen stub direction for an unplaced pin: the pin angle points tip->body in
# the y-up library frame, so outward is angle+180, y-flipped onto the screen.
OUTWARD = {0: (-1, 0), 90: (0, 1), 180: (1, 0), 270: (0, -1)}
MIRROR = {
    None: lambda v: v,
    "x": lambda v: (v[0], -v[1]),
    "y": lambda v: (-v[0], v[1]),
}
ROTATE = {
    0: lambda v: v,
    90: lambda v: (v[1], -v[0]),
    180: lambda v: (-v[0], -v[1]),
    270: lambda v: (-v[1], v[0]),
}

PIN_DIR_HAND_CHECKED = [
    # (pin angle, rot, mirror, expected stub direction)
    (0, 0, None, (-1, 0)),   # pin left of body: stub goes left
    (90, 0, None, (0, 1)),   # pin below body: stub goes down (screen y is down)
    (0, 90, None, (0, 1)),   # left stub turned a quarter CCW: down
    (0, 0, "y", (1, 0)),     # horizontal mirror sends the left stub right
    (270, 0, "x", (0, 1)),   # vertical mirror sends the top stub down
]


def stub_symbol(angle, rot, mirror):
    sym = builder.Symbol("Test:Stub", "U1", "stub", 0, 0, rot=rot, mirror=mirror)
    sym._pins = {"1": {"x": 0.0, "y": 0.0, "angle": float(angle), "name": "1"}}
    return sym


def grid_key(point):
    return (round(point[0], 3), round(point[1], 3))


class TransformMatrixTests(unittest.TestCase):
    def test_rotation_mirror_matrix(self):
        for rot, mirror, expected in TRANSFORM_CASES:
            with self.subTest(rot=rot, mirror=mirror):
                got = builder.transform(*ORIGIN, rot, mirror, *LOCAL)
                self.assertEqual(got, expected)

    def test_symbol_origin_is_a_fixed_point(self):
        for rot, mirror, _ in TRANSFORM_CASES:
            with self.subTest(rot=rot, mirror=mirror):
                got = builder.transform(*ORIGIN, rot, mirror, 0.0, 0.0)
                self.assertEqual(got, ORIGIN)


class PinDirMatrixTests(unittest.TestCase):
    def test_hand_checked_cases(self):
        for angle, rot, mirror, expected in PIN_DIR_HAND_CHECKED:
            with self.subTest(angle=angle, rot=rot, mirror=mirror):
                self.assertEqual(stub_symbol(angle, rot, mirror).pin_dir("1"),
                                 expected)

    def test_full_rotation_mirror_matrix(self):
        for angle, outward in OUTWARD.items():
            for mirror, flip in MIRROR.items():
                for rot, turn in ROTATE.items():
                    with self.subTest(angle=angle, rot=rot, mirror=mirror):
                        expected = turn(flip(outward))
                        got = stub_symbol(angle, rot, mirror).pin_dir("1")
                        self.assertEqual(got, expected)


class SchematicPinnedTests(unittest.TestCase):
    """Pin transform against the committed, netlist-verified schematic."""

    @classmethod
    def setUpClass(cls):
        root = sexp.parse(SCHEMATIC.read_text(encoding="utf-8"))[0]
        cls.placed = []
        cls.anchors = set()
        unit_pins = {}
        for node in root:
            if not (isinstance(node, list) and node):
                continue
            if node[0] == "lib_symbols":
                for sym in node:
                    if isinstance(sym, list) and sym and sym[0] == "symbol":
                        unit_pins[sym[1]] = builder._index_unit_pins(sym)
            elif node[0] == "symbol":
                at = sexp.val(node, "at")
                mirror = sexp.val(node, "mirror")
                unit = sexp.val(node, "unit")
                cls.placed.append(dict(
                    lib_id=sexp.val(node, "lib_id")[0],
                    x=float(at[0]), y=float(at[1]),
                    rot=float(at[2]) if len(at) > 2 else 0.0,
                    mirror=mirror[0] if mirror else None,
                    unit=int(unit[0]) if unit else 1))
            elif node[0] == "wire":
                pts = [c for c in node
                       if isinstance(c, list) and c and c[0] == "pts"][0]
                for xy in pts:
                    if isinstance(xy, list) and xy and xy[0] == "xy":
                        cls.anchors.add(grid_key((float(xy[1]), float(xy[2]))))
            elif node[0] in ("junction", "label"):
                at = sexp.val(node, "at")
                cls.anchors.add(grid_key((float(at[0]), float(at[1]))))
        cls.pins = []
        for sym in cls.placed:
            pins = dict(unit_pins[sym["lib_id"]].get(0, {}))
            pins.update(unit_pins[sym["lib_id"]].get(sym["unit"], {}))
            for pin in pins.values():
                where = builder.transform(sym["x"], sym["y"], sym["rot"],
                                          sym["mirror"], pin["x"], pin["y"])
                cls.pins.append((sym["rot"], where))

    def test_exercises_every_rotation(self):
        self.assertEqual({rot for rot, _ in self.pins}, {0.0, 90.0, 180.0, 270.0})
        self.assertGreaterEqual(len(self.pins), 100)

    def test_every_pin_lands_on_the_drawn_net(self):
        touching = {}
        for _, where in self.pins:
            touching[grid_key(where)] = touching.get(grid_key(where), 0) + 1
        floating = [(rot, where) for rot, where in self.pins
                    if grid_key(where) not in self.anchors
                    and touching[grid_key(where)] < 2]
        self.assertEqual(floating, [])


if __name__ == "__main__":
    unittest.main()
