import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
sys.path.insert(0, str(TOOLS))

import gen_gamut_lut as generator  # noqa: E402


class TestGamutLutMirrors(unittest.TestCase):
    def test_cpp_logical_operators_are_evaluated(self):
        parsed = generator._cpp_float_fn(
            "{ if (x > 0.0f && x != 4.0f && "
            "!(x == 2.0f || x == 3.0f)) return 1.0f; "
            "return 0.0f;", ("x",))
        self.assertEqual(parsed(1.0), 1.0)
        self.assertEqual(parsed(2.0), 0.0)
        self.assertEqual(parsed(4.0), 0.0)
        self.assertEqual(parsed(-1.0), 0.0)

    def test_hex_literal_is_not_treated_as_float_suffixed(self):
        parsed = generator._cpp_float_fn("{ return 0x1f;", ())
        self.assertEqual(parsed(), 31)

    def test_malformed_cpp_expression_reports_mirror_guidance(self):
        source = ("inline float diamond_angle(float y, float x) {\n"
                  "  return x &&;\n"
                  "}\n")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "3dmath.h"
            path.write_text(source, encoding="utf-8")
            error = io.StringIO()
            with contextlib.redirect_stderr(error):
                result = generator.check_diamond_angle_mirror(path)
        self.assertFalse(result)
        self.assertIn("could not parse 3dmath.h diamond_angle()", error.getvalue())
        self.assertIn("re-derive the mirror", error.getvalue())

    def test_operand_swaps_fail_the_oklab_mirror_check(self):
        color_h = (ROOT / "core" / "color" / "color.h").read_text(
            encoding="utf-8")
        mutations = (
            ("0.3963377774f * lab.a + 0.2158037573f * lab.b",
             "0.3963377774f * lab.b + 0.2158037573f * lab.a"),
            ("4.0767416621f * l - 3.3077115913f * m",
             "4.0767416621f * m - 3.3077115913f * l"),
        )
        math_h = ROOT / "core" / "math" / "3dmath.h"
        for original, replacement in mutations:
            with self.subTest(original=original):
                mutated = color_h.replace(original, replacement)
                self.assertNotEqual(mutated, color_h)
                with tempfile.TemporaryDirectory() as directory:
                    path = Path(directory) / "color.h"
                    path.write_text(mutated, encoding="utf-8")
                    with contextlib.redirect_stderr(io.StringIO()):
                        result = generator.check_mirrors(path, math_h)
                self.assertFalse(result)


if __name__ == "__main__":
    unittest.main()
