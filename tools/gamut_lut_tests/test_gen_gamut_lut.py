import contextlib
import io
import re
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
sys.path.insert(0, str(TOOLS))

import gen_gamut_lut as generator  # noqa: E402


class TestGamutLutMirrors(unittest.TestCase):
    def test_provenance_rejects_numeric_header_and_newline_drift(self):
        table = generator.np.zeros((1, 1, 2), dtype=generator.np.uint16)
        expected = generator.render(table)
        variants = (
            expected.replace("0, 0,", "0, 1,", 1),
            expected.replace("#pragma once", "#pragma twice", 1),
            expected.replace("\n", "\r\n"),
        )
        with tempfile.TemporaryDirectory() as directory, \
                mock.patch.object(generator, "build_table", return_value=(table, 0)):
            path = Path(directory) / "gamut_lut.h"
            path.write_bytes(expected.encode("utf-8"))
            self.assertTrue(generator.check_provenance(str(path)))
            for changed in variants:
                with self.subTest(changed=changed[:40]):
                    self.assertNotEqual(changed, expected)
                    path.write_bytes(changed.encode("utf-8"))
                    error = io.StringIO()
                    with contextlib.redirect_stderr(error):
                        self.assertFalse(generator.check_provenance(str(path)))
                    self.assertIn("out of sync", error.getvalue())
                    self.assertIn("Regenerate with:", error.getvalue())

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

    def test_refinement_slack_and_matrix_contract_mutations_are_rejected(self):
        color_path = ROOT / "core/color/color_space.h"
        math_path = ROOT / "core/math/3dmath.h"
        source = color_path.read_text(encoding="utf-8")
        mutations = (
            (r"(GAMUT_SCAN_STEPS\s*=\s*)4", r"\g<1>5", "GAMUT_SCAN_STEPS differs"),
            (r"(GAMUT_BRACKET_STEPS\s*=\s*)3", r"\g<1>4", "GAMUT_BRACKET_STEPS differs"),
            (r"(GATE_LO\s*=\s*)-1e-4f", r"\g<1>-2e-4f", "slack"),
            (r"inline constexpr float OKLAB_TO_LMS_CBRT\[3\]\[3\]", "inline constexpr float RENAMED[3][3]", "missing shared matrix"),
            (r"0.3963377774f,", "", "has 8 coefficients"),
        )
        self.assertTrue(generator.check_mirrors(color_path, math_path))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "color_space.h"
            for pattern, replacement, diagnostic in mutations:
                with self.subTest(diagnostic=diagnostic):
                    changed, count = re.subn(pattern, replacement, source, count=1)
                    self.assertEqual(count, 1)
                    path.write_text(changed, encoding="utf-8")
                    with contextlib.redirect_stderr(io.StringIO()) as errors:
                        self.assertFalse(generator.check_mirrors(path, math_path))
                    self.assertIn(diagnostic, errors.getvalue())

    def test_parseable_angle_semantic_drift_is_rejected(self):
        original = ROOT / "core/math/3dmath.h"
        source = original.read_text(encoding="utf-8")
        self.assertTrue(generator.check_diamond_angle_mirror(original))
        changed, count = re.subn(r"(inline float diamond_angle\([^)]*\)\s*\{)",
                                r"\g<1> return 0.0f;", source, count=1)
        self.assertEqual(count, 1)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "3dmath.h"
            path.write_text(changed, encoding="utf-8")
            with contextlib.redirect_stderr(io.StringIO()) as errors:
                self.assertFalse(generator.check_diamond_angle_mirror(path))
            self.assertIn("mirror drift", errors.getvalue())

    def test_angle_roundtrip_refuses_a_changed_direction(self):
        self.assertTrue(generator.check_angle_roundtrip())
        real = generator.diamond_direction

        def changed(angle):
            x, y = real(angle)
            return y, x

        with mock.patch.object(generator, "diamond_direction", side_effect=changed), \
             contextlib.redirect_stderr(io.StringIO()) as errors:
            self.assertFalse(generator.check_angle_roundtrip())
        self.assertIn("no longer inverts", errors.getvalue())

    def test_operand_swaps_fail_the_oklab_mirror_check(self):
        color_h = (ROOT / "core" / "color" / "color_space.h").read_text(
            encoding="utf-8")
        mutations = (
            ("0.3963377774f, 0.2158037573f",
             "0.2158037573f, 0.3963377774f"),
            ("4.0767416621f, -3.3077115913f",
             "-3.3077115913f, 4.0767416621f"),
        )
        math_h = ROOT / "core" / "math" / "3dmath.h"
        self.assertTrue(generator.check_mirrors(
            ROOT / "core" / "color" / "color_space.h", math_h))
        for original, replacement in mutations:
            with self.subTest(original=original):
                pattern = r"\s+".join(re.escape(token) for token in original.split())
                mutated, count = re.subn(pattern, replacement, color_h, count=1)
                self.assertEqual(count, 1)
                self.assertNotEqual(mutated, color_h)
                with tempfile.TemporaryDirectory() as directory:
                    path = Path(directory) / "color_space.h"
                    path.write_text(mutated, encoding="utf-8")
                    error = io.StringIO()
                    with contextlib.redirect_stderr(error):
                        result = generator.check_mirrors(path, math_h)
                self.assertFalse(result)
                self.assertIn("mirror drift in", error.getvalue())


if __name__ == "__main__":
    unittest.main()
