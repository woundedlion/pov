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
