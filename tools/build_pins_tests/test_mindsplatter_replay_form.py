#!/usr/bin/env python3
"""Host tests for the MindSplatter generated-form gate."""

from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
CHECKER = ROOT / "tests/check_mindsplatter_replay_form.cmake"
FORM = """// Canonical form fixture.
inline const unsigned char HEAVY_SEARCH_V1_STATE[] PROGMEM = {0};
inline const unsigned short HEAVY_SEARCH_V1_FRAMEBUFFER[] PROGMEM = {0};
inline const Corpus HEAVY_SEARCH_V1 = {"msp-heavy-search-v1"};
inline const Corpus *const CORPUS_MANIFEST[] = {&HEAVY_SEARCH_V1};
""" + "// payload padding\n" * 6000


class MindSplatterReplayForm(unittest.TestCase):
    def setUp(self):
        self.cmake = shutil.which("cmake")
        self.assertIsNotNone(self.cmake, "cmake is required for generated-form tests")
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def run_checker(self, generated=FORM, compiler="Clang", version="23.0.0",
                    required=True):
        committed = self.root / "committed.h"
        committed.write_text(FORM, encoding="utf-8")
        output = self.root / "generated.py"
        output.write_text(
            "from pathlib import Path\n"
            f"Path(__file__).write_text({generated!r}, encoding='utf-8')\n",
            encoding="utf-8")
        return subprocess.run([
            self.cmake, f"-DGENERATOR={sys.executable}",
            f"-DCOMMITTED={committed}", f"-DGENERATED={output}",
            f"-DCOMPILER_ID={compiler}", f"-DCOMPILER_VERSION={version}",
            f"-DREQUIRE_COMPILER_MATCH={'ON' if required else 'OFF'}",
            "-DSKIP_CODE=77", "-P", str(CHECKER),
        ], capture_output=True, text=True, check=False)

    def test_clang_form_is_independent_of_compiler_major(self):
        for version in ("22.1.8", "23.0.0"):
            with self.subTest(version=version):
                result = self.run_checker(version=version)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertIn("msp-heavy-search-v1 in canonical form", result.stdout)

    def test_required_unsupported_compiler_is_rejected(self):
        result = self.run_checker(compiler="GNU", version="14.2.0")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("requires Clang, got GNU 14.2.0", result.stderr)

    def test_optional_unsupported_compiler_skips_without_running_generator(self):
        result = self.run_checker(compiler="GNU", required=False)
        self.assertEqual(result.returncode, 77, result.stdout + result.stderr)
        self.assertIn("from pathlib import Path", (self.root / "generated.py").read_text())

    def test_truncated_output_is_rejected(self):
        result = self.run_checker(generated="HEAVY_SEARCH_V1_STATE[]")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("generator emitted only", result.stderr)

    def test_missing_manifest_is_rejected(self):
        result = self.run_checker(generated=FORM.replace("CORPUS_MANIFEST", "OTHER_MANIFEST"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("generator omitted CORPUS_MANIFEST", result.stderr)

    def test_changed_preamble_is_rejected(self):
        result = self.run_checker(generated=FORM.replace("Canonical form", "Changed form"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("preamble differs from generated form", result.stderr)

    def test_changed_revision_is_rejected(self):
        result = self.run_checker(generated=FORM.replace("msp-heavy-search-v1", "msp-heavy-search-v2"))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("replay revision drift", result.stderr)


if __name__ == "__main__":
    unittest.main()
