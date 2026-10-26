"""Regression tests for the tracked-text encoding gate."""

import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


TOOL = Path(__file__).resolve().parents[1] / "encoding_check.py"
SPEC = importlib.util.spec_from_file_location("encoding_check", TOOL)
sys.path.insert(0, str(TOOL.parent))
encoding_check = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(encoding_check)


class EncodingCheckTests(unittest.TestCase):
    def test_reversible_runs_are_rejected(self):
        for text in ("\u2014", "\u2192", "\u00b1", "\u00b0", "\u00d7"):
            with self.subTest(text=text):
                corrupt = text.encode("utf-8").decode("cp1252")
                self.assertEqual(encoding_check.mojibake_lines("clean\n" + corrupt), [2])

    def test_correct_unicode_is_accepted(self):
        self.assertEqual(encoding_check.mojibake_lines("ASCII \u2014 \u2192 \u00b1 \u00b0 \u00d7 caf\u00e9"), [])

    def test_gate_reads_only_tracked_text(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            subprocess.run(["git", "-C", str(root), "init", "-q"], check=True)
            corrupt = "\u2014".encode("utf-8").decode("cp1252")
            (root / "tracked.txt").write_text("\u2014\n", encoding="utf-8")
            (root / "untracked.txt").write_text(corrupt, encoding="utf-8")
            (root / "binary.bin").write_bytes(b"\0" + corrupt.encode("utf-8"))
            subprocess.run(["git", "-C", str(root), "add", "--", "tracked.txt", "binary.bin"], check=True)
            result = subprocess.run([sys.executable, str(TOOL), "--root", str(root)], capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            (root / "tracked.txt").write_text("clean\n" + corrupt, encoding="utf-8")
            result = subprocess.run([sys.executable, str(TOOL), "--root", str(root)], capture_output=True)
            self.assertEqual(result.returncode, 1, result.stderr)
            self.assertIn(b"tracked.txt:2:", result.stderr)
            self.assertNotIn(b"untracked.txt", result.stderr)


if __name__ == "__main__":
    unittest.main()
