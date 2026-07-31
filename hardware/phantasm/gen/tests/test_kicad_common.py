import sys
import tempfile
import unittest
from pathlib import Path

GEN = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN))

import kicad_common  # noqa: E402


class RequireWritableTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.path = Path(directory.name) / "phantasm.kicad_pcb"

    def refusal(self, **kwargs):
        with self.assertRaises(SystemExit) as caught:
            kicad_common.require_writable(self.path, False, **kwargs)
        return str(caught.exception)

    def test_allows_an_absent_target(self):
        self.assertIsNone(kicad_common.require_writable(self.path, False))

    def test_allows_a_forced_overwrite(self):
        self.path.write_text("routed", encoding="utf-8")

        self.assertIsNone(kicad_common.require_writable(self.path, True))

    def test_refuses_an_existing_target(self):
        self.path.write_text("routed", encoding="utf-8")

        message = self.refusal()
        self.assertIn(f"refusing to overwrite {self.path}", message)
        self.assertIn("routing, vias, silk, hand edits", message)
        self.assertIn("--force", message)

    def test_states_the_caller_reason(self):
        self.path.write_text("routed", encoding="utf-8")

        message = self.refusal(reason="It is the fabrication source.")
        self.assertIn("It is the fabrication source.", message)
        self.assertNotIn("routing, vias, silk, hand edits", message)


if __name__ == "__main__":
    unittest.main()
