import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

GEN = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(GEN))

import kicad_common  # noqa: E402


class FindKicadCliTests(unittest.TestCase):
    """Every tool that shells out to KiCad resolves the binary through here."""

    def test_env_override_wins_when_it_exists(self):
        with tempfile.TemporaryDirectory() as directory:
            cli = Path(directory) / "kicad-cli"
            cli.touch()
            with mock.patch.dict(os.environ, {"KICAD_CLI": str(cli)}):
                self.assertEqual(kicad_common.find_kicad_cli(), str(cli))

    def test_falls_back_to_the_path_name(self):
        with mock.patch.dict(os.environ, {"KICAD_CLI": "missing-kicad-cli"}), \
                mock.patch.object(kicad_common.glob, "glob", return_value=[]):
            self.assertEqual(kicad_common.find_kicad_cli(), "kicad-cli")

    def test_prefers_the_newest_install(self):
        installs = [r"C:\Program Files\KiCad\9.0\bin\kicad-cli.exe",
                    r"C:\Program Files\KiCad\10.0\bin\kicad-cli.exe"]
        with mock.patch.dict(os.environ, {}, clear=True), \
                mock.patch.object(kicad_common.glob, "glob",
                                  side_effect=lambda p: installs if "Program" in p
                                  and "x86" not in p else []):
            self.assertEqual(kicad_common.find_kicad_cli(), installs[1])


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
