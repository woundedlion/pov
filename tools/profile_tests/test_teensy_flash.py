"""Board selection for locked firmware uploads."""

import os
import subprocess
import tempfile
import unittest
from pathlib import Path


FLASH = Path(__file__).resolve().parents[1] / "teensy_flash.sh"


class TeensyFlashTests(unittest.TestCase):
    def flash(self, port="COM4", token="owned"):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name, body in {
                "teensy_ports.exe": "printf '%s\\n' 'usb-a Serial COM3 Teensy' 'usb-b Serial COM4 Teensy'",
                "teensy_post_compile.exe": "printf '%s\\n' \"$@\"",
            }.items():
                script = root / name
                script.write_text("#!/bin/bash\n" + body + "\n", encoding="utf-8", newline="\n")
                script.chmod(0o755)
            env = dict(os.environ, HS_TEENSY_TOOLS=root.as_posix(), HS_TEENSY_PORT=port or "")
            if port is None:
                env.pop("HS_TEENSY_PORT", None)
            return subprocess.run(
                ["bash", "-c", 'set -u; . "$1"; _HS_TOKEN=$2; cygpath() { echo "$2"; }; hs_teensy_flash bench',
                 "test", FLASH.as_posix(), token],
                env=env, capture_output=True, text=True, encoding="utf-8",
            )

    def test_selects_claimed_board_usb_location(self):
        result = self.flash()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("-port=usb-b\n", result.stdout)
        self.assertNotIn("-port=usb-a", result.stdout)
        self.assertIn(".pio/build/bench", result.stdout)

    def test_portless_lock_reports_missing_pin(self):
        result = self.flash(port=None)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")
        self.assertIn("did not pin a Teensy port", result.stderr)
        self.assertNotIn("unbound variable", result.stderr)

    def test_missing_board_fails_without_upload(self):
        result = self.flash(port="COM9")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")

    def test_missing_lock_fails_without_upload(self):
        result = self.flash(token="")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")
