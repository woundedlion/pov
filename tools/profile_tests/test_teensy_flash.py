"""Board selection for locked firmware uploads."""

import os
import subprocess
import tempfile
import unittest
from pathlib import Path


FLASH = Path(__file__).resolve().parents[1] / "teensy_flash.sh"
LOCK_SH = FLASH.parent / "device_lock.sh"


class TeensyFlashTests(unittest.TestCase):
    def flash(self, port="COM4", token="owned", absent_calls=(), **lock_env):
        """Flash with a stub loader that omits COM4 on the 1-based absent_calls."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            calls = (root / "calls").as_posix()
            absent = " ".join(map(str, absent_calls))
            for name, body in {
                "teensy_ports.exe": (
                    f'n=$(( $(cat "{calls}" 2>/dev/null || echo 0) + 1 )); echo "$n" > "{calls}"\n'
                    "echo 'usb-a COM3 (Teensy 4.0)'\n"
                    f'case " {absent} " in *" $n "*) ;; *) echo "usb-b COM4 (Teensy 4.0)";; esac'),
                "teensy_post_compile.exe": "printf '%s\n' \"$@\"",
            }.items():
                script = root / name
                script.write_text("#!/bin/bash\n" + body + "\n", encoding="utf-8", newline="\n")
                script.chmod(0o755)
            env = {key: value for key, value in os.environ.items()
                   if not key.startswith(("HS_DEVICE_", "HS_TEENSY_", "HS_TYCMD", "HS_UPLOAD_"))}
            env["HS_UPLOAD_BACKEND"] = "vendor"
            env.update(HS_TEENSY_TOOLS=root.as_posix(), HS_TEENSY_PORT=port or "", **lock_env)
            if port is None:
                env.pop("HS_TEENSY_PORT", None)
            return subprocess.run(
                ["bash", "-c",
                 'set -u; . "$3"; . "$1"; _HS_TOKEN=$2; cygpath() { echo "$2"; }; '
                 's=0; sleep() { s=$((s + 1)); }; '
                 'hs_teensy_flash cache/entry; rc=$?; echo "SLEEPS=$s" >&2; exit "$rc"',
                 "test", FLASH.as_posix(), token, LOCK_SH.as_posix()],
                env=env, capture_output=True, text=True, encoding="utf-8",
            )

    def test_selects_claimed_board_usb_location(self):
        result = self.flash()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("-port=usb-b\n", result.stdout)
        self.assertIn("-portlabel=COM4 (Teensy 4.0)\n", result.stdout)
        self.assertNotIn("-port=usb-a", result.stdout)
        self.assertIn("-path=cache/entry\n", result.stdout)

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
        self.assertIn("no Teensy at COM9", result.stderr)
        self.assertNotIn("unbound variable", result.stderr)

    def test_missing_lock_fails_without_upload(self):
        result = self.flash(token="")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")
        self.assertIn("flash requires a device lock", result.stderr)
        self.assertNotIn("unbound variable", result.stderr)

    def test_waiting_flash_outlasts_a_board_briefly_off_the_list(self):
        result = self.flash(absent_calls=(1, 2), HS_DEVICE_WAIT="3600")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("-port=usb-b\n", result.stdout)
        self.assertIn("SLEEPS=2", result.stderr)

    def test_flash_without_a_wait_fails_fast_on_an_absent_board(self):
        result = self.flash(absent_calls=(1,))
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")
        self.assertIn("no Teensy at COM4", result.stderr)
        self.assertIn("SLEEPS=0", result.stderr)

    def test_waiting_flash_gives_up_after_the_absent_grace(self):
        result = self.flash(absent_calls=range(1, 100), HS_DEVICE_WAIT="3600",
                            HS_DEVICE_ABSENT_GRACE="15")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")
        self.assertIn("no Teensy at COM4", result.stderr)
        self.assertIn("SLEEPS=3", result.stderr)
