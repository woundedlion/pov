"""Shell dispatch preserves lock ownership and explicit backend selection."""

import os
from pathlib import Path
import subprocess
import tempfile
import unittest


TOOLS = Path(__file__).resolve().parents[1]


class UploadDispatchTests(unittest.TestCase):
    def dispatch(self, backend="auto", available=True, explicit=False, upload_result=0,
                 token="owned", lock_dir="fixture-claim", resolve_result=None):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name, content in {
                "teensy_ports.exe": "echo 'usb-b COM4 (Teensy 4.0)'",
                "teensy_post_compile.exe": "echo VENDOR",
            }.items():
                path = root / name
                path.write_text("#!/bin/bash\n" + content + "\n", encoding="utf-8", newline="\n")
                path.chmod(0o755)
            environment = {key: value for key, value in os.environ.items() if not key.startswith("HS_")}
            environment.update(HS_TEENSY_TOOLS=root.as_posix(), HS_TEENSY_PORT="COM4",
                               HS_UPLOAD_BACKEND=backend,
                               FIXTURE_AVAILABLE=str(resolve_result if resolve_result is not None else 0 if available else 2),
                               FIXTURE_UPLOAD_RESULT=str(upload_result), FIXTURE_CLAIM=lock_dir)
            if explicit:
                environment["HS_TYCMD"] = "fixture-tycmd.exe"
            return subprocess.run([
                "bash", "-c", '. "$1/device_lock.sh"; . "$1/teensy_flash.sh"; '
                'fake_python() { if [ "$2" = --resolve ]; then return "$FIXTURE_AVAILABLE"; fi; '
                'printf "%s\\n" "$@"; return "$FIXTURE_UPLOAD_RESULT"; }; '
                '_hs_resolve_python() { _HS_LOCK_PYTHON=fake_python; }; '
                'hs_native_path() { printf "%s\\n" "$1"; }; '
                'cygpath() { printf "%s\\n" "$2"; }; '
                '_HS_TOKEN=$2; _HS_LOCK_DIR=$FIXTURE_CLAIM; hs_teensy_flash fixture-image',
                "test", TOOLS.as_posix(), token,
            ], env=environment, capture_output=True, text=True, timeout=10)

    def test_available_pinned_tool_is_the_automatic_default(self):
        result = self.dispatch()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("teensy_upload.py\n--port\nCOM4\n", result.stdout)
        self.assertIn("--lock-dir\nfixture-claim\n--lock-token\nowned\n", result.stdout)
        self.assertIn("fixture-image/firmware.hex\n", result.stdout)
        self.assertNotIn("VENDOR", result.stdout)

    def test_missing_headless_tool_keeps_vendor_default(self):
        result = self.dispatch(available=False)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), "VENDOR")

    def test_explicit_path_selects_headless_without_fallback(self):
        result = self.dispatch(available=False, explicit=True, upload_result=7)
        self.assertEqual(result.returncode, 7)
        self.assertNotIn("VENDOR", result.stdout)

    def test_explicit_vendor_wins_over_headless_override(self):
        result = self.dispatch(backend="vendor", explicit=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), "VENDOR")

    def test_invalid_installed_headless_tool_does_not_fall_back_to_vendor(self):
        result = self.dispatch(resolve_result=1)
        self.assertEqual(result.returncode, 1)
        self.assertNotIn("VENDOR", result.stdout)

    def test_missing_token_or_claim_prevents_headless_process(self):
        for options, diagnostic in (({"token": ""}, "requires a device lock"),
                                    ({"lock_dir": ""}, "lock directory is missing")):
            with self.subTest(options=options):
                result = self.dispatch(**options)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(result.stdout, "")
                self.assertIn(diagnostic, result.stderr)

    def test_unknown_backend_is_rejected_without_any_upload(self):
        result = self.dispatch(backend="unknown")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")
        self.assertIn("unknown HS_UPLOAD_BACKEND", result.stderr)


if __name__ == "__main__":
    unittest.main()
