#!/usr/bin/env python3
"""Host tests for process-environment omission from build captures."""

import subprocess
import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS))

from profile_envdump import sanitize_envdump  # noqa: E402


class EnvdumpTests(unittest.TestCase):
    def test_preserves_build_fields_and_line_endings(self):
        prefix = "Processing profile\r\n{ 'CCFLAGS': ['-O3'],\r\n  'ENV': "
        suffix = ",\r\n  'CPPDEFINES': [('HS_PROFILE_TARGET', 'Fx')]}\r\nSUCCESS\r\n"
        environment = repr({"PATH": "C:\\tools", "TOKEN": "braces } { and 'quotes' \""})
        for value in (environment, f"environ({environment})", "{}"):
            with self.subTest(value=value):
                self.assertEqual(sanitize_envdump(prefix + value + suffix),
                                 prefix + "{}" + suffix)

    def test_multiline_dictionary(self):
        text = "{ 'CC': <object at 0x123>,\n  'ENV': {'PATH': 'bin',\n          'USER': 'user'},\n  'CXX': 'g++'}"
        self.assertEqual(sanitize_envdump(text),
                         "{ 'CC': <object at 0x123>,\n  'ENV': {},\n  'CXX': 'g++'}")

    def test_rejects_missing_duplicate_and_malformed_fields(self):
        for text in ("", "build failed", "  'ENV': {},\n  'ENV': {}", "  'ENV': []}",
                     "  'ENV': {'PATH': 'unterminated}", "  'ENV': environ({}},",
                     "  'ENV': {'PATH': 5},", "  'ENV': {} trailing",
                     "  'ENV': {'PATH': run()},", "  'ENV': {'PATH': 'x' 'USER': 'y'},"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                sanitize_envdump(text)

    def test_cli_preserves_non_utf8_bytes_outside_environment(self):
        result = subprocess.run([sys.executable, str(TOOLS / "profile_envdump.py")],
                                input=b"header\xff\r\n  'ENV': {'TOKEN': 'secret'},\r\n",
                                capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, b"header\xff\r\n  'ENV': {},\r\n")

    def test_cli_failure_emits_no_capture_or_environment(self):
        result = subprocess.run([sys.executable, str(TOOLS / "profile_envdump.py")],
                                input=b"  'ENV': {'TOKEN': 'secret", capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, b"")
        self.assertNotIn(b"secret", result.stderr)


if __name__ == "__main__":
    unittest.main()
