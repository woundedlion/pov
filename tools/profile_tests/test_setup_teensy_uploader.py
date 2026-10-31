"""Pinned release installation never extracts unverified uploader bytes."""

import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import setup_teensy_uploader as setup  # noqa: E402


def archive(executable):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as package:
        package.writestr("tycmd.exe", executable)
        package.writestr("LICENSE.txt", "fixture license")
        package.writestr("README.md", "fixture documentation")
    return buffer.getvalue()


class UploaderSetupTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.destination = Path(temporary.name) / "installed"
        self.executable = b"independent executable fixture"
        self.release = archive(self.executable)
        self.addCleanup(mock.patch.stopall)
        mock.patch.object(setup, "TYCMD_SHA256", hashlib.sha256(self.executable).hexdigest()).start()
        mock.patch.object(setup, "TYTOOLS_ARCHIVE_SHA256", hashlib.sha256(self.release).hexdigest()).start()
        self.download = mock.patch.object(setup, "urlopen").start()
        self.download.return_value.__enter__.return_value = io.BytesIO(self.release)

    def test_valid_archive_installs_verified_tool_license_and_provenance(self):
        path = setup.install(self.destination)
        self.assertEqual(path, self.destination / "tycmd.exe")
        self.assertEqual(path.read_bytes(), self.executable)
        self.assertEqual((self.destination / "LICENSE.txt").read_text(), "fixture license")
        provenance = json.loads((self.destination / "provenance.json").read_text())
        self.assertEqual(provenance["executable_sha256"], hashlib.sha256(self.executable).hexdigest())
        self.assertEqual(provenance["archive_sha256"], hashlib.sha256(self.release).hexdigest())
        self.assertEqual(provenance["url"], setup.TYTOOLS_URL)
        self.download.assert_called_once_with(setup.TYTOOLS_URL, timeout=30)
        self.assertFalse(any(path.name.startswith(".") for path in self.destination.iterdir()))

    def test_corrupt_archive_creates_no_install_directory(self):
        self.download.return_value.__enter__.return_value = io.BytesIO(b"wrong archive")
        with self.assertRaisesRegex(ValueError, "archive SHA-256 mismatch"):
            setup.install(self.destination)
        self.assertFalse(self.destination.exists())

    def test_wrong_executable_inside_known_archive_is_not_installed(self):
        release = archive(b"wrong executable")
        self.download.return_value.__enter__.return_value = io.BytesIO(release)
        with mock.patch.object(setup, "TYTOOLS_ARCHIVE_SHA256", hashlib.sha256(release).hexdigest()):
            with self.assertRaisesRegex(ValueError, "executable SHA-256 mismatch"):
                setup.install(self.destination)
        self.assertFalse(self.destination.exists())

    def test_existing_verified_tool_does_not_download_again(self):
        setup.install(self.destination)
        self.download.reset_mock()
        self.assertEqual(setup.install(self.destination), self.destination / "tycmd.exe")
        self.download.assert_not_called()

    def test_existing_different_tool_is_not_overwritten(self):
        self.destination.mkdir()
        target = self.destination / "tycmd.exe"
        target.write_bytes(b"another tool")
        with self.assertRaisesRegex(ValueError, "different executable already exists"):
            setup.install(self.destination)
        self.assertEqual(target.read_bytes(), b"another tool")
        self.download.assert_not_called()


if __name__ == "__main__":
    unittest.main()
