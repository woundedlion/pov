"""Pinned headless uploads select only the board held by the device lock."""

import copy
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import teensy_upload as upload  # noqa: E402


BOARDS = [
    {"action": "add", "tag": "6362410-Teensy", "serial": "6362410",
     "model": "Teensy 4.0", "location": "usb-2-1",
     "capabilities": ["unique", "reboot", "serial"], "interfaces": [["Serial", "COM3"]]},
    {"action": "add", "tag": "6362390-Teensy", "serial": "6362390",
     "model": "Teensy 4.0", "location": "usb-2-9",
     "capabilities": ["unique", "reboot", "serial"], "interfaces": [["Serial", "COM4"]]},
]


class BoardSelectionTests(unittest.TestCase):
    def test_com_selection_keeps_both_serial_and_usb_location(self):
        self.assertEqual(upload.select_board(json.dumps(BOARDS), "COM4"),
                         ("6362390-Teensy@usb-2-9", "6362390", "usb-2-9"))
        self.assertEqual(upload.select_board(json.dumps(BOARDS), "com3")[0],
                         "6362410-Teensy@usb-2-1")

    def test_missing_com_does_not_select_another_attached_board(self):
        with self.assertRaisesRegex(upload.BoardAbsent, "no Teensy at COM9"):
            upload.select_board(json.dumps(BOARDS), "COM9")

    def test_duplicate_com_is_ambiguous(self):
        boards = copy.deepcopy(BOARDS)
        boards[0]["interfaces"] = [["Serial", "COM4"]]
        with self.assertRaisesRegex(upload.UploadError, "ambiguous"):
            upload.select_board(json.dumps(boards), "COM4")

    def test_duplicate_complete_identity_is_ambiguous(self):
        boards = copy.deepcopy(BOARDS)
        boards[0].update({key: boards[1][key] for key in ("tag", "serial", "location")})
        with self.assertRaisesRegex(upload.UploadError, "ambiguous serial/USB"):
            upload.select_board(json.dumps(boards), "COM4")

    def test_incomplete_or_untrusted_identity_is_rejected(self):
        for key, value in (("serial", ""), ("serial", "6362390*"), ("serial", 6362390),
                           ("tag", "-Teensy"), ("location", ""), ("location", "@COM4"),
                           ("capabilities", []), ("model", "Teensy 4.1"),
                           ("action", "remove")):
            with self.subTest(key=key, value=value):
                boards = copy.deepcopy(BOARDS)
                boards[1][key] = value
                with self.assertRaisesRegex(upload.UploadError, "serial/USB identity"):
                    upload.select_board(json.dumps(boards), "COM4")

    def test_invalid_or_malformed_json_is_rejected(self):
        for value in ("invalid", "{}", "[null]", '[{"interfaces":{}}]',
                      '[{"interfaces":[["COM4"]]}]'):
            with self.subTest(value=value):
                with self.assertRaises(upload.UploadError):
                    upload.select_board(value, "COM4")


class UploadTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.tool = self.root / "tycmd.exe"
        self.tool.write_bytes(b"independent tool fixture")
        self.firmware = self.root / "firmware.hex"
        self.firmware.write_text(":00000001FF\n", encoding="ascii")
        self.claim = self.root / "claim"
        self.claim.mkdir()
        self.info = self.claim / "info"
        self.info.write_text("token=fixture-owned\nport=COM4\n", encoding="utf-8")
        self.listing = json.dumps(BOARDS)
        self.addCleanup(mock.patch.stopall)
        mock.patch.dict(os.environ, {"HS_TYCMD": str(self.tool)}).start()
        mock.patch.object(upload, "TYCMD_SHA256", hashlib.sha256(self.tool.read_bytes()).hexdigest()).start()
        mock.patch("sys.stdout", new_callable=io.StringIO).start()
        self.runner = mock.patch.object(upload.subprocess, "run", side_effect=self.respond).start()

    def respond(self, command, **_options):
        if command[1] == "--version":
            output = "tycmd.exe 0.9.9\n"
        elif command[1] == "list":
            output = self.listing
        else:
            output = "uploaded fixture\n"
        return subprocess.CompletedProcess(command, 0, output, "")

    def invoke(self, **options):
        upload.upload(self.firmware, "COM4", self.claim, "fixture-owned", **options)

    def test_upload_uses_exact_board_no_gui_or_manual_wait_and_preserves_rtc(self):
        self.invoke()
        self.assertEqual(self.runner.call_args_list[-1].args[0],
                         [str(self.tool.resolve()), "upload", "--board",
                          "6362390-Teensy@usb-2-9", "--rtc", "none", str(self.firmware.resolve())])
        self.assertEqual(self.runner.call_args_list[1].args[0][1:],
                         ["list", "--verbose", "--output", "json"])
        self.assertEqual([call.kwargs["timeout"] for call in self.runner.call_args_list], [10, 10, 60])

    def test_missing_token_wrong_port_or_replaced_claim_never_starts_tool(self):
        for text in ("token=other\nport=COM4\n", "token=fixture-owned\nport=COM3\n",
                     "token=fixture-owned\ntoken=other\nport=COM4\n"):
            with self.subTest(text=text):
                self.info.write_text(text, encoding="utf-8")
                with self.assertRaises(upload.UploadError):
                    self.invoke()
                self.runner.assert_not_called()
        with self.assertRaisesRegex(upload.UploadError, "requires a device lock"):
            upload.upload(self.firmware, "COM4", self.claim, "")
        self.runner.assert_not_called()

    def test_missing_claim_file_or_firmware_never_starts_tool(self):
        self.info.unlink()
        with self.assertRaisesRegex(upload.UploadError, "cannot read device claim"):
            self.invoke()
        self.runner.assert_not_called()
        self.info.write_text("token=fixture-owned\nport=COM4\n", encoding="utf-8")
        self.firmware.unlink()
        with self.assertRaisesRegex(upload.UploadError, "HEX does not exist"):
            self.invoke()
        self.runner.assert_not_called()

    def test_missing_or_ambiguous_device_never_uploads(self):
        for boards in ([], [BOARDS[1], BOARDS[1]]):
            with self.subTest(boards=boards):
                self.listing = json.dumps(boards)
                with self.assertRaises(upload.UploadError):
                    self.invoke()
                self.assertFalse(any(call.args[0][1] == "upload"
                                     for call in self.runner.call_args_list))

    def test_claim_is_rechecked_after_board_selection_before_reboot(self):
        def change_claim(command, **options):
            response = self.respond(command, **options)
            if command[1] == "list":
                self.info.write_text("token=replaced\nport=COM4\n", encoding="utf-8")
            return response

        self.runner.side_effect = change_claim
        with self.assertRaisesRegex(upload.UploadError, "no longer owns"):
            self.invoke()
        self.assertEqual(self.runner.call_count, 2)

    def test_mismatched_executable_digest_or_version_never_lists_or_uploads(self):
        self.tool.write_bytes(b"unexpected executable")
        with self.assertRaisesRegex(upload.UploadError, "SHA-256 mismatch"):
            self.invoke()
        self.runner.assert_not_called()
        self.tool.write_bytes(b"independent tool fixture")
        self.runner.return_value = subprocess.CompletedProcess([], 0, "tycmd.exe 0.9.8", "")
        self.runner.side_effect = None
        with self.assertRaisesRegex(upload.UploadError, "unsupported tycmd version"):
            self.invoke()
        self.assertEqual(self.runner.call_count, 1)

    def test_process_failure_and_timeout_are_not_success(self):
        for failure in (subprocess.CompletedProcess([], 7, "", "fixture upload failure"),
                        subprocess.TimeoutExpired("tycmd", 60), OSError("fixture process missing")):
            with self.subTest(failure=failure):
                def fail_upload(command, **options):
                    if command[1] != "upload":
                        return self.respond(command, **options)
                    if isinstance(failure, Exception):
                        raise failure
                    return failure

                self.runner.side_effect = fail_upload
                with self.assertRaises(upload.UploadError):
                    self.invoke()

    def test_absent_board_wait_is_bounded_and_can_recover(self):
        listings = iter(("[]", self.listing))

        def transient(command, **options):
            result = self.respond(command, **options)
            if command[1] == "list":
                result.stdout = next(listings)
            return result

        self.runner.side_effect = transient
        with mock.patch.object(upload.time, "monotonic", side_effect=(0, 1, 1)), \
                mock.patch.object(upload.time, "sleep") as sleep:
            self.invoke(absent_seconds=5)
        sleep.assert_called_once_with(1)
        self.assertEqual(self.runner.call_args_list[-1].args[0][1], "upload")

    def test_resolve_mode_only_verifies_tool_without_board_access(self):
        self.assertEqual(upload.main(["--resolve"]), 0)
        self.assertEqual(self.runner.call_count, 1)
        self.assertEqual(self.runner.call_args.args[0][1:], ["--version"])

    def test_upload_cli_requires_claim_before_any_process(self):
        with mock.patch("sys.stderr", new_callable=io.StringIO), self.assertRaises(SystemExit) as error:
            upload.main([str(self.firmware), "--port", "COM4"])
        self.assertEqual(error.exception.code, 2)
        self.runner.assert_not_called()


if __name__ == "__main__":
    unittest.main()
