#!/usr/bin/env python3
"""Host tests for the Teensy serial capture (tools/profile_capture.py).

The pyserial transport is stubbed before the module under test imports it;
every test drives that stub.

Run:  python -m unittest discover -s tools/profile_tests
"""

import contextlib
import io
import os
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

TOOLS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(TOOLS))


class _SerialException(Exception):
    pass


def _install_serial_stub():
    """Put a minimal pyserial surface on sys.modules, even over a real one."""
    serial = types.ModuleType("serial")
    serial.SerialException = _SerialException
    serial.Serial = None  # every test patches this
    tools = types.ModuleType("serial.tools")
    ports = types.ModuleType("serial.tools.list_ports")
    ports.comports = lambda: []
    tools.list_ports = ports
    serial.tools = tools
    sys.modules["serial"] = serial
    sys.modules["serial.tools"] = tools
    sys.modules["serial.tools.list_ports"] = ports


_install_serial_stub()

import profile_capture as pc   # noqa: E402

TEENSY = pc.TEENSY_VID
OTHER_VID = 0x0403  # an FTDI adapter — present on a bench, never the target


class FakePort:
    """A `list_ports.comports()` row."""

    def __init__(self, device, vid=TEENSY):
        self.device = device
        self.vid = vid


class FakeClock:
    """Stands in for the `time` module; time passes only where the code blocks."""

    def __init__(self):
        self.now = 0.0
        self.sleeps = []

    def monotonic(self):
        return self.now

    def sleep(self, seconds):
        self.sleeps.append(seconds)
        self.now += seconds


class FakeSerial:
    """An opened port that replays a fixed sequence of readline() results.

    An empty read costs the port's 1 s timeout on `clock`; a line costs 10 ms.
    """

    def __init__(self, port, lines=(), clock=None):
        self.port = port
        self.lines = list(lines)
        self.clock = clock
        self.closed = False

    def readline(self):
        line = self.lines.pop(0) if self.lines else b""
        self.clock.now += 0.01 if line else 1.0
        return line

    def close(self):
        self.closed = True


def ports(*rows):
    """Patch the bus to list `rows`; a runaway poll fails at 100 enumerations."""
    polls = []

    def comports():
        polls.append(None)
        if len(polls) > 100:
            raise AssertionError("port enumeration polled without backing off")
        return list(rows)

    return mock.patch.object(pc.list_ports, "comports", comports)


class TestFindPort(unittest.TestCase):
    """Which board the capture picks off the bus."""

    def test_single_teensy_is_found_without_a_pin(self):
        with ports(FakePort("COM3")):
            self.assertEqual(pc.find_port(), "COM3")

    def test_non_teensy_devices_are_ignored(self):
        with ports(FakePort("COM1", OTHER_VID), FakePort("COM7")):
            self.assertEqual(pc.find_port(), "COM7")

    def test_no_teensy_attached_reports_nothing(self):
        with ports(FakePort("COM1", OTHER_VID)):
            self.assertIsNone(pc.find_port())
        with ports():
            self.assertIsNone(pc.find_port())

    def test_two_teensys_the_pin_selects_the_named_board(self):
        # Enumeration order is not stable across sessions, so the pinned name
        # has to win from either position in the list.
        rows = (FakePort("COM3"), FakePort("COM9"))
        with ports(*rows):
            self.assertEqual(pc.find_port("COM9"), "COM9")
            self.assertEqual(pc.find_port("COM3"), "COM3")
        with ports(*reversed(rows)):
            self.assertEqual(pc.find_port("COM9"), "COM9")
            self.assertEqual(pc.find_port("COM3"), "COM3")

    def test_two_teensys_without_a_pin_takes_the_first_enumerated(self):
        # Unpinned, whichever board enumerated first wins.
        with ports(FakePort("COM3"), FakePort("COM9")):
            self.assertEqual(pc.find_port(), "COM3")
        with ports(FakePort("COM9"), FakePort("COM3")):
            self.assertEqual(pc.find_port(), "COM9")

    def test_a_pin_naming_an_absent_board_matches_nothing(self):
        # It must not silently fall back to the other attached Teensy.
        with ports(FakePort("COM3"), FakePort("COM9")):
            self.assertIsNone(pc.find_port("COM5"))

    def test_a_pin_naming_a_non_teensy_port_matches_nothing(self):
        with ports(FakePort("COM1", OTHER_VID)):
            self.assertIsNone(pc.find_port("COM1"))


class TestOpenPort(unittest.TestCase):
    """The connect window: re-enumeration, contention, and giving up."""

    def setUp(self):
        self.clock = FakeClock()
        self.enter = mock.patch.object(pc, "time", self.clock)
        self.enter.start()
        self.addCleanup(self.enter.stop)

    def _serial(self, *results):
        """Patch serial.Serial with one result (value or exception) per call.

        The last result repeats; a runaway retry loop fails at 100 calls.
        """
        calls = []

        def factory(port, **kwargs):
            calls.append((port, kwargs))
            if len(calls) > 100:
                raise AssertionError("open_port retried without backing off")
            result = results[min(len(calls), len(results)) - 1]
            if isinstance(result, Exception):
                raise result
            return result

        return calls, mock.patch.object(pc.serial, "Serial", factory)

    def test_open_returns_the_port_at_the_configured_baud(self):
        opened = FakeSerial("COM3")
        calls, patched = self._serial(opened)
        with ports(FakePort("COM3")), patched:
            self.assertIs(pc.open_port(3.0), opened)
        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0][0], "COM3")
        self.assertEqual(calls[0][1]["baudrate"], 115200)
        self.assertEqual(calls[0][1]["timeout"], 1)

    def test_open_honours_the_board_pin_with_two_attached(self):
        opened = FakeSerial("COM9")
        calls, patched = self._serial(opened)
        with ports(FakePort("COM3"), FakePort("COM9")), patched:
            self.assertIs(pc.open_port(3.0, "COM9"), opened)
        self.assertEqual(calls[0][0], "COM9")

    def test_open_retries_through_re_enumeration(self):
        opened = FakeSerial("COM3")
        calls, patched = self._serial(_SerialException("not ready"), opened)
        with ports(FakePort("COM3")), patched:
            self.assertIs(pc.open_port(3.0), opened)
        self.assertEqual(len(calls), 2)
        self.assertEqual(self.clock.sleeps, [0.5])

    def test_a_held_port_exhausts_the_window_and_names_the_cause(self):
        # A peer session still holding the port: every attempt raises, so the
        # deadline has to end the retry loop rather than spin forever.
        busy = _SerialException("could not open port COM3: Access is denied.")
        calls, patched = self._serial(busy)
        with ports(FakePort("COM3")), patched:
            with self.assertRaises(SystemExit) as caught:
                pc.open_port(3.0, "COM3")
        message = str(caught.exception)
        self.assertIn("COM3", message)
        self.assertIn("Access is denied.", message)
        self.assertGreaterEqual(len(calls), 2)
        self.assertGreaterEqual(sum(self.clock.sleeps), 3.0 - 0.5)
        self.assertLessEqual(sum(self.clock.sleeps), 3.0)

    def test_no_board_at_all_exits_naming_any(self):
        with ports(), mock.patch.object(pc.serial, "Serial", None):
            with self.assertRaises(SystemExit) as caught:
                pc.open_port(3.0)
        self.assertIn("[any]", str(caught.exception))
        self.assertIn("VID 0x16C0 enumerated", str(caught.exception))
        self.assertIn("within 3 s", str(caught.exception))
        self.assertNotIn("None", str(caught.exception))

    def test_an_absent_pinned_board_exits_naming_the_pin(self):
        # The other Teensy is attached and openable; the pin must still fail.
        with ports(FakePort("COM3")), mock.patch.object(
                pc.serial, "Serial", None):
            with self.assertRaises(SystemExit) as caught:
                pc.open_port(3.0, "COM9")
        self.assertIn("[COM9]", str(caught.exception))

    def test_a_zero_window_never_touches_the_bus(self):
        with ports(FakePort("COM3")), mock.patch.object(
                pc.serial, "Serial", None):
            with self.assertRaises(SystemExit):
                pc.open_port(0.0)
        self.assertEqual(self.clock.sleeps, [])


class TestMain(unittest.TestCase):
    """The capture loop: what lands on disk, and what is refused."""

    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.dir.cleanup)
        self.out = os.path.join(self.dir.name, "nested", "capture.log")
        self.clock = FakeClock()
        patched = mock.patch.object(pc, "time", self.clock)
        patched.start()
        self.addCleanup(patched.stop)

    def _run(self, serial_port, *extra):
        """Drive main() over a fake port; returns what it teed to stdout.

        stdout is captured: a cp1252 Windows console cannot print the
        replacement characters of the undecodable-bytes case.
        """
        argv = ["profile_capture.py", "--out", self.out, "--seconds", "5.0"]
        stdout = io.StringIO()
        with mock.patch.object(sys, "argv", argv + list(extra)), \
                mock.patch.object(pc, "open_port", lambda *a: serial_port), \
                contextlib.redirect_stdout(stdout):
            pc.main()
        return stdout.getvalue()

    def test_capture_writes_the_streamed_lines_and_creates_the_directory(self):
        opened = FakeSerial("COM3", [b"phase 1\r\n", b"", b"phase 2\n"], self.clock)
        printed = self._run(opened)
        self.assertTrue(opened.closed)
        with open(self.out, encoding="utf-8", newline="") as handle:
            self.assertEqual(handle.read(), "phase 1\nphase 2\n")
        # The capture is a tee: the same lines reach the console live.
        self.assertIn("phase 1", printed)
        self.assertIn("phase 2", printed)
        self.assertIn("2 lines", printed)

    def test_undecodable_bytes_do_not_abort_the_capture(self):
        opened = FakeSerial("COM3", [b"\xff\xfe ok\n"], self.clock)
        self._run(opened)
        with open(self.out, encoding="utf-8") as handle:
            self.assertIn("ok", handle.read())

    def test_a_serial_drop_removes_the_partial_capture(self):
        opened = FakeSerial("COM3")
        stderr = io.StringIO()
        with mock.patch.object(opened, "readline", side_effect=[
                b"HS_CHECK failed: breadcrumb\n", _SerialException("USB disconnected")]), \
                contextlib.redirect_stderr(stderr):
            with self.assertRaises(SystemExit) as caught:
                self._run(opened)
        self.assertIn("HS_CHECK failed: breadcrumb", stderr.getvalue())
        self.assertIn("COM3 dropped after 1 lines: USB disconnected",
                      str(caught.exception))
        self.assertFalse(os.path.exists(self.out))
        self.assertTrue(opened.closed)

    def test_a_silent_board_leaves_no_log_behind(self):
        # A wrong or hung image enumerates and streams nothing; an empty log
        # reads downstream as a real capture.
        opened = FakeSerial("COM3", clock=self.clock)
        with self.assertRaises(SystemExit) as caught:
            self._run(opened)
        self.assertIn("streamed nothing", str(caught.exception))
        self.assertFalse(os.path.exists(self.out))
        self.assertTrue(opened.closed)

    def test_a_silent_board_leaves_no_remnant_of_a_prior_capture(self):
        # Opening for write already truncated the previous run's log, so the
        # removal has to happen or a zero-byte file survives as a "capture".
        os.makedirs(os.path.dirname(self.out))
        with open(self.out, "w", encoding="utf-8") as handle:
            handle.write("previous run\n")
        with self.assertRaises(SystemExit):
            self._run(FakeSerial("COM3", clock=self.clock))
        self.assertFalse(os.path.exists(self.out))

    def test_a_nonpositive_duration_is_rejected(self):
        for seconds in ("0", "-5"):
            with self.subTest(seconds=seconds):
                argv = ["profile_capture.py", "--out", self.out,
                        "--seconds", seconds]
                with mock.patch.object(sys, "argv", argv), \
                        mock.patch.object(sys, "stderr", mock.MagicMock()):
                    with self.assertRaises(SystemExit) as caught:
                        pc.main()
                self.assertEqual(caught.exception.code, 2)
                self.assertFalse(os.path.exists(self.out))


class TestBoardPinForwarding(unittest.TestCase):
    """main() hands the board pin to the connect step, from flag or bench env."""

    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.dir.cleanup)
        self.out = os.path.join(self.dir.name, "capture.log")
        self.clock = FakeClock()
        patched = mock.patch.object(pc, "time", self.clock)
        patched.start()
        self.addCleanup(patched.stop)

    def _captured_pin(self, environ, *extra):
        seen = {}

        def fake_open_port(timeout_s, want=None):
            seen["timeout"] = timeout_s
            seen["want"] = want
            return FakeSerial("COM9", [b"line\n"], self.clock)

        argv = ["profile_capture.py", "--out", self.out, "--seconds", "1.0"]
        with mock.patch.dict(os.environ, {}, clear=False), \
                mock.patch.object(sys, "argv", argv + list(extra)), \
                mock.patch.object(pc, "open_port", fake_open_port), \
                contextlib.redirect_stdout(io.StringIO()):
            os.environ.pop("HS_TEENSY_PORT", None)
            os.environ.update(environ)
            pc.main()
        return seen

    def test_the_bench_environment_variable_is_the_default_pin(self):
        seen = self._captured_pin({"HS_TEENSY_PORT": "COM9"})
        self.assertEqual(seen["want"], "COM9")

    def test_the_flag_overrides_the_environment_pin(self):
        seen = self._captured_pin({"HS_TEENSY_PORT": "COM9"},
                                  "--port", "COM3")
        self.assertEqual(seen["want"], "COM3")

    def test_no_pin_anywhere_leaves_the_board_unpinned(self):
        seen = self._captured_pin({})
        self.assertIsNone(seen["want"])

    def test_the_connect_window_is_forwarded(self):
        seen = self._captured_pin({}, "--connect-timeout", "7.5")
        self.assertEqual(seen["timeout"], 7.5)


if __name__ == "__main__":
    unittest.main()
