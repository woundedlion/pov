"""Capture the profile image's USB-serial readout from an attached Teensy.

Companion to the `profile` PlatformIO env (targets/Profile/Profile.ino): after
`pio run -e profile -t upload`, the board re-enumerates and starts streaming
HS_PROFILE cycle-counter dumps. This opens the Teensy's serial port (retrying
through the re-enumeration window), tees every line to stdout and --out, and
exits after --seconds.

Uses pyserial, which PlatformIO installs (`python -m pip install pyserial`
otherwise).
"""

import argparse
import os
import sys
import time

import serial
from serial.tools import list_ports

TEENSY_VID = 0x16C0


def find_port(want=None):
    """Return the device name of the Teensy serial port to capture from.

    With more than one Teensy attached, enumeration order is not stable, so
    the first VID match can be the wrong board. `want` (--port /
    HS_TEENSY_PORT) pins the capture to one device name.
    """
    for p in list_ports.comports():
        if p.vid == TEENSY_VID and (want is None or p.device == want):
            return p.device
    return None


def open_port(timeout_s, want=None):
    """Open the Teensy port, retrying through USB re-enumeration after upload."""
    deadline = time.monotonic() + timeout_s
    last_err = None
    while time.monotonic() < deadline:
        port = find_port(want)
        if port:
            try:
                return serial.Serial(port, baudrate=115200, timeout=1)
            except serial.SerialException as e:  # enumerating; not ready yet
                last_err = e
        time.sleep(0.5)
    which = want or "any"
    raise SystemExit(f"profile_capture: no Teensy serial port [{which}] "
                     f"({last_err})")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seconds", type=float, default=150.0,
                    help="capture duration (default 150 s: covers a full "
                         "DisplacementField noise+balls phase cycle)")
    ap.add_argument("--out", default="build/profile_capture.log",
                    help="file the captured lines are written to")
    ap.add_argument("--connect-timeout", type=float, default=30.0,
                    help="seconds to wait for the port to enumerate")
    ap.add_argument("--port", default=os.environ.get("HS_TEENSY_PORT"),
                    help="capture from this device name (e.g. COM3) instead "
                         "of the first Teensy found; defaults to "
                         "$HS_TEENSY_PORT. Required when the host has more "
                         "than one Teensy attached")
    args = ap.parse_args()

    ser = open_port(args.connect_timeout, args.port)
    print(f"profile_capture: reading {ser.port} for {args.seconds:.0f} s",
          flush=True)
    end = time.monotonic() + args.seconds
    with open(args.out, "w", encoding="utf-8", newline="\n") as f:
        while time.monotonic() < end:
            line = ser.readline()  # 1 s timeout keeps the deadline responsive
            if not line:
                continue
            text = line.decode("utf-8", errors="replace").rstrip("\r\n")
            print(text, flush=True)
            f.write(text + "\n")
    ser.close()
    print(f"profile_capture: wrote {args.out}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
