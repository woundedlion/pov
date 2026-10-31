#!/usr/bin/env python3
"""Upload through tycmd to the exact board held by the shared device lock."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

TYTOOLS_VERSION = "0.9.9"
TYTOOLS_URL = f"https://github.com/Koromix/tytools/releases/download/v{TYTOOLS_VERSION}/tytools_{TYTOOLS_VERSION}_win64.zip"
TYTOOLS_ARCHIVE_SHA256 = "1ff7c945de4789950fb6f2d56bac8523e849b86a44c47b87e2c77d130af559b2"
TYCMD_SHA256 = "cd4d199bd33bce2143ce76bfb721587f3bd0b17a3375c613e86834435fc3a5d0"
DEFAULT_TYCMD = Path.home() / ".holosphere" / "tools" / f"tytools-{TYTOOLS_VERSION}" / "tycmd.exe"


class UploadError(RuntimeError):
    """The claimed board could not be selected or programmed."""


class BoardAbsent(UploadError):
    """The claimed COM interface is temporarily absent."""


class ToolAbsent(UploadError):
    """No headless uploader was installed or selected."""


def verify_lock(directory: Path, token: str, port: str) -> None:
    """Require the same claim token and serial port as the lock owner."""
    if not token:
        raise UploadError("flash requires a device lock")
    try:
        lines = (directory / "info").read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as error:
        raise UploadError(f"cannot read device claim: {error}") from error
    fields = {}
    for line in lines:
        key, separator, value = line.partition("=")
        if separator and key in ("token", "port"):
            if key in fields:
                raise UploadError(f"duplicate {key} in device claim")
            fields[key] = value
    if fields.get("token") != token or fields.get("port", "").upper() != port.upper():
        raise UploadError("device claim no longer owns the requested port")


def run_tool(command: list[str], timeout: float) -> str:
    """Bound each tool process and preserve failures in its diagnostic."""
    try:
        result = subprocess.run(command, capture_output=True, text=True,
                                encoding="utf-8", errors="replace", timeout=timeout,
                                check=False)
    except subprocess.TimeoutExpired as error:
        raise UploadError(f"{Path(command[0]).name} timed out after {timeout:g}s") from error
    except OSError as error:
        raise UploadError(f"cannot execute {command[0]}: {error}") from error
    if result.returncode:
        detail = (result.stderr or result.stdout).strip()
        raise UploadError(f"{Path(command[0]).name} exited {result.returncode}: {detail}")
    return result.stdout


def resolve_tycmd() -> Path:
    """Resolve an explicit override, pinned local install, or PATH tool."""
    candidate = (os.environ.get("HS_TYCMD") or
                 (str(DEFAULT_TYCMD) if DEFAULT_TYCMD.is_file() else shutil.which("tycmd")))
    if not candidate:
        raise ToolAbsent("tycmd not found; run tools/setup_teensy_uploader.py or set HS_TYCMD")
    path = Path(candidate).expanduser().resolve()
    if not path.is_file():
        raise ToolAbsent("tycmd not found; set HS_TYCMD to the pinned tycmd.exe")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != TYCMD_SHA256:
        raise UploadError(f"tycmd SHA-256 mismatch: {path} ({digest})")
    version = run_tool([str(path), "--version"], 10).strip()
    if version not in (f"tycmd {TYTOOLS_VERSION}", f"tycmd.exe {TYTOOLS_VERSION}"):
        raise UploadError(f"unsupported tycmd version: {version}")
    return path


def select_board(output: str, port: str) -> tuple[str, str, str]:
    """Resolve one COM interface to a complete serial/family/USB-location tag."""
    try:
        boards = json.loads(output)
    except (ValueError, TypeError) as error:
        raise UploadError("tycmd returned invalid board JSON") from error
    if not isinstance(boards, list) or any(not isinstance(board, dict) for board in boards):
        raise UploadError("tycmd board JSON is not a board list")
    matches = []
    for board in boards:
        interfaces = board.get("interfaces", [])
        if not isinstance(interfaces, list):
            raise UploadError("tycmd returned malformed board interfaces")
        if any(not isinstance(item, list) or len(item) != 2 for item in interfaces):
            raise UploadError("tycmd returned malformed board interface")
        if any(item[0] == "Serial" and str(item[1]).upper() == port.upper() for item in interfaces):
            matches.append(board)
    if not matches:
        raise BoardAbsent(f"no Teensy at {port}")
    if len(matches) != 1:
        raise UploadError(f"ambiguous Teensy interface at {port}: {len(matches)} boards")
    board = matches[0]
    tag, serial, location = board.get("tag"), board.get("serial"), board.get("location")
    capabilities = board.get("capabilities", [])
    if (not isinstance(serial, str) or not re.fullmatch(r"[1-9][0-9]*", serial)
            or not isinstance(tag, str) or tag != f"{serial}-Teensy"
            or not isinstance(location, str) or not re.fullmatch(r"usb(?:-[0-9]+)+", location)
            or not isinstance(capabilities, list) or "unique" not in capabilities
            or board.get("model") != "Teensy 4.0" or board.get("action") != "add"):
        raise UploadError(f"{port} has no unambiguous Teensy 4.0 serial/USB identity")
    if sum(other.get("tag") == tag and other.get("location") == location for other in boards) != 1:
        raise UploadError(f"ambiguous serial/USB identity at {port}")
    return f"{tag}@{location}", serial, location


def upload(firmware: Path, port: str, directory: Path, token: str,
           absent_seconds: int = 0, timeout: int = 60) -> None:
    """Upload with a stable physical identity and recheck ownership before reboot."""
    if not re.fullmatch(r"COM[1-9][0-9]*", port, re.IGNORECASE):
        raise UploadError("device lock did not pin a COM port")
    verify_lock(directory, token, port)
    if not firmware.is_file() or firmware.suffix.lower() != ".hex":
        raise UploadError(f"firmware HEX does not exist: {firmware}")
    tool = resolve_tycmd()
    deadline = time.monotonic() + absent_seconds
    while True:
        listing = run_tool([str(tool), "list", "--verbose", "--output", "json"], 10)
        try:
            tag, serial, location = select_board(listing, port)
            break
        except BoardAbsent:
            if time.monotonic() >= deadline:
                raise
            verify_lock(directory, token, port)
            time.sleep(min(1, max(0, deadline - time.monotonic())))
    verify_lock(directory, token, port)
    digest = hashlib.sha256(firmware.read_bytes()).hexdigest()
    print(f"uploader: tycmd {TYTOOLS_VERSION} sha256={TYCMD_SHA256} path={tool}", flush=True)
    print(f"device: {port} serial={serial} location={location} board={tag}", flush=True)
    print(f"firmware: {firmware.resolve()} sha256={digest}", flush=True)
    output = run_tool([str(tool), "upload", "--board", tag, "--rtc", "none", str(firmware.resolve())], timeout)
    if output:
        print(output, end="" if output.endswith("\n") else "\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("firmware", type=Path, nargs="?")
    parser.add_argument("--resolve", action="store_true", help="verify and print tycmd; no device access")
    parser.add_argument("--port")
    parser.add_argument("--lock-dir", type=Path)
    parser.add_argument("--lock-token")
    parser.add_argument("--absent-seconds", type=int, default=0)
    parser.add_argument("--timeout", type=int, default=60)
    args = parser.parse_args(argv)
    if args.resolve:
        try:
            print(resolve_tycmd())
        except ToolAbsent as error:
            print(f"upload: {error}", file=sys.stderr)
            return 2
        except UploadError as error:
            print(f"upload: {error}", file=sys.stderr)
            return 1
        return 0
    if args.firmware is None or not args.port or args.lock_dir is None or not args.lock_token:
        parser.error("upload requires firmware, --port, --lock-dir and --lock-token")
    if not 0 <= args.absent_seconds <= 3600 or not 1 <= args.timeout <= 300:
        parser.error("absent-seconds must be 0..3600 and timeout 1..300")
    try:
        upload(args.firmware, args.port, args.lock_dir, args.lock_token,
               args.absent_seconds, args.timeout)
    except (UploadError, OSError) as error:
        print(f"upload: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
