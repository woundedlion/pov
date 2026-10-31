#!/usr/bin/env python3
"""Install the checksum-pinned Windows tycmd executable without changing PATH."""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import sys
import tempfile
from urllib.request import urlopen
import zipfile

from teensy_upload import DEFAULT_TYCMD, TYCMD_SHA256, TYTOOLS_ARCHIVE_SHA256, TYTOOLS_URL, TYTOOLS_VERSION


def install(directory: Path) -> Path:
    """Verify the release archive and executable before installing any bytes."""
    destination = directory / "tycmd.exe"
    if destination.exists():
        if hashlib.sha256(destination.read_bytes()).hexdigest() != TYCMD_SHA256:
            raise ValueError(f"a different executable already exists: {destination}")
        return destination
    with urlopen(TYTOOLS_URL, timeout=30) as response:
        archive = response.read(32 * 1024 * 1024 + 1)
    if hashlib.sha256(archive).hexdigest() != TYTOOLS_ARCHIVE_SHA256:
        raise ValueError("TyTools release archive SHA-256 mismatch")
    with zipfile.ZipFile(io.BytesIO(archive)) as release:
        files = {name: release.read(name) for name in ("tycmd.exe", "LICENSE.txt", "README.md")}
    if hashlib.sha256(files["tycmd.exe"]).hexdigest() != TYCMD_SHA256:
        raise ValueError("tycmd executable SHA-256 mismatch")
    files["provenance.json"] = (json.dumps({
        "version": TYTOOLS_VERSION, "url": TYTOOLS_URL,
        "archive_sha256": TYTOOLS_ARCHIVE_SHA256, "executable_sha256": TYCMD_SHA256,
    }, indent=2) + "\n").encode("utf-8")
    directory.mkdir(parents=True, exist_ok=True)
    for name, data in files.items():
        descriptor, temporary = tempfile.mkstemp(prefix=f".{name}.", dir=directory)
        try:
            with os.fdopen(descriptor, "wb") as output:
                output.write(data)
                output.flush()
                os.fsync(output.fileno())
            os.replace(temporary, directory / name)
        finally:
            Path(temporary).unlink(missing_ok=True)
    return destination


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=DEFAULT_TYCMD.parent)
    args = parser.parse_args(argv)
    try:
        path = install(args.directory)
    except (OSError, ValueError, KeyError, zipfile.BadZipFile) as error:
        print(f"tycmd setup: {error}", file=sys.stderr)
        return 1
    print(path.resolve().as_posix())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
