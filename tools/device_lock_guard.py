#!/usr/bin/env python3
"""Serialize per-board claim creation and token-checked removal."""

import contextlib
import errno
import os
from pathlib import Path
import shutil
import sys
import time

if os.name == "nt":
    import msvcrt
else:
    import fcntl


@contextlib.contextmanager
def guard(directory):
    # The guard inode must persist across all claim generations.
    with Path(f"{directory}.guard").open("a+b") as stream:
        stream.seek(0)
        if os.name == "nt":
            while True:
                try:
                    msvcrt.locking(stream.fileno(), msvcrt.LK_NBLCK, 1)
                    break
                except OSError as error:
                    if error.errno not in (errno.EACCES, errno.EAGAIN):
                        raise
                    time.sleep(0.05)
        else:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            stream.seek(0)
            if os.name == "nt":
                msvcrt.locking(stream.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(stream.fileno(), fcntl.LOCK_UN)


def read_token(directory):
    try:
        for line in (directory / "info").read_text().splitlines():
            if line.startswith("token="):
                return line[6:]
    except OSError:
        pass
    return ""


def update_claim(directory, operation, value):
    directory = Path(directory)
    # Opening the guard file under a missing parent raises ENOENT, which the
    # handler below swallows as a lost race; the caller then reads a
    # misconfigured lock path as a device somebody else is holding.
    if not directory.parent.is_dir():
        print(f"device: lock root {directory.parent} does not exist",
              file=sys.stderr)
        return False
    try:
        with guard(directory):
            if operation == "claim":
                directory.mkdir()
                try:
                    (directory / "info").write_text(value)
                except OSError:
                    shutil.rmtree(directory)
                    raise
            elif operation == "break":
                if not directory.is_dir() or read_token(directory) != value:
                    return False
                shutil.rmtree(directory)
            else:
                raise ValueError(operation)
    except OSError as error:
        if error.errno not in (errno.EACCES, errno.EAGAIN, errno.EEXIST, errno.ENOENT):
            print(f"device: cannot update {directory}: {error}", file=sys.stderr)
        return False
    return True


if __name__ == "__main__":
    operation, directory = sys.argv[1:3]
    value = sys.stdin.read() if operation == "claim" else sys.argv[3]
    try:
        result = update_claim(directory, operation, value)
    except Exception as error:
        print(f"device: lock guard failed: {error}", file=sys.stderr)
        sys.exit(2)
    sys.exit(0 if result else 1)
