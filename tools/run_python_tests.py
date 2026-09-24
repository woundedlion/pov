#!/usr/bin/env python3
"""Run every tracked Python unittest suite, rejecting empty discovery."""

import argparse
from pathlib import Path
import subprocess
import sys
import unittest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--suite", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    root = args.root.resolve()
    if args.suite is not None:
        paths = sorted((root / args.suite).glob("test*.py"))
        if not paths:
            parser.error(f"no test files discovered in {args.suite}")
        failed = False
        for path in paths:
            suite = unittest.defaultTestLoader.discover(str(path.parent), pattern=path.name)
            if suite.countTestCases() == 0:
                parser.error(f"no test cases discovered in {path}")
            result = unittest.TextTestRunner(verbosity=2).run(suite)
            if result.testsRun == len(result.skipped):
                print(f"no unskipped test cases in {path}", file=sys.stderr)
                failed = True
            failed |= not result.wasSuccessful()
        return int(failed)

    paths = subprocess.check_output(
        ["git", "-C", str(root), "ls-files", "-z", "--", "*/test*.py"]
    ).decode("utf-8").split("\0")
    directories = sorted({str(Path(path).parent) for path in paths if path})
    if not directories:
        parser.error("no tracked Python test suites discovered")
    failed = False
    for directory in directories:
        result = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--root", str(root),
             "--suite", directory], cwd=root, check=False)
        failed |= result.returncode != 0
    return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())
