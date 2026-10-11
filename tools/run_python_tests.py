#!/usr/bin/env python3
"""Run every tracked Python unittest suite, rejecting empty discovery."""

import argparse
from pathlib import Path
import subprocess
import sys
import unittest

from git_files import git_nul_list


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
        executed = 0
        for path in paths:
            suite = unittest.defaultTestLoader.discover(str(path.parent), pattern=path.name)
            if suite.countTestCases() == 0:
                parser.error(f"no test cases discovered in {path}")
            result = unittest.TextTestRunner(verbosity=2).run(suite)
            executed += result.testsRun - len(result.skipped)
            failed |= not result.wasSuccessful()
        if executed == 0:
            print(f"no unskipped test cases in {args.suite}", file=sys.stderr)
            failed = True
        return int(failed)

    paths = git_nul_list(root, "ls-files", "-z", "--", "*/test*.py", timeout=None)
    directories = sorted({str(Path(path).parent) for path in paths if Path(path).name.startswith("test")})
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
