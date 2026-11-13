#!/usr/bin/env python3
"""Run every tracked Python unittest suite, rejecting empty discovery."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import os
from pathlib import Path
import subprocess
import sys
import unittest

from git_files import git_nul_list

# Exit status of a --only run whose cases all skipped; a directory fails when
# every one of its files reports it.
ALL_SKIPPED = 3


def run_suite(parser, root, suite, only):
    paths = sorted((root / suite).glob(only or "test*.py"))
    if not paths:
        parser.error(f"no test files discovered in {suite}")
    failed = False
    executed = 0
    for path in paths:
        tests = unittest.defaultTestLoader.discover(str(path.parent), pattern=path.name)
        if tests.countTestCases() == 0:
            parser.error(f"no test cases discovered in {path}")
        result = unittest.TextTestRunner(verbosity=2).run(tests)
        executed += result.testsRun - len(result.skipped)
        failed |= not result.wasSuccessful()
    if failed:
        return 1
    if executed == 0:
        if only:
            return ALL_SKIPPED
        print(f"no unskipped test cases in {suite}", file=sys.stderr)
        return 1
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--suite", type=Path, help="Run only the test*.py files directly in this directory, relative to --root.")
    parser.add_argument("--only", help="With --suite, run only this test file of the directory.")
    parser.add_argument("--exclude", action="append", default=[], type=Path,
                        help="Skip this suite directory, relative to --root; repeatable.")
    args = parser.parse_args()
    root = args.root.resolve()
    if args.suite is not None:
        return run_suite(parser, root, args.suite, args.only)

    paths = git_nul_list(root, "ls-files", "-z", "--", "*/test*.py", timeout=None)
    excluded = {str(directory) for directory in args.exclude}
    files = sorted({(str(Path(path).parent), Path(path).name) for path in paths
                    if Path(path).name.startswith("test")
                    and (root / path).is_file()
                    and str(Path(path).parent) not in excluded})
    if not files:
        parser.error("no tracked Python test suites discovered")

    def run_file(entry):
        directory, name = entry
        return entry, subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--root", str(root),
             "--suite", directory, "--only", name],
            cwd=root, check=False, capture_output=True, text=True)

    failed = False
    unskipped = {directory: False for directory, _ in files}
    # One process per file; more workers than this oversubscribe the
    # subprocess-heavy suites until their timeouts trip.
    workers = int(os.environ.get("HS_PYTHON_TEST_WORKERS", min(8, os.cpu_count() or 1)))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for (directory, _name), result in pool.map(run_file, files):
            sys.stdout.write(result.stdout)
            sys.stderr.write(result.stderr)
            if result.returncode != ALL_SKIPPED:
                unskipped[directory] = True
                failed |= result.returncode != 0
    for directory, ran in unskipped.items():
        if not ran:
            print(f"no unskipped test cases in {directory}", file=sys.stderr)
            failed = True
    return int(failed)


if __name__ == "__main__":
    raise SystemExit(main())
