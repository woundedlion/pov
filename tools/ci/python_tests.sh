#!/usr/bin/env bash
# Run every Python suite and check routed PCB candidate provenance and metadata.
#
# usage: python_tests.sh [run_python_tests.py args...]
set -euo pipefail

python tools/run_python_tests.py "$@"
python hardware/phantasm/gen/release.py check-candidates
python hardware/phantasm/gen/board_metadata.py --check
