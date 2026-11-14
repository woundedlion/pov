#!/usr/bin/env bash
# Build the API reference and enforce documentation coverage. Expects the
# theme applied to Doxyfile.
set -euo pipefail

mkdir -p build/docs
doxygen Doxyfile
python tools/doc_coverage.py
python tools/doc_coverage.py --allowlist tools/ci/doc_coverage_tests_allowlist.txt tests
