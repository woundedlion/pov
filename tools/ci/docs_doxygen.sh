#!/usr/bin/env bash
# Build the API reference and enforce documentation coverage. Expects the
# theme applied to Doxyfile.
set -euo pipefail

mkdir -p build/docs
doxygen Doxyfile
python tools/doc_coverage.py
