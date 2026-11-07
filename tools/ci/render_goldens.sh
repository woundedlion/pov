#!/usr/bin/env bash
# Build and run the native render fingerprint comparison in a configured tree.
#
# usage: render_goldens.sh [build-dir]   default build/tests
set -euo pipefail

dir="${1:-build/tests}"
cmake --build "$dir" --target render_goldens
HS_SKIPS_ARE_ERRORS=1 "$dir/tests/render_goldens" | tee "$dir/render-goldens.log"
