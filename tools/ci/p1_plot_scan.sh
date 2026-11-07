#!/usr/bin/env bash
# Run the plot and scan modules against the physical display profile.
#
# usage: p1_plot_scan.sh [cmake configure args...]
set -euo pipefail

cmake --preset tests -B build/tests-p1 -DHS_INSTALL_GIT_HOOKS=OFF -DHS_NATIVE_DISPLAY_PROFILE=1 "$@"
cmake --build build/tests-p1 --target test_engine
HS_SMOKE_FRAMES=120 HS_SKIPS_ARE_ERRORS=1 HS_BUFFER_FREE_WATCHDOG_US=30000000 \
  build/tests-p1/tests/test_engine plot_scan
