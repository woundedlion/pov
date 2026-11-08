#!/usr/bin/env bash
# Build the shipping firmware targets through the size and layout gate, then
# compile the profiling variants.
#
# usage: teensy_size.sh <scratch-dir>   holds the profiling variants' build trees
set -euo pipefail

scratch="${1:?usage: teensy_size.sh <scratch-dir>}"
python tools/teensy_size_table.py

profile_variant() {
  PLATFORMIO_BUILD_DIR="$scratch/$1" PLATFORMIO_BUILD_FLAGS="$2" pio run -e profile
}
profile_variant composed-profile-build "-D HS_PROFILE_TARGET=GridWarp -D HS_PROBE_BREAKDOWN -D HS_PLOT_COUNTS -D HS_PROFILE_PLOT_RENDER_COUNTS -D HS_PROFILE_PULLBACK_PROJECTION -D HS_PROFILE_PC -D HS_PROFILE_MODULATION_SOURCE=1 -D HS_PROFILE_MODULATION_MECHANIC=0 -D HS_SCAN_METRICS"
profile_variant plot-stalls-build "-D HS_PROFILE_TARGET=MindSplatter -D HS_PROFILE_PLOT_STALLS"
