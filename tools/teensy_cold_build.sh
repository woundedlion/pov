#!/usr/bin/env bash
# Cold `pio run -v` over every environment and profiling variant, teed to
# teensy_build.log for tools/teensy_warnings.py.
#
# A cached translation unit emits no warnings, so the object cache and
# .pio/build are wiped first; budget tens of minutes.
#
# usage: teensy_cold_build.sh [log]
set -euo pipefail

if [ "$#" -gt 1 ]; then
  echo "usage: $0 [log]" >&2
  exit 2
fi
log=${1:-teensy_build.log}
case "$log" in
  /*|[A-Za-z]:[\\/]*) ;;
  *) log="$PWD/$log" ;;
esac
# shellcheck source-path=SCRIPTDIR source=device_lock.sh
. "$(dirname "$0")/device_lock.sh"
TREE=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
cd "$TREE"
hs_exit_traps
acquire_tree_lock || exit "$?"

export LC_ALL=C
rm -rf .pio/build_cache .pio/build
pio run -v 2>&1 | tee "$log"
PLATFORMIO_BUILD_FLAGS='-D HS_PROFILE_TARGET=MindSplatter -D HS_MINDSPLATTER_REPLAY' \
  pio run -v -e profile 2>&1 | tee -a "$log"
PLATFORMIO_BUILD_FLAGS='-D HS_PROFILE_TARGET=MindSplatter -D HS_PROFILE_PC' \
  pio run -v -e profile 2>&1 | tee -a "$log"
PLATFORMIO_BUILD_FLAGS='-D HS_PROFILE_TARGET=GridWarp -D HS_PROBE_BREAKDOWN -D HS_PLOT_COUNTS -D HS_PROFILE_PLOT_RENDER_COUNTS -D HS_PROFILE_PULLBACK_PROJECTION -D HS_PROFILE_PC -D HS_PROFILE_MODULATION_SOURCE=1 -D HS_PROFILE_MODULATION_MECHANIC=0 -D HS_SCAN_METRICS' \
  pio run -v -e profile 2>&1 | tee -a "$log"
PLATFORMIO_BUILD_FLAGS='-D HS_PROFILE_TARGET=MindSplatter -D HS_MINDSPLATTER_REPLAY -D HS_MINDSPLATTER_REPLAY_AB' \
  pio run -v -e profile 2>&1 | tee -a "$log"
PLATFORMIO_BUILD_FLAGS='-D HS_PROFILE_TARGET=MindSplatter -D HS_PROFILE_PLOT_STALLS' \
  pio run -v -e profile 2>&1 | tee -a "$log"
