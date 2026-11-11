#!/usr/bin/env bash
# Build the shipping firmware targets through the size and layout gate, and the
# other environments for compile coverage.
#
# usage: teensy_size.sh [<scratch-dir>]   also compile the profiling variants,
#                                         each in its own build tree there
set -euo pipefail

python tools/teensy_size_table.py

[ "$#" -gt 0 ] || exit 0
scratch=$1
# shellcheck source=tools/profile_variants.sh
. tools/profile_variants.sh
profile_variant() {
  PLATFORMIO_BUILD_DIR="$scratch/$1" PLATFORMIO_BUILD_FLAGS="$2" pio run -e profile
}
for i in "${!HS_PROFILE_VARIANT_NAMES[@]}"; do
  profile_variant "${HS_PROFILE_VARIANT_NAMES[$i]}" "${HS_PROFILE_VARIANT_FLAGS[$i]}"
done
