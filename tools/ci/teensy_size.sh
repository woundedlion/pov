#!/usr/bin/env bash
# Build the shipping firmware targets through the size and layout gate, then
# compile the profiling variants.
#
# usage: teensy_size.sh <scratch-dir>   holds the profiling variants' build trees
set -euo pipefail

scratch="${1:?usage: teensy_size.sh <scratch-dir>}"
# shellcheck source=tools/profile_variants.sh
. tools/profile_variants.sh
python tools/teensy_size_table.py

profile_variant() {
  PLATFORMIO_BUILD_DIR="$scratch/$1" PLATFORMIO_BUILD_FLAGS="$2" pio run -e profile
}
for i in "${!HS_PROFILE_VARIANT_NAMES[@]}"; do
  profile_variant "${HS_PROFILE_VARIANT_NAMES[$i]}" "${HS_PROFILE_VARIANT_FLAGS[$i]}"
done
