#!/usr/bin/env bash
# Regenerate the sRGB decode tables and diff them against the committed header.
#
# usage: srgb_decode_provenance.sh [cmake configure args...]
set -euo pipefail
# shellcheck source=tools/ci/common.sh
. "$(dirname "$0")/common.sh"

cmake --preset tests -DHS_INSTALL_GIT_HOOKS=OFF -DHS_REQUIRE_GENERATORS=ON "$@"
cmake --build --preset tests --target srgb_decode_gen
mkdir -p build/provenance
generated=build/provenance/srgb_decode_lut.generated.h
./build/tests/tests/srgb_decode_gen "$generated"
require_regenerated core/color/srgb_decode_lut.h "$generated" \
  scripts/generate_srgb_decode.cpp \
  "cmake --build --preset tests --target srgb_decode_gen && ./build/tests/tests/srgb_decode_gen"
