#!/usr/bin/env bash
# Regenerate the sRGB decode tables and diff them against the committed header.
#
# usage: srgb_decode_provenance.sh [cmake configure args...]
set -euo pipefail

cmake --preset tests -DHS_INSTALL_GIT_HOOKS=OFF -DHS_REQUIRE_GENERATORS=ON "$@"
cmake --build --preset tests --target srgb_decode_gen
mkdir -p build/provenance
generated=build/provenance/srgb_decode_lut.generated.h
./build/tests/tests/srgb_decode_gen "$generated"
if ! diff -u core/color/srgb_decode_lut.h "$generated" > build/provenance/srgb_decode_lut.diff; then
  head -40 build/provenance/srgb_decode_lut.diff
  echo "::error::core/color/srgb_decode_lut.h is out of sync with scripts/generate_srgb_decode.cpp"
  echo "Regenerate with: cmake --build --preset tests --target srgb_decode_gen && ./build/tests/tests/srgb_decode_gen"
  exit 1
fi
echo "Regenerated tables match core/color/srgb_decode_lut.h in full."
