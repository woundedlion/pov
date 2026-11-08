#!/usr/bin/env bash
# Configure, build and run one optimized native flavour.
#
# usage: optimized.sh <ieee|fastmath|ndebug-death> [cmake configure args...]
# HS_CI_BUILD_DIR overrides the build directory (default build/tests).
set -euo pipefail
# shellcheck source=tools/ci/common.sh
. "$(dirname "$0")/common.sh"

flavor="${1:?usage: optimized.sh <ieee|fastmath|ndebug-death> [cmake args...]}"
shift
include=".*"
exclude="$GENERATED_TESTS"
skips_are_errors=1
float_flags=""
case "$flavor" in
  ieee)
    flags="-O2 -g" ;;
  fastmath)
    flags="-Os -g"
    float_flags="-ffast-math -fno-finite-math-only -DHS_TEST_FAST_MATH=1"
    exclude='^unit_(3dmath|sdf|effects|mindsplatter|shapeshifter_tiles|conway_continuity|partition_seam|mesh|shader_chain|pullback_capture_determinism|pov_sync|color_luts|gamut_lut|tangent_field_constants)$'
    skips_are_errors=0 ;;
  ndebug-death)
    flags="-O2 -g -DNDEBUG"
    include='^unit_death_([0-9]+|census)$' ;;
  *)
    echo "optimized.sh: unknown flavor $flavor" >&2
    exit 2 ;;
esac

dir="${HS_CI_BUILD_DIR:-build/tests}"
cmake --preset tests -B "$dir" -DHS_INSTALL_GIT_HOOKS=OFF -DHS_REQUIRE_GENERATORS=ON \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo "-DCMAKE_CXX_FLAGS_RELWITHDEBINFO=$flags $float_flags" "$@"
cmake --build "$dir"
HS_SMOKE_FRAMES=120 HS_EFFECTS_FULL=1 HS_REQUIRE_EFFECTS_FULL=1 \
  HS_SKIPS_ARE_ERRORS="$skips_are_errors" HS_BUFFER_FREE_WATCHDOG_US=30000000 \
  ctest --test-dir "$dir" --parallel 4 --output-on-failure --no-tests=error \
  -R "$include" -E "$exclude"
