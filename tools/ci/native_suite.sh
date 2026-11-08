#!/usr/bin/env bash
# Configure, build and run the native suite in build/tests at the smoke depth
# every CI leg uses.
#
# usage: native_suite.sh [--exclude-generated] [cmake configure args...]
set -euo pipefail
# shellcheck source=tools/ci/common.sh
. "$(dirname "$0")/common.sh"

exclude=()
if [ "${1:-}" = --exclude-generated ]; then
  exclude=(-E "$GENERATED_TESTS")
  shift
fi

cmake --preset tests -DHS_INSTALL_GIT_HOOKS=OFF -DHS_REQUIRE_GENERATORS=ON "$@"
cmake --build --preset tests
cmake --build --preset tests --target excluded_targets
HS_SMOKE_FRAMES=120 HS_EFFECTS_FULL=0 HS_REQUIRE_EFFECTS_FULL=0 HS_SKIPS_ARE_ERRORS=1 \
  HS_BUFFER_FREE_WATCHDOG_US=30000000 \
  ctest --preset tests --output-on-failure --no-tests=error "${exclude[@]}"
