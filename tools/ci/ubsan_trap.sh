#!/usr/bin/env bash
# Build and run the UBSan trap-mode flavour on Windows, where emsdk provides no
# MSVC ASan runtime.
set -euo pipefail

cmake --preset tests -B build/tests-ubsan -DHS_INSTALL_GIT_HOOKS=OFF -DHS_SANITIZE=ON
cmake --build build/tests-ubsan --target mindsplatter_replay
ctest --test-dir build/tests-ubsan --output-on-failure --no-tests=error -R "^unit_mindsplatter_replay$"
