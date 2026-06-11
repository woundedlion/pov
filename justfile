# Holosphere build tasks — thin wrappers over the canonical CMake presets
# (CMakePresets.json). Run `just` with no arguments (or `just --list`) to see
# the available recipes.
#
# The wasm recipes need the Emscripten toolchain: set the EMSDK env var (run
# emsdk_env once) before invoking them. The `test` recipe is native-only and
# does not need EMSDK.

# Show the available recipes when run with no arguments.
default:
    @just --list

# WASM release build of the simulator module (daydream).
build:
    cmake --preset wasm-release
    cmake --build --preset wasm-release

# WASM debug build (-O0 -g -sASSERTIONS, 64 KB stack).
build-debug:
    cmake --preset wasm-debug
    cmake --build --preset wasm-debug

# Native unit-test suite (Clang) + CTest.
test:
    cmake --preset tests
    cmake --build --preset tests
    ctest --preset tests

# WASM release build + install the module into ../daydream.
install:
    cmake --preset wasm-release
    cmake --build --preset wasm-release-install
