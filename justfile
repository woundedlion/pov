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

# Build Doxygen API reference locally into build/docs/html/.
# Clones doxygen-awesome theme into .doxygen-awesome/ on first run.
# Requires doxygen on PATH.
docs: _doxygen-theme
    cmake -E make_directory build/docs
    doxygen Doxyfile.local

# Clone the doxygen-awesome theme on first run. The clone guard is the one
# shell-specific step, so it's split per-OS: the rest of `docs` (cmake/doxygen)
# is shell-agnostic. Without this the cmd.exe `if not exist` form made `just
# docs` a syntax error under sh on the Linux/macOS CI the file otherwise targets.
[unix]
_doxygen-theme:
    test -d .doxygen-awesome || git clone --depth 1 --branch v2.3.4 https://github.com/jothepro/doxygen-awesome-css.git .doxygen-awesome

[windows]
_doxygen-theme:
    if not exist .doxygen-awesome git clone --depth 1 --branch v2.3.4 https://github.com/jothepro/doxygen-awesome-css.git .doxygen-awesome

# WASM release build + install the module into ../daydream.
install:
    cmake --preset wasm-release
    cmake --build --preset wasm-release-install
