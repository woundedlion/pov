# Holosphere build tasks — thin wrappers over the canonical CMake presets
# (CMakePresets.json). Run `just` with no arguments (or `just --list`) to see
# the available recipes.
#
# The wasm recipes need the Emscripten toolchain: set the EMSDK env var (run
# emsdk_env once) before invoking them. The `test` recipe is native-only and
# does not need EMSDK.

# The [windows] recipes use cmd.exe syntax (copy /y, if not exist, parenthesized
# echo); pin the interpreter so they run under cmd regardless of just's default
# shell (a developer defaulting just to sh/pwsh would otherwise hit a syntax error).
set windows-shell := ["cmd", "/c"]

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

# Headless smoke test of the shipped WASM module (instantiates the built
# module and drives every effect; asserts arena/stack high-water marks). This
# is the CI `wasm` job's runtime gate — run it locally so `just build` is not
# shipping an un-exercised module. Builds first so it runs against fresh output.
smoke: build
    node scripts/wasm_smoke.mjs

# Capture the WebGL effect gallery to docs/screenshots/ (Playwright, headless).
# Needs the sibling daydream checkout served (see README) and the chromium
# browser installed once via `npx playwright install chromium`.
screenshots:
    node scripts/capture_screenshots.mjs

# Native unit-test suite (Clang) + CTest.
test:
    cmake --preset tests
    cmake --build --preset tests
    ctest --preset tests

# Validate tracked Markdown using the same commands as docs CI.
docs-check:
    python -m unittest discover -s tools/docs_check_tests
    python tools/docs_check.py

# Build Doxygen API reference locally into build/docs/html/.
# Clones doxygen-awesome theme into .doxygen-awesome/ on first run and
# synthesizes the gitignored Doxyfile.local (Doxyfile + theme overrides, mirroring
# .github/workflows/docs.yml). Requires doxygen on PATH.
docs: docs-check _doxygen-theme _doxyfile-local
    cmake -E make_directory build/docs
    doxygen Doxyfile.local

# Fetch the exact doxygen-awesome revision used by CI. The clone guard is split
# per-OS; the fetch and checkout also refresh existing clones.
[unix]
_doxygen-theme:
    test -d .doxygen-awesome/.git || git clone --filter=blob:none --no-checkout https://github.com/jothepro/doxygen-awesome-css.git .doxygen-awesome
    git -C .doxygen-awesome fetch --depth 1 origin 568f56cde6ac78b6dfcc14acd380b2e745c301ea
    git -C .doxygen-awesome checkout --detach 568f56cde6ac78b6dfcc14acd380b2e745c301ea

[windows]
_doxygen-theme:
    if not exist .doxygen-awesome\.git git clone --filter=blob:none --no-checkout https://github.com/jothepro/doxygen-awesome-css.git .doxygen-awesome
    git -C .doxygen-awesome fetch --depth 1 origin 568f56cde6ac78b6dfcc14acd380b2e745c301ea
    git -C .doxygen-awesome checkout --detach 568f56cde6ac78b6dfcc14acd380b2e745c301ea

# Synthesize Doxyfile.local = Doxyfile + docs/doxygen-theme.cfg (the same theme
# overrides docs.yml appends). The copy+append is shell-specific, so it's split
# per-OS; the appended content is shared, not duplicated.
[unix]
_doxyfile-local:
    cp Doxyfile Doxyfile.local
    cat docs/doxygen-theme.cfg >> Doxyfile.local

[windows]
_doxyfile-local:
    copy /y Doxyfile Doxyfile.local
    type docs\doxygen-theme.cfg >> Doxyfile.local

# WASM release build + install the module into ../daydream.
install:
    cmake --preset wasm-release
    cmake --build --preset wasm-release-install

# Teensy 4 shipping-image gates + compile profiles (CI parity for a VMicro developer).
# Needs PlatformIO (`pip install platformio`); the Teensy toolchain auto-installs
# on first `pio run`. The contract is "same PASS/FAIL under the headroom'd
# ceilings", NOT byte-identity with the VMicro/bench image (docs/teensy_ci_gate_spec.md §11).
teensy-size:
    pio run -e holosphere -e phantasm -e holosphere_dma -e phantasm8

# Host self-tests for the size/layout gate parser + layout invariants + warning
# ratchet — pure Python, no ARM toolchain (spec §9.1). Mirrors the CI job.
teensy-gate-test:
    python -m unittest discover -s tools/teensy_gate_tests

# Profile one effect on an attached Teensy: build the single-effect profiling
# image (Phantasm shipping flags + HS_PROFILE cycle counters, board = segment 0
# of 4), flash it, then capture the serial readout for `seconds` into
# build/profile_capture.log. Pass any roster effect class name.
[windows]
profile effect="DisplacementField" seconds="150":
    set "PLATFORMIO_BUILD_FLAGS=-D HS_PROFILE_TARGET={{effect}}" && pio run -e profile -t upload
    python tools/profile_capture.py --seconds {{seconds}}

[unix]
profile effect="DisplacementField" seconds="150":
    PLATFORMIO_BUILD_FLAGS="-D HS_PROFILE_TARGET={{effect}}" pio run -e profile -t upload
    python tools/profile_capture.py --seconds {{seconds}}

# Regenerate the PHANTASM PCB outputs into hardware/phantasm/gen/out/ (all
# gitignored) from the COMMITTED board. It never re-runs the schematic/PCB
# generators, which would discard the routing + silk; needs kicad-cli on PATH
# (or set KICAD_CLI to its full path).
# Outputs: Gerbers + Excellon drill, JLCPCB upload zip, assembly BOM + CPL, DRC.
pcb:
    python hardware/phantasm/gen/fab.py
