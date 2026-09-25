# Holosphere build tasks — thin wrappers over the canonical CMake presets
# (CMakePresets.json). Run `just` with no arguments (or `just --list`) to see
# the available recipes.
#
# The wasm recipes need the Emscripten toolchain: set the EMSDK env var (run
# emsdk_env once) before invoking them. The `test` recipe is native-only and
# does not need EMSDK.

# The [windows] recipes use cmd.exe syntax (copy /y, if not exist); pin the interpreter so they run under cmd regardless of just's default
# shell (a developer defaulting just to sh/pwsh would otherwise hit a syntax error).
set windows-shell := ["cmd", "/c"]

# Python interpreter for every recipe that runs one. Stock Linux/macOS ship
# `python3` only; on Windows a `python3` on PATH is usually the Store execution
# alias, which resolves and then refuses to run. HS_PYTHON overrides, the same
# override the .githooks probe loops honour. Exported so the `_doxygen-theme`
# parameter default can reach it from inside a backtick.
export py := env_var_or_default("HS_PYTHON", if os_family() == "windows" { "python" } else { "python3" })
export quoted_py := '"' + py + '"'
python_command := if os_family() == "windows" { "%quoted_py%" } else { quoted_py }

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

# Build and smoke-test the WASM engine.
smoke: build
    {{python_command}} tools/build_pins.py --check-tool node
    node scripts/wasm_smoke.mjs

# Capture the effect gallery with headless Chromium.
screenshots:
    node scripts/capture_screenshots.mjs

# Build and run the native suite over the configured smoke window.
test $HS_SMOKE_FRAMES="120" $HS_SKIPS_ARE_ERRORS="1":
    cmake --preset tests
    cmake --build --preset tests
    cmake --build --preset tests --target excluded_targets
    ctest --preset tests

# Normalize CRLF working copies without changing the index or discarding edits.
normalize-eol:
    bash tools/eol_gate.sh --fix-worktree

# Run the shared local and CI lint checks.
lint:
    {{python_command}} tools/build_pins.py --check-tool just
    bash tools/whitespace_gate.sh
    bash tools/eol_gate.sh
    {{python_command}} tools/build_pins.py --check-tool ruff
    bash tools/ruff_selection_guard.sh
    ruff check --no-cache .
    bash tools/eslint_selection_guard.sh
    npm run lint
    {{python_command}} tools/build_pins.py --check-tool shellcheck
    bash tools/shellcheck_gate.sh
    bash tools/profile_sweep.sh check
    {{python_command}} tools/build_pins.py --check-tool actionlint
    actionlint -shellcheck shellcheck

# Check formatting of all tracked first-party C++ sources.
clang-format:
    {{python_command}} tools/build_pins.py --check-tool clang-format
    bash tools/clang_format_gate.sh

# Check license headers on tracked C/C++ sources; unit tests run via python-test.
license-headers:
    {{python_command}} tools/license_check.py

# Check the committed gamut LUT against its pinned generator.
gamut-lut:
    {{python_command}} tools/build_pins.py --check-tool numpy
    {{python_command}} tools/gen_gamut_lut.py --check

# Build every firmware environment and check first-party warnings.
teensy-warnings:
    {{python_command}} tools/build_pins.py --check-tool platformio
    bash tools/teensy_cold_build.sh teensy_build.log
    {{python_command}} tools/teensy_warnings.py --build-log teensy_build.log

# Regenerate documentation maps and derived reference data.
docs-sync:
    {{python_command}} tools/docs_check.py --sync --auto-checkout

# Validate tracked Markdown, image references, and build pins.
docs-check:
    {{python_command}} tools/docs_check.py --auto-checkout
    {{python_command}} tools/docs_images.py
    {{python_command}} tools/build_pins.py --check

# Build the themed Doxygen API reference.
docs: docs-check _doxygen-theme _doxyfile-local
    {{python_command}} tools/build_pins.py --check-tool doxygen
    cmake -E make_directory build/docs
    doxygen Doxyfile.local
    {{python_command}} tools/docs_images.py --stage

# Fetch the exact doxygen-awesome revision used by CI. The clone guard is split
# per-OS; the fetch and checkout also refresh existing clones. The pin is a
# parameter default, not a justfile-level assignment: only that form defers the
# backtick to this recipe, leaving every python-free recipe runnable without it.
[unix]
_doxygen-theme sha=`"$py" tools/build_pins.py doxygen-awesome`:
    test -d .doxygen-awesome/.git || git clone --filter=blob:none --no-checkout https://github.com/jothepro/doxygen-awesome-css.git .doxygen-awesome
    git -C .doxygen-awesome fetch --depth 1 origin {{sha}}
    git -C .doxygen-awesome checkout --detach {{sha}}

[windows]
_doxygen-theme sha=`%quoted_py% tools/build_pins.py doxygen-awesome`:
    if not exist .doxygen-awesome\.git git clone --filter=blob:none --no-checkout https://github.com/jothepro/doxygen-awesome-css.git .doxygen-awesome
    git -C .doxygen-awesome fetch --depth 1 origin {{sha}}
    git -C .doxygen-awesome checkout --detach {{sha}}

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

# Build, smoke-test, and install WASM into ../daydream.
install: smoke
    cmake --build --preset wasm-release-install
    node scripts/wasm_smoke.mjs ../daydream/generated/holosphere_wasm.js

# Windows only: build and flash the bench image under the per-board device lock.
bench:
    bash tools/upload_one.sh bench

# Build all firmware environments and check size budgets.
teensy-size:
    {{python_command}} tools/build_pins.py --check-tool platformio
    {{python_command}} tools/teensy_size_table.py
    -{{python_command}} tools/teensy_size_trail.py record

# All tracked Python unit suites.
python-test:
    {{python_command}} tools/run_python_tests.py

# Python tests and routed PCB metadata.
teensy-gate-test: python-test
    {{python_command}} hardware/phantasm/gen/board_metadata.py --check

# Windows only: build, flash, and capture one effect under the device lock.
profile effect="DisplacementField" seconds="150" $HS_PROFILE_DEEP="0":
    bash tools/profile_one.sh "{{effect}}" profile "{{seconds}}" 32

# Windows only: repeated render of the frozen production-resolution corpus.
profile-mindsplatter-replay env="profile" seconds="150":
    bash tools/profile_one.sh MindSplatter "{{env}}" "{{seconds}}" 32 -D HS_MINDSPLATTER_REPLAY

# Windows only: same-device candidate/reference visual comparison.
profile-mindsplatter-replay-ab env="profile" seconds="150":
    bash tools/profile_one.sh MindSplatter "{{env}}" "{{seconds}}" 32 -D HS_MINDSPLATTER_REPLAY -D HS_MINDSPLATTER_REPLAY_AB

# Export fabrication files from the committed routed board.
pcb:
    {{python_command}} hardware/phantasm/gen/fab.py

# Validate generated color lookup tables.
color-lut-check:
    cmake -E make_directory build/provenance
    {{python_command}} scripts/generate_luts.py -o build/provenance/color_luts.h
    cmake -E compare_files core/color/color_luts.h build/provenance/color_luts.h

# Validate the generated reaction graph.
reaction-graph-check:
    cmake -E make_directory build/provenance
    {{python_command}} scripts/generate_reaction_graph.py -o build/provenance/reaction_graph.cpp
    cmake -E compare_files core/spatial/reaction_graph.cpp build/provenance/reaction_graph.cpp

# Validate the generated sRGB decoder tables.
srgb-decode-check:
    cmake --preset tests
    cmake --build --preset tests --target srgb_decode_gen
    cmake -E make_directory build/provenance
    cmake -E chdir build/tests/tests ./srgb_decode_gen ../../provenance/srgb_decode_lut.h
    cmake -E compare_files core/color/srgb_decode_lut.h build/provenance/srgb_decode_lut.h

# Validate promoted shader document provenance.
patterns-check:
    node scripts/generate_promoted_shader_documents.mjs --check

# Run the Node tooling suite.
scripts-test:
    npm test

# Validate the committed screenshot gallery.
gallery-check:
    node scripts/check_screenshots.mjs

# Validate the profile archive against the effect roster.
profiles-check:
    node scripts/check_profiles.mjs
