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

# Python interpreter for every recipe that runs one. Stock Linux/macOS ship
# `python3` only; on Windows a `python3` on PATH is usually the Store execution
# alias, which resolves and then refuses to run. HS_PYTHON overrides, the same
# override the .githooks probe loops honour. Exported so the `_doxygen-theme`
# parameter default can reach it from inside a backtick.
export py := env_var_or_default("HS_PYTHON", if os_family() == "windows" { "python" } else { "python3" })

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
# Node is held to the pin that job's setup-node installs.
smoke: build
    {{py}} tools/build_pins.py --check-tool node
    node scripts/wasm_smoke.mjs

# Capture the WebGL effect gallery to docs/screenshots/ (Playwright, headless).
# Needs the sibling daydream checkout served (see README) and the chromium
# browser installed once via `npx playwright install chromium`.
screenshots:
    node scripts/capture_screenshots.mjs

# Native unit-test suite (Clang) + CTest at the smoke window every CI leg drives.
# The 8-frame default arms no preset transition, so the pause, slot-reuse and
# FIFO-expiry paths never run; pass a narrower window for a fast iteration loop.
test $HS_SMOKE_FRAMES="120" $HS_SKIPS_ARE_ERRORS="1":
    cmake --preset tests
    cmake --build --preset tests
    cmake --build --preset tests --target excluded_targets
    ctest --preset tests

# Python, JavaScript and shell lint checks used by CI. ruff's and shellcheck's
# rule sets move between releases, so both binaries on PATH are held to the pins
# the ci.yml lint job installs; the npm linters are locked by package-lock.json.
# The shell set is the same one that job enumerates from the index. Each linter
# is preceded by that job's anti-vacuity probe. The line-ending check runs first
# for the reason it does in that job: a working copy that diverged from its
# eol=lf blob is what the linters below would otherwise read.
# Normalize CRLF working copies without changing the index or discarding edits.
normalize-eol:
    bash tools/eol_gate.sh --fix-worktree

lint:
    bash tools/eol_gate.sh
    {{py}} tools/build_pins.py --check-tool ruff
    bash tools/ruff_selection_guard.sh
    ruff check --no-cache .
    bash tools/eslint_selection_guard.sh
    npm run lint
    {{py}} tools/build_pins.py --check-tool shellcheck
    bash tools/shellcheck_gate.sh
    bash tools/profile_sweep.sh check

# Formatting gate over the whole tracked first-party C++ set: the ci.yml
# clang-format job's invocation. Majors reflow differently, so the
# binary on PATH is held to CI's pin the way ruff is above.
clang-format:
    {{py}} tools/build_pins.py --check-tool clang-format
    bash tools/clang_format_gate.sh

# Every tracked C/C++ source carries the header LICENSE grants it, plus the
# checker's own unit tests -- the ci.yml license-headers job.
license-headers:
    {{py}} tools/license_check.py

# The committed gamut LUT matches what the generator emits, plus the generator's
# own unit tests -- the ci.yml gamut-lut-provenance job. The solve runs 1-2
# minutes. numpy decides the emitted bytes, so the module the interpreter imports
# is held to the pin that job installs, the way ruff is above.
gamut-lut:
    {{py}} tools/build_pins.py --check-tool numpy
    {{py}} tools/gen_gamut_lut.py --check

# First-party warning gate over every platformio.ini environment -- the
# ci.yml teensy-warnings job. The warning set is the pinned toolchain's, which
# the pinned PlatformIO selects. The build is cold, so budget tens of minutes;
# teensy_build.log is gitignored.
teensy-warnings:
    {{py}} tools/build_pins.py --check-tool platformio
    bash tools/teensy_cold_build.sh teensy_build.log
    {{py}} tools/teensy_warnings.py --build-log teensy_build.log

# Regenerate documentation maps and derived reference data.
docs-sync:
    {{py}} tools/docs_check.py --sync --auto-checkout

# Validate tracked Markdown using the same commands as the ci.yml docs-markdown
# job, plus the docs-images job's checker: this recipe runs that checker's unit
# tests, which say nothing about the tracked tree on their own.
docs-check:
    {{py}} tools/docs_check.py --auto-checkout
    {{py}} tools/docs_images.py
    {{py}} tools/build_pins.py --check

# Build Doxygen API reference locally into build/docs/html/.
# Clones doxygen-awesome theme into .doxygen-awesome/ on first run and
# synthesizes the gitignored Doxyfile.local (Doxyfile + theme overrides, mirroring
# .github/workflows/docs.yml). Requires doxygen on PATH at the pinned version:
# warning text and generated markup move between releases.
docs: docs-check _doxygen-theme _doxyfile-local
    {{py}} tools/build_pins.py --check-tool doxygen
    cmake -E make_directory build/docs
    doxygen Doxyfile.local
    {{py}} tools/docs_images.py --stage

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
_doxygen-theme sha=`%py% tools/build_pins.py doxygen-awesome`:
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

# WASM release build + install the module into ../daydream. Gated on `smoke`
# (which builds first), so the module and provenance triple written into the
# simulator tree are always the ones the runtime gate just exercised.
install: smoke
    cmake --build --preset wasm-release-install
    node scripts/wasm_smoke.mjs ../daydream/holosphere_wasm.js

# Build and flash the stationary bench test image (targets/Bench) to an attached
# Teensy. Every board on the rig takes the same image, and it reads with the
# sphere at rest: one colour across the whole canvas, holding on red, green,
# blue and white. Flash `phantasm` to put the show back.
bench:
    pio run -e bench -t upload

# Teensy 4 shipping-image gates + compile profiles (CI parity for a VMicro developer).
# Needs PlatformIO (`pip install platformio`); the Teensy toolchain auto-installs
# on first `pio run`. The contract is "same PASS/FAIL under the headroom'd
# ceilings" tools/teensy_budgets.json sets, NOT byte-identity with the
# VMicro/bench image. Those ceilings were calibrated against the pinned
# PlatformIO, which selects the toolchain, so the pin is checked before building.
# The wrapper builds every platformio.ini environment (a bare `pio run`, the set
# the warning ratchet expects too), streams the pio output, then appends a
# combined per-env FLASH/RAM1/RAM2 table from the teensy_size lines.
#
# The last line is the size trail's producer: it parses the ELFs this build just
# linked into the worktree's pending record, which the post-commit hook stamps
# onto the next commit. Error-suppressed (`-`) on purpose — a missing ELF, no
# python or no git repo leaves the trail alone instead of failing the build.
teensy-size:
    {{py}} tools/build_pins.py --check-tool platformio
    {{py}} tools/teensy_size_table.py
    -{{py}} tools/teensy_size_trail.py record

# All tracked Python unit suites.
python-test:
    {{py}} tools/run_python_tests.py

# Python tests and routed PCB metadata.
teensy-gate-test: python-test
    {{py}} hardware/phantasm/gen/board_metadata.py --check

# Profile one effect on an attached Teensy: build the single-effect profiling
# image (Phantasm shipping flags + HS_PROFILE cycle counters, board = segment 0
# of 4), flash it, then capture the serial readout for `seconds` into
# build/prof/<effect>_ship.log. Pass any roster effect class name.
#
# Delegates to profile_one.sh so every device path runs under the one host-
# global device lock (tools/device_lock.sh) and the same header/stale-build
# checks; flashing around it would clobber a concurrent agent's capture.
#
# deep=1 turns on the HS_PROFILE_DEEP sub-scopes (per-pixel/per-cell counters in
# shared render code) and captures to build/prof/<effect>_ship_deep.log instead,
# leaving the roster log untouched.
profile effect="DisplacementField" seconds="150" deep="0":
    HS_PROFILE_DEEP="{{deep}}" bash tools/profile_one.sh "{{effect}}" profile "{{seconds}}" 32

# Repeated physics-free render of the frozen production-resolution corpus.
profile-mindsplatter-replay env="profile" seconds="150":
    bash tools/profile_one.sh MindSplatter "{{env}}" "{{seconds}}" 32 -D HS_MINDSPLATTER_REPLAY

# Same-device candidate/reference visual comparison; timing includes both.
profile-mindsplatter-replay-ab env="profile" seconds="150":
    bash tools/profile_one.sh MindSplatter "{{env}}" "{{seconds}}" 32 -D HS_MINDSPLATTER_REPLAY -D HS_MINDSPLATTER_REPLAY_AB

# Regenerate the PHANTASM PCB outputs into hardware/phantasm/gen/out/ (all
# gitignored) from the COMMITTED board. It never re-runs the schematic/PCB
# generators, which would discard the routing + silk; needs kicad-cli on PATH
# (or set KICAD_CLI to its full path).
# Outputs: Gerbers + Excellon drill, JLCPCB upload zip, assembly BOM + CPL, DRC.
pcb:
    {{py}} hardware/phantasm/gen/fab.py
