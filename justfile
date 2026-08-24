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
# alias, which resolves and then refuses to run. HS_PYTHON overrides, matching
# the probe order in .githooks/pre-commit. Exported so the `_doxygen-theme`
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

# Native unit-test suite (Clang) + CTest.
test:
    cmake --preset tests
    cmake --build --preset tests
    ctest --preset tests

# Python, JavaScript and shell lint checks used by CI. ruff's and shellcheck's
# rule sets move between releases, so both binaries on PATH are held to the pins
# the ci.yml lint job installs; the npm linters are locked by package-lock.json.
# The shell set is the same one that job enumerates from the index -- the gate
# scripts and the hooks, where a shellcheck-class defect passes a gate silently.
# Each linter is preceded by that job's anti-vacuity probe: a linter handed
# nothing exits 0, so the selected set is asserted non-empty first. eslint
# reports per-file objects, so a report naming no file is an empty selection.
lint:
    {{py}} tools/build_pins.py --check-tool ruff
    bash -c 'tmp=$(mktemp); trap "rm -f -- \"$tmp\"" EXIT; ruff check --no-cache --show-files . > "$tmp"; test -s "$tmp" || { echo "ruff selected no files -- the ruff.toml exclude list or a .gitignore rule is broken"; exit 1; }'
    ruff check --no-cache .
    bash -c 'tmp=$(mktemp); trap "rm -f -- \"$tmp\"" EXIT; npx eslint . --format json > "$tmp" || true; grep -q filePath "$tmp" || { echo "eslint selected no files -- the eslint.config.mjs ignores list is broken"; exit 1; }'
    npm run lint
    {{py}} tools/build_pins.py --check-tool shellcheck
    bash -c "tmp=\$(mktemp); trap 'rm -f -- \"\$tmp\"' EXIT; git ls-files -- '*.sh' '.githooks/*' > \"\$tmp\"; test -s \"\$tmp\" || { echo 'no shell files selected -- the shell path list is broken'; exit 1; }; xargs shellcheck --exclude=SC1091,SC2034 < \"\$tmp\""
    bash tools/profile_sweep.sh check

# Formatting gate over the whole tracked first-party C++ set: the ci.yml
# clang-format job's invocation. Majors reflow differently, so the
# binary on PATH is held to CI's pin the way ruff is above. The exclusion regex
# is pinned against the ci.yml and .githooks/pre-commit copies by
# tools/build_pins.py --check. That job's anti-vacuity assertion comes with it:
# xargs handed an empty list runs nothing and exits 0.
clang-format:
    {{py}} tools/build_pins.py --check-tool clang-format
    bash -c "tmp=\$(mktemp); trap 'rm -f -- \"\$tmp\"' EXIT; git ls-files -- '*.h' '*.hpp' '*.cpp' '*.cc' '*.inl' | grep -vE '(^|/)core/vendor/|(^|/)core/color/color_luts\.h$|(^|/)core/color/gamut_lut\.h$|(^|/)core/spatial/reaction_graph\.cpp$|(^|/)tests/mindsplatter_replay_corpus\.h$' > \"\$tmp\" || true; test -s \"\$tmp\" || { echo 'no files selected -- the pathspec or exclusion regex is broken'; exit 1; }; xargs clang-format --dry-run --Werror --style=file < \"\$tmp\""

# Every tracked C/C++ source carries the header LICENSE grants it, plus the
# checker's own unit tests -- the ci.yml license-headers job.
license-headers:
    bash tools/require_test_files.sh "tools/license_check_tests/test*.py"
    {{py}} -m unittest discover -s tools/license_check_tests
    {{py}} tools/license_check.py

# First-party warning gate over every platformio.ini environment -- the
# ci.yml teensy-warnings job. The warning set is the pinned toolchain's, which
# the pinned PlatformIO selects. The build must be COLD (a cached TU emits no
# warnings), so the object cache and .pio/build go first and the whole firmware
# tree recompiles; budget tens of minutes. teensy_build.log is gitignored.
teensy-warnings:
    {{py}} tools/build_pins.py --check-tool platformio
    bash -c "set -o pipefail; rm -rf .pio/build_cache .pio/build && pio run -v 2>&1 | tee teensy_build.log"
    {{py}} tools/teensy_warnings.py --build-log teensy_build.log --baseline /dev/null

# The README's `tree daydream` fence draws the sibling checkout's tracked tree;
# docs_check.py can only validate it against a --checkout root (ci.yml checks the
# sibling out for exactly that). Point it at ../daydream when that checkout is
# there, so a local run gates the same claim CI does; without it the fence is
# accepted unvalidated, which the checker requires be asked for explicitly.
# path_exists resolves against the justfile directory, so the recipe works from
# any working directory.
daydream_checkout := if path_exists("../daydream") == "true" {
    "--checkout daydream=../daydream"
} else { "--skip-checkout daydream" }

# Validate tracked Markdown using the same commands as the ci.yml docs-markdown job.
docs-check:
    bash tools/require_test_files.sh 'tools/docs_check_tests/test*.py'
    bash tools/require_test_files.sh 'tools/docs_images_tests/test*.py'
    {{py}} -m unittest discover -s tools/docs_check_tests
    {{py}} -m unittest discover -s tools/docs_images_tests
    {{py}} tools/docs_check.py {{daydream_checkout}}
    {{py}} tools/build_pins.py --check
    bash tools/require_test_files.sh 'tools/build_pins_tests/test*.py'
    {{py}} -m unittest discover -s tools/build_pins_tests

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

# Teensy 4 shipping-image gates + compile profiles (CI parity for a VMicro developer).
# Needs PlatformIO (`pip install platformio`); the Teensy toolchain auto-installs
# on first `pio run`. The contract is "same PASS/FAIL under the headroom'd
# ceilings" tools/teensy_budgets.json sets, NOT byte-identity with the
# VMicro/bench image. Those ceilings were calibrated against the pinned
# PlatformIO, which selects the toolchain, so the pin is checked before building.
# The wrapper streams the pio output, then appends a combined per-env
# FLASH/RAM1/RAM2 table from the teensy_size lines.
teensy-size:
    {{py}} tools/build_pins.py --check-tool platformio
    {{py}} tools/teensy_size_table.py holosphere phantasm holosphere_dma phantasm8 profile profile_o3

# Host self-tests behind the Teensy toolchain: size/layout gate parser + layout
# invariants + warning ratchet, the PlatformIO build hook, the
# git-hook contract tests, the profile log parser, the
# relax-bake generator and the routed PCB metadata — pure Python, no ARM
# toolchain. Mirrors the ci.yml teensy-gate-tests job, including its
# non-empty discovery guards and the cross-check that every test-suite
# directory is named by the workflow.
teensy-gate-test:
    bash tools/check_test_dir_pins.sh
    bash tools/require_test_files.sh "tools/teensy_gate_tests/test*.py"
    {{py}} -m unittest discover -s tools/teensy_gate_tests -v
    bash tools/require_test_files.sh "tools/coverage_tests/test*.py"
    {{py}} -m unittest discover -s tools/coverage_tests -v
    bash tools/require_test_files.sh "tools/teensy_hook_tests/test*.py"
    {{py}} -m unittest discover -s tools/teensy_hook_tests -v
    bash tools/require_test_files.sh "tools/githook_tests/test*.py"
    {{py}} -m unittest discover -s tools/githook_tests -v
    bash tools/require_test_files.sh "tools/profile_tests/test*.py"
    {{py}} -m unittest discover -s tools/profile_tests -v
    bash tools/require_test_files.sh "tools/relax_bake_tests/test*.py"
    {{py}} -m unittest discover -s tools/relax_bake_tests -v
    bash tools/require_test_files.sh "hardware/phantasm/gen/tests/test*.py"
    {{py}} -m unittest discover -s hardware/phantasm/gen/tests -v
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
