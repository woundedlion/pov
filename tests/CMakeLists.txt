# Holosphere native unit tests.
#
# Canonical flow is the `tests` CMake preset (see CMakePresets.json), which uses
# cmake/toolchain-native-clang.cmake to locate Clang and, on Windows, the SDK
# rc.exe automatically — no Visual Studio Developer Prompt or manual flags:
#
#   cmake --preset tests
#   cmake --build --preset tests
#   ctest --preset tests
#
# The `just test` recipe is a thin wrapper over exactly those three commands.
#
# The root CMakeLists.txt (a WASM/Emscripten build) does add_subdirectory(tests);
# the unit tests are native-only, so we no-op under Emscripten to leave the WASM
# build untouched. (A standalone `cmake -S tests -B build_test` configure also
# works for CI, but the preset is the supported path.)
if(EMSCRIPTEN)
  return()
endif()

# Matches the root floor: the native-clang toolchain sets CMAKE_LINKER_TYPE
# (added in CMake 3.29), which is silently ignored on older CMake. This applies
# to the standalone `cmake -S tests` path too, so gate it here as well.
cmake_minimum_required(VERSION 3.29)

# Allow standalone configuration (cmake -S tests) as well as being pulled in by
# a native root build.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(holosphere_tests CXX)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Clang < 19 defaults to the pre-P0522 template-template matching rules and
# rejects an effect with a defaulted extra template parameter (e.g.
# IslamicStars<W, H, SegueT = Crossfade>) as an argument for the harness's
# `template <int, int> class E` parameters. GCC and Clang >= 19 apply the
# C++17 relaxed matching by default; opt older Clang (the CI / deploy-gate
# apt clang-18) in explicitly. Clang 20 removed the flag, so keep the gate.
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 19)
  add_compile_options(-frelaxed-template-template-args)
endif()

set(HS_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/..)

# The native suite is the one build where diagnostics can run: the device/WASM
# build compiles the arena guards out, so ASan/UBSan here are the only automated
# cover for the arena-lifetime hazards the death harness can't reach. -Werror is
# on so suite hygiene is self-enforcing — the fail-fast doctrine applied to the
# build itself; a new warning is a hard failure, not noise scrolled past. Every
# executable defined below goes through this, including the budget gates and the
# non-CTest compile guards; per-target optimization/math flags stay per-target.
function(_hs_warnings target)
  target_compile_options(${target} PRIVATE -Wall -Wextra -Werror)
endfunction()

add_executable(run_tests
  ${CMAKE_CURRENT_SOURCE_DIR}/run_tests.cpp
  ${HS_ROOT}/core/engine/memory.cpp
  ${HS_ROOT}/core/engine/reaction_graph.cpp)

# HS_TEST_BUILD widens the test-only inline/arena budgets (the 64-bit host
# inflates every embedded pointer past the 32-bit device footprints); the device
# footprint is unchanged. See core/animation/animation.h and core/engine/memory.h.
#
# HS_SCAN_METRICS compiles the per-pixel g_scan_metrics counters back in (off in
# the device/WASM release build, where they would be a per-pixel global store on
# the hot loop). The Face distance-oracle test reads exact_hits to tell culled
# samples from evaluated ones, so the suite needs them live. See
# core/engine/platform.h (HS_SCAN_METRIC) and core/render/sdf.h Face::distance.
target_compile_definitions(run_tests PRIVATE HS_TEST_BUILD HS_SCAN_METRICS)

# Headers are included as "tests/..." and "core/..." from the repo root.
target_include_directories(run_tests PRIVATE ${HS_ROOT} ${HS_ROOT}/core)

_hs_warnings(run_tests)

# Opt-in sanitizers. Default OFF so the pre-commit gate and the default
# `cmake --preset tests` stay fast (ASan ~2x); enable with -DHS_SANITIZE=ON.
# The "sanitizers" CI job (.github/workflows/ci.yml) builds with it on every
# push and pull request. The flavour depends on the runtime the toolchain
# ships:
#   * Linux/full-LLVM: ASan + UBSan with diagnostics (-fno-sanitize-recover
#     makes the first error a hard failure rather than logged-and-continue).
#   * Windows emsdk clang: ships no ASan runtime for the MSVC target, so use
#     runtime-free UBSan in trap mode — a UB hit becomes an immediate trap,
#     consistent with the engine's always-on fail-fast doctrine.
option(HS_SANITIZE "Build the unit tests with sanitizers (ASan+UBSan, or UBSan-trap on Windows)" OFF)

# Applies the sanitizer flavour to one target; a no-op when HS_SANITIZE is OFF.
# Every CTest-registered executable the sanitizers job runs must call this, or it
# rides that job unsanitized. The stack_measure/arena_measure budget gates are
# the deliberate exception: ASan redzones and instrumentation inflate the frame
# sizes and codegen they measure against the device build.
function(_hs_sanitize target)
  if(NOT HS_SANITIZE)
    return()
  endif()
  # Suppress only signed-integer-overflow inside vendored FastNoiseLite (its
  # hashing wraps signed ints by design); first-party code keeps full coverage.
  set(_ignore "-fsanitize-ignorelist=${CMAKE_CURRENT_SOURCE_DIR}/ubsan-ignorelist.txt")
  if(WIN32)
    target_compile_options(${target} PRIVATE
      -fsanitize=undefined -fsanitize-trap=undefined ${_ignore}
      -fno-omit-frame-pointer -g)
  else()
    target_compile_options(${target} PRIVATE
      -fsanitize=address,undefined -fno-sanitize-recover=all ${_ignore}
      -fno-omit-frame-pointer -g)
    target_link_options(${target} PRIVATE -fsanitize=address,undefined)
  endif()
endfunction()

if(HS_SANITIZE)
  if(WIN32)
    message(STATUS "Unit tests: UBSan (trap mode) enabled — no ASan runtime on emsdk clang")
  else()
    message(STATUS "Unit tests: ASan + UBSan enabled (HS_SANITIZE=ON)")
  endif()
endif()

_hs_sanitize(run_tests)

# Opt-in Clang source-based coverage for the main suite binary. Default OFF: the
# instrumentation changes codegen and slows the run, so no gating job may carry
# it. The "code-coverage" CI job turns it on and publishes an llvm-cov summary —
# advisory only, there is no threshold gate. Needs clang_rt.profile, which the
# Windows emsdk clang does not ship, so this is a Linux-toolchain option.
option(HS_COVERAGE "Instrument run_tests for llvm-cov source-based coverage" OFF)
if(HS_COVERAGE)
  target_compile_options(run_tests PRIVATE -fprofile-instr-generate -fcoverage-mapping)
  target_link_options(run_tests PRIVATE -fprofile-instr-generate)
endif()

enable_testing()

# One CTest test per module so CI shards (and local `ctest -R <name>`) can
# target subsets without re-running everything.
set(_hs_test_modules
    3dmath concepts memory spatial scb sdf conway conway_morph conway_continuity
    partition_seam
    conway_soak opchain_probe opchain_arena_survey hankin geometry spherical_field
    mesh solids
    reaction_graph color palettes easing_waves platform filter plot_scan canvas
    scan mesh_raster transformers noise generators animation effects
    shapeshifter_oracle dma_core hd107s
    dma_controller
    pov_segmented pov_single pov_sync param_marshal wasm_predicates util led
    presets styles shading death)
# The suite has a real hang mode (a Canvas constructor that spin-waits for a
# buffer that never frees — see test_effects.h), and CTest's 1500 s default would
# burn that once per shard. The whole roster runs in ~2 min unsanitized, but
# ASan + UBSan cost the slowest module (effects) ~2.7x — 227 s plain against
# 620 s sanitized — so a single budget for both builds is either useless plain
# or red sanitized.
set(_hs_test_timeout 600)
if(HS_SANITIZE)
  set(_hs_test_timeout 1800)
endif()
foreach(_mod IN LISTS _hs_test_modules)
    add_test(NAME unit_${_mod} COMMAND run_tests ${_mod})
    set_tests_properties(unit_${_mod} PROPERTIES TIMEOUT ${_hs_test_timeout})
endforeach()

# Pin _hs_test_modules against run_tests' HS_TEST_MODULE_LIST: run_tests
# --check-modules exits non-zero if the CMake list and the roster diverge in
# either direction, so a module added to one but not the other can never run
# silently in CI sharding.
add_test(NAME unit_module_roster
    COMMAND run_tests --check-modules ${_hs_test_modules})

# Pin run_tests.cpp's per-module include block against the roster, catching an
# include left behind after a module is removed (which compiles silently), and
# every tests/test_*.h on disk against that include block, catching a test file
# nothing ever compiles.
add_test(NAME unit_module_includes COMMAND ${CMAKE_COMMAND}
    -DSRC=${CMAKE_CURRENT_SOURCE_DIR}/run_tests.cpp
    -DTESTS_DIR=${CMAKE_CURRENT_SOURCE_DIR}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_includes.cmake)

# The three generated-artifact gates below re-run a Python generator and diff it
# against the committed header. Default OFF keeps a bare local checkout green:
# a missing python3/numpy exits 77 and CTest reports the gate as skipped. That
# skip is invisible in a green run, so the CI job that provisions the
# interpreter passes -DHS_REQUIRE_GENERATORS=ON, which turns the same missing
# python3/numpy into a hard failure — a broken provisioning step can no longer
# silently retire all three gates.
option(HS_REQUIRE_GENERATORS
    "Fail (rather than skip) the generated-artifact gates when python3/numpy is missing" OFF)

# Pin core/color/color_luts.h against its generator.
find_program(HS_PYTHON NAMES python3 python)
add_test(NAME unit_color_luts COMMAND ${CMAKE_COMMAND}
    -DPYTHON_EXE=${HS_PYTHON}
    -DGENERATOR=${HS_ROOT}/scripts/generate_luts.py
    -DCOMMITTED=${HS_ROOT}/core/color/color_luts.h
    -DGENERATED=${CMAKE_CURRENT_BINARY_DIR}/color_luts.generated.h
    -DSKIP_CODE=77 -DREQUIRE_PYTHON=${HS_REQUIRE_GENERATORS}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_color_luts.cmake)
set_tests_properties(unit_color_luts PROPERTIES SKIP_RETURN_CODE 77)

# Run the generator's own --check self-validation (monotonicity + sRGB
# round-trip), which CI ran only as a shell step.
add_test(NAME unit_color_luts_check COMMAND ${CMAKE_COMMAND}
    -DPYTHON_EXE=${HS_PYTHON}
    -DGENERATOR=${HS_ROOT}/scripts/generate_luts.py
    -DSKIP_CODE=77 -DREQUIRE_PYTHON=${HS_REQUIRE_GENERATORS}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_color_luts_selfcheck.cmake)
set_tests_properties(unit_color_luts_check PROPERTIES SKIP_RETURN_CODE 77)

# Pin core/color/gamut_lut.h against its generator, and the generator's copies of
# the OKLab matrices and the gamut slack against core/color/color.h — drift there
# invalidates the brackets the per-pixel bisection is confined to. A full
# regeneration is ~1 min of numpy.
add_test(NAME unit_gamut_lut COMMAND ${CMAKE_COMMAND}
    -DPYTHON_EXE=${HS_PYTHON}
    -DGENERATOR=${HS_ROOT}/tools/gen_gamut_lut.py
    -DSKIP_CODE=77 -DREQUIRE_PYTHON=${HS_REQUIRE_GENERATORS}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_gamut_lut.cmake)
set_tests_properties(unit_gamut_lut PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 900)

# Pin hardware/pov_segment_map.json against hardware/pov_segment_map.h. The JSON
# is the firmware segment->canvas reference the simulator's cross-check (daydream
# tests/segment_crosscheck.test.js) reads, and the WASM install copies it into the
# daydream checkout alongside the module. Emitting it from a TU compiled against
# the header — rather than re-porting the mapping in JS — means a segment
# convention change here cannot leave the web tiling silently mis-partitioned.
add_executable(pov_segment_map_gen ${HS_ROOT}/tools/pov_segment_map_export.cpp)
target_include_directories(pov_segment_map_gen PRIVATE ${HS_ROOT})
_hs_warnings(pov_segment_map_gen)
_hs_sanitize(pov_segment_map_gen)
add_test(NAME unit_pov_segment_map_golden COMMAND ${CMAKE_COMMAND}
    -DEXPORTER=$<TARGET_FILE:pov_segment_map_gen>
    -DCOMMITTED=${HS_ROOT}/hardware/pov_segment_map.json
    -DGENERATED=${CMAKE_CURRENT_BINARY_DIR}/pov_segment_map.generated.json
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_pov_segment_map.cmake)

# The run_tests build sets NO math flags, so its NaN-clamp
# assertions run under default IEEE semantics where the hs::clamp NaN->hi
# contract holds trivially — the at-risk codegen the WASM release build actually
# ships (-ffast-math, which implies the contract-breaking -ffinite-math-only,
# rescued only by re-applying -fno-finite-math-only) is never exercised. This
# separate executable recompiles the same clamp-before-cast assertions under
# those exact release math flags so the contract is verified against the shipping
# codegen. It complements (does not replace) the compile-time #error guard in
# core/engine/platform.h that fails any build where -fno-finite-math-only is lost.
add_executable(fastmath_clamp_check
  ${CMAKE_CURRENT_SOURCE_DIR}/fastmath_clamp_check.cpp
  ${HS_ROOT}/core/engine/memory.cpp
  ${HS_ROOT}/core/engine/reaction_graph.cpp)
target_compile_definitions(fastmath_clamp_check PRIVATE HS_TEST_BUILD)
target_include_directories(fastmath_clamp_check PRIVATE ${HS_ROOT} ${HS_ROOT}/core)
# The WASM release math flags (CMakeLists.txt). -fno-finite-math-only MUST follow
# -ffast-math: it re-enables the NaN/Inf semantics the contract needs. If a future
# edit drops it, core/engine/platform.h's #error makes this TU fail to compile — exactly
# the regression signal we want.
target_compile_options(fastmath_clamp_check PRIVATE
  -O3 -ffast-math -fno-finite-math-only)
_hs_warnings(fastmath_clamp_check)
_hs_sanitize(fastmath_clamp_check)
add_test(NAME unit_fastmath_clamp COMMAND fastmath_clamp_check)

# Relax-bake re-derivation gate. The committed payloads in
# core/mesh/relax_bakes_generated.h are host-generated and loaded verbatim on
# both host and device. This harness recompiles every bake-bearing recipe with
# HS_RELAX_BAKE_VERIFY, which re-runs exactly bake.iterations smoothing steps
# and asserts the result reproduces the committed bits — so a recipe edit that
# silently drifts the geometry (even one that leaves topology unchanged) fails
# the suite instead of shipping a stale asset. It MUST be its own TU: VERIFY
# rewires SolidBuilder::relax_baked() globally, which the normal run_tests build
# needs in its load-the-bits form.
add_executable(relax_bake_verify ${HS_ROOT}/tools/relax_bake_harness.cpp
  ${HS_ROOT}/core/engine/memory.cpp ${HS_ROOT}/core/engine/reaction_graph.cpp)
target_compile_definitions(relax_bake_verify PRIVATE HS_TEST_BUILD HS_RELAX_BAKE_VERIFY)
# IEEE float semantics whatever the build's global flags are: the committed
# payload is a host-IEEE-captured constant, and live fast-math arithmetic does
# not reproduce its bits. The device loads the payload, it never re-derives it.
target_compile_options(relax_bake_verify PRIVATE -fno-fast-math)
target_include_directories(relax_bake_verify PRIVATE ${HS_ROOT} ${HS_ROOT}/core)
_hs_warnings(relax_bake_verify)
_hs_sanitize(relax_bake_verify)
add_test(NAME unit_relax_bake_verify COMMAND relax_bake_verify)

# The generator half of the same harness (HS_RELAX_BAKE_EXTRACT dumps the
# payloads that tools/relax_bakes.py emits into the header). Built here so a
# change that breaks generation is caught in CI; run it to regenerate:
#   relax_bake_gen | python tools/relax_bakes.py emit --stdin
add_executable(relax_bake_gen ${HS_ROOT}/tools/relax_bake_harness.cpp
  ${HS_ROOT}/core/engine/memory.cpp ${HS_ROOT}/core/engine/reaction_graph.cpp)
target_compile_definitions(relax_bake_gen PRIVATE HS_TEST_BUILD HS_RELAX_BAKE_EXTRACT)
# Same IEEE pin as the verifier above: the two halves must share float semantics
# or generation emits payloads the gate rejects.
target_compile_options(relax_bake_gen PRIVATE -fno-fast-math)
target_include_directories(relax_bake_gen PRIVATE ${HS_ROOT} ${HS_ROOT}/core)
_hs_warnings(relax_bake_gen)

# The south-pole Y-clip renormalization in the rasterizer (Screen::AntiAlias's
# bilinear-tap fold, and the H_VIRT handling threaded through geometry/scan/plot/
# sdf) only does real work when hs::H_OFFSET > 0 — the hardware value, fixed to 0
# on every host/WASM build, so the device's most subtle hot-loop path otherwise
# ships with zero automated coverage on the one platform with no debugger or
# console. This executable recompiles the same engine with HS_TEST_H_OFFSET=3 so
# the whole pipeline carries the device offset and the renorm runs against an
# energy-conservation oracle. It MUST be its own TU: an offset-3 and an offset-0
# PhiLUT<H>/TrigLUT<W,H> instantiation would clash under ODR (see core/engine/platform.h
# and tests/test_h_offset_renorm.h).
add_executable(h_offset_renorm_check
  ${CMAKE_CURRENT_SOURCE_DIR}/h_offset_renorm_check.cpp
  ${HS_ROOT}/core/engine/memory.cpp
  ${HS_ROOT}/core/engine/reaction_graph.cpp)
target_compile_definitions(h_offset_renorm_check PRIVATE HS_TEST_BUILD HS_TEST_H_OFFSET=3)
target_include_directories(h_offset_renorm_check PRIVATE ${HS_ROOT} ${HS_ROOT}/core)
_hs_warnings(h_offset_renorm_check)
_hs_sanitize(h_offset_renorm_check)
add_test(NAME unit_h_offset_renorm COMMAND h_offset_renorm_check)

# Host stack high-water-mark gate across every effect. Built -Os + the device
# math flags so codegen frame sizes track the Phantasm size build as closely as a
# host compiler can; x86-64 frames (8-byte pointers) make it a conservative upper
# bound for the 4-byte-pointer device render/init path (see stack_measure.cpp).
# Run as a CTest: it exits non-zero if the worst effect exceeds its stack budget,
# catching a regression before it can overflow the device's MPU-guarded stack.
add_executable(stack_measure
  ${CMAKE_CURRENT_SOURCE_DIR}/stack_measure.cpp
  ${HS_ROOT}/core/engine/memory.cpp
  ${HS_ROOT}/core/engine/reaction_graph.cpp)
target_compile_definitions(stack_measure PRIVATE HS_TEST_BUILD)
target_include_directories(stack_measure PRIVATE ${HS_ROOT} ${HS_ROOT}/core)
target_compile_options(stack_measure PRIVATE -Os -ffast-math -fno-finite-math-only)
_hs_warnings(stack_measure)
add_test(NAME unit_stack_budget COMMAND stack_measure)

# Arena high-water gate across every effect. Reports per-effect global-arena
# usage to size GLOBAL_ARENA_SIZE, and exits non-zero if the worst single-effect
# total exceeds the device budget — catching an arena regression before it can
# OOM-trap on hardware (see arena_measure.cpp header). Host figures are a
# conservative upper bound vs the 32-bit device, so a host pass implies device.
add_executable(arena_measure
  ${CMAKE_CURRENT_SOURCE_DIR}/arena_measure.cpp
  ${HS_ROOT}/core/engine/memory.cpp
  ${HS_ROOT}/core/engine/reaction_graph.cpp)
target_compile_definitions(arena_measure PRIVATE HS_TEST_BUILD)
target_include_directories(arena_measure PRIVATE ${HS_ROOT} ${HS_ROOT}/core)
target_compile_options(arena_measure PRIVATE -Os -ffast-math -fno-finite-math-only)
_hs_warnings(arena_measure)
add_test(NAME unit_arena_budget COMMAND arena_measure)

# Render micro-benchmark, built at two opt levels to bracket the size config's
# per-pixel codegen cost (perf_bench.cpp). Not CTests — wall time is noisy. The
# O3 variant stays in the default build as a compile guard so perf_bench.cpp
# can't bit-rot silently; the Os variant is built on demand.
foreach(_opt Os O3)
  set(_perf_exclude EXCLUDE_FROM_ALL)
  if(_opt STREQUAL "O3")
    set(_perf_exclude "")
  endif()
  add_executable(perf_bench_${_opt} ${_perf_exclude}
    ${CMAKE_CURRENT_SOURCE_DIR}/perf_bench.cpp
    ${HS_ROOT}/core/engine/memory.cpp ${HS_ROOT}/core/engine/reaction_graph.cpp)
  target_compile_definitions(perf_bench_${_opt} PRIVATE HS_TEST_BUILD)
  target_include_directories(perf_bench_${_opt} PRIVATE ${HS_ROOT} ${HS_ROOT}/core)
  target_compile_options(perf_bench_${_opt} PRIVATE -${_opt} -ffast-math -fno-finite-math-only)
  _hs_warnings(perf_bench_${_opt})
endforeach()

# Drivers for the HS_AA_AUDIT mesh-scan coverage instrumentation, whose hooks
# are compiled into core/render/scan.h and core/render/sdf.h and whose recipes
# are documented in docs/islamicstars_next_levers.md and
# docs/probe_path_open_items.md. Not CTests — aa_audit prints a report with no
# pass/fail and aa_search is a brute-force parameter sweep. They stay in the
# default build as compile guards against bit-rot, like perf_bench.
foreach(_aa audit search)
  add_executable(aa_${_aa}
    ${CMAKE_CURRENT_SOURCE_DIR}/aa_${_aa}_main.cpp
    ${HS_ROOT}/core/engine/memory.cpp ${HS_ROOT}/core/engine/reaction_graph.cpp)
  target_compile_definitions(aa_${_aa} PRIVATE HS_TEST_BUILD HS_AA_AUDIT)
  target_include_directories(aa_${_aa} PRIVATE ${HS_ROOT} ${HS_ROOT}/core)
  target_compile_options(aa_${_aa} PRIVATE -O2 -ffast-math -fno-finite-math-only)
  _hs_warnings(aa_${_aa})
endforeach()

# Divergence trace for perturbations of the feedback colour path
# (feedback_divergence.cpp). Not a CTest — comparing two code paths needs two
# builds, so it dumps frames in one build and compares dumps in another. Stays
# in the default build as a compile guard against bit-rot.
add_executable(feedback_divergence
  ${CMAKE_CURRENT_SOURCE_DIR}/feedback_divergence.cpp
  ${HS_ROOT}/core/engine/memory.cpp)
target_compile_definitions(feedback_divergence PRIVATE HS_TEST_BUILD)
target_include_directories(feedback_divergence PRIVATE ${HS_ROOT} ${HS_ROOT}/core)
target_compile_options(feedback_divergence PRIVATE -O2 -ffast-math -fno-finite-math-only)
_hs_warnings(feedback_divergence)

# Generator of record for core/color/srgb_decode_lut.h. EXCLUDE_FROM_ALL: run
# with no argument it rewrites a committed header, so the sharded CI job only
# compiles it; the srgb-decode-provenance job runs it to a temp path (argv[1])
# and diffs. unit_color pins the table itself.
add_executable(srgb_decode_gen EXCLUDE_FROM_ALL
  ${HS_ROOT}/scripts/generate_srgb_decode.cpp)
target_compile_definitions(srgb_decode_gen PRIVATE HS_TEST_BUILD)
target_include_directories(srgb_decode_gen PRIVATE ${HS_ROOT} ${HS_ROOT}/core)
_hs_warnings(srgb_decode_gen)
