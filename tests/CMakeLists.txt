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

set(HS_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/..)

add_executable(run_tests
  ${CMAKE_CURRENT_SOURCE_DIR}/run_tests.cpp
  ${HS_ROOT}/core/memory.cpp
  ${HS_ROOT}/core/reaction_graph.cpp)

# HS_TEST_BUILD widens the test-only inline/arena budgets (MSVC STL's
# std::function is larger than libc++/inplace_function on the device); the
# device footprint is unchanged. See core/animation.h and core/memory.h.
#
# HS_SCAN_METRICS compiles the per-pixel g_scan_metrics counters back in (off in
# the device/WASM release build, where they would be a per-pixel global store on
# the hot loop). The Face LUT test reads lut_hits/exact_hits to decide which
# code path each sample exercised, so the suite needs them live. See
# core/platform.h (HS_SCAN_METRIC) and core/sdf.h Face::distance.
target_compile_definitions(run_tests PRIVATE HS_TEST_BUILD HS_SCAN_METRICS)

# Headers are included as "tests/..." and "core/..." from the repo root.
target_include_directories(run_tests PRIVATE ${HS_ROOT} ${HS_ROOT}/core)

# The native suite is the one build where diagnostics can run: the device/WASM
# build compiles the arena guards out, so ASan/UBSan here are the only automated
# cover for the arena-lifetime hazards the death harness can't reach. -Werror is
# on so suite hygiene is self-enforcing — the fail-fast doctrine applied to the
# build itself; a new warning is a hard failure, not noise scrolled past.
target_compile_options(run_tests PRIVATE -Wall -Wextra -Werror)

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
if(HS_SANITIZE)
  # Suppress only signed-integer-overflow inside vendored FastNoiseLite (its
  # hashing wraps signed ints by design); first-party code keeps full coverage.
  set(_hs_ubsan_ignore "-fsanitize-ignorelist=${CMAKE_CURRENT_SOURCE_DIR}/ubsan-ignorelist.txt")
  if(WIN32)
    target_compile_options(run_tests PRIVATE
      -fsanitize=undefined -fsanitize-trap=undefined ${_hs_ubsan_ignore}
      -fno-omit-frame-pointer -g)
    message(STATUS "Unit tests: UBSan (trap mode) enabled — no ASan runtime on emsdk clang")
  else()
    target_compile_options(run_tests PRIVATE
      -fsanitize=address,undefined -fno-sanitize-recover=all ${_hs_ubsan_ignore}
      -fno-omit-frame-pointer -g)
    target_link_options(run_tests PRIVATE -fsanitize=address,undefined)
    message(STATUS "Unit tests: ASan + UBSan enabled (HS_SANITIZE=ON)")
  endif()
endif()

enable_testing()

# One CTest test per module so CI shards (and local `ctest -R <name>`) can
# target subsets without re-running everything.
set(_hs_test_modules
    3dmath concepts memory spatial scb sdf conway hankin geometry mesh solids
    reaction_graph color palettes easing_waves platform filter plot_scan canvas
    scan mesh_raster transformers noise generators animation effects dma_core hd107s
    pov_segmented pov_single pov_sync param_marshal wasm_predicates util
    presets styles death)
foreach(_mod IN LISTS _hs_test_modules)
    add_test(NAME unit_${_mod} COMMAND run_tests ${_mod})
endforeach()

# Pin _hs_test_modules against run_tests' HS_TEST_MODULE_LIST: run_tests
# --check-modules exits non-zero if the CMake list and the roster diverge in
# either direction, so a module added to one but not the other can never run
# silently in CI sharding.
add_test(NAME unit_module_roster
    COMMAND run_tests --check-modules ${_hs_test_modules})

# Pin core/color_luts.h against its generator; skips (exit 77) without Python.
find_program(HS_PYTHON NAMES python3 python)
add_test(NAME unit_color_luts COMMAND ${CMAKE_COMMAND}
    -DPYTHON_EXE=${HS_PYTHON}
    -DGENERATOR=${HS_ROOT}/scripts/generate_luts.py
    -DCOMMITTED=${HS_ROOT}/core/color_luts.h
    -DGENERATED=${CMAKE_CURRENT_BINARY_DIR}/color_luts.generated.h
    -DSKIP_CODE=77
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_color_luts.cmake)
set_tests_properties(unit_color_luts PROPERTIES SKIP_RETURN_CODE 77)

# Run the generator's own --check self-validation (monotonicity + sRGB
# round-trip), which CI ran only as a shell step. Absent without Python.
if(HS_PYTHON)
  add_test(NAME unit_color_luts_check
      COMMAND ${HS_PYTHON} ${HS_ROOT}/scripts/generate_luts.py --check)
endif()

# Finding 369: the run_tests build sets NO math flags, so its NaN-clamp
# assertions run under default IEEE semantics where the hs::clamp NaN->hi
# contract holds trivially — the at-risk codegen the WASM release build actually
# ships (-ffast-math, which implies the contract-breaking -ffinite-math-only,
# rescued only by re-applying -fno-finite-math-only) is never exercised. This
# separate executable recompiles the same clamp-before-cast assertions under
# those exact release math flags so the contract is verified against the shipping
# codegen. It complements (does not replace) the compile-time #error guard in
# core/platform.h that fails any build where -fno-finite-math-only is lost.
add_executable(fastmath_clamp_check
  ${CMAKE_CURRENT_SOURCE_DIR}/fastmath_clamp_check.cpp
  ${HS_ROOT}/core/memory.cpp
  ${HS_ROOT}/core/reaction_graph.cpp)
target_compile_definitions(fastmath_clamp_check PRIVATE HS_TEST_BUILD)
target_include_directories(fastmath_clamp_check PRIVATE ${HS_ROOT} ${HS_ROOT}/core)
# The WASM release math flags (CMakeLists.txt). -fno-finite-math-only MUST follow
# -ffast-math: it re-enables the NaN/Inf semantics the contract needs. If a future
# edit drops it, core/platform.h's #error makes this TU fail to compile — exactly
# the regression signal we want.
target_compile_options(fastmath_clamp_check PRIVATE
  -O3 -ffast-math -fno-finite-math-only)
add_test(NAME unit_fastmath_clamp COMMAND fastmath_clamp_check)

# The south-pole Y-clip renormalization in the rasterizer (Screen::AntiAlias's
# bilinear-tap fold, and the H_VIRT handling threaded through geometry/scan/plot/
# sdf) only does real work when hs::H_OFFSET > 0 — the hardware value, fixed to 0
# on every host/WASM build, so the device's most subtle hot-loop path otherwise
# ships with zero automated coverage on the one platform with no debugger or
# console. This executable recompiles the same engine with HS_TEST_H_OFFSET=3 so
# the whole pipeline carries the device offset and the renorm runs against an
# energy-conservation oracle. It MUST be its own TU: an offset-3 and an offset-0
# PhiLUT<H>/TrigLUT<W,H> instantiation would clash under ODR (see core/platform.h
# and tests/test_h_offset_renorm.h).
add_executable(h_offset_renorm_check
  ${CMAKE_CURRENT_SOURCE_DIR}/h_offset_renorm_check.cpp
  ${HS_ROOT}/core/memory.cpp
  ${HS_ROOT}/core/reaction_graph.cpp)
target_compile_definitions(h_offset_renorm_check PRIVATE HS_TEST_BUILD HS_TEST_H_OFFSET=3)
target_include_directories(h_offset_renorm_check PRIVATE ${HS_ROOT} ${HS_ROOT}/core)
# No -Werror here (matching fastmath_clamp_check): this single-purpose TU pulls in
# the full engine barrel but uses only the AntiAlias/TrigLUT slice, so the barrel's
# many static helpers would trip -Wunused-function. Engine-header hygiene is
# enforced by the main run_tests target, which compiles the same headers -Werror.
add_test(NAME unit_h_offset_renorm COMMAND h_offset_renorm_check)

# Host stack high-water-mark gate across every effect. Built -Os + the device
# math flags so codegen frame sizes track the Phantasm size build as closely as a
# host compiler can; x86-64 frames (8-byte pointers) make it a conservative upper
# bound for the 4-byte-pointer device render/init path (see stack_measure.cpp).
# Run as a CTest: it exits non-zero if the worst effect exceeds its stack budget,
# catching a regression before it can overflow the device's MPU-guarded stack.
add_executable(stack_measure
  ${CMAKE_CURRENT_SOURCE_DIR}/stack_measure.cpp
  ${HS_ROOT}/core/memory.cpp
  ${HS_ROOT}/core/reaction_graph.cpp)
target_compile_definitions(stack_measure PRIVATE HS_TEST_BUILD)
target_include_directories(stack_measure PRIVATE ${HS_ROOT} ${HS_ROOT}/core)
target_compile_options(stack_measure PRIVATE -Os -ffast-math -fno-finite-math-only)
add_test(NAME unit_stack_budget COMMAND stack_measure)

# Arena high-water probe: reports per-effect global-arena usage to size
# GLOBAL_ARENA_SIZE. Not a CTest (informational, like perf_bench).
add_executable(arena_measure
  ${CMAKE_CURRENT_SOURCE_DIR}/arena_measure.cpp
  ${HS_ROOT}/core/memory.cpp
  ${HS_ROOT}/core/reaction_graph.cpp)
target_compile_definitions(arena_measure PRIVATE HS_TEST_BUILD)
target_include_directories(arena_measure PRIVATE ${HS_ROOT} ${HS_ROOT}/core)
target_compile_options(arena_measure PRIVATE -Os -ffast-math -fno-finite-math-only)

# Render micro-benchmark, built at two opt levels to bracket the size config's
# per-pixel codegen cost (perf_bench.cpp). Not CTests — wall time is noisy.
foreach(_opt Os O3)
  add_executable(perf_bench_${_opt}
    ${CMAKE_CURRENT_SOURCE_DIR}/perf_bench.cpp
    ${HS_ROOT}/core/memory.cpp ${HS_ROOT}/core/reaction_graph.cpp)
  target_compile_definitions(perf_bench_${_opt} PRIVATE HS_TEST_BUILD)
  target_include_directories(perf_bench_${_opt} PRIVATE ${HS_ROOT} ${HS_ROOT}/core)
  target_compile_options(perf_bench_${_opt} PRIVATE -${_opt} -ffast-math -fno-finite-math-only)
endforeach()
