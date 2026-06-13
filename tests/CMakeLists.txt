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
target_compile_definitions(run_tests PRIVATE HS_TEST_BUILD)

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
    3dmath memory spatial scb sdf conway hankin geometry mesh solids
    reaction_graph color easing_waves platform filter plot_scan canvas
    scan mesh_raster transformers generators animation effects hd107s
    pov_segmented pov_single pov_sync param_marshal util presets styles death)
foreach(_mod IN LISTS _hs_test_modules)
    add_test(NAME unit_${_mod} COMMAND run_tests ${_mod})
endforeach()

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
