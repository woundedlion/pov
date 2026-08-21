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

# cmake/toolchain-native-clang.cmake falls back to an unversioned `clang` from
# PATH, so the major is whatever the developer happens to have installed. Clang
# 18 is the oldest the relaxed-template-template opt-in below was written for;
# nothing older is built anywhere, so refuse it by version rather than let it
# fail as a compile error.
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 18)
  message(FATAL_ERROR
    "Clang ${CMAKE_CXX_COMPILER_VERSION} is below the 18 this suite is built "
    "against. Point CMAKE_CXX_COMPILER at a newer clang++, or let "
    "cmake/toolchain-native-clang.cmake resolve the emsdk's.")
endif()

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

add_library(hs_host_abi INTERFACE)
target_compile_definitions(hs_host_abi INTERFACE
  HS_TIMELINE_MAX_ANIM_BYTES=256)

add_library(hs_host_storage INTERFACE)
target_link_libraries(hs_host_storage INTERFACE hs_host_abi)
target_compile_definitions(hs_host_storage INTERFACE
  HS_GLOBAL_ARENA_BYTES=8388608)

add_library(hs_host_effect_runtime INTERFACE)
target_link_libraries(hs_host_effect_runtime INTERFACE hs_host_storage)
target_compile_definitions(hs_host_effect_runtime INTERFACE
  HS_EXTERNAL_PARAM_STORAGE=1)

add_library(hs_test_whitebox INTERFACE)
target_link_libraries(hs_test_whitebox INTERFACE hs_host_effect_runtime)
target_compile_definitions(hs_test_whitebox INTERFACE
  HS_ENABLE_TEST_HOOKS=1
  HS_ENABLE_TEST_ORACLES=1
  HS_ENABLE_SHADER_WORKBENCH_DYNAMIC_BACKEND=1
  HS_ENABLE_CHAIN_INTERPRETER=1
  HS_ENABLE_EFFECT_CONTROL_API=1)

add_library(hs_test_support INTERFACE)
target_link_libraries(hs_test_support INTERFACE hs_host_storage)
target_compile_definitions(hs_test_support INTERFACE
  HS_ENABLE_EFFECT_REGISTRY=1
  HS_ENABLE_PARAM_GUI_BRIDGE=1
  HS_ENABLE_TEST_HOOKS=1
  HS_ENABLE_TEST_ORACLES=1
  HS_ENABLE_SHADER_WORKBENCH_DYNAMIC_BACKEND=1
  HS_ENABLE_CHAIN_INTERPRETER=1
  HS_ENABLE_STRUCTURAL_AUDITS=1
  HS_ENABLE_EFFECT_CONTROL_API=1)

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
  ${HS_ROOT}/core/engine/static_storage.cpp
  ${HS_ROOT}/core/spatial/reaction_graph.cpp)

# The host configuration widens inline/arena budgets because 64-bit pointers
# inflate several objects past their 32-bit device footprints.
#
# HS_SCAN_METRICS compiles the per-pixel g_scan_metrics counters back in (off in
# the device/WASM release build, where they would be a per-pixel global store on
# the hot loop). The Face distance-oracle test reads exact_hits to tell culled
# samples from evaluated ones, so the suite needs them live. See
# core/platform/platform.h (HS_SCAN_METRIC) and core/render/sdf.h Face::distance.
target_compile_definitions(run_tests PRIVATE HS_SCAN_METRICS
  "HS_SHADER_CHAIN_CATALOG_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/data/shader_chain_catalog.json\"")
target_link_libraries(run_tests PRIVATE hs_test_support)

# Census of the engine's HS_CHECK sites, so the death harness can print how
# much of the fail-fast surface its cases actually pin. Generated, never
# committed: a written-down guard list drifts away from the guards.
include(${CMAKE_CURRENT_SOURCE_DIR}/count_guard_sites.cmake)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/death_guard_sites.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/death_guard_sites.h @ONLY)

# Headers are included as "tests/..." and "core/..." from the repo root; the
# binary dir carries the generated guard census.
target_include_directories(run_tests PRIVATE ${HS_ROOT} ${HS_ROOT}/core
                           ${CMAKE_CURRENT_BINARY_DIR})

_hs_warnings(run_tests)

# Opt-in sanitizers. Default OFF so the pre-commit gate and the default
# `cmake --preset tests` stay fast (ASan ~2x); enable with -DHS_SANITIZE=ON.
# Two CI jobs (.github/workflows/ci.yml) build with it on every push and pull
# request, one per flavour, because the flavour depends on the runtime the
# toolchain ships:
#   * Linux/full-LLVM ("sanitizers"): ASan + UBSan with diagnostics
#     (-fno-sanitize-recover makes the first error a hard failure rather than
#     logged-and-continue).
#   * Windows emsdk clang ("windows-tests", one instrumented target): ships no
#     ASan runtime for the MSVC target, so use runtime-free UBSan in trap mode —
#     a UB hit becomes an immediate trap, consistent with the engine's always-on
#     fail-fast doctrine.
option(HS_SANITIZE "Build the unit tests with sanitizers (ASan+UBSan, or UBSan-trap on Windows)" OFF)

# ThreadSanitizer. The canvas and pov_segmented modules drive the Effect double
# buffer from two real threads (a render thread against a display-ISR stand-in),
# which is the engine's only genuine concurrency; nothing else in the suite would
# gain from TSan. It cannot share a binary with ASan, so it is a separate option
# and a separate CI job ("thread-sanitizer" in .github/workflows/ci.yml).
# Linux/macOS only: clang ships no TSan runtime for the Windows target.
option(HS_TSAN "Build the unit tests with ThreadSanitizer (excludes HS_SANITIZE)" OFF)
if(HS_TSAN)
  if(HS_SANITIZE)
    message(FATAL_ERROR
        "HS_TSAN and HS_SANITIZE are mutually exclusive — ASan and TSan cannot "
        "instrument the same binary. Configure two build directories.")
  endif()
  if(WIN32)
    message(FATAL_ERROR
        "HS_TSAN: clang ships no ThreadSanitizer runtime for the Windows "
        "target. Use a Linux or macOS toolchain.")
  endif()
  message(STATUS "Unit tests: ThreadSanitizer enabled (HS_TSAN=ON)")
endif()

# Applies the sanitizer flavour to one target; a no-op when both HS_SANITIZE and
# HS_TSAN are OFF. Every CTest-registered executable the sanitizers or
# thread-sanitizer job runs must call this, or it rides that job uninstrumented.
# The stack_measure/arena_measure budget gates are the deliberate exception:
# ASan redzones and instrumentation inflate the frame sizes and codegen they
# measure against the device build.
function(_hs_sanitize target)
  if(HS_TSAN)
    target_compile_options(${target} PRIVATE
      -fsanitize=thread -fno-omit-frame-pointer -g)
    target_link_options(${target} PRIVATE -fsanitize=thread)
    return()
  endif()
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

add_executable(mindsplatter_replay_gen EXCLUDE_FROM_ALL
  ${HS_ROOT}/tools/mindsplatter_replay_gen.cpp
  ${HS_ROOT}/core/engine/memory.cpp
  ${HS_ROOT}/core/engine/static_storage.cpp
  ${HS_ROOT}/core/spatial/reaction_graph.cpp)
target_compile_definitions(mindsplatter_replay_gen PRIVATE
  HS_PROFILE_MINDSPLATTER_COUNTS)
target_link_libraries(mindsplatter_replay_gen PRIVATE hs_test_whitebox)
target_include_directories(mindsplatter_replay_gen PRIVATE ${HS_ROOT} ${HS_ROOT}/core)
target_compile_options(mindsplatter_replay_gen PRIVATE
  -O3 -ffast-math -fno-finite-math-only)
_hs_warnings(mindsplatter_replay_gen)

add_executable(mindsplatter_palette_gen EXCLUDE_FROM_ALL
  ${HS_ROOT}/tools/mindsplatter_palette_gen.cpp
  ${HS_ROOT}/core/engine/memory.cpp
  ${HS_ROOT}/core/engine/static_storage.cpp
  ${HS_ROOT}/core/spatial/reaction_graph.cpp)
target_include_directories(mindsplatter_palette_gen PRIVATE ${HS_ROOT} ${HS_ROOT}/core)
target_link_libraries(mindsplatter_palette_gen PRIVATE hs_host_abi)
_hs_warnings(mindsplatter_palette_gen)

# Tolerance pin for the bank mindsplatter_palette_gen emits. Its own TU rather
# than a run_tests module: the committed bank is 393,216 B of static data no
# other test reads.
add_executable(mindsplatter_palette_check
  ${CMAKE_CURRENT_SOURCE_DIR}/mindsplatter_palette_check.cpp
  ${HS_ROOT}/core/engine/memory.cpp
  ${HS_ROOT}/core/engine/static_storage.cpp
  ${HS_ROOT}/core/spatial/reaction_graph.cpp)
target_include_directories(mindsplatter_palette_check PRIVATE ${HS_ROOT} ${HS_ROOT}/core)
target_link_libraries(mindsplatter_palette_check PRIVATE hs_host_abi)
_hs_warnings(mindsplatter_palette_check)
_hs_sanitize(mindsplatter_palette_check)
add_test(NAME unit_mindsplatter_palette COMMAND mindsplatter_palette_check)

add_executable(mindsplatter_replay
  ${HS_ROOT}/tools/mindsplatter_replay_main.cpp
  ${HS_ROOT}/core/engine/memory.cpp
  ${HS_ROOT}/core/engine/static_storage.cpp
  ${HS_ROOT}/core/spatial/reaction_graph.cpp)
target_link_libraries(mindsplatter_replay PRIVATE hs_test_whitebox)
target_include_directories(mindsplatter_replay PRIVATE ${HS_ROOT} ${HS_ROOT}/core)
target_compile_options(mindsplatter_replay PRIVATE
  -O3 -ffast-math -fno-finite-math-only)
_hs_warnings(mindsplatter_replay)
_hs_sanitize(mindsplatter_replay)
add_test(NAME unit_mindsplatter_replay COMMAND mindsplatter_replay)

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
    reaction_graph color palettes easing_waves interpolate platform profiling
    pullback shader_chain filter
    plot_scan canvas
    scan mesh_raster transformers noise noise_field projections
    animation effects
    effects_smoke effect_factory shader_workbench curl_lattice facet_grid
    composed_effect
    shapeshifter_oracle shapeshifter_tiles dma_core hd107s
    dma_controller
    pov_segmented pov_single pov_sync param_marshal wasm_predicates led
    presets styles shading death)
# The suite has a real hang mode (a Canvas constructor that spin-waits for a
# buffer that never frees — see test_effects.h), and CTest's 1500 s default would
# burn that once per shard. The whole roster runs in ~2 min unsanitized, but
# ASan + UBSan cost the slowest module (effects) ~2.7x — 227 s plain against
# 620 s sanitized — so a single budget for both builds is either useless plain
# or red sanitized.
set(_hs_test_timeout 600)
if(HS_SANITIZE OR HS_TSAN)
  set(_hs_test_timeout 1800)
endif()
foreach(_mod IN LISTS _hs_test_modules)
    add_test(NAME unit_${_mod} COMMAND run_tests ${_mod})
    set_tests_properties(unit_${_mod} PROPERTIES TIMEOUT ${_hs_test_timeout})
endforeach()
if(NOT HS_SANITIZE AND NOT HS_TSAN)
    set_tests_properties(unit_effects_smoke PROPERTIES TIMEOUT 1200)
endif()

# Pin _hs_test_modules against run_tests' HS_TEST_MODULE_LIST: run_tests
# --check-modules exits non-zero if the CMake list and the roster diverge in
# either direction, so a module added to one but not the other can never run
# silently in CI sharding.
add_test(NAME unit_module_roster
    COMMAND run_tests --check-modules ${_hs_test_modules})

# Pin run_tests.cpp's per-module include block against the roster, catching an
# include left behind after a module is removed (which compiles silently), and
# every tests/test_*.h on disk against that include block, catching a test file
# nothing ever compiles. The exempt list is pinned too: an entry must still be
# included from tests/ or tools/.
add_test(NAME unit_module_includes COMMAND ${CMAKE_COMMAND}
    -DSRC=${CMAKE_CURRENT_SOURCE_DIR}/run_tests.cpp
    -DTESTS_DIR=${CMAKE_CURRENT_SOURCE_DIR}
    -DTOOLS_DIR=${HS_ROOT}/tools
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_includes.cmake)

# Require every test case defined in a tests/*.h header to have a call site: in
# the same header for module cases, or anywhere under tests/ and tools/ for
# shared helpers.
add_test(NAME unit_case_calls COMMAND ${CMAKE_COMMAND}
    -DTESTS_DIR=${CMAKE_CURRENT_SOURCE_DIR}
    -DTOOLS_DIR=${HS_ROOT}/tools
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_case_calls.cmake)

# Pin the device inline-storage audit's hand-written type pack against the
# AnimationBase subclasses actually declared under core/. The per-type
# static_assert in Timeline::add_get only fires for types a build add()s, so an
# animation missing from the pack is measured by nothing.
add_test(NAME unit_anim_pack COMMAND ${CMAKE_COMMAND}
    -DCORE_DIR=${HS_ROOT}/core
    -DANIM_HEADER=${HS_ROOT}/core/animation/animation.h
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_anim_pack.cmake)

# Pin the in-tree patches to vendored FastNoiseLite.h against the patch record
# in core/vendor/FastNoiseLite_config.h. A version bump that drops them still
# compiles: the size guard goes inert and the sampler loses selective -O3, both
# silently.
add_test(NAME unit_vendor_patches COMMAND ${CMAKE_COMMAND}
    -DHEADER=${HS_ROOT}/core/vendor/FastNoiseLite.h
    -DCONFIG=${HS_ROOT}/core/vendor/FastNoiseLite_config.h
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_vendor_patches.cmake)

# Pin the Phantasm linker script's fork base and flash-routing changes. It is a
# second vendored fork -- of the stock Teensy 4 imxrt1062.ld -- that nothing in
# the build re-derives, and a lost routing change reads as a Teensy size-gate
# red far from its cause, on the hosts that have PlatformIO at all.
add_test(NAME unit_linker_patches COMMAND ${CMAKE_COMMAND}
    -DSCRIPT=${HS_ROOT}/tools/phantasm.ld
    -DPLATFORMIO_INI=${HS_ROOT}/platformio.ini
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_linker_patches.cmake)

# The three generated-artifact gates below re-run a Python generator and diff it
# against the committed header. Default OFF keeps a bare local checkout green:
# a missing python3/numpy exits 77 and CTest reports the gate as skipped. That
# skip is invisible in a green run, so the CI job that provisions the
# interpreter passes -DHS_REQUIRE_GENERATORS=ON, which turns the same missing
# python3/numpy into a hard failure — a broken provisioning step can no longer
# silently retire all three gates.
option(HS_REQUIRE_GENERATORS
    "Fail (rather than skip) the generated-artifact gates when python3/numpy is missing" OFF)

# NO_CACHE: a cached NOTFOUND would outlive installing an interpreter and leave
# the four gates below permanently skipped, which a green run does not show.
find_program(HS_PYTHON NAMES python3 python NO_CACHE)

if(NOT HS_PYTHON AND HS_REQUIRE_GENERATORS)
  message(FATAL_ERROR "Python is required to validate pullback manifests")
endif()

# The pullback manifest is the one generated artifact the suite COMPILES
# against (tests/test_shader_workbench.h): without the header nothing builds, so the
# manifest-derived assertions compile out under HS_HAVE_PULLBACK_MANIFEST=0 and
# the two manifest CTests lose the targets they run, leaving the rest of the
# suite runnable in a bare checkout. They still register, as 77-skips over the
# stand-in script, so a Python-less configure shrinks the pass count rather than
# the test count -- `100% tests passed out of 76` reads exactly as green as `out
# of 78`. The retirement warns at configure time, and the CI job that provisions
# the interpreter passes -DHS_REQUIRE_GENERATORS=ON, so it cannot happen there.
if(HS_PYTHON)
target_compile_definitions(run_tests PRIVATE HS_HAVE_PULLBACK_MANIFEST=1)

set(PULLBACK_MANIFEST_DIR ${HS_ROOT}/tests/data/pullback)
set(PULLBACK_MANIFEST_HEADER
    ${CMAKE_CURRENT_BINARY_DIR}/pullback_manifest.generated.h)
add_custom_command(
  OUTPUT ${PULLBACK_MANIFEST_HEADER}
  COMMAND ${HS_PYTHON} ${HS_ROOT}/tools/generate_pullback_manifest_header.py
          --manifest-dir ${PULLBACK_MANIFEST_DIR}
          --output ${PULLBACK_MANIFEST_HEADER}
  DEPENDS
    ${HS_ROOT}/tools/generate_pullback_manifest_header.py
    ${PULLBACK_MANIFEST_DIR}/schema.json
    ${PULLBACK_MANIFEST_DIR}/programs.json
    ${PULLBACK_MANIFEST_DIR}/peirce_fast_square.json
    ${PULLBACK_MANIFEST_DIR}/hue_rotation_noise_luts.json
  VERBATIM)
add_custom_target(pullback_manifest_header DEPENDS ${PULLBACK_MANIFEST_HEADER})
add_dependencies(run_tests pullback_manifest_header)

add_executable(pullback_manifest_check
  ${CMAKE_CURRENT_SOURCE_DIR}/pullback_manifest_check.cpp
  ${PULLBACK_MANIFEST_HEADER})
add_dependencies(pullback_manifest_check pullback_manifest_header)
target_include_directories(pullback_manifest_check PRIVATE
                           ${CMAKE_CURRENT_BINARY_DIR} ${HS_ROOT})
_hs_warnings(pullback_manifest_check)
_hs_sanitize(pullback_manifest_check)
add_test(NAME unit_pullback_manifest COMMAND pullback_manifest_check)

set(HS_PULLBACK_CAPTURE_SOURCE_ROOT "${HS_ROOT}" CACHE PATH
    "Checkout whose ShaderWorkbench sources the pullback capture producer builds")
add_executable(pullback_capture_native
  ${HS_ROOT}/tools/pullback_capture_native.cpp
  ${HS_PULLBACK_CAPTURE_SOURCE_ROOT}/core/engine/memory.cpp
  ${HS_PULLBACK_CAPTURE_SOURCE_ROOT}/core/engine/static_storage.cpp
  ${HS_PULLBACK_CAPTURE_SOURCE_ROOT}/core/spatial/reaction_graph.cpp)
target_include_directories(pullback_capture_native PRIVATE
  ${HS_PULLBACK_CAPTURE_SOURCE_ROOT}
  ${HS_PULLBACK_CAPTURE_SOURCE_ROOT}/core)
target_compile_definitions(pullback_capture_native PRIVATE
  HS_ENABLE_EFFECT_REGISTRY=1
  HS_ENABLE_PARAM_GUI_BRIDGE=1
  HS_ENABLE_TEST_HOOKS=1
  HS_ENABLE_TEST_ORACLES=1
  HS_ENABLE_SHADER_WORKBENCH_DYNAMIC_BACKEND=0
  HS_ENABLE_CHAIN_INTERPRETER=0
  HS_ENABLE_EFFECT_CONTROL_API=1)
target_link_libraries(pullback_capture_native PRIVATE hs_host_storage)
_hs_warnings(pullback_capture_native)
_hs_sanitize(pullback_capture_native)
add_test(NAME unit_pullback_capture_determinism
  COMMAND ${CMAKE_COMMAND}
    -DPRODUCER=$<TARGET_FILE:pullback_capture_native>
    -DOUTPUT_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DPYTHON=${HS_PYTHON}
    -DCAPTURE_SCRIPT=${HS_ROOT}/tools/pullback_capture.py
    -DMANIFEST_DIR=${PULLBACK_MANIFEST_DIR}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_pullback_capture.cmake)
else()
target_compile_definitions(run_tests PRIVATE HS_HAVE_PULLBACK_MANIFEST=0)
foreach(_pullback_test IN ITEMS unit_pullback_manifest unit_pullback_capture_determinism)
  add_test(NAME ${_pullback_test} COMMAND ${CMAKE_COMMAND}
      -DTEST_NAME=${_pullback_test} -DSKIP_CODE=77
      -P ${CMAKE_CURRENT_SOURCE_DIR}/skip_without_python.cmake)
  set_tests_properties(${_pullback_test} PROPERTIES SKIP_RETURN_CODE 77)
endforeach()
message(WARNING
  "No python3/python found: the pullback manifest header cannot be generated, "
  "so unit_pullback_manifest and unit_pullback_capture_determinism skip and "
  "the manifest assertions in unit_shader_workbench compile out. "
  "Install an interpreter, or configure with -DHS_REQUIRE_GENERATORS=ON to "
  "make this a hard error.")
endif()

# Pin core/color/color_luts.h against its generator.
add_test(NAME unit_color_luts COMMAND ${CMAKE_COMMAND}
    -DPYTHON_EXE=${HS_PYTHON}
    -DGENERATOR=${HS_ROOT}/scripts/generate_luts.py
    -DCOMMITTED=${HS_ROOT}/core/color/color_luts.h
    -DGENERATED=${CMAKE_CURRENT_BINARY_DIR}/color_luts.generated.h
    -DSKIP_CODE=77 -DREQUIRE_PYTHON=${HS_REQUIRE_GENERATORS}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_color_luts.cmake)
set_tests_properties(unit_color_luts PROPERTIES SKIP_RETURN_CODE 77)

# Host self-tests for the per-commit Teensy size trail (tools/teensy_size_trail.py):
# the ELF32 section-header parser the pre-commit recorder reads firmware images
# with, and the regression classifier. Pure stdlib Python — the ci.yml
# teensy-gate-tests job runs the same file through its unittest discovery.
add_test(NAME unit_teensy_size_trail COMMAND ${CMAKE_COMMAND}
    -DPYTHON_EXE=${HS_PYTHON}
    -DHS_ROOT=${HS_ROOT}
    -DSKIP_CODE=77 -DREQUIRE_PYTHON=${HS_REQUIRE_GENERATORS}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_teensy_size_trail.cmake)
set_tests_properties(unit_teensy_size_trail PROPERTIES SKIP_RETURN_CODE 77)

# Run the generator's own --check self-validation (monotonicity + sRGB
# round-trip), which CI ran only as a shell step.
add_test(NAME unit_color_luts_check COMMAND ${CMAKE_COMMAND}
    -DPYTHON_EXE=${HS_PYTHON}
    -DGENERATOR=${HS_ROOT}/scripts/generate_luts.py
    -DSKIP_CODE=77 -DREQUIRE_PYTHON=${HS_REQUIRE_GENERATORS}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_color_luts_selfcheck.cmake)
set_tests_properties(unit_color_luts_check PROPERTIES SKIP_RETURN_CODE 77)

# Pin core/color/gamut_lut.h against its generator, the generator's copies of the
# OKLab matrices and the gamut slack against core/color/color.h, and its diamond
# angle against core/math/3dmath.h — drift there invalidates the brackets the
# per-pixel bisection is confined to, or mis-indexes them. A full regeneration is
# ~1 min of numpy.
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
# core/platform/platform.h that fails any build where -fno-finite-math-only is lost.
add_executable(fastmath_clamp_check
  ${CMAKE_CURRENT_SOURCE_DIR}/fastmath_clamp_check.cpp
  ${HS_ROOT}/core/engine/memory.cpp
  ${HS_ROOT}/core/engine/static_storage.cpp
  ${HS_ROOT}/core/spatial/reaction_graph.cpp)
target_include_directories(fastmath_clamp_check PRIVATE ${HS_ROOT} ${HS_ROOT}/core)
target_link_libraries(fastmath_clamp_check PRIVATE hs_host_abi)
# The WASM release math flags (CMakeLists.txt). -fno-finite-math-only MUST follow
# -ffast-math: it re-enables the NaN/Inf semantics the contract needs. If a future
# edit drops it, core/platform/platform.h's #error makes this TU fail to compile — exactly
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
  ${HS_ROOT}/core/engine/memory.cpp
  ${HS_ROOT}/core/engine/static_storage.cpp
  ${HS_ROOT}/core/spatial/reaction_graph.cpp)
target_compile_definitions(relax_bake_verify PRIVATE HS_RELAX_BAKE_VERIFY)
target_link_libraries(relax_bake_verify PRIVATE hs_host_abi)
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
  ${HS_ROOT}/core/engine/memory.cpp
  ${HS_ROOT}/core/engine/static_storage.cpp
  ${HS_ROOT}/core/spatial/reaction_graph.cpp)
target_compile_definitions(relax_bake_gen PRIVATE HS_RELAX_BAKE_EXTRACT)
target_link_libraries(relax_bake_gen PRIVATE hs_host_abi)
# Same IEEE pin as the verifier above: the two halves must share float semantics
# or generation emits payloads the gate rejects.
target_compile_options(relax_bake_gen PRIVATE -fno-fast-math)
target_include_directories(relax_bake_gen PRIVATE ${HS_ROOT} ${HS_ROOT}/core)
_hs_warnings(relax_bake_gen)
_hs_sanitize(relax_bake_gen)

# Pin the asset's FORM against the emitter: unit_relax_bake_verify above proves
# the committed payloads re-derive bit-exact, but nothing re-ran the emitter, so
# a drift in the banner or the chunking would surface only as a large spurious
# reformat inside the next legitimate regeneration.
add_test(NAME unit_relax_bake_form COMMAND ${CMAKE_COMMAND}
    -DPYTHON_EXE=${HS_PYTHON}
    -DSCRIPT=${HS_ROOT}/tools/relax_bakes.py
    -DHARNESS=$<TARGET_FILE:relax_bake_gen>
    -DDUMP=${CMAKE_CURRENT_BINARY_DIR}/relax_bakes.dump
    -DSKIP_CODE=77 -DREQUIRE_PYTHON=${HS_REQUIRE_GENERATORS}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_relax_bakes.cmake)
set_tests_properties(unit_relax_bake_form PROPERTIES SKIP_RETURN_CODE 77)

# The south-pole Y-clip renormalization in the rasterizer (Screen::AntiAlias's
# bilinear-tap fold, and the H_VIRT handling threaded through geometry/scan/plot/
# sdf) only does real work when hs::H_OFFSET > 0 — the hardware value, fixed to 0
# on every host/WASM build, so the device's most subtle hot-loop path otherwise
# ships with zero automated coverage on the one platform with no debugger or
# console. This executable recompiles the same engine with HS_TEST_H_OFFSET=3 so
# the whole pipeline carries the device offset and the renorm runs against an
# energy-conservation oracle. It MUST be its own TU: an offset-3 and an offset-0
# PhiLUT<H>/TrigLUT<W,H> instantiation would clash under ODR (see core/platform/platform.h
# and tests/test_h_offset_renorm.h).
add_executable(h_offset_renorm_check
  ${CMAKE_CURRENT_SOURCE_DIR}/h_offset_renorm_check.cpp
  ${HS_ROOT}/core/engine/memory.cpp
  ${HS_ROOT}/core/engine/static_storage.cpp
  ${HS_ROOT}/core/spatial/reaction_graph.cpp)
target_compile_definitions(h_offset_renorm_check PRIVATE HS_TEST_H_OFFSET=3)
target_link_libraries(h_offset_renorm_check PRIVATE hs_host_abi)
target_include_directories(h_offset_renorm_check PRIVATE ${HS_ROOT} ${HS_ROOT}/core)
_hs_warnings(h_offset_renorm_check)
_hs_sanitize(h_offset_renorm_check)
add_test(NAME unit_h_offset_renorm COMMAND h_offset_renorm_check)

# Host stack high-water-mark gate across every effect. Built -Os + the device
# math flags so codegen frame sizes track the Phantasm size build as closely as a
# host compiler can. stack_measure.cpp heap-constructs effects like the firmware
# factory so their object storage is not charged to the measured call stack.
# Run as a CTest: it exits non-zero if the worst effect exceeds its stack budget,
# catching a regression before it can overflow the device's MPU-guarded stack.
#
# The budget is the device's own DTCM stack reservation, read out of
# tools/teensy_budgets.json (phantasm.regions.ram1.free_min_bytes) — the same
# figure tools/teensy_gate.py holds free for locals and derives the ITCM code
# ceiling from. Reading it here means changing the device reservation moves this
# gate with it, instead of leaving a stale copy behind. The file is JSONC, so
# line comments are stripped before string(JSON).
file(READ ${HS_ROOT}/tools/teensy_budgets.json _hs_teensy_budgets_jsonc)
string(REGEX REPLACE "//[^\n]*" "" _hs_teensy_budgets_json
  "${_hs_teensy_budgets_jsonc}")
string(JSON HS_DEVICE_STACK_FLOOR_BYTES ERROR_VARIABLE _hs_stack_floor_error
  GET "${_hs_teensy_budgets_json}" phantasm regions ram1 free_min_bytes)
if(_hs_stack_floor_error OR NOT HS_DEVICE_STACK_FLOOR_BYTES MATCHES "^[0-9]+$")
  message(FATAL_ERROR
    "tools/teensy_budgets.json phantasm.regions.ram1.free_min_bytes is missing "
    "or not a number (${_hs_stack_floor_error}); stack_measure has no device "
    "stack reservation to gate against.")
endif()

add_executable(stack_measure
  ${CMAKE_CURRENT_SOURCE_DIR}/stack_measure.cpp
  ${HS_ROOT}/core/engine/memory.cpp
  ${HS_ROOT}/core/engine/static_storage.cpp
  ${HS_ROOT}/core/spatial/reaction_graph.cpp)
target_compile_definitions(stack_measure PRIVATE
  HS_DEVICE_STACK_FLOOR_BYTES=${HS_DEVICE_STACK_FLOOR_BYTES})
target_link_libraries(stack_measure PRIVATE hs_host_effect_runtime)
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
  ${HS_ROOT}/core/engine/static_storage.cpp
  ${HS_ROOT}/core/spatial/reaction_graph.cpp)
target_link_libraries(arena_measure PRIVATE hs_host_effect_runtime)
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
    ${HS_ROOT}/core/engine/memory.cpp
    ${HS_ROOT}/core/engine/static_storage.cpp
    ${HS_ROOT}/core/spatial/reaction_graph.cpp)
  target_link_libraries(perf_bench_${_opt} PRIVATE hs_host_effect_runtime)
  target_include_directories(perf_bench_${_opt} PRIVATE ${HS_ROOT} ${HS_ROOT}/core)
  target_compile_options(perf_bench_${_opt} PRIVATE -${_opt} -ffast-math -fno-finite-math-only)
  _hs_warnings(perf_bench_${_opt})
endforeach()

# Divergence trace for perturbations of the feedback colour path
# (feedback_divergence.cpp). Not a CTest — comparing two code paths needs two
# builds, so it dumps frames in one build and compares dumps in another. Stays
# in the default build as a compile guard against bit-rot.
add_executable(feedback_divergence
  ${CMAKE_CURRENT_SOURCE_DIR}/feedback_divergence.cpp
  ${HS_ROOT}/core/engine/memory.cpp
  ${HS_ROOT}/core/engine/static_storage.cpp)
target_include_directories(feedback_divergence PRIVATE ${HS_ROOT} ${HS_ROOT}/core)
target_link_libraries(feedback_divergence PRIVATE hs_host_abi)
target_compile_options(feedback_divergence PRIVATE -O2 -ffast-math -fno-finite-math-only)
_hs_warnings(feedback_divergence)

# Generator of record for core/color/srgb_decode_lut.h. EXCLUDE_FROM_ALL: run
# with no argument it rewrites a committed header, so the sharded CI job only
# compiles it; the srgb-decode-provenance job runs it to a temp path (argv[1])
# and diffs. unit_color pins the table itself.
add_executable(srgb_decode_gen EXCLUDE_FROM_ALL
  ${HS_ROOT}/scripts/generate_srgb_decode.cpp)
target_include_directories(srgb_decode_gen PRIVATE ${HS_ROOT} ${HS_ROOT}/core)
_hs_warnings(srgb_decode_gen)

# Compile coverage for every EXCLUDE_FROM_ALL target above. Those targets are
# out of the default build (their run rewrites a committed artifact, or is too
# noisy to gate), so an API rename breaks their sources while every job stays
# green. CI builds this one target instead of naming them; the list is read back
# off the buildsystem, so a new excluded target joins it with no workflow edit.
get_property(_hs_dir_targets DIRECTORY PROPERTY BUILDSYSTEM_TARGETS)
set(_hs_excluded)
foreach(_hs_target IN LISTS _hs_dir_targets)
  get_target_property(_hs_is_excluded ${_hs_target} EXCLUDE_FROM_ALL)
  if(_hs_is_excluded)
    list(APPEND _hs_excluded ${_hs_target})
  endif()
endforeach()
if(NOT _hs_excluded)
  message(FATAL_ERROR
    "no EXCLUDE_FROM_ALL targets found — excluded_targets would gate nothing; "
    "drop it along with the CI steps that build it")
endif()
message(STATUS "EXCLUDE_FROM_ALL compile coverage: ${_hs_excluded}")
add_custom_target(excluded_targets)
add_dependencies(excluded_targets ${_hs_excluded})
