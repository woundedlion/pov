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
# Emscripten configurations are excluded.
if(EMSCRIPTEN)
  return()
endif()

# The native-clang toolchain sets CMAKE_LINKER_TYPE (CMake 3.29+), which older
# CMake silently ignores.
cmake_minimum_required(VERSION 3.29)

# Allow standalone configuration (cmake -S tests) as well as being pulled in by
# a native root build.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(holosphere_tests CXX)
  set(HS_PULLBACK_CAPTURE_SOURCE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/.." CACHE PATH
      "Checkout whose sources a pullback capture build compiles")
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
endif()

# Python is required; optional generator dependencies may skip locally.
option(HS_REQUIRE_GENERATORS
    "Fail rather than skip when optional generator dependencies are missing" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# The native-clang toolchain may fall back to an unversioned `clang` from PATH;
# refuse anything older than 18.
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 18)
  message(FATAL_ERROR
    "Clang ${CMAKE_CXX_COMPILER_VERSION} is below the 18 this suite is built "
    "against. Point CMAKE_CXX_COMPILER at a newer clang++, or let "
    "cmake/toolchain-native-clang.cmake resolve the emsdk's.")
endif()

# Clang < 19 needs relaxed template-template matching opted in for the
# harness's `template <int, int> class E` parameters; Clang 20 removed the flag.
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 19)
  add_compile_options(-frelaxed-template-template-args)
endif()

set(HS_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/..)
set(HS_ENGINE_SOURCES
  ${HS_ROOT}/core/memory.cpp
  ${HS_ROOT}/core/engine/static_storage.cpp
  ${HS_ROOT}/core/spatial/reaction_graph.cpp)
set(HS_ENGINE_INCLUDES ${HS_ROOT} ${HS_ROOT}/core)

# The host configuration widens inline/arena budgets because 64-bit pointers
# inflate several objects past their 32-bit device footprints.
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
  HS_ENABLE_CHAIN_INTERPRETER=1
  HS_ENABLE_EFFECT_CONTROL_API=1)

add_library(hs_test_support INTERFACE)
target_link_libraries(hs_test_support INTERFACE hs_host_storage)
target_compile_definitions(hs_test_support INTERFACE
  HS_ENABLE_EFFECT_REGISTRY=1
  HS_ENABLE_PARAM_GUI_BRIDGE=1
  HS_ENABLE_TEST_HOOKS=1
  HS_ENABLE_TEST_ORACLES=1
  HS_ENABLE_CHAIN_INTERPRETER=1
  HS_ENABLE_STRUCTURAL_AUDITS=1
  HS_ENABLE_EFFECT_CONTROL_API=1)

# Host diagnostics and warning flags; per-target optimization and math flags
# remain local.
function(_hs_warnings target)
  target_compile_options(${target} PRIVATE -Wall -Wextra -Werror)
endfunction()

add_executable(run_tests
  ${CMAKE_CURRENT_SOURCE_DIR}/run_tests.cpp
  ${HS_ENGINE_SOURCES})

# HS_SCAN_METRICS compiles in the per-pixel g_scan_metrics counters
# (HS_SCAN_METRIC, core/platform/profiling.h).
target_compile_definitions(run_tests PRIVATE HS_SCAN_METRICS
  "HS_SHADER_CHAIN_CATALOG_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/data/shader_chain_catalog.json\""
  "HS_PROMOTED_PATTERNS_DIR=\"${HS_ROOT}/patterns\"")
target_link_libraries(run_tests PRIVATE hs_test_support)

set(_color_header_sources)
file(GLOB _color_headers CONFIGURE_DEPENDS "${HS_ROOT}/core/color/*.h")
list(FILTER _color_headers EXCLUDE REGEX "_luts?\\.h$")
foreach(_header_path IN LISTS _color_headers)
  get_filename_component(_header "${_header_path}" NAME_WE)
  set(_source "${CMAKE_CURRENT_BINARY_DIR}/color_header_${_header}.cpp")
  file(CONFIGURE OUTPUT "${_source}" CONTENT "#include \"color/${_header}.h\"\n")
  list(APPEND _color_header_sources "${_source}")
endforeach()
add_library(color_header_checks OBJECT ${_color_header_sources})
target_include_directories(color_header_checks PRIVATE ${HS_ROOT}/core)
_hs_warnings(color_header_checks)
add_dependencies(run_tests color_header_checks)


set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/guard_directories.txt"
  "${HS_ROOT}/tools/teensy_budgets.json"
  "${HS_ROOT}/core/platform/build_features.h")

# Generated census of the engine's HS_CHECK sites, for the death harness's
# coverage report.
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/guard_directories.txt" HS_GUARD_DIRS)
set(_guard_sources)
foreach(_dir IN LISTS HS_GUARD_DIRS)
  file(GLOB_RECURSE _found CONFIGURE_DEPENDS
       "${HS_ROOT}/${_dir}/*.h" "${HS_ROOT}/${_dir}/*.cpp"
       "${HS_ROOT}/${_dir}/*.ino")
  list(APPEND _guard_sources ${_found})
endforeach()
set(_guard_output "${CMAKE_CURRENT_BINARY_DIR}/death_guard_sites.h")
add_custom_command(OUTPUT "${_guard_output}"
  COMMAND "${CMAKE_COMMAND}" "-DHS_ROOT=${HS_ROOT}"
          "-DHS_GUARD_OUTPUT=${_guard_output}"
          "-DHS_GUARD_DIRS=${HS_GUARD_DIRS}"
          -P "${CMAKE_CURRENT_SOURCE_DIR}/count_guard_sites.cmake"
  DEPENDS ${_guard_sources} guard_directories.txt count_guard_sites.cmake death_guard_sites.h.in
  VERBATIM)
add_custom_target(death_guard_census DEPENDS "${_guard_output}")
add_dependencies(run_tests death_guard_census)
add_test(NAME unit_death_pins COMMAND "${Python3_EXECUTABLE}"
         "${CMAKE_CURRENT_SOURCE_DIR}/check_death_pins.py" "${HS_ROOT}")

# Headers are included as "tests/..." and "core/..." from the repo root; the
# binary dir carries the generated guard census.
target_include_directories(run_tests PRIVATE ${HS_ENGINE_INCLUDES}
                           ${CMAKE_CURRENT_BINARY_DIR})

_hs_warnings(run_tests)

# Opt-in sanitizers for the native suite; enable with -DHS_SANITIZE=ON.
#   * Linux/full-LLVM: ASan + UBSan; the first error is a hard failure.
#   * Windows emsdk clang: no ASan runtime for the MSVC target, so runtime-free
#     UBSan in trap mode.
option(HS_SANITIZE "Build the unit tests with sanitizers (ASan+UBSan, or UBSan-trap on Windows)" OFF)

# ThreadSanitizer, for the Effect double buffer driven from two real threads.
# Cannot share a binary with ASan. Linux/macOS only: clang ships no TSan runtime
# for the Windows target.
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
# HS_TSAN are OFF. Participating test targets opt in by calling this helper.
# Stack and arena budget measurements exclude sanitizer instrumentation.
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
  # hashing wraps signed ints).
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

add_executable(mindsplatter_replay_gen
  ${HS_ROOT}/tools/mindsplatter_replay_gen.cpp
  ${HS_ENGINE_SOURCES})
target_compile_definitions(mindsplatter_replay_gen PRIVATE
  HS_PROFILE_PLOT_RENDER_COUNTS)
target_link_libraries(mindsplatter_replay_gen PRIVATE hs_test_whitebox)
target_include_directories(mindsplatter_replay_gen PRIVATE ${HS_ENGINE_INCLUDES})
target_compile_options(mindsplatter_replay_gen PRIVATE
  -O3 -ffast-math -fno-finite-math-only)
_hs_warnings(mindsplatter_replay_gen)
_hs_sanitize(mindsplatter_replay_gen)

add_test(NAME unit_mindsplatter_replay_form COMMAND ${CMAKE_COMMAND}
  -DPYTHON=${Python3_EXECUTABLE}
  -DGENERATOR=$<TARGET_FILE:mindsplatter_replay_gen>
  -DCOMMITTED=${HS_ROOT}/tests/mindsplatter_replay_corpus.h
  -DGENERATED=${CMAKE_CURRENT_BINARY_DIR}/mindsplatter_replay_corpus.generated.h
  -DCOMPILER_ID=${CMAKE_CXX_COMPILER_ID}
  -DCOMPILER_VERSION=${CMAKE_CXX_COMPILER_VERSION}
  -DREQUIRE_COMPILER_MATCH=${HS_REQUIRE_GENERATORS}
  -DSKIP_CODE=77
  -P ${CMAKE_CURRENT_SOURCE_DIR}/check_mindsplatter_replay_form.cmake)
set_tests_properties(unit_mindsplatter_replay_form PROPERTIES SKIP_RETURN_CODE 77)
add_test(NAME unit_mindsplatter_replay_refresh COMMAND ${CMAKE_COMMAND}
  -DREFRESH_ONLY=ON
  -DGENERATOR=$<TARGET_FILE:mindsplatter_replay_gen>
  -DCOMMITTED=${HS_ROOT}/tests/mindsplatter_replay_corpus.h
  -DGENERATED=${CMAKE_CURRENT_BINARY_DIR}/mindsplatter_replay_corpus.refreshed.h
  -P ${CMAKE_CURRENT_SOURCE_DIR}/check_mindsplatter_replay_form.cmake)

add_test(NAME unit_mindsplatter_replay_revision_missing COMMAND ${CMAKE_COMMAND}
  -P ${CMAKE_CURRENT_SOURCE_DIR}/check_mindsplatter_replay_revision_missing.cmake)
add_test(NAME unit_mindsplatter_replay_revision_drift COMMAND ${CMAKE_COMMAND}
  -P ${CMAKE_CURRENT_SOURCE_DIR}/check_mindsplatter_replay_revision_drift.cmake)

add_executable(mindsplatter_palette_gen
  ${HS_ROOT}/tools/mindsplatter_palette_gen.cpp
  ${HS_ENGINE_SOURCES})
target_include_directories(mindsplatter_palette_gen PRIVATE ${HS_ENGINE_INCLUDES})
target_link_libraries(mindsplatter_palette_gen PRIVATE hs_host_abi)
_hs_warnings(mindsplatter_palette_gen)
_hs_sanitize(mindsplatter_palette_gen)

add_custom_target(regenerate_mindsplatter_palette
  COMMAND $<TARGET_FILE:mindsplatter_palette_gen>
          ${HS_ROOT}/core/color/mindsplatter_palette_luts.h
  DEPENDS mindsplatter_palette_gen
  VERBATIM)

add_test(NAME unit_mindsplatter_palette_form COMMAND ${CMAKE_COMMAND}
  -DGENERATOR=$<TARGET_FILE:mindsplatter_palette_gen>
  -DCOMMITTED=${HS_ROOT}/core/color/mindsplatter_palette_luts.h
  -DGENERATED=${CMAKE_CURRENT_BINARY_DIR}/mindsplatter_palette_luts.generated.h
  -P ${CMAKE_CURRENT_SOURCE_DIR}/check_mindsplatter_palette_form.cmake)

# Standalone tolerance check for the palette data, with its own error report
# and provenance gate.
add_executable(mindsplatter_palette_check
  ${CMAKE_CURRENT_SOURCE_DIR}/mindsplatter_palette_check.cpp
  ${HS_ENGINE_SOURCES})
target_include_directories(mindsplatter_palette_check PRIVATE ${HS_ENGINE_INCLUDES})
target_link_libraries(mindsplatter_palette_check PRIVATE hs_host_abi)
_hs_warnings(mindsplatter_palette_check)
_hs_sanitize(mindsplatter_palette_check)
add_test(NAME unit_mindsplatter_palette COMMAND mindsplatter_palette_check)

add_executable(mindsplatter_replay
  ${HS_ROOT}/tools/mindsplatter_replay_main.cpp
  ${HS_ENGINE_SOURCES})
target_link_libraries(mindsplatter_replay PRIVATE hs_test_whitebox)
target_include_directories(mindsplatter_replay PRIVATE ${HS_ENGINE_INCLUDES})
target_compile_options(mindsplatter_replay PRIVATE
  -O3 -ffast-math -fno-finite-math-only)
_hs_warnings(mindsplatter_replay)
_hs_sanitize(mindsplatter_replay)
add_test(NAME unit_mindsplatter_replay COMMAND mindsplatter_replay)

# Opt-in Clang source-based coverage for run_tests. The instrumentation changes
# codegen. Needs clang_rt.profile, which the Windows emsdk clang does not ship.
option(HS_COVERAGE "Instrument run_tests for llvm-cov source-based coverage" OFF)
if(HS_COVERAGE)
  target_compile_options(run_tests PRIVATE -fprofile-instr-generate -fcoverage-mapping)
  target_link_options(run_tests PRIVATE -fprofile-instr-generate)
endif()

enable_testing()

# One CTest test per module so CI shards (and local `ctest -R <name>`) can
# target subsets without re-running everything.
set(_hs_test_modules
    3dmath concepts memory spatial scb sdf ray conway conway_morph conway_continuity
    partition_seam
    conway_soak opchain_probe opchain_arena_survey hankin hyper_lattice ray_events ray_demonstrators sdf_patterns geometry
    spherical_field spherical_harmonics
    mesh solids
    reaction_graph color palettes easing_waves interpolate platform profiling
    pullback shader_chain filter
    plot_scan canvas
    scan mesh_raster transformers noise noise_field projections
    animation effects mindsplatter
    effects_smoke effect_factory lattice_melt kaleidoscope_smooth
    composed_effect
    shapeshifter_oracle shapeshifter_tiles dma_core hd107s
    dma_controller
    pov_segmented pov_single pov_sync param_marshal wasm_predicates led
    presets styles shading death)
# Bound each module independently.
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

add_test(NAME unit_death_child_dispatch COMMAND ${CMAKE_COMMAND}
    -DTEST_BIN=$<TARGET_FILE:run_tests>
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_death_child_dispatch.cmake)

# run_tests --check-modules fails if _hs_test_modules and HS_TEST_MODULE_LIST
# diverge in either direction.
add_test(NAME unit_module_roster
    COMMAND run_tests --check-modules ${_hs_test_modules})

# Pin run_tests.cpp's per-module include block against the roster, and every
# tests/*.h and tests/*.hpp against that include block.
add_test(NAME unit_module_includes COMMAND ${CMAKE_COMMAND}
    -DSRC=${CMAKE_CURRENT_SOURCE_DIR}/run_tests.cpp
    -DTESTS_DIR=${CMAKE_CURRENT_SOURCE_DIR}
    -DTOOLS_DIR=${HS_ROOT}/tools
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_includes.cmake)

add_test(NAME unit_module_entry_controls COMMAND ${CMAKE_COMMAND}
    -DGATE=${CMAKE_CURRENT_SOURCE_DIR}/check_includes.cmake
    -DWORK=${CMAKE_CURRENT_BINARY_DIR}/module_entry_controls
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_includes_selftest.cmake)

add_test(NAME unit_header_sections COMMAND ${CMAKE_COMMAND}
    -DWORK=${CMAKE_CURRENT_BINARY_DIR}/header_sections_fixture
    -DCALL_GATE=${CMAKE_CURRENT_SOURCE_DIR}/check_case_calls.cmake
    -DINCLUDE_GATE=${CMAKE_CURRENT_SOURCE_DIR}/check_includes.cmake
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_header_sections_selftest.cmake)

# Check reachability and assertion reach for column-0 void/bool/int/size_t
# cases in tests/*.h and tests/*.hpp, including run_*_cases drivers. Death
# case_* definitions in test_death.h need reachability but not assertion reach.
# Shared helpers resolve across tests/ and tools/; member definitions are
# excluded.
add_test(NAME unit_case_calls COMMAND ${CMAKE_COMMAND}
    -DCACHE_FILE=${CMAKE_CURRENT_BINARY_DIR}/case_calls.sha256
    -DTESTS_DIR=${CMAKE_CURRENT_SOURCE_DIR}
    -DTOOLS_DIR=${HS_ROOT}/tools
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_case_calls.cmake)

# Pin the device inline-storage audit's hand-written type pack against the
# AnimationBase subclasses declared under core/.
add_test(NAME unit_anim_pack COMMAND ${CMAKE_COMMAND}
    -DCORE_DIR=${HS_ROOT}/core
    -DANIM_HEADER=${HS_ROOT}/core/animation/animation.h
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_anim_pack.cmake)
add_test(NAME unit_anim_pack_declarations COMMAND ${CMAKE_COMMAND}
    -DFIXTURE_DIR=${CMAKE_CURRENT_BINARY_DIR}/anim_pack_fixture
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_anim_pack.cmake)

# Pin the in-tree patches to vendored FastNoiseLite.h against the patch record
# in core/vendor/FastNoiseLite_config.h.
add_test(NAME unit_vendor_patches COMMAND ${CMAKE_COMMAND}
    -DHEADER=${HS_ROOT}/core/vendor/FastNoiseLite.h
    -DCONFIG=${HS_ROOT}/core/vendor/FastNoiseLite_config.h
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_vendor_patches.cmake)

# Pin the Phantasm linker script (a fork of the stock Teensy 4 imxrt1062.ld)
# fork base and flash-routing changes.
add_test(NAME unit_linker_patches COMMAND ${CMAKE_COMMAND}
    -DSCRIPT=${HS_ROOT}/tools/phantasm.ld
    -DPLATFORMIO_INI=${HS_ROOT}/platformio.ini
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_linker_patches.cmake)

add_test(NAME unit_mirrored_prune COMMAND ${CMAKE_COMMAND}
    -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}/mirrored_prune_fixture
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_mirrored_prune.cmake)

set(HS_PYTHON "${Python3_EXECUTABLE}")
add_test(NAME unit_chain_capture_fixtures COMMAND "${Python3_EXECUTABLE}"
  "${CMAKE_CURRENT_SOURCE_DIR}/../tools/test_chain_capture_fixtures.py")

set(PULLBACK_MANIFEST_DIR ${HS_ROOT}/tests/data/pullback)
set(PULLBACK_MANIFEST_HEADER
    ${CMAKE_CURRENT_BINARY_DIR}/pullback_manifest.generated.h)
file(GLOB PULLBACK_MANIFEST_FILES CONFIGURE_DEPENDS
     ${PULLBACK_MANIFEST_DIR}/*.json)
add_custom_command(
  OUTPUT ${PULLBACK_MANIFEST_HEADER}
  COMMAND ${HS_PYTHON} ${HS_ROOT}/tools/generate_pullback_manifest_header.py
          --manifest-dir ${PULLBACK_MANIFEST_DIR}
          --output ${PULLBACK_MANIFEST_HEADER}
  DEPENDS
    ${HS_ROOT}/tools/generate_pullback_manifest_header.py
    ${HS_ROOT}/tools/pullback_operations.def
    ${PULLBACK_MANIFEST_FILES}
  VERBATIM)
add_custom_target(pullback_manifest_header DEPENDS ${PULLBACK_MANIFEST_HEADER})
add_dependencies(run_tests pullback_manifest_header)

add_executable(pullback_manifest_check
  ${CMAKE_CURRENT_SOURCE_DIR}/pullback_manifest_check.cpp
  ${PULLBACK_MANIFEST_HEADER})
add_dependencies(pullback_manifest_check pullback_manifest_header)
target_include_directories(pullback_manifest_check PRIVATE
                           ${CMAKE_CURRENT_BINARY_DIR} ${HS_ENGINE_INCLUDES})
_hs_warnings(pullback_manifest_check)
_hs_sanitize(pullback_manifest_check)
add_test(NAME unit_pullback_manifest COMMAND pullback_manifest_check)


set(HS_CAPTURE_ENGINE_SOURCES
  ${HS_PULLBACK_CAPTURE_SOURCE_ROOT}/core/memory.cpp
  ${HS_PULLBACK_CAPTURE_SOURCE_ROOT}/core/engine/static_storage.cpp
  ${HS_PULLBACK_CAPTURE_SOURCE_ROOT}/core/spatial/reaction_graph.cpp)
add_executable(pullback_capture_native
  ${HS_ROOT}/tools/pullback_capture_native.cpp
  ${HS_CAPTURE_ENGINE_SOURCES})
include(${HS_ROOT}/cmake/chain_capture_fixtures.cmake)
hs_chain_capture_fixtures(pullback_capture_native "${HS_ROOT}")
target_include_directories(pullback_capture_native PRIVATE
  ${HS_PULLBACK_CAPTURE_SOURCE_ROOT}
  ${HS_PULLBACK_CAPTURE_SOURCE_ROOT}/core)
target_compile_definitions(pullback_capture_native PRIVATE
  HS_ENABLE_EFFECT_REGISTRY=1
  HS_ENABLE_PARAM_GUI_BRIDGE=1
  HS_ENABLE_TEST_HOOKS=1
  HS_ENABLE_TEST_ORACLES=1
  HS_ENABLE_CHAIN_INTERPRETER=1
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
# Pin core/color/color_luts.h against its generator.
add_test(NAME unit_color_luts COMMAND ${CMAKE_COMMAND}
    -DPYTHON_EXE=${HS_PYTHON}
    -DGENERATOR=${HS_ROOT}/scripts/generate_luts.py
    -DCOMMITTED=${HS_ROOT}/core/color/color_luts.h
    -DGENERATED=${CMAKE_CURRENT_BINARY_DIR}/color_luts.generated.h
    -DSKIP_CODE=77 -DREQUIRE_GENERATORS=${HS_REQUIRE_GENERATORS}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_color_luts.cmake)
set_tests_properties(unit_color_luts PROPERTIES SKIP_RETURN_CODE 77)

# Host self-tests for tools/teensy_size_trail.py: the ELF32 section-header
# parser and the regression classifier.
add_test(NAME unit_teensy_size_trail COMMAND ${CMAKE_COMMAND}
    -DPYTHON_EXE=${HS_PYTHON}
    -DHS_ROOT=${HS_ROOT}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_teensy_size_trail.cmake)

# Pin core/color/gamut_lut.h against its generator, the generator's copies of the
# OKLab matrices and the gamut slack against core/color/color_space.h, and its
# diamond angle against core/math/3dmath.h.
add_test(NAME unit_gamut_lut COMMAND ${CMAKE_COMMAND}
    -DPYTHON_EXE=${HS_PYTHON}
    -DGENERATOR=${HS_ROOT}/tools/gen_gamut_lut.py
    -DCACHE_FILE=${CMAKE_CURRENT_BINARY_DIR}/gamut_lut_check.sha256
    -DSKIP_CODE=77 -DREQUIRE_GENERATORS=${HS_REQUIRE_GENERATORS}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_gamut_lut.cmake)
set_tests_properties(unit_gamut_lut PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 900)

# Pin hardware/pov_segment_map.json against hardware/pov_segment_map.h, emitted
# from a TU compiled against the header.
add_executable(pov_segment_map_gen ${HS_ROOT}/tools/pov_segment_map_export.cpp)
target_include_directories(pov_segment_map_gen PRIVATE ${HS_ROOT})
_hs_warnings(pov_segment_map_gen)
_hs_sanitize(pov_segment_map_gen)
add_test(NAME unit_pov_segment_map_golden COMMAND ${CMAKE_COMMAND}
    -DEXPORTER=$<TARGET_FILE:pov_segment_map_gen>
    -DCOMMITTED=${HS_ROOT}/hardware/pov_segment_map.json
    -DGENERATED=${CMAKE_CURRENT_BINARY_DIR}/pov_segment_map.generated.json
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_pov_segment_map.cmake)

# Exercise the release math contract at -O3. -fno-finite-math-only keeps
# NaN-aware clamps valid under -ffast-math.
add_executable(fastmath_clamp_check
  ${CMAKE_CURRENT_SOURCE_DIR}/fastmath_clamp_check.cpp
  ${HS_ENGINE_SOURCES})
target_include_directories(fastmath_clamp_check PRIVATE ${HS_ENGINE_INCLUDES})
target_link_libraries(fastmath_clamp_check PRIVATE hs_host_abi)
target_compile_definitions(fastmath_clamp_check PRIVATE HS_ENABLE_TEST_ORACLES=1)
# The WASM release math flags. -fno-finite-math-only MUST follow -ffast-math;
# core/platform/platform.h #errors without it.
target_compile_options(fastmath_clamp_check PRIVATE
  -O3 -ffast-math -fno-finite-math-only)
_hs_warnings(fastmath_clamp_check)
_hs_sanitize(fastmath_clamp_check)
add_test(NAME unit_fastmath_clamp COMMAND fastmath_clamp_check)

# Relax-bake re-derivation gate: recompiles every bake-bearing recipe with
# HS_RELAX_BAKE_VERIFY and asserts smoothing reproduces the committed bits in
# core/mesh/relax_bakes_generated.h. Its own TU: VERIFY rewires
# SolidBuilder::relax_baked() globally.
add_executable(relax_bake_verify ${HS_ROOT}/tools/relax_bake_harness.cpp
  ${HS_ENGINE_SOURCES})
target_compile_definitions(relax_bake_verify PRIVATE HS_RELAX_BAKE_VERIFY)
target_link_libraries(relax_bake_verify PRIVATE hs_host_abi)
# IEEE float semantics: the committed payload is a host-IEEE-captured constant.
target_compile_options(relax_bake_verify PRIVATE -fno-fast-math)
target_include_directories(relax_bake_verify PRIVATE ${HS_ENGINE_INCLUDES})
_hs_warnings(relax_bake_verify)
_hs_sanitize(relax_bake_verify)
add_test(NAME unit_relax_bake_verify COMMAND relax_bake_verify)

# The generator half of the same harness (HS_RELAX_BAKE_EXTRACT dumps the
# payloads that tools/relax_bakes.py emits into the header). To regenerate:
#   relax_bake_gen | python tools/relax_bakes.py emit --stdin
add_executable(relax_bake_gen ${HS_ROOT}/tools/relax_bake_harness.cpp
  ${HS_ENGINE_SOURCES})
target_compile_definitions(relax_bake_gen PRIVATE HS_RELAX_BAKE_EXTRACT)
target_link_libraries(relax_bake_gen PRIVATE hs_host_abi)
# Same IEEE pin as relax_bake_verify: the two halves must share float semantics.
target_compile_options(relax_bake_gen PRIVATE -fno-fast-math)
target_include_directories(relax_bake_gen PRIVATE ${HS_ENGINE_INCLUDES})
_hs_warnings(relax_bake_gen)
_hs_sanitize(relax_bake_gen)

# Pin the asset's form (banner, chunking) against the emitter.
add_test(NAME unit_relax_bake_form COMMAND ${CMAKE_COMMAND}
    -DPYTHON_EXE=${HS_PYTHON}
    -DSCRIPT=${HS_ROOT}/tools/relax_bakes.py
    -DHARNESS=$<TARGET_FILE:relax_bake_gen>
    -DDUMP=${CMAKE_CURRENT_BINARY_DIR}/relax_bakes.dump
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_relax_bakes.cmake)

# Legacy offset mapping needs a separate TU: offset-3 and offset-0
# PhiLUT/TrigLUT instantiations have incompatible definitions.
add_executable(h_offset_renorm_check
  ${CMAKE_CURRENT_SOURCE_DIR}/h_offset_renorm_check.cpp
  ${HS_ENGINE_SOURCES})
target_compile_definitions(h_offset_renorm_check PRIVATE HS_TEST_H_OFFSET=3)
# hs_host_storage, not hs_host_abi: the Scan and Plot cases draw through the
# global scratch arena.
target_link_libraries(h_offset_renorm_check PRIVATE hs_host_storage)
target_include_directories(h_offset_renorm_check PRIVATE ${HS_ENGINE_INCLUDES})
_hs_warnings(h_offset_renorm_check)
_hs_sanitize(h_offset_renorm_check)
add_test(NAME unit_h_offset_renorm COMMAND h_offset_renorm_check)

add_executable(h_offset_full_resolution_check
  ${CMAKE_CURRENT_SOURCE_DIR}/h_offset_renorm_check.cpp
  ${HS_ENGINE_SOURCES})
target_compile_definitions(h_offset_full_resolution_check PRIVATE
  HS_TEST_H_OFFSET=3 HS_OFFSET_FULL_RESOLUTION=1)
target_link_libraries(h_offset_full_resolution_check PRIVATE hs_host_storage)
target_include_directories(h_offset_full_resolution_check PRIVATE ${HS_ENGINE_INCLUDES})
_hs_warnings(h_offset_full_resolution_check)
_hs_sanitize(h_offset_full_resolution_check)
add_test(NAME unit_h_offset_full_resolution COMMAND h_offset_full_resolution_check)


add_executable(display_geometry_check
  ${CMAKE_CURRENT_SOURCE_DIR}/display_geometry_check.cpp
  ${HS_ENGINE_SOURCES})
target_compile_definitions(display_geometry_check PRIVATE HS_DISPLAY_PROFILE=1)
target_link_libraries(display_geometry_check PRIVATE hs_host_storage)
target_include_directories(display_geometry_check PRIVATE ${HS_ENGINE_INCLUDES})
_hs_warnings(display_geometry_check)
_hs_sanitize(display_geometry_check)
add_test(NAME unit_display_geometry COMMAND display_geometry_check)

add_executable(runtime_display_geometry_check
  ${CMAKE_CURRENT_SOURCE_DIR}/runtime_display_geometry_check.cpp
  ${HS_ENGINE_SOURCES})
target_compile_definitions(runtime_display_geometry_check PRIVATE
  HS_RUNTIME_DISPLAY_GEOMETRY=1)
target_link_libraries(runtime_display_geometry_check PRIVATE hs_host_storage)
target_include_directories(runtime_display_geometry_check PRIVATE ${HS_ENGINE_INCLUDES})
_hs_warnings(runtime_display_geometry_check)
_hs_sanitize(runtime_display_geometry_check)
add_test(NAME unit_runtime_display_geometry COMMAND runtime_display_geometry_check)

# Host stack high-water-mark gate across every effect, built -Os + the device
# math flags to track device frame sizes. Fails if the worst effect exceeds the
# device DTCM stack reservation, read from tools/teensy_budgets.json
# (phantasm.regions.ram1.free_min_bytes). The file is JSONC, so line comments
# are stripped before string(JSON).
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
  ${HS_ENGINE_SOURCES})
target_compile_definitions(stack_measure PRIVATE
  HS_DEVICE_STACK_FLOOR_BYTES=${HS_DEVICE_STACK_FLOOR_BYTES})
target_link_libraries(stack_measure PRIVATE hs_host_effect_runtime)
target_include_directories(stack_measure PRIVATE ${HS_ENGINE_INCLUDES})
target_compile_options(stack_measure PRIVATE -Os -ffast-math -fno-finite-math-only)
_hs_warnings(stack_measure)
add_test(NAME unit_stack_budget COMMAND stack_measure)

# Arena high-water gate across every effect: exits non-zero if the worst
# single-effect partition charge exceeds the device budget.
add_executable(arena_measure
  ${CMAKE_CURRENT_SOURCE_DIR}/arena_measure.cpp
  ${HS_ENGINE_SOURCES})
target_link_libraries(arena_measure PRIVATE hs_host_effect_runtime)
file(READ ${HS_ROOT}/core/platform/build_features.h _hs_build_features)
if(NOT _hs_build_features MATCHES
   "#if defined\\(ARDUINO\\)\r?\n#define HS_INLINE_PARAM_CAPACITY ([0-9]+)\r?\n")
  message(FATAL_ERROR "core/platform/build_features.h no longer states the "
    "device HS_INLINE_PARAM_CAPACITY that arena_measure models")
endif()
target_compile_definitions(arena_measure PRIVATE
  HS_INLINE_PARAM_CAPACITY=${CMAKE_MATCH_1})
target_include_directories(arena_measure PRIVATE ${HS_ENGINE_INCLUDES})
target_compile_options(arena_measure PRIVATE -Os -ffast-math -fno-finite-math-only)
_hs_warnings(arena_measure)
add_test(NAME unit_arena_budget COMMAND arena_measure)

# Render micro-benchmark at two opt levels (perf_bench.cpp). Not CTests: wall
# time is noisy. The O3 variant builds by default as a compile guard.
foreach(_opt Os O3)
  set(_perf_exclude EXCLUDE_FROM_ALL)
  if(_opt STREQUAL "O3")
    set(_perf_exclude "")
  endif()
  add_executable(perf_bench_${_opt} ${_perf_exclude}
    ${CMAKE_CURRENT_SOURCE_DIR}/perf_bench.cpp
    ${HS_ENGINE_SOURCES})
  target_link_libraries(perf_bench_${_opt} PRIVATE hs_host_effect_runtime)
  target_include_directories(perf_bench_${_opt} PRIVATE ${HS_ENGINE_INCLUDES})
  target_compile_options(perf_bench_${_opt} PRIVATE -${_opt} -ffast-math -fno-finite-math-only)
  _hs_warnings(perf_bench_${_opt})
endforeach()

# Divergence trace for perturbations of the feedback colour path
# (feedback_divergence.cpp). Cross-build comparisons use separately produced
# dumps; CTest exercises a one-frame dump and its reader.
add_executable(feedback_divergence
  ${CMAKE_CURRENT_SOURCE_DIR}/feedback_divergence.cpp
  ${HS_ROOT}/core/memory.cpp
  ${HS_ROOT}/core/engine/static_storage.cpp)
target_include_directories(feedback_divergence PRIVATE ${HS_ENGINE_INCLUDES})
target_link_libraries(feedback_divergence PRIVATE hs_host_abi)
target_compile_options(feedback_divergence PRIVATE -O2 -ffast-math -fno-finite-math-only)
_hs_warnings(feedback_divergence)
_hs_sanitize(feedback_divergence)
add_test(NAME unit_feedback_dump_a COMMAND feedback_divergence --dump
  "${CMAKE_CURRENT_BINARY_DIR}/feedback-smoke-a.bin" 1)
add_test(NAME unit_feedback_dump_b COMMAND feedback_divergence --dump
  "${CMAKE_CURRENT_BINARY_DIR}/feedback-smoke-b.bin" 1)
set_tests_properties(unit_feedback_dump_a unit_feedback_dump_b PROPERTIES
  FIXTURES_SETUP feedback_dump)
add_test(NAME unit_feedback_compare COMMAND feedback_divergence --check-comparator
  "${CMAKE_CURRENT_BINARY_DIR}/feedback-smoke-a.bin"
  "${CMAKE_CURRENT_BINARY_DIR}/feedback-smoke-b.bin")
set_tests_properties(unit_feedback_compare PROPERTIES
  FIXTURES_REQUIRED feedback_dump)

add_executable(feedback_divergence_io
  ${CMAKE_CURRENT_SOURCE_DIR}/feedback_divergence_io.cpp
  ${HS_ROOT}/core/memory.cpp
  ${HS_ROOT}/core/engine/static_storage.cpp)
target_include_directories(feedback_divergence_io PRIVATE ${HS_ENGINE_INCLUDES})
target_link_libraries(feedback_divergence_io PRIVATE hs_host_abi)
target_compile_options(feedback_divergence_io PRIVATE -O2 -ffast-math -fno-finite-math-only)
_hs_warnings(feedback_divergence_io)
_hs_sanitize(feedback_divergence_io)
add_test(NAME unit_feedback_io COMMAND feedback_divergence_io
  "${CMAKE_CURRENT_BINARY_DIR}/feedback-io.bin")

# Generator of record for core/color/srgb_decode_lut.h. EXCLUDE_FROM_ALL: run
# with no argument it rewrites the committed header; argv[1] redirects output.
add_executable(srgb_decode_gen EXCLUDE_FROM_ALL
  ${HS_ROOT}/scripts/generate_srgb_decode.cpp)
target_include_directories(srgb_decode_gen PRIVATE ${HS_ENGINE_INCLUDES})
_hs_warnings(srgb_decode_gen)


# Generator of record for tests/data/shader_chain_catalog.json; excluded from
# the default build.
add_executable(shader_chain_catalog_gen EXCLUDE_FROM_ALL
  ${CMAKE_CURRENT_SOURCE_DIR}/shader_chain_catalog_gen.cpp
  ${HS_ENGINE_SOURCES})
target_include_directories(shader_chain_catalog_gen PRIVATE ${HS_ENGINE_INCLUDES})
target_link_libraries(shader_chain_catalog_gen PRIVATE hs_test_support)
_hs_warnings(shader_chain_catalog_gen)

add_custom_target(regenerate_shader_chain_catalog
  COMMAND $<TARGET_FILE:shader_chain_catalog_gen>
          ${CMAKE_CURRENT_SOURCE_DIR}/data/shader_chain_catalog.json
  DEPENDS shader_chain_catalog_gen
  VERBATIM)

add_executable(noise_golden_gen ${HS_ROOT}/tools/noise_golden_export.cpp)
target_include_directories(noise_golden_gen PRIVATE ${HS_ENGINE_INCLUDES})
target_link_libraries(noise_golden_gen PRIVATE hs_host_abi)
_hs_warnings(noise_golden_gen)

add_executable(sdf_periodic_benchmark periodic_surface_benchmark.cpp)
target_include_directories(sdf_periodic_benchmark PRIVATE ${HS_ENGINE_INCLUDES})
target_compile_options(sdf_periodic_benchmark PRIVATE -O3 -ffast-math -fno-finite-math-only)
_hs_warnings(sdf_periodic_benchmark)

file(CONFIGURE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/spherical_experiment_check.cpp"
  CONTENT "#include \"core/engine/engine.h\"\n#include \"tools/profile_spherical_experiment.h\"\ntemplate class SphericalExperiment<96, 20>;\n")
add_library(spherical_experiment_check OBJECT EXCLUDE_FROM_ALL
  "${CMAKE_CURRENT_BINARY_DIR}/spherical_experiment_check.cpp")
target_include_directories(spherical_experiment_check PRIVATE ${HS_ENGINE_INCLUDES})
target_link_libraries(spherical_experiment_check PRIVATE hs_test_support)
_hs_warnings(spherical_experiment_check)

# Apply the native display profile to compiled targets in this directory.
get_property(_hs_dir_targets DIRECTORY PROPERTY BUILDSYSTEM_TARGETS)
set(_hs_excluded)
foreach(_hs_target IN LISTS _hs_dir_targets)
  get_target_property(_hs_target_type ${_hs_target} TYPE)
  if(NOT _hs_target_type STREQUAL "UTILITY"
      AND NOT _hs_target_type STREQUAL "INTERFACE_LIBRARY"
      AND NOT _hs_target STREQUAL "display_geometry_check")
    target_compile_definitions(${_hs_target} PRIVATE HS_DISPLAY_PROFILE=0)
  endif()
endforeach()

# Compile coverage for non-UTILITY EXCLUDE_FROM_ALL targets.
foreach(_hs_target IN LISTS _hs_dir_targets)
  get_target_property(_hs_is_excluded ${_hs_target} EXCLUDE_FROM_ALL)
  get_target_property(_hs_target_type ${_hs_target} TYPE)
  if(_hs_is_excluded AND NOT _hs_target_type STREQUAL "UTILITY")
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

add_test(NAME unit_excluded_targets COMMAND ${CMAKE_COMMAND}
  "-DCOVERED_TARGETS=$<TARGET_PROPERTY:excluded_targets,MANUALLY_ADDED_DEPENDENCIES>"
  -P ${CMAKE_CURRENT_SOURCE_DIR}/check_excluded_targets.cmake)

add_test(NAME unit_composed_edge_fade_controls COMMAND ${CMAKE_COMMAND}
  "-DCOMPILER=${CMAKE_CXX_COMPILER}" "-DROOT=${HS_ROOT}"
  "-DDEFINES=$<TARGET_PROPERTY:hs_host_abi,INTERFACE_COMPILE_DEFINITIONS>;$<TARGET_PROPERTY:hs_host_storage,INTERFACE_COMPILE_DEFINITIONS>"
  -P ${CMAKE_CURRENT_SOURCE_DIR}/check_composed_edge_fade.cmake)

# Backstop TIMEOUT for every registered test without an explicit one. Applied
# last, over the directory's final test list.
get_property(_hs_registered_tests DIRECTORY PROPERTY TESTS)
foreach(_hs_test IN LISTS _hs_registered_tests)
  get_test_property(${_hs_test} TIMEOUT _hs_existing_timeout)
  if(NOT _hs_existing_timeout)
    set_tests_properties(${_hs_test} PROPERTIES TIMEOUT ${_hs_test_timeout})
  endif()
endforeach()
