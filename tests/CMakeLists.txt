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
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
endif()

# Python is required; optional generator dependencies may skip locally.
option(HS_REQUIRE_GENERATORS
    "Fail rather than skip when optional generator dependencies are missing" OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# The native-clang toolchain may fall back to an unversioned `clang` from PATH;
# refuse anything older than 18.
if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  message(FATAL_ERROR
    "The native test suite requires upstream LLVM Clang 18 or newer; point "
    "CMAKE_CXX_COMPILER at LLVM clang++, or use the emsdk's compiler.")
endif()
if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 18)
  message(FATAL_ERROR
    "Clang ${CMAKE_CXX_COMPILER_VERSION} is below the 18 this suite is built "
    "against. Point CMAKE_CXX_COMPILER at a newer clang++, or let "
    "cmake/toolchain-native-clang.cmake resolve the emsdk's.")
endif()

# Clang < 19 needs relaxed template-template matching opted in for the
# harness's `template <int, int> class E` parameters; Clang 20 removed the flag.
if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 19)
  add_compile_options(-frelaxed-template-template-args)
endif()

set(HS_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/..)
include(${HS_ROOT}/cmake/engine_sources.cmake)
set(HS_ENGINE_INCLUDES ${HS_ROOT})

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

add_library(hs_test_support INTERFACE)
target_link_libraries(hs_test_support INTERFACE hs_host_storage)
include(${HS_ROOT}/cmake/test_capabilities.cmake)
target_compile_definitions(hs_test_support INTERFACE
  ${HS_TEST_CAPABILITY_DEFINITIONS}
  HS_ENABLE_STRUCTURAL_AUDITS=1)

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
  "HS_ENGINE_CATALOG_PATH=\"${HS_ROOT}/scripts/engine_catalog.json\""
  "HS_PROMOTED_PATTERNS_DIR=\"${HS_ROOT}/patterns\"")
target_link_libraries(run_tests PRIVATE hs_test_support)

set(_dma_led_source "${CMAKE_CURRENT_BINARY_DIR}/dma_led_guard_check.cpp")
file(CONFIGURE OUTPUT "${_dma_led_source}" CONTENT [=[
#include "core/platform/led.h"
int main() {
  if (correction_guard_live()) return 1;
  { NoColorCorrection guard; if (!correction_guard_live()) return 2; }
  if (correction_guard_live()) return 3;
  { NoTempCorrection guard; if (!correction_guard_live()) return 4; }
  return correction_guard_live() ? 5 : 0;
}
]=])
add_executable(dma_led_guard_check "${_dma_led_source}" ${HS_ENGINE_SOURCES})
target_include_directories(dma_led_guard_check PRIVATE ${HS_ROOT})
target_compile_definitions(dma_led_guard_check PRIVATE USE_DMA_LEDS)
target_link_libraries(dma_led_guard_check PRIVATE hs_host_storage)
_hs_warnings(dma_led_guard_check)
add_test(NAME unit_dma_led_guards COMMAND dma_led_guard_check)

set(_color_header_sources)
file(GLOB _color_headers CONFIGURE_DEPENDS "${HS_ROOT}/core/color/*.h")
list(FILTER _color_headers EXCLUDE REGEX "/(color_luts|gamut_lut|srgb_decode_lut)\\.h$")
foreach(_header_path IN LISTS _color_headers)
  get_filename_component(_header "${_header_path}" NAME_WE)
  set(_source "${CMAKE_CURRENT_BINARY_DIR}/color_header_${_header}.cpp")
  set(_content "#include \"core/color/${_header}.h\"\n")
  if(_header MATCHES "^(modulation|modulation_types|hue_rotation_lut)$")
    string(APPEND _content "#ifdef FASTNOISELITE_H\n#error generic color header imports noise generator\n#endif\n")
  endif()
  file(CONFIGURE OUTPUT "${_source}" CONTENT "${_content}")
  list(APPEND _color_header_sources "${_source}")
endforeach()
add_library(color_header_checks OBJECT ${_color_header_sources})
target_include_directories(color_header_checks PRIVATE ${HS_ROOT})
_hs_warnings(color_header_checks)
add_dependencies(run_tests color_header_checks)

set(_math_header_sources)
file(GLOB _math_headers CONFIGURE_DEPENDS "${HS_ROOT}/core/math/*.h")
foreach(_header_path IN LISTS _math_headers)
  get_filename_component(_header "${_header_path}" NAME_WE)
  set(_source "${CMAKE_CURRENT_BINARY_DIR}/math_header_${_header}.cpp")
  set(_content "#include \"core/math/${_header}.h\"\n")
  if(_header STREQUAL "math_constants")
    string(APPEND _content "#if defined(HS_MATH_SCALAR_NOISE_H) || defined(FASTNOISELITE_H)\n#error math constants header imports noise\n#endif\n")
  endif()
  file(CONFIGURE OUTPUT "${_source}" CONTENT "${_content}")
  list(APPEND _math_header_sources "${_source}")
endforeach()
add_library(math_header_checks OBJECT ${_math_header_sources})
target_include_directories(math_header_checks PRIVATE ${HS_ROOT})
_hs_warnings(math_header_checks)
add_dependencies(run_tests math_header_checks)

set(_pattern_header_sources)
file(GLOB _pattern_headers CONFIGURE_DEPENDS "${HS_ROOT}/core/pattern/*.h")
foreach(_header_path IN LISTS _pattern_headers)
  get_filename_component(_header "${_header_path}" NAME_WE)
  set(_source "${CMAKE_CURRENT_BINARY_DIR}/pattern_header_${_header}.cpp")
  set(_content "#include \"core/pattern/${_header}.h\"\n")
  if(NOT _header STREQUAL "noise")
    string(APPEND _content "#ifdef FASTNOISELITE_H\n#error analytic pattern header imports noise generator\n#endif\n")
  endif()
  file(CONFIGURE OUTPUT "${_source}" CONTENT "${_content}")
  list(APPEND _pattern_header_sources "${_source}")
endforeach()
add_library(pattern_header_checks OBJECT ${_pattern_header_sources})
target_include_directories(pattern_header_checks PRIVATE ${HS_ROOT})
_hs_warnings(pattern_header_checks)
add_dependencies(run_tests pattern_header_checks)

add_executable(pattern_qualification ${HS_ROOT}/tools/pattern_qualification.cpp)
target_include_directories(pattern_qualification PRIVATE ${HS_ROOT})
_hs_warnings(pattern_qualification)
add_test(NAME unit_pattern_qualification COMMAND pattern_qualification)


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
  DEPENDS ${_guard_sources} guard_directories.txt count_guard_sites.cmake header_sections.cmake death_guard_sites.h.in
  VERBATIM)
add_custom_target(death_guard_census DEPENDS "${_guard_output}")
add_dependencies(run_tests death_guard_census)
add_test(NAME unit_death_pins COMMAND "${Python3_EXECUTABLE}"
         "${CMAKE_CURRENT_SOURCE_DIR}/check_death_pins.py" "${HS_ROOT}")

add_test(NAME unit_include_paths COMMAND "${Python3_EXECUTABLE}"
         "${HS_ROOT}/tools/include_paths.py" --root "${HS_ROOT}")

# Headers are included as "tests/..." and "core/..." from the repo root; the
# binary dir carries the generated guard census.
target_include_directories(run_tests PRIVATE ${HS_ENGINE_INCLUDES}
                           ${CMAKE_CURRENT_BINARY_DIR})

_hs_warnings(run_tests)

add_executable(case_context_check case_context_check.cpp)
target_include_directories(case_context_check PRIVATE ${HS_ROOT})
_hs_warnings(case_context_check)
add_test(NAME unit_case_context COMMAND ${CMAKE_COMMAND}
    -DTEST_BIN=$<TARGET_FILE:case_context_check>
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_case_context.cmake)

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
      -fsanitize=thread -fno-omit-frame-pointer ${HS_DEBUG_INFO_FLAG})
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
      -fno-omit-frame-pointer ${HS_DEBUG_INFO_FLAG})
  else()
    target_compile_options(${target} PRIVATE
      -fsanitize=address,undefined -fno-sanitize-recover=all ${_ignore}
      -fno-omit-frame-pointer ${HS_DEBUG_INFO_FLAG})
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
_hs_sanitize(case_context_check)

include(${CMAKE_CURRENT_SOURCE_DIR}/test_targets.cmake)

add_executable(mindsplatter_palette_gen
  ${HS_ROOT}/tools/mindsplatter_palette_gen.cpp
  ${HS_ENGINE_SOURCES})
target_include_directories(mindsplatter_palette_gen PRIVATE ${HS_ENGINE_INCLUDES})
target_link_libraries(mindsplatter_palette_gen PRIVATE hs_host_abi)
_hs_warnings(mindsplatter_palette_gen)
_hs_sanitize(mindsplatter_palette_gen)

add_custom_target(regenerate_mindsplatter_palette
  COMMAND $<TARGET_FILE:mindsplatter_palette_gen>
          ${HS_ROOT}/effects/generated/mind_splatter_palette_bank.h
  DEPENDS mindsplatter_palette_gen
  VERBATIM)

add_test(NAME unit_mindsplatter_palette_form COMMAND ${CMAKE_COMMAND}
  -DGENERATOR=$<TARGET_FILE:mindsplatter_palette_gen>
  -DCOMMITTED=${HS_ROOT}/effects/generated/mind_splatter_palette_bank.h
  -DGENERATED=${CMAKE_CURRENT_BINARY_DIR}/mind_splatter_palette_bank.generated.h
  -P ${CMAKE_CURRENT_SOURCE_DIR}/check_mindsplatter_palette_form.cmake)

# Standalone tolerance check for the palette data, with its own error report.
add_executable(mindsplatter_palette_check
  ${CMAKE_CURRENT_SOURCE_DIR}/mindsplatter_palette_check.cpp
  ${HS_ENGINE_SOURCES})
target_include_directories(mindsplatter_palette_check PRIVATE ${HS_ENGINE_INCLUDES})
target_link_libraries(mindsplatter_palette_check PRIVATE hs_host_abi)
_hs_warnings(mindsplatter_palette_check)
_hs_sanitize(mindsplatter_palette_check)
add_test(NAME unit_mindsplatter_palette COMMAND mindsplatter_palette_check)

# Opt-in Clang source-based coverage for the shared native suite. Instrumentation changes
# codegen. Needs clang_rt.profile, which the Windows emsdk clang does not ship.
option(HS_COVERAGE "Instrument the native test runners, case objects and their engine objects for llvm-cov source-based coverage" OFF)
if(HS_COVERAGE)
  foreach(_target IN LISTS _hs_module_targets _hs_case_targets)
    target_compile_options(${_target} PRIVATE -fprofile-instr-generate -fcoverage-mapping)
    target_link_options(${_target} PRIVATE -fprofile-instr-generate)
  endforeach()
endif()

enable_testing()

# One CTest test per module so CI shards (and local `ctest -R <name>`) can
# target subsets without re-running everything.
# Bound each module independently.
set(_hs_test_timeout 600)
if(HS_SANITIZE OR HS_TSAN)
  set(_hs_test_timeout 1800)
endif()
foreach(_mod IN LISTS _hs_test_modules)
    if(_mod STREQUAL "death")
        continue()
    endif()
    add_test(NAME unit_${_mod} COMMAND run_tests ${_mod})
    set_tests_properties(unit_${_mod} PROPERTIES TIMEOUT ${_hs_test_timeout})
endforeach()

add_test(NAME unit_shapeshifter_footprint COMMAND test_shapeshifter shapeshifter_oracle)
set_tests_properties(unit_shapeshifter_footprint PROPERTIES
    TIMEOUT ${_hs_test_timeout} ENVIRONMENT "HS_SURFACE_DARK_DIAGNOSTICS=1;CI=")

add_test(NAME unit_shapeshifter_lobes COMMAND test_shapeshifter shapeshifter_oracle)
set_tests_properties(unit_shapeshifter_lobes PROPERTIES
    TIMEOUT ${_hs_test_timeout} ENVIRONMENT "HS_SURFACE_LOBE_DIAGNOSTICS=1;CI=")

set_tests_properties(unit_canvas PROPERTIES
    ENVIRONMENT "HS_BUFFER_FREE_WATCHDOG_US=30000000")

# The death module runs as unit_death_<k> shards, each spawning the cases whose
# table index % _hs_death_shards is k and writing a record, plus
# unit_death_census, which merges the records, requires every enabled case to
# have run exactly once, and reports guard coverage. `run_tests death` with no
# HS_DEATH_SHARD runs the whole module in one process.
set(_hs_death_shards 8)
math(EXPR _hs_death_last_shard "${_hs_death_shards} - 1")
foreach(_shard RANGE ${_hs_death_last_shard})
    add_test(NAME unit_death_${_shard} COMMAND run_tests death)
    set_tests_properties(unit_death_${_shard} PROPERTIES
        TIMEOUT ${_hs_test_timeout}
        FIXTURES_SETUP death_shards
        ENVIRONMENT "HS_DEATH_SHARD=${_shard}/${_hs_death_shards};HS_DEATH_SHARD_DIR=${CMAKE_CURRENT_BINARY_DIR}")
endforeach()
add_test(NAME unit_death_census COMMAND run_tests death)
set_tests_properties(unit_death_census PROPERTIES
    TIMEOUT ${_hs_test_timeout}
    FIXTURES_REQUIRED death_shards
    ENVIRONMENT "HS_DEATH_SHARD=census/${_hs_death_shards};HS_DEATH_SHARD_DIR=${CMAKE_CURRENT_BINARY_DIR}")
if(NOT HS_SANITIZE AND NOT HS_TSAN)
    set_tests_properties(unit_effects_smoke PROPERTIES TIMEOUT 1200)
endif()

# Under ASan and UBSan the production-resolution effect sweep runs as roster
# shards beside the rest of the suite instead of one serial pass after it.
if(HS_SANITIZE)
    set(_hs_smoke_shards 2)
    math(EXPR _hs_smoke_last_shard "${_hs_smoke_shards} - 1")
    foreach(_shard RANGE ${_hs_smoke_last_shard})
        add_test(NAME unit_effects_smoke_full_${_shard} COMMAND run_tests effects_smoke)
        set_tests_properties(unit_effects_smoke_full_${_shard} PROPERTIES
            TIMEOUT ${_hs_test_timeout} COST 240
            ENVIRONMENT "HS_EFFECTS_FULL=1;HS_REQUIRE_EFFECTS_FULL=1;HS_EFFECTS_SHARD=${_shard}/${_hs_smoke_shards}")
    endforeach()
endif()

add_test(NAME unit_death_child_dispatch COMMAND ${CMAKE_COMMAND}
    -DTEST_BIN=$<TARGET_FILE:run_tests>
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_death_child_dispatch.cmake)

# run_tests --check-modules fails if _hs_test_modules and HS_TEST_MODULE_LIST
# diverge in either direction.
add_test(NAME unit_module_roster
    COMMAND run_tests --check-modules ${_hs_test_modules})

# Pin the compiled groups' module includes against the dispatcher roster.
add_test(NAME unit_module_includes COMMAND ${CMAKE_COMMAND}
    "-DGROUP_SOURCES=${_hs_test_group_sources}"
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

add_test(NAME unit_effect_case_parity COMMAND ${CMAKE_COMMAND}
    -DTESTS_DIR=${CMAKE_CURRENT_SOURCE_DIR}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_effect_case_parity.cmake)
add_test(NAME unit_effect_case_parity_controls COMMAND ${CMAKE_COMMAND}
    -DGATE=${CMAKE_CURRENT_SOURCE_DIR}/check_effect_case_parity.cmake
    -DWORK=${CMAKE_CURRENT_BINARY_DIR}/effect_case_parity_fixture
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_effect_case_parity_selftest.cmake)

# Check reachability and assertion reach for column-0 void/bool/int/size_t
# cases in tests/*.h and tests/*.hpp, including run_*_cases drivers. Death
# case_* registration is checked by check_death_pins.py.
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
add_test(NAME unit_linker_patches_controls COMMAND ${CMAKE_COMMAND}
    -DSCRIPT=${HS_ROOT}/tools/phantasm.ld
    -DPLATFORMIO_INI=${HS_ROOT}/platformio.ini
    -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}/linker_patches_controls
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_linker_patches_selftest.cmake)

add_test(NAME unit_mirrored_prune COMMAND ${CMAKE_COMMAND}
    -DTEST_ROOT=${CMAKE_CURRENT_BINARY_DIR}/mirrored_prune_fixture
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_mirrored_prune.cmake)

add_test(NAME unit_chain_capture_fixtures COMMAND "${Python3_EXECUTABLE}"
  "${CMAKE_CURRENT_SOURCE_DIR}/../tools/profile_tests/test_chain_capture_fixtures.py")

set(PULLBACK_MANIFEST_DIR ${HS_ROOT}/tests/data/pullback)
set(PULLBACK_MANIFEST_HEADER
    ${CMAKE_CURRENT_BINARY_DIR}/pullback_manifest.generated.h)
file(GLOB PULLBACK_MANIFEST_FILES CONFIGURE_DEPENDS
     ${PULLBACK_MANIFEST_DIR}/*.json)
add_custom_command(
  OUTPUT ${PULLBACK_MANIFEST_HEADER}
  COMMAND ${Python3_EXECUTABLE} ${HS_ROOT}/tools/generate_pullback_manifest_header.py
          --manifest-dir ${PULLBACK_MANIFEST_DIR}
          --output ${PULLBACK_MANIFEST_HEADER}
  DEPENDS
    ${HS_ROOT}/tools/generate_pullback_manifest_header.py
    ${HS_ROOT}/tools/pullback_operations.def
    ${PULLBACK_MANIFEST_FILES}
  VERBATIM)
add_custom_target(pullback_manifest_header DEPENDS ${PULLBACK_MANIFEST_HEADER})

add_executable(pullback_manifest_check
  ${CMAKE_CURRENT_SOURCE_DIR}/pullback_manifest_check.cpp
  ${PULLBACK_MANIFEST_HEADER})
add_dependencies(pullback_manifest_check pullback_manifest_header)
target_include_directories(pullback_manifest_check PRIVATE
                           ${CMAKE_CURRENT_BINARY_DIR} ${HS_ENGINE_INCLUDES})
_hs_warnings(pullback_manifest_check)
_hs_sanitize(pullback_manifest_check)
add_test(NAME unit_pullback_manifest COMMAND pullback_manifest_check)


add_executable(pullback_capture_native
  ${HS_ROOT}/tools/pullback_capture_native.cpp
  ${HS_ENGINE_SOURCES})
include(${HS_ROOT}/cmake/chain_capture_fixtures.cmake)
hs_chain_capture_fixtures(pullback_capture_native "${HS_ROOT}")
target_include_directories(pullback_capture_native PRIVATE
  ${HS_ROOT})
hs_capture_feature_defines(pullback_capture_native)
target_link_libraries(pullback_capture_native PRIVATE hs_host_storage)
_hs_warnings(pullback_capture_native)
_hs_sanitize(pullback_capture_native)
add_test(NAME unit_pullback_capture_determinism
  COMMAND ${CMAKE_COMMAND}
    -DPRODUCER=$<TARGET_FILE:pullback_capture_native>
    -DOUTPUT_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DPYTHON=${Python3_EXECUTABLE}
    -DCAPTURE_SCRIPT=${HS_ROOT}/tools/pullback_capture.py
    -DMANIFEST_DIR=${PULLBACK_MANIFEST_DIR}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_pullback_capture.cmake)
# Pin core/color/color_luts.h against its generator.
add_test(NAME unit_color_luts COMMAND ${CMAKE_COMMAND}
    -DPYTHON_EXE=${Python3_EXECUTABLE}
    -DGENERATOR=${HS_ROOT}/scripts/generate_luts.py
    -DCOMMITTED=${HS_ROOT}/core/color/color_luts.h
    -DGENERATED=${CMAKE_CURRENT_BINARY_DIR}/color_luts.generated.h
    -DSKIP_CODE=77 -DREQUIRE_GENERATORS=${HS_REQUIRE_GENERATORS}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_color_luts.cmake)
set_tests_properties(unit_color_luts PROPERTIES SKIP_RETURN_CODE 77)

# Host self-tests for tools/teensy_size_trail.py: the ELF32 section-header
# parser and the regression classifier.
add_test(NAME unit_teensy_size_trail COMMAND ${CMAKE_COMMAND}
    -DPYTHON_EXE=${Python3_EXECUTABLE}
    -DHS_ROOT=${HS_ROOT}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_teensy_size_trail.cmake)

# Pin core/color/gamut_lut.h against its generator, the generator's copies of the
# OKLab matrices and the gamut slack against core/color/color_space.h, and its
# diamond angle against core/math/fast_math.h.
add_test(NAME unit_gamut_lut COMMAND ${CMAKE_COMMAND}
    -DPYTHON_EXE=${Python3_EXECUTABLE}
    -DGENERATOR=${HS_ROOT}/tools/gen_gamut_lut.py
    -DCACHE_FILE=${CMAKE_CURRENT_BINARY_DIR}/gamut_lut_check.sha256
    -DSKIP_CODE=77 -DREQUIRE_GENERATORS=${HS_REQUIRE_GENERATORS}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_gamut_lut.cmake)
set_tests_properties(unit_gamut_lut PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 900)

add_test(NAME unit_tangent_field_constants
    COMMAND ${CMAKE_COMMAND}
    -DPYTHON_EXE=${Python3_EXECUTABLE}
    -DGENERATOR=${HS_ROOT}/tools/gen_tangent_field_constants.py
    -DSKIP_CODE=77 -DREQUIRE_GENERATORS=${HS_REQUIRE_GENERATORS}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_tangent_field_constants.cmake)
set_tests_properties(unit_tangent_field_constants PROPERTIES SKIP_RETURN_CODE 77 TIMEOUT 300)

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
target_compile_definitions(fastmath_clamp_check PRIVATE
  HS_ENABLE_TEST_ORACLES=1 HS_TEST_FAST_MATH=1)
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
    -DPYTHON_EXE=${Python3_EXECUTABLE}
    -DSCRIPT=${HS_ROOT}/tools/relax_bakes.py
    -DHARNESS=$<TARGET_FILE:relax_bake_gen>
    -DDUMP=${CMAKE_CURRENT_BINARY_DIR}/relax_bakes.dump
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_relax_bakes.cmake)

include(${CMAKE_CURRENT_SOURCE_DIR}/contract_checks.cmake)

# Backstop TIMEOUT for every registered test without an explicit one. Applied
# last, over the directory's final test list.
get_property(_hs_registered_tests DIRECTORY PROPERTY TESTS)
foreach(_hs_test IN LISTS _hs_registered_tests)
  get_test_property(${_hs_test} TIMEOUT _hs_existing_timeout)
  if(NOT _hs_existing_timeout)
    set_tests_properties(${_hs_test} PROPERTIES TIMEOUT ${_hs_test_timeout})
  endif()
endforeach()

# Link startup into every native executable, including helper tools.
foreach(_hs_target IN LISTS _hs_dir_targets)
  get_target_property(_hs_target_type ${_hs_target} TYPE)
  if(_hs_target_type STREQUAL "EXECUTABLE")
    target_link_libraries(${_hs_target} PRIVATE hs_test_process)
  endif()
endforeach()

include(${CMAKE_CURRENT_SOURCE_DIR}/engine_libraries.cmake)
hs_share_test_engines()

# ctest starts higher-COST tests first, so the longest ones do not finish last.
set_tests_properties(unit_effects PROPERTIES COST 180)
set_tests_properties(unit_effects_smoke unit_choreography_hooks
    unit_shapeshifter_lobes unit_shapeshifter_footprint PROPERTIES COST 80)
set_tests_properties(unit_case_calls unit_conway_soak unit_shapeshifter_oracle
    unit_stack_budget unit_arena_budget unit_pullback PROPERTIES COST 50)

include(${CMAKE_CURRENT_SOURCE_DIR}/header_contracts.cmake)
