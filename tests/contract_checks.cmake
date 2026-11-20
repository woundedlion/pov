# Native engine and compile-contract checks.

# One executable per engine compile context. Each MODULES entry registers
# unit_<module>, which passes the module name as the executable's argument;
# without MODULES, unit_<name minus _check> runs it with no argument.
function(_hs_engine_check name source)
  cmake_parse_arguments(PARSE_ARGV 2 check "" "" "LIBS;DEFINES;MODULES;OPTIONS")
  if(NOT IS_ABSOLUTE "${source}")
    set(source "${CMAKE_CURRENT_SOURCE_DIR}/${source}")
  endif()
  add_executable(${name} "${source}" ${HS_ENGINE_SOURCES})
  target_compile_definitions(${name} PRIVATE ${check_DEFINES})
  target_compile_options(${name} PRIVATE ${check_OPTIONS})
  target_link_libraries(${name} PRIVATE ${check_LIBS})
  target_include_directories(${name} PRIVATE ${HS_ENGINE_INCLUDES})
  _hs_warnings(${name})
  _hs_sanitize(${name})
  if(check_MODULES)
    foreach(module IN LISTS check_MODULES)
      add_test(NAME unit_${module} COMMAND ${name} ${module})
    endforeach()
  else()
    string(REGEX REPLACE "_check$" "" test_suffix "${name}")
    add_test(NAME unit_${test_suffix} COMMAND ${name})
  endif()
endfunction()

_hs_engine_check(dma_led_guards_check "${_dma_led_source}"
  LIBS hs_host_storage
  DEFINES USE_DMA_LEDS)

# NaN-aware clamps under the release math flags.
_hs_engine_check(fastmath_clamp_check fastmath_clamp_check.cpp
  LIBS hs_host_abi
  DEFINES HS_ENABLE_TEST_ORACLES=1 HS_TEST_FAST_MATH=1
  OPTIONS -O3 -ffast-math -fno-finite-math-only)

# Re-derive committed relax bakes with IEEE float semantics.
_hs_engine_check(relax_bake_verify "${HS_ROOT}/tools/relax_bake_harness.cpp"
  LIBS hs_host_abi
  DEFINES HS_RELAX_BAKE_VERIFY
  OPTIONS -fno-fast-math)

# Legacy offset mapping needs a separate TU: offset-3 and offset-0
# PhiLUT/TrigLUT instantiations have incompatible definitions.
_hs_engine_check(h_offset_check h_offset_check.cpp
  LIBS hs_host_storage
  DEFINES HS_TEST_H_OFFSET=3 HS_ENABLE_TEST_HOOKS=1)

# Global operator new replacement conflicts with the sanitizer runtimes'.
if(NOT HS_SANITIZE AND NOT HS_TSAN)
  add_executable(frame_heap_check
    ${CMAKE_CURRENT_SOURCE_DIR}/frame_heap_check.cpp
    ${HS_ENGINE_SOURCES})
  target_link_libraries(frame_heap_check PRIVATE hs_host_effect_runtime)
  target_include_directories(frame_heap_check PRIVATE ${HS_ENGINE_INCLUDES})
  _hs_warnings(frame_heap_check)
  add_test(NAME unit_frame_heap COMMAND frame_heap_check)
endif()


_hs_engine_check(display_geometry_check display_geometry_check.cpp
  LIBS hs_host_storage
  DEFINES HS_DISPLAY_PROFILE=1 HS_ENABLE_TEST_HOOKS=1)

_hs_engine_check(runtime_display_geometry_check runtime_display_geometry_check.cpp
  LIBS hs_host_effect_runtime
  DEFINES HS_RUNTIME_DISPLAY_GEOMETRY=1 HS_ENABLE_TEST_HOOKS=1
  MODULES runtime_display_geometry ashcloud_camera)

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
  target_compile_definitions(perf_bench_${_opt} PRIVATE HS_ENABLE_TEMPORAL_SWITCH=1)
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

add_test(NAME unit_feedback_cli COMMAND "${CMAKE_COMMAND}"
  "-DBIN=$<TARGET_FILE:feedback_divergence>"
  "-DCAPTURE=${CMAKE_CURRENT_BINARY_DIR}/feedback-smoke-a.bin"
  -P "${CMAKE_CURRENT_SOURCE_DIR}/check_feedback_cli.cmake")
set_tests_properties(unit_feedback_cli PROPERTIES FIXTURES_REQUIRED feedback_dump)

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

add_executable(noise_golden_gen ${HS_ROOT}/tools/noise_golden_export.cpp)
target_include_directories(noise_golden_gen PRIVATE ${HS_ENGINE_INCLUDES})
target_link_libraries(noise_golden_gen PRIVATE hs_host_abi)
_hs_warnings(noise_golden_gen)

# Measures query budgets for Raycast::surface_search on periodic surfaces.
add_executable(sdf_periodic_benchmark EXCLUDE_FROM_ALL periodic_surface_benchmark.cpp)
target_include_directories(sdf_periodic_benchmark PRIVATE ${HS_ENGINE_INCLUDES})
target_compile_options(sdf_periodic_benchmark PRIVATE -O3 -ffast-math -fno-finite-math-only)
_hs_warnings(sdf_periodic_benchmark)

# Apply the native display profile to compiled targets in this directory.
set(HS_NATIVE_DISPLAY_PROFILE 0 CACHE STRING "Native test latitude profile (0 or 1)")
set_property(CACHE HS_NATIVE_DISPLAY_PROFILE PROPERTY STRINGS 0 1)
if(NOT HS_NATIVE_DISPLAY_PROFILE MATCHES "^[01]$")
  message(FATAL_ERROR "HS_NATIVE_DISPLAY_PROFILE must be 0 or 1")
endif()
get_property(_hs_dir_targets DIRECTORY PROPERTY BUILDSYSTEM_TARGETS)
set(_hs_excluded)
foreach(_hs_target IN LISTS _hs_dir_targets)
  get_target_property(_hs_target_type ${_hs_target} TYPE)
  if(NOT _hs_target_type STREQUAL "UTILITY"
      AND NOT _hs_target_type STREQUAL "INTERFACE_LIBRARY"
      AND NOT _hs_target STREQUAL "display_geometry_check")
    target_compile_definitions(${_hs_target} PRIVATE HS_DISPLAY_PROFILE=${HS_NATIVE_DISPLAY_PROFILE})
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

function(_hs_compile_contract_test name group)
  add_test(NAME ${name} COMMAND ${CMAKE_COMMAND}
    "-DCOMPILER=${CMAKE_CXX_COMPILER}" "-DROOT=${HS_ROOT}" "-DGROUP=${group}"
    "-DDEFINES=$<TARGET_PROPERTY:hs_host_abi,INTERFACE_COMPILE_DEFINITIONS>;$<TARGET_PROPERTY:hs_host_storage,INTERFACE_COMPILE_DEFINITIONS>"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/check_compile_contracts.cmake)
endfunction()

_hs_compile_contract_test(unit_compile_contract_companions companions)

_hs_compile_contract_test(unit_composed_edge_fade_controls edge-fade)

_hs_compile_contract_test(unit_choreography_hooks hooks)

_hs_compile_contract_test(unit_composed_metadata metadata)

_hs_compile_contract_test(unit_ring_column_pipeline ring-columns)
