# Syntax-only compile contracts. GROUP selects one table; COMPILER, ROOT and
# the DEFINES list give the compile context. The companions group compiles the
# valid configurations of every fixture, sharing one translation unit wherever
# the configurations' macros are disjoint; each refusal compiles alone.
cmake_minimum_required(VERSION 3.29)

set(EDGE_FADE composed_edge_fade_check.cpp)
set(HOOKS choreography_hooks_check.cpp)
set(METADATA composed_metadata_check.cpp)
set(RING ring_column_pipeline_check.cpp)

# Compile the ROOT/tests fixtures in sources as one translation unit with the
# extra flags in ARGN. expected is PASS for a configuration that must compile,
# otherwise a regex the refusal must match; label names the contract in
# failure messages.
function(hs_check_compile_fixture sources label expected)
  set(flags)
  foreach(definition IN LISTS DEFINES)
    list(APPEND flags "-D${definition}")
  endforeach()
  list(POP_BACK sources main)
  foreach(source IN LISTS sources)
    list(APPEND flags -include "${ROOT}/tests/${source}")
  endforeach()
  execute_process(COMMAND "${COMPILER}" -std=gnu++20 -fsyntax-only
    "-I${ROOT}" ${flags} ${ARGN}
    "${ROOT}/tests/${main}"
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
  if(expected STREQUAL "PASS")
    if(NOT result EQUAL 0)
      message(FATAL_ERROR "valid ${label} companion rejected: ${ARGN}\n${output}${error}")
    endif()
  elseif(result EQUAL 0)
    message(FATAL_ERROR "invalid ${label} fixture accepted: ${ARGN}")
  elseif(NOT "${output}${error}" MATCHES "${expected}")
    message(FATAL_ERROR "wrong ${label} refusal: ${output}${error}")
  endif()
endfunction()

function(check_edge_fade expected)
  hs_check_compile_fixture(${EDGE_FADE} edge-fade "${expected}" ${ARGN})
endfunction()

function(check_hooks expected)
  hs_check_compile_fixture(${HOOKS} "hook contract" "${expected}" ${ARGN})
endfunction()

function(check_metadata expected)
  hs_check_compile_fixture(${METADATA} metadata "${expected}" ${ARGN})
endfunction()

if(GROUP STREQUAL "companions")
  hs_check_compile_fixture("${EDGE_FADE};${HOOKS};${METADATA};${RING}"
    "fixture" PASS)
  hs_check_compile_fixture("${EDGE_FADE};${HOOKS};${RING}"
    "edge-fade, hook contract or ring_columns" PASS
    -DHS_EDGE_FADE_NO_DISTANCE -DHS_EDGE_FADE_DISABLED
    -DHS_HOOK_ROW_FIELDS -DHS_RING_COLUMN_DIRECT)
  hs_check_compile_fixture("${EDGE_FADE};${HOOKS}"
    "edge-fade or hook contract" PASS
    -DHS_EDGE_FADE_AUTHORED -DHS_HOOK_PRESET_ROW)
  check_edge_fade(PASS -DHS_EDGE_FADE_AUTHORED -DHS_EDGE_FADE_NO_DISTANCE
    -DHS_EDGE_FADE_DISABLED)
  message(STATUS "compile contract companions passed")
elseif(GROUP STREQUAL "edge-fade")
  check_edge_fade("edge-fade coverage requires projection edge distance"
    -DHS_EDGE_FADE_NO_DISTANCE)
  check_edge_fade("edge-fade coverage requires projection edge distance"
    -DHS_EDGE_FADE_AUTHORED -DHS_EDGE_FADE_NO_DISTANCE)
  message(STATUS "composed edge-fade controls passed")
elseif(GROUP STREQUAL "hooks")
  check_hooks("preset row is outside its parameter fields"
    -DHS_HOOK_ROW_FIELDS -DHS_HOOK_ROW_OUT_OF_RANGE)
  check_hooks("preset.index. requires PRESET_IDS"
    -DHS_HOOK_PRESET_ROW -DHS_HOOK_PRESET_ROW_NO_IDS)
  check_hooks("choreography hooks: wrong hook signature" -DHS_HOOK_BOUND_TYPO)
  check_hooks("choreography hooks: wrong hook signature" -DHS_HOOK_BAD_RETURN)
  check_hooks("choreography hooks: wrong hook signature" -DHS_HOOK_BAD_ARGUMENT)
  check_hooks("choreography hooks: bound inherited default"
    -DHS_HOOK_INHERITED_DEFAULT)
  check_hooks("choreography hooks: bound inherited default"
    -DHS_HOOK_VALIDATION_DEFAULT)
  check_hooks("choreography hooks: duplicate binding" -DHS_HOOK_DUPLICATE)
  check_hooks("choreography hooks: declaration owner mismatch"
    -DHS_HOOK_FACTORY_TYPO)
  check_hooks("effect lacks the departure.s hook" -DHS_HOOK_FADE_NO_OPACITY)
  check_hooks("effect lacks the departure.s hook" -DHS_HOOK_LERP_NO_BLEND)
  check_hooks("departure outlasts PRESET_DWELL_FRAMES"
    -DHS_HOOK_DEPARTURE_OUTLASTS)
  check_hooks("discrete fields require MIDPOINT or SNAP curves"
    -DHS_HOOK_DISCRETE_LERP)
  check_hooks("PRESET_IDS and PRESETS must name the same presets"
    -DHS_HOOK_IDS_SIZE_MISMATCH)
  message(STATUS "choreography hook contracts passed")
elseif(GROUP STREQUAL "metadata")
  check_metadata("declare topology through typed pipeline entries"
    -DHS_METADATA_REDUNDANT_SCALAR)
  check_metadata("unsupported projection coverage policy"
    -DHS_METADATA_UNKNOWN_COVERAGE)
  check_metadata("projection stage is required" -DHS_METADATA_MISSING_PROJECT)
  check_metadata("displacement modulation requires pipeline path tracking"
    -DHS_METADATA_UNTRACKED_DISPLACEMENT)
  check_metadata("composed modulation presets must retain source and mechanic"
    -DHS_METADATA_MODULATION_SOURCE_CHANGE)
  check_metadata("composed modulation presets must retain source and mechanic"
    -DHS_METADATA_MODULATION_MECHANIC_CHANGE)
  message(STATUS "composed metadata contracts passed")
elseif(GROUP STREQUAL "ring-columns")
  hs_check_compile_fixture(${RING} ring_columns
    "the ring column walk requires a screen-only pipeline"
    -DHS_RING_COLUMN_WORLD)
  message(STATUS "ring column pipeline contracts passed")
else()
  message(FATAL_ERROR "unknown compile contract GROUP '${GROUP}'")
endif()
