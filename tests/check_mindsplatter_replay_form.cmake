cmake_minimum_required(VERSION 3.29)

function(check_generated_form)
  file(SIZE "${GENERATED}" _generated_size)
  if(_generated_size LESS 100000)
    message(FATAL_ERROR
      "MindSplatter replay generator emitted only ${_generated_size} bytes")
  endif()

  foreach(_source IN ITEMS _committed _generated)
    string(FIND "${${_source}}" "HEAVY_SEARCH_V1_STATE[]" _state_offset)
    if(_state_offset LESS 0)
      message(FATAL_ERROR "MindSplatter replay corpus omits the state declaration")
    endif()
    string(SUBSTRING "${${_source}}" 0 ${_state_offset} ${_source}_preamble)
  endforeach()
  if(NOT _committed_preamble STREQUAL _generated_preamble)
    message(FATAL_ERROR "MindSplatter replay corpus preamble differs from generated form; "
      "regenerate with cmake --build --preset tests --target regenerate_mindsplatter_replay")
  endif()

  foreach(_source IN ITEMS _committed _generated)
    string(REGEX MATCH "\"seed=[^\"]+\"" ${_source}_provenance "${${_source}}")
    if(${_source}_provenance STREQUAL "")
      message(FATAL_ERROR "MindSplatter replay corpus omits search provenance")
    endif()
  endforeach()
  if(NOT _committed_provenance STREQUAL _generated_provenance)
    message(FATAL_ERROR "MindSplatter replay corpus provenance drift")
  endif()

  # Form, revision and search provenance are portable; payload low bits depend on libm.
  foreach(_symbol IN ITEMS
      "HEAVY_SEARCH_V1_STATE"
      "HEAVY_SEARCH_V1_FRAMEBUFFER"
      "Corpus HEAVY_SEARCH_V1")
    string(FIND "${_generated}" "${_symbol}" _symbol_offset)
    if(_symbol_offset EQUAL -1)
      message(FATAL_ERROR
        "MindSplatter replay generator omitted ${_symbol} from ${GENERATED}")
    endif()
  endforeach()
endfunction()

if(REFRESH_ONLY)
  execute_process(COMMAND "${GENERATOR}" "${GENERATED}" --refresh-framebuffer
    COMMAND_ERROR_IS_FATAL ANY)
  file(READ "${COMMITTED}" _committed)
  file(READ "${GENERATED}" _generated)
  set(_state_pattern "HEAVY_SEARCH_V1_STATE\\[\\] PROGMEM = \\{[^}]+\\}")
  string(REGEX MATCH "${_state_pattern}" _old_state "${_committed}")
  string(REGEX MATCH "${_state_pattern}" _new_state "${_generated}")
  if(_old_state STREQUAL "" OR NOT _old_state STREQUAL _new_state)
    message(FATAL_ERROR "MindSplatter framebuffer refresh changed the captured state")
  endif()
  check_generated_form()
  include("${CMAKE_CURRENT_LIST_DIR}/check_mindsplatter_replay_revision.cmake")
  check_mindsplatter_replay_revision("${_committed}" "${_generated}")
  message(STATUS "MindSplatter framebuffer refresh preserves captured state and revision")
  return()
endif()

if(NOT COMPILER_ID STREQUAL "Clang")
  if(REQUIRE_COMPILER_MATCH)
    message(FATAL_ERROR
      "MindSplatter replay form pin requires Clang, got "
      "${COMPILER_ID} ${COMPILER_VERSION}")
  endif()
  message(STATUS
    "MindSplatter replay form pin requires Clang; skipping "
    "${COMPILER_ID} ${COMPILER_VERSION}")
  cmake_language(EXIT ${SKIP_CODE})
endif()

if(NOT EXISTS "${GENERATOR}")
  message(FATAL_ERROR
    "MindSplatter replay form pin: generator not built: ${GENERATOR}")
endif()
if(NOT EXISTS "${COMMITTED}")
  message(FATAL_ERROR
    "MindSplatter replay form pin: committed corpus missing: ${COMMITTED}")
endif()

execute_process(
  COMMAND "${GENERATOR}" "${GENERATED}"
  RESULT_VARIABLE _rc
  ERROR_VARIABLE _err)
if(NOT _rc EQUAL 0)
  message(FATAL_ERROR "mindsplatter_replay_gen failed (${_rc}):\n${_err}")
endif()

file(READ "${COMMITTED}" _committed)
file(READ "${GENERATED}" _generated)
check_generated_form()

include("${CMAKE_CURRENT_LIST_DIR}/check_mindsplatter_replay_revision.cmake")
check_mindsplatter_replay_revision("${_committed}" "${_generated}")

message(STATUS
  "MindSplatter replay generator emitted ${MINDSPLATTER_REPLAY_REVISION} in canonical form")
