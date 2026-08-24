cmake_minimum_required(VERSION 3.29)

string(REGEX MATCH "^[0-9]+" _compiler_major "${COMPILER_VERSION}")
if(NOT COMPILER_ID STREQUAL "Clang" OR NOT _compiler_major EQUAL 22)
  if(REQUIRE_COMPILER_MATCH)
    message(FATAL_ERROR
      "MindSplatter replay form pin requires Clang 22, got "
      "${COMPILER_ID} ${COMPILER_VERSION}")
  endif()
  message(STATUS
    "MindSplatter replay form pin requires Clang 22; skipping "
    "${COMPILER_ID} ${COMPILER_VERSION}")
  cmake_language(EXIT ${SKIP_CODE})
endif()

if(NOT EXISTS "${GENERATOR}")
  message(FATAL_ERROR
    "MindSplatter replay form pin: generator not built: ${GENERATOR}")
endif()

execute_process(
  COMMAND "${GENERATOR}" "${GENERATED}"
  RESULT_VARIABLE _rc
  ERROR_VARIABLE _err)
if(NOT _rc EQUAL 0)
  message(FATAL_ERROR "mindsplatter_replay_gen failed (${_rc}):\n${_err}")
endif()

execute_process(
  COMMAND "${CMAKE_COMMAND}" -E compare_files "${COMMITTED}" "${GENERATED}"
  RESULT_VARIABLE _rc)
if(NOT _rc EQUAL 0)
  message(FATAL_ERROR
    "mindsplatter_replay_corpus.h differs from mindsplatter_replay_gen; "
    "diff it against ${GENERATED}")
endif()

message(STATUS "MindSplatter replay corpus matches the generator")
