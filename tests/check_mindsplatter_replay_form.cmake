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

foreach(_symbol IN ITEMS
    "HEAVY_SEARCH_V1_STATE"
    "HEAVY_SEARCH_V1_FRAMEBUFFER"
    "Corpus HEAVY_SEARCH_V1"
    "CORPUS_MANIFEST")
  string(FIND "${_generated}" "${_symbol}" _symbol_offset)
  if(_symbol_offset EQUAL -1)
    message(FATAL_ERROR
      "MindSplatter replay generator omitted ${_symbol} from ${GENERATED}")
  endif()
endforeach()

foreach(_identity_pattern IN ITEMS
    "heavy_search_v1_p[0-9]+_f[0-9]+"
    "msp-heavy-search-v[0-9]+")
  string(REGEX MATCH "${_identity_pattern}" _committed_identity "${_committed}")
  string(REGEX MATCH "${_identity_pattern}" _generated_identity "${_generated}")
  if(NOT _committed_identity STREQUAL _generated_identity)
    message(FATAL_ERROR
      "MindSplatter replay identity drift for ${_identity_pattern}: committed "
      "${_committed_identity}, generated ${_generated_identity}")
  endif()
endforeach()

message(STATUS
  "MindSplatter replay generator emitted ${_generated_identity} in canonical form")
