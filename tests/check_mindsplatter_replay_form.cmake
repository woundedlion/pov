cmake_minimum_required(VERSION 3.29)

execute_process(
  COMMAND "${PYTHON}" "${CMAKE_CURRENT_LIST_DIR}/../tools/build_pins.py" clang
  OUTPUT_VARIABLE _required_clang_major
  OUTPUT_STRIP_TRAILING_WHITESPACE
  COMMAND_ERROR_IS_FATAL ANY)
if(NOT _required_clang_major MATCHES "^[0-9]+$")
  message(FATAL_ERROR "Invalid Clang major from tools/build_pins.py")
endif()

string(REGEX MATCH "^[0-9]+" _compiler_major "${COMPILER_VERSION}")
if(NOT COMPILER_ID STREQUAL "Clang" OR NOT _compiler_major EQUAL _required_clang_major)
  if(REQUIRE_COMPILER_MATCH)
    message(FATAL_ERROR
      "MindSplatter replay form pin requires Clang ${_required_clang_major}, got "
      "${COMPILER_ID} ${COMPILER_VERSION}")
  endif()
  message(STATUS
    "MindSplatter replay form pin requires Clang ${_required_clang_major}; skipping "
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
file(SIZE "${GENERATED}" _generated_size)
if(_generated_size LESS 100000)
  message(FATAL_ERROR
    "MindSplatter replay generator emitted only ${_generated_size} bytes")
endif()

# The replay's RandomWalk orientation uses libm under -ffast-math. CPU-specific
# low-bit results amplify over the particle simulation, while instrumentation
# changes can select a different worst-workload frame. Only form and revision
# are portable. unit_mindsplatter_replay validates the committed payload against
# the live renderer.
foreach(_symbol IN ITEMS
    "HEAVY_SEARCH_V1_STATE"
    "HEAVY_SEARCH_V1_FRAMEBUFFER"
    "Corpus HEAVY_SEARCH_V1"
    "CORPUS_MANIFEST[] = {&HEAVY_SEARCH_V1}")
  string(FIND "${_generated}" "${_symbol}" _symbol_offset)
  if(_symbol_offset EQUAL -1)
    message(FATAL_ERROR
      "MindSplatter replay generator omitted ${_symbol} from ${GENERATED}")
  endif()
endforeach()

include("${CMAKE_CURRENT_LIST_DIR}/check_mindsplatter_replay_revision.cmake")
check_mindsplatter_replay_revision("${_committed}" "${_generated}")

message(STATUS
  "MindSplatter replay generator emitted ${MINDSPLATTER_REPLAY_REVISION} in canonical form")
