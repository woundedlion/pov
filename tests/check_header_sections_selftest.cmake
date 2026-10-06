cmake_minimum_required(VERSION 3.29)
get_filename_component(_gate_dir "${CALL_GATE}" DIRECTORY)
file(REMOVE_RECURSE "${WORK}")
file(MAKE_DIRECTORY "${WORK}/tests/alpha" "${WORK}/tools")
file(COPY "${_gate_dir}/header_sections.cmake" DESTINATION "${WORK}/tests")
file(WRITE "${WORK}/tests/off_roster_headers.cmake" "set(HS_OFF_ROSTER_HEADER_NAMES)\n")
set(_header [=[#include "tests/alpha/guards.h"
inline int run_alpha_tests() {
  test_guard();
  return 0;
}
]=])
set(_part [=[inline void test_guard() {
  HS_EXPECT_TRUE(true);
}
]=])
file(WRITE "${WORK}/tests/test_alpha.h" "${_header}")
file(WRITE "${WORK}/tests/alpha/guards.h" "${_part}")
file(WRITE "${WORK}/run_tests.cpp" [=[#include "tests/test_alpha.h"
#define HS_TEST_MODULE_LIST(X)
X("alpha", hs_test::alpha::run_alpha_tests, false)
#define HS_TEST_MODULE_ENTRY
]=])

function(check_gate gate should_pass diagnostic)
  execute_process(COMMAND "${CMAKE_COMMAND}"
    "-DSRC=${WORK}/run_tests.cpp" "-DTESTS_DIR=${WORK}/tests"
    "-DTOOLS_DIR=${WORK}/tools" -P "${gate}"
    RESULT_VARIABLE _result OUTPUT_VARIABLE _output ERROR_VARIABLE _error)
  if(should_pass AND NOT _result EQUAL 0)
    message(FATAL_ERROR "valid section rejected: ${_output}${_error}")
  elseif(NOT should_pass AND _result EQUAL 0)
    message(FATAL_ERROR "invalid section accepted: ${diagnostic}")
  endif()
  if(NOT "${_output}${_error}" MATCHES "${diagnostic}")
    message(FATAL_ERROR "unexpected section verdict: ${_output}${_error}")
  endif()
endfunction()

check_gate("${INCLUDE_GATE}" TRUE "includes match the roster")
check_gate("${CALL_GATE}" TRUE "all 1 case definitions")
string(REPLACE "  test_guard();\n" "" _missing_call "${_header}")
file(WRITE "${WORK}/tests/test_alpha.h" "${_missing_call}")
check_gate("${CALL_GATE}" FALSE "test case defined but never reached")
file(WRITE "${WORK}/tests/test_alpha.h" "${_header}")
string(REPLACE "  HS_EXPECT_TRUE(true);\n" "" _missing_assertion "${_part}")
file(WRITE "${WORK}/tests/alpha/guards.h" "${_missing_assertion}")
check_gate("${CALL_GATE}" FALSE "test cases reach no assertion")
file(WRITE "${WORK}/tests/alpha/guards.h" "${_part}")
string(REPLACE "#include \"tests/alpha/guards.h\"\n" "" _missing_include "${_header}")
file(WRITE "${WORK}/tests/test_alpha.h" "${_missing_include}")
check_gate("${INCLUDE_GATE}" FALSE "owner must include its section once")
file(WRITE "${WORK}/tests/test_alpha.h" "${_header}${_header}")
check_gate("${INCLUDE_GATE}" FALSE "owner must include its section once")
file(WRITE "${WORK}/tests/test_alpha.h" "${_header}")
file(MAKE_DIRECTORY "${WORK}/core/alpha")
file(WRITE "${WORK}/core/alpha/guards.h" "${_part}")
check_gate("${INCLUDE_GATE}" FALSE "section header shadows an engine include")
check_gate("${CALL_GATE}" FALSE "section header shadows an engine include")
file(REMOVE "${WORK}/core/alpha/guards.h")
check_gate("${CALL_GATE}" TRUE "all 1 case definitions")
message(STATUS "test section controls passed")
