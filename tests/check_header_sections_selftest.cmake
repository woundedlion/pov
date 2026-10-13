cmake_minimum_required(VERSION 3.29)
include("${CMAKE_CURRENT_LIST_DIR}/gate_selftest.cmake")
hs_gate_workspace("${CALL_GATE}")
file(MAKE_DIRECTORY "${WORK}/tests/alpha")
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
  hs_run_gate("${gate}" _result _output)
  if(should_pass AND NOT _result EQUAL 0)
    message(FATAL_ERROR "valid section rejected: ${_output}")
  elseif(NOT should_pass AND _result EQUAL 0)
    message(FATAL_ERROR "invalid section accepted: ${diagnostic}")
  endif()
  if(NOT "${_output}" MATCHES "${diagnostic}")
    message(FATAL_ERROR "unexpected section verdict: ${_output}")
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

hs_gate_workspace("${CALL_GATE}")
set(_local_helpers [=[inline void capture_frame() {
}
inline void test_guard() {
  capture_frame();
}
inline int run_alpha_tests() {
  test_guard();
  return 0;
}
]=])
set(_asserting_helpers [=[inline void capture_frame() {
  HS_EXPECT_TRUE(true);
}
]=])
file(WRITE "${WORK}/tests/test_alpha.h" "\n${_local_helpers}")
file(WRITE "${WORK}/tests/pixel_helpers.h" "\n${_asserting_helpers}")
check_gate("${CALL_GATE}" FALSE "test cases reach no assertion: test_alpha.h:test_guard")
string(REPLACE "inline void capture_frame() {\n}" "inline void capture_frame() {\n  HS_EXPECT_TRUE(true);\n}" _asserting_local "${_local_helpers}")
file(WRITE "${WORK}/tests/test_alpha.h" "\n${_asserting_local}")
check_gate("${CALL_GATE}" TRUE "all 1 case definitions")

string(REPLACE "inline void capture_frame() {\n}\n" "" _cross_header "${_local_helpers}")
file(WRITE "${WORK}/tests/test_alpha.h" "\n${_cross_header}")
check_gate("${CALL_GATE}" TRUE "all 1 case definitions")
file(WRITE "${WORK}/tests/replay_helpers.h" "\ninline void capture_frame() {\n}\n")
check_gate("${CALL_GATE}" FALSE "test cases reach no assertion: test_alpha.h:test_guard")
file(REMOVE "${WORK}/tests/replay_helpers.h")

set(_overload [=[inline void capture_frame(int value) {
  HS_EXPECT_TRUE(value);
}
]=])
file(WRITE "${WORK}/tests/test_alpha.h" "\n${_local_helpers}${_overload}")
check_gate("${CALL_GATE}" TRUE "all 1 case definitions")

hs_gate_workspace("${CALL_GATE}")
set(_helper_calls [=[inline void test_guard() {
  HS_EXPECT_TRUE(true);
}
inline void unused_helper() {
  test_guard();
}
inline int run_alpha_tests() {
  return 0;
}
]=])
file(WRITE "${WORK}/tests/test_alpha.h" "\n${_helper_calls}")
check_gate("${CALL_GATE}" FALSE "test case defined but never reached: test_alpha.h:test_guard")
string(REPLACE "  return 0;" "  unused_helper();\n  return 0;" _called_helper "${_helper_calls}")
file(WRITE "${WORK}/tests/test_alpha.h" "\n${_called_helper}")
check_gate("${CALL_GATE}" TRUE "all 1 case definitions")

string(REPLACE "inline void unused_helper()" "inline StrapSweepStats unused_helper()" _typed_helper "${_helper_calls}")
file(WRITE "${WORK}/tests/test_alpha.h" "\n${_typed_helper}")
check_gate("${CALL_GATE}" FALSE "test case defined but never reached: test_alpha.h:test_guard")
string(REPLACE "  return 0;" "  unused_helper();\n  return 0;" _called_typed_helper "${_typed_helper}")
file(WRITE "${WORK}/tests/test_alpha.h" "\n${_called_typed_helper}")
check_gate("${CALL_GATE}" TRUE "all 1 case definitions")
file(WRITE "${WORK}/tests/test_alpha.h" "\n${_helper_calls}\ninline void declared_helper();\n#define CASE_DISPATCH unused_helper\n")
check_gate("${CALL_GATE}" TRUE "all 1 case definitions")
message(STATUS "test section controls passed")
