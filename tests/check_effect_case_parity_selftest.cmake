cmake_minimum_required(VERSION 3.29)
include("${CMAKE_CURRENT_LIST_DIR}/gate_selftest.cmake")
if(NOT DEFINED WORK OR NOT DEFINED GATE)
  message(FATAL_ERROR "WORK and GATE are required")
endif()
file(MAKE_DIRECTORY "${WORK}")
set(umbrella [=[inline int run_effects_tests() {
  if (dynamo_only) {
    run_effect_case_debug();
    return fixture.result();
  }
  run_effect_case_alpha();
  if (effects_full_suite()) {
    run_effect_case_beta();
  }
  return fixture.result();
}
]=])
set(focused [=[void run_effect_case_alpha() {}
void run_effect_case_beta() {}
void run_effect_case_debug() {}
int run_alpha_tests() {
  run_effect_case_alpha();
  if (effects_full_suite()) {
    run_effect_case_beta();
  }
  return fixture.result();
}
]=])
file(WRITE "${WORK}/test_effects.h" "${umbrella}")

function(check_fixture source expected diagnostic)
  file(WRITE "${WORK}/run_effect_cases_alpha.cpp" "${source}")
  execute_process(COMMAND "${CMAKE_COMMAND}" "-DTESTS_DIR=${WORK}" -P "${GATE}"
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
  hs_expect_gate_verdict("${result}" "${output}${error}" "${expected}" "${diagnostic}")
endfunction()

check_fixture("${focused}" TRUE "")
string(REPLACE "  run_effect_case_alpha();\n" "" missing "${focused}")
check_fixture("${missing}" FALSE "membership/tier mismatch")
string(REPLACE "  if (effects_full_suite()) {\n    run_effect_case_beta();\n  }\n"
  "  run_effect_case_beta();\n" demoted "${focused}")
check_fixture("${demoted}" FALSE "membership/tier mismatch")
