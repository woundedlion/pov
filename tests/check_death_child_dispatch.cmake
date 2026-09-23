# Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
# Licensed under the PolyForm Noncommercial License 1.0.0

if(NOT DEFINED TEST_BIN)
  message(FATAL_ERROR "TEST_BIN is required")
endif()

set(child_env
  --unset=HS_EFFECTS_FULL --unset=HS_REQUIRE_EFFECTS_FULL
  CI=1 HS_SKIPS_ARE_ERRORS=1 HS_SMOKE_FRAMES=120
  HS_BUFFER_FREE_WATCHDOG_US=30000000
  HS_DEATH_CHILD=harness HS_DEATH_CASE=__spawn_check__)

function(expect_dispatch expected diagnostic)
  cmake_parse_arguments(CASE "" "" "ENV;ARGS" ${ARGN})
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env ${child_env} ${CASE_ENV}
            "${TEST_BIN}" ${CASE_ARGS}
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error
    TIMEOUT 10)
  if(NOT "${result}" STREQUAL "${expected}")
    message(FATAL_ERROR "dispatch returned ${result}, expected ${expected}: ${output}${error}")
  endif()
  if(NOT "${diagnostic}" STREQUAL "" AND NOT "${error}" MATCHES "${diagnostic}")
    message(FATAL_ERROR "missing diagnostic '${diagnostic}': ${output}${error}")
  endif()
endfunction()

expect_dispatch(0 "")
expect_dispatch(0 "" ENV HS_EFFECTS_FULL=1 HS_REQUIRE_EFFECTS_FULL=1)
expect_dispatch(1 "HS_BUFFER_FREE_WATCHDOG_US" ENV HS_BUFFER_FREE_WATCHDOG_US=0)
expect_dispatch(2 "invalid death-child invocation" ARGS death)
expect_dispatch(2 "unknown module" ENV HS_DEATH_CHILD= ARGS __unknown_module__)
expect_dispatch(1 "HS_EFFECTS_FULL and HS_REQUIRE_EFFECTS_FULL" ENV HS_DEATH_CHILD=)
expect_dispatch(1 "HS_EFFECTS_FULL and HS_REQUIRE_EFFECTS_FULL" ENV HS_DEATH_CHILD= ARGS effects)
