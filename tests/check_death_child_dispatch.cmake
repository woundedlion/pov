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
  cmake_parse_arguments(CASE "QUIET_STDOUT" "" "ENV;ARGS" ${ARGN})
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
  if(CASE_QUIET_STDOUT AND NOT "${output}" STREQUAL "")
    message(FATAL_ERROR "unexpected module output before rejection: ${output}")
  endif()
endfunction()

expect_dispatch(0 "")
expect_dispatch(0 "" ENV HS_EFFECTS_FULL=1 HS_REQUIRE_EFFECTS_FULL=1)
expect_dispatch(1 "HS_BUFFER_FREE_WATCHDOG_US" ENV HS_BUFFER_FREE_WATCHDOG_US=0)
expect_dispatch(1 "HS_SKIPS_ARE_ERRORS" ENV HS_SKIPS_ARE_ERRORS=)
foreach(value IN ITEMS "-1" "+30000000" "30000000us" "4294967296")
  expect_dispatch(1 "HS_BUFFER_FREE_WATCHDOG_US" ENV "HS_BUFFER_FREE_WATCHDOG_US=${value}")
endforeach()
expect_dispatch(1 "HS_SKIPS_ARE_ERRORS" ENV HS_SKIPS_ARE_ERRORS=true)
foreach(lever IN ITEMS HS_EFFECTS_FULL HS_REQUIRE_EFFECTS_FULL)
  expect_dispatch(1 "HS_EFFECTS_FULL and HS_REQUIRE_EFFECTS_FULL" ENV
    HS_DEATH_CHILD= HS_EFFECTS_FULL=1 HS_REQUIRE_EFFECTS_FULL=1 "${lever}=true" ARGS effects)
endforeach()
foreach(value IN ITEMS "8" "1200/2" "120e1" "2147483648")
  expect_dispatch(1 "HS_SMOKE_FRAMES" ENV "HS_SMOKE_FRAMES=${value}")
endforeach()
expect_dispatch(1 "QUICK tier" ENV HS_DEATH_CHILD= HS_EFFECTS_FULL=0 HS_REQUIRE_EFFECTS_FULL=1 ARGS effects)
expect_dispatch(2 "invalid death-child invocation" ARGS death)
expect_dispatch(2 "unknown module" ENV HS_DEATH_CHILD= ARGS __unknown_module__)
expect_dispatch(2 "unknown module" QUIET_STDOUT ENV HS_DEATH_CHILD=
  ARGS wasm_predicates __unknown_module__)
expect_dispatch(1 "HS_EFFECTS_FULL and HS_REQUIRE_EFFECTS_FULL" ENV HS_DEATH_CHILD=)
expect_dispatch(1 "HS_EFFECTS_FULL and HS_REQUIRE_EFFECTS_FULL" ENV HS_DEATH_CHILD= ARGS effects)

foreach(lever IN ITEMS HS_ISLAMICSTARS_ONLY HS_MESH_CAPTURE_ONLY)
  expect_dispatch(0 "" ENV "${lever}=1")
  expect_dispatch(1 "${lever} narrows the run because it is set"
    ENV HS_DEATH_CHILD= HS_EFFECTS_FULL=1 HS_REQUIRE_EFFECTS_FULL=1 "${lever}=1"
    ARGS effects)
endforeach()

foreach(lever IN ITEMS HS_RING_TEMPORAL_SCENARIO HS_SURFACE_LOBE_DIAGNOSTICS
                       HS_SURFACE_DARK_DIAGNOSTICS)
  foreach(value IN ITEMS "" "1")
    expect_dispatch(1 "${lever} narrows the run because it is set"
      ENV HS_DEATH_CHILD= "${lever}=${value}" ARGS shapeshifter_oracle)
    expect_dispatch(1 "${lever} narrows the run because it is set"
      ENV HS_DEATH_CHILD= HS_EFFECTS_FULL=1 HS_REQUIRE_EFFECTS_FULL=1 "${lever}=${value}" ARGS effects)
    expect_dispatch(0 "" ENV CI= "${lever}=${value}")
  endforeach()
endforeach()

get_filename_component(bin_dir "${TEST_BIN}" DIRECTORY)
get_filename_component(runner "${TEST_BIN}" NAME)
set(shard_dir "${CMAKE_CURRENT_BINARY_DIR}/dispatch-shard")
file(MAKE_DIRECTORY "${shard_dir}")
set(executable "${TEST_BIN}")
set(path_env)
if(UNIX)
  set(executable "${runner}")
  set(path_env "PATH=${bin_dir}:$ENV{PATH}")
endif()
execute_process(
  COMMAND "${CMAKE_COMMAND}" -E env ${child_env} ${path_env} HS_DEATH_CHILD= HS_DEATH_CASE=
    HS_DEATH_SHARD=1023/1024 "HS_DEATH_SHARD_DIR=${shard_dir}" "${executable}" death
  WORKING_DIRECTORY "${shard_dir}"
  RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error TIMEOUT 60)
if(NOT result EQUAL 0)
  message(FATAL_ERROR "death shard failed: ${output}${error}")
endif()
file(READ "${shard_dir}/death_shard_1023.txt" record)
if(NOT record MATCHES "^hs-death-shard 1023 1024 ([0-9]+)")
  message(FATAL_ERROR "missing shard case count: ${record}")
endif()
file(WRITE "${shard_dir}/death_shard_0.txt" "hs-death-shard 0 1 ${CMAKE_MATCH_1}\n0 -1\n")
execute_process(
  COMMAND "${CMAKE_COMMAND}" -E env ${child_env} HS_DEATH_CHILD= HS_DEATH_CASE=
    HS_DEATH_SHARD=census/1 "HS_DEATH_SHARD_DIR=${shard_dir}" "${TEST_BIN}" death
  RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error TIMEOUT 15)
if(result EQUAL 0 OR NOT output MATCHES "has a bad row: 0 -1")
  message(FATAL_ERROR "census missed the invalid source line: ${output}${error}")
endif()
