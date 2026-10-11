# Scratch-tree harness for the run_tests.cpp roster gate self-tests.

# Reset WORK to an empty tests/ and tools/ tree carrying gate's section
# expander and an empty off-roster list.
function(hs_gate_workspace gate)
  get_filename_component(_gate_dir "${gate}" DIRECTORY)
  file(REMOVE_RECURSE "${WORK}")
  file(MAKE_DIRECTORY "${WORK}/tests" "${WORK}/tools")
  file(COPY "${_gate_dir}/header_sections.cmake" DESTINATION "${WORK}/tests")
  file(WRITE "${WORK}/tests/off_roster_headers.cmake" "set(HS_OFF_ROSTER_HEADER_NAMES)\n")
endfunction()

# Run gate over WORK/run_tests.cpp; store its exit code in result_var and its
# combined stdout and stderr in output_var.
function(hs_run_gate gate result_var output_var)
  execute_process(COMMAND "${CMAKE_COMMAND}"
    "-DSRC=${WORK}/run_tests.cpp" "-DTESTS_DIR=${WORK}/tests"
    "-DTOOLS_DIR=${WORK}/tools" -P "${gate}"
    RESULT_VARIABLE _result OUTPUT_VARIABLE _output ERROR_VARIABLE _error)
  set(${result_var} "${_result}" PARENT_SCOPE)
  set(${output_var} "${_output}${_error}" PARENT_SCOPE)
endfunction()
