# Run the tangent-field constants provenance check when numpy is available.
cmake_minimum_required(VERSION 3.29)

execute_process(
  COMMAND "${PYTHON_EXE}" -c "import numpy"
  RESULT_VARIABLE _numpy_rc
  OUTPUT_QUIET ERROR_QUIET)
if(NOT _numpy_rc EQUAL 0)
  if(REQUIRE_GENERATORS)
    message(FATAL_ERROR "tangent_field_constants pin: no numpy, and HS_REQUIRE_GENERATORS is ON")
  endif()
  message(STATUS "tangent_field_constants pin: no numpy; skipping")
  cmake_language(EXIT ${SKIP_CODE})
endif()

execute_process(
  COMMAND "${PYTHON_EXE}" "${GENERATOR}" --check
  RESULT_VARIABLE _rc
  ERROR_VARIABLE _err)
if(NOT _rc EQUAL 0)
  message(FATAL_ERROR "gen_tangent_field_constants.py --check failed (${_rc}):\n${_err}")
endif()
message(STATUS "tangent_field_constants pin: generator --check passed")
