# Run the size-trail host self-tests (tools/teensy_gate_tests/test_size_trail.py).
# Skips with SKIP_CODE when no Python is available, or fails outright under
# REQUIRE_PYTHON (CI, which provisions the interpreter).
# -D args: PYTHON_EXE, HS_ROOT, SKIP_CODE, REQUIRE_PYTHON.

if(NOT PYTHON_EXE OR NOT EXISTS "${PYTHON_EXE}")
  if(REQUIRE_PYTHON)
    message(FATAL_ERROR "teensy size trail: no Python interpreter, and HS_REQUIRE_GENERATORS is ON")
  endif()
  message(STATUS "teensy size trail: no Python interpreter; skipping")
  cmake_language(EXIT ${SKIP_CODE})
endif()

execute_process(
  COMMAND "${PYTHON_EXE}" -m unittest discover
          -s "${HS_ROOT}/tools/teensy_gate_tests" -p test_size_trail.py -v
  WORKING_DIRECTORY "${HS_ROOT}"
  RESULT_VARIABLE _rc
  ERROR_VARIABLE _err)
if(NOT _rc EQUAL 0)
  message(FATAL_ERROR "test_size_trail.py failed (${_rc}):\n${_err}")
endif()

message(STATUS "teensy size trail: self-tests passed")
