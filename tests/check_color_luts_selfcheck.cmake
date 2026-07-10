# Run scripts/generate_luts.py --check (the generator's monotonicity + sRGB
# round-trip self-validation). Skips with SKIP_CODE when no Python is available.
# -D args: PYTHON_EXE, GENERATOR, SKIP_CODE.

if(NOT PYTHON_EXE OR NOT EXISTS "${PYTHON_EXE}")
  message(STATUS "color_luts self-check: no Python interpreter; skipping")
  cmake_language(EXIT ${SKIP_CODE})
endif()

execute_process(
  COMMAND "${PYTHON_EXE}" "${GENERATOR}" --check
  RESULT_VARIABLE _rc
  ERROR_VARIABLE _err)
if(NOT _rc EQUAL 0)
  message(FATAL_ERROR "generate_luts.py --check failed (${_rc}):\n${_err}")
endif()

message(STATUS "color_luts self-check: generator --check passed")
