# Run scripts/generate_luts.py --check (the generator's monotonicity + sRGB
# round-trip self-validation). Skips with SKIP_CODE when no Python is available,
# or fails outright under REQUIRE_PYTHON (CI, which provisions the interpreter).
# -D args: PYTHON_EXE, GENERATOR, SKIP_CODE, REQUIRE_PYTHON.

# Script mode inherits no policies from the project, so every policy would
# otherwise default to OLD, and the cmake_language(EXIT) below is a 3.29
# feature. Matches the top-level CMakeLists.
cmake_minimum_required(VERSION 3.29)

if(NOT PYTHON_EXE OR NOT EXISTS "${PYTHON_EXE}")
  if(REQUIRE_PYTHON)
    message(FATAL_ERROR "color_luts self-check: no Python interpreter, and HS_REQUIRE_GENERATORS is ON")
  endif()
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
