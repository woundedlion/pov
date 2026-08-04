# Run tools/gen_gamut_lut.py --check: pins the OKLab matrices and gamut slack
# mirrored in the generator against core/color/color.h, then regenerates the
# table and diffs it against the committed core/color/gamut_lut.h in full.
# Skips with SKIP_CODE when Python or numpy is unavailable, or fails outright
# under REQUIRE_PYTHON (CI, which provisions both).
# -D args: PYTHON_EXE, GENERATOR, SKIP_CODE, REQUIRE_PYTHON.

# Script mode inherits no policies from the project, so every policy would
# otherwise default to OLD, and the cmake_language(EXIT) below is a 3.29
# feature. Matches the top-level CMakeLists.
cmake_minimum_required(VERSION 3.29)

if(NOT PYTHON_EXE OR NOT EXISTS "${PYTHON_EXE}")
  if(REQUIRE_PYTHON)
    message(FATAL_ERROR "gamut_lut pin: no Python interpreter, and HS_REQUIRE_GENERATORS is ON")
  endif()
  message(STATUS "gamut_lut pin: no Python interpreter; skipping")
  cmake_language(EXIT ${SKIP_CODE})
endif()

execute_process(
  COMMAND "${PYTHON_EXE}" -c "import numpy"
  RESULT_VARIABLE _numpy_rc
  OUTPUT_QUIET ERROR_QUIET)
if(NOT _numpy_rc EQUAL 0)
  if(REQUIRE_PYTHON)
    message(FATAL_ERROR "gamut_lut pin: no numpy, and HS_REQUIRE_GENERATORS is ON")
  endif()
  message(STATUS "gamut_lut pin: no numpy; skipping")
  cmake_language(EXIT ${SKIP_CODE})
endif()

execute_process(
  COMMAND "${PYTHON_EXE}" "${GENERATOR}" --check
  RESULT_VARIABLE _rc
  ERROR_VARIABLE _err)
if(NOT _rc EQUAL 0)
  message(FATAL_ERROR "gen_gamut_lut.py --check failed (${_rc}):\n${_err}")
endif()

message(STATUS "gamut_lut pin: generator --check passed")
