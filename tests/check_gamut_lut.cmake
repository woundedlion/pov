# Run tools/gen_gamut_lut.py --check: pins the OKLab matrices and gamut slack
# mirrored in the generator against core/color/color.h, then regenerates the
# table and diffs its numeric tokens against the committed core/color/gamut_lut.h.
# Skips with SKIP_CODE when Python or numpy is unavailable.
# -D args: PYTHON_EXE, GENERATOR, SKIP_CODE.

if(NOT PYTHON_EXE OR NOT EXISTS "${PYTHON_EXE}")
  message(STATUS "gamut_lut pin: no Python interpreter; skipping")
  cmake_language(EXIT ${SKIP_CODE})
endif()

execute_process(
  COMMAND "${PYTHON_EXE}" -c "import numpy"
  RESULT_VARIABLE _numpy_rc
  OUTPUT_QUIET ERROR_QUIET)
if(NOT _numpy_rc EQUAL 0)
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
