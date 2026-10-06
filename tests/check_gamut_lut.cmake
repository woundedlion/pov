# Run tools/gen_gamut_lut.py --check: pins the OKLab matrices and gamut slack
# mirrored in the generator against core/color/color_space.h and its diamond angle
# against core/math/3dmath.h, then regenerates the table and diffs it against
# the committed core/color/gamut_lut.h in full.
# Skips with SKIP_CODE when numpy is unavailable, or fails under
# REQUIRE_GENERATORS.
# -D args: PYTHON_EXE, GENERATOR, SKIP_CODE, REQUIRE_GENERATORS.
# CACHE_FILE is optional: skip when hashed inputs match the last passing run.

# Match the top-level CMake policy version in script mode.
cmake_minimum_required(VERSION 3.29)

execute_process(
  COMMAND "${PYTHON_EXE}" -c "import numpy, sys; print(sys.version); print(numpy.__version__)"
  RESULT_VARIABLE _numpy_rc
  OUTPUT_VARIABLE _runtime ERROR_QUIET)
if(NOT _numpy_rc EQUAL 0)
  if(REQUIRE_GENERATORS)
    message(FATAL_ERROR "gamut_lut pin: no numpy, and HS_REQUIRE_GENERATORS is ON")
  endif()
  message(STATUS "gamut_lut pin: no numpy; skipping")
  cmake_language(EXIT ${SKIP_CODE})
endif()

# CI always regenerates; local runs reuse only a successful identical check.
if(CACHE_FILE AND NOT REQUIRE_GENERATORS)
  get_filename_component(_tools "${GENERATOR}" DIRECTORY)
  get_filename_component(_root "${_tools}" DIRECTORY)
  set(_inputs "${GENERATOR}" "${CMAKE_CURRENT_LIST_FILE}"
      "${_root}/core/color/color_space.h" "${_root}/core/math/3dmath.h"
      "${_root}/core/color/gamut_lut.h")
  set(_signature "${PYTHON_EXE};${_runtime};${CMAKE_HOST_SYSTEM}")
  foreach(_input IN LISTS _inputs)
    file(SHA256 "${_input}" _hash)
    string(APPEND _signature ";${_hash}")
  endforeach()
  string(SHA256 _signature "${_signature}")
  if(EXISTS "${CACHE_FILE}")
    file(READ "${CACHE_FILE}" _cached)
    if(_cached STREQUAL _signature)
      message(STATUS "gamut_lut pin: identical inputs previously regenerated successfully")
      return()
    endif()
  endif()
endif()

execute_process(
  COMMAND "${PYTHON_EXE}" "${GENERATOR}" --check
  RESULT_VARIABLE _rc
  ERROR_VARIABLE _err)
if(NOT _rc EQUAL 0)
  message(FATAL_ERROR "gen_gamut_lut.py --check failed (${_rc}):\n${_err}")
endif()

if(CACHE_FILE AND NOT REQUIRE_GENERATORS)
  file(WRITE "${CACHE_FILE}" "${_signature}")
endif()
message(STATUS "gamut_lut pin: generator --check passed")
