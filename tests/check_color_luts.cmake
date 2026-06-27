# Regenerate core/color_luts.h via scripts/generate_luts.py and compare numeric
# tokens against the committed file (whitespace/clang-format wrapping ignored).
# Skips with SKIP_CODE when no Python is available.
# -D args: PYTHON_EXE, GENERATOR, COMMITTED, GENERATED, SKIP_CODE.

if(NOT PYTHON_EXE OR NOT EXISTS "${PYTHON_EXE}")
  message(STATUS "color_luts pin: no Python interpreter; skipping")
  cmake_language(EXIT ${SKIP_CODE})
endif()

set(_generated "${GENERATED}")
execute_process(
  COMMAND "${PYTHON_EXE}" "${GENERATOR}"
  OUTPUT_FILE "${_generated}"
  RESULT_VARIABLE _rc
  ERROR_VARIABLE _err)
if(NOT _rc EQUAL 0)
  message(FATAL_ERROR "generate_luts.py failed (${_rc}):\n${_err}")
endif()

# Extract numeric tokens in order from each file, matching the CI grep -oE '[0-9]+'.
function(_numeric_tokens path out_var)
  file(READ "${path}" _text)
  string(REGEX MATCHALL "[0-9]+" _toks "${_text}")
  set(${out_var} "${_toks}" PARENT_SCOPE)
endfunction()

_numeric_tokens("${_generated}" _gen_toks)
_numeric_tokens("${COMMITTED}" _com_toks)

if(NOT _gen_toks STREQUAL _com_toks)
  message(FATAL_ERROR
    "core/color_luts.h is out of sync with scripts/generate_luts.py.\n"
    "Regenerate with: python scripts/generate_luts.py > core/color_luts.h "
    "&& clang-format -i core/color_luts.h")
endif()

message(STATUS "color_luts pin: numeric tokens match the generator")
