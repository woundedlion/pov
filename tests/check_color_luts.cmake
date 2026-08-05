# Regenerate core/color/color_luts.h via scripts/generate_luts.py and compare the
# FULL text against the committed file, so array names, element types, the
# flash-section marker, the include and value signs gate alongside the numbers.
# Counterpart of the lut-provenance job in .github/workflows/ci.yml.
# Skips with SKIP_CODE when no Python is available, or fails outright under
# REQUIRE_PYTHON (CI, which provisions the interpreter).
# -D args: PYTHON_EXE, GENERATOR, COMMITTED, GENERATED, SKIP_CODE, REQUIRE_PYTHON.

# Script mode inherits no policies from the project, so every policy would
# otherwise default to OLD, and the cmake_language(EXIT) below is a 3.29
# feature. Matches the top-level CMakeLists.
cmake_minimum_required(VERSION 3.29)

if(NOT PYTHON_EXE OR NOT EXISTS "${PYTHON_EXE}")
  if(REQUIRE_PYTHON)
    message(FATAL_ERROR "color_luts pin: no Python interpreter, and HS_REQUIRE_GENERATORS is ON")
  endif()
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

# Whole file with whitespace runs collapsed. CI normalizes both sides through
# clang-format instead; collapsing whitespace absorbs the same reflow (and a CRLF
# checkout) without requiring clang-format, and changes no other character.
function(_normalized_text path out_var)
  file(READ "${path}" _text)
  string(REGEX REPLACE "[ \t\r\n]+" " " _text "${_text}")
  string(STRIP "${_text}" _text)
  set(${out_var} "${_text}" PARENT_SCOPE)
endfunction()

_normalized_text("${_generated}" _gen_text)
_normalized_text("${COMMITTED}" _com_text)

if(NOT _gen_text STREQUAL _com_text)
  message(FATAL_ERROR
    "core/color/color_luts.h is out of sync with scripts/generate_luts.py.\n"
    "Diff it against the regenerated header: ${_generated}\n"
    "Regenerate with: python scripts/generate_luts.py > core/color/color_luts.h "
    "&& clang-format -i core/color/color_luts.h")
endif()

message(STATUS "color_luts pin: header text matches the generator")
