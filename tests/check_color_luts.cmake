# Regenerate core/color/color_luts.h via scripts/generate_luts.py and compare the
# whitespace-normalized text against the committed file.
# The generator requires clang-format (CLANG_FORMAT or the one on PATH). Skips
# with SKIP_CODE when it is unavailable, or fails under REQUIRE_GENERATORS.
# -D args: PYTHON_EXE, GENERATOR, COMMITTED, GENERATED, SKIP_CODE, REQUIRE_GENERATORS.

# Match the top-level CMake policy version in script mode.
cmake_minimum_required(VERSION 3.29)

if("$ENV{CLANG_FORMAT}" STREQUAL "")
  find_program(_clang_format clang-format)
else()
  set(_clang_format "$ENV{CLANG_FORMAT}")
endif()
execute_process(
  COMMAND "${PYTHON_EXE}" "${CMAKE_CURRENT_LIST_DIR}/../tools/build_pins.py" clang-format
  OUTPUT_VARIABLE _required_format_version
  OUTPUT_STRIP_TRAILING_WHITESPACE
  COMMAND_ERROR_IS_FATAL ANY)
string(REGEX MATCH "^[0-9]+" _required_format_major "${_required_format_version}")

if(_clang_format)
  execute_process(COMMAND "${_clang_format}" --version
    OUTPUT_VARIABLE _format_version RESULT_VARIABLE _format_rc)
  if(NOT _format_rc EQUAL 0 OR NOT _format_version MATCHES "version ${_required_format_major}\\.")
    if(REQUIRE_GENERATORS)
      message(FATAL_ERROR "color_luts pin: clang-format ${_required_format_major} is required")
    endif()
    message(STATUS "color_luts pin: clang-format major mismatch; skipping")
    cmake_language(EXIT ${SKIP_CODE})
  endif()
endif()
if(NOT _clang_format)
  if(REQUIRE_GENERATORS)
    message(FATAL_ERROR "color_luts pin: no clang-format, and HS_REQUIRE_GENERATORS is ON")
  endif()
  message(STATUS "color_luts pin: no clang-format; skipping")
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

# Collapse whitespace to normalize CRLF and residual clang-format reflow.
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
    "Regenerate with: python scripts/generate_luts.py -o core/color/color_luts.h")
endif()

message(STATUS "color_luts pin: header text matches the generator")
