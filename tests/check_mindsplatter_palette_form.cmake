cmake_minimum_required(VERSION 3.29)

execute_process(
  COMMAND "${GENERATOR}" "${GENERATED}"
  RESULT_VARIABLE _rc
  ERROR_VARIABLE _err)
if(NOT _rc EQUAL 0)
  message(FATAL_ERROR "mindsplatter_palette_gen failed (${_rc}):\n${_err}")
endif()

file(READ "${COMMITTED}" _committed)
file(READ "${GENERATED}" _generated)
set(_pixel_pattern "Pixel\\([0-9]+, [0-9]+, [0-9]+\\)")
string(REGEX REPLACE "${_pixel_pattern}" "Pixel(VALUE)" _committed_form
                     "${_committed}")
string(REGEX REPLACE "${_pixel_pattern}" "Pixel(VALUE)" _generated_form
                     "${_generated}")
if(NOT _committed_form STREQUAL _generated_form)
  message(FATAL_ERROR
          "triadic_palette_luts.h form differs from mindsplatter_palette_gen")
endif()

message(STATUS "MindSplatter palette header form matches the generator")
