cmake_minimum_required(VERSION 3.29)

if(NOT TEST_ROOT MATCHES "mirrored_screenshot_prune_fixture$")
  message(FATAL_ERROR "Unsafe or missing TEST_ROOT: ${TEST_ROOT}")
endif()

file(REMOVE_RECURSE "${TEST_ROOT}")
set(_source "${TEST_ROOT}/source")
set(_daydream "${TEST_ROOT}/daydream")
file(MAKE_DIRECTORY "${_source}/nested")
file(MAKE_DIRECTORY "${_daydream}/docs/screenshots/nested")
file(WRITE "${_source}/keep.png" "source keep")
file(WRITE "${_source}/nested/keep.png" "source nested keep")
file(WRITE "${_daydream}/daydream.js" "marker")
file(WRITE "${_daydream}/index.html" "unrelated")
file(WRITE "${_daydream}/docs/screenshots/keep.png" "installed keep")
file(WRITE "${_daydream}/docs/screenshots/stale.png" "installed stale")
file(WRITE "${_daydream}/docs/screenshots/nested/keep.png" "installed nested keep")
file(WRITE "${_daydream}/docs/screenshots/nested/stale.png" "installed nested stale")
file(WRITE "${_daydream}/docs/screenshots/notes.txt" "unrelated")

set(HS_MIRROR_SOURCE "${_source}")
set(HS_DAYDREAM_DIR "${_daydream}")
include("${PRUNE_SCRIPT}")

foreach(_kept IN ITEMS
    "${_daydream}/daydream.js"
    "${_daydream}/index.html"
    "${_daydream}/docs/screenshots/keep.png"
    "${_daydream}/docs/screenshots/nested/keep.png"
    "${_daydream}/docs/screenshots/notes.txt")
  if(NOT EXISTS "${_kept}")
    message(FATAL_ERROR "Prune removed retained file: ${_kept}")
  endif()
endforeach()
foreach(_stale IN ITEMS
    "${_daydream}/docs/screenshots/stale.png"
    "${_daydream}/docs/screenshots/nested/stale.png")
  if(EXISTS "${_stale}")
    message(FATAL_ERROR "Prune retained stale screenshot: ${_stale}")
  endif()
endforeach()

set(_not_daydream "${TEST_ROOT}/not-daydream")
file(MAKE_DIRECTORY "${_not_daydream}/docs/screenshots")
file(WRITE "${_not_daydream}/docs/screenshots/stale.png" "must survive")
execute_process(
  COMMAND "${CMAKE_COMMAND}"
    "-DHS_MIRROR_SOURCE=${_source}"
    "-DHS_DAYDREAM_DIR=${_not_daydream}"
    -P "${PRUNE_SCRIPT}"
  RESULT_VARIABLE _invalid_result
  OUTPUT_QUIET
  ERROR_QUIET)
if(_invalid_result EQUAL 0)
  message(FATAL_ERROR "Prune accepted a destination without daydream.js")
endif()
if(NOT EXISTS "${_not_daydream}/docs/screenshots/stale.png")
  message(FATAL_ERROR "Prune modified an invalid destination")
endif()

file(REMOVE_RECURSE "${TEST_ROOT}")
