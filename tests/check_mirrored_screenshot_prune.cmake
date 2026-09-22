cmake_minimum_required(VERSION 3.29)

if(NOT TEST_ROOT MATCHES "mirrored_screenshot_prune_fixture$")
  message(FATAL_ERROR "Unsafe or missing TEST_ROOT: ${TEST_ROOT}")
endif()

file(REMOVE_RECURSE "${TEST_ROOT}")
set(_source "${TEST_ROOT}/engine/docs/screenshots")
file(MAKE_DIRECTORY "${TEST_ROOT}/engine")
file(WRITE "${TEST_ROOT}/engine/CMakeLists.txt" "marker")
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

set(_source "${TEST_ROOT}/engine/patterns")
file(MAKE_DIRECTORY "${_source}")
set(_patterns "${_daydream}/shader/patterns")
file(MAKE_DIRECTORY "${_patterns}/v1")
file(WRITE "${_patterns}/v1/example.shader.json" "daydream migration fixture")
file(WRITE "${_source}/keep.shader.json" "source")
file(WRITE "${_patterns}/keep.shader.json" "keep")
file(WRITE "${_patterns}/stale.shader.json" "stale")
file(WRITE "${_patterns}/shaderball_migration.json" "stale migration")
file(WRITE "${_patterns}/notes.json" "unrelated")
set(HS_MIRROR_SOURCE "${_source}")
set(HS_DAYDREAM_DIR "${_daydream}")
include("${CMAKE_CURRENT_LIST_DIR}/../cmake/prune_mirrored_patterns.cmake")
if(EXISTS "${_patterns}/stale.shader.json" OR
   EXISTS "${_patterns}/shaderball_migration.json")
  message(FATAL_ERROR "Prune retained stale mirrored patterns")
endif()
if(NOT EXISTS "${_patterns}/keep.shader.json" OR
   NOT EXISTS "${_patterns}/notes.json" OR
   NOT EXISTS "${_patterns}/v1/example.shader.json")
  message(FATAL_ERROR "Prune removed retained pattern files")
endif()

foreach(_kind IN ITEMS screenshots patterns)
  if(_kind STREQUAL "screenshots")
    set(_empty "${TEST_ROOT}/engine/docs/empty")
    set(_retained "${_daydream}/docs/screenshots/keep.png")
  else()
    set(_empty "${TEST_ROOT}/engine/empty")
    set(_retained "${_patterns}/keep.shader.json")
  endif()
  file(MAKE_DIRECTORY "${_empty}")
  execute_process(COMMAND "${CMAKE_COMMAND}"
    "-DHS_MIRROR_SOURCE=${_empty}" "-DHS_DAYDREAM_DIR=${_daydream}"
    -P "${CMAKE_CURRENT_LIST_DIR}/../cmake/prune_mirrored_${_kind}.cmake"
    RESULT_VARIABLE _empty_result OUTPUT_QUIET ERROR_QUIET)
  if(_empty_result EQUAL 0 OR NOT EXISTS "${_retained}")
    message(FATAL_ERROR "Empty ${_kind} source was not safely rejected")
  endif()
endforeach()

file(REMOVE_RECURSE "${TEST_ROOT}")
