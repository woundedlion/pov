cmake_minimum_required(VERSION 3.29)
include("${CMAKE_CURRENT_LIST_DIR}/compile_fixture.cmake")
hs_check_compile_fixture(ring_column_pipeline_check.cpp ring_columns PASS)
hs_check_compile_fixture(ring_column_pipeline_check.cpp ring_columns PASS
  -DHS_RING_COLUMN_DIRECT)
hs_check_compile_fixture(ring_column_pipeline_check.cpp ring_columns
  "the ring column walk requires a screen-only pipeline" -DHS_RING_COLUMN_WORLD)
message(STATUS "ring column pipeline contracts passed")
