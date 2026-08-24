# Re-emit the segment->canvas mapping from hardware/pov_segment_map.h (via the
# pov_segment_map_gen exporter, which compiles against that header) and compare
# the FULL text against the committed hardware/pov_segment_map.json. That golden
# is the firmware reference the simulator's tests/segment_crosscheck.test.js
# reads, so any convention change in the header must land in the JSON — and
# therefore in the web tiling's cross-check — instead of leaving a stale JS
# port green.
# -D args: EXPORTER, COMMITTED, GENERATED.

cmake_minimum_required(VERSION 3.29)

if(NOT EXISTS "${EXPORTER}")
  message(FATAL_ERROR "pov_segment_map pin: exporter not built: ${EXPORTER}")
endif()

execute_process(
  COMMAND "${EXPORTER}"
  OUTPUT_FILE "${GENERATED}"
  RESULT_VARIABLE _rc
  ERROR_VARIABLE _err)
if(NOT _rc EQUAL 0)
  message(FATAL_ERROR "pov_segment_map_gen failed (${_rc}):\n${_err}")
endif()

# Exact bytes pin the readable layout and canonical LF line endings.
execute_process(
  COMMAND "${CMAKE_COMMAND}" -E compare_files "${GENERATED}" "${COMMITTED}"
  RESULT_VARIABLE _compare_rc)

if(NOT _compare_rc EQUAL 0)
  message(FATAL_ERROR
    "hardware/pov_segment_map.json is out of sync with "
    "hardware/pov_segment_map.h.\n"
    "Diff it against the regenerated golden: ${GENERATED}\n"
    "Regenerate with: cmake --build --preset tests --target pov_segment_map_gen "
    "&& build/tests/tests/pov_segment_map_gen > hardware/pov_segment_map.json\n"
    "Then re-run the WASM install so the daydream checkout's copy travels with "
    "the module, and commit both.")
endif()

message(STATUS "pov_segment_map pin: golden JSON matches the header")
