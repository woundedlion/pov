# Emit the pullback operation stream from the manifest, replay it through the
# native capture producer TWICE at each supported resolution, and require the
# two captures to be byte-identical: the cross-checker compares a base and a
# candidate capture, so a producer that is not replay-deterministic reports
# every run as a difference and the comparison means nothing. The captures are
# then handed back to tools/pullback_capture.py --backend-audit, which scores
# them against the manifest oracles, so a run that is stably deterministic but
# stably wrong fails here too.
# -D args: PYTHON, CAPTURE_SCRIPT, MANIFEST_DIR, PRODUCER, OUTPUT_DIR.

# Script mode inherits no policies from the project, so every policy would
# otherwise default to OLD. Matches the top-level CMakeLists.
cmake_minimum_required(VERSION 3.29)

set(OPERATIONS "${OUTPUT_DIR}/pullback_capture_operations.bin")
execute_process(
  COMMAND "${PYTHON}" "${CAPTURE_SCRIPT}"
          --manifest-dir "${MANIFEST_DIR}" --operations-output "${OPERATIONS}"
  RESULT_VARIABLE OPERATIONS_RESULT)
if(NOT OPERATIONS_RESULT EQUAL 0)
  message(FATAL_ERROR "pullback operation generation failed: ${OPERATIONS_RESULT}")
endif()
foreach(RESOLUTION IN ITEMS 96x20 288x144)
  set(FIRST "${OUTPUT_DIR}/pullback_capture_${RESOLUTION}_first.bin")
  set(SECOND "${OUTPUT_DIR}/pullback_capture_${RESOLUTION}_second.bin")
  execute_process(
    COMMAND "${PRODUCER}" --resolution "${RESOLUTION}" --operations "${OPERATIONS}" "${FIRST}"
    RESULT_VARIABLE FIRST_RESULT)
  if(NOT FIRST_RESULT EQUAL 0)
    message(FATAL_ERROR "first ${RESOLUTION} pullback capture failed: ${FIRST_RESULT}")
  endif()
  execute_process(
    COMMAND "${PRODUCER}" --resolution "${RESOLUTION}" --operations "${OPERATIONS}" "${SECOND}"
    RESULT_VARIABLE SECOND_RESULT)
  if(NOT SECOND_RESULT EQUAL 0)
    message(FATAL_ERROR "second ${RESOLUTION} pullback capture failed: ${SECOND_RESULT}")
  endif()
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E compare_files "${FIRST}" "${SECOND}"
    RESULT_VARIABLE COMPARE_RESULT)
  if(NOT COMPARE_RESULT EQUAL 0)
    message(FATAL_ERROR "${RESOLUTION} pullback capture replay is not deterministic")
  endif()
endforeach()
execute_process(
  COMMAND "${PYTHON}" "${CAPTURE_SCRIPT}"
          --manifest-dir "${MANIFEST_DIR}"
          --configuration native-debug
          --backend-audit
          "${OUTPUT_DIR}/pullback_capture_96x20_first.bin"
          "${OUTPUT_DIR}/pullback_capture_288x144_first.bin"
  RESULT_VARIABLE AUDIT_RESULT)
if(NOT AUDIT_RESULT EQUAL 0)
  message(FATAL_ERROR "pullback oracle audit failed: ${AUDIT_RESULT}")
endif()
