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
