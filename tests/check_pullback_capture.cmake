set(FIRST "${OUTPUT_DIR}/pullback_capture_first.bin")
set(SECOND "${OUTPUT_DIR}/pullback_capture_second.bin")
execute_process(
  COMMAND "${PRODUCER}" --resolution 96x20 "${FIRST}"
  RESULT_VARIABLE FIRST_RESULT)
if(NOT FIRST_RESULT EQUAL 0)
  message(FATAL_ERROR "first pullback capture failed: ${FIRST_RESULT}")
endif()
execute_process(
  COMMAND "${PRODUCER}" --resolution 96x20 "${SECOND}"
  RESULT_VARIABLE SECOND_RESULT)
if(NOT SECOND_RESULT EQUAL 0)
  message(FATAL_ERROR "second pullback capture failed: ${SECOND_RESULT}")
endif()
execute_process(
  COMMAND "${CMAKE_COMMAND}" -E compare_files "${FIRST}" "${SECOND}"
  RESULT_VARIABLE COMPARE_RESULT)
if(NOT COMPARE_RESULT EQUAL 0)
  message(FATAL_ERROR "pullback capture replay is not deterministic")
endif()
