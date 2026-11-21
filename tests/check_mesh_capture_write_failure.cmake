# Required Notice: Copyright 2026 Gabriel Levy. All rights reserved.
# Licensed under the PolyForm Noncommercial License 1.0.0

function(expect_capture_failure output diagnostic)
  execute_process(
    COMMAND "${CMAKE_COMMAND}" -E env --unset=HS_FACE_AUDIT_PIX
      --unset=HS_FACE_AUDIT_IMAGES CI= HS_MESH_CAPTURE_ONLY=3
      HS_REPLAY_PRESET=6 HS_REPLAY_TRANS_SPEED=4 HS_REPLAY_FRAMES=1
      "${output}=/dev/full" "${TEST_BIN}" effects
    RESULT_VARIABLE result OUTPUT_VARIABLE stdout ERROR_VARIABLE stderr
    TIMEOUT 60)
  if(NOT "${result}" STREQUAL "1" OR
     NOT "${stdout}${stderr}" MATCHES "${diagnostic} audit capture must write and close")
    message(FATAL_ERROR "${output} returned ${result}: ${stdout}${stderr}")
  endif()
endfunction()

expect_capture_failure(HS_FACE_AUDIT_PIX pixel)
expect_capture_failure(HS_FACE_AUDIT_IMAGES image)
