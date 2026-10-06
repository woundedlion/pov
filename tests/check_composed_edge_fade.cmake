cmake_minimum_required(VERSION 3.29)
function(check_fixture should_pass)
  execute_process(COMMAND "${COMPILER}" -std=gnu++20 -fsyntax-only
    "-I${ROOT}" "-I${ROOT}/core" -DHS_TIMELINE_MAX_ANIM_BYTES=256
    -DHS_GLOBAL_ARENA_BYTES=8388608 ${ARGN}
    "${ROOT}/tests/composed_edge_fade_check.cpp"
    RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error)
  if(should_pass AND NOT result EQUAL 0)
    message(FATAL_ERROR "valid edge-fade companion rejected: ${output}${error}")
  elseif(NOT should_pass)
    if(result EQUAL 0)
      message(FATAL_ERROR "edge fade without edge distance accepted")
    endif()
    if(NOT "${output}${error}" MATCHES "edge-fade coverage requires projection edge distance")
      message(FATAL_ERROR "wrong edge-fade refusal: ${output}${error}")
    endif()
  endif()
endfunction()
check_fixture(TRUE)
check_fixture(TRUE -DHS_EDGE_FADE_NO_DISTANCE -DHS_EDGE_FADE_DISABLED)
check_fixture(FALSE -DHS_EDGE_FADE_NO_DISTANCE)
message(STATUS "composed edge-fade controls passed")
