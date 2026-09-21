execute_process(
  COMMAND "${CMAKE_COMMAND}"
    -P "${CMAKE_CURRENT_LIST_DIR}/check_mindsplatter_replay_revision_drift_case.cmake"
  RESULT_VARIABLE _rc
  ERROR_VARIABLE _err)
if(_rc EQUAL 0)
  message(FATAL_ERROR "drifted replay revisions were accepted")
endif()
string(FIND "${_err}" "MindSplatter replay revision drift" _diagnostic)
if(_diagnostic EQUAL -1)
  message(FATAL_ERROR "unexpected drift diagnostic:\n${_err}")
endif()
message(STATUS "MindSplatter replay revision drift: named rejection confirmed")
