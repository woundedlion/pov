cmake_minimum_required(VERSION 3.29)
execute_process(COMMAND "${TEST_BIN}"
  RESULT_VARIABLE _result OUTPUT_VARIABLE _output ERROR_VARIABLE _error)
if(NOT _result EQUAL 0)
  message(FATAL_ERROR "case-context checks failed: ${_output}${_error}")
endif()
if(NOT _output MATCHES "FAIL \\[operator\\(\\)\\].*\\{lambda_case \\| sample=7\\}")
  message(FATAL_ERROR "lambda failure lost its case or sample context: ${_output}")
endif()
