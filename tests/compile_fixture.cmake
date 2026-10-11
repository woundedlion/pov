# Syntax-only compile fixtures against COMPILER, ROOT and the DEFINES list.

# Compile ROOT/tests/<source> with the extra flags in ARGN. expected is PASS
# for a fixture that must compile, otherwise a regex the refusal must match;
# label names the contract in failure messages.
function(hs_check_compile_fixture source label expected)
  set(_define_flags)
  foreach(_definition IN LISTS DEFINES)
    list(APPEND _define_flags "-D${_definition}")
  endforeach()
  execute_process(COMMAND "${COMPILER}" -std=gnu++20 -fsyntax-only
    "-I${ROOT}" ${_define_flags} ${ARGN}
    "${ROOT}/tests/${source}"
    RESULT_VARIABLE _result OUTPUT_VARIABLE _output ERROR_VARIABLE _error)
  if(expected STREQUAL "PASS")
    if(NOT _result EQUAL 0)
      message(FATAL_ERROR "valid ${label} companion rejected: ${_output}${_error}")
    endif()
  elseif(_result EQUAL 0)
    message(FATAL_ERROR "invalid ${label} fixture accepted: ${ARGN}")
  elseif(NOT "${_output}${_error}" MATCHES "${expected}")
    message(FATAL_ERROR "wrong ${label} refusal: ${_output}${_error}")
  endif()
endfunction()
