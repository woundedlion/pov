# Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
# Licensed under the PolyForm Noncommercial License 1.0.0

set(base "${CAPTURE}.cli")
set(alias "${base}.link")
file(COPY_FILE "${CAPTURE}" "${base}")
file(SHA256 "${base}" original)
file(CREATE_LINK "${base}" "${alias}")

function(expect_rejection)
  execute_process(COMMAND "${BIN}" ${ARGV} RESULT_VARIABLE status)
  file(SHA256 "${base}" actual)
  if(NOT status STREQUAL "1" OR NOT original STREQUAL actual)
    message(FATAL_ERROR "Expected rejection with capture preserved: ${ARGV}")
  endif()
endfunction()

expect_rejection(--check-comparator "${base}" "${base}")
expect_rejection(--check-comparator "${base}" "${alias}")
foreach(frames IN ITEMS "12x" "0" "-1" "2147483648")
  expect_rejection(--dump "${base}" "${frames}")
endforeach()
file(REMOVE "${alias}" "${base}")
