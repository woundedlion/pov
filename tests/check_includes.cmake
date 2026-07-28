# Pin run_tests.cpp's per-module include block against its HS_TEST_MODULE_LIST
# roster. Every roster row references its run_*_tests function, so a missing
# include is already a compile error and the include set is a superset of the
# roster headers. Equal counts therefore prove there is no orphaned include left
# behind after a module was removed. Convention: every tests/test_*.h include in
# run_tests.cpp is a roster module.
#
# The count comparison alone cannot see a tests/test_*.h that exists on disk but
# is included nowhere, so the second half of this script walks the directory and
# requires every test header outside NON_MODULE_HEADERS to be included by name.
# -D args: SRC (path to run_tests.cpp), TESTS_DIR (path to tests/).

file(READ "${SRC}" _text)

string(REGEX MATCHALL "#include \"tests/test_[A-Za-z0-9_]+\\.h\"" _includes "${_text}")
string(REGEX MATCHALL "X\\(\"[A-Za-z0-9_]+\"" _rows "${_text}")

list(LENGTH _includes _ninc)
list(LENGTH _rows _nrow)

if(NOT _ninc EQUAL _nrow)
  message(FATAL_ERROR
    "run_tests.cpp: ${_ninc} test-module includes vs ${_nrow} HS_TEST_MODULE_LIST "
    "rows. An include without a matching roster row compiles dead test source "
    "silently; add the row or drop the orphaned include.")
endif()

# tests/test_*.h files that are deliberately not modules: shared harness code
# pulled in by the module headers, and the separate-TU offset-3 suite compiled
# by h_offset_renorm_check.cpp.
set(NON_MODULE_HEADERS test_fixture.h test_harness.h test_h_offset_renorm.h)

file(GLOB _headers RELATIVE "${TESTS_DIR}" "${TESTS_DIR}/test_*.h")
set(_orphans "")
foreach(_hdr IN LISTS _headers)
  if(_hdr IN_LIST NON_MODULE_HEADERS)
    continue()
  endif()
  if(NOT _text MATCHES "#include \"tests/${_hdr}\"")
    list(APPEND _orphans ${_hdr})
  endif()
endforeach()

if(_orphans)
  string(REPLACE ";" ", " _orphan_list "${_orphans}")
  message(FATAL_ERROR
    "run_tests.cpp includes no such header: ${_orphan_list}. A test header on "
    "disk that nothing includes is never compiled or run; add the include and "
    "its roster row, delete the file, or list it in NON_MODULE_HEADERS.")
endif()

list(LENGTH _headers _nhdr)
message(STATUS
  "run_tests include pin: ${_ninc} includes match the roster; ${_nhdr} test "
  "headers on disk all accounted for")
