# Pin run_tests.cpp's per-module include block against its HS_TEST_MODULE_LIST
# roster. Every roster row references its run_*_tests function, so a missing
# include is already a compile error and the include set is a superset of the
# roster headers. Equal counts therefore prove there is no orphaned include left
# behind after a module was removed. Convention: every test header include in
# run_tests.cpp is a roster module.
#
# The count comparison alone cannot see a test header that exists on disk but
# is included nowhere, so the second half of this script walks the directory and
# requires every test header outside NON_MODULE_HEADERS to be included by name.
# The exemption itself is then held to the same bar: an entry must be included
# by some other source under tests/ or tools/, so listing a header there cannot
# buy it out of ever being compiled.
# -D args: SRC (path to run_tests.cpp), TESTS_DIR (path to tests/),
# TOOLS_DIR (path to tools/).

# Script mode inherits no policies from the project, so without this every
# policy defaults to OLD and the IN_LIST below is a hard "Unknown arguments"
# error (CMP0057). Matches the top-level CMakeLists.
cmake_minimum_required(VERSION 3.29)

file(READ "${SRC}" _text)

string(REGEX MATCHALL "#include \"tests/test_[A-Za-z0-9_]+\\.h(pp)?\"" _includes "${_text}")
string(REGEX MATCHALL "X\\(\"[A-Za-z0-9_]+\"" _rows "${_text}")

list(LENGTH _includes _ninc)
list(LENGTH _rows _nrow)

if(NOT _ninc EQUAL _nrow)
  message(FATAL_ERROR
    "run_tests.cpp: ${_ninc} test-module includes vs ${_nrow} HS_TEST_MODULE_LIST "
    "rows. An include without a matching roster row compiles dead test source "
    "silently; add the row or drop the orphaned include.")
endif()

# Headers that are deliberately not roster modules.
set(NON_MODULE_HEADERS
  aa_audit.h
  mesh_test_util.h
  mindsplatter_replay_corpus.h
  mindsplatter_replay_metrics.h
  mindsplatter_whitebox.h
  test_fixture.h
  test_harness.h
  test_h_offset_renorm.h
  test_pole_wrap.h)

file(GLOB_RECURSE _headers RELATIVE "${TESTS_DIR}"
  "${TESTS_DIR}/*.h"
  "${TESTS_DIR}/*.hpp")
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

# Collect every tests/ header included from a source under tests/ or tools/,
# ignoring a file's include of itself.
file(GLOB_RECURSE _sources
  "${TESTS_DIR}/*.h" "${TESTS_DIR}/*.hpp" "${TESTS_DIR}/*.cpp"
  "${TOOLS_DIR}/*.h" "${TOOLS_DIR}/*.hpp" "${TOOLS_DIR}/*.cpp")
set(_included_headers "")
foreach(_src IN LISTS _sources)
  file(READ "${_src}" _src_text)
  string(REGEX MATCHALL "#include \"tests/[A-Za-z0-9_/]+\\.h(pp)?\""
    _src_includes "${_src_text}")
  foreach(_inc IN LISTS _src_includes)
    string(REGEX REPLACE "^#include \"tests/(.*)\"$" "\\1" _target "${_inc}")
    if(NOT _src STREQUAL "${TESTS_DIR}/${_target}")
      list(APPEND _included_headers "${_target}")
    endif()
  endforeach()
endforeach()

set(_dead_exempt "")
foreach(_hdr IN LISTS NON_MODULE_HEADERS)
  if(NOT _hdr IN_LIST _included_headers)
    list(APPEND _dead_exempt ${_hdr})
  endif()
endforeach()

if(_dead_exempt)
  string(REPLACE ";" ", " _dead_list "${_dead_exempt}")
  message(FATAL_ERROR
    "NON_MODULE_HEADERS entry nothing includes: ${_dead_list}. The exemption "
    "only covers headers pulled in by another source under tests/ or tools/; "
    "an entry no one includes is never compiled, which is what this gate "
    "exists to catch. Include it, delete the file, or drop the entry.")
endif()

list(LENGTH _headers _nhdr)
message(STATUS
  "run_tests include pin: ${_ninc} includes match the roster; ${_nhdr} test "
  "directory headers all accounted for")
