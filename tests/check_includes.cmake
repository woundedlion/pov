# Pin the direct module include count against HS_TEST_MODULE_LIST, require
# every top-level test header outside off_roster_headers.cmake to be included by
# run_tests.cpp, and require every off-roster header to be included by some
# other source under tests/ or tools/.
# -D args: SRC (path to run_tests.cpp), TESTS_DIR (path to tests/),
# TOOLS_DIR (path to tools/).

# Script mode defaults every policy to OLD; IN_LIST needs CMP0057.
cmake_minimum_required(VERSION 3.29)

include("${TESTS_DIR}/header_sections.cmake")
hs_check_test_sections()

file(READ "${SRC}" _text)
string(REGEX REPLACE "/\\*[^*]*\\*+([^/*][^*]*\\*+)*/" "\n" _text "${_text}")
string(REGEX REPLACE "//[^\n]*" "" _text "${_text}")

string(REGEX MATCHALL "#include \"tests/test_[A-Za-z0-9_]+\\.h(pp)?\"" _includes "${_text}")

# Count roster rows inside the HS_TEST_MODULE_LIST block only.
string(FIND "${_text}" "#define HS_TEST_MODULE_LIST(X)" _begin)
string(FIND "${_text}" "#define HS_TEST_MODULE_ENTRY" _end)
if(_begin LESS 0 OR _end LESS _begin)
  message(FATAL_ERROR
    "cannot find the HS_TEST_MODULE_LIST block in ${SRC}")
endif()
math(EXPR _span "${_end} - ${_begin}")
string(SUBSTRING "${_text}" ${_begin} ${_span} _roster)
string(REGEX MATCHALL "X\\(\"[A-Za-z0-9_]+\"" _rows "${_roster}")

list(LENGTH _includes _ninc)
list(LENGTH _rows _nrow)

if(NOT _ninc EQUAL _nrow)
  message(FATAL_ERROR
    "run_tests.cpp: ${_ninc} test-module includes vs ${_nrow} HS_TEST_MODULE_LIST "
    "rows. An include without a matching roster row compiles dead test source "
    "silently; add the row or drop the orphaned include.")
endif()

string(REGEX MATCHALL "X\\(\"[A-Za-z0-9_]+\",[ \t\r\n\\]*[A-Za-z0-9_:]+"
  _function_rows "${_roster}")
set(_roster_functions "")
foreach(_row IN LISTS _function_rows)
  string(REGEX REPLACE ".*::(run_[A-Za-z0-9_]+_tests)$" "\\1" _function "${_row}")
  if(_function IN_LIST _roster_functions)
    message(FATAL_ERROR "duplicate module entry point in roster: ${_function}")
  endif()
  list(APPEND _roster_functions "${_function}")
endforeach()
list(LENGTH _roster_functions _nfunctions)
if(NOT _nfunctions EQUAL _nrow)
  message(FATAL_ERROR "cannot resolve all roster entry points")
endif()
set(_defined_entries "")
foreach(_inc IN LISTS _includes)
  string(REGEX REPLACE "#include \"tests/(.*)\"" "\\1" _hdr "${_inc}")
  file(READ "${TESTS_DIR}/${_hdr}" _header_text)
  string(REGEX REPLACE "/\\*[^*]*\\*+([^/*][^*]*\\*+)*/" "\n" _header_text "${_header_text}")
  string(REGEX REPLACE "//[^\n]*" "" _header_text "${_header_text}")
  string(REGEX MATCHALL "int[ \t\r\n]+run_[A-Za-z0-9_]+_tests[ \t\r\n]*\\("
    _entry_defs "${_header_text}")
  foreach(_def IN LISTS _entry_defs)
    string(REGEX REPLACE "int[ \t\r\n]+(run_[A-Za-z0-9_]+_tests)[ \t\r\n]*\\(" "\\1" _entry "${_def}")
    if(NOT _entry IN_LIST _roster_functions)
      message(FATAL_ERROR "module entry point is absent from roster: ${_hdr}:${_entry}")
    endif()
    list(APPEND _defined_entries "${_entry}")
  endforeach()
endforeach()
foreach(_function IN LISTS _roster_functions)
  if(NOT _function IN_LIST _defined_entries)
    message(FATAL_ERROR "roster entry point has no module definition: ${_function}")
  endif()
endforeach()

# Headers that are not roster modules.
include("${TESTS_DIR}/off_roster_headers.cmake")
set(NON_MODULE_HEADERS ${HS_OFF_ROSTER_HEADER_NAMES})

file(GLOB_RECURSE _headers RELATIVE "${TESTS_DIR}"
  "${TESTS_DIR}/*.h"
  "${TESTS_DIR}/*.hpp")
set(_orphans "")
foreach(_hdr IN LISTS _headers)
  if(_hdr MATCHES "^[A-Za-z0-9_]+/" OR _hdr IN_LIST NON_MODULE_HEADERS)
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
    "its roster row, delete the file, or list it in "
    "tests/off_roster_headers.cmake.")
endif()

# Collect every tests/ header included from a source under tests/ or tools/,
# ignoring a file's include of itself.
file(GLOB_RECURSE _sources
  "${TESTS_DIR}/*.h" "${TESTS_DIR}/*.hpp" "${TESTS_DIR}/*.cpp"
  "${TOOLS_DIR}/*.h" "${TOOLS_DIR}/*.hpp" "${TOOLS_DIR}/*.cpp")
set(_included_headers "")
foreach(_src IN LISTS _sources)
  file(READ "${_src}" _src_text)
  string(REGEX REPLACE "/\\*[^*]*\\*+([^/*][^*]*\\*+)*/" "\n" _src_text "${_src_text}")
  string(REGEX REPLACE "//[^\n]*" "" _src_text "${_src_text}")
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
    "off_roster_headers.cmake entry nothing includes: ${_dead_list}. The "
    "exemption only covers headers pulled in by another source under "
    "tests/ or tools/; "
    "an entry no one includes is never compiled, which is what this gate "
    "exists to catch. Include it, delete the file, or drop the entry.")
endif()

list(LENGTH _headers _nhdr)
message(STATUS
  "run_tests include pin: ${_ninc} includes match the roster; ${_nhdr} test "
  "directory headers all accounted for")
