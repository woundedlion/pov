# Require every test case defined in a tests/*.h header to be referenced
# somewhere else: in the same header for module cases, or anywhere under tests/
# for the shared assertion helpers, which are defined in one header and called
# from others. Also pin how many cases each roster module defines.
# Cases are invoked by hand-written calls from the module's run_*_tests(), so
# dropping a call leaves a definition that still compiles, and deleting a case
# together with its call leaves no trace at all. The per-module assertion floor
# in run_tests.cpp catches either only when the lost assertions outweigh the
# noise: a loop-amplified module runs millions of assertions across a few hundred
# sites, so a concurrent edit that widens a sweep can lift the module back over
# its floor. This check counts sites, not assertions, and is therefore immune to
# trip-count drift; the count is exact, not a floor.
#
# It sees definitions of the form `void test_*(` / `void check_*(` at the start
# of a line (optionally `inline`/`static`), which is how every case in the tree
# is written. A module's header is the one defining the run_*_tests() entry point
# its roster row names. Headers no roster row reaches — helper headers included
# mid-module, and entry points only a standalone tool binary runs — are pinned by
# HELPER_CASE_COUNTS below instead, so every header on disk is counted by exactly
# one pin and a case deleted together with its call is always visible somewhere.
#
# Comment spans are stripped before both the definition scan and the reference
# count: the tree cross-references case names in prose, and an unstripped
# `@details` mention would supply the second reference the real call no longer
# does.
# -D args: TESTS_DIR (path to tests/).

cmake_minimum_required(VERSION 3.29)

file(GLOB _headers "${TESTS_DIR}/*.h")

# Exact case-site count of every header the run_tests.cpp roster does not reach.
# Same contract as the roster's second column: the number is exact, not a floor.
set(HELPER_CASE_COUNTS
  "aa_audit.h=0"
  "mesh_test_util.h=3"
  "mindsplatter_replay_corpus.h=0"
  "mindsplatter_replay_metrics.h=0"
  "mindsplatter_whitebox.h=0"
  "pixel_test_util.h=0"
  "test_fixture.h=0"
  "test_generative_palette.h=25"
  "test_h_offset_renorm.h=5"
  "test_harness.h=0"
  "test_pole_wrap.h=3"
  "vec_test_util.h=0")

# Whole-tree code text, used as the fallback reference scope for cross-header
# calls. Comments are stripped: the helper headers are named in prose by several
# modules, which would otherwise read as call sites.
set(_corpus "")
foreach(_hdr IN LISTS _headers)
  file(READ "${_hdr}" _hdr_text)
  string(REGEX REPLACE "/\\*[^*]*\\*+([^/*][^*]*\\*+)*/" "" _hdr_text
    "${_hdr_text}")
  string(REGEX REPLACE "//[^\n]*" "" _hdr_text "${_hdr_text}")
  string(APPEND _corpus "${_hdr_text}")
endforeach()

set(_uncalled "")
set(_sites 0)
set(_entry_fns "")   # run_*_tests() each header defines, if any
set(_entry_files "") # its header name, index-aligned with _entry_fns
set(_entry_counts "") # its case-site count, index-aligned with _entry_fns
set(_all_files "")    # every header name
set(_all_counts "")   # its case-site count, index-aligned with _all_files
foreach(_hdr IN LISTS _headers)
  file(READ "${_hdr}" _text)
  get_filename_component(_name "${_hdr}" NAME)
  # Blank string bodies, then drop block and line comment spans. Strings go
  # first so a `//` or `/*` inside a message cannot open a comment span; the
  # blanked body is bounded by its own line, an over-eager comment span is not.
  string(REGEX REPLACE "\"([^\"\\\\\n]|\\\\.)*\"" "\"\"" _text "${_text}")
  string(REGEX REPLACE "/\\*[^*]*\\*+([^/*][^*]*\\*+)*/" "\n" _text "${_text}")
  string(REGEX REPLACE "//[^\n]*" "" _text "${_text}")
  string(REGEX MATCHALL "\n[ \t]*(inline )?(static )?void (test|check)_[A-Za-z0-9_]+\\("
    _defs "${_text}")
  set(_seen "")
  set(_count 0)
  foreach(_def IN LISTS _defs)
    string(REGEX REPLACE ".*void ([A-Za-z0-9_]+)\\(" "\\1" _case "${_def}")
    if(_case IN_LIST _seen)
      continue()
    endif()
    list(APPEND _seen ${_case})
    math(EXPR _sites "${_sites} + 1")
    math(EXPR _count "${_count} + 1")
    string(REGEX MATCHALL "[^A-Za-z0-9_]${_case}[^A-Za-z0-9_]" _refs "${_text}")
    list(LENGTH _refs _nref)
    if(_nref LESS 2)
      string(REGEX MATCHALL "[^A-Za-z0-9_]${_case}[^A-Za-z0-9_]" _refs
        "${_corpus}")
      list(LENGTH _refs _nref)
    endif()
    if(_nref LESS 2)
      list(APPEND _uncalled "${_name}:${_case}")
    endif()
  endforeach()
  list(APPEND _all_files "${_name}")
  list(APPEND _all_counts "${_count}")
  string(REGEX MATCH "\n[ \t]*(inline )?(static )?int run_[A-Za-z0-9_]+_tests\\("
    _entry "${_text}")
  if(_entry)
    string(REGEX REPLACE ".*int (run_[A-Za-z0-9_]+_tests)\\(" "\\1" _entry
      "${_entry}")
    list(APPEND _entry_fns "${_entry}")
    list(APPEND _entry_files "${_name}")
    list(APPEND _entry_counts "${_count}")
  endif()
endforeach()

if(_uncalled)
  string(REPLACE ";" ", " _uncalled_list "${_uncalled}")
  message(FATAL_ERROR
    "test case defined but never called: ${_uncalled_list}. A case its module's "
    "run_*_tests() no longer calls compiles clean and asserts nothing; restore "
    "the call or delete the case.")
endif()

# Roster rows are X("<name>", <case sites>, <entry point>, <assertion floor>).
file(READ "${TESTS_DIR}/run_tests.cpp" _run_tests)
string(FIND "${_run_tests}" "#define HS_TEST_MODULE_LIST(X)" _begin)
string(FIND "${_run_tests}" "#define HS_TEST_MODULE_ENTRY" _end)
if(_begin LESS 0 OR _end LESS _begin)
  message(FATAL_ERROR
    "cannot find the HS_TEST_MODULE_LIST block in ${TESTS_DIR}/run_tests.cpp")
endif()
math(EXPR _span "${_end} - ${_begin}")
string(SUBSTRING "${_run_tests}" ${_begin} ${_span} _roster)

string(REGEX MATCHALL "X\\(\"[A-Za-z0-9_]+\"" _heads "${_roster}")
string(REGEX MATCHALL
  "X\\(\"[A-Za-z0-9_]+\",[ \t\r\n\\\\]*[0-9]+,[ \t\r\n\\\\]*[A-Za-z0-9_:]*run_[A-Za-z0-9_]+_tests"
  _rows "${_roster}")
list(LENGTH _heads _nheads)
list(LENGTH _rows _nrows)
if(NOT _nheads EQUAL _nrows)
  message(FATAL_ERROR
    "HS_TEST_MODULE_LIST has ${_nheads} rows but only ${_nrows} parse: every "
    "row must read X(\"<name>\", <case sites>, <entry point>, <assertion "
    "floor>).")
endif()

set(_drift "")
set(_roster_files "") # headers a roster row pins, so the rest need a helper pin
foreach(_row IN LISTS _rows)
  string(REGEX REPLACE "^X\\(\"([A-Za-z0-9_]+)\".*" "\\1" _mod "${_row}")
  string(REGEX REPLACE "^X\\(\"[A-Za-z0-9_]+\",[ \t\r\n\\\\]*([0-9]+),.*" "\\1"
    _pinned "${_row}")
  string(REGEX REPLACE ".*[^A-Za-z0-9_]" "" _fn "${_row}")
  list(FIND _entry_fns "${_fn}" _idx)
  if(_idx LESS 0)
    message(FATAL_ERROR
      "roster module '${_mod}' names ${_fn}, which no header matched by this "
      "check defines; the case-site count cannot be verified.")
  endif()
  list(GET _entry_counts ${_idx} _found)
  list(GET _entry_files ${_idx} _file)
  list(APPEND _roster_files "${_file}")
  if(NOT _pinned EQUAL _found)
    list(APPEND _drift
      "${_mod} (${_file}): expected ${_pinned}, found ${_found}")
  endif()
endforeach()

# Everything the roster does not reach is pinned by HELPER_CASE_COUNTS.
set(_unpinned "")
set(_index 0)
foreach(_file IN LISTS _all_files)
  list(GET _all_counts ${_index} _found)
  math(EXPR _index "${_index} + 1")
  if(_file IN_LIST _roster_files)
    continue()
  endif()
  set(_pinned "")
  foreach(_entry IN LISTS HELPER_CASE_COUNTS)
    if(_entry MATCHES "^${_file}=([0-9]+)$")
      set(_pinned "${CMAKE_MATCH_1}")
    endif()
  endforeach()
  if(_pinned STREQUAL "")
    list(APPEND _unpinned "${_file}")
  elseif(NOT _pinned EQUAL _found)
    list(APPEND _drift "${_file}: expected ${_pinned}, found ${_found}")
  endif()
endforeach()

if(_unpinned)
  string(REPLACE ";" ", " _unpinned_list "${_unpinned}")
  message(FATAL_ERROR
    "test header no case-site pin covers: ${_unpinned_list}. A header no roster "
    "row reaches must be listed in HELPER_CASE_COUNTS in "
    "tests/check_case_calls.cmake with its exact case count, or a case deleted "
    "together with its call leaves no trace.")
endif()

if(_drift)
  string(REPLACE ";" "; " _drift_list "${_drift}")
  message(FATAL_ERROR
    "case-site count drift: ${_drift_list}. The second column of each "
    "HS_TEST_MODULE_LIST row in tests/run_tests.cpp — and the HELPER_CASE_COUNTS "
    "entry in this script for every header outside the roster — is the exact "
    "number of `void test_*(`/`void check_*(` definitions in that header. If you "
    "added or intentionally removed cases, set the pin to the 'found' value "
    "above; otherwise a case was deleted — restore it.")
endif()

list(LENGTH HELPER_CASE_COUNTS _nhelpers)
message(STATUS
  "test case call pin: all ${_sites} case definitions across the tests "
  "directory are called, and ${_nrows} roster modules plus ${_nhelpers} "
  "off-roster headers match their case-site counts")
