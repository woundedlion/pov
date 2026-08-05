# Require every test case defined in a tests/*.h header to be referenced
# somewhere else: in the same header for module cases, or anywhere under tests/
# for the shared assertion helpers, which are defined in one header and called
# from others. Cases are invoked by hand-written calls from the module's
# run_*_tests(), so dropping a call leaves a definition that still compiles.
# The per-module assertion floor in run_tests.cpp catches that only
# when the lost assertions outweigh the noise: a loop-amplified module runs
# millions of assertions across a few hundred sites, so a concurrent edit that
# widens a sweep can lift the module back over its floor. This check counts
# sites, not assertions, and is therefore immune to trip-count drift.
#
# It sees definitions of the form `void test_*(` / `void check_*(` at the start
# of a line (optionally `inline`/`static`), which is how every case in the tree
# is written; a case removed together with its definition is left to code review.
#
# Comment spans are stripped before both the definition scan and the reference
# count: the tree cross-references case names in prose, and an unstripped
# `@details` mention would supply the second reference the real call no longer
# does.
# -D args: TESTS_DIR (path to tests/).

cmake_minimum_required(VERSION 3.29)

file(GLOB _headers "${TESTS_DIR}/*.h")

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
  foreach(_def IN LISTS _defs)
    string(REGEX REPLACE ".*void ([A-Za-z0-9_]+)\\(" "\\1" _case "${_def}")
    if(_case IN_LIST _seen)
      continue()
    endif()
    list(APPEND _seen ${_case})
    math(EXPR _sites "${_sites} + 1")
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
endforeach()

if(_uncalled)
  string(REPLACE ";" ", " _uncalled_list "${_uncalled}")
  message(FATAL_ERROR
    "test case defined but never called: ${_uncalled_list}. A case its module's "
    "run_*_tests() no longer calls compiles clean and asserts nothing; restore "
    "the call or delete the case.")
endif()

message(STATUS
  "test case call pin: all ${_sites} case definitions across the tests "
  "directory are called")
