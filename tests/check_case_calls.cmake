# Require every test case defined in a tests/*.h or tests/*.hpp header to be
# REACHABLE from the module entry point its header defines: run_*_tests() calls
# it, or a case run_*_tests() reaches does. Off-roster helpers defined in one
# header and called from others resolve against the whole tests/ + tools/
# corpus instead.
# Cases are invoked by hand-written calls from the module's run_*_tests(), so
# dropping a call leaves a definition that still compiles — and a cluster of
# dropped cases that still call each other satisfies a plain "referenced
# somewhere" count.
#
# It sees definitions of the form `void test_*(` / `void check_*(` /
# `void case_*(` / `void verify_*(` / `void expect_*(`, plus the named
# roster-wide sweep drivers and the `int run_*_tests(` entry points, at the
# start of a line (optionally `inline`/`static` in either order, behind an
# optional single-line `template <...>` head), which is how every case in the
# tree is written; the head must start at column 0, so an indented member
# function sharing a case-name prefix is not a case. `case_*` names the death
# cases the death module's table drives, and `verify_*`/`expect_*` the
# per-item bodies a case loops a sweep over, which hold the assertions the
# loop amplifies. The name may sit on the next line, which is where
# clang-format puts it when the signature wraps. A
# module's header is the one defining the run_*_tests() entry point its roster
# row names. Headers no roster row reaches — helper headers included mid-module
# and entry points only a standalone tool binary runs — are listed in
# off_roster_headers.cmake. That list is shared with check_includes.cmake, which
# exempts the same headers from the include pin.
#
# A sweep driver one roster header defines and another calls is named in
# HS_CROSS_FILE_CASES below; its reference resolves against the whole-tree
# corpus, like an off-roster helper's.
#
# The call graph comes from splitting each header at its definition heads. A
# definition's own body runs to the first column-0 `}` — every case in the tree
# is a top-level function, so that brace is its close. Whatever follows, up to
# the next definition head, is file scope: dispatch tables, macros and the
# helper functions this scan does not name. Those references seed the closure
# alongside the entry point's, so a case only a table drives counts as reached.
#
# Comment spans and string bodies are stripped from both scan texts first: the
# tree cross-references case names in prose, and an unstripped `@details`
# mention or diagnostic message would supply the reference the real call no
# longer does.
# -D args: TESTS_DIR (path to tests/), TOOLS_DIR (path to tools/).
# Both globs match check_includes.cmake's, so the two gates see one corpus.

cmake_minimum_required(VERSION 3.29)

file(GLOB_RECURSE _headers "${TESTS_DIR}/*.h" "${TESTS_DIR}/*.hpp")

include("${TESTS_DIR}/off_roster_headers.cmake")

# Roster-wide sweep drivers: test_effects.h defines them, test_effects_smoke.h
# runs them over HS_EFFECT_LIST.
set(HS_CROSS_FILE_CASES smoke_one determinism_one clip_clear_parity_one)

# Case names the definition scan accepts.
set(_case_names "(test|check|case|verify|expect)_[A-Za-z0-9_]+")
set(_case_names "${_case_names}|(smoke|determinism|clip_clear_parity)_one")
set(_entry_name "run_[A-Za-z0-9_]*_tests")
set(_def_head "\n(template[ \t]*<[^\n]*>[ \t]*)?((inline|static)[ \t]+)*")
set(_def_head "${_def_head}(void[ \t\r\n]+(${_case_names})")
set(_def_head "${_def_head}|int[ \t\r\n]+(${_entry_name}))\\(")
# Closes the loop over a header's definition heads on the file's own tail.
set(_end_marker "\nvoid test_hs_end_of_header(")

# Whole-tree code text — every header and driver .cpp under tests/, plus the
# sources under tools/, which is where the standalone tool binaries call the
# helpers from — used as the reference scope for the off-roster helpers only.
# A roster module's cases are reached from its own run_*_tests(), so widening
# their scope this far would let a same-named token anywhere under tests/ stand
# in for the call this gate exists to enforce.
# String bodies and comments are stripped, in that order and for the same
# reasons as the per-header pass below: the helper headers are named in prose by
# several modules and in diagnostic messages, either of which would otherwise
# read as a call site.
file(GLOB_RECURSE _driver_srcs "${TESTS_DIR}/*.cpp"
  "${TOOLS_DIR}/*.h" "${TOOLS_DIR}/*.hpp" "${TOOLS_DIR}/*.cpp")
set(_corpus "")
foreach(_file IN LISTS _headers _driver_srcs)
  file(READ "${_file}" _file_text)
  string(REGEX REPLACE "\"([^\"\\\\\n]|\\\\.)*\"" "\"\"" _file_text
    "${_file_text}")
  string(REGEX REPLACE "/\\*[^*]*\\*+([^/*][^*]*\\*+)*/" "\n" _file_text
    "${_file_text}")
  string(REGEX REPLACE "//[^\n]*" "" _file_text "${_file_text}")
  string(APPEND _corpus "${_file_text}")
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
  # Only an off-roster helper is called from another file; a roster module's
  # cases must be reached from the entry point inside their own header.
  set(_cross_file FALSE)
  if(_name IN_LIST HS_OFF_ROSTER_HEADER_NAMES)
    set(_cross_file TRUE)
  endif()

  string(REGEX MATCHALL "${_def_head}" _defs "${_text}")

  set(_seen "")
  set(_entries "")
  foreach(_def IN LISTS _defs)
    string(REGEX REPLACE ".*(void|int)[ \t\r\n]+([A-Za-z0-9_]+)\\(" "\\2" _case
      "${_def}")
    if(_case IN_LIST _seen)
      continue()
    endif()
    list(APPEND _seen ${_case})
    set(_isdef_${_case} 1)
    set(_refs_${_case} "")
    if(_case MATCHES "^${_entry_name}$")
      list(APPEND _entries ${_case})
    else()
      math(EXPR _sites "${_sites} + 1")
    endif()
  endforeach()

  # Walk the heads in order, attributing each span to the definition the
  # previous head opened. The sentinel head closes the walk on the file's tail.
  string(APPEND _text "${_end_marker}")
  list(APPEND _defs "${_end_marker}")
  set(_roots "")
  set(_owner "")
  set(_rest "${_text}")
  foreach(_def IN LISTS _defs)
    string(FIND "${_rest}" "${_def}" _at)
    if(_at LESS 0)
      continue()
    endif()
    string(SUBSTRING "${_rest}" 0 ${_at} _span)
    string(LENGTH "${_def}" _dlen)
    math(EXPR _skip "${_at} + ${_dlen}")
    string(SUBSTRING "${_rest}" ${_skip} -1 _rest)
    set(_body "")
    if(NOT _owner STREQUAL "")
      string(FIND "${_span}" "\n}" _close)
      if(_close LESS 0)
        set(_body "${_span}")
        set(_span "")
      else()
        string(SUBSTRING "${_span}" 0 ${_close} _body)
        math(EXPR _after "${_close} + 2")
        string(SUBSTRING "${_span}" ${_after} -1 _span)
      endif()
    endif()
    string(REGEX MATCHALL "[A-Za-z0-9_]+" _toks "${_body}")
    list(REMOVE_DUPLICATES _toks)
    foreach(_t IN LISTS _toks)
      if(DEFINED _isdef_${_t})
        list(APPEND _refs_${_owner} "${_t}")
      endif()
    endforeach()
    string(REGEX MATCHALL "[A-Za-z0-9_]+" _toks "${_span}")
    list(REMOVE_DUPLICATES _toks)
    foreach(_t IN LISTS _toks)
      if(DEFINED _isdef_${_t})
        list(APPEND _roots "${_t}")
      endif()
    endforeach()
    string(REGEX REPLACE ".*(void|int)[ \t\r\n]+([A-Za-z0-9_]+)\\(" "\\2" _owner
      "${_def}")
  endforeach()

  # Closure from the entry points and everything file scope names.
  set(_reachable ${_entries})
  set(_frontier "${_roots}")
  foreach(_e IN LISTS _entries)
    list(APPEND _frontier ${_refs_${_e}})
  endforeach()
  list(LENGTH _frontier _flen)
  while(_flen GREATER 0)
    set(_next "")
    foreach(_n IN LISTS _frontier)
      if(_n IN_LIST _reachable)
        continue()
      endif()
      list(APPEND _reachable "${_n}")
      list(APPEND _next ${_refs_${_n}})
    endforeach()
    set(_frontier "${_next}")
    list(LENGTH _frontier _flen)
  endwhile()

  foreach(_case IN LISTS _seen)
    unset(_isdef_${_case})
    unset(_refs_${_case})
    if(_case IN_LIST _entries)
      continue()
    endif()
    if(_case IN_LIST _reachable)
      continue()
    endif()
    # An off-roster helper and a named cross-file sweep driver resolve against
    # the whole corpus: their caller is in another file. Every declaration of
    # the case is dropped first, so a forward declaration cannot stand in for
    # the call.
    if(_cross_file OR _case IN_LIST HS_CROSS_FILE_CASES)
      string(REGEX REPLACE "void[ \t\r\n]+${_case}[ \t\r\n]*\\(" "" _scan
        "${_corpus}")
      string(REGEX MATCHALL "[^A-Za-z0-9_]${_case}[^A-Za-z0-9_]" _refs "${_scan}")
      list(LENGTH _refs _nref)
      if(_nref GREATER 0)
        continue()
      endif()
    endif()
    list(APPEND _uncalled "${_name}:${_case}")
  endforeach()
endforeach()

if(_uncalled)
  string(REPLACE ";" ", " _uncalled_list "${_uncalled}")
  message(FATAL_ERROR
    "test case defined but never reached: ${_uncalled_list}. A case its "
    "module's run_*_tests() no longer reaches compiles clean and asserts "
    "nothing; restore the call or delete the case.")
endif()

message(STATUS
  "test case call check: all ${_sites} case definitions across the tests "
  "directory are reachable from their module entry point")
