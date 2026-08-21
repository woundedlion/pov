# Require every test case defined in a tests/*.h or tests/*.hpp header to be
# referenced somewhere else: in the same header for module cases, or anywhere
# under tests/ or tools/ for off-roster helpers defined in one header and called
# from others.
# Cases are invoked by hand-written calls from the module's run_*_tests(), so
# dropping a call leaves a definition that still compiles.
#
# It sees definitions of the form `void test_*(` / `void check_*(` /
# `void case_*(` / `void verify_*(` / `void expect_*(`, plus the named
# roster-wide sweep drivers, at the start of a line (optionally
# `inline`/`static` in either order, behind an optional single-line
# `template <...>` head), which is how every case in the tree is written;
# `case_*` names the death cases the death module's table drives, and
# `verify_*`/`expect_*` the per-item bodies a case loops a sweep over, which
# hold the assertions the loop amplifies. The name may sit on the next line,
# which is where clang-format puts it when the signature wraps. A module's
# header is the one defining the run_*_tests() entry point its roster row names.
# Headers no roster row reaches — helper headers included mid-module and entry
# points only a standalone tool binary runs — are listed in
# off_roster_headers.cmake. That list is shared with check_includes.cmake, which
# exempts the same headers from the include pin.
#
# A sweep driver one roster header defines and another calls is named in
# HS_CROSS_FILE_CASES below; its reference resolves against the whole-tree
# corpus, like an off-roster helper's.
#
# A reference only counts when it is not itself a declaration: `void <case>(`
# spans are removed before the count, so a forward declaration plus the
# definition no longer reads as "referenced twice" and a case whose call was
# dropped stays visible. What remains is a call or a table entry taking the
# case's address, which is how the death module drives its cases.
#
# Comment spans and string bodies are stripped from both scan texts before that
# count: the tree cross-references case names in prose, and an unstripped
# `@details` mention or diagnostic message would supply the reference the real
# call no longer does.
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

# Whole-tree code text — every header and driver .cpp under tests/, plus the
# sources under tools/, which is where the standalone tool binaries call the
# helpers from — used as the reference scope for the off-roster helpers only.
# A roster module's cases are called from its own run_*_tests(), so widening
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
  # cases must be referenced inside their own header.
  set(_cross_file FALSE)
  if(_name IN_LIST HS_OFF_ROSTER_HEADER_NAMES)
    set(_cross_file TRUE)
  endif()
  string(REGEX MATCHALL
    "\n[ \t]*(template[ \t]*<[^\n]*>[ \t]*)?((inline|static)[ \t]+)*void[ \t\r\n]+(${_case_names})\\("
    _defs "${_text}")
  set(_seen "")
  foreach(_def IN LISTS _defs)
    string(REGEX REPLACE ".*void[ \t\r\n]+([A-Za-z0-9_]+)\\(" "\\1" _case
      "${_def}")
    if(_case IN_LIST _seen)
      continue()
    endif()
    list(APPEND _seen ${_case})
    math(EXPR _sites "${_sites} + 1")
    # Drop every declaration of the case — its definition and any forward
    # declaration — then require one reference to survive. Counting raw
    # occurrences instead would let a forward declaration stand in for the call.
    string(REGEX REPLACE "void[ \t\r\n]+${_case}[ \t\r\n]*\\(" "" _scan
      "${_text}")
    string(REGEX MATCHALL "[^A-Za-z0-9_]${_case}[^A-Za-z0-9_]" _refs "${_scan}")
    list(LENGTH _refs _nref)
    if(_nref LESS 1 AND (_cross_file OR _case IN_LIST HS_CROSS_FILE_CASES))
      string(REGEX REPLACE "void[ \t\r\n]+${_case}[ \t\r\n]*\\(" "" _scan
        "${_corpus}")
      string(REGEX MATCHALL "[^A-Za-z0-9_]${_case}[^A-Za-z0-9_]" _refs
        "${_scan}")
      list(LENGTH _refs _nref)
    endif()
    if(_nref LESS 1)
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
  "test case call check: all ${_sites} case definitions across the tests "
  "directory are referenced")
