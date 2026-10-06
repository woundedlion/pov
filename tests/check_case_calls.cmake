# Require every test case defined in a tests/*.h or tests/*.hpp header to be
# reachable from its header's run_*_tests() entry point, and to reach an
# HS_EXPECT assertion or static_assert directly or through helpers (helper
# traversal also recognizes ChainPeaks returns). run_*_cases are case drivers;
# run_*_tests are module entry points.
#
# A case is a void/bool/int/size_t definition named test_*/check_*/case_*/
# verify_*/expect_*, run_*_cases or a named sweep driver, with its head at
# column 0 (optionally `inline`/`static`, behind an optional single-line
# `template <...>`); the name may wrap onto the next line. Indented member
# definitions are not cases. Headers in off_roster_headers.cmake and the
# HS_CROSS_FILE_CASES drivers resolve against the whole tests/ + tools/ corpus.
#
# Each header is split at its definition heads; a body runs to the first
# column-0 `}`, and file-scope text after it (dispatch tables, macros) seeds the
# reachability closure alongside the entry point. Comments and string bodies are
# stripped first so prose and diagnostics cannot supply a reference.
# -D args: TESTS_DIR (path to tests/), TOOLS_DIR (path to tools/).
# CACHE_FILE is optional: skip when hashed inputs match the last passing run.

cmake_minimum_required(VERSION 3.29)

include("${TESTS_DIR}/header_sections.cmake")
hs_check_test_sections()

file(GLOB_RECURSE _headers "${TESTS_DIR}/*.h" "${TESTS_DIR}/*.hpp")

include("${TESTS_DIR}/off_roster_headers.cmake")

# Roster-wide sweep drivers called from a header other than their own.
set(HS_CROSS_FILE_CASES smoke_one determinism_one clip_clear_parity_one)

# Case names the definition scan accepts.
set(_case_names "(test|check|case|verify|expect)_[A-Za-z0-9_]+|run_[A-Za-z0-9_]*_cases")
set(_case_names "${_case_names}|(smoke|determinism|clip_clear_parity)_one")
set(_entry_name "run_[A-Za-z0-9_]*_tests")
set(_def_head "\n(template[ \t]*<[^\n]*>[ \t]*)?((inline|static)[ \t]+)*")
set(_def_head "${_def_head}((void|bool|int|size_t)[ \t\r\n]+(${_case_names})")
set(_def_head "${_def_head}|int[ \t\r\n]+(${_entry_name}))\\(")
# Closes the loop over a header's definition heads on the file's own tail.
set(_end_marker "\nvoid test_hs_end_of_header(")

# Whole-tree code text (tests/ headers and .cpp, tools/ sources), stripped of
# strings and comments: the reference scope for cross-file cases only.
file(GLOB_RECURSE _driver_srcs "${TESTS_DIR}/*.cpp"
  "${TOOLS_DIR}/*.h" "${TOOLS_DIR}/*.hpp" "${TOOLS_DIR}/*.cpp")
if(DEFINED CACHE_FILE)
  set(_cache_inputs ${_headers} ${_driver_srcs}
    "${CMAKE_CURRENT_LIST_FILE}" "${TESTS_DIR}/off_roster_headers.cmake"
    "${TESTS_DIR}/header_sections.cmake")
  list(SORT _cache_inputs)
  set(_cache_material "")
  foreach(_input IN LISTS _cache_inputs)
    file(SHA256 "${_input}" _hash)
    string(APPEND _cache_material "${_input}:${_hash}\n")
  endforeach()
  string(SHA256 _cache_key "${_cache_material}")
  if(EXISTS "${CACHE_FILE}")
    file(READ "${CACHE_FILE}" _cached_key)
    if(_cached_key STREQUAL _cache_key)
      message(STATUS "test case call check: unchanged validated corpus")
      return()
    endif()
  endif()
endif()

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

# Assertion reachability includes helper functions across test headers.
set(_assertion_functions "")
foreach(_hdr IN LISTS _headers)
  hs_read_test_header("${_hdr}" _assertion_text)
  string(REGEX REPLACE "\"([^\"\\\\\n]|\\\\.)*\"" "\"\"" _assertion_text "${_assertion_text}")
  string(REGEX REPLACE "/\\*[^*]*\\*+([^/*][^*]*\\*+)*/" "\n" _assertion_text "${_assertion_text}")
  string(REGEX REPLACE "//[^\n]*" "" _assertion_text "${_assertion_text}")
  string(REGEX MATCHALL "\n(template[ \t]*<[^\n]*>[ \t]*)?((inline|static|constexpr)[ \t]+)*(void|bool|int|size_t|ChainPeaks)[ \t\r\n]+[A-Za-z0-9_]+\\(" _heads "${_assertion_text}")
  set(_assertion_rest "${_assertion_text}")
  foreach(_head IN LISTS _heads)
    string(REGEX REPLACE ".*(void|bool|int|size_t|ChainPeaks)[ \t\r\n]+([A-Za-z0-9_]+)\\(" "\\2" _fn "${_head}")
    string(FIND "${_assertion_rest}" "${_head}" _start)
    string(LENGTH "${_head}" _head_len)
    math(EXPR _start "${_start} + ${_head_len}")
    string(SUBSTRING "${_assertion_rest}" ${_start} -1 _tail)
    set(_assertion_rest "${_tail}")
    string(REGEX MATCH "\n(template[ \t]*<[^\n]*>[ \t]*)?((inline|static|constexpr)[ \t]+)*(void|bool|int|size_t|ChainPeaks)[ \t\r\n]+[A-Za-z0-9_]+\\(" _next_head "${_tail}")
    if(NOT _next_head STREQUAL "")
      string(FIND "${_tail}" "${_next_head}" _next_start)
      string(SUBSTRING "${_tail}" 0 ${_next_start} _tail)
    endif()
    string(FIND "${_tail}" "{" _open)
    string(FIND "${_tail}" ";" _semicolon)
    if(_open LESS 0 OR (_semicolon GREATER_EQUAL 0 AND _semicolon LESS _open))
      continue()
    endif()
    string(FIND "${_tail}" "\n}" _close)
    if(_close LESS 0)
      set(_body "${_tail}")
    else()
      string(SUBSTRING "${_tail}" 0 ${_close} _body)
    endif()
    string(REGEX MATCHALL "[A-Za-z0-9_]+" _tokens "${_body}")
    list(APPEND _assertion_refs_${_fn} ${_tokens})
    list(APPEND _assertion_functions "${_fn}")
    if(_body MATCHES "(HS_EXPECT[A-Z_]*|static_assert)[ \t\r\n]*\\(")
      set(_assertion_reached_${_fn} TRUE)
    endif()
  endforeach()
endforeach()
list(REMOVE_DUPLICATES _assertion_functions)
set(_changed TRUE)
while(_changed)
  set(_changed FALSE)
  foreach(_fn IN LISTS _assertion_functions)
    if(_assertion_reached_${_fn})
      continue()
    endif()
    foreach(_ref IN LISTS _assertion_refs_${_fn})
      if(_assertion_reached_${_ref})
        set(_assertion_reached_${_fn} TRUE)
        set(_changed TRUE)
        break()
      endif()
    endforeach()
  endforeach()
endwhile()

set(_uncalled "")
set(_unasserted "")
set(_sites 0)
foreach(_hdr IN LISTS _headers)
  file(RELATIVE_PATH _relative_header "${TESTS_DIR}" "${_hdr}")
  if(_relative_header MATCHES "^[A-Za-z0-9_]+/")
    continue()
  endif()
  hs_read_test_header("${_hdr}" _text)
  get_filename_component(_name "${_hdr}" NAME)
  # Strings first, so a `//` or `/*` inside a message cannot open a comment span.
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
    string(REGEX REPLACE ".*(void|bool|int|size_t)[ \t\r\n]+([A-Za-z0-9_]+)\\(" "\\2" _case
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
    string(REGEX REPLACE ".*(void|bool|int|size_t)[ \t\r\n]+([A-Za-z0-9_]+)\\(" "\\2" _owner
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
    if(NOT (_name STREQUAL "test_death.h" AND _case MATCHES "^case_") AND
       NOT _case IN_LIST _entries AND NOT _assertion_reached_${_case})
      list(APPEND _unasserted "${_name}:${_case}")
    endif()
    unset(_isdef_${_case})
    unset(_refs_${_case})
    if(_case IN_LIST _entries)
      continue()
    endif()
    if(_case IN_LIST _reachable)
      continue()
    endif()
    # Declarations are dropped first so a forward declaration is not a call.
    if(_cross_file OR _case IN_LIST HS_CROSS_FILE_CASES)
      string(REGEX REPLACE "(void|bool|int|size_t)[ \t\r\n]+${_case}[ \t\r\n]*\\(" "" _scan
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

if(_unasserted)
  message(FATAL_ERROR "test cases reach no assertion: ${_unasserted}")
endif()

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

if(DEFINED CACHE_FILE)
  file(WRITE "${CACHE_FILE}" "${_cache_key}")
endif()
