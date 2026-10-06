# Counts the engine's fail-fast sites (HS_CHECK, HS_AUDIT_CHECK and direct
# check_fail calls) per repository-relative source path, so the death harness
# can print what fraction of them a death case pins.
# Run at build time with HS_ROOT and HS_GUARD_OUTPUT set.
# HS_GUARD_SITE_ROWS and HS_GUARD_SITE_TOTAL expand death_guard_sites.h.in.
# Comments, string bodies and #define lines are stripped before counting.
#
if(NOT DEFINED HS_GUARD_DIRS)
  file(STRINGS "${CMAKE_CURRENT_LIST_DIR}/guard_directories.txt" HS_GUARD_DIRS)
endif()
set(_guard_files "")
foreach(_dir IN LISTS HS_GUARD_DIRS)
  file(GLOB_RECURSE _found
       "${HS_ROOT}/${_dir}/*.h" "${HS_ROOT}/${_dir}/*.cpp"
       "${HS_ROOT}/${_dir}/*.ino")
  list(APPEND _guard_files ${_found})
endforeach()
list(SORT _guard_files)

set(_guard_names "")
set(_guard_counts "")
set(HS_GUARD_SITE_TOTAL 0)
foreach(_file IN LISTS _guard_files)
  file(READ "${_file}" _text)
  # Strings first, so a `//` or `/*` inside a literal cannot open a comment span
  # and swallow the guards that follow it.
  string(REGEX REPLACE "\"([^\"\\\\\n]|\\\\.)*\"" "\"\"" _text "${_text}")
  string(REGEX REPLACE "/\\*[^*]*\\*+([^/*][^*]*\\*+)*/" "" _text "${_text}")
  string(REGEX REPLACE "//[^\n]*" "" _text "${_text}")
  string(REGEX REPLACE "#[ \t]*define[^\n]*" "" _text "${_text}")
  file(RELATIVE_PATH _name "${HS_ROOT}" "${_file}")
  string(REGEX MATCHALL "HS_(AUDIT_)?CHECK\\(" _hits "${_text}")
  list(LENGTH _hits _n)
  # platform.h defines the reporter; its own mentions are not sites.
  if(NOT _name STREQUAL "core/platform/platform.h")
    string(REGEX MATCHALL "(hs::)?check_fail\\(" _direct "${_text}")
    list(LENGTH _direct _n_direct)
    math(EXPR _n "${_n} + ${_n_direct}")
  endif()
  if(_n EQUAL 0)
    continue()
  endif()
  math(EXPR HS_GUARD_SITE_TOTAL "${HS_GUARD_SITE_TOTAL} + ${_n}")
  list(APPEND _guard_names "${_name}")
  list(APPEND _guard_counts "${_n}")
endforeach()

set(HS_GUARD_SITE_ROWS "")
set(_row 0)
foreach(_name IN LISTS _guard_names)
  list(GET _guard_counts ${_row} _count)
  math(EXPR _row "${_row} + 1")
  string(APPEND HS_GUARD_SITE_ROWS "    {\"${_name}\", ${_count}},\n")
endforeach()

list(LENGTH _guard_names _guard_file_count)
message(STATUS
  "death-harness guard census: ${HS_GUARD_SITE_TOTAL} fail-fast sites across "
  "${_guard_file_count} files")

configure_file("${HS_ROOT}/tests/death_guard_sites.h.in"
               "${HS_GUARD_OUTPUT}" @ONLY)
