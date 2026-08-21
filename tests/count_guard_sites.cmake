# Counts the engine's fail-fast sites — HS_CHECK and the HS_AUDIT_CHECK form the
# test build enables — so the death harness can print what fraction of them a
# death case actually pins.
#
# The count is derived here rather than written down anywhere: a hand-kept
# second list of guards is exactly the thing that drifts away from the guards.
# Included from tests/CMakeLists.txt with HS_ROOT set; sets
# HS_GUARD_SITE_ROWS (the initializer body of the generated table) and
# HS_GUARD_SITE_TOTAL, which death_guard_sites.h.in expands.
#
# Counting is per file BASENAME, because that is the granularity check_fail()
# logs and the death table pins. Comment spans are stripped first — the tree
# discusses HS_CHECK in prose — and so is the head of every #define line, which
# is where the macro and its test-build alias are written rather than used.
#
# CONFIGURE_DEPENDS on both the glob (a new guard-bearing file) and every file
# it found (a guard added to an existing one) so the count cannot go stale
# behind an incremental build.

set(_guard_dirs core effects workbench hardware targets tools)
set(_guard_files "")
foreach(_dir IN LISTS _guard_dirs)
  file(GLOB_RECURSE _found CONFIGURE_DEPENDS
       "${HS_ROOT}/${_dir}/*.h" "${HS_ROOT}/${_dir}/*.cpp"
       "${HS_ROOT}/${_dir}/*.ino")
  list(APPEND _guard_files ${_found})
endforeach()
list(SORT _guard_files)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${_guard_files})

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
  get_filename_component(_name "${_file}" NAME)
  string(REGEX MATCHALL "HS_(AUDIT_)?CHECK\\(" _hits "${_text}")
  list(LENGTH _hits _n)
  # A few traps call the reporter directly, where the macro's expression form
  # would not satisfy a [[noreturn]] tail. platform.h is where the macro and the
  # reporter are written rather than used, so its own mentions are not sites.
  if(NOT _name STREQUAL "platform.h")
    string(REGEX MATCHALL "hs::check_fail\\(" _direct "${_text}")
    list(LENGTH _direct _n_direct)
    math(EXPR _n "${_n} + ${_n_direct}")
  endif()
  if(_n EQUAL 0)
    continue()
  endif()
  math(EXPR HS_GUARD_SITE_TOTAL "${HS_GUARD_SITE_TOTAL} + ${_n}")
  list(FIND _guard_names "${_name}" _idx)
  if(_idx LESS 0)
    list(APPEND _guard_names "${_name}")
    list(APPEND _guard_counts "${_n}")
  else()
    list(GET _guard_counts ${_idx} _prev)
    math(EXPR _sum "${_prev} + ${_n}")
    list(REMOVE_AT _guard_counts ${_idx})
    list(INSERT _guard_counts ${_idx} "${_sum}")
  endif()
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
  "death-harness guard census: ${HS_GUARD_SITE_TOTAL} HS_CHECK sites across "
  "${_guard_file_count} files")
