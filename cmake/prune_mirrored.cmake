# Prune obsolete engine-owned files from a validated simulator mirror.
# Required: HS_MIRROR_SOURCE, HS_DAYDREAM_DIR, HS_MIRROR_KIND,
# HS_MIRROR_ENGINE_MARKER, HS_MIRROR_SUBDIRECTORY, HS_MIRROR_GLOB (GLOB or
# GLOB_RECURSE), HS_MIRROR_PATTERNS. HS_MIRROR_EXCLUDE names an optional source.

if(NOT IS_DIRECTORY "${HS_MIRROR_SOURCE}")
  message(FATAL_ERROR "Mirror source is not a directory: ${HS_MIRROR_SOURCE}")
endif()
if(NOT EXISTS "${HS_MIRROR_SOURCE}/${HS_MIRROR_ENGINE_MARKER}")
  message(FATAL_ERROR "Mirror source is not an engine checkout: ${HS_MIRROR_SOURCE}")
endif()
if(NOT EXISTS "${HS_DAYDREAM_DIR}/src/app/daydream.js")
  message(FATAL_ERROR
    "Mirror destination is not a daydream checkout: ${HS_DAYDREAM_DIR}")
endif()
if(NOT HS_MIRROR_GLOB STREQUAL "GLOB" AND NOT HS_MIRROR_GLOB STREQUAL "GLOB_RECURSE")
  message(FATAL_ERROR "Invalid mirror glob mode: ${HS_MIRROR_GLOB}")
endif()
if(NOT HS_MIRROR_SUBDIRECTORY OR NOT HS_MIRROR_PATTERNS)
  message(FATAL_ERROR "Mirror subdirectory and patterns are required")
endif()

set(_hs_mirror_destination "${HS_DAYDREAM_DIR}/${HS_MIRROR_SUBDIRECTORY}")
if(NOT IS_DIRECTORY "${_hs_mirror_destination}")
  return()
endif()

set(_hs_source_globs)
set(_hs_installed_globs)
foreach(_hs_pattern IN LISTS HS_MIRROR_PATTERNS)
  list(APPEND _hs_source_globs "${HS_MIRROR_SOURCE}/${_hs_pattern}")
  list(APPEND _hs_installed_globs "${_hs_mirror_destination}/${_hs_pattern}")
endforeach()
file(${HS_MIRROR_GLOB} _hs_source_files LIST_DIRECTORIES FALSE
  RELATIVE "${HS_MIRROR_SOURCE}" ${_hs_source_globs})
list(REMOVE_ITEM _hs_source_files "${HS_MIRROR_EXCLUDE}")
if(NOT _hs_source_files)
  message(FATAL_ERROR "Mirror source contains no ${HS_MIRROR_KIND}: ${HS_MIRROR_SOURCE}")
endif()
file(${HS_MIRROR_GLOB} _hs_installed_files LIST_DIRECTORIES FALSE
  RELATIVE "${_hs_mirror_destination}" ${_hs_installed_globs})

foreach(_hs_relative_path IN LISTS _hs_installed_files)
  list(FIND _hs_source_files "${_hs_relative_path}" _hs_source_index)
  if(_hs_source_index EQUAL -1)
    file(REMOVE "${_hs_mirror_destination}/${_hs_relative_path}")
    message(STATUS "Removed stale mirrored file: ${_hs_relative_path}")
  endif()
endforeach()
