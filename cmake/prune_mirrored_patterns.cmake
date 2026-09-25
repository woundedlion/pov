# Remove obsolete top-level engine-owned pattern assets.
# Required variables: HS_MIRROR_SOURCE, HS_DAYDREAM_DIR.

if(NOT IS_DIRECTORY "${HS_MIRROR_SOURCE}")
  message(FATAL_ERROR
    "Pattern mirror source is not a directory: ${HS_MIRROR_SOURCE}")
endif()
if(NOT EXISTS "${HS_MIRROR_SOURCE}/../CMakeLists.txt")
  message(FATAL_ERROR "Mirror source is not an engine checkout: ${HS_MIRROR_SOURCE}")
endif()
if(NOT EXISTS "${HS_DAYDREAM_DIR}/src/app/daydream.js")
  message(FATAL_ERROR
    "Pattern mirror destination is not a daydream checkout: "
    "${HS_DAYDREAM_DIR}")
endif()

set(_hs_mirror_destination "${HS_DAYDREAM_DIR}/generated/shader/patterns")
if(NOT IS_DIRECTORY "${_hs_mirror_destination}")
  return()
endif()

file(GLOB _hs_source_patterns
  LIST_DIRECTORIES FALSE
  RELATIVE "${HS_MIRROR_SOURCE}"
  "${HS_MIRROR_SOURCE}/*.shader.json"
  "${HS_MIRROR_SOURCE}/shaderball_migration.json")
if(NOT _hs_source_patterns)
  message(FATAL_ERROR "Mirror source contains no patterns: ${HS_MIRROR_SOURCE}")
endif()
file(GLOB _hs_installed_patterns
  LIST_DIRECTORIES FALSE
  RELATIVE "${_hs_mirror_destination}"
  "${_hs_mirror_destination}/*.shader.json"
  "${_hs_mirror_destination}/shaderball_migration.json")

foreach(_hs_relative_path IN LISTS _hs_installed_patterns)
  list(FIND _hs_source_patterns "${_hs_relative_path}" _hs_source_index)
  if(_hs_source_index EQUAL -1)
    file(REMOVE "${_hs_mirror_destination}/${_hs_relative_path}")
    message(STATUS "Removed stale mirrored pattern: ${_hs_relative_path}")
  endif()
endforeach()
