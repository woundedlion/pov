# Remove installed screenshots that no longer exist in the engine gallery.
# Required variables: HS_MIRROR_SOURCE, HS_DAYDREAM_DIR.

if(NOT IS_DIRECTORY "${HS_MIRROR_SOURCE}")
  message(FATAL_ERROR
    "Screenshot mirror source is not a directory: ${HS_MIRROR_SOURCE}")
endif()
if(NOT EXISTS "${HS_DAYDREAM_DIR}/daydream.js")
  message(FATAL_ERROR
    "Screenshot mirror destination is not a daydream checkout: "
    "${HS_DAYDREAM_DIR}")
endif()

set(_hs_mirror_destination "${HS_DAYDREAM_DIR}/docs/screenshots")
if(NOT IS_DIRECTORY "${_hs_mirror_destination}")
  return()
endif()

file(GLOB_RECURSE _hs_source_pngs
  LIST_DIRECTORIES FALSE
  RELATIVE "${HS_MIRROR_SOURCE}"
  "${HS_MIRROR_SOURCE}/*.png")
file(GLOB_RECURSE _hs_installed_pngs
  LIST_DIRECTORIES FALSE
  RELATIVE "${_hs_mirror_destination}"
  "${_hs_mirror_destination}/*.png")

foreach(_hs_relative_path IN LISTS _hs_installed_pngs)
  list(FIND _hs_source_pngs "${_hs_relative_path}" _hs_source_index)
  if(_hs_source_index EQUAL -1)
    file(REMOVE "${_hs_mirror_destination}/${_hs_relative_path}")
    message(STATUS "Removed stale mirrored screenshot: ${_hs_relative_path}")
  endif()
endforeach()
