include_guard(GLOBAL)

get_filename_component(_hs_engine_root "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
set(_hs_engine_source_list "${CMAKE_CURRENT_LIST_DIR}/engine_sources.ini")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${_hs_engine_source_list}")
set(_hs_engine_source_regex "^[ \t]+\\+<([^>]+)>[ \t]*$")
file(STRINGS "${_hs_engine_source_list}" _hs_engine_source_lines
  REGEX "${_hs_engine_source_regex}")
set(HS_ENGINE_SOURCES)
foreach(_hs_line IN LISTS _hs_engine_source_lines)
  string(REGEX REPLACE "${_hs_engine_source_regex}" "\\1" _hs_source "${_hs_line}")
  list(APPEND HS_ENGINE_SOURCES "${_hs_engine_root}/${_hs_source}")
endforeach()
if(NOT HS_ENGINE_SOURCES)
  message(FATAL_ERROR "${_hs_engine_source_list} lists no engine sources")
endif()

set(HS_WASM_ARENA_BYTES 524288)
