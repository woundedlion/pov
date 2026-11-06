include_guard(GLOBAL)

get_filename_component(_hs_engine_root "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
set(HS_ENGINE_SOURCES
  ${_hs_engine_root}/core/memory.cpp
  ${_hs_engine_root}/core/engine/static_storage.cpp
  ${_hs_engine_root}/core/spatial/reaction_graph.cpp
  ${_hs_engine_root}/core/render/scan/mesh_capture.cpp)

set(HS_WASM_ARENA_BYTES 524288)
