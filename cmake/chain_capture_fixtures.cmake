function(hs_chain_capture_fixtures target source_root)
  set(header "${CMAKE_CURRENT_BINARY_DIR}/chain_capture_fixtures.h")
  add_custom_command(OUTPUT "${header}"
    COMMAND "${Python3_EXECUTABLE}"
      "${source_root}/tools/gen_chain_capture_fixtures.py"
      "${source_root}/tests/data/chain_capture_fixtures.jsonl" "${header}"
    DEPENDS "${source_root}/tools/gen_chain_capture_fixtures.py"
      "${source_root}/tests/data/chain_capture_fixtures.jsonl"
    VERBATIM)
  target_sources(${target} PRIVATE "${header}")
  target_include_directories(${target} PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
endfunction()
