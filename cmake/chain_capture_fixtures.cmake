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

# Engine feature set shared by the native and wasm pullback capture producers.
function(hs_capture_feature_defines target)
  target_compile_definitions(${target} PRIVATE
    HS_ENABLE_EFFECT_REGISTRY=1
    HS_ENABLE_PARAM_GUI_BRIDGE=1
    HS_ENABLE_TEST_HOOKS=1
    HS_ENABLE_TEST_ORACLES=1
    HS_ENABLE_CHAIN_INTERPRETER=1
    HS_ENABLE_EFFECT_CONTROL_API=1)
endfunction()
