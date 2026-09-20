cmake_minimum_required(VERSION 3.29)

foreach(_compiled IN ITEMS perf_bench_Os srgb_decode_gen shader_chain_catalog_gen)
  if(NOT _compiled IN_LIST COVERED_TARGETS)
    message(FATAL_ERROR "excluded_targets lost compile coverage for ${_compiled}")
  endif()
endforeach()

foreach(_utility IN ITEMS
    regenerate_mindsplatter_palette regenerate_shader_chain_catalog
    pullback_manifest_header)
  if(_utility IN_LIST COVERED_TARGETS)
    message(FATAL_ERROR "excluded_targets runs utility target ${_utility}")
  endif()
endforeach()
