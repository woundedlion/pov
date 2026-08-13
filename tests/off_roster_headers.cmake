# The test headers no HS_TEST_MODULE_LIST row reaches, with the exact case-site
# count of each: helper headers included mid-module, and entry points only a
# standalone tool binary runs.
#
# ONE list, two gates. check_includes.cmake reads the names as the set exempt
# from "every test header is a run_tests.cpp include"; check_case_calls.cmake
# reads name and count as the case-site pin for every header the roster does not
# cover. Kept apart, adding a helper to one and not the other either loses case
# counting or raises a confusing orphan error, and nothing would notice.
#
# Same contract as the roster's second column: the number is exact, not a floor.
# An entry here is not an exemption from being compiled — check_includes.cmake
# still requires some source under tests/ or tools/ to include it.

set(HS_OFF_ROSTER_HEADERS
  "color_test_util.h=0"
  "mesh_test_util.h=3"
  "mindsplatter_replay_corpus.h=0"
  "mindsplatter_replay_metrics.h=0"
  "mindsplatter_whitebox.h=0"
  "pixel_test_util.h=0"
  "test_fixture.h=0"
  "test_generative_palette.h=26"
  "test_h_offset_renorm.h=5"
  "test_harness.h=0"
  "test_pole_wrap.h=3"
  "vec_test_util.h=0")

# Names only, for the gate that has no use for the counts.
set(HS_OFF_ROSTER_HEADER_NAMES "")
foreach(_off_roster_entry IN LISTS HS_OFF_ROSTER_HEADERS)
  if(NOT _off_roster_entry MATCHES "^([A-Za-z0-9_.]+)=([0-9]+)$")
    message(FATAL_ERROR
      "off_roster_headers.cmake: '${_off_roster_entry}' is not <header>=<case "
      "sites>.")
  endif()
  list(APPEND HS_OFF_ROSTER_HEADER_NAMES "${CMAKE_MATCH_1}")
endforeach()
