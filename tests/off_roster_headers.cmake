# Test headers no HS_TEST_MODULE_LIST row reaches: helper headers included
# mid-module, and entry points only a standalone tool binary runs.
#
# ONE list, two gates. check_includes.cmake reads the names as the set exempt
# from "every test header is a run_tests.cpp include"; check_case_calls.cmake
# uses them to allow references from other files. An entry is not an exemption
# from being compiled: check_includes.cmake still requires some source under
# tests/ or tools/ to include it.

set(HS_OFF_ROSTER_HEADER_NAMES
  "color_test_util.h"
  "mesh_test_util.h"
  "mindsplatter_replay_corpus.h"
  "mindsplatter_replay_metrics.h"
  "mindsplatter_whitebox.h"
  "pixel_test_util.h"
  "test_fixture.h"
  "test_generative_palette.h"
  "test_h_offset_renorm.h"
  "test_harness.h"
  "test_pole_wrap.h"
  "vec_test_util.h")
