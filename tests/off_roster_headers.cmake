# Test headers no HS_TEST_MODULE_LIST row reaches: helper headers included
# mid-module, and entry points a standalone tool binary runs. An entry must
# still be included by some source under tests/ or tools/.

set(HS_OFF_ROSTER_HEADER_NAMES
  "color_test_util.h"
  "composed_frame_fixture.h"
  "composed_chain_fixture.h"
  "conway_test_util.h"
  "fd_capture_util.h"
  "mesh_test_util.h"
  "mindsplatter_replay_corpus.h"
  "mindsplatter_replay_metrics.h"
  "mindsplatter_whitebox.h"
  "pixel_test_util.h"
  "pole_geometry_test_util.h"
  "pov_tiling_test_util.h"
  "test_fixture.h"
  "test_generative_palette.h"
  "test_h_offset_renorm.h"
  "test_harness.h"
  "test_lattice_trace.h"
  "test_pole_wrap.h"
  "vec_test_util.h"
  "volume_reference.h")
