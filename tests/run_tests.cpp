#include <cstdio>

// Pull in the engine barrel first, exactly as a real target does. geometry.h
// (included via this barrel) defines `template <int CAP = 4> class Orientation`
// with the default CAP, so filter.h and the effects can use the no-arg form
// Orientation<>. Including it here (rather than letting individual test headers
// include geometry.h directly) guarantees correct ordering for the whole
// translation unit.
#include "core/effects_engine.h"

#include "tests/test_3dmath.h"
#include "tests/test_memory.h"
#include "tests/test_spatial.h"
#include "tests/test_static_circular_buffer.h"
#include "tests/test_sdf.h"
#include "tests/test_conway.h"
#include "tests/test_hankin.h"
#include "tests/test_geometry.h"
#include "tests/test_mesh.h"
#include "tests/test_solids.h"
#include "tests/test_reaction_graph.h"
#include "tests/test_color.h"
#include "tests/test_easing_waves.h"
#include "tests/test_platform.h"
#include "tests/test_filter.h"
#include "tests/test_plot_scan.h"
#include "tests/test_canvas.h"
#include "tests/test_scan.h"
#include "tests/test_mesh_raster.h"
#include "tests/test_transformers.h"
#include "tests/test_generators.h"
#include "tests/test_animation.h"
#include "tests/test_effects.h"
#include "tests/test_hd107s_frame.h"
#include "tests/test_pov_segmented.h"
#include "tests/test_pov_single.h"
#include "tests/test_pov_sync.h"
#include "tests/test_param_marshal.h"
#include "tests/test_util.h"
#include "tests/test_presets.h"
#include "tests/test_styles.h"
#include "tests/test_death.h"

int main(int argc, char **argv) {
  std::setvbuf(stdout, nullptr, _IONBF, 0);

  hs_test::death::self_exe() = (argc > 0) ? argv[0] : nullptr;

  // Child death-case dispatch: when HS_DEATH_CASE is set, run ONLY that single
  // trap-triggering case and exit — never the full suite (which would re-spawn
  // children recursively). The case is expected to __builtin_trap(); reaching
  // the return here means it did NOT trap, so we exit 0 and let the parent flag
  // the missing trap.
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
  if (const char *dc = std::getenv("HS_DEATH_CASE")) {
#pragma clang diagnostic pop
    if (dc[0] != '\0') {
      hs_test::death::run_child_case(dc);
      return 0;
    }
  }

  int failures = 0;
  failures += hs_test::math3d::run_3dmath_tests();
  failures += hs_test::mem::run_memory_tests();
  failures += hs_test::spatial::run_spatial_tests();
  failures += hs_test::scb::run_static_circular_buffer_tests();
  failures += hs_test::sdf::run_sdf_tests();
  failures += hs_test::conway_tests::run_conway_tests();
  failures += hs_test::hankin_tests::run_hankin_tests();
  failures += hs_test::geometry::run_geometry_tests();
  failures += hs_test::mesh_tests::run_mesh_tests();
  failures += hs_test::solids_tests::run_solids_tests();
  failures += hs_test::reaction_graph_tests::run_reaction_graph_tests();
  failures += hs_test::color_tests::run_color_tests();
  failures += hs_test::easing_waves_tests::run_easing_waves_tests();
  failures += hs_test::platform_tests::run_platform_tests();
  failures += hs_test::filter_tests::run_filter_tests();
  failures += hs_test::plot_scan_tests::run_plot_scan_tests();
  failures += hs_test::canvas_tests::run_canvas_tests();
  failures += hs_test::scan_tests::run_scan_tests();
  failures += hs_test::mesh_raster_tests::run_mesh_raster_tests();
  failures += hs_test::transformers_tests::run_transformers_tests();
  failures += hs_test::generators_tests::run_generators_tests();
  failures += hs_test::animation_tests::run_animation_tests();
  failures += hs_test::effects_tests::run_effects_tests();
  failures += hs_test::hd107s_tests::run_hd107s_tests();
  failures += hs_test::pov_segmented_tests::run_pov_segmented_tests();
  failures += hs_test::pov_single_tests::run_pov_single_tests();
  failures += hs_test::pov_sync_tests::run_pov_sync_tests();
  failures += hs_test::param_marshal_tests::run_param_marshal_tests();
  failures += hs_test::util_tests::run_util_tests();
  failures += hs_test::presets_tests::run_presets_tests();
  failures += hs_test::styles_tests::run_styles_tests();
  failures += hs_test::death::run_death_tests();
  // Collapse to 0/1: a process exit status is only 8 bits on POSIX, so
  // returning a raw count would wrap (e.g. 256 failures -> 0 -> green CI).
  return failures ? 1 : 0;
}
