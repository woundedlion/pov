/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "tests/test_effects.h"
#include "tests/test_shader_chain.h"
#include "tests/test_composed_effect.h"
#include "tests/test_plot_parity.h"
#include "tests/test_hyper_lattice.h"
int main() {
  hs_test::ModuleFixture fixture("render_goldens");
  hs_test::plot_parity_tests::test_platform_fastmath_golden_reference();
  hs_test::plot_parity_tests::
      test_render_fingerprint_rejects_channel_and_spatial_changes();
  hs_test::effects_tests::test_bz_untimed_frames_are_pinned();
  hs_test::effects_tests::test_gs_untimed_frames_are_pinned();
  hs_test::effects_tests::test_sh_instantaneous_pixel_golden();
  hs_test::shader_chain_tests::test_shader_chain_instantaneous_pixel_golden();
  hs_test::composed_effect_tests::test_mobius_grid_instantaneous_pixel_golden();
  hs_test::plot_parity_tests::
      test_filter_pipeline_effects_match_recorded_frames();
  hs_test::effects_tests::test_voronoi_instantaneous_render_pins();
  hs_test::hyper_lattice_tests::test_instantaneous_render_pins();
  return fixture.result();
}
