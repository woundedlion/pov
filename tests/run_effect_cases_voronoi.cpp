/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Voronoi test entry points. */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "effects/Voronoi.h"
#include "tests/effects/fixture.h"
#include "tests/effects/case_runners.h"

namespace hs_test::effects_tests {

#include "tests/effects/voronoi_cases.h"

void run_effect_case_voronoi_axes_use_uniform_sampler() {
  HS_RUN(run_case, test_voronoi_axes_use_uniform_sampler);
}

void run_effect_case_voronoi_segment_render_matches_full_frame() {
  HS_RUN(run_case, test_voronoi_segment_render_matches_full_frame);
}

void run_effect_case_voronoi_instantaneous_render_pins() {
  HS_RUN(run_case, test_voronoi_instantaneous_render_pins);
}

void run_effect_case_voronoi_endpoint_pairs_match_brute_force() {
  HS_RUN(run_case, test_voronoi_endpoint_pairs_match_brute_force);
}

void run_effect_case_voronoi_temporal_toggle() {
  HS_RUN(run_case, test_voronoi_temporal_toggle);
}

void run_effect_case_voronoi_reseed_resets_history() {
  HS_RUN(run_case, test_voronoi_reseed_resets_history);
}

void run_effect_case_voronoi_retained_corners_match_search() {
  HS_RUN(run_case, test_voronoi_retained_corners_match_search);
}

void run_effect_case_voronoi_union_candidates_cover_nearest() {
  HS_RUN(run_case, test_voronoi_union_candidates_cover_nearest);
}

void run_effect_case_voronoi_column_time_matches_brute_force() {
  HS_RUN(run_case, test_voronoi_column_time_matches_brute_force);
}

} // namespace hs_test::effects_tests

namespace hs_test::module_runners {

/** @brief Runs the voronoi effect cases.
 * @return Module failure count. */
int run_voronoi_tests() {
  using namespace effects_tests;
  ModuleFixture fixture("effects");
  run_effect_case_voronoi_axes_use_uniform_sampler();
  run_effect_case_voronoi_segment_render_matches_full_frame();
  run_effect_case_voronoi_instantaneous_render_pins();
  run_effect_case_voronoi_endpoint_pairs_match_brute_force();
  run_effect_case_voronoi_temporal_toggle();
  run_effect_case_voronoi_reseed_resets_history();
  run_effect_case_voronoi_retained_corners_match_search();
  if (effects_full_suite()) {
    run_effect_case_voronoi_union_candidates_cover_nearest();
    run_effect_case_voronoi_column_time_matches_brute_force();
  }
  return fixture.result();
}

} // namespace hs_test::module_runners
