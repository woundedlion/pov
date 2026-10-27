/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Dynamo test entry points. */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "effects/Dynamo.h"
#include "tests/effects/fixture.h"
#include "tests/effects/case_runners.h"

namespace hs_test::effects_tests {

#include "tests/effects/dynamo_cases.h"
#include "tests/effects/dynamo_temporal.h"

/** @brief Runs a seeded effect case.
 * @param test Case entry point. */
static void run_case(auto test) {
  hs::random().seed(1337u);
  test();
}

void run_effect_case_dynamo_column_time_oracle() {
  HS_RUN(run_case, test_dynamo_column_time_oracle);
}

void run_effect_case_dynamo_column_time_first_frame_and_reset() {
  HS_RUN(run_case, test_dynamo_column_time_first_frame_and_reset);
}

void run_effect_case_dynamo_emitted_points_counts_ring_seeds() {
  HS_RUN(run_case, test_dynamo_emitted_points_counts_ring_seeds);
}

void run_effect_case_dynamo_trail_ceiling_bounds_the_ring() {
  HS_RUN(run_case, test_dynamo_trail_ceiling_bounds_the_ring);
}

void run_effect_case_dynamo_white_box_check_gap_reduction_uses_shortest_path() {
  HS_RUN(DynamoWhiteBox::check_gap_reduction_uses_shortest_path);
}

void run_effect_case_dynamo_white_box_check_overlapping_wipes_stay_in_range() {
  HS_RUN(DynamoWhiteBox::check_overlapping_wipes_stay_in_range);
}

void run_effect_case_dynamo_white_box_check_wipe_band_colors() {
  HS_RUN(DynamoWhiteBox::check_wipe_band_colors);
}

} // namespace hs_test::effects_tests

namespace hs_test::module_runners {

/** @brief Runs the dynamo effect cases.
 * @return Module failure count. */
int run_dynamo_tests() {
  using namespace effects_tests;
  ModuleFixture fixture("effects");
  run_effect_case_dynamo_column_time_oracle();
  run_effect_case_dynamo_column_time_first_frame_and_reset();
  run_effect_case_dynamo_white_box_check_gap_reduction_uses_shortest_path();
  run_effect_case_dynamo_white_box_check_overlapping_wipes_stay_in_range();
  run_effect_case_dynamo_white_box_check_wipe_band_colors();
  run_effect_case_dynamo_emitted_points_counts_ring_seeds();
  if (effects_full_suite()) {
    run_effect_case_dynamo_trail_ceiling_bounds_the_ring();
  }
  return fixture.result();
}

} // namespace hs_test::module_runners
