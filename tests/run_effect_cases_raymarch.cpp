/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Raymarch test entry points. */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "effects/Raymarch.h"
#include "tests/effects/fixture.h"
#include "tests/effects/case_runners.h"

namespace hs_test::effects_tests {

#include "tests/effects/raymarch_cases.h"

/** @brief Runs a seeded effect case.
 * @param test Case entry point. */
static void run_case(auto test) {
  hs::random().seed(1337u);
  test();
}

void run_effect_case_raymarch_volume_random_walks_are_independent() {
  HS_RUN(run_case, test_raymarch_volume_random_walks_are_independent);
}

void run_effect_case_raymarch_preset_and_placement_solids() {
  HS_RUN(run_case, test_raymarch_preset_and_placement_solids);
}

void run_effect_case_raymarch_surface_frame_uv() {
  HS_RUN(run_case, test_raymarch_surface_frame_uv);
}

void run_effect_case_raymarch_lighting_and_palette_at_column_time() {
  HS_RUN(run_case, test_raymarch_lighting_and_palette_at_column_time);
}

void run_effect_case_raymarch_hue_noise_at_column_time() {
  HS_RUN(run_case, test_raymarch_hue_noise_at_column_time);
}

void run_effect_case_raymarch_constexpr_sqrt_converges() {
  HS_RUN(run_case, test_raymarch_constexpr_sqrt_converges);
}

void run_effect_case_raymarch_unit_bounds_contains_twisted_tube() {
  HS_RUN(run_case, test_raymarch_unit_bounds_contains_twisted_tube);
}

} // namespace hs_test::effects_tests

namespace hs_test::module_runners {

/** @brief Runs the raymarch effect cases.
 * @return Module failure count. */
int run_raymarch_tests() {
  using namespace effects_tests;
  ModuleFixture fixture("effects");
  run_effect_case_raymarch_volume_random_walks_are_independent();
  run_effect_case_raymarch_preset_and_placement_solids();
  run_effect_case_raymarch_surface_frame_uv();
  run_effect_case_raymarch_lighting_and_palette_at_column_time();
  run_effect_case_raymarch_hue_noise_at_column_time();
  run_effect_case_raymarch_constexpr_sqrt_converges();
  if (effects_full_suite()) {
    run_effect_case_raymarch_unit_bounds_contains_twisted_tube();
  }
  return fixture.result();
}

} // namespace hs_test::module_runners
