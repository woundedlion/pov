/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Spherical harmonics test entry points. */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "effects/SphericalHarmonics.h"
#include "tests/effects/fixture.h"
#include "tests/effects/case_runners.h"

namespace hs_test::effects_tests {

#include "tests/effects/spherical_harmonics.h"

/** @brief Runs a seeded effect case.
 * @param test Case entry point. */
static void run_case(auto test) {
  hs::random().seed(1337u);
  test();
}

void run_effect_case_sh_preset_mode_mapping() {
  HS_RUN(run_case, test_sh_preset_mode_mapping);
}

void run_effect_case_sh_manual_preset_replaces_inflight_morph() {
  HS_RUN(run_case, test_sh_manual_preset_replaces_inflight_morph);
}

void run_effect_case_sh_field_blend_endpoints_and_frame() {
  HS_RUN(run_case, test_sh_field_blend_endpoints_and_frame);
}

void run_effect_case_sh_field_stays_inside_unit_range() {
  HS_RUN(run_case, test_sh_field_stays_inside_unit_range);
}

void run_effect_case_sh_polarity_split_and_ao_shaping() {
  HS_RUN(run_case, test_sh_polarity_split_and_ao_shaping);
}

void run_effect_case_sh_cycle_covers_every_mode() {
  HS_RUN(run_case, test_sh_cycle_covers_every_mode);
}

void run_effect_case_sh_instantaneous_pixel_golden() {
  HS_RUN(run_case, test_sh_instantaneous_pixel_golden);
}

void run_effect_case_sh_columns_follow_sweep_time() {
  HS_RUN(run_case, test_sh_columns_follow_sweep_time);
}

void run_effect_case_sh_segmented_columns_match_single_engine() {
  HS_RUN(run_case, test_sh_segmented_columns_match_single_engine);
}

void run_effect_case_sh_manual_preset_resets_morph_history() {
  HS_RUN(run_case, test_sh_manual_preset_resets_morph_history);
}

void run_effect_case_sh_fidelity_toggle_keeps_history() {
  HS_RUN(run_case, test_sh_fidelity_toggle_keeps_history);
}

void run_effect_case_sh_frame_sampler_rollover_and_sign() {
  HS_RUN(run_case, test_sh_frame_sampler_rollover_and_sign);
}

void run_effect_case_sh_pullback_matches_direct_colorize() {
  HS_RUN(run_case, test_sh_pullback_matches_direct_colorize);
}

void run_effect_case_sh_morph_chain_rearms() {
  HS_RUN(run_case, test_sh_morph_chain_rearms);
}

} // namespace hs_test::effects_tests

namespace hs_test::module_runners {

/** @brief Runs the spherical harmonics effect cases.
 * @return Module failure count. */
int run_spherical_harmonics_effect_tests() {
  using namespace effects_tests;
  ModuleFixture fixture("effects");
  run_effect_case_sh_preset_mode_mapping();
  run_effect_case_sh_manual_preset_replaces_inflight_morph();
  run_effect_case_sh_field_blend_endpoints_and_frame();
  run_effect_case_sh_field_stays_inside_unit_range();
  run_effect_case_sh_polarity_split_and_ao_shaping();
  run_effect_case_sh_cycle_covers_every_mode();
  run_effect_case_sh_instantaneous_pixel_golden();
  run_effect_case_sh_columns_follow_sweep_time();
  run_effect_case_sh_segmented_columns_match_single_engine();
  run_effect_case_sh_manual_preset_resets_morph_history();
  run_effect_case_sh_fidelity_toggle_keeps_history();
  run_effect_case_sh_frame_sampler_rollover_and_sign();
  if (effects_full_suite()) {
    run_effect_case_sh_pullback_matches_direct_colorize();
    run_effect_case_sh_morph_chain_rearms();
  }
  return fixture.result();
}

} // namespace hs_test::module_runners
