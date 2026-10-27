/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Galaxies test entry points. */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "effects/Galaxies.h"
#include "tests/effects/fixture.h"
#include "tests/effects/case_runners.h"

namespace hs_test::effects_tests {

#include "tests/effects/galaxies_tests.h"

/** @brief Runs a seeded effect case.
 * @param test Case entry point. */
static void run_case(auto test) {
  hs::random().seed(1337u);
  test();
}

void run_effect_case_galaxies_spawn_along_spiral_with_orbital_velocity() {
  HS_RUN(run_case, test_galaxies_spawn_along_spiral_with_orbital_velocity);
}

void run_effect_case_galaxies_alternate_orbits_and_arm_spin() {
  HS_RUN(run_case, test_galaxies_alternate_orbits_and_arm_spin);
}

void run_effect_case_galaxies_arm_contrast_and_age_fade() {
  HS_RUN(run_case, test_galaxies_arm_contrast_and_age_fade);
}

void run_effect_case_galaxies_arm_sharpness() {
  HS_RUN(run_case, test_galaxies_arm_sharpness);
}

void run_effect_case_galaxies_star_color_mix() {
  HS_RUN(run_case, test_galaxies_star_color_mix);
}

void run_effect_case_galaxies_core_fade() {
  HS_RUN(run_case, test_galaxies_core_fade);
}

void run_effect_case_galaxies_white_star_stays_white() {
  HS_RUN(run_case, test_galaxies_white_star_stays_white);
}

void run_effect_case_galaxies_emission_rate_and_capacity() {
  HS_RUN(run_case, test_galaxies_emission_rate_and_capacity);
}

void run_effect_case_galaxies_live_particle_count() {
  HS_RUN(run_case, test_galaxies_live_particle_count);
}

void run_effect_case_galaxies_black_hole_core_toggle() {
  HS_RUN(run_case, test_galaxies_black_hole_core_toggle);
}

void run_effect_case_galaxies_packed_orbits_match_float() {
  HS_RUN(run_case, test_galaxies_packed_orbits_match_float);
}

} // namespace hs_test::effects_tests

namespace hs_test::module_runners {

/** @brief Runs the galaxies effect cases.
 * @return Module failure count. */
int run_galaxies_tests() {
  using namespace effects_tests;
  ModuleFixture fixture("effects");
  run_effect_case_galaxies_spawn_along_spiral_with_orbital_velocity();
  run_effect_case_galaxies_alternate_orbits_and_arm_spin();
  run_effect_case_galaxies_arm_contrast_and_age_fade();
  run_effect_case_galaxies_arm_sharpness();
  run_effect_case_galaxies_star_color_mix();
  run_effect_case_galaxies_core_fade();
  run_effect_case_galaxies_white_star_stays_white();
  run_effect_case_galaxies_emission_rate_and_capacity();
  run_effect_case_galaxies_live_particle_count();
  run_effect_case_galaxies_black_hole_core_toggle();
  run_effect_case_galaxies_packed_orbits_match_float();
  return fixture.result();
}

} // namespace hs_test::module_runners
