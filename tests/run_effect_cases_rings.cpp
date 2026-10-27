/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Rings test entry points. */
#include "core/animation/orientation.h"
#include "core/engine/engine.h"
#include "effects/RingShower.h"
#include "effects/Thrusters.h"
#include "tests/effects/fixture.h"
#include "tests/effects/case_runners.h"

namespace hs_test::effects_tests {

#include "tests/effects/ring_whitebox.h"
#include "tests/effects/rings_temporal.h"

/** @brief Runs a seeded effect case.
 * @param test Case entry point. */
static void run_case(auto test) {
  hs::random().seed(1337u);
  test();
}

void run_effect_case_ring_shower_white_box_check_temporal_rings() {
  HS_RUN(run_case, RingShowerWhiteBox::check_temporal_rings);
}

void run_effect_case_thrusters_white_box_check_temporal_rings() {
  HS_RUN(run_case, ThrustersWhiteBox::check_temporal_rings);
}

void run_effect_case_thrusters_white_box_check_warp_endpoints() {
  HS_RUN(ThrustersWhiteBox::check_warp_endpoints);
}

void run_effect_case_thrusters_white_box_check_fire_spawns_opposed_pair() {
  HS_RUN(ThrustersWhiteBox::check_fire_spawns_opposed_pair);
}

void run_effect_case_thrusters_white_box_check_collapsed_ring_falls_back_to_a_spin_axis() {
  HS_RUN(ThrustersWhiteBox::check_collapsed_ring_falls_back_to_a_spin_axis);
}

void run_effect_case_thrusters_white_box_check_fifo_evicts_the_oldest_pair() {
  HS_RUN(ThrustersWhiteBox::check_fifo_evicts_the_oldest_pair);
}

void run_effect_case_thrusters_white_box_check_expired_slots_retire_by_pair() {
  HS_RUN(ThrustersWhiteBox::check_expired_slots_retire_by_pair);
}

void run_effect_case_ring_shower_white_box_check_radius_endpoints() {
  HS_RUN(RingShowerWhiteBox::check_radius_endpoints);
}

} // namespace hs_test::effects_tests

namespace hs_test::module_runners {

/** @brief Runs the rings effect cases.
 * @return Module failure count. */
int run_rings_tests() {
  using namespace effects_tests;
  ModuleFixture fixture("effects");
  run_effect_case_ring_shower_white_box_check_temporal_rings();
  run_effect_case_thrusters_white_box_check_temporal_rings();
  run_effect_case_thrusters_white_box_check_warp_endpoints();
  run_effect_case_thrusters_white_box_check_fire_spawns_opposed_pair();
  run_effect_case_thrusters_white_box_check_collapsed_ring_falls_back_to_a_spin_axis();
  run_effect_case_thrusters_white_box_check_fifo_evicts_the_oldest_pair();
  run_effect_case_thrusters_white_box_check_expired_slots_retire_by_pair();
  run_effect_case_ring_shower_white_box_check_radius_endpoints();
  return fixture.result();
}

} // namespace hs_test::module_runners
