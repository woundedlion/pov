/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file pattern_qualification.cpp
 * @brief Replays sphere cube seam and scalar-error qualification sweeps. */

#define HS_PATTERN_QUALIFICATION_REPORT 1
#include <cstdio>
#include "tests/pattern/sphere_lut.h"
#include "tests/test_harness.h"

/** @brief Prints measured scalar bounds and returns failed-check status.
 * @return Zero when the qualification and seam checks pass. */
int main() {
  hs_test::pattern_tests::test_sphere_bakes_match_reference();
  hs_test::pattern_tests::test_sphere_lut_all_edges();
  hs_test::pattern_tests::test_sphere_noise_bake_qualification();
  std::printf("checks_passed,%d\nchecks_failed,%d\n", hs_test::stats().passed,
              hs_test::stats().failed);
  return hs_test::stats().failed != 0;
}
