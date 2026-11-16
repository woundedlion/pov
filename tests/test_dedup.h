/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file test_dedup.h
 * @brief Shared hue-knot evaluation through eager and lazy row preparation. */

#include "core/color/hue_rotation_row.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include <array>

namespace hs_test::dedup_tests {

/** @brief Eager and lazy hue rows reproduce the same bound hue domain. */
inline void test_hue_row_preparation_agrees() {
  constexpr int CELLS = 16;
  for (const Pixel color :
       {Pixel(55000, 12000, 7000), Pixel(1000, 30000, 60000),
        Pixel(21000, 21000, 21000)}) {
    const auto base = make_hue_rotate_base(Color4(color, 1.0f));
    for (const float domain : {0.35f, 1.0f}) {
      for (const bool wraps : {false, true}) {
        std::array<Pixel, CELLS + 1> eager_knots{}, lazy_knots{};
        HueRotationRow<CELLS> eager, lazy;
        eager.knots = eager_knots.data();
        lazy.knots = lazy_knots.data();
        eager.bind(base, domain, wraps);
        lazy.bind(base, domain, wraps);
        eager.prepare();
        for (int i = 0; i <= CELLS; ++i) {
          const Pixel expected =
              hue_rotate(base, domain * (static_cast<float>(i) / CELLS)).color;
          HS_EXPECT_EQ(eager_knots[i], expected);
        }
        for (int step = -8; step <= 72; ++step) {
          const float amount = domain * (static_cast<float>(step) / 32.0f);
          HS_EXPECT_EQ(lazy.sample_lazy(amount), eager.sample(amount));
          HS_EXPECT_EQ(lazy.sample_lazy(amount), eager.sample(amount));
        }
      }
    }
  }
}

/** @brief Runs shared reconstruction tests.
 * @return Module failure count. */
inline int run_dedup_tests() {
  ModuleFixture fixture("dedup");
  HS_RUN(test_hue_row_preparation_agrees);
  return fixture.result();
}

} // namespace hs_test::dedup_tests
