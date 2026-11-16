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
#include <limits>

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

/** @brief Hue-knot prefixes cover endpoint interpolation, clamps and wrapping. */
inline void test_hue_row_knots_through() {
  constexpr int CELLS = 16;
  const auto base =
      make_hue_rotate_base(Color4(Pixel(55000, 12000, 7000), 1.0f));
  for (bool cyclic : {false, true}) {
    HueRotationRow<CELLS> row;
    row.bind(base, 0.5f, cyclic);
    HS_EXPECT_EQ(row.knots_through(-0.1f), 2);
    HS_EXPECT_EQ(row.knots_through(0.0f), 2);
    HS_EXPECT_EQ(row.knots_through(0.03125f), 3);
    HS_EXPECT_EQ(row.knots_through(0.0625f), 4);
    HS_EXPECT_EQ(row.knots_through(0.1f), 5);
    HS_EXPECT_EQ(row.knots_through(0.5f), 17);
    HS_EXPECT_EQ(row.knots_through(1.0f), 17);
    for (const float maximum : {0.0f, 0.03125f, 0.1f, 0.49f, 0.5f}) {
      std::array<Pixel, CELLS + 1> partial_knots{}, full_knots{};
      HueRotationRow<CELLS> full;
      full.knots = full_knots.data();
      full.bind(base, row.domain, cyclic);
      full.prepare();
      row.knots = partial_knots.data();
      const int count = row.knots_through(maximum);
      row.prepare(count);
      for (int knot = count; knot <= CELLS; ++knot)
        HS_EXPECT_EQ(partial_knots[knot], Pixel{});
      for (int step = 0; step <= 100; ++step) {
        const float amount = maximum * (static_cast<float>(step) / 100);
        HS_EXPECT_EQ(row.sample(amount), full.sample(amount));
      }
    }
  }
  HueRotationRow<CELLS> row;
  for (const float domain :
       {0.0f, -1.0f, std::numeric_limits<float>::quiet_NaN()}) {
    row.bind(base, domain, false);
    HS_EXPECT_EQ(row.knots_through(0.1f), 17);
  }
}

/** @brief Runs shared reconstruction tests.
 * @return Module failure count. */
inline int run_dedup_tests() {
  ModuleFixture fixture("dedup");
  HS_RUN(test_hue_row_preparation_agrees);
  HS_RUN(test_hue_row_knots_through);
  return fixture.result();
}

} // namespace hs_test::dedup_tests
