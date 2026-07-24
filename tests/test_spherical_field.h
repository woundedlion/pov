/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <array>

#include "core/math/spherical_field.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace spherical_field {

struct Pair {
  float a;
  float b;
};

inline Pair pair_lerp(const Pair &a, const Pair &b, float t) {
  return {hs::lerp(a.a, b.a, t), hs::lerp(a.b, b.b, t)};
}

inline void test_constexpr_layout_counts() {
  constexpr hs::SphericalFieldLayout<288, 144, 0> HOST(4);
  constexpr hs::SphericalFieldLayout<288, 144, 3> DEVICE(4);
  constexpr hs::SphericalFieldLayout<288, 144, 3> FEEDBACK(4, 4, 1, 72);
  static_assert(HOST.ring_count() == 37);
  static_assert(HOST.sample_count() == 1640);
  static_assert(DEVICE.sample_count() == 1699);
  static_assert(FEEDBACK.ring_count() == 40);
  static_assert(FEEDBACK.sample_count() == 1687);
  HS_EXPECT_EQ(HOST.ring(HOST.ring_count() - 1).y, 143);
  HS_EXPECT_EQ(DEVICE.ring(DEVICE.ring_count() - 1).samples, 5);
  HS_EXPECT_EQ(FEEDBACK.ring(0).samples, 1);
  HS_EXPECT_EQ(FEEDBACK.ring(1).samples, 4);
  HS_EXPECT_EQ(FEEDBACK.ring(FEEDBACK.ring_count() - 1).samples, 5);
}

inline void test_offsets_are_contiguous() {
  constexpr hs::SphericalFieldLayout<96, 49, 0> layout(3);
  int end = 0;
  int previous_y = -1;
  for (int i = 0; i < layout.ring_count(); ++i) {
    const auto ring = layout.ring(i);
    HS_EXPECT_EQ(ring.offset, end);
    HS_EXPECT_GT(ring.y, previous_y);
    HS_EXPECT_GT(ring.samples, 0);
    end += ring.samples;
    previous_y = ring.y;
  }
  HS_EXPECT_EQ(end, layout.sample_count());
  HS_EXPECT_EQ(previous_y, 48);
}

inline void test_metric_spacing_is_uniform() {
  constexpr hs::SphericalFieldLayout<288, 144, 3> layout(4);
  constexpr float LATITUDE_STEP = 4.0f * PI_F / 146.0f;
  for (int i = 2; i + 2 < layout.ring_count(); ++i) {
    const auto ring = layout.ring(i);
    const float phi = ring.y * PI_F / 146.0f;
    const float longitude_step = 2.0f * PI_F * std::sin(phi) / ring.samples;
    const float ratio = longitude_step / LATITUDE_STEP;
    HS_EXPECT_GT(ratio, 0.96f);
    HS_EXPECT_LT(ratio, 1.05f);
  }
}

inline void test_unequal_rings_have_independent_longitude_mix() {
  constexpr hs::SphericalFieldLayout<288, 144, 3> layout(4);
  const auto lower = layout.ring(1);
  const auto upper = layout.ring(2);
  HS_EXPECT_TRUE(lower.samples != upper.samples);
  const auto a = layout.longitude(lower, 17.0f);
  const auto b = layout.longitude(upper, 17.0f);
  HS_EXPECT_TRUE(std::abs(a.mix - b.mix) > 0.05f);
}

inline void test_custom_value_interpolation_wraps_seam() {
  constexpr hs::SphericalFieldLayout<64, 33, 0> layout(4);
  std::array<Pair, layout.sample_count()> values{};
  for (int r = 0; r < layout.ring_count(); ++r) {
    const auto ring = layout.ring(r);
    for (int i = 0; i < ring.samples; ++i) {
      const float theta = 2.0f * PI_F * i / ring.samples;
      values[ring.offset + i] = {std::cos(theta), std::sin(theta)};
    }
  }
  hs::SphericalField<Pair, 64, 33, 0> field(values.data(), layout);
  const Pair left = field.interpolate(0.01f, 12.5f, pair_lerp);
  const Pair right = field.interpolate(63.99f, 12.5f, pair_lerp);
  HS_EXPECT_NEAR(left.a, right.a, 0.002f);
  HS_EXPECT_NEAR(left.b, right.b, 0.002f);
}

inline void test_populate_band_preserves_other_samples() {
  constexpr hs::SphericalFieldLayout<64, 33, 0> layout(4);
  constexpr Pair SENTINEL{-99.0f, -99.0f};
  std::array<Pair, layout.sample_count()> values;
  values.fill(SENTINEL);
  hs::SphericalField<Pair, 64, 33, 0> field(values.data(), layout);
  field.populate(2, 3, [](const Vector &v, const auto &point) {
    return Pair{dot(v, v), point.y};
  });
  const auto first = layout.ring(2);
  const auto last = layout.ring(3);
  for (int i = 0; i < layout.sample_count(); ++i) {
    if (i >= first.offset && i < last.offset + last.samples) {
      HS_EXPECT_NEAR(values[i].a, 1.0f, 0.002f);
    } else {
      HS_EXPECT_EQ(values[i].a, SENTINEL.a);
    }
  }
}

inline int run_spherical_field_tests() {
  hs_test::ModuleFixture fixture("spherical_field");
  test_constexpr_layout_counts();
  test_offsets_are_contiguous();
  test_metric_spacing_is_uniform();
  test_unequal_rings_have_independent_longitude_mix();
  test_custom_value_interpolation_wraps_seam();
  test_populate_band_preserves_other_samples();
  return fixture.result();
}

} // namespace spherical_field
} // namespace hs_test
