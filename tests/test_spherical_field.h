/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include <limits>

#include "core/math/spherical_field.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace spherical_field {

struct Pair {
  float a;
  float b;
};

struct IntAccumulator {
  int sum = 0;
  void add(int value) { sum += value; }
  void remove(int value) { sum -= value; }
  int average(int width) const { return sum / width; }
};

struct Rgb {
  int r;
  int g;
  int b;

  constexpr Rgb(int r = 0, int g = 0, int b = 0) : r(r), g(g), b(b) {}
};

inline void test_constexpr_layout_counts() {
  constexpr hs::SphericalFieldLayout<288, 144, 0> HOST(4);
  constexpr hs::SphericalFieldLayout<288, 144, 3> DEVICE(4);
  constexpr hs::SphericalFieldLayout<288, 144, 3> FEEDBACK(4, 4, 1, 72);
  static_assert(HOST.ring_count() == 37);
  static_assert(HOST.sample_count() == 1640);
  static_assert(DEVICE.sample_count() == 1699);
  static_assert(FEEDBACK.ring_count() == 40);
  static_assert(FEEDBACK.sample_count() == 2028);
  HS_EXPECT_EQ(HOST.ring(HOST.ring_count() - 1).y, 143);
  HS_EXPECT_EQ(DEVICE.ring(DEVICE.ring_count() - 1).samples, 5);
  HS_EXPECT_EQ(FEEDBACK.ring(0).samples, 72);
  HS_EXPECT_EQ(FEEDBACK.ring(FEEDBACK.ring_count() - 1).samples, 72);
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

inline void test_longitude_stays_in_ring_at_negative_seam() {
  constexpr hs::SphericalFieldLayout<64, 33, 0> layout(4);
  const auto ring = layout.ring(4);
  // A tiny negative x wraps to the top of the domain.
  const auto seam = layout.longitude(ring, -1e-8f);
  HS_EXPECT_GE(seam.left, ring.offset);
  HS_EXPECT_LT(seam.left, ring.offset + ring.samples);
  HS_EXPECT_GE(seam.right, ring.offset);
  HS_EXPECT_LT(seam.right, ring.offset + ring.samples);
}

inline void test_longitude_saturates_out_of_domain_input() {
  constexpr hs::SphericalFieldLayout<64, 33, 0> layout(4);
  const auto ring = layout.ring(4);
  // volatile keeps the non-finite inputs out of constant folding, so the cast
  // runs at runtime where an unclamped NaN would index far outside the ring.
  volatile float source[]{std::numeric_limits<float>::quiet_NaN(),
                          std::numeric_limits<float>::infinity(),
                          -std::numeric_limits<float>::infinity()};
  for (volatile float &x : source) {
    const auto seam = layout.longitude(ring, x);
    HS_EXPECT_EQ(seam.left, ring.offset + ring.samples - 1);
    HS_EXPECT_EQ(seam.right, ring.offset);
    HS_EXPECT_EQ(seam.mix, 1.0f);
  }

  // Finite longitudes wrap exactly however far out of [0, W) they start.
  const auto expected = layout.longitude(ring, 17.0f);
  for (float x : {17.0f + 262144.0f, 17.0f - 262144.0f}) {
    const auto wrapped = layout.longitude(ring, x);
    HS_EXPECT_EQ(wrapped.left, expected.left);
    HS_EXPECT_EQ(wrapped.right, expected.right);
    HS_EXPECT_EQ(wrapped.mix, expected.mix);
  }
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

inline void test_populate_recurrence_matches_exact_trig() {
  constexpr hs::SphericalFieldLayout<288, 144, 3> layout(4);
  static std::array<Vector, layout.sample_count()> values;
  hs::SphericalField<Vector, 288, 144, 3> field(values.data(), layout);
  // The incremental rotation walks the whole ring from the meridian, so the
  // longest ring is where any drift in the recurrence accumulates.
  field.populate(0, layout.ring_count() - 1,
                 [](const Vector &v, const auto &) { return v; });
  for (int i = 0; i < layout.ring_count(); ++i) {
    const auto ring = layout.ring(i);
    for (int sample = 0; sample < ring.samples; ++sample) {
      const Vector &stepped = values[ring.offset + sample];
      const Vector exact = layout.sample_vector(ring, sample);
      HS_EXPECT_NEAR(stepped.x, exact.x, 2e-5f);
      HS_EXPECT_NEAR(stepped.y, exact.y, 2e-5f);
      HS_EXPECT_NEAR(stepped.z, exact.z, 2e-5f);
      HS_EXPECT_NEAR(stepped.length(), 1.0f, 2e-5f);
    }
  }
}

inline void test_longitude_filter_tracks_spherical_width() {
  constexpr hs::SphericalFieldLayout<288, 144, 3> layout(4, 4, 1, 72);
  // Every longitude of the pole row lands on one point, so it saturates at the
  // widest odd footprint 288 columns admit.
  HS_EXPECT_EQ(layout.longitude_filter_width(0), 287);
  HS_EXPECT_EQ(layout.longitude_filter_width(1), 47);
  HS_EXPECT_EQ(layout.longitude_filter_width(2), 25);
  HS_EXPECT_EQ(layout.longitude_filter_width(3), 17);
  HS_EXPECT_EQ(layout.longitude_filter_width(4), 13);
  HS_EXPECT_EQ(layout.longitude_filter_width(72), 3);
  HS_EXPECT_EQ(layout.longitude_filter_width(143), 17);
}

inline void test_longitude_filter_and_sampler_wrap_poles() {
  constexpr int W = 64;
  constexpr hs::SphericalFieldLayout<W, 34, 3> layout(4, 4, 1, W / 4);
  std::array<int, W> source{};
  std::array<int, W> filtered{};
  source[W - 1] = 1100;
  layout.reconstruct_longitude_row<IntAccumulator>(
      source.data(), 1, [&](int x, int value) { filtered[x] = value; });
  HS_EXPECT_GT(filtered[0], 0);
  HS_EXPECT_EQ(filtered[W / 2], 0);

  const int scalar_poles[]{-1};
  const int sample = layout.sample_bilinear(
      7.0f, -1.0f, scalar_poles, 0, [](int x, int y) { return y * 100 + x; },
      [](int p00, int p10, int p01, int p11, float fx, float fy) {
        const float lower =
            hs::lerp(static_cast<float>(p00), static_cast<float>(p10), fx);
        const float upper =
            hs::lerp(static_cast<float>(p01), static_cast<float>(p11), fx);
        return static_cast<int>(hs::lerp(lower, upper, fy));
      });
  HS_EXPECT_EQ(sample, 100 + 7 + W / 2);

  std::array<Rgb, W * 34> rgb{};
  for (int y = 0; y < 34; ++y)
    for (int x = 0; x < W; ++x)
      rgb[y * W + x] = Rgb(y * 100 + x, y * 100 + x + 1, y * 100 + x + 2);
  float r, g, b;
  const Rgb poles[]{Rgb(-1, -1, -1)};
  layout.sample_bilinear_rgb(rgb.data(), poles, 7.0f, -1.0f, r, g, b);
  HS_EXPECT_EQ(r, 100 + 7 + W / 2);
  HS_EXPECT_EQ(g, 101 + 7 + W / 2);
  HS_EXPECT_EQ(b, 102 + 7 + W / 2);
}

inline void test_sampler_collapses_south_pole_without_virtual_rows() {
  constexpr int W = 64;
  constexpr int H = 34;
  constexpr hs::SphericalFieldLayout<W, H, 0> layout(4, 4, 4, W / 4);
  std::array<Rgb, W * H> rgb{};
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x)
      rgb[y * W + x] = Rgb(y * 100 + x, y * 100 + x + 1, y * 100 + x + 2);
  const Rgb poles[]{Rgb(-1, -1, -1), Rgb(7000, 7001, 7002)};
  static_assert(decltype(layout)::POLE_COUNT == 2);
  float r, g, b;
  layout.sample_bilinear_rgb(rgb.data(), poles, 7.0f, static_cast<float>(H - 1),
                             r, g, b);
  HS_EXPECT_EQ(r, poles[1].r);
  HS_EXPECT_EQ(g, poles[1].g);
  HS_EXPECT_EQ(b, poles[1].b);

  // The row above the pole still loads its own samples.
  layout.sample_bilinear_rgb(rgb.data(), poles, 7.0f, static_cast<float>(H - 2),
                             r, g, b);
  HS_EXPECT_EQ(r, (H - 2) * 100 + 7);

  const int scalar_poles[]{-1, 7000};
  const int sample = layout.sample_bilinear(
      7.0f, static_cast<float>(H - 1), scalar_poles, 0,
      [](int x, int y) { return y * 100 + x; },
      [](int p00, int p10, int p01, int p11, float fx, float fy) {
        const float lower =
            hs::lerp(static_cast<float>(p00), static_cast<float>(p10), fx);
        const float upper =
            hs::lerp(static_cast<float>(p01), static_cast<float>(p11), fx);
        return static_cast<int>(hs::lerp(lower, upper, fy));
      });
  HS_EXPECT_EQ(sample, scalar_poles[1]);
}

inline void test_sampler_wraps_south_pole_with_virtual_rows() {
  constexpr int W = 64;
  constexpr int H = 34;
  constexpr hs::SphericalFieldLayout<W, H, 3> layout(4, 4, 1, W / 4);
  // South pole sits at virtual row H + HOffset - 1 = 36, not H - 1.
  int col = 7;
  int row = 40;
  HS_EXPECT_TRUE(layout.wrap_sample(col, row));
  HS_EXPECT_EQ(row, 32);
  HS_EXPECT_EQ(col, 7 + W / 2);

  col = 7;
  row = 37;
  HS_EXPECT_FALSE(layout.wrap_sample(col, row));

  const int scalar_poles[]{-1};
  const int sample = layout.sample_bilinear(
      7.0f, 40.0f, scalar_poles, 0, [](int x, int y) { return y * 100 + x; },
      [](int p00, int p10, int p01, int p11, float fx, float fy) {
        const float lower =
            hs::lerp(static_cast<float>(p00), static_cast<float>(p10), fx);
        const float upper =
            hs::lerp(static_cast<float>(p01), static_cast<float>(p11), fx);
        return static_cast<int>(hs::lerp(lower, upper, fy));
      });
  HS_EXPECT_EQ(sample, 3200 + 7 + W / 2);
}

inline int run_spherical_field_tests() {
  hs_test::ModuleFixture fixture("spherical_field");
  test_constexpr_layout_counts();
  test_offsets_are_contiguous();
  test_metric_spacing_is_uniform();
  test_unequal_rings_have_independent_longitude_mix();
  test_longitude_stays_in_ring_at_negative_seam();
  test_longitude_saturates_out_of_domain_input();
  test_populate_band_preserves_other_samples();
  test_populate_recurrence_matches_exact_trig();
  test_longitude_filter_tracks_spherical_width();
  test_longitude_filter_and_sampler_wrap_poles();
  test_sampler_wraps_south_pole_with_virtual_rows();
  test_sampler_collapses_south_pole_without_virtual_rows();
  return fixture.result();
}

} // namespace spherical_field
} // namespace hs_test
