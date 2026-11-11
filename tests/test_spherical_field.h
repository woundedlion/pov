/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * Spherical field topology, pole handling, reconstruction, and interpolation.
 */
#pragma once

#include <array>
#include <limits>

#include "core/math/spherical_field.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace spherical_field_tests {

inline constexpr auto bilinear_int = [](int p00, int p10, int p01, int p11,
                                        float fx, float fy) {
  const float lower =
      hs::lerp(static_cast<float>(p00), static_cast<float>(p10), fx);
  const float upper =
      hs::lerp(static_cast<float>(p01), static_cast<float>(p11), fx);
  return static_cast<int>(hs::lerp(lower, upper, fy));
};

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

/**
 * @brief Rings sit on every spacing-th row plus both endpoints, and infill
 * bands add one full-resolution ring per row.
 */
inline void test_constexpr_layout_counts() {
  constexpr int H = 144, SPACING = 4, NORTH = 4, SOUTH = 1, EQUATOR = 72;
  constexpr hs::SphericalFieldLayout<288, H, 0> HOST(SPACING);
  constexpr hs::SphericalFieldLayout<288, H, 3> FEEDBACK(SPACING, NORTH, SOUTH,
                                                         EQUATOR);
  static_assert(HOST.ring_count() == (H - 1 + SPACING - 1) / SPACING + 1);
  HS_EXPECT_EQ(HOST.ring(HOST.ring_count() - 1).y, H - 1);
  int count = 0;
  for (auto ring = FEEDBACK.ring(0);; ring = FEEDBACK.next_ring(ring)) {
    ++count;
    const bool infill = ring.y < NORTH || ring.y >= H - SOUTH;
    HS_EXPECT_TRUE(infill || ring.y % SPACING == 0);
    if (infill)
      HS_EXPECT_EQ(ring.samples, EQUATOR);
    else
      HS_EXPECT_LE(ring.samples, EQUATOR);
    if (ring.y == H - 1)
      break;
  }
  HS_EXPECT_EQ(count, FEEDBACK.ring_count());
  HS_EXPECT_EQ(count, NORTH + (H - SOUTH - 1) / SPACING + SOUTH);

  // A spacing past 4*(H+HOffset-1) rounds the derived longitude count to zero.
  constexpr hs::SphericalFieldLayout<288, 144, 3> SPARSE(4096);
  static_assert(SPARSE.ring_count() == 2);
  static_assert(SPARSE.sample_count() == 2);
  HS_EXPECT_EQ(SPARSE.ring(0).samples, 1);
  HS_EXPECT_EQ(SPARSE.longitude(SPARSE.ring(1), 287.9f).left,
               SPARSE.ring(1).offset);
}

/** @brief Offsets are contiguous. */
inline void test_offsets_are_contiguous() {
  constexpr hs::SphericalFieldLayout<96, 49, 0> LAYOUT(3);
  int end = 0;
  int previous_y = -1;
  for (int i = 0; i < LAYOUT.ring_count(); ++i) {
    const auto ring = LAYOUT.ring(i);
    HS_EXPECT_EQ(ring.offset, end);
    HS_EXPECT_GT(ring.y, previous_y);
    HS_EXPECT_GT(ring.samples, 0);
    end += ring.samples;
    previous_y = ring.y;
  }
  HS_EXPECT_EQ(end, LAYOUT.sample_count());
  HS_EXPECT_EQ(previous_y, 48);
}

/** @brief Metric spacing is uniform. */
inline void test_metric_spacing_is_uniform() {
  constexpr hs::SphericalFieldLayout<288, 144, 3> LAYOUT(4);
  constexpr float LATITUDE_STEP = 4.0f * math::PI_F / 146.0f;
  for (int i = 2; i + 2 < LAYOUT.ring_count(); ++i) {
    const auto ring = LAYOUT.ring(i);
    const float phi = ring.y * math::PI_F / 146.0f;
    const float longitude_step =
        2.0f * math::PI_F * std::sin(phi) / ring.samples;
    const float ratio = longitude_step / LATITUDE_STEP;
    HS_EXPECT_GT(ratio, 0.96f);
    HS_EXPECT_LT(ratio, 1.05f);
  }
}

/** @brief Unequal rings have independent longitude mix. */
inline void test_unequal_rings_have_independent_longitude_mix() {
  constexpr hs::SphericalFieldLayout<288, 144, 3> LAYOUT(4);
  const auto lower = LAYOUT.ring(1);
  const auto upper = LAYOUT.ring(2);
  HS_EXPECT_TRUE(lower.samples != upper.samples);
  const auto a = LAYOUT.longitude(lower, 17.0f);
  const auto b = LAYOUT.longitude(upper, 17.0f);
  HS_EXPECT_TRUE(std::abs(a.mix - b.mix) > 0.05f);
}

/** @brief Longitude stays in ring at negative seam. */
inline void test_longitude_stays_in_ring_at_negative_seam() {
  constexpr hs::SphericalFieldLayout<64, 33, 0> LAYOUT(4);
  const auto ring = LAYOUT.ring(4);
  // A tiny negative x wraps to the top of the domain.
  const auto seam = LAYOUT.longitude(ring, -1e-8f);
  HS_EXPECT_GE(seam.left, ring.offset);
  HS_EXPECT_LT(seam.left, ring.offset + ring.samples);
  HS_EXPECT_GE(seam.right, ring.offset);
  HS_EXPECT_LT(seam.right, ring.offset + ring.samples);
}

/** @brief Longitude saturates out of domain input. */
inline void test_longitude_saturates_out_of_domain_input() {
  constexpr hs::SphericalFieldLayout<64, 33, 0> LAYOUT(4);
  const auto ring = LAYOUT.ring(4);
  // volatile keeps the non-finite inputs out of constant folding.
  volatile float source[]{std::numeric_limits<float>::quiet_NaN(),
                          std::numeric_limits<float>::infinity(),
                          -std::numeric_limits<float>::infinity()};
  for (volatile float &x : source) {
    const auto seam = LAYOUT.longitude(ring, x);
    HS_EXPECT_EQ(seam.left, ring.offset + ring.samples - 1);
    HS_EXPECT_EQ(seam.right, ring.offset);
    HS_EXPECT_EQ(seam.mix, 1.0f);
  }

  // Finite longitudes wrap exactly however far out of [0, W) they start.
  const auto expected = LAYOUT.longitude(ring, 17.0f);
  for (float x : {17.0f + 262144.0f, 17.0f - 262144.0f}) {
    const auto wrapped = LAYOUT.longitude(ring, x);
    HS_EXPECT_EQ(wrapped.left, expected.left);
    HS_EXPECT_EQ(wrapped.right, expected.right);
    HS_EXPECT_EQ(wrapped.mix, expected.mix);
  }
}

/** @brief Populate band preserves other samples. */
inline void test_populate_band_preserves_other_samples() {
  constexpr hs::SphericalFieldLayout<64, 33, 0> LAYOUT(4);
  constexpr Pair SENTINEL{-99.0f, -99.0f};
  std::array<Pair, LAYOUT.sample_count()> values;
  values.fill(SENTINEL);
  hs::SphericalField<Pair, 64, 33, 0> field(values.data(), LAYOUT);
  field.populate(2, 3, [](const math::Vector &v, const auto &point) {
    return Pair{math::dot(v, v), point.y};
  });
  const auto first = LAYOUT.ring(2);
  const auto last = LAYOUT.ring(3);
  for (int i = 0; i < LAYOUT.sample_count(); ++i) {
    if (i >= first.offset && i < last.offset + last.samples) {
      HS_EXPECT_NEAR(values[i].a, 1.0f, 0.002f);
      HS_EXPECT_EQ(values[i].b, i < last.offset ? first.y : last.y);
    } else {
      HS_EXPECT_EQ(values[i].a, SENTINEL.a);
      HS_EXPECT_EQ(values[i].b, SENTINEL.b);
    }
  }
}

/** @brief Populate recurrence matches exact trig. */
inline void test_populate_recurrence_matches_exact_trig() {
  constexpr hs::SphericalFieldLayout<288, 144, 3> LAYOUT(4);
  static std::array<math::Vector, LAYOUT.sample_count()> values;
  hs::SphericalField<math::Vector, 288, 144, 3> field(values.data(), LAYOUT);
  field.populate(0, LAYOUT.ring_count() - 1,
                 [](const math::Vector &v, const auto &) { return v; });
  for (int i = 0; i < LAYOUT.ring_count(); ++i) {
    const auto ring = LAYOUT.ring(i);
    for (int sample = 0; sample < ring.samples; ++sample) {
      const math::Vector &stepped = values[ring.offset + sample];
      const math::Vector exact = LAYOUT.sample_vector(ring, sample);
      HS_EXPECT_NEAR(stepped.x, exact.x, 2e-5f);
      HS_EXPECT_NEAR(stepped.y, exact.y, 2e-5f);
      HS_EXPECT_NEAR(stepped.z, exact.z, 2e-5f);
      HS_EXPECT_NEAR(stepped.length(), 1.0f, 2e-5f);
    }
  }
}

/** @brief Longitude filter tracks spherical width. */
inline void test_longitude_filter_tracks_spherical_width() {
  constexpr hs::SphericalFieldLayout<288, 144, 3> LAYOUT(4, 4, 1, 72);
  // Every longitude of the pole row lands on one point, so it saturates at the
  // widest odd footprint 288 columns admit.
  HS_EXPECT_EQ(LAYOUT.longitude_filter_width(0), 287);
  for (int y = 1; y < 144; ++y) {
    const int width = LAYOUT.longitude_filter_width(y);
    const float equatorial_pixels = 1.0f / LAYOUT.latitude_sine(y);
    HS_EXPECT_EQ(width % 2, 1);
    HS_EXPECT_GE(float(width), equatorial_pixels);
    HS_EXPECT_LT(float(width - 2), equatorial_pixels);
  }
}

/** @brief Arm-split longitude filter keeps each half's footprint inside it,
 *  centred on its column. */
inline void test_longitude_arms_filter_within_halves() {
  constexpr int W = 64, HALF = W / 2;
  constexpr hs::SphericalFieldLayout<W, 34, 3> LAYOUT(4, 4, 1, W / 4);
  std::array<int, W> source{};
  std::array<int, W> arms{};
  std::array<int, W> wrapped{};
  source[HALF] = 1100;
  source[HALF / 2] = 900;
  LAYOUT.reconstruct_longitude_arms<IntAccumulator>(
      source.data(), 1, [&](int x, int value) { arms[x] = value; });
  LAYOUT.reconstruct_longitude_row<IntAccumulator>(
      source.data(), 1, [&](int x, int value) { wrapped[x] = value; });
  const int radius = LAYOUT.longitude_filter_width(1) / 2;
  HS_EXPECT_GT(radius, 0);
  HS_EXPECT_GT(wrapped[HALF - 1], 0);
  HS_EXPECT_EQ(arms[HALF - 1], 0);
  HS_EXPECT_EQ(arms[HALF], 1100);
  for (int x = 0; x < W; ++x) {
    const int column = x % HALF;
    if (std::min(column, HALF - 1 - column) >= radius)
      HS_EXPECT_EQ(arms[x], wrapped[x]);
  }
}

/** @brief Longitude filter and sampler wrap poles. */
inline void test_longitude_filter_and_sampler_wrap_poles() {
  constexpr int W = 64;
  constexpr hs::SphericalFieldLayout<W, 34, 3> LAYOUT(4, 4, 1, W / 4);
  std::array<int, W> source{};
  std::array<int, W> filtered{};
  source[W - 1] = 1100;
  LAYOUT.reconstruct_longitude_row<IntAccumulator>(
      source.data(), 1, [&](int x, int value) { filtered[x] = value; });
  HS_EXPECT_GT(filtered[0], 0);
  HS_EXPECT_EQ(filtered[W / 2], 0);

  // Every longitude of the pole row is one point, so the row reconstructs flat.
  LAYOUT.reconstruct_longitude_row<IntAccumulator>(
      source.data(), 0, [&](int x, int value) { filtered[x] = value; });
  for (int x = 0; x < W; ++x)
    HS_EXPECT_EQ(filtered[x], 1100 / W);

  const int SCALAR_POLES[]{-1};
  const int sample = LAYOUT.sample_bilinear(
      7.0f, -1.0f, SCALAR_POLES, 0, [](int x, int y) { return y * 100 + x; },
      bilinear_int);
  HS_EXPECT_EQ(sample, 100 + 7 + W / 2);

  std::array<Rgb, W * 34> rgb{};
  for (int y = 0; y < 34; ++y)
    for (int x = 0; x < W; ++x)
      rgb[y * W + x] = Rgb(y * 100 + x, y * 100 + x + 1, y * 100 + x + 2);
  float r, g, b;
  const Rgb POLES[]{Rgb(-1, -1, -1)};
  LAYOUT.sample_bilinear_rgb(rgb.data(), POLES, 7.0f, -1.0f, r, g, b);
  HS_EXPECT_EQ(r, 100 + 7 + W / 2);
  HS_EXPECT_EQ(g, 101 + 7 + W / 2);
  HS_EXPECT_EQ(b, 102 + 7 + W / 2);
}

/** @brief Sampler collapses south pole without virtual rows. */
inline void test_sampler_collapses_south_pole_without_virtual_rows() {
  constexpr int W = 64;
  constexpr int H = 34;
  constexpr hs::SphericalFieldLayout<W, H, 0> LAYOUT(4, 4, 4, W / 4);
  std::array<Rgb, W * H> rgb{};
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x)
      rgb[y * W + x] = Rgb(y * 100 + x, y * 100 + x + 1, y * 100 + x + 2);
  const Rgb POLES[]{Rgb(-1, -1, -1), Rgb(7000, 7001, 7002)};
  static_assert(decltype(LAYOUT)::POLE_COUNT == 2);
  float r, g, b;
  LAYOUT.sample_bilinear_rgb(rgb.data(), POLES, 7.0f, static_cast<float>(H - 1),
                             r, g, b);
  HS_EXPECT_EQ(r, POLES[1].r);
  HS_EXPECT_EQ(g, POLES[1].g);
  HS_EXPECT_EQ(b, POLES[1].b);

  // The row above the pole still loads its own samples.
  LAYOUT.sample_bilinear_rgb(rgb.data(), POLES, 7.0f, static_cast<float>(H - 2),
                             r, g, b);
  HS_EXPECT_EQ(r, (H - 2) * 100 + 7);

  const int SCALAR_POLES[]{-1, 7000};
  const int sample = LAYOUT.sample_bilinear(
      7.0f, static_cast<float>(H - 1), SCALAR_POLES, 0,
      [](int x, int y) { return y * 100 + x; }, bilinear_int);
  HS_EXPECT_EQ(sample, SCALAR_POLES[1]);
}

/** @brief Sampler wraps south pole with virtual rows. */
inline void test_sampler_wraps_south_pole_with_virtual_rows() {
  constexpr int W = 64;
  constexpr int H = 34;
  constexpr hs::SphericalFieldLayout<W, H, 3> LAYOUT(4, 4, 1, W / 4);
  // South pole sits at virtual row H + HOffset - 1 = 36, not H - 1.
  int col = 7;
  int row = 40;
  HS_EXPECT_TRUE(LAYOUT.wrap_sample(col, row));
  HS_EXPECT_EQ(row, 32);
  HS_EXPECT_EQ(col, 7 + W / 2);

  col = 7;
  row = 37;
  HS_EXPECT_FALSE(LAYOUT.wrap_sample(col, row));

  const int SCALAR_POLES[]{-1};
  const int sample = LAYOUT.sample_bilinear(
      7.0f, 40.0f, SCALAR_POLES, 0, [](int x, int y) { return y * 100 + x; },
      bilinear_int);
  HS_EXPECT_EQ(sample, 3200 + 7 + W / 2);
}

// Drives the generic sampler with the RGB sampler's own load and blend.
/** @brief Samples RGB channels independently with the generic sampler. */
template <int W, int H, int HOffset>
inline void generic_sample_rgb(
    const hs::SphericalFieldLayout<W, H, HOffset> &layout, const Rgb *source,
    const Rgb (&poles)[hs::SphericalFieldLayout<W, H, HOffset>::POLE_COUNT],
    float x, float y, float &r, float &g, float &b) {
  const Rgb outside;
  layout.sample_bilinear(
      x, y, poles, outside,
      [source](int sample_x, int sample_y) {
        return source[sample_y * W + sample_x];
      },
      [&](const Rgb &p00, const Rgb &p10, const Rgb &p01, const Rgb &p11,
          float fx, float fy) {
        const float w00 = (1.0f - fx) * (1.0f - fy);
        const float w10 = fx * (1.0f - fy);
        const float w01 = (1.0f - fx) * fy;
        const float w11 = fx * fy;
        r = p00.r * w00 + p10.r * w10 + p01.r * w01 + p11.r * w11;
        g = p00.g * w00 + p10.g * w10 + p01.g * w01 + p11.g * w11;
        b = p00.b * w00 + p10.b * w10 + p01.b * w01 + p11.b * w11;
      });
}

/** @brief Compares packed RGB sampling with independent channel sampling. */
template <int HOffset, int PoleCount>
inline void expect_rgb_sampler_matches_generic() {
  constexpr int W = 64;
  constexpr int H = 34;
  constexpr hs::SphericalFieldLayout<W, H, HOffset> LAYOUT(4, 4, 1, W / 4);
  static_assert(decltype(LAYOUT)::POLE_COUNT == PoleCount);
  std::array<Rgb, W * H> rgb{};
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x)
      rgb[y * W + x] = Rgb(y * 100 + x, y * 100 + x + 1, y * 100 + x + 2);
  Rgb poles[PoleCount];
  poles[0] = Rgb(-1, -2, -3);
  if constexpr (PoleCount == 2)
    poles[1] = Rgb(7000, 7001, 7002);

  for (int step = -4; step <= 2 * (H + HOffset); ++step) {
    const float y = step * 0.5f;
    for (const float x : {-0.5f, 0.0f, 7.25f, W - 0.5f, W + 0.75f}) {
      float r, g, b;
      LAYOUT.sample_bilinear_rgb(rgb.data(), poles, x, y, r, g, b);
      float er, eg, eb;
      generic_sample_rgb(LAYOUT, rgb.data(), poles, x, y, er, eg, eb);
      HS_EXPECT_NEAR(r, er, 1e-3f);
      HS_EXPECT_NEAR(g, eg, 1e-3f);
      HS_EXPECT_NEAR(b, eb, 1e-3f);
    }
  }
}

/** @brief Rgb sampler matches generic sampler. */
inline void test_rgb_sampler_matches_generic_sampler() {
  expect_rgb_sampler_matches_generic<0, 2>();
  expect_rgb_sampler_matches_generic<3, 1>();
}

/** @brief U16 rgb sampler error bound. */
inline void test_u16_rgb_sampler_error_bound() {
  struct U16Rgb {
    uint16_t r, g, b;
  };
  const hs::SphericalFieldLayout<4, 4, 0> LAYOUT(1);
  U16Rgb packed[16];
  Rgb floating[16];
  for (int i = 0; i < 16; ++i) {
    const uint16_t v = ((i % 4 + i / 4) % 2) ? 65535 : 0;
    packed[i] = {v, static_cast<uint16_t>(65535 - v), v};
    floating[i] = Rgb(v, 65535 - v, v);
  }
  const U16Rgb PACKED_POLES[2]{};
  const Rgb FLOATING_POLES[2]{};
  for (int ix = 0; ix < 100; ++ix)
    for (int iy = 0; iy < 100; ++iy) {
      float actual[3], expected[3];
      const float x = ix * .01007f;
      const float y = 1 + iy * .01001f;
      LAYOUT.sample_bilinear_rgb(packed, PACKED_POLES, x, y, actual[0],
                                 actual[1], actual[2]);
      LAYOUT.sample_bilinear_rgb(floating, FLOATING_POLES, x, y, expected[0],
                                 expected[1], expected[2]);
      for (int c = 0; c < 3; ++c)
        HS_EXPECT_NEAR(actual[c], expected[c], 8.0f);
    }
}

/** @brief Integer coordinates select the expected longitude and ring brackets. */
inline void test_integer_longitude_and_ring_bracketing() {
  constexpr hs::SphericalFieldLayout<288, 144, 3> LAYOUT(4, 4, 1, 72);
  for (int index = 0; index < LAYOUT.ring_count(); ++index) {
    const auto ring = LAYOUT.ring(index);
    for (int x = 0; x < 288; ++x) {
      const auto integer = LAYOUT.longitude_bounded(ring, x);
      const auto fractional = LAYOUT.longitude(ring, static_cast<float>(x));
      HS_EXPECT_EQ(integer.left, fractional.left);
      HS_EXPECT_EQ(integer.right, fractional.right);
      HS_EXPECT_NEAR(integer.mix, fractional.mix, 1e-5f);
    }
    HS_EXPECT_EQ(LAYOUT.ring_index_at_or_before(ring.y), index);
    HS_EXPECT_EQ(LAYOUT.ring_index_at_or_after(ring.y), index);
    const auto exact = LAYOUT.row(ring.y);
    HS_EXPECT_EQ(exact.lower.y, ring.y);
    HS_EXPECT_EQ(exact.mix, 0.f);
    if (index + 1 < LAYOUT.ring_count()) {
      const auto next = LAYOUT.ring(index + 1);
      const float mid = (ring.y + next.y) * .5f;
      HS_EXPECT_EQ(LAYOUT.ring_index_at_or_before(mid), index);
      HS_EXPECT_EQ(LAYOUT.ring_index_at_or_after(mid), index + 1);
      const auto between = LAYOUT.row(mid);
      HS_EXPECT_EQ(between.lower.y, ring.y);
      HS_EXPECT_EQ(between.upper.y, next.y);
      HS_EXPECT_EQ(between.mix, .5f);
    }
  }
  HS_EXPECT_EQ(LAYOUT.ring_index_at_or_after(143), LAYOUT.ring_count() - 1);
}

/** @brief Runs spherical field tests. */
inline int run_spherical_field_tests() {
  hs_test::ModuleFixture fixture("spherical_field");
  HS_RUN(test_constexpr_layout_counts);
  HS_RUN(test_integer_longitude_and_ring_bracketing);
  HS_RUN(test_offsets_are_contiguous);
  HS_RUN(test_metric_spacing_is_uniform);
  HS_RUN(test_unequal_rings_have_independent_longitude_mix);
  HS_RUN(test_longitude_stays_in_ring_at_negative_seam);
  HS_RUN(test_longitude_saturates_out_of_domain_input);
  HS_RUN(test_populate_band_preserves_other_samples);
  HS_RUN(test_populate_recurrence_matches_exact_trig);
  HS_RUN(test_longitude_filter_tracks_spherical_width);
  HS_RUN(test_longitude_filter_and_sampler_wrap_poles);
  HS_RUN(test_longitude_arms_filter_within_halves);
  HS_RUN(test_sampler_wraps_south_pole_with_virtual_rows);
  HS_RUN(test_sampler_collapses_south_pole_without_virtual_rows);
  HS_RUN(test_rgb_sampler_matches_generic_sampler);
  HS_RUN(test_u16_rgb_sampler_error_bound);
  return fixture.result();
}

} // namespace spherical_field_tests
} // namespace hs_test
