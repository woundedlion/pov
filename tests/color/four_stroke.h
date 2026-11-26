/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file four_stroke.h
 * @brief Four-stroke timing, alignment, contour and coverage reference tests.
 */

#include "core/render/four_stroke_contour.h"
#include "tests/test_harness.h"
#include <numbers>

#include <vector>

namespace hs_test::color_tests {

/** @brief Each hue pair hands its roles through a four-beat stroke cycle and
 *         the beat sequence repeats with its own length. */
inline void test_four_stroke_clock() {
  using namespace FourStroke;
  constexpr uint64_t BEAT_US = 62500;
  constexpr size_t MAX_BEATS = 64;
  size_t beat_count = 1;
  while (beat_count < MAX_BEATS && &colors(beat_count * BEAT_US) != &colors(0))
    ++beat_count;
  HS_EXPECT_LT(beat_count, MAX_BEATS);
  HS_EXPECT_EQ(beat_count % 4, size_t{0});
  HS_EXPECT_GE(beat_count, size_t{8});
  if (beat_count >= MAX_BEATS || beat_count % 4 != 0)
    return;
  std::vector<Pixel> chromatic;
  for (size_t p = 0; p < beat_count / 4; ++p) {
    const FourStroke::Colors &first = colors(4 * p * BEAT_US);
    const FourStroke::Colors &rise = colors((4 * p + 1) * BEAT_US);
    const FourStroke::Colors &second = colors((4 * p + 2) * BEAT_US);
    const FourStroke::Colors &fall = colors((4 * p + 3) * BEAT_US);
    HS_EXPECT_TRUE(first.fill == rise.leading);
    HS_EXPECT_TRUE(first.fill == fall.trailing);
    HS_EXPECT_TRUE(second.fill == rise.trailing);
    HS_EXPECT_TRUE(second.fill == fall.leading);
    HS_EXPECT_FALSE(first.leading == first.trailing);
    HS_EXPECT_TRUE(first.leading == second.trailing);
    HS_EXPECT_TRUE(first.trailing == second.leading);
    chromatic.push_back(first.fill);
    chromatic.push_back(second.fill);
  }
  for (size_t i = 0; i < chromatic.size(); ++i)
    for (size_t j = i + 1; j < chromatic.size(); ++j)
      HS_EXPECT_FALSE(chromatic[i] == chromatic[j]);
  for (size_t i = 0; i < beat_count; ++i) {
    const FourStroke::Colors &beat = colors(i * BEAT_US);
    HS_EXPECT_TRUE(&beat == &colors((i + 1) * BEAT_US - 1));
    HS_EXPECT_TRUE(&beat == &colors((i + beat_count) * BEAT_US));
  }
  HS_EXPECT_EQ(column_time(62498, 62502, 1, 8, false), 62499u);
  HS_EXPECT_EQ(column_time(62498, 62502, 6, 8, false), 62500u);
  HS_EXPECT_EQ(column_time(62498, 62502, 1, 8, true), 62502u);
  HS_EXPECT_EQ(column_time(0, UINT64_MAX, 3, 8, false),
               UINT64_C(13835058055282163711));
}

/** @brief Edge polarity, width and angular contrast follow the nearest feature. */
inline void test_four_stroke_roles() {
  using namespace FourStroke;
  Controls controls;
  const FourStroke::Colors beat = {
      {1000, 2000, 3000}, {31000, 52000, 63000}, {400, 500, 600}};
  const Pixel base(7000, 8000, 9000);
  for (float degrees : {17.25f, 37.0f, 123.4f, 271.5f}) {
    auto dir = heading(degrees);
    for (int deg = 0; deg < 360; ++deg) {
      const double phi = deg * std::numbers::pi / 180;
      const float s =
          static_cast<float>(std::cos(phi) * dir.x + std::sin(phi) * dir.y);
      const double expected =
          std::cos((deg - degrees) * std::numbers::pi / 180);
      const double weight = std::max((std::abs(expected) - 0.05) / 0.95, 0.0);
      const auto pixel = shade(base, {true, 0, s}, beat, controls);
      HS_EXPECT_NEAR(pixel.r, 1000 + (expected > 0 ? 30000 : -600) * weight,
                     1.0);
    }
  }
  const auto square_strip = [](float x, float y, Direction shift) {
    return translated_strips(x, y, shift, [](float sx, float sy) {
      return sx >= 0 && sx < 10 && sy >= 0 && sy < 10;
    });
  };
  HS_EXPECT_EQ(square_strip(9.5f, 5, {1, 0}).alignment, 1);
  HS_EXPECT_EQ(square_strip(0.5f, 5, {1, 0}).alignment, -1);
  HS_EXPECT_EQ(square_strip(5, 0.1f, {1, 0}).alignment, 0);
  HS_EXPECT_EQ(square_strip(5, 9.9f, {1, 0}).alignment, 0);
  HS_EXPECT_EQ(square_strip(9.5f, 5, {-1, 0}).alignment, -1);
  HS_EXPECT_EQ(square_strip(5, 9.5f, {0, 1}).alignment, 1);
  HS_EXPECT_EQ(square_strip(0.5f, 5, {10, 0}).alignment, 1);
  HS_EXPECT_EQ(square_strip(9.5f, 5, {0, 0}).alignment, 0);
  HS_EXPECT_FALSE(square_strip(-0.1f, 5, {1, 0}).inside);

  HS_EXPECT_NEAR(interpolate_heading(359, 1, 0.5f), 0, 1e-5);
  HS_EXPECT_NEAR(interpolate_heading(180, 0, 0.5f), 270, 1e-5);
  Direction dir{};
  HS_EXPECT_TRUE(!normalize_direction(0, 0, dir));
  HS_EXPECT_TRUE(!normalize_direction(INFINITY, 1, dir));
  HS_EXPECT_TRUE(normalize_direction(std::numeric_limits<float>::max(),
                                     std::numeric_limits<float>::max(), dir));
  HS_EXPECT_NEAR(dir.x, 0.70710678f, 1e-6f);

  auto square = [&](float angle) {
    return shade(base, {true, 0, heading(angle).x}, beat, controls);
  };
  HS_EXPECT_TRUE(!(square(17) == square(37)));
  controls.mode = Mode::REFERENCE_POLARITY;
  HS_EXPECT_TRUE(square(17) == square(37));
  HS_EXPECT_TRUE(square(37) == beat.leading);
  HS_EXPECT_TRUE(square(217) == beat.trailing);
  for (int width : {1, 2}) {
    controls.width = width;
    HS_EXPECT_TRUE(shade(base, {true, static_cast<float>(width), 1}, beat,
                         controls) == beat.fill);
    HS_EXPECT_TRUE(shade(base, {true, width - 0.5f, 1}, beat, controls) ==
                   beat.leading);
    HS_EXPECT_TRUE(shade(base, {true, 0, -0.05f}, beat, controls) == beat.fill);
  }
  HS_EXPECT_TRUE(shade(base, {false, 0, 1}, beat, controls) == base);
  controls.strength = 0;
  HS_EXPECT_TRUE(shade(base, {true, 0, 1}, beat, controls) == base);
  controls.strength = 0.5f;
  HS_EXPECT_NEAR(shade(base, {true, 0, 1}, beat, controls).r, 19000, 1);
  Boundary nearest{true};
  consider(nearest, 0.5f, -1);
  consider(nearest, 0.5f, 0.2f);
  consider(nearest, 0.5f, 0.4f);
  consider(nearest, 0.5f, -0.9f);
  HS_EXPECT_EQ(nearest.alignment, 0.4f);
  consider(nearest, 0.49f, 0);
  HS_EXPECT_EQ(nearest.alignment, 0);
}

/** @brief The physical direction metric handles anisotropy, radial arrows and singular frames. */
inline void test_four_stroke_metric() {
  using namespace FourStroke;
  const math::Vector p{0, 0, 1}, dx{0.01f, 0, 0}, dy{0, 1, 0};
  const math::Vector arrow{0.04f, std::sqrt(1 - 0.04f * 0.04f), 0};
  HS_EXPECT_NEAR(tangent_alignment(p, dx, dy, arrow, {1, 0}), 0.04f, 1e-6f);
  HS_EXPECT_EQ(tangent_alignment(p, dx, dy, p, {1, 0}), 0);
  HS_EXPECT_EQ(tangent_alignment(p, {0, 0, 0}, dy, arrow, {1, 0}), 0);
  HS_EXPECT_EQ(tangent_alignment(p, {1e-6f, 0, 0}, dy, arrow, {1, 0}), 0);
  HS_EXPECT_NEAR(tangent_alignment(p, {2e-6f, 0, 0}, dy, arrow, {1, 0}), 0.04f,
                 1e-6f);
  const math::Vector skew{0.5f, 1, 0};
  HS_EXPECT_NEAR(tangent_alignment(p, dx, skew, {1, 0, 0}, {1, 0}),
                 1 / std::sqrt(1.25), 1e-5);
}

/** @brief Reconstructed boundaries preserve depth through a monotone field transfer. */
inline void test_four_stroke_contour() {
  using namespace FourStroke;
  ContourTile tile;
  auto east = [](Direction n) { return n.x; };
  for (bool cubic : {false, true}) {
    for (int offset : {-3, -2, -1, 0}) {
      int probes = 0;
      tile.prepare([&](int x, int) {
        ++probes;
        const float f = x + offset;
        return cubic ? f * f * f : f;
      });
      HS_EXPECT_EQ(probes, 49);
      for (float x : OFFSETS) {
        auto result = tile.query(x, 0, east);
        HS_EXPECT_TRUE(result.inside == (x + offset <= 0));
        if (offset < -2) {
          HS_EXPECT_TRUE(result.depth >= 2);
        } else {
          HS_EXPECT_NEAR(result.depth, std::abs(x + offset), 1e-6f);
          HS_EXPECT_NEAR(result.alignment, 1.0f, 1e-6f);
        }
      }
    }
  }
  tile.prepare([](int x, int y) { return ((x + y) % 2 == 0) ? -1.0f : 1.0f; });
  HS_EXPECT_TRUE(!tile.query(0.5f, 0.5f, east).inside);
  HS_EXPECT_TRUE(tile.query(0.125f, 0.125f, east).inside);
  tile.prepare([](int x, int y) { return ((x + y) % 2 == 0) ? -2.0f : 1.0f; });
  HS_EXPECT_TRUE(tile.query(0.5f, 0.5f, east).inside);
  tile.prepare([](int x, int) { return x == -1 || x == 1 ? -1.0f : 1.0f; });
  const auto between = tile.query(0, 0, east);
  HS_EXPECT_TRUE(!between.inside);
  HS_EXPECT_EQ(between.depth, 0.5f);
  HS_EXPECT_NEAR(between.alignment, 1.0f, 1e-6f);
  tile.prepare([](int x, int) { return x < 0 ? NAN : -1.0f; });
  HS_EXPECT_TRUE(std::isinf(tile.query(0, 0, east).depth));
}

/** @brief Subpixel resolve averages premultiplied RGB and alpha over the same footprint. */
inline void test_four_stroke_coverage() {
  using namespace FourStroke;
  const auto mixed = resolve_pixel([](float x, float y) {
    if (x > 0 && y > 0)
      return Color4({16000, 32000, 8000}, 0.5f);
    if (x < 0 && y > 0)
      return Color4({32000, 8000, 16000}, 1);
    return Color4({32000, 64000, 32000}, 0);
  });
  HS_EXPECT_EQ(mixed.alpha, 0.375f);
  HS_EXPECT_TRUE((mixed.premultiplied == Pixel(10000, 6000, 5000)));
  HS_EXPECT_TRUE(
      (over(mixed, {8000, 16000, 24000}) == Pixel(15000, 16000, 20000)));
  HS_EXPECT_TRUE(over({}, GRAY) == Pixel(CRGB(128, 128, 128)));
}

} // namespace hs_test::color_tests
