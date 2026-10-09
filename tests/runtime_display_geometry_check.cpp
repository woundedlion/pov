/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#include "core/engine/engine.h"
#include "core/render/filter/screen_blur.h"
#include "core/math/spherical_field.h"
#include "effects/MeshFeedback.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include "tests/pole_geometry_test_util.h"
#include <limits>

namespace {
constexpr int FEEDBACK_FRAMES = 4;
template <int W, int H> void check_geometry(float north, float south) {
  using Geometry = math::DisplayGeometry<H>;
  using Field = hs::SphericalFieldLayout<W, H>;
  math::init_geometry_luts<W, H>();
  HS_EXPECT_NEAR(math::y_to_phi<H>(0), north, 1e-6f);
  HS_EXPECT_NEAR(math::y_to_phi<H>(H - 1), south, 1e-6f);
  HS_EXPECT_NEAR(math::RADIANS_PER_ROW<H>, (south - north) / (H - 1), 1e-6f);
  HS_EXPECT_EQ(Field::HAS_NORTH_POLE, north == 0.0f);
  HS_EXPECT_EQ(Field::HAS_SOUTH_POLE, south == math::PI_F);
  constexpr int DOWNSAMPLE = ::Feedback::Style{}.downsample;
  constexpr Field CACHE_FIELD(DOWNSAMPLE, DOWNSAMPLE, DOWNSAMPLE,
                              W / DOWNSAMPLE);
  constexpr size_t CACHE_CELLS = (W / DOWNSAMPLE) * CACHE_FIELD.ring_count();
  constexpr size_t EXPECTED_STORAGE =
      CACHE_CELLS * (4 * sizeof(int16_t) + sizeof(typename Field::Coordinates));
  static_assert(HS_RUNTIME_DISPLAY_GEOMETRY);
  HS_EXPECT_EQ((Filter::Pixel::Feedback<W, H>::STORAGE_BYTES),
               EXPECTED_STORAGE);
  const Field field(4);
  for (int row : {0, H - 1}) {
    const bool pole = row == 0 ? Field::HAS_NORTH_POLE : Field::HAS_SOUTH_POLE;
    if (!pole)
      for (float alpha : {0.5f, 1.0f})
        hs_test::pole_geometry::check_feedback_ring_centroid<W, H>(row, 0, 1,
                                                                   alpha);
  }
  const float poles[Field::POLE_STORAGE_COUNT] = {30.0f, 40.0f};
  for (int row : {0, H - 1}) {
    const float phi = row == 0 ? north : south;
    const bool pole = row == 0 ? Field::HAS_NORTH_POLE : Field::HAS_SOUTH_POLE;
    const auto point = math::pixel_to_vector<W, H>(0, row);
    HS_EXPECT_NEAR(point.y, std::cos(phi), 1e-6f);
    HS_EXPECT_NEAR(point.x, std::sin(phi), 1e-6f);
    HS_EXPECT_NEAR(Geometry::phi_to_row(phi), static_cast<float>(row), 3e-5f);
    const float sampled = field.sample_bilinear(
        3.0f, static_cast<float>(row), poles, 0.0f,
        [](int, int) { return 10.0f; },
        [](float a, float, float, float, float, float) { return a; });
    const float expected_pole =
        row == 0 || !Field::HAS_NORTH_POLE ? 30.0f : 40.0f;
    HS_EXPECT_NEAR(sampled, pole ? expected_pole : 10.0f, 1e-6f);
    const auto taps =
        Filter::Screen::splat_taps<W, H>(3.5f, row == 0 ? -0.5f : H - 0.5f);
    HS_EXPECT_NEAR(taps.v00 + taps.v10 + taps.v01 + taps.v11,
                   pole ? 1.0f : 0.5f, 1e-6f);
    Filter::Screen::Blur<W, H> blur(1.0f);
    float energy = 0.0f;
    blur.plot(3.0f, static_cast<float>(row), Pixel(1, 1, 1), 0.0f, 1.0f,
              [&](float, float, const Pixel &, float, float alpha) {
                energy += alpha;
              });
    HS_EXPECT_NEAR(energy, pole ? 1.0f : 0.75f, 1e-6f);
  }
  hs_test::reset_globals();
  MeshFeedback<W, H> feedback;
  feedback.init();
  int lit = 0;
  int lit_endpoint_rows[2] = {0, 0};
  for (int frame = 0; frame < FEEDBACK_FRAMES; ++frame) {
    feedback.draw_frame();
    feedback.advance_display();
  }
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x) {
      const Pixel &p = feedback.get_pixel(x, y);
      if (!(p.r | p.g | p.b))
        continue;
      ++lit;
      if (y == 0)
        ++lit_endpoint_rows[0];
      if (y == H - 1)
        ++lit_endpoint_rows[1];
    }
  HS_EXPECT_GT(lit, 0);
  if (!Field::HAS_NORTH_POLE)
    HS_EXPECT_GT(lit_endpoint_rows[0], 0);
  if (!Field::HAS_SOUTH_POLE)
    HS_EXPECT_GT(lit_endpoint_rows[1], 0);
}
} // namespace

int main() {
  hs_test::ModuleFixture fixture("runtime display geometry");
  check_geometry<96, 48>(0.0f, math::PI_F);
  for (const auto endpoints :
       {std::pair{0.02f, 0.98f}, std::pair{0.0f, 0.8f}, std::pair{0.2f, 1.0f},
        std::pair{0.25f, 0.75f}, std::pair{0.0f, 1.0f}}) {
    const float north = endpoints.first * math::PI_F;
    const float south = endpoints.second * math::PI_F;
    HS_EXPECT_TRUE(math::set_display_geometry(north, south));
    check_geometry<96, 48>(north, south);
    check_geometry<288, 144>(north, south);
  }
  constexpr float NORTH = 0.2f;
  constexpr float SOUTH = 2.8f;
  HS_EXPECT_TRUE(math::set_display_geometry(NORTH, SOUTH));
  for (const float invalid :
       {-0.01f, std::nextafter(math::PI_F * 0.25f, math::PI_F), SOUTH,
        math::PI_F, std::numeric_limits<float>::quiet_NaN(),
        std::numeric_limits<float>::infinity()}) {
    HS_EXPECT_FALSE(math::set_display_geometry(invalid, SOUTH));
    HS_EXPECT_EQ(math::DISPLAY_NORTH_PHI, NORTH);
    HS_EXPECT_EQ(math::DISPLAY_SOUTH_PHI, SOUTH);
  }
  for (const float invalid :
       {-0.01f, NORTH, std::nextafter(math::PI_F * 0.75f, 0.0f),
        math::PI_F + 0.01f, std::numeric_limits<float>::quiet_NaN(),
        std::numeric_limits<float>::infinity()}) {
    HS_EXPECT_FALSE(math::set_display_geometry(NORTH, invalid));
    HS_EXPECT_EQ(math::DISPLAY_NORTH_PHI, NORTH);
    HS_EXPECT_EQ(math::DISPLAY_SOUTH_PHI, SOUTH);
  }
  return fixture.result() ? 1 : 0;
}
