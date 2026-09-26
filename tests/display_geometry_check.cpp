/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#include "core/engine/engine.h"
#include "core/math/display_geometry.h"
#include "core/render/filter/pixel_feedback.h"
#include "core/render/filter/splat.h"
#include "core/render/plot.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace {
constexpr int W = 96;
constexpr int H = 144;
using Geometry = math::DisplayGeometry<H>;

void test_coordinates() {
  static_assert(!Geometry::HAS_NORTH_POLE && !Geometry::HAS_SOUTH_POLE);
  HS_EXPECT_GT(Geometry::NORTH_PHI, 0.0f);
  HS_EXPECT_LT(Geometry::SOUTH_PHI, math::PI_F);
  HS_EXPECT_NEAR(Geometry::NORTH_PHI + Geometry::SOUTH_PHI, math::PI_F, 1e-6f);
  for (float y : {0.0f, 0.25f, 31.75f, 71.5f, 142.75f, 143.0f}) {
    HS_EXPECT_NEAR(Geometry::phi_to_row(Geometry::row_to_phi(y)), y, 2e-5f);
    HS_EXPECT_NEAR(math::y_to_phi<H>(y), Geometry::row_to_phi(y), 1e-6f);
  }
  for (int y : {0, H - 1}) {
    const auto a = math::pixel_to_vector<W, H>(0, y);
    const auto b = math::pixel_to_vector<W, H>(W / 2, y);
    HS_EXPECT_GT(std::abs(a.x - b.x), 0.01f);
    HS_EXPECT_NEAR(a.y, b.y, 1e-6f);
    HS_EXPECT_NEAR((math::TrigLUT<W, H>::cos_phi[y]),
                   std::cos(Geometry::row_to_phi(y)), 1e-6f);
  }
  const math::LatitudeGeometry asymmetric(H, 0.1f, math::PI_F - 0.2f);
  HS_EXPECT_NEAR(asymmetric.row_to_phi(0), 0.1f, 1e-6f);
  HS_EXPECT_NEAR(asymmetric.row_to_phi(H - 1), math::PI_F - 0.2f, 1e-6f);
  HS_EXPECT_LT(asymmetric.phi_to_row(0), 0.0f);
  HS_EXPECT_GT(asymmetric.phi_to_row(math::PI_F), H - 1.0f);
  for (float y : {0.0f, 0.25f, 72.5f, 143.0f})
    HS_EXPECT_NEAR(asymmetric.phi_to_row(asymmetric.row_to_phi(y)), y, 2e-5f);
}

void test_pole_crossing() {
  for (float y : {-0.25f, H - 0.75f}) {
    float x = 7.0f;
    HS_EXPECT_TRUE(!(math::pole_wrap<W, H>(x, y)));
    HS_EXPECT_NEAR(x, 7.0f, 1e-6f);
  }
  for (float pole : {Geometry::NORTH_POLE_ROW, Geometry::SOUTH_POLE_ROW}) {
    const float destination = pole < 0 ? 2.25f : H - 3.25f;
    float y = 2.0f * pole - destination;
    float x = 7.0f;
    HS_EXPECT_TRUE((math::pole_wrap<W, H>(x, y)));
    HS_EXPECT_NEAR(y, destination, 3e-5f);
    HS_EXPECT_NEAR(x, 7.0f + W / 2, 1e-6f);
  }
}

void test_splat_coverage() {
  Filter::Screen::AntiAlias<W, H> aa;
  for (float y : {-1.25f, -0.75f, -0.5f, -0.25f, 0.0f, H - 1.0f, H - 0.75f,
                  H - 0.5f, H - 0.25f, H + 0.25f}) {
    float energy = 0;
    aa.plot(12.25f, y, Pixel(100, 100, 100), 0.0f, 1.0f,
            [&](float, float row, const Pixel &, float, float alpha) {
              HS_EXPECT_GE(row, 0.0f);
              HS_EXPECT_LT(row, static_cast<float>(H));
              energy += alpha;
            });
    float expected = 0;
    if (y > -1.0f && y < 0.0f)
      expected = math::quintic_kernel(y + 1.0f);
    else if (y >= 0.0f && y <= H - 1.0f)
      expected = 1;
    else if (y > H - 1.0f && y < H)
      expected = 1.0f - math::quintic_kernel(y - (H - 1));
    HS_EXPECT_NEAR(energy, expected, 1e-6f);
  }
}

void test_blur_coverage() {
  for (float factor : {0.0f, 0.5f, 1.0f}) {
    Filter::Screen::Blur<W, H> blur(factor);
    for (float y : {-2.0f, -1.25f, -1.0f, -0.75f, -0.25f, 0.0f, 0.75f, 2.0f}) {
      const int center = static_cast<int>(std::round(y));
      const float expected = center < -1    ? 0.0f
                             : center == -1 ? factor * 0.25f
                             : center == 0  ? 1.0f - factor * 0.25f
                                            : 1.0f;
      for (float row : {y, H - 1.0f - y}) {
        float energy = 0;
        blur.plot(12.25f, row, Pixel(100, 100, 100), 0.0f, 1.0f,
                  [&](float, float tap_row, const Pixel &, float, float alpha) {
                    HS_EXPECT_GE(tap_row, 0.0f);
                    HS_EXPECT_LT(tap_row, static_cast<float>(H));
                    energy += alpha;
                  });
        HS_EXPECT_NEAR(energy, expected, 1e-6f);
      }
    }
  }
}

int plot_stroke(float phi) {
  hs_test::StubEffect effect(W, H);
  Pipeline<W, H> pipeline;
  Fragment first, last;
  first.pos = {std::sin(phi), std::cos(phi), 0.0f};
  last.pos = {std::sin(phi) * std::cos(0.4f), std::cos(phi),
              std::sin(phi) * std::sin(0.4f)};
  {
    Canvas canvas(effect);
    Plot::Line::draw<W, H>(pipeline, canvas, first, last,
                           [](const math::Vector &, Fragment &fragment) {
                             fragment.color =
                                 Color4(Pixel(60000, 60000, 60000), 1.0f);
                           });
  }
  effect.advance_display();
  int lit = 0;
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x) {
      const Pixel &pixel = effect.get_pixel(x, y);
      if ((pixel.r | pixel.g | pixel.b) != 0)
        ++lit;
    }
  return lit;
}

void test_render_caps() {
  HS_EXPECT_EQ(plot_stroke(Geometry::NORTH_PHI * 0.5f), 0);
  HS_EXPECT_EQ(plot_stroke((Geometry::SOUTH_PHI + math::PI_F) * 0.5f), 0);
  HS_EXPECT_GT(plot_stroke(Geometry::row_to_phi(2.0f)), 0);
  HS_EXPECT_GT(plot_stroke(Geometry::row_to_phi(H - 3.0f)), 0);
  const auto north = SDF::phi_bounds_to_rows<H>(0, Geometry::NORTH_PHI * 0.5f);
  const auto south = SDF::phi_bounds_to_rows<H>(
      (Geometry::SOUTH_PHI + math::PI_F) * 0.5f, math::PI_F);
  HS_EXPECT_GT(north.y_min, north.y_max);
  HS_EXPECT_GT(south.y_min, south.y_max);
  const auto band = SDF::phi_bounds_to_rows<H>(Geometry::row_to_phi(10.25f),
                                               Geometry::row_to_phi(12.75f));
  HS_EXPECT_EQ(band.y_min, 10);
  HS_EXPECT_EQ(band.y_max, 13);
}

math::Vector rotate_longitude(const math::Vector &v,
                              const ::Feedback::Style &) {
  constexpr float ANGLE = 0.6f;
  return {std::cos(ANGLE) * v.x - std::sin(ANGLE) * v.z, v.y,
          std::sin(ANGLE) * v.x + std::cos(ANGLE) * v.z};
}

void test_feedback_endpoint_rings() {
  static_assert(hs::SphericalFieldLayout<W, H>::POLE_COUNT == 0);
  constexpr Pixel BRIGHT(12000, 30000, 50000);
  constexpr int SOURCE_X = W / 3;
  for (int row : {0, H - 1}) {
    hs_test::StubEffect effect(W, H);
    ::Feedback::Style style{};
    style.space_fn = &rotate_longitude;
    style.noise = nullptr;
    style.fade = 1.0f;
    style.downsample = 4;
    Pipeline<W, H, Filter::Pixel::Feedback<W, H>> pipeline{
        Filter::Pixel::Feedback<W, H>(style)};
    {
      Canvas canvas(effect);
      canvas(SOURCE_X, row) = BRIGHT;
    }
    effect.advance_display();
    {
      Canvas canvas(effect);
      (void)pipeline.begin_frame(canvas, 1.0f);
    }
    effect.advance_display();
    double mass = 0, moment = 0;
    int lit = 0;
    for (int x = 0; x < W; ++x) {
      const double brightness = effect.get_pixel(x, row).b;
      if (brightness > 0)
        ++lit;
      double dx = x - SOURCE_X;
      if (dx > W * 0.5)
        dx -= W;
      else if (dx < -W * 0.5)
        dx += W;
      mass += brightness;
      moment += brightness * dx;
    }
    HS_EXPECT_GT(mass, BRIGHT.b * 0.5);
    HS_EXPECT_GT(lit, 0);
    HS_EXPECT_LT(lit, W);
    HS_EXPECT_NEAR(moment / mass, -0.6 * W / (2 * math::PI_F), 0.5);
  }
}
} // namespace

int main() {
  hs_test::ModuleFixture fixture("physical display geometry");
  test_coordinates();
  test_pole_crossing();
  test_splat_coverage();
  test_blur_coverage();
  test_feedback_endpoint_rings();
  test_render_caps();
  return fixture.result() ? 1 : 0;
}
