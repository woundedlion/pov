/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <cmath>

#include "core/render/plot.h"
#include "core/render/filter.h"
#include "core/render/filter/pixel_feedback.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test::pole_geometry {
/** @brief Count lit pixels after a short constant-colatitude stroke. */
template <int W, int H> inline int plot_stroke(float phi) {
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

inline math::Vector rotate_longitude(const math::Vector &v,
                                     const ::Feedback::Style &) {
  constexpr float ANGLE = 0.6f;
  return {std::cos(ANGLE) * v.x - std::sin(ANGLE) * v.z, v.y,
          std::sin(ANGLE) * v.x + std::cos(ANGLE) * v.z};
}

/** @brief Check the longitude centroid of a warped endpoint ring. */
template <int W, int H>
inline void check_feedback_ring_centroid(int row, int lit_threshold = 0) {
  constexpr Pixel BRIGHT(12000, 30000, 50000);
  constexpr int SOURCE_X = W / 3;
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
    if (brightness > lit_threshold)
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
} // namespace hs_test::pole_geometry
