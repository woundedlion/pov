/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Pole geometry rendering test helpers.
 */
#pragma once

#include <cmath>

#include "core/render/plot.h"
#include "core/render/filter.h"
#include "core/render/filter/pixel_feedback.h"
#include "effects/MeshFeedback.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"
#include "tests/pixel_test_util.h"

namespace hs_test::pole_geometry {
/** Longitude warp per feedback frame, in radians. */
inline constexpr double LONGITUDE_STEP = 0.6;

/** @brief Checks MeshFeedback's arena peaks against the device partitions.
 * @tparam W,H Canvas dimensions.
 * @tparam Inspect Post-smoke effect inspector.
 * @param inspect Receives the effect after the smoke frames. */
template <int W, int H,
          typename Inspect = decltype([](const MeshFeedback<W, H> &) {})>
inline void check_feedback_device_budget(Inspect inspect = {}) {
  reset_globals();
  MeshFeedback<W, H> effect;
  effect.init();
  for (int frame = 0; frame < smoke_frames(); ++frame) {
    pin_frame_clock(frame);
    effect.draw_frame();
    effect.advance_display();
  }
  inspect(effect);
  std::printf("  [MeshFeedback %dx%d] arenas %zu/%zu/%zu\n", W, H,
              persistent_arena.get_high_water_mark(),
              scratch_arena_a.get_high_water_mark(),
              scratch_arena_b.get_high_water_mark());
  constexpr ArenaSplit SPLIT = MeshFeedback<W, H>::arena_split();
  HS_EXPECT_LE(persistent_arena.get_high_water_mark(),
               SPLIT.device_persistent());
  HS_EXPECT_LE(scratch_arena_a.get_high_water_mark(), SPLIT.scratch_a);
  HS_EXPECT_LE(scratch_arena_b.get_high_water_mark(), SPLIT.scratch_b);
}

/**
 * @brief Count lit pixels after a geodesic between points at colatitude phi.
 * @details Endpoints are 0.4 rad apart in azimuth; the arc bulges poleward.
 */
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
  return static_cast<int>(count_lit_region<W, H>(effect));
}

/** @brief Applies the test longitude warp as a rotation about the polar axis. */
inline math::Vector rotate_longitude(const math::Vector &v,
                                     const ::Feedback::Style &) {
  constexpr float ANGLE = static_cast<float>(LONGITUDE_STEP);
  return {std::cos(ANGLE) * v.x - std::sin(ANGLE) * v.z, v.y,
          std::sin(ANGLE) * v.x + std::cos(ANGLE) * v.z};
}

/** @brief Check the longitude centroid of a warped endpoint ring. */
template <int W, int H>
inline void check_feedback_ring_centroid(int row, int lit_threshold = 0,
                                         int downsample = 4,
                                         float alpha = 1.0f) {
  HS_CONTEXT("resolution", W, H);
  HS_CONTEXT("row / alpha-milli", row, static_cast<int>(alpha * 1000));
  for (bool temporal : {false, true}) {
    HS_CONTEXT("downsample / temporal", downsample, static_cast<int>(temporal));
    const ScopedTemporalSampling sampling(temporal);
    const size_t SCRATCH_B = scratch_arena_b.get_capacity();
    const size_t REQUIRED =
        (Filter::Pixel::Feedback<W, H>::uncached_scratch_bytes(downsample) +
         Filter::Pixel::FeedbackFadeTable::BYTES + 127) &
        ~size_t{63};
    const ScopedArenaSplit split(GLOBAL_ARENA_SIZE - REQUIRED - SCRATCH_B,
                                 REQUIRED, SCRATCH_B);
    constexpr Pixel BRIGHT(12000, 30000, 50000);
    constexpr int SOURCE_X = W / 3;
    hs_test::StubEffect effect(W, H);
    ::Feedback::Style style{};
    style.space_fn = &rotate_longitude;
    style.noise = nullptr;
    style.fade = 1.0f;
    style.downsample = downsample;
    Pipeline<W, H, Filter::Pixel::Feedback<W, H>> pipeline{
        Filter::Pixel::Feedback<W, H>(style)};
    {
      Canvas canvas(effect);
      canvas(SOURCE_X, row) = BRIGHT;
    }
    effect.advance_display();
    {
      Canvas canvas(effect);
      (void)pipeline.begin_frame(canvas, alpha);
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
    HS_EXPECT_GT(mass, BRIGHT.b * alpha * 0.5);
    HS_EXPECT_GT(lit, 0);
    HS_EXPECT_LT(lit, W);
    constexpr double STEP = LONGITUDE_STEP * W / (2 * math::PI_F);
    // A within-arm read solves dx = -STEP * (1 + 2 * dx / W).
    const double EXPECTED = temporal ? -STEP / (1.0 + 2.0 * STEP / W) : -STEP;
    HS_EXPECT_NEAR(moment / mass, EXPECTED, 0.5);
  }
}
} // namespace hs_test::pole_geometry
