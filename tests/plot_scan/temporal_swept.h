/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @brief Swept heights and half-planes from frame-end captures bound every
 *         point of a rotating edge's arcs over the frame. */
inline void test_swept_bounds_cover_moving_edges() {
  int bulging = 0;
  for (const float lat : {-0.5f, 0.1f, 0.5f, 0.8f}) {
    for (const float span : {0.2f, 1.0f}) {
      for (const float angle : {0.02f, 0.25f}) {
        HS_CONTEXT("lat, span x100", std::llround(lat * 100.0f),
                   std::llround(span * 100.0f));
        HS_CONTEXT("angle x100", std::llround(angle * 100.0f));
        const float r = sqrtf(1.0f - lat * lat);
        const math::Vector a0(r * cosf(-0.5f * span), lat,
                              r * sinf(-0.5f * span));
        const math::Vector b0(r * cosf(0.5f * span), lat,
                              r * sinf(0.5f * span));
        const math::Vector axis = math::Vector(0.6f, 0.2f, 0.77f).normalized();
        const auto at = [&](const math::Vector &p, float t) {
          return math::rotate(p, math::make_rotation(axis, angle * t));
        };
        const math::Vector captures[] = {a0, at(a0, 1.0f), b0, at(b0, 1.0f)};
        const float pad =
            Spacetime::PathBound::rotation(angle).chord_error() + 1e-6f;
        float low = 0.0f, high = 0.0f;
        HS_EXPECT_TRUE(Spacetime::swept_heights(captures, 4, pad, low, high));
        const math::Vector normals[] = {
            math::Vector(1.0f, 0.0f, 0.0f), math::Vector(-0.3f, 0.0f, 1.0f),
            math::Vector(0.0f, 1.0f, 0.0f), math::Vector(0.2f, -1.0f, 0.1f)};
        bool beyond[4];
        for (int n = 0; n < 4; ++n)
          beyond[n] = Spacetime::swept_beyond(captures, 4, pad, normals[n]);
        float top = -1.0f;
        for (int i = 0; i <= 64; ++i) {
          const float t = static_cast<float>(i) / 64.0f;
          const math::Vector a = at(a0, t), b = at(b0, t);
          for (int j = 0; j <= 64; ++j) {
            const float s = static_cast<float>(j) / 64.0f;
            const math::Vector p = (a * (1.0f - s) + b * s).normalized();
            top = fmaxf(top, p.y);
            HS_EXPECT_GE(p.y, low);
            HS_EXPECT_LE(p.y, high);
            for (int n = 0; n < 4; ++n)
              if (beyond[n])
                HS_EXPECT_GT(math::dot(p, normals[n]), 0.0f);
          }
        }
        float highest_capture = -1.0f;
        for (const auto &c : captures)
          highest_capture = fmaxf(highest_capture, c.y);
        if (top > highest_capture + pad)
          ++bulging;
      }
    }
  }
  HS_EXPECT_GT(bulging, 0);
}

/** @brief Captures spread past a quarter turn leave the heights unbounded. */
inline void test_swept_heights_refuse_wide_spreads() {
  const math::Vector captures[] = {
      math::X_AXIS, math::Vector(-0.1f, 0.0f, 1.0f).normalized()};
  float low = 0.0f, high = 0.0f;
  HS_EXPECT_FALSE(Spacetime::swept_heights(captures, 2, 0.0f, low, high));
  HS_EXPECT_EQ(low, -1.0f);
  HS_EXPECT_EQ(high, 1.0f);
}
