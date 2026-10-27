/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Lens numeric cases. */

/**
 * @brief Verifies glitch_lens maps unit directions to unit directions and
 *        stays smooth through its equatorial fold.
 * @details Pins |lens(v)| == 1 across a spread of directions, its doubled-
 *          latitude topology, and matching one-sided derivatives at the fold.
 */
inline void test_glitch_lens_unit_norm() {
  for (float polar : {0.4f, 0.9f, 2.0f}) {
    for (float azimuth : {-2.1f, -0.4f, 0.7f, 1.8f}) {
      const math::Vector input(sinf(polar) * cosf(azimuth), cosf(polar),
                               sinf(polar) * sinf(azimuth));
      const math::Vector mapped = lenses::glitch_lens(input);
      HS_EXPECT_VEC(mapped,
                    math::Vector(sinf(2 * polar) * cosf(3 * azimuth),
                                 cosf(2 * polar),
                                 sinf(2 * polar) * sinf(3 * azimuth)),
                    2e-6f);
    }
  }
  const math::Vector dirs[] = {math::Vector(1, 0, 0),
                               math::Vector(0, 0, 1),
                               math::Vector(-1, 0, 0),
                               math::Vector(0, 0, -1),
                               math::Vector(1, 1, 1).normalized(),
                               math::Vector(-1, 2, -3).normalized(),
                               math::Vector(3, -1, 2).normalized(),
                               math::Vector(0.2f, -0.9f, 0.4f).normalized(),
                               math::Vector(-0.7f, 0.1f, 0.7f).normalized()};
  for (const math::Vector &v : dirs) {
    HS_EXPECT_NEAR(lenses::glitch_lens(v).length(), 1.0f, 1e-3f);
  }

  const math::Vector equator_x = lenses::glitch_lens(math::Vector(1, 0, 0));
  HS_EXPECT_VEC(equator_x, math::Vector(0, -1, 0), 1e-6f);

  const math::Vector north = lenses::glitch_lens(math::Vector(0, 1, 0));
  HS_EXPECT_NEAR(north.y, 1.0f, 1e-6f);
  const math::Vector south = lenses::glitch_lens(math::Vector(0, -1, 0));
  HS_EXPECT_NEAR(south.y, 1.0f, 1e-6f);

  constexpr float EPSILON = 1e-4f;
  const float radial = sqrtf(1.0f - EPSILON * EPSILON);
  const math::Vector center =
      lenses::glitch_lens(math::Vector(0.8f, 0.0f, 0.6f));
  const math::Vector above =
      lenses::glitch_lens(math::Vector(0.8f * radial, EPSILON, 0.6f * radial));
  const math::Vector below =
      lenses::glitch_lens(math::Vector(0.8f * radial, -EPSILON, 0.6f * radial));
  const math::Vector forward = above - center;
  const math::Vector backward = center - below;
  HS_EXPECT_VEC(forward, backward, 1e-6f);
}
