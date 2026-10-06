/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Path::get_point
// ============================================================================

/** @brief Verifies adjacent segments share their boundary sample and interpolate linearly. */
inline void test_path_adjacent_segments_fill_exact_capacity() {
  Path<5> path;
  path.append_segment([](float t) { return math::Vector(t, 0, 0); }, 2.0f, 2,
                      math::ease_linear);
  path.append_segment([](float t) { return math::Vector(2.0f + t, 0, 0); },
                      2.0f, 2, math::ease_linear);
  for (int i = 0; i < 5; ++i)
    HS_EXPECT_NEAR(path.get_point(i * 0.25f).x, static_cast<float>(i), 1e-6f);
}

/**
 * @brief Verifies an empty Path returns the origin for any t (the no-points
 * guard).
 */
inline void test_path_empty_returns_origin() {
  Path<32> p;
  math::Vector v = p.get_point(0.0f);
  HS_EXPECT_NEAR(v.x, 0.0f, 0.0f);
  HS_EXPECT_NEAR(v.y, 0.0f, 0.0f);
  HS_EXPECT_NEAR(v.z, 0.0f, 0.0f);
  math::Vector v1 = p.get_point(1.0f);
  HS_EXPECT_NEAR(v1.x, 0.0f, 0.0f);
}

/**
 * @brief Verifies get_point hits the exact endpoints at t=0/1, clamps to back()
 * past t=1 and to front() below t=0, and interpolates linearly between samples.
 */
inline void test_path_endpoints_and_clamp() {
  Path<32> p;
  p.append_segment(
      [](float s) {
        return math::Vector(1.0f + s, 2.0f + 2.0f * s, 3.0f * s - 2.0f);
      },
      1.0f, 4, math::ease_linear);

  math::Vector start = p.get_point(0.0f);
  HS_EXPECT_NEAR(start.x, 1.0f, 1e-5f);
  HS_EXPECT_NEAR(start.y, 2.0f, 1e-5f);
  HS_EXPECT_NEAR(start.z, -2.0f, 1e-5f);

  math::Vector end = p.get_point(1.0f);
  HS_EXPECT_NEAR(end.x, 2.0f, 1e-5f);
  HS_EXPECT_NEAR(end.y, 4.0f, 1e-5f);
  HS_EXPECT_NEAR(end.z, 1.0f, 1e-5f);

  math::Vector over = p.get_point(2.0f);
  HS_EXPECT_NEAR(over.x, end.x, 1e-5f);
  HS_EXPECT_NEAR(over.y, end.y, 1e-5f);
  HS_EXPECT_NEAR(over.z, end.z, 1e-5f);

  // Only the clamp keeps this out of a negative float->size_t cast.
  math::Vector under = p.get_point(-1.0f);
  HS_EXPECT_NEAR(under.x, start.x, 1e-5f);
  HS_EXPECT_NEAR(under.y, start.y, 1e-5f);
  HS_EXPECT_NEAR(under.z, start.z, 1e-5f);

  math::Vector mid = p.get_point(0.5f);
  HS_EXPECT_NEAR(mid.x, 1.5f, 1e-5f);
  HS_EXPECT_NEAR(mid.y, 3.0f, 1e-5f);
  HS_EXPECT_NEAR(mid.z, -0.5f, 1e-5f);
}

/**
 * @brief Verifies collapse() reduces the path to its last point; get_point then
 * returns it for any t.
 */
inline void test_path_collapse_keeps_last() {
  Path<32> p;
  p.append_segment(
      [](float s) {
        return math::Vector(1.0f + s, 2.0f + 2.0f * s, 3.0f * s - 2.0f);
      },
      1.0f, 4, math::ease_linear);
  math::Vector last_before = p.get_point(1.0f);
  p.collapse();
  math::Vector a = p.get_point(0.0f);
  math::Vector b = p.get_point(1.0f);
  HS_EXPECT_NEAR(a.x, last_before.x, 1e-5f);
  HS_EXPECT_NEAR(a.y, last_before.y, 1e-5f);
  HS_EXPECT_NEAR(a.z, last_before.z, 1e-5f);
  HS_EXPECT_NEAR(b.x, last_before.x, 1e-5f);
  HS_EXPECT_NEAR(b.y, last_before.y, 1e-5f);
  HS_EXPECT_NEAR(b.z, last_before.z, 1e-5f);
}
