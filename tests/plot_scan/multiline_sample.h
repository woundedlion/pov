/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Plot::Multiline::sample
// ============================================================================

/**
 * @brief Verifies Multiline::sample arc-length-parameterizes the polyline: v0
 *        progress is monotone 0..1, v2 carries the vertex index, and v1
 *        cumulative arc length sums the per-edge angles.
 */
inline void test_multiline_sample_arclength_param() {
  ScratchScope sc(plot_arena());

  Fragments verts;
  verts.bind(plot_arena(), 4);
  Fragment v;
  v.pos = math::Vector(1, 0, 0);
  verts.push_back(v);
  v.pos = math::Vector(0, 1, 0);
  verts.push_back(v);
  v.pos = math::Vector(-1, 0, 0);
  verts.push_back(v);

  Fragments points;
  points.bind(plot_arena(), 8);
  Plot::Multiline::sample(points, verts, /*closed=*/false);

  HS_EXPECT_SIZE_OR_RETURN(points, (size_t)3);

  HS_EXPECT_NEAR(points[0].v0, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(points.back().v0, 1.0f, 1e-4f);
  float last = -1.0f;
  for (size_t i = 0; i < points.size(); ++i) {
    HS_EXPECT_GE(points[i].v0, last - 1e-6f);
    last = points[i].v0;
    HS_EXPECT_NEAR(points[i].v2, static_cast<float>(i), 1e-6f);
  }
  // v1 cumulative arc length sums the two equal 90deg hops.
  HS_EXPECT_NEAR(points.back().v1, math::PI_F, 1e-3f);

  Fragments closed_points;
  closed_points.bind(plot_arena(), 8);
  Fragment seam =
      Plot::Multiline::sample(closed_points, verts, /*closed=*/true);
  HS_EXPECT_SIZE_OR_RETURN(closed_points, verts.size());
  HS_EXPECT_NEAR(seam.pos.x, closed_points[0].pos.x, 1e-6f);
  HS_EXPECT_NEAR(seam.pos.y, closed_points[0].pos.y, 1e-6f);
  HS_EXPECT_NEAR(seam.v0, 1.0f, 1e-6f);
  HS_EXPECT_NEAR(seam.v1, 2.0f * math::PI_F, 1e-3f);
  HS_EXPECT_NEAR(seam.v2, 3.0f, 1e-6f);
}
