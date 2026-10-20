/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @file temporal_copies.h
 * @brief Selected-copy timed transforms match spatial fan-out. */

/** @brief Checks a selected copy against the corresponding ordinary emission.
 * @param stage Spatial replicator under test. */
template <typename Stage> void check_timed_spatial_copies(Stage &stage) {
  StubEffect effect(288, 144);
  Canvas canvas(effect);
  const math::Vector source = math::Vector(0.3f, 0.5f, 0.7f).normalized();
  const math::Basis basis{math::X_AXIS, math::Y_AXIS, math::Z_AXIS};
  const Pixel color(73, 91, 129);
  std::vector<math::Vector> expected;
  stage.plot(source, color, 7.0f, 0.4f,
             [&](const math::Vector &v, const Pixel &, float, float) {
               expected.push_back(v);
             });
  HS_EXPECT_EQ(expected.size(), static_cast<size_t>(stage.copy_count()));
  auto equal = [](const math::Vector &a, const math::Vector &b) {
    HS_EXPECT_EQ(a.x, b.x);
    HS_EXPECT_EQ(a.y, b.y);
    HS_EXPECT_EQ(a.z, b.z);
  };
  for (int copy = 0; copy < stage.copy_count(); ++copy) {
    stage.select_copy(copy);
    for (float t : {0.0f, 0.25f, 1.0f}) {
      equal(stage.at_time(canvas, source, t), expected[copy]);
      int emitted = 0;
      stage.plot_at(canvas, source, color, 7.0f, 0.4f, t,
                    [&](const math::Vector &v, const Pixel &pixel, float age,
                        float alpha) {
                      equal(v, expected[copy]);
                      HS_EXPECT_EQ(pixel.r, color.r);
                      HS_EXPECT_EQ(pixel.g, color.g);
                      HS_EXPECT_EQ(pixel.b, color.b);
                      HS_EXPECT_EQ(age, 7.0f);
                      HS_EXPECT_EQ(alpha, 0.4f);
                      ++emitted;
                    });
      HS_EXPECT_EQ(emitted, 1);
    }
    for (bool planar : {false, true}) {
      int culled = 0;
      const bool visible = stage.cull_edge_timed(
          source, source, planar ? &basis : nullptr,
          [&](const math::Vector &a, const math::Vector &b,
              const math::Basis *transformed) {
            equal(a, expected[copy]);
            equal(b, expected[copy]);
            HS_EXPECT_EQ(transformed != nullptr, planar);
            if (transformed) {
              equal(transformed->u, stage.at_time(canvas, basis.u, 0.4f));
              equal(transformed->v, stage.at_time(canvas, basis.v, 0.4f));
              equal(transformed->w, stage.at_time(canvas, basis.w, 0.4f));
            }
            ++culled;
            return false;
          });
      HS_EXPECT_FALSE(visible);
      HS_EXPECT_EQ(culled, 1);
    }
  }
}

/** @brief Timed Y-axis copies retain ordinary chained rotation arithmetic. */
inline void test_timed_replicate_selected_copy() {
  Filter::World::Replicate<288> stage(288);
  check_timed_spatial_copies(stage);
  stage.set_count(3);
  check_timed_spatial_copies(stage);
}

/** @brief Timed vertex copies use the selected stored rotation. */
inline void test_timed_vertex_replicate_selected_copy() {
  const std::array<math::Vector, 4> vertices{math::X_AXIS, math::Y_AXIS,
                                             math::Z_AXIS, -math::X_AXIS};
  Filter::World::VertexReplicate<4> stage(vertices);
  check_timed_spatial_copies(stage);
}
