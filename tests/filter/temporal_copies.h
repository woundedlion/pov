/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @file temporal_copies.h
 * @brief Selected-copy timed transforms match spatial fan-out. */

/** @brief Checks a selected copy against the corresponding ordinary emission.
 * @param stage Spatial replicator under test.
 * @param reference Independent double-precision copy transform. */
template <typename Stage, typename Reference>
void check_timed_spatial_copies(Stage &stage, Reference reference) {
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
  for (int copy = 0; copy < stage.copy_count(); ++copy) {
    HS_CONTEXT("copy", copy);
    auto check_reference = [&](const math::Vector &got,
                               const math::Vector &input) {
      const auto independent = reference(input, copy);
      HS_EXPECT_NEAR(got.x, independent[0], 2e-6);
      HS_EXPECT_NEAR(got.y, independent[1], 2e-6);
      HS_EXPECT_NEAR(got.z, independent[2], 2e-6);
      if (copy == 0) {
        HS_EXPECT_EQ(got.x, input.x);
        HS_EXPECT_EQ(got.y, input.y);
        HS_EXPECT_EQ(got.z, input.z);
      }
    };
    check_reference(expected[copy], source);
    stage.select_copy(copy);
    for (float t : {0.0f, 0.25f, 1.0f}) {
      const auto selected = stage.at_time(canvas, source, t);
      check_reference(selected, source);
      int emitted = 0;
      stage.plot_at(canvas, source, color, 7.0f, 0.4f, t,
                    [&](const math::Vector &v, const Pixel &pixel, float age,
                        float alpha) {
                      check_reference(v, source);
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
            check_reference(a, source);
            check_reference(b, source);
            HS_EXPECT_EQ(transformed != nullptr, planar);
            if (transformed) {
              check_reference(transformed->u, basis.u);
              check_reference(transformed->v, basis.v);
              check_reference(transformed->w, basis.w);
            }
            ++culled;
            return false;
          });
      HS_EXPECT_FALSE(visible);
      HS_EXPECT_EQ(culled, 1);
    }
  }
}

/** @brief Timed Y-axis copies select the analytic spatial rotation. */
inline void test_timed_replicate_selected_copy() {
  Filter::World::Replicate<288> stage(288);
  auto reference = [&](const math::Vector &v, int copy) {
    const double ANGLE = 2.0 * double(math::PI_F) * copy / stage.copy_count();
    return std::array<double, 3>{v.x * std::cos(ANGLE) + v.z * std::sin(ANGLE),
                                 double(v.y),
                                 v.z * std::cos(ANGLE) - v.x * std::sin(ANGLE)};
  };
  check_timed_spatial_copies(stage, reference);
  stage.set_count(3);
  check_timed_spatial_copies(stage, reference);
}

/** @brief Timed vertex copies use the selected stored rotation. */
inline void test_timed_vertex_replicate_selected_copy() {
  const std::array<math::Vector, 4> vertices{math::X_AXIS, math::Y_AXIS,
                                             math::Z_AXIS, -math::X_AXIS};
  Filter::World::VertexReplicate<4> stage(vertices);
  check_timed_spatial_copies(stage, [](const math::Vector &v, int copy) {
    const std::array<std::array<double, 3>, 4> copies{{{v.x, v.y, v.z},
                                                       {-v.y, v.x, v.z},
                                                       {-v.z, v.y, v.x},
                                                       {-v.x, -v.y, v.z}}};
    return copies[copy];
  });
}
