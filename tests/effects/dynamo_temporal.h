/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @file dynamo_temporal.h
 * @brief Dynamo column-time strands against independently frozen copies. */

/** @brief Independent world-space reference for Dynamo's live and replayed strands. */
struct DynamoTemporalOracle {
  /** @brief Opaque white palette making fresh samples visible. */
  struct WhitePalette {
    /** @brief Constant source color.
     * @return Opaque white. */
    Color4 get(float) const { return Color4(Pixel(65535, 65535, 65535), 1.0f); }
  };

  /** @brief Live source history records each edge endpoint before fan-out. */
  static void check_live_history_fanout() {
    constexpr int W = 288, H = 144;
    reset_effect_globals();
    Dynamo<W, H> effect;
    effect.init();
    effect.baked_palettes[0].rebake(WhitePalette{});
    effect.orientation.set(math::Quaternion());
    effect.orientation.push(math::make_rotation(math::Y_AXIS, 0.018f));
    for (size_t i = 0; i < effect.NUM_NODES; ++i) {
      effect.nodes[i].x = 20;
      effect.nodes[i].y =
          H * (0.3f + 0.4f * float(i) / float(effect.NUM_NODES - 1));
    }
    Spacetime::geometry_counters = {};
    {
      Canvas canvas(effect);
      effect.draw_nodes_timed(canvas, 0.0f);
    }
    const auto &trails = effect.filters.template get<
        Filter::World::Trails<Dynamo<W, H>::TRAIL_CAPACITY>>();
    const uint32_t edge_samples = uint32_t(2 * (effect.NUM_NODES - 1));
    HS_EXPECT_EQ(trails.records(), edge_samples + 1);
    HS_EXPECT_EQ(Spacetime::geometry_counters.column_samples, 6 * edge_samples);
    effect.advance_display();
  }

  /** @brief Rotates a point with double-precision Rodrigues arithmetic.
   * @param v Source point.
   * @param axis Unit rotation axis.
   * @param angle Rotation angle, radians.
   * @return Rotated point. */
  static math::Vector turn(const math::Vector &v, const math::Vector &axis,
                           double angle) {
    const double sine = std::sin(angle), cosine = std::cos(angle);
    const double projection =
        double(v.x) * axis.x + double(v.y) * axis.y + double(v.z) * axis.z;
    const auto component = [&](double value, double cross, double a) {
      return float(value * cosine + cross * sine +
                   a * projection * (1.0 - cosine));
    };
    return {
        component(v.x, double(axis.y) * v.z - double(axis.z) * v.y, axis.x),
        component(v.y, double(axis.z) * v.x - double(axis.x) * v.z, axis.y),
        component(v.z, double(axis.x) * v.y - double(axis.y) * v.x, axis.z)};
  }

  /** @brief Renders the timed draw or independently placed frozen geometry.
   * @tparam W,H Display dimensions.
   * @param timed Selects production column-time rendering.
   * @param fraction Reference frame fraction.
   * @param displacement Reference displacement in physical pixels.
   * @param start Source strand's starting column.
   * @param spin Frame rotation, radians.
   * @param replay Selects stored samples instead of live edges.
   * @param fresh Holds the frame at its endpoint, as on initialization or reset.
   * @param clipped Forces a segment clip for the reference comparison.
   * @param reset Draws a moving frame before resetting orientation history.
   * @return Linear framebuffer luminance. */
  template <int W, int H>
  static std::vector<float> render(bool timed, float fraction,
                                   float displacement, float start, float spin,
                                   bool replay, bool fresh = false,
                                   bool clipped = false, bool reset = false) {
    reset_effect_globals();
    Dynamo<W, H> effect;
    effect.init();
    effect.baked_palettes[0].rebake(WhitePalette{});
    if (clipped) {
      const auto map = pov::segment_map(1, H * 2, 4);
      const auto clip = pov::segment_clip(map, true, H * 2, 4, W);
      effect.set_clip(clip.y0, clip.y1, clip.x0, clip.x1);
    }
    const auto axis = math::Vector(0.15f, 1.0f, 0.1f).normalized();
    const auto opening = math::make_rotation(math::X_AXIS, 0.27f);
    const auto finish = math::make_rotation(axis, spin) * opening;
    effect.orientation.set(fresh ? finish : opening);
    if (!fresh)
      effect.orientation.push(finish);
    auto &trails = effect.filters.template get<
        Filter::World::Trails<Dynamo<W, H>::TRAIL_CAPACITY>>();
    trails.set_lifetime(4);
    for (size_t i = 0; i < effect.NUM_NODES; ++i) {
      effect.nodes[i].x =
          math::wrap(int(start + 8.0f * sinf(float(i) * 0.13f)), W);
      effect.nodes[i].y =
          H * (0.25f + 0.5f * float(i) / float(effect.NUM_NODES - 1));
    }
    if (reset) {
      effect.orientation.set(opening);
      effect.orientation.push(math::make_rotation(axis, -0.2f) * opening);
      {
        Canvas canvas(effect);
        effect.draw_nodes_timed(canvas, 0.0f);
      }
      effect.advance_display();
      effect.orientation.set(finish);
    }
    {
      Canvas canvas(effect);
      ScratchScope storage_a(scratch_arena_a);
      ScratchScope storage_b(scratch_arena_b);
      auto *source = scratch_arena_a.allocate_n<math::Vector>(effect.NUM_NODES);
      for (size_t i = 0; i < effect.NUM_NODES; ++i) {
        source[i] = math::pixel_to_vector<W, H>(float(effect.nodes[i].x),
                                                effect.nodes[i].y);
        if (replay) {
          trails.plot(source[i], Pixel{}, 0.0f, 0.0f,
                      [](const math::Vector &, const Pixel &, float, float) {});
          source[i] = math::normalized_or(
              math::Snorm3::encode(source[i]).decode(), math::Y_AXIS);
        }
      }
      if (timed) {
        if (replay)
          effect.flush_timed(canvas);
        else
          effect.draw_nodes_timed(canvas, 0.0f);
      } else {
        auto *vertices =
            scratch_arena_a.allocate_n<math::Vector>(effect.NUM_NODES);
        ArenaVector<Plot::Mesh::Edge> edges(scratch_arena_a,
                                            effect.NUM_NODES - 1);
        for (size_t i = 1; i < effect.NUM_NODES; ++i)
          edges.push_back({uint16_t(i - 1), uint16_t(i)});
        auto &sink = effect.filters.next.next.next;
        const Color4 color = effect.color(source[0], 0.0f);
        const auto shader = [&](const math::Vector &, Fragment &f) {
          f.color = color;
          f.color.alpha *= 0.5f;
        };
        Plot::GeometryCallback geometry;
        geometry.rigid_lattice = true;
        for (int copy = 0; copy < effect.STRAND_COPIES; ++copy) {
          for (size_t i = 0; i < effect.NUM_NODES; ++i) {
            auto v = turn(source[i], math::Y_AXIS,
                          double(math::TWO_PI_F) * copy / effect.STRAND_COPIES);
            v = turn(v, math::X_AXIS, 0.27);
            v = turn(v, axis, double(spin) * (fresh ? 1.0 : fraction));
            constexpr double DIAGONAL = 0.7071067811865475;
            v = turn(v, math::Y_AXIS,
                     -double(displacement) * DIAGONAL *
                         math::RADIANS_PER_COLUMN<W>);
            vertices[i] = hs_test::displace_reference_rows<H>(
                v, float(displacement * DIAGONAL));
          }
          if (!replay)
            Plot::Mesh::draw_column_edges<W, H>(
                sink, canvas,
                Plot::Mesh::VertexSpan{{vertices, effect.NUM_NODES}}, edges,
                shader, geometry, 0.0f);
          const size_t count = replay ? effect.NUM_NODES : 1;
          for (size_t i = 0; i < count; ++i) {
            const auto point = math::vector_to_pixel<W, H>(vertices[i]);
            sink.plot(canvas, point.x, point.y, color.color, 0.0f,
                      color.alpha * (replay ? 1.0f : 0.5f));
          }
        }
      }
    }
    effect.advance_display();
    return capture_luminance<W, H>(effect);
  }

  /** @brief Checks moving copies across both arm cuts at one resolution.
   * @tparam W,H Display dimensions. */
  template <int W, int H> static void check_resolution() {
    for (bool replay : {false, true})
      for (float start : {0.0f, float(W / 2)})
        for (float spin : {-0.018f, 0.018f}) {
          auto render_scene = [=](bool timed, float fraction,
                                  float displacement) {
            return render<W, H>(timed, fraction, displacement, start, spin,
                                replay);
          };
          char name[64];
          std::snprintf(name, sizeof(name), "Dynamo %dx%d %s", W, H,
                        replay ? "replay" : "live");
          hs_test::seam_compare_render<W, H>(render_scene, name, 1.0f, true);
        }
  }
};

/** @brief Live and replayed copies match every column's frozen geometry. */
inline void test_dynamo_column_time_oracle() {
  DynamoTemporalOracle::check_live_history_fanout();
  DynamoTemporalOracle::check_resolution<288, 144>();
  DynamoTemporalOracle::check_resolution<96, 20>();
  for (bool replay : {false, true}) {
    auto render = [=](bool timed, float fraction, float displacement) {
      return DynamoTemporalOracle::render<288, 144>(
          timed, fraction, displacement, 0.0f, 0.018f, replay, false, true);
    };
    hs_test::seam_compare_render<288, 144>(render,
                                           replay ? "Dynamo forced-clip replay"
                                                  : "Dynamo forced-clip live");
  }
  Spacetime::geometry_counters = {};
  DynamoTemporalOracle::render<288, 144>(true, 1.0f, 0.0f, 0.0f, 0.018f, false,
                                         false, true);
  const auto &work = Spacetime::geometry_counters;
  HS_EXPECT_GT(work.column_history_visits, 0u);
  std::printf("  [Dynamo forced clip] callback evaluations=%u samples=%u "
              "edge_setups=%u culled=%u\n",
              work.endpoint_evaluations + work.midpoint_evaluations,
              work.column_samples, work.column_edges_prepared,
              work.column_edges_culled);
}

/** @brief Initialization and orientation reset hold the strand at the current endpoint. */
inline void test_dynamo_column_time_first_frame_and_reset() {
  for (bool replay : {false, true}) {
    const auto expected = DynamoTemporalOracle::render<96, 20>(
        false, 0.0f, 0.0f, 0.0f, 0.018f, replay, true);
    const auto actual = DynamoTemporalOracle::render<96, 20>(
        true, 0.0f, 0.0f, 0.0f, 0.018f, replay, true);
    float error = 0.0f;
    int lit = 0;
    for (size_t i = 0; i < actual.size(); ++i) {
      error += fabsf(actual[i] - expected[i]);
      lit += actual[i] > 0.0f;
    }
    std::printf("  Dynamo first frame replay=%d lit=%d\n", int(replay), lit);
    HS_EXPECT_GT(lit, 60);
    HS_EXPECT_LT(error, 0.01f);
  }
  const auto expected = DynamoTemporalOracle::render<96, 20>(
      false, 0.0f, 0.0f, 0.0f, 0.018f, false, true);
  const auto actual = DynamoTemporalOracle::render<96, 20>(
      true, 0.0f, 0.0f, 0.0f, 0.018f, false, true, false, true);
  float error = 0.0f;
  int lit = 0;
  for (size_t i = 0; i < actual.size(); ++i) {
    error += fabsf(actual[i] - expected[i]);
    lit += actual[i] > 0.0f;
  }
  std::printf("  Dynamo orientation reset lit=%d error=%g\n", lit, error);
  HS_EXPECT_GT(lit, 60);
  HS_EXPECT_LT(error, 0.01f);
}
