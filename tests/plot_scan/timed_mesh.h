/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @brief Minimal wireframe vertex storage for timed-mesh tests. */
struct TimedMeshFixture {
  std::array<math::Vector, 3>
      vertices; /**< Endpoints of two connected edges. */
};

/** @brief The hoisted trail gate retains edges whose motion reaches the clip. */
inline void test_trail_gate_uses_motion_region() {
  constexpr int W = 288, H = 144;
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipeline{
      Filter::Screen::AntiAlias<W, H>()};
  const ClipRegion clip{36, 72, 0, W, 1, W, H};
  const auto xc = clip.x_clip();
  ScratchScope guard(scratch_arena_a);
  Fragments points;
  points.bind(scratch_arena_a, 2);
  for (float x : {100.0f, 101.0f}) {
    Fragment f;
    f.pos = math::pixel_to_vector<W, H>(x, 32.0f);
    points.push_back(f);
  }
  uint8_t bits[1];
  const bool still =
      Plot::gate_trail_edges<W, H>(pipeline, clip, xc, points, bits);
  HS_EXPECT(!still, "endpoint edge is outside the clip");
  const bool moving =
      Plot::gate_trail_edges<W, H>(pipeline, clip, xc, points, bits, 0.1f);
  HS_EXPECT(moving, "swept edge reaches the clip");
  HS_EXPECT_EQ(bits[0], Plot::EDGE_VISIBLE);
}

/** @brief A column mesh evaluates each shared vertex once at the frame's
 * start, and once more at its middle when its chord bound fails. */
inline void test_timed_mesh_endpoint_cache() {
  constexpr int W = 288, H = 144;
  Temporal::set_temporal_sampling_enabled(true);
  for (float second : {0.0f, 8.0f}) {
    hs_test::StubEffect effect(W, H);
    Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipeline;
    Canvas canvas(effect);
    ScratchScope a(scratch_arena_a), b(scratch_arena_b);
    const TimedMeshFixture mesh{{math::pixel_to_vector<W, H>(40.0f, 60.0f),
                                 math::pixel_to_vector<W, H>(50.0f, 65.0f),
                                 math::pixel_to_vector<W, H>(60.0f, 70.0f)}};
    ArenaVector<Plot::Mesh::Edge> edges(scratch_arena_b, 2);
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    size_t calls = 0;
    auto at = [&](size_t index, float t) {
      ++calls;
      return math::rotate(
          mesh.vertices[index],
          math::make_rotation(math::Y_AXIS, (t - 1.0f) * 0.01f));
    };
    auto shader = [](const math::Vector &, Fragment &f) {
      f.color = Color4(Pixel(65535, 65535, 65535), 0.5f);
    };
    Plot::Mesh::draw_timed<W, H>(pipeline, canvas, mesh, edges, shader,
                                 {at, {second, second}}, 0.01f);
    HS_EXPECT_EQ(calls, second == 0.0f ? size_t(3) : size_t(6));
  }
}

/** @brief Column-walk edges retain geometry, even coverage and single
 * endpoints. */
inline void test_column_edge_quality() {
  constexpr int W = 288, H = 144;
  Temporal::set_temporal_sampling_enabled(true);
  for (int scene = 0; scene < 3; ++scene) {
    hs_test::StubEffect effect(W, H);
    Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipeline;
    uint32_t starts = 0, ends = 0;
    float geometry_error = 0.0f;
    Temporal::geometry_counters = {};
    {
      Canvas canvas(effect);
      ScratchScope storage(scratch_arena_a), cache(scratch_arena_b);
      const bool vertical = scene == 2;
      const float start = scene == 1 ? 254.0f : 110.0f;
      TimedMeshFixture mesh{};
      for (int i = 0; i < 3; ++i) {
        const int vertex = std::min(i, 1);
        mesh.vertices[i] =
            vertical
                ? math::pixel_to_vector<W, H>(40.3f, 20.0f + 100.0f * vertex)
                : math::Vector::from_spherical((start + 60.0f * vertex) *
                                                   math::RADIANS_PER_COLUMN<W>,
                                               math::PI_F * 0.5f);
      }
      ArenaVector<Plot::Mesh::Edge> edges(scratch_arena_b, 1);
      edges.push_back({0, 1});
      const auto normal =
          math::cross(mesh.vertices[0], mesh.vertices[1]).normalized();
      auto shader = [&](const math::Vector &p, Fragment &f, float) {
        starts += f.v0 == 0.0f;
        ends += f.v0 == 1.0f;
        geometry_error = fmaxf(geometry_error, fabsf(math::dot(normal, p)) *
                                                   math::ROWS_PER_RADIAN<H>);
        f.color = Color4(Pixel(65535, 65535, 65535), 0.5f);
      };
      Plot::Mesh::draw_timed<W, H>(pipeline, canvas, mesh, edges, shader, {},
                                   0.0f);
    }
    effect.advance_display();
    float sum = 0.0f, squares = 0.0f, low = 10.0f, high = 0.0f;
    const int count = scene == 2 ? 90 : 50;
    for (int i = 0; i < count; ++i) {
      float brightness = 0.0f;
      if (scene == 2) {
        for (int x = 0; x < W; ++x)
          brightness += effect.get_pixel(x, 25 + i).r / 65535.0f;
      } else {
        const int x = (scene == 1 ? 259 + i : 115 + i) % W;
        float centroid = 0.0f;
        for (int y = 0; y < H; ++y) {
          const float value = effect.get_pixel(x, y).r / 65535.0f;
          brightness += value;
          centroid += value * y;
        }
        HS_EXPECT_LE(
            fabsf(centroid / brightness - math::phi_to_y<H>(math::PI_F * 0.5f)),
            0.2f);
      }
      sum += brightness;
      squares += brightness * brightness;
      low = fminf(low, brightness);
      high = fmaxf(high, brightness);
    }
    const float mean = sum / count;
    const float variance = fmaxf(0.0f, squares / count - mean * mean);
    std::printf(
        "column quality scene=%d mean=%.6f cv=%.6f range=%.6f..%.6f geometry=%.6f samples=%u endpoints=%u/%u\n",
        scene, mean, sqrtf(variance) / mean, low, high, geometry_error,
        Temporal::geometry_counters.column_generated_samples, starts, ends);
    HS_EXPECT_LE(geometry_error, 0.2f);
    HS_EXPECT_GT(low, 0.48f);
    HS_EXPECT_LT(high, 0.60f);
    HS_EXPECT_LE(sqrtf(variance) / mean, 0.05f);
    HS_EXPECT_EQ(starts, 1u);
    HS_EXPECT_EQ(ends, 1u);
    HS_EXPECT_LT(Temporal::geometry_counters.column_generated_samples, 128u);
  }
}

/** @brief Meridian fractions enclose the column walk sampler's arc parameter. */
inline void test_column_crossing_parameter_bound() {
  constexpr float LENGTH = 3.13f;
  const math::Vector a(1.0f, 0.0f, 0.0f), tangent(0.0f, 0.0f, 1.0f);
  float sine, cosine;
  math::poly_sincosf_0_pi(LENGTH, sine, cosine);
  const math::Vector b(cosine, 0.0f, sine);
  float worst = 0.0f;
  for (int i = 1; i < 16384; ++i) {
    const float parameter = LENGTH * i / 16384.0f;
    math::poly_sincosf_0_pi(parameter, sine, cosine);
    float fraction;
    HS_EXPECT(Plot::column_edge_crossing(a, b, {cosine, 0.0f, sine}, LENGTH,
                                         tangent, fraction),
              "interior meridian crosses arc");
    const float error = fabsf(fraction * LENGTH - parameter);
    worst = fmaxf(worst, error);
    HS_EXPECT_LE(error, Plot::COLUMN_ARC_PARAMETER_PAD);
  }
  std::printf("column parameter maximum error %.9f rad\n", worst);
}

/** @brief Mesh batching rejects hidden caps and retains swept-in geometry. */
inline void test_column_mesh_batch_cull() {
  constexpr int W = 288, H = 144;
  Temporal::set_temporal_sampling_enabled(true);
  for (int scene = 0; scene < 3; ++scene) {
    const float column = scene == 0 ? 180.0f : scene == 1 ? 50.0f : 86.0f;
    const float motion = scene == 2 ? 0.75f : 0.01f;
    hs_test::StubEffect effect(W, H);
    effect.set_clip(36, 72, 0, H);
    Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipeline{
        Filter::Screen::AntiAlias<W, H>()};
    size_t calls = 0;
    Temporal::geometry_counters = {};
    {
      Canvas canvas(effect);
      ScratchScope a(scratch_arena_a), b(scratch_arena_b);
      TimedMeshFixture mesh{
          {math::pixel_to_vector<W, H>(column, 60.0f),
           math::pixel_to_vector<W, H>(column + 2.0f, 65.0f),
           math::pixel_to_vector<W, H>(column + 4.0f, 70.0f)}};
      ArenaVector<Plot::Mesh::Edge> edges(scratch_arena_b, 2);
      edges.push_back({0, 1});
      edges.push_back({1, 2});
      auto at = [&](size_t i, float t) {
        ++calls;
        const auto start = math::rotate(
            mesh.vertices[i], math::make_rotation(math::Y_AXIS, motion));
        return math::nlerp_unit(start, mesh.vertices[i], t);
      };
      auto shader = [](const math::Vector &, Fragment &f) {
        f.color = Color4(Pixel(65535, 32768, 16384), 0.5f);
      };
      Plot::Mesh::draw_timed<W, H>(pipeline, canvas, mesh, edges, shader,
                                   {at, {}}, motion);
    }
    effect.advance_display();
    uint64_t energy = 0;
    for (int y = 0; y < H; ++y)
      for (int x = 36; x < 72; ++x)
        energy += effect.get_pixel(x, y).r;
    HS_EXPECT_EQ(Temporal::geometry_counters.column_mesh_culls, 1u);
    HS_EXPECT_EQ(Temporal::geometry_counters.column_mesh_rejected,
                 scene == 0 ? 1u : 0u);
    HS_EXPECT_EQ(Temporal::geometry_counters.column_draw_setups,
                 scene == 0 ? 0u : 1u);
    HS_EXPECT_EQ(calls, scene == 0 ? size_t(0) : size_t(3));
    if (scene == 0)
      HS_EXPECT_EQ(energy, uint64_t(0));
    else
      HS_EXPECT_GT(energy, uint64_t(0));
  }
}

/** @brief Single-time splats engage within the projected motion budget and in
 * the interior columns of faster edges; a timed splat reuses its own column's
 * taps, evaluates at most its neighbor's and reruns no shader. */
inline void test_column_single_splat_motion_gate() {
  constexpr int W = 288, H = 144;
  Temporal::set_temporal_sampling_enabled(true);
  for (float motion : {0.001f, 0.8f}) {
    hs_test::StubEffect effect(W, H);
    Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipeline;
    Canvas canvas(effect);
    ScratchScope a(scratch_arena_a), b(scratch_arena_b);
    TimedMeshFixture mesh{{math::pixel_to_vector<W, H>(40.0f, 60.0f),
                           math::pixel_to_vector<W, H>(50.0f, 65.0f),
                           math::pixel_to_vector<W, H>(60.0f, 70.0f)}};
    ArenaVector<Plot::Mesh::Edge> edges(scratch_arena_b, 2);
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    auto at = [&](size_t i, float t) {
      return math::rotate(
          mesh.vertices[i],
          math::make_rotation(math::Y_AXIS, motion * (1.0f - t)));
    };
    uint32_t shaders = 0;
    auto shader = [&](const math::Vector &, Fragment &f) {
      ++shaders;
      f.color = Color4(Pixel(65535, 32768, 16384), 0.5f);
    };
    Temporal::geometry_counters = {};
    Plot::Mesh::draw_timed<W, H>(pipeline, canvas, mesh, edges, shader,
                                 {at, {}}, motion);
    const auto &counts = Temporal::geometry_counters;
    HS_EXPECT_GT(counts.column_samples, 0u);
    HS_EXPECT_EQ(shaders, counts.column_shader_calls);
    if (motion < 0.01f) {
      HS_EXPECT_EQ(counts.column_single_samples, counts.column_samples);
      HS_EXPECT_EQ(counts.column_tap_samples, 0u);
      HS_EXPECT_EQ(counts.column_projections, counts.column_samples);
    } else {
      HS_EXPECT_GT(counts.column_single_samples, 0u);
      HS_EXPECT_GT(counts.column_tap_samples, 0u);
      HS_EXPECT_EQ(counts.column_single_samples + counts.column_tap_samples,
                   counts.column_samples);
      HS_EXPECT_LE(counts.boundary_evaluations, counts.column_tap_samples);
    }
  }
}

/** @brief Stationary timed meshes skip trajectory caches and column preparation. */
inline void test_stationary_column_mesh_skips_geometry_cache() {
  constexpr int W = 288, H = 144;
  Temporal::set_temporal_sampling_enabled(true);
  hs_test::StubEffect effect(W, H);
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipeline;
  Canvas canvas(effect);
  ScratchScope a(scratch_arena_a), b(scratch_arena_b);
  TimedMeshFixture mesh{{math::pixel_to_vector<W, H>(40.0f, 60.0f),
                         math::pixel_to_vector<W, H>(50.0f, 65.0f),
                         math::pixel_to_vector<W, H>(60.0f, 70.0f)}};
  ArenaVector<Plot::Mesh::Edge> edges(scratch_arena_b, 2);
  edges.push_back({0, 1});
  edges.push_back({1, 2});
  uint32_t callbacks = 0, shaders = 0;
  auto at = [&](size_t i, float) {
    ++callbacks;
    return mesh.vertices[i];
  };
  auto shader = [&](const math::Vector &, Fragment &f, float t) {
    ++shaders;
    f.color = Color4(Pixel(65535, 32768, 16384), 0.5f * t);
  };
  Temporal::geometry_counters = {};
  Plot::Mesh::draw_timed<W, H>(pipeline, canvas, mesh, edges, shader, {at, {}},
                               0.0f);
  const auto &counts = Temporal::geometry_counters;
  HS_EXPECT_EQ(callbacks, 0u);
  HS_EXPECT_EQ(counts.endpoint_evaluations, 0u);
  HS_EXPECT_EQ(counts.column_draw_setups, 1u);
  HS_EXPECT_EQ(counts.column_crossings, 0u);
  HS_EXPECT_GT(shaders, 0u);
}

/** @brief A chained closed ring plots each vertex once, as its edge's start,
 *         still and moving. */
inline void test_column_walk_chained_ring_plots_vertices_once() {
  constexpr int W = 288, H = 144, N = 12;
  Temporal::set_temporal_sampling_enabled(true);
  struct Ring {
    std::array<math::Vector, N> vertices; /**< Closed ring vertices. */
  } ring;
  const math::Vector centre = math::pixel_to_vector<W, H>(100.0f, 60.0f);
  const math::Vector u = math::cross(centre, math::Y_AXIS).normalized();
  const math::Vector v = math::cross(centre, u);
  for (int i = 0; i < N; ++i) {
    const float turn = math::TWO_PI_F * float(i) / N;
    ring.vertices[i] =
        (centre * cosf(0.2f) + (u * cosf(turn) + v * sinf(turn)) * sinf(0.2f))
            .normalized();
  }
  for (float motion : {0.0f, 0.05f}) {
    HS_CONTEXT("motion x100", std::llround(motion * 100.0f));
    uint32_t samples[2] = {};
    for (bool chained : {false, true}) {
      hs_test::StubEffect effect(W, H);
      Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipeline;
      Canvas canvas(effect);
      ScratchScope a(scratch_arena_a), b(scratch_arena_b);
      ArenaVector<Plot::Mesh::Edge> edges(scratch_arena_b, N);
      for (int i = 0; i < N; ++i)
        edges.push_back(
            {static_cast<uint16_t>(i), static_cast<uint16_t>((i + 1) % N)});
      auto at = [&](size_t i, float t) {
        return math::rotate(
            ring.vertices[i],
            math::make_rotation(math::Y_AXIS, motion * (1.0f - t)));
      };
      auto shader = [&](const math::Vector &, Fragment &f, float) {
        f.color = Color4(Pixel(65535, 32768, 16384), 0.5f);
      };
      Plot::GeometryCallback geometry;
      geometry.evaluate = at;
      geometry.path = Temporal::PathBound::rotation(motion);
      geometry.chained = chained;
      Temporal::geometry_counters = {};
      Plot::Mesh::draw_timed<W, H>(pipeline, canvas, ring, edges, shader,
                                   geometry, motion);
      HS_EXPECT_EQ(Temporal::geometry_counters.column_draw_setups, 1u);
      samples[chained] = Temporal::geometry_counters.column_samples;
    }
    HS_EXPECT_GT(samples[1], 0u);
    HS_EXPECT_EQ(samples[0] - samples[1], uint32_t(N));
  }
}

/** @brief Refined basis normalization stays unit across its input range. */
inline void test_column_basis_normalization() {
  for (int i = 0; i <= 1024; ++i) {
    const float scale = 0.5f + float(i) / 1024.0f;
    const auto value = math::Vector(0.36f, 0.48f, 0.8f) * scale;
    const auto unit = Plot::ColumnEdgeFrame::normalize(value);
    HS_EXPECT_LE(fabsf(math::dot(unit, unit) - 1.0f), 3e-6f);
  }
}

/** @brief Basis crossings select the forward meridian in either arc direction. */
inline void test_column_basis_crossing_orientation() {
  const float tangent = sqrtf(0.75f);
  float fraction;
  for (float sign : {-1.0f, 1.0f}) {
    HS_EXPECT(Plot::column_basis_crossing(-sign * 0.5f, sign * 0.5f,
                                          sign * tangent, tangent, 0.5f,
                                          3.0f / math::PI_F, fraction),
              "forward crossing exists");
    HS_EXPECT_LE(fabsf(fraction - 0.5f), 2e-6f);
    HS_EXPECT(!Plot::column_basis_crossing(-sign * 0.5f, sign * 0.5f,
                                           sign * tangent, -tangent, -0.5f,
                                           3.0f / math::PI_F, fraction),
              "back meridian is rejected");
  }
}

/** @brief Edges crossing an arm boundary keep each pass's columns covered. */
inline void test_column_walk_arm_boundary_coverage() {
  constexpr int W = 288, H = 144;
  const bool enabled = Temporal::temporal_sampling_enabled();
  Temporal::set_temporal_sampling_enabled(true);
  for (float motion : {-3.0f, 3.0f}) {
    const float radians = motion * math::RADIANS_PER_COLUMN<W>;
    auto render = [&](bool timed, float fraction, float displacement) {
      hs_test::StubEffect effect(W, H);
      Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipeline;
      {
        Canvas canvas(effect);
        ScratchScope a(scratch_arena_a), b(scratch_arena_b);
        TimedMeshFixture mesh{{math::pixel_to_vector<W, H>(132.0f, 60.0f),
                               math::pixel_to_vector<W, H>(158.0f, 66.0f),
                               math::X_AXIS}};
        ArenaVector<Plot::Mesh::Edge> edges(scratch_arena_b, 1);
        edges.push_back({0, 1});
        const auto shader = [](const math::Vector &, Fragment &f) {
          f.color = Color4(Pixel(65535, 65535, 65535), 0.5f);
        };
        auto at = [&](size_t i, float t) {
          return hs_test::displace_reference_rows<H>(
              math::rotate(
                  mesh.vertices[i],
                  math::make_rotation(
                      math::Y_AXIS, radians * (1.0f - (timed ? t : fraction)))),
              displacement);
        };
        if (timed) {
          Plot::Mesh::draw_timed<W, H>(pipeline, canvas, mesh, edges, shader,
                                       {at, {}}, fabsf(radians));
        } else {
          TimedMeshFixture frozen{{at(0, 0.0f), at(1, 0.0f), math::X_AXIS}};
          Plot::Mesh::draw_timed<W, H>(pipeline, canvas, frozen, edges, shader,
                                       {}, 0.0f);
        }
      }
      effect.advance_display();
      std::vector<float> pixels(W * H);
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
          pixels[y * W + x] = effect.get_pixel(x, y).r / 65535.0f;
      return pixels;
    };
    hs_test::seam_compare_render<W, H>(render, motion < 0.0f
                                                   ? "arm boundary trailing"
                                                   : "arm boundary leading");
  }
  Temporal::set_temporal_sampling_enabled(enabled);
}

/** @brief Shrinking edges keep their terminal sample in the column walk. */
inline void test_column_walk_terminal_on_shrinking_edge() {
  constexpr int W = 288, H = 144;
  const bool enabled = Temporal::temporal_sampling_enabled();
  Temporal::set_temporal_sampling_enabled(true);
  for (int start = 0; start < 16; ++start) {
    hs_test::StubEffect effect(W, H);
    Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipeline;
    uint32_t ends = 0;
    {
      Canvas canvas(effect);
      ScratchScope storage(scratch_arena_a), cache(scratch_arena_b);
      TimedMeshFixture mesh{};
      for (int i = 0; i < 3; ++i)
        mesh.vertices[i] = math::pixel_to_vector<W, H>(
            i == 0 ? 40.0f + 0.37f * start : 61.0f, i == 0 ? 60.0f : 60.25f);
      ArenaVector<Plot::Mesh::Edge> edges(scratch_arena_b, 1);
      edges.push_back({0, 1});
      constexpr float SHRINK = 0.6f * math::RADIANS_PER_COLUMN<W>;
      auto at = [&](size_t i, float t) {
        return i == 0 ? mesh.vertices[0]
                      : math::rotate(mesh.vertices[i],
                                     math::make_rotation(math::Y_AXIS,
                                                         -SHRINK * (1.0f - t)));
      };
      auto shader = [&](const math::Vector &, Fragment &f, float) {
        ends += f.v0 == 1.0f;
        f.color = Color4(Pixel(65535, 65535, 65535), 0.5f);
      };
      Plot::Mesh::draw_timed<W, H>(pipeline, canvas, mesh, edges, shader,
                                   {at, {}}, SHRINK);
    }
    HS_CONTEXT("shrinking edge start", start);
    HS_EXPECT_EQ(ends, 1u);
  }
  Temporal::set_temporal_sampling_enabled(enabled);
}

/** @brief The swept edge cull never rejects an edge whose arcs over the frame
 * reach the draw's band. */
inline void test_column_edge_cull_is_conservative() {
  constexpr int W = 288, H = 144;
  uint32_t rejected = 0, checked = 0;
  hs::random().seed(20261021u);
  const auto unit = []() {
    return math::Vector(hs::rand_f(-1.0f, 1.0f), hs::rand_f(-1.0f, 1.0f),
                        hs::rand_f(-1.0f, 1.0f))
        .normalized();
  };
  for (const auto band :
       {std::array<int, 4>{0, 72, 0, 144}, std::array<int, 4>{0, 72, 144, 288},
        std::array<int, 4>{36, 72, 100, 136}, std::array<int, 4>{0, H, 30, 230},
        std::array<int, 4>{0, 30, 0, W}, std::array<int, 4>{110, H, 0, 144}}) {
    ClipRegion clip;
    clip.w = W;
    clip.h = H;
    clip.y_start = band[0];
    clip.y_end = band[1];
    clip.x_start = band[2];
    clip.x_end = band[3];
    const float top = cosf(
        math::DisplayGeometry<H>::row_to_phi(float(std::max(0, band[0] - 2))));
    const float bottom = cosf(
        math::DisplayGeometry<H>::row_to_phi(float(std::min(H, band[1] + 2))));
    const auto columns = clip.x_clip();
    const Plot::ColumnEdgeCull<W> cull(columns, top, bottom);
    for (float radians : {0.0f, 0.01f, 0.06f}) {
      for (int trial = 0; trial < 4000; ++trial) {
        const auto a = unit();
        const auto b = math::rotate(
            a, math::make_rotation(math::cross(a, unit()).normalized(),
                                   hs::rand_f(0.0f, 0.7f)));
        const auto a1 = math::rotate(
            a, math::make_rotation(unit(), hs::rand_f(0.0f, radians)));
        const auto b1 = math::rotate(
            b, math::make_rotation(unit(), hs::rand_f(0.0f, radians)));
        const bool rejects = cull.rejects({a, a1, b, b1}, 1e-4f);
        rejected += rejects;
        bool reaches = false;
        for (int move = 0; move <= 8 && !reaches; ++move) {
          const float u = float(move) / 8.0f;
          const auto ma = math::nlerp_unit(a, a1, u);
          const auto mb = math::nlerp_unit(b, b1, u);
          for (int k = 0; k <= 32 && !reaches; ++k) {
            const auto p = math::slerp(ma, mb, float(k) / 32.0f);
            const float column = math::vector_to_theta<W>(p);
            const int owner = math::fast_wrap(int(floorf(column + 0.5f)), W);
            reaches = p.y <= top && p.y >= bottom &&
                      (!columns.active || !columns.clipped(owner) ||
                       !columns.clipped(math::fast_wrap(owner + 1, W)) ||
                       !columns.clipped(math::fast_wrap(owner - 1, W)));
          }
        }
        checked += reaches;
        if (reaches)
          HS_EXPECT(!rejects, "a reachable edge is kept");
      }
    }
  }
  std::printf("  [column edge cull] rejected %u, reachable %u\n", rejected,
              checked);
  HS_EXPECT_GT(rejected, 1000u);
}

/** @brief A canvas's column cull keeps every padded polyline, or the padded
 *         cap around it, that can reach the walk's rows and band. */
inline void test_column_edge_cull_keeps_padded_point_sets() {
  constexpr int W = 288, H = 144;
  uint32_t rejected = 0, cap_rejected = 0, reachable = 0;
  hs::random().seed(20261022u);
  const auto unit = []() {
    return math::Vector(hs::rand_f(-1.0f, 1.0f), hs::rand_f(-1.0f, 1.0f),
                        hs::rand_f(-1.0f, 1.0f))
        .normalized();
  };
  for (const auto band :
       {std::array<int, 4>{36, 72, 0, 144}, std::array<int, 4>{0, 72, 144, 288},
        std::array<int, 4>{110, H, 30, 230}}) {
    hs_test::StubEffect effect(W, H);
    effect.set_clip(band[0], band[1], band[2], band[3]);
    Canvas canvas(effect);
    const auto cull = Plot::ColumnEdgeCull<W>::of<H>(canvas);
    const float top = Plot::ColumnEdgeCull<W>::top_of<H>(canvas);
    const float bottom = Plot::ColumnEdgeCull<W>::bottom_of<H>(canvas);
    HS_EXPECT_NEAR(top,
                   cosf(math::DisplayGeometry<H>::row_to_phi(
                       float(std::max(0, band[0] - 2)))),
                   1e-6f);
    const auto &columns = canvas.render_x_clip();
    for (float pad : {0.0f, 0.02f, 0.06f}) {
      for (int trial = 0; trial < 1500; ++trial) {
        const int count = 2 + trial % 8;
        math::Vector points[9];
        points[0] = unit();
        for (int i = 1; i < count; ++i)
          points[i] =
              math::rotate(points[i - 1],
                           math::make_rotation(
                               math::cross(points[i - 1], unit()).normalized(),
                               hs::rand_f(0.0f, 0.08f)));
        float reach = 1.0f;
        for (int i = 0; i < count; ++i)
          reach = fminf(reach, math::dot(points[count / 2], points[i]));
        const bool rejects = cull.rejects(points, count, pad);
        const bool cap_rejects =
            cull.rejects(&points[count / 2], 1,
                         pad + sqrtf(fmaxf(0.0f, 2.0f - 2.0f * reach)));
        rejected += rejects;
        cap_rejected += cap_rejects;
        bool reaches = false;
        for (int i = 0; i + 1 < count && !reaches; ++i)
          for (int k = 0; k <= 16 && !reaches; ++k) {
            const math::Vector at =
                math::slerp(points[i], points[i + 1], float(k) / 16.0f);
            const math::Vector east =
                math::cross(math::Y_AXIS, at).magnitude() > 1e-4f
                    ? math::cross(math::Y_AXIS, at).normalized()
                    : math::X_AXIS;
            const math::Vector north = math::cross(at, east);
            for (int d = 0; d <= 8 && !reaches; ++d) {
              const float angle = math::TWO_PI_F * float(d) / 8.0f;
              const math::Vector p =
                  d == 8
                      ? at
                      : (at + (east * cosf(angle) + north * sinf(angle)) * pad)
                            .normalized();
              const float column = math::vector_to_theta<W>(p);
              const int owner = math::fast_wrap(int(floorf(column + 0.5f)), W);
              reaches = p.y <= top && p.y >= bottom &&
                        (!columns.active || !columns.clipped(owner) ||
                         !columns.clipped(math::fast_wrap(owner + 1, W)) ||
                         !columns.clipped(math::fast_wrap(owner - 1, W)));
            }
          }
        reachable += reaches;
        if (reaches) {
          HS_EXPECT(!rejects, "a reachable polyline is kept");
          HS_EXPECT(!cap_rejects, "a reachable polyline's cap is kept");
        }
      }
    }
  }
  std::printf("  [column point-set cull] rejected %u, cap %u, reachable %u\n",
              rejected, cap_rejected, reachable);
  HS_EXPECT_GT(rejected, 1000u);
  HS_EXPECT_GT(cap_rejected, 1000u);
}

/** @brief An edge's steady lattice is its walk lattice when the count holds
 *         at both frame ends, and lattice points lie on the arc where the
 *         walk places them. */
inline void test_column_steady_lattice() {
  constexpr int W = 288, H = 144;
  using Frame = Plot::ColumnEdgeFrame;
  hs::random().seed(20261023u);
  int steady = 0, unsteady = 0, beyond = 0;
  for (int trial = 0; trial < 3000; ++trial) {
    const math::Vector a = math::pixel_to_vector<W, H>(
        hs::rand_f(0.0f, float(W)), hs::rand_f(8.0f, 136.0f));
    const float heading = hs::rand_f(0.0f, math::TWO_PI_F);
    const float length = hs::rand_f(0.002f, 0.12f);
    const math::Vector east = math::cross(math::Y_AXIS, a).normalized();
    const math::Vector north = math::cross(a, east);
    const math::Vector b = math::rotate(
        a, math::make_rotation(east * sinf(heading) - north * cosf(heading),
                               length));
    const auto motion = math::make_rotation(
        math::Vector(0.3f, 1.0f, 0.2f).normalized(), hs::rand_f(-0.03f, 0.03f));
    const math::Vector starts[2] = {math::rotate(a, motion),
                                    math::rotate(b, motion)};
    const math::Vector ends[2] = {a, b};
    Plot::GeodesicEdgeSampler arcs[2];
    const int got = Frame::steady_lattice<W, H>(starts, ends, 4, arcs);
    const int start_count = Frame::lattice(
        Frame::pieces_of<W, H>(Frame::arc(starts[0], starts[1])));
    const int end_count =
        Frame::lattice(Frame::pieces_of<W, H>(Frame::arc(ends[0], ends[1])));
    const int want = start_count == end_count && end_count <= 4 ? end_count : 0;
    HS_EXPECT_EQ(got, want);
    steady += want != 0;
    unsteady += start_count != end_count;
    beyond += end_count > 4;
    if (got == 0)
      continue;
    for (int k = 1; k < got; ++k) {
      const float s = float(k) / float(got);
      for (int end = 0; end < 2; ++end) {
        const math::Vector p = Frame::lattice_point(arcs[end], s);
        HS_EXPECT_NEAR(math::dot(p, p), 1.0f, 2e-6f);
        HS_EXPECT_NEAR(math::angle_between(p, arcs[end].v1),
                       s * arcs[end].total_dist, 2e-4f);
        const math::Vector &far = end ? ends[1] : starts[1];
        HS_EXPECT_NEAR(math::angle_between(p, far),
                       (1.0f - s) * arcs[end].total_dist, 2e-4f);
      }
    }
  }
  std::printf("  [steady lattice] steady %d, unsteady %d, beyond %d\n", steady,
              unsteady, beyond);
  HS_EXPECT_GT(steady, 100);
  HS_EXPECT_GT(unsteady, 10);
  HS_EXPECT_GT(beyond, 100);
}

/** @brief A per-vertex path bound replaces the shared one when the column
 *         walk picks each vertex's interpolant. */
inline void test_column_capture_reads_vertex_paths() {
  constexpr int W = 288, H = 144;
  Temporal::set_temporal_sampling_enabled(true);
  hs_test::StubEffect effect(W, H);
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipeline;
  Canvas canvas(effect);
  ScratchScope a(scratch_arena_a), b(scratch_arena_b);
  TimedMeshFixture mesh{{math::pixel_to_vector<W, H>(40.0f, 60.0f),
                         math::pixel_to_vector<W, H>(50.0f, 65.0f),
                         math::pixel_to_vector<W, H>(60.0f, 70.0f)}};
  ArenaVector<Plot::Mesh::Edge> edges(scratch_arena_b, 2);
  edges.push_back({0, 1});
  edges.push_back({1, 2});
  auto at = [&](size_t i, float t) {
    return math::rotate(mesh.vertices[i],
                        math::make_rotation(math::Y_AXIS, 0.02f * (1.0f - t)));
  };
  auto shader = [&](const math::Vector &, Fragment &f, float) {
    f.color = Color4(Pixel(65535, 32768, 16384), 0.5f);
  };
  for (const bool curved : {false, true}) {
    auto vertex_path = [&](size_t i) {
      return curved && i == 1 ? Temporal::PathBound{1.0f, 1.0f}
                              : Temporal::PathBound{};
    };
    Plot::GeometryCallback geometry;
    geometry.evaluate = at;
    geometry.path = Temporal::PathBound{1.0f, 1.0f};
    geometry.vertex_path = vertex_path;
    Temporal::geometry_counters = {};
    Plot::Mesh::draw_timed<W, H>(pipeline, canvas, mesh, edges, shader,
                                 geometry, 0.02f);
    HS_EXPECT_EQ(Temporal::geometry_counters.midpoint_evaluations,
                 curved ? 1u : 0u);
  }
}

/** @brief Each arm pass's walk starts no later than its first lit column and
 * gives up no earlier, with little slack, and passes nothing lights get no
 * walk. */
inline void test_column_plan_covers_lit_columns() {
  constexpr int W = 288, PASS = W / 2;
  hs::random().seed(20261022u);
  const auto unit = []() {
    return math::Vector(hs::rand_f(-1.0f, 1.0f), hs::rand_f(-1.0f, 1.0f),
                        hs::rand_f(-1.0f, 1.0f))
        .normalized();
  };
  const auto basis = [](const math::Vector &a, const math::Vector &b) {
    const auto span = Plot::make_geodesic_edge_span(a, b);
    return Plot::GeodesicEdgeSampler{
        a, span.have_axis ? math::cross(span.axis, a) : math::Vector{},
        span.total};
  };
  uint32_t lit_passes = 0, empty_passes = 0, skipped = 0, slack_count = 0;
  float slack = 0.0f;
  for (int trial = 0; trial < 3000; ++trial) {
    auto a0 = unit();
    if (trial % 4 == 0)
      a0 = math::Vector(a0.x * 0.08f, a0.y < 0.0f ? -1.0f : 1.0f, a0.z * 0.08f)
               .normalized();
    const auto b0 = math::rotate(
        a0, math::make_rotation(math::cross(a0, unit()).normalized(),
                                hs::rand_f(0.01f, 1.5f)));
    const auto axis = unit();
    const float motion = hs::rand_f(0.0f, 0.15f);
    const auto a1 = math::rotate(a0, math::make_rotation(axis, motion));
    const auto b1 = math::rotate(
        b0, math::make_rotation(axis, motion * hs::rand_f(0.5f, 1.5f)));
    Plot::ColumnEdgeFrame frame;
    frame.bases[0] = basis(a0, b0);
    frame.bases[1] = basis(a1, b1);
    frame.ends[0] = b0;
    frame.ends[1] = b1;
    frame.deviation = 0.0f;
    bool lit[W];
    for (int x = 0; x < W; ++x) {
      const float t = Temporal::column_fraction<W>(float(x));
      const float a = math::vector_to_theta<W>(math::nlerp_unit(a0, a1, t));
      const float span = Temporal::centered_remainder(
          math::vector_to_theta<W>(math::nlerp_unit(b0, b1, t)) - a, float(W));
      const float low = span < 0.0f ? a + span : a;
      lit[x] =
          Temporal::centered_remainder(float(x) - low, float(W)) >= -0.5f &&
          Temporal::centered_remainder(float(x) - low, float(W)) <=
              fabsf(span) + 0.5f;
    }
    for (const int direction : {1, -1}) {
      Plot::ColumnSegment segments[2];
      const int count = frame.plan<W>(direction, segments);
      for (int pass = 0; pass < 2; ++pass) {
        const auto offset = [&](int x) {
          return direction > 0 ? x - pass * PASS : pass * PASS + PASS - 1 - x;
        };
        int first = -1;
        for (int k = 0; k < PASS && first < 0; ++k) {
          const int x =
              direction > 0 ? pass * PASS + k : pass * PASS + PASS - 1 - k;
          if (lit[x])
            first = k;
        }
        const Plot::ColumnSegment *walk = nullptr;
        for (int i = 0; i < count; ++i)
          if (segments[i].start / PASS == pass)
            walk = &segments[i];
        if (first < 0) {
          ++empty_passes;
          skipped += walk == nullptr;
          continue;
        }
        ++lit_passes;
        HS_EXPECT(walk != nullptr, "a lit pass gets a walk");
        if (walk == nullptr)
          continue;
        const int start = offset(walk->start);
        HS_EXPECT(start <= first, "the walk starts at or before the lit run");
        HS_EXPECT(first - start < walk->deadline,
                  "the walk reaches the lit run before giving up");
        HS_EXPECT_EQ(start + walk->count, PASS);
        if (walk->deadline < PASS) {
          slack += float(first - start);
          ++slack_count;
        }
      }
    }
  }
  const float mean = slack / float(std::max(slack_count, 1u));
  std::printf("  [column plan] lit passes %u (bounded %u), empty %u skipped "
              "%u, mean leading slack %.2f\n",
              lit_passes, slack_count, empty_passes, skipped, mean);
  HS_EXPECT_GT(lit_passes, 2000u);
  HS_EXPECT_GT(slack_count, lit_passes * 3 / 4);
  HS_EXPECT_GT(skipped, empty_passes * 9 / 10);
  HS_EXPECT(mean < 2.5f, "walks start close to their lit runs");
}
