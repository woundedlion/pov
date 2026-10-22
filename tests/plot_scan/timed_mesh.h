/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @brief Minimal wireframe vertex storage for timed-mesh tests. */
struct TimedMeshFixture {
  std::array<math::Vector, 3>
      vertices; /**< Endpoints of two connected edges. */
};

/** @brief An axial callback retains its small column bound near the pole. */
inline void test_timed_axial_callback_bound() {
  constexpr int W = 288, H = 144;
  hs_test::StubEffect effect(W, H);
  Filter::Screen::DirectAntiAliasSink<W, H> sink;
  Canvas canvas(effect);
  sink.prepare(canvas);
  ScratchScope guard(scratch_arena_a);
  Fragments points;
  points.bind(scratch_arena_a, 2);
  for (float x : {40.0f, 50.0f}) {
    Fragment f;
    f.pos = math::pixel_to_vector<W, H>(x, 0.5f);
    points.push_back(f);
  }
  auto at = [&](size_t index, float t) {
    return math::rotate(points[index].pos,
                        math::make_rotation(math::Y_AXIS, (t - 1.0f) * 0.02f));
  };
  auto shader = [](const math::Vector &, Fragment &f) {
    f.color = Color4(Pixel(65535, 65535, 65535), 0.5f);
  };
  Temporal::geometry_counters = {};
  Plot::rasterize<W, H, {.timed = true}>(sink, canvas, points, shader, {},
                                         {at, false, {}, true}, 0.02f);
  HS_EXPECT_GT(Temporal::geometry_counters.end_probes, 0u);
  HS_EXPECT_EQ(Temporal::geometry_counters.held_anchors, 0u);
}

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

/** @brief Mesh callbacks retain vertex indices, edge registers and still pixels. */
inline void test_timed_mesh_vertex_callback_mapping() {
  constexpr int W = 288, H = 144;
  for (bool enabled : {false, true}) {
    Temporal::set_temporal_sampling_enabled(enabled);
    for (bool linear : {false, true}) {
      auto render = [&](bool mesh_entry) {
        hs_test::StubEffect effect(W, H);
        Filter::Screen::DirectAntiAliasSink<W, H> sink;
        {
          Canvas canvas(effect);
          sink.prepare(canvas);
          ScratchScope a(scratch_arena_a), b(scratch_arena_b);
          const TimedMeshFixture mesh{
              {math::pixel_to_vector<W, H>(30.0f, 60.0f),
               math::pixel_to_vector<W, H>(155.0f, 75.0f),
               math::pixel_to_vector<W, H>(125.0f, 65.0f)}};
          ArenaVector<Plot::Mesh::Edge> edges(scratch_arena_b, 2);
          edges.push_back({2, 0});
          edges.push_back({0, 1});
          auto at = [&](size_t index, float t) {
            return math::rotate(
                mesh.vertices[index],
                math::make_rotation(math::Y_AXIS, (t - 1.0f) * 0.04f));
          };
          auto shader = [](const math::Vector &, Fragment &f, float t) {
            f.color = Color4(Pixel(65535, 32768, 16384),
                             (0.2f + f.v2 * 0.1f) * (0.5f + 0.5f * t));
          };
          if (mesh_entry) {
            Plot::Mesh::draw_timed<W, H>(sink, canvas, mesh, edges, shader,
                                         {at, linear, {}}, 0.04f);
          } else {
            for (size_t edge = 0; edge < edges.size(); ++edge) {
              ScratchScope edge_scope(scratch_arena_a);
              Fragments points;
              points.bind(scratch_arena_a, 2);
              Fragment first, last;
              first.pos = mesh.vertices[edges[edge].u];
              last.pos = mesh.vertices[edges[edge].v];
              Plot::Line::sample(
                  points, first, last,
                  Plot::make_geodesic_edge_span(first.pos, last.pos));
              for (auto &point : points)
                point.v2 = static_cast<float>(edge);
              auto edge_at = [&](size_t index, float t) {
                return at(index == 0 ? edges[edge].u : edges[edge].v, t);
              };
              Plot::rasterize<W, H, {.timed = true}>(
                  sink, canvas, points, shader, {}, {edge_at, linear, {}},
                  0.04f);
            }
          }
        }
        effect.advance_display();
        std::vector<Pixel> pixels;
        for (int y = 0; y < H; ++y)
          for (int x = 0; x < W; ++x)
            pixels.push_back(effect.get_pixel(x, y));
        return pixels;
      };
      const auto mesh = render(true), reference = render(false);
      size_t differences = 0;
      for (size_t i = 0; i < mesh.size(); ++i)
        differences += mesh[i].r != reference[i].r ||
                       mesh[i].g != reference[i].g ||
                       mesh[i].b != reference[i].b;
      HS_EXPECT_EQ(differences, size_t{0});
    }
  }
  Temporal::set_temporal_sampling_enabled(true);
}

/** @brief Source mesh edges preserve identity and shade at column time. */
inline void test_timed_mesh_preserves_edge_register_and_time() {
  constexpr int W = 288, H = 144;
  hs_test::StubEffect effect(W, H);
  Filter::Screen::DirectAntiAliasSink<W, H> sink;
  Canvas canvas(effect);
  sink.prepare(canvas);
  ScratchScope scratch(scratch_arena_b);
  ArenaVector<Plot::Mesh::Edge> edges;
  edges.bind(scratch_arena_b, 2);
  edges.push_back({0, 1});
  edges.push_back({1, 2});
  const TimedMeshFixture mesh{{math::pixel_to_vector<W, H>(40.0f, 65.0f),
                               math::pixel_to_vector<W, H>(80.0f, 65.0f),
                               math::pixel_to_vector<W, H>(100.0f, 60.0f)}};
  std::array<int, 2> shaded{};
  auto shader = [&](const math::Vector &point, Fragment &fragment, float time) {
    const int edge = static_cast<int>(fragment.v2);
    HS_EXPECT_GE(edge, 0);
    HS_EXPECT_LT(edge, 2);
    if (edge >= 0 && edge < 2)
      ++shaded[edge];
    const auto pixel = math::vector_to_pixel<W, H>(point);
    HS_EXPECT_NEAR(time, Temporal::column_fraction<W>(pixel.x), 0.004f);
    fragment.color = Color4(Pixel(65535, 65535, 65535), 0.4f);
  };
  Plot::Mesh::draw_timed<W, H>(sink, canvas, mesh, edges, shader);
  HS_EXPECT_GT(shaded[0], 5);
  HS_EXPECT_GT(shaded[1], 5);
}

/** @brief Renders static source geometry through an ordinary or timed entry.
 * @param timed Selects the timed API.
 * @param mesh_mode Selects the mesh entry.
 * @param clipped Enables a segment clip.
 * @param replicated Emits three static spatial copies.
 * @return Frame pixels after the draw. */
inline std::vector<Pixel> render_still_timed_route(bool timed, bool mesh_mode,
                                                   bool clipped,
                                                   bool replicated) {
  constexpr int W = 288, H = 144;
  hs_test::StubEffect effect(W, H);
  if (clipped)
    effect.set_clip(36, 72, 0, 144);
  Pipeline<W, H, Filter::World::Replicate<W>, Filter::Screen::AntiAlias<W, H>>
      pipeline(Filter::World::Replicate<W>(replicated ? 3 : 1),
               Filter::Screen::AntiAlias<W, H>{});
  {
    Canvas canvas(effect);
    ScratchScope a(scratch_arena_a), b(scratch_arena_b);
    Fragments points;
    points.bind(scratch_arena_a, 3);
    const TimedMeshFixture mesh{{math::pixel_to_vector<W, H>(-12.f, 65.f),
                                 math::pixel_to_vector<W, H>(12.f, 65.f),
                                 math::pixel_to_vector<W, H>(65.f, 15.f)}};
    for (const auto &v : mesh.vertices) {
      Fragment f;
      f.pos = v;
      points.push_back(f);
    }
    ArenaVector<Plot::Mesh::Edge> edges;
    edges.bind(scratch_arena_b, 2);
    edges.push_back({0, 1});
    edges.push_back({1, 2});
    Temporal::geometry_counters = {};
    auto shader = white_shader(.4f);
    if (mesh_mode) {
      if (timed)
        Plot::Mesh::draw_timed<W, H>(pipeline, canvas, mesh, edges, shader);
      else
        Plot::Mesh::draw<W, H>(pipeline, canvas, mesh, edges,
                               FragmentShaderFn(shader));
    } else {
      if (timed)
        Plot::rasterize<W, H, {.timed = true}>(pipeline, canvas, points,
                                               shader);
      else
        Plot::rasterize<W, H>(pipeline, canvas, points, shader);
    }
    if (timed)
      HS_EXPECT_EQ(Temporal::geometry_counters.end_probes, 0u);
  }
  effect.advance_display();
  std::vector<Pixel> pixels;
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x)
      pixels.push_back(effect.get_pixel(x, y));
  return pixels;
}

/** @brief Static timed entries retain ordinary pixels and skip capture work. */
inline void test_timed_static_route_matches_ordinary() {
  for (bool replicated : {false, true})
    for (bool mesh : {false, true})
      for (bool clip : {false, true}) {
        const auto a = render_still_timed_route(true, mesh, clip, replicated),
                   b = render_still_timed_route(false, mesh, clip, replicated);
        size_t different = 0;
        for (size_t i = 0; i < a.size(); ++i)
          if (a[i].r != b[i].r || a[i].g != b[i].g || a[i].b != b[i].b)
            ++different;
        HS_EXPECT_EQ(different, 0u);
      }
}

/** @brief Opposite-longitude polar edges retain both arm branches. */
inline void test_timed_polar_mesh_branch_coverage() {
  constexpr int W = 288, H = 144;
  for (float pole : {-1.0f, 1.0f}) {
    hs_test::StubEffect effect(W, H);
    math::Orientation<> orientation(math::make_rotation(math::X_AXIS, 0.27f));
    const auto axis = math::Vector(0.2f, 1.0f, 0.1f).normalized();
    orientation.push(math::make_rotation(axis, 0.012f) * orientation.get());
    Pipeline<W, H, Filter::World::Orient, Filter::Screen::AntiAlias<W, H>>
        pipeline{Filter::World::Orient(orientation),
                 Filter::Screen::AntiAlias<W, H>{}};
    pipeline.init_storage(persistent_arena);
    {
      Canvas canvas(effect);
      ScratchScope scratch(scratch_arena_b);
      ArenaVector<Plot::Mesh::Edge> edges;
      edges.bind(scratch_arena_b, 1);
      edges.push_back({0, 1});
      const TimedMeshFixture mesh{
          {math::Vector(0.0f, pole * 1.61803398875f, 1.0f).normalized(),
           math::Vector(0.0f, pole * 1.61803398875f, -1.0f).normalized(),
           math::X_AXIS}};
      Plot::Mesh::draw_timed<W, H>(pipeline, canvas, mesh, edges,
                                   white_shader(0.4f));
    }
    effect.advance_display();
    uint64_t energy = 0;
    for (int y = 0; y < H; ++y) {
      const auto pixel = effect.get_pixel(216, y);
      energy += pixel.r + pixel.g + pixel.b;
    }
    HS_CONTEXT("pole=%f", pole);
    HS_EXPECT_GT(energy, uint64_t{0});
  }
}

/** @brief Zero-bound K2 callbacks preserve literal endpoint pixels. */
inline void test_timed_stationary_callback_pixel_parity() {
  constexpr int W = 288, H = 144;
  for (bool clipped : {false, true}) {
    const auto render = [&](bool timed) {
      hs_test::StubEffect effect(W, H);
      if (clipped)
        effect.set_clip(36, 72, 0, H);
      math::Orientation<> orientation(math::make_rotation(math::X_AXIS, 0.27f));
      Pipeline<W, H, Filter::World::Orient, Filter::Screen::AntiAlias<W, H>>
          pipeline{Filter::World::Orient(orientation),
                   Filter::Screen::AntiAlias<W, H>{}};
      pipeline.init_storage(persistent_arena);
      {
        Canvas canvas(effect);
        ScratchScope scratch(scratch_arena_a);
        Fragments points;
        points.bind(scratch_arena_a, 3);
        for (const auto &position :
             {math::pixel_to_vector<W, H>(-12.0f, 65.0f),
              math::pixel_to_vector<W, H>(12.0f, 65.0f),
              math::pixel_to_vector<W, H>(65.0f, 15.0f)}) {
          Fragment point;
          point.pos = position;
          points.push_back(point);
        }
        auto geometry = [&](size_t index, float) { return points[index].pos; };
        if (timed)
          Plot::rasterize<W, H, {.timed = true}>(
              pipeline, canvas, points, white_shader(0.4f), {},
              Plot::GeometryCallback::interpolating(geometry));
        else
          Plot::rasterize<W, H>(pipeline, canvas, points, white_shader(0.4f));
      }
      effect.advance_display();
      std::vector<Pixel> pixels;
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
          pixels.push_back(effect.get_pixel(x, y));
      return pixels;
    };
    const auto timed = render(true), ordinary = render(false);
    size_t differences = 0;
    for (size_t i = 0; i < timed.size(); ++i)
      if (timed[i].r != ordinary[i].r || timed[i].g != ordinary[i].g ||
          timed[i].b != ordinary[i].b)
        ++differences;
    HS_EXPECT_EQ(differences, size_t{0});
  }
}

/** @brief Cadence transitions retain the frozen column's coverage. */
inline void test_timed_mesh_cadence_transition() {
  constexpr int W = 288, H = 144, COLUMN = 228;
  const bool enabled = Temporal::temporal_sampling_enabled();
  Temporal::set_temporal_sampling_enabled(true);
  const auto render = [](bool timed, float displacement) {
    hs_test::StubEffect effect(W, H);
    Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipeline;
    const auto start = math::make_rotation(math::X_AXIS, 0.27f);
    const auto finish =
        math::make_rotation(math::Vector(0.2f, 1.0f, 0.1f).normalized(),
                            0.012f) *
        start;
    const auto at = [&](size_t index, float t) {
      const auto base = index == 0 ? -math::Y_AXIS : math::X_AXIS;
      const auto tangent = math::tangent_axis(base);
      const float phase =
          (0.2f + 0.004f * t) * 2.0f * math::PI_F + (index == 0 ? 0.3f : 0.0f);
      const auto offset = tangent * math::fast_cosf(phase) +
                          math::cross(base, tangent) * math::fast_sinf(phase);
      return math::rotate((base + offset * 0.2f).normalized(),
                          math::slerp(start, finish, t));
    };
    {
      Canvas canvas(effect);
      ScratchScope scratch(scratch_arena_b);
      ArenaVector<Plot::Mesh::Edge> edges(scratch_arena_b, 1);
      edges.push_back({0, 1});
      const float fraction = timed ? 1.0f : Temporal::sweep_fraction(COLUMN, W);
      const TimedMeshFixture mesh{
          {hs_test::displace_reference_rows<H>(at(0, fraction), displacement),
           hs_test::displace_reference_rows<H>(at(1, fraction), displacement),
           math::X_AXIS}};
      const auto shader = [](const math::Vector &, Fragment &f) {
        f.color = Color4(Pixel(65535, 65535, 65535), 0.5f * (1.0f - f.v0));
      };
      if (timed)
        Plot::Mesh::draw_timed<W, H>(pipeline, canvas, mesh, edges, shader,
                                     Plot::GeometryCallback::exact(at),
                                     0.0172f);
      else
        Plot::Mesh::draw<W, H>(pipeline, canvas, mesh, edges, shader);
    }
    effect.advance_display();
    std::array<float, H> column{};
    for (int y = 0; y < H; ++y)
      column[y] = effect.get_pixel(COLUMN, y).r / 65535.0f;
    return column;
  };
  const auto timed = render(true, 0), frozen = render(false, 0),
             shifted = render(false, 0.2f);
  float error = 0.0f, displacement = 0.0f;
  for (int y = 0; y < H; ++y) {
    error += fabsf(timed[y] - frozen[y]);
    displacement += fabsf(shifted[y] - frozen[y]);
  }
  std::printf("  mesh cadence column: error=%.9f displacement=%.9f\n", error,
              displacement);
  HS_EXPECT_GT(displacement, 0.0f);
  HS_EXPECT(error <= displacement + 1e-5f,
            "captured coverage stays within the displacement reference");
  Temporal::set_temporal_sampling_enabled(enabled);
}

/** @brief Stationary long edges shade at projected column time. */
inline void test_stationary_timed_mesh_shader_time() {
  constexpr int W = 288, H = 144;
  const bool enabled = Temporal::temporal_sampling_enabled();
  Temporal::set_temporal_sampling_enabled(true);
  const auto render = [](bool timed, float fraction, float displacement) {
    hs_test::StubEffect effect(W, H);
    Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipeline;
    {
      Canvas canvas(effect);
      ScratchScope scratch(scratch_arena_b);
      ArenaVector<Plot::Mesh::Edge> edges(scratch_arena_b, 1);
      edges.push_back({0, 1});
      const TimedMeshFixture mesh{
          {hs_test::displace_reference_rows<H>(
               math::pixel_to_vector<W, H>(190, 135), displacement),
           hs_test::displace_reference_rows<H>(
               math::pixel_to_vector<W, H>(285, 70), displacement),
           math::X_AXIS}};
      const auto shader = [](const math::Vector &, Fragment &f, float time) {
        f.color = Color4(Pixel(65535, 65535, 65535), 0.5f * time);
      };
      if (timed)
        Plot::Mesh::draw_timed<W, H>(pipeline, canvas, mesh, edges, shader);
      else {
        const auto frozen = [&](const math::Vector &p, Fragment &f) {
          shader(p, f, fraction);
        };
        Plot::Mesh::draw<W, H>(pipeline, canvas, mesh, edges, frozen);
      }
    }
    effect.advance_display();
    std::vector<float> pixels(W * H);
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x)
        pixels[y * W + x] = effect.get_pixel(x, y).r / 65535.0f;
    return pixels;
  };
  hs_test::seam_compare_render<W, H>(render, "stationary mesh fade");
  Temporal::set_temporal_sampling_enabled(enabled);
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
                                 {at, false, {}, false, {second, second}, true},
                                 0.01f);
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
      Plot::Mesh::draw_timed<W, H>(pipeline, canvas, mesh, edges, shader,
                                   {{}, false, {}, false, {}, true}, 0.0f);
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
    HS_EXPECT_EQ(Temporal::geometry_counters.legacy_edge_walks, 0u);
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
                                   {at, false, {}, false, {}, true}, motion);
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
                                 {at, true, {}, false, {}, true}, motion);
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

/** @brief Resolved cadence preserves endpoint ownership at an interval change. */
inline void test_timed_cadence_resolved_phase() {
  float first[] = {0.0f, 0.4f, 0.8f, 1.2f};
  float last[] = {0.0f, 0.6f, 1.2f, 1.8f};
  const Plot::TimedCadence cadence{first, last, 3};
  const auto before = cadence.at_time(0.25f);
  const auto crossing = cadence.at_time(0.5f);
  HS_EXPECT_EQ(before.intervals, size_t(3));
  HS_EXPECT_EQ(crossing.intervals, size_t(2));
  bool exists;
  HS_EXPECT_EQ(crossing.fraction_at(2, true, exists), 1.0f);
  HS_EXPECT(exists, "terminal sample exists at the shorter cadence endpoint");
  crossing.fraction_at(3, true, exists);
  HS_EXPECT(!exists, "no duplicate terminal sample after the cadence change");
  Plot::TimedCadence empty;
  HS_EXPECT_EQ(empty.at_time(0.5f).fraction_at(0, true, exists), 1.0f);
  HS_EXPECT(exists, "empty cadence retains its terminal sample");
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
  Plot::Mesh::draw_timed<W, H>(pipeline, canvas, mesh, edges, shader,
                               {at, true, {}, false, {}, true}, 0.0f);
  const auto &counts = Temporal::geometry_counters;
  HS_EXPECT_EQ(callbacks, 0u);
  HS_EXPECT_EQ(counts.endpoint_evaluations, 0u);
  HS_EXPECT_EQ(counts.column_draw_setups, 1u);
  HS_EXPECT_EQ(counts.legacy_edge_walks, 0u);
  HS_EXPECT_EQ(counts.column_crossings, 0u);
  HS_EXPECT_EQ(counts.capture_solves, 0u);
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
      geometry.column_edges = true;
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
                                       {at, false, {}, false, {}, true},
                                       fabsf(radians));
        } else {
          TimedMeshFixture frozen{{at(0, 0.0f), at(1, 0.0f), math::X_AXIS}};
          Plot::Mesh::draw_timed<W, H>(pipeline, canvas, frozen, edges, shader,
                                       {{}, false, {}, false, {}, true}, 0.0f);
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
                                   {at, false, {}, false, {}, true}, SHRINK);
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
