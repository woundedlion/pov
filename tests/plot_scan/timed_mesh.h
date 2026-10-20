/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @brief Minimal wireframe vertex storage for timed-mesh tests. */
struct TimedMeshFixture {
  std::array<math::Vector, 3>
      vertices; /**< Endpoints of two connected edges. */
};

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
