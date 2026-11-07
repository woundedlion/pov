/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @brief Minimal wireframe vertex storage for timed-mesh tests. */
struct TimedMeshFixture {
  std::array<math::Vector, 3>
      vertices; /**< Endpoints of two connected edges. */
};

/** @brief Timed edge registers collapse rounded-zero arcs to their origin. */
inline void test_timed_edge_endpoint_registers() {
  reset_globals();
  ScratchScope scope(scratch_arena_a);
  Fragments points;
  points.bind(scratch_arena_a, 2);
  const math::Vector near(1.0f, 1e-8f, 0.0f);
  Plot::prepare_timed_edge(points, math::X_AXIS, near, 7, 0.0f);
  HS_EXPECT_EQ(points.size(), 2u);
  HS_EXPECT_EQ(points[0].pos, math::X_AXIS);
  HS_EXPECT_EQ(points[1].pos, math::X_AXIS);
  HS_EXPECT_EQ(points[0].v0, 0.0f);
  HS_EXPECT_EQ(points[1].v0, 0.0f);
  HS_EXPECT_EQ(points[0].v2, 7.0f);
  HS_EXPECT_EQ(points[1].v2, 7.0f);
  HS_EXPECT_EQ(points[1].v1, 0.0f);

  Plot::prepare_timed_edge(points, math::X_AXIS, math::Z_AXIS, 3,
                           math::PI_F / 2.0f);
  HS_EXPECT_EQ(points.size(), 2u);
  HS_EXPECT_EQ(points[0].pos, math::X_AXIS);
  HS_EXPECT_EQ(points[1].pos, math::Z_AXIS);
  HS_EXPECT_EQ(points[0].v0, 0.0f);
  HS_EXPECT_EQ(points[1].v0, 1.0f);
  HS_EXPECT_EQ(points[0].v2, 3.0f);
  HS_EXPECT_EQ(points[1].v2, 3.0f);
  HS_EXPECT_EQ(points[1].v1, math::PI_F / 2.0f);
}

/** @brief Creates three vertices for two connected timed-mesh edges. */
template <int W, int H> inline TimedMeshFixture chain_mesh() {
  return {{math::pixel_to_vector<W, H>(40.0f, 60.0f),
           math::pixel_to_vector<W, H>(50.0f, 65.0f),
           math::pixel_to_vector<W, H>(60.0f, 70.0f)}};
}

/** @brief Creates one timed-mesh edge in the supplied arena. */
inline ArenaVector<Plot::Mesh::Edge> single_edge(Arena &arena) {
  ArenaVector<Plot::Mesh::Edge> edges(arena, 1);
  edges.push_back({0, 1});
  return edges;
}

/** @brief Creates the two connected edges in the supplied arena. */
inline ArenaVector<Plot::Mesh::Edge> chain_edges(Arena &arena) {
  ArenaVector<Plot::Mesh::Edge> edges(arena, 2);
  edges.push_back({0, 1});
  edges.push_back({1, 2});
  return edges;
}

/** @brief Pole-floored long edges reach their far endpoint within the budget. */
inline void test_column_long_polar_edge_finishes() {
  constexpr int W = 288, H = 144;
  const ScopedTemporalSampling sampling(true);
  const math::Vector MID(sinf(0.02f), cosf(0.02f), 0);
  const float HALF_LENGTH = math::PI_F / 4;
  const math::Vector A =
      MID * cosf(HALF_LENGTH) - math::Z_AXIS * sinf(HALF_LENGTH);
  const math::Vector B =
      MID * cosf(HALF_LENGTH) + math::Z_AXIS * sinf(HALF_LENGTH);
  for (bool reversed : {false, true}) {
    hs_test::StubEffect effect(W, H);
    Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipeline;
    TimedMeshFixture mesh{{reversed ? B : A, reversed ? A : B, {}}};
    const auto ARC =
        Plot::EdgeWorldsheet::arc(mesh.vertices[0], mesh.vertices[1]);
    HS_EXPECT_GT((Plot::EdgeWorldsheet::pieces_of<W, H>(ARC)),
                 2 * Plot::rasterize_step_budget<W>());
    Spacetime::geometry_counters = {};
    hs::g_scan_metrics = {};
    uint32_t ends = 0;
    {
      Canvas canvas(effect);
      ScratchScope a(scratch_arena_a), b(scratch_arena_b);
      auto edges = single_edge(scratch_arena_b);
      auto at = [&](size_t vertex, float) { return mesh.vertices[vertex]; };
      auto shade = [&](const math::Vector &, Fragment &f, float) {
        ends += f.v0 == 1.0f;
        f.color = Color4(Pixel(65535, 65535, 65535), 0.5f);
      };
      const Plot::GeometryCallback geometry{at, {}, false};
      Plot::Mesh::draw_timed<W, H>(pipeline, canvas, mesh, edges, shade,
                                   geometry);
    }
    effect.advance_display();
    const auto PIXEL = math::vector_to_pixel<W, H>(mesh.vertices[1]);
    uint64_t far_energy = 0;
    for (int y = int(PIXEL.y) - 1; y <= int(PIXEL.y) + 2; ++y)
      for (int x = int(PIXEL.x) - 1; x <= int(PIXEL.x) + 2; ++x)
        far_energy += effect.get_pixel((x + W) % W, std::clamp(y, 0, H - 1)).r;
    HS_EXPECT_GT(ends, 0u);
    HS_EXPECT_GT(far_energy, uint64_t{0});
    HS_EXPECT_LT(Spacetime::geometry_counters.column_generated_samples,
                 2 * Plot::rasterize_step_budget<W>());
    HS_EXPECT_EQ(hs::g_scan_metrics.plot_backstop_hits, uint64_t{0});
    std::printf("  long polar edge reversed=%d samples=%u far_energy=%llu\n",
                reversed, Spacetime::geometry_counters.column_generated_samples,
                static_cast<unsigned long long>(far_energy));
  }
}

/** @brief Polar arcs cover their independently intersected meridians. */
inline void test_column_polar_arc_coverage() {
  constexpr int W = 288, H = 144;
  const ScopedTemporalSampling temporal_sampling(!Spacetime::flat_surface());
  using Vec = std::array<double, 3>;
  const auto cross = [](const Vec &a, const Vec &b) -> Vec {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]};
  };
  const auto dot = [](const Vec &a, const Vec &b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  };
  const auto angle = [&](const Vec &a, const Vec &b) {
    const auto normal = cross(a, b);
    return std::atan2(std::sqrt(dot(normal, normal)), dot(a, b));
  };
  for (float hemisphere : {-1.0f, 1.0f})
    for (bool reversed : {false, true})
      for (bool timed : {false, true}) {
        hs_test::StubEffect effect(W, H);
        Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipeline;
        TimedMeshFixture mesh{{
            math::Vector(-0.359562576f, hemisphere * 0.882092655f,
                         -0.304347575f),
            math::Vector(0.157168522f, hemisphere * 0.987457812f, 0.015009786f),
            math::Vector{},
        }};
        if (reversed)
          std::swap(mesh.vertices[0], mesh.vertices[1]);
        Spacetime::set_flat_surface(!timed);
        {
          Canvas canvas(effect);
          ScratchScope a(scratch_arena_a), b(scratch_arena_b);
          auto edges = single_edge(scratch_arena_b);
          auto shader = [](const math::Vector &, Fragment &f, float) {
            f.color = Color4(Pixel(65535, 65535, 65535), 0.2f);
          };
          auto at = [&](size_t vertex, float) { return mesh.vertices[vertex]; };
          const Plot::GeometryCallback geometry{at, {}, false};
          Plot::Mesh::draw_timed<W, H>(pipeline, canvas, mesh, edges, shader,
                                       geometry);
        }
        effect.advance_display();
        const Vec a{mesh.vertices[0].x, mesh.vertices[0].y, mesh.vertices[0].z};
        const Vec b{mesh.vertices[1].x, mesh.vertices[1].y, mesh.vertices[1].z};
        const auto normal = cross(a, b);
        const double length = angle(a, b);
        const double pi = std::acos(-1.0);
        int checked = 0, missing = 0;
        for (int column = 0; column < W; ++column) {
          const double theta = column * 2.0 * pi / W;
          const double sine = std::sin(theta), cosine = std::cos(theta);
          const double horizontal = normal[0] * cosine + normal[2] * sine;
          const Vec point{std::abs(normal[1]) * cosine,
                          -std::copysign(1.0, normal[1]) * horizontal,
                          std::abs(normal[1]) * sine};
          if (angle(a, point) + angle(point, b) > length + 1e-7)
            continue;
          const double phi =
              std::atan2(std::hypot(point[0], point[2]), point[1]);
          const double row = (phi - math::DisplayGeometry<H>::NORTH_PHI) *
                             math::ROWS_PER_RADIAN<H>;
          const double pole_row = hemisphere > 0.0f ? row : H - 1 - row;
          if (pole_row < 2.0 || pole_row > 11.0)
            continue;
          ++checked;
          uint64_t coverage = 0;
          for (int y = std::max(0, int(std::floor(row)) - 1);
               y <= std::min(H - 1, int(std::ceil(row)) + 1); ++y)
            coverage += effect.get_pixel(column, y).r;
          missing += coverage == 0;
        }
        HS_CONTEXT("polar arc hemisphere", static_cast<int>(hemisphere));
        HS_CONTEXT("polar arc reversed", reversed);
        HS_CONTEXT("polar arc timed", timed);
        HS_EXPECT_GT(checked, 32);
        HS_EXPECT_EQ(missing, 0);
      }
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

/** @brief A column mesh evaluates each shared vertex once at the frame's
 * start, and once more at its middle when its chord bound fails. */
inline void test_timed_mesh_endpoint_cache() {
  constexpr int W = 288, H = 144;
  Spacetime::set_flat_surface(false);
  for (float second : {0.0f, 8.0f}) {
    hs_test::StubEffect effect(W, H);
    Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipeline;
    Canvas canvas(effect);
    ScratchScope a(scratch_arena_a), b(scratch_arena_b);
    const TimedMeshFixture mesh = chain_mesh<W, H>();
    auto edges = chain_edges(scratch_arena_b);
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
  Spacetime::set_flat_surface(false);
  for (int scene = 0; scene < 3; ++scene) {
    hs_test::StubEffect effect(W, H);
    Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipeline;
    uint32_t starts = 0, ends = 0;
    float geometry_error = 0.0f;
    Spacetime::geometry_counters = {};
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
      auto edges = single_edge(scratch_arena_b);
      const auto normal =
          math::cross(mesh.vertices[0], mesh.vertices[1]).normalized();
      auto shader = [&](const math::Vector &p, Fragment &f, float) {
        starts += f.v0 == 0.0f;
        ends += f.v0 == 1.0f;
        geometry_error =
            fold_worst(geometry_error,
                       fabsf(math::dot(normal, p)) * math::ROWS_PER_RADIAN<H>);
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
        Spacetime::geometry_counters.column_generated_samples, starts, ends);
    HS_EXPECT_LE(geometry_error, 0.2f);
    HS_EXPECT_GT(low, 0.48f);
    HS_EXPECT_LT(high, 0.60f);
    HS_EXPECT_LE(sqrtf(variance) / mean, 0.05f);
    HS_EXPECT_EQ(starts, 1u);
    HS_EXPECT_EQ(ends, 1u);
    HS_EXPECT_LT(Spacetime::geometry_counters.column_generated_samples, 128u);
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
  Spacetime::set_flat_surface(false);
  for (int scene = 0; scene < 3; ++scene) {
    const float column = scene == 0 ? 180.0f : scene == 1 ? 50.0f : 96.0f;
    const float motion = scene == 2 ? 1.2f : 0.01f;
    hs_test::StubEffect effect(W, H);
    effect.set_clip(0, H, 36, 72);
    Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipeline{
        Filter::Screen::AntiAlias<W, H>()};
    size_t calls = 0;
    Spacetime::geometry_counters = {};
    {
      Canvas canvas(effect);
      ScratchScope a(scratch_arena_a), b(scratch_arena_b);
      TimedMeshFixture mesh{
          {math::pixel_to_vector<W, H>(column, 60.0f),
           math::pixel_to_vector<W, H>(column + 2.0f, 65.0f),
           math::pixel_to_vector<W, H>(column + 4.0f, 70.0f)}};
      auto edges = chain_edges(scratch_arena_b);
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
    HS_EXPECT_EQ(Spacetime::geometry_counters.column_mesh_culls, 1u);
    HS_EXPECT_EQ(Spacetime::geometry_counters.column_mesh_rejected,
                 scene == 0 ? 1u : 0u);
    HS_EXPECT_EQ(Spacetime::geometry_counters.column_draw_setups,
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
  Spacetime::set_flat_surface(false);
  for (float motion : {0.001f, 0.8f}) {
    hs_test::StubEffect effect(W, H);
    Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipeline;
    Canvas canvas(effect);
    ScratchScope a(scratch_arena_a), b(scratch_arena_b);
    TimedMeshFixture mesh = chain_mesh<W, H>();
    auto edges = chain_edges(scratch_arena_b);
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
    Spacetime::geometry_counters = {};
    Plot::Mesh::draw_timed<W, H>(pipeline, canvas, mesh, edges, shader,
                                 {at, {}}, motion);
    const auto &counts = Spacetime::geometry_counters;
    HS_EXPECT_GT(counts.column_samples, 0u);
    HS_EXPECT_EQ(shaders, counts.column_shader_calls);
    HS_EXPECT_EQ(shaders, counts.column_samples);
    HS_EXPECT_EQ(counts.column_projections, counts.column_samples);
    if (motion < 0.01f) {
      HS_EXPECT_EQ(counts.column_single_samples, counts.column_samples);
      HS_EXPECT_EQ(counts.column_tap_samples, 0u);
    } else {
      HS_EXPECT_GT(counts.column_single_samples, 0u);
      HS_EXPECT_GT(counts.column_tap_samples, 0u);
      HS_EXPECT_EQ(counts.column_single_samples + counts.column_tap_samples,
                   counts.column_samples);
    }
  }
}

/** @brief Stationary timed meshes skip trajectory caches and column preparation. */
inline void test_stationary_column_mesh_skips_geometry_cache() {
  constexpr int W = 288, H = 144;
  Spacetime::set_flat_surface(false);
  hs_test::StubEffect effect(W, H);
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipeline;
  Canvas canvas(effect);
  ScratchScope a(scratch_arena_a), b(scratch_arena_b);
  TimedMeshFixture mesh = chain_mesh<W, H>();
  auto edges = chain_edges(scratch_arena_b);
  uint32_t callbacks = 0, shaders = 0;
  auto at = [&](size_t i, float) {
    ++callbacks;
    return mesh.vertices[i];
  };
  auto shader = [&](const math::Vector &, Fragment &f, float t) {
    ++shaders;
    f.color = Color4(Pixel(65535, 32768, 16384), 0.5f * t);
  };
  Spacetime::geometry_counters = {};
  Plot::Mesh::draw_timed<W, H>(pipeline, canvas, mesh, edges, shader, {at, {}},
                               0.0f);
  const auto &counts = Spacetime::geometry_counters;
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
  Spacetime::set_flat_surface(false);
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
      geometry.path = Spacetime::PathBound::rotation(motion);
      geometry.chained = chained;
      Spacetime::geometry_counters = {};
      Plot::Mesh::draw_timed<W, H>(pipeline, canvas, ring, edges, shader,
                                   geometry, motion);
      HS_EXPECT_EQ(Spacetime::geometry_counters.column_draw_setups, 1u);
      samples[chained] = Spacetime::geometry_counters.column_samples;
    }
    HS_EXPECT_GT(samples[1], 0u);
    HS_EXPECT_EQ(samples[0] - samples[1], uint32_t(N));
  }
}

/** @brief A still edge across an arm seam splats each tap column in the
 *         color shaded at that column's time, on both sides of the seam. */
inline void test_column_walk_seam_taps_shade_at_their_times() {
  constexpr int W = 288, H = 144;
  const ScopedTemporalSampling temporal_sampling(true);
  const auto shade = [](float time) { return 10000.0f + 20000.0f * time; };
  for (float seam : {0.0f, W / 2.0f}) {
    HS_CONTEXT("seam", static_cast<int>(seam));
    hs_test::StubEffect effect(W, H);
    Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipeline;
    {
      Canvas canvas(effect);
      ScratchScope a(scratch_arena_a), b(scratch_arena_b);
      auto edges = single_edge(scratch_arena_b);
      const TimedMeshFixture mesh{
          {math::pixel_to_vector<W, H>(seam - 9.3f, 60.0f),
           math::pixel_to_vector<W, H>(seam + 8.6f, 71.0f), math::X_AXIS}};
      const auto shader = [&](const math::Vector &, Fragment &f, float time) {
        f.color = Color4(
            Pixel(static_cast<uint16_t>(shade(time) + 0.5f), 2000, 3000), 1.0f);
      };
      Plot::Mesh::draw_timed<W, H>(pipeline, canvas, mesh, edges, shader);
    }
    effect.advance_display();
    for (int dx = -2; dx <= 1; ++dx) {
      const int x = math::fast_wrap(static_cast<int>(seam) + dx, W);
      HS_CONTEXT("column", x);
      int column_checked = 0;
      const float expected =
          shade(Spacetime::SamplingSurface<W>::time(static_cast<float>(x))) /
          2000.0f;
      for (int y = 0; y < H; ++y) {
        const Pixel pixel = effect.get_pixel(x, y);
        if (pixel.g < 400)
          continue;
        HS_EXPECT_NEAR(float(pixel.r) / float(pixel.g), expected,
                       0.01f * expected);
        ++column_checked;
      }
      HS_EXPECT_GT(column_checked, 0);
    }
  }
}

/** @brief Refined basis normalization stays unit across its input range. */
inline void test_column_basis_normalization() {
  for (int i = 0; i <= 1024; ++i) {
    const float scale = 0.5f + float(i) / 1024.0f;
    const auto value = math::Vector(0.36f, 0.48f, 0.8f) * scale;
    const auto unit = Plot::EdgeWorldsheet::normalize(value);
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
  const ScopedTemporalSampling temporal_sampling(true);
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
        auto edges = single_edge(scratch_arena_b);
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
}

/** @brief Shrinking edges keep their terminal sample in the column walk. */
inline void test_column_walk_terminal_on_shrinking_edge() {
  constexpr int W = 288, H = 144;
  const ScopedTemporalSampling temporal_sampling(true);
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
      auto edges = single_edge(scratch_arena_b);
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
}

/** @brief Directional spacing certificates preserve the full predicate. */
inline void test_column_directional_spacing_certificate() {
  constexpr int W = 288, H = 144;
  using Frame = Plot::EdgeWorldsheet;
  hs::random().seed(20261027u);
  Plot::g_column_lattice_checks = Plot::g_column_lattice_certificates = 0;
  for (int trial = 0; trial < 6000; ++trial) {
    const auto a = hs_test::rand_unit();
    const auto axis = math::cross(a, hs_test::rand_unit()).normalized();
    const float scale[] = {0.00001f, 0.002f, 0.02f, 0.2f, 1.0f, 2.8f};
    const float length = scale[trial % 6] * hs::rand_f(0.75f, 1.1f);
    const auto b = math::rotate(a, math::make_rotation(axis, length));
    const auto basis = Frame::arc(a, b);
    const float pieces = Frame::pieces_of<W, H>(basis);
    for (float factor : {0.25f, 0.5f, 0.9999f, 1.0f, 1.0001f, 2.0f, 4.0f}) {
      HS_CONTEXT("spacing trial", trial);
      HS_EXPECT_EQ(
          (Frame::needs_local_lattice<W, H>(basis, pieces * factor)),
          (Frame::needs_local_lattice<W, H, false>(basis, pieces * factor)));
    }
    int low = 1, high = 65536;
    while (low < high) {
      const int middle = (low + high) / 2;
      if (Frame::needs_local_lattice<W, H, false>(basis, float(middle)))
        low = middle + 1;
      else
        high = middle;
    }
    for (int count : {std::max(1, low - 1), low, low + 1})
      for (float rounding : {-0.4999f, 0.0f, 0.4999f}) {
        const float near = float(count) + rounding;
        HS_CONTEXT("critical spacing trial", trial);
        HS_EXPECT_EQ((Frame::needs_local_lattice<W, H>(basis, near)),
                     (Frame::needs_local_lattice<W, H, false>(basis, near)));
      }
  }
  std::printf("  [directional spacing] certificates %u / %u\n",
              Plot::g_column_lattice_certificates,
              Plot::g_column_lattice_checks);
  HS_EXPECT_GT(Plot::g_column_lattice_certificates, 100u);
  HS_EXPECT_GT(Plot::g_column_lattice_checks, 10000u);
}

/** @brief An edge's steady lattice is its walk lattice when the count holds
 *         at both frame ends, and lattice points lie on the arc where the
 *         walk places them. */
inline void test_column_steady_lattice() {
  constexpr int W = 288, H = 144;
  using Frame = Plot::EdgeWorldsheet;
  hs::random().seed(20261023u);
  int steady = 0, unsteady = 0, beyond = 0, transported_trials = 0;
  for (int trial = 0; trial < 3000; ++trial) {
    const float column = hs::rand_f(0.0f, float(W));
    const float row = hs::rand_f(8.0f, 136.0f);
    const math::Vector a = math::pixel_to_vector<W, H>(column, row);
    const float heading = hs::rand_f(0.0f, math::TWO_PI_F);
    const float length = hs::rand_f(0.002f, 0.12f);
    const math::Vector west = math::cross(math::Y_AXIS, a).normalized();
    const math::Vector north = math::cross(a, west);
    const math::Vector b = math::rotate(
        a, math::make_rotation(west * sinf(heading) - north * cosf(heading),
                               length));
    const auto motion = math::make_rotation(
        math::Vector(0.3f, 1.0f, 0.2f).normalized(), hs::rand_f(-0.03f, 0.03f));
    const math::Vector starts[2] = {math::rotate(a, motion),
                                    math::rotate(b, motion)};
    const math::Vector ends[2] = {a, b};
    Plot::GeodesicEdgeSampler arcs[2];
    const int got = Frame::steady_lattice<W, H>(starts, ends, 4, arcs);
    const math::RotationMatrix end_to_start(motion);
    Plot::GeodesicEdgeSampler carried[2];
    const int transported =
        Frame::steady_lattice<W, H>(starts, ends, 4, carried, &end_to_start);
    if (transported > 0) {
      ++transported_trials;
      const math::Vector expected =
          math::rotate(Frame::lattice_point(arcs[1], 0.37f), motion);
      const math::Vector actual = Frame::lattice_point(carried[0], 0.37f);
      HS_EXPECT_EQ(transported,
                   Frame::lattice(Frame::pieces_of<W, H>(carried[0])));
      HS_EXPECT_LT((actual - expected).magnitude(), 0.0002f);
    }
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
  HS_EXPECT_GT(transported_trials, 1500);
  HS_EXPECT_GT(steady, 100);
  HS_EXPECT_GT(unsteady, 10);
  HS_EXPECT_GT(beyond, 100);
}

/** @brief Carried steady-edge samples match their point-route rendering. */
inline void test_plot_steady_edge_matches_point_route() {
  constexpr int W = 288, H = 144;
  using Frame = Plot::EdgeWorldsheet;
  for (const float column : {0.4f, 71.2f, 143.4f, 215.3f}) {
    const math::Vector ends[2] = {
        math::pixel_to_vector<W, H>(column, 45.3f),
        math::pixel_to_vector<W, H>(column + 2.1f, 46.2f)};
    const auto rotation =
        math::make_rotation(math::Vector(0.3f, 1.0f, 0.2f).normalized(), 0.02f);
    const math::RotationMatrix transform(rotation);
    const math::Vector starts[2] = {math::rotate(ends[0], rotation),
                                    math::rotate(ends[1], rotation)};
    auto render = [&](bool carried) {
      hs_test::StubEffect effect(W, H);
      Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> sink;
      {
        Canvas canvas(effect);
        Plot::PointRoute<W, H> route;
        route.prepare(canvas, 0.02f);
        auto shade = [](float, float) {
          return Color4(Pixel(65535, 32768, 16384), 0.3f);
        };
        if (carried) {
          HS_EXPECT_TRUE((Plot::Mesh::plot_steady_edge<W, H>(
              route, sink, canvas, starts, ends, 16, shade, &transform)));
        } else {
          Plot::GeodesicEdgeSampler arcs[2];
          const int count = Frame::steady_lattice<W, H>(starts, ends, 16, arcs);
          HS_EXPECT_GT(count, 1);
          for (int k = 0; k < count; ++k) {
            const auto end =
                k == 0 ? ends[0]
                       : Frame::lattice_point(arcs[1], float(k) / count);
            route.plot(
                sink, canvas, end, [&] { return math::rotate(end, rotation); },
                [&](float t) { return shade(float(k) / count, t); });
          }
        }
      }
      effect.advance_display();
      std::vector<Pixel> pixels;
      capture_frame<W, H>(effect, pixels);
      return pixels;
    };
    const auto reference = render(false), candidate = render(true);
    uint32_t energy = 0;
    for (size_t i = 0; i < reference.size(); ++i) {
      energy += reference[i].r;
      HS_EXPECT_LE(abs(int(reference[i].r) - int(candidate[i].r)), 32);
    }
    HS_EXPECT_GT(energy, 1000u);
  }
}

/** @brief Rigid world poses preserve each shaded source arc coordinate. */
inline void test_column_source_basis_transport() {
  constexpr int W = 288, H = 144;
  Spacetime::set_flat_surface(false);
  hs_test::StubEffect effect(W, H);
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipeline;
  Canvas canvas(effect);
  ScratchScope storage(scratch_arena_a), cache(scratch_arena_b);
  const TimedMeshFixture mesh{{math::pixel_to_vector<W, H>(45.3f, 34.8f),
                               math::pixel_to_vector<W, H>(73.4f, 62.1f),
                               math::X_AXIS}};
  const math::Vector axis = math::Vector(0.2f, 0.7f, -0.4f).normalized();
  const math::Quaternion history[] = {math::make_rotation(axis, 0.31f),
                                      math::make_rotation(axis, 0.39f)};
  Spacetime::PoseTable poses;
  poses.bind(scratch_arena_a, W);
  poses.build(history, 2);
  auto edges = single_edge(scratch_arena_b);
  const auto normal =
      math::cross(mesh.vertices[0], mesh.vertices[1]).normalized();
  const float length = acosf(math::dot(mesh.vertices[0], mesh.vertices[1]));
  uint32_t samples = 0;
  auto shader = [&](const math::Vector &p, Fragment &f, float) {
    ++samples;
    HS_EXPECT_NEAR(math::dot(p, normal), 0.0f, 0.0001f);
    HS_EXPECT_NEAR(math::dot(p, mesh.vertices[0]), cosf(length * f.v0),
                   0.0001f);
    f.color = Color4(Pixel(65535, 32768, 16384), 0.3f);
  };
  Plot::Mesh::draw_column_edges<W, H>(pipeline, canvas, mesh, edges, shader, {},
                                      0.0f, &poses, 0.08f);
  HS_EXPECT_GT(samples, 20u);
}

/** @brief Records ages and projected positions before screen antialiasing. */
struct ColumnAgeProbe : Filter::Is2D {
  /** @brief One projected sample before AA taps. */
  struct Sample {
    math::PixelCoords position; ///< Projected sample center.
    float age;                  ///< Interpolated source age.
  };
  std::vector<Sample> &samples; ///< Borrowed sample log.
  /** @brief Binds the sample log. */
  explicit ColumnAgeProbe(std::vector<Sample> &samples) : samples(samples) {}
  /** @brief Records and forwards one screen sample.
   * @tparam Pass Downstream stage.
   * @param x,y Projected position.
   * @param color,age,alpha Sample registers.
   * @param pass Downstream continuation. */
  template <typename Pass>
  void plot(float x, float y, const Pixel &color, float age, float alpha,
            Pass &&pass) {
    samples.push_back({{x, y}, age});
    pass(x, y, color, age, alpha);
  }
};

/** @brief Edge-uniform color preserves rendered pixels and interpolated ages. */
inline void test_column_uniform_color_preserves_age() {
  constexpr int W = 288, H = 144;
  Spacetime::set_flat_surface(false);
  auto run = []<bool CUT>() {
    std::vector<Pixel> reference;
    std::vector<ColumnAgeProbe::Sample> reference_samples;
    uint32_t reference_calls = 0;
    const math::Vector vertices[] = {math::pixel_to_vector<W, H>(42.3f, 41.4f),
                                     math::pixel_to_vector<W, H>(83.2f, 69.8f)};
    const auto rotation = math::make_rotation(math::Y_AXIS, math::PI_F);
    for (bool uniform : {false, true}) {
      hs_test::StubEffect effect(W, H);
      std::vector<ColumnAgeProbe::Sample> samples;
      auto pipeline = [&] {
        if constexpr (CUT)
          return Pipeline<W, H, Filter::Screen::AntiAlias<W, H>,
                          ColumnAgeProbe>{Filter::Screen::AntiAlias<W, H>(),
                                          ColumnAgeProbe(samples)};
        else
          return Pipeline<W, H, ColumnAgeProbe,
                          Filter::Screen::AntiAlias<W, H>>{
              ColumnAgeProbe(samples), Filter::Screen::AntiAlias<W, H>()};
      }();
      uint32_t calls = 0;
      {
        Canvas canvas(effect);
        ScratchScope a(scratch_arena_a), b(scratch_arena_b);
        const Plot::MeshEdge pair{0, 1};
        Spacetime::Worldlines<> motion;
        motion.bind(scratch_arena_b, 2);
        motion.begin_frame(2, false);
        for (int i = 0; i < 2; ++i)
          motion.capture_end(i, vertices[i]);
        motion.begin_frame(2, true);
        for (int i = 0; i < 2; ++i)
          motion.capture_end(i, vertices[i]);
        Spacetime::PoseTable poses;
        poses.bind(scratch_arena_a, W);
        const math::Quaternion history[] = {rotation, rotation};
        poses.build(history, 2);
        Fragments points;
        points.bind(scratch_arena_a, 2);
        auto prepare = [](Fragments &out, const math::Vector &from,
                          const math::Vector &to, size_t, float span) {
          out.clear();
          Fragment first, last;
          first.pos = from;
          last.pos = to;
          first.age = 3.0f;
          last.age = 7.0f;
          last.v0 = 1.0f;
          last.v1 = span;
          out.push_back(first);
          out.push_back(last);
        };
        auto shader = [&](const math::Vector &p, Fragment &f, float) {
          ++calls;
          if constexpr (!CUT)
            HS_EXPECT_VEC(f.pos, p, 0.0f);
          f.color = Color4(!CUT && p.z > 0.0f ? Pixel(65535, 32768, 16384)
                                              : Pixel(16384, 32768, 65535),
                           0.3f);
        };
        const Plot::EdgeBatch batch{.count = 1,
                                    .prepare = prepare,
                                    .vertices = vertices,
                                    .motion = motion,
                                    .pairs = &pair,
                                    .chained = false,
                                    .rigid_lattice = false,
                                    .uniform_color = uniform};
        PipelineRef erased(pipeline, Spacetime::Timed{});
        Plot::rasterize_column_edges<W, H, {.timed = true}>(
            erased, canvas, points, TimedShaderRef(shader), batch, &poses, 0.0f,
            false, 0.0f);
      }
      effect.advance_display();
      std::vector<Pixel> pixels;
      uint64_t energy = 0, blue_energy = 0;
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x) {
          const auto pixel = effect.get_pixel(x, y);
          pixels.push_back(pixel);
          energy += pixel.r;
          blue_energy += pixel.b;
        }
      HS_EXPECT_GT(energy, uint64_t(10000));
      HS_EXPECT_GT(blue_energy, 2 * energy);
      HS_EXPECT_GT(samples.size(), size_t(20));
      float largest_age = 0.0f;
      for (size_t i = 0; i < samples.size(); ++i) {
        if constexpr (!CUT) {
          const float fraction = float(i) / float(samples.size() - 1);
          HS_EXPECT_NEAR(samples[i].age, 3.0f + 4.0f * fraction, 0.00001f);
        } else {
          HS_EXPECT_GE(samples[i].age, 3.0f);
          HS_EXPECT_LE(samples[i].age, 7.0f);
        }
        largest_age = fmaxf(largest_age, samples[i].age);
      }
      HS_EXPECT_GT(largest_age, 6.0f);
      if (!uniform) {
        reference = std::move(pixels);
        reference_samples = std::move(samples);
        reference_calls = calls;
        continue;
      }
      HS_EXPECT_EQ(calls, 1u);
      HS_EXPECT_GT(reference_calls, 20u);
      HS_EXPECT_EQ(samples.size(), reference_samples.size());
      for (size_t i = 0; i < std::min(samples.size(), reference_samples.size());
           ++i)
        HS_EXPECT_EQ(samples[i].age, reference_samples[i].age);
      HS_EXPECT_EQ(pixels.size(), reference.size());
      for (size_t i = 0; i < reference.size(); ++i) {
        HS_EXPECT_EQ(pixels[i].r, reference[i].r);
        HS_EXPECT_EQ(pixels[i].g, reference[i].g);
        HS_EXPECT_EQ(pixels[i].b, reference[i].b);
      }
    }
  };
  run.template operator()<false>();
  run.template operator()<true>();
}

/** @brief A uniform seam color leaves far-tap material registers interpolated. */
inline void test_column_uniform_color_seam_registers() {
  ScratchScope storage(scratch_arena_a);
  Fragments points;
  points.bind(scratch_arena_a, 2);
  Fragment first, last;
  first.v1 = 10.0f;
  first.v2 = 3.0f;
  first.v3 = 4.0f;
  first.age = 3.0f;
  last.v0 = 1.0f;
  last.v1 = 20.0f;
  last.v2 = 7.0f;
  last.v3 = 8.0f;
  last.age = 7.0f;
  last.size = 3.0f;
  points.push_back(first);
  points.push_back(last);
  Fragment owner;
  owner.pos = math::Vector(0.2f, 0.3f, 0.4f);
  owner.v0 = 101.0f;
  owner.v1 = 102.0f;
  owner.v2 = 103.0f;
  owner.v3 = 104.0f;
  owner.age = 105.0f;
  owner.size = 106.0f;
  owner.color = Color4(Pixel(1000, 2000, 3000), 0.3f);
  const Plot::ColumnTapSample sample{{143.8f, 50.2f}, 0.998f, 0.25f};
  for (const float far : {0.0f, sample.t}) {
    uint32_t calls = 0;
    auto shade = [&](const math::Vector &p, Fragment &f, float time) {
      ++calls;
      HS_EXPECT_VEC(p, owner.pos, 0.0f);
      HS_EXPECT_EQ(time, far);
      f.color = owner.color;
    };
    const TimedShaderRef shader(shade);
    const auto normal = Plot::column_seam_registers<{.timed = true}, false>(
        points, sample, owner, far, shader);
    const auto uniform = Plot::column_seam_registers<{.timed = true}, true>(
        points, sample, owner, far, shader);
    HS_EXPECT_EQ(calls, far == sample.t ? 0u : 1u);
    for (const Fragment *result : {&normal, &uniform}) {
      const bool same_time = far == sample.t;
      HS_EXPECT_VEC(result->pos, same_time ? owner.pos : math::Vector{}, 0.0f);
      HS_EXPECT_EQ(result->v0, same_time ? 101.0f : 0.25f);
      HS_EXPECT_EQ(result->v1, same_time ? 102.0f : 12.5f);
      HS_EXPECT_EQ(result->v2, same_time ? 103.0f : 4.0f);
      HS_EXPECT_EQ(result->v3, same_time ? 104.0f : 5.0f);
      HS_EXPECT_EQ(result->age, same_time ? 105.0f : 4.0f);
      HS_EXPECT_EQ(result->size, same_time ? 106.0f : 1.5f);
      HS_EXPECT_EQ(result->color.color.r, 1000u);
      HS_EXPECT_EQ(result->color.alpha, 0.3f);
    }
  }
}

/** @brief A per-vertex path bound replaces the shared one when the column
 *         walk picks each vertex's interpolant. */
inline void test_column_capture_reads_vertex_paths() {
  constexpr int W = 288, H = 144;
  Spacetime::set_flat_surface(false);
  hs_test::StubEffect effect(W, H);
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipeline;
  Canvas canvas(effect);
  ScratchScope a(scratch_arena_a), b(scratch_arena_b);
  TimedMeshFixture mesh = chain_mesh<W, H>();
  auto edges = chain_edges(scratch_arena_b);
  auto at = [&](size_t i, float t) {
    return math::rotate(mesh.vertices[i],
                        math::make_rotation(math::Y_AXIS, 0.02f * (1.0f - t)));
  };
  auto shader = [&](const math::Vector &, Fragment &f, float) {
    f.color = Color4(Pixel(65535, 32768, 16384), 0.5f);
  };
  for (const bool curved : {false, true}) {
    auto vertex_path = [&](size_t i) {
      return curved && i == 1 ? Spacetime::PathBound{1.0f, 1.0f}
                              : Spacetime::PathBound{};
    };
    Plot::GeometryCallback geometry;
    geometry.evaluate = at;
    geometry.path = Spacetime::PathBound{1.0f, 1.0f};
    geometry.vertex_path = vertex_path;
    Spacetime::geometry_counters = {};
    Plot::Mesh::draw_timed<W, H>(pipeline, canvas, mesh, edges, shader,
                                 geometry, 0.02f);
    HS_EXPECT_EQ(Spacetime::geometry_counters.midpoint_evaluations,
                 curved ? 1u : 0u);
  }
}

/** @brief Each arm pass's walk starts no later than its first lit column and
 * gives up no earlier, with little slack, and most passes nothing lights get no
 * walk. */
inline void test_column_plan_covers_lit_columns() {
  constexpr int W = 288, PASS = W / 2;
  hs::random().seed(20261022u);
  const auto basis = [](const math::Vector &a, const math::Vector &b) {
    const auto span = Plot::make_geodesic_edge_span(a, b);
    return Plot::GeodesicEdgeSampler{
        a, span.have_axis ? math::cross(span.axis, a) : math::Vector{},
        span.total};
  };
  uint32_t lit_passes = 0, empty_passes = 0, skipped = 0, slack_count = 0;
  uint32_t trimmed = 0;
  float slack = 0.0f;
  for (int trial = 0; trial < 3000; ++trial) {
    auto a0 = hs_test::rand_unit();
    if (trial % 4 == 0)
      a0 = math::Vector(a0.x * 0.08f, a0.y < 0.0f ? -1.0f : 1.0f, a0.z * 0.08f)
               .normalized();
    const auto edge_axis = math::cross(a0, hs_test::rand_unit()).normalized();
    const float edge_angle = hs::rand_f(0.01f, 1.5f);
    const auto b0 =
        math::rotate(a0, math::make_rotation(edge_axis, edge_angle));
    const auto axis = hs_test::rand_unit();
    const float motion = hs::rand_f(0.0f, 0.15f);
    const auto a1 = math::rotate(a0, math::make_rotation(axis, motion));
    const auto b1 = math::rotate(
        b0, math::make_rotation(axis, motion * hs::rand_f(0.5f, 1.5f)));
    Plot::EdgeWorldsheet frame;
    frame.bases[0] = basis(a0, b0);
    frame.bases[1] = basis(a1, b1);
    frame.ends[0] = b0;
    frame.ends[1] = b1;
    frame.deviation = 0.0f;
    bool lit[W];
    for (int x = 0; x < W; ++x) {
      const float t = Spacetime::SamplingSurface<W>::time(float(x));
      const float a = math::vector_to_column<W>(math::nlerp_unit(a0, a1, t));
      const float span = math::centered_remainder(
          math::vector_to_column<W>(math::nlerp_unit(b0, b1, t)) - a, float(W));
      const float low = span < 0.0f ? a + span : a;
      lit[x] = math::centered_remainder(float(x) - low, float(W)) >= -0.5f &&
               math::centered_remainder(float(x) - low, float(W)) <=
                   fabsf(span) + 0.5f;
    }
    for (const bool bounded : {false, true}) {
      for (const int direction : {1, -1}) {
        Plot::detail::ColumnSegment segments[2];
        const int count = frame.plan<W>(direction, segments, bounded);
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
          const Plot::detail::ColumnSegment *walk = nullptr;
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
          if (!bounded)
            HS_EXPECT_EQ(start + walk->count, PASS);
          else {
            trimmed += start + walk->count < PASS;
            for (int x = pass * PASS; x < (pass + 1) * PASS; ++x)
              if (lit[x])
                HS_EXPECT_LT(offset(x), start + walk->count);
          }
          if (walk->deadline < PASS) {
            slack += float(first - start);
            ++slack_count;
          }
        }
      }
    }
  }
  const float mean = slack / float(std::max(slack_count, 1u));
  std::printf(
      "  [column plan] lit passes %u (with deadline %u), empty %u skipped "
      "%u, mean leading slack %.2f\n",
      lit_passes, slack_count, empty_passes, skipped, mean);
  HS_EXPECT_GT(trimmed, 100u);
  HS_EXPECT_GT(lit_passes, 2000u);
  HS_EXPECT_GT(slack_count, lit_passes * 3 / 4);
  HS_EXPECT_GT(skipped, empty_passes * 9 / 10);
  HS_EXPECT(mean < 2.5f, "walks start close to their lit runs");
}
