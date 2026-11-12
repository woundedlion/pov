/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @brief Chained edges share the last slot; unchained edges reserve both slots. */
inline void test_column_chunk_capacity_boundary() {
  constexpr int W = 72, H = 36;
  const ScopedTemporalSampling sampling(false);
  for (bool chained : {true, false}) {
    StubEffect effect(W, H);
    Canvas canvas(effect);
    Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipeline;
    Filter::Screen::DirectAntiAliasSink<W, H> sink;
    sink.prepare(canvas);
    Plot::PointRoute<W, H> route;
    route.prepare(canvas, 0.0f);
    ScratchScope storage(scratch_arena_a);
    Plot::Mesh::ColumnChunk<4, float> chunk(scratch_arena_a);
    size_t flushes = 0, stores = 0;
    bool overflow = false;
    const auto at = [&](size_t slot, float) { return chunk.position(slot); };
    const Plot::GeometryCallback geometry{at, {}, true};
    const auto shade = [](float) { return Color4(Pixel(65535, 0, 0), 1.0f); };
    const auto setup = Plot::Mesh::column_draw_setup<W, H>(pipeline, canvas);
    const auto flush = [&] {
      ++flushes;
      if (!overflow)
        chunk.flush<W, H>(
            pipeline, canvas, setup,
            [&](size_t slot) { return chunk.attribute<0>(slot); },
            [](size_t, float, float) {
              return Color4(Pixel(65535, 0, 0), 1.0f);
            },
            geometry, 0.0f);
    };
    const auto store = [&](size_t slot, size_t source) {
      HS_EXPECT_LT(slot, size_t{4});
      if (slot >= 4) {
        overflow = true;
        return;
      }
      if (!chained && source == 3)
        HS_EXPECT_EQ(flushes, size_t{1});
      ++stores;
      chunk.position(slot) =
          math::pixel_to_vector<W, H>(10.0f + float(source), 18.0f);
      chunk.attribute<0>(slot) = float(source);
    };
    chunk.begin_chain(0, flush);
    for (uint16_t edge = 0; edge < 3; ++edge) {
      const uint16_t head = !chained && edge == 2 ? 3 : edge;
      const Plot::Mesh::Edge source{head, uint16_t(head + 1)};
      const math::Vector ends[2] = {
          math::pixel_to_vector<W, H>(10.0f + head, 18.0f),
          math::pixel_to_vector<W, H>(11.0f + head, 18.0f)};
      chunk.plot_edge<W, H>(route, sink, canvas, source, ends, ends, 0, shade,
                            store, flush);
    }
    HS_EXPECT_FALSE(overflow);
    HS_EXPECT_EQ(flushes, chained ? size_t{0} : size_t{1});
    HS_EXPECT_EQ(stores, chained ? size_t{4} : size_t{5});
    flush();
    HS_EXPECT_EQ(flushes, chained ? size_t{1} : size_t{2});
  }
}

/** @brief Cap-ring intervals contain every knot inside the cap. */
inline void test_cap_ring_arc_contains_cap_knots() {
  constexpr int N = 96;
  for (int trial = 0; trial < 500; ++trial) {
    const float theta = math::PI_F * float(trial % 43) / 42.0f;
    const float polar = math::PI_F * float((trial * 17) % 47) / 46.0f;
    const float azimuth = 2.0f * math::PI_F * float((trial * 23) % 53) / 53.0f;
    const float radius = math::PI_F * float((trial * 13) % 59) / 58.0f;
    const float cosine = cosf(radius);
    const auto arc = Plot::cap_ring_arc(cosf(theta), sinf(theta), cosf(polar),
                                        sinf(polar), azimuth, cosine, N, 1e-5f);
    for (int k = 0; k < N; ++k) {
      const float a = 2.0f * math::PI_F * k / N;
      const float product = cosf(theta) * cosf(polar) +
                            sinf(theta) * sinf(polar) * cosf(a - azimuth);
      if (product > cosine) {
        HS_EXPECT_FALSE(arc.empty);
        HS_EXPECT_LT((k - arc.first + N) % N, arc.count);
      }
    }
  }
}

/** @brief Swept bump cap-ring arcs contain every sampled footprint knot. */
inline void test_swept_bump_arc_contains_footprint_knots() {
  constexpr double PI_EXACT = 3.14159265358979323846;
  constexpr int MAX_KNOTS = 288, TIME_STEPS = 128;
  using Point = std::array<double, 3>;
  const auto dot = [](const Point &a, const Point &b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  };
  const auto normalized = [&](Point p) {
    const double length = std::sqrt(dot(p, p));
    for (double &component : p)
      component /= length;
    return p;
  };
  const auto cross = [](const Point &a, const Point &b) {
    return Point{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                 a[0] * b[1] - a[1] * b[0]};
  };
  const auto point = [](const math::Vector &p) {
    return Point{double(p.x), double(p.y), double(p.z)};
  };
  const auto vector = [](const Point &p) {
    return math::Vector(float(p[0]), float(p[1]), float(p[2]));
  };
  uint32_t seed = 0x73c5e21bu;
  const auto unit = [&]() {
    seed = seed * 1664525u + 1013904223u;
    return double(seed >> 8) / 16777216.0;
  };
  const auto direction = [&]() {
    const double x = 2.0 * unit() - 1.0;
    const double y = 2.0 * unit() - 1.0;
    const double z = 2.0 * unit() - 1.0;
    return normalized(Point{x, y, z});
  };
  int hits = 0, fast_hits = 0, antipodal_hits = 0;
  int partial_arcs = 0, empty_arcs = 0, whole_arcs = 0;
  for (int trial = 0; trial < 160; ++trial) {
    const Point from = direction();
    Point tangent = cross(from, direction());
    if (dot(tangent, tangent) < 1e-10)
      tangent = cross(from, Point{0.0, 1.0, 0.0});
    tangent = normalized(tangent);
    double separation = PI_EXACT * unit();
    if (trial == 0)
      separation = 0.0;
    else if (trial < 8)
      separation = PI_EXACT - std::pow(10.0, -trial);
    Point to;
    for (int component = 0; component < 3; ++component)
      to[component] = from[component] * std::cos(separation) +
                      tangent[component] * std::sin(separation);
    Animation::BumpParams params;
    params.previous_center = vector(from);
    params.center = vector(normalized(to));
    params.radius = float(0.02 + 0.9 * unit());
    params.previous_envelope = float(0.1 + 0.9 * unit());
    params.envelope = float(0.1 + 0.9 * unit());
    const auto basis =
        math::make_basis(math::Quaternion(), vector(direction()));
    Animation::BumpPath path{};
    Animation::BumpCover cover{};
    HS_EXPECT_TRUE(Animation::make_bump_path(params, {}, basis.v, path, cover));
    const Point first = point(params.previous_center),
                last = point(params.center);
    const Point axis = point(basis.v), u = point(basis.u), w = point(basis.w);
    const float axial = math::dot(basis.v, cover.centre);
    const float cu = math::dot(basis.u, cover.centre);
    const float cw = math::dot(basis.w, cover.centre);
    const float radial = std::sqrt(cu * cu + cw * cw);
    float azimuth = std::atan2(cw, cu);
    if (azimuth < 0.0f)
      azimuth += math::TWO_PI_F;
    const int n = trial % 3 == 0 ? 48 : trial % 3 == 1 ? 96 : MAX_KNOTS;
    for (int ring = 0; ring < 3; ++ring) {
      const Point endpoint = normalized(ring == 1 ? first : last);
      float theta =
          ring == 0
              ? float(PI_EXACT * unit())
              : float(std::acos(std::clamp(dot(axis, endpoint), -1.0, 1.0)));
      if (ring == 0 && trial % 11 == 0)
        theta = trial % 2 ? 0.0f : math::PI_F;
      const auto arc =
          Plot::cap_ring_arc(std::cos(theta), std::sin(theta), axial, radial,
                             azimuth, path.cos_cover, n, 2e-4f);
      empty_arcs += arc.empty;
      whole_arcs += !arc.empty && arc.count == n;
      partial_arcs += !arc.empty && arc.count < n;
      std::array<Point, MAX_KNOTS> knots;
      for (int k = 0; k < n; ++k) {
        const double a = 2.0 * PI_EXACT * k / n;
        Point p;
        for (int component = 0; component < 3; ++component)
          p[component] =
              axis[component] * std::cos(double(theta)) +
              (u[component] * std::cos(a) + w[component] * std::sin(a)) *
                  std::sin(double(theta));
        knots[k] = normalized(p);
      }
      for (int step = 0; step <= TIME_STEPS; ++step) {
        const double t = double(step) / TIME_STEPS;
        Point centre;
        for (int component = 0; component < 3; ++component)
          centre[component] =
              first[component] + (last[component] - first[component]) * t;
        if (dot(centre, centre) < 1e-20)
          continue;
        centre = normalized(centre);
        const double reach =
            double(params.radius) *
            (params.previous_envelope +
             (double(params.envelope) - params.previous_envelope) * t);
        const double cosine = std::cos(reach);
        for (int k = 0; k < n; ++k) {
          if (dot(knots[k], centre) <= cosine)
            continue;
          ++hits;
          fast_hits += separation > 0.5;
          antipodal_hits += trial > 0 && trial < 8;
          HS_EXPECT_FALSE(arc.empty);
          HS_EXPECT_LT((k - arc.first + n) % n, arc.count);
        }
      }
    }
  }
  HS_EXPECT_GT(hits, 0);
  HS_EXPECT_GT(fast_hits, 0);
  HS_EXPECT_GT(antipodal_hits, 0);
  HS_EXPECT_GT(partial_arcs, 0);
  HS_EXPECT_GT(empty_arcs, 0);
  HS_EXPECT_GT(whole_arcs, 0);
}
