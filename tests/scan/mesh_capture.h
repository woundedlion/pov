/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @file mesh_capture.h
 * @brief Mesh deformation capture and column-window replay tests. */

/**
 * @brief Renders a mesh draw under a held column basis into a fresh stub
 *        frame, optionally clipped to a column band.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 * @param clip_begin First column the clip admits.
 * @param clip_end One past the last column the clip admits.
 * @param draw Callable `(Pipeline<W, H> &, Canvas &)`.
 * @return Displayed frame.
 */
template <int W, int H, typename DrawFn>
inline std::vector<Pixel> render_window_frame(int clip_begin, int clip_end,
                                              DrawFn &&draw) {
  hs_test::StubEffect fx(W, H);
  fx.set_clip_x(clip_begin, clip_end);
  Pipeline<W, H> pipeline;
  {
    Canvas canvas(fx);
    auto basis = std::make_unique<Scan::ColumnBasisTable<W, H>>();
    const Scan::ColumnBasis time =
        basis->prepare(math::Quaternion(), math::Quaternion(), true);
    Scan::ColumnBasisScope column_scope(canvas, time);
    draw(pipeline, canvas);
  }
  fx.advance_display();
  std::vector<Pixel> out;
  capture_frame<W, H>(fx, out);
  return out;
}

/** @brief Capture fractions keep the branch nearest each window centre, across
 *  the arm boundary too. */
inline void test_mesh_capture_unwrap() {
  HS_EXPECT_NEAR(Scan::unwrap_capture(0.51f, 0.25f), 0.51f, 1e-6f);
  HS_EXPECT_NEAR(Scan::unwrap_capture(0.49f, 0.75f), 0.49f, 1e-6f);
  HS_EXPECT_NEAR(Scan::unwrap_capture(0.99f, 0.25f), -0.01f, 1e-6f);
  HS_EXPECT_NEAR(Scan::unwrap_capture(0.01f, 0.75f), 1.01f, 1e-6f);
}

/** @brief Window interval rejection and run clipping match column membership. */
inline void test_capture_window_runs() {
  constexpr int W = 16;
  for (int begin = 0; begin < W / 2; ++begin)
    for (int end = begin + 1; end <= W / 2; ++end) {
      const Spacetime::TimeWindow window{begin, end, true};
      for (int a = -W; a < 2 * W; ++a)
        for (int length = 0; length <= W; ++length) {
          const std::array<SDF::Interval, 1> intervals{
              SDF::Interval{float(a), float(a + length)}};
          bool expected = false;
          for (int k = 0; k < std::max(length, 1); ++k) {
            const int x = ((a + k) % W + W) % W;
            expected |= x % (W / 2) >= begin && x % (W / 2) < end;
          }
          HS_EXPECT_EQ(Scan::intervals_meet_window<W>(intervals, window),
                       expected);
        }
      for (int a = 0; a < W; ++a)
        for (int b = a + 1; b <= W; ++b) {
          Scan::ColumnRun runs[3];
          runs[0] = {a, b};
          const size_t count = Scan::window_runs<W, 3>(runs, 1, window);
          for (int x = 0; x < W; ++x) {
            int visits = 0;
            for (size_t i = 0; i < count; ++i)
              visits += x >= runs[i].first && x < runs[i].second;
            HS_EXPECT_EQ(visits, int(x >= a && x < b && x % (W / 2) >= begin &&
                                     x % (W / 2) < end));
          }
        }
    }
}

/** @brief Whole-face window rejection retains both arms and inactive windows. */
inline void test_capture_window_faces() {
  constexpr int W = 96, H = 48;
  const ClipRegion clip{0, H, 0, W, 0, W, H};
  const uint16_t indices[] = {0, 1, 2};
  ScratchScope guard(scratch_arena_a);
  auto *storage = Scan::new_face_storage(scratch_arena_a);
  for (float arm : {0.0f, math::PI_F}) {
    math::Vector vertices[3];
    for (int i = 0; i < 3; ++i) {
      const float theta = arm + 2.2f + 0.1f * i;
      const float phi = i == 1 ? 1.3f : 1.0f;
      vertices[i] = {sinf(phi) * cosf(theta), cosf(phi),
                     sinf(phi) * sinf(theta)};
    }
    const SDF::Face face(vertices, indices, *storage, math::LatitudeGeometry(H),
                         H, &clip, nullptr, math::TWO_PI_F / W);
    const auto &shape = face.prepared();
    HS_EXPECT_LE(shape.y_min, shape.y_max);
    HS_EXPECT_FALSE(Scan::face_meets_window(
        shape, Spacetime::TimeWindow{0, W / 4, true}, clip));
    HS_EXPECT_TRUE(Scan::face_meets_window(
        shape, Spacetime::TimeWindow{W / 4, W / 2, true}, clip));
    HS_EXPECT_TRUE(Scan::face_meets_window(
        shape, Spacetime::TimeWindow{0, W / 4, false}, clip));
  }
}

/**
 * @brief Vertex arc bounds cover an independent angular reference, and
 *        held vertices retain zero motion flags and reach.
 */
inline void test_vertex_arc_and_reach_flags() {
  const double U = std::numeric_limits<float>::epsilon() / 2.0;
  // Dot/norm roundoff and one Newton-refined reciprocal square root.
  const double COSINE_ROUNDING = 16 * U / (1 - 16 * U);
  const double ROUNDING_ARC = std::acos(1 - COSINE_ROUNDING);
  const double PIXELS_PER_RADIAN = 288 / (2 * std::acos(-1.0));
  HS_EXPECT_LT(ROUNDING_ARC * PIXELS_PER_RADIAN, Scan::HELD_VERTEX_PIXELS);
  for (const float scale : {1.0f, 1.5f, 0.65f, 0.35f}) {
    for (const float angle :
         {0.0f, 0.01f, 0.3f, 0.9f, 1.0f, 1.2f, 2.0f, 3.0f}) {
      HS_CONTEXT("angle milliradians", static_cast<int>(angle * 1000.0f));
      const math::Vector a(2 * scale, 0, 0),
          b(scale * cosf(angle), scale * sinf(angle), 0);
      const double ANGLE = hs_test::small_angle_between(a, b);
      const float arc = Scan::vertex_arc(a, b);
      HS_EXPECT_GE(arc, ANGLE - 1e-6);
      const double UPPER =
          ANGLE == 0 ? ROUNDING_ARC * 1.015 + 1e-6 : ANGLE * 1.015 + 1e-6;
      HS_EXPECT_LE(arc, UPPER);
    }
  }
  const math::Vector still_previous[] = {
      {2, 0, 0}, {3, 0, 0}, {1.3f, 0, 0}, {0.7f, 0, 0}};
  math::Vector still_current[4];
  for (int i = 0; i < 4; ++i) {
    still_current[i] = still_previous[i] * 0.5f;
    HS_EXPECT_EQ(
        hs_test::small_angle_between(still_previous[i], still_current[i]), 0.0);
    HS_EXPECT_LT(Scan::vertex_arc(still_previous[i], still_current[i]) *
                     PIXELS_PER_RADIAN,
                 Scan::HELD_VERTEX_PIXELS);
  }
  uint8_t still_motion[4] = {255, 255, 255, 255};
  Scan::flag_moving_vertices<288>(still_previous, still_current, 4,
                                  still_motion);
  uint8_t still_reach[4] = {255, 255, 255, 255};
  Scan::flag_vertex_reaches<288>(still_previous, still_current, 4, still_reach);
  for (int i = 0; i < 4; ++i) {
    HS_EXPECT_EQ(still_motion[i], 0);
    HS_EXPECT_EQ(still_reach[i], 0);
    HS_EXPECT_EQ(Scan::vertex_reach(still_reach[i]), 0.0f);
  }
  // A reach level stands for at least the reach, and at most a third of an
  // octave more.
  for (const float reach : {1e-6f, 1e-3f, 0.003f, 0.02f, 0.4f, 1.0f}) {
    HS_CONTEXT("reach microradians", static_cast<int>(reach * 1000000.0f));
    const uint8_t flags = Scan::vertex_flags(reach);
    HS_EXPECT_EQ(flags & Scan::VERTEX_MOVES, Scan::VERTEX_MOVES);
    HS_EXPECT_GE(Scan::vertex_reach(flags), reach);
    HS_EXPECT_LE(Scan::vertex_reach(flags),
                 fmaxf(1.26f * reach, Scan::VERTEX_REACH_UNIT));
  }
  HS_EXPECT_EQ(Scan::vertex_reach(0), 0.0f);
  // A reach past the levels stands as unbounded.
  HS_EXPECT_GE(Scan::vertex_reach(Scan::vertex_flags(2.0f)), 1e30f);
  // A face reaches the pole by its fastest vertex from its nearest one.
  {
    const auto at = [](float horizontal, float azimuth) {
      return math::Vector(horizontal * cosf(azimuth),
                          sqrtf(1.0f - horizontal * horizontal),
                          horizontal * sinf(azimuth));
    };
    const math::Vector vertices[3] = {at(0.01f, 0.3f), at(0.3f, 0.6f),
                                      at(0.3f, 0.9f)};
    const float angles[3] = {0.3f / math::PI_F, 0.6f / math::PI_F,
                             0.9f / math::PI_F};
    uint8_t flags[3] = {Scan::vertex_flags(0.001f), Scan::vertex_flags(0.05f),
                        Scan::vertex_flags(0.05f)};
    const auto interval = [&](float sweep) {
      float reach = 0.0f, nearest = 1.0f;
      for (int k = 0; k < 3; ++k) {
        reach = fmaxf(reach, Scan::vertex_reach(flags[k]));
        nearest = fminf(nearest, sqrtf(vertices[k].x * vertices[k].x +
                                       vertices[k].z * vertices[k].z));
      }
      return Scan::face_capture_interval(angles, 3, reach, nearest, sweep);
    };
    HS_EXPECT_LT(interval(0.0f).lo, 0.0f);
    flags[1] = flags[2] = Scan::vertex_flags(0.001f);
    HS_EXPECT_GE(interval(0.0f).lo, 0.0f);
    // The basis's own sweep counts toward the reach.
    HS_EXPECT_LT(interval(0.01f).lo, 0.0f);
  }
  const math::Vector placed = math::Vector(0.1f, 1.0f, 0.0f).normalized();
  const math::Vector current[2] = {placed, placed};
  const math::Vector previous[2] = {
      math::rotate(placed, math::make_rotation(math::Vector(0, 0, 1), 0.2f)),
      placed};
  // The second vertex, picked but still, holds.
  uint8_t moving[2] = {Scan::VERTEX_MOVES, Scan::VERTEX_MOVES};
  Scan::flag_vertex_reaches<288>(previous, current, 2, moving);
  HS_EXPECT_GE(Scan::vertex_reach(moving[0]), 0.2f);
  HS_EXPECT_LE(Scan::vertex_reach(moving[0]), 0.26f);
  HS_EXPECT_EQ(moving[1], 0);
}

/** @brief Captures wrap at the arm seam and follow the equator. */
inline void test_mesh_capture_vertex() {
  constexpr int W = 288, H = 144;
  auto basis = std::make_unique<Scan::ColumnBasisTable<W, H>>();
  const auto time =
      basis->prepare(math::Quaternion(), math::Quaternion(), true);
  HS_EXPECT_EQ(
      Scan::capture_fraction(math::Vector(1, 0, -1e-9f), time, 0.0f, W), 0.0f);
  const float start = 0.4f * math::PI_F, step = 0.04f;
  int calls = 0;
  const auto vertex = [&](float u) {
    ++calls;
    return math::Vector(cosf(start + step * u), 0, sinf(start + step * u));
  };
  const uint16_t end = Scan::vertex_column(vertex(1.0f), time, 1.0f, W),
                 begin = Scan::vertex_column(vertex(0.0f), time, 0.0f, W);
  calls = 0;
  const auto captured = Scan::capture_vertex(begin, end, 0.25f, W, vertex);
  const float solved = start / (math::PI_F - step);
  const auto expected = vertex(solved);
  HS_EXPECT_LT(math::distance_squared(captured, expected), 1e-8f);
  HS_EXPECT_EQ(calls, 2);
  HS_EXPECT_NEAR(math::dot(captured, captured), 1.0f, 1e-6f);
  // A fast vertex turning at a steady rate is captured where its column
  // meets the column's time.
  const float fast = 0.4f;
  const auto swing = [&](float u) {
    return math::Vector(cosf(start + fast * u), 0, sinf(start + fast * u));
  };
  const auto fast_capture = Scan::capture_vertex(
      Scan::vertex_column(swing(0.0f), time, 0.0f, W),
      Scan::vertex_column(swing(1.0f), time, 1.0f, W), 0.25f, W, swing);
  HS_EXPECT_LT(
      math::distance_squared(fast_capture, swing(start / (math::PI_F - fast))),
      1e-7f);
}

/**
 * @brief Capture-branch centre of a vertex in a scan of one arm.
 * @param end The vertex at the current endpoint.
 * @param arm The arm scanned, 0 or 1, or -1 for the pass over both.
 * @return Arm fraction the capture's branch lies within a quarter frame of:
 *         the vertex's own time on its arm or in the pass over both, else its
 *         time continued across the arm boundary nearer it.
 */
inline float capture_centre(const math::Vector &end, int arm) {
  const float turns = atan2f(end.z, end.x) / math::PI_F;
  const bool early = turns - floorf(turns) < 0.5f;
  if (arm >= 0 && (turns < 0.0f) != (arm == 1))
    return early ? 1.25f : -0.25f;
  return early ? 0.25f : 0.75f;
}

/**
 * @brief Builds a one-face mesh from vertex positions.
 * @param mesh Receives the face.
 * @param count Number of vertices.
 * @param position Callable (size_t i) giving vertex i.
 */
template <typename PositionFn>
inline void build_single_face(MeshState &mesh, size_t count,
                              PositionFn &&position) {
  mesh.vertices.bind(scratch_arena_a, count);
  mesh.face_counts.bind(scratch_arena_a, 1);
  mesh.face_offsets.bind(scratch_arena_a, 1);
  mesh.faces.bind(scratch_arena_a, count);
  for (uint16_t i = 0; i < count; ++i) {
    mesh.vertices.push_back(position(i));
    mesh.faces.push_back(i);
  }
  mesh.face_counts.push_back(static_cast<uint8_t>(count));
  mesh.face_offsets.push_back(0);
}

/** @brief Faces clear of the arm boundaries scan once, unwindowed, each
 *  vertex captured once on its own side and shared by the faces around it. */
inline void test_mesh_capture_reuse() {
  constexpr int W = 288, H = 144;
  for (float fraction : {0.1f, 0.5f, 0.9f}) {
    ScratchScope guard(scratch_arena_a);
    MeshState mesh;
    const auto position = [&](size_t i, float u) {
      const float angle =
          (fraction + 0.03f * (float(i) - 1.5f)) * math::PI_F + 0.02f * (u - 1);
      return math::Vector(cosf(angle), i % 2 ? -0.1f : 0.1f, sinf(angle))
          .normalized();
    };
    mesh.vertices.bind(scratch_arena_a, 4);
    mesh.face_counts.bind(scratch_arena_a, 2);
    mesh.face_offsets.bind(scratch_arena_a, 2);
    mesh.faces.bind(scratch_arena_a, 6);
    for (uint16_t i = 0; i < 4; ++i)
      mesh.vertices.push_back(position(i, 1));
    for (uint16_t i : {0, 1, 2, 2, 1, 3})
      mesh.faces.push_back(i);
    for (uint16_t f = 0; f < 2; ++f) {
      mesh.face_counts.push_back(3);
      mesh.face_offsets.push_back(static_cast<uint16_t>(3 * f));
    }
    const uint8_t moving[4] = {1, 1, 1, 1};
    const uint16_t start = static_cast<uint16_t>((fraction - 0.06f) * (W / 2));
    const uint16_t columns[4] = {start, 20, start, 20};
    uint16_t starts[4]{};
    int calls = 0, scans = 0;
    render_window_frame<W, H>(0, W, [&](Pipeline<W, H> &, Canvas &canvas) {
      for (size_t i = 0; i < 4; ++i)
        starts[i] = Scan::vertex_column(position(i, 0), *canvas.column_basis(),
                                        0.0f, W);
      Scan::deform_mesh<W, H>(
          canvas, mesh, true, columns, moving, starts, scratch_arena_b, 0.02f,
          [&](size_t i, float u) {
            ++calls;
            return position(i, u);
          },
          [&] {
            ++scans;
            HS_EXPECT_FALSE(canvas.time_window().active);
            for (size_t i = 0; i < 4; ++i) {
              const math::Vector end = position(i, 1);
              const auto &time = *canvas.column_basis();
              const auto expected = Scan::capture_vertex(
                  starts[i], Scan::vertex_column(end, time, 1.0f, W),
                  capture_centre(end, 0), W,
                  [&](float u) { return position(i, u); });
              HS_EXPECT_TRUE(mesh.vertices[i] == expected);
            }
          });
    });
    HS_EXPECT_EQ(scans, 1);
    HS_EXPECT_EQ(calls, 4);
  }
}

/**
 * @brief A face across an arm boundary scans once per arm, clipped to that
 *        arm, each vertex captured at that arm's time; a segment clip
 *        holding one arm scans that arm alone.
 */
inline void test_mesh_seam_face_copies() {
  constexpr int W = 288, H = 144;
  constexpr float FRACTIONS[] = {0.96f, 1.005f, 1.04f};
  constexpr float HEIGHTS[] = {0.1f, -0.1f, 0.1f};
  const auto position = [&](size_t i, float u) {
    const float angle = FRACTIONS[i] * math::PI_F + 0.02f * (u - 1);
    return math::Vector(cosf(angle), HEIGHTS[i], sinf(angle)).normalized();
  };
  // One column of render margin reaches past the clip's end.
  for (const int clip_end : {W, W / 2 - 1}) {
    HS_CONTEXT("clip end", clip_end);
    ScratchScope guard(scratch_arena_a);
    MeshState mesh;
    build_single_face(mesh, 3, [&](size_t i) { return position(i, 1); });
    uint8_t moving[3] = {1, 1, 1};
    uint16_t starts[3]{};
    uint16_t columns[2]{};
    int calls = 0;
    std::vector<int> sides;
    render_window_frame<W, H>(
        0, clip_end, [&](Pipeline<W, H> &, Canvas &canvas) {
          for (size_t i = 0; i < 3; ++i)
            starts[i] = Scan::vertex_column(position(i, 0),
                                            *canvas.column_basis(), 0.0f, W);
          HS_EXPECT_EQ(
              Scan::pick_moving_faces<W>(canvas, mesh, moving, 0.03f, columns),
              size_t(1));
          Scan::deform_mesh<W, H>(
              canvas, mesh, true, columns, moving, starts, scratch_arena_b,
              0.02f,
              [&](size_t i, float u) {
                ++calls;
                return position(i, u);
              },
              [&] {
                const auto &window = canvas.time_window();
                HS_EXPECT_TRUE(window.active);
                HS_EXPECT_EQ(window.end - window.begin, W / 2);
                const int side = window.arm(W);
                sides.push_back(side);
                for (size_t i = 0; i < 3; ++i) {
                  const math::Vector end = position(i, 1);
                  const auto &time = *canvas.column_basis();
                  const auto expected = Scan::capture_vertex(
                      starts[i], Scan::vertex_column(end, time, 1.0f, W),
                      capture_centre(end, side), W,
                      [&](float u) { return position(i, u); });
                  HS_EXPECT_TRUE(mesh.vertices[i] == expected);
                }
              });
        });
    if (clip_end == W) {
      HS_EXPECT_TRUE(sides == std::vector<int>({0, 1}));
      HS_EXPECT_EQ(calls, 6);
    } else {
      HS_EXPECT_TRUE(sides == std::vector<int>({0}));
      HS_EXPECT_EQ(calls, 3);
    }
  }
}

/** @brief A wide face's copies, one per arm, each hold the face's one time
 *  branch at that arm's times. */
inline void test_mesh_wide_face_copies() {
  constexpr int W = 288, H = 144;
  constexpr float STEP = 0.01f;
  constexpr float FRACTIONS[] = {0.40f, 0.45f, 0.85f, 1.30f,
                                 1.30f, 0.85f, 0.45f, 0.40f};
  ScratchScope guard(scratch_arena_a);
  MeshState mesh;
  const auto position = [&](size_t i, float u) {
    const float angle = FRACTIONS[i] * math::PI_F + STEP * (u - 1.0f);
    return math::Vector(cosf(angle), i < 4 ? 0.2f : -0.2f, sinf(angle))
        .normalized();
  };
  build_single_face(mesh, 8, [&](size_t i) { return position(i, 1); });
  const uint8_t moving[8] = {1, 1, 1, 1, 1, 1, 1, 1};
  uint16_t starts[8]{};
  const uint16_t columns[2] = {0, W - 1};
  std::array<int, W> owners{};
  int copies = 0;
  render_window_frame<W, H>(0, W, [&](Pipeline<W, H> &, Canvas &canvas) {
    for (size_t i = 0; i < 8; ++i)
      starts[i] =
          Scan::vertex_column(position(i, 0), *canvas.column_basis(), 0.0f, W);
    Scan::deform_mesh<W, H>(
        canvas, mesh, true, columns, moving, starts, scratch_arena_b, STEP,
        position, [&] {
          const auto &window = canvas.time_window();
          const auto *copy = canvas.column_basis()->mesh.face_copy;
          HS_EXPECT_TRUE(copy != nullptr);
          if (!copy)
            return;
          for (int x = window.begin; x < window.end; ++x)
            for (const int column : {x, x + W / 2})
              if (column >= 0 && column < W)
                ++owners[column];
          ++copies;
          const int shift = window.arm(W) == 1 ? -1 : 0;
          HS_EXPECT_EQ(copy->count, size_t(8));
          for (size_t m = 0; m < copy->count; ++m) {
            const size_t i = copy->sources[m];
            const float u = (FRACTIONS[i] + shift - STEP / math::PI_F) /
                            (1.0f - STEP / math::PI_F);
            HS_EXPECT_LT(
                math::distance_squared(copy->vertices[m], position(i, u)),
                1e-8f);
          }
        });
  });
  HS_EXPECT_EQ(copies, 2);
  for (int x = 0; x < W; ++x) {
    HS_CONTEXT("column", x);
    HS_EXPECT_EQ(owners[x], 1);
  }
}

/** @brief Builds adjacent triangular faces at the requested longitude.
 * @param mesh Destination mesh, backed by scratch arena A.
 * @param rotation Longitude rotation in radians. */
inline void build_shared_edge_mesh(MeshState &mesh, float rotation) {
  mesh.vertices.bind(scratch_arena_a, 5);
  for (const auto &v : {math::Vector(1, .3f, -.2f), math::Vector(1, -.2f, -.1f),
                        math::Vector(1, .2f, .1f), math::Vector(1, -.3f, .2f),
                        math::Vector(1, .5f, .3f)})
    mesh.vertices.push_back(math::rotate(
        v.normalized(), math::make_rotation(math::Vector(0, 1, 0), rotation)));
  mesh.face_counts.bind(scratch_arena_a, 3);
  mesh.face_offsets.bind(scratch_arena_a, 3);
  mesh.faces.bind(scratch_arena_a, 9);
  for (int f = 0; f < 3; ++f) {
    mesh.face_counts.push_back(3);
    mesh.face_offsets.push_back(3 * f);
  }
  for (uint16_t i : {0, 1, 2, 2, 1, 3, 0, 2, 4})
    mesh.faces.push_back(i);
}

/**
 * @brief Adjacent moving and held faces use coincident shared-edge endpoints.
 * @param rotation Longitude rotating the fixture across the arm boundaries.
 * @details A face scanned alone is read from its copy.
 */
inline void check_deformed_shared_edges(float rotation) {
  constexpr int W = 288, H = 144;
  ScratchScope guard(scratch_arena_a);
  MeshState mesh;
  build_shared_edge_mesh(mesh, rotation);
  const std::vector<math::Vector> current(mesh.vertices.begin(),
                                          mesh.vertices.end());
  const uint8_t moving[5] = {1, 0, 0, 0, 1};
  std::vector<uint8_t> picked(moving, moving + 5);
  uint16_t columns[6]{};
  const auto vertex = [&](size_t i, float u) {
    return math::rotate(
        current[i], math::make_rotation(math::Vector(0, 1, 0),
                                        (u - 1) * (moving[i] ? .12f : .001f)));
  };
  uint16_t starts[5]{};
  int scans = 0, held = 0, shared = 0, captured = 0;
  render_window_frame<W, H>(0, W, [&](Pipeline<W, H> &, Canvas &canvas) {
    for (size_t i = 0; i < 5; ++i)
      starts[i] =
          Scan::vertex_column(vertex(i, 0), *canvas.column_basis(), 0.0f, W);
    HS_EXPECT_EQ(
        Scan::pick_moving_faces<W>(canvas, mesh, picked.data(), .13f, columns),
        size_t(2));
    Scan::deform_mesh<W, H>(
        canvas, mesh, true, columns, moving, starts, scratch_arena_b, .12f,
        vertex, [&] {
          ++scans;
          const auto &window = canvas.time_window();
          const auto *selected = canvas.column_basis()->mesh.selected_faces;
          const uint8_t drawn = window.active ? Scan::detail::FacePass::ARM_COPY
                                              : Scan::detail::FacePass::SHARED;
          const int arm = window.active ? window.arm(W) : -1;
          const auto *copy = canvas.column_basis()->mesh.face_copy;
          const auto at = [&](size_t i) {
            if (copy)
              for (size_t m = 0; m < copy->count; ++m)
                if (copy->sources[m] == i)
                  return copy->vertices[m];
            return mesh.vertices[i];
          };
          for (size_t f = 0; f < 3; ++f) {
            if (selected[f] != drawn)
              continue;
            held += f == 1;
            for (uint16_t i : {1, 2}) {
              HS_EXPECT_TRUE(at(i) == current[i]);
              ++shared;
            }
          }
          for (size_t f : {size_t(0), size_t(2)}) {
            if (selected[f] != drawn)
              continue;
            const size_t i = f == 0 ? 0 : 4;
            const auto &time = *canvas.column_basis();
            const auto expected = Scan::capture_vertex(
                starts[i], Scan::vertex_column(current[i], time, 1.0f, W),
                capture_centre(current[i], arm), W,
                [&](float u) { return vertex(i, u); });
            // Each branch's capture solves from the stored endpoint, not
            // from the pool's other capture.
            HS_EXPECT_TRUE(at(i) == expected);
            ++captured;
          }
        });
  });
  HS_EXPECT_EQ(held, 1);
  HS_EXPECT_GE(scans, 1);
  HS_EXPECT_GE(shared, 6);
  HS_EXPECT_GE(captured, 2);
}

/** @brief Shared edges coincide across ordinary and arm-boundary faces. */
inline void test_mesh_deformation_shared_edges() {
  for (float rotation :
       {0.0f, .5f * math::PI_F, math::PI_F, 1.5f * math::PI_F, .3f})
    check_deformed_shared_edges(rotation);
}

/** @brief Columns from recorded cap radii match columns measured afresh. */
inline void test_recorded_face_caps_match_measured_columns() {
  constexpr int W = 288, H = 144;
  ScratchScope guard(scratch_arena_a);
  MeshState mesh;
  build_shared_edge_mesh(mesh, 0.0f);
  render_window_frame<W, H>(0, W, [&](Pipeline<W, H> &, Canvas &canvas) {
    uint8_t moving[5] = {1, 1, 1, 1, 1};
    uint16_t recorded[6]{}, measured[6]{}, picked[6]{};
    float radii[3]{};
    HS_EXPECT_EQ(
        Scan::pick_moving_faces<W>(canvas, mesh, moving, .13f, picked, radii),
        size_t(3));
    std::copy_n(picked, 6, recorded);
    std::copy_n(picked, 6, measured);
    Scan::moving_face_columns<W>(canvas, mesh, .02f, recorded, radii);
    Scan::moving_face_columns<W>(canvas, mesh, .02f, measured);
    for (int k = 0; k < 6; ++k)
      HS_EXPECT_EQ(recorded[k], measured[k]);
    for (int f = 0; f < 3; ++f) {
      HS_EXPECT_LT(measured[2 * f + 1], picked[2 * f + 1]);
      HS_EXPECT_GT(radii[f], 0.0f);
    }
  });
}

/** @brief Moving faces that no placement brings into the clip's rows hold. */
inline void test_moving_faces_outside_clip_rows_hold() {
  constexpr int W = 288, H = 144;
  ScratchScope guard(scratch_arena_a);
  MeshState mesh;
  mesh.vertices.bind(scratch_arena_a, 9);
  // North of the clip's last row, just south of it, and far south.
  for (const float y : {.8f, -.02f, -.8f})
    for (const auto &v :
         {math::Vector(1, y, -.1f), math::Vector(1, y - .05f, 0),
          math::Vector(1, y, .1f)})
      mesh.vertices.push_back(v.normalized());
  mesh.face_counts.bind(scratch_arena_a, 3);
  mesh.face_offsets.bind(scratch_arena_a, 3);
  mesh.faces.bind(scratch_arena_a, 9);
  for (uint16_t f = 0; f < 3; ++f) {
    mesh.face_counts.push_back(3);
    mesh.face_offsets.push_back(static_cast<uint16_t>(3 * f));
    for (uint16_t k = 0; k < 3; ++k)
      mesh.faces.push_back(static_cast<uint16_t>(3 * f + k));
  }
  auto pick = [&](int rows, uint8_t *moving, uint16_t *columns) {
    hs_test::StubEffect fx(W, H);
    fx.set_clip(0, rows, 0, W);
    Canvas canvas(fx);
    auto basis = std::make_unique<Scan::ColumnBasisTable<W, H>>();
    const Scan::ColumnBasis time =
        basis->prepare(math::Quaternion(), math::Quaternion(), true);
    Scan::ColumnBasisScope column_scope(canvas, time);
    const size_t picked =
        Scan::pick_moving_faces<W>(canvas, mesh, moving, .02f, columns);
    Scan::moving_face_columns<W>(canvas, mesh, .02f, columns);
    return picked;
  };
  uint8_t moving[9];
  uint16_t columns[6];
  std::fill_n(moving, 9, uint8_t(1));
  HS_EXPECT_EQ(pick(H, moving, columns), size_t(3));
  std::fill_n(moving, 9, uint8_t(1));
  HS_EXPECT_EQ(pick(H / 2, moving, columns), size_t(2));
  HS_EXPECT_NE(columns[0], Scan::HELD_FACE_COLUMN);
  HS_EXPECT_NE(columns[2], Scan::HELD_FACE_COLUMN);
  HS_EXPECT_EQ(columns[4], Scan::HELD_FACE_COLUMN);
  for (int i = 6; i < 9; ++i)
    HS_EXPECT_EQ(moving[i], 0);
}

/** @brief Temporal sampling deforms ordinary motion and holds tiny steps and jumps. */
inline void test_mesh_motion_activation() {
  constexpr int W = 288;
  constexpr float PIXEL = math::TWO_PI_F / W;
  HS_EXPECT_FALSE(Scan::mesh_motion_active<W>(0.0f));
  HS_EXPECT_FALSE(Scan::mesh_motion_active<W>(0.49f * PIXEL));
  HS_EXPECT_TRUE(Scan::mesh_motion_active<W>(0.51f * PIXEL));
  HS_EXPECT_TRUE(Scan::mesh_motion_active<W>(7.9f * PIXEL));
  HS_EXPECT_FALSE(Scan::mesh_motion_active<W>(8.1f * PIXEL));
  Spacetime::set_flat_surface(true);
  HS_EXPECT_FALSE(Scan::mesh_motion_active<W>(3.4f * PIXEL));
  Spacetime::set_flat_surface(false);
}

/// Largest channel difference between a swept pole face and its per-column
/// oracle: under a tenth of an antialiased edge pixel, from the oracle's
/// faces projecting about their own centroids and the sweep's quadratic
/// vertex paths.
inline constexpr int POLAR_SWEEP_BOUND = 6144;

/** @brief Moving polar-face draws and independent column-time references. */
template <typename Position> struct PolarFaceCase {
  static constexpr int W = 288, H = 144;
  ScratchScope guard{scratch_arena_a};
  MeshState mesh;
  size_t count;
  float reach;
  BasisMotion motion;
  Position position;
  std::vector<uint8_t> moving;
  std::vector<uint16_t> starts;
  std::array<uint16_t, 2> columns{};

  PolarFaceCase(size_t count, float reach, const math::Quaternion &to,
                Position position)
      : count(count), reach(reach), motion{{}, to}, position(position),
        moving(count, 1), starts(count) {
    build_single_face(mesh, count, [&](size_t i) { return position(i, 1); });
  }

  template <typename Draw> std::vector<Pixel> render(Draw draw) {
    return render_column_basis<W, H>(
        motion,
        [&](Pipeline<W, H> &pipeline, Canvas &canvas, const auto &) {
          draw(pipeline, canvas);
        },
        {}, true);
  }

  void scan(Pipeline<W, H> &pipeline, Canvas &canvas) {
    float scale = 0.0f;
    const auto setup = [&](size_t, float size) { scale = 1.0f / size; };
    const auto shader = [&](const math::Vector &, Fragment &f) {
      const float depth = hs::clamp(-f.v1 * scale, 0.0f, 1.0f);
      f.color = Color4(
          Pixel(static_cast<uint16_t>(2000 + 60000 * depth), 20000, 30000), 1);
    };
    Scan::Mesh::draw_specialized<W, H>(pipeline, canvas, mesh, shader,
                                       scratch_arena_b, nullptr, setup);
  }

  template <typename Hook> std::vector<Pixel> draw(Hook before_scan) {
    Scan::mesh_capture_counters = {};
    return render([&](Pipeline<W, H> &pipeline, Canvas &canvas) {
      for (size_t i = 0; i < count; ++i)
        starts[i] = Scan::vertex_column(position(i, 0), *canvas.column_basis(),
                                        0.0f, W);
      Scan::pick_moving_faces<W>(canvas, mesh, moving.data(), reach,
                                 columns.data());
      for (size_t i = 0; i < count; ++i)
        if (moving[i])
          moving[i] = Scan::vertex_flags(reach);
      Scan::deform_mesh<W, H>(canvas, mesh, true, columns.data(), moving.data(),
                              starts.data(), scratch_arena_b, reach, position,
                              [&] {
                                before_scan(canvas);
                                scan(pipeline, canvas);
                              });
    });
  }

  std::vector<Pixel> frozen(float fraction) {
    for (size_t i = 0; i < count; ++i)
      mesh.vertices[i] = position(i, fraction);
    return render([&](Pipeline<W, H> &pipeline, Canvas &canvas) {
      scan(pipeline, canvas);
    });
  }

  std::vector<Pixel> reference() {
    return column_time_oracle<W, H>([&](float u) { return frozen(u); });
  }
};

/**
 * @brief A face around a pole draws each column pair as the face stood at
 *        the pair's sweep time.
 * @param points Star points; the face alternates tip and notch vertices.
 * @param notch The notch vertices' polar angle; the tips' is 0.4 radians.
 * @param north Whether the face winds around the north pole.
 * @param turn Camera turn over the frame (radians).
 * @param grow Whether the face grows from a point over the frame, so that
 *        its first build culls and every column pair builds it afresh.
 * @param tilt Arc the face tips off the pole by, about the z axis; a tilted
 *        face reaches the pole (Scan::face_capture_interval()) without
 *        winding around it, so that only the column pairs near it scan.
 * @details Compares the deformed scan with frames of the face held at each
 * column pair's time, every column pair taken from its own time's frame.
 */
inline void check_polar_face_sweep(size_t points, float notch, bool north,
                                   float turn, bool grow = false,
                                   float tilt = 0.0f) {
  constexpr int W = 288;
  constexpr float SPIN = 0.2f, TIP = 0.4f, WOBBLE = 3.0f * math::PI_F / 180.0f;
  HS_CONTEXT("star points", static_cast<long long>(points));
  HS_CONTEXT("north", north ? 1 : 0);
  const size_t count = 2 * points;
  const auto position = [&](size_t i, float u) {
    const float azimuth =
        math::TWO_PI_F * static_cast<float>(i) / count + SPIN * (u - 1.0f);
    const float polar =
        grow ? (i % 2 ? notch : TIP) * u
             : (i % 2 ? notch : TIP) +
                   WOBBLE * (u - 1.0f) * (i % 4 < 2 ? 1.0f : -1.0f);
    const float y = cosf(polar);
    return math::rotate(math::Vector(sinf(polar) * cosf(azimuth),
                                     north ? y : -y,
                                     sinf(polar) * sinf(azimuth)),
                        math::make_rotation(math::Vector(0, 0, 1), tilt));
  };
  const math::Quaternion to =
      math::make_rotation(math::Vector(0.3f, 0.2f, 0.9f).normalized(), turn);
  const float reach = grow ? TIP + SPIN * sinf(TIP) : SPIN * sinf(TIP) + WOBBLE;
  PolarFaceCase test(count, reach, to, position);
  const std::vector<Pixel> drawn = test.draw([](Canvas &canvas) {
    HS_EXPECT_TRUE(canvas.column_basis()->mesh.face_sweep != nullptr);
  });
  const Scan::MeshCaptureCounters work = Scan::mesh_capture_counters;
  HS_EXPECT_EQ(work.polar_faces, size_t(1));
  HS_EXPECT_EQ(work.swept_faces, size_t(1));
  if (tilt != 0.0f) {
    HS_EXPECT_LT(work.sweep_columns, size_t(W / 2));
    HS_EXPECT_GT(work.sweep_columns, size_t(0));
  } else {
    // A face grown from a point culls at the first column.
    HS_EXPECT_EQ(work.sweep_columns, size_t(grow ? W / 2 - 1 : W / 2));
  }
  // Classification measures each moving vertex's motion at the frame's start
  // and middle, and the sweep's paths evaluate it there again.
  HS_EXPECT_EQ(work.evaluations, 4 * count);
  const std::vector<Pixel> expected = test.reference();
  HS_EXPECT_GT(
      expect_frames_within("polar sweep", drawn, expected, POLAR_SWEEP_BOUND),
      tilt != 0.0f ? W : W * 4);
}

/** @brief Faces around or reaching either pole draw each column pair at its
 *  own time: convex, starred, turning convex during the frame, growing from a
 *  point or tipped off the pole, under a still or a turning camera. */
inline void test_mesh_polar_face_sweep() {
  check_polar_face_sweep(6, 0.12f, true, 0.0f);
  check_polar_face_sweep(6, 0.12f, false, 0.02f);
  check_polar_face_sweep(3, 0.05f, true, 0.02f);
  check_polar_face_sweep(4, 0.2f, false, 0.0f);
  check_polar_face_sweep(3, 0.12f, true, 0.0f);
  check_polar_face_sweep(4, 0.2f, true, 0.02f, true);
  check_polar_face_sweep(6, 0.12f, true, 0.02f, false, 0.485f);
}

/**
 * @brief A face reaching a pole, too slow for its single placement to
 *        stray a pixel, draws in each half-arm window's scan, its vertices
 *        captured near the window's time.
 * @param points Star points; the face alternates tip and notch vertices.
 * @param notch The notch vertices' polar angle; the tips' is 0.2 radians.
 * @param north Whether the face winds around the north pole.
 * @param turn Camera turn over the frame (radians).
 * @details Compares the deformed scan and the face held at the endpoint with
 * frames of the face at each column pair's time.
 */
inline void check_polar_face_windows(size_t points, float notch, bool north,
                                     float turn) {
  constexpr float SPIN = 0.05f, TIP = 0.2f, WOBBLE = 0.5f * math::PI_F / 180.0f;
  HS_CONTEXT("star points", static_cast<long long>(points));
  HS_CONTEXT("north", north ? 1 : 0);
  const size_t count = 2 * points;
  const auto position = [&](size_t i, float u) {
    const float azimuth =
        math::TWO_PI_F * static_cast<float>(i) / count + SPIN * (u - 1.0f);
    const float polar = (i % 2 ? notch : TIP) +
                        WOBBLE * (u - 1.0f) * (i % 4 < 2 ? 1.0f : -1.0f);
    const float y = cosf(polar);
    return math::Vector(sinf(polar) * cosf(azimuth), north ? y : -y,
                        sinf(polar) * sinf(azimuth));
  };
  const math::Quaternion to =
      math::make_rotation(math::Vector(0.3f, 0.2f, 0.9f).normalized(), turn);
  const float reach = SPIN * sinf(TIP) + WOBBLE;
  PolarFaceCase test(count, reach, to, position);
  const std::vector<Pixel> drawn = test.draw([](Canvas &) {});
  const Scan::MeshCaptureCounters work = Scan::mesh_capture_counters;
  HS_EXPECT_EQ(work.polar_faces, size_t(1));
  HS_EXPECT_EQ(work.swept_faces, size_t(0));
  HS_EXPECT_EQ(work.face_copies, size_t(0));
  HS_EXPECT_GE(work.scans, size_t(1));
  HS_EXPECT_LE(work.scans, size_t(2));
  // Classification measures each moving vertex's motion at the frame's start
  // and middle.
  HS_EXPECT_EQ(work.evaluations, work.captures + 2 * count);
  HS_EXPECT_LE(work.captures, work.scans * count);
  const std::vector<Pixel> held = test.frozen(1.0f);
  const std::vector<Pixel> expected = test.reference();
  const auto error = [&](const std::vector<Pixel> &frame) {
    double sum = 0.0;
    for (size_t i = 0; i < expected.size(); ++i)
      sum += std::abs(int(frame[i].r) - expected[i].r) +
             std::abs(int(frame[i].g) - expected[i].g) +
             std::abs(int(frame[i].b) - expected[i].b);
    return sum;
  };
  const double drawn_error = error(drawn), held_error = error(held);
  std::printf("  [polar face %zu] error %.0f vs held %.0f (ratio %.3f)\n",
              points, drawn_error, held_error, drawn_error / held_error);
  HS_EXPECT_GT(held_error, 0.0);
  HS_EXPECT_LT(drawn_error, held_error);
}

/** @brief Faces around or reaching either pole draw per half-arm window,
 *  closer to their column times than held at the endpoint: convex or starred,
 *  under a still or a turning camera. */
inline void test_mesh_polar_face_windows() {
  check_polar_face_windows(6, 0.12f, true, 0.0f);
  check_polar_face_windows(6, 0.12f, false, 0.02f);
  check_polar_face_windows(3, 0.05f, true, 0.02f);
  check_polar_face_windows(4, 0.2f, false, 0.0f);
}
