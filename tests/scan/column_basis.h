/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Column-time rotated basis: scans against shapes rebuilt at column time
// ============================================================================

/** @brief Orientation endpoints of one oracle frame. */
struct BasisMotion {
  math::Quaternion previous; /**< Orientation at the previous endpoint. */
  math::Quaternion current;  /**< Orientation at the current endpoint. */
};

/** @brief Display clip of one oracle render. */
struct BasisClip {
  int x0 = 0;      /**< First display column. */
  int x1 = 0;      /**< One past the last display column; 0 = full width. */
  int margin = -1; /**< Render margin, or -1 to keep the default. */
};

/**
 * @brief Renders a draw at the current endpoint through the rotated basis.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 * @param motion Frame endpoints.
 * @param draw Callable `(Pipeline<W, H> &, Canvas &, const Quaternion &)`
 *        drawing the shapes at an orientation.
 * @param clip Display clip.
 * @return Displayed frame.
 */
template <int W, int H, typename DrawFn>
inline std::vector<Pixel> render_column_basis(const BasisMotion &motion,
                                              DrawFn &&draw,
                                              const BasisClip &clip = {}) {
  hs_test::StubEffect fx(W, H);
  if (clip.x1 > 0)
    fx.set_clip(0, H, clip.x0, clip.x1);
  if (clip.margin >= 0)
    fx.set_margin(clip.margin);
  auto basis = std::make_unique<Scan::ColumnBasis<W, H>>();
  const Temporal::ColumnTime time =
      basis->prepare(motion.previous, motion.current);
  HS_EXPECT_TRUE(time.azimuth != nullptr);
  Pipeline<W, H> pipeline;
  {
    Canvas canvas(fx);
    Scan::ColumnTimeScope scope(canvas, time);
    draw(pipeline, canvas, motion.current);
  }
  fx.advance_display();
  std::vector<Pixel> out;
  capture_frame<W, H>(fx, out);
  return out;
}

/**
 * @brief Renders each column from the shapes rebuilt at its sweep time.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 * @param motion Frame endpoints, slerped independently per column.
 * @param draw Callable drawing the shapes at an orientation.
 * @return Frame whose columns x and x + W / 2 come from the fixed-meridian
 *         scan at `slerp(previous, current, (x % (W / 2)) / (W / 2))`.
 */
template <int W, int H, typename DrawFn>
inline std::vector<Pixel> render_column_oracle(const BasisMotion &motion,
                                               DrawFn &&draw) {
  std::vector<Pixel> out(static_cast<size_t>(W) * H);
  std::vector<Pixel> frame;
  for (int column = 0; column < W / 2; ++column) {
    const float time = static_cast<float>(column) / (W / 2);
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H> pipeline;
    {
      Canvas canvas(fx);
      draw(pipeline, canvas,
           math::slerp(motion.previous, motion.current, time).normalized());
    }
    fx.advance_display();
    capture_frame<W, H>(fx, frame);
    for (int y = 0; y < H; ++y)
      for (int x : {column, column + W / 2})
        out[static_cast<size_t>(y) * W + x] =
            frame[static_cast<size_t>(y) * W + x];
  }
  return out;
}

/** @brief Channel ceiling of a faint fringe pixel. */
inline constexpr int BASIS_FRINGE_CHANNEL = 1024;

/**
 * @brief Checks two frames within the composed factored-column tolerance.
 * @param label Shape named in failure context.
 * @param actual Rendered frame.
 * @param expected Oracle frame.
 * @return Number of lit oracle pixels, so a vacuous comparison is visible.
 * @details A faint pixel where the oracle is black is exempt: the oracle's
 * fixed-meridian bounds drop AA fringe that the swept bounds keep.
 */
inline int expect_factored_tolerance(const char *label,
                                     const std::vector<Pixel> &actual,
                                     const std::vector<Pixel> &expected) {
  HS_CONTEXT(label);
  int maximum_error = 0;
  uint64_t error_sum = 0;
  uint64_t channels = 0;
  int lit = 0;
  for (size_t i = 0; i < expected.size(); ++i) {
    lit += !is_black(expected[i]);
    const Pixel &a = actual[i];
    const Pixel &b = expected[i];
    const bool fringe =
        is_black(b) && !is_black(a) && a.r <= BASIS_FRINGE_CHANNEL &&
        a.g <= BASIS_FRINGE_CHANNEL && a.b <= BASIS_FRINGE_CHANNEL;
    for (int error : {std::abs(int(a.r) - b.r), std::abs(int(a.g) - b.g),
                      std::abs(int(a.b) - b.b)}) {
      if (fringe)
        continue;
      maximum_error = std::max(maximum_error, error);
      error_sum += error;
      ++channels;
    }
  }
  HS_EXPECT_TRUE(maximum_error <= 16);
  HS_EXPECT_TRUE(static_cast<double>(error_sum) / channels < 0.05);
  return lit;
}

/** @brief A frame rotating 0.35 rad about a tilted axis. */
inline BasisMotion tilted_motion(float angle = 0.35f) {
  const math::Quaternion previous =
      math::make_rotation(math::Vector(0.3f, 0.8f, -0.5f).normalized(), 0.7f);
  return {previous, math::make_rotation(
                        math::Vector(0.6f, -0.2f, 0.77f).normalized(), angle) *
                        previous};
}

/** @brief Constant-colour shader shared by the oracle shapes. */
inline void basis_flat_shader(const math::Vector &, Fragment &f) {
  f.color = Color4(Pixel(52000, 21000, 9000), 1.0f);
}

/**
 * @brief Builds a scannable octahedron, face offsets included.
 * @param mesh Destination mesh.
 * @param arena Arena backing its arrays.
 */
inline void build_basis_octahedron(MeshState &mesh, Arena &arena) {
  build_meshstate_solid<Solids::Octahedron>(mesh, arena);
  mesh.face_offsets.bind(arena, mesh.face_counts.size());
  int offset = 0;
  for (const uint8_t count : mesh.face_counts) {
    mesh.face_offsets.push_back(static_cast<uint16_t>(offset));
    offset += count;
  }
}

/**
 * @brief Rings, stars, meshes and fused ring scans through the rotated basis
 *        match the shapes rebuilt at each column's sweep time.
 */
inline void test_column_basis_matches_rotated_shapes() {
  constexpr int W = 96, H = 48;
  const ScopedPoleLod lod(0.0f);
  const BasisMotion motion = tilted_motion();

  const auto ring = [](Pipeline<W, H> &pipeline, Canvas &canvas,
                       const math::Quaternion &q) {
    Scan::Ring::draw<W, H, false>(pipeline, canvas,
                                  math::make_basis(q, math::Vector(0, 0, 1)),
                                  0.6f, 0.04f, basis_flat_shader);
  };
  HS_EXPECT_TRUE(
      expect_factored_tolerance("ring", render_column_basis<W, H>(motion, ring),
                                render_column_oracle<W, H>(motion, ring)) > 40);

  const auto star = [](Pipeline<W, H> &pipeline, Canvas &canvas,
                       const math::Quaternion &q) {
    Scan::Star::draw<W, H, false>(pipeline, canvas,
                                  math::make_basis(q, math::Vector(1, 0, 0)),
                                  0.2f, 5, basis_flat_shader);
  };
  HS_EXPECT_TRUE(
      expect_factored_tolerance("star", render_column_basis<W, H>(motion, star),
                                render_column_oracle<W, H>(motion, star)) > 10);

  ScratchScope scope(scratch_arena_a);
  MeshState mesh;
  build_basis_octahedron(mesh, scratch_arena_a);
  std::vector<math::Vector> base(mesh.vertices.begin(), mesh.vertices.end());
  const auto solid = [&](Pipeline<W, H> &pipeline, Canvas &canvas,
                         const math::Quaternion &q) {
    for (size_t i = 0; i < base.size(); ++i)
      mesh.vertices[i] = math::rotate(base[i], q);
    Scan::Mesh::draw<W, H>(
        pipeline, canvas, mesh,
        [](const math::Vector &, Fragment &f) {
          const auto face = static_cast<uint16_t>(f.v2);
          f.color = Color4(
              Pixel(static_cast<uint16_t>(6000 + 7000 * face), 30000, 12000),
              1.0f);
        },
        scratch_arena_b);
  };
  HS_EXPECT_TRUE(expect_factored_tolerance(
                     "mesh", render_column_basis<W, H>(motion, solid),
                     render_column_oracle<W, H>(motion, solid)) > 1000);

  constexpr int N = 3;
  const auto group = [](Pipeline<W, H> &pipeline, Canvas &canvas,
                        const math::Quaternion &q) {
    SDF::Ring shapes[N];
    for (int s = 0; s < N; ++s)
      shapes[s] = SDF::Ring(
          math::make_basis(q * math::make_rotation(math::X_AXIS, 0.15f * s),
                           math::Vector(0.77f, 0.0f, -0.6f).normalized()),
          0.35f, 0.05f);
    Scan::RingGroup::draw<W, H>(
        pipeline, canvas, shapes, N,
        [](int s, const math::Vector &, Fragment &f) {
          f.color = Color4(
              Pixel(static_cast<uint16_t>(20000 + 15000 * s), 9000, 40000),
              0.8f);
        });
  };
  HS_EXPECT_TRUE(expect_factored_tolerance(
                     "ring group", render_column_basis<W, H>(motion, group),
                     render_column_oracle<W, H>(motion, group)) > 100);

  constexpr int RINGS = 3, LUT_N = 32;
  static float knots[RINGS][LUT_N + 1];
  build_ring_knots<RINGS, LUT_N>(knots);
  const auto stack = [&](Pipeline<W, H> &pipeline, Canvas &canvas,
                         const math::Quaternion &q) {
    const math::Basis basis = math::make_basis(q, math::Vector(0, 0, 1));
    alignas(SDF::DistortedRing) unsigned char
        storage[RINGS * sizeof(SDF::DistortedRing)];
    auto *shapes = reinterpret_cast<SDF::DistortedRing *>(storage);
    int8_t slots[RINGS];
    for (int i = 0; i < RINGS; ++i) {
      new (&shapes[i]) SDF::DistortedRing(basis, 0.2f + 0.1f * i, 0.06f,
                                          knots[i], LUT_N, 0.0f, nullptr);
      slots[i] = static_cast<int8_t>(i);
    }
    static Scan::DistortedRingStack::CandidateTable<W, H> table;
    Scan::DistortedRingStack::draw<W, H>(
        pipeline, canvas, RINGS, shapes, slots, RINGS, table,
        [](int s, const math::Vector &, Fragment &f) {
          f.color = Color4(Pixel(static_cast<uint16_t>(15000 * (s + 1)),
                                 static_cast<uint16_t>(f.v0 * 8000.0f), 7000),
                           f.v2);
        });
  };
  const BasisMotion swing{math::Quaternion(),
                          math::make_rotation(math::X_AXIS, 0.8f)};
  for (const BasisMotion &frame :
       {swing, BasisMotion{swing.current, swing.previous}})
    HS_EXPECT_TRUE(expect_factored_tolerance(
                       "ring stack", render_column_basis<W, H>(frame, stack),
                       render_column_oracle<W, H>(frame, stack)) > 100);
}

/**
 * @brief Clipped segments of a rotated-basis scan reproduce the single-engine
 *        frame, across both arm boundaries, and motionless and fidelity-off
 *        frames keep the fixed meridian.
 */
inline void test_column_basis_segments_and_endpoints() {
  constexpr int W = 96, H = 48;
  const ScopedPoleLod lod(0.0f);
  const BasisMotion motion = tilted_motion();
  const auto ring = [](Pipeline<W, H> &pipeline, Canvas &canvas,
                       const math::Quaternion &q) {
    Scan::Ring::draw<W, H, false>(
        pipeline, canvas, math::make_basis(q, math::Vector(0.2f, 0.1f, 1)),
        0.9f, 0.05f, basis_flat_shader);
  };
  const std::vector<Pixel> whole = render_column_basis<W, H>(motion, ring);
  int compared = 0;
  for (const BasisClip clip : {BasisClip{0, 20, 4}, BasisClip{20, 52, 4},
                               BasisClip{52, 96, 4}, BasisClip{40, 56, 0}}) {
    const std::vector<Pixel> part =
        render_column_basis<W, H>(motion, ring, clip);
    for (int y = 0; y < H; ++y)
      for (int x = clip.x0; x < clip.x1; ++x) {
        const size_t i = static_cast<size_t>(y) * W + x;
        HS_EXPECT_TRUE(part[i] == whole[i]);
        compared += !is_black(whole[i]);
      }
  }
  HS_EXPECT_TRUE(compared > 100);

  const auto fixed = [&](const math::Quaternion &q) {
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H> pipeline;
    {
      Canvas canvas(fx);
      ring(pipeline, canvas, q);
    }
    fx.advance_display();
    std::vector<Pixel> out;
    capture_frame<W, H>(fx, out);
    return out;
  };
  auto basis = std::make_unique<Scan::ColumnBasis<W, H>>();
  Temporal::set_temporal_sampling_enabled(false);
  const Temporal::ColumnTime off =
      basis->prepare(motion.previous, motion.current);
  Temporal::set_temporal_sampling_enabled(true);
  HS_EXPECT_TRUE(off.azimuth == nullptr);
  HS_EXPECT_TRUE(basis->prepare(motion.current, motion.current).azimuth ==
                 nullptr);
  std::vector<Pixel> drawn;
  {
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H> pipeline;
    {
      Canvas canvas(fx);
      Scan::ColumnTimeScope scope(canvas, off);
      ring(pipeline, canvas, motion.current);
    }
    fx.advance_display();
    capture_frame<W, H>(fx, drawn);
  }
  HS_EXPECT_TRUE(drawn == fixed(motion.current));
  HS_EXPECT_TRUE(whole != fixed(motion.current));
}

/**
 * @brief Builds a fan of small triangles around a centre direction.
 * @param mesh Destination mesh.
 * @param arena Arena backing its arrays.
 * @param centre Fan centre, unit length.
 * @param radius Angular radius of the rim (radians).
 */
inline void build_basis_fan(MeshState &mesh, Arena &arena,
                            const math::Vector &centre, float radius) {
  constexpr int RIM = 6;
  const math::Basis frame = math::make_basis(math::Quaternion(), centre);
  mesh.vertices.bind(arena, RIM + 1);
  mesh.vertices.push_back(centre);
  for (int i = 0; i < RIM; ++i) {
    const float azimuth = math::TWO_PI_F * i / RIM;
    mesh.vertices.push_back(
        (centre * cosf(radius) +
         (frame.u * cosf(azimuth) + frame.w * sinf(azimuth)) * sinf(radius))
            .normalized());
  }
  mesh.face_counts.bind(arena, RIM);
  mesh.faces.bind(arena, 3 * RIM);
  mesh.face_offsets.bind(arena, RIM);
  for (int i = 0; i < RIM; ++i) {
    mesh.face_counts.push_back(3);
    mesh.face_offsets.push_back(static_cast<uint16_t>(3 * i));
    mesh.faces.push_back(0);
    mesh.faces.push_back(static_cast<uint16_t>(1 + i));
    mesh.faces.push_back(static_cast<uint16_t>(1 + (i + 1) % RIM));
  }
}

/**
 * @brief Pole-LOD blocks under a rotated basis settle only across rotations
 *        their slack covers.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 */
template <int W, int H> inline void check_column_basis_blocks() {
  HS_CONTEXT("width", static_cast<long long>(W));
  const ScopedPoleLod lod(0.4f);
  const math::Quaternion previous = math::make_rotation(math::X_AXIS, -0.3f);
  const BasisMotion motion{previous, math::Quaternion()};
  ScratchScope scope(scratch_arena_a);
  MeshState mesh;
  build_basis_fan(mesh, scratch_arena_a,
                  math::Vector(0.05f, 1.0f, 0.08f).normalized(), 0.12f);
  std::vector<math::Vector> base(mesh.vertices.begin(), mesh.vertices.end());
  const auto fan = [&](Pipeline<W, H> &pipeline, Canvas &canvas,
                       const math::Quaternion &q) {
    for (size_t i = 0; i < base.size(); ++i)
      mesh.vertices[i] = math::rotate(base[i], q);
    Scan::Mesh::draw<W, H>(
        pipeline, canvas, mesh,
        [](const math::Vector &, Fragment &f) {
          const auto face = static_cast<uint16_t>(f.v2);
          f.color = Color4(
              Pixel(static_cast<uint16_t>(6000 + 7000 * face), 30000, 12000),
              1.0f);
        },
        scratch_arena_b);
  };
  const std::vector<Pixel> decimated = render_column_basis<W, H>(motion, fan);
  std::vector<Pixel> exact;
  {
    const ScopedPoleLod off(0.0f);
    exact = render_column_oracle<W, H>(motion, fan);
  }
  HS_EXPECT_TRUE(expect_factored_tolerance("pole-LOD mesh", decimated, exact) >
                 100);
}

/**
 * @brief Pole-LOD blocks under a rotated basis settle only across rotations
 *        their slack covers, at both test widths.
 */
inline void test_column_basis_blocks_cover_their_rotation() {
  check_column_basis_blocks<96, 48>();
  check_column_basis_blocks<288, 144>();
}

/** @brief SDF::Star bounded only by its swept() proxy. */
struct SweptOnlyStar : SDF::Star {
  /**
   * @brief Wraps a star.
   * @param star Star whose basis outlives the wrapper.
   */
  explicit SweptOnlyStar(const SDF::Star &star) : SDF::Star(star) {}
  /**
   * @brief The star's swept proxy, kept in the wrapper type.
   * @param arc Sweep arc (radians).
   * @return Wrapped swept proxy.
   */
  SweptOnlyStar swept(float arc) const {
    return SweptOnlyStar(SDF::Star::swept(arc));
  }
  /** @brief Hides SDF::Star::moved, so scans fall back to swept(). */
  void moved() const = delete;
};

/** @brief SDF::Ring bounded only by its swept() proxy. */
struct SweptOnlyRing : SDF::Ring {
  /**
   * @brief Wraps a ring.
   * @param ring Ring to wrap.
   */
  explicit SweptOnlyRing(const SDF::Ring &ring) : SDF::Ring(ring) {}
  /**
   * @brief The ring's swept proxy, kept in the wrapper type.
   * @param arc Sweep arc (radians).
   * @return Wrapped swept proxy.
   */
  SweptOnlyRing swept(float arc) const {
    return SweptOnlyRing(SDF::Ring::swept(arc));
  }
  /** @brief Hides SDF::Ring::moved, so scans fall back to swept(). */
  void moved() const = delete;
};

/**
 * @brief Small stars and discs bounded by caps re-centred at their columns'
 *        sweep time draw the frame their whole-frame swept bounds draw.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 * @param angle Frame rotation (radians); below 0.1 most caps must re-centre
 *        and move by a small share of it.
 */
template <int W, int H> inline void check_lit_caps(float angle) {
  HS_CONTEXT(W == 96 ? "lit caps 96" : "lit caps 288");
  const ScopedPoleLod lod(0.0f);
  const BasisMotion motion = tilted_motion(angle);
  constexpr int SHAPES = 160;
  std::vector<math::Basis> bases;
  for (int i = 0; i < SHAPES; ++i)
    bases.push_back(math::make_basis(math::Quaternion(),
                                     math::fib_spiral(SHAPES, 0.5f, i)));
  const auto draw = [&](bool moved) {
    return render_column_basis<W, H>(motion, [&](Pipeline<W, H> &pipeline,
                                                 Canvas &canvas,
                                                 const math::Quaternion &) {
      for (int i = 0; i < SHAPES; ++i) {
        const auto star_basis = math::get_antipode(bases[i], 0.04f);
        const SDF::Star star(star_basis.first, star_basis.second, 5, 0.0f);
        const SDF::Ring disc(bases[(i * 7) % SHAPES], 0.0f, 0.03f);
        if (moved) {
          Scan::rasterize<W, H, false>(pipeline, canvas, star,
                                       basis_flat_shader);
          Scan::rasterize<W, H, false>(pipeline, canvas, disc,
                                       basis_flat_shader);
        } else {
          Scan::rasterize<W, H, false>(pipeline, canvas, SweptOnlyStar(star),
                                       basis_flat_shader);
          Scan::rasterize<W, H, false>(pipeline, canvas, SweptOnlyRing(disc),
                                       basis_flat_shader);
        }
      }
    });
  };
  const std::vector<Pixel> lit = draw(true);
  const std::vector<Pixel> swept = draw(false);
  int differing = 0, drawn = 0;
  for (size_t i = 0; i < lit.size(); ++i) {
    differing += !(lit[i] == swept[i]);
    drawn += !is_black(swept[i]);
  }
  HS_EXPECT_EQ(differing, 0);
  HS_EXPECT_TRUE(drawn > 4 * SHAPES);

  auto basis = std::make_unique<Scan::ColumnBasis<W, H>>();
  const Temporal::ColumnTime time =
      basis->prepare(motion.previous, motion.current);
  int recentred = 0;
  float arc_sum = 0.0f;
  for (const math::Basis &b : bases) {
    Scan::LitCap cap;
    if (!Scan::lit_cap<W, H>(b.v, 0.03f, time, cap))
      continue;
    ++recentred;
    arc_sum += cap.arc;
    HS_EXPECT_LT(cap.arc, time.sweep);
  }
  if (angle < 0.1f) {
    HS_EXPECT_TRUE(recentred > SHAPES * 3 / 4);
    HS_EXPECT_LT(arc_sum / recentred, 0.15f * time.sweep);
  }
}

/**
 * @brief Re-centred caps keep every pixel the swept bounds draw, at both test
 *        widths.
 */
inline void test_lit_caps_keep_swept_frames() {
  for (const float angle : {0.35f, 0.05f}) {
    check_lit_caps<96, 48>(angle);
    check_lit_caps<288, 144>(angle);
  }
}

/**
 * @brief A disc with a birth fraction lights only the columns swept at or
 *        after it, on both arms.
 */
inline void test_point_birth_skips_earlier_columns() {
  constexpr int W = 96, H = 48;
  const math::Vector centre = math::Vector(0.0f, 0.2f, 1.0f).normalized();
  const auto draw = [&](const Temporal::ColumnTime *time) {
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H> pipeline;
    {
      Canvas canvas(fx);
      if (time)
        canvas.set_column_time(time);
      Scan::Point::draw<W, H>(pipeline, canvas, centre, 0.5f,
                              basis_flat_shader);
      Scan::Point::draw<W, H>(pipeline, canvas, -centre, 0.5f,
                              basis_flat_shader);
    }
    fx.advance_display();
    std::vector<Pixel> out;
    capture_frame<W, H>(fx, out);
    return out;
  };
  const std::vector<Pixel> whole = draw(nullptr);
  Temporal::ColumnTime born;
  born.birth = 0.4f;
  const std::vector<Pixel> cut = draw(&born);
  int before = 0, after = 0;
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x) {
      const size_t i = static_cast<size_t>(y) * W + x;
      if (Temporal::sweep_fraction(x, W) < born.birth) {
        before += !is_black(whole[i]);
        HS_EXPECT_TRUE(is_black(cut[i]));
      } else {
        after += !is_black(whole[i]);
        HS_EXPECT_TRUE(cut[i] == whole[i]);
      }
    }
  HS_EXPECT_TRUE(before > 20);
  HS_EXPECT_TRUE(after > 20);
  HS_EXPECT_EQ(Temporal::first_column_at(0.4f, W), 20);
  HS_EXPECT_EQ(Temporal::first_column_at(0.0f, W), 0);
  HS_EXPECT_EQ(Temporal::first_column_at(1.0f, W), W / 2);
}

/**
 * @brief A volume marched through per-column placements matches the volume
 *        placed at each column's sweep time.
 */
inline void test_volume_columns_match_rotated_placements() {
  constexpr int W = 96, H = 48;
  const SDF::Torus torus{0.25f, 0.08f};
  const math::Vector anchor = math::Vector(0.3f, 0.4f, 0.85f).normalized();
  const BasisMotion camera{math::Quaternion(),
                           math::make_rotation(math::Y_AXIS, 0.6f)};
  const auto spin = [](float u) {
    return math::make_rotation(math::Vector(1.0f, 0.2f, 0.1f).normalized(),
                               1.2f * u);
  };
  constexpr float BOUNDS = 0.36f;
  const auto shade = [](const math::Vector &, Fragment &f) {
    f.color = Color4(Pixel(41000, 23000, 7000), 1.0f);
  };
  const auto render = [&](auto &&draw) {
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H> pipeline;
    {
      Canvas canvas(fx);
      draw(pipeline, canvas);
    }
    fx.advance_display();
    std::vector<Pixel> out;
    capture_frame<W, H>(fx, out);
    return out;
  };
  const auto place = [&](float u) {
    const math::Quaternion q =
        math::slerp(camera.previous, camera.current, u).normalized();
    return std::pair{math::rotate(anchor, q), q * spin(u)};
  };

  auto columns = std::make_unique<Scan::VolumeColumns<W>>();
  for (int c = 0; c < W / 2; ++c) {
    const auto [centre, q] = place(Temporal::sweep_fraction(c, W));
    columns->set(c, centre, q);
  }
  const auto [centre, q] = place(1.0f);
  const std::vector<Pixel> swept = render([&](auto &pipeline, Canvas &canvas) {
    Scan::TransformedVolume volume(torus, centre, q);
    Scan::Volume::draw<W, H>(pipeline, canvas, centre, BOUNDS, volume, shade,
                             24, 0.02f, columns.get());
  });

  std::vector<Pixel> expected(static_cast<size_t>(W) * H);
  for (int c = 0; c < W / 2; ++c) {
    const auto [at, rotation] = place(Temporal::sweep_fraction(c, W));
    const std::vector<Pixel> frame =
        render([&](auto &pipeline, Canvas &canvas) {
          Scan::TransformedVolume volume(torus, at, rotation);
          Scan::Volume::draw<W, H>(pipeline, canvas, at, BOUNDS, volume, shade,
                                   24, 0.02f);
        });
    for (int y = 0; y < H; ++y)
      for (int x : {c, c + W / 2})
        expected[static_cast<size_t>(y) * W + x] =
            frame[static_cast<size_t>(y) * W + x];
  }
  HS_EXPECT_TRUE(expect_factored_tolerance("volume", swept, expected) > 30);
}

/**
 * @brief A ring group with previous axes matches its rings rebuilt about the
 *        axes interpolated to each column's sweep time.
 */
inline void test_ring_group_previous_axes_match_column_rings() {
  constexpr int W = 96, H = 48;
  constexpr int N = 3;
  const ScopedPoleLod lod(0.0f);
  math::Vector previous[N], current[N];
  for (int s = 0; s < N; ++s) {
    current[s] = math::rotate(math::Vector(0.77f, 0.0f, -0.6f).normalized(),
                              math::make_rotation(math::X_AXIS, 0.3f * s));
    previous[s] = math::rotate(
        current[s],
        math::make_rotation(math::Vector(0.6f, -0.2f, 0.77f).normalized(),
                            0.25f + 0.1f * s));
  }
  const auto draw = [&](const math::Vector *axes, const math::Vector *before) {
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H> pipeline;
    {
      Canvas canvas(fx);
      SDF::Ring shapes[N];
      for (int s = 0; s < N; ++s)
        shapes[s] = SDF::Ring(math::make_basis(math::Quaternion(), axes[s]),
                              0.4f, 0.05f);
      Scan::RingGroup::draw<W, H>(
          pipeline, canvas, shapes, N,
          [](int s, const math::Vector &, Fragment &f) {
            f.color = Color4(
                Pixel(static_cast<uint16_t>(20000 + 15000 * s), 9000, 40000),
                0.8f);
          },
          false, before);
    }
    fx.advance_display();
    std::vector<Pixel> out;
    capture_frame<W, H>(fx, out);
    return out;
  };
  const std::vector<Pixel> swept = draw(current, previous);
  std::vector<Pixel> expected(static_cast<size_t>(W) * H);
  for (int c = 0; c < W / 2; ++c) {
    const float u = static_cast<float>(c) / (W / 2);
    math::Vector axes[N];
    for (int s = 0; s < N; ++s)
      axes[s] = (previous[s] + (current[s] - previous[s]) * u).normalized();
    const std::vector<Pixel> frame = draw(axes, nullptr);
    for (int y = 0; y < H; ++y)
      for (int x : {c, c + W / 2})
        expected[static_cast<size_t>(y) * W + x] =
            frame[static_cast<size_t>(y) * W + x];
  }
  HS_EXPECT_TRUE(expect_factored_tolerance("previous axes", swept, expected) >
                 100);
}

/**
 * @brief An object-frame walk hands each pixel its center and sub-samples
 *        rotated into the object frame at its column's sweep time, on both
 *        arms, and visits only the columns the render-x window admits.
 */
inline void test_object_column_basis_matches_rotated_meridian() {
  constexpr int W = 32, H = 12, HALF = W / 2;
  const BasisMotion motion = tilted_motion();
  auto basis = std::make_unique<Scan::ObjectColumnBasis<W, H>>();
  basis->prepare(motion.previous, motion.current);
  if (!math::TrigLUT<W, H>::initialized)
    math::TrigLUT<W, H>::init();
  for (const ColumnWindow window : {ColumnWindow{}, ColumnWindow{3, 9, true}}) {
    HS_CONTEXT("window", window.begin, window.end);
    hs_test::StubEffect fx(W, H);
    std::vector<int> visits(static_cast<size_t>(W) * H, 0);
    float worst = 0.0f;
    {
      Canvas canvas(fx);
      std::optional<RenderXWindow> scope;
      if (window.active)
        scope.emplace(canvas, window.begin, window.end);
      int y = -1, last = W;
      Scan::Shader::walk_object_grid<W, H>(
          canvas, *basis,
          [&](const math::Vector &center, const Scan::ObjectColumn &column,
              const auto &grid, int x) {
            y += x <= last;
            last = x;
            ++visits[static_cast<size_t>(y) * W + x];
            const math::Quaternion inverse =
                math::slerp(motion.previous, motion.current,
                            Temporal::sweep_fraction(x, W))
                    .normalized()
                    .conjugate();
            const float sp = math::TrigLUT<W, H>::sin_phi[y];
            const float cp = math::TrigLUT<W, H>::cos_phi[y];
            const math::Vector meridian(sp * math::TrigLUT<W, H>::cos_theta(x),
                                        cp,
                                        sp * math::TrigLUT<W, H>::sin_theta[x]);
            const auto gap = [](const math::Vector &a, const math::Vector &b) {
              return std::max({std::fabs(a.x - b.x), std::fabs(a.y - b.y),
                               std::fabs(a.z - b.z)});
            };
            worst =
                std::max(worst, gap(center, math::rotate(meridian, inverse)));
            for (int i = 0; i < 4; ++i)
              worst =
                  std::max(worst, gap(grid.at(column, i),
                                      math::rotate(grid.at(x, i), inverse)));
            return Pixel(0, 0, 0);
          });
    }
    HS_EXPECT_LT(worst, 2e-6f);
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x) {
        HS_CONTEXT("px", x, y);
        const bool admitted = !window.active || (x % HALF >= window.begin &&
                                                 x % HALF < window.end);
        HS_EXPECT_EQ(visits[static_cast<size_t>(y) * W + x], admitted ? 1 : 0);
      }
  }
}
