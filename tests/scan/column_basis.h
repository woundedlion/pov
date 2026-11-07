/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Column-time rotated basis: scans against shapes rebuilt at column time
// ============================================================================

template <typename P, typename Shader>
concept VolumeColumnPipeline =
    requires(P &pipeline, Canvas &canvas, const math::Vector &center,
             const Scan::TransformedVolume<SDF::Torus> &shape, Shader shader,
             const Scan::VolumeColumns<96> &columns) {
      Scan::Volume::draw<96, 48>(pipeline, canvas, center, 1.0f, shape, shader,
                                 15, 0.01f, columns);
    };

using MovingVolumePipeline = Pipeline<96, 48, Filter::World::Replicate<96>>;
static_assert(VolumeColumnPipeline<Pipeline<96, 48>, FragmentShaderFn>);
static_assert(VolumeColumnPipeline<Pipeline<96, 48>, TimedShaderFn>);
static_assert(!VolumeColumnPipeline<MovingVolumePipeline, FragmentShaderFn>);
static_assert(!VolumeColumnPipeline<MovingVolumePipeline, TimedShaderFn>);
static_assert(requires(MovingVolumePipeline &pipeline, Canvas &canvas,
                       const math::Vector &center,
                       const Scan::TransformedVolume<SDF::Torus> &shape,
                       FragmentShaderFn shader) {
  Scan::Volume::draw<96, 48>(pipeline, canvas, center, 1.0f, shape, shader);
});

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
 * @param state_moves Whether geometry changes over the frame.
 * @return Displayed frame.
 */
template <int W, int H, typename DrawFn>
inline std::vector<Pixel>
render_column_basis(const BasisMotion &motion, DrawFn &&draw,
                    const BasisClip &clip = {}, bool state_moves = false) {
  hs_test::StubEffect fx(W, H);
  if (clip.x1 > 0)
    fx.set_clip(0, H, clip.x0, clip.x1);
  if (clip.margin >= 0)
    fx.set_margin(clip.margin);
  auto basis = std::make_unique<Scan::ColumnBasisTable<W, H>>();
  const Scan::ColumnBasis time =
      basis->prepare(motion.previous, motion.current, state_moves);
  HS_EXPECT_TRUE(time.azimuth != nullptr);
  Pipeline<W, H> pipeline;
  {
    Canvas canvas(fx);
    Scan::ColumnBasisScope scope(canvas, time);
    draw(pipeline, canvas, motion.current);
  }
  fx.advance_display();
  std::vector<Pixel> out;
  capture_frame<W, H>(fx, out);
  return out;
}

/**
 * @brief Renders each column pair from a fixed-meridian draw of the state at
 *        the pair's sweep time.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 * @param draw Callable `(Pipeline<W, H> &, Canvas &, float u)`.
 * @return Frame whose columns x and x + W / 2 come from @p draw at
 *         `u = (x % (W / 2)) / (W / 2)`.
 */
template <int W, int H, typename DrawFn>
inline std::vector<Pixel> render_state_oracle(DrawFn &&draw) {
  return column_time_oracle<W, H>([&](float u) {
    return render_stub_frame<W, H>(
        [&](Pipeline<W, H> &pipeline, Canvas &canvas) {
          draw(pipeline, canvas, u);
        });
  });
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
  return render_state_oracle<W, H>(
      [&](auto &pipeline, Canvas &canvas, float u) {
        draw(pipeline, canvas,
             math::slerp(motion.previous, motion.current, u).normalized());
      });
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
  HS_EXPECT_LE(maximum_error, 16);
  HS_EXPECT_LT(static_cast<double>(error_sum) / channels, 0.05);
  return lit;
}

/** @brief A frame rotating @p angle about a tilted axis.
 * @param angle Rotation angle in radians. */
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
  HS_EXPECT_GT(
      expect_factored_tolerance("ring", render_column_basis<W, H>(motion, ring),
                                render_column_oracle<W, H>(motion, ring)),
      40);

  const auto star = [](Pipeline<W, H> &pipeline, Canvas &canvas,
                       const math::Quaternion &q) {
    Scan::Star::draw<W, H, false>(pipeline, canvas,
                                  math::make_basis(q, math::Vector(1, 0, 0)),
                                  0.2f, 5, basis_flat_shader);
  };
  HS_EXPECT_GT(
      expect_factored_tolerance("star", render_column_basis<W, H>(motion, star),
                                render_column_oracle<W, H>(motion, star)),
      10);

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
  HS_EXPECT_GT(expect_factored_tolerance(
                   "mesh", render_column_basis<W, H>(motion, solid),
                   render_column_oracle<W, H>(motion, solid)),
               1000);

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
  HS_EXPECT_GT(expect_factored_tolerance(
                   "ring group", render_column_basis<W, H>(motion, group),
                   render_column_oracle<W, H>(motion, group)),
               100);

  constexpr int RINGS = 3, LUT_N = 32;
  static float knots[RINGS][LUT_N + 1];
  build_ring_knots<RINGS, LUT_N>(knots);
  const auto stack = [&](Pipeline<W, H> &pipeline, Canvas &canvas,
                         const math::Quaternion &q) {
    const math::Basis basis = math::make_basis(q, math::Vector(0, 0, 1));
    DistortedRingArray<RINGS> rings;
    auto *shapes = rings.data();
    int8_t slots[RINGS];
    for (int i = 0; i < RINGS; ++i) {
      rings.emplace(i, basis, 0.2f + 0.1f * i, 0.06f, knots[i], LUT_N, 0.0f,
                    nullptr);
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
    HS_EXPECT_GT(expect_factored_tolerance(
                     "ring stack", render_column_basis<W, H>(frame, stack),
                     render_column_oracle<W, H>(frame, stack)),
                 100);
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
        HS_EXPECT_EQ(part[i], whole[i]);
        compared += !is_black(whole[i]);
      }
  }
  HS_EXPECT_GT(compared, 100);

  const auto fixed = [&](const math::Quaternion &q) {
    return render_stub_frame<W, H>(
        [&](Pipeline<W, H> &pipeline, Canvas &canvas) {
          ring(pipeline, canvas, q);
        });
  };
  auto basis = std::make_unique<Scan::ColumnBasisTable<W, H>>();
  const Scan::ColumnBasis off = [&] {
    const ScopedTemporalSampling sampling(false);
    return basis->prepare(motion.previous, motion.current);
  }();
  HS_EXPECT_TRUE(off.azimuth == nullptr);
  HS_EXPECT_TRUE(basis->prepare(motion.current, motion.current).azimuth ==
                 nullptr);
  const std::vector<Pixel> drawn =
      render_stub_frame<W, H>([&](auto &pipeline, Canvas &canvas) {
        Scan::ColumnBasisScope scope(canvas, off);
        ring(pipeline, canvas, motion.current);
      });
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
  HS_EXPECT_GT(expect_factored_tolerance("pole-LOD mesh", decimated, exact),
               100);

  const math::Vector centre = math::Vector(0.05f, 1.0f, 0.08f).normalized();
  const math::Vector dot = math::Vector(-0.2f, 1.0f, 0.1f).normalized();
  const auto polar = [&](Pipeline<W, H> &pipeline, Canvas &canvas,
                         const math::Quaternion &q) {
    Scan::Star::draw<W, H, false>(pipeline, canvas, math::make_basis(q, centre),
                                  0.1f, 5, basis_flat_shader);
    Scan::Point::draw<W, H>(pipeline, canvas, math::rotate(dot, q), 0.05f,
                            basis_flat_shader);
  };
  const std::vector<Pixel> decimated_shapes =
      render_column_basis<W, H>(motion, polar);
  std::vector<Pixel> exact_shapes;
  {
    const ScopedPoleLod off(0.0f);
    exact_shapes = render_column_oracle<W, H>(motion, polar);
  }
  HS_EXPECT_GT(expect_factored_tolerance("pole-LOD shapes", decimated_shapes,
                                         exact_shapes),
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
 * @param angle Frame rotation (radians); every cap re-centres within half of
 *        it, and below 0.1 two in three within a quarter.
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
  HS_EXPECT_GT(drawn, 4 * SHAPES);

  auto basis = std::make_unique<Scan::ColumnBasisTable<W, H>>();
  const Scan::ColumnBasis time =
      basis->prepare(motion.previous, motion.current);
  int tight = 0;
  for (const math::Basis &b : bases) {
    Scan::LitCap cap;
    HS_EXPECT_TRUE((Scan::lit_cap<W, H>(b.v, 0.03f, time, cap)));
    HS_EXPECT_LE(cap.arc, 0.5f * time.sweep_arc + 0.01f);
    tight += cap.arc < 0.25f * time.sweep_arc;
  }
  if (angle < 0.1f)
    HS_EXPECT_GT(tight, SHAPES * 2 / 3);
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
 * @brief A column-time scan over a shape's swept() proxy visits every pixel
 *        the shape lights at its column's rotation, and fewer than whole rows.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 * @tparam ShapeT SDF::SweptBounded shape.
 * @param label Shape named in failure context.
 * @param shape Shape at the current endpoint.
 * @param time Column-time state with a rotated basis.
 */
template <int W, int H, typename ShapeT>
inline void check_swept_proxy_covers(const char *label, const ShapeT &shape,
                                     const Scan::ColumnBasis &time) {
  HS_CONTEXT(label);
  const ScopedPoleLod lod(0.0f);
  const auto proxy = shape.swept(time.sweep_arc);
  std::vector<uint8_t> visited(static_cast<size_t>(W) * H, 0);
  const SDF::Bounds rows = proxy.template get_vertical_bounds<H>();
  Scan::scan_region_column_time<W, H>(
      std::max(0, rows.y_min), std::min(H - 1, rows.y_max),
      [&](int y, auto &&out) {
        return proxy.template get_horizontal_intervals<W, H>(y, out);
      },
      [&](int x, int y, const math::Vector &, int) {
        visited[static_cast<size_t>(y) * W + x] = 1;
        return 1;
      },
      {}, time);
  constexpr float THRESHOLD = ShapeT::is_solid ? 2.0f * math::PI_F / W : 0.0f;
  int lit = 0, missed = 0, scanned = 0;
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x) {
      const bool seen = visited[static_cast<size_t>(y) * W + x] != 0;
      scanned += seen;
      const math::Vector p = Scan::rotated_direction<W>(
          time.azimuth, time.up, x, math::TrigLUT<W, H>::sin_phi[y],
          math::TrigLUT<W, H>::cos_phi[y]);
      SDF::DistanceResult result;
      shape.template distance<false>(p, result);
      if (result.dist < THRESHOLD) {
        ++lit;
        missed += !seen;
      }
    }
  const SDF::Bounds whole = Scan::swept_rows<H>(
      shape.template get_vertical_bounds<H>(), time.sweep_arc);
  HS_EXPECT_EQ(missed, 0);
  HS_EXPECT_GT(lit, 20);
  HS_EXPECT_LT(scanned, (whole.y_max - whole.y_min + 1) * W);
}

/**
 * @brief Distorted rings, lines and CSG compositions bound their column-time
 *        scans by swept proxies that keep every lit pixel.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 * @param motion Frame endpoints, turning at least 0.25 rad.
 */
template <int W, int H>
inline void check_swept_proxies(const BasisMotion &motion) {
  HS_CONTEXT(W == 96 ? "swept proxies 96" : "swept proxies 288");
  if (!math::TrigLUT<W, H>::initialized)
    math::TrigLUT<W, H>::init();
  auto basis = std::make_unique<Scan::ColumnBasisTable<W, H>>();
  const Scan::ColumnBasis time =
      basis->prepare(motion.previous, motion.current);
  HS_EXPECT_GT(time.sweep_arc, 0.25f);

  const math::Basis tilted = math::make_basis(
      math::Quaternion(), math::Vector(0.4f, 0.5f, 0.75f).normalized());
  const ScalarFn wave = [](float t) {
    return 0.08f * sinf(3.0f * math::TWO_PI_F * t);
  };
  const SDF::DistortedRing wavy(tilted, 0.5f, 0.03f, wave, 0.08f, 0.0f);
  check_swept_proxy_covers<W, H>("distorted ring", wavy, time);

  float knots[24];
  for (int k = 0; k < 24; ++k)
    knots[k] = 0.06f * sinf(2.0f * math::TWO_PI_F * k / 24.0f) - 0.02f;
  SDF::KnotPrefilter prefilter;
  const SDF::DistortedRing knotted(tilted, 0.7f, 0.025f, knots, 24, 0.0f,
                                   prefilter);
  check_swept_proxy_covers<W, H>("knot ring", knotted, time);

  const SDF::FlatDistortedRing flat(tilted, 0.3f, 0.03f);
  check_swept_proxy_covers<W, H>("flat ring", flat, time);

  const SDF::Line line(math::Vector(0.8f, 0.3f, 0.2f).normalized(),
                       math::Vector(0.1f, 0.4f, 0.9f).normalized(), 0.05f);
  check_swept_proxy_covers<W, H>("line", line, time);

  const SDF::Ring ring(tilted, 0.4f, 0.03f);
  const math::Basis star_basis = math::make_basis(
      math::Quaternion(), math::Vector(0.9f, 0.25f, 0.3f).normalized());
  const SDF::Star star(star_basis, 0.25f, 5, 0.0f);
  const SDF::Star small(star_basis, 0.12f, 4, 0.0f);
  const SDF::PlanarPolygon square(star_basis, 0.2f, 4, 0.0f);
  const SDF::PlanarPolygon apart(tilted, 0.15f, 4, 0.0f);
  check_swept_proxy_covers<W, H>("union", SDF::Union(ring, line), time);
  check_swept_proxy_covers<W, H>("smooth union",
                                 SDF::SmoothUnion(star, apart, 0.05f), time);
  check_swept_proxy_covers<W, H>("intersection",
                                 SDF::Intersection(small, square), time);
  check_swept_proxy_covers<W, H>("subtract", SDF::Subtract(star, small), time);
  const SDF::Union<SDF::Ring, SDF::Line> inner(ring, line);
  check_swept_proxy_covers<W, H>("nested union", SDF::Union(inner, flat), time);
  const math::Basis petal = math::make_basis(
      math::Quaternion(), math::Vector(0.9f, 0.35f, 0.0f).normalized());
  const SDF::Star seed(petal, 0.1f, 5, 0.0f);
  check_swept_proxy_covers<W, H>("repeat", SDF::AngularRepeat(seed, 6), time);
}

/**
 * @brief Swept proxies keep every lit pixel of distorted rings, lines and CSG
 *        compositions, at both test widths.
 */
inline void test_swept_proxies_keep_lit_pixels() {
  const BasisMotion tilt{math::Quaternion(),
                         math::make_rotation(math::X_AXIS, 0.3f)};
  for (const BasisMotion &motion : {tilted_motion(0.3f), tilt}) {
    check_swept_proxies<96, 48>(motion);
    check_swept_proxies<288, 144>(motion);
  }
}

/**
 * @brief Spacetime::rotate_by() turns a direction as the exact rotation does,
 *        short or long, and a drifting shape's column entries follow it.
 */
inline void test_rotate_by_matches_exact_rotation() {
  const math::Vector axis = math::Vector(0.3f, -0.7f, 0.5f).normalized();
  const math::Vector p = math::Vector(0.8f, 0.1f, -0.6f).normalized();
  for (const float angle : {0.0f, 0.1f, 0.25f, 0.26f, 0.9f, 2.0f, 3.1f}) {
    HS_CONTEXT("angle milliradians", static_cast<int>(angle * 1000.0f));
    HS_EXPECT_VEC(Spacetime::rotate_by(p, axis * angle),
                  math::rotate(p, math::make_rotation(axis, angle)), 2e-6f);
  }
  constexpr int W = 96, H = 48;
  auto basis = std::make_unique<Scan::ColumnBasisTable<W, H>>();
  const BasisMotion motion = tilted_motion(0.2f);
  Scan::ColumnBasis time = basis->prepare(motion.previous, motion.current);
  time.drift = axis * 1.4f;
  ScratchScope scope(scratch_arena_a);
  Scan::DriftColumns<W> drifted(scratch_arena_a, time);
  for (int c = 0; c < W / 2; c += 7) {
    drifted.reach(c);
    const math::Quaternion rest = math::make_rotation(
        axis, 1.4f * (1.0f - Spacetime::SamplingSurface<W>::time(c)));
    HS_EXPECT_VEC(drifted.view.azimuth[c], math::rotate(time.azimuth[c], rest),
                  2e-6f);
    HS_EXPECT_VEC(drifted.view.up[c], math::rotate(time.up[c], rest), 2e-6f);
  }
}

/**
 * @brief A cap that the swept clip test rules out lights no pixel of
 *        it under the frame's rotation.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 */
template <int W, int H> inline void check_swept_cap_clip() {
  HS_CONTEXT(W == 96 ? "cap clip 96" : "cap clip 288");
  const ScopedPoleLod lod(0.0f);
  const BasisMotion motion = tilted_motion(0.2f);
  auto basis = std::make_unique<Scan::ColumnBasisTable<W, H>>();
  const Scan::ColumnBasis time =
      basis->prepare(motion.previous, motion.current);
  struct Band {
    int y0, y1, x0, x1;
  };
  const Band BANDS[] = {{0, H / 2, 0, W / 2},
                        {0, H / 2, W / 2, W},
                        {H / 2, H, 0, W / 2},
                        {H / 2, H, W / 2, W},
                        {H / 3, H / 3 + 2, W / 3, W / 3 + 3}};
  constexpr int SHAPES = 48;
  int ruled_out = 0;
  for (int i = 0; i < SHAPES; ++i) {
    HS_CONTEXT("shape", i);
    const math::Basis b =
        math::make_basis(math::Quaternion(), math::fib_spiral(SHAPES, 0.5f, i));
    const auto star_basis = math::get_antipode(b, 0.06f);
    const SDF::Star star(star_basis.first, star_basis.second, 5, 0.0f);
    const SDF::Ring disc(b, 0.0f, 0.05f);
    const auto check = [&](const auto &shape) {
      math::Vector axis;
      float radius = 0.0f;
      HS_EXPECT_TRUE(shape.bounding_cap(axis, radius));
      const std::vector<Pixel> lit = render_column_basis<W, H>(
          motion, [&](Pipeline<W, H> &pipeline, Canvas &canvas,
                      const math::Quaternion &) {
            Scan::rasterize<W, H, false>(pipeline, canvas, shape,
                                         basis_flat_shader);
          });
      for (const Band &band : BANDS) {
        hs_test::StubEffect fx(W, H);
        fx.set_clip(band.y0, band.y1, band.x0, band.x1);
        const ClipRegion &cr = fx.clip();
        const auto xc = cr.x_clip();
        bool reaches = false;
        for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y)
          for (int x = 0; x < W; ++x)
            reaches |= !xc.clipped(x) &&
                       !is_black(lit[static_cast<size_t>(y) * W + x]);
        const bool meets = Plot::swept_cap_may_touch_clip<W, H>(
            cr, axis, radius, time.sweep_arc);
        const auto cap =
            Plot::make_swept_cap_clip<W, H>(cr, radius, time.sweep_arc);
        HS_EXPECT_EQ(meets, Plot::cap_may_touch_clip<H>(cap, axis));
        const auto center = Plot::make_cap_center(cr, axis);
        HS_EXPECT_EQ(meets, (Plot::swept_cap_may_touch_clip<W, H>(
                                cr, center, radius, time.sweep_arc)));
        if (reaches)
          HS_EXPECT_TRUE(meets);
        ruled_out += !meets;
      }
    };
    check(star);
    check(disc);
  }
  HS_EXPECT_GT(ruled_out, SHAPES);
}

/** @brief Moved shape intervals cover their current center azimuth. */
inline void test_moved_proxy_intervals_follow_centre() {
  constexpr int W = 288, H = 144;
  const auto check = [](const auto &shape) {
    const auto proxy =
        shape.moved(math::Vector(0.0f, cosf(0.5f), sinf(0.5f)), 0.0f);
    bool covered = false;
    for (int y = 0; y < H; ++y)
      proxy.template get_horizontal_intervals<W, H>(y, [&](float lo, float hi) {
        covered |= lo <= W / 4.0f && hi >= W / 4.0f;
      });
    HS_EXPECT_TRUE(covered);
  };
  const auto basis = math::make_basis(math::Quaternion(), math::Y_AXIS);
  check(SDF::Ring(basis, 0.2f, 0.03f));
  const auto star_basis = math::get_antipode(basis, 0.2f);
  check(SDF::Star(star_basis.first, star_basis.second, 5, 0.0f));
}

/** @brief Swept cap culling keeps every clip a swept star or disc lights. */
inline void test_swept_cap_clip_keeps_lit_clips() {
  check_swept_cap_clip<96, 20>();
  check_swept_cap_clip<96, 48>();
  check_swept_cap_clip<288, 144>();
}

/**
 * @brief A disc with a birth fraction lights only the columns swept at or
 *        after it, on both arms.
 */
inline void test_point_birth_skips_earlier_columns() {
  constexpr int W = 96, H = 48;
  const math::Vector centre = math::Vector(0.0f, 0.2f, 1.0f).normalized();
  const auto draw = [&](const Scan::ColumnBasis *time) {
    return render_stub_frame<W, H>([&](auto &pipeline, Canvas &canvas) {
      if (time)
        canvas.set_column_basis(time);
      Scan::Point::draw<W, H>(pipeline, canvas, centre, 0.5f,
                              basis_flat_shader);
      Scan::Point::draw<W, H>(pipeline, canvas, -centre, 0.5f,
                              basis_flat_shader);
    });
  };
  const std::vector<Pixel> whole = draw(nullptr);
  Scan::ColumnBasis born;
  born.birth = 0.4f;
  const std::vector<Pixel> cut = draw(&born);
  int before = 0, after = 0;
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x) {
      const size_t i = static_cast<size_t>(y) * W + x;
      if (Spacetime::SamplingSurface<W>::time(x) < born.birth) {
        before += !is_black(whole[i]);
        HS_EXPECT_TRUE(is_black(cut[i]));
      } else {
        after += !is_black(whole[i]);
        HS_EXPECT_EQ(cut[i], whole[i]);
      }
    }
  HS_EXPECT_GT(before, 20);
  HS_EXPECT_GT(after, 20);
  HS_EXPECT_EQ(Spacetime::SamplingSurface<W>::first_column_at(0.4f), 20);
  HS_EXPECT_EQ(Spacetime::SamplingSurface<W>::first_column_at(0.0f), 0);
  HS_EXPECT_EQ(Spacetime::SamplingSurface<W>::first_column_at(1.0f), W / 2);
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
  const auto place = [&](float u) {
    const math::Quaternion q =
        math::slerp(camera.previous, camera.current, u).normalized();
    return std::pair{math::rotate(anchor, q), q * spin(u)};
  };

  auto columns = std::make_unique<Scan::VolumeColumns<W>>();
  for (int c = 0; c < W / 2; ++c) {
    const auto [centre, q] = place(Spacetime::SamplingSurface<W>::time(c));
    columns->set(c, centre, q);
  }
  const auto [centre, q] = place(1.0f);
  const std::vector<Pixel> swept =
      render_stub_frame<W, H>([&](auto &pipeline, Canvas &canvas) {
        Scan::TransformedVolume volume(torus, centre, q);
        Scan::Volume::draw<W, H>(pipeline, canvas, centre, BOUNDS, volume,
                                 shade, 24, 0.02f, *columns);
      });

  const auto expected =
      render_state_oracle<W, H>([&](auto &pipeline, Canvas &canvas, float u) {
        const auto [at, rotation] = place(u);
        Scan::TransformedVolume volume(torus, at, rotation);
        Scan::Volume::draw<W, H>(pipeline, canvas, at, BOUNDS, volume, shade,
                                 24, 0.02f);
      });
  HS_EXPECT_GT(expect_factored_tolerance("volume", swept, expected), 30);
}

/**
 * @brief A ring group with previous axes matches its rings rebuilt about the
 *        axes at each column's sweep time, by the interpolant each slot's turn
 *        bound picks: the normalized chord, or the great-circle arc.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 * @param turns Per-slot angle between the previous and current axes, slow
 *        slots taking the chord and fast ones the arc.
 */
template <int W, int H>
inline void check_ring_group_previous_axes(std::initializer_list<float> turns) {
  constexpr int N = 3;
  HS_CONTEXT(W == 96 ? "previous axes 96" : "previous axes 288");
  const ScopedPoleLod lod(0.0f);
  math::Vector previous[N], current[N];
  int s = 0;
  for (const float turn : turns) {
    current[s] = math::rotate(math::Vector(0.77f, 0.0f, -0.6f).normalized(),
                              math::make_rotation(math::X_AXIS, 0.3f * s));
    previous[s] = math::rotate(
        current[s], math::make_rotation(
                        math::Vector(0.6f, -0.2f, 0.77f).normalized(), turn));
    ++s;
  }
  const auto draw = [&](const math::Vector *axes, const math::Vector *before) {
    return render_stub_frame<W, H>([&](auto &pipeline, Canvas &canvas) {
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
    });
  };
  const std::vector<Pixel> swept = draw(current, previous);
  int arcs = 0;
  const auto expected = column_time_oracle<W, H>([&](float u) {
    math::Vector axes[N];
    for (int s = 0; s < N; ++s) {
      const math::Vector normal = math::cross(previous[s], current[s]);
      const float angle = atan2f(sqrtf(math::dot(normal, normal)),
                                 math::dot(previous[s], current[s]));
      const bool arc = Spacetime::choose_turn_interpolant(
                           angle, Scan::RingGroup::MAX_AXIS_ERROR_PX *
                                      math::RADIANS_PER_COLUMN<W>) ==
                       Spacetime::Interpolant::ARC;
      arcs += arc && u == 0.0f;
      axes[s] =
          arc ? math::rotate(previous[s], math::make_rotation(
                                              normal.normalized(), u * angle))
              : (previous[s] + (current[s] - previous[s]) * u).normalized();
    }
    return draw(axes, nullptr);
  });
  HS_EXPECT_GT(expect_factored_tolerance("previous axes", swept, expected),
               100);
  HS_EXPECT_GT(arcs, 0);
  HS_EXPECT_LT(arcs, N);
}

/**
 * @brief The normalized chord of a constant-rate turn stays within
 *        PathBound::normalized_chord_error() of the turn up to a quarter turn.
 */
inline void test_normalized_chord_error_bounds_the_turn() {
  for (int i = 1; i <= 64; ++i) {
    const float angle = 0.5f * math::PI_F * static_cast<float>(i) / 64.0f;
    const double bound =
        Spacetime::PathBound::rotation(angle).normalized_chord_error();
    double worst = 0.0;
    for (int k = 0; k <= 256; ++k) {
      const double u = k / 256.0;
      const double lead =
          std::atan2(u * std::sin(angle), 1.0 - u + u * std::cos(angle));
      worst = std::max(worst, std::fabs(lead - u * angle));
    }
    HS_CONTEXT("angle milliradians", static_cast<int>(angle * 1000.0f));
    HS_EXPECT_LE(worst, bound);
    HS_EXPECT_GE(worst, 0.8 * bound);
  }
}

/**
 * @brief Ring group slots follow their axes' great circles at slow and fast
 *        turns.
 */
inline void test_ring_group_previous_axes_match_column_rings() {
  check_ring_group_previous_axes<96, 48>({0.25f, 0.35f, 0.45f});
  check_ring_group_previous_axes<288, 144>({0.02f, 0.2f, 1.2f});
}

/** @brief Column row spans contain every sample inside the cosine band. */
inline void test_ring_group_column_rows_cover_band() {
  constexpr int W = 288, H = 144;
  math::TrigLUT<W, H>::init();
  hs::Pcg32 rng(8291u);
  for (int i = 0; i < 2000; ++i) {
    const float a = i == 0 ? 0.0f : hs_test::rand_uniform(rng, -1.0f, 1.0f);
    const float b = i == 0 ? 0.0f : hs_test::rand_uniform(rng, -1.0f, 1.0f);
    const float low = hs_test::rand_uniform(rng, -1.1f, 0.9f);
    const float high = low + hs_test::rand_uniform(rng, 0.001f, 1.0f);
    const auto band = Scan::RingGroup::column_band(a, b, low, high);
    const bool far_arm = (i & 1) != 0;
    const float signed_a = far_arm ? -a : a;
    const auto rows =
        Scan::RingGroup::column_rows<W, H>(a, b, low, high, band, far_arm);
    HS_EXPECT_LE(rows.count, 2);
    for (int span = 0; span < rows.count; ++span) {
      HS_EXPECT_GE(rows.lo[span], 0);
      HS_EXPECT_LT(rows.hi[span], H);
      HS_EXPECT_LE(rows.lo[span], rows.hi[span]);
      for (int y : {rows.lo[span], rows.hi[span]}) {
        const float d = signed_a * math::TrigLUT<W, H>::sin_phi[y] +
                        b * math::TrigLUT<W, H>::cos_phi[y];
        HS_EXPECT_GE(d, low);
        HS_EXPECT_LE(d, high);
      }
      if (span)
        HS_EXPECT_GT(rows.lo[span], rows.hi[span - 1] + 1);
    }
    for (int y = 0; y < H; ++y) {
      const float d = signed_a * math::TrigLUT<W, H>::sin_phi[y] +
                      b * math::TrigLUT<W, H>::cos_phi[y];
      if (d < low || d > high)
        continue;
      bool covered = false;
      for (int span = 0; span < rows.count; ++span)
        covered |= y >= rows.lo[span] && y <= rows.hi[span];
      HS_EXPECT_TRUE(covered);
    }
  }
}

/** @brief Windowed column row spans cover the window as the full spans do. */
inline void test_ring_group_column_rows_window_matches_full() {
  constexpr int W = 288, H = 144;
  math::TrigLUT<W, H>::init();
  hs::Pcg32 rng(5113u);
  int dropped = 0;
  for (int i = 0; i < 2000; ++i) {
    const float a = hs_test::rand_uniform(rng, -1.0f, 1.0f);
    const float b = hs_test::rand_uniform(rng, -1.0f, 1.0f);
    const float low = hs_test::rand_uniform(rng, -1.1f, 0.9f);
    const float high = low + hs_test::rand_uniform(rng, 0.001f, 0.2f);
    const auto band = Scan::RingGroup::column_band(a, b, low, high);
    const bool far_arm = (i & 1) != 0;
    const int y_lo = static_cast<int>(rng() % H);
    const int y_hi =
        y_lo + static_cast<int>(rng() % static_cast<uint32_t>(H - y_lo));
    const auto full =
        Scan::RingGroup::column_rows<W, H>(a, b, low, high, band, far_arm);
    const auto window = Scan::RingGroup::column_rows<W, H>(
        a, b, low, high, band, far_arm, y_lo, y_hi);
    HS_EXPECT_LE(window.count, full.count);
    dropped += full.count - window.count;
    const auto covers = [](const Scan::RingGroup::ColumnRows &rows, int y) {
      bool covered = false;
      for (int span = 0; span < rows.count; ++span)
        covered |= y >= rows.lo[span] && y <= rows.hi[span];
      return covered;
    };
    for (int y = y_lo; y <= y_hi; ++y)
      HS_EXPECT_EQ(covers(window, y), covers(full, y));
  }
  HS_EXPECT_GT(dropped, 0);
}

/**
 * @brief A column basis prepared from an orientation follows its sub-frame
 *        history: each column's entries turn the fixed meridian by the
 *        current pose times the inverse of the history's pose at the
 *        column's sweep time, within the sweep bound.
 */
inline void test_column_basis_follows_history_path() {
  constexpr int W = 288, H = 144, HALF = W / 2;
  math::Orientation<> orientation(math::make_rotation(math::X_AXIS, 0.3f));
  orientation.push(math::make_rotation(math::Y_AXIS, 0.15f) *
                   orientation.get());
  orientation.push(math::make_rotation(math::Z_AXIS, 0.15f) *
                   orientation.get());
  orientation.push(math::make_rotation(math::X_AXIS, -0.15f) *
                   orientation.get());
  ScratchScope scope(scratch_arena_a);
  const Scan::ColumnBasis time =
      Scan::prepare_column_basis<W, H>(scratch_arena_a, orientation);
  HS_EXPECT_TRUE(time.azimuth != nullptr);
  if (!time.azimuth)
    return;
  if (!math::TrigLUT<W, H>::initialized)
    math::TrigLUT<W, H>::init();
  float worst = 0.0f, endpoint = 0.0f, farthest = 0.0f;
  for (int c = 0; c < HALF; ++c) {
    const float u = Spacetime::SamplingSurface<W>::time(c);
    const math::Vector meridian(math::TrigLUT<W, H>::cos_theta(c), 0.0f,
                                math::TrigLUT<W, H>::sin_theta[c]);
    const math::Quaternion delta =
        orientation.get() * orientation.sample(u).conjugate();
    const math::Vector up = math::rotate(math::Y_AXIS, delta);
    worst = fmaxf(worst, (time.up[c] - up).magnitude());
    worst = fmaxf(
        worst, (time.azimuth[c] - math::rotate(meridian, delta)).magnitude());
    const math::Quaternion slerped =
        orientation.get() *
        math::slerp(orientation.get(0), orientation.get(), u).conjugate();
    endpoint =
        fmaxf(endpoint, (math::rotate(math::Y_AXIS, slerped) - up).magnitude());
    farthest =
        fmaxf(farthest,
              2.0f * asinf(fminf(sqrtf(math::dot(delta.v, delta.v)), 1.0f)));
  }
  HS_EXPECT_LT(worst, 2e-6f);
  HS_EXPECT_GT(endpoint, 100.0f * 2e-6f);
  HS_EXPECT_LE(farthest, time.sweep_arc);
}

/**
 * @brief An object-frame walk hands each pixel its center and sub-samples
 *        rotated into the object frame at its column's sweep time, on both
 *        arms, and visits only the columns the render-x window admits.
 */
inline void test_object_column_basis_matches_rotated_meridian() {
  constexpr int W = 32, H = 12, HALF = W / 2;
  const BasisMotion motion = tilted_motion();
  auto basis = std::make_unique<Scan::ObjectBasisTable<W, H>>();
  basis->prepare(motion.previous, motion.current);
  if (!math::TrigLUT<W, H>::initialized)
    math::TrigLUT<W, H>::init();
  for (const Spacetime::TimeWindow window :
       {Spacetime::TimeWindow{}, Spacetime::TimeWindow{3, 9, true}}) {
    HS_CONTEXT("window", window.begin, window.end);
    hs_test::StubEffect fx(W, H);
    std::vector<int> visits(static_cast<size_t>(W) * H, 0);
    float worst = 0.0f;
    {
      Canvas canvas(fx);
      std::optional<TimeWindowScope> scope;
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
                            Spacetime::SamplingSurface<W>::time(x))
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

/**
 * @brief A fused ring stack scanned under a render-x window draws exactly the
 *        whole frame's pixels in the window's columns on both arms and none
 *        elsewhere, at the fixed meridian and through the rotated basis.
 */
inline void test_ring_stack_honours_render_window() {
  constexpr int W = 96, H = 48, HALF = W / 2;
  constexpr int RINGS = 3, LUT_N = 32;
  static float knots[RINGS][LUT_N + 1];
  build_ring_knots<RINGS, LUT_N>(knots);
  const BasisMotion motion = tilted_motion();
  auto basis = std::make_unique<Scan::ColumnBasisTable<W, H>>();
  const Scan::ColumnBasis swept =
      basis->prepare(motion.previous, motion.current);
  const auto render = [&](bool rotated, const Spacetime::TimeWindow &window) {
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H> pipeline;
    {
      Canvas canvas(fx);
      const Scan::ColumnBasis fixed;
      Scan::ColumnBasisScope time(canvas, rotated ? swept : fixed);
      std::optional<TimeWindowScope> scope;
      if (window.active)
        scope.emplace(canvas, window.begin, window.end);
      const math::Basis frame =
          math::make_basis(motion.current, math::Vector(0.3f, 0.2f, 1.0f));
      DistortedRingArray<RINGS> rings;
      auto *shapes = rings.data();
      int8_t slots[RINGS];
      for (int i = 0; i < RINGS; ++i) {
        rings.emplace(i, frame, 0.4f + 0.3f * i, 0.06f, knots[i], LUT_N, 0.0f,
                      nullptr);
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
    }
    fx.advance_display();
    std::vector<Pixel> out;
    capture_frame<W, H>(fx, out);
    return out;
  };
  for (const bool rotated : {false, true}) {
    HS_CONTEXT("rotated", rotated);
    const std::vector<Pixel> whole = render(rotated, Spacetime::TimeWindow{});
    for (const Spacetime::TimeWindow window :
         {Spacetime::TimeWindow{0, HALF / 2, true},
          Spacetime::TimeWindow{HALF / 2, HALF, true},
          Spacetime::TimeWindow{5, 17, true}}) {
      HS_CONTEXT("window", window.begin, window.end);
      const std::vector<Pixel> part = render(rotated, window);
      int lit = 0;
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x) {
          const size_t i = static_cast<size_t>(y) * W + x;
          const bool admitted =
              x % HALF >= window.begin && x % HALF < window.end;
          HS_EXPECT_EQ(part[i], (admitted ? whole[i] : Pixel(0, 0, 0)));
          lit += admitted && !is_black(whole[i]);
        }
      HS_EXPECT_GT(lit, 0);
    }
  }
}

/**
 * @brief Plot::cap_may_touch_window() admits every cap holding a pixel of the
 *        window's columns on either arm, and rejects some that hold none.
 */
inline void test_cap_may_touch_window_covers_window_pixels() {
  constexpr int W = 96, H = 48, HALF = W / 2;
  if (!math::TrigLUT<W, H>::initialized)
    math::TrigLUT<W, H>::init();
  int rejected = 0;
  for (const Spacetime::TimeWindow window :
       {Spacetime::TimeWindow{0, HALF / 2, true},
        Spacetime::TimeWindow{HALF / 2, HALF, true},
        Spacetime::TimeWindow{7, 9, true}}) {
    HS_CONTEXT("window", window.begin, window.end);
    for (int trial = 0; trial < 13 * 17 * 3; ++trial) {
      const float colatitude = math::PI_F * (trial % 13 + 0.5f) / 13.0f;
      const float longitude = math::TWO_PI_F * (trial / 13 % 17 + 0.3f) / 17.0f;
      const math::Vector dir(sinf(colatitude) * cosf(longitude),
                             cosf(colatitude),
                             sinf(colatitude) * sinf(longitude));
      const float reach = 0.05f + 0.2f * (trial / (13 * 17));
      bool holds = false;
      for (int y = 0; y < H && !holds; ++y)
        for (int x = 0; x < W && !holds; ++x) {
          if (x % HALF < window.begin || x % HALF >= window.end)
            continue;
          const float sp = math::TrigLUT<W, H>::sin_phi[y];
          const math::Vector pixel(sp * math::TrigLUT<W, H>::cos_theta(x),
                                   math::TrigLUT<W, H>::cos_phi[y],
                                   sp * math::TrigLUT<W, H>::sin_theta[x]);
          holds = math::dot(pixel, dir) > cosf(reach);
        }
      const bool admitted =
          Plot::cap_may_touch_window<W>(window, dir, reach, sinf(reach));
      HS_CONTEXT("trial", trial);
      if (holds)
        HS_EXPECT_TRUE(admitted);
      rejected += !admitted;
    }
  }
  HS_EXPECT_GT(rejected, 100);
}

// ============================================================================
// Column-time state: shaders, drift and forced tables at each column's time
// ============================================================================

/**
 * @brief Colour of the oracle shaders at a frame fraction.
 * @param u Frame fraction.
 * @return Red rising with @p u, alpha from 0.4 to 1.
 */
inline Color4 column_tint(float u) {
  return Color4(
      Pixel(static_cast<uint16_t>(2000.0f + 60000.0f * u), 18000, 9000),
      0.4f + 0.6f * u);
}

/**
 * @brief Checks that Scan::rasterize runs a time-aware shader at each
 *        column's sweep time.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 * @param lod Pole-LOD aggressiveness held for the check.
 */
template <int W, int H> inline void check_timed_shader(float lod) {
  HS_CONTEXT("width", static_cast<long long>(W));
  const ScopedPoleLod scoped(lod);
  const math::Vector discs[] = {math::Vector(1.0f, 0.1f, 0.05f).normalized(),
                                math::Vector(-1.0f, -0.2f, 0.05f).normalized(),
                                math::Vector(0.1f, 1.0f, 0.1f).normalized()};
  const auto antipode = math::get_antipode(
      math::make_basis(math::Quaternion(), math::Vector(0.05f, -1.0f, 0.1f)),
      0.2f);
  const SDF::Star star(antipode.first, antipode.second, 5, 0.0f, false);
  const auto timed = [&](Pipeline<W, H> &pipeline, Canvas &canvas) {
    const auto shader = [](const math::Vector &, Fragment &f, float u) {
      f.color = column_tint(u);
    };
    for (const math::Vector &centre : discs)
      Scan::Point::draw<W, H>(pipeline, canvas, centre, 0.3f, shader);
    Scan::rasterize<W, H>(pipeline, canvas, star, TimedShaderFn(shader));
  };
  const auto held = [&](Pipeline<W, H> &pipeline, Canvas &canvas, float u) {
    const auto shader = [u](const math::Vector &, Fragment &f) {
      f.color = column_tint(u);
    };
    for (const math::Vector &centre : discs)
      Scan::Point::draw<W, H>(pipeline, canvas, centre, 0.3f, shader);
    Scan::rasterize<W, H>(pipeline, canvas, star, FragmentShaderFn(shader));
  };
  const std::vector<Pixel> drawn = render_stub_frame<W, H>(timed);
  const std::vector<Pixel> expected = render_state_oracle<W, H>(held);
  int lit = 0, mismatched = 0;
  for (size_t i = 0; i < expected.size(); ++i) {
    lit += !is_black(expected[i]);
    mismatched += !(drawn[i] == expected[i]);
  }
  HS_EXPECT_EQ(mismatched, 0);
  HS_EXPECT_GT(lit, 300);
  const size_t seam = static_cast<size_t>(H / 2) * W + W / 2;
  HS_EXPECT_GT(drawn[seam - 1].r, drawn[seam].r);

  const std::vector<Pixel> off = [&] {
    const ScopedTemporalSampling sampling(false);
    return render_stub_frame<W, H>(timed);
  }();
  const std::vector<Pixel> endpoint =
      render_stub_frame<W, H>([&](Pipeline<W, H> &pipeline, Canvas &canvas) {
        held(pipeline, canvas, 1.0f);
      });
  HS_EXPECT_TRUE(off == endpoint);
}

/**
 * @brief Scan::rasterize runs a time-aware shader at each column's sweep
 *        time, per column through pole-LOD blocks and across both arm
 *        boundaries, and at the current endpoint with temporal sampling off.
 */
inline void test_timed_shader_reads_column_sweep_time() {
  check_timed_shader<96, 48>(0.0f);
  check_timed_shader<288, 144>(0.6f);
}

/**
 * @brief Checks a face palette shader's column-time motion against the
 *        shader held at each column's sweep time.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 * @param lod Pole-LOD aggressiveness held for the check.
 */
template <int W, int H> inline void check_face_palette_motion(float lod) {
  HS_CONTEXT("width", static_cast<long long>(W));
  const ScopedPoleLod scoped(lod);
  constexpr float BLEND_FROM = 0.15f, BLEND_TO = 0.85f;
  constexpr float RIM_FROM = 0.35f, RIM_TO = 0.9f;
  constexpr float ALPHA_FROM = 0.25f, ALPHA_TO = 1.0f;
  ScratchScope scope(scratch_arena_a);
  BakedPaletteStorage from, to, rim;
  from.bake(scratch_arena_a, Palettes::UNDERSEA);
  to.bake(scratch_arena_a, Palettes::MANGO_PEEL);
  rim.bake(scratch_arena_a, Palettes::RICH_SUNSET);
  MeshState mesh;
  build_basis_octahedron(mesh, scratch_arena_a);
  const auto lerp = [](float a, float b, float u) { return a + (b - a) * u; };
  const auto draw = [&](Pipeline<W, H> &pipeline, Canvas &canvas, float blend,
                        float previous_blend, float rim_weight,
                        float previous_rim, float alpha, float previous_alpha) {
    ScratchScope blends(scratch_arena_b);
    const BakedPalette ramp =
        bake_palette_blend(scratch_arena_b, from, to, blend);
    FacePaletteShader shader;
    auto setup = [&](size_t, float size) {
      shader.set_palette(&ramp);
      if (blend != previous_blend)
        shader.set_blend(&from.view(), &to.view(), blend, previous_blend);
      shader.set_counterpart(&rim.view(), rim_weight, previous_rim);
      shader.set_alpha(alpha, previous_alpha);
      shader.scale = size > math::TOLERANCE ? 1.0f / size : 0.0f;
    };
    Scan::Mesh::draw_specialized<W, H>(pipeline, canvas, mesh, shader,
                                       scratch_arena_b, nullptr, setup);
  };
  const std::vector<Pixel> moving =
      render_stub_frame<W, H>([&](Pipeline<W, H> &pipeline, Canvas &canvas) {
        auto basis = std::make_unique<Scan::ColumnBasisTable<W, H>>();
        const Scan::ColumnBasis time =
            basis->prepare(math::Quaternion(), math::Quaternion(), true);
        Scan::ColumnBasisScope column_scope(canvas, time);
        draw(pipeline, canvas, BLEND_TO, BLEND_FROM, RIM_TO, RIM_FROM, ALPHA_TO,
             ALPHA_FROM);
      });
  const std::vector<Pixel> expected = render_state_oracle<W, H>(
      [&](Pipeline<W, H> &pipeline, Canvas &canvas, float u) {
        const float blend = lerp(BLEND_FROM, BLEND_TO, u);
        const float rim_weight = lerp(RIM_FROM, RIM_TO, u);
        const float alpha = lerp(ALPHA_FROM, ALPHA_TO, u);
        draw(pipeline, canvas, blend, blend, rim_weight, rim_weight, alpha,
             alpha);
      });
  HS_EXPECT_GT(expect_frames_within("face palette shader", moving, expected, 2),
               1000);
}

/**
 * @brief A face palette shader under column time moves its alpha, crossfade
 *        and counterpart weight from their previous endpoints to the current
 *        ones across each arm, per column through pole-LOD blocks.
 */
inline void test_face_palette_shader_moves_at_column_time() {
  constexpr size_t SCRATCH_A = DEFAULT_SCRATCH_A_SIZE;
  constexpr size_t SCRATCH_B = DEFAULT_SCRATCH_B_SIZE +
                               sizeof(SDF::FaceStorage) -
                               sizeof(SDF::FaceScratchBuffer);
  hs_test::ScopedArenaSplit split(GLOBAL_ARENA_SIZE - SCRATCH_A - SCRATCH_B,
                                  SCRATCH_A, SCRATCH_B);
  check_face_palette_motion<96, 48>(0.0f);
  check_face_palette_motion<288, 144>(0.6f);
}

/**
 * @brief A volume's time-aware shader receives each hit's column sweep time.
 */
inline void test_volume_timed_shader_reads_column_sweep_time() {
  constexpr int W = 96, H = 48;
  const SDF::Torus torus{0.3f, 0.1f};
  const math::Vector centre = math::Vector(-1.0f, 0.2f, 0.1f).normalized();
  const math::Quaternion q =
      math::make_rotation(math::Vector(0.3f, 1.0f, 0.2f).normalized(), 0.8f);
  constexpr float BOUNDS = 0.42f;
  auto columns = std::make_unique<Scan::VolumeColumns<W>>();
  for (int c = 0; c < W / 2; ++c)
    columns->set(c, centre, q);
  const std::vector<Pixel> drawn =
      render_stub_frame<W, H>([&](Pipeline<W, H> &pipeline, Canvas &canvas) {
        Scan::TransformedVolume volume(torus, centre, q);
        Scan::Volume::draw<W, H>(
            pipeline, canvas, centre, BOUNDS, volume,
            [](const math::Vector &, Fragment &f, float u) {
              f.color = Color4(column_tint(u).color, 1.0f);
            },
            24, 0.02f, *columns);
      });
  const std::vector<Pixel> expected = render_state_oracle<W, H>(
      [&](Pipeline<W, H> &pipeline, Canvas &canvas, float u) {
        Scan::TransformedVolume volume(torus, centre, q);
        Scan::Volume::draw<W, H>(
            pipeline, canvas, centre, BOUNDS, volume,
            [u](const math::Vector &, Fragment &f) {
              f.color = Color4(column_tint(u).color, 1.0f);
            },
            24, 0.02f, *columns);
      });
  HS_EXPECT_GT(expect_factored_tolerance("timed volume", drawn, expected), 50);
}

/**
 * @brief surface_column inverts surface_time on both arms.
 */
inline void test_surface_column_inverts_surface_time() {
  for (int width = 4; width <= 512; width += 2) {
    HS_CONTEXT("width", static_cast<long long>(width));
    for (int x = 0; x < width; ++x)
      HS_EXPECT_EQ(
          Spacetime::surface_column(Spacetime::surface_time(x, width), width),
          x % (width / 2));
  }
}

/**
 * @brief Moving per-column state prepares a table without rotation, and only
 *        with temporal sampling on.
 */
inline void test_column_basis_prepares_for_moving_state() {
  constexpr int W = 96, H = 48;
  const math::Quaternion q =
      math::make_rotation(math::Vector(0.2f, 0.9f, 0.4f).normalized(), 0.6f);
  auto basis = std::make_unique<Scan::ColumnBasisTable<W, H>>();
  HS_EXPECT_TRUE(basis->prepare(q, q).azimuth == nullptr);
  const Scan::ColumnBasis forced = basis->prepare(q, q, true);
  HS_EXPECT_TRUE(forced.azimuth != nullptr);
  HS_EXPECT_TRUE(forced.sweep_arc > 0.0f && forced.sweep_arc < 1e-3f);
  {
    const ScopedTemporalSampling sampling(false);
    HS_EXPECT_TRUE(basis->prepare(q, q, true).azimuth == nullptr);
  }
  ScratchScope scope(scratch_arena_a);
  const Scan::ColumnBasis still =
      Scan::prepare_column_basis<W, H>(scratch_arena_a, q, q);
  const Scan::ColumnBasis moving =
      Scan::prepare_column_basis<W, H>(scratch_arena_a, q, q, true);
  HS_EXPECT_TRUE(still.azimuth == nullptr);
  HS_EXPECT_TRUE(moving.azimuth != nullptr);
}

/// Largest channel error between a column-time scan and its per-column
/// oracle under a multi-interval pose history.
inline constexpr int HISTORY_ORACLE_CHANNEL_ERROR = 64;

/**
 * @brief Re-centred caps keep every pixel a per-column oracle lights under a
 *        pose history of unequal or reversing interval turns.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 * @param label History named in failure context.
 * @param turns Turn of each history pose about one tilted axis (radians),
 *        oldest first.
 */
template <int W, int H>
inline void check_history_lit_caps(const char *label,
                                   const std::array<float, 4> &turns) {
  HS_CONTEXT(label, W);
  const ScopedPoleLod lod(0.0f);
  const math::Vector axis = math::Vector(0.6f, -0.2f, 0.77f).normalized();
  math::Orientation<> orientation(math::make_rotation(axis, turns[0]));
  for (size_t i = 1; i < turns.size(); ++i)
    orientation.push(math::make_rotation(axis, turns[i]));
  ScratchScope scope(scratch_arena_a);
  const Scan::ColumnBasis time =
      Scan::prepare_column_basis<W, H>(scratch_arena_a, orientation);
  HS_EXPECT_TRUE(time.azimuth != nullptr);
  if (!time.azimuth)
    return;

  constexpr int HALF = W / 2;
  float travel = 0.0f;
  for (int c = 0; c < HALF; ++c)
    for (int m = 0; m < HALF; ++m) {
      const math::Quaternion step =
          orientation.sample(Spacetime::surface_time(c, W)) *
          orientation.sample(Spacetime::surface_time(m, W)).conjugate();
      const float turn =
          2.0f * asinf(fminf(sqrtf(math::dot(step.v, step.v)), 1.0f));
      travel =
          fmaxf(travel, turn - time.sweep_rate *
                                   static_cast<float>(std::abs(c - m)) / HALF);
    }
  HS_EXPECT_LE(travel, (Scan::ColumnBasisTable<W, H>::SWEEP_PAD));

  constexpr int SHAPES = 160;
  std::vector<math::Basis> bases;
  int moved = 0;
  for (int i = 0; i < SHAPES; ++i) {
    bases.push_back(math::make_basis(math::Quaternion(),
                                     math::fib_spiral(SHAPES, 0.5f, i)));
    Scan::LitCap cap;
    moved += Scan::lit_cap<W, H>(bases.back().v, 0.03f, time, cap);
  }
  HS_EXPECT_GT(moved, SHAPES / 2);
  const auto draw = [&](Pipeline<W, H> &pipeline, Canvas &canvas,
                        const math::Quaternion &placement, bool lit) {
    const auto place = [&](const math::Basis &b) {
      return math::Basis{math::rotate(b.u, placement),
                         math::rotate(b.v, placement),
                         math::rotate(b.w, placement)};
    };
    for (int i = 0; i < SHAPES; ++i) {
      const auto star_basis = math::get_antipode(place(bases[i]), 0.04f);
      const SDF::Star star(star_basis.first, star_basis.second, 5, 0.0f);
      const SDF::Ring disc(place(bases[(i * 7) % SHAPES]), 0.0f, 0.03f);
      if (lit) {
        Scan::rasterize<W, H, false>(pipeline, canvas, star, basis_flat_shader);
        Scan::rasterize<W, H, false>(pipeline, canvas, disc, basis_flat_shader);
      } else {
        Scan::rasterize<W, H, false>(pipeline, canvas, SweptOnlyStar(star),
                                     basis_flat_shader);
        Scan::rasterize<W, H, false>(pipeline, canvas, SweptOnlyRing(disc),
                                     basis_flat_shader);
      }
    }
  };
  const auto scanned = [&](bool lit) {
    return render_stub_frame<W, H>(
        [&](Pipeline<W, H> &pipeline, Canvas &canvas) {
          Scan::ColumnBasisScope column_scope(canvas, time);
          draw(pipeline, canvas, math::Quaternion(), lit);
        });
  };
  const std::vector<Pixel> actual = scanned(true);
  const std::vector<Pixel> swept = scanned(false);
  const std::vector<Pixel> expected = render_state_oracle<W, H>(
      [&](Pipeline<W, H> &pipeline, Canvas &canvas, float u) {
        draw(pipeline, canvas,
             (orientation.sample(u) * orientation.get().conjugate())
                 .normalized(),
             true);
      });
  int differing = 0, lit = 0, worst = 0;
  for (size_t i = 0; i < expected.size(); ++i) {
    differing += !(actual[i] == swept[i]);
    lit += !is_black(expected[i]);
    const Pixel &a = actual[i];
    const Pixel &b = expected[i];
    worst = std::max<int>(worst, max_channel_gap(a, b));
  }
  HS_EXPECT_EQ(differing, 0);
  HS_EXPECT_LE(worst, HISTORY_ORACLE_CHANNEL_ERROR);
  HS_EXPECT_GT(lit, 4 * SHAPES);
}

/**
 * @brief Compact shapes keep their oracle pixels under pose histories whose
 *        interval turns differ or reverse; an endpoint pair's rate bound is
 *        its arc.
 */
inline void test_lit_caps_follow_history_rate() {
  const BasisMotion motion = tilted_motion();
  auto basis = std::make_unique<Scan::ColumnBasisTable<96, 48>>();
  const Scan::ColumnBasis endpoints =
      basis->prepare(motion.previous, motion.current);
  HS_EXPECT_EQ(endpoints.sweep_rate, endpoints.sweep_arc);
  for (const auto &[label, turns] :
       {std::pair{"variable rate", std::array{0.0f, 0.6f, 0.6f, 0.6f}},
        std::pair{"reversing", std::array{0.0f, 0.9f, 0.9f, 0.3f}}}) {
    check_history_lit_caps<96, 48>(label, turns);
    check_history_lit_caps<288, 144>(label, turns);
  }
}

/**
 * @brief Drifting discs and stars match the shapes rotated to where they
 *        stood at each column's sweep time, across the arm boundary, for
 *        short and long drifts.
 */
inline void test_drift_matches_shape_at_column_centre() {
  constexpr int W = 96, H = 48;
  const ScopedPoleLod lod(0.0f);
  const math::Vector centres[] = {
      math::Vector(-1.0f, 0.15f, 0.1f).normalized(),
      math::Vector(0.3f, -0.4f, 0.85f).normalized()};
  const math::Vector axis = math::Vector(0.2f, 1.0f, 0.1f).normalized();
  for (const float angle : {3.0f * math::RADIANS_PER_COLUMN<W>, 0.7f}) {
    HS_CONTEXT("drift milliradians", static_cast<int>(angle * 1000.0f));
    const auto drifting = [&](auto &&shape) {
      return render_stub_frame<W, H>(
          [&](Pipeline<W, H> &pipeline, Canvas &canvas) {
            auto basis = std::make_unique<Scan::ColumnBasisTable<W, H>>();
            Scan::ColumnBasis time =
                basis->prepare(math::Quaternion(), math::Quaternion(), true);
            Scan::ColumnBasisScope column_scope(canvas, time);
            time.drift = axis * angle;
            for (const math::Vector &c : centres)
              shape(pipeline, canvas, math::make_basis(math::Quaternion(), c));
          });
    };
    const auto oracle = [&](auto &&shape) {
      return render_state_oracle<W, H>(
          [&](Pipeline<W, H> &pipeline, Canvas &canvas, float u) {
            const math::Quaternion back =
                math::make_rotation(axis, -angle * (1.0f - u));
            for (const math::Vector &c : centres) {
              const math::Basis basis = math::make_basis(math::Quaternion(), c);
              shape(pipeline, canvas,
                    math::Basis{math::rotate(basis.u, back),
                                math::rotate(basis.v, back),
                                math::rotate(basis.w, back)});
            }
          });
    };
    const auto disc = [](Pipeline<W, H> &pipeline, Canvas &canvas,
                         const math::Basis &basis) {
      Scan::Point::draw<W, H>(pipeline, canvas, basis.v, 0.2f,
                              basis_flat_shader);
    };
    const auto star = [](Pipeline<W, H> &pipeline, Canvas &canvas,
                         const math::Basis &basis) {
      Scan::Star::draw<W, H, false>(pipeline, canvas, basis, 0.15f, 5,
                                    basis_flat_shader);
    };
    HS_EXPECT_GT(expect_factored_tolerance("drifting disc", drifting(disc),
                                           oracle(disc)),
                 40);
    HS_EXPECT_GT(expect_factored_tolerance("drifting star", drifting(star),
                                           oracle(star)),
                 40);
  }
  const float angle = 3.0f * math::RADIANS_PER_COLUMN<W>;
  HS_EXPECT_NEAR(
      math::dot(Spacetime::rotate_by(centres[1], axis * angle),
                math::rotate(centres[1], math::make_rotation(axis, angle))),
      1.0f, 1e-6f);
}

/**
 * @brief A star whose drift barely changes across the columns that see it
 *        settles to one turn that matches its per-column drift scan; one
 *        seen across an arm boundary does not settle.
 */
inline void test_settled_drift_matches_column_drift() {
  constexpr int W = 96, H = 48;
  const auto settle_shape = &Scan::settle_drift<W, H>;
  const ScopedPoleLod lod(0.0f);
  const BasisMotion motion = tilted_motion(0.05f);
  const math::Vector axis = math::Vector(0.2f, 1.0f, 0.1f).normalized();
  const math::Vector drift = axis * (0.3f * math::RADIANS_PER_COLUMN<W>);
  const math::Vector centre = math::Vector(0.3f, -0.4f, 0.85f).normalized();
  const float RADIUS = 0.15f;
  const auto star = [&](Pipeline<W, H> &pipeline, Canvas &canvas,
                        const math::Basis &basis) {
    Scan::Star::draw<W, H, false>(pipeline, canvas, basis, RADIUS, 5,
                                  basis_flat_shader);
  };
  const auto render = [&](bool settle) {
    return render_stub_frame<W, H>(
        [&](Pipeline<W, H> &pipeline, Canvas &canvas) {
          auto basis = std::make_unique<Scan::ColumnBasisTable<W, H>>();
          Scan::ColumnBasis time =
              basis->prepare(motion.previous, motion.current, true);
          Scan::ColumnBasisScope column_scope(canvas, time);
          math::Basis shape = math::make_basis(math::Quaternion(), centre);
          if (settle)
            HS_EXPECT_TRUE(
                settle_shape(time, drift, RADIUS * (math::PI_F / 2.0f), shape));
          else
            time.drift = drift;
          star(pipeline, canvas, shape);
        });
  };
  // Within DRIFT_SETTLE_ARC a star's edge moves under a tenth of a pixel, so
  // its coverage, and each channel, by under a tenth of full scale.
  const std::vector<Pixel> settled = render(true);
  const std::vector<Pixel> drifting = render(false);
  int lit = 0, worst = 0;
  for (size_t i = 0; i < settled.size(); ++i) {
    lit += !is_black(drifting[i]);
    worst = std::max<int>(worst, max_channel_gap(settled[i], drifting[i]));
  }
  std::printf("  [settled drift] lit=%d worst channel error=%d\n", lit, worst);
  HS_EXPECT_GT(lit, 20);
  HS_EXPECT_LE(worst, 52000 / 10);

  auto basis = std::make_unique<Scan::ColumnBasisTable<W, H>>();
  const Scan::ColumnBasis time =
      basis->prepare(motion.previous, motion.current, true);
  math::Basis seam = math::make_basis(
      math::Quaternion(), math::Vector(-1.0f, 0.1f, 0.02f).normalized());
  const math::Basis held = seam;
  HS_EXPECT_FALSE(
      settle_shape(time, drift, RADIUS * (math::PI_F / 2.0f), seam));
  HS_EXPECT_VEC(seam.u, held.u, 0.0f);
  HS_EXPECT_VEC(seam.v, held.v, 0.0f);
  HS_EXPECT_VEC(seam.w, held.w, 0.0f);
  math::Basis fast = math::make_basis(math::Quaternion(), centre);
  HS_EXPECT_FALSE(
      settle_shape(time, drift * 20.0f, RADIUS * (math::PI_F / 2.0f), fast));
}
