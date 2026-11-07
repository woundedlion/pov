/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

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
    Scan::rasterize<W, H>(pipeline, canvas, star, TimedShaderRef(shader));
  };
  const auto held = [&](Pipeline<W, H> &pipeline, Canvas &canvas, float u) {
    const auto shader = [u](const math::Vector &, Fragment &f) {
      f.color = column_tint(u);
    };
    for (const math::Vector &centre : discs)
      Scan::Point::draw<W, H>(pipeline, canvas, centre, 0.3f, shader);
    Scan::rasterize<W, H>(pipeline, canvas, star, FragmentShaderRef(shader));
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
  const Scan::ColumnBasis fixed =
      basis->prepare(math::Quaternion(), math::Quaternion(), true);
  const float THETA = 15 * math::TWO_PI_F / W;
  const float PHI = 0.4f;
  const math::Vector ODD_CENTRE{cosf(THETA) * sinf(PHI), cosf(PHI),
                                sinf(THETA) * sinf(PHI)};
  math::Basis odd = math::make_basis(math::Quaternion(), ODD_CENTRE);
  const math::Basis UNMOVED = odd;
  const math::Vector ODD_DRIFT = math::perpendicular_axis(ODD_CENTRE) * 0.033f;
  HS_EXPECT_FALSE(settle_shape(fixed, ODD_DRIFT, 0.045f, odd));
  HS_EXPECT_VEC(odd.u, UNMOVED.u, 0.0f);
  HS_EXPECT_VEC(odd.v, UNMOVED.v, 0.0f);
  HS_EXPECT_VEC(odd.w, UNMOVED.w, 0.0f);
}
