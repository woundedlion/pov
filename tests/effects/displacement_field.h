/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// DisplacementField phase machine, ball spans, span and hue-table bakes.
// ---------------------------------------------------------------------------

/**
 * @brief White-box accessor for DisplacementField's phase machine, ball pool,
 * span bake, and hue-table bake.
 */
struct DisplacementFieldWhiteBox {
  template <int W, int H>
  static void enter_balls(DisplacementField<W, H> &effect) {
    effect.enter_balls();
    effect.spawn_ball();
  }

  template <int W, int H>
  static void fill_ball_pool(DisplacementField<W, H> &effect) {
    effect.params.ball_speed_min = effect.params.ball_speed_max = 3.0f;
    effect.enter_balls();
    for (int i = 0; i <= effect.MAX_BALLS; ++i)
      effect.spawn_ball();
    HS_EXPECT_EQ(effect.balls.active_count(), effect.MAX_BALLS);
    HS_EXPECT_TRUE(effect.logged_pool_full);
    effect.ball_phase_left = 0;
  }

  template <int W, int H>
  static bool in_noise(const DisplacementField<W, H> &effect) {
    return effect.phase == DisplacementField<W, H>::Phase::NOISE;
  }

  template <int W, int H, size_t N>
  static void
  check_ball_spans(DisplacementField<W, H> &effect,
                   const std::array<Animation::BumpParams, N> &balls,
                   const math::Basis &basis, float theta) {
    for (size_t i = 0; i < N; ++i) {
      const auto &ball = balls[i];
      effect.ball_params[i] = &ball;
      effect.ball_local[i] = static_cast<int>(i);
      const float cv = math::dot(basis.v, ball.center);
      const float cu = math::dot(basis.u, ball.center);
      const float cw = math::dot(basis.w, ball.center);
      effect.ball_cv[i] = cv;
      effect.ball_rho[i] = sqrtf(cu * cu + cw * cw);
      effect.ball_colat[i] = math::fast_acos(hs::clamp(cv, -1.0f, 1.0f));
      const float azimuth = atan2f(cw, cu);
      effect.ball_azimuth[i] = azimuth < 0 ? azimuth + math::TWO_PI_F : azimuth;
    }
    std::array<float, W> shifts{};
    const float step = math::TWO_PI_F / W;
    effect.bake_ball_spans(basis, theta, cosf(theta), sinf(theta), cosf(step),
                           sinf(step), W, effect.CHUNK_MASK, N, shifts.data());
    int displaced = 0;
    for (int x = 0; x < W; ++x) {
      const float expected =
          effect.ball_field(effect.knot_pos[x], effect.ball_local, N, theta) +
          effect.noise_field.field(effect.knot_pos[x]);
      HS_EXPECT_NEAR(shifts[x], expected, 1e-6f);
      displaced += std::abs(expected) > 1e-5f;
    }
    HS_EXPECT_GT(displaced, 0);
  }

  template <int W, int H>
  static void prepare_hue_table(DisplacementField<W, H> &effect,
                                const Color4 &color, float domain) {
    effect.prepare_hue_table(make_hue_rotate_base(color), domain);
  }

  template <int W, int H>
  static Pixel sample_hue_table(const DisplacementField<W, H> &effect,
                                float amount, float domain, bool cyclic) {
    return effect.sample_hue_table(amount, domain, cyclic);
  }

  template <int W, int H>
  static Pixel sample_hue_table_cached(DisplacementField<W, H> &effect,
                                       float amount, float domain, bool cyclic,
                                       const HueRotateBase &base,
                                       uint64_t *valid) {
    return effect.sample_hue_table_cached(amount, domain, cyclic, base, valid);
  }

  template <int W, int H>
  static int hue_table_size(const DisplacementField<W, H> &) {
    return DisplacementField<W, H>::HUE_TABLE_SIZE;
  }

  template <int W, int H>
  static Pixel hue_table_value(const DisplacementField<W, H> &effect,
                               int index) {
    return effect.hue_table[index];
  }

  template <int W, int H>
  static void clear_hue_table(DisplacementField<W, H> &effect) {
    for (int i = 0; i <= DisplacementField<W, H>::HUE_TABLE_SIZE; ++i)
      effect.hue_table[i] = Pixel(0, 0, 0);
  }

  template <int W, int H>
  static void configure_noise(DisplacementField<W, H> &effect, float rings,
                              float hue_scale) {
    effect.params.num_rings = rings;
    effect.params.hue_scale = hue_scale;
    effect.master_gain = 1.0f;
  }

  /** @brief Drives the ring stack to its widest registered footprint. */
  template <int W, int H>
  static void configure_max_footprint(DisplacementField<W, H> &effect) {
    using Effect = DisplacementField<W, H>;
    effect.params.num_rings = static_cast<float>(Effect::RING_SLOTS);
    effect.params.thickness = 6.0f * Effect::THICKNESS_PX;
    effect.params.ball_amp = 0.8f;
    effect.params.noise_amp = 0.8f;
  }

  template <int W, int H>
  static void set_force_exact_hue(DisplacementField<W, H> &effect, bool exact) {
    effect.force_exact_hue = exact;
  }

  template <int W, int H>
  static int hue_table_uses(const DisplacementField<W, H> &effect) {
    return effect.hue_table_uses;
  }

  template <int W, int H>
  static Pixel hue_lut_value(const DisplacementField<W, H> &effect, int index) {
    // Slot 0 of the pooled bake: the first drawn ring's hue LUT.
    return effect.hue_pool[index];
  }

  template <int W, int H>
  static Color4 current_ring_color(const DisplacementField<W, H> &effect,
                                   float ring_t) {
    return effect.palette.get(math::wrap_t(ring_t + effect.color_spin));
  }

  template <int W, int H>
  static int baked_lut_samples(const DisplacementField<W, H> &effect) {
    return static_cast<int>(effect.slot_lut_nf[0]);
  }

  template <int W, int H>
  static int noise_lut_samples(const DisplacementField<W, H> &,
                               float feature_scale, float sin_theta) {
    return DisplacementField<W, H>::noise_lut_samples(feature_scale, sin_theta,
                                                      true);
  }

  template <int W, int H>
  static int octave_grid(const DisplacementField<W, H> &) {
    return DisplacementField<W, H>::OCTAVE_GRID;
  }

  /** @brief Bakes one fully visible ring through the octave-grid noise path. */
  template <int W, int H>
  static void bake_noise_octaves(DisplacementField<W, H> &effect,
                                 const Animation::NoiseProductParams &np,
                                 const math::Basis &basis, float theta,
                                 int lut_n, float *slut) {
    const float dphi = 2.0f * math::PI_F / lut_n;
    effect.bake_noise_octaves(np, basis, theta, cosf(theta), sinf(theta),
                              cosf(dphi), sinf(dphi), lut_n,
                              DisplacementField<W, H>::CHUNK_MASK, 0, slut);
  }
};

/**
 * @brief Verifies lazy table endpoints and outputs match eager construction.
 */
inline void test_displacement_field_lazy_hue_table_matches_eager() {
  reset_effect_globals();
  DisplacementField<SMALL_W, SMALL_H> effect;
  effect.init();
  GenerativePalette palette(PaletteRecipes::profile(PaletteDomain::MIRROR,
                                                    PaletteHarmony::ANALOGOUS,
                                                    AxisCurve::CONSTANT, 0.0f));
  const Color4 color = palette.get(0.37f);
  const HueRotateBase base = make_hue_rotate_base(color);
  struct TableCase {
    float domain;
    float max_amount;
    bool cyclic;
  };
  const TableCase cases[] = {{0.4f, 0.4f, false},
                             {-0.4f, -0.4f, false},
                             {1.0f, 3.0f, true},
                             {-1.0f, -3.0f, true}};
  constexpr int SAMPLE_COUNT = 1024;
  const int table_size = DisplacementFieldWhiteBox::hue_table_size(effect);

  for (const TableCase &table_case : cases) {
    DisplacementFieldWhiteBox::prepare_hue_table(effect, color,
                                                 table_case.domain);
    std::vector<Pixel> endpoints(table_size + 1);
    for (int i = 0; i <= table_size; ++i)
      endpoints[i] = DisplacementFieldWhiteBox::hue_table_value(effect, i);
    std::vector<Pixel> expected(SAMPLE_COUNT + 1);
    for (int i = 0; i <= SAMPLE_COUNT; ++i) {
      const float amount = table_case.max_amount * i / SAMPLE_COUNT;
      expected[i] = DisplacementFieldWhiteBox::sample_hue_table(
          effect, amount, table_case.domain, table_case.cyclic);
    }

    DisplacementFieldWhiteBox::clear_hue_table(effect);
    std::vector<uint64_t> valid((table_size + 64) / 64);
    for (int i = 0; i <= SAMPLE_COUNT; ++i) {
      const float amount = table_case.max_amount * i / SAMPLE_COUNT;
      Pixel actual = DisplacementFieldWhiteBox::sample_hue_table_cached(
          effect, amount, table_case.domain, table_case.cyclic, base,
          valid.data());
      HS_EXPECT_EQ(actual.r, expected[i].r);
      HS_EXPECT_EQ(actual.g, expected[i].g);
      HS_EXPECT_EQ(actual.b, expected[i].b);
    }
    for (int i = 0; i <= table_size; ++i)
      HS_EXPECT_TRUE(valid[i >> 6] & (uint64_t{1} << (i & 63)));
    for (int i = 0; i <= table_size; ++i) {
      Pixel actual = DisplacementFieldWhiteBox::hue_table_value(effect, i);
      HS_EXPECT_EQ(actual.r, endpoints[i].r);
      HS_EXPECT_EQ(actual.g, endpoints[i].g);
      HS_EXPECT_EQ(actual.b, endpoints[i].b);
    }
  }
}

/**
 * @brief Bounds dynamic and periodic hue tables over effect palette colors.
 * @details Each bound is the sweep's measured worst case plus headroom. The
 * sRGB8 bounds are looser because the encode is non-linear.
 */
inline void test_displacement_field_hue_table_fidelity() {
  reset_effect_globals();
  DisplacementField<SMALL_W, SMALL_H> effect;
  effect.init();

  constexpr float INV16 = 1.0f / 65535.0f;
  float default_delta_e = 0.0f;
  float cyclic_delta_e = 0.0f;
  int default_srgb_delta = 0;
  int cyclic_srgb_delta = 0;
  for (uint32_t seed = 0; seed < 12; ++seed) {
    GenerativePalette palette(PaletteRecipes::profile(
        PaletteDomain::MIRROR, PaletteHarmony::ANALOGOUS, AxisCurve::CONSTANT,
        PaletteRecipes::hue_turns(seed * 17u)));
    for (int color_index = 0; color_index < 48; ++color_index) {
      Color4 base = palette.get((color_index + 0.5f) / 48.0f);
      HueRotateBase exact_base = make_hue_rotate_base(base);
      for (int mode = 0; mode < 2; ++mode) {
        const bool cyclic = mode == 1;
        const float domain = cyclic ? 1.0f : 0.4f;
        const float max_amount = cyclic ? 3.0f : domain;
        DisplacementFieldWhiteBox::prepare_hue_table(effect, base, domain);
        for (int i = 0; i <= 1024; ++i) {
          const float amount = max_amount * i / 1024.0f;
          Pixel exact = hue_rotate(exact_base, amount).color;
          Pixel approx = DisplacementFieldWhiteBox::sample_hue_table(
              effect, amount, domain, cyclic);
          OKLab exact_lab = linear_rgb_to_oklab(
              exact.r * INV16, exact.g * INV16, exact.b * INV16);
          OKLab approx_lab = linear_rgb_to_oklab(
              approx.r * INV16, approx.g * INV16, approx.b * INV16);
          const float dl = exact_lab.L - approx_lab.L;
          const float da = exact_lab.a - approx_lab.a;
          const float db = exact_lab.b - approx_lab.b;
          const float delta_e = sqrtf(dl * dl + da * da + db * db);
          const int srgb_delta = std::max(
              std::abs(
                  static_cast<int>(linear_float_to_srgb8(exact.r * INV16)) -
                  linear_float_to_srgb8(approx.r * INV16)),
              std::max(std::abs(static_cast<int>(
                                    linear_float_to_srgb8(exact.g * INV16)) -
                                linear_float_to_srgb8(approx.g * INV16)),
                       std::abs(static_cast<int>(
                                    linear_float_to_srgb8(exact.b * INV16)) -
                                linear_float_to_srgb8(approx.b * INV16))));
          if (cyclic) {
            cyclic_delta_e = hs_test::fold_worst(cyclic_delta_e, delta_e);
            cyclic_srgb_delta = std::max(cyclic_srgb_delta, srgb_delta);
          } else {
            default_delta_e = hs_test::fold_worst(default_delta_e, delta_e);
            default_srgb_delta = std::max(default_srgb_delta, srgb_delta);
          }
        }
      }
    }
  }
  std::printf("  [hue table] default deltaE=%g sRGB8=%d; cyclic deltaE=%g "
              "sRGB8=%d\n",
              default_delta_e, default_srgb_delta, cyclic_delta_e,
              cyclic_srgb_delta);
  HS_EXPECT_LE(default_delta_e, 0.002f);
  HS_EXPECT_LE(default_srgb_delta, 10);
  HS_EXPECT_LE(cyclic_delta_e, 0.006f);
  HS_EXPECT_LE(cyclic_srgb_delta, 21);
}

struct DisplacementHueFrame {
  std::vector<Pixel> pixels;
  int table_uses;
};

/**
 * @brief Renders a deterministic DisplacementField frame with table control.
 * @param exact Whether to bypass the hue table.
 * @return Final RGB16 framebuffer and table-use count.
 */
inline DisplacementHueFrame render_displacement_hue_frame(bool exact) {
  constexpr int W = 192;
  constexpr int H = 40;
  constexpr int FRAMES = 64;
  reset_effect_globals();
  hs::set_mock_time(0, 0);
  DisplacementField<W, H> effect;
  effect.init();
  DisplacementFieldWhiteBox::set_force_exact_hue(effect, exact);
  for (int frame = 0; frame < FRAMES; ++frame) {
    hs::set_mock_time(frame * FRAME_MS, frame * FRAME_US);
    effect.draw_frame();
    effect.advance_display();
  }
  DisplacementHueFrame result{
      {}, DisplacementFieldWhiteBox::hue_table_uses(effect)};
  result.pixels.resize(W * H);
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x)
      result.pixels[y * W + x] = effect.get_pixel(x, y);
  hs::clear_mock_time();
  return result;
}

/**
 * @brief Bounds the table's final-frame delta against exact hue rotation.
 */
inline void test_displacement_field_hue_table_frame_fidelity() {
  DisplacementHueFrame exact = render_displacement_hue_frame(true);
  DisplacementHueFrame table = render_displacement_hue_frame(false);
  HS_EXPECT_EQ(exact.table_uses, 0);
  HS_EXPECT_GT(table.table_uses, 0);
  constexpr float INV16 = 1.0f / 65535.0f;
  float max_delta_e = 0.0f;
  int max_srgb_delta = 0;
  int changed = 0;
  size_t lit = 0;
  for (size_t i = 0; i < exact.pixels.size(); ++i) {
    lit += !is_black(exact.pixels[i]) && !is_black(table.pixels[i]);
    if (exact.pixels[i].r != table.pixels[i].r ||
        exact.pixels[i].g != table.pixels[i].g ||
        exact.pixels[i].b != table.pixels[i].b)
      ++changed;
    OKLab a = linear_rgb_to_oklab(exact.pixels[i].r * INV16,
                                  exact.pixels[i].g * INV16,
                                  exact.pixels[i].b * INV16);
    OKLab b = linear_rgb_to_oklab(table.pixels[i].r * INV16,
                                  table.pixels[i].g * INV16,
                                  table.pixels[i].b * INV16);
    const float dl = a.L - b.L;
    const float da = a.a - b.a;
    const float db = a.b - b.b;
    max_delta_e =
        hs_test::fold_worst(max_delta_e, sqrtf(dl * dl + da * da + db * db));
    max_srgb_delta = std::max(
        max_srgb_delta,
        std::abs(
            static_cast<int>(linear_float_to_srgb8(exact.pixels[i].r * INV16)) -
            linear_float_to_srgb8(table.pixels[i].r * INV16)));
    max_srgb_delta = std::max(
        max_srgb_delta,
        std::abs(
            static_cast<int>(linear_float_to_srgb8(exact.pixels[i].g * INV16)) -
            linear_float_to_srgb8(table.pixels[i].g * INV16)));
    max_srgb_delta = std::max(
        max_srgb_delta,
        std::abs(
            static_cast<int>(linear_float_to_srgb8(exact.pixels[i].b * INV16)) -
            linear_float_to_srgb8(table.pixels[i].b * INV16)));
  }
  std::printf("  [hue frame] changed=%d/%zu uses=%d deltaE=%g sRGB8=%d\n",
              changed, exact.pixels.size(), table.table_uses, max_delta_e,
              max_srgb_delta);
  HS_EXPECT_GT(lit, exact.pixels.size() / 100);
  HS_EXPECT_LE(max_delta_e, 0.01f);
  HS_EXPECT_LE(max_srgb_delta, 3);
}

/**
 * @brief Bounds DisplacementField's octave-grid noise bake against the exact
 *        noise field.
 * @details The bake evaluates each noise octave on every other knot and
 *          fills the rest with a Catmull-Rom spline; every knot stays within
 *          half a canvas column of the exact field.
 */
inline void test_displacement_field_octave_bake_tracks_noise() {
  constexpr int W = DEFAULT_W;
  constexpr int H = DEFAULT_H;
  reset_effect_globals();
  DisplacementField<W, H> effect;
  effect.init();
  Animation::NoiseProductParams np;
  np.noise.SetSeed(1234);
  np.amplitude = 0.2f;
  np.scale1 = 1.5f;
  np.scale2 = 3.0f;
  const int grid = DisplacementFieldWhiteBox::octave_grid(effect);
  std::vector<float> slut(W + 1);
  float worst = 0.0f;
  for (int trial = 0; trial < 12; ++trial) {
    np.time = 0.37f * trial;
    const math::Basis basis = math::make_basis(
        math::make_rotation(math::Vector(0.3f, 0.5f, 0.8f).normalized(),
                            0.9f * trial),
        math::X_AXIS);
    for (int ring = 0; ring < 24; ++ring) {
      const float theta = math::PI_F * (ring + 1) / 25.0f;
      const int lut_n = DisplacementFieldWhiteBox::noise_lut_samples(
          effect, np.scale1 + np.scale2, sinf(theta));
      DisplacementFieldWhiteBox::bake_noise_octaves(effect, np, basis, theta,
                                                    lut_n, slut.data());
      for (int x = 0; x < lut_n; ++x) {
        const float a = 2.0f * math::PI_F * x / lut_n;
        const math::Vector p =
            (basis.v * cosf(theta)) +
            ((basis.u * cosf(a)) + (basis.w * sinf(a))) * sinf(theta);
        const float exact = noise_product_field(p, np);
        const float err = std::fabs(slut[x] - exact);
        if (x % grid == 0) {
          HS_CONTEXT("trial / knot", trial, (ring << 16) | x);
          HS_EXPECT_NEAR(slut[x], exact, 1e-4f);
        }
        worst = hs_test::fold_worst(worst, err);
      }
    }
  }
  const float half_column = 0.5f * math::RADIANS_PER_COLUMN<W>;
  std::printf("  [octave bake] worst %g rad (half column %g)\n", worst,
              half_column);
  HS_EXPECT_LE(worst, half_column);
}

/**
 * @brief Verifies zero hue scale reproduces one exact zero rotation per ring.
 */
inline void test_displacement_field_zero_hue_scale_is_exact() {
  constexpr int W = 256;
  constexpr int H = 40;
  reset_effect_globals();
  hs::set_mock_time(0, 0);
  DisplacementField<W, H> effect;
  effect.init();
  DisplacementFieldWhiteBox::configure_noise(effect, 1.0f, 0.0f);
  effect.draw_frame();

  Color4 base = DisplacementFieldWhiteBox::current_ring_color(effect, 0.5f);
  Pixel expected = hue_rotate(make_hue_rotate_base(base), 0.0f).color;
  const int lut_n = DisplacementFieldWhiteBox::baked_lut_samples(effect);
  for (int i = 0; i <= lut_n; ++i) {
    Pixel actual = DisplacementFieldWhiteBox::hue_lut_value(effect, i);
    HS_EXPECT_EQ(actual.r, expected.r);
    HS_EXPECT_EQ(actual.g, expected.g);
    HS_EXPECT_EQ(actual.b, expected.b);
  }
  effect.advance_display();
  hs::clear_mock_time();
}

/**
 * @brief Verifies DisplacementField's clipped render tiles the full render:
 *        under a quadrant clip, every display-region pixel matches the
 *        full-canvas render within one 16-bit channel step, frame by frame.
 * @details Identical seeds and mock clock per run isolate the clip-only paths.
 */
inline void test_displacement_field_clip_tiles_full() {
  struct Quad {
    int x0, x1, y0, y1;
  };
  const Quad quads[] = {{0, DEFAULT_W / 2, 0, DEFAULT_H / 2},
                        {DEFAULT_W / 2, DEFAULT_W, DEFAULT_H / 2, DEFAULT_H}};
  const int frames = 60;
  // The shorter window still sits inside the same displacement phase.
  const int widest_frames = 24;

  size_t lit = 0;
  auto capture_region = [&](bool clip, const Quad &q, bool widest) {
    reset_effect_globals();
    hs::set_mock_time(0, 0);
    DisplacementField<DEFAULT_W, DEFAULT_H> fx;
    fx.init();
    if (widest)
      DisplacementFieldWhiteBox::configure_max_footprint(fx);
    DisplacementFieldWhiteBox::enter_balls(fx);
    if (clip)
      fx.set_clip(q.y0, q.y1, q.x0, q.x1);
    std::vector<Pixel> pixels;
    for (int f = 0; f < (widest ? widest_frames : frames); ++f) {
      hs::set_mock_time(static_cast<unsigned long>(f) * FRAME_MS,
                        static_cast<unsigned long>(f) * FRAME_US);
      fx.draw_frame();
      fx.advance_display();
      for (int y = q.y0; y < q.y1; ++y)
        for (int x = q.x0; x < q.x1; ++x) {
          const Pixel p = fx.get_pixel(x, y);
          if (!clip) {
            if (p.r | p.g | p.b)
              ++lit;
          }
          pixels.push_back(p);
        }
    }
    hs::clear_mock_time();
    return pixels;
  };

  for (bool widest : {false, true})
    for (const Quad &q : quads) {
      HS_CONTEXT("clip pair", widest, q.x0);
      lit = 0;
      const auto full = capture_region(false, q, widest);
      const auto clipped = capture_region(true, q, widest);
      HS_EXPECT_GT(lit, size_t{0});
      HS_EXPECT_EQ(full.size(), clipped.size());
      if (full.size() != clipped.size())
        continue;
      int max_error = 0;
      for (size_t i = 0; i < full.size(); ++i)
        max_error =
            std::max({max_error, std::abs(int(full[i].r) - clipped[i].r),
                      std::abs(int(full[i].g) - clipped[i].g),
                      std::abs(int(full[i].b) - clipped[i].b)});
      HS_EXPECT_LE(max_error, 1);
    }
}

inline void test_displacement_field_ball_spans_and_lifecycle() {
  reset_effect_globals();
  hs::set_mock_time(0, 0);
  DisplacementField<SMALL_W, SMALL_H> effect;
  effect.init();
  DisplacementFieldWhiteBox::fill_ball_pool(effect);
  for (int frame = 0; frame < 40; ++frame) {
    effect.draw_frame();
    effect.advance_display();
  }
  HS_EXPECT_TRUE(DisplacementFieldWhiteBox::in_noise(effect));

  const math::Basis basis = math::make_basis(math::Quaternion(), math::X_AXIS);
  for (float theta : {0.15f, 1.0f, 1.5f, 2.9f}) {
    std::array<Animation::BumpParams, 3> balls;
    for (size_t i = 0; i < balls.size(); ++i) {
      const float polar = theta + 0.08f;
      const float azimuth = i == 0 ? 0.0f : i == 1 ? -0.03f : 2.0f;
      auto &ball = balls[i];
      ball.center =
          basis.v * cosf(polar) +
          (basis.u * cosf(azimuth) + basis.w * sinf(azimuth)) * sinf(polar);
      ball.axis = basis.v;
      ball.radius = 0.3f;
      ball.amplitude = 0.8f;
      ball.envelope = 1.0f;
      ball.sync();
    }
    DisplacementFieldWhiteBox::check_ball_spans(effect, balls, basis, theta);
  }
  hs::clear_mock_time();
}
