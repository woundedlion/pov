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
  static bool history_ready(const DisplacementField<W, H> &effect) {
    return effect.draw_history.initialized;
  }
  /** @brief Holds the stack still with narrow, rapidly flowing noise.
   * @param effect Effect before its frame draw. */
  template <int W, int H>
  static void boundary_stencil_pose(DisplacementField<W, H> &effect) {
    effect.orientation.set(math::Quaternion());
    effect.params.noise_amp = 0.01f;
    effect.params.flow_speed = 0.15f;
    effect.params.thickness = 0.4f * math::RADIANS_PER_COLUMN<W>;
  }

  /** @brief Selects the single pass, window oracle, or per-column reference.
   * @param effect Initialized effect.
   * @param mode Zero for single pass, one for window bakes, two for columns. */
  template <int W, int H>
  static void sweep_mode(DisplacementField<W, H> &effect, int mode) {
    effect.force_window_bakes = mode != 0;
    effect.force_column_reference = mode == 2;
  }

  /** @brief Checks one setup per non-flat ring and reports the last frame's work.
   * @param effect Effect after a default swept frame. */
  template <int W, int H>
  static void check_single_pass(const DisplacementField<W, H> &effect) {
    HS_EXPECT_EQ(effect.stack_flushes, 1);
    HS_EXPECT_LE(effect.primary_bakes, effect.params.num_rings);
    HS_EXPECT_GT(effect.duplicate_knots, 0);
    std::printf(
        "  [single pass] primary=%d alternate=%d octave-extra=%d rings=%d\n",
        effect.primary_knots, effect.duplicate_knots,
        effect.alternate_noise_samples, effect.primary_bakes);
  }

  /** @brief Turns the stack so its rings circle the arm-boundary plane's
   *  normal, rocking it slightly each frame, with every ring drawn.
   * @param effect Effect before its frame draw.
   * @param frame Frame number. */
  template <int W, int H>
  static void arm_plane_pose(DisplacementField<W, H> &effect, int frame) {
    effect.orientation.set(math::make_rotation(
        math::Y_AXIS, -0.5f * math::PI_F + 0.004f * static_cast<float>(frame)));
    effect.params.num_rings = DisplacementField<W, H>::RING_SLOTS;
    effect.params.noise_amp = 0.6f;
  }

  /** @brief Stack scans of the last swept frame.
   * @param effect Effect after a swept frame.
   * @return The frame's stack scans. */
  template <int W, int H>
  static int stack_flushes(const DisplacementField<W, H> &effect) {
    return effect.stack_flushes;
  }

  /** @brief The frame turn a swept draw derives from two endpoints.
   * @param effect Effect whose noise field is live.
   * @param previous Stack orientation at the previous endpoint.
   * @param current Stack orientation at the current endpoint.
   * @return The turn carrying the previous endpoint's points to the current
   *         one's. */
  template <int W, int H>
  static math::Vector sweep_turn(DisplacementField<W, H> &effect,
                                 const math::Quaternion &previous,
                                 const math::Quaternion &current) {
    effect.draw_history.capture({.pose = previous});
    typename DisplacementField<W, H>::Sweep sweep;
    HS_EXPECT_TRUE(effect.prepare_sweep(previous, current, 1.0f, 0, sweep));
    return sweep.turn;
  }

  /** @brief Denies preset advancement at one phase departure boundary.
   * @param effect Initialized effect with no live balls.
   * @param leaving_noise Selects the Noise departure rather than Balls. */
  template <int W, int H>
  static void deny_phase_advance(DisplacementField<W, H> &effect,
                                 bool leaving_noise) {
    effect.preset_count = 0;
    effect.preset_index =
        leaving_noise ? effect.NOISE_PRESET : effect.BALLS_PRESET;
    effect.noise_hold = 0;
    effect.master_gain = 0.0f;
    effect.ball_phase_left = 0;
  }

  /** @brief Forces the slowest default fall and a spawn on the window's last frame.
   * @param effect Initialized effect. */
  template <int W, int H>
  static void force_last_slow_ball(DisplacementField<W, H> &effect) {
    effect.params.ball_speed_max = effect.params.ball_speed_min;
    if (effect.ball_phase_left == 1)
      effect.spawn_cooldown = 0.0f;
  }

  /** @brief Requests preview progress without an active frame canvas.
   * @param effect Initialized effect with a free ball slot. */
  template <int W, int H>
  static void spawn_preview_without_canvas(DisplacementField<W, H> &effect) {
    effect.spawn_ball(0.5f);
  }

  template <int W, int H>
  static void start_noise_fade(DisplacementField<W, H> &effect, bool fade_out) {
    if (fade_out) {
      effect.master_gain = 1.0f;
      effect.noise_hold = 1;
    } else {
      effect.enter_noise();
    }
  }

  template <int W, int H>
  static float gain(const DisplacementField<W, H> &effect) {
    return effect.master_gain;
  }

  template <int W, int H>
  static float noise_amplitude(const DisplacementField<W, H> &effect) {
    return effect.noise_field.template_params.amplitude;
  }

  template <int W, int H>
  static int hold(const DisplacementField<W, H> &effect) {
    return effect.noise_hold;
  }

  template <int W, int H> static constexpr int fade_frames() {
    return DisplacementField<W, H>::NOISE_FADE_FRAMES;
  }

  template <int W, int H>
  static void enter_balls(DisplacementField<W, H> &effect) {
    effect.preset_index = effect.BALLS_PRESET;
    effect.enter_balls();
    effect.spawn_ball();
  }

  template <int W, int H>
  static void fill_ball_pool(DisplacementField<W, H> &effect) {
    effect.params.ball_speed_min = effect.params.ball_speed_max = 3.0f;
    effect.preset_index = effect.BALLS_PRESET;
    effect.enter_balls();
    for (int i = 0; i <= effect.MAX_BALLS; ++i)
      effect.spawn_ball();
    HS_EXPECT_EQ(effect.balls.active_count(), effect.MAX_BALLS);
    HS_EXPECT_TRUE(effect.logged_pool_full);
    effect.ball_phase_left = 0;
  }

  template <int W, int H>
  static int active_balls(const DisplacementField<W, H> &effect) {
    return effect.balls.active_count();
  }

  template <int W, int H>
  static std::vector<math::Vector> knot_positions(float theta, int knots) {
    const math::Basis basis =
        math::make_basis(math::Quaternion(), math::X_AXIS);
    const float delta = 2.0f * math::PI_F / knots;
    typename DisplacementField<W, H>::RingKnotWalk walk{
        basis, cosf(theta), sinf(theta), cosf(delta), sinf(delta)};
    std::vector<math::Vector> points;
    for (int knot = 0; knot < knots; ++knot) {
      points.push_back(walk.position());
      walk.advance();
    }
    return points;
  }

  template <int W, int H>
  static float first_ball_height(const DisplacementField<W, H> &effect) {
    return effect.balls.active_params(0).center.y;
  }

  template <int W, int H>
  static int ball_window(const DisplacementField<W, H> &effect) {
    return effect.ball_phase_left;
  }

  template <int W, int H>
  static float ball_bound(const DisplacementField<W, H> &effect) {
    return effect.balls.field_bound();
  }

  template <int W, int H>
  static float ball_amplitude(const DisplacementField<W, H> &effect) {
    return effect.balls.active_params(0).amplitude;
  }

  template <int W, int H>
  static void expire_balls(DisplacementField<W, H> &effect) {
    effect.ball_phase_left = 0;
  }

  template <int W, int H>
  static float color_spin(const DisplacementField<W, H> &effect) {
    return effect.color_spin;
  }

  template <int W, int H>
  static bool in_noise(const DisplacementField<W, H> &effect) {
    return effect.getPresetIndex() == effect.NOISE_PRESET;
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
    effect.bake_ball_spans(basis,
                           {cosf(theta), sinf(theta), cosf(step), sinf(step)},
                           theta, W, effect.CHUNK_MASK, N, shifts.data());
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
  static void bind_hue_table(DisplacementField<W, H> &effect,
                             const HueRotateBase &base, float domain,
                             bool cyclic) {
    effect.hue_table.bind(base, domain, cyclic);
  }

  template <int W, int H>
  static void prepare_hue_table(DisplacementField<W, H> &effect) {
    effect.hue_table.prepare();
  }

  template <int W, int H>
  static Pixel sample_hue_table(const DisplacementField<W, H> &effect,
                                float amount) {
    return effect.hue_table.sample(amount);
  }

  template <int W, int H>
  static Pixel sample_hue_table_lazy(DisplacementField<W, H> &effect,
                                     float amount) {
    return effect.hue_table.sample_lazy(amount);
  }

  template <int W, int H>
  static bool hue_knot_valid(const DisplacementField<W, H> &effect, int index) {
    return (effect.hue_table.valid[index >> 6] >> (index & 63)) & 1u;
  }

  template <int W, int H>
  static int hue_table_size(const DisplacementField<W, H> &) {
    return DisplacementField<W, H>::HUE_TABLE_SIZE;
  }

  template <int W, int H>
  static Pixel hue_table_value(const DisplacementField<W, H> &effect,
                               int index) {
    return effect.hue_table.knots[index];
  }

  template <int W, int H>
  static void clear_hue_table(DisplacementField<W, H> &effect) {
    for (int i = 0; i <= DisplacementField<W, H>::HUE_TABLE_SIZE; ++i)
      effect.hue_table.knots[i] = Pixel(0, 0, 0);
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
    effect.params.ball_min = effect.params.ball_max = 1.0f;
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
    return effect.rings.hue_row(0)[index];
  }

  template <int W, int H>
  static auto &ring_pool(DisplacementField<W, H> &effect) {
    return effect.rings;
  }

  template <int W, int H>
  static int ring_slots(const DisplacementField<W, H> &) {
    return DisplacementField<W, H>::RingPool::SLOTS;
  }

  template <int W, int H>
  static size_t footprint_bytes(const DisplacementField<W, H> &) {
    return DisplacementField<W, H>::FOOTPRINT_BYTES;
  }

  template <int W, int H>
  static Color4 current_ring_color(const DisplacementField<W, H> &effect,
                                   float ring_t) {
    return effect.palette.get(math::wrap_t(ring_t + effect.color_spin));
  }

  template <int W, int H>
  static int baked_lut_samples(const DisplacementField<W, H> &effect) {
    return static_cast<int>(effect.rings.lut_columns(0));
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
    effect.bake_noise_octaves(
        np, basis, {cosf(theta), sinf(theta), cosf(dphi), sinf(dphi)}, theta,
        lut_n, DisplacementField<W, H>::CHUNK_MASK, 0, slut);
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
  const HueRotateBase base = make_hue_rotate_base(palette.get(0.37f));
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
    DisplacementFieldWhiteBox::bind_hue_table(effect, base, table_case.domain,
                                              table_case.cyclic);
    DisplacementFieldWhiteBox::prepare_hue_table(effect);
    std::vector<Pixel> endpoints(table_size + 1);
    for (int i = 0; i <= table_size; ++i)
      endpoints[i] = DisplacementFieldWhiteBox::hue_table_value(effect, i);
    std::vector<Pixel> expected(SAMPLE_COUNT + 1);
    for (int i = 0; i <= SAMPLE_COUNT; ++i) {
      const float amount = table_case.max_amount * i / SAMPLE_COUNT;
      expected[i] = DisplacementFieldWhiteBox::sample_hue_table(effect, amount);
    }

    DisplacementFieldWhiteBox::clear_hue_table(effect);
    DisplacementFieldWhiteBox::bind_hue_table(effect, base, table_case.domain,
                                              table_case.cyclic);
    for (int i = 0; i <= SAMPLE_COUNT; ++i) {
      const float amount = table_case.max_amount * i / SAMPLE_COUNT;
      Pixel actual =
          DisplacementFieldWhiteBox::sample_hue_table_lazy(effect, amount);
      HS_EXPECT_EQ(actual, expected[i]);
    }
    for (int i = 0; i <= table_size; ++i)
      HS_EXPECT_TRUE(DisplacementFieldWhiteBox::hue_knot_valid(effect, i));
    for (int i = 0; i <= table_size; ++i) {
      Pixel actual = DisplacementFieldWhiteBox::hue_table_value(effect, i);
      HS_EXPECT_EQ(actual, endpoints[i]);
    }
  }
}

/**
 * @brief Verifies a bind() for a new ring rebakes a knot the previous ring's
 *        lazy samples already baked.
 */
inline void test_displacement_field_hue_table_bind_clears_knots() {
  reset_effect_globals();
  DisplacementField<SMALL_W, SMALL_H> effect;
  effect.init();
  GenerativePalette palette(PaletteRecipes::profile(PaletteDomain::MIRROR,
                                                    PaletteHarmony::ANALOGOUS,
                                                    AxisCurve::CONSTANT, 0.0f));
  const HueRotateBase first = make_hue_rotate_base(palette.get(0.1f));
  const HueRotateBase second = make_hue_rotate_base(palette.get(0.6f));
  const float domain = 0.4f;
  const float amount = 0.0f;

  DisplacementFieldWhiteBox::bind_hue_table(effect, first, domain, false);
  const Pixel first_sample =
      DisplacementFieldWhiteBox::sample_hue_table_lazy(effect, amount);
  DisplacementFieldWhiteBox::bind_hue_table(effect, second, domain, false);
  HS_EXPECT_FALSE(DisplacementFieldWhiteBox::hue_knot_valid(effect, 0));
  const Pixel second_sample =
      DisplacementFieldWhiteBox::sample_hue_table_lazy(effect, amount);

  const Pixel expected = hue_rotate(second, amount).color;
  HS_EXPECT_EQ(second_sample, expected);
  HS_EXPECT_TRUE(first_sample != second_sample);
}

/** @brief Reciprocal of the RGB16 channel full scale. */
constexpr float HUE_INV16 = 1.0f / 65535.0f;

/** @brief OKLab Euclidean distance between two linear RGB16 pixels. */
inline float oklab_delta_e(const Pixel &a, const Pixel &b) {
  const OKLab lab_a =
      linear_rgb_to_oklab(a.r * HUE_INV16, a.g * HUE_INV16, a.b * HUE_INV16);
  const OKLab lab_b =
      linear_rgb_to_oklab(b.r * HUE_INV16, b.g * HUE_INV16, b.b * HUE_INV16);
  const float dl = lab_a.L - lab_b.L;
  const float da = lab_a.a - lab_b.a;
  const float db = lab_a.b - lab_b.b;
  return sqrtf(dl * dl + da * da + db * db);
}

/** @brief Largest per-channel gap between two linear RGB16 pixels in sRGB8. */
inline int srgb8_channel_gap(const Pixel &a, const Pixel &b) {
  const auto gap = [](uint16_t x, uint16_t y) {
    return std::abs(static_cast<int>(linear_float_to_srgb8(x * HUE_INV16)) -
                    linear_float_to_srgb8(y * HUE_INV16));
  };
  return std::max({gap(a.r, b.r), gap(a.g, b.g), gap(a.b, b.b)});
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
        DisplacementFieldWhiteBox::bind_hue_table(effect, exact_base, domain,
                                                  cyclic);
        DisplacementFieldWhiteBox::prepare_hue_table(effect);
        for (int i = 0; i <= 1024; ++i) {
          const float amount = max_amount * i / 1024.0f;
          Pixel exact = hue_rotate(exact_base, amount).color;
          Pixel approx =
              DisplacementFieldWhiteBox::sample_hue_table(effect, amount);
          const float delta_e = oklab_delta_e(exact, approx);
          const int srgb_delta = srgb8_channel_gap(exact, approx);
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
  pin_frame_clock(0);
  DisplacementField<W, H> effect;
  effect.init();
  DisplacementFieldWhiteBox::set_force_exact_hue(effect, exact);
  for (int frame = 0; frame < FRAMES; ++frame) {
    pin_frame_clock(frame);
    effect.draw_frame();
    effect.advance_display();
  }
  DisplacementHueFrame result{
      {}, DisplacementFieldWhiteBox::hue_table_uses(effect)};
  capture_frame<W, H>(effect, result.pixels);
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
  float max_delta_e = 0.0f;
  int max_srgb_delta = 0;
  int changed = 0;
  size_t lit = 0;
  for (size_t i = 0; i < exact.pixels.size(); ++i) {
    lit += !is_black(exact.pixels[i]) && !is_black(table.pixels[i]);
    if (exact.pixels[i] != table.pixels[i])
      ++changed;
    max_delta_e = hs_test::fold_worst(
        max_delta_e, oklab_delta_e(exact.pixels[i], table.pixels[i]));
    max_srgb_delta = std::max(
        max_srgb_delta, srgb8_channel_gap(exact.pixels[i], table.pixels[i]));
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
  pin_frame_clock(0);
  DisplacementField<W, H> effect;
  effect.init();
  DisplacementFieldWhiteBox::configure_noise(effect, 1.0f, 0.0f);
  effect.draw_frame();

  Color4 base = DisplacementFieldWhiteBox::current_ring_color(effect, 0.5f);
  Pixel expected = hue_rotate(make_hue_rotate_base(base), 0.0f).color;
  const int lut_n = DisplacementFieldWhiteBox::baked_lut_samples(effect);
  for (int i = 0; i <= lut_n; ++i) {
    Pixel actual = DisplacementFieldWhiteBox::hue_lut_value(effect, i);
    HS_EXPECT_EQ(actual, expected);
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
    pin_frame_clock(0);
    DisplacementField<DEFAULT_W, DEFAULT_H> fx;
    fx.init();
    if (widest)
      DisplacementFieldWhiteBox::configure_max_footprint(fx);
    DisplacementFieldWhiteBox::enter_balls(fx);
    if (clip)
      fx.set_clip(q.y0, q.y1, q.x0, q.x1);
    std::vector<Pixel> pixels;
    for (int f = 0; f < (widest ? widest_frames : frames); ++f) {
      pin_frame_clock(f);
      fx.draw_frame();
      fx.advance_display();
      for (int y = q.y0; y < q.y1; ++y)
        for (int x = q.x0; x < q.x1; ++x) {
          const Pixel p = fx.get_pixel(x, y);
          if (!clip) {
            if (!is_black(p))
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
            std::max<int>(max_error, max_channel_gap(full[i], clipped[i]));
      HS_EXPECT_LE(max_error, 1);
    }
}

inline void test_displacement_field_ball_spans_and_lifecycle() {
  reset_effect_globals();
  pin_frame_clock(0);
  DisplacementField<SMALL_W, SMALL_H> effect;
  effect.init();
  DisplacementFieldWhiteBox::fill_ball_pool(effect);
  hs_test::run_frames(effect, 40);
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

/**
 * @brief Verifies next() is idempotent until a commit, and a ring culled after
 *        next() leaves the pool unchanged.
 */
inline void test_displacement_field_ring_pool_abandoned_slot() {
  using WB = DisplacementFieldWhiteBox;
  reset_effect_globals();
  DisplacementField<SMALL_W, SMALL_H> effect;
  effect.init();
  auto &pool = WB::ring_pool(effect);
  using Pool = std::remove_reference_t<decltype(pool)>;
  const math::Basis basis = math::make_basis(math::Quaternion(), math::X_AXIS);
  constexpr int LUT_N = 16;

  pool.begin_frame();
  const auto first = pool.next();
  const auto again = pool.next();
  HS_EXPECT_TRUE(first.shift_row == again.shift_row);
  HS_EXPECT_TRUE(first.hue_row == again.hue_row);
  HS_EXPECT_EQ(pool.size(), 0);

  for (int x = 0; x <= LUT_N; ++x)
    first.shift_row[x] = 0.0f;
  pool.commit(2, 0.5f, LUT_N, basis, 0.5f, 0.03f, first.shift_row, LUT_N, 0.0f,
              nullptr);
  HS_EXPECT_EQ(pool.size(), 1);
  HS_EXPECT_EQ(pool.slot_of(2), 0);
  HS_EXPECT_EQ(pool.slot_of(0), Pool::CULLED);
  HS_EXPECT_EQ(pool.slot_of(1), Pool::CULLED);
  HS_EXPECT_EQ(pool.lut_columns(0), static_cast<float>(LUT_N));
  HS_EXPECT_EQ(pool.frag_alpha(0), 0.5f);
  HS_EXPECT_TRUE(pool.hue_row(0) == first.hue_row);
  const auto second = pool.next();
  HS_EXPECT_TRUE(second.shift_row != first.shift_row);
  HS_EXPECT_TRUE(second.hue_row != first.hue_row);
  pool.release();
  HS_EXPECT_EQ(pool.size(), 0);
}

/**
 * @brief Verifies the pool commits up to capacity and release() destroys
 *        exactly the committed rings.
 */
inline void test_displacement_field_ring_pool_capacity_release() {
  using WB = DisplacementFieldWhiteBox;
  reset_effect_globals();
  DisplacementField<SMALL_W, SMALL_H> effect;
  effect.init();
  auto &pool = WB::ring_pool(effect);
  const int slots = WB::ring_slots(effect);
  const math::Basis basis = math::make_basis(math::Quaternion(), math::X_AXIS);

  for (int frame = 0; frame < 2; ++frame) {
    pool.begin_frame();
    const int committed = frame == 0 ? slots : 3;
    for (int i = 0; i < committed; ++i)
      pool.commit(i, 1.0f, 16, basis, 0.5f, 0.03f,
                  ScalarFn([](float) { return 0.0f; }), 0.0f, 0.0f);
    HS_EXPECT_EQ(pool.size(), committed);
    for (int i = 0; i < committed; ++i)
      HS_EXPECT_EQ(pool.slot_of(i), i);
    const int destroyed_before = pool.destroyed;
    pool.release();
    HS_EXPECT_EQ(pool.destroyed - destroyed_before, committed);
    HS_EXPECT_EQ(pool.size(), 0);
  }
}

/** @brief Verifies FOOTPRINT_BYTES covers every persistent allocation. */
inline void test_displacement_field_footprint_covers_init() {
  reset_effect_globals();
  const size_t before = persistent_arena.get_offset();
  DisplacementField<DEFAULT_W, DEFAULT_H> effect;
  effect.init();
  const size_t used = persistent_arena.get_offset() - before;
  HS_EXPECT_GT(used, size_t{0});
  HS_EXPECT_LE(used, DisplacementFieldWhiteBox::footprint_bytes(effect));
}

/** @brief Pause holds noise departure while its entry fade remains visible. */
inline void test_displacement_field_pause_holds_phase_only() {
  using WB = DisplacementFieldWhiteBox;
  for (bool fade_out : {false, true}) {
    reset_effect_globals();
    DisplacementField<SMALL_W, SMALL_H> effect;
    effect.init();
    if (fade_out)
      hs_test::run_frames(effect, WB::fade_frames<SMALL_W, SMALL_H>() + 1);
    WB::start_noise_fade(effect, fade_out);
    if (fade_out)
      hs_test::run_frames(effect, WB::fade_frames<SMALL_W, SMALL_H>() / 2);
    effect.setAnimationsPaused(true);
    const int HOLD = WB::hold(effect);
    const float GAIN = WB::gain(effect);
    if (fade_out) {
      HS_EXPECT_GT(GAIN, 0.0f);
      HS_EXPECT_LT(GAIN, 1.0f);
    }
    for (int frame = 0; frame < WB::fade_frames<SMALL_W, SMALL_H>() + 2;
         ++frame) {
      effect.draw_frame();
      effect.advance_display();
    }
    HS_EXPECT_TRUE(WB::in_noise(effect));
    HS_EXPECT_EQ(WB::hold(effect), HOLD);
    HS_EXPECT_NEAR(WB::gain(effect), fade_out ? GAIN : 1.0f, 1e-6f);

    if (fade_out) {
      effect.setAnimationsPaused(false);
      hs_test::run_frames(effect, WB::fade_frames<SMALL_W, SMALL_H>() + 2);
      HS_EXPECT_FALSE(WB::in_noise(effect));
    }
  }
}

/** @brief Both phases support authored defaults, navigation, pause and tuning. */
inline void test_displacement_field_presets() {
  using FX = DisplacementField<SMALL_W, SMALL_H>;
  using WB = DisplacementFieldWhiteBox;
  reset_effect_globals();
  FX effect;
  effect.init();
  HS_EXPECT_EQ(effect.getPresetCount(), size_t{2});
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{0});
  HS_EXPECT_EQ(FX::PRESET_IDS[0], "noise");
  HS_EXPECT_EQ(FX::PRESET_IDS[1], "balls");
  HS_EXPECT_FALSE(effect.selectPreset(2));
  HS_EXPECT_FALSE(effect.animations_paused());
  HS_EXPECT_TRUE(effect.previousPreset());
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});
  HS_EXPECT_TRUE(effect.animations_paused());
  hs_test::run_frames(effect, 2);
  HS_EXPECT_GE(WB::active_balls(effect), 8);
  HS_EXPECT_GT(WB::ball_bound(effect), 0.0f);
  HS_EXPECT_EQ(effect.updateParameter("Ball Amp", 0.6f),
               ParamSetResult::APPLIED);
  hs_test::run_frames(effect, 2);
  HS_EXPECT_NEAR(WB::ball_amplitude(effect), 2.4f, 1e-6f);
  const float SPIN = WB::color_spin(effect);
  hs_test::run_frames(effect, 10);
  HS_EXPECT_TRUE(WB::color_spin(effect) != SPIN);
  HS_EXPECT_GE(WB::active_balls(effect), 8);
  HS_EXPECT_TRUE(effect.selectPreset(1));
  hs_test::run_frames(effect, 2);
  HS_EXPECT_GE(WB::active_balls(effect), 8);
  HS_EXPECT_NEAR(effect.getParameters().find("Ball Amp")->get(),
                 FX::PRESETS[1].ball_amp, 1e-6f);
  for (int repeat = 0; repeat < 100; ++repeat) {
    HS_EXPECT_TRUE(effect.selectPreset(repeat % 2));
  }
  hs_test::run_frames(effect, 2);
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});
  HS_EXPECT_GE(WB::active_balls(effect), 8);
  HS_EXPECT_TRUE(effect.nextPreset());
  hs_test::run_frames(effect, 2);
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{0});
  HS_EXPECT_EQ(WB::active_balls(effect), 0);
  HS_EXPECT_EQ(WB::gain(effect), 1.0f);
  HS_EXPECT_TRUE(effect.synchronizePreset(1));
  HS_EXPECT_TRUE(effect.animations_paused());
  hs_test::run_frames(effect, 2);
  HS_EXPECT_GE(WB::active_balls(effect), 8);
  effect.setAnimationsPaused(false);
  HS_EXPECT_TRUE(effect.synchronizePreset(0));
  HS_EXPECT_FALSE(effect.animations_paused());
  hs_test::run_frames(effect, 800);
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});
  WB::expire_balls(effect);
  hs_test::run_frames(effect, 200);
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{0});
  HS_EXPECT_TRUE(WB::in_noise(effect));
}

/** @brief Editing a departing noise preset restores its tuned field and dwell. */
inline void test_displacement_field_edit_cancels_departure() {
  using WB = DisplacementFieldWhiteBox;
  reset_effect_globals();
  DisplacementField<SMALL_W, SMALL_H> effect;
  effect.init();
  hs_test::run_frames(effect, WB::fade_frames<SMALL_W, SMALL_H>() + 1);
  WB::start_noise_fade(effect, true);
  hs_test::run_frames(effect, WB::fade_frames<SMALL_W, SMALL_H>() / 2);
  HS_EXPECT_LT(WB::gain(effect), 1.0f);
  effect.setAnimationsPaused(true);
  HS_EXPECT_EQ(effect.updateParameter("Noise Amp", 0.6f),
               ParamSetResult::APPLIED);
  hs_test::run_frames(effect, 200);
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{0});
  HS_EXPECT_EQ(WB::gain(effect), 1.0f);
  HS_EXPECT_NEAR(WB::noise_amplitude(effect), 0.6f, 1e-6f);
  effect.setAnimationsPaused(false);
  hs_test::run_frames(effect, 10);
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{0});
}

/** @brief A paused Balls preset renders its tuned deformation. */
inline void test_displacement_field_paused_balls_deform_render() {
  const auto render = [](float amplitude) {
    reset_effect_globals();
    pin_frame_clock(0);
    DisplacementField<DEFAULT_W, DEFAULT_H> effect;
    effect.init();
    HS_EXPECT_TRUE(effect.selectPreset(1));
    HS_EXPECT_EQ(effect.updateParameter("Ball Amp", amplitude),
                 ParamSetResult::APPLIED);
    for (int frame = 0; frame < 32; ++frame) {
      pin_frame_clock(frame);
      effect.draw_frame();
      effect.advance_display();
    }
    std::vector<Pixel> pixels;
    capture_frame<DEFAULT_W, DEFAULT_H>(effect, pixels);
    hs::clear_mock_time();
    return pixels;
  };
  const auto plain = render(0.0f);
  const auto deformed = render(0.6f);
  HS_EXPECT_TRUE(plain != deformed);
}

/** @brief Preset pause holds the Balls window while falls and spawning continue. */
inline void test_displacement_field_paused_balls_keep_moving() {
  using WB = DisplacementFieldWhiteBox;
  reset_effect_globals();
  DisplacementField<SMALL_W, SMALL_H> effect;
  effect.init();
  HS_EXPECT_TRUE(effect.selectPreset(1));
  hs_test::run_frames(effect, 2);
  const float HEIGHT = WB::first_ball_height(effect);
  const int COUNT = WB::active_balls(effect);
  const int WINDOW = WB::ball_window(effect);
  hs_test::run_frames(effect, 10);
  HS_EXPECT_LT(WB::first_ball_height(effect), HEIGHT);
  HS_EXPECT_GT(WB::active_balls(effect), COUNT);
  HS_EXPECT_EQ(WB::ball_window(effect), WINDOW);
  HS_EXPECT_EQ(effect.updateParameter("Ball Amp", 0.6f),
               ParamSetResult::APPLIED);
  hs_test::run_frames(effect, 1200);
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});
  HS_EXPECT_EQ(WB::ball_window(effect), WINDOW);
  HS_EXPECT_GT(WB::active_balls(effect), 0);
  HS_EXPECT_NEAR(WB::ball_amplitude(effect), 2.4f, 1e-6f);
  WB::expire_balls(effect);
  hs_test::run_frames(effect, 200);
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});
  HS_EXPECT_GT(WB::active_balls(effect), 0);
  effect.setAnimationsPaused(false);
  hs_test::run_frames(effect, 200);
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{0});
}

/** @brief Knot recurrence preserves unit samples beside near-coincident balls. */
inline void test_displacement_field_knot_walk_unit_near_ball() {
  constexpr float THETA = 3.0f;
  constexpr float OFFSET = 0.0015f;
  const math::Basis basis = math::make_basis(math::Quaternion(), math::X_AXIS);
  for (int knots : {32, 144, 288}) {
    const auto points =
        DisplacementFieldWhiteBox::knot_positions<288, 144>(THETA, knots);
    for (const auto &point : points) {
      HS_EXPECT_NEAR(math::dot(point, point), 1.0f, 3e-7f);
      const math::Vector tangent = (point - basis.v * cosf(THETA)).normalized();
      Animation::BumpParams ball;
      ball.axis = basis.v;
      ball.center =
          (basis.v * cosf(THETA + OFFSET) + tangent * sinf(THETA + OFFSET))
              .normalized();
      ball.radius = 0.1f;
      ball.amplitude = 2.4f;
      ball.envelope = 1.0f;
      ball.sync();
      const float shift = bump_field_with_y(point, ball, -OFFSET);
      HS_EXPECT_TRUE(std::isfinite(shift));
      HS_EXPECT_LE(std::fabs(shift), ball.field_bound());
    }
  }
}

/** @brief Startup parameter writes preserve the paused preset's noise entry ramp. */
inline void test_displacement_field_startup_writes_preserve_entry() {
  using WB = DisplacementFieldWhiteBox;
  for (float amplitude : {0.2f, 0.6f}) {
    reset_effect_globals();
    DisplacementField<SMALL_W, SMALL_H> effect;
    effect.init();
    HS_EXPECT_EQ(WB::gain(effect), 0.0f);
    HS_EXPECT_EQ(effect.updateParameter("Noise Amp", amplitude),
                 ParamSetResult::APPLIED);
    effect.setAnimationsPaused(true);
    HS_EXPECT_EQ(WB::gain(effect), 0.0f);
    HS_EXPECT_EQ(effect.getPresetIndex(), size_t{0});
    hs_test::run_frames(effect, WB::fade_frames<SMALL_W, SMALL_H>() / 2);
    HS_EXPECT_GT(WB::gain(effect), 0.0f);
    HS_EXPECT_LT(WB::gain(effect), 1.0f);
    HS_EXPECT_EQ(effect.getPresetIndex(), size_t{0});
    hs_test::run_frames(effect, WB::fade_frames<SMALL_W, SMALL_H>() + 2);
    HS_EXPECT_EQ(WB::gain(effect), 1.0f);
    HS_EXPECT_NEAR(WB::noise_amplitude(effect), amplitude, 1e-6f);
    HS_EXPECT_EQ(effect.getPresetIndex(), size_t{0});
  }
}

/** @brief The default cycle budget covers the last ball draining before Noise returns. */
inline void test_displacement_field_default_cycle_frames() {
  reset_effect_globals();
  using DF = DisplacementField<SMALL_W, SMALL_H>;
  DF effect;
  effect.init();
  bool saw_balls = false;
  bool saw_drain = false;
  bool returned = false;
  for (size_t frame = 0; frame < DF::cycle_frames(); ++frame) {
    DisplacementFieldWhiteBox::force_last_slow_ball(effect);
    effect.draw_frame();
    effect.advance_display();
    if (!DisplacementFieldWhiteBox::in_noise(effect)) {
      saw_balls = true;
      if (DisplacementFieldWhiteBox::ball_window(effect) == 0) {
        saw_drain |= DisplacementFieldWhiteBox::active_balls(effect) > 0;
      }
    } else if (saw_balls) {
      HS_EXPECT_EQ(DisplacementFieldWhiteBox::active_balls(effect), 0);
      returned = true;
      break;
    }
  }
  HS_EXPECT_TRUE(saw_balls);
  HS_EXPECT_TRUE(saw_drain);
  HS_EXPECT_TRUE(returned);
}

/** @brief A frame's displayed endpoints and its column-time draw. */
struct DisplacementFieldEndpoints {
  std::vector<Pixel> previous; ///< Previous frame without temporal sampling.
  std::vector<Pixel> current;  ///< This frame without temporal sampling.
  std::vector<Pixel> swept;    ///< This frame with temporal sampling.
};

/**
 * @brief Draws one frame of DisplacementField with and without temporal
 *        sampling, from the same start.
 * @tparam W Render width.
 * @tparam H Render height.
 * @param frames Frames drawn; the last is the one compared.
 * @param setup Called on the initialized effect before the first frame.
 * @return The frame's endpoints and its column-time draw.
 */
template <int W, int H, typename Setup>
inline DisplacementFieldEndpoints displacement_field_endpoints(int frames,
                                                               Setup &&setup) {
  DisplacementFieldEndpoints out;
  for (bool enabled : {false, true}) {
    reset_effect_globals();
    hs::random().seed(1337u);
    pin_frame_clock(0);
    Spacetime::set_flat_surface(!enabled);
    DisplacementField<W, H> effect;
    effect.init();
    setup(effect);
    for (int f = 0; f < frames; ++f) {
      pin_frame_clock(f);
      effect.draw_frame();
      effect.advance_display();
      if (!enabled && f == frames - 2)
        capture_frame<W, H>(effect, out.previous);
      if (f == frames - 1)
        capture_frame<W, H>(effect, enabled ? out.swept : out.current);
    }
  }
  Spacetime::set_flat_surface(false);
  hs::clear_mock_time();
  return out;
}

/**
 * @brief Largest channel difference of two frames over some columns.
 * @tparam W Frame width.
 * @tparam H Frame height.
 * @param a First frame.
 * @param b Second frame.
 * @param columns Columns compared.
 * @return The largest difference, in 16-bit channel units.
 */
template <int W, int H>
inline int column_difference(const std::vector<Pixel> &a,
                             const std::vector<Pixel> &b,
                             std::initializer_list<int> columns) {
  int worst = 0;
  for (int x : columns)
    for (int y = 0; y < H; ++y) {
      const Pixel &p = a[static_cast<size_t>(y) * W + x];
      const Pixel &q = b[static_cast<size_t>(y) * W + x];
      for (int e : {std::abs(int(p.r) - q.r), std::abs(int(p.g) - q.g),
                    std::abs(int(p.b) - q.b)})
        worst = std::max(worst, e);
    }
  return worst;
}

/** @brief Largest difference, in 8-bit steps, between a column-time draw's
 *  arm-boundary columns and the frame endpoint they are lit at. */
inline constexpr int DF_ENDPOINT_STEPS = 2;
/** @brief Smallest endpoint jump a case must show at its first columns, in
 *  8-bit steps, so a held field fails it. */
inline constexpr int DF_ENDPOINT_JUMP_STEPS = 8;

/**
 * @brief Checks one case's arm-boundary columns against its endpoints.
 * @param label Case name.
 * @param frames The case's endpoints and column-time draw.
 */
template <int W, int H>
inline void
check_displacement_field_endpoints(const char *label,
                                   const DisplacementFieldEndpoints &frames) {
  HS_CONTEXT(label);
  constexpr int STEP = 257;
  const int jump =
      column_difference<W, H>(frames.previous, frames.current, {0, W / 2});
  const int head =
      column_difference<W, H>(frames.swept, frames.previous, {0, W / 2});
  const int tail =
      column_difference<W, H>(frames.swept, frames.current, {W / 2 - 1, W - 1});
  std::printf("  [%s] jump %.2f head %.2f tail %.2f steps\n", label,
              static_cast<double>(jump) / STEP,
              static_cast<double>(head) / STEP,
              static_cast<double>(tail) / STEP);
  HS_EXPECT_GE(jump, DF_ENDPOINT_JUMP_STEPS * STEP);
  HS_EXPECT_LE(head, DF_ENDPOINT_STEPS * STEP);
  HS_EXPECT_LE(tail, DF_ENDPOINT_STEPS * STEP);
}

/**
 * @brief A column-time frame's arm-boundary columns meet the frame endpoints
 *        they are lit at: the first columns show the previous frame and the
 *        last the current one, while the noise flows, balls fall and the
 *        stack turns through both fields.
 */
inline void test_displacement_field_arm_boundaries_meet_endpoints() {
  constexpr int W = SMALL_W, H = SMALL_H;
  check_displacement_field_endpoints<W, H>(
      "noise", displacement_field_endpoints<W, H>(200, [](auto &) {}));
  check_displacement_field_endpoints<W, H>(
      "balls", displacement_field_endpoints<W, H>(40, [](auto &effect) {
        HS_EXPECT_TRUE(effect.selectPreset(1));
        HS_EXPECT_EQ(effect.updateParameter("Ball Amp", 0.6f),
                     ParamSetResult::APPLIED);
      }));
}

/** @brief Single-pass residuals do not exceed the two-window residuals. */
inline void test_displacement_field_single_pass_reference() {
  constexpr int W = 288, H = 144;
  for (bool balls : {false, true}) {
    HS_CONTEXT(balls ? "balls" : "noise");
    std::vector<Pixel> images[3];
    const int frames = balls ? 40 : 200;
    for (int mode = 0; mode < 3; ++mode) {
      reset_effect_globals();
      hs::random().seed(1337u);
      pin_frame_clock(0);
      Spacetime::set_flat_surface(false);
      DisplacementField<W, H> effect;
      effect.init();
      if (balls)
        HS_EXPECT_TRUE(effect.selectPreset(1));
      for (int f = 0; f < frames; ++f) {
        DisplacementFieldWhiteBox::sweep_mode(
            effect, mode == 2 && f != frames - 1 ? 1 : mode);
        pin_frame_clock(f);
        effect.draw_frame();
        effect.advance_display();
      }
      if (mode == 0)
        DisplacementFieldWhiteBox::check_single_pass(effect);
      capture_frame<W, H>(effect, images[mode]);
    }
    std::array<int, 4> worst[2]{};
    for (int y = 0; y < H; ++y) {
      for (int x = 0; x < W; ++x) {
        const size_t index = static_cast<size_t>(y) * W + x;
        auto difference = [&](int a, int b) {
          const Pixel &p = images[a][index], &q = images[b][index];
          return std::max({std::abs(int(p.r) - q.r), std::abs(int(p.g) - q.g),
                           std::abs(int(p.b) - q.b)});
        };
        HS_EXPECT_LE(difference(0, 1), 1);
        for (int mode = 0; mode < 2; ++mode) {
          const int error = difference(mode, 2);
          worst[mode][0] = std::max(worst[mode][0], error);
          if ((y >= 1 && y <= 12) || (y >= 125 && y <= 142))
            worst[mode][1] = std::max(worst[mode][1], error);
          if (y >= 13 && y <= 124)
            worst[mode][2] = std::max(worst[mode][2], error);
          if (x == 0 || x == W / 2 - 1 || x == W / 2 || x == W - 1)
            worst[mode][3] = std::max(worst[mode][3], error);
        }
      }
    }
    for (int region = 0; region < 4; ++region) {
      HS_EXPECT_LE(worst[0][region], worst[1][region]);
      std::printf("  [reference region %d] window=%.5f single=%.5f steps\n",
                  region, worst[1][region] / 257.0, worst[0][region] / 257.0);
    }
  }
  hs::clear_mock_time();
}

/**
 * @brief A swept draw's frame turn carries the previous endpoint's points to
 *        the current one's however far the stack turns.
 */
inline void test_displacement_field_turn_follows_large_steps() {
  constexpr int W = SMALL_W, H = SMALL_H;
  reset_effect_globals();
  hs::random().seed(1337u);
  Spacetime::set_flat_surface(false);
  DisplacementField<W, H> effect;
  effect.init();
  for (int f = 0; f < 60; ++f) {
    pin_frame_clock(f);
    effect.draw_frame();
    effect.advance_display();
  }
  const math::Quaternion previous =
      math::make_rotation(math::Vector(0.3f, 0.8f, -0.5f).normalized(), 0.4f);
  for (const float angle : {0.1f, 0.6f, 2.0f, 3.0f}) {
    HS_CONTEXT("angle", angle);
    const math::Quaternion current =
        math::make_rotation(math::Vector(0.6f, -0.2f, 0.77f).normalized(),
                            angle) *
        previous;
    const math::Vector turn =
        DisplacementFieldWhiteBox::sweep_turn(effect, previous, current);
    HS_EXPECT_NEAR(sqrtf(math::dot(turn, turn)), angle, 1e-4f);
    for (const math::Vector &local :
         {math::X_AXIS, math::Y_AXIS, math::Vector(0.5f, -0.4f, 0.77f)}) {
      HS_EXPECT_VEC(Spacetime::rotate_by(math::rotate(local, previous), turn),
                    math::rotate(local, current), 2e-5f);
    }
  }
  hs::clear_mock_time();
}

/**
 * @brief A frame whose boundary knots overfill the alternate buffer scans
 *        its stack in successive groups, matching the two-window bakes.
 */
inline void test_displacement_field_alternate_groups_match_windows() {
  constexpr int W = 288, H = 144, FRAMES = 6;
  std::vector<Pixel> images[2];
  int flushes = 0;
  for (int mode = 0; mode < 2; ++mode) {
    reset_effect_globals();
    hs::random().seed(1337u);
    pin_frame_clock(0);
    Spacetime::set_flat_surface(false);
    DisplacementField<W, H> effect;
    effect.init();
    for (int f = 0; f < 200 + FRAMES; ++f) {
      DisplacementFieldWhiteBox::sweep_mode(effect, f < 200 ? 0 : mode);
      if (f >= 200)
        DisplacementFieldWhiteBox::arm_plane_pose(effect, f);
      pin_frame_clock(f);
      effect.draw_frame();
      effect.advance_display();
    }
    if (mode == 0)
      flushes = DisplacementFieldWhiteBox::stack_flushes(effect);
    capture_frame<W, H>(effect, images[mode]);
  }
  std::printf("  [alternate groups] scans=%d\n", flushes);
  HS_EXPECT_GT(flushes, 1);
  int worst = 0, lit = 0;
  for (size_t i = 0; i < images[0].size(); ++i) {
    const Pixel &p = images[0][i], &q = images[1][i];
    lit += !is_black(q);
    worst = std::max({worst, std::abs(int(p.r) - q.r), std::abs(int(p.g) - q.g),
                      std::abs(int(p.b) - q.b)});
  }
  HS_EXPECT_LE(worst, 1);
  HS_EXPECT_GT(lit, 1000);
  hs::clear_mock_time();
}

/** @brief Boundary copies cover neighboring spline and segment reads. */
inline void test_displacement_field_boundary_stencil() {
  constexpr int W = SMALL_W, H = SMALL_H, FRAMES = 80;
  std::vector<Pixel> primary;
  primary.reserve(FRAMES * W * H);
  int worst = 0;
  for (int mode = 0; mode < 2; ++mode) {
    reset_effect_globals();
    hs::random().seed(1337u);
    pin_frame_clock(0);
    Spacetime::set_flat_surface(false);
    DisplacementField<W, H> effect;
    effect.init();
    DisplacementFieldWhiteBox::sweep_mode(effect, mode);
    size_t index = 0;
    for (int f = 0; f < FRAMES; ++f) {
      DisplacementFieldWhiteBox::boundary_stencil_pose(effect);
      pin_frame_clock(f);
      effect.draw_frame();
      effect.advance_display();
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x) {
          const Pixel pixel = effect.get_pixel(x, y);
          if (mode == 0)
            primary.push_back(pixel);
          else {
            const Pixel &other = primary[index++];
            worst = std::max({worst, std::abs(int(pixel.r) - other.r),
                              std::abs(int(pixel.g) - other.g),
                              std::abs(int(pixel.b) - other.b)});
          }
        }
    }
  }
  HS_EXPECT_LE(worst, 1);
  std::printf("  [boundary stencil] worst=%d channel units\n", worst);
  hs::clear_mock_time();
}

/** @brief Manual preset cuts discard DisplacementField draw history. */
inline void test_displacement_field_retained_history_resets() {
  reset_effect_globals();
  DisplacementField<SMALL_W, SMALL_H> effect;
  effect.init();
  HS_EXPECT_FALSE(DisplacementFieldWhiteBox::history_ready(effect));
  run_frames(effect, 4);
  HS_EXPECT_TRUE(DisplacementFieldWhiteBox::history_ready(effect));
  HS_EXPECT_TRUE(effect.selectPreset(1));
  HS_EXPECT_FALSE(DisplacementFieldWhiteBox::history_ready(effect));
}
