/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Raymarch white-box cases. */

/** @brief White-box accessor for Raymarch's torus proportions. */
struct RaymarchWhiteBox {
  using RM = Raymarch<DEFAULT_W, DEFAULT_H>;
  static constexpr int MAX_POINTS = RM::MAX_POINTS;
  static constexpr float MAJOR = RM::MAJOR_K;
  static constexpr float MINOR = RM::MINOR_K;
  static constexpr float TWIST = RM::TWIST_K;
  static constexpr int TWIST_MAX = static_cast<int>(RM::TWIST_MAX);
  static constexpr float VIS = RM::VIS_K;
  static constexpr float BOUNDS = RM::UNIT_BOUNDS;

  template <int W, int H>
  static int volume_count(const Raymarch<W, H> &effect) {
    return effect.active_count;
  }

  template <int W, int H>
  static math::Quaternion volume_spin(const Raymarch<W, H> &effect, int index) {
    return effect.volume_spins[index].orientation.get();
  }

  template <int W, int H>
  static int animation_count(const Raymarch<W, H> &effect) {
    return effect.timeline.event_count();
  }

  template <int W, int H> static void refresh_points(Raymarch<W, H> &effect) {
    effect.refresh_points();
  }
  template <int W, int H> static void build_points(Raymarch<W, H> &effect) {
    effect.build_points();
  }
  template <int W, int H>
  static void set_base_solid(Raymarch<W, H> &effect,
                             RaymarchPlacementSolid solid) {
    effect.params.base_solid = solid;
  }

  template <int W, int H>
  static RaymarchPlacementSolid
  active_base_solid(const Raymarch<W, H> &effect) {
    return effect.active_base_solid;
  }

  template <int W, int H>
  static std::array<float, 7> surface_frame(const math::Vector &loc, int twist,
                                            float amplitude, float major_r,
                                            float minor_r) {
    SDF::WarpedVolume<SDF::Torus, SDF::Warp::Twist> torus{
        {major_r, minor_r}, {twist, amplitude, major_r}};
    const auto frame = Raymarch<W, H>::surface_frame(torus, loc);
    return {frame.normal.x, frame.normal.y, frame.normal.z, frame.cos_u,
            frame.sin_u,    frame.cos_v,    frame.sin_v};
  }

  /**
   * @brief Draws a frame through draw_fn from given previous endpoints, with
   *        the hue noise at its current bake.
   * @param fx Effect drawn.
   * @param previous_spin Spin phase the last draw used.
   * @param previous_palette Palette scroll the last draw used.
   * @return Displayed frame.
   */
  template <int W, int H>
  static std::vector<Pixel> draw_from(Raymarch<W, H> &fx, float previous_spin,
                                      float previous_palette) {
    fx.drawn_spin_phase = previous_spin;
    fx.drawn_palette_phase = previous_palette;
    fx.spin_drawn = true;
    fx.hue_noise_moved = false;
    {
      Canvas canvas(fx);
      fx.draw_fn(canvas);
    }
    fx.advance_display();
    std::vector<Pixel> out;
    capture_frame<W, H>(fx, out);
    return out;
  }

  /**
   * @brief Draws every torus placed, lit and scrolled as at one frame
   *        fraction, through the fixed-meridian scan.
   * @param fx Effect drawn.
   * @param u Frame fraction.
   * @param previous_spin Spin phase at the previous endpoint.
   * @param palette_step Palette scroll over the frame.
   * @return Displayed frame.
   */
  template <int W, int H>
  static std::vector<Pixel> draw_at(Raymarch<W, H> &fx, float u,
                                    float previous_spin, float palette_step) {
    using RM = Raymarch<W, H>;
    {
      Canvas canvas(fx);
      const auto sample = [u](const math::Quaternion &from,
                              const math::Quaternion &to) {
        return math::slerp(from, to, u);
      };
      const math::Quaternion camera = sample(fx.camera.get(0), fx.camera.get());
      const math::Quaternion spin =
          math::make_rotation(math::X_AXIS, fx.spin_phase * math::TWO_PI_F);
      const math::Quaternion previous_spin_q =
          math::make_rotation(math::X_AXIS, previous_spin * math::TWO_PI_F);
      for (int i = 0; i < fx.active_count; ++i) {
        const float outer = sinf(0.5f * fx.nn_angle[i] * fx.params.fill);
        const float scale = outer / RM::VIS_K;
        const float minor = scale * RM::MINOR_K;
        const float aa_width = minor * fx.params.aa_mult;
        const SDF::WarpedVolume<SDF::Torus, SDF::Warp::Twist> torus{
            {scale * RM::MAJOR_K, minor},
            {static_cast<int>(fx.params.twist + 0.5f), scale * RM::TWIST_K,
             scale * RM::MAJOR_K},
            2.0f * aa_width};
        const math::Quaternion q =
            sample(fx.camera.get(0) * fx.raw_quats[i] *
                       fx.volume_spins[i].orientation.get(0) * previous_spin_q,
                   fx.camera.get() * fx.raw_quats[i] *
                       fx.volume_spins[i].orientation.get() * spin);
        const math::Vector centre = math::rotate(fx.points[i], camera);
        const math::Vector half =
            blinn_phong_half(centre, math::rotate(math::X_AXIS, q));
        const float offset = fx.palette_phase +
                             static_cast<float>(i) / fx.active_count -
                             palette_step * (1.0f - u);
        Scan::TransformedVolume volume(torus, centre, q);
        Scan::Volume::draw<W, H>(
            fx.pipeline, canvas, centre, scale * RM::UNIT_BOUNDS + aa_width,
            volume,
            [&](const math::Vector &loc, Fragment &frag) {
              frag.color = fx.shade_hit(torus, loc, q, centre, half, offset);
            },
            fx.params.max_steps, aa_width);
      }
    }
    fx.advance_display();
    std::vector<Pixel> out;
    capture_frame<W, H>(fx, out);
    return out;
  }

  template <int W, int H> static float palette_phase(const Raymarch<W, H> &fx) {
    return fx.palette_phase;
  }
  template <int W, int H>
  static float drawn_palette_phase(const Raymarch<W, H> &fx) {
    return fx.drawn_palette_phase;
  }
  template <int W, int H> static float spin_phase(const Raymarch<W, H> &fx) {
    return fx.spin_phase;
  }
  template <int W, int H>
  static float hue_noise_phase(const Raymarch<W, H> &fx) {
    return fx.hue_noise_phase;
  }

  /**
   * @brief Rebakes the current hue-noise table at a phase.
   * @param fx Effect rebaked.
   * @param phase Hue-noise phase in turns.
   */
  template <int W, int H>
  static void bake_hue_noise(Raymarch<W, H> &fx, float phase) {
    fx.hue_noise_phase = phase;
    fx.refresh_hue_noise();
  }
};

/** @brief Channel bound of a Raymarch column against its oracle frame. */
inline constexpr int RAYMARCH_COLUMN_BOUND = 2;

/**
 * @brief Each column of a Raymarch frame marches, lights and scrolls the
 *        palette of every torus as at its sweep time, on both arms and across
 *        the palette scroll's wrap; the last draw's palette scroll is the
 *        next frame's previous endpoint.
 */
inline void test_raymarch_lighting_and_palette_at_column_time() {
  using WB = RaymarchWhiteBox;
  constexpr int W = SMALL_W, H = SMALL_H;
  reset_effect_globals();
  Raymarch<W, H> fx;
  fx.init();
  hs_test::run_frames(fx, 6);
  HS_EXPECT_EQ(WB::drawn_palette_phase(fx), WB::palette_phase(fx));

  constexpr float PALETTE_STEP = 0.04f;
  const float previous_spin = WB::spin_phase(fx) - 0.03f;
  const float previous_palette =
      math::wrap_t(WB::palette_phase(fx) - PALETTE_STEP);
  HS_EXPECT_GT(previous_palette, WB::palette_phase(fx));
  const std::vector<Pixel> drawn =
      WB::draw_from(fx, previous_spin, previous_palette);
  const std::vector<Pixel> expected = column_time_oracle<W, H>(
      [&](float u) { return WB::draw_at(fx, u, previous_spin, PALETTE_STEP); });
  HS_EXPECT_TRUE(expect_frames_within("raymarch columns", drawn, expected,
                                      RAYMARCH_COLUMN_BOUND) > 100);
}

/** @brief Largest sRGB8 channel gap a Raymarch arm-boundary column may show
 *  against its hue-noise oracle. */
inline constexpr int RAYMARCH_NOISE_ARM_BOUND = 2;
/** @brief Mean per-pixel sRGB8 gap a Raymarch frame may show against its
 *  hue-noise oracle. */
inline constexpr float RAYMARCH_NOISE_MEAN_BOUND = 2.0f;

/**
 * @brief Each column of a Raymarch frame samples the hue noise at its sweep
 *        time, starting from the bake the previous frame drew.
 */
inline void test_raymarch_hue_noise_at_column_time() {
  using WB = RaymarchWhiteBox;
  constexpr int W = SMALL_W, H = SMALL_H;
  reset_effect_globals();
  Raymarch<W, H> fx;
  fx.init();
  hs_test::run_frames(fx, 6);
  const float previous_spin = WB::spin_phase(fx);
  const float previous_palette = WB::palette_phase(fx);
  const float previous_noise = WB::hue_noise_phase(fx);
  fx.draw_frame();
  fx.advance_display();
  std::vector<Pixel> drawn;
  capture_frame<W, H>(fx, drawn);

  const float noise = WB::hue_noise_phase(fx);
  const float noise_step = std::remainder(noise - previous_noise, 1.0f);
  HS_EXPECT_GT(noise_step, 0.0f);
  const float palette_step =
      std::remainder(WB::palette_phase(fx) - previous_palette, 1.0f);
  const std::vector<Pixel> expected = column_time_oracle<W, H>([&](float u) {
    WB::bake_hue_noise(fx, math::wrap_t(previous_noise + noise_step * u));
    return WB::draw_at(fx, u, previous_spin, palette_step);
  });
  WB::bake_hue_noise(fx, noise);

  int arm_gap = 0;
  int lit = 0;
  double gap_sum = 0.0;
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x) {
      const Pixel &a = drawn[static_cast<size_t>(y) * W + x];
      const Pixel &b = expected[static_cast<size_t>(y) * W + x];
      if (is_black(a) && is_black(b))
        continue;
      const int gap = srgb8_channel_gap(a, b);
      ++lit;
      gap_sum += gap;
      if (x % (W / 2) == 0)
        arm_gap = std::max(arm_gap, gap);
    }
  HS_EXPECT_GT(lit, 100);
  HS_EXPECT_LE(arm_gap, RAYMARCH_NOISE_ARM_BOUND);
  HS_EXPECT_LE(gap_sum / std::max(lit, 1), RAYMARCH_NOISE_MEAN_BOUND);
}

/**
 * @brief Verifies every Raymarch volume receives its own live RandomWalk.
 */
inline void test_raymarch_volume_random_walks_are_independent() {
  reset_effect_globals();
  Raymarch<SMALL_W, SMALL_H> effect;
  effect.init();
  const int count = RaymarchWhiteBox::volume_count(effect);
  HS_EXPECT_EQ(count, 26);
  HS_EXPECT_EQ(RaymarchWhiteBox::animation_count(effect),
               RaymarchWhiteBox::MAX_POINTS + 3);

  effect.draw_frame();
  effect.advance_display();
  for (int i = 0; i < count; ++i)
    for (int j = i + 1; j < count; ++j) {
      const math::Quaternion a = RaymarchWhiteBox::volume_spin(effect, i);
      const math::Quaternion b = RaymarchWhiteBox::volume_spin(effect, j);
      const float dr = a.r - b.r;
      const math::Vector dv = a.v - b.v;
      HS_EXPECT_GT(dr * dr + math::dot(dv, dv), 1e-8f);
    }
}

/** @brief Pins Raymarch's named preset and selectable placement solids. */
inline void test_raymarch_preset_and_placement_solids() {
  using RM = Raymarch<SMALL_W, SMALL_H>;
  reset_effect_globals();
  RM effect;
  effect.init();

  auto value = [&](const char *name) {
    const auto *def = effect.getParameters().find(name);
    HS_EXPECT(def != nullptr, "Raymarch parameter is missing");
    return def ? def->get() : -1.0f;
  };

  HS_EXPECT_EQ(effect.getPresetCount(), 1u);
  HS_EXPECT_EQ(effect.getPresetIndex(), 0u);
  HS_EXPECT_EQ(RM::PRESET_IDS[0], std::string_view("uv-surface-noise"));

  hs_test::expect_params_in_range(effect, "Raymarch preset");

  const auto *base_solid = effect.getParameters().find("Base Solid");
  HS_EXPECT_TRUE(base_solid != nullptr);
  if (!base_solid)
    return;
  HS_EXPECT_TRUE(base_solid->animated);
  HS_EXPECT_TRUE(base_solid->preset);
  HS_EXPECT_EQ(base_solid->option_count,
               static_cast<int>(RM::PLACEMENT_SOLID_COUNT));
  HS_EXPECT_EQ(std::string_view(base_solid->options[0]), "Tetrahedron");
  HS_EXPECT_EQ(std::string_view(base_solid->options[17]),
               "Disdyakis Dodecahedron");
  HS_EXPECT_EQ(std::string_view(base_solid->options[20]),
               "Pentakis Dodecahedron");
  HS_EXPECT_EQ(std::string_view(base_solid->export_options[17]),
               "RaymarchPlacementSolid::DISDYAKIS_DODECAHEDRON");
  HS_EXPECT_EQ(value("Base Solid"), 17.0f);

  static constexpr std::array<int, RM::PLACEMENT_SOLID_COUNT> VERTEX_COUNTS{
      4,  8, 6,  20, 12, 12, 12, 24, 24, 24, 24,
      30, 8, 14, 14, 14, 26, 26, 32, 32, 32};
  for (size_t i = 0; i < RM::PLACEMENT_SOLID_COUNT; ++i) {
    HS_EXPECT_EQ(effect.updateParameter("Base Solid", static_cast<float>(i)),
                 ParamSetResult::APPLIED);
    RaymarchWhiteBox::refresh_points(effect);
    HS_EXPECT_EQ(RaymarchWhiteBox::volume_count(effect), VERTEX_COUNTS[i]);
    HS_EXPECT_EQ(
        static_cast<size_t>(RaymarchWhiteBox::active_base_solid(effect)), i);
  }

  HS_EXPECT_TRUE(effect.selectPreset(0));
  HS_EXPECT_EQ(value("Base Solid"), 17.0f);
  HS_EXPECT_NEAR(value("Hue Shift"), 0.76f, 1e-6f);
  effect.draw_frame();
  HS_EXPECT_EQ(RaymarchWhiteBox::volume_count(effect), 26);
}

/**
 * @brief Pins Raymarch surface coordinates to both periodic torus axes.
 */
inline void test_raymarch_surface_frame_uv() {
  constexpr float MAJOR_R = 0.6f;
  constexpr float MINOR_R = 0.12f;
  constexpr float AMPLITUDE = 0.08f;
  constexpr int TWIST = 3;
  const float u = 0.7f;
  const float v = -1.1f;
  const float cos_u = cosf(u);
  const float sin_u = sinf(u);
  const float cos_v = cosf(v);
  const float sin_v = sinf(v);
  const float radial = MAJOR_R + MINOR_R * cos_v;
  const math::Vector loc(radial * cos_u,
                         AMPLITUDE * sinf(TWIST * u) + MINOR_R * sin_v,
                         radial * sin_u);
  const auto frame = RaymarchWhiteBox::surface_frame<SMALL_W, SMALL_H>(
      loc, TWIST, AMPLITUDE, MAJOR_R, MINOR_R);
  HS_EXPECT_NEAR(frame[3], cos_u, 1e-5f);
  HS_EXPECT_NEAR(frame[4], sin_u, 1e-5f);
  HS_EXPECT_NEAR(frame[5], cos_v, 1e-5f);
  HS_EXPECT_NEAR(frame[6], sin_v, 1e-5f);

  SDF::WarpedVolume<SDF::Torus, SDF::Warp::Twist> torus{
      {MAJOR_R, MINOR_R}, {TWIST, AMPLITUDE, MAJOR_R}};
  const math::Vector expected_normal = torus.normal(loc);
  HS_EXPECT_NEAR(frame[0], expected_normal.x, 1e-5f);
  HS_EXPECT_NEAR(frame[1], expected_normal.y, 1e-5f);
  HS_EXPECT_NEAR(frame[2], expected_normal.z, 1e-5f);
  for (float radius : {0.0f, 0.5f * math::TOLERANCE, math::TOLERANCE,
                       2.0f * math::TOLERANCE}) {
    const auto frame = RaymarchWhiteBox::surface_frame<SMALL_W, SMALL_H>(
        math::Vector(radius, 0.1f, 0.0f), 2, 0.35f, 0.45f, 0.14f);
    for (float value : frame)
      HS_EXPECT_TRUE(std::isfinite(value));
  }
}

/**
 * @brief Pins the constexpr Newton square root behind UNIT_BOUNDS against libm.
 * @details A non-positive radicand short-circuits to 0.
 */
inline void test_raymarch_constexpr_sqrt_converges() {
  using WB = RaymarchWhiteBox;
  double worst = 0.0;
  for (int i = 1; i <= 80; ++i) {
    const float x = 0.05f * static_cast<float>(i); // 0.05 .. 4.0
    worst = hs_test::fold_worst(
        worst, std::fabs(static_cast<double>(math::constexpr_sqrt(x)) -
                         std::sqrt(static_cast<double>(x))));
  }
  HS_EXPECT_LT(worst, 1e-6);

  const float radicand = WB::MAJOR * WB::MAJOR + WB::TWIST * WB::TWIST;
  HS_EXPECT_NEAR(math::constexpr_sqrt(radicand), std::sqrt(radicand), 1e-7);

  // Radicands far outside O(1), where a fixed iteration count stops short.
  for (int i = 0; i < 6; ++i) {
    const float x = std::pow(10.0f, static_cast<float>(4 * i + 2));
    HS_CONTEXT("decade", static_cast<long long>(4 * i + 2));
    HS_EXPECT_NEAR_REL(static_cast<double>(math::constexpr_sqrt(x)),
                       std::sqrt(static_cast<double>(x)), 1e-6);
  }

  static_assert(math::constexpr_sqrt(400.0f) == 20.0f);
  static_assert(math::constexpr_sqrt(1.0e4f) == 100.0f);
  static_assert(math::constexpr_sqrt(0.0f) == 0.0f);
  static_assert(math::constexpr_sqrt(-1.0f) == 0.0f);
  HS_EXPECT_EQ(math::constexpr_sqrt(0.0f), 0.0f);
}

/**
 * @brief Verifies UNIT_BOUNDS really bounds the twisted tube it culls against.
 * @details The twisted torus surface is the tube circle of radius MINOR_K
 *          about the centerline (MAJOR_K, TWIST_K*sin(n*theta)) in the
 *          (xz-radius, y) half plane. Swept over the "Twist" slider's integer
 *          domain, it must match the production SDF's zero set and sit inside
 *          the sphere, and the sphere must be tight.
 */
inline void test_raymarch_unit_bounds_contains_twisted_tube() {
  using WB = RaymarchWhiteBox;
  constexpr double TWO_PI_DBL = 6.283185307179586;
  const double R = WB::MAJOR, r = WB::MINOR, A = WB::TWIST;
  const double bound = WB::BOUNDS;

  double max_off_surface = 0.0; // |SDF| at the sampled surface points
  double max_radius = 0.0;      // farthest surface point from the torus centre
  double min_on_shell = 1e30;   // least SDF value anywhere on the cull sphere

  for (int n = 0; n <= WB::TWIST_MAX; ++n) {
    const SDF::WarpedVolume<SDF::Torus, SDF::Warp::Twist> wv{
        SDF::Torus{static_cast<float>(R), static_cast<float>(r)},
        SDF::Warp::Twist{n, static_cast<float>(A), static_cast<float>(R)}};

    for (int i = 0; i < 256; ++i) {
      const double th = TWO_PI_DBL * i / 256;
      const double ct = std::cos(th), st = std::sin(th);
      const double y_mid = A * std::sin(n * th);
      for (int j = 0; j < 256; ++j) {
        const double a = TWO_PI_DBL * j / 256;
        const double s = R + r * std::cos(a);
        const double y = y_mid + r * std::sin(a);
        const math::Vector p(static_cast<float>(s * ct), static_cast<float>(y),
                             static_cast<float>(s * st));
        max_off_surface = hs_test::fold_worst(
            max_off_surface,
            std::fabs(static_cast<double>(wv.raw_distance(p))));
        max_radius = hs_test::fold_worst(max_radius, std::sqrt(s * s + y * y));
      }
    }

    // Sweep the cull sphere itself: no point of it may be inside the volume.
    for (int i = 0; i < 128; ++i) {
      const double th = TWO_PI_DBL * i / 128;
      for (int j = 0; j <= 64; ++j) {
        const double ph = 0.5 * TWO_PI_DBL * j / 64;
        const double sp = std::sin(ph);
        const math::Vector q(static_cast<float>(bound * sp * std::cos(th)),
                             static_cast<float>(bound * std::cos(ph)),
                             static_cast<float>(bound * sp * std::sin(th)));
        const double shell_distance = static_cast<double>(wv.raw_distance(q));
        if (std::isnan(shell_distance) || shell_distance < min_on_shell)
          min_on_shell = shell_distance;
      }
    }
  }

  HS_EXPECT_LT(max_off_surface, 2e-4); // the samples are the SDF's zero set
  HS_EXPECT_LE(max_radius, bound + 1e-6);
  HS_EXPECT_GT(min_on_shell, -1e-4);
  HS_EXPECT_NEAR(max_radius, bound, 1e-3);

  // The visible outer radius the per-vertex auto-size fits to the neighbour gap
  // is the ring rim, so it must sit inside the cull sphere.
  HS_EXPECT_LT(static_cast<double>(WB::VIS), bound);
}
