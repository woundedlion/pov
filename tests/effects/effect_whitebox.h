/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Per-effect white-box checks without a section of their own.
// ---------------------------------------------------------------------------

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
};

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

  for (const auto &def : effect.getParameters()) {
    HS_EXPECT_GE(def.get(), def.min);
    HS_EXPECT_LE(def.get(), def.max);
  }

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

/**
 * @brief White-box accessor for GnomonicStars' pixel-pitch star radius and its
 *        spiral cache.
 */
struct GnomonicStarsWhiteBox {
  template <int W, int H> static constexpr float radius_px() {
    return GnomonicStars<W, H>::radius_px();
  }
  template <int W, int H> static int max_points(const GnomonicStars<W, H> &) {
    return GnomonicStars<W, H>::MAX_POINTS;
  }
  template <int W, int H>
  static void set_points(GnomonicStars<W, H> &fx, float n) {
    HS_EXPECT_TRUE(fx.updateParameter("Points", n) == ParamSetResult::APPLIED);
  }
  template <int W, int H>
  static int cached_points(const GnomonicStars<W, H> &fx) {
    return fx.cached_points;
  }
  template <int W, int H>
  static math::Vector cache_at(const GnomonicStars<W, H> &fx, int i) {
    return fx.spiral_cache[i];
  }
  template <int W, int H>
  static void poison_cache(GnomonicStars<W, H> &fx, int i,
                           const math::Vector &v) {
    fx.spiral_cache[i] = v;
  }
};

/**
 * @brief Pins GnomonicStars' spiral-cache invalidation: the trig-heavy base
 *        lattice is rebuilt exactly when the clamped point count changes.
 * @details Parameter writes round and clamp the count before the cache sees it.
 */
inline void test_gnomonicstars_spiral_cache_invalidation() {
  using WB = GnomonicStarsWhiteBox;
  reset_effect_globals();
  GnomonicStars<SMALL_W, SMALL_H> fx;
  fx.init();

  auto step = [&](float points) {
    WB::set_points(fx, points);
    fx.draw_frame();
    fx.advance_display();
  };
  auto expect_lattice = [&](int n) {
    for (int i = 0; i < n; ++i) {
      const math::Vector want = math::fib_spiral(n, 0.5f, i);
      const math::Vector got = WB::cache_at(fx, i);
      HS_EXPECT_EQ(got.x, want.x);
      HS_EXPECT_EQ(got.y, want.y);
      HS_EXPECT_EQ(got.z, want.z);
    }
  };

  step(164.0f);
  HS_EXPECT_EQ(WB::cached_points(fx), 164);
  expect_lattice(164);

  // An unchanged count leaves the cache alone, so the poison survives.
  const math::Vector poison(0.0f, 0.0f, 1.0f);
  WB::poison_cache(fx, 7, poison);
  step(164.0f);
  HS_EXPECT_EQ(WB::cached_points(fx), 164);
  HS_EXPECT_EQ(WB::cache_at(fx, 7).z, poison.z);

  // A changed count rebuilds every slot, the poisoned one included.
  step(140.0f);
  HS_EXPECT_EQ(WB::cached_points(fx), 140);
  expect_lattice(140);

  // Out-of-range slider values clamp before they reach the cache.
  step(0.0f);
  HS_EXPECT_EQ(WB::cached_points(fx), 100);
  expect_lattice(100);
  step(static_cast<float>(WB::max_points(fx)) + 500.0f);
  HS_EXPECT_EQ(WB::cached_points(fx), WB::max_points(fx));
}

/** @brief Star radius covers the coarser row or column pitch. */
inline void test_gnomonicstars_radius_px_covers_both_axes() {
  using WB = GnomonicStarsWhiteBox;
  constexpr int W = DEFAULT_W, H = DEFAULT_H;
  const double pixel_pitch =
      std::max(math::TWO_PI_F / W, math::RADIANS_PER_ROW<H>);
  const double small_pixel_pitch =
      std::max(math::TWO_PI_F / SMALL_W, math::RADIANS_PER_ROW<SMALL_H>);
  const math::Basis basis = math::make_basis(math::Quaternion(), math::X_AXIS);

  for (int k : {1, 2, 7}) {
    const SDF::Star shape(basis, k * WB::radius_px<W, H>(), 5, 0.0f);
    HS_EXPECT_NEAR(shape.circumradius, k * pixel_pitch, 1e-6);
    const SDF::Star small(basis, k * WB::radius_px<SMALL_W, SMALL_H>(), 5,
                          0.0f);
    HS_EXPECT_NEAR(small.circumradius, k * small_pixel_pitch, 1e-6);
    HS_EXPECT_GE(small.circumradius / math::RADIANS_PER_ROW<SMALL_H>,
                 k - 1e-5f);
    HS_EXPECT_GE(small.circumradius / (math::TWO_PI_F / SMALL_W), k - 1e-5f);
  }

  constexpr int SPAN_PX = 12;
  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> pipe;
  {
    Canvas c(fx);
    Scan::Star::draw<W, H, false>(
        pipe, c, basis, SPAN_PX * WB::radius_px<W, H>(), /*sides=*/5,
        [](const math::Vector &, Fragment &f) {
          f.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
        });
  }
  fx.advance_display();

  size_t lit = 0;
  double max_arc = 0.0;
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x) {
      const Pixel &p = fx.get_pixel(x, y);
      if (p.r == 0 && p.g == 0 && p.b == 0)
        continue;
      ++lit;
      const math::Vector v = math::pixel_to_vector<W, H>(x, y);
      max_arc =
          std::max(max_arc, static_cast<double>(std::acos(hs::clamp(
                                math::dot(v, math::X_AXIS), -1.0f, 1.0f))));
    }

  HS_EXPECT_GT(lit, (size_t)0);
  // Tip reach uses the coarser pixel pitch, with AA and quantization slack.
  HS_EXPECT_LE(max_arc, (SPAN_PX + 2) * pixel_pitch);
  HS_EXPECT_GE(max_arc, (SPAN_PX - 2) * pixel_pitch);
}

/** @brief White-box accessor for Fishbowl's trail node. */
struct FishbowlWhiteBox {
  using EffectType = Fishbowl<SMALL_W, SMALL_H>;
  static size_t tween_vertices(const EffectType &fx) {
    size_t n = 0;
    deep_tween(fx.node->trail, [&](const math::Quaternion &, float) { ++n; });
    return n;
  }
  static Color4 sample_fire(const EffectType &fx, float t) {
    return fx.fire_palette.get(fx.duty_mod.modify(t));
  }
  static float palette_fill_scale(size_t trail_length) {
    return EffectType::palette_fill_scale(trail_length);
  }
  static float palette_phase_speed(float cycle_speed, float scale_factor,
                                   size_t trail_length) {
    return EffectType::palette_phase_speed(cycle_speed, scale_factor,
                                           trail_length);
  }
  static Color4 shade_trail(const EffectType &fx, float palette_t,
                            float age_t) {
    return fx.shade_trail(palette_t, age_t);
  }
  static bool needs_adaptive_midpoint(const math::Vector &a,
                                      const math::Vector &mid,
                                      const math::Vector &b) {
    return EffectType::needs_adaptive_midpoint(a, mid, b);
  }
};

/** @brief Verifies Fishbowl's sole preset and fire-palette duty window. */
inline void test_fishbowl_preset_and_fire_duty_cycle() {
  reset_effect_globals();
  using WB = FishbowlWhiteBox;
  WB::EffectType fx;
  fx.init();

  const auto snapshot = fx.serialize_parameters();
  auto edited = snapshot;
  edited.params.speed = 1.0f;
  HS_EXPECT_TRUE(fx.restore_parameters(edited));
  HS_EXPECT_EQ(fx.serialize_parameters().params.speed, 1.0f);
  edited.params.speed = std::numeric_limits<float>::quiet_NaN();
  HS_EXPECT_FALSE(fx.restore_parameters(edited));
  edited = snapshot;
  ++edited.schema_version;
  HS_EXPECT_FALSE(fx.restore_parameters(edited));
  HS_EXPECT_TRUE(fx.restore_parameters(snapshot));

  auto value = [&](const char *name) {
    const auto *def = fx.getParameters().find(name);
    HS_EXPECT(def != nullptr, "Fishbowl parameter is missing");
    return def ? def->get() : -1.0f;
  };

  for (const auto &def : fx.getParameters()) {
    HS_EXPECT_GE(def.get(), def.min);
    HS_EXPECT_LE(def.get(), def.max);
  }
  HS_EXPECT_EQ(fx.getPresetCount(), 1u);
  HS_EXPECT_EQ(fx.getPresetIndex(), 0u);

  HS_EXPECT_EQ(fx.updateParameter("Cycle Dur", 111.0f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(fx.updateParameter("Duty Cycle", 0.25f),
               ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(fx.selectPreset(0));
  HS_EXPECT_TRUE(fx.animations_paused());

  HS_EXPECT_EQ(value("Cycle Dur"), snapshot.params.cycle_duration);
  HS_EXPECT_EQ(value("Duty Cycle"), snapshot.params.duty_cycle);
  HS_EXPECT_EQ(fx.updateParameter("Duty Cycle", 0.5f), ParamSetResult::APPLIED);

  const Pixel black(0, 0, 0);
  const Color4 red = WB::sample_fire(fx, 0.20f);
  const Color4 yellow = WB::sample_fire(fx, 0.39f);
  HS_EXPECT_GT(red.color.r, red.color.g);
  HS_EXPECT_GT(red.color.r, red.color.b);
  HS_EXPECT_GT(yellow.color.g, red.color.g);
  HS_EXPECT_EQ(WB::sample_fire(fx, 0.50f).color, black);
  HS_EXPECT_EQ(WB::sample_fire(fx, 0.75f).color, black);
  const float dark_palette_t = 0.75f / value("Scale Factor");
  HS_EXPECT_EQ(WB::shade_trail(fx, dark_palette_t, 0.50f).alpha, 0.0f);

  HS_EXPECT_EQ(fx.updateParameter("Duty Cycle", 0.25f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(WB::sample_fire(fx, 0.30f).color, black);

  HS_EXPECT_EQ(WB::palette_fill_scale(1), 1.0f / WB::EffectType::TRAIL_LENGTH);
  HS_EXPECT_EQ(WB::palette_fill_scale(WB::EffectType::TRAIL_LENGTH / 2),
               static_cast<float>(WB::EffectType::TRAIL_LENGTH / 2) /
                   WB::EffectType::TRAIL_LENGTH);
  HS_EXPECT_EQ(WB::palette_fill_scale(WB::EffectType::TRAIL_LENGTH), 1.0f);
  const float sample_index = 12.0f;
  const float short_coord = (sample_index / 32.0f) * WB::palette_fill_scale(32);
  const float long_coord = (sample_index / 96.0f) * WB::palette_fill_scale(96);
  HS_EXPECT_NEAR(short_coord, long_coord, 1e-7f);

  const float cycle_speed = value("Cycle Speed");
  const float scale_factor = value("Scale Factor");
  const float coordinate_shift = scale_factor / WB::EffectType::TRAIL_LENGTH;
  const float filling_speed = WB::palette_phase_speed(
      cycle_speed, scale_factor, WB::EffectType::TRAIL_LENGTH - 1);
  const float saturated_speed = WB::palette_phase_speed(
      cycle_speed, scale_factor, WB::EffectType::TRAIL_LENGTH);
  HS_EXPECT_NEAR(filling_speed, cycle_speed - coordinate_shift, 1e-7f);
  HS_EXPECT_NEAR(saturated_speed - coordinate_shift, filling_speed, 1e-7f);

  const math::Vector a = math::X_AXIS;
  const math::Vector b = math::Y_AXIS;
  HS_EXPECT_FALSE(WB::needs_adaptive_midpoint(a, math::slerp(a, b, 0.5f), b));
  HS_EXPECT_TRUE(WB::needs_adaptive_midpoint(a, math::Z_AXIS, b));
}

/**
 * @brief Verifies the scratch-A split's predicted worst case really covers a
 *        saturated Fishbowl frame.
 * @details SCRATCH_A_BYTES is sized against a closed-form worst case: the
 *          MAX_FRAGMENTS vertex buffer, the Multiline fragment buffer it binds,
 *          and rasterize's sub-step cache, all live at once.
 */
inline void test_fishbowl_scratch_estimate_covers_peak() {
  reset_effect_globals();
  using WB = FishbowlWhiteBox;
  using EffectType = WB::EffectType;

  EffectType fx;
  fx.init();
  scratch_arena_a.reset_high_water_mark();

  size_t worst_vertices = 0;
  for (int i = 0; i < EffectType::TRAIL_LENGTH + 4; ++i) {
    fx.draw_frame();
    fx.advance_display();
    worst_vertices = std::max(worst_vertices, WB::tween_vertices(fx));
  }

  constexpr size_t PREDICTED = EffectType::SCRATCH_A_ESTIMATE;
  HS_EXPECT_LE(scratch_arena_a.get_high_water_mark(), PREDICTED);
  // The trail fills, so the peak above is a saturated frame and not a warm-up.
  HS_EXPECT_GT(worst_vertices, (size_t)EffectType::TRAIL_LENGTH);
  HS_EXPECT_LE(worst_vertices, (size_t)EffectType::MAX_FRAGMENTS);
}

/** @brief White-box accessor for RingSpin's ring pool. */
struct RingSpinWhiteBox {
  using RS = RingSpin<DEFAULT_W, DEFAULT_H>;

  /** @brief Sub-frame count of ring @p i's current orientation step. */
  static int substeps(const RS &fx, int i) {
    return fx.rings[i].orientation.length();
  }
  /** @brief World-space axis of ring @p i's great circle at sub-frame @p s. */
  static math::Vector axis(const RS &fx, int i, int s) {
    return fx.rings[i].orientation.orient(math::Y_AXIS, s).normalized();
  }
  /** @brief Ring-pool size. */
  static int num_rings() { return RS::NUM_RINGS; }
  /** @brief Half-width in radians of the trail-head stroke. */
  static float head_half_width(const RS &fx) {
    const float pixel_w = math::coarse_pixel_pitch<DEFAULT_W, DEFAULT_H>();
    return 2.0f * pixel_w * fx.params.thickness;
  }
};

/** @brief Strobe configuration preserves a nonblack RingSpin render. */
inline void test_ringspin_strobe_configuration_preserves_rendering() {
  std::vector<Pixel> expected;
  bool lit = false;
  render_capture<RingSpin, SMALL_W, SMALL_H>(expected, 2, nullptr, &lit);
  HS_EXPECT_TRUE(lit);
  reset_effect_globals();
  pin_frame_clock(0);
  {
    RingSpin<SMALL_W, SMALL_H> roster_effect;
    HS_EXPECT_TRUE(roster_effect.strobe_columns());
  }
  reset_effect_globals();
  pin_frame_clock(0);
  RingSpin<SMALL_W, SMALL_H> firmware_effect(false);
  HS_EXPECT_FALSE(firmware_effect.strobe_columns());
  firmware_effect.init();
  for (int frame = 0; frame < 2; ++frame) {
    pin_frame_clock(frame);
    firmware_effect.draw_frame();
    firmware_effect.advance_display();
  }
  HS_EXPECT_EQ(expected.size(), size_t{SMALL_W * SMALL_H});
  for (int y = 0; y < SMALL_H; ++y)
    for (int x = 0; x < SMALL_W; ++x) {
      const Pixel &actual = firmware_effect.get_pixel(x, y);
      const Pixel &reference = expected[y * SMALL_W + x];
      HS_EXPECT_EQ(actual.r, reference.r);
      HS_EXPECT_EQ(actual.g, reference.g);
      HS_EXPECT_EQ(actual.b, reference.b);
    }
}

/**
 * @brief Verifies every lit RingSpin pixel sits on one of its great circles.
 * @details On the first frame, where the trail holds a single orientation
 *          step, every lit pixel must lie within the stroke half-width of some
 *          ring's equator. Alpha 0 must then blank the frame.
 */
inline void test_ringspin_trail_hugs_its_great_circles() {
  reset_effect_globals();
  pin_frame_clock(0);
  using WB = RingSpinWhiteBox;
  WB::RS fx;
  fx.init();
  pin_frame_clock(1);
  fx.draw_frame();
  fx.advance_display();

  // Polar deviation from the equator that the stroke half-width admits.
  const float band = std::sin(WB::head_half_width(fx));
  int lit = 0, off_circle = 0;
  float worst = 0.0f;
  for (int y = 0; y < DEFAULT_H; ++y)
    for (int x = 0; x < DEFAULT_W; ++x) {
      const Pixel p = fx.get_pixel(x, y);
      if (p.r == 0 && p.g == 0 && p.b == 0)
        continue;
      ++lit;
      const math::Vector v = math::pixel_to_vector<DEFAULT_W, DEFAULT_H>(x, y);
      float nearest = 1.0f;
      for (int i = 0; i < WB::num_rings(); ++i)
        for (int s = 0; s < WB::substeps(fx, i); ++s)
          nearest =
              std::min(nearest, std::abs(math::dot(v, WB::axis(fx, i, s))));
      worst = hs_test::fold_worst(worst, nearest);
      if (nearest > band)
        ++off_circle;
    }

  std::printf("  [info] RingSpin: lit=%d off-circle=%d worst |dot|=%.6f band "
              "%.6f\n",
              lit, off_circle, worst, band);
  HS_EXPECT_GT(lit, 0);
  HS_EXPECT_EQ(off_circle, 0);

  HS_EXPECT_EQ(fx.updateParameter("Alpha", 0.0f), ParamSetResult::APPLIED);
  fx.draw_frame();
  fx.advance_display();
  int still_lit = 0;
  for (int y = 0; y < DEFAULT_H; ++y)
    for (int x = 0; x < DEFAULT_W; ++x) {
      const Pixel p = fx.get_pixel(x, y);
      if (p.r != 0 || p.g != 0 || p.b != 0)
        ++still_lit;
    }
  HS_EXPECT_EQ(still_lit, 0);
}

// PetalFlow spawn-gap, hue-cursor and ring-pool bounds on every frame.

/**
 * @brief White-box accessor for PetalFlow's spawn-gap accumulator and hue
 *        cursor.
 */
struct PetalFlowWhiteBox {
  using PF = PetalFlow<DEFAULT_W, DEFAULT_H>;
  static float gap(const PF &pf) { return pf.gap_accumulator; }
  static float start_rho() { return PF::START_RHO; }
  static float speed_max() { return PF::SPEED_MAX; }
  static float density_max() { return PF::DENSITY_MAX; }
  static float youngest_rho(const PF &pf) {
    float youngest = PF::END_RHO;
    for (const auto &ring : pf.rings)
      if (ring.active && ring.rho < youngest)
        youngest = ring.rho;
    return youngest;
  }
  static float next_hue(const PF &pf) { return pf.next_hue; }
  static float live_spacing(const PF &pf) { return pf.spacing(); }
  static float move_dist(const PF &pf) { return pf.move_dist(); }
  static int spawns(const PF &pf) { return pf.spawns; }
  static uint32_t pool_full_drops(const PF &pf) { return pf.pool_full_drops; }
  static int max_rings() { return PF::MAX_RINGS; }
  static int rings_on_path() { return PF::RINGS_ON_PATH; }
  static int active_rings(const PF &pf) {
    int active = 0;
    for (int i = 0; i < PF::MAX_RINGS; ++i)
      if (pf.rings[i].active)
        ++active;
    return active;
  }
};

/**
 * @brief Verifies the spawn-gap accumulator drains every frame, the hue cursor
 *        stays wrapped, and the ring pool absorbs the worst slider corner.
 * @details After any frame the residue must satisfy 0 <= gap < spacing() and
 *          next_hue must stay in [0, 1).
 */
inline void test_petalflow_spawn_gap_bounded() {
  using WB = PetalFlowWhiteBox;
  reset_effect_globals();
  WB::PF pf;
  pf.init();
  HS_EXPECT_NEAR(WB::youngest_rho(pf) - WB::start_rho(), WB::gap(pf), 1e-5f);
  HS_EXPECT_EQ(pf.updateParameter("Speed", WB::speed_max()),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(pf.updateParameter("Density", WB::density_max()),
               ParamSetResult::APPLIED);

  const float spacing = WB::live_spacing(pf);
  // Below this the while-loop emits at most one ring per frame, whatever the
  // frame count.
  HS_EXPECT_GT(WB::move_dist(pf), spacing);

  // Long enough for the path prefilled at the default density to be fully
  // replaced at the tightest one, so the pool bound is read at steady state
  // rather than mid-fill.
  const int frames = smoke_frames() < 128 ? 128 : smoke_frames();
  int worst_active = 0;
  int worst_burst = 0;
  int previous_spawns = WB::spawns(pf);
  for (int f = 0; f < frames; ++f) {
    pf.draw_frame();
    pf.advance_display();
    const float gap = WB::gap(pf);
    HS_EXPECT_GE(gap, 0.0f);
    HS_EXPECT_LT(gap, spacing); // accumulator never runs away past one spacing
    const float hue = WB::next_hue(pf);
    HS_EXPECT_GE(hue, 0.0f);
    HS_EXPECT_LT(hue, 1.0f); // hue cursor stays wrapped
    const int spawns = WB::spawns(pf);
    worst_burst = std::max(worst_burst, spawns - previous_spawns);
    previous_spawns = spawns;
    worst_active = std::max(worst_active, WB::active_rings(pf));
    // every spawn found a free slot
    HS_EXPECT_EQ(WB::pool_full_drops(pf), (uint32_t)0);
  }
  HS_EXPECT_GE(worst_burst, 2); // the multi-spawn branch actually ran
  HS_EXPECT_LE(worst_active, WB::rings_on_path());
  HS_EXPECT_LE(WB::rings_on_path(), WB::max_rings());
}

/**
 * @brief Verifies glitch_lens maps unit directions to unit directions and
 *        stays smooth through its equatorial fold.
 * @details Pins |lens(v)| == 1 across a spread of directions, its doubled-
 *          latitude topology, and matching one-sided derivatives at the fold.
 */
inline void test_glitch_lens_unit_norm() {
  for (float polar : {0.4f, 0.9f, 2.0f}) {
    for (float azimuth : {-2.1f, -0.4f, 0.7f, 1.8f}) {
      const math::Vector input(sinf(polar) * cosf(azimuth), cosf(polar),
                               sinf(polar) * sinf(azimuth));
      const math::Vector mapped = lenses::glitch_lens(input);
      HS_EXPECT_NEAR(mapped.x, sinf(2 * polar) * cosf(3 * azimuth), 2e-6f);
      HS_EXPECT_NEAR(mapped.y, cosf(2 * polar), 2e-6f);
      HS_EXPECT_NEAR(mapped.z, sinf(2 * polar) * sinf(3 * azimuth), 2e-6f);
    }
  }
  const math::Vector dirs[] = {math::Vector(1, 0, 0),
                               math::Vector(0, 0, 1),
                               math::Vector(-1, 0, 0),
                               math::Vector(0, 0, -1),
                               math::Vector(1, 1, 1).normalized(),
                               math::Vector(-1, 2, -3).normalized(),
                               math::Vector(3, -1, 2).normalized(),
                               math::Vector(0.2f, -0.9f, 0.4f).normalized(),
                               math::Vector(-0.7f, 0.1f, 0.7f).normalized()};
  for (const math::Vector &v : dirs) {
    HS_EXPECT_NEAR(lenses::glitch_lens(v).length(), 1.0f, 1e-3f);
  }

  const math::Vector equator_x = lenses::glitch_lens(math::Vector(1, 0, 0));
  HS_EXPECT_NEAR(equator_x.x, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(equator_x.y, -1.0f, 1e-6f);
  HS_EXPECT_NEAR(equator_x.z, 0.0f, 1e-6f);

  const math::Vector north = lenses::glitch_lens(math::Vector(0, 1, 0));
  HS_EXPECT_NEAR(north.y, 1.0f, 1e-6f);
  const math::Vector south = lenses::glitch_lens(math::Vector(0, -1, 0));
  HS_EXPECT_NEAR(south.y, 1.0f, 1e-6f);

  constexpr float EPSILON = 1e-4f;
  const float radial = sqrtf(1.0f - EPSILON * EPSILON);
  const math::Vector center =
      lenses::glitch_lens(math::Vector(0.8f, 0.0f, 0.6f));
  const math::Vector above =
      lenses::glitch_lens(math::Vector(0.8f * radial, EPSILON, 0.6f * radial));
  const math::Vector below =
      lenses::glitch_lens(math::Vector(0.8f * radial, -EPSILON, 0.6f * radial));
  const math::Vector forward = above - center;
  const math::Vector backward = center - below;
  HS_EXPECT_NEAR(forward.x, backward.x, 1e-6f);
  HS_EXPECT_NEAR(forward.y, backward.y, 1e-6f);
  HS_EXPECT_NEAR(forward.z, backward.z, 1e-6f);
}

/** @brief White-box accessor for MobiusRings' singular numeric branches. */
struct MobiusRingsWhiteBox {
  using MR = MobiusRings<DEFAULT_W, DEFAULT_H>;
  static float conformal_coord(float z, float phase) {
    return MR::conformal_coord(z, phase);
  }
  static math::Quaternion counter_rotation(const math::Vector &mid) {
    return MR::counter_rotation(mid);
  }
};

inline void test_mobius_rings_conformal_and_counter_rotation() {
  using WB = MobiusRingsWhiteBox;

  const float zs[] = {-1.0f, -0.999999f, -0.5f, 0.0f, 0.5f, 0.999999f, 1.0f};
  for (float z : zs) {
    for (float phase : {0.0f, 0.25f, 0.6f, 0.95f}) {
      float coord = WB::conformal_coord(z, phase);
      HS_EXPECT_TRUE(std::isfinite(coord));
      HS_EXPECT_GE(coord, 0.0f);
      HS_EXPECT_LE(coord, 1.0f);
    }
  }
  HS_EXPECT_NEAR(WB::conformal_coord(1.0f, 0.3f), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(WB::conformal_coord(-1.0f, 0.7f), 1.0f, 1e-6f);

  HS_EXPECT_EQ(WB::counter_rotation(math::Vector(0.0f, 0.0f, 0.0f)),
               math::Quaternion());

  math::Vector mid(0.3f, -0.7f, 0.4f);
  math::Vector r = math::rotate(mid.normalized(), WB::counter_rotation(mid));
  HS_EXPECT_NEAR(r.x, 0.0f, 1e-3f);
  HS_EXPECT_NEAR(r.y, 0.0f, 1e-3f);
  HS_EXPECT_NEAR(r.z, 1.0f, 1e-3f);
}

/**
 * @brief Checks the full-frame query on four constructed effects.
 * @details MeshFeedback and Dynamo fold any_crosses_segments true; the two
 *          sampled non-crossing effects keep segment clipping.
 */
inline void test_needs_full_frame_gate() {
  // Each effect aliases the same static double buffer, so construct one at a
  // time.
  auto reset = [] { reset_effect_globals(); };
  auto check = [&](Effect &fx, bool expected, const char *name) {
    fx.init();
    if (fx.needs_full_frame() != expected)
      std::printf("  [FAIL] %s needs_full_frame()=%d expected=%d\n", name,
                  fx.needs_full_frame(), expected);
    HS_EXPECT_EQ(fx.needs_full_frame(), expected);
  };

  // Cross-segment stateful effects -> full-frame.
  {
    reset();
    MeshFeedback<DEFAULT_W, DEFAULT_H> fx;
    check(fx, true, "MeshFeedback");
  }
  {
    reset();
    Dynamo<DEFAULT_W, DEFAULT_H> fx;
    check(fx, true, "Dynamo");
  }

  // Representative non-stateful effects -> keep band clipping (default false).
  {
    reset();
    Voronoi<DEFAULT_W, DEFAULT_H> fx;
    check(fx, false, "Voronoi");
  }
  {
    reset();
    RingSpin<DEFAULT_W, DEFAULT_H> fx;
    check(fx, false, "RingSpin");
  }
}

/**
 * @brief Bounds whole-solid generation, classification and rendering scratch
 * against HankinSolids' exported budgets at the device height.
 */
inline void test_hankinsolids_arena_budget_covers_every_solid() {
  constexpr int W = 288, H = 144;
  constexpr size_t SCRATCH_A = HankinSolids<W, H>::SCRATCH_A_BYTES;
  constexpr size_t SCRATCH_B = HankinSolids<W, H>::SCRATCH_B_BYTES;
  constexpr size_t MEASURE = 1024 * 1024; // headroom so a peak never traps here
  constexpr float ANGLE = math::PI_F / 4.0f;

  auto solids = Solids::Collections::get_simple_solids();
  for (size_t idx = 0; idx < solids.size(); ++idx) {
    configure_arenas(GLOBAL_ARENA_SIZE - 2 * MEASURE, MEASURE, MEASURE);

    MeshPaletteBank palette_bank;
    palette_bank.bake_all(persistent_arena);

    // The effect's held graph-walk seed; dodecahedron is the largest Platonic.
    PolyMesh seed;
    hs::generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      seed =
          Solids::finalize_solid(Solids::Platonic::dodecahedron(a, b), target);
    });

    MeshState mesh;
    CompiledHankin hankin;
    hs::generate(persistent_arena, [&](Arena &target, Arena &a, Arena &b) {
      PolyMesh base = Solids::finalize_solid(solids[idx].generate(a, b), a);
      hankin = CompiledHankin();
      MeshOps::compile_hankin(base, hankin, target, a);
      mesh.clear();
      MeshOps::update_hankin(hankin, mesh, target, ANGLE);
    });
    {
      ScratchScope a_guard(scratch_arena_a);
      ScratchScope b_guard(scratch_arena_b);
      MeshOps::classify_faces_by_topology(mesh, scratch_arena_a,
                                          scratch_arena_b, persistent_arena);
    }

    // Render peak: transform into scratch_a, then Scan::Mesh::draw stacks a
    // FaceScratchBuffer on top.
    {
      ScratchScope a_guard(scratch_arena_a);
      math::Orientation<> orientation;
      OrientTransformer camera(orientation);
      MeshState rotated;
      MeshOps::transform(mesh, rotated, scratch_arena_a, camera);
      hs_test::StubEffect fx(W, H);
      Canvas canvas(fx);
      Pipeline<W, H> filters;
      auto frag = [](const math::Vector &, Fragment &f) {
        f.color = Color4(Pixel(1000, 1000, 1000), 1.0f);
      };
      Scan::Mesh::draw<W, H>(filters, canvas, rotated, frag, scratch_arena_a);
    }

    // Morph compaction peak: the CompiledHankin + palette bank survive into
    // scratch_b, the mesh + walk seed into scratch_a, then persistent is
    // reset.
    {
      Persist<CompiledHankin> ph(hankin, scratch_arena_b, persistent_arena);
      Persist<MeshState> pf(mesh, scratch_arena_a, persistent_arena);
      Persist<MeshPaletteBank> pp(palette_bank, scratch_arena_b,
                                  persistent_arena);
      Persist<PolyMesh> ps(seed, scratch_arena_a, persistent_arena);
      persistent_arena.reset();
    }

    const size_t a_peak = scratch_arena_a.get_high_water_mark();
    const size_t b_peak = scratch_arena_b.get_high_water_mark();
    if (a_peak > SCRATCH_A || b_peak > SCRATCH_B)
      std::printf("  HankinSolids arena OVER BUDGET solid[%zu] '%s': "
                  "scratchA=%zu/%zu scratchB=%zu/%zu\n",
                  idx, solids[idx].name, a_peak, SCRATCH_A, b_peak, SCRATCH_B);
    HS_EXPECT_TRUE(a_peak <= SCRATCH_A);
    HS_EXPECT_TRUE(b_peak <= SCRATCH_B);
  }
  configure_arenas_default();
}

template <typename EffectT> inline void check_manual_preset_navigation() {
  constexpr size_t expected_count = EffectT::authored_preset_count();
  HS_EXPECT_GT(expected_count, size_t(1));
  reset_effect_globals();
  EffectT effect;
  effect.init();
  HS_EXPECT_EQ(effect.getPresetCount(), expected_count);
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t(0));
  HS_EXPECT_FALSE(effect.selectPreset(expected_count));
  HS_EXPECT_FALSE(effect.animations_paused());
  HS_EXPECT_TRUE(effect.previousPreset());
  HS_EXPECT_EQ(effect.getPresetIndex(), expected_count - 1);
  HS_EXPECT_TRUE(effect.nextPreset());
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t(0));
  HS_EXPECT_TRUE(effect.animations_paused());
}

inline void test_manual_preset_navigation() {
  check_manual_preset_navigation<MindSplatter<SMALL_W, SMALL_H>>();
  check_manual_preset_navigation<DreamBalls<SMALL_W, SMALL_H>>();
  check_manual_preset_navigation<Comets<SMALL_W, SMALL_H>>();
  check_manual_preset_navigation<MeshFeedback<SMALL_W, SMALL_H>>();
  check_manual_preset_navigation<ShapeShifter<SMALL_W, SMALL_H>>();
}
