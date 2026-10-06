/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_effects.h.

// ---------------------------------------------------------------------------
// Geometry, caches, presets and arena budgets across the effect roster.
// ---------------------------------------------------------------------------

/**
 * @brief White-box accessor for Raymarch's torus proportions (befriended in
 *        effects/Raymarch.h).
 */
struct RaymarchWhiteBox {
  using RM = Raymarch<DEFAULT_W, DEFAULT_H>;
  static constexpr int MAX_POINTS = RM::MAX_POINTS;
  static constexpr float MAJOR = RM::MAJOR_K;
  static constexpr float MINOR = RM::MINOR_K;
  static constexpr float TWIST = RM::TWIST_K;
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
  const math::Quaternion first = RaymarchWhiteBox::volume_spin(effect, 0);
  int distinct = 0;
  for (int i = 1; i < count; ++i) {
    const math::Quaternion q = RaymarchWhiteBox::volume_spin(effect, i);
    const float dr = q.r - first.r;
    const math::Vector dv = q.v - first.v;
    if (dr * dr + math::dot(dv, dv) > 1e-8f)
      ++distinct;
  }
  HS_EXPECT_EQ(distinct, count - 1);
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
 * @details It runs only at compile time, so nothing else would notice it
 *          under-converging, and its consumers are the Raymarch cull-sphere
 *          radius — where a low answer culls real surface — and Voronoi's
 *          coherence-block floor at MAX_SITES. A non-positive radicand
 *          short-circuits to 0.
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
 * @details Raymarch hands scale*UNIT_BOUNDS (+ the AA pad) to Scan::Volume::draw
 *          as the ray-march cull sphere, so any surface outside it is dropped
 *          with no other tell. The twisted torus surface has a closed form —
 *          the tube circle of radius MINOR_K about the centerline
 *          (MAJOR_K, TWIST_K*sin(n*theta)) in the (xz-radius, y) half plane —
 *          so the whole surface is swept over the "Twist" slider's integer
 *          domain, checked against the production SDF's zero set, and required
 *          to sit inside the sphere. The sphere is also required to be tight: a
 *          slack radius is wasted ray steps at every vertex.
 */
inline void test_raymarch_unit_bounds_contains_twisted_tube() {
  using WB = RaymarchWhiteBox;
  constexpr double TWO_PI_DBL = 6.283185307179586;
  const double R = WB::MAJOR, r = WB::MINOR, A = WB::TWIST;
  const double bound = WB::BOUNDS;

  double max_off_surface = 0.0; // |SDF| at the sampled surface points
  double max_radius = 0.0;      // farthest surface point from the torus centre
  double min_on_shell = 1e30;   // least SDF value anywhere on the cull sphere

  // twist_n rounds the "Twist" slider, whose range is 0..8.
  for (int n = 0; n <= 8; ++n) {
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
 *        spiral cache (befriended in effects/GnomonicStars.h).
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

/**
 * @brief White-box accessor for Fishbowl's trail node (befriended in
 *        effects/Fishbowl.h).
 */
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
  HS_EXPECT_EQ(WB::EffectType::MAX_FRAGMENTS,
               2 * WB::EffectType::TRAIL_LENGTH *
                   WB::EffectType::ORIENTATION_SUBSTEPS);
}

/**
 * @brief Verifies the scratch-A split's predicted worst case really covers a
 *        saturated Fishbowl frame.
 * @details SCRATCH_A_BYTES is carved from the global arena against a closed-form
 *          worst case — the MAX_FRAGMENTS vertex buffer, the Multiline fragment
 *          buffer it binds, and rasterize's sub-step cache, all live at once.
 *          The static_assert only checks that estimate against the split, never
 *          against a real frame, so an estimate that under-counts would go
 *          unnoticed until the split shrank. Runs past TRAIL_LENGTH frames so
 *          the trail is full when the peak is read.
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

/**
 * @brief White-box accessor for RingSpin's ring pool (befriended in
 *        effects/RingSpin.h).
 */
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
    const float pixel_w =
        std::max(math::TWO_PI_F / DEFAULT_W, math::RADIANS_PER_ROW<DEFAULT_H>);
    return 2.0f * pixel_w * fx.params.thickness;
  }
};

inline void test_ringspin_strobe_configuration_preserves_rendering() {
  std::vector<Pixel> expected;
  render_capture<RingSpin, 96, 20>(expected, 2);
  reset_effect_globals();
  pin_frame_clock(0);
  {
    RingSpin<96, 20> roster_effect;
    HS_EXPECT_TRUE(roster_effect.strobe_columns());
  }
  reset_effect_globals();
  pin_frame_clock(0);
  RingSpin<96, 20> firmware_effect(false);
  HS_EXPECT_FALSE(firmware_effect.strobe_columns());
  firmware_effect.init();
  for (int frame = 0; frame < 2; ++frame) {
    pin_frame_clock(frame);
    firmware_effect.draw_frame();
    firmware_effect.advance_display();
  }
  HS_EXPECT_EQ(expected.size(), size_t{96 * 20});
  for (int y = 0; y < 20; ++y)
    for (int x = 0; x < 96; ++x) {
      const Pixel &actual = firmware_effect.get_pixel(x, y);
      const Pixel &reference = expected[y * 96 + x];
      HS_EXPECT_EQ(actual.r, reference.r);
      HS_EXPECT_EQ(actual.g, reference.g);
      HS_EXPECT_EQ(actual.b, reference.b);
    }
}

/**
 * @brief Verifies every lit RingSpin pixel sits on one of its great circles.
 * @details RingSpin draws SDF::Ring at radius 1 (a great circle) about each
 *          ring's rotated plane normal, so on the first frame — where the trail
 *          holds a single orientation step — every lit pixel must lie within the
 *          stroke half-width of some ring's equator. This catches a basis,
 *          radius or trail-record regression that still renders a plausible
 *          frame. Alpha 0 must then blank the frame, pinning the trail-colour
 *          early-out.
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

  // Polar deviation from the equator that the stroke half-width admits. The
  // pipeline carries no Screen::AntiAlias, so there is no fringe to allow for;
  // the measured worst case is 0.011 against this 0.035 bound.
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

// White-box bounds for spawn gaps, emit phases and pool indices on every frame.

/**
 * @brief White-box accessor for PetalFlow's spawn-gap accumulator and hue cursor
 *        (befriended in effects/PetalFlow.h).
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
  static float spacing() { return PF::SPACING; }
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
 * @details check_spawn() integrates Speed*RHO_PER_SPEED into gap_accumulator and
 *          drains it by the live spacing per spawn, so after any frame the
 *          residue must satisfy 0 <= gap < spacing() — a runaway (missing drain)
 *          or a negative residue both fail here. Run at Speed_max x Density_max:
 *          the only corner where a frame's travel exceeds the live spacing, so
 *          the while-loop's multi-spawn branch runs, and the corner RINGS_ON_PATH
 *          is derived for, where the pool bound has no margin left — a dropped
 *          spawn there is a real defect, not a rounding allowance. next_hue is
 *          advanced wrap(.,1) per spawn and must stay [0, 1).
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
 * @brief White-box accessor for DisplacementField's hue-table bake.
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
    // Slot 0 of the pooled bake: the first (only, in the one-ring tests)
    // drawn ring's hue LUT.
    return effect.hue_pool[index];
  }

  template <int W, int H>
  static Color4 current_ring_color(const DisplacementField<W, H> &effect,
                                   float ring_t) {
    return effect.palette.get(math::wrap_t(ring_t + effect.color_spin));
  }

  template <int W, int H>
  static int noise_lut_samples(const DisplacementField<W, H> &effect,
                               float sin_theta) {
    return noise_lut_samples(
        effect, effect.params.scale1 + effect.params.scale2, sin_theta);
  }

  /** @brief Knots a noise-displaced ring bakes at a feature scale. */
  template <int W, int H>
  static int noise_lut_samples(const DisplacementField<W, H> &,
                               float feature_scale, float sin_theta) {
    using Effect = DisplacementField<W, H>;
    const int lut_n = hs::clamp(
        static_cast<int>(ceilf(Effect::LUT_SAMPLES_PER_UNIT * 2.0f *
                               math::PI_F * feature_scale * sin_theta)),
        Effect::LUT_MIN_SAMPLES, W);
    return (lut_n + Effect::OCTAVE_GRID - 1) / Effect::OCTAVE_GRID *
           Effect::OCTAVE_GRID;
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
    uint64_t valid[2] = {0, 0};
    for (int i = 0; i <= SAMPLE_COUNT; ++i) {
      const float amount = table_case.max_amount * i / SAMPLE_COUNT;
      Pixel actual = DisplacementFieldWhiteBox::sample_hue_table_cached(
          effect, amount, table_case.domain, table_case.cyclic, base, valid);
      HS_EXPECT_EQ(actual.r, expected[i].r);
      HS_EXPECT_EQ(actual.g, expected[i].g);
      HS_EXPECT_EQ(actual.b, expected[i].b);
    }
    HS_EXPECT_EQ(valid[0], ~uint64_t{0});
    HS_EXPECT_TRUE(valid[1] & uint64_t{1});
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
 * @details Every bound below is the measured worst case over the sweep with
 * headroom: peak deltaE is 0.0015 default / 0.0053 cyclic against 0.002 / 0.006
 * bounds, and the paired peaks in encoded space are 8 / 19 sRGB8 codes against
 * 10 / 21. The sRGB8 pair is the looser gate because the encode is non-linear —
 * the same table-interpolation error spans more 8-bit codes where the transfer
 * curve is steep than deltaE weights it — so it is bounded rather than pinned to
 * the perceptual figure.
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
 *          fills the rest with a Catmull-Rom spline. Knots on both octave
 *          grids reproduce the exact product field up to the recurrence's
 *          knot-position drift, and every knot stays within half a canvas
 *          column of it: the spline error is well inside the linear polyline's
 *          own deviation from the field between knots.
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
  const int lut_n = DisplacementFieldWhiteBox::noise_lut_samples(effect, 1.0f);
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
 * @details Identical seeds and mock clock per run, so a divergence isolates
 *          the clip-only paths (the per-ring cap cull and the azimuth-chunk
 *          bake cull) dropping a reachable fragment or sampling a stale LUT
 *          entry. The test explicitly enters the ball phase, whose
 *          footprints drive both culls. Both culls widen with the ring band,
 *          so each quadrant runs at the registered defaults and again at the
 *          full ring pool with maximum thickness and displacement amplitudes.
 */
inline void test_displacement_field_clip_tiles_full() {
  struct Quad {
    int x0, x1, y0, y1;
  };
  const Quad quads[] = {{0, DEFAULT_W / 2, 0, DEFAULT_H / 2},
                        {DEFAULT_W / 2, DEFAULT_W, DEFAULT_H / 2, DEFAULT_H}};
  const int frames = 60;
  // A widest-footprint frame costs ~7x a default one; the shorter window still
  // sits inside the same displacement phase.
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
 * @brief Verifies ShapeShifter's slider contract and preset-row invariants.
 * @details A preset row is authored artistic data that is retuned freely, so its
 *          magnitudes (count, sides, amplitude, speed) are not pinned: a golden
 *          copy of them reds on every intentional retune and reports nothing. The
 *          invariants a retune must not break are pinned instead — every value
 *          inside the range register_param bound it to, the structural selections
 *          each row draws, the non-preset Alpha slider surviving a selection, the
 *          shape/falloff pairing, and every row landing on a distinct parameter
 *          vector.
 */
inline void test_shapeshifter_preset_defaults() {
  reset_effect_globals();
  ShapeShifter<DEFAULT_W, DEFAULT_H> ss;
  ss.init();

  auto value = [&](const char *name) {
    for (const auto &def : ss.getParameters())
      if (std::strcmp(def.name, name) == 0)
        return def.get();
    HS_EXPECT(false, "ShapeShifter parameter is missing");
    return -1.0f;
  };

  // register_param traps on a default outside its range, but a preset row is
  // assigned straight into params and never passes through it.
  auto expect_in_range = [&](const char *label) {
    HS_CONTEXT(label);
    for (const auto &def : ss.getParameters()) {
      HS_CONTEXT(def.name);
      const float v = def.get();
      HS_EXPECT_TRUE(std::isfinite(v));
      HS_EXPECT_GE(v, def.min);
      HS_EXPECT_LE(v, def.max);
      if (def.option_count > 0) {
        HS_EXPECT_EQ(v, std::floor(v));
        HS_EXPECT_LT(v, static_cast<float>(def.option_count));
      }
    }
  };

  expect_in_range("boot state");
  HS_EXPECT_EQ(value("Alpha"), 1.0f); // boots fully opaque
  HS_EXPECT_EQ(value("Shape"), 3.0f);
  HS_EXPECT_EQ(value("Spacing"), 1.0f);
  HS_EXPECT_EQ(value("Function"),
               static_cast<float>(
                   ShapeShifter<DEFAULT_W, DEFAULT_H>::PhaseFunction::SINE));
  HS_EXPECT_EQ(value("Opposite"), 0.0f);
  HS_EXPECT_EQ(value("Alpha Falloff"), 1.0f);

  const char *expected_export_order[] = {
      "Shape", "Count",    "Sides",         "Function", "Amplitude",
      "Speed", "Opposite", "Alpha Falloff", "Spacing"};
  size_t export_index = 0;
  for (const auto &def : ss.getParameters()) {
    if (!def.preset)
      continue;
    HS_EXPECT(export_index < std::size(expected_export_order),
              "ShapeShifter exports an unexpected parameter");
    if (export_index < std::size(expected_export_order))
      HS_EXPECT_EQ(std::string_view(def.name),
                   std::string_view(expected_export_order[export_index]));
    ++export_index;
  }
  HS_EXPECT_EQ(export_index, std::size(expected_export_order));

  const auto *alpha = ss.getParameters().find("Alpha");
  const auto *shape = ss.getParameters().find("Shape");
  const auto *falloff = ss.getParameters().find("Alpha Falloff");
  const auto *spacing = ss.getParameters().find("Spacing");
  const auto *count = ss.getParameters().find("Count");
  const auto *speed = ss.getParameters().find("Speed");
  HS_EXPECT(alpha != nullptr, "ShapeShifter Alpha parameter is missing");
  HS_EXPECT(shape != nullptr, "ShapeShifter Shape parameter is missing");
  HS_EXPECT(falloff != nullptr,
            "ShapeShifter Alpha Falloff parameter is missing");
  HS_EXPECT(spacing != nullptr, "ShapeShifter Spacing parameter is missing");
  HS_EXPECT(count != nullptr, "ShapeShifter Count parameter is missing");
  HS_EXPECT(speed != nullptr, "ShapeShifter Speed parameter is missing");
  if (count)
    HS_EXPECT_EQ(count->max, 288.0f);
  if (speed)
    HS_EXPECT_EQ(speed->max, 0.16f);
  if (alpha)
    HS_EXPECT_FALSE(alpha->preset);
  if (shape) {
    HS_EXPECT_TRUE(shape->is_enum());
    HS_EXPECT_EQ(std::string_view(shape->export_options[3]),
                 std::string_view("ShapeType::PLANAR_STAR"));
    HS_EXPECT_EQ(std::string_view(shape->export_options[4]),
                 std::string_view("ShapeType::SPHERICAL_STAR"));
  }
  if (falloff) {
    HS_EXPECT_TRUE(falloff->is_enum());
    HS_EXPECT_EQ(std::string_view(falloff->export_options[1]),
                 std::string_view("AlphaFalloff::TOWARD_EQUATOR"));
  }
  if (spacing) {
    HS_EXPECT_TRUE(spacing->is_enum());
    HS_EXPECT_EQ(std::string_view(spacing->export_options[1]),
                 std::string_view("RadiusSpacing::SCREEN_BALANCED"));
  }

  HS_EXPECT_EQ(ss.updateParameter("Alpha", 0.37f), ParamSetResult::APPLIED);

  // Structural selections: which primitive and falloff each row draws. These
  // pin the index -> row mapping the profile reports are keyed by; the
  // magnitudes each row sets are deliberately left free.
  const float expected_shapes[] = {3.0f, 1.0f, 3.0f, 2.0f, 3.0f,
                                   1.0f, 1.0f, 1.0f, 2.0f};
  const float expected_falloffs[] = {1.0f, 0.0f, 1.0f, 0.0f, 1.0f,
                                     0.0f, 0.0f, 0.0f, 0.0f};
  HS_EXPECT_EQ(std::size(expected_shapes), ss.getPresetCount());
  HS_EXPECT_EQ(std::size(expected_falloffs), ss.getPresetCount());
  std::vector<std::vector<float>> rows;
  for (size_t i = 0; i < std::size(expected_shapes); ++i) {
    HS_CONTEXT("preset", static_cast<int>(i));
    ss.profile_select_preset(i);
    expect_in_range("preset row");
    HS_EXPECT_EQ(value("Alpha"), 0.37f); // a preset never writes a non-preset
    HS_EXPECT_TRUE(ss.animations_paused());
    HS_EXPECT_EQ(value("Shape"), expected_shapes[i]);
    HS_EXPECT_EQ(value("Alpha Falloff"), expected_falloffs[i]);
    HS_EXPECT_EQ(value("Shape") == 3.0f, value("Alpha Falloff") == 1.0f);

    std::vector<float> row;
    for (const auto &def : ss.getParameters())
      if (def.preset)
        row.push_back(def.get());
    rows.push_back(row);
  }

  // Two rows that collapse onto the same parameter vector are one preset the
  // cycle visits twice, which no range or selection check above would show.
  for (size_t i = 0; i < rows.size(); ++i)
    for (size_t j = i + 1; j < rows.size(); ++j) {
      HS_CONTEXT("preset pair", static_cast<int>(i), static_cast<int>(j));
      HS_EXPECT(rows[i] != rows[j],
                "each ShapeShifter preset must be distinct");
    }
}

/**
 * @brief Renders every Shape and Function slider selection.
 * @details Each primitive is exercised at radii on both sides of the antipode
 * fold while the four phase functions advance through the same Plot pipeline.
 * Every selection renders one frame from an identical fresh state with the
 * preset timer paused, so the only input that moves between renders is the
 * slider: two selections folding to the same frame means the selection switch
 * did not dispatch on them.
 */
inline void test_shapeshifter_slider_selections_render() {
  using SS = ShapeShifter<SMALL_W, SMALL_H>;

  auto render = [](const char *slider, int selection) {
    reset_effect_globals();
    SS ss;
    ss.init();
    ss.setAnimationsPaused(true);
    HS_EXPECT_EQ(ss.updateParameter(slider, static_cast<float>(selection)),
                 ParamSetResult::APPLIED);
    ss.draw_frame();
    ss.advance_display();

    const uint64_t acc = frame_energy<SMALL_W, SMALL_H>(ss);
    uint64_t fold = hs_test::FNV1A64_BASIS;
    for (int y = 0; y < SMALL_H; ++y)
      for (int x = 0; x < SMALL_W; ++x) {
        const Pixel &pixel = ss.get_pixel(x, y);
        for (uint16_t channel : {pixel.r, pixel.g, pixel.b})
          fold = hs_test::fnv1a64_channel(fold, channel);
      }
    HS_EXPECT_GT(acc, 0u);
    return fold;
  };

  auto sweep = [&](const char *slider, int selections) {
    std::vector<uint64_t> folds;
    for (int selection = 0; selection < selections; ++selection)
      folds.push_back(render(slider, selection));
    for (int i = 0; i < selections; ++i)
      for (int j = i + 1; j < selections; ++j) {
        if (folds[i] == folds[j])
          std::printf("  SHAPESHIFTER %s selections %d and %d render the same "
                      "frame (fold %llu)\n",
                      slider, i, j, static_cast<unsigned long long>(folds[i]));
        HS_EXPECT(folds[i] != folds[j],
                  "each slider selection must render a distinct frame");
      }
  };

  sweep("Shape", SS::NUM_SHAPES);
  sweep("Function", SS::NUM_FUNCTIONS);
}

/**
 * @brief Checks the full-frame query on four constructed effects.
 * @details MeshFeedback and Dynamo fold any_crosses_segments true; the two
 *          sampled non-crossing effects keep segment clipping. The WASM
 *          setClip bridge reads this non-virtual Effect accessor.
 */
inline void test_needs_full_frame_gate() {
  // Each effect aliases the same static double buffer (single-live guard) and
  // reconfigures the shared arenas/timeline in init(), so construct one at a
  // time with the same fresh setup smoke_one uses.
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
 * @brief White-box accessor for Voronoi's seeded sites and adaptive block floor.
 */
struct VoronoiWhiteBox {
  using VO = Voronoi<DEFAULT_W, DEFAULT_H>;
  static constexpr int MAX_SITES = VO::MAX_SITES;

  /** @brief Adaptive block floor at render resolution W x H. */
  template <int W, int H> static constexpr int coherence_block_min() {
    return Voronoi<W, H>::COHERENCE_BLOCK_MIN;
  }

  /** @brief Number of currently seeded sites. */
  template <int W, int H> static size_t site_count(const Voronoi<W, H> &v) {
    return v.sites_buffer.size();
  }
  /** @brief Spin axis of seeded site @p i. */
  template <int W, int H>
  static math::Vector site_axis(const Voronoi<W, H> &v, size_t i) {
    return v.sites_buffer[i].axis;
  }

  /** @brief Position of seeded site @p i. */
  template <int W, int H>
  static math::Vector site_position(const Voronoi<W, H> &v, size_t i) {
    return v.sites_buffer[i].pos;
  }

  /** @brief Installs deterministic sites and index-coded colors. */
  template <int W, int H>
  static void set_sites(Voronoi<W, H> &v, std::span<const math::Vector> sites) {
    v.sites_buffer.clear();
    for (size_t i = 0; i < sites.size(); ++i) {
      v.sites_buffer.push_back(
          {sites[i], math::Vector(1, 0, 0),
           Color4(Pixel(static_cast<uint16_t>(i + 1), 0, 0))});
    }
    v.current_num_sites = static_cast<int>(sites.size());
    v.params.num_sites = static_cast<float>(sites.size());
    v.params.speed = 0.0f;
    v.params.sharpness = 0.0f;
    v.params.border_thickness = 0.0f;
  }
};

/**
 * @brief Verifies Voronoi seeds every spin axis through random_vector().
 */
inline void test_voronoi_axes_use_uniform_sampler() {
  using WB = VoronoiWhiteBox;
  reset_effect_globals();
  Voronoi<SMALL_W, SMALL_H> effect;
  effect.init();

  const size_t sites = WB::site_count(effect);
  HS_EXPECT_GT(sites, 0u);
  hs::random().seed(1337u);
  for (size_t i = 0; i < sites; ++i)
    HS_EXPECT_VEC(WB::site_axis(effect, i), math::random_vector(), 0.0f);
}

/**
 * @brief Renders production Voronoi and compares each pixel with exact nearest.
 * @tparam W,H Render resolution.
 * @param sites Site positions on the unit sphere.
 * @param max_deficit Out: worst dot(p, true nearest) - dot(p, union nearest)
 *        over the mismatched pixels (0 when every pixel matches).
 * @return Fraction of pixels whose union-of-corner-pairs nearest matches the
 *         true nearest site.
 */
template <int W, int H>
inline double voronoi_render_nearest_match(std::span<const math::Vector> sites,
                                           float &max_deficit) {
  using WB = VoronoiWhiteBox;
  reset_effect_globals();
  Voronoi<W, H> effect;
  effect.init();
  WB::set_sites(effect, sites);
  effect.draw_frame();
  effect.advance_display();

  long matched = 0;
  max_deficit = 0.0f;
  for (int y = 0; y < H; ++y) {
    for (int x = 0; x < W; ++x) {
      const math::Vector p = math::pixel_to_vector<W, H>(x, y);
      float exact = -2.0f;
      for (size_t i = 0; i < WB::site_count(effect); ++i)
        exact = std::max(exact, math::dot(p, WB::site_position(effect, i)));

      const Pixel rendered = effect.get_pixel(x, y);
      const size_t rendered_site = rendered.r - 1u;
      const bool encoded = rendered.r > 0 && rendered.g == 0 &&
                           rendered.b == 0 &&
                           rendered_site < WB::site_count(effect);
      const float deficit =
          encoded
              ? exact - math::dot(p, WB::site_position(effect, rendered_site))
              : 4.0f;
      if (deficit <= 0.0f)
        ++matched;
      else
        max_deficit = std::max(max_deficit, deficit);
    }
  }
  return static_cast<double>(matched) / (static_cast<double>(W) * H);
}

/**
 * @brief Pins Voronoi's block-candidate-union coverage across the adaptive
 *        block regime at both render resolutions: the union of a block's four
 *        corner pairs contains the true nearest site at every pixel in the
 *        low-density octahedral case, and at >= 99.9% of pixels (with a
 *        sub-visibility dot deficit on the rest) at the MAX_SITES Fibonacci
 *        spread.
 * @details The dense cases floor the adaptive block at COHERENCE_BLOCK_MIN. A
 *          block edge that outruns the cell pixel extent straddles whole cells
 *          and collapses the match fraction, so the short-canvas resolution —
 *          where MAX_SITES cells are sub-pixel and the floor drops to 1 — is
 *          checked for exact coverage.
 */
inline void test_voronoi_union_candidates_cover_nearest() {
  static_assert(VoronoiWhiteBox::coherence_block_min<SMALL_W, SMALL_H>() == 1);
  static_assert(VoronoiWhiteBox::coherence_block_min<DEFAULT_W, DEFAULT_H>() ==
                4);
  float deficit = 0.0f;

  const math::Vector octahedral[] = {
      math::Vector(1, 0, 0),  math::Vector(-1, 0, 0), math::Vector(0, 1, 0),
      math::Vector(0, -1, 0), math::Vector(0, 0, 1),  math::Vector(0, 0, -1),
  };
  const size_t octa_count = sizeof(octahedral) / sizeof(octahedral[0]);
  const std::span<const math::Vector> sparse(octahedral, octa_count);
  const double octa_match =
      voronoi_render_nearest_match<DEFAULT_W, DEFAULT_H>(sparse, deficit);
  HS_EXPECT_EQ(octa_match, 1.0);
  const double octa_match_dev =
      voronoi_render_nearest_match<SMALL_W, SMALL_H>(sparse, deficit);
  HS_EXPECT_EQ(octa_match_dev, 1.0);

  // Dense regime: seed MAX_SITES on a Fibonacci sphere exactly as
  // Voronoi::seed_sites places them, so the adaptive block floors at
  // COHERENCE_BLOCK_MIN.
  constexpr int N = VoronoiWhiteBox::MAX_SITES;
  static math::Vector fib[N];
  for (int i = 0; i < N; ++i)
    fib[i] = math::fib_spiral(N, /*eps=*/0.5f, i);
  const std::span<const math::Vector> dense(fib, N);
  const double fib_match =
      voronoi_render_nearest_match<DEFAULT_W, DEFAULT_H>(dense, deficit);
  HS_EXPECT_GE(fib_match, 0.999);
  HS_EXPECT_LE(deficit, 0.005f);

  const double fib_match_dev =
      voronoi_render_nearest_match<SMALL_W, SMALL_H>(dense, deficit);
  HS_EXPECT_EQ(fib_match_dev, 1.0);
  HS_EXPECT_EQ(deficit, 0.0f);
}

/**
 * @brief Requires a Voronoi segment band to shade every pixel exactly as the
 *        full-canvas render does.
 * @details The coarse-coherence grid decides per block which sites reach a
 *          pixel's candidate union, so a grid whose phase follows the clip
 *          origin shades the same pixel differently depending on which band
 *          renders it — a discontinuity pinned to the segment seam, since
 *          Voronoi is neither full-frame nor persisting. The bands below start
 *          off a block boundary (the adaptive block is 6 px at the default site
 *          count), which an aligned split would hide.
 */
inline void test_voronoi_segment_render_matches_full_frame() {
  constexpr int W = DEFAULT_W;
  constexpr int H = DEFAULT_H;

  auto render = [](int x0, int x1, int y0, int y1,
                   const std::vector<Pixel> *reference = nullptr) {
    reset_effect_globals();
    Voronoi<W, H> effect;
    effect.init();
    effect.set_margin(3);
    effect.set_clip(y0, y1, x0, x1);
    effect.draw_frame();
    effect.advance_display();
    if (reference)
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
          if (effect.clip().contains_y(y) && effect.clip().contains_x(x))
            HS_EXPECT_EQ(effect.get_pixel(x, y),
                         (*reference)[static_cast<size_t>(y) * W + x]);
    std::vector<Pixel> band;
    band.reserve(static_cast<size_t>(x1 - x0) * (y1 - y0));
    for (int y = y0; y < y1; ++y)
      for (int x = x0; x < x1; ++x)
        band.push_back(effect.get_pixel(x, y));
    return band;
  };

  const std::vector<Pixel> full = render(0, W, 0, H);

  struct Band {
    int x0, x1, y0, y1;
  };
  const Band bands[] = {
      {0, 100, 0, H}, {100, W, 0, H}, {0, W, 50, H}, {37, 205, 11, 93}};

  for (const Band &b : bands) {
    HS_CONTEXT("band", b.x0, b.y0);
    size_t lit = 0;
    const std::vector<Pixel> banded = render(b.x0, b.x1, b.y0, b.y1, &full);
    HS_EXPECT_EQ(banded.size(),
                 static_cast<size_t>(b.x1 - b.x0) * (b.y1 - b.y0));
    size_t different = 0;
    size_t i = 0;
    for (int y = b.y0; y < b.y1; ++y)
      for (int x = b.x0; x < b.x1; ++x, ++i) {
        const Pixel &reference = full[static_cast<size_t>(y) * W + x];
        if (banded[i] != reference)
          ++different;
        if (reference.r | reference.g | reference.b)
          ++lit;
      }
    if (different)
      std::printf("  VORONOI SEAM band x[%d,%d) y[%d,%d): %zu of %zu pixels "
                  "differ from the full-canvas render\n",
                  b.x0, b.x1, b.y0, b.y1, different, banded.size());
    HS_EXPECT_EQ(different, static_cast<size_t>(0));
    HS_EXPECT_GT(lit, size_t{0});
  }
}

/**
 * @brief Bounds whole-solid generation, classification and rendering scratch
 * against HankinSolids' exported budgets at the device height.
 * @details The graph-walk soak separately exercises the shipping OpLeg path.
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
    // FaceScratchBuffer on top (the scratch_a-binding path documented on SCRATCH_A_BYTES).
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
    // reset — the same Persist discipline finish_morph_cycle uses to compact
    // between legs.
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

/**
 * @brief White-box accessor for IslamicStars' private build-chain state
 *        (befriended in effects/IslamicStars.h).
 * @details init() opens on recipe entry 0 (dodecahedron_hk62_ambo_hk62).
 * Its build starts after the 16-frame seed fade-in. The probe pre-sets Trans
 * Speed to shorten the build and reads build_active/solid_idx to pin completion.
 * Build bookkeeping is resolution-independent.
 */
struct IslamicBuildProbe {
  using IS = IslamicStars<SMALL_W, SMALL_H>;
  static_assert(IS::MACRO_TRUNCATE_T == RECONCILE_TRUNCATE_T);
  static void check_build_budget(IS &e, size_t budget) {
    e.device_persistent_budget = budget;
    e.check_build_budget();
  }
  static void invalid_bridge_continuation(IS &e) {
    e.schedule_dual_bridge(IS::BuildContinuation::DUAL_DONE);
  }

  template <int W, int H>
  static void set_trans_speed(IslamicStars<W, H> &e, float v) {
    e.params.trans_speed = v;
  }
  template <int W, int H>
  static bool build_active(const IslamicStars<W, H> &e) {
    return e.build_active;
  }
  template <int W, int H> static int solid_idx(const IslamicStars<W, H> &e) {
    return e.solid_idx;
  }
  template <int W, int H> static int dual_bridges(const IslamicStars<W, H> &e) {
    return e.dual_bridges_built;
  }
  template <int W, int H>
  static size_t persistent_budget(const IslamicStars<W, H> &e) {
    return e.device_persistent_budget;
  }
  template <int W, int H> static int front_slot(IslamicStars<W, H> &e) {
    return e.carousel.front_index();
  }
  template <int W, int H>
  static const uint8_t *slot_palette(const IslamicStars<W, H> &e, int slot) {
    return e.slot_face_palette[slot];
  }
  template <int W, int H>
  static size_t slot_faces(IslamicStars<W, H> &e, int slot) {
    return e.carousel.slot(slot).topology.size();
  }
  static constexpr size_t bridge_scratch_a() {
    return IS::BRIDGE_BUDGET.scratch_a;
  }
  static constexpr size_t bridge_scratch_b() {
    return IS::BRIDGE_BUDGET.scratch_b;
  }
  static constexpr int sprite_fade_frames() { return IS::SPRITE_FADE_FRAMES; }
  /**
   * @brief Spawns @p entry out of band, whatever the registry cycle is on.
   * @details Clears the timeline first: the pending scheduled spawn and the
   * outgoing shape's sprite would otherwise run against the injected shape.
   */
  template <int W, int H>
  static void spawn_entry(IslamicStars<W, H> &e, const Solids::Entry &entry) {
    e.timeline.clear();
    e.spawn_entry(entry);
  }
  template <int W, int H>
  static void set_burst_size(IslamicStars<W, H> &e, int size) {
    e.params.burst_size = size;
  }
  template <int W, int H>
  static int cached_burst_size(const IslamicStars<W, H> &e) {
    return e.burst_size_eff;
  }
  template <int W, int H>
  static int cached_burst_window(const IslamicStars<W, H> &e) {
    return (e.burst_size_eff - 1) * e.ripple_stagger_eff + e.ripple_dur_eff;
  }
  template <int W, int H> static int fire_ripple(IslamicStars<W, H> &e) {
    int count = 0;
    {
      Canvas canvas(e);
      e.ripple(canvas);
      count = e.ripple_gen.active_count();
    }
    e.advance_display();
    return count;
  }
  template <int W, int H>
  static PolyMesh clean_endpoint(IslamicStars<W, H> &e,
                                 const Solids::OpStep &step, Arena &a,
                                 Arena &b) {
    return e.clean_endpoint(step, a, b);
  }
};

/** Whole-solid generator of the needle-ending recipe. */
inline PolyMesh generate_needle_recipe_solid(Arena &a, Arena &b) {
  return Solids::build_recipe(
      TRUNCATED_ICOSAHEDRON_AMBO_RELAX_HK54_NEEDLE_RECIPE, a, b);
}

/** The needle-ending recipe as a spawnable entry. It is not in
 * islamic_registry, so the arena gate below builds it from the recipe constant
 * rather than finding it in the roster. */
inline constexpr Solids::Entry NEEDLE_ENTRY = {
    "truncatedIcosahedron_ambo_relax_hk54_needle", generate_needle_recipe_solid,
    Solids::Category::Complex,
    &TRUNCATED_ICOSAHEDRON_AMBO_RELAX_HK54_NEEDLE_RECIPE};

/**
 * @brief Verifies the first recipe seed uses the Sprite envelope's 16-frame
 *        fade-in and starts its build at the full-opacity boundary.
 */
inline void test_islamicstars_seed_sprite_fade_in() {
  reset_effect_globals();
  IslamicBuildProbe::IS effect;
  effect.init();

  constexpr int FADE_FRAMES = IslamicBuildProbe::sprite_fade_frames();
  HS_EXPECT_EQ(FADE_FRAMES, 16);
  for (int frame = 1; frame < FADE_FRAMES; ++frame) {
    effect.draw_frame();
    effect.advance_display();
    HS_EXPECT_FALSE(IslamicBuildProbe::build_active(effect));
  }

  effect.draw_frame();
  effect.advance_display();
  HS_EXPECT_TRUE(IslamicBuildProbe::build_active(effect));
}

inline void test_islamicstars_burst_size_is_snapshotted_per_spawn() {
  reset_effect_globals();
  IslamicBuildProbe::IS effect;
  effect.init();
  const auto entry = Solids::Collections::get_islamic_solids().front();

  IslamicBuildProbe::set_burst_size(effect, 1);
  IslamicBuildProbe::spawn_entry(effect, entry);
  const int single_window = IslamicBuildProbe::cached_burst_window(effect);

  IslamicBuildProbe::set_burst_size(effect, 4);
  HS_EXPECT_EQ(IslamicBuildProbe::cached_burst_size(effect), 1);
  HS_EXPECT_EQ(IslamicBuildProbe::cached_burst_window(effect), single_window);
  HS_EXPECT_EQ(IslamicBuildProbe::fire_ripple(effect), 1);

  IslamicBuildProbe::spawn_entry(effect, entry);
  HS_EXPECT_EQ(IslamicBuildProbe::cached_burst_size(effect), 4);
  HS_EXPECT_TRUE(IslamicBuildProbe::cached_burst_window(effect) >
                 single_window);
  HS_EXPECT_EQ(IslamicBuildProbe::fire_ripple(effect), 4);
}

/** @brief Small recipe builds finish every bridge and retain landed palettes. */
inline void test_islamicstars_smooth_recipe_completion() {
  struct Case {
    Solids::Op op;
    size_t faces;
    int bridges;
  };
  for (const Case &c :
       {Case{Solids::Op::DUAL, 8, 1}, Case{Solids::Op::NEEDLE, 24, 1},
        Case{Solids::Op::KIS, 24, 2}}) {
    reset_effect_globals();
    const Solids::OpStep steps[] = {{c.op}};
    const Solids::Recipe recipe = Solids::make_recipe(1, steps);
    HS_EXPECT_TRUE(
        std::string_view(Solids::simple_registry[recipe.seed].name) == "cube");
    const Solids::Entry entry = {"cube_bridge", Solids::Platonic::cube,
                                 Solids::Category::Complex, &recipe};
    IslamicBuildProbe::IS effect;
    IslamicBuildProbe::set_trans_speed(effect, 2.0f);
    effect.init();
    IslamicBuildProbe::spawn_entry(effect, entry);
    bool started = false;
    bool finished = false;
    for (int frame = 0; frame < 256; ++frame) {
      effect.draw_frame();
      effect.advance_display();
      HS_EXPECT_LE(persistent_arena.get_offset(),
                   IslamicBuildProbe::persistent_budget(effect));
      if (IslamicBuildProbe::build_active(effect))
        started = true;
      else if (started) {
        finished = true;
        break;
      }
    }
    HS_EXPECT_TRUE(started);
    HS_EXPECT_TRUE(finished);
    HS_EXPECT_EQ(IslamicBuildProbe::dual_bridges(effect), c.bridges);
    const int slot = IslamicBuildProbe::front_slot(effect);
    HS_EXPECT_EQ(IslamicBuildProbe::slot_faces(effect, slot), c.faces);
    const uint8_t *palette = IslamicBuildProbe::slot_palette(effect, slot);
    std::vector<uint8_t> landed(palette, palette + c.faces);
    for (uint8_t value : landed)
      HS_EXPECT_LT(value, MeshPaletteBank::N);
    for (int frame = 0; frame < 4; ++frame) {
      effect.draw_frame();
      effect.advance_display();
      HS_EXPECT_FALSE(IslamicBuildProbe::build_active(effect));
      HS_EXPECT_EQ(IslamicBuildProbe::front_slot(effect), slot);
      HS_EXPECT_TRUE(std::equal(landed.begin(), landed.end(),
                                IslamicBuildProbe::slot_palette(effect, slot)));
    }
    HS_EXPECT_GT((frame_energy<SMALL_W, SMALL_H>(effect)), uint64_t(0));
  }
}

/**
 * @brief Drives IslamicStars across the first registry entry's complete
 *        op-by-op build at max trans speed: the build must activate and
 *        finish without a trap, the built shape's per-face colours must never
 *        change from finish_build through its still/ripple/fade display, and
 *        the shape after it must start cleanly and light pixels.
 */
inline void test_islamicstars_recipe_build_smoke() {
  reset_effect_globals();
  IslamicBuildProbe::IS effect;
  IslamicBuildProbe::set_trans_speed(effect, 8.0f);
  effect.init();

  // The snapshot captures entry 0's first completed recipe build.
  constexpr int MAX_FRAMES = 400;
  int frames = 0;
  int build_frames = 0;
  bool was_building = false;
  int built_slot = -1;
  int snap_solid = -1;
  int changed_after_build = 0;
  int constant_frames = 0;
  std::vector<uint8_t> built_pal;
  while (frames < MAX_FRAMES && IslamicBuildProbe::solid_idx(effect) < 2) {
    effect.draw_frame();
    effect.advance_display();
    ++frames;
    const bool building = IslamicBuildProbe::build_active(effect);
    if (building)
      ++build_frames;
    // Snapshot the built shape's per-face colours the moment its build
    // completes; they must stay byte-identical through its still, ripple, and
    // fade phases (the next spawn retires the shape and may reuse its array).
    if (was_building && !building && built_slot < 0) {
      built_slot = IslamicBuildProbe::front_slot(effect);
      snap_solid = IslamicBuildProbe::solid_idx(effect);
      const uint8_t *pal = IslamicBuildProbe::slot_palette(effect, built_slot);
      built_pal.assign(pal,
                       pal + IslamicBuildProbe::slot_faces(effect, built_slot));
    } else if (built_slot >= 0 &&
               IslamicBuildProbe::solid_idx(effect) == snap_solid) {
      const uint8_t *pal = IslamicBuildProbe::slot_palette(effect, built_slot);
      for (size_t f = 0; f < built_pal.size(); ++f)
        if (pal[f] != built_pal[f])
          ++changed_after_build;
      ++constant_frames;
    }
    was_building = building;
  }
  HS_EXPECT_LT(frames, MAX_FRAMES);
  HS_EXPECT_GT(build_frames, 0);
  HS_EXPECT_TRUE(!IslamicBuildProbe::build_active(effect));
  HS_EXPECT_GE(built_slot, 0);
  HS_EXPECT_GT(built_pal.size(), (size_t)0);
  HS_EXPECT_GT(constant_frames, 0);
  HS_EXPECT_EQ(changed_after_build, 0);

  // The following shape renders lit frames.
  for (int f = 0; f < 12; ++f) {
    effect.draw_frame();
    effect.advance_display();
  }
  const uint64_t acc = frame_energy<SMALL_W, SMALL_H>(effect);
  HS_EXPECT_GT(acc, (uint64_t)0);
}

/**
 * @brief Drives IslamicStars through every registry entry and then through the
 *        needle recipe, pinning the persistent arena against the effect's own
 *        budget.
 * @details The per-chain gate measures a build in isolation, which misses what
 *          the effect actually holds: a build follows whatever shape preceded
 *          it. Cycling the roster exercises builds against their predecessors. The roster runs at 288x144 with Trans Speed 8; the
 *          separate dual-bridge gate covers closing legs that this fast
 *          cadence can omit. An arena overrun traps.
 *          The needle is the heaviest smooth-bridge shape and sets the
 *          scratch_a-heavy split, so it is measured separately.
 */
inline void test_islamicstars_roster_cycle_fits_budget() {
  reset_effect_globals();
  // spawn_entry selects BRIDGE_BUDGET, RECIPE_BUDGET, or GENERATED_BUDGET.
  // Host scratch uses the device caps; persistent usage is checked against the
  // live per-shape device budget. Resplitting rebases scratch high-water marks.
  {
    IslamicStars<288, 144> effect;
    IslamicBuildProbe::set_trans_speed(effect, 8.0f);
    effect.init();

    auto solids = Solids::Collections::get_islamic_solids();
    const int entries = static_cast<int>(solids.size());
    constexpr int MAX_FRAMES = 20000;
    size_t a_peak = 0, b_peak = 0, persist_peak = 0;
    size_t worst_p = 0, worst_p_budget = 1;
    int worst_p_idx = -1; // shape at the worst persistent/budget ratio
    int frames = 0, shapes = 0, builds = 0;
    bool was_building = false;
    int last = IslamicBuildProbe::solid_idx(effect);
    while (frames < MAX_FRAMES && shapes <= entries) {
      effect.draw_frame();
      effect.advance_display();
      ++frames;
      const int cur = IslamicBuildProbe::solid_idx(effect);
      // Palette variety of each finished build: distinct landed palettes on the
      // shape the real leg chain just landed.
      const bool building = IslamicBuildProbe::build_active(effect);
      if (was_building && !building) {
        const int front = IslamicBuildProbe::front_slot(effect);
        const uint8_t *pal = IslamicBuildProbe::slot_palette(effect, front);
        const size_t nf = IslamicBuildProbe::slot_faces(effect, front);
        bool seen[MeshPaletteBank::N] = {};
        int distinct = 0;
        for (size_t f = 0; f < nf; ++f)
          if (pal[f] < MeshPaletteBank::N && !seen[pal[f]]) {
            seen[pal[f]] = true;
            ++distinct;
          }
        std::printf("  [built] %s: %d/%d palettes on %zu faces\n",
                    (cur >= 0 && cur < entries) ? solids[cur].name : "?",
                    distinct, MeshPaletteBank::N, nf);
      }
      was_building = building;
      const size_t p = persistent_arena.get_offset();
      const size_t p_budget = IslamicBuildProbe::persistent_budget(effect);
      a_peak = std::max(a_peak, scratch_arena_a.get_high_water_mark());
      b_peak = std::max(b_peak, scratch_arena_b.get_high_water_mark());
      persist_peak = std::max(persist_peak, p);
      if (p_budget &&
          uint64_t(p) * worst_p_budget > uint64_t(worst_p) * p_budget) {
        worst_p = p;
        worst_p_budget = p_budget;
        worst_p_idx = cur;
      }
      // Per-shape persistent budget (device figure); scratch is trap-enforced.
      HS_EXPECT_LE(p, p_budget);
      if (building)
        ++builds;
      if (cur != last) {
        last = cur;
        ++shapes;
      }
    }

    const char *worst_name = (worst_p_idx >= 0 && worst_p_idx < entries)
                                 ? solids[worst_p_idx].name
                                 : "?";
    std::printf(
        "  [roster] %d shapes over %d frames, %d build frames: scratch_a "
        "peak=%zu B, scratch_b peak=%zu B, persistent peak=%zu B; tightest "
        "persistent %zu/%zu at %s\n",
        shapes, frames, builds, a_peak, b_peak, persist_peak, worst_p,
        worst_p_budget, worst_name);
    HS_EXPECT_GT(shapes, entries - 1);
    HS_EXPECT_GT(builds, 0);
  }

  // The needle build. Trans Speed 2, not the roster's 8: a compressed stage can
  // drop the closing bridge leg before it runs, which is the peak.
  reset_effect_globals();
  size_t na_peak = 0, nb_peak = 0, np_peak = 0;
  int needle_frames = 0;
  bool needle_built = false;
  {
    constexpr int NEEDLE_MAX_FRAMES = 4000;
    IslamicStars<288, 144> effect;
    IslamicBuildProbe::set_trans_speed(effect, 2.0f);
    effect.init();
    IslamicBuildProbe::spawn_entry(effect, NEEDLE_ENTRY);
    bool was_building = false;
    while (needle_frames < NEEDLE_MAX_FRAMES) {
      effect.draw_frame();
      effect.advance_display();
      ++needle_frames;
      const bool building = IslamicBuildProbe::build_active(effect);
      const size_t p = persistent_arena.get_offset();
      na_peak = std::max(na_peak, scratch_arena_a.get_high_water_mark());
      nb_peak = std::max(nb_peak, scratch_arena_b.get_high_water_mark());
      np_peak = std::max(np_peak, p);
      HS_EXPECT_LE(p, IslamicBuildProbe::persistent_budget(effect));
      if (building)
        was_building = true;
      else if (was_building) {
        needle_built = true;
        break;
      }
    }
  }
  std::printf("  [needle] smooth-path peaks over %d frames: scratch_a=%zu/%zu "
              "scratch_b=%zu/%zu persistent=%zu/%zu\n",
              needle_frames, na_peak, IslamicBuildProbe::bridge_scratch_a(),
              nb_peak, IslamicBuildProbe::bridge_scratch_b(), np_peak,
              DEVICE_GLOBAL_ARENA_SIZE - IslamicBuildProbe::bridge_scratch_a() -
                  IslamicBuildProbe::bridge_scratch_b());
  HS_EXPECT_TRUE(needle_built);
  // needle actually reached its scratch_a-heavy split (proves the smooth path
  // ran, not a silently-dropped build).
  HS_EXPECT_GT(na_peak, 120u * 1024u);
}

/**
 * @brief Drives IslamicStars until TARGET_BRIDGES dual bridges complete,
 *        pinning the scratch peaks against the effect's budget.
 * @details The roster gate at Trans Speed 8 compresses each build so far that a
 *          heavy shape's closing dual leg can be dropped before it runs; this
 *          drives at a modest speed so closing bridge legs can complete. The
 *          bridge's leg 3 rebuilds the medial for
 *          its handoff centroids, whose scratch must not co-reside with the
 *          leg's own arrival mesh -- an over-budget leg traps in the host arena.
 */
inline void test_islamicstars_dual_bridge_fits_budget() {
  reset_effect_globals();
  IslamicStars<288, 144> effect;
  IslamicBuildProbe::set_trans_speed(effect, 2.0f);
  effect.init();

  constexpr int TARGET_BRIDGES = 5;
  constexpr int MAX_FRAMES = 40000;
  size_t a_peak = 0, b_peak = 0, persist_peak = 0;
  int frames = 0;
  // Scratch is hard-capped per-shape (spawn_shape's resplit), so a leg over its
  // split traps here; completing without a trap proves every full bridge fit.
  // Persistent is host-inflated, so it is checked per frame against the effect's
  // live per-shape device budget.
  while (frames < MAX_FRAMES &&
         IslamicBuildProbe::dual_bridges(effect) < TARGET_BRIDGES) {
    effect.draw_frame();
    effect.advance_display();
    ++frames;
    a_peak = std::max(a_peak, scratch_arena_a.get_high_water_mark());
    b_peak = std::max(b_peak, scratch_arena_b.get_high_water_mark());
    const size_t p = persistent_arena.get_offset();
    persist_peak = std::max(persist_peak, p);
    HS_EXPECT_LE(p, IslamicBuildProbe::persistent_budget(effect));
  }
  std::printf(
      "  [dual-bridge] %d bridges over %d frames: scratch_a peak=%zu B, "
      "scratch_b peak=%zu B, persistent peak=%zu B\n",
      IslamicBuildProbe::dual_bridges(effect), frames, a_peak, b_peak,
      persist_peak);
  HS_EXPECT_GE(IslamicBuildProbe::dual_bridges(effect), TARGET_BRIDGES);
}

template <typename EffectT>
inline void check_manual_preset_navigation(size_t expected_count) {
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
  check_manual_preset_navigation<MindSplatter<SMALL_W, SMALL_H>>(8);
  check_manual_preset_navigation<DreamBalls<SMALL_W, SMALL_H>>(10);
  check_manual_preset_navigation<Comets<SMALL_W, SMALL_H>>(12);
  check_manual_preset_navigation<MeshFeedback<SMALL_W, SMALL_H>>(12);
  check_manual_preset_navigation<ShapeShifter<SMALL_W, SMALL_H>>(9);
}
