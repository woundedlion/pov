/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief GnomonicStars white-box cases. */

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

  template <int W, int H>
  static bool warp_drawn(const GnomonicStars<W, H> &fx) {
    return fx.warp_drawn;
  }
  template <int W, int H>
  static const math::Vector &drawn_centre(const GnomonicStars<W, H> &fx,
                                          int i) {
    return fx.drawn_centres[i];
  }
  template <int W, int H>
  static const math::MobiusParams &warp(const GnomonicStars<W, H> &fx) {
    return fx.transformer.active_params(0);
  }
  template <int W, int H>
  static const math::Vector &spiral(const GnomonicStars<W, H> &fx, int i) {
    return fx.spiral_cache[i];
  }
  template <int W, int H>
  static math::Vector warped(const GnomonicStars<W, H> &fx, int i) {
    return fx.transformer.transform(fx.spiral_cache[i]);
  }
  template <int W, int H>
  static math::Vector star_drift(const GnomonicStars<W, H> &fx,
                                 const math::Vector &previous,
                                 const math::Vector &warped) {
    return fx.star_drift(previous, warped);
  }
  template <int W, int H>
  static const math::Quaternion &orientation(const GnomonicStars<W, H> &fx) {
    return fx.orientation.get();
  }
  template <int W, int H> static void refresh_spiral(GnomonicStars<W, H> &fx) {
    fx.refresh_spiral();
  }
  template <int W, int H> static constexpr float cull_drift() {
    return GnomonicStars<W, H>::CULL_DRIFT;
  }
};

/**
 * @brief Whether two Mobius parameter sets hold the same coefficients.
 * @param a First set.
 * @param b Second set.
 * @return True when every coefficient matches bit for bit.
 */
inline bool same_mobius(const math::MobiusParams &a,
                        const math::MobiusParams &b) {
  return std::memcmp(&a, &b, sizeof(math::MobiusParams)) == 0;
}

/**
 * @brief Each star drifts from where the last frame drew it to its current
 *        warped centre, at the current orientation, however far it moves; a
 *        rebuilt spiral holds it.
 */
inline void test_gnomonicstars_stars_drift_from_last_draw() {
  using WB = GnomonicStarsWhiteBox;
  constexpr int W = SMALL_W, H = SMALL_H;
  reset_effect_globals();
  GnomonicStars<W, H> fx;
  fx.init();
  HS_EXPECT_FALSE(WB::warp_drawn(fx));
  fx.draw_frame();
  fx.advance_display();
  HS_EXPECT_TRUE(WB::warp_drawn(fx));

  const float cull_drift = WB::cull_drift<W, H>();
  int moving = 0;
  for (int f = 0; f < 6; ++f) {
    const math::MobiusParams previous = WB::warp(fx);
    fx.draw_frame();
    fx.advance_display();
    HS_EXPECT_FALSE(same_mobius(previous, WB::warp(fx)));
    const math::Quaternion &q = WB::orientation(fx);
    for (int i = 0; i < 600; i += 37) {
      HS_EXPECT_VEC(WB::drawn_centre(fx, i), WB::warped(fx, i), 0.0f);
      const math::Vector point = WB::spiral(fx, i);
      const math::Vector drawn =
          math::gnomonic_mobius_transform(point, previous);
      const math::Vector before = math::rotate(drawn, q);
      const math::Vector now = math::rotate(WB::warped(fx, i), q);
      const math::Vector drift = WB::star_drift(fx, drawn, WB::warped(fx, i));
      if (math::dot(drift, drift) == 0.0f)
        continue;
      ++moving;
      HS_EXPECT_NEAR(math::dot(Spacetime::rotate_by(before, drift), now), 1.0f,
                     1e-6f);
    }
  }
  HS_EXPECT_GT(moving, 50);

  const math::Vector drawn = WB::drawn_centre(fx, 5);
  const math::Quaternion &q = WB::orientation(fx);
  for (const float angle : {1.5f * cull_drift, 0.6f, 1.7f, 3.0f}) {
    HS_CONTEXT("angle", angle);
    const math::Vector jump = math::rotate(
        drawn, math::make_rotation(
                   math::cross(drawn, math::Y_AXIS).normalized(), angle));
    const math::Vector drift = WB::star_drift(fx, drawn, jump);
    HS_EXPECT_NEAR(sqrtf(math::dot(drift, drift)), angle, 1e-4f);
    const math::Vector before = math::rotate(drawn, q);
    const math::Vector now = math::rotate(jump, q);
    HS_EXPECT_VEC(Spacetime::rotate_by(before, drift), now, 1e-5f);
  }

  HS_EXPECT_TRUE(fx.updateParameter("Points", 300.0f) ==
                 ParamSetResult::APPLIED);
  WB::refresh_spiral(fx);
  HS_EXPECT_FALSE(WB::warp_drawn(fx));
}

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
      HS_EXPECT_VEC(got, want, 0.0f);
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
      if (is_black(fx.get_pixel(x, y)))
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
