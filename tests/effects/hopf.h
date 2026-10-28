/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief HopfFibration white-box cases. */

/**
 * @brief White-box accessor for HopfFibration's S3-lift + stereographic
 *        projection.
 * @details Sets the private per-frame cache (tumble sines/cosines, fold/flow/
 *          tumble-y phases) and a fiber's base coordinates, then calls the real
 *          projection.
 */
struct HopfWhiteBox {
  using HF = HopfFibration<DEFAULT_W, DEFAULT_H>;

  /** @brief Renders a growing production trail or its frozen column reference.
   * @details The reference holds the pose of @p fraction, turned
   * @p displacement columns about the vertical and @p displacement rows along
   * the meridian of @p fraction's columns, and draws through the same timed
   * route with no motion: its trail nodes are steady single samples, at the
   * point route's tolerance.
   * @param timed Selects the production temporal draw.
   * @param fraction Frozen frame fraction.
   * @param displacement Reference displacement in rows and columns.
   * @param start_column Start of the retained trail.
   * @return Linear framebuffer luminance. */
  static std::vector<float> render_temporal_trail(bool timed, float fraction,
                                                  float displacement,
                                                  float start_column) {
    reset_effect_globals();
    HF fx;
    fx.init();
    for (size_t i = 0; i < HF::ACTUAL_FIBERS; ++i)
      fx.trails[i].clear();
    const auto source = [&](int index) {
      return math::pixel_to_vector<DEFAULT_W, DEFAULT_H>(
          start_column + index * 3.3f, 65.0f + 8.0f * sinf(index * 0.2f));
    };
    for (int i = 0; i < HF::TRAIL_LEN; ++i)
      fx.trails[0].record(source(i));
    fx.retired[0] = source(-1);
    fx.trail_was_full = true;
    fx.history_ready = true;
    const auto end_pose =
        math::make_rotation(math::Y_AXIS, 2 * math::PI_F / 600);
    if (timed) {
      fx.orientation.set(math::Quaternion());
      fx.orientation.push(end_pose);
    } else {
      int column = 0;
      for (int c = 1; c < DEFAULT_W / 2; ++c)
        if (fabsf(Spacetime::SamplingSurface<DEFAULT_W>::time(c) - fraction) <
            fabsf(Spacetime::SamplingSurface<DEFAULT_W>::time(column) -
                  fraction))
          column = c;
      const math::Vector meridian = math::pixel_to_vector<DEFAULT_W, DEFAULT_H>(
          static_cast<float>(column), DEFAULT_H / 2.0f);
      const math::Vector across =
          math::cross(math::Y_AXIS, meridian).normalized();
      fx.orientation.set(
          math::make_rotation(across,
                              displacement * math::RADIANS_PER_ROW<DEFAULT_H>) *
          math::make_rotation(math::Y_AXIS,
                              displacement *
                                  math::RADIANS_PER_COLUMN<DEFAULT_W>) *
          math::slerp(math::Quaternion(), end_pose, fraction));
    }
    {
      Canvas canvas(fx);
      fx.render_trails(canvas);
    }
    fx.advance_display();
    return capture_luminance<DEFAULT_W, DEFAULT_H>(fx);
  }

  static size_t fiber_count() { return HF::ACTUAL_FIBERS; }
  static math::Vector project(HF &fx, size_t i) { return fx.hopf_project(i); }
  static void set_fiber(HF &fx, size_t i, float azimuth, float polar) {
    fx.fibers[i] = math::Spherical(azimuth, polar);
  }
  static void set_cache(HF &fx, float cx, float sx, float cy, float sy,
                        float fold_base, float flow_rad, float ty_rad) {
    fx.cx = cx;
    fx.sx = sx;
    fx.cy = cy;
    fx.sy = sy;
    fx.fold_base = fold_base;
    fx.flow_rad = flow_rad;
    fx.ty_rad = ty_rad;
  }
  static const math::Spherical &fiber(const HF &fx, size_t i) {
    return fx.fibers[i];
  }
  static float phase_step() { return HF::PHASE_STEP; }
  static void set_folding(HF &fx, float v) { fx.params.folding = v; }
  static void set_twist(HF &fx, float v) { fx.params.twist = v; }
  static size_t trim_start(size_t len, float alpha) {
    return HF::trail_trim_start(len, alpha);
  }
};

/**
 * @brief Pins HopfFibration::hopf_project against closed forms under identity
 *        and active 4D tumble.
 * @details Identity tumble with zero folding and twist reduces fiber 0
 *          (beta == 0) to a plain S3 lift, whose stereographic image is the
 *          normalized (q0, q1, q2). Under active tumble, with and without
 *          fold and twist, every fiber off the S3 pole matches a
 *          double-precision reference.
 */
inline void test_hopf_projection_math() {
  reset_effect_globals();
  using WB = HopfWhiteBox;

  HopfFibration<DEFAULT_W, DEFAULT_H> fx;
  fx.init();

  // Identity tumble + no folding/twist: fiber 0 has beta == 0, so q3 == 0 and the
  // projection collapses to the normalized (q0, q1, q2).
  WB::set_folding(fx, 0.0f);
  WB::set_twist(fx, 0.0f);
  WB::set_cache(fx, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f);

  const float polar = 1.2f, azimuth = 0.7f;
  WB::set_fiber(fx, 0, azimuth, polar);
  const math::Vector v = WB::project(fx, 0);

  const float eta = polar * 0.5f;
  const float ce = math::fast_cosf(eta), se = math::fast_sinf(eta);
  const math::Vector expected = math::Vector(ce * math::fast_cosf(azimuth),
                                             ce * math::fast_sinf(azimuth), se)
                                    .normalized();
  HS_EXPECT_VEC(v, expected, 1e-4f);
  HS_EXPECT_NEAR(v.magnitude(), 1.0f, 1e-4f);

  // Active tumble, then tumble with fold and twist, against a double-precision
  // fold -> twist -> S3 -> R_xw/R_yz -> stereographic reference.
  const double ax = 0.9, ay = 0.4, flow = 0.3;
  const double fold_base = std::sin(ax * 0.5) * 0.5;
  WB::set_cache(fx, static_cast<float>(std::cos(ax)),
                static_cast<float>(std::sin(ax)),
                static_cast<float>(std::cos(ay)),
                static_cast<float>(std::sin(ay)), static_cast<float>(fold_base),
                static_cast<float>(flow), static_cast<float>(ay));
  for (const auto [folding, twist] : {std::pair{0.0f, 0.0f}, {0.5f, 2.0f}}) {
    WB::set_folding(fx, folding);
    WB::set_twist(fx, twist);
    int checked = 0;
    for (size_t i = 0; i < WB::fiber_count(); ++i) {
      const math::Vector p = WB::project(fx, i);
      HS_EXPECT_TRUE(finite_vec(p));
      HS_EXPECT_NEAR(p.magnitude(), 1.0f, 1e-3f);

      const math::Spherical &sph = WB::fiber(fx, i);
      const double azimuth0 = sph.theta;
      const double eta =
          sph.phi / 2.0 +
          std::sin(azimuth0 * 2.0 + ay + fold_base) * 0.2 * folding;
      const double azimuth = azimuth0 + eta * twist;
      const double beta = flow + static_cast<double>(i) * WB::phase_step();
      const double q0 = std::cos(eta) * std::cos(azimuth + beta);
      const double q1 = std::cos(eta) * std::sin(azimuth + beta);
      const double q2 = std::sin(eta) * std::cos(beta);
      const double q3 = std::sin(eta) * std::sin(beta);
      const double r0 = q0 * std::cos(ax) - q3 * std::sin(ax);
      const double r3 = q0 * std::sin(ax) + q3 * std::cos(ax);
      const double r1 = q1 * std::cos(ay) - q2 * std::sin(ay);
      const double r2 = q1 * std::sin(ay) + q2 * std::cos(ay);
      const double factor = 1.0 / (1.0 + 0.001 - r3);
      const double x = r0 * factor, y = r1 * factor, z = r2 * factor;
      const double len = std::sqrt(x * x + y * y + z * z);
      // Near the S3 pole the direction is ill-conditioned in q.
      if (std::sqrt(r0 * r0 + r1 * r1 + r2 * r2) < 0.5)
        continue;
      ++checked;
      // Tolerance absorbs the fast_sinf/fast_cosf approximation.
      HS_EXPECT_VEC(p,
                    math::Vector(static_cast<float>(x / len),
                                 static_cast<float>(y / len),
                                 static_cast<float>(z / len)),
                    5e-3f);
    }
    HS_EXPECT_GT(checked, static_cast<int>(WB::fiber_count() / 2));
  }
}

/**
 * @brief Checks HopfFibration's trail trim across supported trail lengths and
 *        visible alphas.
 * @details For alpha >= MIN_VISIBLE_ALPHA, first stays at or below len - 2. The
 *          trim must match its per-sample cutoff and be monotone in alpha.
 */
inline void test_hopf_trail_trim_keeps_a_segment() {
  using WB = HopfWhiteBox;
  using HF = HopfFibration<DEFAULT_W, DEFAULT_H>;
  // Ascending, all at or above MIN_VISIBLE_ALPHA.
  const float alphas[] = {MIN_VISIBLE_ALPHA, 0.005f, 0.01f, 0.05f, 0.4f, 1.0f};
  for (size_t len = 2; len <= static_cast<size_t>(HF::TRAIL_LEN); ++len) {
    size_t previous_first = len;
    for (float alpha : alphas) {
      const size_t first = WB::trim_start(len, alpha);
      HS_EXPECT_LT(first + 1, len);
      HS_EXPECT_LE(first, previous_first);
      // The kept segment meets the cutoff; its predecessor falls below it.
      HS_EXPECT_GE(static_cast<float>(first + 1) / (len - 1) * alpha,
                   MIN_ENCODABLE_ALPHA);
      if (first > 0)
        HS_EXPECT_LT(static_cast<float>(first) / (len - 1) * alpha,
                     MIN_ENCODABLE_ALPHA);
      previous_first = first;
    }
  }
}

/** @brief Hopf orientation, births and retiring-tail shading match frozen columns. */
inline void test_hopf_temporal_trail_oracle() {
  const ScopedTemporalSampling temporal_sampling(true);
  for (float start : {-15.0f, 125.0f}) {
    auto render = [start](bool timed, float fraction, float displacement) {
      return HopfWhiteBox::render_temporal_trail(timed, fraction, displacement,
                                                 start);
    };
    hs_test::seam_compare_render<DEFAULT_W, DEFAULT_H>(render, "Hopf trail");
  }
}
