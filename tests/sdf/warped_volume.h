/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// WarpedVolume + Warp::Twist (domain warp, Lipschitz bound, normal correction)
// ============================================================================

/** @brief Verifies Twist::apply displaces Y by amplitude·sin(twist·θ), θ=atan2(z,x). */
inline void test_twist_apply_displaces_y() {
  SDF::Warp::Twist tw{/*twist=*/1, /*amplitude=*/0.3f, /*R=*/1.0f};

  // θ = atan2(0, 1) = 0 → no displacement.
  math::Vector a(1.0f, 0.5f, 0.0f);
  math::Vector ra = tw.apply(a, tw.make_ctx(a));
  HS_EXPECT_VEC(ra, math::Vector(1.0f, 0.5f, 0.0f), 1e-3f);

  // θ = π/2 → sin(twist·π/2) = 1 → Y drops by amplitude.
  math::Vector b(0.0f, 0.5f, 1.0f);
  math::Vector rb = tw.apply(b, tw.make_ctx(b));
  HS_EXPECT_VEC(rb, math::Vector(0.0f, 0.5f - 0.3f, 1.0f), 5e-4f);
}

/** @brief Verifies Twist::lipschitz is 1 for twist 0 and matches the closed form otherwise. */
inline void test_twist_lipschitz_identity_and_closed_form() {
  SDF::Warp::Twist flat{0, 0.5f, 1.0f};
  HS_EXPECT_NEAR(flat.lipschitz(math::Vector(2, 0, 0),
                                flat.make_ctx(math::Vector(2, 0, 0))),
                 1.0f, 1e-6f);

  // twist=2, amplitude=0.5 at s=2: γ = 0.5, bound = γ/2 + √(1 + γ²/4).
  SDF::Warp::Twist tw{2, 0.5f, 1.0f};
  math::Vector p(2, 0, 0);
  float s = tw.make_ctx(p);
  HS_EXPECT_NEAR(s, 2.0f, 1e-6f);
  float gamma = 0.5f;
  float expected = 0.5f * gamma + std::sqrt(1.0f + 0.25f * gamma * gamma);
  HS_EXPECT_NEAR(tw.lipschitz(p, s), expected, 1e-5f);
}

/** @brief Verifies Twist::bounding_inflation returns the displacement amplitude. */
inline void test_twist_bounding_inflation() {
  SDF::Warp::Twist tw{3, 0.42f, 1.0f};
  HS_EXPECT_NEAR(tw.bounding_inflation(), 0.42f, 1e-6f);
}

/** @brief Compares twist kernels with the unspecialized harmonic recurrence. */
inline void test_twisted_torus_matches_recurrence() {
  hs::Pcg32 rng(0x73145u);
  float worst_distance = 0.0f, worst_normal = 0.0f;
  for (int twist = 0; twist <= 8; ++twist) {
    for (float scale : {0.001f, 0.03f, 0.3f, 1.0f, 4.0f}) {
      for (float precision : {0.0f, 0.014f * scale, 0.14f * scale}) {
        const SDF::WarpedVolume<SDF::Torus, SDF::Warp::Twist> torus{
            {0.45f * scale, 0.14f * scale},
            {twist, 0.35f * scale, 0.45f * scale},
            precision};
        const auto &base = torus.base;
        const auto &warp = torus.warp;
        for (int i = 0; i < 512; ++i) {
          const float theta = rand_uniform(rng, -math::PI_F, math::PI_F);
          const float tube_angle = rand_uniform(rng, -math::PI_F, math::PI_F);
          float radius = base.R + rand_uniform(rng, -2.0f, 2.0f) * base.r;
          if (i < 4)
            radius = static_cast<float>(i) * 0.5f * math::TOLERANCE;
          const math::Vector p(radius * cosf(theta),
                               warp.amplitude * sinf(twist * theta) +
                                   base.r * sinf(tube_angle),
                               radius * sinf(theta));
          const float s = sqrtf(p.x * p.x + p.z * p.z);
          const float inv_s = s > math::TOLERANCE ? 1.0f / s : 0.0f;
          float sin_prev = 0.0f, sin_n = p.z * inv_s;
          float cos_prev = 1.0f, cos_n = p.x * inv_s;
          const float two_cos = 2.0f * p.x * inv_s;
          for (int k = 1; k < twist; ++k) {
            const float sin_next = two_cos * sin_n - sin_prev;
            sin_prev = sin_n;
            sin_n = sin_next;
            const float cos_next = two_cos * cos_n - cos_prev;
            cos_prev = cos_n;
            cos_n = cos_next;
          }
          if (twist == 0 || s <= math::TOLERANCE) {
            sin_n = 0.0f;
            cos_n = 1.0f;
          }
          const float gate = precision > 0.0f ? precision : warp.amplitude;
          const float q = s - base.R;
          const float dy = std::max(fabsf(p.y) - warp.amplitude, 0.0f);
          const float qq = q * q + dy * dy;
          const float threshold = gate + base.r;
          const math::Vector warped(p.x, p.y - warp.amplitude * sin_n, p.z);
          float expected = base.distance(warped);
          if (qq > threshold * threshold)
            expected = sqrtf(qq) - base.r;
          else if (expected > 0.0f)
            expected *=
                warp.lipschitz_inv(inv_s == 0.0f ? warp.two_over_r : inv_s);
          worst_distance = fold_worst(
              worst_distance, fabsf(torus.distance(p) - expected) / scale);
          if (s > math::TOLERANCE) {
            const math::Vector expected_normal = warp.correct_normal_inv(
                p, base.normal_raw(warped, inv_s), inv_s, cos_n);
            const math::Vector difference = torus.normal(p) - expected_normal;
            worst_normal = fold_worst(worst_normal, difference.length());
          }
        }
      }
    }
  }
  HS_EXPECT_LT(worst_distance, 2e-6f);
  HS_EXPECT_LT(worst_normal, 2e-6f);
}

/**
 * @brief Verifies the recurrence siblings agree at the degenerate-axis
 *        threshold itself: an XZ radius of exactly TOLERANCE is the axis for
 *        sin_ntheta/cos_ntheta and sincos_ntheta alike, and one ulp above it
 *        is not.
 */
inline void test_twist_axis_threshold_siblings_agree() {
  const float above = std::nextafter(math::TOLERANCE, 1.0f);
  for (int twist = 0; twist <= 8; ++twist) {
    const SDF::Warp::Twist warp{twist, 0.35f, 0.45f};
    for (float s : {math::TOLERANCE, above}) {
      for (const math::Vector &p :
           {math::Vector(s, 0.1f, 0.0f), math::Vector(0.0f, 0.1f, s),
            math::Vector(-0.6f * s, -0.1f, 0.8f * s)}) {
        const auto both = warp.sincos_ntheta(p, s);
        const auto sin_inv = warp.sin_ntheta_inv(p, s);
        HS_EXPECT_EQ(warp.sin_ntheta(p, s), both.sin_n);
        HS_EXPECT_EQ(warp.cos_ntheta(p, s), both.cos_n);
        HS_EXPECT_EQ(sin_inv.sin_n, both.sin_n);
        if (twist > 0 && s > math::TOLERANCE)
          HS_EXPECT_EQ(sin_inv.lipschitz_arg, 1.0f / s);
        else
          HS_EXPECT_EQ(sin_inv.lipschitz_arg, warp.two_over_r);
      }
    }
  }
}

constexpr double PI_DBL = 3.14159265358979323846;

/**
 * @brief Samples distance from a point to a twisted torus surface.
 * @param p Query point.
 * @param R Major radius.
 * @param r Minor radius.
 * @param n Twist count.
 * @param A Twist amplitude.
 * @param steps Number of theta samples around the surface.
 * @return Minimum unsigned distance over the sampled tube circles.
 * @details The surface is the union over theta of tube circles of radius r
 * centered at (R cos t, A sin(n t), R sin t) in the plane spanned by the radial
 * direction at t and the y axis, so the tube angle is solved in closed form and
 * only theta is sampled.
 */
inline double twisted_torus_distance(const math::Vector &p, double R, double r,
                                     int n, double A, int steps) {
  double best = 1e30;
  for (int i = 0; i < steps; ++i) {
    double t = 2.0 * PI_DBL * i / steps;
    double ct = std::cos(t), st = std::sin(t);
    double dx = p.x - R * ct, dy = p.y - A * std::sin(n * t), dz = p.z - R * st;
    double u = dx * ct + dz * st;
    double w = -dx * st + dz * ct;
    double rad = std::sqrt(u * u + dy * dy) - r;
    best = std::min(best, std::sqrt(rad * rad + w * w));
  }
  return best;
}

/** @brief Pins march distances against an independently sampled surface. */
inline void test_warped_volume_distance_is_sphere_trace_safe() {
  struct Case {
    double R, r, A;
    int n;
  };
  const Case cases[] = {{1.0, 0.3, 0.2, 3},    {1.0, 0.3, 0.0, 3},
                        {1.0, 0.3, 0.2, 0},    {1.0, 0.31, 2.5, 5},
                        {0.45, 0.14, 0.35, 2}, {0.45, 0.14, 0.35, 8},
                        {2.0, 0.05, 0.9, 7}};
  hs::Pcg32 rng(0x51afeu);
  int correction_needed = 0;
  for (const Case &c : cases) {
    SDF::WarpedVolume<SDF::Torus, SDF::Warp::Twist> volume{
        SDF::Torus{static_cast<float>(c.R), static_cast<float>(c.r)},
        SDF::Warp::Twist{c.n, static_cast<float>(c.A),
                         static_cast<float>(c.R)}};
    for (int i = 0; i < 96; ++i) {
      const double theta = rand_uniform(rng, 0.0f, 2.0f * math::PI_F);
      const double phi = rand_uniform(rng, 0.0f, 2.0f * math::PI_F);
      const double radius = c.r * rand_uniform(rng, 1.01f, 1.8f);
      const double radial = c.R + radius * std::cos(phi);
      const math::Vector p(static_cast<float>(radial * std::cos(theta)),
                           static_cast<float>(radius * std::sin(phi) +
                                              c.A * std::sin(c.n * theta)),
                           static_cast<float>(radial * std::sin(theta)));
      const double truth = twisted_torus_distance(p, c.R, c.r, c.n, c.A, 20000);
      const float distance = volume.distance(p);
      HS_EXPECT_TRUE(std::isfinite(distance));
      if (distance > 0)
        HS_EXPECT_LE(distance, truth + 1e-4);
      correction_needed += volume.raw_distance(p) > truth + 1e-4;
    }
  }
  HS_EXPECT_GT(correction_needed, 0);
}

/**
 * @brief Verifies WarpedVolume::bounding_distance never over-estimates the
 *        distance to the warped surface, over randomized points and torus/twist
 *        parameter sets.
 * @details The fast path returns this bound directly.
 */
inline void test_warped_volume_bounding_distance_never_over_estimates() {
  struct Case {
    double R, r, A;
    int n;
  };
  const Case cases[] = {
      {1.0, 0.3, 0.2, 3},    // nominal
      {1.0, 0.3, 0.0, 3},    // zero amplitude
      {1.0, 0.3, 0.2, 0},    // zero twist
      {1.0, 0.31, 2.5, 5},   // amplitude well past the major radius
      {0.45, 0.14, 0.35, 2}, // Raymarch proportions, default twist
      {0.45, 0.14, 0.35, 8}, // Raymarch proportions, max twist
      {2.0, 0.05, 0.9, 7},   // thin tube
      {0.5, 0.45, 0.1, 1},   // fat tube
  };

  int violations = 0;
  hs::Pcg32 rng(0x5eed1337u);
  const auto next = [&rng]() {
    return static_cast<double>(rand_uniform(rng, 0.0f, 1.0f));
  };

  for (const Case &c : cases) {
    SDF::WarpedVolume<SDF::Torus, SDF::Warp::Twist> wv{
        SDF::Torus{static_cast<float>(c.R), static_cast<float>(c.r)},
        SDF::Warp::Twist{c.n, static_cast<float>(c.A),
                         static_cast<float>(c.R)}};
    const double reach = c.R + c.r + c.A + 1.0;

    for (int k = 0; k < 240; ++k) {
      math::Vector p;
      const int kind = k % 4;
      if (kind == 0) {
        p = math::Vector(0.0f, static_cast<float>((next() * 2 - 1) * reach),
                         0.0f);
      } else if (kind == 1 || kind == 2) {
        // On or inside the tube: radius scaled to at most the minor radius.
        const double t = next() * 2 * PI_DBL, ph = next() * 2 * PI_DBL;
        const double rr = c.r * (kind == 1 ? next() * 0.9 : 1.0);
        const double X = c.R + rr * std::cos(ph);
        p = math::Vector(
            static_cast<float>(X * std::cos(t)),
            static_cast<float>(rr * std::sin(ph) + c.A * std::sin(c.n * t)),
            static_cast<float>(X * std::sin(t)));
      } else {
        const double px = (next() * 2 - 1) * reach;
        const double py = (next() * 2 - 1) * reach;
        const double pz = (next() * 2 - 1) * reach;
        p = math::Vector(static_cast<float>(px), static_cast<float>(py),
                         static_cast<float>(pz));
      }

      const double bd = wv.bounding_distance(p);
      const double truth = twisted_torus_distance(p, c.R, c.r, c.n, c.A, 20000);
      if (!std::isfinite(bd) || bd - truth > 1e-4)
        ++violations;
      // The bound must also never exceed the slow path's raw warped distance.
      if (bd - static_cast<double>(wv.raw_distance(p)) > 1e-4)
        ++violations;
    }
  }
  HS_EXPECT_EQ(violations, 0);
}

/**
 * @brief Verifies the Lipschitz-corrected path returns raw/lipschitz on a
 *        near-surface outside point (not on the bounding fast-path).
 */
inline void test_warped_volume_distance_matches_lipschitz_correction() {
  SDF::Torus torus{1.0f, 0.3f};
  SDF::Warp::Twist tw{3, 0.2f, 1.0f};
  SDF::WarpedVolume<SDF::Torus, SDF::Warp::Twist> wv{torus, tw};

  // Just outside the outer rim at θ=0: small positive base distance lands off
  // the bounding fast-path and triggers the Lipschitz divide.
  math::Vector p(1.4f, 0.1f, 0.0f);
  float raw = wv.raw_distance(p);
  HS_EXPECT_TRUE(raw > 0.0f);
  auto ctx = tw.make_ctx(p);
  float lip = tw.lipschitz(p, ctx);
  HS_EXPECT_TRUE(lip > 1.0f);
  HS_EXPECT_NEAR(wv.distance(p), raw / lip, 1e-4f);
}

/**
 * @brief Verifies Twist::correct_normal returns a unit vector, is identity at
 *        twist 0, and reproduces the gradient of the warped field.
 */
inline void test_twist_correct_normal_unit_length() {
  math::Vector base_n = math::Vector(0.6f, 0.8f, 0.0f); // already unit

  SDF::Warp::Twist flat{0, 0.3f, 1.0f};
  math::Vector cf =
      flat.correct_normal(math::Vector(1, 0.2f, 0.5f), base_n,
                          flat.make_ctx(math::Vector(1, 0.2f, 0.5f)));
  HS_EXPECT_VEC(cf, base_n, 1e-6f);

  SDF::Warp::Twist tw{4, 0.25f, 1.0f};
  for (float x = -1.0f; x <= 1.0f; x += 0.5f)
    for (float z = -1.0f; z <= 1.0f; z += 0.5f) {
      math::Vector p(x, 0.3f, z);
      math::Vector c = tw.correct_normal(p, base_n, tw.make_ctx(p));
      HS_EXPECT_NEAR(c.length(), 1.0f, 1e-4f);
    }

  // The analytic chain rule must match central differences of raw_distance,
  // the field the normal is the gradient of.
  const float R = 2.0f, r = 0.5f, A = 0.3f;
  SDF::WarpedVolume<SDF::Torus, SDF::Warp::Twist> wv{SDF::Torus{R, r},
                                                     SDF::Warp::Twist{2, A, R}};
  const float h = 1e-3f;
  for (int i = 0; i < 16; ++i) {
    const float theta = static_cast<float>(i) * 0.53f;
    const float phi = static_cast<float>(i) * 0.91f;
    const float s = R + r * std::cos(phi);
    // A point on the warped surface: torus surface, displaced as apply() does.
    math::Vector p(s * std::cos(theta),
                   r * std::sin(phi) + A * std::sin(2.0f * theta),
                   s * std::sin(theta));
    math::Vector grad((wv.raw_distance(p + math::Vector(h, 0, 0)) -
                       wv.raw_distance(p - math::Vector(h, 0, 0))),
                      (wv.raw_distance(p + math::Vector(0, h, 0)) -
                       wv.raw_distance(p - math::Vector(0, h, 0))),
                      (wv.raw_distance(p + math::Vector(0, 0, h)) -
                       wv.raw_distance(p - math::Vector(0, 0, h))));
    HS_EXPECT_VEC(wv.normal(p), grad.normalized(), 5e-3f);
  }
}
