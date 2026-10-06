/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Scan::Volume / TransformedVolume — orthographic ray-march
// ============================================================================

/**
 * @brief Minimal analytic SDF: a sphere of `radius` centred at the local origin.
 * @details Supplies the distance() half of the Volume shape concept; wrapped in a
 * TransformedVolume to gain ray_to_local()/origin_to_local().
 */
struct SphereSDF {
  float radius; /**< Sphere radius in local units. */
  /**
   * @brief Signed distance to the sphere surface.
   * @param p Query point in local space.
   * @return |p| - radius (negative inside, zero on the surface).
   */
  float distance(const math::Vector &p) const { return p.length() - radius; }
};

/**
 * @brief Analytic SDF of two spheres unioned by min(): a small foreground sphere
 *        floating in front of a larger background sphere along the view axis.
 * @details The foreground silhouette has empty space immediately behind its edge
 * and the background surface a short march deeper, so a grazing ray at that edge
 * self-occludes the background.
 */
struct TwoSphereSDF {
  math::Vector fg_center; /**< Foreground sphere centre (nearer the camera). */
  float fg_radius;        /**< Foreground sphere radius. */
  math::Vector
      bg_center;   /**< Background sphere centre (deeper along the ray). */
  float bg_radius; /**< Background sphere radius. */
  /**
   * @brief Signed distance to the union of the two spheres.
   * @param p Query point in local space.
   * @return The nearer of the two surface distances.
   */
  float distance(const math::Vector &p) const {
    return std::min((p - fg_center).length() - fg_radius,
                    (p - bg_center).length() - bg_radius);
  }
};

/** @brief Counts signed-distance evaluations in the scalar volume oracle. */
template <typename Shape> struct CountedVolume {
  const Shape &shape;
  mutable int samples = 0;
  float distance(const math::Vector &p) const {
    ++samples;
    return shape.distance(p);
  }
};

/** @brief Compares scalar ray state with the vector-accumulating baseline. */
inline void test_volume_scalar_state_differential() {
  float max_distance = 0.0f, max_position = 0.0f, max_coverage = 0.0f;
  float max_probe_position = 0.0f, max_probe_coverage = 0.0f;
  int rays = 0, limited = 0, halos = 0, background_grazes = 0;
  int solid_changes = 0;
  auto compare = [&](const auto &shape, const math::Vector &ro,
                     const math::Vector &vd, float radius, int steps,
                     float aa) {
    CountedVolume counted{shape};
    math::Vector reference_p, actual_p;
    float reference_d = VolumeScalarRegression::trace_closest(
        counted, ro, vd, radius, steps, aa, reference_p);
    limited += counted.samples == steps;
    float actual_d =
        Scan::Volume::trace_closest(shape, ro, vd, radius, steps, aa, actual_p);
    ++rays;
    float reference_alpha =
        Scan::volume_edge_coverage(reference_d, aa * 0.1f, aa);
    float actual_alpha = Scan::volume_edge_coverage(actual_d, aa * 0.1f, aa);
    max_distance =
        hs_test::fold_worst(max_distance, fabsf(reference_d - actual_d));
    if (reference_alpha > 0.0f || actual_alpha > 0.0f)
      max_position =
          hs_test::fold_worst(max_position, (reference_p - actual_p).length());
    max_coverage = hs_test::fold_worst(max_coverage,
                                       fabsf(reference_alpha - actual_alpha));
    if (reference_d > aa * 0.1f && reference_d < aa && actual_d > aa * 0.1f &&
        actual_d < aa) {
      ++halos;
      auto reference_occlusion = VolumeScalarRegression::probe_occluder(
          shape, reference_p, vd, radius, aa * 0.1f, aa);
      auto actual_occlusion = Scan::Volume::probe_occluder(
          shape, actual_p, vd, radius, aa * 0.1f, aa, actual_d);
      solid_changes += reference_occlusion.solid != actual_occlusion.solid;
      if (reference_occlusion.solid || reference_occlusion.soft > 0.0f)
        max_probe_position = hs_test::fold_worst(
            max_probe_position,
            (reference_occlusion.behind - actual_occlusion.behind).length());
      max_probe_coverage = hs_test::fold_worst(
          max_probe_coverage,
          fabsf(reference_occlusion.soft - actual_occlusion.soft));
      background_grazes +=
          !reference_occlusion.solid && reference_occlusion.soft > 0.0f;
    }
  };
  for (int twist : {0, 1, 2, 4, 7, 8}) {
    for (float scale : {0.08f, 0.3f, 1.0f}) {
      SDF::WarpedVolume<SDF::Torus, SDF::Warp::Twist> torus{
          {0.45f * scale, 0.14f * scale},
          {twist, 0.35f * scale, 0.45f * scale}};
      torus.precision = 0.14f * scale;
      for (int angle = 0; angle < 5; ++angle) {
        math::Quaternion q = math::make_rotation(
            math::Vector(0.3f, 1.0f, -0.2f).normalized(), angle * 0.59f);
        math::Vector vd = math::rotate(math::Vector(0, 0, -1), q);
        for (int steps : {1, 14, 18, 40}) {
          for (int y = -16; y <= 16; ++y) {
            for (int x = -16; x <= 16; ++x) {
              math::Vector ro =
                  math::rotate(math::Vector(x * 0.045f * scale,
                                            y * 0.045f * scale, 0.72f * scale),
                               q);
              compare(torus, ro, vd, 0.72f * scale, steps, 0.07f * scale);
            }
          }
        }
      }
    }
  }
  TwoSphereSDF pair{math::Vector(0, 0, 0.2f), 0.18f,
                    math::Vector(0.28f, 0, -0.2f), 0.3f};
  for (int x = -20; x <= 20; ++x)
    for (int y = -20; y <= 20; ++y)
      compare(pair,
              math::Vector(0.0357f + x * 0.0003f, 0.1797f + y * 0.0003f, 1),
              math::Vector(0, 0, -1), 0.6f, 40, 0.01f);
  printf(
      "scalar ray differential: %d rays, %d limited, %d halos, %d background "
      "grazes; max distance %.9g, position %.9g, coverage %.9g, probe "
      "position %.9g, probe coverage %.9g, solid changes %d\n",
      rays, limited, halos, background_grazes, max_distance, max_position,
      max_coverage, max_probe_position, max_probe_coverage, solid_changes);
  HS_EXPECT_GT(limited, 1000);
  HS_EXPECT_GT(halos, 1000);
  HS_EXPECT_GT(background_grazes, 100);
  HS_EXPECT_LT(max_distance, 1e-5f);
  HS_EXPECT_LT(max_position, 1e-4f);
  HS_EXPECT_LT(max_coverage, 1e-4f);
  HS_EXPECT_LT(max_probe_position, 1e-4f);
  HS_EXPECT_LT(max_probe_coverage, 1e-4f);
  HS_EXPECT_EQ(solid_changes, 0);
}

/** @brief Checks silhouettes and occlusion against a dense fixed-step march. */
inline void test_volume_dense_reference() {
  const math::Vector DIRECTION(0, 0, -1);
  const float AA = 0.01f;
  const float THRESHOLD = AA * 0.1f;
  auto compare = [&](const auto &shape, const math::Vector &origin,
                     float radius) {
    auto expected =
        VolumeReference::trace(shape, origin, DIRECTION, radius, AA);
    math::Vector actual;
    float distance = Scan::Volume::trace_closest(shape, origin, DIRECTION,
                                                 radius, 128, AA, actual);
    HS_EXPECT_NEAR(VolumeReference::coverage(expected.distance, THRESHOLD, AA),
                   Scan::volume_edge_coverage(distance, THRESHOLD, AA), 0.015f);
    if (expected.distance > THRESHOLD && expected.distance < AA) {
      HS_EXPECT_NEAR(distance, expected.distance, 0.0003f);
      auto background =
          VolumeReference::behind(shape, expected.position, DIRECTION, radius);
      auto occluder = Scan::Volume::probe_occluder(shape, expected.position,
                                                   DIRECTION, radius, THRESHOLD,
                                                   AA, expected.distance);
      HS_EXPECT_EQ(occluder.solid, background.distance < THRESHOLD);
      if (!occluder.solid)
        HS_EXPECT_NEAR(
            occluder.soft,
            VolumeReference::coverage(background.distance, THRESHOLD, AA),
            0.03f);
    }
  };
  for (int i = -8; i <= 20; ++i) {
    float offset = i * 0.0005f;
    compare(SphereSDF{0.3f}, math::Vector(0.3f + offset, 0, 0.7f), 0.7f);
    compare(SDF::Torus{0.3f, 0.1f}, math::Vector(0.4f + offset, 0, 0.7f), 0.7f);
    TwoSphereSDF pair{math::Vector(0, 0, 0.2f), 0.18f,
                      math::Vector(0.28f, 0, -0.2f), 0.3f};
    compare(pair, math::Vector(0.0357f + offset, 0.1797f, 1), 0.6f);
  }
}

/** @brief Pins closest-sample ownership at a nearly tied silhouette minimum. */
inline void test_volume_trace_nearly_tied_minimum() {
  SDF::WarpedVolume<SDF::Torus, SDF::Warp::Twist> torus{
      {0x1.7cafe8p-3f, 0x1.d9be78p-5f}, {7, 0x1.28170ap-3f, 0x1.7cafe8p-3f}};
  const float AA = 0x1.634edap-4f;
  const float RADIUS = 0x1.852f16p-2f;
  torus.precision = 2.0f * AA;
  const math::Vector ORIGIN(-0x1.c9d22p-2f, -0x1.5bf194p-3f, 0x1.2c3592p-3f);
  const math::Vector DIRECTION(0x1.acc1c2p-2f, 0x1.af3f2ep-1f, -0x1.5ba266p-2f);
  math::Vector expected, actual;
  float expected_d = VolumeScalarRegression::trace_closest(
      torus, ORIGIN, DIRECTION, RADIUS, 17, AA, expected);
  float actual_d = Scan::Volume::trace_closest(torus, ORIGIN, DIRECTION, RADIUS,
                                               17, AA, actual);
#if defined(HS_TEST_FAST_MATH)
  HS_EXPECT_NEAR(expected_d, actual_d, 1e-5f);
  HS_EXPECT_NEAR(expected.x, actual.x, 1e-4f);
  HS_EXPECT_NEAR(expected.y, actual.y, 1e-4f);
  HS_EXPECT_NEAR(expected.z, actual.z, 1e-4f);
#else
  HS_EXPECT_EQ(expected_d, actual_d);
  HS_EXPECT_EQ(expected.x, actual.x);
  HS_EXPECT_EQ(expected.y, actual.y);
  HS_EXPECT_EQ(expected.z, actual.z);
#endif
}

/**
 * @brief Capturing volume sink: records plotted pixel coordinates and alpha.
 * @details Only the integer-coordinate plot() records.
 */
struct VolumeSink {
  /** Pixel coordinates handed to the integer plot(). */
  std::vector<std::pair<int, int>> plotted;
  std::vector<float> alpha; /**< Composited alpha per plotted pixel. */
  void plot(Canvas &, int x, int y, const Pixel &, float, float a) {
    plotted.push_back({x, y});
    alpha.push_back(a);
  }
  void plot(Canvas &, float, float, const Pixel &, float, float) {}
  void plot(Canvas &, const math::Vector &, const Pixel &, float, float) {}
};

/**
 * @brief Verifies TransformedVolume's world<->local contract: the round trip is
 *        the identity, the mapped ray direction stays unit length, and distance()
 *        delegates to the wrapped SDF.
 */
inline void test_transformed_volume_world_local_roundtrip() {
  SphereSDF sphere{0.3f};
  const math::Vector center(0.2f, -0.5f, 0.8f);
  const math::Quaternion q =
      math::make_rotation(math::Vector(0.3f, 1.0f, -0.2f).normalized(), 0.7f);
  using Volume = Scan::TransformedVolume<decltype(sphere)>;
  static_assert(
      std::is_constructible_v<Volume, const decltype(sphere) &,
                              const math::Vector &, const math::Quaternion &>);
  static_assert(
      !std::is_constructible_v<Volume, decltype(sphere) &&,
                               const math::Vector &, const math::Quaternion &>);
  Scan::TransformedVolume vol(sphere, center, q);

  // bounds_center maps to the local origin (the cull precondition).
  math::Vector local_bc = vol.origin_to_local(center);
  HS_EXPECT_NEAR(local_bc.length(), 0.0f, 1e-5f);

  // A local point pushed out to world and back is recovered exactly.
  const math::Vector lp(0.1f, 0.2f, -0.25f);
  math::Vector world = center + math::rotate(lp, q);
  math::Vector back = vol.origin_to_local(world);
  HS_EXPECT_NEAR(back.x, lp.x, 1e-5f);
  HS_EXPECT_NEAR(back.y, lp.y, 1e-5f);
  HS_EXPECT_NEAR(back.z, lp.z, 1e-5f);

  // ray_to_local maps the origin like origin_to_local and keeps a unit direction
  // unit (rigid map, no scale) — the |local_vd| == 1 precondition.
  const math::Vector vd(0.0f, 0.0f, -1.0f);
  auto [lro, lvd] = vol.ray_to_local(world, vd);
  HS_EXPECT_NEAR(lro.x, lp.x, 1e-5f);
  HS_EXPECT_NEAR(lro.y, lp.y, 1e-5f);
  HS_EXPECT_NEAR(lro.z, lp.z, 1e-5f);
  HS_EXPECT_NEAR(lvd.length(), 1.0f, 1e-5f);

  // distance() forwards straight to the wrapped SDF.
  HS_EXPECT_NEAR(vol.distance(lp), sphere.distance(lp), 1e-6f);
}

/**
 * @brief Verifies Volume::draw ray-marches a sphere SDF into a bounded silhouette
 *        whose every shaded fragment's hit registers land on the surface, on the
 *        camera-facing cap.
 * @details The silhouette is non-empty and smaller than the canvas, with no more
 * plots than shades; each hit's frag.pos (closest_local) and frag.size
 * (closest_d) lie within the AA band of the surface; the hit centroid lies on
 * the +Z cap facing the camera (rays travel along -Z).
 */
inline void test_volume_raymarch_silhouette_and_registers() {
  constexpr int W = 96, H = 64;
  const math::Vector center(0.0f, 0.0f, 1.0f); // bounds centre in LED space
  const float bounds_radius = 0.35f;
  const float sphere_r = 0.28f; // < bounds so the SDF fits the cull sphere
  const float aa_width = 0.01f;

  SphereSDF sphere{sphere_r};
  Scan::TransformedVolume vol(sphere, center, math::Quaternion());

  hs_test::StubEffect fx(W, H);
  VolumeSink sink;

  int hits = 0;
  float max_surf_err = 0.0f; // worst |‖pos‖ - radius| over all hits
  float max_reg_d = 0.0f;    // worst |frag.size| (closest_d) over all hits
  math::Vector centroid_sum(0.0f, 0.0f, 0.0f);
  {
    Canvas c(fx);
    Scan::Volume::draw<W, H>(
        sink, c, center, bounds_radius, vol,
        [&](const math::Vector &loc, Fragment &frag) {
          ++hits;
          HS_EXPECT_VEC(frag.pos, loc, 0.0f);
          max_surf_err =
              fold_worst(max_surf_err, std::fabs(frag.pos.length() - sphere_r));
          max_reg_d = fold_worst(max_reg_d, std::fabs(frag.size));
          centroid_sum = centroid_sum + loc;
          frag.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
        },
        /*max_steps=*/24, aa_width);
  }

  // Real silhouette: some fragments, never the whole canvas; plotted ⊆ shaded.
  HS_EXPECT_GT(hits, 0);
  HS_EXPECT_LT((size_t)hits, (size_t)(W * H));
  HS_EXPECT_GT(sink.plotted.size(), (size_t)0);
  HS_EXPECT_LE(sink.plotted.size(), (size_t)hits);

  // Every shaded fragment is a genuine surface hit inside the AA band.
  HS_EXPECT_LE(max_surf_err, aa_width + 1e-3f);
  HS_EXPECT_LE(max_reg_d, aa_width);

  // The hit centroid is on the camera-facing (+Z) cap.
  math::Vector centroid = centroid_sum * (1.0f / static_cast<float>(hits));
  HS_EXPECT_GT(centroid.z, 0.1f);
}

/**
 * @brief Verifies Volume::draw antialiases a self-occlusion edge over the surface
 *        behind it rather than fading the edge to black.
 * @details Around the foreground silhouette Volume::draw plots the background
 * opaque, then blends the foreground over it by the edge coverage (0<α<1) at the
 * same pixel.
 */
inline void test_volume_draw_occluded_edge_blends_over_background() {
  constexpr int W = 96, H = 64;
  const math::Vector center(0.0f, 0.0f, 1.0f);
  const float bounds_radius = 0.50f;
  const float aa_width = 0.01f;

  // Small foreground sphere nearer the camera (+Z local); larger background
  // sphere deeper and wider so it sits behind the whole foreground silhouette.
  // The step budget lets the grazing ray stall on the foreground edge in the AA
  // band, so the occluder probe discovers the surface behind.
  TwoSphereSDF shape{math::Vector(0.0f, 0.0f, 0.20f), 0.18f,
                     math::Vector(0.0f, 0.0f, -0.20f), 0.30f};
  Scan::TransformedVolume vol(shape, center, math::Quaternion());

  hs_test::StubEffect fx(W, H);
  VolumeSink sink;
  {
    Canvas c(fx);
    Scan::Volume::draw<W, H>(
        sink, c, center, bounds_radius, vol,
        [&](const math::Vector &, Fragment &frag) {
          frag.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
        },
        /*max_steps=*/12, aa_width);
  }

  // Solid occlusion and soft corner fill both emit background/foreground pairs.
  // An opaque background distinguishes the solid occluder branch.
  int occ_pairs = 0;
  float bg_alpha = 0.0f, fg_alpha = 1.0f;
  for (size_t i = 1; i < sink.plotted.size(); ++i) {
    if (sink.plotted[i] == sink.plotted[i - 1]) {
      ++occ_pairs;
      if (occ_pairs == 1) {
        bg_alpha = sink.alpha[i - 1];
        fg_alpha = sink.alpha[i];
      }
    }
  }

  HS_EXPECT_GT(occ_pairs, 0);
  // Background laid down opaque; foreground blended over it as a partial edge.
  HS_EXPECT_EQ(bg_alpha, 1.0f);
  HS_EXPECT_GT(fg_alpha, 0.0f);
  HS_EXPECT_LT(fg_alpha, bg_alpha);
}

/**
 * @brief Verifies trace_closest stops at the first silhouette graze instead of
 *        letting an occluded surface behind it steal the closest approach.
 * @details A ray grazes the foreground sphere's edge inside the AA band and then
 * passes solidly through the background sphere. The returned distance is the
 * foreground graze and the returned point sits on the foreground silhouette,
 * independent of the step budget.
 */
inline void test_volume_trace_closest_stops_at_first_graze() {
  const float aa_width = 0.01f;
  TwoSphereSDF shape{math::Vector(0.0f, 0.0f, 0.20f), 0.18f,
                     math::Vector(0.0f, 0.0f, -0.20f), 0.30f};

  // Grazing ray: passes 0.005 outside the foreground silhouette, then through
  // the background sphere's interior.
  const math::Vector ro(0.185f, 0.0f, 1.0f);
  const math::Vector vd(0.0f, 0.0f, -1.0f);

  for (int max_steps : {14, 40}) {
    math::Vector closest_local;
    float closest_d = Scan::Volume::trace_closest(
        shape, ro, vd, 0.5f, max_steps, aa_width, closest_local);
    HS_EXPECT_GT(closest_d, 0.003f);
    HS_EXPECT_LT(closest_d, 0.0075f);
    HS_EXPECT_GT(closest_local.z, 0.1f);
  }
}

/**
 * @brief Verifies probe_occluder reports the grazed background edge's own
 *        closest-approach point, so the corner fill is shaded on the background
 *        surface rather than reusing the foreground fragment.
 * @details The ray passes through the corner where the two spheres' silhouettes
 * cross, inside both AA bands but hitting neither: the trace grazes the
 * foreground sphere, and the probe must come back non-solid with the graze
 * coverage and a `behind` point at the background sphere's edge (negative z),
 * not the foreground graze it was seeded with.
 */
inline void test_volume_probe_occluder_reports_background_graze_point() {
  const float aa_width = 0.01f;
  const float hit_threshold = aa_width * 0.1f;
  TwoSphereSDF shape{math::Vector(0.0f, 0.0f, 0.20f), 0.18f,
                     math::Vector(0.28f, 0.0f, -0.20f), 0.30f};

  // Silhouette-crossing corner, offset outward of both circles by ~0.0033.
  const math::Vector ro(0.0357f, 0.1797f, 1.0f);
  const math::Vector vd(0.0f, 0.0f, -1.0f);

  math::Vector closest_local;
  float closest_d = Scan::Volume::trace_closest(shape, ro, vd, 0.6f, 40,
                                                aa_width, closest_local);
  HS_EXPECT_GT(closest_d, hit_threshold);
  HS_EXPECT_LT(closest_d, aa_width);
  HS_EXPECT_GT(closest_local.z, 0.1f);

  auto occ = Scan::Volume::probe_occluder(shape, closest_local, vd, 0.6f,
                                          hit_threshold, aa_width);
  HS_EXPECT_FALSE(occ.solid);
  // The ray passes 0.0033 outside the background sphere, so the analytic
  // coverage is quintic(1 - (0.0033 - 0.001)/0.009) ~= 0.89; the parabolic
  // refinement must land near it (the coarse stride alone reads ~0.65).
  HS_EXPECT_GT(occ.soft, 0.8f);
  HS_EXPECT_LT(occ.soft, 0.95f);
  HS_EXPECT_LT(occ.behind.z, -0.05f);
}

/**
 * @brief Verifies overrelaxed sphere tracing never steps over a surface.
 * @details Sweeps rays across a twisted torus whose Lipschitz-divided distance
 * badly underestimates the true one and compares each trace against a dense
 * fixed-step scan of the same ray. A ray whose true closest approach lies well
 * inside the AA band reports a hit; no ray reports a hit the dense scan cannot
 * corroborate.
 */
inline void test_volume_trace_closest_overrelax_never_skips_surface() {
  SDF::WarpedVolume<SDF::Torus, SDF::Warp::Twist> torus{{0.45f, 0.14f},
                                                        {2, 0.35f, 0.45f}};
  const float bounds_radius = 0.72f;
  const float aa_width = 0.07f;
  const math::Vector vd(0.0f, 0.0f, -1.0f);

  int covered = 0, hits = 0;
  for (int iy = -14; iy <= 14; ++iy) {
    for (int ix = -14; ix <= 14; ++ix) {
      math::Vector ro(ix * 0.045f, iy * 0.045f, bounds_radius);

      math::Vector closest_local;
      float closest_d = Scan::Volume::trace_closest(
          torus, ro, vd, bounds_radius, 18, aa_width, closest_local);

      // Dense reference: true closest approach along the same ray segment.
      float ref_min = FLT_MAX;
      for (int s = 0; s <= 4000; ++s) {
        math::Vector p(ro.x, ro.y, ro.z - s * (2.0f * bounds_radius / 4000.0f));
        float d = torus.distance(p);
        if (d < ref_min)
          ref_min = d;
      }

      if (closest_d < aa_width) {
        ++hits;
        // A reported hit must correspond to a real approach on the ray.
        HS_EXPECT_LT(ref_min, aa_width);
      }
      // Rays grazing the band edge may legitimately land either side of it;
      // anything this far inside is unambiguous coverage the march owes.
      if (ref_min < aa_width * 0.9f) {
        ++covered;
        HS_EXPECT_LT(closest_d, aa_width);
      }
    }
  }
  HS_EXPECT_GT(covered, 400);
  HS_EXPECT_GT(hits, covered);
}
