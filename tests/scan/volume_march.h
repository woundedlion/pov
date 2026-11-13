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

/** @brief Capped ray samples remain on the ray and bound longer traces. */
inline void test_volume_capped_trace_properties() {
  int rays = 0;
  auto compare = [&](const auto &shape, const math::Vector &ro,
                     const math::Vector &vd, float radius, float aa) {
    math::Vector full_position;
    const float full_distance = Scan::Volume::trace_closest(
        shape, ro, vd, radius, 128, aa, full_position);
    for (int steps : {1, 14, 18, 40}) {
      math::Vector position;
      const float distance = Scan::Volume::trace_closest(shape, ro, vd, radius,
                                                         steps, aa, position);
      HS_EXPECT_NEAR(distance, shape.distance(position), 1e-6f);
      HS_EXPECT_LT(math::cross(position - ro, vd).length(), 1e-6f);
      HS_EXPECT_GE(distance + 1e-6f, full_distance);
      ++rays;
    }
  };
  for (int twist : {0, 1, 2, 4, 7, 8})
    for (float scale : {0.08f, 0.3f, 1.0f}) {
      SDF::WarpedVolume<SDF::Torus, SDF::Warp::Twist> torus{
          {0.45f * scale, 0.14f * scale},
          {twist, 0.35f * scale, 0.45f * scale}};
      torus.precision = 0.14f * scale;
      for (int angle = 0; angle < 5; ++angle) {
        const auto q = math::make_rotation(
            math::Vector(0.3f, 1.0f, -0.2f).normalized(), angle * 0.59f);
        const auto vd = math::rotate(math::Vector(0, 0, -1), q);
        for (int y = -16; y <= 16; ++y)
          for (int x = -16; x <= 16; ++x)
            compare(
                torus,
                math::rotate(math::Vector(x * 0.045f * scale,
                                          y * 0.045f * scale, 0.72f * scale),
                             q),
                vd, 0.72f * scale, 0.07f * scale);
      }
    }
  const TwoSphereSDF PAIR{math::Vector(0, 0, 0.2f), 0.18f,
                          math::Vector(0.28f, 0, -0.2f), 0.3f};
  for (int x = -20; x <= 20; ++x)
    for (int y = -20; y <= 20; ++y)
      compare(PAIR,
              math::Vector(0.0357f + x * 0.0003f, 0.1797f + y * 0.0003f, 1),
              math::Vector(0, 0, -1), 0.6f, 0.01f);
  HS_EXPECT_GT(rays, 100000);
}

/** @brief Checks silhouettes and occlusion against a dense fixed-step march. */
inline void test_volume_dense_reference() {
  const math::Vector DIRECTION(0, 0, -1);
  const float AA = 0.01f;
  const float THRESHOLD = AA * Scan::Volume::SOLID_FRACTION;
  int band_entries = 0;
  int soft_entries = 0;
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
      ++band_entries;
      HS_EXPECT_NEAR(distance, expected.distance, 0.0003f);
      auto background =
          VolumeReference::behind(shape, expected.position, DIRECTION, radius);
      auto occluder = Scan::Volume::probe_occluder(shape, expected.position,
                                                   DIRECTION, radius, THRESHOLD,
                                                   AA, expected.distance);
      soft_entries += occluder.soft > 0.0f;
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
  std::printf("  dense volume oracle band=%d soft=%d\n", band_entries,
              soft_entries);
  HS_EXPECT_GE(band_entries, 50);
  HS_EXPECT_GE(soft_entries, 10);
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
 *        the identity, the direction maps through the inverse rotation, and distance()
 *        delegates to the wrapped SDF.
 */
inline void test_transformed_volume_world_local_roundtrip() {
  SphereSDF sphere{0.3f};
  const math::Vector CENTER(0.2f, -0.5f, 0.8f);
  const math::Quaternion q =
      math::make_rotation(math::Vector(0.3f, 1.0f, -0.2f).normalized(), 0.7f);
  using Volume = Scan::TransformedVolume<decltype(sphere)>;
  static_assert(
      std::is_constructible_v<Volume, const decltype(sphere) &,
                              const math::Vector &, const math::Quaternion &>);
  static_assert(
      !std::is_constructible_v<Volume, decltype(sphere) &&,
                               const math::Vector &, const math::Quaternion &>);
  Scan::TransformedVolume vol(sphere, CENTER, q);

  // bounds_center maps to the local origin (the cull precondition).
  math::Vector local_bc = vol.origin_to_local(CENTER);
  HS_EXPECT_NEAR(local_bc.length(), 0.0f, 1e-5f);

  // A local point pushed out to world and back is recovered exactly.
  const math::Vector LP(0.1f, 0.2f, -0.25f);
  math::Vector world = CENTER + math::rotate(LP, q);
  math::Vector back = vol.origin_to_local(world);
  HS_EXPECT_NEAR(back.x, LP.x, 1e-5f);
  HS_EXPECT_NEAR(back.y, LP.y, 1e-5f);
  HS_EXPECT_NEAR(back.z, LP.z, 1e-5f);

  // The ray direction maps through the inverse rotation.
  const math::Vector VD(0.0f, 0.0f, -1.0f);
  auto [lro, lvd] = vol.ray_to_local(world, VD);
  HS_EXPECT_NEAR(lro.x, LP.x, 1e-5f);
  HS_EXPECT_NEAR(lro.y, LP.y, 1e-5f);
  HS_EXPECT_NEAR(lro.z, LP.z, 1e-5f);
  HS_EXPECT_NEAR(lvd.length(), 1.0f, 1e-5f);
  HS_EXPECT_VEC(lvd, math::rotate(VD, q.conjugate()), 1e-5f);

  // distance() forwards straight to the wrapped SDF.
  HS_EXPECT_NEAR(vol.distance(LP), sphere.distance(LP), 1e-6f);
}

/**
 * @brief Verifies Volume::draw ray-marches a sphere SDF into a bounded silhouette
 *        whose every shaded fragment's hit registers land on the surface, on the
 *        camera-facing cap.
 * @details frag.pos (closest_local) and frag.size (closest_d) lie within the
 * AA band of the surface; rays travel along -Z.
 */
inline void test_volume_raymarch_silhouette_and_registers() {
  constexpr int W = 96, H = 64;
  const math::Vector CENTER(0.0f, 0.0f, 1.0f); // bounds centre in LED space
  const float BOUNDS_RADIUS = 0.35f;
  const float SPHERE_R = 0.28f; // < bounds so the SDF fits the cull sphere
  const float AA_WIDTH = 0.01f;

  SphereSDF sphere{SPHERE_R};
  Scan::TransformedVolume vol(sphere, CENTER, math::Quaternion());

  hs_test::StubEffect fx(W, H);
  VolumeSink sink;

  int hits = 0;
  float max_surf_err = 0.0f; // worst |‖pos‖ - radius| over all hits
  float max_reg_d = 0.0f;    // worst |frag.size| (closest_d) over all hits
  math::Vector centroid_sum(0.0f, 0.0f, 0.0f);
  {
    Canvas c(fx);
    Scan::Volume::draw<W, H>(
        sink, c, CENTER, BOUNDS_RADIUS, vol,
        [&](const math::Vector &loc, Fragment &frag) {
          ++hits;
          HS_EXPECT_VEC(frag.pos, loc, 0.0f);
          max_surf_err =
              fold_worst(max_surf_err, std::fabs(frag.pos.length() - SPHERE_R));
          max_reg_d = fold_worst(max_reg_d, std::fabs(frag.size));
          centroid_sum = centroid_sum + loc;
          frag.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
        },
        /*max_steps=*/24, AA_WIDTH);
  }

  // Real silhouette: some fragments, never the whole canvas; plotted ⊆ shaded.
  HS_EXPECT_GT(hits, 0);
  HS_EXPECT_LT((size_t)hits, (size_t)(W * H));
  HS_EXPECT_GT(sink.plotted.size(), (size_t)0);
  HS_EXPECT_LE(sink.plotted.size(), (size_t)hits);
  const float PROJECTED_DISC_PIXELS =
      math::PI_F * SPHERE_R * SPHERE_R /
      ((math::TWO_PI_F / W) * (math::PI_F / (H - 1)));
  HS_EXPECT_GT(sink.plotted.size(),
               static_cast<size_t>(PROJECTED_DISC_PIXELS * 0.5f));
  for (const auto &[x, y] : sink.plotted) {
    const math::Vector point = math::pixel_to_vector<W, H>(x, y);
    HS_EXPECT_GT(point.z, 0.0f);
    HS_EXPECT_LT(std::hypot(point.x, point.y), SPHERE_R + AA_WIDTH + 1e-5f);
  }

  // Every shaded fragment is a genuine surface hit inside the AA band.
  HS_EXPECT_LE(max_surf_err, AA_WIDTH + 1e-3f);
  HS_EXPECT_LE(max_reg_d, AA_WIDTH);

  // The hit centroid is on the camera-facing (+Z) cap.
  math::Vector centroid = centroid_sum * (1.0f / static_cast<float>(hits));
  HS_EXPECT_GT(centroid.z, 0.1f);
}

/**
 * @brief Verifies Volume::draw antialiases a self-occlusion edge over the surface
 *        behind it rather than fading the edge to black.
 * @details Background coverage precedes partial foreground coverage at each
 * silhouette overlap, including soft corners.
 */
inline void test_volume_draw_occluded_edge_blends_over_background() {
  constexpr int W = 96, H = 64;
  const math::Vector CENTER(0.0f, 0.0f, 1.0f);
  const float BOUNDS_RADIUS = 0.65f;
  const float AA_WIDTH = 0.08f;

  TwoSphereSDF shape{math::Vector(0.0f, 0.0f, 0.20f), 0.18f,
                     math::Vector(0.28f, 0.0f, -0.20f), 0.30f};
  Scan::TransformedVolume vol(shape, CENTER, math::Quaternion());

  hs_test::StubEffect fx(W, H);
  VolumeSink sink;
  {
    Canvas c(fx);
    Scan::Volume::draw<W, H>(
        sink, c, CENTER, BOUNDS_RADIUS, vol,
        [&](const math::Vector &, Fragment &frag) {
          frag.color = Color4(Pixel(60000, 60000, 60000), 1.0f);
        },
        /*max_steps=*/40, AA_WIDTH);
  }

  int occ_pairs = 0;
  int solid_pairs = 0, soft_pairs = 0;
  for (size_t i = 1; i < sink.plotted.size(); ++i) {
    if (sink.plotted[i] == sink.plotted[i - 1]) {
      ++occ_pairs;
      const float bg = sink.alpha[i - 1], fg = sink.alpha[i];
      HS_EXPECT_GT(bg, 0.0f);
      HS_EXPECT_LE(bg, 1.0f);
      HS_EXPECT_GT(fg, 0.0f);
      HS_EXPECT_LT(fg, 1.0f);
      solid_pairs += bg == 1.0f;
      soft_pairs += bg < 1.0f;
    }
  }

  HS_EXPECT_GT(occ_pairs, 0);
  HS_EXPECT_GT(solid_pairs, 0);
  HS_EXPECT_GT(soft_pairs, 0);
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
  const float AA_WIDTH = 0.01f;
  TwoSphereSDF shape{math::Vector(0.0f, 0.0f, 0.20f), 0.18f,
                     math::Vector(0.0f, 0.0f, -0.20f), 0.30f};

  // Grazing ray: passes 0.005 outside the foreground silhouette, then through
  // the background sphere's interior.
  const math::Vector RO(0.185f, 0.0f, 1.0f);
  const math::Vector VD(0.0f, 0.0f, -1.0f);

  for (int max_steps : {14, 40}) {
    math::Vector closest_local;
    float closest_d = Scan::Volume::trace_closest(
        shape, RO, VD, 0.5f, max_steps, AA_WIDTH, closest_local);
    HS_EXPECT_GT(closest_d, 0.003f);
    HS_EXPECT_LT(closest_d, 0.0075f);
    HS_EXPECT_GT(closest_local.z, 0.1f);
  }
}

/**
 * @brief Verifies probe_occluder reports the grazed background edge's own
 *        closest-approach point, so the corner fill is shaded on the background
 *        surface rather than reusing the foreground fragment.
 * @details The ray passes the corner where the two silhouettes cross, inside
 * both AA bands; the probe is non-solid with `behind` at the background
 * sphere's edge (negative z).
 */
inline void test_volume_probe_occluder_reports_background_graze_point() {
  const float AA_WIDTH = 0.01f;
  const float HIT_THRESHOLD = AA_WIDTH * Scan::Volume::SOLID_FRACTION;
  TwoSphereSDF shape{math::Vector(0.0f, 0.0f, 0.20f), 0.18f,
                     math::Vector(0.28f, 0.0f, -0.20f), 0.30f};

  // Silhouette-crossing corner, offset outward of both circles by ~0.0033.
  const math::Vector RO(0.0357f, 0.1797f, 1.0f);
  const math::Vector VD(0.0f, 0.0f, -1.0f);

  math::Vector closest_local;
  float closest_d = Scan::Volume::trace_closest(shape, RO, VD, 0.6f, 40,
                                                AA_WIDTH, closest_local);
  HS_EXPECT_GT(closest_d, HIT_THRESHOLD);
  HS_EXPECT_LT(closest_d, AA_WIDTH);
  HS_EXPECT_GT(closest_local.z, 0.1f);

  auto occ = Scan::Volume::probe_occluder(shape, closest_local, VD, 0.6f,
                                          HIT_THRESHOLD, AA_WIDTH);
  HS_EXPECT_FALSE(occ.solid);
  // The ray passes 0.0033 outside the background sphere, so the analytic
  // coverage is quintic(1 - (0.0033 - 0.001)/0.009) ~= 0.89.
  HS_EXPECT_GT(occ.soft, 0.8f);
  HS_EXPECT_LT(occ.soft, 0.95f);
  HS_EXPECT_LT(occ.behind.z, -0.05f);
}

/**
 * @brief Verifies overrelaxed sphere tracing never steps over a surface.
 * @details Sweeps rays across a twisted torus whose Lipschitz-divided distance
 * badly underestimates the true one and compares each trace against a dense
 * fixed-step scan of the same march-safe distance along the ray. A dense-scan
 * distance well inside the AA band reports a hit; no ray reports a hit the
 * dense scan cannot corroborate.
 */
inline void test_volume_trace_closest_overrelax_never_skips_surface() {
  SDF::WarpedVolume<SDF::Torus, SDF::Warp::Twist> torus{{0.45f, 0.14f},
                                                        {2, 0.35f, 0.45f}};
  const float BOUNDS_RADIUS = 0.72f;
  const float AA_WIDTH = 0.07f;
  const math::Vector VD(0.0f, 0.0f, -1.0f);

  int covered = 0, hits = 0;
  for (int iy = -14; iy <= 14; ++iy) {
    for (int ix = -14; ix <= 14; ++ix) {
      math::Vector ro(ix * 0.045f, iy * 0.045f, BOUNDS_RADIUS);

      math::Vector closest_local;
      float closest_d = Scan::Volume::trace_closest(
          torus, ro, VD, BOUNDS_RADIUS, 18, AA_WIDTH, closest_local);

      float ref_min = FLT_MAX;
      VolumeReference::march(
          torus, ro, VD, BOUNDS_RADIUS,
          [&](VolumeReference::Sample sample) {
            ref_min = std::min(ref_min, sample.distance);
            return true;
          },
          4000);

      if (closest_d < AA_WIDTH) {
        ++hits;
        // The dense scan of the same distance must corroborate the hit.
        HS_EXPECT_LT(ref_min, AA_WIDTH);
      }
      // Rays grazing the band edge may legitimately land either side of it;
      // anything this far inside is unambiguous coverage the march owes.
      if (ref_min < AA_WIDTH * (1.0f - Scan::Volume::SOLID_FRACTION)) {
        ++covered;
        HS_EXPECT_LT(closest_d, AA_WIDTH);
      }
    }
  }
  HS_EXPECT_GT(covered, 400);
  HS_EXPECT_GT(hits, covered);
}
