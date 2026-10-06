/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Union — min of distances
// ============================================================================

/** @brief Verifies Union returns the min of member distances, picking the closest shape. */
inline void test_union_picks_closest_shape() {
  SDF::Line la(math::Vector(1, 0, 0), math::Vector(0, 0, 1), 0.1f);
  SDF::Line lb(math::Vector(-1, 0, 0), math::Vector(0, 0, -1), 0.1f);

  SDF::Union<SDF::Line, SDF::Line> u(la, lb);

  math::Vector mid_a =
      ((math::Vector(1, 0, 0) + math::Vector(0, 0, 1)) * 0.5f).normalized();
  auto r = SDF::distance_of(u, mid_a);
  HS_EXPECT_NEAR(r.dist, -0.1f, 5e-4f);

  math::Vector mid_b =
      ((math::Vector(-1, 0, 0) + math::Vector(0, 0, -1)) * 0.5f).normalized();
  auto r2 = SDF::distance_of(u, mid_b);
  HS_EXPECT_NEAR(r2.dist, -0.1f, 5e-4f);
}

// ============================================================================
// Subtract — max(A, -B)
// ============================================================================

/** @brief Verifies a point inside A but outside B stays inside the difference A - B. */
inline void test_subtract_inside_a_outside_b_remains_inside() {
  SDF::Line la(math::Vector(1, 0, 0), math::Vector(0, 0, 1), 0.2f);
  SDF::Line lb(math::Vector(-1, 0, 0), math::Vector(0, 0, -1), 0.1f);
  SDF::Subtract<SDF::Line, SDF::Line> s(la, lb);

  math::Vector mid_a =
      ((math::Vector(1, 0, 0) + math::Vector(0, 0, 1)) * 0.5f).normalized();
  auto r = SDF::distance_of(s, mid_a);
  HS_EXPECT_TRUE(r.dist < 0.0f);
}

/** @brief Verifies a point inside both A and B becomes outside the difference A - B. */
inline void test_subtract_inside_both_becomes_outside() {
  // Same line for A and B → A - A is empty everywhere.
  SDF::Line la(math::Vector(1, 0, 0), math::Vector(0, 0, 1), 0.1f);
  SDF::Subtract<SDF::Line, SDF::Line> s(la, la);
  math::Vector mid =
      ((math::Vector(1, 0, 0) + math::Vector(0, 0, 1)) * 0.5f).normalized();
  auto r = SDF::distance_of(s, mid);
  // max(dist(A), -dist(B)) = max(-0.1, 0.1) = 0.1 → outside.
  HS_EXPECT_NEAR(r.dist, 0.1f, 1e-3f);
}

/** @brief Verifies the AA size metric stays the minuend's when the subtrahend wins. */
inline void test_subtract_keeps_minuend_size_when_b_wins() {
  math::Vector p(1, 0, 0), q(0, 0, 1);
  SDF::Line la(p, q, 0.1f);
  SDF::Line lb(p, q, 0.4f);
  SDF::Subtract<SDF::Line, SDF::Line> s(la, lb);
  auto r = SDF::distance_of(s, ((p + q) * 0.5f).normalized());
  HS_EXPECT_NEAR(r.dist, 0.4f, 1e-3f);
  HS_EXPECT_NEAR(r.size, 0.1f, 1e-6f);
}

namespace sdf_interval_detail {
/**
 * @brief Mock SDF shape that emits a fixed (possibly unsorted, multi-) interval list.
 */
struct MockIntervalShape {
  static constexpr bool BLENDS_SMOOTHLY = true;
  const std::vector<std::pair<float, float>>
      *ivs; /**< Interval list this mock replays. */
  static constexpr bool is_solid =
      true; /**< Marks the mock as a solid fill shape. */
  /**
   * @brief Emits the stored intervals to the scanline sink.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam Out Interval-sink callable type taking (start, end).
   * @param out Sink invoked once per stored interval.
   * @return Always true: this mock definitively produces intervals.
   */
  template <int W, int H, typename Out>
  bool get_horizontal_intervals(int, Out out) const {
    for (const auto &p : *ivs)
      out(p.first, p.second);
    return true;
  }
  /**
   * @brief Claims every row, matching the unconditional interval emission.
   * @tparam H Canvas height in rows.
   * @return Row bounds spanning the canvas.
   */
  template <int H> SDF::Bounds get_vertical_bounds() const {
    return {0, H - 1};
  }
};

/**
 * @brief Mock that falls back to a full-row scan.
 * @details Returning false means "I cannot produce intervals — caller must scan
 *   the whole row with distance()", NOT that the shape covers the row.
 */
struct MockFullWidthShape {
  static constexpr bool is_solid =
      true; /**< Marks the mock as a solid fill shape. */
  /**
   * @brief Declines to emit intervals, forcing a full-row distance scan.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam Out Interval-sink callable type (unused).
   * @return Always false: caller must scan the whole row with distance().
   */
  template <int W, int H, typename Out>
  bool get_horizontal_intervals(int, Out) const {
    return false;
  }
};

/**
 * @brief Mock that violates the interval protocol by emitting before falling back.
 * @details The protocol requires a false return to emit nothing.
 */
struct MockEmitThenFullWidthShape {
  static constexpr bool is_solid =
      true; /**< Marks the mock as a solid fill shape. */
  /**
   * @brief Emits one span and then requests a full-row distance scan anyway.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam Out Interval-sink callable type taking (start, end).
   * @param out Sink invoked once before the fallback.
   * @return Always false.
   */
  template <int W, int H, typename Out>
  bool get_horizontal_intervals(int, Out out) const {
    out(10.0f, 30.0f);
    return false;
  }
};
} // namespace sdf_interval_detail

/**
 * @brief Verifies a solid B's spans never carve the minuend's scanline emission.
 * @details A child's spans bound its coverage, so every subtrahend is carved
 *   per pixel by max(A, -B), leaving A's spans intact whatever B emits.
 */
inline void test_subtract_solid_b_leaves_the_minuend_uncarved() {
  using P = std::pair<float, float>;
  using Mock = sdf_interval_detail::MockIntervalShape;
  std::vector<P> a_ivs = {{0.0f, 40.0f}, {60.0f, 100.0f}};
  std::vector<P> b_ivs = {{60.0f, 70.0f}, {20.0f, 30.0f}}; // unsorted, multi
  Mock A{&a_ivs}, B{&b_ivs};
  SDF::Subtract<Mock, Mock> s(A, B);

  std::vector<P> out;
  bool ok = s.get_horizontal_intervals<256, 128>(
      0, [&](float st, float en) { out.push_back({st, en}); });
  HS_EXPECT_TRUE(ok);

  HS_EXPECT_SIZE_OR_RETURN(out, 2);
  HS_EXPECT_NEAR(out[0].first, 0.0f, 1e-4f);
  HS_EXPECT_NEAR(out[0].second, 40.0f, 1e-4f);
  HS_EXPECT_NEAR(out[1].first, 60.0f, 1e-4f);
  HS_EXPECT_NEAR(out[1].second, 100.0f, 1e-4f);
}

/**
 * @brief Verifies a star notch inside the subtrahend's bounding cap still gets scanned.
 * @details Star emits its circumscribed disc as one span. The notch point is
 *   outside the star and inside the polygon, so its column must lie in an
 *   emitted span.
 */
inline void test_subtract_star_notch_columns_survive_the_carve() {
  using P = std::pair<float, float>;
  constexpr int W = 256, H = 128;
  // Axis +Z puts the shared center at column W/4, clear of the theta=0 seam.
  const math::Basis b{math::Vector(1, 0, 0), math::Vector(0, 0, 1),
                      math::Vector(0, 1, 0)};
  constexpr float outer = 0.5f; // star tip radius (radians)
  SDF::PlanarPolygon poly(b, 0.7f / (math::PI_F / 2.0f), 5,
                          0.0f); // apothem 0.566 > outer
  SDF::Star star(b, outer / (math::PI_F / 2.0f), 5, 0.0f);
  SDF::Subtract<SDF::PlanarPolygon, SDF::Star> s(poly, star);

  // Half a sector off a tip, just short of the tip radius: a notch.
  const float polar = 0.9f * outer, az = math::PI_F / 5.0f;
  math::Vector p = (b.v * std::cos(polar) +
                    (b.u * std::cos(az) + b.w * std::sin(az)) * std::sin(polar))
                       .normalized();
  HS_EXPECT_TRUE(SDF::distance_of(star, p).dist > 0.0f);
  HS_EXPECT_TRUE(SDF::distance_of(poly, p).dist < 0.0f);
  HS_EXPECT_TRUE(SDF::distance_of(s, p).dist < 0.0f);

  float phi = std::acos(hs::clamp(p.y, -1.0f, 1.0f));
  int y = static_cast<int>(phi * (H + hs::H_OFFSET - 1) / math::PI_F + 0.5f);
  float theta = std::atan2(p.z, p.x);
  if (theta < 0.0f)
    theta += math::TWO_PI_F;
  const float col = theta * W / math::TWO_PI_F;

  std::vector<P> out;
  bool ok = s.get_horizontal_intervals<W, H>(
      y, [&](float st, float en) { out.push_back({st, en}); });
  HS_EXPECT_TRUE(ok);
  bool covered = false;
  for (const auto &iv : out)
    covered = covered || (col >= iv.first && col <= iv.second);
  HS_EXPECT_TRUE(covered);
}

/**
 * @brief Verifies Subtract forwards A's intervals verbatim, including with empty B.
 * @details Subtract never consults B's intervals. No A span is dropped, merged
 *   or reordered.
 */
inline void test_subtract_empty_b_passes_a_through_verbatim() {
  using P = std::pair<float, float>;
  using Mock = sdf_interval_detail::MockIntervalShape;
  std::vector<P> a_ivs = {{50.0f, 60.0f}, {0.0f, 10.0f}}; // unsorted
  std::vector<P> b_ivs = {};
  Mock A{&a_ivs}, B{&b_ivs};
  SDF::Subtract<Mock, Mock> s(A, B);

  std::vector<P> out;
  s.get_horizontal_intervals<256, 128>(
      0, [&](float st, float en) { out.push_back({st, en}); });
  HS_EXPECT_SIZE_OR_RETURN(out, 2);
  HS_EXPECT_NEAR(out[0].first, 50.0f, 1e-4f);
  HS_EXPECT_NEAR(out[0].second, 60.0f, 1e-4f);
  HS_EXPECT_NEAR(out[1].first, 0.0f, 1e-4f);
  HS_EXPECT_NEAR(out[1].second, 10.0f, 1e-4f);
}

/**
 * @brief Verifies a subtrahend that cannot produce intervals costs the minuend nothing.
 * @details B's fallback neither widens the row to a full scan nor erases A.
 */
inline void test_subtract_full_width_b_still_emits_the_minuend() {
  using P = std::pair<float, float>;
  using MockA = sdf_interval_detail::MockIntervalShape;
  using MockB = sdf_interval_detail::MockFullWidthShape;
  std::vector<P> a_ivs = {{0.0f, 100.0f}};
  MockA A{&a_ivs};
  MockB B;
  SDF::Subtract<MockA, MockB> s(A, B);

  std::vector<P> out;
  bool ok = s.get_horizontal_intervals<256, 128>(
      0, [&](float st, float en) { out.push_back({st, en}); });
  HS_EXPECT_TRUE(ok);
  HS_EXPECT_SIZE_OR_RETURN(out, 1);
  HS_EXPECT_NEAR(out[0].first, 0.0f, 1e-4f);
  HS_EXPECT_NEAR(out[0].second, 100.0f, 1e-4f);
}

/** @brief Verifies Subtract forwards an unwrapped minuend band. */
inline void test_subtract_seam_straddle_forwards_minuend() {
  using P = std::pair<float, float>;
  using Mock = sdf_interval_detail::MockIntervalShape;
  std::vector<P> a_ivs = {{-10.0f, 10.0f}};  // seam band, negative frame
  std::vector<P> b_ivs = {{246.0f, 266.0f}}; // same band, [W,2W) frame (W=256)
  Mock A{&a_ivs}, B{&b_ivs};
  SDF::Subtract<Mock, Mock> s(A, B);

  std::vector<P> out;
  bool ok = s.get_horizontal_intervals<256, 128>(
      0, [&](float st, float en) { out.push_back({st, en}); });
  HS_EXPECT_TRUE(ok);
  HS_EXPECT_SIZE_OR_RETURN(out, 1);
  HS_EXPECT_NEAR(out[0].first, -10.0f, 1e-4f);
  HS_EXPECT_NEAR(out[0].second, 10.0f, 1e-4f);
}

/** @brief Verifies Subtract preserves all minuend spans. */
inline void test_subtract_many_arc_preserves_minuend() {
  using P = std::pair<float, float>;
  using Mock = sdf_interval_detail::MockIntervalShape;
  std::vector<P> a_ivs = {{-6.0f, 4.0f},    {250.0f, 262.0f}, // seam straddlers
                          {20.0f, 30.0f},   {40.0f, 50.0f},   {60.0f, 70.0f},
                          {80.0f, 90.0f},   {100.0f, 110.0f}, {120.0f, 130.0f},
                          {140.0f, 150.0f}, {160.0f, 170.0f}, {180.0f, 190.0f},
                          {200.0f, 210.0f}};
  while (a_ivs.size() < SDF::INTERVAL_SPAN_CAP) {
    const float start = static_cast<float>(a_ivs.size()) * 0.25f;
    a_ivs.push_back({start, start + 0.125f});
  }
  std::vector<P> b_ivs = {{45.0f, 55.0f}, {125.0f, 135.0f}};
  Mock A{&a_ivs}, B{&b_ivs};
  SDF::Subtract<Mock, Mock> s(A, B);

  std::vector<P> out;
  bool ok = s.get_horizontal_intervals<256, 128>(
      0, [&](float st, float en) { out.push_back({st, en}); });
  HS_EXPECT_TRUE(ok);

  HS_EXPECT_SIZE_OR_RETURN(out, a_ivs.size());
  for (size_t i = 0; i < out.size(); ++i) {
    HS_EXPECT_NEAR(out[i].first, a_ivs[i].first, 1e-4f);
    HS_EXPECT_NEAR(out[i].second, a_ivs[i].second, 1e-4f);
  }
}

// ============================================================================
// Intersection — max(A, B)
// ============================================================================

/** @brief Verifies Intersection is inside only where both children are inside. */
inline void test_intersection_requires_both_inside() {
  SDF::Line la(math::Vector(1, 0, 0), math::Vector(0, 0, 1), 0.3f);
  SDF::Line lb(math::Vector(1, 0, 0), math::Vector(0, 1, 0), 0.3f);
  SDF::Intersection<SDF::Line, SDF::Line> inter(la, lb);

  // Endpoint a is on both arcs.
  math::Vector a(1, 0, 0);
  auto r = SDF::distance_of(inter, a);
  HS_EXPECT_TRUE(r.dist < 0.0f);

  math::Vector far_pt(-1, 0, 0);
  auto r2 = SDF::distance_of(inter, far_pt);
  HS_EXPECT_TRUE(r2.dist > 0.0f);

  const math::Vector ONLY_A = math::Vector(1, 0, 1).normalized();
  const math::Vector ONLY_B = math::Vector(1, 1, 0).normalized();
  HS_EXPECT_LT(SDF::distance_of(la, ONLY_A).dist, 0.0f);
  HS_EXPECT_GT(SDF::distance_of(lb, ONLY_A).dist, 0.0f);
  HS_EXPECT_GT(SDF::distance_of(inter, ONLY_A).dist, 0.0f);
  HS_EXPECT_GT(SDF::distance_of(la, ONLY_B).dist, 0.0f);
  HS_EXPECT_LT(SDF::distance_of(lb, ONLY_B).dist, 0.0f);
  HS_EXPECT_GT(SDF::distance_of(inter, ONLY_B).dist, 0.0f);
}

/**
 * @brief Verifies an UNSORTED multi-interval child still yields a start-sorted intersection.
 * @details Intersection's merge-sweep advances the child lists in start order,
 *   so it sorts the seam-split lists itself before sweeping them.
 */
inline void test_intersection_unsorted_child_yields_sorted_result() {
  using P = std::pair<float, float>;
  using Mock = sdf_interval_detail::MockIntervalShape;
  std::vector<P> a_ivs = {{0.0f, 100.0f}};
  std::vector<P> b_ivs = {{60.0f, 80.0f}, {20.0f, 40.0f}}; // unsorted, multi
  Mock A{&a_ivs}, B{&b_ivs};
  SDF::Intersection<Mock, Mock> s(A, B);

  std::vector<P> out;
  bool ok = s.get_horizontal_intervals<256, 128>(
      0, [&](float st, float en) { out.push_back({st, en}); });
  HS_EXPECT_TRUE(ok);

  // [0,100] ∩ {[20,40],[60,80]} = [20,40],[60,80], start-sorted.
  HS_EXPECT_SIZE_OR_RETURN(out, 2);
  HS_EXPECT_NEAR(out[0].first, 20.0f, 1e-4f);
  HS_EXPECT_NEAR(out[0].second, 40.0f, 1e-4f);
  HS_EXPECT_NEAR(out[1].first, 60.0f, 1e-4f);
  HS_EXPECT_NEAR(out[1].second, 80.0f, 1e-4f);
}

/**
 * @brief Verifies a full-width child intersected with the other replays the other's intervals.
 * @details When one child falls back to a full-width scan, the intersection is
 *   just the other child's intervals (replayed from the buffer already collected).
 *   Pins the equivalence in both orientations and the both-fall-back full-scan case.
 */
inline void test_intersection_full_width_child_replays_other() {
  using P = std::pair<float, float>;
  using MockI = sdf_interval_detail::MockIntervalShape;
  using MockF = sdf_interval_detail::MockFullWidthShape;
  std::vector<P> ivs = {{60.0f, 80.0f},
                        {20.0f, 40.0f}}; // multi, emission order
  MockI shape{&ivs};
  MockF full;

  {
    SDF::Intersection<MockF, MockI> s(full, shape);
    std::vector<P> out;
    bool ok = s.get_horizontal_intervals<256, 128>(
        0, [&](float st, float en) { out.push_back({st, en}); });
    HS_EXPECT_TRUE(ok);
    HS_EXPECT_SIZE_OR_RETURN(out, 2);
    HS_EXPECT_NEAR(out[0].first, 60.0f, 1e-4f);
    HS_EXPECT_NEAR(out[1].first, 20.0f, 1e-4f);
  }

  // Symmetric: A has intervals, B full-width.
  {
    SDF::Intersection<MockI, MockF> s(shape, full);
    std::vector<P> out;
    bool ok = s.get_horizontal_intervals<256, 128>(
        0, [&](float st, float en) { out.push_back({st, en}); });
    HS_EXPECT_TRUE(ok);
    HS_EXPECT_SIZE_OR_RETURN(out, 2);
    HS_EXPECT_NEAR(out[0].first, 60.0f, 1e-4f);
    HS_EXPECT_NEAR(out[1].first, 20.0f, 1e-4f);
  }

  // Both fall back -> full-scan fallback (return false).
  {
    SDF::Intersection<MockF, MockF> s(full, full);
    std::vector<P> out;
    bool ok = s.get_horizontal_intervals<256, 128>(
        0, [&](float st, float en) { out.push_back({st, en}); });
    HS_EXPECT_FALSE(ok);
    HS_EXPECT_EQ(out.size(), static_cast<size_t>(0));
  }
}

/**
 * @brief Verifies Intersection emits nothing whenever it requests a full-row scan.
 * @details The interval protocol pairs a false return with an empty emission.
 *   Driven by a child that itself violates the protocol.
 */
inline void test_intersection_full_scan_emits_no_spans() {
  using P = std::pair<float, float>;
  using MockE = sdf_interval_detail::MockEmitThenFullWidthShape;
  using MockF = sdf_interval_detail::MockFullWidthShape;
  MockE emitting;
  MockF full;

  {
    SDF::Intersection<MockF, MockE> s(full, emitting);
    std::vector<P> out;
    bool ok = s.get_horizontal_intervals<256, 128>(
        0, [&](float st, float en) { out.push_back({st, en}); });
    HS_EXPECT_FALSE(ok);
    HS_EXPECT_EQ(out.size(), static_cast<size_t>(0));
  }

  {
    SDF::Intersection<MockE, MockF> s(emitting, full);
    std::vector<P> out;
    bool ok = s.get_horizontal_intervals<256, 128>(
        0, [&](float st, float en) { out.push_back({st, en}); });
    HS_EXPECT_FALSE(ok);
    HS_EXPECT_EQ(out.size(), static_cast<size_t>(0));
  }
}

/**
 * @brief Verifies a seam-straddling band shared by both children is intersected across wrap frames.
 * @details A emits the seam band as [-10, 10]; B emits the SAME physical band as
 *   [W-10, W+10]. After normalizing both into [0, W) the bands coincide and the
 *   intersection is the full shared band, split at the seam.
 */
inline void test_intersection_seam_straddle_overlaps_across_wrap_frames() {
  using P = std::pair<float, float>;
  using Mock = sdf_interval_detail::MockIntervalShape;
  std::vector<P> a_ivs = {{-10.0f, 10.0f}};  // seam band, negative frame
  std::vector<P> b_ivs = {{246.0f, 266.0f}}; // same band, [W,2W) frame (W=256)
  Mock A{&a_ivs}, B{&b_ivs};
  SDF::Intersection<Mock, Mock> s(A, B);

  std::vector<P> out;
  bool ok = s.get_horizontal_intervals<256, 128>(
      0, [&](float st, float en) { out.push_back({st, en}); });
  HS_EXPECT_TRUE(ok);
  HS_EXPECT_SIZE_OR_RETURN(out, 2);
  HS_EXPECT_NEAR(out[0].first, 0.0f, 1e-5f);
  HS_EXPECT_NEAR(out[0].second, 10.0f, 1e-5f);
  HS_EXPECT_NEAR(out[1].first, 246.0f, 1e-5f);
  HS_EXPECT_NEAR(out[1].second, 256.0f, 1e-5f);
}

// ============================================================================
// SmoothUnion — blends at the boundary
// ============================================================================

/** @brief Verifies that away from the blend zone, SmoothUnion's distance equals the hard Union's. */
inline void test_smooth_union_matches_union_far_from_boundary() {
  SDF::Line la(math::Vector(1, 0, 0), math::Vector(0, 0, 1), 0.1f);
  SDF::Line lb(math::Vector(-1, 0, 0), math::Vector(0, 0, -1), 0.1f);
  SDF::Union<SDF::Line, SDF::Line> u(la, lb);
  SDF::SmoothUnion<SDF::Line, SDF::Line> su(la, lb, /*k*/ 0.05f);

  math::Vector mid_a =
      ((math::Vector(1, 0, 0) + math::Vector(0, 0, 1)) * 0.5f).normalized();
  auto r_hard = SDF::distance_of(u, mid_a);
  auto r_soft = SDF::distance_of(su, mid_a);
  HS_EXPECT_NEAR(r_hard.dist, r_soft.dist, 1e-3f);
}

/**
 * @brief Verifies the cubic smin blend term inside the blend band.
 * @details Where |dA - dB| < k the smooth distance dips below min(dA, dB) by
 *   the expected m; outside the band it collapses to the hard min.
 */
inline void test_smooth_union_blends_inside_band() {
  SDF::Line la(math::Vector(1, 0, 0), math::Vector(0, 0, 1), 0.1f);
  SDF::Line lb(math::Vector(-1, 0, 0), math::Vector(0, 0, -1), 0.1f);
  const float k = 0.5f;
  SDF::SmoothUnion<SDF::Line, SDF::Line> su(la, lb, k);

  // Both arcs lie in the y=0 plane, so the north pole is equidistant from both
  // (|dA - dB| ≈ 0), maximizing the cubic blend (h == 1, m == k/6).
  math::Vector p(0, 1, 0);
  float dA = SDF::distance_of(la, p).dist;
  float dB = SDF::distance_of(lb, p).dist;
  HS_EXPECT_TRUE(std::abs(dA - dB) < k);

  float h = std::max(k - std::abs(dA - dB), 0.0f) / k;
  float m = h * h * h * k * (1.0f / 6.0f);
  float soft = SDF::distance_of(su, p).dist;
  HS_EXPECT_NEAR(soft, std::min(dA, dB) - m, 1e-4f);
  HS_EXPECT_LT(soft, std::min(dA, dB) - 1e-4f);
  // Independent of the cubic: the smooth min never rises above the hard min,
  // and k/6 is the deepest it can dip below it anywhere.
  HS_EXPECT_LE(soft, std::min(dA, dB) + 1e-5f);
  HS_EXPECT_GE(soft, std::min(dA, dB) - k * (1.0f / 6.0f) - 1e-5f);

  // Outside the band the blend vanishes and collapses to the hard min.
  math::Vector q =
      ((math::Vector(1, 0, 0) + math::Vector(0, 0, 1)) * 0.5f).normalized();
  float qA = SDF::distance_of(la, q).dist;
  float qB = SDF::distance_of(lb, q).dist;
  HS_EXPECT_TRUE(std::abs(qA - qB) >= k);
  HS_EXPECT_NEAR(SDF::distance_of(su, q).dist, std::min(qA, qB), 1e-5f);
}

/**
 * @brief Verifies SmoothUnion solidity follows its children.
 * @details Children must share solidity (enforced by static_assert): two strokes
 *   take the soft falloff path, two solids the hard 1-px silhouette path.
 */
inline void test_smooth_union_solidity_follows_children() {
  static_assert(!SDF::SmoothUnion<SDF::Line, SDF::Line>::is_solid,
                "two strokes -> falloff path");
  static_assert(
      SDF::SmoothUnion<SDF::PlanarPolygon, SDF::PlanarPolygon>::is_solid,
      "two solids -> silhouette path");
}

/**
 * @brief Verifies the blendability trait tracks which shapes clamp to the far
 *   sentinel, and that every combinator blends only when each child does --
 *   including Subtract, whose is_solid tracks the minuend alone.
 * @details Ring, DistortedRing, FlatDistortedRing and Face report dist = 100
 *   outside their reject band, so both children read the sentinel across the
 *   weld and SmoothUnion collapses to Union; SmoothUnion static_asserts the
 *   trait to reject those instantiations at compile time.
 */
inline void test_sentinel_clampers_are_not_blendable() {
  static_assert(!SDF::blends_smoothly<SDF::Ring>);
  static_assert(!SDF::blends_smoothly<SDF::DistortedRing>);
  static_assert(!SDF::blends_smoothly<SDF::FlatDistortedRing>);
  static_assert(!SDF::blends_smoothly<SDF::Face>);
  static_assert(SDF::blends_smoothly<SDF::Line>);
  static_assert(SDF::blends_smoothly<SDF::PlanarPolygon>);
  static_assert(SDF::blends_smoothly<SDF::SphericalPolygon>);
  static_assert(SDF::blends_smoothly<SDF::Star>);
  static_assert(SDF::blends_smoothly<SDF::Flower>);
  static_assert(!SDF::blends_smoothly<SDF::Union<SDF::Line, SDF::Ring>>);
  static_assert(SDF::blends_smoothly<SDF::Union<SDF::Line, SDF::Line>>);
  static_assert(!SDF::blends_smoothly<SDF::AngularRepeat<SDF::Ring>>);
  static_assert(SDF::blends_smoothly<SDF::AngularRepeat<SDF::Line>>);
  static_assert(!SDF::blends_smoothly<SDF::SmoothUnion<SDF::Line, SDF::Ring>>);
  static_assert(SDF::blends_smoothly<SDF::SmoothUnion<SDF::Line, SDF::Line>>);
  static_assert(!SDF::blends_smoothly<SDF::Intersection<SDF::Line, SDF::Ring>>);
  static_assert(SDF::blends_smoothly<SDF::Intersection<SDF::Line, SDF::Line>>);
  static_assert(!SDF::blends_smoothly<SDF::Subtract<SDF::Line, SDF::Ring>>);
  static_assert(!SDF::blends_smoothly<SDF::Subtract<SDF::Ring, SDF::Line>>);
  static_assert(SDF::blends_smoothly<SDF::Subtract<SDF::Line, SDF::Line>>);
}

/**
 * @brief Verifies the CSG combinators reject a temporary child.
 * @details Every combinator holds its children by reference.
 */
inline void test_csg_combinators_reject_temporary_children() {
  using L = SDF::Line;
  using Poly = SDF::PlanarPolygon;
  static_assert(std::is_constructible_v<SDF::Union<L, L>, L &, L &>);
  static_assert(!std::is_constructible_v<SDF::Union<L, L>, L &&, L &>);
  static_assert(!std::is_constructible_v<SDF::Union<L, L>, L &, L &&>);
  static_assert(!std::is_constructible_v<SDF::Union<L, L>, L &&, L &&>);

  static_assert(
      std::is_constructible_v<SDF::SmoothUnion<L, L>, L &, L &, float>);
  static_assert(
      !std::is_constructible_v<SDF::SmoothUnion<L, L>, L &&, L &, float>);
  static_assert(
      !std::is_constructible_v<SDF::SmoothUnion<L, L>, L &, L &&, float>);

  static_assert(std::is_constructible_v<SDF::Subtract<Poly, L>, Poly &, L &>);
  static_assert(!std::is_constructible_v<SDF::Subtract<Poly, L>, Poly &&, L &>);
  static_assert(!std::is_constructible_v<SDF::Subtract<Poly, L>, Poly &, L &&>);

  static_assert(std::is_constructible_v<SDF::Intersection<L, L>, L &, L &>);
  static_assert(!std::is_constructible_v<SDF::Intersection<L, L>, L &&, L &>);
  static_assert(!std::is_constructible_v<SDF::Intersection<L, L>, L &, L &&>);

  // AngularRepeat copies its axis, so only the child rejects a temporary.
  static_assert(std::is_constructible_v<SDF::AngularRepeat<L>, L &, int>);
  static_assert(!std::is_constructible_v<SDF::AngularRepeat<L>, L &&, int>);
  static_assert(std::is_constructible_v<SDF::AngularRepeat<L>, L &, int,
                                        math::Vector &&>);
  static_assert(!std::is_constructible_v<SDF::AngularRepeat<L>, L &&, int,
                                         math::Vector &&>);
}

/**
 * @brief Verifies Union coalesces two overlapping child intervals into one span.
 * @details The children emit overlapping bands (A [0,40], B [30,70]) that must
 *   collapse to a single [0,70].
 */
inline void test_union_merges_overlapping_intervals() {
  using P = std::pair<float, float>;
  using Mock = sdf_interval_detail::MockIntervalShape;
  std::vector<P> a_ivs = {{0.0f, 40.0f}};
  std::vector<P> b_ivs = {{30.0f, 70.0f}};
  Mock A{&a_ivs}, B{&b_ivs};
  SDF::Union<Mock, Mock> u(A, B);

  std::vector<P> out;
  bool ok = u.get_horizontal_intervals<256, 128>(
      0, [&](float st, float en) { out.push_back({st, en}); });
  HS_EXPECT_TRUE(ok);
  HS_EXPECT_SIZE_OR_RETURN(out, 1);
  HS_EXPECT_NEAR(out[0].first, 0.0f, 1e-4f);
  HS_EXPECT_NEAR(out[0].second, 70.0f, 1e-4f);
}

/**
 * @brief Verifies Union welds a seam-straddling span with an overlapping one.
 * @details A [-10,6] straddles θ=0 in the negative frame and overlaps B's
 *   in-frame [2,12]; the overlap must coalesce to a single [-10,12] span.
 */
inline void test_union_seam_straddle_merges_overlapping_intervals() {
  using P = std::pair<float, float>;
  using Mock = sdf_interval_detail::MockIntervalShape;
  std::vector<P> a_ivs = {{-10.0f, 6.0f}};
  std::vector<P> b_ivs = {{2.0f, 12.0f}};
  Mock A{&a_ivs}, B{&b_ivs};
  SDF::Union<Mock, Mock> u(A, B);

  std::vector<P> out;
  bool ok = u.get_horizontal_intervals<256, 128>(
      0, [&](float st, float en) { out.push_back({st, en}); });
  HS_EXPECT_TRUE(ok);
  HS_EXPECT_SIZE_OR_RETURN(out, 1);
  HS_EXPECT_NEAR(out[0].first, -10.0f, 1e-4f);
  HS_EXPECT_NEAR(out[0].second, 12.0f, 1e-4f);
}

/**
 * @brief Verifies three- and four-way nested Unions of real leaves compile and
 *        emit every child arc.
 * @details Nesting depth is gated by sdf_max_spans. Four coaxial rings at
 *   disjoint radii cross an equatorial row in 8 disjoint spans, which is also
 *   the bound the trait reports.
 */
inline void test_nested_union_emits_every_child_arc() {
  using P = std::pair<float, float>;
  constexpr int W = 256, H = 128;
  // Ring axis along +X so the row math has a non-degenerate horizontal
  // projection (a +Y axis full-row scans instead).
  const math::Basis b{math::Vector(0, 1, 0), math::Vector(1, 0, 0),
                      math::Vector(0, 0, 1)};
  SDF::Ring r1(b, 0.4f, 0.05f), r2(b, 0.6f, 0.05f);
  SDF::Ring r3(b, 0.8f, 0.05f), r4(b, 1.0f, 0.05f);

  using U2 = SDF::Union<SDF::Ring, SDF::Ring>;
  using U3 = SDF::Union<U2, SDF::Ring>;
  using U4 = SDF::Union<U2, U2>;
  static_assert(SDF::sdf_max_spans<SDF::Ring>::value == 2);
  static_assert(SDF::sdf_max_spans<U3>::value == 6);
  static_assert(SDF::sdf_max_spans<U4>::value == 8);
  // Under Intersection, each child must fit one IntervalBuffer.
  static_assert(sizeof(SDF::Intersection<U2, U2>) > 0);

  U2 u_lo(r1, r2), u_hi(r3, r4);
  U3 u3(u_lo, r3);
  U4 u4(u_lo, u_hi);

  std::vector<P> out3, out4;
  bool ok3 = u3.get_horizontal_intervals<W, H>(
      H / 2, [&](float st, float en) { out3.push_back({st, en}); });
  bool ok4 = u4.get_horizontal_intervals<W, H>(
      H / 2, [&](float st, float en) { out4.push_back({st, en}); });
  HS_EXPECT_TRUE(ok3);
  HS_EXPECT_TRUE(ok4);

  HS_EXPECT_EQ(out3.size(), static_cast<size_t>(6));
  HS_EXPECT_EQ(out4.size(), static_cast<size_t>(8));
  for (size_t i = 1; i < out4.size(); ++i)
    HS_EXPECT_TRUE(out4[i - 1].second < out4[i].first);
}

/**
 * @brief Verifies SmoothUnion's k-padded union welds a seam-straddling span
 *        with an overlapping one.
 * @details Each child interval is inflated by pad_px = k·W/(2π·sinφ) before the
 *   merge; a band straddling θ=0 and an overlapping in-frame band must coalesce
 *   to one padded span.
 */
inline void test_smooth_union_seam_straddle_merges_padded_intervals() {
  using P = std::pair<float, float>;
  using Mock = sdf_interval_detail::MockIntervalShape;
  constexpr int W = 256, H = 128;
  math::init_geometry_luts<
      W,
      H>();              // fill sin_phi; scan_region does this in production
  const int row = H / 2; // equatorial row: sinφ ≈ 1, pad ≈ k·W/(2π)
  const float k = 0.02f;
  const float sin_phi = math::TrigLUT<W, H>::sin_phi[row];
  const float pad =
      std::min(k * W / (2.0f * math::PI_F) / sin_phi, static_cast<float>(W));
  std::vector<P> a_ivs = {{-10.0f, 6.0f}};
  std::vector<P> b_ivs = {{2.0f, 12.0f}};
  Mock A{&a_ivs}, B{&b_ivs};
  SDF::SmoothUnion<Mock, Mock> su(A, B, k);

  std::vector<P> out;
  bool ok = su.get_horizontal_intervals<W, H>(
      row, [&](float st, float en) { out.push_back({st, en}); });
  HS_EXPECT_TRUE(ok);
  HS_EXPECT_SIZE_OR_RETURN(out, 1);
  HS_EXPECT_NEAR(out[0].first, -10.0f - pad, 1e-4f);
  HS_EXPECT_NEAR(out[0].second, 12.0f + pad, 1e-4f);
  // Independent of the pad formula: the weld must cover both raw spans and
  // widen past their union in both directions.
  HS_EXPECT_LT(out[0].first, -10.0f);
  HS_EXPECT_GT(out[0].second, 12.0f);
}

/**
 * @brief Verifies the weld pad widens toward the poles by the 1/sinφ factor.
 * @details A point interval's emitted span is exactly twice the row pad, so a
 *   near-pole row (small sinφ) must yield a strictly wider span than the equator.
 */
inline void test_smooth_union_pad_widens_toward_pole() {
  using P = std::pair<float, float>;
  using Mock = sdf_interval_detail::MockIntervalShape;
  constexpr int W = 256, H = 128;
  math::init_geometry_luts<
      W,
      H>(); // fill sin_phi; scan_region does this in production
  const float k = 0.05f;
  std::vector<P> ivs = {{100.0f, 100.0f}}; // a point; only the pad sets width
  Mock A{&ivs}, B{&ivs};
  SDF::SmoothUnion<Mock, Mock> su(A, B, k);

  auto span_at = [&](int row) {
    std::vector<P> out;
    const bool HANDLED = su.get_horizontal_intervals<W, H>(
        row, [&](float st, float en) { out.push_back({st, en}); });
    HS_EXPECT_TRUE(HANDLED);
    if (!HANDLED)
      return -1.0f;
    return out.empty() ? 0.0f : out[0].second - out[0].first;
  };
  const int NEAR_POLE = su.pad_rows<H>() + 4;
  const float POLE_SPAN = span_at(NEAR_POLE);
  HS_EXPECT_GT(POLE_SPAN, span_at(H / 2));
  HS_EXPECT_NEAR(POLE_SPAN,
                 2.0f * fminf(k * W / math::TWO_PI_F /
                                  math::TrigLUT<W, H>::sin_phi[NEAR_POLE],
                              static_cast<float>(W)),
                 1e-3f);
}

// ============================================================================
// AngularRepeat — folds azimuth into N sectors
// ============================================================================

/** @brief Verifies AngularRepeat agrees with the base shape in the canonical (zero-angle) sector. */
inline void test_angular_repeat_matches_base_at_zero_angle() {
  SDF::Line ln(math::Vector(1, 0, 0), math::Vector(0.7071f, 0, 0.7071f), 0.05f);
  SDF::AngularRepeat<SDF::Line> rep(ln, /*reps*/ 4, math::Vector(0, 1, 0));

  math::Vector mid =
      ((math::Vector(1, 0, 0) + math::Vector(0.7071f, 0, 0.7071f)) * 0.5f)
          .normalized();
  auto r_base = SDF::distance_of(ln, mid);
  auto r_rep = SDF::distance_of(rep, mid);
  HS_EXPECT_TRUE(r_rep.dist < 0.0f);
  HS_EXPECT_NEAR(r_rep.dist, r_base.dist, 1e-3f);
}

/** @brief Verifies AngularRepeat folds a line in the canonical sector into a folded copy. */
inline void test_angular_repeat_creates_copies() {
  SDF::Line ln(math::Vector(1, 0, 0), math::Vector(0.7071f, 0, 0.7071f), 0.05f);
  SDF::AngularRepeat<SDF::Line> rep(ln, 4, math::Vector(0, 1, 0));

  // Rotate the midpoint by one sector (90° around Y) onto a folded copy.
  math::Vector mid =
      ((math::Vector(1, 0, 0) + math::Vector(0.7071f, 0, 0.7071f)) * 0.5f)
          .normalized();
  math::Quaternion q90 =
      math::make_rotation(math::Vector(0, 1, 0), math::PI_F * 0.5f);
  math::Vector mid_rot = math::rotate(mid, q90);

  auto r = SDF::distance_of(rep, mid_rot);
  HS_EXPECT_TRUE(r.dist < 0.0f);
}

/**
 * @brief Verifies AngularRepeat's child UV (t) is sector-local, not global.
 * @details distance() folds p into one sector before evaluating the child, so
 *   a point and its copy one full sector away share a t.
 */
inline void test_angular_repeat_t_is_sector_local() {
  math::Basis b = equator_basis();
  SDF::Ring ring(b, 1.0f, 0.1f);
  const int reps = 4;
  SDF::AngularRepeat<SDF::Ring> rep(ring, reps, math::Vector(0, 1, 0));

  // A point at 30° azimuth (inside sector 0, off a boundary) and its copy one
  // sector away.
  float az = math::PI_F / 6.0f;
  math::Vector p(cosf(az), 0.0f, sinf(az));
  math::Quaternion q_sector =
      math::make_rotation(math::Vector(0, 1, 0), 2 * math::PI_F / reps);
  math::Vector p2 = math::rotate(p, q_sector);

  // The un-repeated ring sees two global azimuths one sector apart (sign is
  // handedness-dependent, so compare the wrapped gap).
  float t_global_1 = SDF::distance_of(ring, p).t;
  float t_global_2 = SDF::distance_of(ring, p2).t;
  float dg = t_global_2 - t_global_1;
  dg -= floorf(dg);
  HS_EXPECT_NEAR(std::min(dg, 1.0f - dg), 1.0f / reps, 1e-3f);

  // The repeated shape folds both into the same sector → identical sector-local t.
  float t_rep_1 = SDF::distance_of(rep, p).t;
  float t_rep_2 = SDF::distance_of(rep, p2).t;
  HS_EXPECT_NEAR(t_rep_1, t_rep_2, 1e-3f);
  HS_EXPECT_NEAR(t_rep_1, t_global_1, 1e-3f);
}
