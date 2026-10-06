/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/render/plot/cull.h.

/**
 * @brief True when @p P statically declares it has no world cull stage, so a
 *        cull predicate may be evaluated against the raw geometry.
 * @tparam P Pipeline type; types without the has_world_cull member are
 *           conservatively not hoistable.
 */
template <typename P> static consteval bool pipeline_hoistable_cull() {
  if constexpr (requires { P::has_world_cull; })
    return !P::has_world_cull;
  else
    return false;
}

/**
 * @brief True when @p P statically declares it has no world-space stage, so a
 *        caller may plot a point through precomputed screen coordinates.
 * @tparam P Pipeline type; types without the has_world_stage member are
 *           conservatively not hoistable.
 */
template <typename P> static consteval bool pipeline_hoistable_projection() {
  if constexpr (requires { P::has_world_stage; })
    return !P::has_world_stage;
  else
    return false;
}

/**
 * @brief Conservative screen-length test: true only when the geodesic edge
 *        a->b provably spans at most SCREEN_STEP_PX on screen.
 * @tparam W,H Rasterization resolution (pixel grid).
 * @param a Edge start (unit sphere point).
 * @param b Edge end (unit sphere point).
 * @details Tightened form of the fast-path test
 * `total_dist <= screen_step(sample(0))` in multiplies only, via
 * sin(theta)*tangent = b - a*cos(theta) and theta/sin(theta) <= F for
 * theta <= base_step (enforced by the chord cap). True also implies
 * theta >= EPS_GEOMETRIC.
 */
HS_O3_BEGIN
template <int W, int H>
static inline bool edge_fits_one_dot(const math::Vector &a,
                                     const math::Vector &b) {
  constexpr float BASE = (2.0f * math::PI_F) / W;
  constexpr float B2 = BASE * BASE;
  static_assert(B2 < 1.0f, "chord/angle bounds assume base_step < 1 rad");
  constexpr float KX2 = (W / (2.0f * math::PI_F)) * (W / (2.0f * math::PI_F));
  const float KY2 = (math::ROWS_PER_RADIAN<H>)*(math::ROWS_PER_RADIAN<H>);
  constexpr float SPX2 = SCREEN_STEP_PX * SCREEN_STEP_PX;
  // Preserve the fast-path implication under screen_rsqrt's <0.1% undershoot.
  constexpr float SCREEN_RSQRT_MIN2 = 0.999f * 0.999f;
  // chord^2 caps: (2 sin(BASE/2))^2 >= B2*(1 - B2/12) bounds theta <= BASE.
  // The lower cap stays above process_segment's EPS_GEOMETRIC degenerate
  // branch.
  constexpr float CHORD2_MAX = B2 * (1.0f - B2 / 12.0f);
  constexpr float CHORD2_MIN = 4.0e-6f;
  // (theta/sin(theta))^2 <= F2 for theta <= BASE, plus float-rounding slack.
  constexpr float F2 = (1.0001f / ((1.0f - B2 / 6.0f) * (1.0f - B2 / 6.0f)));
  const math::Vector d = b - a;
  const float chord2 = math::dot(d, d);
  if (chord2 > CHORD2_MAX || chord2 < CHORD2_MIN)
    return false;
  const float sin2 = 1.0f - a.y * a.y;
  if (sin2 < 1e-7f)
    return false;
  const float c = math::dot(a, b);
  const float cx = a.x * b.z - a.z * b.x;
  const float ty = b.y - c * a.y;
  return F2 * (KX2 * cx * cx + KY2 * ty * ty * sin2) <=
         SPX2 * SCREEN_RSQRT_MIN2 * sin2 * sin2;
}
HS_O3_END
