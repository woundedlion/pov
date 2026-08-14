/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <cmath>
#include "math/geometry.h"
#include "engine/constants.h"
#include "engine/concepts.h"
#include "engine/util.h"
#include "render/sdf_common.h"

/**
 * @file sdf_rings.h
 * @brief The ring leaves: Ring, DistortedRing and FlatDistortedRing.
 */

namespace SDF {

/**
 * @brief Calculates signed distance to a ring.
 * @details Register semantics: the DistanceResult table (stroke row: Ring).
 */
struct Ring {
  const Basis &basis; /**< Orientation frame (v = ring axis); retained by
                         reference, so it must outlive the shape. */
  float radius;       /**< Ring radius as a fraction of the hemisphere. */
  float thickness;    /**< Half-width of the stroke (radians). */
  float phase;        /**< Azimuth phase offset (radians). */

  Vector normal, u, w; /**< Ring axis and the two in-plane basis vectors. */
  float ny;            /**< y-component of the ring axis. */
  float target_angle,
      center_phi; /**< Centerline polar angle and axis colatitude. */
  float cos_max, cos_min, cos_target,
      inv_sin_target; /**< Precomputed band trig. */

  float r_val;       /**< Horizontal projection length of the axis (for full-row
                         check). */
  float alpha_angle; /**< Azimuth angle of the normal vector in the XZ plane. */
  static constexpr bool is_solid = false; /**< Ring renders as a stroke. */

  /**
   * @brief Builds a ring from its basis, radius, thickness, and phase.
   * @param b Orientation frame (v = ring axis); retained by reference, so it
   *          must outlive the shape.
   * @param r Ring radius as a fraction of the hemisphere.
   * @param th Half-width of the stroke (radians).
   * @param ph Azimuth phase offset (radians).
   */
  Ring(const Basis &b, float r, float th, float ph = 0)
      : basis(b), radius(r), thickness(th), phase(ph) {
    normal = basis.v;
    u = basis.u;
    w = basis.w;
    AxisProjection ap = project_axis(normal);
    ny = ap.ny;

    target_angle = radius * (PI_F / 2.0f);
    center_phi = acosf(std::max(-1.0f, std::min(1.0f, ny)));

    float ang_min = std::max(0.0f, target_angle - thickness);
    float ang_max = std::min(PI_F, target_angle + thickness);
    cos_max = cosf(ang_min);
    cos_min = cosf(ang_max);
    cos_target = cosf(target_angle);

    bool safe_approx =
        (target_angle > POLE_SAFE_MARGIN &&
         target_angle < PI_F - POLE_SAFE_MARGIN &&
         thickness < RING_LINEARIZE_TAN_FRAC * std::abs(tanf(target_angle)));
    inv_sin_target = safe_approx ? (1.0f / sinf(target_angle)) : 0.0f;

    r_val = ap.R_val;
    alpha_angle = ap.alpha_angle;
  }

  /**
   * @brief Deleted constructor from a temporary Basis.
   * @details The ring retains its basis by reference, so binding a temporary
   * would leave every later read of basis dangling.
   */
  Ring(const Basis &&, float, float, float = 0) = delete;

  /**
   * @brief Maps the ring's latitude band to its inclusive scanline row range.
   * @tparam H Canvas height in rows.
   * @return Row bounds covering the ring plus its AA falloff.
   */
  template <int H> Bounds get_vertical_bounds() const {
    PhiBand band = clamp_phi_band(center_phi, target_angle);

    // Deliberately tighter than the true band: the trimmed fringe carries
    // coverage up to quintic_kernel(0.05) = 1.2e-3, a hair over the
    // rasterizer's 1e-3 alpha cut.
    float eff_th = 0.95f * thickness;
    float f_phi_min = std::max(0.0f, band.phi_min - eff_th);
    float f_phi_max = std::min(PI_F, band.phi_max + eff_th);

    return phi_bounds_to_rows<H>(f_phi_min, f_phi_max);
  }

  /**
   * @brief Computes the horizontal scanline intervals for this shape at a given
   * y-coordinate.
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @tparam OutputIt Sink type invoked as out(float start, float end).
   * @param y The vertical pixel coordinate (row index).
   * @param out Output iterator or callback accepting (float start, float end).
   * @return True if intervals were found and reported; false requests a full
   * scan.
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
    float cos_phi = TrigLUT<W, H>::cos_phi[y];
    float sin_phi = TrigLUT<W, H>::sin_phi[y];

    if (needs_full_row_scan(sin_phi))
      return false;

    float denom = r_val * sin_phi;
    emit_annular_band<W>(cos_min, cos_max, ny, cos_phi, denom, alpha_angle,
                         out);
    return true;
  }

  /**
   * @brief Whether row interval math is degenerate and the row must be
   *        full-row scanned.
   * @param sin_phi sin of the row's colatitude.
   * @return True when get_horizontal_intervals would return false at this row.
   */
  bool needs_full_row_scan(float sin_phi) const {
    return r_val < MIN_HORIZONTAL_PROJ ||
           std::abs(r_val * sin_phi) < INTERVAL_DENOM_EPS;
  }

  /**
   * @brief Stroke coverage from a precomputed axis dot, skipping the
   *        DistanceResult round trip.
   * @param d dot(p, normal) for the pixel's unit vector.
   * @return quintic stroke alpha in [0, 1]; 0 outside the band. Same distance
   *         branches and float ops as distance<false>() + process_pixel's
   *         stroke epilogue, so coverage is bit-identical to that path.
   */
  __attribute__((always_inline)) float stroke_alpha(float d) const {
    if (d < cos_min || d > cos_max)
      return 0.0f;
    float dist;
    if (inv_sin_target != 0) {
      dist = std::abs(d - cos_target) * inv_sin_target;
    } else {
      float polar = fast_acos(hs::clamp(d, -1.0f, 1.0f));
      dist = std::abs(polar - target_angle);
    }
    float sd = dist - thickness;
    if (sd >= 0.0f || thickness <= 0.0f)
      return 0.0f;
    return quintic_kernel(-sd / thickness);
  }

  /**
   * @brief Computes signed distance to the ring, writing into res.
   * @tparam ComputeUVs When true, also computes the azimuthal t parameter.
   * @param p Point on sphere (normalized).
   * @param res Output result; dist = signed distance, raw_dist = unsigned
   *        centerline distance, t = azimuth in [0,1) when ComputeUVs.
   */
  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    float d = dot(p, normal);
    if (d < cos_min || d > cos_max) {
      res = DistanceResult(FAR_SENTINEL, 0.0f, FAR_SENTINEL, 0.0f, thickness);
      return;
    }

    float dist = 0;
    if (inv_sin_target != 0) {
      dist = std::abs(d - cos_target) * inv_sin_target;
    } else {
      float polar = fast_acos(hs::clamp(d, -1.0f, 1.0f));
      dist = std::abs(polar - target_angle);
    }

    float t = 0.0f;
    if constexpr (ComputeUVs) {
      t = wrap_t(basis_azimuth(p, u, w, phase) / TWO_PI_F);
    }

    res = DistanceResult(dist - thickness, t, dist, 0.0f, thickness);
  }
};

/**
 * @brief Per-azimuth-chunk knot ranges backing DistortedRing's knot prefilter.
 * @details Caller-owned, so the callback and undisplaced modes carry none of
 * it. One instance per knot-mode ring; must outlive the shape.
 */
struct KnotPrefilter {
  static constexpr int CHUNKS = 32; /**< Azimuth chunks in the prefilter. */
  float lo[CHUNKS];                 /**< Min knot per azimuth chunk. */
  float hi[CHUNKS];                 /**< Max knot per azimuth chunk. */
};

/**
 * @brief Calculates signed distance to a distorted ring.
 * @details Register semantics: the DistanceResult table (stroke row:
 * DistortedRing).
 */
struct DistortedRing {
  const Basis &basis; /**< Orientation frame (v = ring axis); retained by
                         reference, so it must outlive the shape. */
  float radius;       /**< Ring radius as a fraction of the hemisphere. */
  float thickness;    /**< Half-width of the stroke (radians). */
  ScalarFn shift_fn;  /**< Per-azimuth centerline shift, t in [0,1) -> radians;
                         empty in knot mode. */
  const float *knots =
      nullptr;   /**< Optional lut_n + 1 shift knots (entry lut_n repeats entry
                    0); selects exact polyline distance. */
  int lut_n = 0; /**< Knot cell count when knots is set. */
  const KnotPrefilter *prefilter =
      nullptr;          /**< Caller-owned chunk ranges over knots. */
  float max_distortion; /**< Maximum magnitude of the shift (radians). */
  float phase;          /**< Azimuth phase offset (radians). */

  Vector normal, u, w; /**< Ring axis and the two in-plane basis vectors. */
  float ny;            /**< y-component of the ring axis. */
  float target_angle,
      center_phi;      /**< Centerline polar angle and axis colatitude. */
  float max_thickness; /**< thickness + max_distortion (radians). */

  float r_val;       /**< Horizontal projection length of the axis. */
  float alpha_angle; /**< Azimuth of the normal in the XZ plane. */
  float cos_max_limit, cos_min_limit; /**< Cosines of the widened band edges. */
  bool suppress_pole_fill =
      false; /**< Drop the degenerate exact-pole row rather than full-row
                filling it (see get_horizontal_intervals). */
  static constexpr bool is_solid =
      false; /**< Distorted ring renders as a stroke. */

  /**
   * @brief Builds a distorted ring with a per-azimuth centerline shift.
   * @param b Orientation frame (v = ring axis); retained by reference, so it
   *          must outlive the shape.
   * @param r Ring radius as a fraction of the hemisphere.
   * @param th Half-width of the stroke (radians).
   * @param sf Per-azimuth centerline shift function, t in [0,1) -> radians.
   * @param md Maximum magnitude of sf over t in [0,1) (radians). PRECONDITION:
   *           md must be a true upper bound on |sf|. It widens the reject
   * bands, so an underestimate silently culls genuine arcs. Pinned by the cull
   *           tests, not checked here.
   * @param ph Azimuth phase offset (radians).
   */
  DistortedRing(const Basis &b, float r, float th, ScalarFn sf, float md,
                float ph)
      : basis(b), radius(r), thickness(th), shift_fn(sf), max_distortion(md),
        phase(ph) {
    normal = basis.v;
    u = basis.u;
    w = basis.w;
    AxisProjection ap = project_axis(normal);
    ny = ap.ny;
    target_angle = radius * (PI_F / 2.0f);
    center_phi = acosf(std::max(-1.0f, std::min(1.0f, ny)));
    max_thickness = thickness + max_distortion;

    r_val = ap.R_val;
    alpha_angle = ap.alpha_angle;

    float ang_min = std::max(0.0f, target_angle - max_thickness);
    float ang_max = std::min(PI_F, target_angle + max_thickness);
    cos_max_limit = cosf(ang_min);
    cos_min_limit = cosf(ang_max);
  }

  /**
   * @brief Deleted constructor from a temporary Basis.
   * @details The ring retains its basis by reference, so binding a temporary
   * would leave every later read of basis dangling.
   */
  DistortedRing(const Basis &&, float, float, ScalarFn, float, float) = delete;

  /**
   * @brief Builds a distorted ring whose centerline is a shift-knot polyline.
   * @param b Orientation frame (v = ring axis); retained by reference, so it
   *          must outlive the shape.
   * @param r Ring radius as a fraction of the hemisphere.
   * @param th Half-width of the stroke (radians).
   * @param kn n + 1 centerline shifts (radians), one per equal azimuth cell,
   *           entry n repeating entry 0; must outlive the shape. distance()
   *           returns the exact distance to this polyline (within the local
   *           tangent chart), so steep or sharply curved segments render at
   *           full stroke width with no slope approximation.
   * @param n Number of knot cells; at least 1.
   * @param ph Azimuth phase offset (radians).
   * @param pf Prefilter storage filled here; must outlive the shape.
   */
  DistortedRing(const Basis &b, float r, float th, const float *kn, int n,
                float ph, KnotPrefilter &pf)
      : DistortedRing(b, r, th, ScalarFn{}, 0.0f, ph) {
    HS_CHECK(kn != nullptr && n >= 1); // n == 0: knot_v(-1) reads knots[-1]
    knots = kn;
    lut_n = n;
    prefilter = &pf;
    float min_shift = kn[0];
    float max_shift = kn[0];
    // Per-chunk knot ranges for the per-pixel prefilter. A segment registers
    // its endpoints in every chunk its azimuth extent touches (a straddling
    // segment spans two), so the polyline inside a chunk never leaves
    // [pf.lo, pf.hi].
    for (int c = 0; c < KnotPrefilter::CHUNKS; ++c) {
      pf.lo[c] = 1e9f;
      pf.hi[c] = -1e9f;
    }
    for (int k = 0; k < n; ++k) {
      float lo = std::min(kn[k], kn[k + 1]);
      float hi = std::max(kn[k], kn[k + 1]);
      min_shift = std::min(min_shift, lo);
      max_shift = std::max(max_shift, hi);
      int c1 = k * KnotPrefilter::CHUNKS / n;
      int c2 = std::min((k + 1) * KnotPrefilter::CHUNKS / n,
                        KnotPrefilter::CHUNKS - 1);
      for (int c = c1; c <= c2; ++c) {
        pf.lo[c] = std::min(pf.lo[c], lo);
        pf.hi[c] = std::max(pf.hi[c], hi);
      }
    }

    max_distortion = std::max(std::abs(min_shift), std::abs(max_shift));
    max_thickness = thickness + max_distortion;
    float ang_min = hs::clamp(target_angle + min_shift - thickness, 0.0f, PI_F);
    float ang_max = hs::clamp(target_angle + max_shift + thickness, 0.0f, PI_F);
    cos_max_limit = cosf(ang_min);
    cos_min_limit = cosf(ang_max);
  }

  /**
   * @brief Deleted constructor from a temporary Basis.
   * @details The ring retains its basis by reference, so binding a temporary
   * would leave every later read of basis dangling.
   */
  DistortedRing(const Basis &&, float, float, const float *, int, float,
                KnotPrefilter &) = delete;

  /**
   * @brief Maps the distorted ring's widened latitude band to its row range.
   * @tparam H Canvas height in rows.
   * @return Inclusive row bounds covering the ring plus distortion margin.
   */
  template <int H> Bounds get_vertical_bounds() const {
    PhiBand band = clamp_phi_band(center_phi, target_angle);

    float margin = max_thickness + BOUNDS_MARGIN_WIDE;
    float f_phi_min = std::max(0.0f, band.phi_min - margin);
    float f_phi_max = std::min(PI_F, band.phi_max + margin);

    return phi_bounds_to_rows<H>(f_phi_min, f_phi_max);
  }

  /**
   * @brief Emits the widened annular-band intervals for one scanline row.
   * @tparam W Canvas width in columns.
   * @tparam H Canvas height in rows.
   * @tparam OutputIt Sink type invoked as out(float start, float end).
   * @param y The row index.
   * @param out Sink accepting (float start, float end).
   * @return True if the row was handled; false requests a full scan.
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
    float cos_phi = TrigLUT<W, H>::cos_phi[y];
    float sin_phi = TrigLUT<W, H>::sin_phi[y];

    if (r_val < MIN_HORIZONTAL_PROJ)
      return false;

    float denom = r_val * sin_phi;
    if (std::abs(denom) < INTERVAL_DENOM_EPS)
      // r_val cleared MIN_HORIZONTAL_PROJ above, so a vanishing denom is the
      // exact pole row: every column aliases to the one pole point. A displaced
      // stroke reaches the pole at a single azimuth, but the degenerate row
      // math can't recover which column, so the default (return false) full-row
      // scans and fills the whole aliased row -- fine for a lone ring, but a
      // dense ring stack renders it as a solid pole cap. suppress_pole_fill
      // drops the row.
      return suppress_pole_fill;

    emit_annular_band<W>(cos_min_limit, cos_max_limit, ny, cos_phi, denom,
                         alpha_angle, out);
    return true;
  }

  /**
   * @brief Computes signed distance to the distorted ring, writing into res.
   * @tparam ComputeUVs When true, also computes the azimuthal t parameter.
   * @param p Point on sphere (normalized).
   * @param res Output result; dist = signed distance minus thickness, raw_dist
   * = unsigned centerline distance, t = azimuth in [0,1) when ComputeUVs.
   */
  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    float d = dot(p, normal);
    // Early reject: outside bounding annulus
    if (d < cos_min_limit || d > cos_max_limit) {
      res = DistanceResult(FAR_SENTINEL, 0.0f, FAR_SENTINEL, 0.0f, thickness);
      return;
    }
    float polar = fast_acos(hs::clamp(d, -1.0f, 1.0f));

    float t_norm = wrap_t(basis_azimuth(p, u, w, phase) / TWO_PI_F);

    float dist;
    if (knots)
      dist = polyline_distance(t_norm, polar,
                               sqrtf(std::max(1.0f - d * d, POLE_SIN2_FLOOR)));
    else
      dist = std::abs(polar - (target_angle + shift_fn(t_norm)));

    if constexpr (!ComputeUVs)
      t_norm = 0.0f;

    res = DistanceResult(dist - thickness, t_norm, dist, 0.0f, thickness);
  }

  /**
   * @brief distance<true>() from a precomputed pixel frame.
   * @param d Pixel dot ring axis (= dot(p, normal)).
   * @param polar fast_acos(clamp(d, -1, 1)).
   * @param sin_polar sqrtf(max(1 - d * d, POLE_SIN2_FLOOR)).
   * @param t_norm Pixel azimuth in [0, 1), phase applied.
   * @param res Output result, identical to distance<true>() except for
   *        undisplaced knot rings, which take the exact polar distance (the
   *        zero-knot polyline agrees only to within an ulp).
   * @details Same-axis ring stacks share d/polar/sin_polar/t_norm across every
   * ring at a pixel; hoisting them there drops the per-ring dot/acos/atan2
   * recompute.
   */
  HS_O3_FN void distance_from_frame(float d, float polar, float sin_polar,
                                    float t_norm, DistanceResult &res) const {
    if (d < cos_min_limit || d > cos_max_limit) {
      res = DistanceResult(FAR_SENTINEL, 0.0f, FAR_SENTINEL, 0.0f, thickness);
      return;
    }
    float dist;
    if (knots)
      dist = max_distortion > 0.0f ? polyline_distance(t_norm, polar, sin_polar)
                                   : std::abs(polar - target_angle);
    else
      dist = std::abs(polar - (target_angle + shift_fn(t_norm)));
    res = DistanceResult(dist - thickness, t_norm, dist, 0.0f, thickness);
  }

  static constexpr float POLE_SIN2_FLOOR =
      1e-6f; /**< sin^2 floor keeping the chart's azimuth scale positive at the
                poles. */

private:
  static constexpr int MAX_SEARCH_CELLS =
      64; /**< Outward search budget per side; only near-pole chart compression
             approaches it. */

  /**
   * @brief Distance from a pixel to the knot polyline.
   * @param t_norm Pixel azimuth in [0, 1) (phase applied).
   * @param polar Pixel polar angle (radians).
   * @param sin_polar sqrtf(max(1 - d * d, POLE_SIN2_FLOOR)) for the pixel;
   *        hoisted to the caller so a ring stack pays it once per pixel.
   * @return Geodesic distance to the nearest polyline point (radians), exact
   *         within the local tangent chart for distances up to `thickness`, or
   *         a lower bound when the outward search hits its cell budget; past
   *         that reach either a lower bound above `thickness` (prefilter) or
   *         FAR_SENTINEL, never the reach itself.
   * @details Works in the chart (azimuth * sin(polar), polar) centered on the
   * pixel: exact point-to-segment distances, searched outward from the
   * pixel's own cell. A segment o cells away is at least (o - 1) * cell_u
   * away, so the search stops once that gap exceeds the best distance found
   * (or the stroke reach, past which alpha is zero regardless). True distance
   * is 1-Lipschitz in screen position, so the stroke's alpha cannot ripple
   * along the centerline the way slope-corrected vertical-distance estimates
   * do at curvature extrema.
   */
  HS_O3_FN float polyline_distance(float t_norm, float polar,
                                   float sin_polar) const {
    const float base =
        target_angle - polar; // knot m sits at v = base + knots[m]

    // Prefilter: when a chunk's arc exceeds the stroke reach, only the pixel's
    // chunk and its neighbours can hold a within-reach curve point; a pixel
    // whose polar offset clears all three knot ranges by more than thickness
    // skips the segment search (most band pixels, in a displaced ring).
    constexpr int CHUNKS = KnotPrefilter::CHUNKS;
    const float chunk_u = (TWO_PI_F / CHUNKS) * sin_polar;
    if (chunk_u >= thickness) {
      const KnotPrefilter &pf = *prefilter;
      int c = static_cast<int>(t_norm * CHUNKS);
      if (c >= CHUNKS)
        c = CHUNKS - 1;
      int cl = c == 0 ? CHUNKS - 1 : c - 1;
      int cr = c == CHUNKS - 1 ? 0 : c + 1;
      float lo = std::min(pf.lo[cl], std::min(pf.lo[c], pf.lo[cr]));
      float hi = std::max(pf.hi[cl], std::max(pf.hi[c], pf.hi[cr]));
      float gap = std::max(base + lo, -(base + hi));
      if (gap > thickness)
        return gap;
    }

    float x = t_norm * lut_n;
    int j = static_cast<int>(x);
    if (j >= lut_n) // t_norm * lut_n can round up to lut_n at the seam
      j = lut_n - 1;
    float f = x - j;
    const float cell_u = (TWO_PI_F / lut_n) * sin_polar;
    const float cell_u2 = cell_u * cell_u;

    // Chart v of knot m relative to the pixel (pixel at the origin).
    auto knot_v = [&](int m) {
      if (m >= lut_n)
        m -= lut_n;
      else if (m < 0)
        m += lut_n;
      return base + knots[m];
    };
    // Perpendicular foot on the segment rising cell_u wide from (u0, v0) to
    // v1; contributes only when the foot lands inside the segment, so the
    // division is paid roughly once per pixel.
    auto interior_d2 = [&](float u0, float v0, float v1, float &best2) {
      float dv = v1 - v0;
      float numer = -(u0 * cell_u + v0 * dv);
      float len2 = cell_u2 + dv * dv;
      if (numer > 0.0f && numer < len2) {
        float cross = u0 * dv - v0 * cell_u;
        best2 = std::min(best2, cross * cross / len2);
      }
    };

    // Straddling cell first, then knot-by-knot outward on each arm; endpoint
    // distances are shared between adjacent segments so each step loads one
    // new knot.
    float ul = -f * cell_u; // arm frontiers: knots j (left), j + 1 (right)
    float ur = (1.0f - f) * cell_u;
    float vl = knot_v(j);
    float vr = knot_v(j + 1);
    float best2 = std::min(ul * ul + vl * vl, ur * ur + vr * vr);
    interior_d2(ul, vl, vr, best2);

    const float th2 = thickness * thickness;
    float bound2 = std::min(best2, th2);
    const int max_o = std::min(MAX_SEARCH_CELLS, lut_n / 2 + 1);
    for (int o = 1; o <= max_o; ++o) {
      // A segment past a frontier knot can't beat that knot's |u|.
      bool right = ur * ur < bound2;
      bool left = ul * ul < bound2;
      if (!right && !left)
        break;
      if (right) {
        float un = ur + cell_u;
        float vn = knot_v(j + o + 1);
        best2 = std::min(best2, un * un + vn * vn);
        interior_d2(ur, vr, vn, best2);
        ur = un;
        vr = vn;
      }
      if (left) {
        float un = ul - cell_u;
        float vn = knot_v(j - o);
        best2 = std::min(best2, un * un + vn * vn);
        interior_d2(un, vn, vl, best2);
        ul = un;
        vl = vn;
      }
      bound2 = std::min(bound2, best2);
    }
    // Past the stroke reach the search leaves best2 an upper bound only, so
    // report the reject band's far sentinel and skip the sqrt. Returning
    // thickness instead would land dist on exactly 0, which a CSG parent reads
    // as on-surface across the whole bounding annulus. Every unsearched point
    // is at least a frontier |u| away, so folding the frontiers in keeps a
    // budget-truncated search (near-pole chart compression) from sentinelling a
    // pixel the unreached cells could still light.
    float lb2 = std::min(best2, std::min(ul * ul, ur * ur));
    return lb2 >= th2 ? FAR_SENTINEL : sqrtf(lb2);
  }
};

/**
 * @brief Undisplaced DistortedRing geometry with exact polar distance.
 */
struct FlatDistortedRing : private DistortedRing {
  using DistortedRing::get_horizontal_intervals;
  using DistortedRing::get_vertical_bounds;
  using DistortedRing::is_solid;
  using DistortedRing::suppress_pole_fill;
  using DistortedRing::thickness;

  /**
   * @brief Builds an undisplaced ring using exact polar centerline distance.
   * @param b Orientation frame (v = ring axis); retained by reference, so it
   *          must outlive the shape.
   * @param r Ring radius as a fraction of the hemisphere.
   * @param th Half-width of the stroke (radians).
   * @param ph Azimuth phase offset (radians).
   */
  FlatDistortedRing(const Basis &b, float r, float th, float ph = 0.0f)
      : DistortedRing(b, r, th, ScalarFn{}, 0.0f, ph) {}

  /**
   * @brief Deleted constructor from a temporary Basis.
   * @details The ring retains its basis by reference, so binding a temporary
   * would leave every later read of basis dangling.
   */
  FlatDistortedRing(const Basis &&, float, float, float = 0.0f) = delete;

  /**
   * @brief Computes signed distance to the undisplaced ring, writing into res.
   * @tparam ComputeUVs When true, also computes the azimuthal t parameter.
   * @param p Point on sphere (normalized).
   * @param res Output result; dist = signed distance, raw_dist = unsigned
   *        centerline distance, t = azimuth in [0,1) when ComputeUVs.
   */
  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    float d = dot(p, normal);
    if (d < cos_min_limit || d > cos_max_limit) {
      res = DistanceResult(FAR_SENTINEL, 0.0f, FAR_SENTINEL, 0.0f, thickness);
      return;
    }

    float polar = fast_acos(hs::clamp(d, -1.0f, 1.0f));
    float t_norm = 0.0f;
    if constexpr (ComputeUVs) {
      t_norm = wrap_t(basis_azimuth(p, u, w, phase) / TWO_PI_F);
    }

    float dist = std::abs(polar - target_angle);
    res = DistanceResult(dist - thickness, t_norm, dist, 0.0f, thickness);
  }
};

} // namespace SDF
