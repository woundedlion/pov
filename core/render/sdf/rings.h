/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <cmath>
#include "math/geometry.h"
#include "platform/constants.h"
#include "engine/concepts.h"
#include "render/sdf/common.h"

/**
 * @file rings.h
 * @brief The ring leaves: Ring, DistortedRing and FlatDistortedRing.
 */

namespace SDF {

/**
 * @brief Calculates signed distance to a ring.
 * @details Register semantics: the DistanceResult table (stroke row: Ring).
 */
struct Ring {
  float radius;    /**< Ring radius as a fraction of the hemisphere. */
  float thickness; /**< Half-width of the stroke (radians). */
  float phase;     /**< Azimuth phase offset (radians). */

  math::Vector normal, u,
      w;    /**< Ring axis and the two in-plane basis vectors. */
  float ny; /**< y-component of the ring axis. */
  float target_angle,
      center_phi; /**< Centerline polar angle and axis colatitude. */
  float cos_max, cos_min, cos_target,
      inv_sin_target; /**< Precomputed band trig. */

  float r_val;       /**< Horizontal projection length of the axis (for full-row
                         check). */
  float alpha_angle; /**< Azimuth angle of the normal vector in the XZ plane. */
  static constexpr bool is_solid = false; /**< Ring renders as a stroke. */

  Ring() = default;

  /**
   * @brief Builds a ring from its basis, radius, thickness, and phase.
   * @param b Orientation frame (v = ring axis).
   * @param r Ring radius as a fraction of the hemisphere.
   * @param th Half-width of the stroke (radians).
   * @param ph Azimuth phase offset (radians).
   */
  Ring(const math::Basis &b, float r, float th, float ph = 0)
      : radius(r), thickness(th), phase(ph) {
    HS_CHECK(radius >= 0.0f && radius <= 2.0f, "Ring: radius outside [0, 2]");
    // A negative half-width inverts the band, culling every probe.
    HS_CHECK(thickness >= 0.0f, "Ring: negative stroke half-width");
    normal = b.v;
    u = b.u;
    w = b.w;
    AxisProjection ap = project_axis(normal);
    ny = ap.ny;

    target_angle = radius * (math::PI_F / 2.0f);
    center_phi = acosf(std::max(-1.0f, std::min(1.0f, ny)));

    float ang_min = std::max(0.0f, target_angle - thickness);
    float ang_max = std::min(math::PI_F, target_angle + thickness);
    cos_max = cosf(ang_min);
    cos_min = cosf(ang_max);
    if (cos_max == 1.0f)
      cos_max = 2.0f;
    if (cos_min == -1.0f)
      cos_min = -2.0f;
    cos_target = cosf(target_angle);

    bool safe_approx =
        (target_angle > POLE_SAFE_MARGIN &&
         target_angle < math::PI_F - POLE_SAFE_MARGIN &&
         thickness < RING_LINEARIZE_TAN_FRAC * std::abs(tanf(target_angle)));
    inv_sin_target = safe_approx ? (1.0f / sinf(target_angle)) : 0.0f;

    r_val = ap.r_val;
    alpha_angle = ap.alpha_angle;
  }

  /**
   * @brief Maps the ring's latitude band to its inclusive scanline row range.
   * @tparam H Canvas height in rows.
   * @return Row bounds covering the ring plus its AA falloff.
   */
  template <int H> Bounds get_vertical_bounds() const {
    PhiBand band = clamp_phi_band(center_phi, target_angle);

    // Exact distance trims the outer 5% (negligible quintic alpha); the
    // linearized metric retains the full stroke band.
    float eff_th = inv_sin_target != 0.0f ? thickness : 0.95f * thickness;
    float f_phi_min = std::max(0.0f, band.phi_min - eff_th);
    float f_phi_max = std::min(math::PI_F, band.phi_max + eff_th);

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
   * @return True if the row was handled, possibly with no intervals; false
   *         requests a full scan.
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    float cos_phi = math::TrigLUT<W, H>::cos_phi[y];
    float sin_phi = math::TrigLUT<W, H>::sin_phi[y];

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
      float polar = math::fast_acos(hs::clamp(d, -1.0f, 1.0f));
      dist = std::abs(polar - target_angle);
    }
    float sd = dist - thickness;
    if (sd >= 0.0f || thickness <= 0.0f)
      return 0.0f;
    return math::quintic_kernel(-sd / thickness);
  }

  /**
   * @brief Computes signed distance to the ring, writing into res.
   * @tparam ComputeUVs When true, also computes the azimuthal t parameter.
   * @param p Point on sphere (normalized).
   * @param res Output result; dist = signed distance, raw_dist = unsigned
   *        centerline distance, t = azimuth in [0,1) when ComputeUVs.
   */
  template <bool ComputeUVs = true>
  void distance(const math::Vector &p, DistanceResult &res) const {
    float d = math::dot(p, normal);
    if (d < cos_min || d > cos_max) {
      res = DistanceResult(FAR_SENTINEL, 0.0f, FAR_SENTINEL, 0.0f, thickness);
      return;
    }

    float dist = 0;
    if (inv_sin_target != 0) {
      dist = std::abs(d - cos_target) * inv_sin_target;
    } else {
      float polar = math::fast_acos(hs::clamp(d, -1.0f, 1.0f));
      dist = std::abs(polar - target_angle);
    }

    float t = 0.0f;
    if constexpr (ComputeUVs) {
      t = math::wrap_t(basis_azimuth(p, u, w, phase) / math::TWO_PI_F);
    }

    res = DistanceResult(dist - thickness, t, dist, 0.0f, thickness);
  }
};

/**
 * @brief Per-azimuth-chunk knot ranges backing DistortedRing's knot prefilter.
 * @details Caller-owned; one instance per knot-mode ring; must outlive the
 * shape.
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
  const math::Basis &basis; /**< Orientation frame (v = ring axis); retained by
                         reference, so it must outlive the shape. */
  float radius;             /**< Ring radius as a fraction of the hemisphere. */
  float thickness;          /**< Half-width of the stroke (radians). */
  float thickness2;         /**< Squared half-width of the stroke. */
  ScalarFn shift_fn; /**< Per-azimuth centerline shift, t in [0,1) -> radians;
                         empty in knot mode. */
  const float *knots =
      nullptr; /**< Optional lut_n shift knots; entry lut_n is ignored if present.
                    Selects exact polyline distance. */
  int lut_n = 0;                /**< Knot cell count when knots is set. */
  float knot_count = 0.0f;      /**< Knot cell count as a float. */
  float knot_cell_angle = 0.0f; /**< Angular width of one knot cell. */
  const KnotPrefilter *prefilter =
      nullptr;          /**< Caller-owned chunk ranges over knots. */
  float max_distortion; /**< Maximum magnitude of the shift (radians). */
  float phase;          /**< Azimuth phase offset (radians). */

  math::Vector normal, u,
      w;    /**< Ring axis and the two in-plane basis vectors. */
  float ny; /**< y-component of the ring axis. */
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
   *           bands, so an underestimate silently culls genuine arcs. Not
   *           checked here.
   * @param ph Azimuth phase offset (radians).
   * @details distance() measures polar separation at the query azimuth. It is
   *          not the exact distance to the shifted centerline, so steep shifts
   *          can render thinner than the requested stroke; the knot
   *          constructor takes the exact polyline distance instead.
   */
  DistortedRing(const math::Basis &b, float r, float th, ScalarFn sf, float md,
                float ph)
      : DistortedRing(b, r, th, md, ph) {
    HS_CHECK(sf, "DistortedRing: shift_fn must be non-null");
    HS_CHECK(md >= 0.0f, "DistortedRing: negative maximum distortion");
    shift_fn = sf;
  }

  /**
   * @brief Deleted constructor from a temporary Basis.
   * @details The ring retains its basis by reference, so binding a temporary
   * would leave every later read of basis dangling.
   */
  DistortedRing(const math::Basis &&, float, float, ScalarFn, float,
                float) = delete;

protected:
  /**
   * @brief Builds the ring geometry shared by the modes carrying no shift
   *        callback.
   * @param b Orientation frame (v = ring axis); retained by reference, so it
   *          must outlive the shape.
   * @param r Ring radius as a fraction of the hemisphere.
   * @param th Half-width of the stroke (radians).
   * @param md Maximum magnitude of the centerline shift (radians).
   * @param ph Azimuth phase offset (radians).
   */
  DistortedRing(const math::Basis &b, float r, float th, float md, float ph)
      : basis(b), radius(r), thickness(th), thickness2(th * th),
        max_distortion(md), phase(ph) {
    HS_CHECK(radius >= 0.0f && radius <= 2.0f,
             "DistortedRing: radius outside [0, 2]");
    // A negative half-width inverts the band, culling every probe.
    HS_CHECK(thickness >= 0.0f, "DistortedRing: negative stroke half-width");
    normal = basis.v;
    u = basis.u;
    w = basis.w;
    AxisProjection ap = project_axis(normal);
    ny = ap.ny;
    target_angle = radius * (math::PI_F / 2.0f);
    center_phi = acosf(std::max(-1.0f, std::min(1.0f, ny)));
    max_thickness = thickness + max_distortion;

    r_val = ap.r_val;
    alpha_angle = ap.alpha_angle;

    float ang_min = std::max(0.0f, target_angle - max_thickness);
    float ang_max = std::min(math::PI_F, target_angle + max_thickness);
    cos_max_limit = cosf(ang_min);
    cos_min_limit = cosf(ang_max);
    if (cos_max_limit == 1.0f)
      cos_max_limit = 2.0f;
    if (cos_min_limit == -1.0f)
      cos_min_limit = -2.0f;
  }

public:
  /**
   * @brief Builds a distorted ring whose centerline is a shift-knot polyline.
   * @param b Orientation frame (v = ring axis); retained by reference, so it
   *          must outlive the shape.
   * @param r Ring radius as a fraction of the hemisphere.
   * @param th Half-width of the stroke (radians).
   * @param kn n centerline shifts (radians), one per equal azimuth cell;
   *           closure wraps to entry 0; entry n is ignored if present.
   *           Must outlive the shape. distance() returns the exact distance to
   *           this polyline (within the local tangent chart).
   * @param n Number of knot cells; at least 3.
   * @param ph Azimuth phase offset (radians).
   * @param pf Prefilter storage filled here; must outlive the shape. nullptr
   *           skips the per-pixel prefilter, for callers that cull candidate
   *           pixels themselves.
   */
  DistortedRing(const math::Basis &b, float r, float th, const float *kn, int n,
                float ph, KnotPrefilter *pf)
      : DistortedRing(b, r, th, 0.0f, ph) {
    HS_CHECK(kn != nullptr && n >= 3,
             "DistortedRing: knot storage requires at least three knots");
    knots = kn;
    lut_n = n;
    knot_count = static_cast<float>(n);
    knot_cell_angle = math::TWO_PI_F / n;
    prefilter = pf;
    float min_shift = kn[0];
    float max_shift = kn[0];
    if (pf) {
      // Per-chunk knot ranges for the per-pixel prefilter. A segment registers
      // its endpoints in every chunk its azimuth extent touches (a straddling
      // segment spans two), so the polyline inside a chunk never leaves
      // [pf->lo, pf->hi].
      for (int c = 0; c < KnotPrefilter::CHUNKS; ++c) {
        pf->lo[c] = 1e9f;
        pf->hi[c] = -1e9f;
      }
      for (int k = 0; k < n; ++k) {
        const float next = kn[k + 1 == n ? 0 : k + 1];
        float lo = std::min(kn[k], next);
        float hi = std::max(kn[k], next);
        min_shift = std::min(min_shift, lo);
        max_shift = std::max(max_shift, hi);
        int c1 = k * KnotPrefilter::CHUNKS / n;
        int c2 = std::min((k + 1) * KnotPrefilter::CHUNKS / n,
                          KnotPrefilter::CHUNKS - 1);
        for (int c = c1; c <= c2; ++c) {
          pf->lo[c] = std::min(pf->lo[c], lo);
          pf->hi[c] = std::max(pf->hi[c], hi);
        }
      }
    } else {
      for (int k = 1; k < n; ++k) {
        min_shift = std::min(min_shift, kn[k]);
        max_shift = std::max(max_shift, kn[k]);
      }
    }

    max_distortion = std::max(std::abs(min_shift), std::abs(max_shift));
    max_thickness = thickness + max_distortion;
    float ang_min =
        hs::clamp(target_angle + min_shift - thickness, 0.0f, math::PI_F);
    float ang_max =
        hs::clamp(target_angle + max_shift + thickness, 0.0f, math::PI_F);
    cos_max_limit = cosf(ang_min);
    cos_min_limit = cosf(ang_max);
    if (cos_max_limit == 1.0f)
      cos_max_limit = 2.0f;
    if (cos_min_limit == -1.0f)
      cos_min_limit = -2.0f;
  }

  /**
   * @brief Deleted constructor from a temporary Basis.
   * @details The ring retains its basis by reference, so binding a temporary
   * would leave every later read of basis dangling.
   */
  DistortedRing(const math::Basis &&, float, float, const float *, int, float,
                KnotPrefilter *) = delete;

  /** @brief The knot constructor with its prefilter storage attached. */
  DistortedRing(const math::Basis &b, float r, float th, const float *kn, int n,
                float ph, KnotPrefilter &pf)
      : DistortedRing(b, r, th, kn, n, ph, &pf) {}

  /**
   * @brief Deleted constructor from a temporary Basis.
   * @details The ring retains its basis by reference, so binding a temporary
   * would leave every later read of basis dangling.
   */
  DistortedRing(const math::Basis &&, float, float, const float *, int, float,
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
    float f_phi_max = std::min(math::PI_F, band.phi_max + margin);

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
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    float cos_phi = math::TrigLUT<W, H>::cos_phi[y];
    float sin_phi = math::TrigLUT<W, H>::sin_phi[y];

    if (r_val < MIN_HORIZONTAL_PROJ)
      return false;

    float denom = r_val * sin_phi;
    if (std::abs(denom) < INTERVAL_DENOM_EPS)
      // Exact pole row: every column aliases to the one pole point, so the
      // default full-row scan fills the whole row; suppress_pole_fill drops it.
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
  void distance(const math::Vector &p, DistanceResult &res) const {
    float d = math::dot(p, normal);
    // Early reject: outside bounding annulus
    if (d < cos_min_limit || d > cos_max_limit) {
      res = DistanceResult(FAR_SENTINEL, 0.0f, FAR_SENTINEL, 0.0f, thickness);
      return;
    }
    float polar = math::fast_acos(hs::clamp(d, -1.0f, 1.0f));

    float t_norm = math::wrap_t(basis_azimuth(p, u, w, phase) / math::TWO_PI_F);

    float dist;
    if (knots)
      dist = polyline_distance(t_norm, polar,
                               sqrtf(fmaxf(1.0f - d * d, POLE_SIN2_FLOOR)));
    else
      dist = std::abs(polar - (target_angle + shift_fn(t_norm)));

    if constexpr (!ComputeUVs)
      t_norm = 0.0f;

    res = DistanceResult(dist - thickness, t_norm, dist, 0.0f, thickness);
  }

  /**
   * @brief distance<true>() of a knot ring from a precomputed pixel frame.
   * @param d Pixel dot ring axis (= dot(p, normal)).
   * @param polar fast_acos(clamp(d, -1, 1)).
   * @param sin_polar sqrtf(max(1 - d * d, POLE_SIN2_FLOOR)).
   * @param t_norm Pixel azimuth in [0, 1), phase applied.
   * @param res Output result. Wherever the stroke can light (dist < 0) it is
   *        equal to distance<true>() to float rounding; elsewhere dist is only
   *        known to be non-negative.
   * @pre The ring is in knot mode.
   * @details Skips the chunk prefilter.
   */
  __attribute__((always_inline)) void
  distance_from_frame(float d, float polar, float sin_polar, float t_norm,
                      DistanceResult &res) const {
    if (d < cos_min_limit || d > cos_max_limit) {
      res = DistanceResult(FAR_SENTINEL, 0.0f, FAR_SENTINEL, 0.0f, thickness);
      return;
    }
    const float dist = max_distortion > 0.0f
                           ? polyline_window(t_norm, polar, sin_polar)
                           : std::abs(polar - target_angle);
    res = DistanceResult(dist - thickness, t_norm, dist, 0.0f, thickness);
  }

  static constexpr float POLE_SIN2_FLOOR =
      1e-6f; /**< sin^2 floor keeping the chart's azimuth scale positive at the
                poles. */

  static constexpr int MAX_SEARCH_CELLS =
      64; /**< Outward search budget per side; only near-pole chart compression
             approaches it. */

private:
  /**
   * @brief Distance from a pixel to the knot polyline.
   * @param t_norm Pixel azimuth in [0, 1) (phase applied).
   * @param polar Pixel polar angle (radians).
   * @param sin_polar sqrtf(max(1 - d * d, POLE_SIN2_FLOOR)) for the pixel.
   * @return Geodesic distance to the nearest polyline point (radians), exact
   *         within the local tangent chart for distances up to `thickness`, or
   *         a lower bound when the outward search hits its cell budget; past
   *         that reach either a value above `thickness` that may over-estimate
   *         (prefilter, which bounds only the pixel's own chunk and its
   *         neighbours) or FAR_SENTINEL, never the reach itself.
   * @details Works in the chart (azimuth * sin(polar), polar) centered on the
   * pixel: exact point-to-segment distances, searched outward from the
   * pixel's own cell. A segment o cells away is at least (o - 1) * cell_u
   * away, so the search stops once that gap exceeds the best distance found
   * (or the stroke reach, past which alpha is zero regardless).
   */
  HS_O3_FN float polyline_distance(float t_norm, float polar,
                                   float sin_polar) const {
    return polyline_search<true>(t_norm, polar, sin_polar);
  }

  /**
   * @brief polyline_distance()'s body, inlined into each caller.
   * @tparam UsePrefilter Consult the chunk prefilter when one is attached.
   */
  template <bool UsePrefilter>
  __attribute__((always_inline)) float
  polyline_search(float t_norm, float polar, float sin_polar) const {
    const float base =
        target_angle - polar; // knot m sits at v = base + knots[m]

    // Prefilter: when a chunk's arc exceeds the stroke reach, only the pixel's
    // chunk and its neighbours can hold a within-reach curve point; a pixel
    // whose polar offset clears all three knot ranges by more than thickness
    // skips the segment search.
    constexpr int CHUNKS = KnotPrefilter::CHUNKS;
    const float chunk_u = (math::TWO_PI_F / CHUNKS) * sin_polar;
    if (UsePrefilter && prefilter && chunk_u >= thickness) {
      const KnotPrefilter &pf = *prefilter;
      int c = static_cast<int>(t_norm * CHUNKS);
      if (c >= CHUNKS)
        c = CHUNKS - 1;
      int cl = c == 0 ? CHUNKS - 1 : c - 1;
      int cr = c == CHUNKS - 1 ? 0 : c + 1;
      float lo = fminf(pf.lo[cl], fminf(pf.lo[c], pf.lo[cr]));
      float hi = fmaxf(pf.hi[cl], fmaxf(pf.hi[c], pf.hi[cr]));
      float gap = fmaxf(base + lo, -(base + hi));
      if (gap > thickness)
        return gap;
    }

    float x = t_norm * knot_count;
    int j = static_cast<int>(x);
    if (j >= lut_n) // t_norm * lut_n can round up to lut_n at the seam
      j = lut_n - 1;
    float f = x - j;
    const float cell_u = knot_cell_angle * sin_polar;
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
    // v1; contributes only when the foot lands inside the segment.
    auto interior_d2 = [&](float u0, float v0, float v1, float &best2) {
      float dv = v1 - v0;
      float numer = -(u0 * cell_u + v0 * dv);
      float len2 = cell_u2 + dv * dv;
      if (numer > 0.0f && numer < len2) {
        float cross = u0 * dv - v0 * cell_u;
        best2 = fminf(best2, cross * cross / len2);
      }
    };

    // Straddling cell first, then knot-by-knot outward on each arm; endpoint
    // distances are shared between adjacent segments so each step loads one
    // new knot.
    float ul = -f * cell_u; // arm frontiers: knots j (left), j + 1 (right)
    float ur = (1.0f - f) * cell_u;
    float vl = knot_v(j);
    float vr = knot_v(j + 1);
    float best2 = fminf(ul * ul + vl * vl, ur * ur + vr * vr);
    interior_d2(ul, vl, vr, best2);

    const float th2 = thickness2;
    float bound2 = fminf(best2, th2);
    const int sweep_o = lut_n / 2 + 1; // covers every knot on both arms
    const bool budget_capped = MAX_SEARCH_CELLS < sweep_o;
    const int max_o = budget_capped ? MAX_SEARCH_CELLS : sweep_o;
    for (int o = 1; o <= max_o; ++o) {
      // A segment past a frontier knot can't beat that knot's |u|.
      bool right = ur * ur < bound2;
      bool left = ul * ul < bound2;
      if (!right && !left)
        break;
      if (right) {
        float un = ur + cell_u;
        float vn = knot_v(j + o + 1);
        best2 = fminf(best2, un * un + vn * vn);
        interior_d2(ur, vr, vn, best2);
        ur = un;
        vr = vn;
      }
      if (left) {
        float un = ul - cell_u;
        float vn = knot_v(j - o);
        best2 = fminf(best2, un * un + vn * vn);
        interior_d2(un, vn, vl, best2);
        ul = un;
        vl = vn;
      }
      bound2 = fminf(bound2, best2);
    }
    if (best2 < th2)
      return sqrtf(best2);
    // A budget-capped search leaves knots unvisited, each at least a frontier
    // |u| away; report that bound.
    if (budget_capped) {
      float frontier2 = fminf(ul * ul, ur * ur);
      if (frontier2 < th2)
        return sqrtf(frontier2);
    }
    // Past the stroke reach best2 is only an upper bound: report the far
    // sentinel (dist == 0 would read as on-surface to a CSG parent).
    return FAR_SENTINEL;
  }

  /** @brief polyline_search<false>() kept out of the window's inline body. */
  HS_HOT_FLASH_MEMBER float polyline_search_far(float t_norm, float polar,
                                                float sin_polar) const {
    return polyline_search<false>(t_norm, polar, sin_polar);
  }

  /**
   * @brief polyline_search<false>() over a fixed knot window.
   * @details Once a knot cell spans at least half the stroke reach, the
   * outward search never passes the straddling cell's neighbours plus one more
   * cell on a side whose frontier knot sits inside the reach, so that window is
   * evaluated straight through. Narrower cells take the full search. Agrees
   * with polyline_search<false>() to float rounding.
   */
  __attribute__((always_inline)) float
  polyline_window(float t_norm, float polar, float sin_polar) const {
    const float cell_u = knot_cell_angle * sin_polar;
    if (2.0f * cell_u < thickness || lut_n < 8)
      return polyline_search_far(t_norm, polar, sin_polar);
    const float base = target_angle - polar;
    const float x = t_norm * knot_count;
    int j = static_cast<int>(x);
    if (j >= lut_n)
      j = lut_n - 1;
    const float f = x - j;
    const float cell_u2 = cell_u * cell_u;
    const float th2 = thickness2;

    float k[6];
    if (j >= 2 && j + 3 < lut_n) {
      for (int i = 0; i < 6; ++i)
        k[i] = knots[j - 2 + i];
    } else {
      for (int i = 0; i < 6; ++i) {
        int m = j - 2 + i;
        m = m < 0 ? m + lut_n : m >= lut_n ? m - lut_n : m;
        k[i] = knots[m];
      }
    }

    // Best squared distance so far, as num / den.
    float num, den = 1.0f;
    auto interior = [&](float u0, float v0, float v1) {
      const float dv = v1 - v0;
      const float numer = -(u0 * cell_u + v0 * dv);
      const float len2 = cell_u2 + dv * dv;
      const float cross = u0 * dv - v0 * cell_u;
      const float c2 = cross * cross;
      const bool better = numer > 0.0f && numer < len2 && c2 * den < num * len2;
      num = better ? c2 : num;
      den = better ? len2 : den;
    };

    const float ul = -f * cell_u;
    const float ur = (1.0f - f) * cell_u;
    const float vl = base + k[2];
    const float vr = base + k[3];
    const float ur1 = ur + cell_u;
    const float vr1 = base + k[4];
    const float ul1 = ul - cell_u;
    const float vl1 = base + k[1];
    num = fminf(fminf(ul * ul + vl * vl, ur * ur + vr * vr),
                fminf(ur1 * ur1 + vr1 * vr1, ul1 * ul1 + vl1 * vl1));
    interior(ul, vl, vr);
    interior(ur, vr, vr1);
    interior(ul1, vl1, vl);
    // One more cell on a side whose frontier still sits inside the bound.
    if (ur1 * ur1 < th2 && ur1 * ur1 * den < num) {
      const float un = ur1 + cell_u;
      const float vn = base + k[5];
      const float e = un * un + vn * vn;
      if (e * den < num) {
        num = e;
        den = 1.0f;
      }
      interior(ur1, vr1, vn);
    }
    if (ul1 * ul1 < th2 && ul1 * ul1 * den < num) {
      const float un = ul1 - cell_u;
      const float vn = base + k[0];
      const float e = un * un + vn * vn;
      if (e * den < num) {
        num = e;
        den = 1.0f;
      }
      interior(un, vn, vl1);
    }
    if (num < th2 * den)
      return sqrtf(num / den);
    return FAR_SENTINEL;
  }
};

/**
 * @brief Undisplaced DistortedRing geometry with exact polar distance.
 */
struct FlatDistortedRing : private DistortedRing {
  using DistortedRing::get_horizontal_intervals;
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
  FlatDistortedRing(const math::Basis &b, float r, float th, float ph = 0.0f)
      : DistortedRing(b, r, th, 0.0f, ph) {}

  /**
   * @brief Deleted constructor from a temporary Basis.
   * @details The ring retains its basis by reference, so binding a temporary
   * would leave every later read of basis dangling.
   */
  FlatDistortedRing(const math::Basis &&, float, float, float = 0.0f) = delete;

  /**
   * @brief Maps the undisplaced ring's latitude band to its row range.
   * @tparam H Canvas height in rows.
   * @return Inclusive row bounds covering the stroke.
   * @details distance() is the exact polar offset, so the band widened by
   * `thickness` is tight.
   */
  template <int H> Bounds get_vertical_bounds() const {
    PhiBand band = clamp_phi_band(center_phi, target_angle);
    return phi_bounds_to_rows<H>(
        std::max(0.0f, band.phi_min - thickness),
        std::min(math::PI_F, band.phi_max + thickness));
  }

  /**
   * @brief Computes signed distance to the undisplaced ring, writing into res.
   * @tparam ComputeUVs When true, also computes the azimuthal t parameter.
   * @param p Point on sphere (normalized).
   * @param res Output result; dist = signed distance, raw_dist = unsigned
   *        centerline distance, t = azimuth in [0,1) when ComputeUVs.
   */
  template <bool ComputeUVs = true>
  void distance(const math::Vector &p, DistanceResult &res) const {
    float d = math::dot(p, normal);
    if (d < cos_min_limit || d > cos_max_limit) {
      res = DistanceResult(FAR_SENTINEL, 0.0f, FAR_SENTINEL, 0.0f, thickness);
      return;
    }

    float polar = math::fast_acos(hs::clamp(d, -1.0f, 1.0f));
    float t_norm = 0.0f;
    if constexpr (ComputeUVs) {
      t_norm = math::wrap_t(basis_azimuth(p, u, w, phase) / math::TWO_PI_F);
    }

    float dist = std::abs(polar - target_angle);
    res = DistanceResult(dist - thickness, t_norm, dist, 0.0f, thickness);
  }
};

// Leaf roster for the CSG composition contract.
static_assert(SDFShape<Ring>);
static_assert(SDFShape<DistortedRing>);
static_assert(SDFShape<FlatDistortedRing>);

} // namespace SDF
