/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <concepts>
#include <type_traits>
#include "math/3dmath.h"
#include "render/shading.h"
#include "engine/platform.h"

/**
 * @file sdf_volume.h
 * @brief 3D volumetric signed-distance shapes and the domain warps composed
 * with them, marched by Scan::Volume.
 */

namespace SDF {

// ============================================================================
// 3D Volumetric SDF Shapes (for Scan::Volume raymarching)
// ============================================================================
// These shapes are reached only from Scan::Volume::draw's march loop (Raymarch
// is their sole instantiation site), so they share that region's -O3 options —
// an -Os distance() here would be a wall inside the promoted loop.

HS_O3_BEGIN
/**
 * @brief 3D Torus signed distance field.
 *
 * The torus ring lies in the XZ plane with symmetry axis along Y.
 * Major radius R = ring centerline distance from origin.
 * Minor radius r = tube cross-section radius.
 *
 * Unlike the 2D spherical SDF shapes above (which use DistanceResult),
 * 3D volumetric shapes return plain float distances and operate in
 * Cartesian ray-space.
 */
struct Torus {
  float R; /**< Major radius (ring centerline). */
  float r; /**< Minor radius (tube cross-section). */

  /**
   * @brief Signed distance from a point to the torus surface.
   * @param p Query point in Cartesian ray-space.
   * @return Signed distance (negative inside, positive outside).
   */
  float distance(const Vector &p) const {
    float q = sqrtf(p.x * p.x + p.z * p.z) - R;
    return sqrtf(q * q + p.y * p.y) - r;
  }

  /**
   * @brief Outward normal direction at a point, from a precomputed XZ radius.
   * @param p Query point in Cartesian ray-space.
   * @param inv_xz_len 1 / sqrtf(p.x² + p.z²) for `p`, or 0 on the Y axis.
   * @return The outward normal direction, NOT unit length; the caller
   *         normalizes, or folds `p` through a linear map that it then
   *         normalizes.
   */
  Vector normal_raw(const Vector &p, float inv_xz_len) const {
    float scale = R * inv_xz_len;
    return p - Vector(p.x * scale, 0.0f, p.z * scale);
  }

  /**
   * @brief Surface normal at a point near the torus surface.
   * @param p Query point in Cartesian ray-space.
   * @return Unit outward normal.
   */
  Vector normal(const Vector &p) const {
    float xz_len = sqrtf(p.x * p.x + p.z * p.z);
    float inv_xz_len = (xz_len > TOLERANCE) ? 1.0f / xz_len : 0.0f;
    return normal_raw(p, inv_xz_len).normalized();
  }

  /**
   * @brief Populates a Fragment's registers for shading.
   * @param p Query point in Cartesian ray-space.
   * @param frag Output fragment; v0 = ring angle (0-1, for palette lookup),
   *        v1/v2/v3 = surface normal (x, y, z).
   * @note Volumetric register convention (README "Volumetric Path"), distinct
   *        from Scan's v2 stroke-coverage and mesh face-index conventions.
   */
  void populate(const Vector &p, Fragment &frag) const {
    Vector n = normal(p);
    frag.v0 = (fast_atan2(p.z, p.x) + PI_F) / TWO_PI_F;
    frag.v1 = n.x;
    frag.v2 = n.y;
    frag.v3 = n.z;
  }
};
HS_O3_END

/**
 * @brief Domain warp functions for composing with WarpedSDF.
 */
namespace Warp {

HS_O3_BEGIN
/**
 * @brief Oscillates the Y coordinate sinusoidally around the azimuthal
 * angle θ = atan2(z, x).
 *
 * Produces twisted/undulating geometry when composed with a torus.
 * Provides an analytical Lipschitz bound for safe sphere tracing.
 */
struct Twist {
  int twist;           /**< Number of oscillations around the ring. */
  float amplitude;     /**< Vertical displacement magnitude. */
  float R;             /**< Major radius (needed for the Lipschitz bound). */
  float twist_amp;     /**< twist * amplitude, the warp's angular gradient. */
  float twist_amp_abs; /**< |twist * amplitude|, the Lipschitz numerator. */
  float two_over_r;    /**< 2/R, the Lipschitz reciprocal clamp. */

  /**
   * @brief Constructs a twist warp around a torus of major radius R.
   * @param oscillations Number of oscillations around the ring; must be >= 0.
   *        The harmonic recurrence counts up from 1, so a negative count yields
   *        the first harmonic while the Lipschitz bound describes the requested
   *        one.
   * @param displacement Vertical displacement magnitude; must be >= 0. A
   *        negative one makes bounding_inflation() tighten the bound rather
   *        than relax it, so a sphere trace steps through the surface, and the
   *        amplitude < TOLERANCE tests take a no-warp branch that apply()
   *        contradicts.
   * @param major_radius Major radius; must be > 0. The Lipschitz bound scales
   *        by 2/R, so R == 0 yields a non-finite bound on the XZ axis. Guarded
   *        at the cold construction site, not per-call.
   */
  Twist(int oscillations, float displacement, float major_radius)
      : twist(oscillations), amplitude(displacement), R(major_radius),
        twist_amp(static_cast<float>(oscillations) * displacement),
        twist_amp_abs(fabsf(twist_amp)), two_over_r(2.0f / major_radius) {
    HS_CHECK(R > 0.0f);
    HS_CHECK(twist >= 0);
    HS_CHECK(amplitude >= 0.0f);
  }

  /** @brief Precomputed context: s = sqrtf(x² + z²), shared across
   * apply/lipschitz. */
  using Ctx = float;

  /**
   * @brief Precomputes the shared per-point context s = sqrtf(x² + z²).
   * @param p Query point.
   * @return The radial distance s in the XZ plane.
   */
  Ctx make_ctx(const Vector &p) const { return sqrtf(p.x * p.x + p.z * p.z); }

  /**
   * @brief Warps the domain by displacing Y by amplitude * sin(twist * θ).
   * @param p Query point.
   * @param s Shared context from make_ctx(p).
   * @return The warped point.
   */
  Vector apply(const Vector &p, Ctx s) const {
    return Vector(p.x, p.y - amplitude * sin_ntheta(p, s), p.z);
  }

  /**
   * @brief sin(n*theta) at theta = atan2(z, x), without evaluating either.
   * @param p Query point.
   * @param s Precomputed context (radial distance in the XZ plane).
   * @return sin(twist * theta).
   * @details Chebyshev recurrence sin((k+1)t) = 2cos(t)sin(kt) - sin((k-1)t),
   * seeded from (cos t, sin t) = (x/s, z/s), so the cost is one reciprocal
   * rather than an atan2 and a sine. Exact to float rounding, where
   * fast_atan2/fast_sinf each carry approximation error.
   */
  float sin_ntheta(const Vector &p, Ctx s) const {
    return sin_ntheta_inv(p, s).sin_n;
  }

  /** @brief sin(twist*theta) with the Lipschitz argument that seeded it. */
  struct SinInv {
    float sin_n;         /**< sin(twist * theta). */
    float lipschitz_arg; /**< 1/s; 2/R where the recurrence degenerates. Only
                            valid as a lipschitz()/lipschitz_inv() argument —
                            correct_normal_inv() marks the axis with 0 instead
                            and would read 2/R as a finite radius. */
  };

  /**
   * @brief sin(n*theta) and the Lipschitz argument from one recurrence seed.
   * @param p Query point.
   * @param s Precomputed context (radial distance in the XZ plane).
   * @return sin(twist * theta) and 1/s, the latter feeding lipschitz().
   * @details On the degenerate axis lipschitz_arg carries 2/R, the value the
   * Lipschitz clamp would select there anyway, so no caller needs a second
   * branch.
   */
  SinInv sin_ntheta_inv(const Vector &p, Ctx s) const {
    if (twist == 0 || s < TOLERANCE)
      return {0.0f, two_over_r};
    const float inv_s = 1.0f / s;
    const float two_cos = 2.0f * p.x * inv_s;
    float prev = 0.0f, cur = p.z * inv_s;
    for (int k = 1; k < twist; ++k) {
      const float next = two_cos * cur - prev;
      prev = cur;
      cur = next;
    }
    return {cur, inv_s};
  }

  /**
   * @brief cos(n*theta) at theta = atan2(z, x), by the same recurrence.
   * @param p Query point.
   * @param s Precomputed context (radial distance in the XZ plane).
   * @return cos(twist * theta).
   */
  float cos_ntheta(const Vector &p, Ctx s) const {
    if (twist == 0 || s < TOLERANCE)
      return 1.0f;
    const float inv_s = 1.0f / s;
    const float two_cos = 2.0f * p.x * inv_s;
    float prev = 1.0f, cur = p.x * inv_s;
    for (int k = 1; k < twist; ++k) {
      const float next = two_cos * cur - prev;
      prev = cur;
      cur = next;
    }
    return cur;
  }

  /** @brief Both harmonics of twist*theta from one recurrence pass. */
  struct SinCos {
    float sin_n; /**< sin(twist * theta). */
    float cos_n; /**< cos(twist * theta). */
  };

  /**
   * @brief sin(n*theta) and cos(n*theta) together, sharing one recurrence
   *        setup and loop.
   * @param p Query point.
   * @param s Precomputed context (radial distance in the XZ plane).
   * @return Both harmonics, matching sin_ntheta/cos_ntheta bit for bit.
   * @details For the hit path, which needs both: the march path needs only the
   * sine and calls sin_ntheta so it does not pay for the cosine sequence.
   */
  SinCos sincos_ntheta(const Vector &p, Ctx s) const {
    return sincos_ntheta_inv(p, (s > TOLERANCE) ? 1.0f / s : 0.0f);
  }

  /**
   * @brief sin(n*theta) and cos(n*theta) from an already-computed 1/s.
   * @param p Query point.
   * @param inv_s Reciprocal of the XZ radius; 0 marks the degenerate axis.
   * @return Both harmonics, matching sincos_ntheta bit for bit.
   */
  SinCos sincos_ntheta_inv(const Vector &p, float inv_s) const {
    if (twist == 0 || inv_s == 0.0f)
      return {0.0f, 1.0f};
    const float two_cos = 2.0f * p.x * inv_s;
    float sin_prev = 0.0f, sin_cur = p.z * inv_s;
    float cos_prev = 1.0f, cos_cur = p.x * inv_s;
    for (int k = 1; k < twist; ++k) {
      const float sin_next = two_cos * sin_cur - sin_prev;
      sin_prev = sin_cur;
      sin_cur = sin_next;
      const float cos_next = two_cos * cos_cur - cos_prev;
      cos_prev = cos_cur;
      cos_cur = cos_next;
    }
    return {sin_cur, cos_cur};
  }

  /**
   * @brief Analytical Lipschitz constant of the warp at a point.
   * @param s Precomputed context (radial distance in the XZ plane).
   * @return Operator norm of the warp Jacobian (>= 1), with the radial factor
   * clamped at 1/max(s, R/2).
   * @note Precondition: the composed torus keeps r <= R/2, so its surface never
   * reaches s < R/2, where the clamp under-reports the true norm and a sphere
   * trace would step through.
   */
  float lipschitz(const Vector & /*p*/, Ctx s) const {
    return lipschitz(1.0f / std::max(s, R * 0.5f));
  }

  /**
   * @brief Analytical Lipschitz constant from an already-computed 1/s.
   * @param inv_s Reciprocal of the XZ radius, from sin_ntheta_inv().
   * @return Operator norm of the warp Jacobian (>= 1), with the radial factor
   * clamped as described below.
   * @details The warp Jacobian is the shear I - e_y·gᵀ with e_y ⊥ g and
   * |g| = γ; its operator norm (largest singular value) is γ/2 + √(1 + γ²/4).
   * γ uses |twist·amplitude| so the bound stays conservative regardless of
   * sign, and min(1/s, 2/R) is the clamp 1/max(s, R/2).
   */
  float lipschitz(float inv_s) const {
    if (twist == 0)
      return 1.0f;
    const float gamma = twist_amp_abs * std::min(inv_s, two_over_r);
    return 0.5f * gamma + sqrtf(1.0f + 0.25f * gamma * gamma);
  }

  /**
   * @brief Reciprocal of the Lipschitz constant, from an already-computed 1/s.
   * @param inv_s Reciprocal of the XZ radius, from sin_ntheta_inv().
   * @return 1 / lipschitz(inv_s), in (0, 1].
   * @details (√(1+γ²/4) - γ/2)(√(1+γ²/4) + γ/2) = 1 exactly, so the reciprocal
   * needs no divide. The result never exceeds 1, so callers scale by it
   * unconditionally.
   */
  float lipschitz_inv(float inv_s) const {
    if (twist == 0)
      return 1.0f;
    const float gamma = twist_amp_abs * std::min(inv_s, two_over_r);
    return sqrtf(1.0f + 0.25f * gamma * gamma) - 0.5f * gamma;
  }

  /**
   * @brief Maximum possible inflation of the bounding volume.
   * @return The displacement amplitude (radians of XYZ space).
   */
  float bounding_inflation() const { return amplitude; }

  /**
   * @brief Analytical normal correction via the chain rule (once per hit).
   * @param p Query point.
   * @param base_n Unwarped surface normal.
   * @param s Precomputed context (radial distance in the XZ plane).
   * @return The corrected unit normal accounting for the warp.
   */
  Vector correct_normal(const Vector &p, const Vector &base_n, Ctx s) const {
    if (twist == 0 || amplitude < TOLERANCE)
      return base_n;
    return correct_normal(p, base_n, s, cos_ntheta(p, s));
  }

  /**
   * @brief Normal correction from an already-computed cos(twist*theta).
   * @param p Query point.
   * @param base_n Unwarped surface normal; need not be unit length.
   * @param s Precomputed context (radial distance in the XZ plane).
   * @param cos_n cos(twist * theta) at `p`.
   * @return The corrected unit normal.
   * @details The correction is a linear map of base_n, so scaling base_n scales
   * the result and the final normalize cancels it — an unnormalized base normal
   * gives the identical unit result and saves a normalize.
   */
  Vector correct_normal(const Vector &p, const Vector &base_n, Ctx s,
                        float cos_n) const {
    return correct_normal_inv(p, base_n, (s > TOLERANCE) ? 1.0f / s : 0.0f,
                              cos_n);
  }

  /**
   * @brief Normal correction from an already-computed 1/s and cos(twist*theta).
   * @param p Query point.
   * @param base_n Unwarped surface normal; need not be unit length.
   * @param inv_s Reciprocal of the XZ radius; 0 marks the degenerate axis.
   * @param cos_n cos(twist * theta) at `p`.
   * @return The corrected unit normal.
   */
  Vector correct_normal_inv(const Vector &p, const Vector &base_n, float inv_s,
                            float cos_n) const {
    float dh_dtheta = twist_amp * cos_n;
    float inv_s2 = inv_s * inv_s;

    float dh_dx = dh_dtheta * (-p.z) * inv_s2;
    float dh_dz = dh_dtheta * p.x * inv_s2;

    return Vector(base_n.x - base_n.y * dh_dx, base_n.y,
                  base_n.z - base_n.y * dh_dz)
        .normalized();
  }
};
HS_O3_END

} // namespace Warp

/**
 * @brief The domain-warp interface WarpedVolume's generic path calls.
 * @tparam T Candidate warp type.
 * @details correct_normal() is optional and detected separately; the generic
 * path returns the base normal unchanged when it is absent.
 */
template <typename T>
concept VolumeWarp =
    requires(const T &w, const Vector &p, typename T::Ctx ctx) {
      { w.make_ctx(p) } -> std::same_as<typename T::Ctx>;
      { w.apply(p, ctx) } -> std::same_as<Vector>;
      { w.lipschitz(p, ctx) } -> std::same_as<float>;
      { w.bounding_inflation() } -> std::same_as<float>;
    };

/**
 * @brief The volume-SDF interface WarpedVolume's generic path calls.
 * @tparam T Candidate base shape type.
 */
template <typename T>
concept VolumeShape = requires(const T &s, const Vector &p, Fragment &f) {
  { s.distance(p) } -> std::same_as<float>;
  { s.normal(p) } -> std::same_as<Vector>;
  s.populate(p, f);
};

/**
 * @brief Composable domain-warped volume SDF.
 *
 * Applies a Warp to the input domain of any 3D volume SDF, with:
 * - Bounding fast-path (skips warp trig when far from surface)
 * - Analytical Lipschitz correction for safe sphere tracing
 * - Normal correction via the warp's chain rule
 *
 * The base and warp satisfy VolumeShape and VolumeWarp. The TORUS_TWIST
 * specialization additionally calls the reciprocal-sharing members of the
 * concrete Torus/Twist pair it is keyed on, which neither concept describes.
 */
HS_O3_BEGIN
template <typename SDF, typename Warp> struct WarpedVolume {
  static_assert(VolumeShape<SDF>);
  static_assert(VolumeWarp<Warp>);

  SDF base;  /**< The underlying volume SDF being warped. */
  Warp warp; /**< The domain warp applied before the base SDF. */

  /**
   * @brief Smallest distance the caller needs an accurate value for.
   * @details The cheap bound is returned only strictly above this, so a coarse
   * value can never reach a caller's hit test or antialiasing band. Zero means
   * unstated and selects the warp's maximum displacement, which is safe for any
   * caller but forfeits most of the fast path.
   */
  float precision = 0.0f;

  /** True when the base/warp pair admits the tight per-axis bound below. */
  static constexpr bool TORUS_TWIST = std::is_same_v<SDF, ::SDF::Torus> &&
                                      std::is_same_v<Warp, ::SDF::Warp::Twist>;

  /**
   * @brief Cheap lower bound on the warped distance.
   * @param p Query point in Cartesian ray-space.
   * @return A lower bound on the distance to the warped surface.
   */
  float bounding_distance(const Vector &p) const {
    if constexpr (TORUS_TWIST) {
      // Twist moves only y, by at most bounding_inflation(), so the warped
      // surface lies inside the torus swept +-A along y; this is that solid's
      // exact distance, hence a lower bound. Relaxing the y term alone keeps it
      // far tighter than subtracting A from the whole distance.
      const float q = sqrtf(p.x * p.x + p.z * p.z) - base.R;
      const float dy = std::max(fabsf(p.y) - warp.bounding_inflation(), 0.0f);
      return sqrtf(q * q + dy * dy) - base.r;
    } else {
      return base.distance(p) - warp.bounding_inflation();
    }
  }

  /**
   * @brief Raw warped SDF distance with no Lipschitz correction.
   * @param p Query point in Cartesian ray-space.
   * @return The base SDF evaluated at the warped point (use for surface
   * projection).
   */
  float raw_distance(const Vector &p) const {
    auto ctx = warp.make_ctx(p);
    return base.distance(warp.apply(p, ctx));
  }

  /**
   * @brief March-safe distance with bounding fast-path and Lipschitz
   * correction.
   * @param p Query point in Cartesian ray-space.
   * @return A sphere-tracing-safe (under-estimated) distance to the surface.
   */
  float distance(const Vector &p) const {
    const float gate =
        (precision > 0.0f) ? precision : warp.bounding_inflation();
    if constexpr (TORUS_TWIST) {
      // gate + base.r > 0, so bounding_distance(p) > gate is exactly
      // qq > (gate + base.r)²; the sqrt is then only the fast path's result.
      const float s = sqrtf(p.x * p.x + p.z * p.z);
      const float q = s - base.R;
      const float dy = std::max(fabsf(p.y) - warp.bounding_inflation(), 0.0f);
      const float qq = q * q + dy * dy;
      const float t = gate + base.r;
      if (qq > t * t)
        return sqrtf(qq) - base.r;

      const auto h = warp.sin_ntheta_inv(p, s);
      const Vector warped(p.x, p.y - warp.amplitude * h.sin_n, p.z);
      float d = base.distance(warped);
      if (d > 0.0f)
        d *= warp.lipschitz_inv(h.lipschitz_arg);
      return d;
    } else {
      const float bd = base.distance(p) - warp.bounding_inflation();
      if (bd > gate)
        return bd;

      auto ctx = warp.make_ctx(p);
      float d = base.distance(warp.apply(p, ctx));
      if (d > 0.0f) {
        float lip = warp.lipschitz(p, ctx);
        if (lip > 1.0f)
          d /= lip;
      }
      return d;
    }
  }

  /**
   * @brief Surface normal with the warp's chain-rule correction.
   * @param p Query point in Cartesian ray-space.
   * @return The corrected unit surface normal.
   */
  Vector normal(const Vector &p) const {
    auto ctx = warp.make_ctx(p);
    if constexpr (TORUS_TWIST) {
      // Twist displaces only y, so the warped point keeps p's XZ radius `ctx`
      // and the torus normal needs no second sqrt; one recurrence yields both
      // the sine the warp needs and the cosine the correction needs; and the
      // correction normalizes, so the base normal can stay unnormalized.
      // twist == 0 needs no special case: it gives cos_n = 1, hence a zero
      // gradient and an unchanged normal.
      const float inv_s = (ctx > TOLERANCE) ? 1.0f / ctx : 0.0f;
      if (warp.amplitude < TOLERANCE)
        return base.normal_raw(p, inv_s).normalized();
      auto h = warp.sincos_ntheta_inv(p, inv_s);
      Vector warped(p.x, p.y - warp.amplitude * h.sin_n, p.z);
      return warp.correct_normal_inv(p, base.normal_raw(warped, inv_s), inv_s,
                                     h.cos_n);
    } else {
      Vector base_n = base.normal(warp.apply(p, ctx));
      if constexpr (requires { warp.correct_normal(p, base_n, ctx); }) {
        return warp.correct_normal(p, base_n, ctx);
      }
      return base_n;
    }
  }

  /**
   * @brief Populates a Fragment's registers for shading.
   * @param p Query point in Cartesian ray-space.
   * @param frag Output fragment; v0 = ring angle (0-1), v1/v2/v3 = surface
   * normal.
   * @note Volumetric register convention (README "Volumetric Path"), distinct
   *        from Scan's v2 stroke-coverage and mesh face-index conventions.
   */
  void populate(const Vector &p, Fragment &frag) const {
    Vector n = normal(p);
    frag.v0 = (fast_atan2(p.z, p.x) + PI_F) / TWO_PI_F;
    frag.v1 = n.x;
    frag.v2 = n.y;
    frag.v3 = n.z;
  }
};
HS_O3_END

} // namespace SDF
