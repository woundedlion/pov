/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file swept_bound.h
 * @brief Bounds on everything a moving primitive sweeps over a frame, from
 * its vertices' frame-end captures and their paths' deviation bound, for
 * culling.
 */

#include "core/math/3dmath.h"
#include "core/platform/attributes.h"
#include <cmath>

namespace Spacetime {

/**
 * @brief Padded height extremes and pairwise clearance of a frame's swept
 *        primitive's captures.
 * @details Every vertex between its captures, and every great-circle arc
 * between such vertices, lies within the path deviation bound of the
 * normalized convex hull of the captures. A convex combination of unit vectors
 * whose pairwise dot products are at least `d > 0` has length at least
 * `sqrt(d)`, which bounds how far normalizing lifts its height past the
 * captures'. With the pad, `d` is the least pairwise dot less twice the pad,
 * which is the clearance.
 */
struct SweptHull {
  float top;       /**< Greatest capture height plus the pad. */
  float bottom;    /**< Least capture height less the pad. */
  float clearance; /**< Least pairwise capture dot less twice the pad. */

  /** @brief Whether the padded hull clears the origin, so its normalized
   *         heights are bounded.
   * @return True when the clearance is positive. */
  HS_HOT_INLINE bool bounded() const { return clearance > 0.0f; }

  /** @brief Height interval of the directions the swept primitive can reach.
   * @param low Receives the least reachable height.
   * @param high Receives the greatest reachable height.
   * @details Unbounded hulls give [-1, 1]. */
  HS_HOT_INLINE void heights(float &low, float &high) const {
    if (!bounded()) {
      low = -1.0f;
      high = 1.0f;
      return;
    }
    const float bulge = 1.0f / sqrtf(clearance);
    high = fminf(top > 0.0f ? top * bulge : top, 1.0f);
    low = fmaxf(bottom < 0.0f ? bottom * bulge : bottom, -1.0f);
  }
};

/**
 * @brief Measures the swept hull of a frame's primitive.
 * @param captures Unit captures of the primitive's vertices at the frame's
 *        start and end.
 * @param count Number of captures, at least 1.
 * @param pad Bound on any vertex path's distance from its chord between its
 *        captures (Spacetime::Worldline::deviation()).
 * @return The padded extremes and clearance.
 */
HS_HOT_INLINE inline SweptHull swept_hull(const math::Vector *captures,
                                          int count, float pad) {
  float top = captures[0].y, bottom = captures[0].y, least = 1.0f;
  for (int i = 0; i < count; ++i) {
    top = fmaxf(top, captures[i].y);
    bottom = fminf(bottom, captures[i].y);
    for (int j = i + 1; j < count; ++j)
      least = fminf(least, math::dot(captures[i], captures[j]));
  }
  return {top + pad, bottom - pad, least - 2.0f * pad};
}

/**
 * @brief Height interval of the directions a frame's swept primitive can
 *        reach.
 * @param points Unit captures of the primitive's vertices at the frame's
 *        start and end.
 * @param count Number of captures, at least 1.
 * @param pad As swept_hull().
 * @param low Receives the least reachable height.
 * @param high Receives the greatest reachable height.
 * @return False when the captures spread too far apart to bound, leaving
 *         [-1, 1].
 */
HS_HOT_INLINE inline bool swept_heights(const math::Vector *points, int count,
                                        float pad, float &low, float &high) {
  const SweptHull hull = swept_hull(points, count, pad);
  hull.heights(low, high);
  return hull.bounded();
}

/**
 * @brief Whether a frame's swept primitive stays on the side a plane's normal
 *        points to, for a plane through the centre.
 * @param points Unit captures of the primitive's vertices at the frame's
 *        start and end.
 * @param count Number of captures.
 * @param pad As swept_hull().
 * @param normal Plane normal.
 * @return True when every capture lies more than @p pad from the plane on the
 *         side @p normal points to (@p normal need not be unit), so every swept
 *         direction lies on that side.
 * @details The half-space is convex and holds the origin on its boundary, so
 * it holds the captures' normalized convex hull with them.
 */
HS_HOT_INLINE inline bool swept_beyond(const math::Vector *points, int count,
                                       float pad, const math::Vector &normal) {
  const float margin = pad * sqrtf(math::dot(normal, normal));
  for (int i = 0; i < count; ++i)
    if (!(math::dot(points[i], normal) > margin))
      return false;
  return true;
}

} // namespace Spacetime
