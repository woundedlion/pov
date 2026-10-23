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
 * @brief Height interval of the directions a frame's swept primitive can
 *        reach.
 * @param points Unit captures of the primitive's vertices at the frame's
 *        start and end.
 * @param count Number of captures, at least 1.
 * @param pad Bound on any vertex path's distance from its chord between its
 *        captures (Spacetime::Worldline::deviation()).
 * @param low Receives the least reachable height.
 * @param high Receives the greatest reachable height.
 * @return False when the captures spread too far apart to bound, leaving
 *         [-1, 1].
 * @details Every vertex between its captures, and every great-circle arc
 * between such vertices, lies within @p pad of the normalized convex hull of
 * the captures. A convex combination of unit vectors whose pairwise dot
 * products are at least `d > 0` has length at least `sqrt(d)`, which bounds
 * how far normalizing lifts its height past the captures'.
 */
HS_HOT_INLINE inline bool swept_heights(const math::Vector *points, int count,
                                        float pad, float &low, float &high) {
  float top = points[0].y, bottom = points[0].y, least = 1.0f;
  for (int i = 0; i < count; ++i) {
    top = fmaxf(top, points[i].y);
    bottom = fminf(bottom, points[i].y);
    for (int j = i + 1; j < count; ++j)
      least = fminf(least, math::dot(points[i], points[j]));
  }
  top += pad;
  bottom -= pad;
  least -= 2.0f * pad;
  if (!(least > 0.0f)) {
    low = -1.0f;
    high = 1.0f;
    return false;
  }
  const float bulge = 1.0f / sqrtf(least);
  high = fminf(top > 0.0f ? top * bulge : top, 1.0f);
  low = fmaxf(bottom < 0.0f ? bottom * bulge : bottom, -1.0f);
  return true;
}

/**
 * @brief Whether a frame's swept primitive stays on the far side of a plane
 *        through the centre.
 * @param points Unit captures of the primitive's vertices at the frame's
 *        start and end.
 * @param count Number of captures.
 * @param pad As swept_heights().
 * @param normal Plane normal.
 * @return True when every capture lies beyond `pad / |normal|` on the side
 *         @p normal points to, so every swept direction lies on that side.
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
