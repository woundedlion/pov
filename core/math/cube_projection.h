/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file cube_projection.h @brief Dominant-axis cube-face coordinates. */

#include "core/math/vector.h"

namespace math {
/** @brief Cube face and unnormalized local coordinates. */
struct CubeFaceProjection {
  int face;    ///< Face index in +X, -X, +Y, -Y, +Z, -Z order.
  float u;     ///< Horizontal numerator; divide by major for [-1, 1].
  float v;     ///< Vertical numerator; divide by major for [-1, 1].
  float major; ///< Dominant component magnitude.
};

/** @brief Projects a direction onto its dominant cube face.
 * @param direction Finite direction; zero returns a zero major magnitude.
 * @return Raw face coordinates, with X then Y winning equal-magnitude ties. */
HS_HOT_INLINE inline CubeFaceProjection
project_cube_face(const Vector &direction) {
  const float ax = fabsf(direction.x);
  const float ay = fabsf(direction.y);
  const float az = fabsf(direction.z);
  if (ax >= ay && ax >= az)
    return {direction.x >= 0.0f ? 0 : 1,
            direction.x >= 0.0f ? direction.z : -direction.z, direction.y, ax};
  if (ay >= az)
    return {direction.y >= 0.0f ? 2 : 3, direction.x,
            direction.y >= 0.0f ? direction.z : -direction.z, ay};
  return {direction.z >= 0.0f ? 4 : 5,
          direction.z >= 0.0f ? direction.x : -direction.x, direction.y, az};
}
} // namespace math
