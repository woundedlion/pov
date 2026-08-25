/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file 4dmath.h
 * @brief Four-dimensional linear algebra: the Vec4 point, the Mat4 transform,
 *        and the plane rotation that composes a 4D orientation.
 */

#include "math/3dmath.h"

/** @brief Component count of the four-dimensional primitives below. */
inline constexpr int VEC4_DIMENSIONS = 4;

/**
 * @brief A point or direction in four dimensions.
 */
struct Vec4 {
  float v[VEC4_DIMENSIONS];

  constexpr float &operator[](int index) { return v[index]; }
  constexpr float operator[](int index) const { return v[index]; }
};

/**
 * @brief A 4x4 matrix in row-major order, applied to Vec4 on the left.
 */
struct Mat4 {
  float m[VEC4_DIMENSIONS][VEC4_DIMENSIONS]{};

  /** @brief The identity transform. */
  static constexpr Mat4 identity() {
    Mat4 result;
    for (int i = 0; i < VEC4_DIMENSIONS; ++i)
      result.m[i][i] = 1.0f;
    return result;
  }

  /** @brief Transforms a Vec4 by this matrix. */
  __attribute__((always_inline)) Vec4 apply(const Vec4 &input) const {
    return {{m[0][0] * input[0] + m[0][1] * input[1] + m[0][2] * input[2] +
                 m[0][3] * input[3],
             m[1][0] * input[0] + m[1][1] * input[1] + m[1][2] * input[2] +
                 m[1][3] * input[3],
             m[2][0] * input[0] + m[2][1] * input[1] + m[2][2] * input[2] +
                 m[2][3] * input[3],
             m[3][0] * input[0] + m[3][1] * input[1] + m[3][2] * input[2] +
                 m[3][3] * input[3]}};
  }
};

/**
 * @brief Right-multiplies `matrix` by a rotation of `angle` in the (a, b)
 *   coordinate plane.
 * @details Four dimensions have six such planes; composing rotations in them
 *   builds an arbitrary 4D orientation the way axis rotations build a 3D one.
 * @param matrix Transform to rotate in place.
 * @param a First plane axis, in [0, VEC4_DIMENSIONS).
 * @param b Second plane axis, in [0, VEC4_DIMENSIONS), distinct from `a`.
 * @param angle Rotation angle in radians.
 */
HS_COLD static inline void rotate_plane(Mat4 &matrix, int a, int b,
                                        float angle) {
  const float c = fast_cosf(angle);
  const float s = fast_sinf(angle);
  for (int column = 0; column < VEC4_DIMENSIONS; ++column) {
    const float av = matrix.m[a][column];
    const float bv = matrix.m[b][column];
    matrix.m[a][column] = c * av - s * bv;
    matrix.m[b][column] = s * av + c * bv;
  }
}
