/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/pattern/field.h"
#include "core/math/cube_projection.h"
#include <array>
#include <span>
#include <cstdint>
#include <cstring>

/** @file sphere_lut.h
 * @brief Reusable signed scalar cube storage for sphere fields. */

namespace Pattern {
/** @brief View over a baked cube-map sphere scalar field (values in [-1, 1]). */
struct SphereLutView {
  static constexpr uint16_t LAYOUT_VERSION =
      2; ///< Cube layout and quantization version.
  static constexpr int FACE_COUNT = 6;  ///< Cube faces.
  static constexpr int FACE_STEPS = 24; ///< Samples per face edge.
  static constexpr int FACE_SIZE = FACE_STEPS * FACE_STEPS; ///< Per face.
  static constexpr size_t SIZE = FACE_COUNT * FACE_SIZE;    ///< Entries.

  /// Borrowed face-major, row-major table of signed scalars in [-1, 1] times 127.
  const int8_t *data;
  bool active; ///< Whether `data` holds a baked table to sample.
  /// Float table on the same [-127, 127] scale. Takes precedence when set.
  const float *linear_data = nullptr;
};

/**
 * @brief Unit direction through a cube-face sample point.
 * @param face Cube face in [0, 6): +X, -X, +Y, -Y, +Z, -Z.
 * @param u Face-local horizontal coordinate in [-1, 1].
 * @param v Face-local vertical coordinate in [-1, 1].
 * @return The normalized direction.
 */
HS_HOT_INLINE inline SphereDirection sphere_face_coordinate(int face, float u,
                                                            float v) {
  const auto point = [&] {
    switch (face) {
    case 0:
      return math::Vector(1.0f, v, u);
    case 1:
      return math::Vector(-1.0f, v, -u);
    case 2:
      return math::Vector(u, 1.0f, v);
    case 3:
      return math::Vector(u, -1.0f, -v);
    case 4:
      return math::Vector(u, v, 1.0f);
    default:
      return math::Vector(-u, v, -1.0f);
    }
  };
  return SphereDirection::from(point());
}

/** @brief Unit vector through a cube face sample.
 * @param face Cube face index.
 * @param u Horizontal coordinate in [-1, 1].
 * @param v Vertical coordinate in [-1, 1].
 * @return Unit direction. */
inline math::Vector sphere_face_direction(int face, float u, float v) {
  return sphere_face_coordinate(face, u, v).value();
}

/** @brief Bilinear coordinates shared by cube-map scalar samples. */
struct SphereLutSample {
  uint16_t offset;      ///< Lower-left sample in the face-major point table.
  uint16_t cell_offset; ///< Corresponding cell in the packed face grid.
  float x_fraction;     ///< Horizontal interpolation fraction.
  float y_fraction;     ///< Vertical interpolation fraction.

  /** @brief Tests the four signed-byte taps used by this bilinear sample.
   * @param a First cube's face-major data.
   * @param b Second cube's face-major data.
   * @return Whether interpolation uses identical inputs. */
  HS_HOT_INLINE bool same_byte_corners(const int8_t *a, const int8_t *b) const {
    a += offset;
    b += offset;
    constexpr int stride = SphereLutView::FACE_STEPS;
    return std::memcmp(a, b, 2) == 0 &&
           std::memcmp(a + stride, b + stride, 2) == 0;
  }

  /** @brief Interpolates four taps in a face-major table.
   * @param table Float or signed-byte cube data.
   * @return Scalar value in [-1, 1]. */
  template <typename T> HS_HOT_INLINE float sample_data(const T *table) const {
    const auto *row = table + offset;
    const auto *next = row + SphereLutView::FACE_STEPS;
    const auto across = [&](const T *pair) HS_HOT_INLINE {
      // Signed-byte corner differences are exact before conversion to float.
      if constexpr (std::is_same_v<T, int8_t>)
        return static_cast<float>(pair[0]) +
               static_cast<float>(pair[1] - pair[0]) * x_fraction;
      else
        return hs::lerp(pair[0], pair[1], x_fraction);
    };
    const float low = across(row);
    const float high = across(next);
    return hs::lerp(low, high, y_fraction) * (1.0f / 127.0f);
  }
  /** @brief Samples a table using these cube-map coordinates.
   * @param view Prepared scalar cube LUT.
   * @return Scalar value in [-1, 1]. */
  HS_HOT_INLINE float sample(const SphereLutView &view) const {
    return view.linear_data ? sample_data(view.linear_data)
                            : sample_data(view.data);
  }
};

/** @brief Prepares bilinear cube-map coordinates within a face.
 * @param face Cube face index, [0, 6).
 * @param u Face coordinate in [-1, 1].
 * @param v Face coordinate in [-1, 1].
 * @return Cube-map coordinates for a bilinear sample. */
HS_HOT_INLINE inline SphereLutSample prepare_sphere_lut_face(int face, float u,
                                                             float v) {
  constexpr float SCALE =
      0.5f * static_cast<float>(SphereLutView::FACE_STEPS - 1);
  const float x = (u + 1.0f) * SCALE;
  const float y = (v + 1.0f) * SCALE;
  const int x_low =
      std::min(static_cast<int>(x), SphereLutView::FACE_STEPS - 2);
  const int y_low =
      std::min(static_cast<int>(y), SphereLutView::FACE_STEPS - 2);
  constexpr int N = SphereLutView::FACE_STEPS;
  static_assert(SphereLutView::SIZE <= UINT16_MAX);
  return {static_cast<uint16_t>((face * N + y_low) * N + x_low),
          static_cast<uint16_t>((face * (N - 1) + y_low) * (N - 1) + x_low),
          x - x_low, y - y_low};
}

/**
 * @brief Bilinear scalar tap in the cube field's face coordinates.
 * @param view Prepared scalar cube LUT.
 * @param face Cube face index, [0, 6).
 * @param u Face coordinate in [-1, 1].
 * @param v Face coordinate in [-1, 1].
 * @return Scalar value in [-1, 1].
 */
HS_HOT_INLINE inline float sample_sphere_lut_face(const SphereLutView &view,
                                                  int face, float u, float v) {
  return prepare_sphere_lut_face(face, u, v).sample(view);
}

/**
 * @brief Prepares cube-map scalar coordinates along an admitted direction.
 * @param direction Finite direction with positive float squared length;
 *          normalization is not required.
 * @return Cube-map coordinates for a bilinear sample.
 */
HS_O3_FN HS_HOT_INLINE inline SphereLutSample
prepare_sphere_lut_sample_unchecked(const math::Vector &direction) {
  const auto projection = math::project_cube_face(direction);
  const float inverse = 1.0f / projection.major;
  return prepare_sphere_lut_face(projection.face, projection.u * inverse,
                                 projection.v * inverse);
}

/** @brief Admits a raw direction and prepares its cube-map coordinates.
 * @param direction Finite nondegenerate direction; normalization is not required.
 * @return Cube-map coordinates for a bilinear sample. */
HS_O3_FN inline SphereLutSample
prepare_sphere_lut_sample(const math::Vector &direction) {
  HS_CHECK(SphereDirection::valid(direction),
           "sphere lookup requires an admitted direction");
  return prepare_sphere_lut_sample_unchecked(direction);
}

/** @brief Samples the cube-map scalar field along a nondegenerate direction.
 * @param view Prepared scalar cube LUT.
 * @param direction Finite nondegenerate direction; normalization is not required.
 * @return Scalar value in [-1, 1]. */
HS_O3_FN inline float sample_sphere_lut(const SphereLutView &view,
                                        const math::Vector &direction) {
  return prepare_sphere_lut_sample(direction).sample(view);
}

/** @brief Samples a cube field at an admitted direction.
 * @param view Borrowed cube table.
 * @param direction Admitted unit direction.
 * @return Interpolated signed scalar. */
inline float sample_sphere_lut(const SphereLutView &view,
                               const SphereDirection &direction) {
  return prepare_sphere_lut_sample_unchecked(direction.value()).sample(view);
}
/** @brief Prepares shared cube coordinates for an admitted direction.
 * @param direction Admitted unit direction.
 * @return Prepared cube taps. */
inline SphereLutSample
prepare_sphere_lut_sample(const SphereDirection &direction) {
  return prepare_sphere_lut_sample_unchecked(direction.value());
}
/** @brief Symmetric face-lattice coordinate for an edge index.
 * @param index Sample index in [0, FACE_STEPS).
 * @return Face coordinate in [-1, 1], with exact opposite-index symmetry. */
inline float sphere_face_grid_coordinate(int index) {
  constexpr float HALF_STEP = 1.0f / (SphereLutView::FACE_STEPS - 1);
  return (2 * index - (SphereLutView::FACE_STEPS - 1)) * HALF_STEP;
}

/** @brief Copies shared cube edges from their first face in X, Y, Z order. */
inline void
copy_sphere_lut_edges(std::span<int8_t, SphereLutView::SIZE> table) {
  constexpr int N = SphereLutView::FACE_STEPS;
  const auto at = [&](int face, int x, int y) -> int8_t & {
    return table[face * N * N + y * N + x];
  };
  for (int i = 0; i < N; ++i) {
    const int reverse = N - 1 - i;
    at(2, 0, i) = at(1, reverse, N - 1);
    at(2, N - 1, i) = at(0, i, N - 1);
    at(3, 0, i) = at(1, i, 0);
    at(3, N - 1, i) = at(0, reverse, 0);
    at(4, 0, i) = at(1, 0, i);
    at(4, N - 1, i) = at(0, N - 1, i);
    at(5, 0, i) = at(0, 0, i);
    at(5, N - 1, i) = at(1, N - 1, i);
  }
  for (int i = 1; i < N - 1; ++i) {
    at(4, i, 0) = at(3, i, 0);
    at(4, i, N - 1) = at(2, i, N - 1);
    at(5, i, 0) = at(3, N - 1 - i, N - 1);
    at(5, i, N - 1) = at(2, N - 1 - i, 0);
  }
}

/** @brief Bakes a typed source into signed-byte cube storage.
 * @tparam Source Sphere pattern source policy.
 * @param output Destination storage.
 * @param prepared Immutable prepared source inputs. */
template <class Source>
HS_FLASH_INLINE inline void
prepare_sphere_lut(std::span<int8_t, SphereLutView::SIZE> output,
                   const typename Source::Prepared &prepared) {
  constexpr int N = SphereLutView::FACE_STEPS;
  static_assert(N % 2 == 0);
  constexpr bool shifted_coordinates =
      requires { Source::configured_single(prepared); };
  const auto shift = [&](const math::Vector &direction) HS_HOT_INLINE {
    if constexpr (shifted_coordinates)
      return Source::shift_coordinate(direction, prepared);
    else
      return direction;
  };
  const bool single = [&] {
    if constexpr (shifted_coordinates)
      return Source::configured_single(prepared);
    else
      return false;
  }();
  const auto evaluate = [&](const math::Vector &coordinate) {
    if constexpr (shifted_coordinates)
      return Source::sample_coordinate(coordinate, prepared, single);
    else
      return Source::sample(SphereDirection::from_unit(coordinate), prepared);
  };
  for (int axis = 0; axis < 3; ++axis) {
    const int x_start = axis == 0 ? 0 : 1;
    for (int y = axis == 2 ? 1 : 0; y < N / 2; ++y) {
      for (int x = x_start; x < N / 2; ++x) {
        const auto v =
            sphere_face_coordinate(axis * 2, sphere_face_grid_coordinate(x),
                                   sphere_face_grid_coordinate(y))
                .value();
        const auto pos = shift(v);
        const auto neg = shift(-v);
        const auto write = [&](int column, const math::Vector &positive,
                               const math::Vector &negative) HS_HOT_INLINE {
          for (int side = 0; side < 4; ++side) {
            const bool odd = side & 1, bottom = side & 2;
            const math::Vector q(
                (odd && axis != 1 ? negative : positive).x,
                ((axis == 1 ? odd : bottom) ? negative : positive).y,
                (odd != (bottom && axis == 1) ? negative : positive).z);
            const int face = axis * 2 + odd;
            const int index =
                (face * N + (bottom ? N - 1 - y : y)) * N + column;
            const float sample = bounded_scalar(evaluate(q));
            output[index] = static_cast<int8_t>(sample * 127.0f +
                                                (sample < 0.0f ? -0.5f : 0.5f));
          }
        };
        write(x, pos, neg);
        if (axis == 0)
          write(N - 1 - x, {pos.x, pos.y, neg.z}, {neg.x, neg.y, pos.z});
        else
          write(N - 1 - x, {neg.x, pos.y, pos.z}, {pos.x, neg.y, neg.z});
      }
    }
  }
  copy_sphere_lut_edges(output);
}
/** @brief Baked sphere-field adapter for a typed source policy.
 * @tparam Source Sphere pattern source policy. */
template <class Source> struct BakedSphereField {
  using Coordinate = typename Source::Coordinate; ///< Typed source coordinate.
  using Prepared = SphereLutView;                 ///< Borrowed immutable table.
  /** @brief Bakes a source and returns its field view.
  * @param output Owned destination table.
  * @param source Prepared source inputs.
  * @return Borrowed prepared field. */
  static Prepared prepare(std::span<int8_t, SphereLutView::SIZE> output,
                          const typename Source::Prepared &source) {
    prepare_sphere_lut<Source>(output, source);
    return {output.data(), true};
  }
  /** @brief Samples a prepared cube field.
  * @param coordinate Admitted direction.
  * @param prepared Borrowed table.
  * @return Signed scalar sample. */
  static float sample(const Coordinate &coordinate, const Prepared &prepared) {
    return sample_sphere_lut(prepared, coordinate);
  }
};
/** @brief Cache for a complete caller-supplied source and layout identity.
 * @tparam Identity Equality-comparable record of every prepared source input. */
template <class Identity> struct SphereBakeCache {
  Identity identity{};                 ///< Complete source identity.
  const int8_t *destination = nullptr; ///< Storage associated with this bake.
  /** @brief Invalidates retained derived storage. */
  void reset() { destination = nullptr; }
  /** @brief Whether the table matches a complete source identity.
  * @param output Destination address.
  * @param key Complete prepared-input identity.
  * @return True if no rebake is needed. */
  bool holds(const int8_t *output, const Identity &key) const {
    return destination && destination == output && identity == key;
  }
  /** @brief Rebakes after any identity or ownership change.
  * @tparam Source Typed sphere source policy.
  * @param output Destination table.
  * @param key Complete source identity.
  * @param source Prepared source inputs.
  * @return Whether the table was rebuilt. */
  template <class Source>
  bool refresh(std::span<int8_t, SphereLutView::SIZE> output,
               const Identity &key, const typename Source::Prepared &source) {
    if (holds(output.data(), key))
      return false;
    prepare_sphere_lut<Source>(output, source);
    identity = key;
    destination = output.data();
    return true;
  }
};
} // namespace Pattern
