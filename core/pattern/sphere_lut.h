/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/pattern/field.h"
#include <span>
#include <cstdint>

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
  const float *linear_data =
      nullptr; ///< Exact float copies of quantized samples.
};

/**
 * @brief Unit direction through a cube-face sample point.
 * @param face Cube face in [0, 6): +X, -X, +Y, -Y, +Z, -Z.
 * @param u Face-local horizontal coordinate in [-1, 1].
 * @param v Face-local vertical coordinate in [-1, 1].
 * @return The normalized direction.
 */
__attribute__((always_inline)) inline SphereDirection
sphere_face_coordinate(int face, float u, float v) {
  switch (face) {
  case 0:
    return SphereDirection::from(math::Vector(1.0f, v, u));
  case 1:
    return SphereDirection::from(math::Vector(-1.0f, v, -u));
  case 2:
    return SphereDirection::from(math::Vector(u, 1.0f, v));
  case 3:
    return SphereDirection::from(math::Vector(u, -1.0f, -v));
  case 4:
    return SphereDirection::from(math::Vector(u, v, 1.0f));
  default:
    return SphereDirection::from(math::Vector(-u, v, -1.0f));
  }
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
  int offset;       ///< Lower-left tap's offset within the table.
  float x_fraction; ///< Horizontal interpolation fraction.
  float y_fraction; ///< Vertical interpolation fraction.

  /** @brief Samples a table using these cube-map coordinates.
   * @param view Prepared scalar cube LUT.
   * @return Scalar value in [-1, 1]. */
  __attribute__((always_inline)) float sample(const SphereLutView &view) const {
    float a;
    float b;
    float c;
    float d;
    if (view.linear_data) {
      a = view.linear_data[offset];
      b = view.linear_data[offset + 1];
      c = view.linear_data[offset + SphereLutView::FACE_STEPS];
      d = view.linear_data[offset + SphereLutView::FACE_STEPS + 1];
    } else {
      a = static_cast<float>(view.data[offset]);
      b = static_cast<float>(view.data[offset + 1]);
      c = static_cast<float>(view.data[offset + SphereLutView::FACE_STEPS]);
      d = static_cast<float>(view.data[offset + SphereLutView::FACE_STEPS + 1]);
    }
    const float row_low = hs::lerp(a, b, x_fraction);
    const float row_high = hs::lerp(c, d, x_fraction);
    return hs::lerp(row_low, row_high, y_fraction) * (1.0f / 127.0f);
  }
};

/** @brief Prepares bilinear cube-map coordinates within a face.
 * @param face Cube face index, [0, 6).
 * @param u Face coordinate in [-1, 1].
 * @param v Face coordinate in [-1, 1].
 * @return Cube-map coordinates for a bilinear sample. */
__attribute__((always_inline)) inline SphereLutSample
prepare_sphere_lut_face(int face, float u, float v) {
  constexpr float SCALE =
      0.5f * static_cast<float>(SphereLutView::FACE_STEPS - 1);
  const float X = (u + 1.0f) * SCALE;
  const float Y = (v + 1.0f) * SCALE;
  const int X_LOW =
      std::min(static_cast<int>(X), SphereLutView::FACE_STEPS - 2);
  const int Y_LOW =
      std::min(static_cast<int>(Y), SphereLutView::FACE_STEPS - 2);
  return {face * SphereLutView::FACE_SIZE + Y_LOW * SphereLutView::FACE_STEPS +
              X_LOW,
          X - X_LOW, Y - Y_LOW};
}

/**
 * @brief Bilinear scalar tap in the cube field's face coordinates.
 * @param view Prepared scalar cube LUT.
 * @param face Cube face index, [0, 6).
 * @param u Face coordinate in [-1, 1].
 * @param v Face coordinate in [-1, 1].
 * @return Scalar value in [-1, 1].
 */
__attribute__((always_inline)) inline float
sample_sphere_lut_face(const SphereLutView &view, int face, float u, float v) {
  return prepare_sphere_lut_face(face, u, v).sample(view);
}

/**
 * @brief Prepares cube-map scalar coordinates along a nonzero direction.
 * @param direction Finite direction with positive float squared length;
 *          normalization is not required.
 * @return Cube-map coordinates for a bilinear sample.
 */
HS_O3_FN inline SphereLutSample
prepare_sphere_lut_sample(const math::Vector &direction) {
  HS_CHECK(SphereDirection::valid(direction),
           "sphere lookup requires an admitted direction");
  const float ax = fabsf(direction.x);
  const float ay = fabsf(direction.y);
  const float az = fabsf(direction.z);
  int face;
  float u;
  float v;
  if (ax >= ay && ax >= az) {
    const float inverse = 1.0f / ax;
    face = direction.x >= 0.0f ? 0 : 1;
    u = (direction.x >= 0.0f ? direction.z : -direction.z) * inverse;
    v = direction.y * inverse;
  } else if (ay >= az) {
    const float inverse = 1.0f / ay;
    face = direction.y >= 0.0f ? 2 : 3;
    u = direction.x * inverse;
    v = (direction.y >= 0.0f ? direction.z : -direction.z) * inverse;
  } else {
    const float inverse = 1.0f / az;
    face = direction.z >= 0.0f ? 4 : 5;
    u = (direction.z >= 0.0f ? direction.x : -direction.x) * inverse;
    v = direction.y * inverse;
  }

  return prepare_sphere_lut_face(face, u, v);
}

/** @brief Samples the cube-map scalar field along a nonzero direction.
 * @param view Prepared scalar cube LUT.
 * @param direction Finite nonzero direction; normalization is not required.
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
  return sample_sphere_lut(view, direction.value());
}
/** @brief Prepares shared cube coordinates for an admitted direction.
 * @param direction Admitted unit direction.
 * @return Prepared cube taps. */
inline SphereLutSample
prepare_sphere_lut_sample(const SphereDirection &direction) {
  return prepare_sphere_lut_sample(direction.value());
}
/** @brief Symmetric face-lattice coordinate for an edge index.
 * @param index Sample index in [0, FACE_STEPS).
 * @return Face coordinate in [-1, 1], with exact opposite-index symmetry. */
inline float sphere_face_grid_coordinate(int index) {
  constexpr float HALF_STEP = 1.0f / (SphereLutView::FACE_STEPS - 1);
  return (2 * index - (SphereLutView::FACE_STEPS - 1)) * HALF_STEP;
}

/** @brief Bakes a typed source into signed-byte cube storage.
 * @tparam Source Sphere pattern source policy.
 * @param output Destination storage.
 * @param prepared Immutable prepared source inputs. */
template <class Source>
HS_FLASH_INLINE inline void
prepare_sphere_lut(std::span<int8_t, SphereLutView::SIZE> output,
                   const typename Source::Prepared &prepared) {
  for (int face = 0; face < SphereLutView::FACE_COUNT; ++face)
    for (int y = 0; y < SphereLutView::FACE_STEPS; ++y)
      for (int x = 0; x < SphereLutView::FACE_STEPS; ++x) {
        const SphereDirection D =
            sphere_face_coordinate(face, sphere_face_grid_coordinate(x),
                                   sphere_face_grid_coordinate(y));
        const float SAMPLE = bounded_scalar(Source::sample(D, prepared));
        output[face * SphereLutView::FACE_SIZE + y * SphereLutView::FACE_STEPS +
               x] = static_cast<int8_t>(SAMPLE * 127.0f +
                                        (SAMPLE < 0.0f ? -0.5f : 0.5f));
      }
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
  static constexpr uint16_t LAYOUT_VERSION =
      SphereLutView::LAYOUT_VERSION; ///< Cube layout and quantization version.
  Identity identity{};               ///< Complete source identity.
  const int8_t *destination = nullptr; ///< Storage associated with this bake.
  uint16_t layout_version = 0; ///< Baked layout version; zero when invalid.
  /** @brief Invalidates retained derived storage. */
  void reset() {
    destination = nullptr;
    layout_version = 0;
  }
  /** @brief Whether the table matches a complete source identity.
  * @param output Destination address.
  * @param key Complete prepared-input identity.
  * @return True if no rebake is needed. */
  bool holds(const int8_t *output, const Identity &key) const {
    return destination == output && layout_version == LAYOUT_VERSION &&
           identity == key;
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
    layout_version = LAYOUT_VERSION;
    return true;
  }
};
} // namespace Pattern
