/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once
#include <cmath>
#include <cstdint>
#include <cstring>
#include "math/3dmath.h"
#include "math/pixel_mapping.h"
#include "math/spherical_field.h"

/**
 * @file feedback_cap_plane.h
 * @brief Filter::Pixel::FeedbackCapPlane: the pole cap-plane coordinates the
 * feedback filter interpolates polar warp offsets in.
 */

namespace Filter {

namespace Pixel {

/**
 * @brief Cap-plane encoding of polar warp offsets and their conversion back
 * to field coordinates.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details A cap plane lays each direction's angle from a pole along its
 * longitude, so offsets near the pole stay bounded where equirect longitude
 * offsets blow up as 1/sin(phi). Offsets are stored as CapOffset in
 * CAP_SCALE units per radian.
 */
template <int W, int H> struct FeedbackCapPlane {
  /// Spherical field layout of the W x H canvas.
  using SphereField = hs::SphericalFieldLayout<W, H>;

  /** @brief Cap-plane offset units per radian. */
  static constexpr float CAP_SCALE = 8192.0f;

  /** @brief An offset in a pole's cap plane, CAP_SCALE units per radian. */
  struct CapOffset {
    int16_t u; ///< Offset along world x, CAP_SCALE units per radian.
    int16_t v; ///< Offset along world z, CAP_SCALE units per radian.
  };

  /** @brief A point in a pole's cap plane: the angle from that pole, in
   *  radians, laid along the point's longitude. */
  struct CapPoint {
    float u; ///< Component along world x, in radians.
    float v; ///< Component along world z, in radians.
  };

  /** @brief A polar cell's cap-plane offset at its left edge and its change
   *  across the cell, both blended between the cell's two rings. */
  struct CapCell {
    CapPoint left;  ///< Offset at the cell's left column, in radians.
    CapPoint slope; ///< Change from left to right column, in radians.
  };

  /** @brief Cap-plane coordinates of a direction, from the north pole or,
   *  with @p south, from the south pole. */
  __attribute__((noinline)) static CapPoint cap_point(const math::Vector &v,
                                                      bool south) {
    const float horizontal = sqrtf(v.x * v.x + v.z * v.z);
    if (!(horizontal > 1e-9f))
      return {0.0f, 0.0f};
    const float scale =
        math::precise_atan2(horizontal, south ? -v.y : v.y) / horizontal;
    return {v.x * scale, v.z * scale};
  }

  /** @brief Quantizes a cap-plane offset in radians, saturating at the
   *  int16_t range. */
  static CapOffset encode_cap(float u, float v) {
    auto quantize = [](float c) {
      const float scaled = hs::clamp(c * CAP_SCALE, -32767.0f, 32767.0f);
      return static_cast<int16_t>(scaled + (scaled < 0.0f ? -0.5f : 0.5f));
    };
    return {quantize(u), quantize(v)};
  }

  /** @brief Encoded offsets @p a and @p b blended by @p mix, truncated. */
  static __attribute__((always_inline)) CapOffset
  lerp_offset(const CapOffset &a, const CapOffset &b, float mix) {
    return {static_cast<int16_t>(hs::lerp(static_cast<float>(a.u),
                                          static_cast<float>(b.u), mix)),
            static_cast<int16_t>(hs::lerp(static_cast<float>(a.v),
                                          static_cast<float>(b.v), mix))};
  }

  /**
   * @brief Decodes the cell between columns @p cx0 and @p cx1 of two ring
   * rows of encoded offsets, blended by the row weights.
   * @param caps0 Encoded offsets of the upper ring.
   * @param caps1 Encoded offsets of the lower ring.
   * @param cx0 Left column.
   * @param cx1 Right column.
   * @param wy0 Weight of @p caps0.
   * @param wy1 Weight of @p caps1.
   */
  static __attribute__((always_inline)) CapCell
  decode_cell(const CapOffset *caps0, const CapOffset *caps1, int cx0, int cx1,
              float wy0, float wy1) {
    constexpr float INVERSE_CAP = 1.0f / CAP_SCALE;
    const float left_u =
        (caps0[cx0].u * wy0 + caps1[cx0].u * wy1) * INVERSE_CAP;
    const float left_v =
        (caps0[cx0].v * wy0 + caps1[cx0].v * wy1) * INVERSE_CAP;
    return {{left_u, left_v},
            {(caps0[cx1].u * wy0 + caps1[cx1].u * wy1) * INVERSE_CAP - left_u,
             (caps0[cx1].v * wy0 + caps1[cx1].v * wy1) * INVERSE_CAP - left_v}};
  }

  /**
   * @brief Source coordinates of one lane of a polar row.
   * @param cell The lane's decoded cell.
   * @param x The lane's column.
   * @param cap_angle The row's angle from its pole, in radians.
   * @param north Whether the row's pole is the north pole.
   * @param midpoint Whether the lane stands for its column pair.
   * @param fx The lane's fraction across @p cell.
   */
  static __attribute__((always_inline)) typename SphereField::Coordinates
  polar_lane(const CapCell &cell, int x, float cap_angle, bool north,
             bool midpoint, float fx) {
    const CapPoint p = cap_lane(cell, x, cap_angle, midpoint, fx);
    const AtanTerms a = atan_terms(p.v, p.u);
    return {atan_fold(a.numerator / a.denominator, a) *
                (W / (2.0f * math::PI_F)),
            cap_row(p, north)};
  }

private:
  /** @brief One arctangent's ratio terms and the octant folds it needs. */
  struct AtanTerms {
    float numerator;
    float denominator;
    bool steep;
    bool negative_x;
    bool negative_y;
  };

  static __attribute__((always_inline)) AtanTerms atan_terms(float y, float x) {
    uint32_t y_bits, x_bits;
    std::memcpy(&y_bits, &y, sizeof(y_bits));
    std::memcpy(&x_bits, &x, sizeof(x_bits));
    // Sign-cleared bit patterns order like the magnitudes.
    const bool steep = (y_bits & 0x7fffffffu) > (x_bits & 0x7fffffffu);
    const float abs_y = fabsf(y), abs_x = fabsf(x);
    // A floor on the denominator keeps a pair's shared product normal; only
    // the undefined longitude at an exact pole reaches it.
    return {steep ? abs_x : abs_y, fmaxf(steep ? abs_y : abs_x, 1e-6f), steep,
            (x_bits >> 31) != 0, (y_bits >> 31) != 0};
  }

  static __attribute__((always_inline)) float atan_fold(float ratio,
                                                        const AtanTerms &t) {
    float angle = math::atan_unit(ratio);
    if (t.steep)
      angle = 1.57079633f - angle;
    if (t.negative_x)
      angle = 3.14159265f - angle;
    return t.negative_y ? -angle : angle;
  }

  /**
   * @brief A polar lane's target in the cap plane: its own column at the
   *        row's cap angle plus the cell's offset at @p fx.
   * @details A half-resolution lane stands for its column pair, so its
   * column turns half a column past the even column.
   */
  static __attribute__((always_inline)) CapPoint cap_lane(
      const CapCell &cell, int x, float cap_angle, bool midpoint, float fx) {
    constexpr float HALF = math::PI_F / W;
    constexpr float HALF_COS =
        1.0f - HALF * HALF * 0.5f + HALF * HALF * HALF * HALF * (1.0f / 24.0f);
    constexpr float HALF_SIN =
        HALF - HALF * HALF * HALF * (1.0f / 6.0f) +
        HALF * HALF * HALF * HALF * HALF * (1.0f / 120.0f);
    float cosine = math::TrigLUT<W, H>::cos_theta(x);
    float sine = math::TrigLUT<W, H>::sin_theta[x];
    if (midpoint) {
      const float turned = cosine * HALF_COS - sine * HALF_SIN;
      sine = cosine * HALF_SIN + sine * HALF_COS;
      cosine = turned;
    }
    return {cap_angle * cosine + cell.left.u + cell.slope.u * fx,
            cap_angle * sine + cell.left.v + cell.slope.v * fx};
  }

  /** @brief Field row of a cap angle from the north or the south pole. */
  static __attribute__((always_inline)) float cap_row(const CapPoint &p,
                                                      bool north) {
    const float length_sq = fmaxf(p.u * p.u + p.v * p.v, 1e-12f);
    const float angle = length_sq * math::fast_rsqrt(length_sq);
    return SphereField::Geometry::phi_to_row(north ? angle
                                                   : math::PI_F - angle);
  }
};

} // namespace Pixel

} // namespace Filter
