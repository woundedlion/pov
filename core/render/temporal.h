/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file temporal.h
 * @brief Column-time sampling contract shared by every renderer family: the
 * fidelity toggle, the sweep fraction of a column, sweep slabs, the shortest
 * rotation arc between frame endpoints, the column-time state scans read, and
 * capture-time solving for plotted points.
 */

#include "core/animation/orientation.h"
#include "core/math/3dmath.h"
#include "core/platform/attributes.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <span>

namespace Temporal {

#if defined(__EMSCRIPTEN__) || HS_ENABLE_TEST_HOOKS
namespace Detail {
inline bool enabled = true; ///< Simulator sweep-time sampling selection.
}

/** @brief Selects sweep-time sampling for simulator presentation.
 * @param enabled Whether columns use their individual emission times. */
inline void set_temporal_sampling_enabled(bool enabled) {
  Detail::enabled = enabled;
}

/** @brief Whether columns use their individual emission times.
 * @return The simulator's temporal sampling mode. */
inline bool temporal_sampling_enabled() { return Detail::enabled; }
#else
/** @brief Whether columns use their individual emission times.
 * @return True for physical display rendering. */
constexpr bool temporal_sampling_enabled() { return true; }
#endif

/** @brief Column time within a two-arm half-turn.
 * @param column Absolute canvas column.
 * @param width Even canvas width.
 * @return Fraction of one effect frame, in [0, 1). */
constexpr float sweep_fraction(int column, int width) {
  return static_cast<float>(column % (width / 2)) / (width / 2);
}

/** @brief Shortest quaternion arc with frame-invariant trigonometry. */
struct TemporalRotation {
  math::Quaternion from; ///< Sign-adjusted start orientation.
  math::Quaternion to;   ///< End orientation.
  float theta = 0.0f;    ///< Arc angle, or zero for normalized linear sampling.
  float sine = 1.0f;     ///< Sine of the arc angle.

  /** @brief Prepares the shortest arc between unit orientations.
   * @param a Previous unit orientation.
   * @param b Current unit orientation. */
  HS_FLASH_INLINE TemporalRotation(const math::Quaternion &a,
                                   const math::Quaternion &b)
      : from(a), to(b) {
    assert(!(std::fabs(math::dot(a, a) - 1.0f) > math::EPS_UNIT_QUAT_SQ));
    assert(!(std::fabs(math::dot(b, b) - 1.0f) > math::EPS_UNIT_QUAT_SQ));
    float d = math::dot(a, b);
    if (d < 0.0f) {
      from = -from;
      d = -d;
    }
    if (d <= 1.0f - math::TOLERANCE) {
      theta = acosf(hs::clamp(d, -1.0f, 1.0f));
      sine = sinf(theta);
    }
  }

  /** @brief Samples an interior fraction of the prepared arc.
   * @param fraction Interior progress in (0, 1).
   * @return Unit orientation along the shortest arc. */
  __attribute__((always_inline)) math::Quaternion sample(float fraction) const {
    if (theta == 0.0f)
      return (from + fraction * (to - from)).normalized();
    const float A = sinf((1.0f - fraction) * theta) / sine;
    const float B = sinf(fraction * theta) / sine;
    return (A * from + B * to).normalized();
  }
};

/**
 * @brief Column-time state shared by the scans of one effect frame.
 * @details A scan with a basis builds a pixel's direction from its column
 * `c = x % (W / 2)` as `±sin(phi) * azimuth[c] + cos(phi) * up[c]`, negated on
 * the far arm: the fixed-meridian direction rotated from the column's sweep
 * time to the current endpoint, so shapes built at the current endpoint are
 * evaluated as they stood when the column was lit.
 */
struct ColumnTime {
  /// Rotated azimuth directions of the first arm, W / 2 entries; null for the
  /// fixed meridian.
  const math::Vector *azimuth = nullptr;
  /// Per-fraction rotated up directions, W / 2 entries.
  const math::Vector *up = nullptr;
  /// Bound on the arc any direction moves over the frame (radians).
  float sweep = 0.0f;
  /// Earliest frame fraction drawn; columns lit before it are left untouched.
  float birth = 0.0f;
};

/**
 * @brief First half-turn column at or after a frame fraction.
 * @param fraction Frame fraction in [0, 1].
 * @param width Even canvas width.
 * @return Smallest column `c` in [0, width / 2] with
 *         `sweep_fraction(c, width) >= fraction`, or width / 2 when none is.
 */
HS_FLASH_MEMBER inline int first_column_at(float fraction, int width) {
  const int half = width / 2;
  int column = static_cast<int>(
      std::ceil(hs::clamp(fraction, 0.0f, 1.0f) * static_cast<float>(half)));
  while (column > 0 && sweep_fraction(column - 1, width) >= fraction)
    --column;
  while (column < half && sweep_fraction(column, width) < fraction)
    ++column;
  return column;
}

/** @brief Half-turn columns of one sweep slab. */
struct SweepSlab {
  int begin; ///< First half-turn column.
  int end;   ///< One past the last half-turn column.
};

/** @brief Columns of one of `count` slabs that tile the half-turn in sweep
 * order, for a render-x window.
 * @param slab Slab index, in [0, count).
 * @param count Slab count, in [1, width / 2].
 * @param width Even canvas width.
 * @return Slab `slab`, whose columns' sweep fractions lie in
 * `[slab / count, (slab + 1) / count)` up to column rounding. */
constexpr SweepSlab sweep_slab(int slab, int count, int width) {
  return {slab * (width / 2) / count, (slab + 1) * (width / 2) / count};
}

/** @brief Most history nodes capture() accepts. */
inline constexpr int MAX_CAPTURE_NODES = 16;

/** @brief Frame fractions at which a display arm passes a moving point. */
struct Captures {
  static constexpr int CAPACITY = 2;      ///< Most arm passes in one frame.
  std::array<float, CAPACITY> fraction{}; ///< Frame fraction per pass.
  int count = 0;                          ///< Live passes.
};

/** @brief Azimuth of a point in half-turns, so that a column's sweep
 * fraction is this value modulo one.
 * @param p Point on the sphere.
 * @return Azimuth divided by pi, in [-1, 1]. */
__attribute__((always_inline)) inline float half_turns(const math::Vector &p) {
  return math::fast_atan2(p.z, p.x) / math::PI_F;
}

/**
 * @brief Solves `u = (theta(u) / pi) mod 1` for a point sampled at evenly
 *        spaced frame fractions.
 * @details Azimuth is linear between nodes. The arm-relative position
 * `theta / pi - u` must fall strictly from node to node; each integer it
 * crosses is one pass, located by binary search over the nodes.
 * @param nodes The point at fractions `i / (size - 1)`; at least two and at
 *        most `MAX_CAPTURE_NODES`, unit length.
 * @return Zero, one or two passes with u in [0, 1). A rise between nodes,
 * from a swing near a pole or motion with the sweep faster than the arm, or a
 * fall of two or more, returns one pass at u = 1.
 */
HS_FLASH_MEMBER inline Captures capture(std::span<const math::Vector> nodes) {
  assert(nodes.size() >= 2 && nodes.size() <= MAX_CAPTURE_NODES);
  const int last = static_cast<int>(nodes.size()) - 1;
  std::array<float, MAX_CAPTURE_NODES> relative;
  float azimuth = half_turns(nodes[0]);
  float unwrapped = azimuth;
  relative[0] = unwrapped;
  Captures out;
  for (int i = 1; i <= last; ++i) {
    const float next = half_turns(nodes[i]);
    unwrapped += std::remainder(next - azimuth, 2.0f);
    azimuth = next;
    relative[i] = unwrapped - static_cast<float>(i) / last;
    if (!(relative[i] < relative[i - 1])) {
      out.fraction[0] = 1.0f;
      out.count = 1;
      return out;
    }
  }
  if (!(relative[0] - relative[last] < Captures::CAPACITY)) {
    out.fraction[0] = 1.0f;
    out.count = 1;
    return out;
  }
  const auto begin = relative.begin();
  const auto end = begin + last + 1;
  for (float n = std::floor(relative[0]); n > relative[last]; n -= 1.0f) {
    const int j = static_cast<int>(
        std::partition_point(begin, end, [n](float r) { return r >= n; }) -
        begin);
    const int i = j - 1;
    const float t = (relative[i] - n) / (relative[i] - relative[j]);
    out.fraction[out.count++] = (static_cast<float>(i) + t) / last;
  }
  return out;
}

/**
 * @brief Arm passes of a point carried along an orientation's sub-frame path.
 * @tparam CAP Orientation capacity.
 * @param orientation History of at least two frames spanning one frame.
 * @param v Unrotated point, unit length.
 * @return capture() over @p v rotated by each history frame.
 */
template <int CAP>
HS_FLASH_MEMBER Captures capture(const math::Orientation<CAP> &orientation,
                                 const math::Vector &v) {
  static_assert(CAP <= MAX_CAPTURE_NODES);
  std::array<math::Vector, CAP> nodes;
  for (int i = 0; i < orientation.length(); ++i)
    nodes[i] = math::rotate(v, orientation.get(i));
  return capture(std::span<const math::Vector>(
      nodes.data(), static_cast<size_t>(orientation.length())));
}

} // namespace Temporal
