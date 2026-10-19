/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include "core/animation/orientation.h"
#include "core/math/3dmath.h"
#include "core/memory.h"
#include "core/render/canvas.h"
#include "core/render/sdf.h"
#include "core/render/temporal.h"

/**
 * @file column_basis.h
 * @brief Column-time rigid rotation for the scan loops: the per-column pixel
 * basis, its frame installation, and swept row bounds.
 */

namespace Scan {

/**
 * @brief Column-time state a scan through a pipeline honours.
 * @param pipeline Destination pipeline, typed or erased.
 * @param canvas Destination canvas.
 * @return The canvas's installed state, or null when none is installed or the
 *         pipeline moves world positions after the scan.
 */
template <typename PipelineT>
__attribute__((always_inline)) inline const Temporal::ColumnTime *
column_time(const PipelineT &pipeline, const Canvas &canvas) {
  if constexpr (requires { pipeline.world_transform_is_identity; }) {
    if (!pipeline.world_transform_is_identity)
      return nullptr;
  }
  return canvas.column_time();
}

/**
 * @brief Unit direction a scan evaluates at a pixel.
 * @tparam W Canvas width in columns.
 * @param azimuth Rotated azimuth directions of the first arm's W / 2 columns,
 *        or null for the fixed meridian.
 * @param up Rotated up directions, W / 2 entries; read only with @p azimuth.
 * @param x Column.
 * @param sin_phi Sine of the row's colatitude.
 * @param cos_phi Cosine of the row's colatitude.
 * @param cos_theta Column cosine table.
 * @param sin_theta Column sine table.
 * @return The fixed-meridian direction, or that direction rotated from the
 *         column's sweep time to the current endpoint.
 */
template <int W>
__attribute__((always_inline)) inline math::Vector
column_direction(const math::Vector *azimuth, const math::Vector *up, int x,
                 float sin_phi, float cos_phi, const float *cos_theta,
                 const float *sin_theta) {
  if (!azimuth)
    return math::Vector(sin_phi * cos_theta[x], cos_phi,
                        sin_phi * sin_theta[x]);
  const bool far_arm = x >= W / 2;
  const int column = far_arm ? x - W / 2 : x;
  const float s = far_arm ? -sin_phi : sin_phi;
  const math::Vector &e = azimuth[column];
  const math::Vector &v = up[column];
  return math::Vector(s * e.x + cos_phi * v.x, s * e.y + cos_phi * v.y,
                      s * e.z + cos_phi * v.z);
}

/**
 * @brief Widens a row band to cover a shape swept by up to an arc.
 * @tparam H Canvas height in rows.
 * @param bounds Rows covering the shape at the current endpoint.
 * @param sweep Largest arc any point moves (radians); 0 leaves the band.
 * @return Rows covering every point within @p sweep of the band.
 */
template <int H>
inline SDF::Bounds swept_rows(const SDF::Bounds &bounds, float sweep) {
  if (!(sweep > 0.0f) || bounds.y_min > bounds.y_max)
    return bounds;
  const math::LatitudeGeometry geometry(H);
  return SDF::phi_bounds_to_rows<H>(
      geometry.row_to_phi(static_cast<float>(bounds.y_min)) - sweep,
      geometry.row_to_phi(static_cast<float>(bounds.y_max)) + sweep);
}

/**
 * @brief Per-column pixel basis of a rigid rotation over one effect frame.
 * @tparam W Canvas width in columns.
 * @tparam H Canvas height in rows.
 * @details Shapes are built at the current endpoint. Column `x`, lit at
 * fraction `u = sweep_fraction(x, W)`, evaluates them along
 * `current * conj(slerp(previous, current, u))` applied to its fixed-meridian
 * direction, which is where the shape stood when the column was lit. Opposing
 * columns share a fraction, so one entry serves both arms. The table is
 * frame-scoped: an effect allocates it in a scratch arena for its draw.
 */
template <int W, int H> struct ColumnBasis {
  static_assert(W % 4 == 0, "column basis reads the quarter-turn trig table");
  static constexpr int HALF = W / 2; ///< Distinct sweep fractions.

  std::array<math::Vector, HALF> azimuth; ///< Rotated first-arm azimuths.
  std::array<math::Vector, HALF> up;      ///< Rotated up per sweep fraction.

  /**
   * @brief Prepares the frame's basis from its orientation endpoints.
   * @param previous Orientation at the previous endpoint, unit length.
   * @param current Orientation at the current endpoint, unit length.
   * @return Column-time state viewing this table, or the fixed meridian with
   *         temporal sampling disabled or without motion.
   * @details An orientation's history supplies both endpoints; one without
   * history, as after construction or a snap, has no motion.
   */
  HS_FLASH_MEMBER Temporal::ColumnTime
  prepare(const math::Quaternion &previous, const math::Quaternion &current) {
    Temporal::ColumnTime time;
    if (!moving(previous, current))
      return time;
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    const float *cos_theta = math::TrigLUT<W, H>::cos_theta_data();
    const float *sin_theta = math::TrigLUT<W, H>::sin_theta.data();
    const Temporal::TemporalRotation path(previous, current);
    for (int i = 0; i < HALF; ++i) {
      const math::Quaternion delta =
          current * path.sample(Temporal::sweep_fraction(i, W)).conjugate();
      up[i] = math::rotate(math::Vector(0.0f, 1.0f, 0.0f), delta);
      azimuth[i] =
          math::rotate(math::Vector(cos_theta[i], 0.0f, sin_theta[i]), delta);
    }
    const math::Quaternion step = current * path.from.conjugate();
    const float half_sine = sqrtf(math::dot(step.v, step.v));
    time.azimuth = azimuth.data();
    time.up = up.data();
    time.sweep = 2.0f * asinf(fminf(half_sine, 1.0f)) + SWEEP_PAD;
    return time;
  }

  /// Arc added to the sweep bound for the slerp's rounding.
  static constexpr float SWEEP_PAD = 1e-4f;

  /**
   * @brief Whether a frame's endpoints call for a rotated basis.
   * @param previous Orientation at the previous endpoint.
   * @param current Orientation at the current endpoint.
   * @return True with temporal sampling enabled and distinct endpoints.
   */
  static bool moving(const math::Quaternion &previous,
                     const math::Quaternion &current) {
    return Temporal::temporal_sampling_enabled() &&
           !(previous.r == current.r && previous.v.x == current.v.x &&
             previous.v.y == current.v.y && previous.v.z == current.v.z);
  }
};

/**
 * @brief Prepares a frame's column basis from its orientation endpoints.
 * @tparam W Canvas width in columns.
 * @tparam H Canvas height in rows.
 * @param arena Arena holding the table when the frame moves; the table lives
 *        until the caller's scope on it closes.
 * @param previous Orientation at the previous endpoint.
 * @param current Orientation at the current endpoint.
 * @return Column-time state for the frame's scans.
 */
template <int W, int H>
inline Temporal::ColumnTime
prepare_column_basis(Arena &arena, const math::Quaternion &previous,
                     const math::Quaternion &current) {
  if (!ColumnBasis<W, H>::moving(previous, current))
    return {};
  return arena.make_default_initialized<ColumnBasis<W, H>>()->prepare(previous,
                                                                      current);
}

/**
 * @brief Prepares an orientation's column basis for this frame.
 * @tparam W Canvas width in columns.
 * @tparam H Canvas height in rows.
 * @tparam CAP Orientation history capacity.
 * @param arena Arena holding the table when the frame moves; the table lives
 *        until the caller's scope on it closes.
 * @param orientation Orientation whose history spans the frame: index 0 is
 *        the previous endpoint and the last index the current one.
 * @return Column-time state for the frame's scans.
 */
template <int W, int H, int CAP>
inline Temporal::ColumnTime
prepare_column_basis(Arena &arena, const math::Orientation<CAP> &orientation) {
  return prepare_column_basis<W, H>(arena, orientation.get(0),
                                    orientation.get());
}

/**
 * @brief Installs column-time state on a canvas for a scope.
 * @details Restores the previously installed state on exit, so scopes nest.
 */
class ColumnTimeScope {
public:
  /**
   * @brief Installs @p state.
   * @param canvas Canvas whose later scans read the state.
   * @param state State to install; must outlive the scope.
   */
  ColumnTimeScope(Canvas &canvas, const Temporal::ColumnTime &state)
      : canvas(canvas), saved(canvas.column_time()) {
    canvas.set_column_time(&state);
  }
  ~ColumnTimeScope() { canvas.set_column_time(saved); }
  ColumnTimeScope(const ColumnTimeScope &) = delete;
  ColumnTimeScope &operator=(const ColumnTimeScope &) = delete;

private:
  Canvas &canvas;                    ///< Canvas carrying the state.
  const Temporal::ColumnTime *saved; ///< State installed before the scope.
};

} // namespace Scan
