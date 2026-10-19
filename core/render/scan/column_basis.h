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
 * basis, its object-space counterpart, its frame installation, and swept row
 * bounds.
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
 * @brief Direction a column's rotated basis gives a pixel.
 * @tparam W Canvas width in columns.
 * @param azimuth Rotated azimuth directions of the first arm's W / 2 columns.
 * @param up Rotated up directions, W / 2 entries.
 * @param x Column.
 * @param sin_phi Sine of the row's colatitude.
 * @param cos_phi Cosine of the row's colatitude.
 * @return The fixed-meridian direction rotated from the column's sweep time to
 *         the current endpoint.
 */
template <int W>
__attribute__((always_inline)) inline math::Vector
rotated_direction(const math::Vector *azimuth, const math::Vector *up, int x,
                  float sin_phi, float cos_phi) {
  const bool far_arm = x >= W / 2;
  const int column = far_arm ? x - W / 2 : x;
  const float s = far_arm ? -sin_phi : sin_phi;
  const math::Vector &e = azimuth[column];
  const math::Vector &v = up[column];
  return math::Vector(s * e.x + cos_phi * v.x, s * e.y + cos_phi * v.y,
                      s * e.z + cos_phi * v.z);
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
  return rotated_direction<W>(azimuth, up, x, sin_phi, cos_phi);
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
   * @param state_moves Whether other per-column state of the frame moves, so
   *        its scans take their column-time variants without rotation.
   * @return Column-time state viewing this table, or the fixed meridian with
   *         temporal sampling disabled or without motion.
   * @details An orientation's history supplies both endpoints; one without
   * history, as after construction or a snap, has no motion.
   */
  HS_FLASH_MEMBER Temporal::ColumnTime prepare(const math::Quaternion &previous,
                                               const math::Quaternion &current,
                                               bool state_moves = false) {
    Temporal::ColumnTime time;
    if (!needed(previous, current, state_moves))
      return time;
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    const float *cos_theta = math::TrigLUT<W, H>::cos_theta_data();
    const float *sin_theta = math::TrigLUT<W, H>::sin_theta.data();
    const Temporal::TemporalRotation path(previous, current);
    path.for_each_sweep_sample(
        W, [&](int i, const math::Quaternion &at) HS_HOT_INLINE {
          const math::Quaternion delta = current * at.conjugate();
          up[i] = math::rotate(math::Vector(0.0f, 1.0f, 0.0f), delta);
          azimuth[i] = math::rotate(
              math::Vector(cos_theta[i], 0.0f, sin_theta[i]), delta);
        });
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

  /**
   * @brief Whether a frame calls for a basis table.
   * @param previous Orientation at the previous endpoint.
   * @param current Orientation at the current endpoint.
   * @param state_moves Whether other per-column state of the frame moves.
   * @return True with temporal sampling enabled and either motion.
   */
  static bool needed(const math::Quaternion &previous,
                     const math::Quaternion &current, bool state_moves) {
    return moving(previous, current) ||
           (state_moves && Temporal::temporal_sampling_enabled());
  }
};

/** @brief A compact shape's bounding cap where the columns that see it lit it.
 */
struct LitCap {
  math::Vector axis; ///< Cap axis turned to the middle seeing column's time.
  float arc; ///< Largest arc the cap moves over the seeing columns' times.
};

/**
 * @brief First-arm columns that can see a compact shape over one frame.
 * @tparam W Canvas width in columns.
 * @param axis Axis of the cap bounding the shape at the current endpoint.
 * @param radius Angular radius of that cap (radians).
 * @param sweep Largest arc the cap moves over the frame (radians).
 * @param first Receives the first seeing column, in [0, W / 2).
 * @param last Receives the last seeing column, in [@p first, W / 2).
 * @details Columns outside the range cannot see the shape on either arm. A
 * cap that holds a pole or crosses an arm boundary is seen by every column.
 */
template <int W>
inline void lit_columns(const math::Vector &axis, float radius, float sweep,
                        int &first, int &last) {
  constexpr float COLUMN = math::TWO_PI_F / W;
  constexpr int HALF = W / 2;
  // Two columns of reach cover the AA fringe and the column rounding; the
  // slack covers fast_atan2's and fast_acos's errors.
  const float reach = radius + sweep + 2.0f * COLUMN + SDF::MOVED_AXIS_SLACK;
  const float r_val = sqrtf(axis.x * axis.x + axis.z * axis.z);
  first = 0;
  last = HALF - 1;
  // reach < r_val <= 1 keeps the cap clear of both poles; asin(reach / r_val)
  // bounds the cap's azimuth half-width asin(sin(reach) / r_val).
  if (reach < r_val) {
    const float half =
        (0.5f * math::PI_F - math::fast_acos(reach / r_val)) / COLUMN;
    const float center = math::fast_atan2(axis.z, axis.x) / COLUMN;
    const int lo = static_cast<int>(floorf(center - half));
    const int hi = static_cast<int>(ceilf(center + half));
    const int arm_start =
        (lo >= 0 ? lo / HALF : -((HALF - 1 - lo) / HALF)) * HALF;
    if (hi < arm_start + HALF) {
      first = lo - arm_start;
      last = hi - arm_start;
    }
  }
}

/**
 * @brief Re-centres a compact shape's bounding cap at the sweep time of the
 *        columns that can see it.
 * @tparam W Canvas width in columns.
 * @tparam H Canvas height in rows.
 * @param axis Axis of the cap bounding the shape at the current endpoint.
 * @param radius Angular radius of that cap (radians).
 * @param time Column-time state with a rotated basis.
 * @param lit Receives the re-centred cap.
 * @return False when the re-centred cap would pass a hemisphere; the caller
 *         then sweeps the cap through the whole frame.
 * @details Column `x` sees the shape turned by the inverse of the rotation
 * that maps its fixed-meridian frame to its table entry; between two columns
 * that turn changes by at most `time.sweep` times their fraction difference.
 * The swept cap bounds the columns that can see the shape, which
 * lit_columns() finds. When they stay on one arm, clear of the poles, the cap
 * turned to their middle column's time and widened by that arc covers every
 * column's view; otherwise the cap turned to mid-frame and widened by half the
 * frame's turn does. Only the middle column's table entry is read.
 */
template <int W, int H>
inline bool lit_cap(const math::Vector &axis, float radius,
                    const Temporal::ColumnTime &time, LitCap &lit) {
  constexpr int HALF = W / 2;
  int first, last;
  lit_columns<W>(axis, radius, time.sweep, first, last);
  const int m = (first + last) / 2;
  const float arc =
      time.sweep * static_cast<float>(std::max(m - first, last - m)) / HALF +
      ColumnBasis<W, H>::SWEEP_PAD;
  if (!(radius + arc < 0.5f * math::PI_F))
    return false;
  if (!math::TrigLUT<W, H>::initialized)
    math::TrigLUT<W, H>::init();
  const math::Vector &e = time.azimuth[m];
  const math::Vector &v = time.up[m];
  const math::Vector t = math::cross(e, v);
  const float a = math::dot(axis, e);
  const float b = math::dot(axis, v);
  const float c = math::dot(axis, t);
  const float ct = math::TrigLUT<W, H>::cos_theta_data()[m];
  const float st = math::TrigLUT<W, H>::sin_theta[m];
  lit.axis = math::Vector(a * ct - c * st, b, a * st + c * ct);
  lit.arc = arc;
  return true;
}

/**
 * @brief Column-time state viewing a drifting shape's own basis table.
 * @tparam W Canvas width in columns.
 * @param arena Arena receiving the table; it lives until the caller's scope
 *        on the arena closes.
 * @param time Column-time state with a rotated basis and a drift.
 * @return State without a drift whose basis entries fill_drift_columns()
 *         supplies and whose sweep bounds the drifted basis's turn.
 */
template <int W>
HS_FLASH_MEMBER Temporal::ColumnTime
begin_drift_columns(Arena &arena, const Temporal::ColumnTime &time) {
  math::Vector *table = arena.allocate_n<math::Vector>(W);
  Temporal::ColumnTime view;
  view.azimuth = table;
  view.up = table + W / 2;
  view.sweep = time.sweep + sqrtf(math::dot(time.drift, time.drift));
  view.birth = time.birth;
  return view;
}

/**
 * @brief Fills a drifting shape's basis entries for a range of columns.
 * @tparam W Canvas width in columns.
 * @param view State begin_drift_columns() returned for @p time.
 * @param time Column-time state with a rotated basis and a drift.
 * @param first First first-arm column to fill.
 * @param last Last first-arm column to fill, inclusive.
 * @details Column `c`, lit at `u`, takes its basis entry rotated by
 * `drift * (1 - u)`, which rotates its pixel directions by the shape's motion
 * still to come.
 */
template <int W>
HS_FLASH_MEMBER void fill_drift_columns(const Temporal::ColumnTime &view,
                                        const Temporal::ColumnTime &time,
                                        int first, int last) {
  auto *azimuth = const_cast<math::Vector *>(view.azimuth);
  auto *up = const_cast<math::Vector *>(view.up);
  for (int c = first; c <= last; ++c) {
    const math::Vector rotation =
        time.drift * (1.0f - Temporal::sweep_fraction(c, W));
    azimuth[c] = Temporal::rotate_short(time.azimuth[c], rotation);
    up[c] = Temporal::rotate_short(time.up[c], rotation);
  }
}

/**
 * @brief lit_cap() for a drifting shape, filling the basis entries of every
 *        column its bounds reach.
 * @tparam W Canvas width in columns.
 * @tparam H Canvas height in rows.
 * @param axis Axis of the cap bounding the shape at the current endpoint.
 * @param radius Angular radius of that cap (radians).
 * @param view State begin_drift_columns() returned for @p time.
 * @param time Column-time state with a rotated basis and a drift.
 * @param lit Receives the re-centred cap.
 * @return lit_cap()'s result under @p view; on false every entry is filled.
 * @details The entries of the columns that can see the shape are filled
 * first, which holds the one lit_cap() reads, then those of the columns that
 * can see the re-centred cap.
 */
template <int W, int H>
HS_FLASH_MEMBER bool drift_lit_cap(const math::Vector &axis, float radius,
                                   const Temporal::ColumnTime &view,
                                   const Temporal::ColumnTime &time,
                                   LitCap &lit) {
  int first, last;
  lit_columns<W>(axis, radius, view.sweep, first, last);
  fill_drift_columns<W>(view, time, first, last);
  if (!lit_cap<W, H>(axis, radius, view, lit)) {
    fill_drift_columns<W>(view, time, 0, W / 2 - 1);
    return false;
  }
  int lit_first, lit_last;
  lit_columns<W>(lit.axis, radius + lit.arc, 0.0f, lit_first, lit_last);
  if (lit_first < first)
    fill_drift_columns<W>(view, time, lit_first, first - 1);
  if (lit_last > last)
    fill_drift_columns<W>(view, time, last + 1, lit_last);
  return true;
}

/**
 * @brief Object-space pixel frame of one column at its sweep time.
 * @details The column's pixel at colatitude phi looks along
 * `sin(phi) * radial + cos(phi) * up`; turning the azimuth by d replaces
 * `radial` with `cos(d) * radial + sin(d) * tangent`.
 */
struct ObjectColumn {
  math::Vector radial;  ///< Object-space direction of the column's azimuth.
  math::Vector tangent; ///< Object-space direction of increasing azimuth.
  math::Vector up;      ///< Object-space north.

  /**
   * @brief Object-space direction of the column's pixel at a colatitude.
   * @param sin_phi Sine of the colatitude.
   * @param cos_phi Cosine of the colatitude.
   * @return Unit direction.
   */
  HS_HOT_INLINE math::Vector direction(float sin_phi, float cos_phi) const {
    return math::Vector(sin_phi * radial.x + cos_phi * up.x,
                        sin_phi * radial.y + cos_phi * up.y,
                        sin_phi * radial.z + cos_phi * up.z);
  }
};

/**
 * @brief Per-column object-space pixel frames of an orientation over one
 *        effect frame.
 * @tparam W Canvas width in columns.
 * @tparam H Canvas height in rows.
 * @details An object placed in the world by `R(u)`, slerped between the
 * frame's endpoints, is seen by column `x`, lit at `u = sweep_fraction(x, W)`,
 * along its fixed-meridian frame rotated by `conj(R(u))`. A shader that
 * evaluates the object in its own space reads the pixel directions from this
 * table and needs no oriented copy of the object. Opposing columns share a
 * fraction, so one entry serves both arms. The table is frame-scoped.
 */
template <int W, int H> struct ObjectColumnBasis {
  static_assert(W % 4 == 0, "column basis reads the quarter-turn trig table");
  static constexpr int HALF = W / 2; ///< Distinct sweep fractions.

  std::array<math::Vector, HALF> radial;  ///< First-arm object azimuths.
  std::array<math::Vector, HALF> tangent; ///< First-arm azimuth tangents.
  std::array<math::Vector, HALF> up;      ///< Object north per fraction.

  /**
   * @brief Prepares the frame's object frames from its orientation endpoints.
   * @param previous Object-to-world orientation at the previous endpoint.
   * @param current Object-to-world orientation at the current endpoint.
   * @details Equal endpoints give every column the current frame.
   */
  HS_FLASH_MEMBER void prepare(const math::Quaternion &previous,
                               const math::Quaternion &current) {
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    const float *cos_theta = math::TrigLUT<W, H>::cos_theta_data();
    const float *sin_theta = math::TrigLUT<W, H>::sin_theta.data();
    const Temporal::TemporalRotation path(previous, current);
    path.for_each_sweep_sample(
        W, [&](int i, const math::Quaternion &at) HS_HOT_INLINE {
          const math::Quaternion inverse = at.conjugate();
          radial[i] = math::rotate(
              math::Vector(cos_theta[i], 0.0f, sin_theta[i]), inverse);
          tangent[i] = math::rotate(
              math::Vector(-sin_theta[i], 0.0f, cos_theta[i]), inverse);
          up[i] = math::rotate(math::Vector(0.0f, 1.0f, 0.0f), inverse);
        });
  }

  /**
   * @brief Object-space frame of a column.
   * @param x Column in [0, W).
   * @return The column's frame; the far arm's azimuth and tangent are the
   *         near arm's, negated.
   */
  HS_HOT_INLINE ObjectColumn column(int x) const {
    if (x < HALF)
      return {radial[x], tangent[x], up[x]};
    return {-radial[x - HALF], -tangent[x - HALF], up[x - HALF]};
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
 * @param state_moves Whether other per-column state of the frame moves; the
 *        table is then prepared without rotation too.
 * @return Column-time state for the frame's scans.
 */
template <int W, int H>
inline Temporal::ColumnTime
prepare_column_basis(Arena &arena, const math::Quaternion &previous,
                     const math::Quaternion &current,
                     bool state_moves = false) {
  if (!ColumnBasis<W, H>::needed(previous, current, state_moves))
    return {};
  return arena.make_default_initialized<ColumnBasis<W, H>>()->prepare(
      previous, current, state_moves);
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
 * @param state_moves Whether other per-column state of the frame moves.
 * @return Column-time state for the frame's scans.
 */
template <int W, int H, int CAP>
inline Temporal::ColumnTime
prepare_column_basis(Arena &arena, const math::Orientation<CAP> &orientation,
                     bool state_moves = false) {
  return prepare_column_basis<W, H>(arena, orientation.get(0),
                                    orientation.get(), state_moves);
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
