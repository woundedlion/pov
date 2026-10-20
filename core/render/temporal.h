/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file temporal.h
 * @brief Column-time sampling contract shared by every renderer family: the
 * fidelity toggle, the sweep fraction of a column, slab-wise simulation, the
 * shortest rotation arc between frame endpoints, the column-time state scans
 * read, shaders that read a column's sweep time, and capture-time solving for
 * plotted points.
 */

#include "core/animation/orientation.h"
#include "core/math/3dmath.h"
#include "core/math/pixel_mapping.h"
#include "core/platform/attributes.h"
#include "core/render/canvas.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include <span>
#include <type_traits>
#include <utility>

struct Pixel;

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

  /** @brief Visits the arc at each first-arm column's sweep fraction.
   * @tparam Fn Callable `(int column, const math::Quaternion &q)`.
   * @param width Even canvas width.
   * @param fn Receives each column `c` in [0, width / 2) with the arc at
   *        `sweep_fraction(c, width)`, as sample() gives it up to rounding.
   * @details Sines of evenly spaced angles obey
   * `sin((k + 1) d) = 2 cos(d) sin(k d) - sin((k - 1) d)`, so the sweep costs
   * three transcendental calls rather than two per column; the recurrence
   * runs in double, as its rounding grows with the column count. */
  template <typename Fn>
  HS_HOT_INLINE void for_each_sweep_sample(int width, Fn &&fn) const {
    const int half = width / 2;
    if (theta == 0.0f) {
      for (int c = 0; c < half; ++c)
        fn(c, (from + sweep_fraction(c, width) * (to - from)).normalized());
      return;
    }
    const double step = static_cast<double>(theta) / half;
    const double twice_cos = 2.0 * cos(step);
    const double inverse_sine = 1.0 / static_cast<double>(sine);
    // rising = sin(c * step) and falling = sin((half - c) * step), scaled.
    double rising = 0.0, rising_next = sin(step) * inverse_sine;
    double falling = 1.0;
    double falling_next = sin(static_cast<double>(theta) - step) * inverse_sine;
    for (int c = 0; c < half; ++c) {
      fn(c,
         (static_cast<float>(falling) * from + static_cast<float>(rising) * to)
             .normalized());
      const double up = twice_cos * rising_next - rising;
      rising = rising_next;
      rising_next = up;
      const double down = twice_cos * falling_next - falling;
      falling = falling_next;
      falling_next = down;
    }
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
  /// Rotation vector, axis times angle, carrying the scanned shape from where
  /// it stood at the previous endpoint to the current one. A Scan::rasterize
  /// column lit at fraction `u` evaluates the shape along its direction
  /// rotated by `drift * (1 - u)`; at most MAX_DRIFT_ANGLE long.
  math::Vector drift = math::Vector(0.0f, 0.0f, 0.0f);
  /// Per face of the mesh a deformation scan draws, which of its passes
  /// draws the face: MESH_FACE_HELD for the pass without a render-x window,
  /// MESH_FACE_WINDOW for the open window's; null outside such a scan.
  const uint8_t *selected_faces = nullptr;
  /// Per-row azimuth pads shared by the mesh's capture-window passes.
  const float *face_pads = nullptr;
};

/// ColumnTime::selected_faces entry of a face drawn once over the frame.
inline constexpr uint8_t MESH_FACE_HELD = 0;
/// ColumnTime::selected_faces entry of a moving face the open window draws.
inline constexpr uint8_t MESH_FACE_WINDOW = 1;
/// ColumnTime::selected_faces entry of a moving face the open window skips.
inline constexpr uint8_t MESH_FACE_SKIPPED = 2;

/// Longest drift a scan evaluates, radians.
inline constexpr float MAX_DRIFT_ANGLE = 0.25f;

/**
 * @brief Rotates a direction by a short rotation vector.
 * @param p Direction to rotate.
 * @param rotation Axis times angle, at most MAX_DRIFT_ANGLE long.
 * @return @p p rotated, by Rodrigues' formula with series coefficients.
 */
__attribute__((always_inline)) inline math::Vector
rotate_short(const math::Vector &p, const math::Vector &rotation) {
  const float angle_sq = math::dot(rotation, rotation);
  // Series for sin(a) / a and (1 - cos(a)) / a^2.
  const float sine = 1.0f - angle_sq * (1.0f / 6.0f - angle_sq / 120.0f);
  const float versine = 0.5f - angle_sq * (1.0f / 24.0f - angle_sq / 720.0f);
  const math::Vector across = math::cross(rotation, p);
  return p + across * sine + math::cross(rotation, across) * versine;
}

/** @brief First-arm column lit at a sweep fraction.
 * @param fraction A column's sweep fraction, as sweep_fraction() returns it.
 * @param width Even canvas width.
 * @return The column in [0, width / 2) lit at @p fraction. */
constexpr int sweep_column(float fraction, int width) {
  return static_cast<int>(fraction * static_cast<float>(width / 2) + 0.5f);
}

/**
 * @brief Whether a scan shader reads the sweep time of the pixel it shades.
 * @tparam ShaderT Shader type.
 * @tparam Args The scan's shader arguments.
 * @details A shader declares it by its signature: one invocable with a
 * trailing `float u` after the scan's arguments receives, in a column-time
 * scan, the sweep fraction of the shaded column, in [0, 1), and evaluates its
 * time-varying state there. The fixed-meridian scans call the plain form.
 */
template <typename ShaderT, typename... Args>
inline constexpr bool reads_column_time =
    std::is_invocable_v<ShaderT &, Args..., float>;

/**
 * @brief Runs a scan shader at a column's sweep time when it reads one.
 * @tparam W Canvas width in columns.
 * @tparam ShaderT Shader type.
 * @tparam Args The scan's shader arguments.
 * @param shader Shader to run.
 * @param x Shaded canvas column.
 * @param args The scan's shader arguments.
 */
template <int W, typename ShaderT, typename... Args>
__attribute__((always_inline)) inline void
shade_at_column(ShaderT &shader, int x, Args &&...args) {
  if constexpr (reads_column_time<ShaderT, Args...>)
    shader(std::forward<Args>(args)..., sweep_fraction(x, W));
  else
    shader(std::forward<Args>(args)...);
}

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

/**
 * @brief Advances a simulation and shades it in sweep slabs, each slab after
 *        the substeps that end its share of the frame.
 * @tparam SUBSTEPS Substeps per frame.
 * @tparam SLABS Slab count, a divisor of `SUBSTEPS`.
 * @param canvas Target canvas; each slab shades under a RenderXWindow.
 * @param advance Callable (int substeps) advancing the simulation.
 * @param shade Callable () shading the current state into the open window.
 * @details With temporal sampling disabled, runs every substep and then one
 * frame-wide shade.
 */
template <int SUBSTEPS, int SLABS, typename Advance, typename Shade>
HS_HOT_INLINE void simulate_in_slabs(Canvas &canvas, Advance &&advance,
                                     Shade &&shade) {
  static_assert(SLABS > 0 && SUBSTEPS % SLABS == 0);
  if (!temporal_sampling_enabled()) {
    advance(SUBSTEPS);
    shade();
    return;
  }
  for (int slab = 0; slab < SLABS; ++slab) {
    advance(SUBSTEPS / SLABS);
    const SweepSlab columns = sweep_slab(slab, SLABS, canvas.width());
    RenderXWindow window(canvas, columns.begin, columns.end);
    shade();
  }
}

/** @brief Most history nodes capture() accepts. */
inline constexpr int MAX_CAPTURE_NODES = 16;

/** @brief Frame fractions at which a display arm passes a moving point. */
struct Captures {
  static constexpr int CAPACITY = 2;      ///< Most arm passes in one frame.
  std::array<float, CAPACITY> fraction{}; ///< Frame fraction per pass.
  int count = 0;                          ///< Live passes.
};

/** @brief Tagged construction of a pipeline with explicit-time entry points. */
struct Timed {};

/** @brief Reduces a bounded displacement with nearest-even period ties.
 * @param value Displacement with an integer-representable period quotient.
 * @param period Positive period.
 * @return Remainder in [-period/2, period/2]. */
HS_HOT_INLINE inline float centered_remainder(float value, float period) {
  const int quotient = static_cast<int>(value / period);
  float result = value - static_cast<float>(quotient) * period;
  const float half = period * 0.5f;
  if (result > half || (result == half && (quotient & 1)))
    result -= period;
  else if (result < -half || (result == -half && (quotient & 1)))
    result += period;
  return result;
}

/** @brief Maximum dense tap-position error retained by a single-time splat. */
inline constexpr float MAX_SINGLE_TIME_TAP_ERROR = 0.003f;

/** @brief Pure boundary-column geometry and shading evaluator. */
struct BoundaryEvaluator {
  void *context;             ///< Borrowed evaluator state.
  bool sample_exists = true; ///< Material exists at the original splat time.
  bool force_timing = false; ///< Evaluate interior taps at their own times.
  /** Evaluates projected position, color, age and opacity; false before birth. */
  bool (*evaluate)(void *, float, math::PixelCoords &, Pixel &, float &,
                   float &);

  /** @brief Borrows a pure evaluator for one synchronous plot call.
   * @param fn Callable receiving time and output registers. */
  template <typename Fn>
  explicit BoundaryEvaluator(Fn &fn) : context(&fn), evaluate(invoke<Fn>) {}

private:
  /** @brief Calls the borrowed evaluator from flash.
   * @tparam Fn Pure callable type.
   * @param ctx Borrowed callable.
   * @param t Frame fraction.
   * @param p,c,age,alpha Output registers.
   * @return Whether the material is born. */
  template <typename Fn>
  HS_HOT_FLASH_MEMBER static bool invoke(void *ctx, float t,
                                         math::PixelCoords &p, Pixel &c,
                                         float &age, float &alpha) {
    return (*static_cast<Fn *>(ctx))(t, p, c, age, alpha);
  }
};

/** @brief Predicted unwrapped landing and pre-write acceptance window. */
struct ColumnPrediction {
  float column;    /**< Predicted unwrapped column. */
  float tolerance; /**< Largest accepted prediction error in columns. */
  float min_column =
      -std::numeric_limits<float>::infinity();               /**< Part start. */
  float max_column = std::numeric_limits<float>::infinity(); /**< Part end. */

  const BoundaryEvaluator *boundary = nullptr; ///< Optional tap-time evaluator.
  bool record_sample = true; ///< Whether this sample exists at its own time.

  /** @brief Unwraps a projected column near the prediction.
   * @param landed Wrapped projected column.
   * @param width Canvas width.
   * @return Nearest unwrapped column. */
  HS_HOT_INLINE float unwrap(float landed, float width) const {
    const float turns = (column - landed) / width;
    const int nearest = static_cast<int>(turns + (turns < 0.0f ? -0.5f : 0.5f));
    return landed + static_cast<float>(nearest) * width;
  }

  /** @brief Tests the time-prediction error bound.
   * @param landed Unwrapped projected column.
   * @return Whether the point needs its single correction. */
  HS_HOT_INLINE bool needs_correction(float landed) const {
    return fabsf(landed - column) > tolerance;
  }

  /** @brief Tests prediction and part membership before emission.
   * @param landed Unwrapped projected column.
   * @return Whether the point may write framebuffer taps. */
  HS_HOT_INLINE bool accepts(float landed) const {
    return !needs_correction(landed) && landed >= min_column &&
           landed <= max_column;
  }
};

/** @brief Work performed by geometry-callback capture and rasterization. */
struct GeometryCounters {
  uint32_t common_time_edges =
      0; /**< Edges using a common-time spacing metric. */
  uint32_t interior_timed_samples =
      0; /**< Samples requesting own-time tap evaluation. */
  uint32_t interior_single_time_samples =
      0; /**< Dense samples below the tap-motion limit. */
  uint32_t dense_cadence_edges =
      0; /**< Wider edges whose endpoint sample counts differ. */
  uint32_t dense_wide_edges =
      0; /**< Wider edges passing both curvature tests. */
  uint32_t dense_curved_edges = 0; /**< Wider edges rejected by curvature. */
  uint32_t dense_midpoint_probes =
      0;                            /**< Projected curvature-test midpoints. */
  uint32_t dense_spacing_edges = 0; /**< Edges beyond dense spacing. */
  uint32_t dense_polar_edges = 0; /**< Close edges rejected by polar bounds. */
  uint32_t dense_balanced_candidates =
      0; /**< Draws retaining the balanced edge sampler. */
  uint32_t dense_reject_scratch =
      0; /**< Candidates exceeding available scratch B. */
  uint32_t dense_sparse_edges =
      0; /**< Sparse edges inside partially dense polylines. */
  uint32_t dense_run_edges = 0;  /**< Edges emitted as dense runs. */
  uint32_t dense_candidates = 0; /**< Interpolating rigid density candidates. */
  uint32_t dense_candidate_controls = 0; /**< Controls examined for density. */
  uint32_t dense_reject_spacing =
      0; /**< Candidates rejected by control spacing. */
  uint32_t dense_reject_polar =
      0;                   /**< Candidates rejected by polar motion bounds. */
  uint32_t dense_runs = 0; /**< Polylines dispatched to dense sampling. */
  uint32_t dense_controls = 0;    /**< Dense original controls processed. */
  uint32_t dense_gap_samples = 0; /**< Interpolated dense material samples. */
  uint32_t dense_culls = 0; /**< Whole-polyline dense visibility decisions. */
  uint32_t dense_endpoint_probes = 0; /**< Dense endpoint projections. */

  uint32_t end_probes = 0;   /**< Probes at the current endpoint. */
  uint32_t solve_probes = 0; /**< Additional probes during capture solving. */
  uint32_t rigid_fallback_edges =
      0; /**< Edges exceeding adaptive capacity or precision. */
  uint32_t control_evaluations = 0; /**< Rigid world/control evaluations. */
  uint32_t control_corrections =
      0; /**< Corrected control capture predictions. */
  uint32_t rigid_samples =
      0; /**< Samples interpolated from transformed controls. */
  uint32_t boundary_evaluations = 0; /**< Kept boundary columns re-evaluated. */
  uint32_t corrections = 0; /**< Re-evaluations after a missed prediction. */
  uint32_t anchors = 0;     /**< Points solved for their arm passes. */
  uint32_t refined_controls = 0;     /**< Rigid controls added by refinement. */
  uint32_t callback_evaluations = 0; /**< Geometry callback invocations. */
  uint32_t parts = 0;                /**< Walks between arm-boundary cuts. */
  uint32_t held_anchors = 0; /**< Capture solves held at the column limit. */
  uint32_t unsettled_passes = 0; /**< Passes exhausting their probe budget. */
};

/** @brief Host instrumentation for the geometry-callback route. */
inline GeometryCounters geometry_counters;

/** @brief Increments a geometry counter in instrumented builds.
 * @tparam field Counter to increment.
 * @param amount Work items performed. */
template <uint32_t GeometryCounters::*field>
HS_HOT_INLINE inline void count_geometry(uint32_t amount = 1) {
#if HS_ENABLE_TEST_HOOKS
  geometry_counters.*field += amount;
#else
  (void)amount;
#endif
}

/** @brief Maximum solved motion in columns, one eighth of the canvas width. */
template <int W> inline constexpr float MAX_CAPTURE_COLUMNS = W / 8.0f;

/** @brief Motion bound supplied by geometry or a world stage. */
struct MotionBound {
  float radians = 0.0f; /**< Angular distance from the endpoint. */
  bool azimuth = false; /**< The bound is directly on azimuth. */

  /** @brief Converts the bound to columns at an endpoint row.
   * @tparam W,H Canvas dimensions.
   * @param row Projected endpoint row.
   * @return Maximum column displacement. */
  template <int W, int H> HS_HOT_INLINE float columns(float row) const {
    if (azimuth || radians == 0.0f)
      return radians * (W / (2.0f * math::PI_F));
    const float sine = fabsf(sinf(math::y_to_phi<H>(row)));
    if (radians >= math::PI_F / 2.0f || sine <= sinf(radians))
      return W / 4.0f;
    return asinf(sinf(radians) / sine) * (W / (2.0f * math::PI_F));
  }
};

/** @brief Continuous sweep fraction of an unwrapped column.
 * @tparam W Canvas width.
 * @param column Column coordinate.
 * @return Fraction in [0, 1). */
template <int W> HS_HOT_INLINE inline float column_fraction(float column) {
  const float half_turn = column / (W / 2.0f);
  const int whole = static_cast<int>(half_turn);
  const float fraction = half_turn - static_cast<float>(whole);
  return fraction < 0.0f ? fraction + 1.0f : fraction;
}

/** @brief Solves one arm branch, clamping extensions to frame endpoints.
 * @tparam W Canvas width.
 * @param probe Projected position at a frame fraction.
 * @param end Endpoint projection.
 * @param columns Maximum column displacement.
 * @param linear_pair Whether geometry interpolates two positions.
 * @param branch Unwrapped arm branch.
 * @return Capture fraction in [0, 1]. */
template <int W, typename Probe>
HS_HOT_INLINE float solve_capture_branch(Probe &&probe, math::PixelCoords end,
                                         float columns, bool linear_pair,
                                         int branch) {
  const float fraction = column_fraction<W>(end.x);
  const float half = W / 2.0f;
  const float boundary = fminf(fraction, 1.0f - fraction) * half;
  if (columns >= MAX_CAPTURE_COLUMNS<W> ||
      (columns < 3.8f && boundary > columns + 1.0f &&
       branch == static_cast<int>(floorf(end.x / half))))
    return fraction;
  count_geometry<&GeometryCounters::solve_probes>();
  const math::PixelCoords first = probe(0.0f);
  const float x0 =
      end.x + std::remainder(first.x - end.x, static_cast<float>(W));
  float lo = 0.0f, hi = 1.0f;
  float flo = x0 / half - branch, fhi = end.x / half - 1.0f - branch;
  if (flo <= 0.0f)
    return 0.0f;
  if (fhi >= 0.0f)
    return 1.0f;
  float t = flo / (flo - fhi);
  if (!(linear_pair && columns <= 4.0f)) {
    for (int iteration = 0; iteration < 4; ++iteration) {
      count_geometry<&GeometryCounters::solve_probes>();
      const float x =
          end.x + std::remainder(probe(t).x - end.x, static_cast<float>(W));
      const float residual = x / half - t - branch;
      if (fabsf(residual) * half <= 0.5f)
        break;
      if (iteration == 3) {
        count_geometry<&GeometryCounters::unsettled_passes>();
        break;
      }
      if (residual > 0.0f) {
        lo = t;
        flo = residual;
      } else {
        hi = t;
        fhi = residual;
      }
      t = lo + (hi - lo) * flo / (flo - fhi);
    }
  }
  return t;
}

/** @brief Solves each arm branch from a previously projected endpoint.
 * @tparam W Canvas width.
 * @param probe Projected position at a frame fraction.
 * @param end Endpoint projection at t = 1.
 * @param columns Maximum azimuth displacement in columns.
 * @param linear_pair Geometry uses a two-position interpolating callback.
 * @return Zero, one or two capture fractions. */
template <int W, typename Probe>
HS_HOT_INLINE Captures solve_capture(Probe &&probe, math::PixelCoords end,
                                     float columns, bool linear_pair) {
  count_geometry<&GeometryCounters::anchors>();
  Captures out;
  const float fraction = column_fraction<W>(end.x);
  const float boundary = fminf(fraction, 1.0f - fraction) * (W / 2.0f);
  if (columns >= MAX_CAPTURE_COLUMNS<W> ||
      (columns < 3.8f && boundary > columns + 1.0f)) {
    if (columns >= MAX_CAPTURE_COLUMNS<W>)
      count_geometry<&GeometryCounters::held_anchors>();
    out.count = 1;
    out.fraction[0] = fraction;
    return out;
  }
  count_geometry<&GeometryCounters::solve_probes>();
  const float start =
      end.x + std::remainder(probe(0.0f).x - end.x, static_cast<float>(W));
  const float r0 = start / (W / 2.0f);
  const float r1 = end.x / (W / 2.0f) - 1.0f;
  for (float n = floorf(r0); n > r1 && out.count < Captures::CAPACITY;
       n -= 1.0f) {
    float lo = 0.0f, hi = 1.0f, flo = r0 - n, fhi = r1 - n;
    float t = flo / (flo - fhi);
    if (!(linear_pair && columns <= 4.0f)) {
      for (int iteration = 0; iteration < 4; ++iteration) {
        count_geometry<&GeometryCounters::solve_probes>();
        const float x =
            end.x + std::remainder(probe(t).x - end.x, static_cast<float>(W));
        const float residual = x / (W / 2.0f) - t - n;
        if (fabsf(residual) * (W / 2.0f) <= 0.5f)
          break;
        if (iteration == 3) {
          count_geometry<&GeometryCounters::unsettled_passes>();
          break;
        }
        if (residual > 0.0f) {
          lo = t;
          flo = residual;
        } else {
          hi = t;
          fhi = residual;
        }
        t = lo + (hi - lo) * flo / (flo - fhi);
      }
    }
    out.fraction[out.count++] = t;
  }
  return out;
}

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
