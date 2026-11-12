/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <type_traits>
#include "core/render/sdf.h"
#include "core/render/shading.h"
#include "core/color/color.h"
#include "core/render/filter/pipeline.h"
#include "core/containers/static_circular_buffer.h"
#include "core/render/canvas.h"
#include "core/render/spacetime/worldlines.h"
#include "core/platform/platform.h"
#include "core/render/scan/raster.h"

/**
 * @file shapes.h
 * @brief The SDF-backed scan primitives: rings, polygons, circles, points, stars and flowers.
 */

namespace Scan {

/** @brief Optional DistortedRing draw settings; defaults draw a plain ring. */
struct RingDrawOptions {
  /** Angular phase offset in radians. */
  float phase = 0.0f;
  /** Render the bounding box for debugging. */
  bool debug_bb = false;
  /** Drop the degenerate exact-pole row instead of full-row filling it. */
  bool suppress_pole_fill = false;
};

/**
 * @brief Draws a ring whose radius is modulated around the circumference by
 *        shift_fn.
 */
struct DistortedRing {
  /**
   * @brief Rasterizes an undisplaced ring with exact polar centerline distance.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the ring plane.
   * @param radius Ring radius as a fraction of the hemisphere.
   * @param thickness Ring stroke half-width (radians).
   * @param fragment_shader Shader invoked per covered pixel.
   * @param options Phase, debug bounding box and pole-fill settings.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw_flat(ScanRef pipeline, Canvas &canvas,
                        const math::Basis &basis, float radius, float thickness,
                        FragmentShaderRef fragment_shader,
                        const RingDrawOptions &options = {}) {
    SDF::FlatDistortedRing shape(basis, radius, thickness, options.phase);
    shape.suppress_pole_fill = options.suppress_pole_fill;
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      options.debug_bb);
  }

  /**
   * @brief Rasterizes a circumference-modulated ring.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the ring plane.
   * @param radius Ring radius as a fraction of the hemisphere.
   * @param thickness Ring stroke half-width (radians).
   * @param shift_fn Scalar modulation function over the circumference.
   * @param amplitude Modulation amplitude (radians); must upper-bound
   *        |shift_fn|.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param options Phase, debug bounding box and pole-fill settings.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(ScanRef pipeline, Canvas &canvas, const math::Basis &basis,
                   float radius, float thickness, ScalarFn shift_fn,
                   float amplitude, FragmentShaderRef fragment_shader,
                   const RingDrawOptions &options = {}) {
    SDF::DistortedRing shape(basis, radius, thickness, shift_fn, amplitude,
                             options.phase);
    shape.suppress_pole_fill = options.suppress_pole_fill;
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      options.debug_bb);
  }

  /**
   * @brief Rasterizes a ring whose centerline is a shift-knot polyline with
   *        exact stroke distance (see SDF::DistortedRing's knot overload).
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the ring plane.
   * @param radius Ring radius as a fraction of the hemisphere.
   * @param thickness Ring stroke half-width (radians).
   * @param knots lut_n centerline shifts; closure wraps to entry 0.
   *        Storage must outlive the call.
   * @param lut_n Number of knot cells; at least 3.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param options Phase, debug bounding box and pole-fill settings.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(ScanRef pipeline, Canvas &canvas, const math::Basis &basis,
                   float radius, float thickness, const float *knots, int lut_n,
                   FragmentShaderRef fragment_shader,
                   const RingDrawOptions &options = {}) {
    SDF::KnotPrefilter prefilter;
    SDF::DistortedRing shape(basis, radius, thickness, knots, lut_n,
                             options.phase, prefilter);
    shape.suppress_pole_fill = options.suppress_pole_fill;
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      options.debug_bb);
  }
};

/**
 * @brief Draws a flat (tangent-plane) regular polygon projected onto the
 *        sphere.
 */
struct PlanarPolygon {
  /**
   * @brief Rasterizes a tangent-plane regular polygon.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the polygon plane.
   * @param radius Polygon circumradius as a fraction of the hemisphere, in (0, 2).
   * @param sides Number of polygon sides.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param phase Angular phase offset in radians. Added to the query azimuth,
   *        so the shape turns in the opposite sense to Plot's phase for the
   *        same basis and radius.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(ScanRef pipeline, Canvas &canvas, const math::Basis &basis,
                   float radius, int sides, FragmentShaderRef fragment_shader,
                   float phase = 0, bool debug_bb = false) {
    auto res = math::get_antipode(basis, radius);
    SDF::PlanarPolygon shape(res.first, res.second, sides, phase,
                             radius > 1.0f);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }

  /**
   * @brief Rasterizes a constant-color tangent-plane regular polygon.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the polygon plane.
   * @param radius Circumradius as a fraction of the hemisphere, in (0, 2).
   * @param sides Number of polygon sides.
   * @param color Fill color and alpha.
   * @param phase Angular phase offset in radians. Added to the query azimuth,
   *        so the shape turns in the opposite sense to Plot's phase for the
   *        same basis and radius.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, typename PipelineT>
  static void draw_solid(PipelineT &pipeline, Canvas &canvas,
                         const math::Basis &basis, float radius, int sides,
                         const Color4 &color, float phase = 0,
                         bool debug_bb = false) {
    auto res = math::get_antipode(basis, radius);
    SDF::PlanarPolygon shape(res.first, res.second, sides, phase,
                             radius > 1.0f);
    Scan::rasterize_solid<W, H>(pipeline, canvas, shape, color, debug_bb);
  }
};

/**
 * @brief Draws a great-circle line segment of given thickness between two
 *        vectors.
 */
struct Line {
  /**
   * @brief Rasterizes a great-circle segment between two unit vectors.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param v1 First endpoint as a world-space unit vector.
   * @param v2 Second endpoint as a world-space unit vector.
   * @param thickness Stroke half-width (radians).
   * @param fragment_shader Shader invoked per covered pixel.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H>
  static void draw(ScanRef pipeline, Canvas &canvas, const math::Vector &v1,
                   const math::Vector &v2, float thickness,
                   FragmentShaderRef fragment_shader, bool debug_bb = false) {
    SDF::Line shape(v1, v2, thickness);
    Scan::rasterize<W, H>(pipeline, canvas, shape, fragment_shader, debug_bb);
  }
};

/**
 * @brief Draws a ring stroke using SDF rasterization.
 */
struct Ring {
  /**
   * @brief Draws a ring stroke from an orientation basis.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the ring plane.
   * @param radius Ring radius as a fraction of the hemisphere, in [0, 2].
   * @param thickness Ring stroke half-width (radians).
   * @param fragment_shader Shader invoked per covered pixel.
   * @param phase Angular phase offset in radians.
   * @param debug_bb When true, renders the bounding box for debugging.
   * @details A radius past 1 keeps the caller's azimuth origin and handedness.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(ScanRef pipeline, Canvas &canvas, const math::Basis &basis,
                   float radius, float thickness,
                   FragmentShaderRef fragment_shader, float phase = 0,
                   bool debug_bb = false) {
    SDF::Ring shape(basis, radius, thickness, phase);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }

  /**
   * @brief Draws a ring stroke from a plane-normal vector.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param normal Plane normal as a world-space vector.
   * @param radius Ring radius as a fraction of the hemisphere.
   * @param thickness Ring stroke half-width (radians).
   * @param fragment_shader Shader invoked per covered pixel.
   * @param phase Angular phase offset in radians.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(ScanRef pipeline, Canvas &canvas, const math::Vector &normal,
                   float radius, float thickness,
                   FragmentShaderRef fragment_shader, float phase = 0,
                   bool debug_bb = false) {
    math::Basis basis = math::make_basis(math::Quaternion(), normal);
    draw<W, H, ComputeUVs>(pipeline, canvas, basis, radius, thickness,
                           fragment_shader, phase, debug_bb);
  }
};

/**
 * @brief Fused single-pass rasterizer for a small group of rings.
 */
HS_O3_BEGIN
struct RingGroup {
  /** Maximum number of rings in one group. */
  static constexpr int MAX_RINGS = 8;

  /** Largest distance, in columns, a slot's interpolated axis may stray
   *  from its great circle. */
  static constexpr float MAX_AXIS_ERROR_PX = 0.025f;

  /** Inclusive row spans where one column meets the covering band. */
  struct ColumnRows {
    int16_t lo[2]{}; /**< First rows. */
    int16_t hi[2]{}; /**< Last rows. */
    int count = 0;   /**< Number of disjoint spans. */
  };

  /** Angular solutions shared by the two arms at one sweep fraction. */
  struct ColumnBand {
    float center; /**< First-arm sinusoid phase. */
    float inner;  /**< Inner angular offset; greater than outer when empty. */
    float outer;  /**< Outer angular offset. */
  };

  /**
   * @brief Prepares a column's angular band intersections.
   * @param a First-arm coefficient of sin(phi).
   * @param b Coefficient of cos(phi).
   * @param low Lower cosine bound.
   * @param high Upper cosine bound.
   * @return Angular solutions about the first-arm phase.
   */
  HS_HOT_INLINE static ColumnBand column_band(float a, float b, float low,
                                              float high) {
    const float radius = sqrtf(a * a + b * b);
    if (low > radius || high < -radius)
      return {0.0f, 1.0f, 0.0f};
    if (radius < 1e-7f)
      return {0.0f, 0.0f, math::PI_F};
    const float inverse = 1.0f / radius;
    return {math::fast_atan2(a * inverse, b * inverse),
            math::fast_acos(hs::clamp(high * inverse, -1.0f, 1.0f)),
            math::fast_acos(hs::clamp(low * inverse, -1.0f, 1.0f))};
  }

  /**
   * @brief Fills the covering band's shared angular solutions once per draw.
   * @param bands Half-width table receiving angular solutions.
   * @param terms Per-column, per-slot axis terms.
   * @param half Number of distinct sweep fractions.
   * @param n Number of slots.
   * @param mid Covering slot index.
   * @param low Lower cosine bound.
   * @param high Upper cosine bound.
   */
  HS_HOT_FLASH_MEMBER static void fill_column_bands(ColumnBand *bands,
                                                    const float *terms,
                                                    int half, int n, int mid,
                                                    float low, float high) {
    for (int c = 0; c < half; ++c) {
      const float *axis = terms + 2 * (n * c + mid);
      bands[c] = column_band(axis[0], axis[1], low, high);
    }
  }

  /**
   * @brief Solves a column's sinusoidal axis dot against a cosine band.
   * @tparam W Canvas width.
   * @tparam H Canvas height.
   * @param a First-arm coefficient of sin(phi).
   * @param b Coefficient of cos(phi).
   * @param low Lower cosine bound.
   * @param high Upper cosine bound.
   * @param band Shared angular solutions for this sweep fraction.
   * @param far_arm Whether to reflect the first-arm phase.
   * @param y_lo First row the caller reads.
   * @param y_hi Last row the caller reads.
   * @return Up to two disjoint row spans with in-band endpoints. Within
   *         [@p y_lo, @p y_hi] their union is the band's whole-canvas spans';
   *         spans that cannot reach that window are dropped.
   * @pre The canvas trigonometric table is initialized.
   */
  template <int W, int H>
  HS_HOT_FLASH_MEMBER static ColumnRows
  column_rows(float a, float b, float low, float high, const ColumnBand &band,
              bool far_arm, int y_lo = 0, int y_hi = H - 1) {
    ColumnRows out;
    if (band.inner > band.outer)
      return out;
    const float center = far_arm ? -band.center : band.center;
    const float inner = band.inner, outer = band.outer;
    a = far_arm ? -a : a;
    constexpr float ANGLE_PAD = 4e-3f;
    const auto inside = [&](int y) HS_HOT_INLINE {
      const float d = a * math::TrigLUT<W, H>::sin_phi[y] +
                      b * math::TrigLUT<W, H>::cos_phi[y];
      return d >= low && d <= high;
    };
    SDF::Bounds candidates[6];
    int count = 0;
    for (int turn = -1; turn <= 1; ++turn) {
      const float origin = center + turn * math::TWO_PI_F;
      for (int side : {-1, 1}) {
        const float start = origin + (side < 0 ? -outer : inner) - ANGLE_PAD;
        const float end = origin + (side < 0 ? -inner : outer) + ANGLE_PAD;
        if (end < 0.0f || start > math::PI_F)
          continue;
        auto rows = SDF::phi_bounds_to_rows<H>(
            std::max(0.0f, start - 1e-5f), std::min(math::PI_F, end + 1e-5f));
        if (rows.y_max < y_lo || rows.y_min > y_hi)
          continue;
        while (rows.y_min <= rows.y_max && !inside(rows.y_min))
          ++rows.y_min;
        while (rows.y_min <= rows.y_max && !inside(rows.y_max))
          --rows.y_max;
        if (rows.y_min > rows.y_max)
          continue;
        int at = count++;
        while (at > 0 && candidates[at - 1].y_min > rows.y_min) {
          candidates[at] = candidates[at - 1];
          --at;
        }
        candidates[at] = rows;
      }
    }
    for (int i = 0; i < count; ++i) {
      const auto rows = candidates[i];
      if (out.count && rows.y_min <= out.hi[out.count - 1] + 1) {
        out.hi[out.count - 1] =
            std::max<int>(out.hi[out.count - 1], rows.y_max);
      } else {
        out.lo[out.count] = rows.y_min;
        out.hi[out.count++] = rows.y_max;
      }
    }
    return out;
  }

  /**
   * @brief Rasterizes every ring of a group in one scan over the union band.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @tparam RingShaderT Per-ring shader: shader(int slot, const Vector &p,
   *         Fragment &f). One Fragment serves the whole scan and only color,
   *         pos, v2 (stroke coverage), size and age are refreshed per pixel —
   *         no UVs, no raw distance. v0, v1 and v3 are NaN in debug builds and
   *         retain defaults or previous values in release builds. The shader
   *         must not read inputs it did not set.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param shapes Ring shapes in draw order.
   * @param n Number of shapes, in [1, MAX_RINGS].
   * @param shader Per-ring fragment shader (see RingShaderT).
   * @param debug_bb When true, or under canvas.debug(), falls back to per-ring
   *        rasterizes, which scan each ring's own row intervals and fill
   *        v0/v1/v3 per pixel.
   * @param previous_normals Optional per-slot axes at the previous frame
   *        endpoint. With temporal sampling on, column `x` evaluates each
   *        slot about its axis at `surface_time(x, W)` along the great
   *        circle from its previous axis to its current one, and each slot's
   *        bounds grow by the angle between the two.
   * @details Row intervals come from one covering ring: the middle member (n/2)
   * inflated by the group's maximum plane/radius deviation plus thickness.
   * Members evaluate per pixel in ascending slot order. The canvas's
   * column-time basis is honoured. This draw ignores the birth fraction and
   * the render-x window of an open TimeWindowScope. At
   * pole_lod_aggressiveness 0 and one destination pixel per plot, the only
   * divergence from rasterizing the rings one by one is AA-tail pixels a
   * member's own interval clip would drop.
   */
  template <int W, int H, typename PipelineT, typename RingShaderT>
  static void draw(PipelineT &pipeline, Canvas &canvas, const SDF::Ring *shapes,
                   int n, RingShaderT &&shader, bool debug_bb = false,
                   const math::Vector *previous_normals = nullptr) {
    check_canvas_dims<W, H>(canvas);
    check_pipeline_prepared(pipeline, canvas);
    HS_CHECK(n >= 1 && n <= MAX_RINGS,
             "ring group size must be in [1, MAX_RINGS]");
    if (debug_bb || canvas.debug()) {
      for (int s = 0; s < n; ++s) {
        auto slot_shader = [&](const math::Vector &p, Fragment &f) {
          shader(s, p, f);
        };
        if constexpr (std::is_constructible_v<ScanRef, PipelineT &>) {
          ScanRef ref(pipeline);
          Scan::rasterize<W, H>(ref, canvas, shapes[s], slot_shader, true);
        } else {
          Scan::rasterize<W, H, true>(pipeline, canvas, shapes[s], slot_shader,
                                      true);
        }
      }
      return;
    }

    const Scan::ColumnBasis *time = column_basis(pipeline, canvas);
    const bool rotated = time && time->azimuth;
    const bool moving_slots = previous_normals && !Spacetime::flat_surface();
    if (rotated || moving_slots) {
      draw_column_time<W, H>(pipeline, canvas, shapes, n, shader,
                             rotated ? time : nullptr,
                             moving_slots ? previous_normals : nullptr);
      return;
    }
    draw_fixed<W, H>(pipeline, canvas, shapes, n, shader);
  }

  /**
   * @brief RingGroup::draw at the fixed meridian with still slots, placed in
   *        FLASH.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @tparam RingShaderT Per-ring shader.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param shapes Ring shapes in draw order.
   * @param n Number of shapes.
   * @param shader Per-ring fragment shader.
   */
  template <int W, int H, typename PipelineT, typename RingShaderT>
  HS_HOT_FLASH_MEMBER static void
  draw_fixed(PipelineT &pipeline, Canvas &canvas, const SDF::Ring *shapes,
             int n, RingShaderT &shader) {
    draw_body<W, H, false, true>(pipeline, canvas, shapes, n, shader, nullptr,
                                 nullptr);
  }

  /**
   * @brief RingGroup::draw under a rotated basis or moving slot axes.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @tparam RingShaderT Per-ring shader.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param shapes Ring shapes in draw order.
   * @param n Number of shapes.
   * @param shader Per-ring fragment shader.
   * @param time Column-time state with a rotated basis, or null.
   * @param previous_normals Per-slot axes at the previous endpoint, or null.
   */
  template <int W, int H, typename PipelineT, typename RingShaderT>
  HS_NOINLINE_NOCLONE static void
  draw_column_time(PipelineT &pipeline, Canvas &canvas, const SDF::Ring *shapes,
                   int n, RingShaderT &shader, const Scan::ColumnBasis *time,
                   const math::Vector *previous_normals) {
    draw_body<W, H, true, false>(pipeline, canvas, shapes, n, shader, time,
                                 previous_normals);
  }

  /**
   * @brief Fills each column's ring-axis terms for the column-time scan.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @param terms 2 * n * W / 2 entries receiving, per first-arm column and
   *        slot, the slot's axis at the column's sweep time dotted with the
   *        column's azimuth and with its up direction.
   * @param shapes Ring shapes in draw order.
   * @param n Number of shapes.
   * @param azimuth Rotated first-arm azimuths, or null for the fixed meridian.
   * @param up Rotated up directions; read only with @p azimuth.
   * @param previous_normals Per-slot axes at the previous endpoint, or null.
   * @details A pixel at colatitude phi in column x has axis dot
   * `+-sin(phi) * a + cos(phi) * b`, negated on the far arm. A slot takes its
   * normalized chord, unless turn_fast_axes() moves it along its great
   * circle.
   */
  template <int W, int H>
  HS_NOINLINE_NOCLONE static void
  fill_column_axes(float *terms, const SDF::Ring *shapes, int n,
                   const math::Vector *azimuth, const math::Vector *up,
                   const math::Vector *previous_normals) {
    const float *cos_theta = math::TrigLUT<W, H>::cos_theta_data();
    const float *sin_theta = math::TrigLUT<W, H>::sin_theta.data();
    // Slot s at fraction u has axis (start + u * motion) / |start + u * motion|.
    math::Vector start[MAX_RINGS], motion[MAX_RINGS];
    float length_terms[MAX_RINGS][3];
    for (int s = 0; s < n; ++s) {
      start[s] = previous_normals ? previous_normals[s] : shapes[s].normal;
      motion[s] = shapes[s].normal - start[s];
      length_terms[s][0] = math::dot(start[s], start[s]);
      length_terms[s][1] = 2.0f * math::dot(start[s], motion[s]);
      length_terms[s][2] = math::dot(motion[s], motion[s]);
    }
    for (int c = 0; c < W / 2; ++c) {
      const math::Vector e =
          azimuth ? azimuth[c] : math::Vector(cos_theta[c], 0.0f, sin_theta[c]);
      const math::Vector v = azimuth ? up[c] : math::Vector(0.0f, 1.0f, 0.0f);
      const float u = Spacetime::SamplingSurface<W>::time(c);
      float *out = terms + 2 * n * c;
      for (int s = 0; s < n; ++s) {
        const float length_sq =
            length_terms[s][0] +
            u * (length_terms[s][1] + u * length_terms[s][2]);
        float inverse = math::fast_rsqrt(length_sq);
        // a third Newton step: the ring's cosine band needs full precision
        inverse *= 1.5f - 0.5f * length_sq * inverse * inverse;
        out[2 * s] =
            (math::dot(e, start[s]) + u * math::dot(e, motion[s])) * inverse;
        out[2 * s + 1] =
            (math::dot(v, start[s]) + u * math::dot(v, motion[s])) * inverse;
      }
    }
    turn_fast_axes<W, H>(terms, shapes, n, azimuth, up, previous_normals);
  }

  /**
   * @brief Rewrites the axis terms of the slots too fast for their chord with
   *        their great-circle turn.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @param terms Axis terms fill_column_axes() wrote.
   * @param shapes Ring shapes in draw order.
   * @param n Number of shapes.
   * @param azimuth Rotated first-arm azimuths, or null for the fixed meridian.
   * @param up Rotated up directions; read only with @p azimuth.
   * @param previous_normals Per-slot axes at the previous endpoint, or null
   *        for slots that hold still.
   * @details Each slot's interpolant comes from its turn's analytic bound: the
   * normalized chord while it stays within MAX_AXIS_ERROR_PX of the
   * constant-rate turn, else the arc, turned per column by a rotation
   * recurrence.
   */
  template <int W, int H>
  HS_FLASH_MEMBER static void
  turn_fast_axes(float *terms, const SDF::Ring *shapes, int n,
                 const math::Vector *azimuth, const math::Vector *up,
                 const math::Vector *previous_normals) {
    const float *cos_theta = math::TrigLUT<W, H>::cos_theta_data();
    const float *sin_theta = math::TrigLUT<W, H>::sin_theta.data();
    constexpr float TOLERANCE = MAX_AXIS_ERROR_PX * math::RADIANS_PER_COLUMN<W>;
    for (int s = 0; previous_normals && s < n; ++s) {
      const math::Vector &start = previous_normals[s];
      const math::Vector &end = shapes[s].normal;
      // The turn is at most pi / 2 times its chord.
      const math::Vector chord = end - start;
      if (Spacetime::choose_turn_interpolant(
              0.5f * math::PI_F * sqrtf(math::dot(chord, chord)), TOLERANCE) ==
          Spacetime::Interpolant::CHORD)
        continue;
      const math::Vector normal = math::cross(start, end);
      const float sine = sqrtf(math::dot(normal, normal));
      const float angle = atan2f(sine, math::dot(start, end));
      if (!(sine > 0.0f) ||
          Spacetime::choose_turn_interpolant(angle, TOLERANCE) !=
              Spacetime::Interpolant::ARC)
        continue;
      const math::Vector tangent = math::cross(normal, start) * (1.0f / sine);
      // Double precision keeps the recurrence's drift over the half-turn of
      // columns far below the stroke's coverage slope.
      const double step = static_cast<double>(angle) / (W / 2);
      const double step_cos = cos(step), step_sin = sin(step);
      double turn_cos = 1.0, turn_sin = 0.0;
      for (int c = 0; c < W / 2; ++c) {
        const math::Vector e =
            azimuth ? azimuth[c]
                    : math::Vector(cos_theta[c], 0.0f, sin_theta[c]);
        const math::Vector v = azimuth ? up[c] : math::Vector(0.0f, 1.0f, 0.0f);
        const math::Vector axis = start * static_cast<float>(turn_cos) +
                                  tangent * static_cast<float>(turn_sin);
        float *out = terms + 2 * (n * c + s);
        out[0] = math::dot(e, axis);
        out[1] = math::dot(v, axis);
        const double next = turn_cos * step_cos - turn_sin * step_sin;
        turn_sin = turn_sin * step_cos + turn_cos * step_sin;
        turn_cos = next;
      }
    }
  }

  /**
   * @brief Fused scan body shared by the fixed and column-time draws.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam COLUMN_TIME Whether @p time and @p previous_normals are read.
   * @tparam INLINE_WALK Whether the row walk inlines into the body, so a
   *         FLASH variant takes it along.
   * @tparam PipelineT Plotting pipeline type.
   * @tparam RingShaderT Per-ring shader.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param shapes Ring shapes in draw order.
   * @param n Number of shapes.
   * @param shader Per-ring fragment shader.
   * @param time Column-time state with a rotated basis, or null.
   * @param previous_normals Per-slot axes at the previous endpoint, or null.
   */
  template <int W, int H, bool COLUMN_TIME, bool INLINE_WALK,
            typename PipelineT, typename RingShaderT>
  HS_HOT_INLINE static void
  draw_body(PipelineT &pipeline, Canvas &canvas, const SDF::Ring *shapes, int n,
            RingShaderT &shader, const Scan::ColumnBasis *time,
            const math::Vector *previous_normals) {
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();

    [[maybe_unused]] const math::Vector *azimuth = nullptr;
    [[maybe_unused]] const math::Vector *up = nullptr;
    [[maybe_unused]] float sweep = 0.0f;
    [[maybe_unused]] float slot_sweep[MAX_RINGS] = {};
    [[maybe_unused]] float max_slot_sweep = 0.0f;
    ScratchScope scratch(scratch_arena_b);
    [[maybe_unused]] float *axis_terms = nullptr;
    [[maybe_unused]] ColumnBand *column_bands = nullptr;
    if constexpr (COLUMN_TIME) {
      if (time) {
        azimuth = time->azimuth;
        up = time->up;
        sweep = time->sweep_arc;
      }
      if (previous_normals) {
        for (int s = 0; s < n; ++s) {
          slot_sweep[s] = math::fast_acos(hs::clamp(
                              math::dot(shapes[s].normal, previous_normals[s]),
                              -1.0f, 1.0f)) +
                          1e-3f;
          max_slot_sweep = std::max(max_slot_sweep, slot_sweep[s]);
        }
      }
    }
    int sy_min[MAX_RINGS], sy_max[MAX_RINGS];
    int y_lo = H, y_hi = -1;
    for (int s = 0; s < n; ++s) {
      auto b = shapes[s].template get_vertical_bounds<H>();
      if constexpr (COLUMN_TIME)
        b = swept_rows<H>(b, sweep + slot_sweep[s]);
      sy_min[s] = b.y_min;
      sy_max[s] = b.y_max;
      y_lo = std::min(y_lo, b.y_min);
      y_hi = std::max(y_hi, b.y_max);
    }
    const auto &cr = source_clip<W, H>(pipeline, canvas);
    y_lo = std::max(y_lo, cr.render_y_start());
    y_hi = std::min(y_hi, cr.render_y_end() - 1);
    if (y_lo > y_hi)
      return;

    // Covering ring: the middle member inflated by the worst centerline
    // deviation plus thickness. The 1e-3 absorbs fast_acos error.
    const int mid = n / 2;
    auto member_reach = [&](const math::Vector &mid_axis,
                            const math::Vector &member_axis,
                            int s) HS_HOT_INLINE {
      const float dev =
          math::fast_acos(
              hs::clamp(math::dot(mid_axis, member_axis), -1.0f, 1.0f)) +
          std::abs(shapes[s].target_angle - shapes[mid].target_angle) + 1e-3f;
      return shapes[s].thickness + dev;
    };
    float pad_th = shapes[mid].thickness;
    for (int s = 0; s < n; ++s) {
      if (s == mid)
        continue;
      pad_th = std::max(pad_th,
                        member_reach(shapes[mid].normal, shapes[s].normal, s));
    }
    SDF::Ring cover(
        math::Basis{shapes[mid].u, shapes[mid].normal, shapes[mid].w},
        shapes[mid].radius, pad_th);
    // A column's lit pixels lie in the middle member's band, about its axis at
    // the column's time, widened by every member's deviation from it at
    // either endpoint; the lerped axes' deviations stay within the larger,
    // up to their normalisation.
    [[maybe_unused]] float cover_cos_min = -2.0f, cover_cos_max = 2.0f;
    if constexpr (COLUMN_TIME) {
      float pad = pad_th;
      if (previous_normals) {
        for (int s = 0; s < n; ++s) {
          if (s == mid)
            continue;
          pad = std::max(
              pad, member_reach(previous_normals[mid], previous_normals[s], s));
        }
        pad = pad * (1.0f + max_slot_sweep * max_slot_sweep) + 1e-3f;
      }
      SDF::band_edge_cosines(shapes[mid].target_angle - pad,
                             shapes[mid].target_angle + pad, cover_cos_max,
                             cover_cos_min);
    }
    if constexpr (COLUMN_TIME) {
      axis_terms = scratch_arena_b.allocate_n<float>(2 * n * (W / 2));
      fill_column_axes<W, H>(axis_terms, shapes, n, azimuth, up,
                             previous_normals);
      column_bands = scratch_arena_b.allocate_n<ColumnBand>(W / 2);
      fill_column_bands(column_bands, axis_terms, W / 2, n, mid, cover_cos_min,
                        cover_cos_max);
    }

    Fragment frag;

    // Walk the covering ring's row runs (shared emit_row_runs).
    const float *cos_theta = math::TrigLUT<W, H>::cos_theta_data();
    const float *sin_theta = math::TrigLUT<W, H>::sin_theta.data();
    const auto xc = cr.x_clip();
    StaticCircularBuffer<SDF::Interval, 4> intervals;
    StaticCircularBuffer<SDF::Interval, 8> norm;
    static_assert(decltype(intervals)::CAPACITY >=
                      SDF::sdf_max_spans<SDF::Ring>::value,
                  "intervals must hold the covering Ring's per-row emission");
    static_assert(decltype(norm)::CAPACITY == 2 * decltype(intervals)::CAPACITY,
                  "norm must hold 2 spans per input interval (seam split)");
    int active[MAX_RINGS];
    int n_active = 0;

    auto pixel_run_body = [&](int x1, int x2, int y, float sp,
                              float cp) HS_HOT_INLINE {
      for (int x = x1; x < x2; ++x) {
        [[maybe_unused]] float s_phi = sp;
        [[maybe_unused]] const float *terms = nullptr;
        math::Vector p;
        if constexpr (COLUMN_TIME) {
          const bool far_arm = x >= W / 2;
          s_phi = far_arm ? -sp : sp;
          terms = axis_terms + 2 * n * (far_arm ? x - W / 2 : x);
          const float d_mid = s_phi * terms[2 * mid] + cp * terms[2 * mid + 1];
          if (d_mid < cover_cos_min || d_mid > cover_cos_max)
            continue;
        } else {
          p = math::Vector(sp * cos_theta[x], cp, sp * sin_theta[x]);
        }
        for (int i = 0; i < (COLUMN_TIME ? n : n_active); ++i) {
          const int s = COLUMN_TIME ? i : active[i];
          if constexpr (COLUMN_TIME)
            if (y < sy_min[s] || y > sy_max[s])
              continue;
          float d;
          if constexpr (COLUMN_TIME)
            d = s_phi * terms[2 * s] + cp * terms[2 * s + 1];
          else
            d = math::dot(p, shapes[s].normal);
          const float alpha = shapes[s].stroke_alpha(d);
          if (alpha <= MIN_ALPHA)
            continue;
          if constexpr (COLUMN_TIME)
            p = column_direction<W>(azimuth, up, x, sp, cp, cos_theta,
                                    sin_theta);
#ifndef NDEBUG
          frag.poison_inputs();
#endif
          frag.color = Color4(0, 0, 0, 0);
          frag.pos = p;
          frag.v2 = alpha;
          frag.size = shapes[s].thickness;
          frag.age = 0;
          shader(s, p, frag);
          if (frag.color.alpha > MIN_ALPHA)
            plot_admitted(pipeline, canvas, x, y, frag.color.color, frag.age,
                          frag.color.alpha * alpha);
        }
      }
    };
    // One out-of-line copy of the run body serves every emit site.
    [[maybe_unused]] auto pixel_run =
        [&](int x1, int x2, int y, float sp, float cp)
            HS_NOINLINE_NOCLONE { pixel_run_body(x1, x2, y, sp, cp); };
    if constexpr (COLUMN_TIME) {
      walk_clip_columns<W>(xc, [&](int x) HS_HOT_INLINE {
        const bool far_arm = x >= W / 2;
        const int c = far_arm ? x - W / 2 : x;
        const float *terms = axis_terms + 2 * n * c;
        const auto rows = column_rows<W, H>(
            terms[2 * mid], terms[2 * mid + 1], cover_cos_min, cover_cos_max,
            column_bands[c], far_arm, y_lo, y_hi);
        for (int span = 0; span < rows.count; ++span) {
          const int end = std::min<int>(y_hi, rows.hi[span]);
          for (int y = std::max<int>(y_lo, rows.lo[span]); y <= end; ++y)
            pixel_run_body(x, x + 1, y, math::TrigLUT<W, H>::sin_phi[y],
                           math::TrigLUT<W, H>::cos_phi[y]);
        }
      });
      return;
    }
    for (int y = y_lo; y <= y_hi; ++y) {
      // Per-slot domains match solo rasterization; slot order preserves blending.
      n_active = 0;
      for (int s = 0; s < n; ++s)
        if (y >= sy_min[s] && y <= sy_max[s])
          active[n_active++] = s;
      if (n_active == 0)
        continue;
      float sp = math::TrigLUT<W, H>::sin_phi[y];
      float cp = math::TrigLUT<W, H>::cos_phi[y];
      auto walk_row = [&](auto &emit) HS_HOT_INLINE {
        if (cover.needs_full_row_scan(sp)) {
          clip_run(0, W, xc, emit);
          return;
        }
        intervals.clear();
        cover.template get_horizontal_intervals<W, H>(
            y,
            [&](float t1, float t2) { SDF::push_interval(intervals, t1, t2); });
        // The full-row case is taken above, so the covering ring always
        // answers.
        emit_row_runs<W>(true, intervals, norm, xc, emit);
      };
      if constexpr (INLINE_WALK) {
        auto emit = [&](int x1, int x2)
                        HS_HOT_INLINE { pixel_run_body(x1, x2, y, sp, cp); };
        walk_row(emit);
      } else {
        auto emit = [&](int x1, int x2) { pixel_run(x1, x2, y, sp, cp); };
        walk_row(emit);
      }
    }
  }
};
HS_O3_END

/**
 * @brief Draws a disc as a zero-radius ring whose stroke spans it.
 * @details Stroke coverage falls off quintically from 1 at the center to 0 at
 * the rim; it modulates the plotted alpha and reaches the shader as the
 * fragment's register 2, so the shader picks the radial style on top of it.
 */
struct Circle {
  /**
   * @brief Draws a disc from an orientation basis.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the circle plane.
   * @param radius Circle radius as a fraction of the hemisphere.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(ScanRef pipeline, Canvas &canvas, const math::Basis &basis,
                   float radius, FragmentShaderRef fragment_shader,
                   bool debug_bb = false) {
    float th = radius * (math::PI_F / 2.0f);
    Ring::draw<W, H, ComputeUVs>(pipeline, canvas, basis, 0.0f, th,
                                 fragment_shader, 0, debug_bb);
  }

  /**
   * @brief Draws a disc from a plane-normal vector.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param normal Plane normal as a world-space vector.
   * @param radius Circle radius as a fraction of the hemisphere.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(ScanRef pipeline, Canvas &canvas, const math::Vector &normal,
                   float radius, FragmentShaderRef fragment_shader,
                   bool debug_bb = false) {
    math::Basis basis = math::make_basis(math::Quaternion(), normal);
    draw<W, H, ComputeUVs>(pipeline, canvas, basis, radius, fragment_shader,
                           debug_bb);
  }
};

/**
 * @brief Draws a dot as a zero-radius ring of the given thickness.
 * @details Stroke coverage falls off quintically from 1 at the center to 0 at
 * the rim, giving the soft glow effects shade through the fragment's
 * register 2.
 */
struct Point {
  /**
   * @brief Draws a dot centered on a unit vector.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param p Point center as a world-space unit vector.
   * @param thickness Point radius as a stroke half-width (radians).
   * @param fragment_shader Shader invoked per covered pixel.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H>
  static void draw(ScanRef pipeline, Canvas &canvas, const math::Vector &p,
                   float thickness, FragmentShaderRef fragment_shader,
                   bool debug_bb = false) {
    math::Basis basis = math::make_basis(math::Quaternion(), p);
    Ring::draw<W, H>(pipeline, canvas, basis, 0.0f, thickness, fragment_shader,
                     0.0f, debug_bb);
  }

  /**
   * @brief Draws a dot whose shader reads its column's sweep time.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param p Point center as a world-space unit vector.
   * @param thickness Point radius as a stroke half-width (radians).
   * @param fragment_shader Shader invoked per covered pixel with its column's
   *        sweep fraction.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H>
  static void draw(ScanRef pipeline, Canvas &canvas, const math::Vector &p,
                   float thickness, TimedShaderRef fragment_shader,
                   bool debug_bb = false) {
    const SDF::Ring shape(math::make_basis(math::Quaternion(), p), 0.0f,
                          thickness, 0.0f);
    Scan::rasterize<W, H>(pipeline, canvas, shape, fragment_shader, debug_bb);
  }
};

/**
 * @brief Draws a solid star shape.
 */
struct Star {
  /**
   * @brief Rasterizes a solid star.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the star plane.
   * @param radius Star circumradius as a fraction of the hemisphere, in (0, 2).
   * @param sides Number of star points.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param phase Angular phase offset in radians. Added to the query azimuth,
   *        so the shape turns in the opposite sense to Plot's phase for the
   *        same basis and radius.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(ScanRef pipeline, Canvas &canvas, const math::Basis &basis,
                   float radius, int sides, FragmentShaderRef fragment_shader,
                   float phase = 0, bool debug_bb = false) {
    auto res = math::get_antipode(basis, radius);
    SDF::Star shape(res.first, res.second, sides, phase, radius > 1.0f);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }

  /**
   * @brief Rasterizes a constant-color solid star.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the star plane.
   * @param radius Circumradius as a fraction of the hemisphere, in (0, 2).
   * @param sides Number of star points.
   * @param color Fill color and alpha.
   * @param phase Angular phase offset in radians. Added to the query azimuth,
   *        so the shape turns in the opposite sense to Plot's phase for the
   *        same basis and radius.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, typename PipelineT>
  static void draw_solid(PipelineT &pipeline, Canvas &canvas,
                         const math::Basis &basis, float radius, int sides,
                         const Color4 &color, float phase = 0,
                         bool debug_bb = false) {
    auto res = math::get_antipode(basis, radius);
    SDF::Star shape(res.first, res.second, sides, phase, radius > 1.0f);
    Scan::rasterize_solid<W, H>(pipeline, canvas, shape, color, debug_bb);
  }
};

/**
 * @brief Draws a solid flower shape.
 */
struct Flower {
  /**
   * @brief Rasterizes a solid flower.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the flower plane.
   * @param radius Flower circumradius as a fraction of the hemisphere, in (0, 2).
   * @param sides Number of flower petals.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param phase Angular phase offset in radians. Added to the query azimuth,
   *        so the shape turns in the opposite sense to Plot's phase for the
   *        same basis and radius.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(ScanRef pipeline, Canvas &canvas, const math::Basis &basis,
                   float radius, int sides, FragmentShaderRef fragment_shader,
                   float phase = 0, bool debug_bb = false) {
    auto res = math::get_antipode(basis, radius);
    SDF::Flower shape(res.first, res.second, sides, phase, radius > 1.0f);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }

  /**
   * @brief Rasterizes a constant-color solid flower.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the flower plane.
   * @param radius Circumradius as a fraction of the hemisphere, in (0, 2).
   * @param sides Number of flower petals.
   * @param color Fill color and alpha.
   * @param phase Angular phase offset in radians. Added to the query azimuth,
   *        so the shape turns in the opposite sense to Plot's phase for the
   *        same basis and radius.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, typename PipelineT>
  static void draw_solid(PipelineT &pipeline, Canvas &canvas,
                         const math::Basis &basis, float radius, int sides,
                         const Color4 &color, float phase = 0,
                         bool debug_bb = false) {
    auto res = math::get_antipode(basis, radius);
    SDF::Flower shape(res.first, res.second, sides, phase, radius > 1.0f);
    Scan::rasterize_solid<W, H>(pipeline, canvas, shape, color, debug_bb);
  }
};

/**
 * @brief Draws a solid spherical polygon.
 * @details Both entry points add half a sector to the caller's phase, so phase
 * 0 puts a vertex on the basis u-axis.
 */
struct SphericalPolygon {
  /**
   * @brief Rasterizes a solid spherical polygon.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam ComputeUVs Whether to compute UV coordinates during distance eval.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the polygon.
   * @param radius Polygon circumradius as a fraction of the hemisphere, in (0, 2).
   * @param sides Number of polygon sides.
   * @param fragment_shader Shader invoked per covered pixel.
   * @param phase Angular phase offset in radians. Added to the query azimuth,
   *        so the shape turns in the opposite sense to Plot's phase for the
   *        same basis and radius.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(ScanRef pipeline, Canvas &canvas, const math::Basis &basis,
                   float radius, int sides, FragmentShaderRef fragment_shader,
                   float phase = 0, bool debug_bb = false) {
    auto res = math::get_antipode(basis, radius);
    float offset = math::PI_F / sides;

    SDF::SphericalPolygon shape(res.first, res.second, sides, phase + offset,
                                radius > 1.0f);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }

  /**
   * @brief Rasterizes a constant-color solid spherical polygon.
   * @tparam SineDistance Use edge-plane distance for the AA band.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PipelineT Plotting pipeline type.
   * @param pipeline Plotting pipeline receiving the final colors.
   * @param canvas Destination canvas.
   * @param basis Orientation basis of the polygon.
   * @param radius Circumradius as a fraction of the hemisphere, in (0, 2).
   * @param sides Number of polygon sides.
   * @param color Fill color and alpha.
   * @param phase Angular phase offset in radians. Added to the query azimuth,
   *        so the shape turns in the opposite sense to Plot's phase for the
   *        same basis and radius.
   * @param debug_bb When true, renders the bounding box for debugging.
   */
  template <int W, int H, bool SineDistance = false, typename PipelineT>
  static void draw_solid(PipelineT &pipeline, Canvas &canvas,
                         const math::Basis &basis, float radius, int sides,
                         const Color4 &color, float phase = 0,
                         bool debug_bb = false) {
    auto res = math::get_antipode(basis, radius);
    float offset = math::PI_F / sides;
    SDF::SphericalPolygon shape(res.first, res.second, sides, phase + offset,
                                radius > 1.0f);
    Scan::rasterize_solid<W, H, SineDistance>(pipeline, canvas, shape, color,
                                              debug_bb);
  }
};

} // namespace Scan
