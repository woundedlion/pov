/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file
 * @brief Scan admission, shape certificates and shader validation. */

#include "core/math/scalar_kernels.h"
#include <algorithm>
#include <optional>
#include <type_traits>
#include "core/render/sdf.h"
#include "core/render/shading.h"
#include "core/render/render_policy.h"
#include "core/color/color.h"
#include "core/render/filter/pipeline.h"
#include "core/containers/static_circular_buffer.h"
#include "core/render/canvas.h"
#include "core/platform/platform.h"
#include "core/render/scan/column_basis.h"
#include "core/render/spacetime.h"

namespace Scan {

/** @brief Plots an admitted sample through the sink's unchecked entry when available.
 * @tparam PipelineT Plotting sink type.
 * @param pipeline Plotting sink.
 * @param canvas Destination canvas.
 * @param x Admitted column.
 * @param y Admitted row.
 * @param color Sample color.
 * @param age Sample age.
 * @param alpha Sample coverage. */
template <typename PipelineT>
HS_HOT_INLINE inline void plot_admitted(PipelineT &pipeline, Canvas &canvas,
                                        int x, int y, const Pixel &color,
                                        float age, float alpha) {
  if constexpr (requires {
                  pipeline.plot_in_bounds(canvas, x, y, color, age, alpha);
                })
    pipeline.plot_in_bounds(canvas, x, y, color, age, alpha);
  else
    pipeline.plot(canvas, x, y, color, age, alpha);
}

/** @brief Loads shader registers from one distance query.
 * @param fragment Receives initialized shader inputs.
 * @param result Distance query result.
 * @param position Surface point.
 * @param v2 Coverage register value. */
HS_HOT_INLINE inline void load_fragment(Fragment &fragment,
                                        const SDF::DistanceResult &result,
                                        const math::Vector &position,
                                        float v2) {
  fragment.color = Color4(0, 0, 0, 0);
  fragment.pos = position;
  fragment.v0 = result.t;
  fragment.v1 = result.raw_dist;
  fragment.v2 = v2;
  fragment.v3 = result.aux;
  fragment.size = result.size;
  fragment.age = 0;
}

/** @brief Coverage/alpha at or below which a fragment is dropped unplotted. */
inline constexpr float MIN_ALPHA = 0.001f;

/**
 * @brief Columns one shade may cover at a row's colatitude.
 * @param sin_phi Sine of the row's colatitude; sign is ignored.
 * @return Run length in columns, 1 when decimation is off or unwarranted.
 * @details `Render::pole_lod_aggressiveness / sin(phi)`, clamped to Render::POLE_LOD_MAX_RUN.
 *          Returns 1 at aggressiveness 0.
 */
inline int pole_lod_run(float sin_phi) {
  const float lod = Render::pole_lod_aggressiveness;
  if (lod <= 0.0f)
    return 1;
  const float a = sin_phi < 0.0f ? -sin_phi : sin_phi;
  if (a < 1e-6f)
    return Render::POLE_LOD_MAX_RUN;
  const int run = static_cast<int>(lod / a);
  if (run <= 1)
    return 1;
  return run < Render::POLE_LOD_MAX_RUN ? run : Render::POLE_LOD_MAX_RUN;
}

/**
 * @brief First canvas-aligned block anchor at or after a column.
 * @param x Column a walk starts at; non-negative.
 * @param stride Columns per block.
 * @return The anchor. A walk advances it by stride from here, so block heads
 *         are selected without a per-column divide.
 */
__attribute__((always_inline)) inline int pole_lod_block_anchor(int x,
                                                                int stride) {
  return x + (stride - x % stride) % stride;
}

/**
 * @brief Latitude-parallel arc from a pole-LOD block's first column to its
 *        last.
 * @tparam W Canvas width in pixels.
 * @param run Columns in the block; 1 yields no arc.
 * @param sin_phi Sine of the row's colatitude.
 * @return The block's longitude span foreshortened by sin(phi), in radians.
 */
template <int W>
__attribute__((always_inline)) inline float pole_lod_arc(int run,
                                                         float sin_phi) {
  return static_cast<float>(run - 1) * (2.0f * math::PI_F / W) * sin_phi;
}

/**
 * @brief Factor from an angular step to the units a shape's distance() reports.
 * @details A block probe's slack is an arc, so a shape reporting in another
 * unit scales it. A shape without an overload does not compile.
 * @return 1 for SDF::Ring (no rescaling).
 */
__attribute__((always_inline)) inline float report_stretch(const SDF::Ring &) {
  return 1.0f;
}
/**
 * @brief Report stretch of an SDF::DistortedRing.
 * @return 1 (no rescaling).
 */
__attribute__((always_inline)) inline float
report_stretch(const SDF::DistortedRing &) {
  return 1.0f;
}
/**
 * @brief Report stretch of an SDF::FlatDistortedRing.
 * @return 1 (no rescaling).
 */
__attribute__((always_inline)) inline float
report_stretch(const SDF::FlatDistortedRing &) {
  return 1.0f;
}
/**
 * @brief Report stretch of an SDF::PlanarPolygon.
 * @return 1 (no rescaling).
 */
__attribute__((always_inline)) inline float
report_stretch(const SDF::PlanarPolygon &) {
  return 1.0f;
}
/**
 * @brief Report stretch of an SDF::SphericalPolygon.
 * @return 1 (no rescaling).
 */
__attribute__((always_inline)) inline float
report_stretch(const SDF::SphericalPolygon &) {
  return 1.0f;
}
/**
 * @brief Report stretch of an SDF::Star.
 * @return 1 (no rescaling).
 */
__attribute__((always_inline)) inline float report_stretch(const SDF::Star &) {
  return 1.0f;
}
/**
 * @brief Report stretch of an SDF::Flower.
 * @return 1 (no rescaling).
 */
__attribute__((always_inline)) inline float
report_stretch(const SDF::Flower &) {
  return 1.0f;
}
/**
 * @brief Report stretch of an SDF::Line.
 * @return 1 (no rescaling).
 */
__attribute__((always_inline)) inline float report_stretch(const SDF::Line &) {
  return 1.0f;
}
/** @brief SDF::PreparedFace reports gnomonic-plane distance, whose angular-step stretch
 *         grows as 1 + r^2 out to the farthest probe it reports.
 *  @param shape Face whose report radius bounds the stretch.
 *  @return 1 + `SDF::PreparedFace::report_radius_sq()`. */
__attribute__((always_inline)) inline float
report_stretch(const SDF::PreparedFace &shape) {
  return 1.0f + shape.report_radius_sq();
}
/** @brief Reports projection stretch through a checked face handle.
 * @param shape Prepared face handle.
 * @return Prepared geometry's angular-step stretch. */
inline float report_stretch(const SDF::Face &shape) {
  return report_stretch(shape.prepared());
}
// Forward-declared so a nested composite's recursive call sees every overload.
template <template <typename, typename> class Csg, typename A, typename B>
  requires(std::is_same_v<Csg<A, B>, SDF::Union<A, B>> ||
           std::is_same_v<Csg<A, B>, SDF::SmoothUnion<A, B>> ||
           std::is_same_v<Csg<A, B>, SDF::Intersection<A, B>> ||
           std::is_same_v<Csg<A, B>, SDF::Subtract<A, B>>)
inline float report_stretch(const Csg<A, B> &shape);
template <typename Shape>
inline float report_stretch(const SDF::AngularRepeat<Shape> &shape);

/**
 * @brief Report stretch of a binary CSG shape: the larger operand stretch.
 * @tparam Csg Union, SmoothUnion, Intersection or Subtract.
 * @tparam A First operand shape.
 * @tparam B Second operand shape.
 * @param shape Composite whose operands are queried.
 * @return max(report_stretch(shape.a), report_stretch(shape.b)).
 */
template <template <typename, typename> class Csg, typename A, typename B>
  requires(std::is_same_v<Csg<A, B>, SDF::Union<A, B>> ||
           std::is_same_v<Csg<A, B>, SDF::SmoothUnion<A, B>> ||
           std::is_same_v<Csg<A, B>, SDF::Intersection<A, B>> ||
           std::is_same_v<Csg<A, B>, SDF::Subtract<A, B>>)
__attribute__((always_inline)) inline float
report_stretch(const Csg<A, B> &shape) {
  return fmaxf(report_stretch(shape.a), report_stretch(shape.b));
}
/**
 * @brief Report stretch of an SDF::AngularRepeat: that of its repeated shape.
 * @tparam Shape Repeated shape.
 * @param shape Repeat whose inner shape is queried.
 * @return report_stretch(shape.shape).
 */
template <typename Shape>
__attribute__((always_inline)) inline float
report_stretch(const SDF::AngularRepeat<Shape> &shape) {
  return report_stretch(shape.shape);
}

/**
 * @brief Factor from a pole-LOD block's angular slack to the units a shape's
 *        distance() reports.
 * @details The shape's change-per-arc factor (SDF::arc_stretch) times its
 * report stretch (report_stretch). A probe farther from the surface than the
 * scaled arc cannot change side anywhere in the block. Only block_slack_scale
 * constructs one.
 */
class SlackScale {
public:
  /**
   * @brief Scales an arc.
   * @param arc Angular slack, radians.
   * @return Slack in the shape's report units.
   */
  HS_HOT_INLINE float operator()(float arc) const {
    return arc * arc_factor * report;
  }
  /**
   * @brief Scales the sum of a block's arc and the basis turn across it.
   * @param arc Latitude-parallel block arc, radians.
   * @param turn Extra arc the evaluated directions turn, radians.
   * @return Slack in the shape's report units.
   */
  HS_HOT_INLINE float operator()(float arc, float turn) const {
    return (arc * arc_factor + turn * arc_factor) * report;
  }

private:
  float arc_factor; ///< SDF::arc_stretch of the shape.
  float report;     ///< report_stretch of the shape.

  /**
   * @brief Captures both factors.
   * @param arc_factor SDF::arc_stretch of the shape.
   * @param report report_stretch of the shape.
   */
  HS_HOT_INLINE SlackScale(float arc_factor, float report)
      : arc_factor(arc_factor), report(report) {}

  template <typename ShapeT>
  friend SlackScale block_slack_scale(const ShapeT &shape);
};

/**
 * @brief Pole-LOD slack factor of a shape.
 * @tparam ShapeT Shape the probe reads distance() from.
 * @param shape Shape whose report stretch scales the slack.
 * @return The shape's slack factor.
 */
template <typename ShapeT>
HS_HOT_INLINE SlackScale block_slack_scale(const ShapeT &shape) {
  return SlackScale(SDF::arc_stretch<std::remove_cvref_t<ShapeT>>,
                    report_stretch(shape));
}

/**
 * @brief Whether a walk over a shape may settle a whole pole-LOD block from one
 *        probe of it.
 * @tparam ShapeT Shape the probe reads distance() from.
 * @details False with the knob compiled out, false for a shape with no reject
 * band (`SDF::reject_margin` <= 0), and false for a shape whose distance() has
 * no finite change-per-arc factor: no slack makes one probe vouch for its
 * neighbours, so the walk stays per column.
 */
template <typename ShapeT>
inline constexpr bool pole_lod_blocks =
    Render::POLE_LOD_ENABLED &&
    SDF::reject_margin<std::remove_cvref_t<ShapeT>> > 0.0f &&
    SDF::arc_stretch<std::remove_cvref_t<ShapeT>> < SDF::ARC_STRETCH_UNBOUNDED;

/**
 * @brief Whether a probe's clearance report bounds a whole pole-LOD block.
 * @tparam ShapeT Shape the probe read distance() from.
 * @param threshold Clearance the walk tests the probe against.
 * @param block_slack Extra clearance demanded of a block probe.
 * @return True when a report at or above threshold + block_slack holds for
 *         every column in the block.
 * @details Past its reject band (SDF::reject_margin) a shape may report
 * FAR_SENTINEL, which bounds nothing; the block test holds only where the
 * widened threshold sits inside the static margin.
 */
template <typename ShapeT>
__attribute__((always_inline)) inline bool
probe_bounds_block(float threshold, float block_slack) {
  return threshold + block_slack <=
         SDF::reject_margin<std::remove_cvref_t<ShapeT>>;
}

/**
 * @brief Extra clearance a block probe must carry, in the shape's report units.
 * @tparam W Canvas width in columns.
 * @tparam ShapeT Shape the probe reads distance() from.
 * @param run Columns in the block; 1 yields no slack.
 * @param p_y Probe's y coordinate on the unit sphere (cos of the colatitude).
 * @param shape Shape the probe reads distance() from.
 * @return Slack in the shape's report units.
 */
template <int W, typename ShapeT>
__attribute__((always_inline)) inline float
pole_lod_block_slack(int run, float p_y, const ShapeT &shape) {
  const float sin_phi = sqrtf(fmaxf(0.0f, 1.0f - p_y * p_y));
  return block_slack_scale(shape)(pole_lod_arc<W>(run, sin_phi));
}

/**
 * @brief Whether one probe settles a whole pole-LOD block on one side.
 * @tparam ShapeT Shape the probe read distance() from.
 * @param clearance Probe's distance from the surface on the side under test;
 *        negate a solid's report to test the interior side.
 * @param threshold Clearance the walk tests the probe against.
 * @param block_slack Extra clearance demanded of a block probe.
 * @return True when the probe holds for every column in the block.
 */
template <typename ShapeT>
__attribute__((always_inline)) inline bool
pole_lod_block_settles(float clearance, float threshold, float block_slack) {
  return clearance >= threshold + block_slack &&
         probe_bounds_block<ShapeT>(threshold, block_slack);
}

/** @brief How one probe settles a pole-LOD block. */
enum class BlockSettle : uint8_t {
  PER_COLUMN, ///< Unsettled: the walk probes each column.
  CLEAR,      ///< Every column is clear: skip the block.
  INTERIOR,   ///< Every column is solid interior: splat the probe's result.
};

/**
 * @brief Settles a pole-LOD block from its first column's probe.
 * @tparam Solid Whether the shape is solid, so its interior can settle too.
 * @tparam ShapeT Shape the probe read distance() from.
 * @param d Probe's signed distance report.
 * @param threshold Clearance the walk tests the probe against.
 * @param block_slack Prepared clearance for the offered block.
 * @param shape Shape the probe read.
 * @param probe Optional prepared distance and row certificate.
 * @return The side the probe settles the block on, or PER_COLUMN.
 */
template <bool Solid, typename ShapeT, typename Probe = SDF::DirectDistance>
__attribute__((always_inline)) inline BlockSettle
pole_lod_settle_with_slack(float d, float threshold, float block_slack,
                           [[maybe_unused]] const ShapeT &shape,
                           const Probe &probe = {}) {
  const auto settles = [&](float clearance) {
    if constexpr (requires {
                    probe.block_allowed;
                    probe.threshold;
                  })
      return probe.block_allowed && clearance >= probe.threshold;
    else
      return pole_lod_block_settles<ShapeT>(clearance, threshold, block_slack);
  };
  if (settles(d))
    return BlockSettle::CLEAR;
  if (Solid && settles(-d))
    return BlockSettle::INTERIOR;
  return BlockSettle::PER_COLUMN;
}

/**
 * @brief Settles a pole-LOD block from its first column's probe.
 * @tparam W Canvas width in pixels.
 * @tparam Solid Whether the shape is solid, so its interior can settle too.
 * @tparam ShapeT Shape the probe read distance() from.
 * @param d Probe's signed distance report.
 * @param threshold Clearance the walk tests the probe against.
 * @param run Block width in columns.
 * @param p_y Probe's y coordinate.
 * @param shape Shape the probe read.
 * @return The side the probe settles the block on, or PER_COLUMN.
 */
template <int W, bool Solid, typename ShapeT>
__attribute__((always_inline)) inline BlockSettle
pole_lod_settle_block(float d, float threshold, int run, float p_y,
                      const ShapeT &shape) {
  return pole_lod_settle_with_slack<Solid>(
      d, threshold, pole_lod_block_slack<W>(run, p_y, shape), shape);
}

/**
 * @brief Validates that the canvas is the size the scan was instantiated for.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @param canvas Destination canvas.
 * @details Canvas dimensions must match the trig LUT extents and the
 * pipeline's wrap period and framebuffer stride.
 */
template <int W, int H>
HS_NOINLINE_NOCLONE inline void check_canvas_dims(const Canvas &canvas) {
  HS_CHECK(canvas.width() == W && canvas.height() == H,
           "canvas size differs from the scan's W/H");
}

/**
 * @brief Validates that a direct-raster pipeline was prepared for this canvas.
 * @param pipeline Plotting pipeline the draw writes through.
 * @param canvas Destination canvas.
 * @details A direct-raster sink writes through a cached framebuffer base; the
 * canvas double-buffers, so a stale base is the buffer the display is scanning
 * out. Compiles away for a pipeline without the hook.
 */
template <typename PipelineT>
inline void check_pipeline_prepared(PipelineT &pipeline, Canvas &canvas) {
  if constexpr (requires { pipeline.prepared_for(canvas); })
    HS_CHECK(pipeline.prepared_for(canvas),
             "direct raster pipeline not prepared for this canvas");
}

/**
 * @brief Source footprint for a scan whose pipeline may move world positions.
 * @param pipeline Destination pipeline, typed or erased.
 * @param canvas Destination canvas.
 * @return Full source domain for world transforms, otherwise the render clip.
 */
template <int W, int H, typename PipelineT>
inline const ClipRegion &source_clip(const PipelineT &pipeline,
                                     const Canvas &canvas) {
  if constexpr (requires { pipeline.world_transform_is_identity; }) {
    if (!pipeline.world_transform_is_identity) {
      static constexpr ClipRegion FULL{0, H, 0, W, 0, W, H};
      return FULL;
    }
  }
  return canvas.clip();
}

/**
 * @brief Validates that a type-erased fragment shader refers to a callable.
 * @param fragment_shader Shader the draw invokes per pixel.
 * @details Traps with a scan-specific diagnostic before invoking the shader.
 */
HS_NOINLINE_NOCLONE inline void
check_fragment_shader(FragmentShaderRef fragment_shader) {
  HS_CHECK(fragment_shader, "Scan: fragment_shader must be non-null");
}

/**
 * @brief Validates that a type-erased time-aware shader refers to a callable.
 * @param fragment_shader Shader the draw invokes per pixel.
 */
HS_FLASH_MEMBER inline void
check_fragment_shader(TimedShaderRef fragment_shader) {
  HS_CHECK(fragment_shader,
           "Scan: time-aware fragment_shader must be non-null");
}

} // namespace Scan
