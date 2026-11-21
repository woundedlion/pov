/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file four_stroke.h
 * @brief Boundary-aware forward shading on a fixed-time output footprint. */

#include "core/render/four_stroke_contour.h"
#include "core/render/scan/raster.h"
#include "core/render/scan/four_stroke_star.h"

namespace Scan {

/** @brief Invokes a fragment stage with the resolved contour registers.
 * @details v0 retains producer progress, v1 is signed contour distance in output pixels
 * (negative inside), v2 is tap coverage, and v3 is signed outward alignment.
 * The base material has already consumed its shape-specific registers.
 * @param frag Base material outputs and position.
 * @param edge Producer's resolved geometry sample.
 * @param shader Fragment stage consuming the contour register convention. */
inline void shade_boundary_fragment(Fragment &frag, FourStroke::Boundary edge,
                                    const auto &shader) {
  frag.v1 = edge.inside ? -edge.depth : edge.depth;
  frag.v2 = edge.inside ? 1.0f : 0.0f;
  frag.v3 = edge.alignment;
  shader(frag.pos, frag);
}

/** @brief Fragment material consuming the contour register convention. */
struct FourStrokeShader {
  const FourStroke::Controls &controls; ///< Draw-wide material controls.
  const FourStroke::Colors &beat;       ///< Prepared colors at the column time.

  /** @brief Modulates the wrapped fragment's straight RGB.
   * @param position Producer's surface position.
   * @param frag Wrapped material and contour registers. */
  void operator()([[maybe_unused]] const math::Vector &position,
                  Fragment &frag) const {
    frag.color.color = FourStroke::shade(
        frag.color.color, {frag.v2 > 0, -frag.v1, frag.v3}, beat, controls);
  }
};

/** @brief Pixel resolved from one visible owner per covered tap. */
struct StrokePixel {
  FourStroke::Resolved rgba; ///< Premultiplied linear color and covered alpha.
  float age;                 ///< Alpha-weighted fragment age.
};

/** @brief Evaluates a visible fragment material on the reference footprint.
 * @param x,y Output pixel center.
 * @param time Presentation time shared by its taps.
 * @param beat Prepared role colors.
 * @param controls Four-stroke material settings.
 * @param query Prepared query of center-relative offsets and normal alignment.
 * @param alignment Alignment at absolute tap coordinates and presentation time.
 * @param material Loads and shades the wrapped fragment at absolute coordinates.
 * @param boundary_shader Optional replacement consuming the same contour registers.
 * @return Premultiplied color, coverage and age. */
inline StrokePixel four_stroke_pixel(int x, int y, uint64_t time,
                                     const FourStroke::Colors &beat,
                                     const FourStroke::Controls &controls,
                                     const auto &query, const auto &alignment,
                                     const auto &material,
                                     FragmentShaderRef boundary_shader = {}) {
  float age = 0;
  auto resolved = FourStroke::resolve_pixel([&](float dx, float dy) -> Color4 {
    const float sx = x + dx, sy = y + dy;
    const auto edge = query(dx, dy, [&](FourStroke::Direction n) {
      return alignment(sx, sy, time, n);
    });
    if (!edge.inside)
      return {};
    Fragment frag;
    material(sx, sy, time, frag);
    if (boundary_shader)
      shade_boundary_fragment(frag, edge, boundary_shader);
    else
      shade_boundary_fragment(frag, edge, FourStrokeShader{controls, beat});
    age += frag.age * frag.color.alpha;
    return frag.color;
  });
  return {resolved, resolved.alpha > 0 ? age / (16 * resolved.alpha) : 0};
}

/** @brief Binds an output-pixel heading to a draw.
 * @param direction Unit heading or normalized instance override.
 * @return Signed alignment function. */
inline auto four_stroke_output_direction(FourStroke::Direction direction) {
  return [direction](float, float, uint64_t, FourStroke::Direction n) {
    return n.x * direction.x + n.y * direction.y;
  };
}

/** @brief Binds a physical arrow in view coordinates to a draw.
 * @tparam W,H Output dimensions.
 * @param arrow Returns the sphere or rotated object arrow at the sample time.
 * @return Physical tangent alignment using calibrated latitude spacing. */
template <int W, int H> auto four_stroke_tangent_direction(const auto &arrow) {
  return [&arrow](float x, float y, uint64_t time, FourStroke::Direction n) {
    constexpr float STEP = 2 * math::PI_F / W;
    const float phi = math::DisplayGeometry<H>::row_to_phi(y);
    const float st = std::sin(x * STEP), ct = std::cos(x * STEP);
    const float sp = std::sin(phi), cp = std::cos(phi);
    const math::Vector p{sp * ct, cp, sp * st};
    const math::Vector dx{-sp * st * STEP, 0, sp * ct * STEP};
    const math::Vector dy = math::Vector{cp * ct, -sp, cp * st} *
                            math::DisplayGeometry<H>::RADIANS_PER_ROW;
    return FourStroke::tangent_alignment(p, dx, dy, arrow(p, time), n);
  };
}

/** @brief Projects a prepared spherical SDF into its sampled output contour.
 * @tparam W,H Output dimensions.
 * @param pipeline Prepared destination pipeline.
 * @param canvas Destination canvas.
 * @param select_column Fixes the center column for all footprint samples.
 * @tparam Shape Final signed-distance producer type.
 * @tparam Shader Wrapped material shader type.
 * @tparam ColumnFn Fixed-time geometry context selector.
 * @param shape Final signed-distance producer, including any CSG.
 * @param shader Original material receiving the producer's fragment registers.
 * @param view Fixed-time mapping from logical pixel coordinates to the shape's frame.
 * @param alignment Signed alignment for a unit outward normal at a tap.
 * @param controls Four-stroke material controls.
 * @param t0,t1 Integer presentation endpoints.
 * @param boundary_shader Optional replacement consuming the same contour registers. */
template <int W, int H, typename Shape, typename Shader, typename ColumnFn>
void rasterize_four_stroke_sdf(auto &pipeline, Canvas &canvas,
                               const Shape &shape, const Shader &shader,
                               const auto &view, const auto &alignment,
                               const FourStroke::Controls &controls,
                               uint64_t t0, uint64_t t1,
                               FragmentShaderRef boundary_shader,
                               const ColumnFn &select_column) {
  if (controls.strength == 0) {
    rasterize<W, H>(pipeline, canvas, shape, shader);
    return;
  }
  check_pipeline_prepared(pipeline, canvas);
  const auto field = [&](int x, int y, uint64_t time) {
    SDF::DistanceResult hit;
    shape.template distance<false>(view(x, y, time), hit);
    return hit.dist;
  };
  const auto shade = [&](float x, float y, uint64_t time, Fragment &f) {
    const auto p = view(x, y, time);
    SDF::DistanceResult hit;
    shape.template distance<true>(p, hit);
    load_fragment(f, hit, p, 0);
    shader(p, f);
  };
  const auto emit = [&](int x, int y, const StrokePixel &result) {
    const float alpha = result.rgba.alpha;
    if (alpha <= 0)
      return;
    const Pixel &c = result.rgba.premultiplied;
    const Pixel straight{round_linear_channel(c.r / alpha),
                         round_linear_channel(c.g / alpha),
                         round_linear_channel(c.b / alpha)};
    pipeline.plot(canvas, x, y, straight, result.age, alpha);
  };
  if constexpr (std::is_same_v<Shape, SDF::Star>) {
    if (!canvas.column_basis()) {
      ScratchScope storage(scratch_arena_b);
      auto &contour = *scratch_arena_b.make<ProjectedStar<W, H>>();
      if (contour.prepare(shape)) {
        const auto &clip = canvas.clip();
        const int y0 = std::max(clip.render_y_start(), contour.begin_row());
        const int y1 = std::min(clip.render_y_end(), contour.end_row());
        walk_window_columns<W>(
            canvas.render_x_clip(), canvas.time_window(), [&](int x) {
              const float ux = contour.unwrap(float(x));
              if (!contour.admits_x(ux))
                return;
              select_column(x);
              const uint64_t time = FourStroke::column_time(
                  t0, t1, x, W, Spacetime::flat_surface());
              const auto &beat = FourStroke::colors(time);
              for (int y = y0; y < y1; ++y) {
                const auto query = [&](float dx, float dy, const auto &align) {
                  return contour.query(ux + dx, y + dy, align);
                };
                emit(x, y,
                     four_stroke_pixel(x, y, time, beat, controls, query,
                                       alignment, shade, boundary_shader));
              }
            });
        return;
      }
    }
  }
  math::Vector axis{};
  float cap_radius = 0;
  bool bounded = false;
  if constexpr (requires { shape.bounding_cap(axis, cap_radius); })
    bounded = shape.bounding_cap(axis, cap_radius);
  const float halo =
      2 * math::PI_F / W + math::DisplayGeometry<H>::RADIANS_PER_ROW;
  const float min_dot =
      std::cos(std::min(math::PI_F, cap_radius + halo + 0.001f));
  const auto &clip = canvas.clip();
  walk_window_columns<W>(
      canvas.render_x_clip(), canvas.time_window(), [&](int x) {
        select_column(x);
        const uint64_t time =
            FourStroke::column_time(t0, t1, x, W, Spacetime::flat_surface());
        const auto &beat = FourStroke::colors(time);
        for (int y = clip.render_y_start(); y < clip.render_y_end(); ++y) {
          if (bounded && math::dot(view(x, y, time), axis) < min_dot)
            continue;
          float center[3][3];
          bool covered = false;
          for (int dy = -1; dy <= 1; ++dy)
            for (int dx = -1; dx <= 1; ++dx) {
              const float f = field(math::wrap(x + dx, W), y + dy, time);
              center[dy + 1][dx + 1] = f;
              covered |= f <= 0;
            }
          if (!covered)
            continue;
          FourStroke::ContourTile tile;
          bool solid = true;
          tile.prepare([&](int dx, int dy) {
            const float f = std::abs(dx) <= 1 && std::abs(dy) <= 1
                                ? center[dy + 1][dx + 1]
                                : field(math::wrap(x + dx, W), y + dy, time);
            solid &= f <= 0;
            return f;
          });
          const auto query = [&](float dx, float dy, auto align) {
            if (solid)
              return FourStroke::Boundary{true, INFINITY, 0};
            return tile.query(dx, dy, align);
          };
          emit(x, y,
               four_stroke_pixel(x, y, time, beat, controls, query, alignment,
                                 shade, boundary_shader));
        }
      });
}

/** @brief Four-stroke binding shared by area and path draws. */
struct FourStrokeMaterial {
  FourStroke::Controls controls; ///< Direction, contrast, width and strength.
  FragmentShaderRef
      base; ///< Wrapped shader, invoked on producer-specific inputs.
  FragmentShaderRef
      boundary{}; ///< Optional replacement for the standard register shader.
  std::optional<math::Vector>
      arrow{}; ///< Physical direction in the view frame.
};

/** @brief Rasterizes an SDF with a boundary-aware fragment material.
 * @tparam W,H Output dimensions.
 * @param pipeline Prepared plotting sink.
 * @param canvas Destination canvas.
 * @param shape Signed-distance producer.
 * @param material Base material and the Four-Stroke register stage. */
template <int W, int H, typename Shape>
void rasterize(auto &pipeline, Canvas &canvas, const Shape &shape,
               FourStrokeMaterial material) {
  const auto interval = canvas.presentation_interval();
  const auto output = four_stroke_output_direction(
      FourStroke::heading(material.controls.heading));
  const auto arrow = [&](const auto &, uint64_t) { return *material.arrow; };
  const auto tangent = four_stroke_tangent_direction<W, H>(arrow);
  const auto alignment = [&](float x, float y, uint64_t time,
                             FourStroke::Direction n) {
    return material.arrow ? tangent(x, y, time, n) : output(x, y, time, n);
  };
  int center = 0;
  const auto view = [&](float x, float y, uint64_t) {
    const float phi = math::DisplayGeometry<H>::row_to_phi(y);
    const auto *basis = canvas.column_basis();
    if (!basis || !basis->azimuth || Spacetime::flat_surface())
      return math::Vector(math::Spherical(x * (2 * math::PI_F / W), phi));
    const int k = center % (W / 2);
    const auto radial = center < W / 2 ? basis->azimuth[k] : -basis->azimuth[k];
    const auto &up = basis->up[k];
    const auto tangent = math::cross(radial, up);
    const float d = (x - center) * (2 * math::PI_F / W);
    return (radial * std::cos(d) + tangent * std::sin(d)) * std::sin(phi) +
           up * std::cos(phi);
  };
  const auto save_column = [&](int x) { center = x; };
  rasterize_four_stroke_sdf<W, H>(pipeline, canvas, shape, material.base, view,
                                  alignment, material.controls, interval.t0,
                                  interval.t1, material.boundary, save_column);
}

} // namespace Scan
