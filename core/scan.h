/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#ifndef HOLOSPHERE_CORE_SCAN_H_
#define HOLOSPHERE_CORE_SCAN_H_

#include <utility>
#include <functional>
#include "sdf.h"
#include "color.h"
#include "filter.h"
#include "static_circular_buffer.h"
#include "canvas.h"
#include "platform.h"

/**
 * @brief The Scan struct contains volumetric (raster) drawing primitives.
 *
 * General Register Mapping for Scan Primitives:
 *  v0: Normalized parameter t (0-1) or angle
 *  v1: Raw Distance or Supplementary Value
 */
namespace Scan {

/**
 * @brief Processes a single pixel for rasterization.
 */
template <int W, int H, bool ComputeUVs = true,
          typename PipelineT = PipelineRef>
static void process_pixel(int x, int y, const Vector &p, PipelineT &pipeline,
                          Canvas &canvas, const auto &shape,
                          FragmentShaderFn fragment_shader, bool debug_bb,
                          SDF::DistanceResult &result_scratch,
                          Fragment &frag_scratch) {
  uint32_t t_sdf = HS_OS_CYCLES();
  shape.template distance<ComputeUVs>(p, result_scratch);
  hs::g_scan_metrics.sdf_dist += (HS_OS_CYCLES() - t_sdf);

  float d = result_scratch.dist;
  float pixel_width = 2.0f * PI_F / W;
  bool is_solid = shape.is_solid;
  float threshold = is_solid ? pixel_width : 0.0f;

  if (debug_bb || d < threshold) {
    float alpha = 1.0f;

    constexpr bool solid = std::remove_cvref_t<decltype(shape)>::is_solid;
    if constexpr (solid) {
      if (d <= -pixel_width) {
        // FAST PATH: Pixel is fully inside the shape. Skip AA
        alpha = 1.0f;
      } else {
        // SLOW PATH: Pixel is on the boundary. Calculate AA
        float t_aa = 0.5f - d / (2.0f * pixel_width);
        alpha = quintic_kernel(std::max(0.0f, std::min(1.0f, t_aa)));
      }
    } else {
      // Stroke Falloff: Opacity fades over the entire thickness
      if (shape.thickness > 0) {
        alpha = quintic_kernel(-d / shape.thickness);
      }
    }

    if (!debug_bb && alpha <= 0.001f)
      return;

    frag_scratch.pos = p;
    frag_scratch.v0 = result_scratch.t;
    frag_scratch.v1 = result_scratch.raw_dist;
    frag_scratch.v2 = 0.0f;
    frag_scratch.v3 = result_scratch.aux;
    frag_scratch.size = result_scratch.size;
    frag_scratch.age = 0;

    uint32_t t_frag = HS_OS_CYCLES();
    fragment_shader(p, frag_scratch);
    hs::g_scan_metrics.frag_shader += (HS_OS_CYCLES() - t_frag);

    if (debug_bb) {
      frag_scratch.color.color = frag_scratch.color.color.lerp16(
          Pixel(65535, 65535, 65535), 65535 / 2);
      frag_scratch.color.alpha = 1.0f;
      alpha = 1.0f;
    }

    if (frag_scratch.color.alpha > 0.001f) {
      uint32_t t0 = HS_OS_CYCLES();
      pipeline.plot(canvas, x, y, frag_scratch.color.color, frag_scratch.age,
                    frag_scratch.color.alpha * alpha);
      hs::g_scan_metrics.plot += (HS_OS_CYCLES() - t0);
    }
  }
}

/**
 * @brief Shared pixel iteration utility for bounded spherical regions.
 * Iterates y in [y_min, y_max], collects float intervals per row via
 * get_intervals, wraps x coordinates, and calls pixel_fn(wx, y, p) per pixel.
 *
 * @param get_intervals  (int y, auto &out) -> bool. Pushes {float,float}
 *                       intervals via out(start, end). Returns true if
 *                       intervals were produced, false for full-row scan.
 * @param pixel_fn       (int wx, int y, const Vector &p) -> void.
 */
template <int W, int H, typename IntervalFn, typename PixelFn>
static void scan_region(int y_min, int y_max, IntervalFn &&get_intervals,
                        PixelFn &&pixel_fn) {
  if (!TrigLUT<W, H>::initialized)
    TrigLUT<W, H>::init();

  StaticCircularBuffer<std::pair<float, float>, 32> intervals;

  for (int y = y_min; y <= y_max; ++y) {
    bool handled = get_intervals(
        y, [&](float t1, float t2) { intervals.push_back({t1, t2}); });

    if (handled && !intervals.is_empty()) {
      float current_end = -FLT_MAX;
      for (const auto &iv : intervals) {
        if (iv.second <= current_end)
          continue;

        float start = std::max(iv.first, current_end);
        float end = iv.second;
        current_end = end;

        if (end - start >= W) {
          for (int x = 0; x < W; ++x)
            pixel_fn(x, y, pixel_to_vector<W, H>(x, y));
          continue;
        }

        int x1 = static_cast<int>(floorf(start));
        int x2 = static_cast<int>(ceilf(end));
        if (x1 == x2)
          x2++;

        // Modulo-free coordinate wrap
        int wx = x1 % W;
        if (wx < 0)
          wx += W;

        for (int x = x1; x < x2; ++x) {
          pixel_fn(wx, y, pixel_to_vector<W, H>(wx, y));
          wx++;
          if (wx >= W)
            wx -= W;
        }
      }
      intervals.clear();
    } else if (!handled) {
      for (int x = 0; x < W; ++x)
        pixel_fn(x, y, pixel_to_vector<W, H>(x, y));
    }
  }
}

/**
 * @brief Computes bounding sphere y-range and per-row x intervals.
 */
template <int W, int H> struct BoundingSphere {
  int y_min, y_max;
  float center_theta; // longitude of center in pixel units
  float angular_radius;

  BoundingSphere(const Vector &center, float bounds_radius)
      : angular_radius(asinf(std::min(bounds_radius, 1.0f))) {
    PixelCoords center_px = vector_to_pixel<W, H>(center);
    center_theta = center_px.x;
    constexpr int H_VIRT = H + hs::H_OFFSET;
    float center_phi = (center_px.y * PI_F) / (H_VIRT - 1);
    y_min =
        std::max(0, static_cast<int>(phi_to_y<H>(center_phi - angular_radius)));
    y_max = std::min(H - 1, static_cast<int>(ceilf(
                                phi_to_y<H>(center_phi + angular_radius))));
  }

  /// Push a single interval for row y based on longitude span at that latitude.
  template <typename OutFn> bool get_intervals(int y, OutFn &&out) const {
    float phi = y_to_phi<H>(y);
    float sin_phi = sinf(phi);
    float theta_span;
    if (sin_phi < angular_radius / PI_F) {
      // Near pole: span exceeds half the row, scan all columns
      theta_span = static_cast<float>(W);
    } else {
      theta_span = angular_radius / sin_phi * W / (2.0f * PI_F);
    }
    int x_half = std::min(W / 2, static_cast<int>(ceilf(theta_span)) + 1);
    out(center_theta - x_half, center_theta + x_half);
    return true;
  }
};

/**
 * @brief Main rasterization routine for SDF shapes.
 * Scans the bounding box, computes intervals, and executes the shader for valid
 * pixels.
 */
template <int W, int H, bool ComputeUVs = true,
          typename PipelineT = PipelineRef>
static void rasterize(PipelineT &pipeline, Canvas &canvas, const auto &shape,
                      FragmentShaderFn fragment_shader, bool debug_bb = false) {
  bool effective_debug = debug_bb || canvas.debug();
  auto bounds = shape.template get_vertical_bounds<H>();

  // Tier 2: Clamp SDF bounding box to clip region
  const auto &cr = canvas.clip();
  int y_lo =
      bounds.y_min > cr.render_y_start() ? bounds.y_min : cr.render_y_start();
  int y_hi = bounds.y_max < cr.render_y_end() - 1 ? bounds.y_max
                                                  : cr.render_y_end() - 1;
  if (y_lo > y_hi)
    return;

  SDF::DistanceResult result_scratch;
  Fragment frag_scratch;

  scan_region<W, H>(
      y_lo, y_hi,
      [&](int y, auto &&out) {
        return shape.template get_horizontal_intervals<W, H>(y, out);
      },
      [&](int wx, int y, const Vector &p) {
        process_pixel<W, H, ComputeUVs>(wx, y, p, pipeline, canvas, shape,
                                        fragment_shader, effective_debug,
                                        result_scratch, frag_scratch);
      });
}

struct DistortedRing {
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, float thickness, ScalarFn shift_fn,
                   float amplitude, FragmentShaderFn fragment_shader,
                   float phase = 0, bool debug_bb = false) {
    SDF::DistortedRing shape(basis, radius, thickness, shift_fn, amplitude,
                             phase);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }
};

struct PlanarPolygon {
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, int sides, FragmentShaderFn fragment_shader,
                   float phase = 0, bool debug_bb = false) {
    auto res = get_antipode(basis, radius);
    float thickness = res.second * (PI_F / 2.0f);

    SDF::Polygon shape(res.first, res.second, thickness, sides, phase,
                       H + hs::H_OFFSET, H);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }
};

struct Line {
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Vector &v1,
                   const Vector &v2, float thickness,
                   FragmentShaderFn fragment_shader, bool debug_bb = false) {
    SDF::Line shape(v1, v2, thickness);
    Scan::rasterize<W, H>(pipeline, canvas, shape, fragment_shader, debug_bb);
  }
};

struct Ring {
  /**
   * @brief Draws a solid ring using SDF rasterization.
   */
  template <int W, int H, bool ComputeUVs = true,
            typename PipelineT = PipelineRef>
  static void draw(PipelineT &pipeline, Canvas &canvas, const Basis &basis,
                   float radius, float thickness,
                   FragmentShaderFn fragment_shader, float phase = 0,
                   bool debug_bb = false) {
    auto res = get_antipode(basis, radius);
    SDF::Ring shape(res.first, res.second, thickness, phase);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }

  // Overload for Vector normal inputs
  template <int W, int H, bool ComputeUVs = true,
            typename PipelineT = PipelineRef>
  static void draw(PipelineT &pipeline, Canvas &canvas, const Vector &normal,
                   float radius, float thickness,
                   FragmentShaderFn fragment_shader, float phase = 0,
                   bool debug_bb = false) {
    Basis basis = make_basis(Quaternion(), normal);
    draw<W, H, ComputeUVs>(pipeline, canvas, basis, radius, thickness,
                           fragment_shader, phase, debug_bb);
  }
};

struct Circle {
  /**
   * @brief Draws a solid circle (filled ring).
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, FragmentShaderFn fragment_shader,
                   bool debug_bb = false) {
    // Circle is a Ring with inner radius 0
    float th = radius * (PI_F / 2.0f);
    Ring::draw<W, H, ComputeUVs>(pipeline, canvas, basis, 0.0f, th,
                                 fragment_shader, 0, debug_bb);
  }

  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Vector &normal,
                   float radius, FragmentShaderFn fragment_shader,
                   bool debug_bb = false) {
    Basis basis = make_basis(Quaternion(), normal);
    draw<W, H, ComputeUVs>(pipeline, canvas, basis, radius, fragment_shader,
                           debug_bb);
  }
};

struct Point {
  /**
   * @brief Draws a solid point (thick circle).
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Vector &p,
                   float thickness, FragmentShaderFn fragment_shader,
                   bool debug_bb = false) {
    // Scan.Point uses Scan.Ring with radius 0
    Basis basis = make_basis(Quaternion(), p);
    Ring::draw<W, H>(pipeline, canvas, basis, 0.0f, thickness, fragment_shader,
                     0.0f, debug_bb);
  }
};

struct Star {
  /**
   * @brief Draws a solid star shape.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, int sides, FragmentShaderFn fragment_shader,
                   float phase = 0, bool debug_bb = false) {
    auto res = get_antipode(basis, radius);
    SDF::Star shape(res.first, res.second, sides, phase, H, H);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }
};

struct Flower {
  /**
   * @brief Draws a solid flower shape.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, int sides, FragmentShaderFn fragment_shader,
                   float phase = 0, bool debug_bb = false) {
    auto res = get_antipode(basis, radius);
    SDF::Flower shape(res.first, res.second, sides, phase, H + hs::H_OFFSET, H);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }
};

struct SphericalPolygon {
  /**
   * @brief Draws a solid spherical polygon.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, int sides, FragmentShaderFn fragment_shader,
                   float phase = 0, bool debug_bb = false) {
    auto res = get_antipode(basis, radius);
    float offset = PI_F / sides;

    SDF::SphericalPolygon shape(res.first, res.second, sides, phase + offset,
                                H + hs::H_OFFSET, H);
    Scan::rasterize<W, H>(pipeline, canvas, shape, fragment_shader, debug_bb);
  }
};

struct HarmonicBlob {
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, int l, int m,
                   float amplitude, const Quaternion &orientation,
                   HarmonicWaveFn harmonic_fn, FragmentShaderFn fragment_shader,
                   bool debug_bb = false) {
    SDF::HarmonicBlob shape(l, m, amplitude, orientation, harmonic_fn);
    Scan::rasterize<W, H>(pipeline, canvas, shape, fragment_shader, debug_bb);
  }
};

struct Mesh {
  // MeshState overload
  template <int W, int H, typename PipelineT = PipelineRef>
  static void draw(PipelineT &pipeline, Canvas &canvas, const MeshState &mesh,
                   FragmentShaderFn fragment_shader, Arena &scratch_arena,
                   bool debug_bb = false) {
    ScratchScope scope(scratch_arena);
    auto *scratch =
        static_cast<SDF::FaceScratchBuffer *>(scratch_arena.allocate(
            sizeof(SDF::FaceScratchBuffer), alignof(SDF::FaceScratchBuffer)));

    const uint8_t *fc = mesh.get_face_counts_data();
    size_t num_f = mesh.get_face_counts_size();
    const uint16_t *fi = mesh.get_faces_data();
    const uint16_t *fo = mesh.get_face_offsets_data();

    for (size_t i = 0; i < num_f; ++i) {
      size_t count = fc[i];

      std::span<const Vector> verts(mesh.vertices.data(), mesh.vertices.size());
      std::span<const uint16_t> indices(fi + fo[i], count);

      uint32_t t_face = HS_OS_CYCLES();
      SDF::Face shape(verts, indices, 0.0f, *scratch, H + hs::H_OFFSET, H);
      hs::g_scan_metrics.face_setup += (HS_OS_CYCLES() - t_face);

      auto wrapper = [&](const Vector &p, Fragment &f_in) {
        f_in.v2 = static_cast<float>(i);
        fragment_shader(p, f_in);
      };

      uint32_t t_rast = HS_OS_CYCLES();
      Scan::rasterize<W, H, true>(pipeline, canvas, shape, wrapper, debug_bb);
      hs::g_scan_metrics.scan_loop += (HS_OS_CYCLES() - t_rast);
    }
  }
};
/**
 * @brief Full-screen per-pixel shader with SAMPLES× SSAA.
 *
 * Accepts a single callable ShaderFn(const Vector &v) -> Color4
 * that maps a world-space unit vector to a final color.
 * The utility calls it SAMPLES× per pixel at sub-pixel offsets and averages.
 */
struct Shader {
  template <int W, int H, int SAMPLES = 1, typename ShaderFn>
  static void draw(Canvas &canvas, ShaderFn &&shader) {
    constexpr float h_virt_minus_1 = static_cast<float>(H + hs::H_OFFSET - 1);
    constexpr float w_float = static_cast<float>(W);

    if constexpr (SAMPLES == 1) {
      // 1× — single sample at pixel center (no SSAA)
      for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
          Vector v = pixel_to_vector<W, H>(x, y);
          Color4 sample = shader(v);
          canvas(x, y) = sample.color * sample.alpha;
        }
      }
    } else {
      constexpr float inv_samples = 1.0f / SAMPLES;
      // Generalized rotated-grid sub-pixel offsets
      struct SampleOffsets {
        float x[SAMPLES];
        float y[SAMPLES];
      };
      constexpr float eps = 0.5f;
      constexpr SampleOffsets offsets = [&]() constexpr {
        SampleOffsets o{};
        // Rotated grid: distribute samples at equal angular spacing
        for (int i = 0; i < SAMPLES; ++i) {
          // Map index to quadrant-based pattern
          int qx = (i & 1) ? -1 : 1;
          int qy = (i & 2) ? -1 : 1;
          o.x[i] = eps * qx;
          o.y[i] = eps * qy;
        }
        return o;
      }();

      for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
          Color4 accum(Pixel(0, 0, 0), 0.0f);

          for (int i = 0; i < SAMPLES; ++i) {
            float px = static_cast<float>(x) + offsets.x[i];
            float py = static_cast<float>(y) + offsets.y[i];

            float theta = (px * 2.0f * PI_F) / w_float;
            float phi = (py * PI_F) / h_virt_minus_1;
            float sin_phi = sinf(phi);
            Vector v(sin_phi * cosf(theta), cosf(phi), sin_phi * sinf(theta));

            Color4 sample = shader(v);
            sample *= inv_samples;
            accum += sample;
          }

          canvas(x, y) = accum.color * accum.alpha;
        }
      }
    }
  }

  /**
   * @brief Full-screen per-pixel shader with SAMPLES× SSAA.
   *
   * Separates expensive per-pixel work (vertex_shader, called once at pixel
   * center) from per-sub-sample evaluation (fragment_shader, called SAMPLES×).
   */
  template <int W, int H, int SAMPLES = 4>
  static void draw(Canvas &canvas, FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader) {
    constexpr float h_virt_minus_1 = static_cast<float>(H + hs::H_OFFSET - 1);
    constexpr float w_float = static_cast<float>(W);

    Fragment frag_base;

    if constexpr (SAMPLES == 1) {
      const auto &cr = canvas.clip();
      for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y) {
        for (int x = 0; x < W; ++x) {
          Vector center_v = pixel_to_vector<W, H>(x, y);
          frag_base.pos = center_v;
          vertex_shader(frag_base);
          fragment_shader(center_v, frag_base);
          canvas(x, y) = frag_base.color.color;
        }
      }
    } else {
      constexpr float inv_samples = 1.0f / SAMPLES;
      struct SampleOffsets {
        float x[SAMPLES];
        float y[SAMPLES];
      };
      constexpr float eps = 0.5f;
      constexpr SampleOffsets offsets = [&]() constexpr {
        SampleOffsets o{};
        for (int i = 0; i < SAMPLES; ++i) {
          int qx = (i & 1) ? -1 : 1;
          int qy = (i & 2) ? -1 : 1;
          o.x[i] = eps * qx;
          o.y[i] = eps * qy;
        }
        return o;
      }();

      const auto &cr = canvas.clip();
      for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y) {
        for (int x = 0; x < W; ++x) {
          Vector center_v = pixel_to_vector<W, H>(x, y);

          frag_base.pos = center_v;
          vertex_shader(frag_base);

          Color4 accum(Pixel(0, 0, 0), 0.0f);

          for (int i = 0; i < SAMPLES; ++i) {
            float px = static_cast<float>(x) + offsets.x[i];
            float py = static_cast<float>(y) + offsets.y[i];
            float theta = (px * 2.0f * PI_F) / w_float;
            float phi = (py * PI_F) / h_virt_minus_1;
            float sin_phi = sinf(phi);
            Vector v(sin_phi * cosf(theta), cosf(phi), sin_phi * sinf(theta));

            Fragment sub_frag = frag_base;
            sub_frag.pos = v;

            fragment_shader(v, sub_frag);

            Color4 sample = sub_frag.color;
            sample *= inv_samples;
            accum += sample;
          }

          canvas(x, y) = accum.color;
        }
      }
    }
  }
};

/**
 * @brief Generic wrapper that places an SDF in world space via a center
 * point and a rotation quaternion. Satisfies the Volume::draw shape concept.
 *
 * The quaternion q maps local→world: world_p = center + rotate(local_p, q)
 * ray_to_local uses q.inverse() to map world→local.
 */
template <typename SDF> struct TransformedVolume {
  const SDF &sdf;
  Vector center;
  Quaternion q_inv; // precomputed inverse (world→local)

  TransformedVolume(const SDF &sdf, const Vector &center, const Quaternion &q)
      : sdf(sdf), center(center), q_inv(q.inverse()) {}

  /// Transform ray origin and direction from world space to local space.
  std::pair<Vector, Vector> ray_to_local(const Vector &ro,
                                         const Vector &vd) const {
    return {rotate(ro - center, q_inv), rotate(vd, q_inv)};
  }

  /// Delegate to the underlying SDF in local space.
  float distance(const Vector &local_p) const { return sdf.distance(local_p); }
};

/**
 * @brief Raymarch volume renderer with orthographic projection.
 *
 * The render loop is internal: callers provide a Shape with a
 * `float distance(const Vector&) const` method (evaluated per march step)
 * and a fragment shader (evaluated once per hit to populate Fragment
 * registers for shading).
 *
 * Coordinate-space contract:
 *   - `view_dir` is the normalized direction all rays travel (camera → scene).
 *   - Ray origins are computed via orthographic projection: each pixel's
 *     position is projected onto the plane perpendicular to `view_dir`,
 *     then offset backward along `view_dir`.
 *   - `bounds_center` and `view_dir` must both be in physical LED space.
 *   - Filter::World::Orient rotates the *output* position passed to the
 *     canvas, not the ray.
 */
struct Volume {
  /**
   * Shape concept:
   *   std::pair<Vector, Vector> ray_to_local(const Vector &ro, const Vector
   * &vd) const; float distance(const Vector &local_point) const;
   */
  template <int W, int H, typename Shape>
  static void
  draw(PipelineRef pipeline, Canvas &canvas, const Vector &bounds_center,
       float bounds_radius, const Vector &view_dir, const Shape &shape,
       FragmentShaderFn frag_fn, int max_steps = 15, float aa_width = 0.01f) {

    // Normalized view direction
    float vd_len = sqrtf(view_dir.x * view_dir.x + view_dir.y * view_dir.y +
                         view_dir.z * view_dir.z);
    float vd_inv = (vd_len > TOLERANCE) ? 1.0f / vd_len : 1.0f;
    Vector vd(view_dir.x * vd_inv, view_dir.y * vd_inv, view_dir.z * vd_inv);

    // Ray must start behind the farthest extent of the shape
    float start_offset = 1.0f + bounds_radius;

    // Precompute bounds_center projected onto the view plane (⊥ vd)
    float bc_dot_vd = bounds_center.x * vd.x + bounds_center.y * vd.y +
                      bounds_center.z * vd.z;
    Vector bc_proj(bounds_center.x - bc_dot_vd * vd.x,
                   bounds_center.y - bc_dot_vd * vd.y,
                   bounds_center.z - bc_dot_vd * vd.z);
    float bounds_r2 = bounds_radius * bounds_radius;

    // Precompute local-space view direction (shared across all pixels)
    auto [dummy_ro, local_vd] = shape.ray_to_local(Vector(0, 0, 0), vd);

    BoundingSphere<W, H> bounds(bounds_center, bounds_radius);

    // Tier 2: Clamp volume bounds to clip region
    const auto &cr = canvas.clip();
    int vol_y_lo =
        bounds.y_min > cr.render_y_start() ? bounds.y_min : cr.render_y_start();
    int vol_y_hi = bounds.y_max < cr.render_y_end() - 1 ? bounds.y_max
                                                        : cr.render_y_end() - 1;
    if (vol_y_lo > vol_y_hi)
      return;

    scan_region<W, H>(
        vol_y_lo, vol_y_hi,
        [&](int y, auto &&out) { return bounds.get_intervals(y, out); },
        [&](int wx, int y, const Vector &p) {
          // Back-face cull
          float facing = p.x * vd.x + p.y * vd.y + p.z * vd.z;
          if (facing >= 0.0f)
            return;

          // Orthographic ray-sphere cull
          float pp_x = p.x - facing * vd.x;
          float pp_y = p.y - facing * vd.y;
          float pp_z = p.z - facing * vd.z;
          float dx = pp_x - bc_proj.x;
          float dy = pp_y - bc_proj.y;
          float dz = pp_z - bc_proj.z;
          if (dx * dx + dy * dy + dz * dz > bounds_r2)
            return;

          // Orthographic ray origin: outside the unit sphere
          Vector ro(pp_x - vd.x * start_offset, pp_y - vd.y * start_offset,
                    pp_z - vd.z * start_offset);

          // Transform ray to local space once per pixel
          auto [local_ro, unused_ld] = shape.ray_to_local(ro, vd);
          Vector local_p = local_ro;
          Vector closest_local = local_ro;
          float closest_d = 999.0f;

          // --- Sphere tracing in local space ---
          for (int i = 0; i < max_steps; ++i) {
            // Early out: ray has exited the back of the bounding sphere
            if (local_p.x * local_vd.x + local_p.y * local_vd.y +
                    local_p.z * local_vd.z >
                bounds_radius)
              break;

            float d = shape.distance(local_p);

            if (d < closest_d) {
              closest_d = d;
              closest_local = local_p;
            }

            if (d < -aa_width)
              break;

            float advance = std::max(d * 0.9f, 1e-5f);
            local_p = Vector(local_p.x + local_vd.x * advance,
                             local_p.y + local_vd.y * advance,
                             local_p.z + local_vd.z * advance);
          }

          if (closest_d >= aa_width)
            return;

          // --- Fragment shading ---
          Fragment frag;
          frag.pos = closest_local;
          frag.size = closest_d;
          frag_fn(closest_local, frag);

          // One-sided AA with quintic kernel
          float hit_threshold = aa_width * 0.1f;
          float edge_alpha;

          if (closest_d <= hit_threshold) {
            // FAST PATH: Solid hit. No probe needed.
            edge_alpha = 1.0f;
          } else {
            // SLOW PATH: We are in the fuzzy AA border.
            // Calculate standard AA alpha...
            edge_alpha = quintic_kernel(1.0f - (closest_d - hit_threshold) /
                                                   (aa_width - hit_threshold));

            // ...and fire the occlusion probe to see if this is a false halo!
            Vector probe = local_p;
            float probed = 0.0f;
            for (int i = 0; i < 4; ++i) {
              float pd = shape.distance(probe);
              if (pd < hit_threshold) {
                edge_alpha = 1.0f; // Punched through the halo, solid hit!
                break;
              }
              float step = std::max(pd * 0.9f, bounds_radius * 0.15f);
              probe = Vector(probe.x + local_vd.x * step,
                             probe.y + local_vd.y * step,
                             probe.z + local_vd.z * step);
              probed += step;
              if (probed > bounds_radius)
                break;
            }
          }

          if (frag.color.alpha * edge_alpha > 0.001f) {
            pipeline.plot(canvas, p, frag.color.color, 0.0f,
                          frag.color.alpha * edge_alpha);
          }
        });
  }
};

} // namespace Scan
#endif // HOLOSPHERE_CORE_SCAN_H_
