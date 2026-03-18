/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <utility>
#include <functional>
#include "sdf.h"
#include "color.h"
#include "filter.h"
#include "static_circular_buffer.h"
#include "canvas.h"

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
template <int W, int H, bool ComputeUVs = true>
static void process_pixel(int x, int y, const Vector &p, PipelineRef pipeline,
                          Canvas &canvas, const auto &shape,
                          FragmentShaderFn fragment_shader, bool debug_bb,
                          SDF::DistanceResult &result_scratch,
                          Fragment &frag_scratch) {
  shape.template distance<ComputeUVs>(p, result_scratch);
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

    fragment_shader(p, frag_scratch);

    if (debug_bb) {
      frag_scratch.color.color = frag_scratch.color.color.lerp16(
          Pixel(65535, 65535, 65535), 65535 / 2);
      frag_scratch.color.alpha = 1.0f;
      alpha = 1.0f;
    }

    if (frag_scratch.color.alpha > 0.001f) {
      pipeline.plot(canvas, x, y, frag_scratch.color.color, frag_scratch.age,
                    frag_scratch.color.alpha * alpha);
    }
  }
}

/**
 * @brief Main rasterization routine for SDF shapes.
 * Scans the bounding box, computes intervals, and executes the shader for valid
 * pixels.
 */
template <int W, int H, bool ComputeUVs = true>
static void rasterize(PipelineRef pipeline, Canvas &canvas, const auto &shape,
                      FragmentShaderFn fragment_shader, bool debug_bb = false) {
  bool effective_debug = debug_bb || canvas.debug();
  auto bounds = shape.template get_vertical_bounds<H>();

  if (!TrigLUT<W, H>::initialized)
    TrigLUT<W, H>::init();

  StaticCircularBuffer<std::pair<float, float>, 32> intervals;
  SDF::DistanceResult result_scratch;
  Fragment frag_scratch;

  for (int y = bounds.y_min; y <= bounds.y_max; ++y) {

    bool handled = shape.template get_horizontal_intervals<W, H>(
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
          for (int x = 0; x < W; ++x) {
            Vector p = pixel_to_vector<W, H>(x, y);
            process_pixel<W, H, ComputeUVs>(x, y, p, pipeline, canvas, shape,
                                            fragment_shader, effective_debug,
                                            result_scratch, frag_scratch);
          }
          continue;
        }

        int x1 = static_cast<int>(floorf(start));
        int x2 = static_cast<int>(ceilf(end));
        if (x1 == x2)
          x2++;

        // FAST MODULO LOGIC
        int wx = x1 % W;
        if (wx < 0)
          wx += W;

        for (int x = x1; x < x2; ++x) {
          Vector p = pixel_to_vector<W, H>(wx, y);
          process_pixel<W, H, ComputeUVs>(wx, y, p, pipeline, canvas, shape,
                                          fragment_shader, effective_debug,
                                          result_scratch, frag_scratch);

          // Modulo-free coordinate wrap!
          wx++;
          if (wx >= W)
            wx -= W;
        }
      }
      intervals.clear();
    } else {
      if (!handled) {
        for (int x = 0; x < W; ++x) {
          Vector p = pixel_to_vector<W, H>(x, y);
          process_pixel<W, H, ComputeUVs>(x, y, p, pipeline, canvas, shape,
                                          fragment_shader, effective_debug,
                                          result_scratch, frag_scratch);
        }
      }
    }
  }
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
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, float thickness,
                   FragmentShaderFn fragment_shader, float phase = 0,
                   bool debug_bb = false) {
    auto res = get_antipode(basis, radius);
    SDF::Ring shape(res.first, res.second, thickness, phase);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }

  // Overload for Vector normal inputs
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Vector &normal,
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
    SDF::Star<W> shape(res.first, res.second, sides, phase, H, H);
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
  template <int W, int H, typename MeshT>
  static void draw(PipelineRef pipeline, Canvas &canvas, const MeshT &mesh,
                   FragmentShaderFn fragment_shader, Arena &scratch_arena,
                   bool debug_bb = false) {
    ScratchScope scope(scratch_arena);
    auto *scratch =
        static_cast<SDF::FaceScratchBuffer *>(scratch_arena.allocate(
            sizeof(SDF::FaceScratchBuffer), alignof(SDF::FaceScratchBuffer)));

    size_t idx_offset = 0;
    for (size_t i = 0; i < mesh.num_faces; ++i) {
      size_t count = mesh.face_counts[i];
      std::span<const Vector> verts(mesh.vertices.data(), mesh.num_vertices);
      std::span<const uint16_t> indices(&mesh.faces[idx_offset], count);

      SDF::Face shape(verts, indices, 0.0f, *scratch, H + hs::H_OFFSET, H);
      idx_offset += count;

      auto wrapper = [&](const Vector &p, Fragment &f_in) {
        f_in.v2 = static_cast<float>(i);
        fragment_shader(p, f_in);
      };
      Scan::rasterize<W, H, true>(pipeline, canvas, shape, wrapper, debug_bb);
    }
  }

  // Overload for MeshState
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const MeshState &mesh,
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

      SDF::Face shape(verts, indices, 0.0f, *scratch, H + hs::H_OFFSET, H);

      auto wrapper = [&](const Vector &p, Fragment &f_in) {
        f_in.v2 = static_cast<float>(i);
        fragment_shader(p, f_in);
      };

      Scan::rasterize<W, H, true>(pipeline, canvas, shape, wrapper, debug_bb);
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
  template <int W, int H, int SAMPLES = 4, typename ShaderFn>
  static void draw(Canvas &canvas, ShaderFn &&shader) {
    constexpr float h_virt_minus_1 = static_cast<float>(H + hs::H_OFFSET - 1);
    constexpr float w_float = static_cast<float>(W);

    if constexpr (SAMPLES == 1) {
      // 1× — single sample at pixel center (no SSAA)
      for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
          float px = static_cast<float>(x);
          float py = static_cast<float>(y);
          float theta = (px * 2.0f * PI_F) / w_float;
          float phi = (py * PI_F) / h_virt_minus_1;
          float sin_phi = sinf(phi);
          Vector v(sin_phi * cosf(theta), cosf(phi), sin_phi * sinf(theta));
          Color4 sample = shader(v);
          canvas(x, y) = sample.color * sample.alpha;
        }
      }
    } else {
      // SAMPLES× SSAA — sub-pixel samples averaged
      constexpr float inv_samples = 1.0f / SAMPLES;
      constexpr float eps = 0.5f;
      constexpr float offsets_x[4] = {eps, -eps, eps, -eps};
      constexpr float offsets_y[4] = {eps, eps, -eps, -eps};

      for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
          Color4 accum(Pixel(0, 0, 0), 0.0f);

          for (int i = 0; i < SAMPLES; ++i) {
            float px = static_cast<float>(x) + offsets_x[i];
            float py = static_cast<float>(y) + offsets_y[i];

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
      for (int y = 0; y < H; ++y) {
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
      constexpr float eps = 0.5f;
      constexpr float offsets_x[4] = {eps, -eps, eps, -eps};
      constexpr float offsets_y[4] = {eps, eps, -eps, -eps};

      for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
          Vector center_v = pixel_to_vector<W, H>(x, y);

          frag_base.pos = center_v;
          vertex_shader(frag_base);

          Color4 accum(Pixel(0, 0, 0), 0.0f);

          for (int i = 0; i < SAMPLES; ++i) {
            float px = static_cast<float>(x) + offsets_x[i];
            float py = static_cast<float>(y) + offsets_y[i];
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
 * @brief Raymarch volume renderer with orthographic projection.
 *
 * The render loop is internal: callers provide an SDF distance function
 * (evaluated per march step) and a fragment shader (evaluated once per hit
 * to populate Fragment registers for shading).
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
  template <int W, int H>
  static void
  draw(PipelineRef pipeline, Canvas &canvas, const Vector &bounds_center,
       float bounds_radius, const Vector &view_dir, SDFDistanceFn sdf_fn,
       SDFFragmentFn frag_fn, int max_steps = 15, float aa_width = 0.01f) {
    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();

    // Normalized view direction
    float vd_len = sqrtf(view_dir.x * view_dir.x + view_dir.y * view_dir.y +
                         view_dir.z * view_dir.z);
    float vd_inv = (vd_len > TOLERANCE) ? 1.0f / vd_len : 1.0f;
    Vector vd(view_dir.x * vd_inv, view_dir.y * vd_inv, view_dir.z * vd_inv);

    // Ray must start behind the farthest extent of the torus:
    // vertex is at distance 1.0 (on unit sphere), torus extends bounds_radius
    float start_offset = 1.0f + bounds_radius;

    // Precompute bounds_center projected onto the view plane (⊥ vd)
    float bc_dot_vd = bounds_center.x * vd.x + bounds_center.y * vd.y +
                      bounds_center.z * vd.z;
    Vector bc_proj(bounds_center.x - bc_dot_vd * vd.x,
                   bounds_center.y - bc_dot_vd * vd.y,
                   bounds_center.z - bc_dot_vd * vd.z);
    float bounds_r2 = bounds_radius * bounds_radius;

    for (int y = 0; y < H; ++y) {
      for (int wx = 0; wx < W; ++wx) {
        Vector p = pixel_to_vector<W, H>(wx, y);

        // Back-face cull
        float facing = p.x * vd.x + p.y * vd.y + p.z * vd.z;
        if (facing >= 0.0f)
          continue;

        // Orthographic ray-sphere cull
        float pp_x = p.x - facing * vd.x;
        float pp_y = p.y - facing * vd.y;
        float pp_z = p.z - facing * vd.z;
        float dx = pp_x - bc_proj.x;
        float dy = pp_y - bc_proj.y;
        float dz = pp_z - bc_proj.z;
        if (dx * dx + dy * dy + dz * dz > bounds_r2)
          continue;

        // Orthographic ray origin: outside the unit sphere
        Vector ro(pp_x - vd.x * start_offset, pp_y - vd.y * start_offset,
                  pp_z - vd.z * start_offset);

        // --- Sphere tracing ---
        Vector march_p = ro;
        float closest_d = 999.0f;
        Vector closest_p = ro;

        for (int i = 0; i < max_steps; ++i) {
          float d = sdf_fn(march_p);

          if (d < closest_d) {
            closest_d = d;
            closest_p = march_p;
          }

          // Stop if deep inside surface
          if (d < -aa_width)
            break;

          // SDF-guided step with very small floor to avoid quantization
          float advance = std::max(d * 0.9f, 1e-5f);
          march_p =
              Vector(march_p.x + vd.x * advance, march_p.y + vd.y * advance,
                     march_p.z + vd.z * advance);
        }

        if (closest_d >= aa_width)
          continue;

        // --- Fragment shading ---
        Fragment frag;
        frag.pos = closest_p;
        frag.size = closest_d;
        frag_fn(closest_p, frag);

        // One-sided AA with quintic kernel (matches other Scan shapes)
        float hit_threshold = aa_width * 0.1f;
        float edge_alpha;
        if (closest_d <= hit_threshold) {
          edge_alpha = 1.0f;
        } else {
          edge_alpha = quintic_kernel(
              1.0f - (closest_d - hit_threshold) / (aa_width - hit_threshold));
        }

        if (frag.color.alpha * edge_alpha > 0.001f) {
          pipeline.plot(canvas, p, frag.color.color, 0.0f,
                        frag.color.alpha * edge_alpha);
        }
      }
    }
  }
};

} // namespace Scan
