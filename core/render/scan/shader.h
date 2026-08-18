/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <type_traits>
#include <utility>
#include "render/sdf.h"
#include "render/shading.h"
#include "mesh/mesh_class_types.h"
#include "color/color.h"
#include "render/filter.h"
#include "engine/static_circular_buffer.h"
#include "render/canvas.h"
#include "engine/platform.h"
#include "render/scan/raster.h"

/**
 * @file shader.h
 * @brief Scan::Shader: full-screen per-pixel shaders with SSAA.
 */

namespace Scan {

/**
 * @brief Full-screen per-pixel shaders with SAMPLES× SSAA.
 *
 * Three entry points, in increasing order of caller control:
 * - draw(canvas, shader): one callable ShaderFn(const Vector &v) -> Color4,
 *   invoked SAMPLES× per pixel at sub-pixel offsets and averaged.
 * - draw(canvas, fragment_shader, vertex_shader): splits per-pixel setup
 *   (vertex_shader, once at the pixel center) from per-sub-sample evaluation.
 *   Both callables are required; a null one traps.
 * - draw_grid(canvas, vertex_shader, pixel_shader): hands the seeded fragment
 *   and the row's SsaaGrid to pixel_shader, which owns the sampling and returns
 *   the finished pixel.
 *
 * @details Every entry point assigns the finished premultiplied color to the
 * canvas rather than plotting it, so the destination is overwritten: alpha < 1
 * darkens the pixel instead of blending with what is under it, and no
 * plot-time filter stage (Filter::World / Filter::Pixel) sees it. These entry
 * points take no pipeline; an effect needing the filter chain must plot
 * through it itself.
 */
struct Shader {
  // --- Shared SSAA helpers (used by every entry point) -----------------------
  /**
   * @brief Per-draw sub-pixel trig for the 2×2 SSAA sample grid, derived from
   *        the resident engine trig LUT.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @details Each sub-sample sits a constant ±0.25 px from an integer pixel,
   * i.e. a constant angular offset from that pixel's theta/phi under the same
   * parameterization as pixel_to_vector (theta = 2*pi*px/W,
   * phi = py*pi/(H_VIRT-1)). So sin/cos at a sub-sample follow from the
   * integer-pixel TrigLUT<W,H> tables by one angle-addition, keeping the
   * samples LUT-consistent with the non-SSAA path and needing no per-column
   * storage. The two rotation angles (d_theta, d_phi) are resolution constants;
   * their sin/cos are built once per draw (libm is not constexpr) into a
   * handful of floats.
   */
  template <int W, int H> struct SsaaGrid {
    /** @brief Sub-pixel samples per pixel supplied by the 2×2 grid. */
    static constexpr int SAMPLES = 4;

    float sin_phi[2]; /**< Current row's sin(phi) at y+0.25 [0] / y-0.25 [1]. */
    float cos_phi[2]; /**< Current row's cos(phi) at y+0.25 [0] / y-0.25 [1]. */
    float cos_dtheta; /**< cos of the ±0.25 px column rotation. */
    float sin_dtheta; /**< sin of the ±0.25 px column rotation. */
    float cos_dphi;   /**< cos of the ±0.25 px row rotation. */
    float sin_dphi;   /**< sin of the ±0.25 px row rotation. */

    SsaaGrid() {
      // d_theta = 0.25 px * (2*pi/W); d_phi = 0.25 px * (pi/(H_VIRT-1)).
      constexpr float d_theta = 0.5f * PI_F / static_cast<float>(W);
      constexpr float h_virt_minus_1 = static_cast<float>(H + hs::H_OFFSET - 1);
      constexpr float d_phi = 0.25f * PI_F / h_virt_minus_1;
      cos_dtheta = cosf(d_theta);
      sin_dtheta = sinf(d_theta);
      cos_dphi = cosf(d_phi);
      sin_dphi = sinf(d_phi);
    }

    /** @brief Loads the two phi trig pairs for pixel row y from the LUT. */
    void set_row(int y) {
      const float sy = TrigLUT<W, H>::sin_phi[y];
      const float cy = TrigLUT<W, H>::cos_phi[y];
      // Row 0 = y+0.25, row 1 = y-0.25 (the 2×2 grid's centered ±0.25 offsets).
      sin_phi[0] = sy * cos_dphi + cy * sin_dphi;
      cos_phi[0] = cy * cos_dphi - sy * sin_dphi;
      sin_phi[1] = sy * cos_dphi - cy * sin_dphi;
      cos_phi[1] = cy * cos_dphi + sy * sin_dphi;
    }

    /**
     * @brief World-space unit vector for sample i of pixel x in the current
     *        row (see set_row).
     * @param x Pixel column.
     * @param i Sample index in [0, SAMPLES); the low bit selects the column
     * offset (±0.25 px) and bit 1 the row offset (±0.25 px).
     */
    Vector at(int x, int i) const {
      const float st = TrigLUT<W, H>::sin_theta[x];
      const float ct = TrigLUT<W, H>::cos_theta(x);
      // Column 0 (i&1==0) = x+0.25, column 1 = x-0.25.
      const float s = (i & 1) ? -sin_dtheta : sin_dtheta;
      const float sin_theta = st * cos_dtheta + ct * s;
      const float cos_theta = ct * cos_dtheta - st * s;
      const float sp = sin_phi[(i >> 1) & 1];
      return Vector(sp * cos_theta, cos_phi[(i >> 1) & 1], sp * sin_theta);
    }
  };

  /**
   * @brief Validates the LUT-domain invariant shared by every entry point.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @param cr Clip region whose bounds are checked against the LUT extents.
   * @details Checked once per draw, not per pixel: every (x,y) the loops feed to
   * pixel_to_vector indexes the trig LUTs within bounds. The row bound is H, not
   * the LUT's H_VIRT (H plus the pole offset): those rows also subscript the
   * canvas, which holds H of them, so the tighter bound is the binding one.
   */
  template <int W, int H> static void check_lut_domain(const ClipRegion &cr) {
    HS_CHECK(cr.x_start >= 0 && cr.x_end <= W && cr.render_y_start() >= 0 &&
             cr.render_y_end() <= H);
  }
  // --------------------------------------------------------------------------

  /**
   * @brief Full-screen per-pixel shader with SAMPLES× SSAA from a single
   *        callable.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam SAMPLES Number of sub-pixel samples per pixel (1 disables SSAA).
   * @tparam ShaderFn Callable ShaderFn(const Vector &v) -> Color4.
   * @param canvas Destination canvas.
   * @param shader Maps a world-space unit vector to a final color; invoked
   *               SAMPLES× per pixel at sub-pixel offsets and averaged.
   */
private:
  template <int W, int H, int SAMPLES, typename ShaderFn>
  HS_O3_FN __attribute__((always_inline)) static void
  draw_typed(Canvas &canvas, ShaderFn &&shader) {
    // The sample-offset table has four distinct sub-pixel positions; only 1 and
    // the 2x2 grid (4) are valid.
    static_assert(SAMPLES == 1 || SAMPLES == 4,
                  "Scan::Shader SSAA supports only SAMPLES == 1 or 4");
    check_canvas_dims<W, H>(canvas);
    const auto &cr = canvas.clip();
    check_lut_domain<W, H>(cr);
    const auto xc = cr.x_clip();

    if constexpr (SAMPLES == 1) {
      for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y) {
        walk_clip_columns<W>(xc, [&](int x) {
          Vector v = pixel_to_vector<W, H>(x, y);
          Color4 sample = shader(v);
          canvas(x, y) = sample.color * sample.alpha;
        });
      }
    } else {
      constexpr float inv_samples = 1.0f / SAMPLES;
      if (!TrigLUT<W, H>::initialized)
        TrigLUT<W, H>::init();
      SsaaGrid<W, H> grid;

      for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y) {
        grid.set_row(y);
        walk_clip_columns<W>(xc, [&](int x) {
          // Premultiplied SSAA: accumulate each sample's coverage-weighted color
          // (color * alpha / N), matching the SAMPLES==1 path.
          Pixel accum(0, 0, 0);

          for (int i = 0; i < SAMPLES; ++i) {
            Color4 sample = shader(grid.at(x, i));
            accum += sample.color * (sample.alpha * inv_samples);
          }

          canvas(x, y) = accum;
        });
      }
    }
  }

public:
  template <int W, int H, int SAMPLES = 1, typename ShaderFn>
  HS_O3_FN static void draw(Canvas &canvas, ShaderFn &&shader) {
    draw_typed<W, H, SAMPLES>(canvas, static_cast<ShaderFn &&>(shader));
  }

  /**
   * @brief Direct typed shader draw whose traversal executes from cached flash.
   * @details The callable remains statically bound and is inlined into this
   * instantiation; only its code placement differs from draw().
   */
  template <int W, int H, int SAMPLES = 1, typename ShaderFn>
  HS_HOT_FLASH_MEMBER static void draw_cached(Canvas &canvas,
                                              ShaderFn &&shader) {
    draw_typed<W, H, SAMPLES>(canvas, static_cast<ShaderFn &&>(shader));
  }

  /**
   * @brief Full-screen per-pixel shader with SAMPLES× SSAA and a split
   *        vertex/fragment shader.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam SAMPLES Number of sub-pixel samples per pixel (1 disables SSAA).
   * @param canvas Destination canvas.
   * @param fragment_shader Per-sub-sample shader, called SAMPLES× per pixel.
   * @param vertex_shader Per-pixel shader, called once at the pixel center.
   * @details Splits expensive per-pixel work (vertex_shader, once at pixel
   * center) from per-sub-sample evaluation (fragment_shader, SAMPLES×).
   *
   * @note SAMPLES defaults to 1 (no SSAA), matching the single-callback overload.
   */
  template <int W, int H, int SAMPLES = 1>
  static void draw(Canvas &canvas, FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader) {
    // Only 1 and 4 are supported (see the single-callback overload).
    static_assert(SAMPLES == 1 || SAMPLES == 4,
                  "Scan::Shader SSAA supports only SAMPLES == 1 or 4");
    // Cold (once per draw), not per-pixel: trap null shaders here so they fail
    // deterministically instead of calling a null thunk under NDEBUG.
    HS_CHECK(vertex_shader,
             "Scan::Shader::draw requires a non-null vertex_shader");
    HS_CHECK(fragment_shader,
             "Scan::Shader::draw requires a non-null fragment_shader");
    check_canvas_dims<W, H>(canvas);
    // frag_base is per pixel, not per draw: each pixel starts from a default
    // Fragment, so a vertex shader writing only some registers (v0-v3/size/age/
    // color) can't inherit the previous pixel's values.
    if constexpr (SAMPLES == 1) {
      const auto &cr = canvas.clip();
      check_lut_domain<W, H>(cr);
      const auto xc = cr.x_clip();
      for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y) {
        walk_clip_columns<W>(xc, [&](int x) {
          Vector center_v = pixel_to_vector<W, H>(x, y);
          Fragment frag_base;
          frag_base.pos = center_v;
          vertex_shader(frag_base);
          fragment_shader(center_v, frag_base);
          // Premultiply by alpha, matching the single-callback overload and the
          // process_pixel/Volume contract.
          canvas(x, y) = frag_base.color.color * frag_base.color.alpha;
        });
      }
    } else {
      constexpr float inv_samples = 1.0f / SAMPLES;
      if (!TrigLUT<W, H>::initialized)
        TrigLUT<W, H>::init();
      SsaaGrid<W, H> grid;

      const auto &cr = canvas.clip();
      check_lut_domain<W, H>(cr);
      const auto xc = cr.x_clip();
      for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y) {
        grid.set_row(y);
        walk_clip_columns<W>(xc, [&](int x) {
          Vector center_v = pixel_to_vector<W, H>(x, y);

          Fragment frag_base;
          frag_base.pos = center_v;
          vertex_shader(frag_base);

          // Premultiplied SSAA: accumulate coverage-weighted color directly (see
          // the single-callback overload).
          Pixel accum(0, 0, 0);

          for (int i = 0; i < SAMPLES; ++i) {
            Vector v = grid.at(x, i);

            Fragment sub_frag = frag_base;
            sub_frag.pos = v;

            fragment_shader(v, sub_frag);

            Color4 sample = sub_frag.color;
            accum += sample.color * (sample.alpha * inv_samples);
          }

          canvas(x, y) = accum;
        });
      }
    }
  }

  /**
   * @brief Full-screen draw handing the effect the whole per-pixel body.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam VertexFn Callable VertexFn(Fragment&) — per-pixel seed.
   * @tparam PixelFn Callable PixelFn(Fragment&, const SsaaGrid<W,H>&, int x)
   *         -> Pixel — computes the finished (already premultiplied) pixel.
   * @param canvas Destination canvas.
   * @param vertex_shader Per-pixel shader, called once at the pixel center.
   * @param pixel_shader Owns the pixel: receives the seeded fragment, the
   *        sub-pixel SSAA grid for the current row, and the pixel column, and
   *        returns the final (premultiplied) pixel. A template, not a
   *        type-erased FunctionRef: the whole body inlines into this loop, and
   *        the effect can hoist work its sub-samples share.
   * @details Same outer scaffolding as the SSAA draw() overloads (clip,
   * LUT-domain check, trig-LUT init, per-row SsaaGrid); the
   * per-pixel work is delegated whole so the caller controls the sampling.
   */
  template <int W, int H, typename VertexFn, typename PixelFn>
  HS_O3_FN static void draw_grid(Canvas &canvas, VertexFn &&vertex_shader,
                                 PixelFn &&pixel_shader) {
    check_canvas_dims<W, H>(canvas);
    const auto &cr = canvas.clip();
    check_lut_domain<W, H>(cr);
    const auto xc = cr.x_clip();
    if (!TrigLUT<W, H>::initialized)
      TrigLUT<W, H>::init();
    SsaaGrid<W, H> grid;
    for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y) {
      grid.set_row(y);
      walk_clip_columns<W>(xc, [&](int x) {
        Fragment frag_base;
        frag_base.pos = pixel_to_vector<W, H>(x, y);
        vertex_shader(frag_base);
        canvas(x, y) = pixel_shader(frag_base, grid, x);
      });
    }
  }
};

} // namespace Scan
