/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include <type_traits>
#include "engine/memory.h"
#include "render/shading.h"
#include "color/color.h"
#include "render/canvas.h"
#include "platform/platform.h"
#include "render/scan/raster.h"

/**
 * @file shader.h
 * @brief Scan::Shader: full-screen per-pixel shaders with SSAA.
 */

namespace Scan {

/**
 * @brief Full-screen per-pixel shaders with SAMPLES× SSAA.
 *
 * Entry points:
 * - draw(canvas, shader): one callable ShaderFn(const Vector &v) -> Color4
 *   or premultiplied Pixel, invoked SAMPLES× per pixel at sub-pixel offsets
 *   and averaged.
 * - draw_cached(canvas, shader): the same typed draw with its traversal placed
 *   in cached flash.
 * - draw(canvas, fragment_shader, vertex_shader): splits per-pixel setup
 *   (vertex_shader, once at the pixel center) from per-sub-sample evaluation.
 *   Both callables are required; a null one traps.
 * - draw_block_coherent(canvas, block, positions, scratch, classify, shade):
 *   shades from the union of a canvas-anchored block's corner candidates.
 * - draw_grid(canvas, vertex_shader, pixel_shader): hands the seeded fragment
 *   and the row's SsaaGrid to pixel_shader, which owns the sampling and returns
 *   the finished pixel.
 *
 * @details Every entry point assigns the finished premultiplied color to the
 * canvas rather than plotting it, so the destination is overwritten: alpha < 1
 * darkens the pixel instead of blending with what is under it, and no
 * plot-time filter stage (Filter::World / Filter::Screen / Filter::Pixel) sees it. These entry
 * points take no pipeline; an effect needing the filter chain must plot
 * through it itself. They take no debug flag and do not read canvas.debug():
 * every pixel is covered, so there is no scan bound for the bounding-box tint
 * to mark.
 */
struct Shader {
  // --- Shared SSAA helpers (used by every entry point) -----------------------
  /**
   * @brief Per-draw sub-pixel trig for the 2×2 SSAA sample grid, derived from
   *        the resident engine trig LUT.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @details Samples sit at constant quarter-row and quarter-column angular
   * offsets from the calibrated integer-pixel lookup tables.
   */
  template <int W, int H> struct SsaaGrid {
    static constexpr int WIDTH = W;  /**< Column lookup-table extent. */
    static constexpr int HEIGHT = H; /**< Row lookup-table extent. */
    /** @brief Sub-pixel samples per pixel supplied by the 2×2 grid. */
    static constexpr int SAMPLES = 4;

    float sin_phi[2]; /**< Current row's sin(phi) at y+0.25 [0] / y-0.25 [1]. */
    float cos_phi[2]; /**< Current row's cos(phi) at y+0.25 [0] / y-0.25 [1]. */
    float cos_dtheta; /**< cos of the ±0.25 px column rotation. */
    float sin_dtheta; /**< sin of the ±0.25 px column rotation. */
    float cos_dphi;   /**< cos of the ±0.25 px row rotation. */
    float sin_dphi;   /**< sin of the ±0.25 px row rotation. */

    SsaaGrid() {
      constexpr float d_theta = 0.5f * math::PI_F / static_cast<float>(W);
      const float d_phi = 0.25f * math::RADIANS_PER_ROW<H>;
      cos_dtheta = cosf(d_theta);
      sin_dtheta = sinf(d_theta);
      cos_dphi = cosf(d_phi);
      sin_dphi = sinf(d_phi);
    }

    /** @brief Loads the two phi trig pairs for pixel row y from the LUT. */
    void set_row(int y) {
      const float sy = math::TrigLUT<W, H>::sin_phi[y];
      const float cy = math::TrigLUT<W, H>::cos_phi[y];
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
    math::Vector at(int x, int i) const {
      const float st = math::TrigLUT<W, H>::sin_theta[x];
      const float ct = math::TrigLUT<W, H>::cos_theta(x);
      // Column 0 (i&1==0) = x+0.25, column 1 = x-0.25.
      const float s = (i & 1) ? -sin_dtheta : sin_dtheta;
      const float sin_theta = st * cos_dtheta + ct * s;
      const float cos_theta = ct * cos_dtheta - st * s;
      const float sp = sin_phi[(i >> 1) & 1];
      return math::Vector(sp * cos_theta, cos_phi[(i >> 1) & 1],
                          sp * sin_theta);
    }
  };

  /**
   * @brief Validates the LUT-domain invariant shared by every entry point.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @param cr Clip region whose bounds are checked against the LUT extents.
   * @details Clip bounds must fit both the canvas and row lookup tables.
   */
  template <int W, int H> static void check_lut_domain(const ClipRegion &cr) {
    HS_CHECK(cr.x_start >= 0 && cr.x_end <= W && cr.render_y_start() >= 0 &&
                 cr.render_y_end() <= H,
             "scan clip region outside the canvas and trig LUT domain");
  }
  // --------------------------------------------------------------------------

private:
  // Shared traversal for typed shader entry points.
  template <int W, int H, int SAMPLES, typename ShaderFn,
            typename RowFn = std::nullptr_t>
  HS_O3_FN __attribute__((always_inline)) static void
  draw_typed(Canvas &canvas, ShaderFn &&shader, RowFn &&begin_row = nullptr) {
    // SsaaGrid::at() provides the four +/-0.25-pixel positions.
    static_assert(SAMPLES == 1 || SAMPLES == 4,
                  "Scan::Shader SSAA supports only SAMPLES == 1 or 4");
    auto sample = [&](const math::Vector &v, int x, int y) {
      if constexpr (std::is_invocable_v<ShaderFn &, const math::Vector &, int,
                                        int>)
        return shader(v, x, y);
      else
        return shader(v);
    };
    constexpr bool PREMULTIPLIED =
        std::is_same_v<std::decay_t<decltype(sample(
                           std::declval<const math::Vector &>(), 0, 0))>,
                       Pixel>;
    check_canvas_dims<W, H>(canvas);
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    const auto &cr = canvas.clip();
    check_lut_domain<W, H>(cr);
    const auto xc = cr.x_clip();

    if constexpr (SAMPLES == 1) {
      for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y) {
        if constexpr (!std::is_same_v<RowFn, std::nullptr_t>)
          begin_row(y);
        const float sp = math::TrigLUT<W, H>::sin_phi[y];
        const float cp = math::TrigLUT<W, H>::cos_phi[y];
        walk_clip_columns<W>(xc, [&](int x) {
          math::Vector v =
              math::Vector(sp * math::TrigLUT<W, H>::cos_theta(x), cp,
                           sp * math::TrigLUT<W, H>::sin_theta[x]);
          if constexpr (PREMULTIPLIED)
            canvas(x, y) = sample(v, x, y);
          else {
            Color4 color = sample(v, x, y);
            canvas(x, y) = color.color * color.alpha;
          }
        });
      }
    } else {
      constexpr float inv_samples = 1.0f / SAMPLES;
      SsaaGrid<W, H> grid;

      for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y) {
        if constexpr (!std::is_same_v<RowFn, std::nullptr_t>)
          begin_row(y);
        grid.set_row(y);
        walk_clip_columns<W>(xc, [&](int x) {
          // Premultiplied SSAA: accumulate each sample's coverage-weighted color
          // (color * alpha / N), matching the SAMPLES==1 path.
          Pixel accum(0, 0, 0);

          for (int i = 0; i < SAMPLES; ++i) {
            if constexpr (PREMULTIPLIED)
              accum += sample(grid.at(x, i), x, y) * inv_samples;
            else {
              Color4 color = sample(grid.at(x, i), x, y);
              accum += color.color * (color.alpha * inv_samples);
            }
          }

          canvas(x, y) = accum;
        });
      }
    }
  }

public:
  /**
   * @brief Full-screen per-pixel shader with SAMPLES× SSAA from a single
   *        callable.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam SAMPLES Number of sub-pixel samples per pixel (1 disables SSAA).
   * @tparam ShaderFn Callable shader(const Vector &v), optionally with integer
   *         pixel coordinates as shader(v, x, y); returns Color4 or a Pixel
   *         already premultiplied by its alpha.
   * @param canvas Destination canvas.
   * @param shader Maps a world-space unit vector to a final color; invoked
   *               SAMPLES× per pixel at sub-pixel offsets and averaged.
   */
  template <int W, int H, int SAMPLES = 1, typename ShaderFn>
  HS_O3_FN static void draw(Canvas &canvas, ShaderFn &&shader) {
    draw_typed<W, H, SAMPLES>(canvas, static_cast<ShaderFn &&>(shader));
  }

  /**
   * @brief Direct typed shader draw whose traversal executes from cached flash.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam SAMPLES Number of sub-pixel samples per pixel (1 disables SSAA).
   * @tparam ShaderFn Callable shader(const Vector &v), optionally with integer
   *         pixel coordinates as shader(v, x, y); returns Color4 or a Pixel
   *         already premultiplied by its alpha.
   * @param canvas Destination canvas.
   * @param shader Maps a world-space unit vector to a final color; invoked
   *               SAMPLES× per pixel at sub-pixel offsets and averaged.
   *               May also accept integer pixel coordinates after the vector.
   * @param begin_row Optional callable (int y), invoked once before each row.
   * @details The callable remains statically bound and is inlined into this
   * instantiation; only its code placement differs from draw().
   */
  template <int W, int H, int SAMPLES = 1, typename ShaderFn,
            typename RowFn = std::nullptr_t>
  HS_HOT_FLASH_MEMBER static void draw_cached(Canvas &canvas, ShaderFn &&shader,
                                              RowFn &&begin_row = nullptr) {
    draw_typed<W, H, SAMPLES>(canvas, static_cast<ShaderFn &&>(shader),
                              static_cast<RowFn &&>(begin_row));
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
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    // frag_base is per pixel, not per draw: each pixel starts from a default
    // Fragment, so a vertex shader writing only some registers (v0-v3/size/age/
    // color) can't inherit the previous pixel's values.
    if constexpr (SAMPLES == 1) {
      const auto &cr = canvas.clip();
      check_lut_domain<W, H>(cr);
      const auto xc = cr.x_clip();
      for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y) {
        const float sp = math::TrigLUT<W, H>::sin_phi[y];
        const float cp = math::TrigLUT<W, H>::cos_phi[y];
        walk_clip_columns<W>(xc, [&](int x) {
          math::Vector center_v =
              math::Vector(sp * math::TrigLUT<W, H>::cos_theta(x), cp,
                           sp * math::TrigLUT<W, H>::sin_theta[x]);
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
      SsaaGrid<W, H> grid;

      const auto &cr = canvas.clip();
      check_lut_domain<W, H>(cr);
      const auto xc = cr.x_clip();
      for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y) {
        const float sp = math::TrigLUT<W, H>::sin_phi[y];
        const float cp = math::TrigLUT<W, H>::cos_phi[y];
        grid.set_row(y);
        walk_clip_columns<W>(xc, [&](int x) {
          math::Vector center_v =
              math::Vector(sp * math::TrigLUT<W, H>::cos_theta(x), cp,
                           sp * math::TrigLUT<W, H>::sin_theta[x]);

          Fragment frag_base;
          frag_base.pos = center_v;
          vertex_shader(frag_base);

          // Premultiplied SSAA: accumulate coverage-weighted color directly (see
          // the single-callback overload).
          Pixel accum(0, 0, 0);

          for (int i = 0; i < SAMPLES; ++i) {
            math::Vector v = grid.at(x, i);

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
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    const auto &cr = canvas.clip();
    check_lut_domain<W, H>(cr);
    const auto xc = cr.x_clip();
    SsaaGrid<W, H> grid;
    for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y) {
      const float sp = math::TrigLUT<W, H>::sin_phi[y];
      const float cp = math::TrigLUT<W, H>::cos_phi[y];
      grid.set_row(y);
      walk_clip_columns<W>(xc, [&](int x) {
        Fragment frag_base;
        frag_base.pos = math::Vector(sp * math::TrigLUT<W, H>::cos_theta(x), cp,
                                     sp * math::TrigLUT<W, H>::sin_theta[x]);
        vertex_shader(frag_base);
        canvas(x, y) = pixel_shader(frag_base, grid, x);
      });
    }
  }
  /** @brief K site indices classified at one block-grid corner. */
  template <size_t K> using BlockCell = std::array<uint16_t, K>;

  /** @brief Deduped union of four corners, with contiguous shading positions. */
  template <size_t K> struct BlockCandidates {
    static_assert(K > 0 && 4 * K <= UINT8_MAX,
                  "block candidates require 1..63 sites per corner");
    math::Vector pos[4 * K];
    uint16_t idx[4 * K];
    uint8_t n;
  };

  /** @brief Maximum corner-grid and candidate-row bytes at a minimum block size.
   * @details Excludes caller-owned site positions and any arena alignment pad.
   */
  template <int W, int H, size_t K, int MIN_BLOCK>
  static constexpr size_t block_coherent_scratch_bytes() {
    static_assert(W > 0 && H > 0 && MIN_BLOCK > 0);
    constexpr size_t COLS = (W - 1) / MIN_BLOCK + 2;
    constexpr size_t ROWS = (H - 1) / MIN_BLOCK + 2;
    return COLS * ROWS * sizeof(BlockCell<K>) +
           (COLS - 1) * sizeof(BlockCandidates<K>);
  }

  /**
   * @brief Shade pixels from their block's deduped corner candidates.
   * @param canvas Destination canvas.
   * @param block Positive block edge in pixels; corners anchor to the canvas.
   * @param positions Site positions indexed by every returned corner index.
   * @param scratch Arena for the corner grid and one candidate row. The caller
   * owns its scope; block_coherent_scratch_bytes() bounds these allocations.
   * @param classify_corner Callable (const Vector&) returning BlockCell<K>.
   * @param shade_candidates Callable (const Vector&, const BlockCandidates<K>&)
   * returning Color4 or a premultiplied Pixel.
   * @details Sites absent from all four corners are omitted. Wrapped column
   * bands skip their interior gap; candidate order follows corner order.
   */
  template <int W, int H, size_t K, typename ClassifyFn, typename ShadeFn>
  HS_O3_FN __attribute__((always_inline)) static void
  draw_block_coherent(Canvas &canvas, int block, const math::Vector *positions,
                      Arena &scratch, ClassifyFn &&classify_corner,
                      ShadeFn &&shade_candidates) {
    HS_CHECK(block > 0, "block-coherent shader requires a positive block size");
    check_canvas_dims<W, H>(canvas);
    const auto &cr = canvas.clip();
    check_lut_domain<W, H>(cr);
    const auto columns = cr.x_clip();
    const int x0 = columns.active && !columns.wrap ? columns.rs : 0;
    const int x1 = columns.active && !columns.wrap ? columns.re : W;
    const int y0 = cr.render_y_start();
    const int y1 = cr.render_y_end();
    if (x1 <= x0 || y1 <= y0)
      return;
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    const int B = block;
    // Canvas-anchored corners give each pixel the same candidates in every band.
    const int gx0 = (x0 / B) * B;
    const int gy0 = (y0 / B) * B;
    const int nbx = (x1 - 1 - gx0) / B + 2; // corner columns spanning [x0, x1)
    const int nby = (y1 - 1 - gy0) / B + 2; // corner rows spanning    [y0, y1)
    const int GAP_BEGIN =
        columns.active && columns.wrap ? (columns.re + B - 1) / B : nbx;
    const int GAP_END = columns.active && columns.wrap ? columns.rs / B : nbx;
    BlockCell<K> *cells = scratch.allocate_n<BlockCell<K>>(nbx * nby);
    // End corners clamp to the last canvas pixel, independent of the clip.
    auto corner_x = [&](int j) { return std::min(gx0 + j * B, W - 1); };
    auto corner_y = [&](int k) { return std::min(gy0 + k * B, H - 1); };

    for (int k = 0; k < nby; ++k)
      for (int j = 0; j < nbx; ++j) {
        if (j > GAP_BEGIN && j < GAP_END)
          continue;
        cells[k * nbx + j] = classify_corner(
            math::pixel_to_vector<W, H>(corner_x(j), corner_y(k)));
      }

    // One candidate set per block column, rebuilt on each block-row change.
    // Positions are copied in so the per-pixel scan runs over contiguous data.
    const int nblk = nbx - 1;
    BlockCandidates<K> *cands = scratch.allocate_n<BlockCandidates<K>>(nblk);
    auto build_candidate_row = [&](int ky) {
      for (int jx = 0; jx < nblk; ++jx) {
        if (jx >= GAP_BEGIN && jx < GAP_END)
          continue;
        BlockCandidates<K> &cs = cands[jx];
        cs.n = 0;
        auto add = [&](uint16_t s) {
          for (uint8_t i = 0; i < cs.n; ++i)
            if (cs.idx[i] == s)
              return;
          cs.idx[cs.n] = s;
          cs.pos[cs.n] = positions[s];
          ++cs.n;
        };
        for (const BlockCell<K> *c :
             {&cells[ky * nbx + jx], &cells[ky * nbx + jx + 1],
              &cells[(ky + 1) * nbx + jx], &cells[(ky + 1) * nbx + jx + 1]}) {
          for (uint16_t site : *c)
            add(site);
        }
      }
    };

    int last_ky = -1;
    draw_cached<W, H>(
        canvas,
        [&](const math::Vector &p, int x, int) {
          const BlockCandidates<K> &cs = cands[(x - gx0) / B];
          return shade_candidates(p, cs);
        },
        [&](int y) {
          const int ky = (y - gy0) / B;
          if (ky != last_ky) {
            build_candidate_row(ky);
            last_ky = ky;
          }
        });
  }
};

} // namespace Scan
