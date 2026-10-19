/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include <type_traits>
#include "core/memory.h"
#include "core/render/shading.h"
#include "core/color/color.h"
#include "core/render/canvas.h"
#include "core/platform/platform.h"
#include "core/render/scan/raster.h"

/**
 * @file shader.h
 * @brief Scan::Shader: full-screen per-pixel shaders with SSAA.
 */

namespace Scan {

/**
 * @brief Full-screen per-pixel shaders with SAMPLES× SSAA.
 * @details Every entry point assigns the finished premultiplied color to the
 * canvas rather than plotting it: alpha < 1 darkens the pixel instead of
 * blending, and no plot-time filter stage (Filter::World / Filter::Screen /
 * Filter::Pixel) sees it. Entry points take no pipeline and do not read
 * canvas.debug().
 */
struct Shader {
  // --- Shared SSAA helpers ---------------------------------------------------
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
      constexpr float D_THETA = 0.5f * math::PI_F / static_cast<float>(W);
      const float D_PHI = 0.25f * math::RADIANS_PER_ROW<H>;
      cos_dtheta = cosf(D_THETA);
      sin_dtheta = sinf(D_THETA);
      cos_dphi = cosf(D_PHI);
      sin_dphi = sinf(D_PHI);
    }

    /**
     * @brief Loads the two phi trig pairs for pixel row y from the LUT.
     * @param y Pixel row.
     */
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
     * @return Unit sample direction.
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

    /**
     * @brief Object-space unit vector for sample i of a pixel in the current
     *        row, in its column's object frame.
     * @param column The pixel's column frame.
     * @param i Sample index, as for at(int, int).
     * @return Unit sample direction.
     */
    HS_HOT_INLINE math::Vector at(const ObjectColumn &column, int i) const {
      const float sp = sin_phi[(i >> 1) & 1];
      const float cp = cos_phi[(i >> 1) & 1];
      const float a = sp * cos_dtheta;
      const float t = (i & 1) ? -sp * sin_dtheta : sp * sin_dtheta;
      return math::Vector(
          a * column.radial.x + t * column.tangent.x + cp * column.up.x,
          a * column.radial.y + t * column.tangent.y + cp * column.up.y,
          a * column.radial.z + t * column.tangent.z + cp * column.up.z);
    }
  };

  /**
   * @brief Validates the clip region against the LUT domain.
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
      constexpr float INV_SAMPLES = 1.0f / SAMPLES;
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
              accum += sample(grid.at(x, i), x, y) * INV_SAMPLES;
            else {
              Color4 color = sample(grid.at(x, i), x, y);
              accum += color.color * (color.alpha * INV_SAMPLES);
            }
          }

          canvas(x, y) = accum;
        });
      }
    }
  }

public:
  /** @brief Shades columns with one prepared callable per column.
   * @tparam W Canvas width.
   * @tparam H Canvas height.
   * @tparam Prepare Callable taking an absolute column and returning its shader.
   * @param canvas Destination canvas.
   * @param prepare Prepares a shader accepting a unit sphere direction or
   * sine and cosine of the row polar angle, and returning Color4 or a Pixel
   * already premultiplied by its alpha. */
  template <int W, int H, typename Prepare>
  HS_HOT_FLASH_MEMBER static void draw_columns(Canvas &canvas,
                                               Prepare &&prepare) {
    check_canvas_dims<W, H>(canvas);
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    const auto &clip = canvas.clip();
    check_lut_domain<W, H>(clip);
    walk_clip_columns<W>(clip.x_clip(), [&](int x) {
      auto shader = prepare(x);
      const float COS_THETA = math::TrigLUT<W, H>::cos_theta(x);
      const float SIN_THETA = math::TrigLUT<W, H>::sin_theta[x];
      for (int y = clip.render_y_start(); y < clip.render_y_end(); ++y) {
        const float SIN_PHI = math::TrigLUT<W, H>::sin_phi[y];
        const auto COLOR = [&]() {
          if constexpr (std::is_invocable_v<decltype(shader) &, float, float>)
            return shader(SIN_PHI, math::TrigLUT<W, H>::cos_phi[y]);
          else {
            const math::Vector VIEW(SIN_PHI * COS_THETA,
                                    math::TrigLUT<W, H>::cos_phi[y],
                                    SIN_PHI * SIN_THETA);
            return shader(VIEW);
          }
        }();
        if constexpr (std::is_same_v<std::decay_t<decltype(COLOR)>, Pixel>)
          canvas(x, y) = COLOR;
        else
          canvas(x, y) = COLOR.color * COLOR.alpha;
      }
    });
  }

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
   */
  template <int W, int H, int SAMPLES = 1>
  static void draw(Canvas &canvas, FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader) {
    static_assert(SAMPLES == 1 || SAMPLES == 4,
                  "Scan::Shader SSAA supports only SAMPLES == 1 or 4");
    HS_CHECK(vertex_shader,
             "Scan::Shader::draw requires a non-null vertex_shader");
    HS_CHECK(fragment_shader,
             "Scan::Shader::draw requires a non-null fragment_shader");
    auto shade = [&](Fragment &frag_base, const auto &grid, int x) {
      if constexpr (SAMPLES == 1) {
        fragment_shader(frag_base.pos, frag_base);
        return frag_base.color.color * frag_base.color.alpha;
      } else {
        constexpr float INV_SAMPLES = 1.0f / SAMPLES;
        Pixel accum(0, 0, 0);
        for (int i = 0; i < SAMPLES; ++i) {
          const math::Vector v = grid.at(x, i);
          Fragment sub_frag = frag_base;
          sub_frag.pos = v;
          fragment_shader(v, sub_frag);
          accum += sub_frag.color.color * (sub_frag.color.alpha * INV_SAMPLES);
        }
        return accum;
      }
    };
    draw_grid<W, H>(canvas, vertex_shader, shade);
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
   *        returns the final (premultiplied) pixel.
   * @details Writes only the columns the canvas's RenderXWindow admits.
   */
  template <int W, int H, typename VertexFn, typename PixelFn>
  HS_O3_FN static void draw_grid(Canvas &canvas, VertexFn &&vertex_shader,
                                 PixelFn &&pixel_shader) {
    walk_grid<W, H>(canvas,
                    [&](const math::Vector &center, const auto &grid, int x) {
                      Fragment frag_base;
                      frag_base.pos = center;
                      vertex_shader(frag_base);
                      return pixel_shader(frag_base, grid, x);
                    });
  }

  /**
   * @brief Full-screen SSAA-grid traversal without fragment construction.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PixelFn Callable (const Vector&, const SsaaGrid<W,H>&, int) -> Pixel.
   * @param canvas Destination canvas.
   * @param pixel_shader Computes the premultiplied pixel from its center,
   *        current row's subpixel grid, and column.
   * @details Writes only the columns the canvas's RenderXWindow admits.
   */
  template <int W, int H, typename PixelFn>
  HS_O3_FN __attribute__((always_inline)) static void
  walk_grid(Canvas &canvas, PixelFn &&pixel_shader) {
    check_canvas_dims<W, H>(canvas);
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    const auto &cr = canvas.clip();
    check_lut_domain<W, H>(cr);
    const ColumnWindow window = canvas.render_x_window();
    const auto xc = cr.x_clip();
    SsaaGrid<W, H> grid;
    for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y) {
      const float sp = math::TrigLUT<W, H>::sin_phi[y];
      const float cp = math::TrigLUT<W, H>::cos_phi[y];
      grid.set_row(y);
      walk_window_columns<W>(xc, window, [&](int x) {
        const math::Vector center(sp * math::TrigLUT<W, H>::cos_theta(x), cp,
                                  sp * math::TrigLUT<W, H>::sin_theta[x]);
        canvas(x, y) = pixel_shader(center, grid, x);
      });
    }
  }

  /**
   * @brief walk_grid() in the object frame of a turning object: each column's
   *        pixels look along its frame at the column's sweep time.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam PixelFn Callable (const Vector &, const ObjectColumn &,
   *         const SsaaGrid<W,H> &, int) -> Pixel.
   * @param canvas Destination canvas.
   * @param basis Per-column object frames of the frame's orientation.
   * @param pixel_shader Computes the premultiplied pixel from its object-space
   *        center, its column frame (`grid.at(column, i)` gives the object-space
   *        samples), the current row's subpixel grid, and the column.
   * @details Writes only the columns the canvas's RenderXWindow admits.
   */
  template <int W, int H, typename PixelFn>
  HS_O3_FN HS_HOT_INLINE static void
  walk_object_grid(Canvas &canvas, const ObjectColumnBasis<W, H> &basis,
                   PixelFn &&pixel_shader) {
    check_canvas_dims<W, H>(canvas);
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    const auto &cr = canvas.clip();
    check_lut_domain<W, H>(cr);
    const ColumnWindow window = canvas.render_x_window();
    const auto xc = cr.x_clip();
    SsaaGrid<W, H> grid;
    for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y) {
      const float sp = math::TrigLUT<W, H>::sin_phi[y];
      const float cp = math::TrigLUT<W, H>::cos_phi[y];
      grid.set_row(y);
      walk_window_columns<W>(xc, window, [&](int x) HS_HOT_INLINE {
        const ObjectColumn column = basis.column(x);
        canvas(x, y) = pixel_shader(column.direction(sp, cp), column, grid, x);
      });
    }
  }
  /** @brief K site indices classified at one block-grid corner. */
  template <size_t K> using BlockCell = std::array<uint16_t, K>;

  /** @brief Deduped union of four corners, with contiguous shading positions. */
  template <size_t K> struct BlockCandidates {
    static_assert(K > 0 && 4 * K <= UINT8_MAX,
                  "block candidates require 1..63 sites per corner");
    math::Vector pos[4 * K]; ///< Site positions, parallel to `idx`.
    uint16_t idx[4 * K];     ///< Distinct site indices; first `n` valid.
    uint8_t n;               ///< Number of valid entries.
  };

  /** @brief A block's four corner cells and their deduped site union. */
  template <size_t K> struct BlockSites {
    static_assert(K > 0 && 4 * K <= UINT8_MAX,
                  "block sites require 1..63 sites per corner");
    /// Corners in (x, y), (x + 1, y), (x, y + 1), (x + 1, y + 1) order.
    std::array<const BlockCell<K> *, 4> corners;
    uint16_t idx[4 * K]; ///< Distinct site indices; first `n` valid.
    uint8_t n;           ///< Number of valid entries.
  };

private:
  /** @brief Canvas-anchored block corners classified over the render clip. */
  template <size_t K> struct BlockGrid {
    int block;           ///< Block edge in pixels.
    int gx0;             ///< Column of corner 0.
    int gy0;             ///< Row of corner 0.
    int nbx;             ///< Corner columns.
    int nby;             ///< Corner rows.
    int gap_begin;       ///< First block column a wrapped clip skips.
    int gap_end;         ///< Block column after the wrapped clip's gap.
    BlockCell<K> *cells; ///< Corner cells, row-major; null when empty.

    /** @brief Whether block column @p jx lies in the wrapped clip's gap.
     *  @param jx Block column.
     *  @return True when the clip skips the block. */
    bool skips(int jx) const { return jx >= gap_begin && jx < gap_end; }
  };

  /**
   * @brief Classifies the block corners covering the render clip.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam K Sites classified per corner.
   * @param cr Render clip.
   * @param block Positive block edge in pixels.
   * @param scratch Arena for the corner cells.
   * @param classify_corner Callable (const Vector&) returning BlockCell<K>.
   * @return The classified grid; `cells` is null when the clip is empty.
   */
  template <int W, int H, size_t K, typename ClassifyFn>
  __attribute__((always_inline)) static BlockGrid<K>
  classify_block_grid(const ClipRegion &cr, int block, Arena &scratch,
                      ClassifyFn &classify_corner) {
    HS_CHECK(block > 0, "block-coherent shader requires a positive block size");
    const auto columns = cr.x_clip();
    const int x0 = columns.active && !columns.wrap ? columns.rs : 0;
    const int x1 = columns.active && !columns.wrap ? columns.re : W;
    const int y0 = cr.render_y_start();
    const int y1 = cr.render_y_end();
    BlockGrid<K> grid{};
    if (x1 <= x0 || y1 <= y0)
      return grid;
    if (!math::TrigLUT<W, H>::initialized)
      math::TrigLUT<W, H>::init();
    const int B = block;
    grid.block = B;
    // Canvas-anchored corners give each pixel the same candidates in every band.
    grid.gx0 = (x0 / B) * B;
    grid.gy0 = (y0 / B) * B;
    grid.nbx = (x1 - 1 - grid.gx0) / B + 2; // corner columns spanning [x0, x1)
    grid.nby = (y1 - 1 - grid.gy0) / B + 2; // corner rows spanning    [y0, y1)
    grid.gap_begin =
        columns.active && columns.wrap ? (columns.re + B - 1) / B : grid.nbx;
    grid.gap_end = columns.active && columns.wrap ? columns.rs / B : grid.nbx;
    grid.cells = scratch.allocate_n<BlockCell<K>>(grid.nbx * grid.nby);
    // End corners clamp to the last canvas pixel, independent of the clip.
    auto corner_x = [&](int j) { return std::min(grid.gx0 + j * B, W - 1); };
    auto corner_y = [&](int k) { return std::min(grid.gy0 + k * B, H - 1); };

    for (int k = 0; k < grid.nby; ++k)
      for (int j = 0; j < grid.nbx; ++j) {
        if (j > grid.gap_begin && j < grid.gap_end)
          continue;
        grid.cells[k * grid.nbx + j] = classify_corner(
            math::pixel_to_vector<W, H>(corner_x(j), corner_y(k)));
      }
    return grid;
  }

public:
  /** @brief Maximum corner-grid and candidate-row bytes at a minimum block size.
   * @details Excludes caller-owned site positions and any arena alignment pad.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam K Sites classified per corner.
   * @tparam MIN_BLOCK Minimum block edge, in pixels.
   * @return Scratch byte count.
   */
  template <int W, int H, size_t K, int MIN_BLOCK>
  static constexpr size_t block_coherent_scratch_bytes() {
    static_assert(W > 0 && H > 0 && MIN_BLOCK > 0);
    constexpr size_t COLS = (W - 1) / MIN_BLOCK + 2;
    constexpr size_t ROWS = (H - 1) / MIN_BLOCK + 2;
    return COLS * ROWS * sizeof(BlockCell<K>) +
           (COLS - 1) * sizeof(BlockCandidates<K>);
  }

  /** @brief Maximum draw_block_columns() scratch at a minimum block size.
   * @details Excludes any arena alignment pad.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam K Sites classified per corner.
   * @tparam MIN_BLOCK Minimum block edge, in pixels.
   * @return Scratch byte count.
   */
  template <int W, int H, size_t K, int MIN_BLOCK>
  static constexpr size_t block_columns_scratch_bytes() {
    static_assert(W > 0 && H > 0 && MIN_BLOCK > 0);
    constexpr size_t COLS = (W - 1) / MIN_BLOCK + 2;
    constexpr size_t ROWS = (H - 1) / MIN_BLOCK + 2;
    return COLS * ROWS * sizeof(BlockCell<K>) +
           (ROWS - 1) * sizeof(BlockSites<K>);
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
    check_canvas_dims<W, H>(canvas);
    const auto &cr = canvas.clip();
    check_lut_domain<W, H>(cr);
    const BlockGrid<K> grid =
        classify_block_grid<W, H, K>(cr, block, scratch, classify_corner);
    if (!grid.cells)
      return;
    const int B = grid.block;
    const int gx0 = grid.gx0;
    const int gy0 = grid.gy0;
    const int nbx = grid.nbx;
    const BlockCell<K> *cells = grid.cells;

    // One candidate set per block column, rebuilt on each block-row change.
    const int nblk = nbx - 1;
    BlockCandidates<K> *cands = scratch.allocate_n<BlockCandidates<K>>(nblk);
    auto build_candidate_row = [&](int ky) {
      for (int jx = 0; jx < nblk; ++jx) {
        if (grid.skips(jx))
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

  /**
   * @brief Shades block-coherent candidates column by column.
   * @tparam W Canvas width in pixels.
   * @tparam H Canvas height in pixels.
   * @tparam K Sites classified per corner.
   * @param canvas Destination canvas.
   * @param block Positive block edge in pixels; corners anchor to the canvas.
   * @param scratch Arena for the corner grid and one block column's site
   * unions. The caller owns its scope.
   * @param classify_corner Callable (const Vector&) returning BlockCell<K>.
   * @param prepare_column Callable (int x) returning a column object whose
   * `block(const BlockSites<K>&)` returns the block's shader, a callable
   * (const Vector&) returning Color4 or a premultiplied Pixel.
   * @details Corners and unions are those of draw_block_coherent(). Each
   * clipped column is prepared once, then its blocks are prepared and shaded
   * top to bottom.
   */
  template <int W, int H, size_t K, typename ClassifyFn, typename PrepareFn>
  HS_O3_FN __attribute__((always_inline)) static void
  draw_block_columns(Canvas &canvas, int block, Arena &scratch,
                     ClassifyFn &&classify_corner, PrepareFn &&prepare_column) {
    check_canvas_dims<W, H>(canvas);
    const auto &cr = canvas.clip();
    check_lut_domain<W, H>(cr);
    const BlockGrid<K> grid =
        classify_block_grid<W, H, K>(cr, block, scratch, classify_corner);
    if (!grid.cells)
      return;
    const int B = grid.block;
    const int y0 = cr.render_y_start();
    const int y1 = cr.render_y_end();
    const int nblk = grid.nby - 1;
    BlockSites<K> *sites = scratch.allocate_n<BlockSites<K>>(nblk);
    auto build_block_column = [&](int jx) {
      for (int ky = 0; ky < nblk; ++ky) {
        BlockSites<K> &set = sites[ky];
        const BlockCell<K> *row = grid.cells + ky * grid.nbx + jx;
        set.corners = {row, row + 1, row + grid.nbx, row + grid.nbx + 1};
        set.n = 0;
        for (const BlockCell<K> *corner : set.corners)
          for (uint16_t site : *corner) {
            bool seen = false;
            for (uint8_t i = 0; i < set.n; ++i)
              seen |= set.idx[i] == site;
            if (!seen)
              set.idx[set.n++] = site;
          }
      }
    };

    int built = -1;
    walk_clip_columns<W>(cr.x_clip(), [&](int x) {
      const int jx = (x - grid.gx0) / B;
      if (jx != built) {
        build_block_column(jx);
        built = jx;
      }
      auto column = prepare_column(x);
      const float COS_THETA = math::TrigLUT<W, H>::cos_theta(x);
      const float SIN_THETA = math::TrigLUT<W, H>::sin_theta[x];
      for (int ky = 0; ky < nblk; ++ky) {
        const int ROW_BEGIN = std::max(y0, grid.gy0 + ky * B);
        const int ROW_END = std::min(y1, grid.gy0 + (ky + 1) * B);
        auto shade = column.block(sites[ky]);
        for (int y = ROW_BEGIN; y < ROW_END; ++y) {
          const float SIN_PHI = math::TrigLUT<W, H>::sin_phi[y];
          const auto COLOR = shade(math::Vector(SIN_PHI * COS_THETA,
                                                math::TrigLUT<W, H>::cos_phi[y],
                                                SIN_PHI * SIN_THETA));
          if constexpr (std::is_same_v<std::decay_t<decltype(COLOR)>, Pixel>)
            canvas(x, y) = COLOR;
          else
            canvas(x, y) = COLOR.color * COLOR.alpha;
        }
      }
    });
  }
};

} // namespace Scan
