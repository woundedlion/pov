/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file band_split.h
 * @brief Conservative clip bands and chart-straight polyline splitting. */

#include "core/memory.h"
#include "core/math/display_geometry.h"
#include "core/math/pixel_mapping.h"
#include "core/render/clip.h"
#include "core/render/plot/cull.h"
#include "core/render/plot/raster.h"

namespace Plot {

/**
 * @brief Conservative screen rows and unwrapped columns for splats reaching
 *        a clip region, and the reach test for chart-straight pieces.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H> struct ClipBand {
  float row_lo = 0.0f;   ///< Lowest row, one above the render region.
  float row_hi = 0.0f;   ///< Exclusive end row of the render region.
  float x_center = 0.0f; ///< Band centre column, unwrapped past W.
  float x_half = 0.0f;   ///< Half-width before chart padding, in columns.
  bool x_active = false; ///< Whether the clip limits columns.

  /** @brief Rows and columns a splat and the fast-trig pixel map reach past a
   *  bound. */
  static constexpr float PAD = 2.0f;

  /** @brief The band a sample must splat into to touch @p clip's render
   *  region.
   *  @param clip Clip region to cover.
   *  @return The band for @p clip. */
  HS_FLASH_MEMBER static ClipBand of(const ClipRegion &clip) {
    ClipBand band;
    const ClipRegion::XClip x_clip = clip.x_clip();
    const ColumnBand columns = ColumnBand::of(x_clip);
    band.row_lo = static_cast<float>(clip.render_y_start() - 1);
    band.row_hi = static_cast<float>(clip.render_y_end());
    band.x_active = x_clip.active;
    if (!band.x_active)
      return band;
    band.x_center = columns.centre(W);
    band.x_half = columns.clip_half(W) - PAD;
    return band;
  }

  /**
   * @brief Whether a chart-straight piece between two projected points can
   *        splat into the band.
   * @param pa Pixel coordinates of the piece's start.
   * @param pb Pixel coordinates of the piece's end.
   * @param arc The piece's chart length, in radians.
   * @return False only when the piece cannot touch the band.
   * @details The azimuthal-equidistant chart-to-sphere map is 1-Lipschitz, so
   * every point of the piece lies within colatitude (phi_a + phi_b +- arc) / 2
   * and, above the smallest colatitude sine, longitude (lambda_a + lambda_b) /
   * 2 +- arc / (2 sin phi).
   */
  __attribute__((always_inline)) bool may_reach(const math::PixelCoords &pa,
                                                const math::PixelCoords &pb,
                                                float arc) const {
    const float half_rows = 0.5f * arc * math::ROWS_PER_RADIAN<H>;
    const float mid_row = 0.5f * (pa.y + pb.y);
    const float lo = mid_row - half_rows;
    const float hi = mid_row + half_rows;
    if (hi + PAD < row_lo || lo - PAD >= row_hi)
      return false;
    if (!x_active)
      return true;
    const float phi_lo = math::DisplayGeometry<H>::row_to_phi(lo);
    const float phi_hi = math::DisplayGeometry<H>::row_to_phi(hi);
    if (phi_lo <= 0.0f || phi_hi >= math::PI_F)
      return true;
    const float sin_min =
        fminf(math::fast_sinf(phi_lo), math::fast_sinf(phi_hi));
    const float span_x = shift_column_near<W>(pb.x - pa.x, 0.0f);
    // 0.98 absorbs fast_sinf's error.
    const float half_cols =
        0.5f * arc * (W / (2.0f * math::PI_F)) / (0.98f * sin_min);
    const float centered = pa.x + 0.5f * span_x - x_center + 0.5f * W;
    const float offset =
        fabsf(centered - W * floorf(centered * (1.0f / W)) - 0.5f * W);
    return offset <= half_cols + x_half + PAD;
  }
};

/**
 * @brief Cuts a polyline's chart-straight edges into equal chart pieces
 *        and flags the runs of pieces that cannot reach a clip band, so
 *        Plot::rasterize never simulates them.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Points are added only where a run begins or ends. A run's start
 * depends on the clip, so clipped strokes match the unclipped one to a fraction
 * of a pixel, not bit for bit. Edges touching the chart antipode stay whole.
 */
template <int W, int H> class PlanarBandSplit {
public:
  /** @brief Arena bytes init_storage() takes for up to @p max_points points.
   *  @param max_points Point capacity.
   *  @return Byte count. */
  static constexpr size_t storage_bytes(int max_points) {
    return static_cast<size_t>(max_points) + alignof(uint8_t);
  }

  /** @brief Most points split() emits for @p edges edges of @p pieces each.
   *  @param edges Input edge count.
   *  @param pieces Pieces per edge.
   *  @return Point count. */
  static constexpr int64_t max_points(int edges, int pieces) {
    return static_cast<int64_t>(edges) * pieces + 1;
  }

  /** @brief Binds flag storage for polylines of up to @p max_points points.
   *  @param arena Arena the flags are allocated from.
   *  @param max_points Point capacity; at least 2. */
  HS_FLASH_MEMBER void init_storage(Arena &arena, int max_points) {
    HS_CHECK(max_points >= 2, "PlanarBandSplit: max_points %d < 2", max_points);
    capacity = max_points;
    flag_storage = arena.allocate_n<uint8_t>(max_points);
  }

  /**
   * @brief Splits @p ring into @p out and flags each output edge.
   * @param out Empty destination receiving positions only; registers are dropped. Bound for
   * max_points(edges, pieces) points.
   * @param ring Open or closed polyline with @p edges + 1 points.
   * @param edges Edge count.
   * @param pieces Chart pieces per edge.
   * @param planar_basis Azimuthal-equidistant chart the edges are straight in.
   * @param band Clip band the pieces are tested against.
   * @return Rasterize edge flags for @p out, one per edge.
   */
  HS_HOT_FLASH_MEMBER std::span<const uint8_t>
  split(Fragments &out, const Fragments &ring, int edges, int pieces,
        const math::Basis &planar_basis, const ClipBand<W, H> &band) {
    HS_CHECK(out.empty(), "PlanarBandSplit::split: out must be empty");
    HS_CHECK(edges >= 1 && pieces >= 1 &&
                 ring.size() > static_cast<size_t>(edges),
             "PlanarBandSplit: %d edges of %d pieces over %d ring points",
             edges, pieces, static_cast<int>(ring.size()));
    HS_CHECK(max_points(edges, pieces) <= capacity,
             "PlanarBandSplit: %d edges of %d pieces exceed capacity %d", edges,
             pieces, capacity);
    Fragment point;
    point.pos = ring[0].pos;
    out.push_back(point);
    auto push = [&](const math::Vector &position, bool visible)
                    __attribute__((always_inline)) {
                      flag_storage[out.size() - 1] =
                          visible ? RasterOptions::EDGE_VISIBLE : uint8_t{0};
                      point.pos = position;
                      out.push_back(point);
                    };
    const float inv_pieces = 1.0f / static_cast<float>(pieces);
    for (int edge = 0; edge < edges; ++edge) {
      const math::Vector &a = ring[edge].pos;
      const math::Vector &b = ring[edge + 1].pos;
      if (at_planar_antipode(a, planar_basis) ||
          at_planar_antipode(b, planar_basis)) {
        push(b, true);
        continue;
      }
      const auto pa = azimuthal_project(a, planar_basis);
      const auto pb = azimuthal_project(b, planar_basis);
      const float dx = pb.x - pa.x;
      const float dy = pb.y - pa.y;
      const float piece_arc = sqrtf(dx * dx + dy * dy) * inv_pieces;
      math::Vector piece_start = a;
      math::PixelCoords piece_start_px = math::vector_to_pixel<W, H>(a);
      bool run_visible = false;
      for (int j = 1; j <= pieces; ++j) {
        math::Vector piece_end = b;
        if (j < pieces) {
          const float t = static_cast<float>(j) * inv_pieces;
          piece_end = newton_unit(
              azimuthal_unproject(pa.x + dx * t, pa.y + dy * t, planar_basis));
        }
        const math::PixelCoords piece_end_px =
            math::vector_to_pixel<W, H>(piece_end);
        const bool visible =
            band.may_reach(piece_start_px, piece_end_px, piece_arc);
        if (j > 1 && visible != run_visible)
          push(piece_start, run_visible);
        run_visible = visible;
        piece_start = piece_end;
        piece_start_px = piece_end_px;
      }
      push(b, run_visible);
    }
    return {flag_storage, out.size() - 1};
  }

private:
  uint8_t *flag_storage = nullptr;
  int capacity = 0;
};

} // namespace Plot
