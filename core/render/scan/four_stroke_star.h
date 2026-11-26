/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file four_stroke_star.h
 * @brief Projected exposed edges of a stationary spherical star. */

#include "core/render/four_stroke.h"
#include "core/render/sdf/shapes.h"

namespace Scan {

/** @brief Piecewise-linear output contour of a star in its azimuthal chart.
 * @tparam W,H Output dimensions. */
template <int W, int H> class ProjectedStar {
  using Point = FourStroke::Direction;
  static constexpr int CAPACITY = 80;
  struct Edge {
    Point a, b;
    Point normal;
    float inv_length;
  };
  std::array<Edge, CAPACITY> edges;
  int count = 0;
  bool valid = false;
  float anchor = 0, orientation = 0;
  float x_min = INFINITY, x_max = -INFINITY, y_min = INFINITY,
        y_max = -INFINITY;

  Point project(const SDF::Star &star, Point local) const {
    const float r = std::hypot(local.x, local.y);
    const float s = r > 0 ? std::sin(r) / r : 1;
    const auto p = star.basis.v * std::cos(r) +
                   (star.basis.u * local.x + star.basis.w * local.y) * s;
    const float x = std::atan2(p.z, p.x) * (W / math::TWO_PI_F);
    return {anchor + std::remainder(x - anchor, float(W)),
            math::DisplayGeometry<H>::phi_to_row(
                std::atan2(std::hypot(p.x, p.z), p.y))};
  }

  bool append(const SDF::Star &star, Point a, Point b, Point pa, Point pb,
              int depth) {
    const Point m{(a.x + b.x) * 0.5f, (a.y + b.y) * 0.5f};
    const auto pm = project(star, m);
    const float ex = pm.x - (pa.x + pb.x) * 0.5f,
                ey = pm.y - (pa.y + pb.y) * 0.5f;
    if (ex * ex + ey * ey > 1.0f / (64 * 64)) {
      return depth < 7 && append(star, a, m, pa, pm, depth + 1) &&
             append(star, m, b, pm, pb, depth + 1);
    }
    if (count == CAPACITY)
      return false;
    const float dx = pb.x - pa.x, dy = pb.y - pa.y;
    const float len = std::hypot(dx, dy);
    if (len == 0)
      return true;
    edges[count++] = {pa, pb, {dy / len, -dx / len}, 1 / len};
    orientation += (pa.x - anchor) * pb.y - (pb.x - anchor) * pa.y;
    x_min = std::min(x_min, pa.x);
    x_max = std::max(x_max, pa.x);
    y_min = std::min(y_min, pa.y);
    y_max = std::max(y_max, pa.y);
    return true;
  }

public:
  /** @brief Builds a contiguous longitude chart for the star's visible perimeter.
   * @param star Solid SDF star.
   * @return False for a pole-crossing, inverted or oversized contour. */
  bool prepare(const SDF::Star &star) {
    valid = false;
    count = 0;
    orientation = 0;
    x_min = y_min = INFINITY;
    x_max = y_max = -INFINITY;
    const float phi = std::acos(std::clamp(star.basis.v.y, -1.0f, 1.0f));
    if (star.sign < 0 || star.sides > CAPACITY / 2 ||
        star.circumradius >= std::min(phi, math::PI_F - phi))
      return false;
    anchor = std::atan2(star.basis.v.z, star.basis.v.x) * (W / math::TWO_PI_F);
    const auto vertex = [&](int k) {
      const float r =
          star.circumradius * (k % 2 ? SDF::STAR_INNER_RATIO : 1.0f);
      const float theta = k * (math::PI_F / star.sides) - star.phase;
      return Point{r * std::cos(theta), r * std::sin(theta)};
    };
    Point a = vertex(0), pa = project(star, a);
    const Point first = a, first_pixel = pa;
    for (int k = 1; k <= 2 * star.sides; ++k) {
      const bool closes = k == 2 * star.sides;
      const Point b = closes ? first : vertex(k),
                  pb = closes ? first_pixel : project(star, b);
      if (!append(star, a, b, pa, pb, 0))
        return false;
      a = b;
      pa = pb;
    }
    if (orientation < 0)
      for (int i = 0; i < count; ++i) {
        edges[i].normal.x *= -1;
        edges[i].normal.y *= -1;
      }
    valid = count >= 3;
    return valid;
  }

  /** @brief Whether the latest preparation produced a usable projected contour.
   * @return True for a complete exposed outline. */
  bool ready() const { return valid; }
  /** @brief First potentially covered column in the unwrapped chart.
   * @return Inclusive column limit. */
  int begin_column() const { return static_cast<int>(std::ceil(x_min - 0.5f)); }
  /** @brief End of the contour's unwrapped column range.
   * @return Exclusive column limit. */
  int end_column() const {
    return static_cast<int>(std::floor(x_max + 0.5f)) + 1;
  }

  /** @brief Maps a pixel column to the contour's unwrapped chart.
   * @param x Output column.
   * @return Unwrapped column. */
  float unwrap(float x) const {
    return anchor + std::remainder(x - anchor, float(W));
  }
  /** @brief Tests a pixel footprint against the contour's horizontal extent.
   * @param x Unwrapped column center.
   * @return Whether its subpixel taps may touch the contour. */
  bool admits_x(float x) const {
    return x >= x_min - 0.5f && x <= x_max + 0.5f;
  }
  /** @brief First potentially covered row.
   * @return Inclusive row limit. */
  int begin_row() const { return static_cast<int>(std::ceil(y_min - 0.5f)); }
  /** @brief End of the covered row range.
   * @return Exclusive row limit. */
  int end_row() const { return static_cast<int>(std::floor(y_max + 0.5f)) + 1; }

  /** @brief Ordered edge indices for the samples in one raster row. */
  struct RowEdges {
    std::array<uint8_t, CAPACITY>
        indices;  ///< Contour indices in winding order.
    uint8_t size; ///< Number of indices in the row's vertical band.
  };

  /** @brief Selects the edges that meet a vertical sample band.
   * @param top,bottom Inclusive row coordinates of the band.
   * @return Ordered candidates for coverage or local contour distance. */
  RowEdges row_edges(float top, float bottom) const {
    RowEdges result{};
    for (int i = 0; i < count; ++i) {
      const auto &edge = edges[i];
      if (fmaxf(edge.a.y, edge.b.y) >= top &&
          fminf(edge.a.y, edge.b.y) <= bottom)
        result.indices[result.size++] = static_cast<uint8_t>(i);
    }
    return result;
  }

  /** @brief Tests coverage without resolving contour distance or direction.
   * @param x,y Sample coordinates in the contour's unwrapped chart.
   * @param row Candidates spanning the sample's row.
   * @return True inside or on the exposed outline. */
  bool contains(float x, float y, const RowEdges &row) const {
    bool inside = false;
    for (int i = 0; i < row.size; ++i) {
      const auto &edge = edges[row.indices[i]];
      const float dx = edge.b.x - edge.a.x, dy = edge.b.y - edge.a.y;
      const float px = x - edge.a.x, py = y - edge.a.y;
      if (y < fminf(edge.a.y, edge.b.y) || y > fmaxf(edge.a.y, edge.b.y))
        continue;
      const float cross = px * dy - py * dx;
      if (cross == 0 && px * dx + py * dy >= 0 &&
          (x - edge.b.x) * dx + (y - edge.b.y) * dy <= 0)
        return true;
      if ((edge.a.y > y) != (edge.b.y > y) && x < edge.a.x + py * dx / dy)
        inside = !inside;
    }
    return inside;
  }

  /** @brief Classifies one sample against the exposed segments.
   * @param x,y Unwrapped sample coordinates.
   * @param align Alignment of the physical or output heading with a unit normal.
   * @param alignment_width Maximum depth that needs directional alignment.
   * @param row Optional candidates covering the sample and alignment width.
   * @return Membership and closest exposed contour in output pixels. */
  FourStroke::Boundary query(float x, float y, const auto &align,
                             float alignment_width = INFINITY,
                             const RowEdges *row = nullptr) const {
    FourStroke::Boundary result;
    float nearest = INFINITY;
    const float limit = alignment_width * alignment_width;
    int winner = 0;
    bool tied = false;
    const int row_size = row ? row->size : count;
    for (int offset = 0; offset < row_size; ++offset) {
      const int i = row ? row->indices[offset] : offset;
      const auto &e = edges[i];
      const float dx = e.b.x - e.a.x, dy = e.b.y - e.a.y;
      if ((e.a.y > y) != (e.b.y > y) && x < e.a.x + (y - e.a.y) * dx / dy)
        result.inside = !result.inside;
      const float u = hs::clamp(((x - e.a.x) * dx + (y - e.a.y) * dy) *
                                    e.inv_length * e.inv_length,
                                0.0f, 1.0f);
      const float px = x - (u == 1 ? e.b.x : e.a.x + u * dx),
                  py = y - (u == 1 ? e.b.y : e.a.y + u * dy);
      const float sq = px * px + py * py;
      if (sq < nearest) {
        nearest = sq;
        winner = i;
        tied = false;
      } else if (sq == nearest && sq < limit) {
        if (!tied) {
          result.depth = nearest;
          result.alignment = align(edges[winner].normal);
        }
        tied = true;
        FourStroke::consider(result, sq, align(e.normal));
      }
    }
    result.inside |= nearest == 0;
    if (result.inside && nearest < limit && !tied)
      result.alignment = align(edges[winner].normal);
    result.depth = std::sqrt(nearest);
    return result;
  }
};
} // namespace Scan
