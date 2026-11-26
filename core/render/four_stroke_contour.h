/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file four_stroke_contour.h
 * @brief Marching-squares contour and subpixel membership on a fixed-time output lattice.
 */

#include "core/platform/platform.h"
#include "core/render/four_stroke.h"

namespace FourStroke {

/// A three-cell halo shared by all subpoints of one output pixel.
class ContourTile {
public:
  /** @brief Samples the signed field on the output-pixel lattice.
   * @tparam Field Signed scalar or signed binary occupancy query.
   * @param field Callable taking integer column and row offsets from the center.
   *              Nonpositive values are inside. A nonfinite value marks a missing sample.
   * @param radius Lattice radius, from 1 to 3. Queries must stay in these cells. */
  template <typename Field> void prepare(Field &&field, int radius = 3) {
    for (int y = -radius; y <= radius; ++y)
      for (int x = -radius; x <= radius; ++x)
        values[y + 3][x + 3] = field(x, y);
#if HS_ENABLE_STRUCTURAL_AUDITS
    prepared_radius = radius;
#endif
  }

  /** @brief Classifies a subpixel against the reconstructed exposed segments.
   * @tparam MembershipOnly Evaluate only the owner's coverage.
   * @tparam Align Normal-to-direction alignment query.
   * @tparam ValidCell Chart continuity predicate.
   * @param dx Pixel-relative column coordinate in [-0.5, 0.5].
   * @param dy Pixel-relative row coordinate in [-0.5, 0.5].
   * @param align Returns signed alignment for an outward unit normal.
   * @param valid_cell Accepts the integer origin of cells with valid interpolation.
   * @return Reconstructed membership and nearest exposed feature in the tile. */
  template <bool MembershipOnly = false, typename Align, typename ValidCell>
  Boundary query(float dx, float dy, Align &&align,
                 ValidCell &&valid_cell) const {
    Boundary nearest;
    const int owner_x = static_cast<int>(std::floor(dx));
    const int owner_y = static_cast<int>(std::floor(dy));
    const int x_start = MembershipOnly ? owner_x + 3 : 0;
    const int y_start = MembershipOnly ? owner_y + 3 : 0;
    const int x_end = MembershipOnly ? x_start + 1 : 6;
    const int y_end = MembershipOnly ? y_start + 1 : 6;
    for (int y = y_start; y < y_end; ++y) {
      for (int x = x_start; x < x_end; ++x) {
        if (!valid_cell(x - 3, y - 3))
          continue;
        HS_AUDIT_CHECK(x - 3 >= -prepared_radius && x - 2 <= prepared_radius &&
                           y - 3 >= -prepared_radius &&
                           y - 2 <= prepared_radius,
                       "ContourTile: query cell outside the prepared radius");
        const float field[4] = {values[y][x], values[y][x + 1],
                                values[y + 1][x + 1], values[y + 1][x]};
        if (!(std::isfinite(field[0]) && std::isfinite(field[1]) &&
              std::isfinite(field[2]) && std::isfinite(field[3])))
          continue;
        const bool owner = owner_x == x - 3 && owner_y == y - 3;
        const float px = dx - (x - 3), py = dy - (y - 3);
        const int mask = (field[0] <= 0) | ((field[1] <= 0) << 1) |
                         ((field[2] <= 0) << 2) | ((field[3] <= 0) << 3);
        if (mask == 0 || mask == 15) {
          if (owner)
            nearest.inside = mask == 15;
          continue;
        }

        constexpr float CX[] = {0, 1, 1, 0}, CY[] = {0, 0, 1, 1};
        Direction crossing[4];
        int exits[2], entries[2];
        int n_exits = 0, n_entries = 0;
        for (int i = 0; i < 4; ++i) {
          const int next = (i + 1) % 4;
          const bool inside = field[i] <= 0;
          if (inside == (field[next] <= 0))
            continue;
          const float t =
              static_cast<float>(static_cast<double>(field[i]) /
                                 (static_cast<double>(field[i]) - field[next]));
          crossing[i] = {CX[i] + t * (CX[next] - CX[i]),
                         CY[i] + t * (CY[next] - CY[i])};
          if (inside)
            exits[n_exits++] = i;
          else
            entries[n_entries++] = i;
        }
        const bool connected =
            n_exits == 2 &&
            static_cast<double>(field[0]) + field[1] + field[2] + field[3] < 0;
        bool inside = connected;
        for (int i = 0; i < n_exits; ++i) {
          const int exit = exits[i];
          const int entry =
              n_exits == 1 ? entries[0] : (exit + (connected ? 1 : 3)) % 4;
          const auto a = crossing[exit], b = crossing[entry];
          const float vx = b.x - a.x, vy = b.y - a.y;
          const float rx = px - a.x, ry = py - a.y;
          const float length_sq = vx * vx + vy * vy;
          const bool half_plane =
              length_sq == 0 ? rx == 0 && ry == 0 : rx * vy - ry * vx <= 0;
          inside = connected ? inside && half_plane : inside || half_plane;
          if constexpr (!MembershipOnly) {
            if (length_sq == 0)
              continue;
            const float t =
                hs::clamp((rx * vx + ry * vy) / length_sq, 0.0f, 1.0f);
            const float ex = rx - t * vx, ey = ry - t * vy;
            const float distance = std::sqrt(ex * ex + ey * ey);
            if (distance <= nearest.depth) {
              const float length = std::sqrt(length_sq);
              consider(nearest, distance, align({vy / length, -vx / length}));
            }
          }
        }
        if (owner)
          nearest.inside = inside;
      }
    }
    return nearest;
  }

  /** @brief Classifies a subpixel on a continuous chart.
   * @details Visits every cell of the halo, so the tile must be prepared
   *          with radius 3.
   * @tparam Align Normal-to-direction alignment query.
   * @param dx Pixel-relative column coordinate in [-0.5, 0.5].
   * @param dy Pixel-relative row coordinate in [-0.5, 0.5].
   * @param align Returns signed alignment for an outward unit normal.
   * @return Reconstructed membership and nearest exposed feature. */
  template <typename Align>
  Boundary query(float dx, float dy, Align &&align) const {
    return query(dx, dy, align, [](int, int) { return true; });
  }

  /** @brief Tests coverage within the reconstructed halo.
   * @details The cell containing (dx, dy) must lie within the prepared radius.
   * @param dx,dy Pixel-relative coordinates.
   * @return True inside the interpolated contour. */
  bool contains(float dx, float dy) const {
    if (dx < -3 || dx >= 3 || dy < -3 || dy >= 3)
      return false;
    return query<true>(
               dx, dy, [](Direction) { return 0.f; },
               [](int, int) { return true; })
        .inside;
  }

private:
  float values[7][7]; ///< Signed samples centered on the output pixel.
#if HS_ENABLE_STRUCTURAL_AUDITS
  int prepared_radius = 0; ///< Lattice radius of the last prepare().
#endif
};

} // namespace FourStroke
