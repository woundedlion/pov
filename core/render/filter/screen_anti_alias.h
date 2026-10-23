/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "core/render/filter/splat.h"
#include "core/color/color.h"

/**
 * @file screen_anti_alias.h
 * @brief Filter::Screen::AntiAlias: splats a sub-pixel sample across its four
 * nearest pixel neighbors.
 */

namespace Filter {
namespace Screen {

HS_O3_BEGIN

/**
 * @brief Applies 2D anti-aliasing to sub-pixel coordinates.
 * @details Distributes intensity to the 4 nearest neighbors using a quintic kernel.
 */
template <int W, int H> class AntiAlias : public Is2D {
public:
  /** @brief Tap-position error of splatting a linearly moving sample at one time.
 * @param x Column of the splat.
 * @param dx,dy Column and row motion over the frame.
 * @return Bound, in pixels, on how far a tap's own-time position lies from
 *         the splat's: the motion over the farther tap column's lag. */
  HS_HOT_INLINE static inline float tap_motion_error(float x, float dx,
                                                     float dy) {
    const float fraction = x - floorf(x);
    return (fabsf(dx) + fabsf(dy)) *
           (fmaxf(fraction, 1.0f - fraction) / (W / 2.0f));
  }

  static constexpr bool has_boundary_splat =
      true; ///< Provides a separate boundary splat.
  /** @brief Splats a timed sample with arm-boundary tap ownership.
   * @param x Column coordinate.
   * @param y Row coordinate.
   * @param c Source color.
   * @param age Age register.
   * @param alpha Source opacity.
   * @param time Frame fraction.
   * @param pass Downstream tap continuation.
   */
  template <typename PassFnT>
  HS_O3_FN HS_NOINLINE_NOCLONE void plot_at(float x, float y, const ::Pixel &c,
                                            float age, float alpha, float time,
                                            PassFnT &&pass) {
    const bool boundary = !Spacetime::flat_surface() &&
                          fabsf(math::centered_remainder(x, W / 2.0f)) <= 1.0f;
    if (boundary) {
      plot_boundary_static(x, y, c, age, alpha, time, pass);
      return;
    }
    plot(x, y, c, age, alpha, pass);
  }

  /** @brief Emits a sample with arm-pass tap ownership.
   * @tparam PassFnT Downstream tap continuation type.
   * @param x,y Projected sample coordinates.
   * @param c,age,alpha Shading registers.
   * @param time Owning arm fraction.
   * @param pass Downstream tap continuation.
   */
  template <typename PassFnT>
  HS_HOT_INLINE void plot_boundary_at(float x, float y, const ::Pixel &c,
                                      float age, float alpha, float time,
                                      PassFnT &pass) {
    plot_boundary_static(x, y, c, age, alpha, time, pass);
  }

  /** @brief Visits a splat's tap columns lit within half a frame of its
   *         own time.
   * @tparam VisitFnT Callable `(int column, float fraction)`.
   * @param x Projected column.
   * @param time Owning arm fraction.
   * @param visit Receives each kept wrapped column and its sweep fraction. */
  template <typename VisitFnT>
  HS_HOT_INLINE static void for_each_owned_column(float x, float time,
                                                  VisitFnT &&visit) {
    const int x0 = math::fast_wrap(static_cast<int>(floorf(x)), W);
    const int columns[] = {x0, math::fast_wrap(x0 + 1, W)};
    for (int column : columns) {
      const float at =
          Spacetime::SamplingSurface<W>::time(static_cast<float>(column));
      if (fabsf(at - time) <= 0.5f)
        visit(column, at);
    }
  }

  /** @brief Visits the tap columns of a splat beside an arm seam, each with
   *         the sweep fraction it is lit at.
   * @tparam VisitFnT Callable `(int column, float fraction)`.
   * @param x Projected column, within a pixel of a seam.
   * @param visit Receives each distinct wrapped column and its fraction.
   * @details The splat's own two columns are joined by the two columns
   * astride the nearest seam, lit a frame apart, so both arm passes hold the
   * taps they own. */
  template <typename VisitFnT>
  HS_HOT_INLINE static void for_each_seam_column(float x, VisitFnT &&visit) {
    const int x0 = math::fast_wrap(static_cast<int>(floorf(x)), W);
    const int arm = (int(floorf(x / (W / 2.0f) + 0.5f)) & 1) * (W / 2);
    const int columns[] = {x0, math::fast_wrap(x0 + 1, W),
                           math::fast_wrap(arm - 1, W), arm};
    for (int i = 0; i < 4; ++i)
      if (i < 2 || (columns[i] != columns[0] && columns[i] != columns[1]))
        visit(columns[i], Spacetime::SamplingSurface<W>::time(
                              static_cast<float>(columns[i])));
  }

  /** @brief Emits the taps of a splat at (x, y) that fall in one column.
   * @tparam PassFnT Downstream tap continuation type.
   * @param x,y Projected coordinates at the column time.
   * @param c,age,alpha Shading registers.
   * @param column Wrapped destination column.
   * @param pass Downstream tap continuation. */
  template <typename PassFnT>
  HS_HOT_INLINE static void column_taps(float x, float y, const ::Pixel &c,
                                        float age, float alpha, int column,
                                        PassFnT &pass) {
    assert(age >= 0.0f && alpha >= 0.0f);
    const SplatTaps t = splat_taps<W, H>(x, y);
    if constexpr (W > 1) {
      if (column != t.x0 && column != t.x1)
        return;
      const float TOP = column == t.x0 ? t.v00 : t.v10;
      const float BOTTOM = column == t.x0 ? t.v01 : t.v11;
      if (t.y0_physical && TOP > SPLAT_TAP_CUTOFF)
        pass(static_cast<float>(column), static_cast<float>(t.y0), c, age,
             alpha * TOP);
      if (t.y1_physical && BOTTOM > SPLAT_TAP_CUTOFF)
        pass(static_cast<float>(column), static_cast<float>(t.y1), c, age,
             alpha * BOTTOM);
      return;
    }
    if (t.y0_physical && t.x0 == column && t.v00 > SPLAT_TAP_CUTOFF)
      pass(static_cast<float>(t.x0), static_cast<float>(t.y0), c, age,
           alpha * t.v00);
    if (t.y0_physical && t.x1 == column && t.v10 > SPLAT_TAP_CUTOFF)
      pass(static_cast<float>(t.x1), static_cast<float>(t.y0), c, age,
           alpha * t.v10);
    if (t.y1_physical && t.x0 == column && t.v01 > SPLAT_TAP_CUTOFF)
      pass(static_cast<float>(t.x0), static_cast<float>(t.y1), c, age,
           alpha * t.v01);
    if (t.y1_physical && t.x1 == column && t.v11 > SPLAT_TAP_CUTOFF)
      pass(static_cast<float>(t.x1), static_cast<float>(t.y1), c, age,
           alpha * t.v11);
  }

  /** @brief column_taps() out of line.
   * @tparam PassFnT Downstream tap continuation type.
   * @param x,y Projected coordinates at the column time.
   * @param c,age,alpha Shading registers.
   * @param column Wrapped destination column.
   * @param pass Downstream tap continuation. */
  template <typename PassFnT>
  HS_O3_FN HS_NOINLINE_NOCLONE void
  plot_owned_column(float x, float y, const ::Pixel &c, float age, float alpha,
                    int column, PassFnT &pass) {
    column_taps(x, y, c, age, alpha, column, pass);
  }

private:
  /** @brief Filters a stationary splat by arm-pass ownership.
   * @tparam PassFnT Downstream tap continuation type.
   * @param x,y Projected sample coordinates.
   * @param c,age,alpha Shading registers.
   * @param time Owning arm fraction.
   * @param pass Downstream tap continuation. */
  template <typename PassFnT>
  HS_HOT_FLASH_MEMBER void
  plot_boundary_static(float x, float y, const ::Pixel &c, float age,
                       float alpha, float time, PassFnT &pass) {
    plot(x, y, c, age, alpha,
         [&](float nx, float ny, const ::Pixel &nc, float nage,
             float nalpha) HS_HOT_INLINE {
           if (fabsf(Spacetime::SamplingSurface<W>::time(nx) - time) <= 0.5f)
             pass(nx, ny, nc, nage, nalpha);
         });
  }

public:
  /// Splat taps land up to one pixel away.
  static constexpr int segment_margin = 1;
  /// Taps land on integer pixel coordinates.
  static constexpr bool emits_pixel_centers = true;
  /** @brief The whole splat is the sub-pixel fraction: at an integer
   * coordinate the four taps collapse to one at full weight. */
  static constexpr bool requires_subpixel_input = true;

  /**
   * @brief Splats a sub-pixel sample across its four nearest pixel neighbors.
   * @tparam PassFnT Downstream 2D callback type.
   * @param x Sub-pixel column coordinate.
   * @param y Sub-pixel row coordinate.
   * @param c Source color, forwarded to each tap.
   * @param age Temporal age channel (frames), forwarded unchanged.
   * @param alpha Blend alpha in [0, 1]; scaled per tap by its quintic-eased splat weight.
   * @param pass Downstream 2D callback receiving each weighted tap.
   */
  template <typename PassFnT>
  HS_HOT_INLINE void plot(float x, float y, const ::Pixel &c, float age,
                          float alpha, PassFnT &&pass) {
    assert(age >= 0.0f && alpha >= 0.0f);
    const SplatTaps t = splat_taps<W, H>(x, y);

    if (t.y0_physical && t.v00 > SPLAT_TAP_CUTOFF)
      pass(static_cast<float>(t.x0), static_cast<float>(t.y0), c, age,
           alpha * t.v00);
    if (t.y0_physical && t.v10 > SPLAT_TAP_CUTOFF)
      pass(static_cast<float>(t.x1), static_cast<float>(t.y0), c, age,
           alpha * t.v10);
    if (t.y1_physical && t.v01 > SPLAT_TAP_CUTOFF)
      pass(static_cast<float>(t.x0), static_cast<float>(t.y1), c, age,
           alpha * t.v01);
    if (t.y1_physical && t.v11 > SPLAT_TAP_CUTOFF)
      pass(static_cast<float>(t.x1), static_cast<float>(t.y1), c, age,
           alpha * t.v11);
  }

  /** @brief Splats a sample straight into a terminal that writes taps.
   * @tparam TerminalT Terminal pipeline.
   * @param terminal Destination terminal.
   * @param cv Target canvas.
   * @param x,y Projected sample coordinates.
   * @param c,age,alpha Shading registers.
   * @details Writes the taps, coverage and order plot() passes on. */
  template <typename TerminalT>
  HS_HOT_INLINE static void splat(TerminalT &terminal, Canvas &cv, float x,
                                  float y, const ::Pixel &c, float age,
                                  float alpha) {
    assert(age >= 0.0f && alpha >= 0.0f);
    (void)age;
    terminal.plot_taps(cv, splat_taps<W, H>(x, y), c, alpha, SPLAT_TAP_CUTOFF);
  }
};
HS_O3_END

} // namespace Screen
} // namespace Filter
