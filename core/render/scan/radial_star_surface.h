/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file radial_star_surface.h
 * @brief Continuous twisted radial star family drawn as one surface.
 */

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>

#include "core/animation/orientation.h"
#include "core/color/baked_palette.h"
#include "core/memory.h"
#include "core/platform/attributes.h"
#include "core/platform/profiling.h"
#include "core/render/canvas.h"
#include "core/render/plot/contour_layout.h"
#include "core/render/scan/column_basis.h"
#include "core/render/scan/radial_surface_candidate_mask.h"
#include "core/render/scan/radial_surface_renderer.h"
#include "core/render/spacetime/surface.h"

namespace Scan {

/** @brief One frame's inputs to a RadialStarSurface. */
struct SurfaceFrame {
  float phase;      ///< Waveform phase at the frame's end, turns.
  float twist_step; ///< Waveform phase advance over the frame, turns.
  float amplitude;  ///< Peak contour rotation, radians.
  bool opposite;    ///< Toggles the rotation sense past the equator.
  float opacity;    ///< Frame opacity.
  const BakedPalette &palette;            ///< Palette sampled by radius.
  const math::Orientation<> &orientation; ///< Surface orientation.
};

/**
 * @brief Draws the continuum of sine-twisted screen-balanced planar stars as
 *        a first-arrival radial surface around both chart poles.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @tparam SIDES Star lobe count.
 * @details Holds no persistent storage; every table lives in scratch_arena_a
 * and scratch_arena_b for the duration of draw().
 */
template <int W, int H, int SIDES> class RadialStarSurface {
public:
  /** @brief Most radial knot intervals per family. */
  static constexpr int MAX_INTERVALS = 240;
  /** @brief Kernel family queried at each pixel. */
  using Family = RadialSurface<MAX_INTERVALS, 16>;

  /** @brief Peak scratch_arena_a bytes draw() takes.
   *  @return Byte count. */
  static constexpr size_t scratch_a_bytes() {
    return 2 * Family::storage_bytes(SIDES, MAX_INTERVALS) +
           W * sizeof(math::Vector) + 32;
  }

  /** @brief Peak scratch_arena_b bytes draw() takes.
   *  @return Byte count. */
  static constexpr size_t scratch_b_bytes() {
    return 2 * RadialSurfaceCandidateMask<>::STORAGE_BYTES +
           (MAX_INTERVALS + 1) * (sizeof(typename Family::Knot) +
                                  3 * sizeof(float) + sizeof(SurfaceColor)) +
           257 * sizeof(uint16_t) + 32;
  }

  /** @brief Draws the surface.
   * @tparam COUNT Whether to collect coverage statistics.
   * @tparam MASK Whether to use a conservative candidate mask.
   * @param canvas Destination canvas.
   * @param frame Phase, twist, opacity, palette and orientation.
   * @param timed Whether to evaluate geometry at each column time.
   * @param stats Optional coverage counters.
   * @param row_displacement Vertical scan offset in rows.
   */
  template <bool COUNT = false, bool MASK = true>
  HS_HOT_FLASH_MEMBER static void
  draw(Canvas &canvas, const SurfaceFrame &frame, bool timed,
       SurfaceAreaStats *stats = nullptr, float row_displacement = 0) {
    const int intervals =
        std::clamp(int(ceilf(fabsf(frame.amplitude) * 24)), 24, MAX_INTERVALS);
    ScratchScope scope_a(scratch_arena_a), scope_b(scratch_arena_b);
    const auto columns = timed ? prepare_column_basis<W, H>(
                                     scratch_arena_a, frame.orientation, true)
                               : ColumnBasis{};
    const auto basis = math::make_basis(frame.orientation.get(), math::X_AXIS);
    Family near, far;
    near.bind(scratch_arena_a, SIDES, Plot::STAR_INNER_RATIO, intervals);
    far.bind(scratch_arena_a, SIDES, Plot::STAR_INNER_RATIO, intervals);
    RadialSurfaceCandidateMask<> near_mask, far_mask;
    if constexpr (MASK) {
      near_mask.bind(scratch_arena_b);
      far_mask.bind(scratch_arena_b);
    }
    const float footprint =
        .25f * math::PI_F *
            sqrtf(math::RADIANS_PER_COLUMN<W> * math::RADIANS_PER_COLUMN<W> +
                  math::RADIANS_PER_ROW<H> * math::RADIANS_PER_ROW<H>) +
        1e-5f;
    auto *knots =
        scratch_arena_b.allocate_n<typename Family::Knot>(intervals + 1);
    auto *radii = scratch_arena_b.allocate_n<float>(intervals + 1);
    auto *sines = scratch_arena_b.allocate_n<float>(intervals + 1);
    auto *cosines = scratch_arena_b.allocate_n<float>(intervals + 1);
    auto *colors = scratch_arena_b.allocate_n<SurfaceColor>(intervals + 1);
    auto *first = scratch_arena_b.allocate_n<uint16_t>(257);
    const float opacity = frame.opacity;
    for (int i = 0; i <= intervals; ++i) {
      const float q = .5f * i / intervals;
      radii[i] = math::PI_F * Plot::screen_balanced_radius_t(q);
      sines[i] = sinf(math::TWO_PI_F * q);
      cosines[i] = cosf(math::TWO_PI_F * q);
      const auto color = frame.palette.get(1 - q);
      colors[i] = {color.color.r * color.alpha * opacity,
                   color.color.g * color.alpha * opacity,
                   color.color.b * color.alpha * opacity};
    }
    int index = 0;
    for (int i = 0; i <= 256; ++i) {
      const float radius = i * (.5f * math::PI_F / 256);
      while (index < intervals - 1 && radii[index + 1] < radius)
        ++index;
      first[i] = uint16_t(index);
    }
    float sine = 0, cosine = 1, current_time = 1;
    float near_mask_time = -2, far_mask_time = -2;
    constexpr float MASK_TIME_RANGE = 16.0f / W;
    const float phase_motion =
        timed ? fabsf(frame.amplitude * math::TWO_PI_F * frame.twist_step) *
                    MASK_TIME_RANGE
              : 0;
    const float complex_motion = .5f * math::PI_F * phase_motion;
    const float relative_phase_motion = phase_motion * math::PI_F / intervals;
    auto prepare = [&](int column) HS_HOT_INLINE {
      const float time =
          timed ? Spacetime::SamplingSurface<W>::time(column) : 1;
      current_time = time;
      const float angle =
          math::TWO_PI_F * (frame.phase - frame.twist_step * (1 - time));
      sine = sinf(angle);
      cosine = cosf(angle);
    };
    auto ensure = [&](int side) HS_HOT_FLASH_MEMBER {
      HS_PROFILE(ss_surface_knots);
      for (int i = 0; i <= intervals; ++i) {
        const float wave = side ? sine * cosines[i] - cosine * sines[i]
                                : sine * cosines[i] + cosine * sines[i];
        const float turn =
            (side && !frame.opposite ? -1 : 1) * frame.amplitude * wave;
        float angle = turn - floorf(turn / math::TWO_PI_F) * math::TWO_PI_F;
        const float sign = angle > math::PI_F ? -1.0f : 1.0f;
        if (angle > math::PI_F)
          angle -= math::PI_F;
        float turn_sine, turn_cosine;
        math::poly_sincosf_0_pi(angle, turn_sine, turn_cosine);
        knots[i] = {radii[i] * turn_cosine * sign, radii[i] * turn_sine * sign,
                    turn};
      }
      (side ? far : near)
          .set_knots(
              std::span<const typename Family::Knot>(knots, intervals + 1));
      if constexpr (MASK) {
        float &mask_time = side ? far_mask_time : near_mask_time;
        if (fabsf(current_time - mask_time) > MASK_TIME_RANGE) {
          HS_PROFILE(ss_surface_mask);
          const auto work = (side ? far_mask : near_mask)
                                .template prepare<COUNT>(
                                    side ? far : near, footprint,
                                    complex_motion, relative_phase_motion);
          mask_time = current_time;
          if constexpr (COUNT)
            if (stats) {
              ++stats->mask_prepares;
              stats->mask_candidates += work.candidates;
              stats->mask_pieces += work.pieces;
              stats->mask_row_writes += work.row_writes;
            }
          (void)work;
        }
      }
    };
    auto lower_bound = [&](float radius) HS_HOT_INLINE {
      return int(first[std::min(256, int(radius * (512 / math::PI_F)))]);
    };
    RadialSurfaceArea<Family, COUNT, MASK> area(near, near_mask, far_mask);
    auto shade = [&](const typename Family::Hit &hit) HS_HOT_INLINE {
      const auto a = colors[hit.interval], b = colors[hit.interval + 1];
      const float t = hit.fraction;
      return SurfaceColor{a.r + (b.r - a.r) * t, a.g + (b.g - a.g) * t,
                          a.b + (b.b - a.b) * t};
    };
    draw_radial_surface<W, H>(canvas, near, far, area, basis, columns, prepare,
                              ensure, lower_bound, shade, stats,
                              row_displacement);
  }
};

} // namespace Scan
