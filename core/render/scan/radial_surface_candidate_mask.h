/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file radial_surface_candidate_mask.h
 * @brief Conservative periodic candidate-interval masks for radial surfaces. */

#include <type_traits>
#include "core/memory.h"
#include "core/render/scan/radial_surface_boundary.h"
#include "core/render/scan/radial_surface_math.h"

namespace Scan {

/** @brief Conservative polar cells requiring first-arrival edge AA.
 * @tparam RADIAL_BINS Squared-radius cells.
 * @tparam ANGULAR_BINS Angular bits per repeated star period.
 * @tparam MAX_PIECES Maximum subdivisions per candidate segment. */
template <int RADIAL_BINS = 64, int ANGULAR_BINS = 32, int MAX_PIECES = 4>
class RadialSurfaceCandidateMask {
  static_assert(RADIAL_BINS > 0 && ANGULAR_BINS > 0 && ANGULAR_BINS <= 64);
  static_assert(MAX_PIECES > 0);
  using Bits = std::conditional_t<(ANGULAR_BINS <= 32), uint32_t, uint64_t>;
  using Storage = uint32_t;
  using Boundaries = RadialSurfaceBoundaryCandidates;
  using Point = typename Boundaries::Point;

public:
  /** @brief Preparation operation census. */
  struct Stats {
    uint32_t candidates = 0; ///< Canonical candidate segments.
    uint32_t pieces = 0;     ///< Padded segment pieces.
    uint32_t row_writes = 0; ///< Candidate-cell updates.
  };

  static constexpr size_t STORAGE_COUNT = RADIAL_BINS * ANGULAR_BINS;
  static constexpr size_t STORAGE_BYTES = sizeof(Storage) * STORAGE_COUNT;

  /** @brief Allocates reusable mask rows.
   * @param arena Storage owner. */
  HS_FLASH_MEMBER void bind(Arena &arena) {
    rows = arena.allocate_n<Storage>(STORAGE_COUNT);
  }

  /** @brief Marks every possible jump within the supplied pixel radius.
   * @tparam COUNT Enables preparation work counters.
   * @tparam Surface Prepared regular alternating-radius polygon family.
   * @param surface Similarity knots and unit geometry.
   * @param footprint_radius Maximum chart distance from center to pixel.
   * @param complex_motion Maximum displacement of every knot multiplier.
   * @param adjacent_phase_motion Optional adjacent phase-difference drift;
   * nonnegative values require fixed knot radii throughout the cache lifetime.
   * @return Bounded preparation work counts.
   * @details False cells certify that no family boundary candidate reaches
   * the footprint. Hidden boundaries remain marked. */
  template <bool COUNT = true, class Surface>
  HS_HOT_FLASH_MEMBER Stats prepare(const Surface &surface,
                                    float footprint_radius,
                                    float complex_motion = 0,
                                    float adjacent_phase_motion = -1) {
    HS_CHECK(rows != nullptr, "RadialSurfaceCandidateMask: unbound mask");
    HS_CHECK(footprint_radius >= 0,
             "RadialSurfaceCandidateMask: negative footprint");
    HS_CHECK(complex_motion >= 0,
             "RadialSurfaceCandidateMask: negative motion");
    std::fill_n(rows, STORAGE_COUNT, Storage(0));
    group_shift = 0;
    while ((surface.interval_count() >> group_shift) >= 32)
      ++group_shift;
    float radius2 = 0;
    for (int i = 0; i <= surface.interval_count(); ++i) {
      const auto &k = surface.knot(i);
      radius2 = fmaxf(radius2, k.real * k.real + k.imag * k.imag);
    }
    const float RADIUS = sqrtf(radius2);
    padding = footprint_radius + complex_motion +
              64 * FLT_EPSILON * (1 + RADIUS + complex_motion);
    float derivative_motion = 0;
    for (int edge = 0; edge < 2; ++edge) {
      const auto &m = surface.edge_metrics(edge);
      const float K = fmaxf(fabsf(m.dot), fabsf(m.dot + m.squared_length));
      const float ERROR =
          adjacent_phase_motion >= 0
              ? radius2 * adjacent_phase_motion * (K + m.determinant)
              : (2 * RADIUS * complex_motion +
                 complex_motion * complex_motion) *
                    (K + 2 * m.determinant);
      derivative_motion = fmaxf(derivative_motion, ERROR);
    }
    if (derivative_motion > 0)
      derivative_motion += 64 * FLT_EPSILON * (1 + derivative_motion);
    extent2 = (RADIUS + padding) * (RADIUS + padding);
    HS_CHECK(extent2 > 0, "RadialSurfaceCandidateMask: empty radius");
    radial_scale = RADIAL_BINS / extent2;
    angular_scale = ANGULAR_BINS * (surface.edge_count() / 2) / math::TWO_PI_F;
    inverse_step = 1 / fmaxf(2 * padding, RADIUS / RADIAL_BINS);
    Stats stats;
    Boundaries::visit_candidates_detailed(
        surface,
        [&](const auto &candidate) HS_HOT_FLASH_MEMBER {
          if constexpr (COUNT)
            ++stats.candidates;
          source_bit = uint32_t(1)
                       << (candidate.origin_interval >> group_shift);
          const Point A{candidate.a.x, candidate.a.y};
          const Point B{candidate.b.x, candidate.b.y};
          const Point D{B.x - A.x, B.y - A.y};
          const float LENGTH = fabsf(D.x) + fabsf(D.y);
          const bool TIP = !candidate.knot &&
                           candidate.origin_interval < surface.interval_count();
          const int PIECES =
              TIP ? 1
                  : std::clamp(int(LENGTH * inverse_step) + 1, 1, MAX_PIECES);
          if (PIECES == 1) {
            mark<COUNT>(A, B, stats, candidate.angles_known || TIP,
                        candidate.start_angle, candidate.finish_angle);
          } else {
            Point previous = A;
            const float INV = 1.0f / PIECES;
            for (int i = 1; i <= PIECES; ++i) {
              const float T = i * INV;
              const Point NEXT{i == PIECES ? B.x : A.x + T * D.x,
                               i == PIECES ? B.y : A.y + T * D.y};
              mark<COUNT>(previous, NEXT, stats, false, 0, 0);
              previous = NEXT;
            }
          }
        },
        derivative_motion);
    return stats;
  }

  /** @brief Checks a center against the prepared conservative mask.
   * @param rho Chart radius, nonnegative.
   * @param theta Chart argument, radians.
   * @return Whether edge AA may be necessary. */
  HS_HOT_INLINE bool near_boundary(float rho, float theta) const {
    const float R2 = rho * rho;
    if (R2 >= extent2)
      return false;
    const int RADIAL = std::clamp(int(R2 * radial_scale), 0, RADIAL_BINS - 1);
    const int ANGLE = wrap(floor_int(theta * angular_scale));
    return rows[RADIAL * ANGULAR_BINS + ANGLE] != 0;
  }

  /** @brief Returns candidate interval groups for a pixel center.
   * @param rho Chart radius.
   * @param theta Chart argument.
   * @return Source interval bits, including the terminal endpoint. */
  HS_HOT_INLINE uint32_t interval_bits(float rho, float theta) const {
    const float R2 = rho * rho;
    if (R2 >= extent2)
      return 0;
    const int RADIAL = std::clamp(int(R2 * radial_scale), 0, RADIAL_BINS - 1);
    const int ANGLE = wrap(floor_int(theta * angular_scale));
    return rows[RADIAL * ANGULAR_BINS + ANGLE];
  }

  /** @return Shift mapping a source interval to its candidate bit. */
  HS_HOT_INLINE int interval_shift() const { return group_shift; }

private:
  static HS_HOT_INLINE int floor_int(float value) {
    const int TRUNC = int(value);
    return TRUNC - (float(TRUNC) > value);
  }
  static HS_HOT_INLINE int wrap(int value) {
    const int REMAINDER = value % ANGULAR_BINS;
    return REMAINDER < 0 ? REMAINDER + ANGULAR_BINS : REMAINDER;
  }
  static HS_HOT_INLINE float square(Point p) { return p.x * p.x + p.y * p.y; }
  static HS_HOT_INLINE float coordinate_lower(float a, float b) {
    return (a <= 0 && b >= 0) || (b <= 0 && a >= 0) ? 0
                                                    : fminf(fabsf(a), fabsf(b));
  }
  template <bool COUNT>
  HS_HOT_INLINE void mark(Point a, Point b, Stats &stats, bool angles_known,
                          float start_angle, float finish_angle) {
    if constexpr (COUNT)
      ++stats.pieces;
    const Point D{b.x - a.x, b.y - a.y};
    const float LENGTH2 = square(D);
    const float T =
        LENGTH2 > 0 ? std::clamp(-(a.x * D.x + a.y * D.y) / LENGTH2, 0.0f, 1.0f)
                    : 0;
    const Point NEAREST{a.x + T * D.x, a.y + T * D.y};
    const float MIN2 = square(NEAREST), MAX2 = fmaxf(square(a), square(b));
    const float UPPER = fmaxf(fabsf(a.x) + fabsf(a.y), fabsf(b.x) + fabsf(b.y));
    const float PAD2 = 2 * UPPER * padding + padding * padding;
    const int FIRST_ROW =
        std::clamp(floor_int((MIN2 - PAD2) * radial_scale), 0, RADIAL_BINS - 1);
    const int LAST_ROW =
        std::clamp(floor_int((MAX2 + PAD2) * radial_scale), 0, RADIAL_BINS - 1);
    constexpr Bits ALL = ~Bits(0) >> (sizeof(Bits) * 8 - ANGULAR_BINS);
    Bits bits = ALL;
    const float LOWER =
        fmaxf(coordinate_lower(a.x, b.x), coordinate_lower(a.y, b.y));
    if (LOWER > padding) {
      float first = start_angle, last = finish_angle;
      if (!angles_known) {
        first = radial_atan2(a.y, a.x);
        float delta = radial_atan2(b.y, b.x) - first;
        if (delta > math::PI_F)
          delta -= math::TWO_PI_F;
        if (delta < -math::PI_F)
          delta += math::TWO_PI_F;
        last = first + delta;
      }
      const float PAD = padding / (LOWER - padding) + 8 * FLT_EPSILON;
      const int FIRST = floor_int((fminf(first, last) - PAD) * angular_scale);
      const int LAST = floor_int((fmaxf(first, last) + PAD) * angular_scale);
      const int BIT_COUNT = LAST - FIRST + 1;
      if (BIT_COUNT < ANGULAR_BINS) {
        const int SHIFT = wrap(FIRST);
        const Bits RANGE = (Bits(1) << BIT_COUNT) - 1;
        bits =
            SHIFT == 0
                ? RANGE
                : ((RANGE << SHIFT) | (RANGE >> (ANGULAR_BINS - SHIFT))) & ALL;
      }
    }
    for (int row = FIRST_ROW; row <= LAST_ROW; ++row) {
      Bits pending = bits;
      while (pending) {
        const int ANGLE = sizeof(Bits) <= 4
                              ? __builtin_ctz(uint32_t(pending))
                              : __builtin_ctzll(uint64_t(pending));
        rows[row * ANGULAR_BINS + ANGLE] |= source_bit;
        pending &= pending - 1;
        if constexpr (COUNT)
          ++stats.row_writes;
      }
    }
  }

  Storage *rows = nullptr;
  uint32_t source_bit = 0;
  int group_shift = 0;
  float padding = 0, extent2 = 0, radial_scale = 0, angular_scale = 0,
        inverse_step = 0;
};

} // namespace Scan
