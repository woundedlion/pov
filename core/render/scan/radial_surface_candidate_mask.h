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

  /** @brief Mask cells per hemisphere. */
  static constexpr size_t STORAGE_COUNT = RADIAL_BINS * ANGULAR_BINS;
  /** @brief Bytes of mask cells per hemisphere. */
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
    const float radius = sqrtf(radius2);
    padding = footprint_radius + complex_motion +
              64 * FLT_EPSILON * (1 + radius + complex_motion);
    float derivative_motion = 0, relative_derivative_motion = 0;
    for (int edge = 0; edge < 2; ++edge) {
      const auto &m = surface.edge_metrics(edge);
      const float k = fmaxf(fabsf(m.dot), fabsf(m.dot + m.squared_length));
      const float error = adjacent_phase_motion >= 0
                              ? adjacent_phase_motion * (k + m.determinant)
                              : (2 * radius * complex_motion +
                                 complex_motion * complex_motion) *
                                    (k + 2 * m.determinant);
      if (adjacent_phase_motion >= 0)
        relative_derivative_motion = fmaxf(relative_derivative_motion, error);
      else
        derivative_motion = fmaxf(derivative_motion, error);
    }
    if (derivative_motion > 0)
      derivative_motion += 64 * FLT_EPSILON * (1 + derivative_motion);
    if (relative_derivative_motion > 0)
      relative_derivative_motion +=
          64 * FLT_EPSILON * (1 + relative_derivative_motion);
    extent2 = (radius + padding) * (radius + padding);
    radial_scale = RADIAL_BINS / extent2;
    angular_scale = ANGULAR_BINS * (surface.edge_count() / 2) / math::TWO_PI_F;
    inverse_step = 1 / fmaxf(2 * padding, radius / RADIAL_BINS);
    Stats stats;
    Boundaries::visit_candidates(
        surface,
        [&](const auto &candidate) HS_HOT_FLASH_MEMBER {
          if constexpr (COUNT)
            ++stats.candidates;
          source_bit = uint32_t(1)
                       << (candidate.origin_interval >> group_shift);
          const Point a{candidate.a.x, candidate.a.y};
          const Point b{candidate.b.x, candidate.b.y};
          const Point d{b.x - a.x, b.y - a.y};
          const float length = fabsf(d.x) + fabsf(d.y);
          const bool tip = !candidate.knot &&
                           candidate.origin_interval < surface.interval_count();
          const int pieces =
              tip ? 1
                  : std::clamp(int(length * inverse_step) + 1, 1, MAX_PIECES);
          if (pieces == 1) {
            mark<COUNT>(a, b, stats, candidate.angles_known || tip,
                        candidate.start_angle, candidate.finish_angle);
          } else {
            Point previous = a;
            const float inv = 1.0f / pieces;
            for (int i = 1; i <= pieces; ++i) {
              const float t = i * inv;
              const Point next{i == pieces ? b.x : a.x + t * d.x,
                               i == pieces ? b.y : a.y + t * d.y};
              mark<COUNT>(previous, next, stats, false, 0, 0);
              previous = next;
            }
          }
        },
        derivative_motion, relative_derivative_motion);
    return stats;
  }

  /** @brief Checks a center against the prepared conservative mask.
   * @param rho Chart radius, nonnegative.
   * @param theta Chart argument, radians.
   * @return Whether edge AA may be necessary. */
  HS_HOT_INLINE bool near_boundary(float rho, float theta) const {
    return interval_bits(rho, theta) != 0;
  }

  /** @brief Returns candidate interval groups for a pixel center.
   * @param rho Chart radius.
   * @param theta Chart argument.
   * @return Source interval bits, including the terminal endpoint. */
  HS_HOT_INLINE uint32_t interval_bits(float rho, float theta) const {
    const float r2 = rho * rho;
    if (r2 >= extent2)
      return 0;
    const int radial = std::clamp(int(r2 * radial_scale), 0, RADIAL_BINS - 1);
    const int angle = radial_wrap_index(radial_floor_int(theta * angular_scale),
                                        ANGULAR_BINS);
    return rows[radial * ANGULAR_BINS + angle];
  }

  /** @return Shift mapping a source interval to its candidate bit. */
  HS_HOT_INLINE int interval_shift() const { return group_shift; }

private:
  /** @brief Squared planar distance from the origin. */
  static HS_HOT_INLINE float square(Point p) { return p.x * p.x + p.y * p.y; }
  /** @brief Marks candidate bins touched by a planar boundary segment. */
  template <bool COUNT>
  HS_HOT_INLINE void mark(Point a, Point b, Stats &stats, bool angles_known,
                          float start_angle, float finish_angle) {
    if constexpr (COUNT)
      ++stats.pieces;
    const Point d{b.x - a.x, b.y - a.y};
    const float length2 = square(d);
    const float t =
        length2 > 0 ? std::clamp(-(a.x * d.x + a.y * d.y) / length2, 0.0f, 1.0f)
                    : 0;
    const Point nearest{a.x + t * d.x, a.y + t * d.y};
    const float min2 = square(nearest), max2 = fmaxf(square(a), square(b));
    const float upper = sqrtf(max2);
    const float pad2 = 2 * upper * padding + padding * padding;
    const int first_row = std::clamp(
        radial_floor_int((min2 - pad2) * radial_scale), 0, RADIAL_BINS - 1);
    const int last_row = std::clamp(
        radial_floor_int((max2 + pad2) * radial_scale), 0, RADIAL_BINS - 1);
    constexpr Bits ALL = ~Bits(0) >> (sizeof(Bits) * 8 - ANGULAR_BINS);
    float first = start_angle, last = finish_angle;
    const float padding2 = padding * padding;
    const bool bounded =
        min2 > padding2 + 32 * FLT_EPSILON * fmaxf(min2, padding2);
    if (bounded && !angles_known) {
      first = radial_atan2(a.y, a.x);
      float delta = radial_atan2(b.y, b.x) - first;
      if (delta > math::PI_F)
        delta -= math::TWO_PI_F;
      if (delta < -math::PI_F)
        delta += math::TWO_PI_F;
      last = first + delta;
    }
    const float pad =
        bounded ? padding / sqrtf(min2 - padding * padding) + 8 * FLT_EPSILON
                : 0;
    const float radial_delta = square(b) - square(a);
    const float angular_rate =
        max2 > 0 ? fabsf(a.x * d.y - a.y * d.x) / max2 : 0;
    for (int row = first_row; row <= last_row; ++row) {
      Bits bits = ALL;
      if (bounded) {
        float low = 0, high = 1;
        const float lower_norm = row / radial_scale - pad2 - square(a);
        const float upper_norm =
            (row + 1) / radial_scale + pad2 + .25f * length2 - square(a);
        if (radial_delta > 0) {
          low = fmaxf(low, lower_norm / radial_delta);
          high = fminf(high, upper_norm / radial_delta);
        } else if (radial_delta < 0) {
          low = fmaxf(low, upper_norm / radial_delta);
          high = fminf(high, lower_norm / radial_delta);
        }
        if (low > high)
          continue;
        const float begin = first <= last ? first + angular_rate * low
                                          : last + angular_rate * (1 - high);
        const float end = first <= last ? last - angular_rate * (1 - high)
                                        : first - angular_rate * low;
        const int first_value = radial_floor_int((begin - pad) * angular_scale);
        const int last_value = radial_floor_int((end + pad) * angular_scale);
        const int bit_count = last_value - first_value + 1;
        if (bit_count < ANGULAR_BINS) {
          const int shift = radial_wrap_index(first_value, ANGULAR_BINS);
          const Bits range = (Bits(1) << bit_count) - 1;
          bits = shift == 0
                     ? range
                     : ((range << shift) | (range >> (ANGULAR_BINS - shift))) &
                           ALL;
        }
      }
      Bits pending = bits;
      while (pending) {
        const int angle = sizeof(Bits) <= 4
                              ? __builtin_ctz(uint32_t(pending))
                              : __builtin_ctzll(uint64_t(pending));
        rows[row * ANGULAR_BINS + angle] |= source_bit;
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
