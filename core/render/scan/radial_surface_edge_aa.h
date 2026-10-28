/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file radial_surface_edge_aa.h
 * @brief Two-color area coverage at radial first-arrival discontinuities. */

#include "core/render/scan/radial_surface_math.h"
#include "core/render/scan/radial_surface_types.h"

namespace Scan {

/** @brief First-arrival sampling with one local discontinuity per pixel.
 * @tparam Surface Bounded similar-polygon radial surface. */
template <class Surface> class RadialSurfaceEdgeAA {
public:
  /** @brief Planar chart point. */
  using Point = typename Surface::Point;

  /** @brief Area of a unit square on the positive side of a line.
   * @param a,b,c Coefficients of a*x+b*y+c >= 0.
   * @return Exact straight-line coverage, apart from floating rounding. */
  static HS_HOT_INLINE float coverage(float a, float b, float c) {
    const float A = fabsf(a), B = fabsf(b), H = .5f * (A + B);
    if (c >= H)
      return 1;
    if (c <= -H)
      return 0;
    const float LARGE = fmaxf(A, B), SMALL = fminf(A, B);
    if (fabsf(c) <= .5f * (LARGE - SMALL) || SMALL == 0)
      return .5f + c / LARGE;
    const float T = H - fabsf(c);
    const float CORNER = T * T / (2 * A * B);
    return c > 0 ? 1 - CORNER : CORNER;
  }

  /** @brief Shades a center hit and the opposite side of the nearest jump.
   * @tparam COUNT Enables optional work counters.
   * @tparam Shade Callable mapping a surface hit to linear color.
   * @param surface Prepared radial surface.
   * @param p,dx,dy Affine chart footprint.
   * @param theta Center argument.
   * @param first Conservative first radial interval for the footprint.
   * @param shade First-arrival color mapping.
   * @param stats Optional work counters.
   * @param prepared_radius,prepared_extent Optional precomputed chart bounds.
   * @param candidate_bits Conservative candidate groups, or all bits to disable.
   * @param candidate_shift Shift mapping source intervals, including the terminal
   * endpoint, to bit indices below 32. Ignored when candidate_bits is all bits.
   * @return Pixel color with straight-line area coverage. */
  template <bool COUNT = true, class Shade>
  HS_HOT_INLINE SurfaceColor
  sample(const Surface &surface, Point p, Point dx, Point dy, float theta,
         int first, Shade &&shade, SurfaceAreaStats *stats = nullptr,
         float prepared_radius = -1, float prepared_extent = -1,
         uint32_t candidate_bits = UINT32_MAX, int candidate_shift = 0) const {
    Detector detector(surface, p, dx, dy, theta, prepared_radius,
                      prepared_extent, candidate_bits, candidate_shift);
    typename Surface::Stats work;
    const auto hit = [&]() HS_HOT_INLINE {
      if constexpr (COUNT)
        if (stats)
          return surface.template sample<true, true, false>(p, theta, first,
                                                            &work, &detector);
      return surface.template sample<false, true, false>(p, theta, first,
                                                         nullptr, &detector);
    }();
    uint32_t tail = 0;
    if (hit.found()) {
      for (int i = std::max(hit.interval, detector.last + 1);
           i < surface.interval_count(); ++i) {
        if (!detector.any_from(i) || detector.covered(surface, i))
          break;
        if (!detector.active(i)) {
          detector.last = i;
          continue;
        }
        if (i % Surface::BLOCK_INTERVALS == 0) {
          const int NEXT = surface.skip_outside(detector.query, i);
          if (NEXT != i) {
            i = NEXT - 1;
            continue;
          }
        }
        detector.interval(surface, i);
        ++tail;
      }
    }
    if (!hit.found() || detector.last + 1 == surface.interval_count())
      detector.endpoint(surface, surface.interval_count());
    SurfaceColor color = hit.found() ? shade(hit) : SurfaceColor{};
    if constexpr (COUNT)
      if (stats) {
        ++stats->pixels;
        ++stats->fragments;
        stats->point_intervals += work.intervals;
        stats->intervals += work.intervals + tail;
        stats->max_point_intervals =
            std::max(stats->max_point_intervals, int(work.intervals));
        stats->max_intervals =
            std::max(stats->max_intervals, int(work.intervals + tail));
      }
    if (!detector.found)
      return color;
    if (hidden_only(surface, detector, hit)) {
      if constexpr (COUNT)
        if (stats)
          ++stats->hidden_pixels;
      return color;
    }
    float front_weight = coverage(detector.a, detector.b, detector.c);
    bool center_front = detector.c >= 0;
    if (detector.paired && !detector.terminal_back) {
      front_weight *=
          coverage(detector.second_a, detector.second_b, detector.second_c);
      center_front &= detector.second_c >= 0;
    }
    SurfaceColor front, back;
    int back_interval = hit.interval;
    uint32_t extra_intervals = 0;
    if (!center_front && detector.tip && surface.strictly_increasing_radius()) {
      typename Surface::Hit birth_hit;
      birth_hit.interval =
          std::clamp(int(detector.birth), 0, surface.interval_count() - 1);
      birth_hit.fraction =
          std::clamp(detector.birth - birth_hit.interval, 0.0f, 1.0f);
      front = shade(birth_hit);
      back = color;
    } else {
      Point other;
      if (detector.paired && !detector.terminal_back && !center_front) {
        const float DET =
            detector.a * detector.second_b - detector.b * detector.second_a;
        const float C1 =
            detector.c - .01f * (fabsf(detector.a) + fabsf(detector.b));
        const float C2 = detector.second_c - .01f * (fabsf(detector.second_a) +
                                                     fabsf(detector.second_b));
        const float U = (detector.b * C2 - detector.second_b * C1) / DET;
        const float V = (detector.second_a * C1 - detector.a * C2) / DET;
        other = {p.x + U * dx.x + V * dy.x, p.y + U * dx.y + V * dy.y};
      } else {
        const float U = copysignf(1.0f, detector.a),
                    V = copysignf(1.0f, detector.b);
        const float ALONG = center_front ? -.01f : .01f;
        other = {detector.boundary.x + ALONG * (U * dx.x + V * dy.x),
                 detector.boundary.y + ALONG * (U * dx.y + V * dy.y)};
      }
      typename Surface::Stats outside_work;
      const auto outside_hit = [&]() HS_HOT_INLINE {
        if constexpr (COUNT)
          if (stats)
            return surface.template sample<true, true, false>(
                other, radial_atan2(other.y, other.x), std::max(0, first - 1),
                &outside_work);
        return surface.template sample<false, true, false>(
            other, radial_atan2(other.y, other.x), std::max(0, first - 1));
      }();
      const SurfaceColor outside =
          outside_hit.found() ? shade(outside_hit) : SurfaceColor{};
      if (center_front) {
        front = color;
        back = outside;
        back_interval = outside_hit.interval;
      } else {
        front = outside;
        back = color;
      }
      if constexpr (COUNT)
        extra_intervals = outside_work.intervals;
    }
    float back_weight = 1 - front_weight;
    if (detector.terminal_back &&
        (back_interval < 0 || back_interval + 1 == surface.interval_count())) {
      typename Surface::Hit terminal_hit;
      terminal_hit.interval = surface.interval_count() - 1;
      terminal_hit.fraction = 1;
      back = shade(terminal_hit);
      back_weight *=
          coverage(detector.second_a, detector.second_b, detector.second_c);
    }
    if constexpr (COUNT)
      if (stats) {
        ++stats->edge_pixels;
        ++stats->fragments;
        stats->max_intervals = std::max(
            stats->max_intervals, int(work.intervals + tail + extra_intervals));
        stats->outside_intervals += extra_intervals;
        stats->intervals += extra_intervals;
        stats->max_parts = 2;
        stats->clips += detector.paired;
      }
    return {front.r * front_weight + back.r * back_weight,
            front.g * front_weight + back.g * back_weight,
            front.b * front_weight + back.b * back_weight};
  }

private:
  struct Detector;
  /** @brief Certifies all intersecting candidates behind one complete star.
   * @param surface Prepared radial family.
   * @param detector Intersecting candidates and affine footprint.
   * @param hit Center first arrival.
   * @return True only when a covering star precedes every candidate. */
  static HS_HOT_FLASH_MEMBER bool
  hidden_only(const Surface &surface, const Detector &detector,
              const typename Surface::Hit &hit) {
    if (!hit.found() || detector.first_candidate <= hit.interval + hit.fraction)
      return false;
    const int LAST =
        std::min(int(detector.first_candidate), surface.interval_count() - 1);
    for (int i = hit.interval; i <= LAST; ++i) {
      auto candidate = hit;
      candidate.interval = i;
      if (i != hit.interval)
        candidate.fraction = 0;
      const float WITNESS = detector.witness(surface, candidate);
      if (WITNESS >= 0 && WITNESS <= detector.first_candidate)
        return true;
    }
    return false;
  }

  static HS_HOT_INLINE float cross(Point a, Point b) {
    return a.x * b.y - a.y * b.x;
  }
  static HS_HOT_INLINE float dot(Point a, Point b) {
    return a.x * b.x + a.y * b.y;
  }
  static HS_HOT_INLINE Point multiply(Point a, Point b) {
    return {a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x};
  }
  static HS_HOT_INLINE int floor_int(float x) {
    const int T = int(x);
    return T - (float(T) > x);
  }
  static HS_HOT_INLINE int wrap(int x, int count) {
    const int R = x % count;
    return R < 0 ? R + count : R;
  }
  static HS_HOT_INLINE bool clip(float value, float slope, float &low,
                                 float &high) {
    if (slope > 0)
      high = fminf(high, -value / slope);
    else if (slope < 0)
      low = fmaxf(low, -value / slope);
    else if (value > 0)
      return false;
    return low <= high;
  }
  struct Detector {
    Point p, dx, dy, shape0, shape1, boundary{};
    typename Surface::QueryBounds query;
    float theta, extent, angular, period, inverse_sector, h0, h1, radius,
        edge_reach;
    float a = 0, b = 0, c = 0, best_distance = 0, birth = 0;
    float second_a = 0, second_b = 0, second_c = 0;
    bool found = false, all = false, tip = false, paired = false,
         terminal_back = false;
    int last = -1, sides, cached_knot = -1, cached_sector = 0;
    float cached_value = 0;
    float first_candidate =
        INFINITY; ///< Lower bound on intersecting candidate births.
    uint32_t candidate_bits;
    int candidate_shift;
    HS_HOT_INLINE Detector(const Surface &surface, Point point, Point x,
                           Point y, float angle, float prepared_radius = -1,
                           float prepared_extent = -1,
                           uint32_t bits = UINT32_MAX, int shift = 0)
        : p(point), dx(x), dy(y), theta(angle), candidate_bits(bits),
          candidate_shift(shift) {
      sides = surface.edge_count() / 2;
      period = math::TWO_PI_F / sides;
      inverse_sector = 2.0f * sides / math::TWO_PI_F;
      extent = (prepared_extent >= 0
                    ? prepared_extent
                    : .5f * (sqrtf(dot(dx, dx)) + sqrtf(dot(dy, dy)))) +
               1e-6f;
      const float R = prepared_radius >= 0 ? prepared_radius : sqrtf(dot(p, p));
      radius = R;
      angular = R > extent
                    ? fminf(math::PI_F, extent / (R - extent) + 8 * FLT_EPSILON)
                    : math::PI_F;
      all = angular >= period * .5f;
      const auto &M0 = surface.edge_metrics(0), &M1 = surface.edge_metrics(1);
      shape0 = {M0.dot, M0.squared_length};
      shape1 = {M1.dot, M1.squared_length};
      h0 = M0.determinant;
      h1 = M1.determinant;
      edge_reach = extent * surface.unit_edge_length();
      query = {p, theta, radius, extent, angular};
    }
    HS_HOT_INLINE bool active(int index) const {
      return candidate_bits == UINT32_MAX ||
             ((candidate_bits >> (index >> candidate_shift)) & 1u) != 0;
    }
    HS_HOT_INLINE bool any_from(int index) const {
      return candidate_bits == UINT32_MAX ||
             (candidate_bits >> (index >> candidate_shift)) != 0;
    }
    HS_HOT_INLINE bool compatible(float x, float y, float z, float u, float v,
                                  float w) const {
      const float DET = x * v - y * u;
      if (DET * DET <= .0975f * (x * x + y * y) * (u * u + v * v))
        return false;
      return 2 * fabsf(y * w - v * z) <= fabsf(DET) &&
             2 * fabsf(u * z - x * w) <= fabsf(DET);
    }
    HS_HOT_FLASH_MEMBER void record(Point A, Point B, float first_birth,
                                    float last_birth, float front_sign = 1,
                                    bool is_tip = false,
                                    bool is_terminal = false) {
      const Point E{B.x - A.x, B.y - A.y}, OFFSET{p.x - A.x, p.y - A.y};
      const float LENGTH = dot(E, E), PROJECTION = dot(OFFSET, E);
      if (LENGTH < 1e-20f)
        return;
      if (PROJECTION < 0 && dot(OFFSET, OFFSET) > extent * extent)
        return;
      if (PROJECTION > LENGTH) {
        const Point Q{p.x - B.x, p.y - B.y};
        if (dot(Q, Q) > extent * extent)
          return;
      }
      const float LA = cross(E, dx), LB = cross(E, dy), LC = cross(E, OFFSET);
      const float SUPPORT = .5f * (fabsf(LA) + fabsf(LB));
      if (SUPPORT <= 0 || fabsf(LC) >= SUPPORT)
        return;
      const float ALONG = .5f * (fabsf(dot(E, dx)) + fabsf(dot(E, dy)));
      if (PROJECTION + ALONG < 0 || PROJECTION - ALONG > LENGTH)
        return;
      if (PROJECTION - ALONG < 0 || PROJECTION + ALONG > LENGTH) {
        const float DET = cross(dx, dy), SIGN = copysignf(1.0f, DET);
        const float HALF = .5f * fabsf(DET);
        const float U = -SIGN * cross(OFFSET, dy),
                    V = -SIGN * cross(dx, OFFSET);
        const float DU = SIGN * cross(E, dy), DV = SIGN * cross(dx, E);
        float entry = 0, leave = 1;
        if (!clip(U - HALF, DU, entry, leave) ||
            !clip(-U - HALF, -DU, entry, leave) ||
            !clip(V - HALF, DV, entry, leave) ||
            !clip(-V - HALF, -DV, entry, leave))
          return;
      }
      const float EARLIEST = fminf(first_birth, last_birth);
      const float ORDER_GUARD = 64 * FLT_EPSILON * (1 + fabsf(EARLIEST));
      first_candidate = fminf(first_candidate, EARLIEST - ORDER_GUARD);
      float distance;
      if (PROJECTION <= 0)
        distance = dot(OFFSET, OFFSET);
      else if (PROJECTION >= LENGTH) {
        const Point Q{p.x - B.x, p.y - B.y};
        distance = dot(Q, Q);
      } else
        distance = LC * LC / LENGTH;
      const float T = std::clamp(PROJECTION / LENGTH, 0.0f, 1.0f);
      const float NEW_A = front_sign * LA, NEW_B = front_sign * LB,
                  NEW_C = front_sign * LC;
      const float NEW_BIRTH = first_birth + T * (last_birth - first_birth);
      if (is_terminal && found && NEW_BIRTH > birth + .75f) {
        second_a = NEW_A;
        second_b = NEW_B;
        second_c = NEW_C;
        paired = true;
        terminal_back = true;
        return;
      }
      const bool MATCH = found && fabsf(NEW_BIRTH - birth) <= .75f &&
                         compatible(a, b, c, NEW_A, NEW_B, NEW_C);
      if (found && distance >= best_distance) {
        if (MATCH) {
          second_a = NEW_A;
          second_b = NEW_B;
          second_c = NEW_C;
          paired = true;
        }
        return;
      }
      if (MATCH) {
        second_a = a;
        second_b = b;
        second_c = c;
        paired = true;
      } else if (!(paired && fabsf(NEW_BIRTH - birth) <= .25f &&
                   compatible(NEW_A, NEW_B, NEW_C, second_a, second_b,
                              second_c)))
        paired = false;
      boundary = {A.x + T * E.x, A.y + T * E.y};
      a = NEW_A;
      b = NEW_B;
      c = NEW_C;
      birth = NEW_BIRTH;
      tip = is_tip;
      best_distance = distance;
      found = true;
    }
    HS_HOT_INLINE void endpoint(const Surface &surface, int index) {
      if (!active(index))
        return;
      const auto &k = surface.knot(index);
      const Point C{k.real, k.imag};
      const int LOW_SECTOR =
          all ? 0 : floor_int((theta - angular - k.phase) * inverse_sector);
      const int HIGH_SECTOR =
          all ? 2 * sides - 1
              : floor_int((theta + angular - k.phase) * inverse_sector);
      for (int sector = LOW_SECTOR; sector <= HIGH_SECTOR; ++sector) {
        const int I = wrap(sector, 2 * sides);
        float low = 0, high = 1;
        if (index > 0) {
          const auto &previous = surface.knot(index - 1);
          const Point D{C.x - previous.real, C.y - previous.imag};
          const Point SHAPE = (I & 1) ? shape1 : shape0;
          const float H = (I & 1) ? h1 : h0;
          const float SPIN = cross(C, D);
          const float VALUE = SPIN * SHAPE.x - H * dot(C, D),
                      SLOPE = SPIN * SHAPE.y;
          const float TOL =
              64 * FLT_EPSILON * (fabsf(VALUE) + fabsf(SLOPE) + 1e-8f);
          if (!clip(VALUE - TOL, SLOPE, low, high))
            continue;
        }
        const Point A = multiply(C, surface.vertex(I)),
                    E = multiply(C, surface.edge(I));
        record({A.x + low * E.x, A.y + low * E.y},
               {A.x + high * E.x, A.y + high * E.y}, float(index), float(index),
               1, false, index == surface.interval_count());
      }
    }
    HS_HOT_INLINE float bound(Point c) const {
      const float X = fabsf(c.x), Y = fabsf(c.y);
      return edge_reach * (fmaxf(X, Y) + .5f * fminf(X, Y)) + 32 * FLT_EPSILON;
    }
    HS_HOT_INLINE void center_predicate(const Surface &surface, int i,
                                        int &sector, float &value) const {
      if (cached_knot == i) {
        sector = cached_sector;
        value = cached_value;
        return;
      }
      const auto &k = surface.knot(i);
      const Point C{k.real, k.imag};
      sector = floor_int((theta - k.phase) * inverse_sector);
      value = cross(p, multiply(C, surface.edge(wrap(sector, 2 * sides)))) -
              h0 * dot(C, C);
    }
    HS_HOT_INLINE bool covered(const Surface &surface, int i) const {
      int sector;
      float value;
      center_predicate(surface, i, sector, value);
      return surface.contains_footprint(query, dx, dy, surface.knot(i), sector,
                                        value);
    }
    /** @brief Finds a whole-footprint witness at an endpoint or outer-tip event.
     * @param surface Prepared radial family.
     * @param hit Earliest allowed parameter in the tested interval.
     * @return Witness parameter, or -1 when containment is uncertified. */
    HS_HOT_INLINE float witness(const Surface &surface,
                                const typename Surface::Hit &hit) const {
      if (!hit.found())
        return -1;
      const int I = hit.interval;
      if (covered(surface, I + 1))
        return float(I + 1);
      const auto &ka = surface.knot(I), &kb = surface.knot(I + 1);
      const int START = floor_int((theta - ka.phase) * inverse_sector);
      const int END = floor_int((theta - kb.phase) * inverse_sector);
      if (START == END)
        return -1;
      const int CORNER = std::max(START, END);
      if (CORNER & 1)
        return -1;
      const Point V = surface.vertex(wrap(CORNER, 2 * sides));
      const Point C{ka.real, ka.imag}, D{kb.real - ka.real, kb.imag - ka.imag};
      const float DENOMINATOR = cross(p, multiply(D, V));
      if (DENOMINATOR == 0)
        return -1;
      const float T = -cross(p, multiply(C, V)) / DENOMINATOR;
      if (T < hit.fraction || T > 1)
        return -1;
      const typename Surface::Knot K{C.x + T * D.x, C.y + T * D.y,
                                     theta - CORNER / inverse_sector};
      return surface.contains_footprint(query, dx, dy, K) ? I + T : -1;
    }

    HS_HOT_INLINE void interval(const Surface &surface, int i) {
      int start, end;
      float first, last_value;
      center_predicate(surface, i, start, first);
      center_predicate(surface, i + 1, end, last_value);
      interval(surface, i, start, first, end, last_value);
    }
    HS_HOT_INLINE
    void interval(const Surface &surface, int i, int start_sector,
                  float start_value, int end_sector, float end_value) {
      last = i;
      cached_knot = i + 1;
      cached_sector = end_sector;
      cached_value = end_value;
      if (!active(i))
        return;
      (void)start_sector;
      const auto &ka = surface.knot(i), &kb = surface.knot(i + 1);
      const Point C{ka.real, ka.imag}, D{kb.real - ka.real, kb.imag - ka.imag};
      if (i == 0)
        endpoint(surface, 0);
      const float SPIN = cross(C, D), RADIAL = dot(C, D);
      const float REACH = radius + extent, INNER = fmaxf(0.0f, radius - extent);
      const float RA = dot(C, C);
      const Point CB{kb.real, kb.imag};
      const float RB = dot(CB, CB);
      const bool RADIAL_TIP =
          RADIAL >= 0 ? (RA <= REACH * REACH && RB >= INNER * INNER) : true;
      if (RADIAL_TIP) {
        const int LOW_SECTOR =
            all ? -1
                : floor_int((theta - angular - fmaxf(ka.phase, kb.phase)) *
                            inverse_sector);
        const int HIGH_SECTOR =
            all ? 2 * sides - 1
                : floor_int((theta + angular - fminf(ka.phase, kb.phase)) *
                            inverse_sector);
        int corner = LOW_SECTOR + 1;
        if (corner & 1)
          ++corner;
        if (corner <= HIGH_SECTOR) {
          const float LENGTH = dot(D, D);
          const auto BEFORE = surface.edge(SPIN >= 0 ? 0 : 2 * sides - 1),
                     AFTER = surface.edge(SPIN >= 0 ? 2 * sides - 1 : 0);
          const Point L{SPIN * BEFORE.x - h0 * RADIAL, -h0 * LENGTH};
          const Point R{SPIN * AFTER.x - h0 * RADIAL, -h0 * LENGTH};
          const float TOL =
              64 * FLT_EPSILON *
              (fabsf(L.x) + fabsf(L.y) + fabsf(R.x) + fabsf(R.y) + 1e-8f);
          float low = 0, high = 1;
          if (fminf(L.x, L.x + L.y) <= TOL && fmaxf(R.x, R.x + R.y) >= -TOL &&
              clip(L.x - TOL, L.y, low, high) &&
              clip(-R.x - TOL, -R.y, low, high)) {
            const Point A{C.x + low * D.x, C.y + low * D.y},
                B{C.x + high * D.x, C.y + high * D.y};
            for (; corner <= HIGH_SECTOR; corner += 2) {
              const Point TIP = surface.vertex(wrap(corner, 2 * sides));
              record(multiply(A, TIP), multiply(B, TIP), i + low, i + high,
                     SPIN >= 0 ? 1.0f : -1.0f, true);
            }
          }
        }
      }
      if (i == 0 || fabsf(start_value) > bound(C))
        return;
      const auto &previous = surface.knot(i - 1);
      const Point LEFT{C.x - previous.real, C.y - previous.imag};
      const float LS = cross(C, LEFT), LR = dot(C, LEFT);
      const int LOW_SECTOR =
          all ? 0 : floor_int((theta - angular - ka.phase) * inverse_sector);
      const int HIGH_SECTOR =
          all ? 2 * sides - 1
              : floor_int((theta + angular - ka.phase) * inverse_sector);
      for (int sector = LOW_SECTOR; sector <= HIGH_SECTOR; ++sector) {
        const int INDEX = wrap(sector, 2 * sides), edge = INDEX & 1;
        const Point SHAPE = edge == 0 ? shape0 : shape1;
        const float H = edge == 0 ? h0 : h1;
        const Point L{LS * SHAPE.x - H * LR, LS * SHAPE.y};
        const Point R{SPIN * SHAPE.x - H * RADIAL, SPIN * SHAPE.y};
        const float TOL =
            64 * FLT_EPSILON *
            (fabsf(L.x) + fabsf(L.y) + fabsf(R.x) + fabsf(R.y) + 1e-8f);
        if (fminf(L.x, L.x + L.y) > TOL || fmaxf(R.x, R.x + R.y) < -TOL)
          continue;
        float low = 0, high = 1;
        if (!clip(L.x - TOL, L.y, low, high) ||
            !clip(-R.x - TOL, -R.y, low, high))
          continue;
        const Point A = multiply(C, surface.vertex(INDEX)),
                    E = multiply(C, surface.edge(INDEX));
        record({A.x + low * E.x, A.y + low * E.y},
               {A.x + high * E.x, A.y + high * E.y}, float(i), float(i));
      }
    }
  };
};

} // namespace Scan
