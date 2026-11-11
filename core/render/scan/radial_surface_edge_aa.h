/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file radial_surface_edge_aa.h
 * @brief Two-color area coverage at radial first-arrival discontinuities. */

#include <type_traits>
#include "core/render/scan/radial_surface_boundary.h"
#include "core/render/scan/radial_surface_math.h"
#include "core/render/scan/radial_surface_types.h"

namespace Scan {

/** @brief First-arrival sampling with one local discontinuity per pixel.
 * @tparam Surface Bounded similar-polygon radial surface. */
template <class Surface> class RadialSurfaceEdgeAA {
public:
  /** @brief Planar chart point. */
  using Point = typename Surface::Point;
  using Prepared =
      SurfaceEdgePrepared; ///< Optional chart and candidate bounds.
#if HS_ENABLE_TEST_HOOKS
  static constexpr int TIP_PATH = 1;      ///< Increasing-radius tip path.
  static constexpr int CORNER_PATH = 2;   ///< Paired corner path.
  static constexpr int OFFSET_PATH = 4;   ///< Offset boundary path.
  static constexpr int TERMINAL_PATH = 8; ///< Terminal boundary flag.
  /** @brief An AA result and its exact first-arrival shader. */
  struct Observation {
    const Surface &surface; ///< Prepared family.
    Point p;                ///< Footprint centre.
    Point dx;               ///< Column derivative.
    Point dy;               ///< Row derivative.
    SurfaceColor color;     ///< Blended output.
    int path;               ///< Tip/corner/offset bits and terminal flag.
    void *shader;           ///< Shader object for this sample.
    SurfaceColor (*shade)(void *,
                          const typename Surface::Hit &); ///< Shader adapter.
  };
  static inline void (*observer)(const Observation &, void *) =
      nullptr;                                 ///< Test-only result visitor.
  static inline void *observer_data = nullptr; ///< Visitor context.
#endif

  /** @brief Area of a unit square on the positive side of a line.
   * @param a,b,c Coefficients of a*x+b*y+c >= 0.
   * @return Exact straight-line coverage, apart from floating rounding. */
  static HS_HOT_FLASH_MEMBER float coverage(float a, float b, float c) {
    const float a_value = fabsf(a), b_value = fabsf(b),
                h = .5f * (a_value + b_value);
    if (c >= h)
      return 1;
    if (c <= -h)
      return 0;
    const float large = fmaxf(a_value, b_value),
                small = fminf(a_value, b_value);
    if (fabsf(c) <= .5f * (large - small) || small == 0)
      return .5f + c / large;
    const float t = h - fabsf(c);
    const float corner = t * t / (2 * a_value * b_value);
    return c > 0 ? 1 - corner : corner;
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
   * @param prepared Optional chart bounds and candidate interval lookup.
   * @return Pixel color with straight-line edge coverage. A paired corner uses
   * the product of its two edge coverages. */
  template <bool COUNT = true, class Shade>
  HS_HOT_FLASH_MEMBER SurfaceColor sample(const Surface &surface, Point p,
                                          Point dx, Point dy, float theta,
                                          int first, Shade &&shade,
                                          SurfaceAreaStats *stats = nullptr,
                                          const Prepared &prepared = {}) const {
    auto query = [&](Point q, float angle, int lower,
                     typename Surface::Stats &work) HS_HOT_INLINE {
      if constexpr (COUNT)
        if (stats)
          return surface.template sample<true, true, false>(q, angle, lower,
                                                            &work);
      return surface.template sample<false, true, false>(q, angle, lower);
    };
    typename Surface::Stats work;
    const auto hit =
        query(p, theta,
              prepared.point_first >= 0 ? prepared.point_first : first, work);
    Detector detector(surface, p, dx, dy, theta, prepared);
    uint32_t tail = 0;
    for (int i = first; i < surface.interval_count(); ++i) {
      if (!detector.any_from(i))
        break;
      if (hit.found() &&
          (i > hit.interval || (i == hit.interval && hit.fraction == 0)) &&
          detector.covered(surface, i))
        break;
      if (!detector.active(i)) {
        detector.last = i;
        continue;
      }
      if (i % Surface::BLOCK_INTERVALS == 0) {
        const int next = surface.skip_outside(detector.query, i);
        if (next != i) {
          i = next - 1;
          continue;
        }
      }
      detector.interval(surface, i);
      ++tail;
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
#if HS_ENABLE_TEST_HOOKS
    int sample_path =
        (!center_front && detector.tip && surface.strictly_increasing_radius())
            ? TIP_PATH
        : (detector.paired && !detector.terminal_back && !center_front)
            ? CORNER_PATH
            : OFFSET_PATH;
    if (detector.terminal_back)
      sample_path |= TERMINAL_PATH;
#endif
    SurfaceColor front, back;
    int back_interval = hit.interval;
    uint32_t extra_intervals = 0;
    if (!center_front && detector.tip && surface.strictly_increasing_radius()) {
      typename Surface::Hit birth_hit;
      birth_hit.interval =
          std::clamp(int(detector.birth), 0, surface.interval_count() - 1);
      birth_hit.fraction =
          hs::clamp(detector.birth - birth_hit.interval, 0.0f, 1.0f);
      front = shade_edge(shade, birth_hit);
      back = color;
    } else {
      Point other;
      if (detector.paired && !detector.terminal_back && !center_front) {
        if (!detector.corner(SIDE_OFFSET, SIDE_OFFSET, other))
          other = p;
      } else {
        const float u = copysignf(1.0f, detector.a),
                    v = copysignf(1.0f, detector.b);
        const float along = center_front ? -SIDE_OFFSET : SIDE_OFFSET;
        other = detector.side(detector.boundary, along * u, along * v);
      }
      typename Surface::Stats outside_work;
      const auto outside_hit = query(other, radial_atan2(other.y, other.x),
                                     std::max(0, first - 1), outside_work);
      const SurfaceColor outside =
          outside_hit.found() ? shade_edge(shade, outside_hit) : SurfaceColor{};
      if (center_front) {
        front = color;
        back = outside;
        back_interval = outside_hit.interval;
      } else {
        front = outside;
        back = color;
      }
      if constexpr (COUNT)
        extra_intervals += outside_work.intervals;
    }
    float back_weight = 1 - front_weight;
    if (detector.terminal_back &&
        (back_interval < 0 || back_interval + 1 == surface.interval_count())) {
      typename Surface::Stats terminal_work, background_work;
      Point outside;
      if (!detector.corner(-SIDE_OFFSET, -SIDE_OFFSET, outside))
        outside =
            detector.side(detector.second_boundary,
                          -SIDE_OFFSET * copysignf(1.0f, detector.second_a),
                          -SIDE_OFFSET * copysignf(1.0f, detector.second_b));
      const auto background_hit =
          query(outside, radial_atan2(outside.y, outside.x),
                std::max(0, first - 1), background_work);
      if (!background_hit.found()) {
        Point inside;
        if (!detector.corner(-SIDE_OFFSET, SIDE_OFFSET, inside))
          inside =
              detector.side(detector.second_boundary,
                            SIDE_OFFSET * copysignf(1.0f, detector.second_a),
                            SIDE_OFFSET * copysignf(1.0f, detector.second_b));
        const auto terminal_hit =
            query(inside, radial_atan2(inside.y, inside.x),
                  std::max(0, first - 1), terminal_work);
        if (terminal_hit.found()) {
          back = shade_edge(shade, terminal_hit);
          back_weight *=
              coverage(detector.second_a, detector.second_b, detector.second_c);
        }
      }
      if constexpr (COUNT)
        extra_intervals += terminal_work.intervals + background_work.intervals;
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
    const SurfaceColor result{front.r * front_weight + back.r * back_weight,
                              front.g * front_weight + back.g * back_weight,
                              front.b * front_weight + back.b * back_weight};
#if HS_ENABLE_TEST_HOOKS
    if (observer) {
      auto adapter = [](void *shader, const typename Surface::Hit &hit) {
        return (*static_cast<std::remove_reference_t<Shade> *>(shader))(hit);
      };
      observer({surface, p, dx, dy, result, sample_path,
                const_cast<void *>(static_cast<const void *>(&shade)), adapter},
               observer_data);
    }
#endif
    return result;
  }

private:
  /** Offset across a boundary, in pixel units. */
  static constexpr float SIDE_OFFSET = .01f;
  /** Maximum birth separation for one sheet, in source intervals. */
  static constexpr float SAME_SHEET_BIRTHS = .75f;
  /** Maximum birth separation for paired sheets, in source intervals. */
  static constexpr float PAIRED_SHEET_BIRTHS = .25f;
  /** Minimum squared sine between compatible boundary lines. */
  static constexpr float MIN_LINE_SINE_SQ = .0975f;

  /** @brief Shades one first arrival on an edge sample.
   * @tparam Shade First-arrival color callable.
   * @param shade First-arrival shader.
   * @param hit Contained edge first arrival.
   * @return Linear edge color. */
  template <class Shade>
  static HS_HOT_FLASH_MEMBER SurfaceColor
  shade_edge(Shade &shade, const typename Surface::Hit &hit) {
    return shade(hit);
  }
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
    const int last =
        std::min(int(detector.first_candidate), surface.interval_count() - 1);
    for (int i = hit.interval; i <= last; ++i) {
      auto candidate = hit;
      candidate.interval = i;
      if (i != hit.interval)
        candidate.fraction = 0;
      const float witness = detector.witness(surface, candidate);
      if (witness >= 0 && witness <= detector.first_candidate)
        return true;
    }
    return false;
  }

  static HS_HOT_INLINE float dot(Point a, Point b) {
    return a.x * b.x + a.y * b.y;
  }
  struct Detector {
    Point p, dx, dy, shape0, shape1, boundary{}, second_boundary{};
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
                           Point y, float angle, const Prepared &prepared)
        : p(point), dx(x), dy(y), theta(angle),
          candidate_bits(prepared.candidate_bits),
          candidate_shift(prepared.candidate_shift) {
      sides = surface.edge_count() / 2;
      period = math::TWO_PI_F / sides;
      inverse_sector = 2.0f * sides / math::TWO_PI_F;
      extent = (prepared.extent >= 0
                    ? prepared.extent
                    : .5f * (sqrtf(dot(dx, dx)) + sqrtf(dot(dy, dy)))) +
               1e-6f;
      const float r = prepared.radius >= 0 ? prepared.radius : sqrtf(dot(p, p));
      radius = r;
      angular = r > extent
                    ? fminf(math::PI_F, extent / (r - extent) + 8 * FLT_EPSILON)
                    : math::PI_F;
      all = angular >= period * .5f;
      const auto &m0 = surface.edge_metrics(0), &m1 = surface.edge_metrics(1);
      shape0 = {m0.dot, m0.squared_length};
      shape1 = {m1.dot, m1.squared_length};
      h0 = m0.determinant;
      h1 = m1.determinant;
      edge_reach = extent * surface.unit_edge_length();
      query = {p, theta, radius, extent, angular};
    }
    /** @brief Maps contained pixel coordinates to the affine chart.
     * @param u,v Pixel coordinates.
     * @return Point in the pixel footprint. */
    HS_HOT_INLINE Point point(float u, float v) const {
      u = hs::clamp(u, -.5f, .5f);
      v = hs::clamp(v, -.5f, .5f);
      return {p.x + u * dx.x + v * dy.x, p.y + u * dx.y + v * dy.y};
    }
    /** @brief Certifies a point strictly inside the affine footprint.
     * @param q Chart point.
     * @param determinant Nonzero footprint determinant.
     * @return Whether cross-product bounds certify containment. */
    HS_HOT_INLINE bool contains(Point q, float determinant) const {
      const Point offset{q.x - p.x, q.y - p.y};
      const float x = fabsf(dx.x) + fabsf(dx.y);
      const float y = fabsf(dy.x) + fabsf(dy.y);
      const float error =
          16 * FLT_EPSILON *
          (x * y + (fabsf(offset.x) + fabsf(offset.y)) * fmaxf(x, y));
      const float half = .5f * fabsf(determinant) - error;
      return fabsf(radial_cross(offset, dy)) < half &&
             fabsf(radial_cross(dx, offset)) < half;
    }
    /** @brief Offsets a boundary point within the affine footprint.
     * @param q Boundary point.
     * @param u,v Pixel-coordinate offset.
     * @return Contained side sample. */
    HS_HOT_FLASH_MEMBER Point side(Point q, float u, float v) const {
      const Point offset{q.x - p.x, q.y - p.y};
      const float det = radial_cross(dx, dy);
      if (det == 0) {
        const Point e = dot(dx, dx) > dot(dy, dy) ? dx : dy;
        const float length = dot(e, e);
        if (length == 0)
          return p;
        const float x = dot(dx, e), y = dot(dy, e);
        const float h = .5f * (fabsf(x) + fabsf(y));
        const float t =
            hs::clamp(dot(offset, e) + u * x + v * y, -h, h) / length;
        return {p.x + t * e.x, p.y + t * e.y};
      }
      const float u_value = radial_cross(offset, dy) / det + u;
      const float v_value = radial_cross(dx, offset) / det + v;
      if (fabsf(u_value) <= .5f && fabsf(v_value) <= .5f)
        return {q.x + u * dx.x + v * dy.x, q.y + u * dx.y + v * dy.y};
      return point(u_value, v_value);
    }
    /** @brief Samples prescribed sides of two intersecting footprint lines.
     * @param first_side,second_side Signed offsets in pixel coordinates.
     * @param result Contained sample when the lines meet inside the footprint.
     * @return Whether the line intersection lies in the footprint. */
    HS_HOT_FLASH_MEMBER bool corner(float first_side, float second_side,
                                    Point &result) const {
      const float det = a * second_b - b * second_a;
      if (det == 0)
        return false;
      const float u0 = (b * second_c - second_b * c) / det;
      const float v0 = (second_a * c - a * second_c) / det;
      if (fabsf(u0) > .5f || fabsf(v0) > .5f)
        return false;
      const float c1 = c - first_side * (fabsf(a) + fabsf(b));
      const float c2 =
          second_c - second_side * (fabsf(second_a) + fabsf(second_b));
      const float du = (b * c2 - second_b * c1) / det - u0;
      const float dv = (second_a * c1 - a * c2) / det - v0;
      float along = 1;
      if (du != 0)
        along = fminf(along, (.5f - copysignf(1.0f, du) * u0) / fabsf(du));
      if (dv != 0)
        along = fminf(along, (.5f - copysignf(1.0f, dv) * v0) / fabsf(dv));
      result = point(u0 + along * du, v0 + along * dv);
      return true;
    }
    /** @brief Whether the source interval has an active candidate bit. */
    HS_HOT_INLINE bool active(int index) const {
      return candidate_bits == UINT32_MAX ||
             ((candidate_bits >> (index >> candidate_shift)) & 1u) != 0;
    }
    /** @brief Whether any candidate bit remains from the source interval. */
    HS_HOT_INLINE bool any_from(int index) const {
      return candidate_bits == UINT32_MAX ||
             (candidate_bits >> (index >> candidate_shift)) != 0;
    }
    /** @brief Whether boundary lines cross inside the pixel at sufficient angle. */
    HS_HOT_INLINE bool compatible(float x, float y, float z, float u, float v,
                                  float w) const {
      const float det = x * v - y * u;
      if (det * det <= MIN_LINE_SINE_SQ * (x * x + y * y) * (u * u + v * v))
        return false;
      return 2 * fabsf(y * w - v * z) <= fabsf(det) &&
             2 * fabsf(u * z - x * w) <= fabsf(det);
    }
    /** @brief Retains nearest compatible boundary lines in the pixel footprint. */
    HS_HOT_FLASH_MEMBER void record(Point A, Point B, float first_birth,
                                    float last_birth, float front_sign = 1,
                                    bool is_tip = false,
                                    bool is_terminal = false) {
      const Point e{B.x - A.x, B.y - A.y}, offset{p.x - A.x, p.y - A.y};
      const float length = dot(e, e), projection = dot(offset, e);
      if (length < 1e-20f)
        return;
      if (projection < 0 && dot(offset, offset) > extent * extent)
        return;
      if (projection > length) {
        const Point q{p.x - B.x, p.y - B.y};
        if (dot(q, q) > extent * extent)
          return;
      }
      const float la = radial_cross(e, dx), lb = radial_cross(e, dy),
                  lc = radial_cross(e, offset);
      const float support = .5f * (fabsf(la) + fabsf(lb));
      if (support <= 0 || fabsf(lc) >= support)
        return;
      const float along = .5f * (fabsf(dot(e, dx)) + fabsf(dot(e, dy)));
      if (projection + along < 0 || projection - along > length)
        return;
      float entry = 0, leave = 1;
      const float det = radial_cross(dx, dy);
      if (det == 0) {
        const Point axis = dot(dx, dx) > dot(dy, dy) ? dx : dy;
        const float denominator = radial_cross(e, axis);
        if (denominator == 0)
          return;
        const float t = radial_cross(offset, axis) / denominator;
        if (t < 0 || t > 1)
          return;
        const Point q{t * e.x - offset.x, t * e.y - offset.y};
        const float half = .5f * (fabsf(dot(dx, axis)) + fabsf(dot(dy, axis)));
        if (fabsf(dot(q, axis)) > half)
          return;
        entry = leave = t;
      } else if (const float closest =
                     hs::clamp(projection / length, 0.0f, 1.0f);
                 !contains({A.x + closest * e.x, A.y + closest * e.y}, det)) {
        const float sign = copysignf(1.0f, det);
        const float half = .5f * fabsf(det);
        const float u = -sign * radial_cross(offset, dy),
                    v = -sign * radial_cross(dx, offset);
        const float du = sign * radial_cross(e, dy),
                    dv = sign * radial_cross(dx, e);
        if (!RadialSurfaceBoundaryCandidates::clip_parameter(u - half, du,
                                                             entry, leave) ||
            !RadialSurfaceBoundaryCandidates::clip_parameter(-u - half, -du,
                                                             entry, leave) ||
            !RadialSurfaceBoundaryCandidates::clip_parameter(v - half, dv,
                                                             entry, leave) ||
            !RadialSurfaceBoundaryCandidates::clip_parameter(-v - half, -dv,
                                                             entry, leave))
          return;
      }
      const float earliest = fminf(first_birth, last_birth);
      const float order_guard = 64 * FLT_EPSILON * (1 + fabsf(earliest));
      first_candidate = fminf(first_candidate, earliest - order_guard);
      float distance;
      if (projection <= 0)
        distance = dot(offset, offset);
      else if (projection >= length) {
        const Point q{p.x - B.x, p.y - B.y};
        distance = dot(q, q);
      } else
        distance = lc * lc / length;
      const float t = hs::clamp(projection / length, entry, leave);
      const float new_a = front_sign * la, new_b = front_sign * lb,
                  new_c = front_sign * lc;
      const float new_birth = first_birth + t * (last_birth - first_birth);
      if (is_terminal && found && new_birth > birth + SAME_SHEET_BIRTHS) {
        second_boundary = {A.x + t * e.x, A.y + t * e.y};
        second_a = new_a;
        second_b = new_b;
        second_c = new_c;
        paired = true;
        terminal_back = true;
        return;
      }
      const bool match = found &&
                         fabsf(new_birth - birth) <= SAME_SHEET_BIRTHS &&
                         compatible(a, b, c, new_a, new_b, new_c);
      if (found && distance >= best_distance) {
        if (match) {
          second_a = new_a;
          second_b = new_b;
          second_c = new_c;
          paired = true;
        }
        return;
      }
      if (match) {
        second_a = a;
        second_b = b;
        second_c = c;
        paired = true;
      } else if (!(paired && fabsf(new_birth - birth) <= PAIRED_SHEET_BIRTHS &&
                   compatible(new_a, new_b, new_c, second_a, second_b,
                              second_c)))
        paired = false;
      boundary = {A.x + t * e.x, A.y + t * e.y};
      a = new_a;
      b = new_b;
      c = new_c;
      birth = new_birth;
      tip = is_tip;
      best_distance = distance;
      found = true;
    }
    HS_HOT_FLASH_MEMBER void endpoint(const Surface &surface, int index) {
      if (!active(index))
        return;
      const auto &k = surface.knot(index);
      const Point c{k.real, k.imag};
      const int low_sector =
          all ? 0
              : radial_floor_int((theta - angular - k.phase) * inverse_sector);
      const int high_sector =
          all ? 2 * sides - 1
              : radial_floor_int((theta + angular - k.phase) * inverse_sector);
      for (int sector = low_sector; sector <= high_sector; ++sector) {
        const int i = radial_wrap_index(sector, 2 * sides);
        float low = 0, high = 1;
        if (index > 0) {
          const auto &previous = surface.knot(index - 1);
          const Point d{c.x - previous.real, c.y - previous.imag};
          const Point shape = (i & 1) ? shape1 : shape0;
          const float h = (i & 1) ? h1 : h0;
          const float spin = radial_cross(c, d);
          if (!RadialSurfaceBoundaryCandidates::endpoint_range(
                  spin, dot(c, d), shape.x, shape.y, h, 0.0f, low, high))
            continue;
        }
        const Point a = radial_multiply(c, surface.vertex(i)),
                    e = radial_multiply(c, surface.edge(i));
        record({a.x + low * e.x, a.y + low * e.y},
               {a.x + high * e.x, a.y + high * e.y}, float(index), float(index),
               1, false, index == surface.interval_count());
      }
    }
    HS_HOT_INLINE float bound(Point c) const {
      const float x = fabsf(c.x), y = fabsf(c.y);
      return edge_reach * (fmaxf(x, y) + .5f * fminf(x, y)) + 32 * FLT_EPSILON;
    }
    HS_HOT_INLINE void center_predicate(const Surface &surface, int i,
                                        int &sector, float &value) const {
      if (cached_knot == i) {
        sector = cached_sector;
        value = cached_value;
        return;
      }
      const auto &k = surface.knot(i);
      const Point c{k.real, k.imag};
      sector = radial_floor_int((theta - k.phase) * inverse_sector);
      value = radial_cross(p, radial_multiply(c, surface.edge(radial_wrap_index(
                                                     sector, 2 * sides)))) -
              h0 * dot(c, c);
    }
    HS_HOT_FLASH_MEMBER bool covered(const Surface &surface, int i) const {
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
      const int i = hit.interval;
      if (covered(surface, i + 1))
        return float(i + 1);
      const auto &ka = surface.knot(i), &kb = surface.knot(i + 1);
      const int start = radial_floor_int((theta - ka.phase) * inverse_sector);
      const int end = radial_floor_int((theta - kb.phase) * inverse_sector);
      if (start == end)
        return -1;
      const int corner = std::max(start, end);
      if (corner & 1)
        return -1;
      const Point v = surface.vertex(radial_wrap_index(corner, 2 * sides));
      const Point c{ka.real, ka.imag}, d{kb.real - ka.real, kb.imag - ka.imag};
      const float denominator = radial_cross(p, radial_multiply(d, v));
      if (denominator == 0)
        return -1;
      const float t = -radial_cross(p, radial_multiply(c, v)) / denominator;
      if (t < hit.fraction || t > 1)
        return -1;
      const typename Surface::Knot k{c.x + t * d.x, c.y + t * d.y,
                                     theta - corner / inverse_sector};
      return surface.contains_footprint(query, dx, dy, k) ? i + t : -1;
    }

    HS_HOT_INLINE void interval(const Surface &surface, int i) {
      int start_sector, end_sector;
      float start_value, end_value;
      center_predicate(surface, i, start_sector, start_value);
      center_predicate(surface, i + 1, end_sector, end_value);
      last = i;
      cached_knot = i + 1;
      cached_sector = end_sector;
      cached_value = end_value;
      const auto &ka = surface.knot(i), &kb = surface.knot(i + 1);
      const Point c{ka.real, ka.imag}, d{kb.real - ka.real, kb.imag - ka.imag};
      if (i == 0)
        endpoint(surface, 0);
      const float spin = radial_cross(c, d), radial = dot(c, d);
      const float reach = radius + extent, inner = fmaxf(0.0f, radius - extent);
      const float ra = dot(c, c);
      const Point cb{kb.real, kb.imag};
      const float rb = dot(cb, cb);
      const bool radial_tip =
          radial >= 0 ? (ra <= reach * reach && rb >= inner * inner) : true;
      if (radial_tip) {
        const int low_sector =
            all ? -1
                : radial_floor_int(
                      (theta - angular - fmaxf(ka.phase, kb.phase)) *
                      inverse_sector);
        const int high_sector =
            all ? 2 * sides - 1
                : radial_floor_int(
                      (theta + angular - fminf(ka.phase, kb.phase)) *
                      inverse_sector);
        int corner = low_sector + 1;
        if (corner & 1)
          ++corner;
        if (corner <= high_sector) {
          const float length = dot(d, d);
          const auto before = surface.edge(spin >= 0 ? 0 : 2 * sides - 1),
                     after = surface.edge(spin >= 0 ? 2 * sides - 1 : 0);
          float low = 0, high = 1;
          if (RadialSurfaceBoundaryCandidates::tip_range(spin, radial, length,
                                                         h0, before.x, after.x,
                                                         0.0f, low, high)) {
            const Point a{c.x + low * d.x, c.y + low * d.y},
                b{c.x + high * d.x, c.y + high * d.y};
            for (; corner <= high_sector; corner += 2) {
              const Point tip =
                  surface.vertex(radial_wrap_index(corner, 2 * sides));
              record(radial_multiply(a, tip), radial_multiply(b, tip), i + low,
                     i + high, spin >= 0 ? 1.0f : -1.0f, true);
            }
          }
        }
      }
      if (i == 0 || fabsf(start_value) > bound(c))
        return;
      const auto &previous = surface.knot(i - 1);
      const Point left{c.x - previous.real, c.y - previous.imag};
      const float ls = radial_cross(c, left), lr = dot(c, left);
      const int low_sector =
          all ? 0
              : radial_floor_int((theta - angular - ka.phase) * inverse_sector);
      const int high_sector =
          all ? 2 * sides - 1
              : radial_floor_int((theta + angular - ka.phase) * inverse_sector);
      for (int sector = low_sector; sector <= high_sector; ++sector) {
        const int index = radial_wrap_index(sector, 2 * sides),
                  edge = index & 1;
        const Point shape = edge == 0 ? shape0 : shape1;
        const float h = edge == 0 ? h0 : h1;
        float low = 0, high = 1;
        if (!RadialSurfaceBoundaryCandidates::knot_range(
                ls, lr, spin, radial, shape.x, shape.y, h, 0.0f, low, high))
          continue;
        const Point a = radial_multiply(c, surface.vertex(index)),
                    e = radial_multiply(c, surface.edge(index));
        record({a.x + low * e.x, a.y + low * e.y},
               {a.x + high * e.x, a.y + high * e.y}, float(i), float(i));
      }
    }
  };
};

} // namespace Scan
