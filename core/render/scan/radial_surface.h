/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file radial_surface.h
 * @brief First-arrival queries on linearly swept similar planar stars. */

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <span>
#include <type_traits>
#include "core/math/3dmath.h"
#include "core/render/scan/radial_surface_math.h"
#include "core/memory.h"

namespace Scan {

/** @brief An arena-backed sweep of scaled and rotated alternating-radius stars.
 * @tparam MAX_INTERVALS Fixed upper bound on a query's interval visits.
 * @tparam MAX_SIDES Fixed upper bound on the number of outer tips.
 * @tparam BLOCK_SIZE Intervals grouped by each conservative rejection bound. */
template <int MAX_INTERVALS, int MAX_SIDES, int BLOCK_SIZE = 2>
class RadialSurface {
  static_assert(MAX_INTERVALS > 0 && MAX_SIDES >= 3 && BLOCK_SIZE > 0);

public:
  /** @brief A vector in the planar chart. */
  struct Point {
    float x = 0.0f; ///< First chart coordinate.
    float y = 0.0f; ///< Second chart coordinate.
  };

  /** @brief Similarity multiplier at one sweep knot. */
  struct Knot {
    float real;  ///< Real component of radius times exp(i phase).
    float imag;  ///< Imaginary component of radius times exp(i phase).
    float phase; ///< Unwrapped argument of the multiplier, radians.
  };

  /** @brief Unit-profile edge constants. */
  struct EdgeMetrics {
    float dot;            ///< dot(A, E) for edge start A and edge vector E.
    float squared_length; ///< |E|^2.
    float determinant;    ///< radial_cross(A, E).
  };

  /** @brief First boundary sheet reached by a chart point. */
  struct Hit {
    int interval = -1;     ///< Interval index; -1 means uncovered.
    float fraction = 0.0f; ///< Fraction within the interval.
    int edge = 0;          ///< Unit-star edge containing the hit.
    Point normal{};        ///< Unnormalized outward edge normal.
    float derivative = 0;  ///< Edge predicate's interval derivative.
    Point gradient{};      ///< Derivative of interval fraction by chart point.

    /** @return Whether the point belongs to the swept solid. */
    bool found() const { return interval >= 0; }
  };

  /** @brief Optional query work counters. */
  struct Stats {
    uint32_t intervals = 0;     ///< Visited radial intervals.
    uint32_t events = 0;        ///< Evaluated sector-crossing events.
    uint32_t roots = 0;         ///< Closed-form quadratic or linear solves.
    uint32_t coarse_checks = 0; ///< Tested radial blocks.
    uint32_t skipped = 0;       ///< Fine intervals rejected by block bounds.
  };

  /** @brief A chart disk and its conservative angular extent. */
  struct QueryBounds {
    Point point;          ///< Disk center.
    float theta;          ///< Argument of the center.
    float radius;         ///< Distance of the center from the chart origin.
    float extent;         ///< Enclosing disk radius.
    float angular_extent; ///< Upper bound on the disk's argument variation.
  };

  /** @brief Intervals per certified skip block. */
  static constexpr int BLOCK_INTERVALS = BLOCK_SIZE;
  /** @brief Empty scan observer, removed at compile time. */
  struct NoObserver {};

  /** @brief Arena bound for one reusable surface.
   * @param sides Outer tip count.
   * @param intervals Sweep interval count.
   * @return Bytes including allocation alignment slack. */
  static constexpr size_t storage_bytes(int sides, int intervals) {
    return 4 * sides * sizeof(Point) +
           (intervals + 1) * (sizeof(Knot) + sizeof(float)) +
           intervals * sizeof(Piece) +
           ((intervals + BLOCK_SIZE - 1) / BLOCK_SIZE) * sizeof(Block) + 24;
  }

  /** @brief Allocates reusable sweep storage and prepares the unit star.
   * @param arena Storage owner.
   * @param sides Number of outer tips, at least 3.
   * @param ratio Inner-to-outer vertex radius ratio, in (0,cos(pi/sides)].
   * @param intervals Number of sweep intervals. */
  HS_FLASH_MEMBER void bind(Arena &arena, int sides, float ratio,
                            int intervals) {
    HS_CHECK(sides >= 3 && sides <= MAX_SIDES, "RadialSurface: side capacity");
    HS_CHECK(intervals >= 1 && intervals <= MAX_INTERVALS,
             "RadialSurface: interval capacity");
    HS_CHECK(ratio > 0.0f && ratio <= cosf(math::PI_F / sides),
             "RadialSurface: inner radius ratio");
    num_edges = 2 * sides;
    num_intervals = intervals;
    inner = ratio;
    sector = math::PI_F / sides;
    inverse_sector = sides / math::PI_F;
    determinant = ratio * sinf(sector);
    vertices = arena.allocate_n<Point>(size_t(num_edges));
    edges = arena.allocate_n<Point>(size_t(num_edges));
    knots = arena.allocate_n<Knot>(size_t(intervals + 1));
    norms = arena.allocate_n<float>(size_t(intervals + 1));
    pieces = arena.allocate_n<Piece>(size_t(intervals));
    blocks = arena.allocate_n<Block>(
        size_t((intervals + BLOCK_INTERVALS - 1) / BLOCK_INTERVALS));
    edge_length = 0.0f;
    for (int i = 0; i < num_edges; ++i) {
      const float radius = (i & 1) ? inner : 1.0f;
      vertices[i] = {radius * cosf(i * sector), radius * sinf(i * sector)};
    }
    for (int i = 0; i < num_edges; ++i) {
      const Point next = vertices[(i + 1) % num_edges];
      edges[i] = {next.x - vertices[i].x, next.y - vertices[i].y};
      edge_length = fmaxf(edge_length, sqrtf(edges[i].x * edges[i].x +
                                             edges[i].y * edges[i].y));
    }
    auto metrics = [&](int i) {
      const Point a = vertices[i], e = edges[i];
      return EdgeMetrics{a.x * e.x + a.y * e.y, e.x * e.x + e.y * e.y,
                         a.x * e.y - a.y * e.x};
    };
    first_edge_metrics = metrics(0);
    second_edge_metrics = metrics(1);
  }

  /** @brief Replaces the sweep without allocating.
   * @param source Similarity knots; radius is positive after the first knot.
   * @details Adjacent arguments differ by less than one edge sector. Arguments
   * must match the nonzero multipliers, with a continuous unwrapping. */
  HS_FLASH_MEMBER void set_knots(std::span<const Knot> source) {
    HS_CHECK(source.size() == size_t(num_intervals + 1),
             "RadialSurface: knot count");
    increasing_radius = true;
    for (int i = 0; i <= num_intervals; ++i) {
      knots[i] = source[size_t(i)];
      norms[i] = knots[i].real * knots[i].real + knots[i].imag * knots[i].imag;
      HS_CHECK(i == 0 || norms[i] > 0.0f,
               "RadialSurface: noninitial zero multiplier");
      if (i != 0) {
        HS_CHECK(fabsf(knots[i].phase - knots[i - 1].phase) < sector,
                 "RadialSurface: interval crosses multiple sectors");
        const float dx = knots[i].real - knots[i - 1].real;
        const float dy = knots[i].imag - knots[i - 1].imag;
        pieces[i - 1] = {{dx, dy},
                         determinant * (dx * dx + dy * dy),
                         2.0f * determinant *
                             (knots[i - 1].real * dx + knots[i - 1].imag * dy)};
        increasing_radius &=
            pieces[i - 1].radial_linear >= 0 && pieces[i - 1].quadratic > 0;
      }
    }
    for (int first = 0; first < num_intervals; first += BLOCK_INTERVALS) {
      const int end = std::min(first + BLOCK_INTERVALS, num_intervals);
      const float inv = 1.0f / (end - first);
      const Knot &start = knots[first], &finish = knots[end];
      float low = start.phase, high = start.phase, deviation2 = 0;
      float radius2 = norms[first];
      for (int i = first + 1; i <= end; ++i) {
        low = fminf(low, knots[i].phase);
        high = fmaxf(high, knots[i].phase);
        radius2 = fmaxf(radius2, norms[i]);
        const float t = (i - first) * inv;
        const float dx =
            knots[i].real - (start.real + t * (finish.real - start.real));
        const float dy =
            knots[i].imag - (start.imag + t * (finish.imag - start.imag));
        deviation2 = fmaxf(deviation2, dx * dx + dy * dy);
      }
      const float radius = sqrtf(radius2);
      const float coarse_radius = sqrtf(fmaxf(norms[first], norms[end]));
      const float deviation = sqrtf(deviation2) + 16.0f * FLT_EPSILON * radius;
      blocks[first / BLOCK_INTERVALS] = {low,
                                         high,
                                         radius * (1.0f + 16.0f * FLT_EPSILON),
                                         edge_length * deviation,
                                         determinant * deviation *
                                             (2.0f * coarse_radius + deviation),
                                         edge_length * coarse_radius};
    }
  }

  /** @brief Allocates and prepares a complete sweep.
   * @param arena Storage owner.
   * @param sides Number of outer tips.
   * @param ratio Inner-to-outer radius ratio.
   * @param source Similarity knots. */
  HS_FLASH_MEMBER void prepare(Arena &arena, int sides, float ratio,
                               std::span<const Knot> source) {
    bind(arena, sides, ratio, int(source.size()) - 1);
    set_knots(source);
  }

  /** @return Number of prepared radial intervals. */
  int interval_count() const { return num_intervals; }
  /** @return Number of edges in each star. */
  int edge_count() const { return num_edges; }
  /** @return Whether multiplier radius increases strictly throughout the sweep. */
  bool strictly_increasing_radius() const { return increasing_radius; }
  /** @return Maximum unit edge length. */
  float unit_edge_length() const { return edge_length; }
  /** @param parity Alternating edge parity.
   * @return Cached unit-profile edge constants. */
  const EdgeMetrics &edge_metrics(int parity) const {
    return parity ? second_edge_metrics : first_edge_metrics;
  }
  /** @return Inner-to-outer unit-star vertex radius ratio. */
  float inner_ratio() const { return inner; }
  /** @param i Knot index.
   * @return Similarity multiplier at the knot. */
  const Knot &knot(int i) const { return knots[i]; }
  /** @param i Unit-star edge index.
   * @return Edge vector from vertex i to vertex i+1. */
  Point edge(int i) const { return edges[i]; }
  /** @param i Unit-star vertex index.
   * @return Alternating-radius unit-star vertex. */
  Point vertex(int i) const { return vertices[i]; }

  /** @brief Prepares a disk bound shared by consecutive block checks.
   * @param p Chart center.
   * @param theta Argument of p.
   * @param extent Enclosing footprint radius.
   * @return Radius and angular extent for conservative block rejection. */
  static HS_HOT_INLINE QueryBounds query_bounds(Point p, float theta,
                                                float extent) {
    const float radius = sqrtf(p.x * p.x + p.y * p.y);
    const float angle =
        extent == 0.0f
            ? 0.0f
            : (radius > extent ? fminf(math::PI_F, extent / (radius - extent))
                               : math::PI_F);
    return {p, theta, radius, extent, angle};
  }

  /** @brief Skips the remainder of a block only when its whole sweep misses.
   * @param query Chart disk with conservative radius and angular bounds.
   * @param interval Current fine interval.
   * @return Block end when certified outside; otherwise interval. */
  HS_HOT_INLINE int skip_outside(const QueryBounds &query, int interval) const {
    if (interval < 0 || interval >= num_intervals)
      return interval;
    const int first = interval / BLOCK_INTERVALS * BLOCK_INTERVALS;
    const int end = std::min(first + BLOCK_INTERVALS, num_intervals);
    const Block &block = blocks[first / BLOCK_INTERVALS];
    if (query.radius - query.extent > block.radius)
      return end;
    const int first_sector =
        sector_index(query.theta - query.angular_extent - block.high_phase);
    const int last_sector =
        sector_index(query.theta + query.angular_extent - block.low_phase);
    if (first_sector != last_sector)
      return interval;
    const Point edge = edges[radial_wrap_index(first_sector, num_edges)];
    const float a = radial_cross(query.point, multiply(knots[first], edge)) -
                    determinant * norms[first];
    const float b = radial_cross(query.point, multiply(knots[end], edge)) -
                    determinant * norms[end];
    const float rounding = 32.0f * FLT_EPSILON *
                           (1.0f + query.radius * block.reach +
                            determinant * block.radius * block.radius);
    const float margin = (query.radius + query.extent) * block.point_error +
                         block.constant_error + query.extent * block.reach +
                         rounding;
    return fminf(a, b) > margin ? end : interval;
  }

  /** @brief Tests a chart footprint against a prepared radial block.
   * @param p Chart center.
   * @param footprint_radius Enclosing footprint radius.
   * @param interval Current fine interval.
   * @return Certified block end, or interval. */
  HS_HOT_INLINE int skip_outside(Point p, float footprint_radius,
                                 int interval) const {
    return skip_outside(query_bounds(p, atan2f(p.y, p.x), footprint_radius),
                        interval);
  }

  /** @brief Certifies a complete affine pixel footprint inside one star.
   * @param p Footprint center.
   * @param dx,dy Full pixel-coordinate derivatives.
   * @param theta Argument of the footprint center.
   * @param knot_index Star knot to test.
   * @return True only when all touched sectors certify containment. */
  HS_HOT_INLINE bool contains_footprint_at_knot(Point p, Point dx, Point dy,
                                                float theta,
                                                int knot_index) const {
    if (knot_index < 0 || knot_index > num_intervals || norms[knot_index] == 0)
      return false;
    const float extent = .5f * (sqrtf(dx.x * dx.x + dx.y * dx.y) +
                                sqrtf(dy.x * dy.x + dy.y * dy.y)) +
                         1e-6f;
    const auto query = query_bounds(p, theta, extent);
    return contains_footprint(query, dx, dy, knots[knot_index]);
  }

  /** @brief Certifies an affine footprint inside a similarity star.
   * @param query Conservative spatial and angular footprint bounds.
   * @param dx,dy Full pixel-coordinate derivatives.
   * @param c Similarity multiplier and its unwrapped argument.
   * @return True only if every touched edge half-plane contains the footprint. */
  HS_HOT_INLINE bool contains_footprint(const QueryBounds &query, Point dx,
                                        Point dy, const Knot &c) const {
    const int center = sector_index(query.theta - c.phase);
    const Point edge = multiply(c, edges[radial_wrap_index(center, num_edges)]);
    const float value = radial_cross(query.point, edge) -
                        determinant * (c.real * c.real + c.imag * c.imag);
    return contains_footprint(query, dx, dy, c, center, value);
  }

  /** @brief Certifies containment using a cached center edge predicate.
   * @param query Conservative footprint bounds.
   * @param dx,dy Full pixel-coordinate derivatives.
   * @param c Star similarity multiplier.
   * @param center Center angular sector.
   * @param value Center containment predicate in that sector.
   * @return True only when all touched sectors certify containment. */
  HS_HOT_INLINE bool contains_footprint(const QueryBounds &query, Point dx,
                                        Point dy, const Knot &c, int center,
                                        float value) const {
    if (value >= 0)
      return false;
    const Point c_value{c.real, c.imag};
    const float constant =
        determinant * (c_value.x * c_value.x + c_value.y * c_value.y);
    auto inside = [&](int index, float predicate) HS_HOT_INLINE {
      const Point edge =
          radial_multiply(c_value, edges[radial_wrap_index(index, num_edges)]);
      const float support =
          .5f * (fabsf(radial_cross(dx, edge)) + fabsf(radial_cross(dy, edge)));
      const float guard =
          32 * FLT_EPSILON * (1 + fabsf(predicate) + support + constant);
      return predicate + support < -guard;
    };
    if (!inside(center, value))
      return false;
    const int low_sector =
        sector_index(query.theta - query.angular_extent - c.phase);
    const int high_sector =
        std::min(sector_index(query.theta + query.angular_extent - c.phase),
                 low_sector + num_edges - 1);
    for (int i = low_sector; i <= high_sector; ++i)
      if (i != center &&
          !inside(i, radial_cross(
                         query.point,
                         radial_multiply(
                             c_value, edges[radial_wrap_index(i, num_edges)])) -
                         constant))
        return false;
    return true;
  }

  /** @brief Finds the first visible boundary sheet in a bounded radial scan.
   * @tparam COUNT Enables optional work counters.
   * @tparam SKIP Enables conservative block rejection.
   * @tparam DERIVATIVES Computes normal, derivative, and gradient fields.
   * @tparam Observer Optional interval observer with prepared footprint query bounds.
   * @param p Chart point.
   * @param theta Argument of p, radians.
   * @param first_interval Conservative first possible interval.
   * @param stats Optional counter destination when COUNT is true.
   * @param observer Receives interval sectors and endpoint predicates.
   * @return Earliest hit, or a result with interval -1.
   * @details A point already inside the first supplied star returns its start.
   * The caller must not skip an interval that could contain an earlier hit.
   * For a sweep growing from zero scale, an unordered root fraction selects
   * the interval endpoint. */
  template <bool COUNT = false, bool SKIP = true, bool DERIVATIVES = true,
            class Observer = NoObserver>
  HS_HOT_FLASH_MEMBER Hit sample(Point p, float theta, int first_interval = 0,
                                 Stats *stats = nullptr,
                                 Observer *observer = nullptr) const {
    if (first_interval < 0)
      first_interval = 0;
    if (first_interval >= num_intervals)
      return {};
    if (p.x == 0.0f && p.y == 0.0f)
      return {first_interval, 0.0f};
    QueryBounds query;
    if constexpr (std::is_same_v<Observer, NoObserver>)
      query = query_bounds(p, theta, 0.0f);
    else
      query = observer->query;

    int start_sector = sector_index(theta - knots[first_interval].phase);
    int start_edge = radial_wrap_index(start_sector, num_edges);
    const Point start_edge_value =
        multiply(knots[first_interval], edges[start_edge]);
    float start_value =
        radial_cross(p, start_edge_value) - determinant * norms[first_interval];
    if (norms[first_interval] > 0.0f && start_value <= 0.0f)
      return finish<DERIVATIVES>(p, first_interval, 0.0f, start_edge);

    int checked_block = -1;
    for (int i = first_interval; i < num_intervals; ++i) {
      if constexpr (SKIP) {
        if (i / BLOCK_INTERVALS != checked_block) {
          checked_block = i / BLOCK_INTERVALS;
          count<COUNT>(stats, &Stats::coarse_checks);
          const int next = skip_outside(query, i);
          if (next != i) {
            if constexpr (COUNT)
              if (stats)
                stats->skipped += uint32_t(next - i);
            i = next;
            if (i >= num_intervals)
              break;
            start_sector = sector_index(theta - knots[i].phase);
            start_edge = radial_wrap_index(start_sector, num_edges);
            start_value =
                radial_cross(p, multiply(knots[i], edges[start_edge])) -
                determinant * norms[i];
            --i;
            continue;
          }
        }
      }
      count<COUNT>(stats, &Stats::intervals);
      const int end_sector = sector_index(theta - knots[i + 1].phase);
      const int unwrapped_end_edge = start_edge + end_sector - start_sector;
      const int end_edge =
          unwrapped_end_edge < 0            ? unwrapped_end_edge + num_edges
          : unwrapped_end_edge >= num_edges ? unwrapped_end_edge - num_edges
                                            : unwrapped_end_edge;
      const Point end_vector = multiply(knots[i + 1], edges[end_edge]);
      const float end_value =
          radial_cross(p, end_vector) - determinant * norms[i + 1];
      if constexpr (!std::is_same_v<Observer, NoObserver>)
        observer->interval(*this, i, start_sector, start_value, end_sector,
                           end_value);
      if (norms[i] == 0.0f) {
        if (end_value <= 0.0f) {
          count<COUNT>(stats, &Stats::roots);
          const float t =
              radial_cross(p, end_vector) / (determinant * norms[i + 1]);
          return finish<DERIVATIVES>(p, i, hs::clamp(t, 0.0f, 1.0f), end_edge);
        }
      } else if (start_sector == end_sector) {
        if (end_value <= 0.0f) {
          count<COUNT>(stats, &Stats::roots);
          return solve<DERIVATIVES>(p, i, start_edge, 0.0f, 1.0f, start_value);
        }
      } else {
        const int corner = end_sector > start_sector ? end_edge : start_edge;
        if ((corner & 1) == 0 || end_value <= 0.0f) {
          count<COUNT>(stats, &Stats::events);
          const Point corner_start = multiply(knots[i], vertices[corner]);
          const Point corner_delta =
              radial_multiply(pieces[i].delta, vertices[corner]);
          const float denominator = radial_cross(p, corner_delta);
          const float event = hs::clamp(
              -radial_cross(p, corner_start) / denominator, 0.0f, 1.0f);
          const float event_value = predicate(p, i, start_edge, event);
          if (event_value <= 0.0f) {
            count<COUNT>(stats, &Stats::roots);
            return solve<DERIVATIVES>(p, i, start_edge, 0.0f, event,
                                      start_value);
          }
          if (end_value <= 0.0f) {
            count<COUNT>(stats, &Stats::roots);
            const float end_event_value = predicate(p, i, end_edge, event);
            return solve<DERIVATIVES>(p, i, end_edge, event, 1.0f,
                                      end_event_value);
          }
        }
      }
      start_sector = end_sector;
      start_edge = end_edge;
      start_value = end_value;
    }
    return {};
  }

private:
  /** @brief Conservative displacement from a block's endpoint sweep. */
  struct Block {
    float low_phase;      ///< Minimum fine-knot phase.
    float high_phase;     ///< Maximum fine-knot phase.
    float radius;         ///< Maximum fine multiplier norm.
    float point_error;    ///< Edge length times multiplier deviation.
    float constant_error; ///< Radial predicate displacement bound.
    float reach;          ///< Maximum coarse edge norm.
  };
  /** @brief Coefficients shared by every edge of a radial interval. */
  struct Piece {
    Point delta;         ///< Similarity multiplier difference.
    float quadratic;     ///< Positive magnitude of the quadratic coefficient.
    float radial_linear; ///< Linear coefficient of determinant times radius².
  };

  /** @param a Similarity multiplier.
   * @param b Chart vector.
   * @return Transformed vector. */
  static HS_HOT_INLINE Point multiply(const Knot &a, Point b) {
    return radial_multiply(Point{a.real, a.imag}, b);
  }
  /** @param angle Relative point argument.
   * @return Unwrapped edge sector. */
  HS_HOT_INLINE int sector_index(float angle) const {
    return radial_floor_int(angle * inverse_sector);
  }
  /** @param p Chart point.
   * @param i Interval.
   * @param edge_index Unit-star edge.
   * @param t Interval fraction.
   * @return Positive outside, negative inside the edge half-plane. */
  HS_HOT_INLINE float predicate(Point p, int i, int edge_index, float t) const {
    const Point c{knots[i].real + t * pieces[i].delta.x,
                  knots[i].imag + t * pieces[i].delta.y};
    return radial_cross(p, radial_multiply(c, edges[edge_index])) -
           determinant * (c.x * c.x + c.y * c.y);
  }
  /** @param p Chart point.
   * @param i Interval.
   * @param edge_index Unit-star edge.
   * @param low,high Sector's interval range.
   * @param value Predicate at low, nonnegative.
   * @return Closed-form descending root and sheet derivatives. */
  template <bool DERIVATIVES>
  HS_HOT_INLINE Hit solve(Point p, int i, int edge_index, float low, float high,
                          float value) const {
    if (value <= 0.0f)
      return finish<DERIVATIVES>(p, i, low, edge_index);
    const Piece &piece = pieces[i];
    const float linear =
        radial_cross(p, radial_multiply(piece.delta, edges[edge_index])) -
        piece.radial_linear - 2.0f * piece.quadratic * low;
    const float root = sqrtf(linear * linear + 4.0f * piece.quadratic * value);
    const float t =
        low + (linear > 0.0f ? (linear + root) / (2.0f * piece.quadratic)
                             : 2.0f * value / (root - linear));
    return finish<DERIVATIVES>(p, i, hs::clamp(t, low, high), edge_index);
  }
  /** @param p Chart point.
   * @param i Interval.
   * @param t Hit fraction.
   * @param edge_index Unit-star edge.
   * @return Hit with the implicit sheet's chart derivatives. */
  template <bool DERIVATIVES>
  HS_HOT_INLINE Hit finish(Point p, int i, float t, int edge_index) const {
    if constexpr (!DERIVATIVES)
      return {i, t, edge_index};
    const Piece &piece = pieces[i];
    const Point c{knots[i].real + t * piece.delta.x,
                  knots[i].imag + t * piece.delta.y};
    const Point edge = radial_multiply(c, edges[edge_index]);
    const float derivative =
        radial_cross(p, radial_multiply(piece.delta, edges[edge_index])) -
        piece.radial_linear - 2.0f * piece.quadratic * t;
    const float inv = derivative != 0.0f ? -1.0f / derivative : 0.0f;
    return {i,          t,
            edge_index, {edge.y, -edge.x},
            derivative, {edge.y * inv, -edge.x * inv}};
  }
  /** @tparam COUNT Whether counting is compiled in.
   * @param stats Counter destination.
   * @param member Counter to increment. */
  template <bool COUNT>
  static HS_HOT_INLINE void count(Stats *stats, uint32_t Stats::*member) {
    if constexpr (COUNT)
      if (stats)
        ++(stats->*member);
  }

  int num_edges = 0;           ///< Prepared edge count.
  int num_intervals = 0;       ///< Prepared radial interval count.
  float inner = 0.0f;          ///< Inner-to-outer vertex radius ratio.
  float sector = 0.0f;         ///< One edge's angular extent.
  float inverse_sector = 0.0f; ///< Reciprocal edge sector width.
  float determinant = 0.0f;    ///< Unit edge's signed vertex determinant.
  float edge_length = 0.0f;    ///< Maximum prepared unit edge length.
  Point *vertices = nullptr;   ///< Unit-star vertices.
  Point *edges = nullptr;      ///< Unit-star edge vectors.
  Knot *knots = nullptr;       ///< Prepared similarity multipliers.
  float *norms = nullptr;      ///< Squared multiplier norms.
  Piece *pieces = nullptr;     ///< Shared interval coefficients.
  Block *blocks = nullptr;     ///< Coarse sweep rejection certificates.
  EdgeMetrics first_edge_metrics{}, second_edge_metrics{};
  bool increasing_radius = false;
};

} // namespace Scan
