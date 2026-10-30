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
      const float RADIUS = (i & 1) ? inner : 1.0f;
      vertices[i] = {RADIUS * cosf(i * sector), RADIUS * sinf(i * sector)};
    }
    for (int i = 0; i < num_edges; ++i) {
      const Point NEXT = vertices[(i + 1) % num_edges];
      edges[i] = {NEXT.x - vertices[i].x, NEXT.y - vertices[i].y};
      edge_length = fmaxf(edge_length, sqrtf(edges[i].x * edges[i].x +
                                             edges[i].y * edges[i].y));
    }
    auto metrics = [&](int i) {
      const Point A = vertices[i], E = edges[i];
      return EdgeMetrics{A.x * E.x + A.y * E.y, E.x * E.x + E.y * E.y,
                         A.x * E.y - A.y * E.x};
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
        const float DX = knots[i].real - knots[i - 1].real;
        const float DY = knots[i].imag - knots[i - 1].imag;
        pieces[i - 1] = {{DX, DY},
                         determinant * (DX * DX + DY * DY),
                         2.0f * determinant *
                             (knots[i - 1].real * DX + knots[i - 1].imag * DY)};
        increasing_radius &=
            pieces[i - 1].radial_linear >= 0 && pieces[i - 1].quadratic > 0;
      }
    }
    for (int first = 0; first < num_intervals; first += BLOCK_INTERVALS) {
      const int END = std::min(first + BLOCK_INTERVALS, num_intervals);
      const float INV = 1.0f / (END - first);
      const Knot &START = knots[first], &FINISH = knots[END];
      float low = START.phase, high = START.phase, deviation2 = 0;
      float radius2 = norms[first];
      for (int i = first + 1; i <= END; ++i) {
        low = fminf(low, knots[i].phase);
        high = fmaxf(high, knots[i].phase);
        radius2 = fmaxf(radius2, norms[i]);
        const float T = (i - first) * INV;
        const float DX =
            knots[i].real - (START.real + T * (FINISH.real - START.real));
        const float DY =
            knots[i].imag - (START.imag + T * (FINISH.imag - START.imag));
        deviation2 = fmaxf(deviation2, DX * DX + DY * DY);
      }
      const float RADIUS = sqrtf(radius2);
      const float COARSE_RADIUS = sqrtf(fmaxf(norms[first], norms[END]));
      const float DEVIATION = sqrtf(deviation2) + 16.0f * FLT_EPSILON * RADIUS;
      blocks[first / BLOCK_INTERVALS] = {low,
                                         high,
                                         RADIUS * (1.0f + 16.0f * FLT_EPSILON),
                                         edge_length * DEVIATION,
                                         determinant * DEVIATION *
                                             (2.0f * COARSE_RADIUS + DEVIATION),
                                         edge_length * COARSE_RADIUS};
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
    const float RADIUS = sqrtf(p.x * p.x + p.y * p.y);
    const float ANGLE =
        extent == 0.0f
            ? 0.0f
            : (RADIUS > extent ? fminf(math::PI_F, extent / (RADIUS - extent))
                               : math::PI_F);
    return {p, theta, RADIUS, extent, ANGLE};
  }

  /** @brief Skips the remainder of a block only when its whole sweep misses.
   * @param query Chart disk with conservative radius and angular bounds.
   * @param interval Current fine interval.
   * @return Block end when certified outside; otherwise interval. */
  HS_HOT_INLINE int skip_outside(const QueryBounds &query, int interval) const {
    if (interval < 0 || interval >= num_intervals)
      return interval;
    const int FIRST = interval / BLOCK_INTERVALS * BLOCK_INTERVALS;
    const int END = std::min(FIRST + BLOCK_INTERVALS, num_intervals);
    const Block &BLOCK = blocks[FIRST / BLOCK_INTERVALS];
    if (query.radius - query.extent > BLOCK.radius)
      return END;
    const int FIRST_SECTOR =
        sector_index(query.theta - query.angular_extent - BLOCK.high_phase);
    const int LAST_SECTOR =
        sector_index(query.theta + query.angular_extent - BLOCK.low_phase);
    if (FIRST_SECTOR != LAST_SECTOR)
      return interval;
    const Point EDGE = edges[radial_wrap_index(FIRST_SECTOR, num_edges)];
    const float A = radial_cross(query.point, multiply(knots[FIRST], EDGE)) -
                    determinant * norms[FIRST];
    const float B = radial_cross(query.point, multiply(knots[END], EDGE)) -
                    determinant * norms[END];
    const float ROUNDING = 32.0f * FLT_EPSILON *
                           (1.0f + query.radius * BLOCK.reach +
                            determinant * BLOCK.radius * BLOCK.radius);
    const float MARGIN = (query.radius + query.extent) * BLOCK.point_error +
                         BLOCK.constant_error + query.extent * BLOCK.reach +
                         ROUNDING;
    return fminf(A, B) > MARGIN ? END : interval;
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
    const float EXTENT = .5f * (sqrtf(dx.x * dx.x + dx.y * dx.y) +
                                sqrtf(dy.x * dy.x + dy.y * dy.y)) +
                         1e-6f;
    const auto QUERY = query_bounds(p, theta, EXTENT);
    return contains_footprint(QUERY, dx, dy, knots[knot_index]);
  }

  /** @brief Certifies an affine footprint inside a similarity star.
   * @param query Conservative spatial and angular footprint bounds.
   * @param dx,dy Full pixel-coordinate derivatives.
   * @param c Similarity multiplier and its unwrapped argument.
   * @return True only if every touched edge half-plane contains the footprint. */
  HS_HOT_INLINE bool contains_footprint(const QueryBounds &query, Point dx,
                                        Point dy, const Knot &c) const {
    const int CENTER = sector_index(query.theta - c.phase);
    const Point EDGE = multiply(c, edges[radial_wrap_index(CENTER, num_edges)]);
    const float VALUE = radial_cross(query.point, EDGE) -
                        determinant * (c.real * c.real + c.imag * c.imag);
    return contains_footprint(query, dx, dy, c, CENTER, VALUE);
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
    const Point C{c.real, c.imag};
    const float CONSTANT = determinant * (C.x * C.x + C.y * C.y);
    auto inside = [&](int index, float predicate) HS_HOT_INLINE {
      const Point EDGE =
          radial_multiply(C, edges[radial_wrap_index(index, num_edges)]);
      const float SUPPORT =
          .5f * (fabsf(radial_cross(dx, EDGE)) + fabsf(radial_cross(dy, EDGE)));
      const float GUARD =
          32 * FLT_EPSILON * (1 + fabsf(predicate) + SUPPORT + CONSTANT);
      return predicate + SUPPORT < -GUARD;
    };
    if (!inside(center, value))
      return false;
    const int LOW_SECTOR =
        sector_index(query.theta - query.angular_extent - c.phase);
    const int HIGH_SECTOR =
        std::min(sector_index(query.theta + query.angular_extent - c.phase),
                 LOW_SECTOR + num_edges - 1);
    for (int i = LOW_SECTOR; i <= HIGH_SECTOR; ++i)
      if (i != center &&
          !inside(i,
                  radial_cross(query.point,
                               radial_multiply(
                                   C, edges[radial_wrap_index(i, num_edges)])) -
                      CONSTANT))
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
   * The caller must not skip an interval that could contain an earlier hit. */
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
    QueryBounds QUERY;
    if constexpr (std::is_same_v<Observer, NoObserver>)
      QUERY = query_bounds(p, theta, 0.0f);
    else
      QUERY = observer->query;

    int start_sector = sector_index(theta - knots[first_interval].phase);
    int start_edge = radial_wrap_index(start_sector, num_edges);
    const Point START_EDGE = multiply(knots[first_interval], edges[start_edge]);
    float start_value =
        radial_cross(p, START_EDGE) - determinant * norms[first_interval];
    if (norms[first_interval] > 0.0f && start_value <= 0.0f)
      return finish<DERIVATIVES>(p, first_interval, 0.0f, start_edge);

    int checked_block = -1;
    for (int i = first_interval; i < num_intervals; ++i) {
      if constexpr (SKIP) {
        if (i / BLOCK_INTERVALS != checked_block) {
          checked_block = i / BLOCK_INTERVALS;
          count<COUNT>(stats, &Stats::coarse_checks);
          const int NEXT = skip_outside(QUERY, i);
          if (NEXT != i) {
            if constexpr (COUNT)
              if (stats)
                stats->skipped += uint32_t(NEXT - i);
            i = NEXT;
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
      const int END_SECTOR = sector_index(theta - knots[i + 1].phase);
      int end_edge = start_edge + END_SECTOR - start_sector;
      if (end_edge < 0)
        end_edge += num_edges;
      else if (end_edge >= num_edges)
        end_edge -= num_edges;
      const int END_EDGE = end_edge;
      const Point END_VECTOR = multiply(knots[i + 1], edges[END_EDGE]);
      const float END_VALUE =
          radial_cross(p, END_VECTOR) - determinant * norms[i + 1];
      if constexpr (!std::is_same_v<Observer, NoObserver>)
        observer->interval(*this, i, start_sector, start_value, END_SECTOR,
                           END_VALUE);
      if (norms[i] == 0.0f) {
        if (END_VALUE <= 0.0f) {
          count<COUNT>(stats, &Stats::roots);
          const float T =
              radial_cross(p, END_VECTOR) / (determinant * norms[i + 1]);
          return finish<DERIVATIVES>(p, i, std::clamp(T, 0.0f, 1.0f), END_EDGE);
        }
      } else if (start_sector == END_SECTOR) {
        if (END_VALUE <= 0.0f) {
          count<COUNT>(stats, &Stats::roots);
          return solve<DERIVATIVES>(p, i, start_edge, 0.0f, 1.0f, start_value);
        }
      } else {
        const int CORNER = END_SECTOR > start_sector ? END_EDGE : start_edge;
        if ((CORNER & 1) == 0 || END_VALUE <= 0.0f) {
          count<COUNT>(stats, &Stats::events);
          const Point CORNER_START = multiply(knots[i], vertices[CORNER]);
          const Point CORNER_DELTA =
              radial_multiply(pieces[i].delta, vertices[CORNER]);
          const float DENOMINATOR = radial_cross(p, CORNER_DELTA);
          const float EVENT = std::clamp(
              -radial_cross(p, CORNER_START) / DENOMINATOR, 0.0f, 1.0f);
          const int START_EDGE = start_edge;
          const float EVENT_VALUE = predicate(p, i, START_EDGE, EVENT);
          if (EVENT_VALUE <= 0.0f) {
            count<COUNT>(stats, &Stats::roots);
            return solve<DERIVATIVES>(p, i, START_EDGE, 0.0f, EVENT,
                                      start_value);
          }
          if (END_VALUE <= 0.0f) {
            count<COUNT>(stats, &Stats::roots);
            const float START_VALUE = predicate(p, i, END_EDGE, EVENT);
            return solve<DERIVATIVES>(p, i, END_EDGE, EVENT, 1.0f, START_VALUE);
          }
        }
      }
      start_sector = END_SECTOR;
      start_edge = END_EDGE;
      start_value = END_VALUE;
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
    const Point C{knots[i].real + t * pieces[i].delta.x,
                  knots[i].imag + t * pieces[i].delta.y};
    return radial_cross(p, radial_multiply(C, edges[edge_index])) -
           determinant * (C.x * C.x + C.y * C.y);
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
    const Piece &PIECE = pieces[i];
    const float LINEAR =
        radial_cross(p, radial_multiply(PIECE.delta, edges[edge_index])) -
        PIECE.radial_linear - 2.0f * PIECE.quadratic * low;
    const float ROOT = sqrtf(LINEAR * LINEAR + 4.0f * PIECE.quadratic * value);
    const float DENOMINATOR = ROOT - LINEAR;
    HS_CHECK(DENOMINATOR > 0.0f, "RadialSurface: descending root denominator");
    const float T = low + 2.0f * value / DENOMINATOR;
    return finish<DERIVATIVES>(p, i, std::clamp(T, low, high), edge_index);
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
    const Piece &PIECE = pieces[i];
    const Point C{knots[i].real + t * PIECE.delta.x,
                  knots[i].imag + t * PIECE.delta.y};
    const Point EDGE = radial_multiply(C, edges[edge_index]);
    const float DERIVATIVE =
        radial_cross(p, radial_multiply(PIECE.delta, edges[edge_index])) -
        PIECE.radial_linear - 2.0f * PIECE.quadratic * t;
    const float INV = DERIVATIVE != 0.0f ? -1.0f / DERIVATIVE : 0.0f;
    return {i,          t,
            edge_index, {EDGE.y, -EDGE.x},
            DERIVATIVE, {EDGE.y * INV, -EDGE.x * INV}};
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
