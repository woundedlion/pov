/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file swept_face.h
 * @brief Column-time polygon sweep state and meridian crossings. */

#include "core/render/scan/raster.h"

namespace Scan {

/// Angle past a swept face's edges within which its column pairs shade: a
/// column of antialiasing with slack.
template <int W>
inline constexpr float SWEEP_ROW_REACH = 1.25f * math::TWO_PI_F / W;

/**
 * @brief A face a deformation scan sweeps, standing at one column pair's
 *        sweep time at a time.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Per frame it projects the face's vertices, as SDF::PreparedFace projects
 * them, about the vertex centroid at the frame's middle, at the frame's
 * start, middle and end, and fits each projected vertex, a convex face's
 * half-planes and the face's size as quadratics in time through those
 * three. Per column pair it evaluates the fits by Horner's rule and finds
 * where the pair's meridian, a line in the projection, crosses the edges. A pixel is inside
 * where the meridian, walked from the pole away from the face, has crossed
 * an odd number of edges on the way to it; its distance is SDF::PreparedFace's:
 * signed nearest-edge distance, with convex interior half-plane distances, in
 * radians for a face whose size is. Its rows are those within
 * SWEEP_ROW_REACH of the edges clipped to the meridian's antialiasing strip,
 * on either arm and past either pole, and from the pole on when the face
 * holds it.
 */
template <int W, int H> class SweptFace {
public:
  /**
   * @brief Fits a sweep's face over the frame.
   * @param sweep The face and its vertex paths.
   * @param arena Arena holding the face's fits and per-column state for its
   *        scope.
   */
  HS_HOT_FLASH_MEMBER SweptFace(const Scan::PoleFaceSweep &sweep, Arena &arena)
      : count(static_cast<int>(sweep.count)), north(sweep.north) {
    edges = arena.allocate_n<Edge>(sweep.count);
    planes = arena.allocate_n<Plane>(sweep.count);
    motion = arena.allocate_n<Motion>(sweep.count);
    line2 = arena.allocate_n<float>(sweep.count);
    span = arena.allocate_n<Span>(sweep.count);
    order = arena.allocate_n<uint8_t>(sweep.count);
    for (int k = 0; k < count; ++k)
      order[k] = static_cast<uint8_t>(k);
    crossings[0] = arena.allocate_n<float>(sweep.count);
    crossings[1] = arena.allocate_n<float>(sweep.count);
    math::Vector middle(0.0f, 0.0f, 0.0f);
    for (int k = 0; k < count; ++k)
      middle = middle + sweep.paths[k].middle;
    axis = unit(middle);
    basis_u = math::perpendicular_axis(axis);
    basis_w = unit(math::cross(axis, basis_u));
    // The fitted quantities at the frame's start, middle and end.
    float sizes[3];
    bool convex_at[3];
    for (int j = 0; j < 3; ++j) {
      const float t = 0.5f * static_cast<float>(j);
      math::Vector sum(0.0f, 0.0f, 0.0f);
      for (int k = 0; k < count; ++k) {
        const math::Vector v = sweep.paths[k].at(t);
        sum = sum + v;
        float px, py;
        project(v, px, py);
        edges[k].x = px;
        edges[k].y = py;
        motion[k].x.sample(j, px);
        motion[k].y.sample(j, py);
      }
      float area2 = 0.0f, radius2 = 0.0f;
      stand_edges(area2, radius2);
      convex_at[j] = half_planes(area2);
      for (int k = 0; k < count; ++k) {
        motion[k].nx.sample(j, planes[k].nx);
        motion[k].ny.sample(j, planes[k].ny);
        motion[k].off.sample(j, planes[k].off);
      }
      float cx, cy;
      project(sum, cx, cy);
      sizes[j] = fmaxf(sqrtf(nearest_squared(cx, cy)),
                       sqrtf(radius2) * SDF::MIN_SIZE_RADIUS_RATIO);
    }
    convex = convex_at[0] && convex_at[1] && convex_at[2];
    linear = sizes[1] < LINEAR_SIZE;
    for (int j = 0; j < 3; ++j)
      if (!linear)
        sizes[j] = math::fast_atan2(sizes[j], 1.0f);
    for (int k = 0; k < count; ++k) {
      Motion &m = motion[k];
      m.x.fit();
      m.y.fit();
      m.nx.fit();
      m.ny.fit();
      m.off.fit();
    }
    for (int j = 0; j < 3; ++j)
      size_fit.sample(j, sizes[j]);
    size_fit.fit();
  }

  /** @brief The size shading normalizes by at the column pair column()
   *  stood the face at, as SDF::PreparedFace::size.
   *  @return Plane units or, for a large face, radians. */
  float size() const { return face_size; }

  /**
   * @brief Stands the face at a column pair's sweep time.
   * @param time Frame's column basis.
   * @param x The pair's first column, in [0, W / 2).
   * @param arms Arms to bound: bit 0 the first, bit 1 the second.
   * @param rows Receives each bounded arm's rows, culled where it misses the
   *        face.
   * @return Whether a bounded arm has rows; false where the face collapses,
   *         as a build culls.
   */
  HS_HOT_FLASH_MEMBER bool column(const Scan::ColumnBasis &time, int x,
                                  int arms, SDF::Bounds rows[2]) {
    rows[0] = rows[1] = SDF::BOUNDS_CULLED;
    const float t = Spacetime::SamplingSurface<W>::time(x);
    for (int k = 0; k < count; ++k) {
      edges[k].x = motion[k].x.at(t);
      edges[k].y = motion[k].y.at(t);
    }
    float area2 = 0.0f, radius2 = 0.0f;
    for (int k = 0; k < count; ++k) {
      Edge &e = edges[k];
      const Edge &next = edges[k + 1 == count ? 0 : k + 1];
      e.ex = next.x - e.x;
      e.ey = next.y - e.y;
      const float length2 = e.ex * e.ex + e.ey * e.ey;
      e.inv_length2 = length2 > 0.0f ? 1.0f / length2 : 0.0f;
      area2 += e.x * next.y - next.x * e.y;
      radius2 = fmaxf(radius2, e.x * e.x + e.y * e.y);
      if (convex) {
        planes[k].nx = motion[k].nx.at(t);
        planes[k].ny = motion[k].ny.at(t);
        planes[k].off = motion[k].off.at(t);
      }
    }
    if (fabsf(area2) <= SDF::COLLAPSED_AREA_RATIO * radius2)
      return false;
    face_size = size_fit.at(t);
    azimuth = time.azimuth[x];
    up = time.up[x];
    // The meridian's line in the projection, by its unit normal and offset.
    const math::Vector meridian = math::cross(azimuth, up);
    float la = math::dot(basis_u, meridian), lb = math::dot(basis_w, meridian);
    const float inv_norm = 1.0f / sqrtf(la * la + lb * lb);
    la *= inv_norm;
    lb *= inv_norm;
    const float le = math::dot(axis, meridian) * inv_norm;
    constexpr float REACH = SWEEP_ROW_REACH<W>;
    float lo[2] = {2.0f, 2.0f}, hi[2] = {-2.0f, -2.0f};
    // A projected point's arm and polar cosine along the meridian.
    const auto include = [&](float px, float py) HS_HOT_INLINE {
      const math::Vector q = axis + basis_u * px + basis_w * py;
      const int arm = math::dot(q, azimuth) >= 0.0f ? 0 : 1;
      const float c = math::dot(q, up) / sqrtf(math::dot(q, q));
      lo[arm] = fminf(lo[arm], c);
      hi[arm] = fmaxf(hi[arm], c);
      return std::pair<int, float>{arm, c};
    };
    crossed[0] = crossed[1] = 0;
    // Each edge's extent along the line, by the projection on its direction.
    along_x = -lb;
    along_y = la;
    const float pole_depth = math::dot(axis, up);
    const float pole_along = pole_depth != 0.0f
                                 ? (along_x * math::dot(basis_u, up) +
                                    along_y * math::dot(basis_w, up)) /
                                       pole_depth
                                 : FLT_MAX;
    bool touches_pole = false;
    float from_line = la * edges[0].x + lb * edges[0].y + le;
    float from_along = along_x * edges[0].x + along_y * edges[0].y;
    for (int k = 0; k < count; ++k) {
      Edge &e = edges[k];
      const Edge &next = edges[k + 1 == count ? 0 : k + 1];
      const float to_line = la * next.x + lb * next.y + le;
      const float to_along = along_x * next.x + along_y * next.y;
      span[k] = {fminf(from_along, to_along), fmaxf(from_along, to_along)};
      if (fabsf(from_line) <= REACH)
        include(e.x, e.y);
      float inv_delta = 0.0f;
      // No pixel of the pair, on the meridian's line, comes nearer the edge
      // than its nearer end, unless the edge crosses the line.
      if ((from_line > 0.0f) == (to_line > 0.0f)) {
        const float apart = fminf(fabsf(from_line), fabsf(to_line));
        line2[k] = apart * apart;
      } else {
        line2[k] = 0.0f;
        inv_delta = 1.0f / (from_line - to_line);
        const float s = from_line * inv_delta;
        const auto [arm, c] = include(e.x + e.ex * s, e.y + e.ey * s);
        crossings[arm][crossed[arm]++] = c;
      }
      if (line2[k] <= REACH * REACH) {
        float begin = 0.0f, end = 1.0f;
        const bool from_outside = fabsf(from_line) > REACH;
        const bool to_outside = fabsf(to_line) > REACH;
        if (from_outside || to_outside) {
          if (inv_delta == 0.0f)
            inv_delta = 1.0f / (from_line - to_line);
          if (from_outside) {
            begin = (from_line - copysignf(REACH, from_line)) * inv_delta;
            include(e.x + e.ex * begin, e.y + e.ey * begin);
          }
          if (to_outside) {
            end = (from_line - copysignf(REACH, to_line)) * inv_delta;
            include(e.x + e.ex * end, e.y + e.ey * end);
          }
        }
        const float first = from_along + (to_along - from_along) * begin;
        const float last = from_along + (to_along - from_along) * end;
        touches_pole |= pole_along >= fminf(first, last) &&
                        pole_along <= fmaxf(first, last);
      }
      from_line = to_line;
      from_along = to_along;
    }
    // The edges nearest the line first; the order holds from the last pair,
    // so the insertion sort runs over nearly sorted keys.
    for (int i = 1; i < count; ++i) {
      const uint8_t k = order[i];
      const float key = line2[k];
      int j = i;
      for (; j > 0 && line2[order[j - 1]] > key; --j)
        order[j] = order[j - 1];
      order[j] = k;
    }
    // The pole is covered by the interior or by an edge's antialiasing strip.
    if ((crossed[0] & 1) || touches_pole) {
      for (int arm = 0; arm < 2; ++arm) {
        if (north)
          hi[arm] = 1.0f;
        else
          lo[arm] = -1.0f;
      }
    }
    // Each arm's polar angles, widened by the reach and math::fast_acos's
    // error, and past either pole onto the other arm.
    const float pad = REACH + static_cast<float>(math::FAST_ACOS_MAX_ERROR);
    float first[2], last[2];
    for (int arm = 0; arm < 2; ++arm) {
      const bool empty = lo[arm] > hi[arm];
      first[arm] = empty ? math::PI_F + pad : math::fast_acos(hi[arm]) - pad;
      last[arm] = empty ? -pad : math::fast_acos(lo[arm]) + pad;
    }
    bool any = false;
    for (int arm = 0; arm < 2; ++arm) {
      const int other = 1 - arm;
      float from = first[arm], to = last[arm];
      if (first[other] < 0.0f) {
        from = fminf(from, 0.0f);
        to = fmaxf(to, -first[other]);
      }
      if (last[other] > math::PI_F) {
        from = fminf(from, 2.0f * math::PI_F - last[other]);
        to = fmaxf(to, math::PI_F);
      }
      if (!(arms & (1 << arm)) || from > to)
        continue;
      rows[arm] =
          SDF::phi_bounds_to_rows<H>(fmaxf(from, 0.0f), fminf(to, math::PI_F));
      any = true;
    }
    return any;
  }

  /**
   * @brief The face's signed distance at a pixel of the column pair column()
   *        stood it at.
   * @param p Unit pixel direction.
   * @param res Receives the distance, negative inside, in the distance and
   *        raw-distance registers, and size().
   */
  HS_HOT_INLINE void probe(const math::Vector &p,
                           SDF::DistanceResult &res) const {
    const float depth = math::dot(p, axis);
    if (depth < math::TOLERANCE) {
      res = SDF::DistanceResult(SDF::FAR_SENTINEL, 0.0f, SDF::FAR_SENTINEL,
                                0.0f, face_size);
      return;
    }
    const float inv = 1.0f / depth;
    const float px = math::dot(p, basis_u) * inv,
                py = math::dot(p, basis_w) * inv;
    float signed_d;
    if (convex) {
      signed_d = -FLT_MAX;
      int winning = 0;
      for (int k = 0; k < count; ++k) {
        const Plane &h = planes[k];
        const float d = h.nx * px + h.ny * py + h.off;
        if (d > signed_d)
          winning = k;
        signed_d = fmaxf(signed_d, d);
      }
      if (signed_d > 0.0f) {
        const Edge &edge = edges[winning];
        const float foot = ((px - edge.x) * edge.ex + (py - edge.y) * edge.ey) *
                           edge.inv_length2;
        if (foot < 0.0f || foot > 1.0f)
          signed_d = sqrtf(nearest_squared(px, py));
      }
    } else {
      const float c = math::dot(p, up);
      int passed;
      if (math::dot(p, azimuth) >= 0.0f) {
        passed = 0;
        for (int i = 0; i < crossed[0]; ++i)
          passed += north ? crossings[0][i] < c : crossings[0][i] > c;
      } else {
        passed = crossed[0];
        for (int i = 0; i < crossed[1]; ++i)
          passed += north ? crossings[1][i] > c : crossings[1][i] < c;
      }
      // From the edge nearest the last pixel, then the edges by their
      // distance from the meridian's line until that distance passes the
      // nearest found; an edge whose extent along the line keeps it farther
      // is passed over.
      float best = edge_squared(edges[nearest], px, py);
      const float along = along_x * px + along_y * py;
      for (int i = 0; i < count; ++i) {
        const int k = order[i];
        const float across = line2[k];
        if (across >= best)
          break;
        const float gap =
            fmaxf(fmaxf(span[k].lo - along, along - span[k].hi), 0.0f);
        if (across + gap * gap >= best)
          continue;
        const float d = edge_squared(edges[k], px, py);
        if (d < best) {
          best = d;
          nearest = k;
        }
      }
      const float d = sqrtf(best);
      signed_d = passed & 1 ? -d : d;
    }
    if (!linear)
      signed_d = math::fast_atan2(signed_d, 1.0f);
    res = SDF::DistanceResult(signed_d, 0.0f, signed_d, 0.0f, face_size);
  }

private:
  /// Size, in plane units, below which a face reports plane distances, as
  /// SDF::PreparedFace does.
  static constexpr float LINEAR_SIZE = 0.2f;

  /** @brief A quadratic in time through values at t = 0, 1/2 and 1. */
  struct Quadratic {
    float c0 = 0.0f, c1 = 0.0f, c2 = 0.0f; ///< Coefficients, or samples.

    /// Records the sample at t = j / 2.
    HS_FLASH_INLINE void sample(int j, float value) {
      (j == 0 ? c0 : j == 1 ? c1 : c2) = value;
    }

    /// Turns the three samples into the quadratic through them.
    HS_FLASH_INLINE void fit() {
      const float s0 = c0, s1 = c1, s2 = c2;
      c1 = -3.0f * s0 + 4.0f * s1 - s2;
      c2 = 2.0f * s0 - 4.0f * s1 + 2.0f * s2;
    }

    /// The value at a frame fraction.
    HS_HOT_INLINE float at(float t) const { return c0 + t * (c1 + t * c2); }
  };

  /** @brief One projected vertex and the edge from it to the next. */
  struct Edge {
    float x, y;        ///< The vertex.
    float ex, ey;      ///< To the next.
    float inv_length2; ///< The edge's inverse squared length; 0 collapsed.
  };

  /** @brief An edge's extent along the meridian's line. */
  struct Span {
    float lo, hi; ///< Least and greatest projection on the line's direction.
  };

  /** @brief A convex face's outward half-plane along an edge. */
  struct Plane {
    float nx, ny, off; ///< Unit normal and offset.
  };

  /** @brief A vertex's and its edge's half-plane's fits over the frame. */
  struct Motion {
    Quadratic x, y;        ///< The projected vertex.
    Quadratic nx, ny, off; ///< The half-plane, on a convex face.
  };

  /// A direction scaled to unit length, or zero left zero.
  HS_HOT_INLINE static math::Vector unit(const math::Vector &v) {
    const float length2 = math::dot(v, v);
    return length2 > 0.0f ? v * (1.0f / sqrtf(length2)) : v;
  }

  /// Projects a direction, any length, about the axis as SDF::PreparedFace does.
  HS_FLASH_INLINE void project(const math::Vector &v, float &px,
                               float &py) const {
    float depth = math::dot(v, axis);
    if (fabsf(depth) < math::TOLERANCE)
      depth = copysignf(math::TOLERANCE, depth);
    const float inv = 1.0f / depth;
    px = math::dot(v, basis_u) * inv;
    py = math::dot(v, basis_w) * inv;
  }

  /// Fills the edges from the projected vertices, with twice the polygon's
  /// signed area and its squared circumradius about the axis.
  HS_FLASH_INLINE void stand_edges(float &area2, float &radius2) {
    for (int k = 0; k < count; ++k) {
      Edge &e = edges[k];
      const Edge &next = edges[k + 1 == count ? 0 : k + 1];
      e.ex = next.x - e.x;
      e.ey = next.y - e.y;
      const float length2 = e.ex * e.ex + e.ey * e.ey;
      e.inv_length2 = length2 > 0.0f ? 1.0f / length2 : 0.0f;
      area2 += e.x * next.y - next.x * e.y;
      radius2 = fmaxf(radius2, e.x * e.x + e.y * e.y);
      planes[k] = {0.0f, 0.0f, 0.0f};
    }
  }

  /**
   * @brief Finds whether the standing polygon is convex, as SDF::PreparedFace finds
   *        it, and its half-planes if so.
   * @param area2 Twice its signed area.
   * @return Whether it is.
   */
  HS_FLASH_INLINE bool half_planes(float area2) {
    bool pos = false, neg = false;
    int flips = 0;
    float dir = 0.0f;
    const Edge *previous = &edges[count - 1];
    for (int k = 0; k < count; ++k) {
      const Edge &e = edges[k];
      const float turn = previous->ex * e.ey - previous->ey * e.ex;
      const float scale =
          (previous->ex * previous->ex + previous->ey * previous->ey) *
          (e.ex * e.ex + e.ey * e.ey);
      previous = &e;
      SDF::scan_polygon_turn(turn, scale, pos, neg);
      if (e.ey != 0.0f) {
        flips += dir * e.ey < 0.0f;
        dir = e.ey;
      }
    }
    if ((pos && neg) || flips < 1 || flips > 2)
      return false;
    const float sign = area2 >= 0.0f ? 1.0f : -1.0f;
    float centre_side = -FLT_MAX;
    for (int k = 0; k < count; ++k) {
      const Edge &e = edges[k];
      Plane &h = planes[k];
      const float length2 = e.ex * e.ex + e.ey * e.ey;
      if (length2 < 1e-12f)
        return false;
      const float inv = sign / sqrtf(length2);
      h.nx = e.ey * inv;
      h.ny = -e.ex * inv;
      h.off = -(h.nx * e.x + h.ny * e.y);
      centre_side = fmaxf(centre_side, h.off);
    }
    return centre_side < 0.0f;
  }

  /// Squared distance in the projection plane to an edge.
  HS_HOT_INLINE static float edge_squared(const Edge &e, float px, float py) {
    const float wx = px - e.x, wy = py - e.y;
    const float s =
        hs::clamp((wx * e.ex + wy * e.ey) * e.inv_length2, 0.0f, 1.0f);
    const float bx = wx - e.ex * s, by = wy - e.ey * s;
    return bx * bx + by * by;
  }

  /// Squared distance in the projection plane to the nearest edge.
  HS_FLASH_INLINE float nearest_squared(float px, float py) const {
    float best = FLT_MAX;
    for (int k = 0; k < count; ++k)
      best = fminf(best, edge_squared(edges[k], px, py));
    return best;
  }

  int count;               ///< Number of vertices.
  bool north;              ///< Whether it lies around the north pole.
  bool convex = false;     ///< Whether it probes by its half-planes.
  bool linear = true;      ///< Whether distances stay in plane units.
  Edge *edges;             ///< The vertices and edges.
  Plane *planes;           ///< A convex face's half-planes.
  Motion *motion;          ///< The fits.
  float *line2;            ///< Squared bound on each edge's distance from
                           ///< the column pair's pixels.
  uint8_t *order;          ///< The edges by line2, ascending.
  Span *span;              ///< Each edge's extent along the meridian's line.
  float along_x = 0.0f;    ///< The meridian's line's unit direction, x.
  float along_y = 0.0f;    ///< The meridian's line's unit direction, y.
  float *crossings[2];     ///< Cosines of the meridian's crossings, per arm.
  int crossed[2] = {0, 0}; ///< Crossings per arm.
  math::Vector axis;       ///< Projection axis.
  math::Vector basis_u;    ///< Projection x axis.
  math::Vector basis_w;    ///< Projection y axis.
  math::Vector azimuth;    ///< The pair's azimuth.
  math::Vector up;         ///< The pair's up.
  Quadratic size_fit;      ///< size()'s fit.
  float face_size = 0.0f;  ///< size().
  mutable int nearest = 0; ///< The edge nearest the last pixel probed.
};

} // namespace Scan
