/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <utility>
#include <array>
#include <functional>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <cfloat>
#include <span>
#include "geometry.h"
#include "color.h"
#include "constants.h"
#include "filter.h"
#include "concepts.h"
#include "static_circular_buffer.h"
#include "canvas.h"

namespace SDF {
/** @brief Vertical scanline bounds (min/max Y). */
struct Bounds {
  int y_min, y_max;
};
/** @brief Horizontal scanline interval (start/end X). */
struct Interval {
  int start, end;
};
/** @brief Result of a signed distance query. */
struct DistanceResult {
  float dist;        // Signed distance (negative inside)
  float t;           // Normalized parameter (0-1) or angle
  float raw_dist;    // Unsigned or supplementary distance
  float aux;         // Auxiliary value (e.g. barycentric coordinate)
  float size = 1.0f; // Size metric for normalization

  DistanceResult() = default;
  DistanceResult(float d, float t_val, float rd, float ax, float sz)
      : dist(d), t(t_val), raw_dist(rd), aux(ax), size(sz) {}
};

/**
 * @brief Calculates signed distance to a ring.
 * Returns:
 *  dist: Signed distance (negative inside)
 *  t: Normalized parameter (0-1) corresponding to angle/2PI
 *  raw_dist: Unsigned distance to centerline
 */
struct Ring {
  const Basis &basis;
  float radius;
  float thickness;
  float phase;

  Vector normal, u, w;
  float nx, ny, nz;
  float target_angle, center_phi;
  float cos_max, cos_min, cos_target, inv_sin_target;

  // Optional optimization for "Scan Full Row" check
  float r_val;
  float alpha_angle; /**< The azimuth angle of the normal vector in the XZ
                        plane. */
  static constexpr bool is_solid = false;

  Ring(const Basis &b, float r, float th, float ph = 0)
      : basis(b), radius(r), thickness(th), phase(ph) {
    normal = basis.v;
    u = basis.u;
    w = basis.w;
    nx = normal.i;
    ny = normal.j;
    nz = normal.k;

    target_angle = radius * (PI_F / 2.0f);
    center_phi = acosf(ny);

    float ang_min = std::max(0.0f, target_angle - thickness);
    float ang_max = std::min(PI_F, target_angle + thickness);
    cos_max = cosf(ang_min);
    cos_min = cosf(ang_max);
    cos_target = cosf(target_angle);

    bool safe_approx = (target_angle > 0.05f && target_angle < PI_F - 0.05f);
    inv_sin_target = safe_approx ? (1.0f / sinf(target_angle)) : 0.0f;

    // For getHorizontalBounds
    r_val = sqrtf(nx * nx + nz * nz);
    alpha_angle = atan2f(nz, nx);
  }

  template <int H> Bounds get_vertical_bounds() const {
    constexpr int H_VIRT = H + hs::H_OFFSET;
    float a1 = center_phi - target_angle;
    float a2 = center_phi + target_angle;
    float phi_min = 0, phi_max = PI_F;

    if (a1 > 0) {
      float p1 = a1;
      float p2 = acosf(cosf(a2));
      phi_min = std::min(p1, p2);
    }
    if (a2 < PI_F) {
      float p1 = acosf(cosf(a1));
      float p2 = acosf(cosf(a2));
      phi_max = std::max(p1, p2);
    }

    float f_phi_min = std::max(0.0f, phi_min - thickness);
    float f_phi_max = std::min(PI_F, phi_max + thickness);

    int y_min = std::max(
        0, static_cast<int>(floorf((f_phi_min * (H_VIRT - 1)) / PI_F)));
    int y_max = std::min(
        H - 1, static_cast<int>(ceilf((f_phi_max * (H_VIRT - 1)) / PI_F)));
    return {y_min, y_max};
  }

  /**
   * @brief Computes the horizontal scanline intervals for this shape at a given
   * y-coordinate.
   * @param y The vertical pixel coordinate.
   * @param out Output iterator or callback accepting (float start, float end).
   * @return True if intervals were found and reported.
   */
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    constexpr int H_VIRT = H + hs::H_OFFSET;
    float phi = y_to_phi<H_VIRT>(static_cast<float>(y));
    float cos_phi = cosf(phi);
    float sin_phi = sinf(phi);

    if (r_val < 0.01f)
      return false;

    float denom = r_val * sin_phi;
    if (std::abs(denom) < 0.000001f)
      return false;

    float C_min = (cos_min - ny * cos_phi) / denom;
    float C_max = (cos_max - ny * cos_phi) / denom;
    float min_cos = std::max(-1.0f, C_min);
    float max_cos = std::min(1.0f, C_max);

    if (min_cos > max_cos)
      return true; // Empty row

    float angle_min = acosf(max_cos);
    float angle_max = acosf(min_cos);

    float pixel_width = 2.0f * PI_F / W;
    float safe_threshold = pixel_width;

    if (angle_min <= safe_threshold) {
      // Merge across alpha
      float f_x1 = (alpha_angle - angle_max) * W / (2 * PI_F);
      float f_x2 = (alpha_angle + angle_max) * W / (2 * PI_F);
      out(floorf(f_x1), ceilf(f_x2));
    } else if (angle_max >= PI_F - safe_threshold) {
      // Merge across antipode
      float f_x1 = (alpha_angle + angle_min) * W / (2 * PI_F);
      float f_x2 = (alpha_angle + 2 * PI_F - angle_min) * W / (2 * PI_F);
      out(floorf(f_x1), ceilf(f_x2));
    } else {
      // Merge across alpha
      float f_x1 = (alpha_angle - angle_max) * W / (2 * PI_F);
      float f_x2 = (alpha_angle - angle_min) * W / (2 * PI_F);
      float f_x3 = (alpha_angle + angle_min) * W / (2 * PI_F);
      float f_x4 = (alpha_angle + angle_max) * W / (2 * PI_F);

      out(floorf(f_x1), ceilf(f_x2));
      out(floorf(f_x3), ceilf(f_x4));
    }

    return true;
  }

  /**
   * @brief Computes signed distance to the ring.
   * @param p Point on sphere (normalized).
   * @return DistanceResult {dist, t, raw_dist}.
   */
  DistanceResult distance(const Vector &p) const {
    DistanceResult res;
    distance<true>(p, res);
    return res;
  }

  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    float d = dot(p, normal);
    if (d < cos_min || d > cos_max) {
      res = DistanceResult(100.0f, 0.0f, 100.0f, 0.0f, thickness);
      return;
    }

    float dist = 0;
    if (inv_sin_target != 0) {
      dist = std::abs(d - cos_target) * inv_sin_target;
    } else {
      float polar = acosf(std::max(-1.0f, std::min(1.0f, d)));
      dist = std::abs(polar - target_angle);
    }

    float t = 0.0f;
    if constexpr (ComputeUVs) {
      float dot_u = dot(p, u);
      float dot_w = dot(p, w);
      float azimuth = atan2f(dot_w, dot_u);
      if (azimuth < 0)
        azimuth += 2 * PI_F;
      azimuth += phase;
      t = azimuth / (2 * PI_F);
    }

    res = DistanceResult(dist - thickness, t, dist, 0.0f, thickness);
  }
};

/**
 * @brief Calculates signed distance to a distorted ring.
 * Returns:
 *  dist: Signed distance - thickness
 *  t: Normalized parameter (0-1) corresponding to angle/2PI
 *  raw_dist: Unsigned distance to centerline
 */
struct DistortedRing {
  const Basis &basis;
  float radius;
  float thickness;
  ScalarFn shift_fn;
  float max_distortion;
  float phase;

  Vector normal, u, w;
  float nx, ny, nz;
  float target_angle, center_phi;
  float max_thickness;

  // Optimization
  float r_val;
  float alpha_angle;
  float cos_max_limit, cos_min_limit;
  static constexpr bool is_solid = false;

  DistortedRing(const Basis &b, float r, float th, ScalarFn sf, float md,
                float ph)
      : basis(b), radius(r), thickness(th), shift_fn(sf), max_distortion(md),
        phase(ph) {
    normal = basis.v;
    u = basis.u;
    w = basis.w;
    nx = normal.i;
    ny = normal.j;
    nz = normal.k;
    target_angle = radius * (PI_F / 2.0f);
    center_phi = acosf(ny);
    max_thickness = thickness + max_distortion;

    r_val = sqrtf(nx * nx + nz * nz);
    alpha_angle = atan2f(nz, nx);

    float ang_min = std::max(0.0f, target_angle - max_thickness);
    float ang_max = std::min(PI_F, target_angle + max_thickness);
    cos_max_limit = cosf(ang_min);
    cos_min_limit = cosf(ang_max);
  }

  template <int H> Bounds get_vertical_bounds() const {
    constexpr int H_VIRT = H + hs::H_OFFSET;
    float a1 = center_phi - target_angle;
    float a2 = center_phi + target_angle;
    float phi_min = 0, phi_max = PI_F;

    if (a1 > 0) {
      float p1 = acosf(cosf(a1));
      float p2 = acosf(cosf(a2));
      phi_min = std::min(p1, p2);
    }
    if (a2 < PI_F) {
      float p1 = acosf(cosf(a1));
      float p2 = acosf(cosf(a2));
      phi_max = std::max(p1, p2);
    }

    float margin = max_thickness + 0.1f;
    float f_phi_min = std::max(0.0f, phi_min - margin);
    float f_phi_max = std::min(PI_F, phi_max + margin);

    int y_min = std::max(
        0, static_cast<int>(floorf((f_phi_min * (H_VIRT - 1)) / PI_F)));
    int y_max = std::min(
        H - 1, static_cast<int>(ceilf((f_phi_max * (H_VIRT - 1)) / PI_F)));
    return {y_min, y_max};
  }

  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    constexpr int H_VIRT = H + hs::H_OFFSET;
    float phi = y_to_phi<H_VIRT>(static_cast<float>(y));
    float cos_phi = cosf(phi);
    float sin_phi = sinf(phi);

    if (r_val < 0.01f)
      return false;

    float denom = r_val * sin_phi;
    if (std::abs(denom) < 0.000001f)
      return false;

    float C_min = (cos_min_limit - ny * cos_phi) / denom;
    float C_max = (cos_max_limit - ny * cos_phi) / denom;
    float min_cos = std::max(-1.0f, C_min);
    float max_cos = std::min(1.0f, C_max);

    if (min_cos > max_cos)
      return true; // Empty row

    float angle_min = acosf(max_cos);
    float angle_max = acosf(min_cos);

    float pixel_width = 2.0f * PI_F / W;
    float safe_threshold = pixel_width;

    if (angle_min <= safe_threshold) {
      float f_x1 = (alpha_angle - angle_max) * W / (2 * PI_F);
      float f_x2 = (alpha_angle + angle_max) * W / (2 * PI_F);
      out(floorf(f_x1), ceilf(f_x2));
    } else if (angle_max >= PI_F - safe_threshold) {
      float f_x1 = (alpha_angle + angle_min) * W / (2 * PI_F);
      float f_x2 = (alpha_angle + 2 * PI_F - angle_min) * W / (2 * PI_F);
      out(floorf(f_x1), ceilf(f_x2));
    } else {
      float f_x1 = (alpha_angle - angle_max) * W / (2 * PI_F);
      float f_x2 = (alpha_angle - angle_min) * W / (2 * PI_F);
      float f_x3 = (alpha_angle + angle_min) * W / (2 * PI_F);
      float f_x4 = (alpha_angle + angle_max) * W / (2 * PI_F);
      out(floorf(f_x1), ceilf(f_x2));
      out(floorf(f_x3), ceilf(f_x4));
    }
    return true;
  }

  /**
   * @brief Computes signed distance to the distorted ring.
   * @param p Point on sphere (normalized).
   * @return DistanceResult {dist, t, raw_dist}.
   */
  DistanceResult distance(const Vector &p) const {
    DistanceResult res;
    distance<true>(p, res);
    return res;
  }

  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    float polar = angle_between(p, normal);

    float t_norm = 0.0f;
    float shift = 0.0f;

    if constexpr (ComputeUVs) {
      float dot_u = dot(p, u);
      float dot_w = dot(p, w);
      float azimuth = atan2f(dot_w, dot_u);
      if (azimuth < 0)
        azimuth += 2 * PI_F;

      float t_val = azimuth + phase;
      t_norm = t_val / (2 * PI_F);
      shift = shift_fn(t_norm);
    } else {
      float dot_u = dot(p, u);
      float dot_w = dot(p, w);
      float azimuth = atan2f(dot_w, dot_u);
      if (azimuth < 0)
        azimuth += 2 * PI_F;

      float t_val = azimuth + phase;
      t_norm = t_val / (2 * PI_F);
      shift = shift_fn(t_norm);
    }

    float local_target = target_angle + shift;
    float dist = std::abs(polar - local_target);

    res = DistanceResult(dist - thickness, t_norm, dist, 0.0f, thickness);
  }
};

/**
 * @brief CSG Union operation (A + B).
 * Takes the minimum distance of two shapes.
 */
template <typename A, typename B> struct Union {
  const A &a;
  const B &b;
  float thickness;
  static constexpr bool is_solid = true;

  Union(const A &shapeA, const B &shapeB)
      : a(shapeA), b(shapeB),
        thickness(std::max(shapeA.thickness, shapeB.thickness)) {}

  template <int H> Bounds get_vertical_bounds() const {
    auto b1 = a.template get_vertical_bounds<H>();
    auto b2 = b.template get_vertical_bounds<H>();
    return {std::min(b1.y_min, b2.y_min), std::max(b1.y_max, b2.y_max)};
  }

  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    // TODO: Union intervals are complex; fallback to full scan.
    return false;
  }

  /**
   * @brief Signed distance to Union.
   */
  DistanceResult distance(const Vector &p) const {
    DistanceResult res;
    distance<true>(p, res);
    return res;
  }

  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    a.template distance<ComputeUVs>(p, res);
    DistanceResult resB;
    b.template distance<ComputeUVs>(p, resB);
    if (res.dist < resB.dist)
      return; // Min distance
    res = resB;
  }
};

/**
 * @brief CSG Subtraction operation (A - B).
 * Takes Max(A, -B).
 */
template <typename A, typename B> struct Subtract {
  const A &a;
  const B &b;
  float thickness;
  static constexpr bool is_solid = true;

  Subtract(const A &shapeA, const B &shapeB)
      : a(shapeA), b(shapeB), thickness(shapeA.thickness) {}

  template <int H> Bounds get_vertical_bounds() const {
    return a.template get_vertical_bounds<H>(); // Subtraction is bounded by A
  }

  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    // Subtraction intervals delegate to A (conservative)
    return a.template get_horizontal_intervals<W, H>(y, out);
  }

  /**
   * @brief Signed distance to Subtraction (A - B).
   */
  DistanceResult distance(const Vector &p) const {
    DistanceResult res;
    distance<true>(p, res);
    return res;
  }

  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    a.template distance<ComputeUVs>(p, res);
    DistanceResult resB;
    b.template distance<ComputeUVs>(p, resB);
    // Max(A, -B)
    if (-resB.dist > res.dist) {
      resB.dist = -resB.dist;
      res = resB;
    }
  }
};

/**
 * @brief CSG Intersection operation (A & B).
 * Takes the maximum distance of two shapes.
 */
template <typename A, typename B> struct Intersection {
  const A &a;
  const B &b;
  float thickness;
  static constexpr bool is_solid = true;

  Intersection(const A &shapeA, const B &shapeB)
      : a(shapeA), b(shapeB),
        thickness(std::min(shapeA.thickness, shapeB.thickness)) {}

  template <int H> Bounds get_vertical_bounds() const {
    auto b1 = a.template get_vertical_bounds<H>();
    auto b2 = b.template get_vertical_bounds<H>();
    return {std::max(b1.y_min, b2.y_min), std::min(b1.y_max, b2.y_max)};
  }

  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    StaticCircularBuffer<std::pair<float, float>, 32> intervalsA;
    StaticCircularBuffer<std::pair<float, float>, 32> intervalsB;

    bool hasA = a.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) { intervalsA.push_back({start, end}); });

    bool hasB = b.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) { intervalsB.push_back({start, end}); });

    if (!hasA) {
      return b.template get_horizontal_intervals<W, H>(y, out);
    }
    if (!hasB) {
      return a.template get_horizontal_intervals<W, H>(y, out);
    }

    if (intervalsA.is_empty() || intervalsB.is_empty())
      return true;

    // Intersect sorted intervals
    size_t idxA = 0;
    size_t idxB = 0;

    bool found = false;
    while (idxA < intervalsA.size() && idxB < intervalsB.size()) {
      auto ivA = intervalsA[idxA];
      auto ivB = intervalsB[idxB];

      float start = std::max(ivA.first, ivB.first);
      float end = std::min(ivA.second, ivB.second);

      if (start < end) {
        out(start, end);
        found = true;
      }

      if (ivA.second < ivB.second) {
        idxA++;
      } else {
        idxB++;
      }
    }

    return true;
  }

  /**
   * @brief Signed distance to Intersection.
   */
  DistanceResult distance(const Vector &p) const {
    DistanceResult res;
    distance<true>(p, res);
    return res;
  }

  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    a.template distance<ComputeUVs>(p, res);
    DistanceResult resB;
    b.template distance<ComputeUVs>(p, resB);
    // Max(A, B)
    if (res.dist > resB.dist)
      return;
    res = resB;
  }
};

/**
 * @brief Scratch buffer for Face computations to avoid allocations.
 */
struct FaceScratchBuffer {
  static constexpr int MAX_VERTS = 64;
  std::array<Vector, MAX_VERTS + 1> poly2D; // +1 to avoid modulo
  std::array<Vector, MAX_VERTS> edgeVectors;
  std::array<float, MAX_VERTS> edgeLengthsSq;
  std::array<Vector, MAX_VERTS> planes;
  std::array<std::pair<float, float>, 4> intervals;
  std::array<float, MAX_VERTS> thetas;
  std::array<float, MAX_VERTS> invEdgeLengthsSq;
  std::array<float, MAX_VERTS> invEdgeJ;
};

/**
 * @brief Represents a planar face for SDF rendering.
 * Computes 2D projection and vertical/horizontal bounds for optimization.
 */
struct Face {
  Vector center;
  Vector basisV, basisU, basisW;
  int count;
  float thickness;
  float size;
  float max_r2 = 0.0f;
  float radius = 0.0f;
  float max_dist = 0.0f;
  float max_dist_sq = 0.0f;

  std::span<Vector> poly2D;
  std::span<Vector> edgeVectors;
  std::span<float> edgeLengthsSq;
  std::span<Vector> planes;
  std::span<float> invEdgeLengthsSq;
  std::span<float> invEdgeJ;

  int y_min, y_max;
  std::span<std::pair<float, float>> intervals;
  bool full_width;
  static constexpr bool is_solid = true;

  Face(std::span<const Vector> vertices, std::span<const uint16_t> indices,
       float th, FaceScratchBuffer &scratch, int h_virt, int height)
      : thickness(th), full_width(true) {

    // Early Vertical Exit
    float min_y_val = 2.0f;
    float max_y_val = -2.0f;

    for (int idx : indices) {
      float y = vertices[idx].j;
      if (y < min_y_val)
        min_y_val = y;
      if (y > max_y_val)
        max_y_val = y;
    }

    float min_phi_check = acosf(hs::clamp(max_y_val, -1.0f, 1.0f));
    float max_phi_check = acosf(hs::clamp(min_y_val, -1.0f, 1.0f));
    float margin_check = thickness + 0.05f;

    int y_min_check = std::max(
        0, static_cast<int>(floorf(
               (std::max(0.0f, min_phi_check - margin_check) * (h_virt - 1)) /
               PI_F)));
    int y_max_check = std::min(
        height - 1,
        static_cast<int>(ceilf(
            (std::min(PI_F, max_phi_check + margin_check) * (h_virt - 1)) /
            PI_F)));

    if (y_min_check > y_max_check) {
      y_min = 1;
      y_max = 0;
      return;
    }

    count = indices.size();
    if (count > FaceScratchBuffer::MAX_VERTS)
      count = FaceScratchBuffer::MAX_VERTS;

    center = Vector(0, 0, 0);
    for (int i = 0; i < count; ++i)
      center = center + vertices[indices[i]];
    center.normalize();

    basisV = center;
    Vector ref = (std::abs(dot(center, X_AXIS)) > 0.9f) ? Y_AXIS : X_AXIS;
    basisU = cross(center, ref).normalize();
    basisW = cross(center, basisU).normalize();

    max_r2 = 0.0f;
    for (int i = 0; i < count; ++i) {
      const Vector &v = vertices[indices[i]];
      float d = dot(v, basisV);
      float px = dot(v, basisU) / d;
      float py = dot(v, basisW) / d;

      scratch.poly2D[i] = Vector(px, py, 0);

      float r2 = px * px + py * py;
      if (r2 > max_r2)
        max_r2 = r2;
    }
    radius = sqrtf(max_r2);
    max_dist = radius + 0.1f;
    max_dist_sq = max_dist * max_dist;

    scratch.poly2D[count] = scratch.poly2D[0];
    poly2D = std::span<Vector>(scratch.poly2D.data(), count + 1);

    float min_edge_dist = 1e9f;
    for (int i = 0; i < count; ++i) {
      Vector v1 = scratch.poly2D[i];
      Vector v2 = scratch.poly2D[(i + 1) % count];

      Vector edge = v2 - v1;
      float edge_len_sq = dot(edge, edge);
      float t = 0.0f;
      if (edge_len_sq > 1e-9f) {
        t = dot(-v1, edge) / edge_len_sq;
        t = std::max(0.0f, std::min(1.0f, t));
      }
      Vector closest = v1 + edge * t;
      float d_line = closest.magnitude();

      if (d_line < min_edge_dist)
        min_edge_dist = d_line;
    }
    size = (min_edge_dist > 1e8f) ? 1.0f : min_edge_dist;

    if (size < radius * 0.25f)
      size = radius * 0.25f;

    // Determine if this is a small face eligible for the fast path.
    // For small faces (vertex phi span < 0.3f ~17°), vertex-based bounds
    // are tight enough — skip expensive arc extrema computation.
    float vertex_phi_span = max_phi_check - min_phi_check;
    bool use_fast_bounds = (vertex_phi_span < 0.3f);

    int planes_count = 0;

    if (use_fast_bounds) {
      // FAST PATH: Skip plane normals and arc extrema entirely.
      // Use vertex-based phi bounds from the early check above.
      for (int i = 0; i < count; ++i) {
        Vector edge = scratch.poly2D[(i + 1) % count] - scratch.poly2D[i];
        scratch.edgeVectors[i] = edge;

        float edge_len_sq = dot(edge, edge);
        scratch.edgeLengthsSq[i] = edge_len_sq;
        scratch.invEdgeLengthsSq[i] =
            (edge_len_sq > 1e-12f) ? (1.0f / edge_len_sq) : 0.0f;
        scratch.invEdgeJ[i] =
            (std::abs(edge.j) > 1e-12f) ? (1.0f / edge.j) : 0.0f;

        float theta = atan2f(vertices[indices[i]].k, vertices[indices[i]].i);
        if (theta < 0)
          theta += 2 * PI_F;
        scratch.thetas[i] = theta;
      }

      float margin = thickness + 0.05f;
      y_min = std::max(
          0,
          static_cast<int>(floorf(
              (std::max(0.0f, min_phi_check - margin) * (h_virt - 1)) / PI_F)));
      y_max = std::min(
          height - 1,
          static_cast<int>(ceilf(
              (std::min(PI_F, max_phi_check + margin) * (h_virt - 1)) / PI_F)));
    } else {
      // SLOW PATH: Full arc extrema computation for large faces.
      float min_phi = 100.0f;
      float max_phi = -100.0f;

      for (int i = 0; i < count; ++i) {
        int idx1 = indices[i];
        int idx2 = indices[(i + 1) % count];
        const Vector &v1 = vertices[idx1];
        const Vector &v2 = vertices[idx2];

        Vector edge = scratch.poly2D[(i + 1) % count] - scratch.poly2D[i];
        scratch.edgeVectors[i] = edge;

        float edge_len_sq = dot(edge, edge);
        scratch.edgeLengthsSq[i] = edge_len_sq;
        scratch.invEdgeLengthsSq[i] =
            (edge_len_sq > 1e-12f) ? (1.0f / edge_len_sq) : 0.0f;
        scratch.invEdgeJ[i] =
            (std::abs(edge.j) > 1e-12f) ? (1.0f / edge.j) : 0.0f;

        Vector normal = cross(v1, v2);
        float lenSq = dot(normal, normal);
        if (lenSq > 1e-12f)
          scratch.planes[planes_count++] = normal.normalize();

        float phi_val = acosf(hs::clamp(v1.j, -1.0f, 1.0f));
        if (phi_val < min_phi)
          min_phi = phi_val;
        if (phi_val > max_phi)
          max_phi = phi_val;

        // Arc Extrema Logic
        if (planes_count > 0) {
          const Vector &n = scratch.planes[planes_count - 1];
          float ny = n.j;
          if (std::abs(ny) < 0.99999f) {
            float nx = n.i;
            float nz = n.k;
            float tx = -nx * ny;
            float ty = 1.0f - ny * ny;
            float tz = -nz * ny;
            float tLenSq = tx * tx + ty * ty + tz * tz;
            if (tLenSq > 1e-12f) {
              float invLen = 1.0f / sqrtf(tLenSq);
              float ptx = tx * invLen;
              float pty = ty * invLen;
              float ptz = tz * invLen;

              float cx1 = (v1.j * ptz - v1.k * pty) * nx +
                          (v1.k * ptx - v1.i * ptz) * ny +
                          (v1.i * pty - v1.j * ptx) * nz;
              float cx2 = (pty * v2.k - ptz * v2.j) * nx +
                          (ptz * v2.i - ptx * v2.k) * ny +
                          (ptx * v2.j - pty * v2.i) * nz;

              if (cx1 > 0 && cx2 > 0) {
                float phiTop = acosf(hs::clamp(pty, -1.0f, 1.0f));
                if (phiTop < min_phi)
                  min_phi = phiTop;
              }
              if (cx1 < 0 && cx2 < 0) {
                float phiBot = acosf(hs::clamp(-pty, -1.0f, 1.0f));
                if (phiBot > max_phi)
                  max_phi = phiBot;
              }
            }
          }
        }

        float theta = atan2f(v1.k, v1.i);
        if (theta < 0)
          theta += 2 * PI_F;
        scratch.thetas[i] = theta;
      }

      bool npInside = (planes_count > 0);
      bool spInside = (planes_count > 0);
      for (int pi = 0; pi < planes_count; ++pi) {
        float center_side = dot(center, scratch.planes[pi]);
        // North pole: dot((0,1,0), plane) = plane.j
        if ((scratch.planes[pi].j > 0) != (center_side > 0))
          npInside = false;
        // South pole: dot((0,-1,0), plane) = -plane.j
        if ((-scratch.planes[pi].j > 0) != (center_side > 0))
          spInside = false;
      }
      if (npInside)
        min_phi = 0.0f;
      if (spInside)
        max_phi = PI_F;

      float margin = thickness + 0.05f;
      y_min = std::max(
          0, static_cast<int>(floorf(
                 (std::max(0.0f, min_phi - margin) * (h_virt - 1)) / PI_F)));
      y_max = std::min(
          height - 1,
          static_cast<int>(
              ceilf((std::min(PI_F, max_phi + margin) * (h_virt - 1)) / PI_F)));
    }

    edgeVectors = std::span<Vector>(scratch.edgeVectors.data(), count);
    edgeLengthsSq = std::span<float>(scratch.edgeLengthsSq.data(), count);
    invEdgeLengthsSq = std::span<float>(scratch.invEdgeLengthsSq.data(), count);
    invEdgeJ = std::span<float>(scratch.invEdgeJ.data(), count);
    planes = std::span<Vector>(scratch.planes.data(), planes_count);

    std::sort(scratch.thetas.begin(), scratch.thetas.begin() + count);
    float maxGap = 0;
    float gapStart = 0;
    for (int i = 0; i < count; ++i) {
      float next = (i + 1 < count) ? scratch.thetas[i + 1]
                                   : (scratch.thetas[0] + 2 * PI_F);
      float diff = next - scratch.thetas[i];
      if (diff > maxGap) {
        maxGap = diff;
        gapStart = scratch.thetas[i];
      }
    }

    int interval_count = 0;
    if (maxGap > PI_F) {
      full_width = false;
      float startT = fmodf(gapStart + maxGap, 2 * PI_F);
      float endT = gapStart;

      if (startT <= endT) {
        scratch.intervals[interval_count++] = {startT, endT};
      } else {
        // PUSH IN STRICTLY ASCENDING ORDER
        scratch.intervals[interval_count++] = {0.0f, endT}; // Left side first
        scratch.intervals[interval_count++] = {startT,
                                               2 * PI_F}; // Right side second
      }
    } else {
      full_width = true;
    }
    intervals = std::span<std::pair<float, float>>(scratch.intervals.data(),
                                                   interval_count);

    // Pole containment: if a face encircles a pole, extend vertical bounds
    // to cover it. The fast path skips the slow-path plane-normal check,
    // so this unified check works for both paths.
    // North pole (0, 1, 0)
    if (center.j > 0.01f) {
      float inv_c = 1.0f / center.j;
      float ppx = basisU.j * inv_c;
      float ppy = basisW.j * inv_c;
      bool pole_inside = false;
      for (int i = 0; i < count; ++i) {
        if ((poly2D[i].j > ppy) != (poly2D[i + 1].j > ppy)) {
          float ix = poly2D[i].i +
                     (ppy - poly2D[i].j) * edgeVectors[i].i * invEdgeJ[i];
          if (ppx < ix)
            pole_inside = !pole_inside;
        }
      }
      if (pole_inside) {
        y_min = 0;
        full_width = true;
      }
    }
    // South pole (0, -1, 0)
    if (center.j < -0.01f) {
      float cos_sp = -center.j;
      float inv_c = 1.0f / cos_sp;
      float ppx = -basisU.j * inv_c;
      float ppy = -basisW.j * inv_c;
      bool pole_inside = false;
      for (int i = 0; i < count; ++i) {
        if ((poly2D[i].j > ppy) != (poly2D[i + 1].j > ppy)) {
          float ix = poly2D[i].i +
                     (ppy - poly2D[i].j) * edgeVectors[i].i * invEdgeJ[i];
          if (ppx < ix)
            pole_inside = !pole_inside;
        }
      }
      if (pole_inside) {
        y_max = height - 1;
        full_width = true;
      }
    }
  }

  template <int H> Bounds get_vertical_bounds() const { return {y_min, y_max}; }

  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    if (full_width)
      return false;
    for (const auto &iv : intervals) {
      float f_x1 = iv.first * W / (2 * PI_F);
      float f_x2 = iv.second * W / (2 * PI_F);
      out(floorf(f_x1), ceilf(f_x2));
    }
    return true;
  }

  DistanceResult distance(const Vector &p) const {
    DistanceResult res;
    distance<true>(p, res);
    return res;
  }

  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    float cos_angle = dot(p, center);
    if (cos_angle <= 0.01f) {
      res = DistanceResult(100.0f, 0.0f, 100.0f, 0.0f, size);
      return;
    }

    float inv_cos = 1.0f / cos_angle;
    float px = dot(p, basisU) * inv_cos;
    float py = dot(p, basisW) * inv_cos;

    float pR2 = px * px + py * py;
    if (pR2 > max_dist_sq) {
      res = DistanceResult(100.0f, 0.0f, 100.0f, 0.0f, size);
      return;
    }

    float d = FLT_MAX;
    bool inside = false;

    for (int i = 0; i < count; ++i) {
      const Vector &Vi = poly2D[i];
      const Vector &Vnext = poly2D[i + 1];
      const Vector &edge = edgeVectors[i];

      float wx = px - Vi.i;
      float wy = py - Vi.j;

      float t = (wx * edge.i + wy * edge.j) * invEdgeLengthsSq[i];
      float clampVal = hs::clamp(t, 0.0f, 1.0f);

      float bx = wx - edge.i * clampVal;
      float by = wy - edge.j * clampVal;
      float distSq = bx * bx + by * by;

      if (distSq < d)
        d = distSq;

      if ((Vi.j > py) != (Vnext.j > py)) {
        float intersect_x = Vi.i + (py - Vi.j) * edge.i * invEdgeJ[i];
        if (px < intersect_x) {
          inside = !inside;
        }
      }
    }

    float s = inside ? -1.0f : 1.0f;
    float plane_dist = s * sqrtf(d);
    float dist = plane_dist - thickness;

    res = DistanceResult(dist, 0.0f, plane_dist, 0.0f, size);
  }
};

/**
 * @brief Calculates signed distance to a planar polygon.
 * Returns:
 *  dist: Signed distance from edge (negative inside)
 *  t: Normalized polar distance (polar / thickness)
 *  raw_dist: Polar distance from center
 */
struct Polygon {
  const Basis &basis;
  float thickness;
  int sides;
  float phase;
  float apothem;
  float nx, ny, nz, R_val, alpha_angle;
  int y_min, y_max;
  const bool is_solid = true;

  Polygon(const Basis &b, float r, float th, int s, float ph, int h_virt,
          int height)
      : basis(b), thickness(th), sides(s), phase(ph) {
    apothem = thickness * cosf(PI_F / sides);
    nx = basis.v.i;
    ny = basis.v.j;
    nz = basis.v.k;
    R_val = sqrtf(nx * nx + nz * nz);
    alpha_angle = atan2f(nz, nx);

    float center_phi = acosf(std::max(-1.0f, std::min(1.0f, ny)));
    float margin = thickness + 0.1f;
    y_min = std::max(
        0, static_cast<int>(floorf(
               (std::max(0.0f, center_phi - margin) * (h_virt - 1)) / PI_F)));
    y_max = std::min(
        height - 1,
        static_cast<int>(ceilf(
            (std::min(PI_F, center_phi + margin) * (h_virt - 1)) / PI_F)));
  }

  template <int H> Bounds get_vertical_bounds() const { return {y_min, y_max}; }

  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    float phi = y_to_phi<H>(static_cast<float>(y));
    float cos_phi = cosf(phi);
    float sin_phi = sinf(phi);

    if (R_val < 0.01f)
      return false;

    float ang_high = thickness;
    float D_min = cosf(ang_high);

    float denom = R_val * sin_phi;
    if (std::abs(denom) < 0.000001f)
      return false;

    float C_min = (D_min - ny * cos_phi) / denom;

    if (C_min > 1.0f)
      return true; // Completely outside the cap
    if (C_min < -1.0f)
      return false; // Full scan fallback (simplification)

    float d_alpha = acosf(C_min);
    float f_x1 = (alpha_angle - d_alpha) * W / (2 * PI_F);
    float f_x2 = (alpha_angle + d_alpha) * W / (2 * PI_F);
    out(floorf(f_x1), ceilf(f_x2));
    return true;
  }

  /**
   * @brief Signed distance to Polygon edge.
   */
  DistanceResult distance(const Vector &p) const {
    DistanceResult res;
    distance<true>(p, res);
    return res;
  }

  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    float polar = angle_between(p, basis.v);

    float t_val = 0.0f;
    float dist_edge = 0.0f;

    if constexpr (ComputeUVs) {
      float dot_u = dot(p, basis.u);
      float dot_w = dot(p, basis.w);
      float azimuth = atan2f(dot_w, dot_u);
      if (azimuth < 0)
        azimuth += 2 * PI_F;
      azimuth += phase;

      float sector = 2 * PI_F / sides;
      float local = wrap(azimuth + sector / 2.0f, sector) - sector / 2.0f;

      dist_edge = polar * cosf(local) - apothem;
      t_val = polar / thickness;
    } else {
      // For Polygon, distance depends on angle (azimuth).
      // Like DistortedRing, we cannot skip azimuth calculation if we want
      // correct distance. So ComputeUVs here acts only as "don't need to return
      // valid t".
      float dot_u = dot(p, basis.u);
      float dot_w = dot(p, basis.w);
      float azimuth = atan2f(dot_w, dot_u);
      if (azimuth < 0)
        azimuth += 2 * PI_F;
      azimuth += phase;

      float sector = 2 * PI_F / sides;
      float local = wrap(azimuth + sector / 2.0f, sector) - sector / 2.0f;

      dist_edge = polar * cosf(local) - apothem;
      t_val = 0.0f;
    }

    res = DistanceResult(dist_edge, t_val, polar, 0.0f, apothem);
  }
};

/**
 * @brief Calculates signed distance to a star shape.
 * Returns:
 *  dist: Signed distance from edge (negative inside)
 *  t: Normalized polar distance (polar / thickness)
 *  raw_dist: Polar distance from center
 */
template <int W> struct Star {
  const Basis &basis;
  int sides;
  float phase;
  const bool is_solid = true;

  float nx, ny, planeD;
  float thickness;

  // Scan
  float scanNy, scanNx, scanNz, scanR, scanAlpha;
  int yMin, yMax;

  Star(const Basis &b, float radius, int s, float ph, int h_virt, int height)
      : basis(b), sides(s), phase(ph) {
    float outerRadius = radius * (PI_F / 2.0f);
    float innerRadius = outerRadius * 0.382f;
    float angleStep = PI_F / sides;

    float vT = outerRadius;
    float vVx = innerRadius * cosf(angleStep);
    float vVy = innerRadius * sinf(angleStep);

    float dx = vVx - vT;
    float dy = vVy;
    float len = sqrtf(dx * dx + dy * dy);
    nx = -dy / len;
    ny = dx / len;
    planeD = -(nx * vT);
    thickness = outerRadius;

    // Scan
    scanNx = basis.v.i;
    scanNy = basis.v.j;
    scanNz = basis.v.k;
    scanR = sqrtf(scanNx * scanNx + scanNz * scanNz);
    scanAlpha = atan2f(scanNz, scanNx);

    float centerPhi = acosf(std::max(-1.0f, std::min(1.0f, basis.v.j)));
    float margin = outerRadius + 0.1f;
    yMin = std::max(
        0, static_cast<int>(floorf(
               (std::max(0.0f, centerPhi - margin) * (height - 1)) / PI_F)));
    yMax = std::min(
        height - 1,
        static_cast<int>(
            ceilf((std::min(PI_F, centerPhi + margin) * (height - 1)) / PI_F)));
  }

  template <int H> Bounds get_vertical_bounds() const { return {yMin, yMax}; }

  template <int W_scan, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    // Bounding circle
    float phi = y_to_phi<H>(static_cast<float>(y));
    float cosPhi = cosf(phi);
    float sinPhi = sinf(phi);

    if (scanR < 0.01f)
      return false;

    float pixelWidth = 2.0f * PI_F / W_scan;
    float D_min = cosf(thickness + pixelWidth);
    float denom = scanR * sinPhi;

    if (std::abs(denom) < 0.000001f)
      return false;

    float C_min = (D_min - scanNy * cosPhi) / denom;

    if (C_min > 1.0f)
      return true;
    if (C_min < -1.0f)
      return false;

    float dAlpha = acosf(C_min);
    float t1 = scanAlpha - dAlpha;
    float t2 = scanAlpha + dAlpha;
    float f_x1 = (t1 * W) / (2 * PI_F);
    float f_x2 = (t2 * W) / (2 * PI_F);
    float x1 = floorf(f_x1);
    float x2 = ceilf(f_x2);

    // Check for full line wraparound case equivalent to "x2 - x1 >= W" in JS
    if (x2 - x1 >= W)
      return false;

    out(x1, x2);
    return true;
  }

  DistanceResult distance(const Vector &p) const {
    DistanceResult res;
    distance<true>(p, res);
    return res;
  }

  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    float scanDist = angle_between(p, basis.v);
    float dotU = dot(p, basis.u);
    float dotW = dot(p, basis.w);
    float azimuth = atan2f(dotW, dotU);
    if (azimuth < 0)
      azimuth += 2 * PI_F;

    azimuth += phase;

    float sectorAngle = 2 * PI_F / sides;
    float localAzimuth =
        wrap(azimuth + sectorAngle / 2.0f, sectorAngle) - sectorAngle / 2.0f;
    localAzimuth = std::abs(localAzimuth);

    float px = scanDist * cosf(localAzimuth);
    float py = scanDist * sinf(localAzimuth);

    float distToEdge = px * nx + py * ny + planeD;

    res = DistanceResult(-distToEdge, azimuth / (2 * PI_F), scanDist, 0.0f,
                         thickness);
  }
};

/**
 * @brief Calculates signed distance to a flower shape.
 * Returns:
 *  dist: Signed distance from flower edge
 *  t: Normalized scan distance (scan_dist / thickness)
 *  raw_dist: Scan distance from antipode
 */
struct Flower {
  const Basis &basis;
  int sides;
  float phase;
  float thickness;
  float apothem;
  Vector antipode;
  // Optimization
  float scan_nx, scan_ny, scan_nz, scan_R, scan_alpha;
  int y_min, y_max;
  const bool is_solid = true;

  Flower(const Basis &b, float radius, int s, float ph, int h_virt, int height)
      : basis(b), sides(s), phase(ph) {
    float outer = radius * (PI_F / 2.0f);
    apothem = PI_F - outer;
    thickness = outer;
    antipode = -basis.v;

    scan_nx = antipode.i;
    scan_ny = antipode.j;
    scan_nz = antipode.k;
    scan_R = sqrtf(scan_nx * scan_nx + scan_nz * scan_nz);
    scan_alpha = atan2f(scan_nz, scan_nx);

    float center_phi = acosf(std::max(-1.0f, std::min(1.0f, antipode.j)));
    float margin = thickness + 0.1f;
    y_min = std::max(
        0, static_cast<int>(floorf(
               (std::max(0.0f, center_phi - margin) * (h_virt - 1)) / PI_F)));
    y_max = std::min(
        height - 1,
        static_cast<int>(ceilf(
            (std::min(PI_F, center_phi + margin) * (h_virt - 1)) / PI_F)));
  }

  template <int H> Bounds get_vertical_bounds() const { return {y_min, y_max}; }

  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    float phi = y_to_phi<H>(static_cast<float>(y));
    float cos_phi = cosf(phi);
    float sin_phi = sinf(phi);

    if (scan_R < 0.01f)
      return false;

    float denom = scan_R * sin_phi;
    if (std::abs(denom) < 0.000001f)
      return false;

    float ang_max = thickness;
    float cos_limit = cosf(ang_max);

    float C_min = (cos_limit - scan_ny * cos_phi) / denom;
    float C_max = (1.0f - scan_ny * cos_phi) / denom;

    float min_cos = std::max(-1.0f, C_min);
    float max_cos = std::min(1.0f, C_max);

    if (min_cos > max_cos)
      return true;

    float angle_min = acosf(max_cos);
    float angle_max = acosf(min_cos);

    float f_x1 = (scan_alpha - angle_max) * W / (2 * PI_F);
    float f_x2 = (scan_alpha - angle_min) * W / (2 * PI_F);
    float f_x3 = (scan_alpha + angle_min) * W / (2 * PI_F);
    float f_x4 = (scan_alpha + angle_max) * W / (2 * PI_F);

    out(floorf(f_x1), ceilf(f_x2));
    out(floorf(f_x3), ceilf(f_x4));
    return true;
  }

  DistanceResult distance(const Vector &p) const {
    DistanceResult res;
    distance<true>(p, res);
    return res;
  }

  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    float scan_dist = angle_between(p, antipode);
    float polar = PI_F - scan_dist;

    float t_val = 0.0f;
    float dist_edge = 0.0f;

    if constexpr (ComputeUVs) {
      float dot_u = dot(p, basis.u);
      float dot_w = dot(p, basis.w);
      float azimuth = atan2f(dot_w, dot_u);
      if (azimuth < 0)
        azimuth += 2 * PI_F;
      azimuth += phase;

      float sector = 2 * PI_F / sides;
      float local = wrap(azimuth + sector / 2.0f, sector) - sector / 2.0f;

      dist_edge = polar * cosf(local) - apothem;
      t_val = scan_dist / thickness;
    } else {
      // Flower shape depends on angle.
      float dot_u = dot(p, basis.u);
      float dot_w = dot(p, basis.w);
      float azimuth = atan2f(dot_w, dot_u);
      if (azimuth < 0)
        azimuth += 2 * PI_F;
      azimuth += phase;

      float sector = 2 * PI_F / sides;
      float local = wrap(azimuth + sector / 2.0f, sector) - sector / 2.0f;

      dist_edge = polar * cosf(local) - apothem;
      t_val = 0.0f;
    }
    res = DistanceResult(-dist_edge, t_val, scan_dist, 0.0f, thickness);
  }
};

struct HarmonicBlob {
  int l;
  int m;
  float amplitude;
  Quaternion inv_q;
  HarmonicWaveFn harmonic_fn;
  const bool is_solid = true;

  HarmonicBlob(int l, int m, float amplitude, const Quaternion &orientation,
               HarmonicWaveFn harmonic_fn)
      : l(l), m(m), amplitude(amplitude), harmonic_fn(harmonic_fn) {
    inv_q = orientation.inverse();
  }

  template <int H> Bounds get_vertical_bounds() const { return {0, H - 1}; }

  DistanceResult distance(const Vector &p) const {
    DistanceResult res;
    distance<true>(p, res);
    return res;
  }

  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    // Transform world pixel position to local harmonic space
    Vector v = rotate(p, inv_q);
    float phi = acosf(std::max(-1.0f, std::min(1.0f, v.j)));
    float theta = atan2f(v.k, v.i);
    float harmonic_val = harmonic_fn(l, m, theta, phi);

    float lobe_radius = 1.0f + std::abs(harmonic_val) * amplitude;
    float d = 1.0f - lobe_radius;

    // HarmonicBlob needs theta/phi for shape, so we can't fully skip UVs
    // (theta/phi) because they ARE the shape. But we can skip the return value
    // 't' calculation if 'ComputeUVs' meant 'return texture coords' separately.
    // In this engine, 't' is often used for texture mapping.
    // Here we return 'tanhf(...)' as 't'.
    float t_val = 0.0f;
    if constexpr (ComputeUVs) {
      t_val = tanhf(std::abs(harmonic_val) * amplitude);
    }

    res = DistanceResult(d, t_val, harmonic_val, 0.0f, 1.0f);
  }

  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    // HarmonicBlob is complex, full scan is safest/easiest port
    return false;
  }
};

struct Line {
  Vector a, b;
  float thickness;

  Vector n;
  float len;
  const bool is_solid = false;

  Line(const Vector &start, const Vector &end, float th)
      : a(start), b(end), thickness(th) {
    len = angle_between(a, b);
    if (len < 1e-6f) {
      n = Vector(0, 0, 0);
    } else {
      n = cross(a, b).normalize();
    }
  }

  template <int H> Bounds get_vertical_bounds() const { return {0, H - 1}; }

  DistanceResult distance(const Vector &p) const {
    DistanceResult res;
    distance<true>(p, res);
    return res;
  }

  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    if (len < 1e-6f) {
      float dist = angle_between(p, a);
      res = DistanceResult(dist - thickness, 0.0f, dist, 0.0f, thickness);
      return;
    }

    float d_plane = dot(p, n);
    Vector p_proj_plane = p - n * d_plane;
    float proj_mag = p_proj_plane.magnitude();

    if (proj_mag < 1e-6f) {
      float dA = angle_between(p, a);
      float dB = angle_between(p, b);
      float dist = std::min(dA, dB);
      res = DistanceResult(dist - thickness, 0.0f, dist, 0.0f, thickness);
      return;
    }

    Vector p_proj = p_proj_plane / proj_mag;
    float ang_a = angle_between(a, p_proj);
    float ang_b = angle_between(b, p_proj);

    float dist_seg = 0.0f;
    if (ang_a + ang_b <= len + 1e-4f) {
      dist_seg = asinf(std::abs(d_plane));
    } else {
      float dA = angle_between(p, a);
      float dB = angle_between(p, b);
      dist_seg = std::min(dA, dB);
    }

    res = DistanceResult(dist_seg - thickness, 0.0f, dist_seg, 0.0f, thickness);
  }

  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    return false;
  }
};
} // namespace SDF

/**
 * @brief The Scan struct contains volumetric (raster) drawing primitives.
 *
 * General Register Mapping for Scan Primitives:
 *  v0: Normalized parameter t (0-1) or angle
 *  v1: Raw Distance or Supplementary Value
 */
namespace Scan {

/**
 * @brief Processes a single pixel for rasterization.
 */
template <int W, int H, bool ComputeUVs = true>
static void process_pixel(int x, int y, const Vector &p, PipelineRef pipeline,
                          Canvas &canvas, const auto &shape,
                          FragmentShaderFn fragment_shader, bool debug_bb,
                          SDF::DistanceResult &result_scratch,
                          Fragment &frag_scratch) {
  shape.template distance<ComputeUVs>(p, result_scratch);
  float d = result_scratch.dist;

  float pixel_width = 2.0f * PI_F / W;

  // Determine solidity directly (compiler enforces existence)
  bool is_solid = shape.is_solid;

  // Threshold: Solids allow 1px AA bleed; Strokes (non-solid) must be strictly
  // inside (d <= 0)
  float threshold = is_solid ? pixel_width : 0.0f;

  if (debug_bb || d < threshold) {
    float alpha = 1.0f;

    if (is_solid) {
      if (d <= -pixel_width) {
        // FAST PATH: Pixel is fully inside the shape. Skip all AA math!
        alpha = 1.0f;
      } else {
        // SLOW PATH: Pixel is on the boundary. Calculate AA blending.
        float t_aa = 0.5f - d / (2.0f * pixel_width);
        alpha = quintic_kernel(std::max(0.0f, std::min(1.0f, t_aa)));
      }
    } else {
      // Stroke Falloff: Opacity fades over the entire thickness
      if constexpr (requires { shape.thickness; }) {
        if (shape.thickness > 0) {
          alpha = quintic_kernel(-d / shape.thickness);
        }
      }
    }

    if (!debug_bb && alpha <= 0.001f)
      return;

    frag_scratch.pos = p;
    frag_scratch.v0 = result_scratch.t;
    frag_scratch.v1 = result_scratch.raw_dist;
    frag_scratch.v2 = 0.0f;
    frag_scratch.v3 = result_scratch.aux;
    frag_scratch.size = result_scratch.size;
    frag_scratch.age = 0;

    fragment_shader(p, frag_scratch);

    if (debug_bb) {
      frag_scratch.color.color = frag_scratch.color.color.lerp16(
          Pixel(65535, 65535, 65535), 65535 / 2);
      frag_scratch.color.alpha = 1.0f;
      alpha = 1.0f;
    }

    if (frag_scratch.color.alpha > 0.001f) {
      pipeline.plot(canvas, x, y, frag_scratch.color.color, frag_scratch.age,
                    frag_scratch.color.alpha * alpha, frag_scratch.blend);
    }
  }
}

/**
 * @brief Main rasterization routine for SDF shapes.
 * Scans the bounding box, computes intervals, and executes the shader for valid
 * pixels.
 */
template <int W, int H, bool ComputeUVs = true>
static void rasterize(PipelineRef pipeline, Canvas &canvas, const auto &shape,
                      FragmentShaderFn fragment_shader, bool debug_bb = false) {
  bool effective_debug = debug_bb || canvas.debug();
  auto bounds = shape.template get_vertical_bounds<H>();

  if (!PixelLUT<W, H>::initialized)
    PixelLUT<W, H>::init();

  StaticCircularBuffer<std::pair<float, float>, 32> intervals;
  SDF::DistanceResult result_scratch;
  Fragment frag_scratch;

  for (int y = bounds.y_min; y <= bounds.y_max; ++y) {
    const Vector *row_vectors = &PixelLUT<W, H>::data[y * W];

    bool handled = shape.template get_horizontal_intervals<W, H>(
        y, [&](float t1, float t2) { intervals.push_back({t1, t2}); });

    if (handled && !intervals.is_empty()) {
      float current_end = -FLT_MAX;
      for (const auto &iv : intervals) {
        if (iv.second <= current_end)
          continue;

        float start = std::max(iv.first, current_end);
        float end = iv.second;
        current_end = end;

        if (end - start >= W) {
          for (int x = 0; x < W; ++x) {
            process_pixel<W, H, ComputeUVs>(
                x, y, row_vectors[x], pipeline, canvas, shape, fragment_shader,
                effective_debug, result_scratch, frag_scratch);
          }
          continue;
        }

        int x1 = static_cast<int>(floorf(start));
        int x2 = static_cast<int>(ceilf(end));
        if (x1 == x2)
          x2++;

        // FAST MODULO LOGIC
        int wx = x1 % W;
        if (wx < 0)
          wx += W;

        for (int x = x1; x < x2; ++x) {
          process_pixel<W, H, ComputeUVs>(
              wx, y, row_vectors[wx], pipeline, canvas, shape, fragment_shader,
              effective_debug, result_scratch, frag_scratch);

          // Modulo-free coordinate wrap!
          wx++;
          if (wx >= W)
            wx -= W;
        }
      }
      intervals.clear();
    } else {
      if (!handled) {
        for (int x = 0; x < W; ++x) {
          process_pixel<W, H, ComputeUVs>(
              x, y, row_vectors[x], pipeline, canvas, shape, fragment_shader,
              effective_debug, result_scratch, frag_scratch);
        }
      }
    }
  }
}

struct DistortedRing {
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, float thickness, ScalarFn shift_fn,
                   float amplitude, FragmentShaderFn fragment_shader,
                   float phase = 0, bool debug_bb = false) {
    SDF::DistortedRing shape(basis, radius, thickness, shift_fn, amplitude,
                             phase);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }
};

struct PlanarPolygon {
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, int sides, FragmentShaderFn fragment_shader,
                   float phase = 0, bool debug_bb = false) {
    auto res = get_antipode(basis, radius);
    float thickness = res.second * (PI_F / 2.0f);

    SDF::Polygon shape(res.first, res.second, thickness, sides, phase, H + 3,
                       H);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }
};

struct Line {
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Vector &v1,
                   const Vector &v2, float thickness,
                   FragmentShaderFn fragment_shader, bool debug_bb = false) {
    SDF::Line shape(v1, v2, thickness);
    Scan::rasterize<W, H>(pipeline, canvas, shape, fragment_shader, debug_bb);
  }
};

struct Ring {
  /**
   * @brief Draws a solid ring using SDF rasterization.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, float thickness,
                   FragmentShaderFn fragment_shader, float phase = 0,
                   bool debug_bb = false) {
    auto res = get_antipode(basis, radius);
    SDF::Ring shape(res.first, res.second, thickness, phase);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }

  // Overload for Vector normal inputs
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Vector &normal,
                   float radius, float thickness,
                   FragmentShaderFn fragment_shader, float phase = 0,
                   bool debug_bb = false) {
    Basis basis = make_basis(Quaternion(), normal);
    draw<W, H, ComputeUVs>(pipeline, canvas, basis, radius, thickness,
                           fragment_shader, phase, debug_bb);
  }
};

struct Circle {
  /**
   * @brief Draws a solid circle (filled ring).
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, FragmentShaderFn fragment_shader,
                   bool debug_bb = false) {
    // Circle is a Ring with inner radius 0
    float th = radius * (PI_F / 2.0f);
    Ring::draw<W, H, ComputeUVs>(pipeline, canvas, basis, 0.0f, th,
                                 fragment_shader, 0, debug_bb);
  }

  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Vector &normal,
                   float radius, FragmentShaderFn fragment_shader,
                   bool debug_bb = false) {
    Basis basis = make_basis(Quaternion(), normal);
    draw<W, H, ComputeUVs>(pipeline, canvas, basis, radius, fragment_shader,
                           debug_bb);
  }
};

struct Point {
  /**
   * @brief Draws a solid point (thick circle).
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Vector &p,
                   float thickness, FragmentShaderFn fragment_shader) {
    // Scan.Point uses Scan.Ring with radius 0
    Basis basis = make_basis(Quaternion(), p);
    Ring::draw<W, H>(pipeline, canvas, basis, 0.0f, thickness, fragment_shader);
  }
};

struct Star {
  /**
   * @brief Draws a solid star shape.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, int sides, FragmentShaderFn fragment_shader,
                   float phase = 0, bool debug_bb = false) {
    auto res = get_antipode(basis, radius);
    SDF::Star<W> shape(res.first, res.second, sides, phase, H, H);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }
};

struct Flower {
  /**
   * @brief Draws a solid flower shape.
   */
  template <int W, int H, bool ComputeUVs = true>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, int sides, FragmentShaderFn fragment_shader,
                   float phase = 0, bool debug_bb = false) {
    auto res = get_antipode(basis, radius);
    SDF::Flower shape(res.first, res.second, sides, phase, H + 3, H);
    Scan::rasterize<W, H, ComputeUVs>(pipeline, canvas, shape, fragment_shader,
                                      debug_bb);
  }
};

struct SphericalPolygon {
  /**
   * @brief Draws a solid spherical polygon.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, int sides, FragmentShaderFn fragment_shader,
                   float phase = 0, bool debug_bb = false) {
    // 1. Get Antipode and Basis
    Basis eff_basis = basis;
    float effective_radius = radius;
    if (radius > 1.0f) {
      eff_basis.v = -eff_basis.v;
      eff_basis.u = -eff_basis.u;
      effective_radius = 2.0f - radius;
    }

    // 2. Sample Points
    float offset = PI_F / sides;
    float theta_eq = effective_radius * (PI_F / 2.0f);
    float r = sinf(theta_eq);
    float d = cosf(theta_eq);
    float step = (2.0f * PI_F) / sides;

    std::array<Vector, 32> vertices;

    for (int i = 0; i < sides; ++i) {
      float theta = i * step + phase + offset;
      float cos_t = cosf(theta);
      float sin_t = sinf(theta);

      // pos = u * cosT + w * sinT
      // pos *= r
      // pos += v * d
      Vector pos = eff_basis.u * cos_t + eff_basis.w * sin_t;
      pos = pos * r + eff_basis.v * d;
      pos.normalize();
      vertices[i] = pos;
    }

    // 3. Create Indices (0, 1, ..., sides-1)
    std::array<uint16_t, 32> indices;
    for (uint16_t i = 0; i < sides; ++i) {
      indices[i] = i;
    }

    // 4. Create SDF::Face
    SDF::FaceScratchBuffer scratch;
    SDF::Face shape(std::span<const Vector>(vertices.data(), sides),
                    std::span<const uint16_t>(indices.data(), sides), 0.0f,
                    scratch, H + 3, H);

    Scan::rasterize<W, H>(pipeline, canvas, shape, fragment_shader, debug_bb);
  }
};

struct HarmonicBlob {
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, int l, int m,
                   float amplitude, const Quaternion &orientation,
                   HarmonicWaveFn harmonic_fn, FragmentShaderFn fragment_shader,
                   bool debug_bb = false) {
    SDF::HarmonicBlob shape(l, m, amplitude, orientation, harmonic_fn);
    Scan::rasterize<W, H>(pipeline, canvas, shape, fragment_shader, debug_bb);
  }
};

struct Mesh {
  template <int W, int H, typename MeshT>
  static void draw(PipelineRef pipeline, Canvas &canvas, const MeshT &mesh,
                   FragmentShaderFn fragment_shader, bool debug_bb = false) {
    SDF::FaceScratchBuffer scratch;

    size_t idx_offset = 0;
    for (size_t i = 0; i < mesh.num_faces; ++i) {
      size_t count = mesh.face_counts[i];
      std::span<const Vector> verts(mesh.vertices.data(), mesh.num_vertices);
      std::span<const uint16_t> indices(&mesh.faces[idx_offset], count);

      SDF::Face shape(verts, indices, 0.0f, scratch, H + 3, H);
      idx_offset += count;

      auto wrapper = [&](const Vector &p, Fragment &f_in) {
        f_in.v2 = static_cast<float>(i);
        fragment_shader(p, f_in);
      };
      Scan::rasterize<W, H, true>(pipeline, canvas, shape, wrapper, debug_bb);
    }
  }

  // Overload for MeshState
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const MeshState &mesh,
                   FragmentShaderFn fragment_shader, bool debug_bb = false) {
    SDF::FaceScratchBuffer scratch;

    const uint8_t *fc = mesh.get_face_counts_data();
    size_t num_f = mesh.get_face_counts_size();
    const uint16_t *fi = mesh.get_faces_data();
    const uint16_t *fo = mesh.get_face_offsets_data();

    for (size_t i = 0; i < num_f; ++i) {
      size_t count = fc[i];

      // Create Spans from MeshState
      std::span<const Vector> verts(mesh.vertices.data(), mesh.vertices.size());
      std::span<const uint16_t> indices(fi + fo[i], count);

      SDF::Face shape(verts, indices, 0.0f, scratch, H + 3, H);

      auto wrapper = [&](const Vector &p, Fragment &f_in) {
        f_in.v2 = static_cast<float>(i);
        fragment_shader(p, f_in);
      };

      Scan::rasterize<W, H, true>(pipeline, canvas, shape, wrapper, debug_bb);
    }
  }
};
} // namespace Scan
