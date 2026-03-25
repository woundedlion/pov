/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <array>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <cfloat>
#include <span>
#include "geometry.h"
#include "constants.h"
#include "concepts.h"
#include "static_circular_buffer.h"
#include "util.h"

namespace SDF {

// --- Rasterization constants ------------------------------------------------
/** Margin added to bounding boxes for AA kernel width (small faces). */
static constexpr float BOUNDS_MARGIN = 0.05f;
/** Expanded margin for shapes with variable thickness or large faces. */
static constexpr float BOUNDS_MARGIN_WIDE = 0.1f;
/** Minimum horizontal projection length for interval culling. */
static constexpr float MIN_HORIZONTAL_PROJ = 0.01f;
/** Epsilon for near-zero denominators in interval math. */
static constexpr float INTERVAL_DENOM_EPS = 1e-6f;
/** Threshold for near-pole ring approximation safety. */
static constexpr float POLE_SAFE_MARGIN = 0.05f;
/** Inner/outer radius ratio for star shapes (≈ golden ratio conjugate). */
static constexpr float STAR_INNER_RATIO = 0.382f;
/** Vertex phi-span threshold below which Face uses the fast vertex-only
 *  bounding path instead of the full arc-extrema computation. */
static constexpr float FAST_BOUNDS_PHI_THRESHOLD = 0.3f;
/** Minimum inradius-to-circumradius ratio used to floor Face::size,
 *  preventing degenerate near-zero inradii from collapsing AA. */
static constexpr float MIN_SIZE_RADIUS_RATIO = 0.25f;

/** Fold an angle into [0, π] (equivalent to acosf(cosf(x)) without trig). */
inline float clamp_phi(float x) {
  if (x < 0.0f)
    return -x;
  if (x > PI_F)
    return 2.0f * PI_F - x;
  return x;
}
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
    nx = normal.x;
    ny = normal.y;
    nz = normal.z;

    target_angle = radius * (PI_F / 2.0f);
    center_phi = fast_acos(ny);

    float ang_min = std::max(0.0f, target_angle - thickness);
    float ang_max = std::min(PI_F, target_angle + thickness);
    cos_max = fast_cosf(ang_min);
    cos_min = fast_cosf(ang_max);
    cos_target = fast_cosf(target_angle);

    bool safe_approx = (target_angle > POLE_SAFE_MARGIN &&
                        target_angle < PI_F - POLE_SAFE_MARGIN);
    inv_sin_target = safe_approx ? (1.0f / fast_sinf(target_angle)) : 0.0f;

    // For getHorizontalBounds
    r_val = sqrtf(nx * nx + nz * nz);
    alpha_angle = fast_atan2(nz, nx);
  }

  template <int H> Bounds get_vertical_bounds() const {
    constexpr int H_VIRT = H + hs::H_OFFSET;
    float a1 = center_phi - target_angle;
    float a2 = center_phi + target_angle;
    float phi_min = 0, phi_max = PI_F;

    if (a1 > 0) {
      float p1 = a1;
      float p2 = clamp_phi(a2);
      phi_min = std::min(p1, p2);
    }
    if (a2 < PI_F) {
      float p1 = clamp_phi(a1);
      float p2 = a2;
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
    if (!TrigLUT<W, H>::initialized) TrigLUT<W, H>::init();
    float cos_phi = TrigLUT<W, H>::cos_phi[y];
    float sin_phi = TrigLUT<W, H>::sin_phi[y];

    if (r_val < MIN_HORIZONTAL_PROJ)
      return false;

    float denom = r_val * sin_phi;
    if (std::abs(denom) < INTERVAL_DENOM_EPS)
      return false;

    float C_min = (cos_min - ny * cos_phi) / denom;
    float C_max = (cos_max - ny * cos_phi) / denom;
    float min_cos = std::max(-1.0f, C_min);
    float max_cos = std::min(1.0f, C_max);

    if (min_cos > max_cos)
      return true; // Empty row

    float angle_min = fast_acos(max_cos);
    float angle_max = fast_acos(min_cos);

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
      float polar = fast_acos(hs::clamp(d, -1.0f, 1.0f));
      dist = std::abs(polar - target_angle);
    }

    float t = 0.0f;
    if constexpr (ComputeUVs) {
      float dot_u = dot(p, u);
      float dot_w = dot(p, w);
      float azimuth = fast_atan2(dot_w, dot_u);
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
    nx = normal.x;
    ny = normal.y;
    nz = normal.z;
    target_angle = radius * (PI_F / 2.0f);
    center_phi = fast_acos(ny);
    max_thickness = thickness + max_distortion;

    r_val = sqrtf(nx * nx + nz * nz);
    alpha_angle = fast_atan2(nz, nx);

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
      float p1 = clamp_phi(a1);
      float p2 = clamp_phi(a2);
      phi_min = std::min(p1, p2);
    }
    if (a2 < PI_F) {
      float p1 = clamp_phi(a1);
      float p2 = clamp_phi(a2);
      phi_max = std::max(p1, p2);
    }

    float margin = max_thickness + BOUNDS_MARGIN_WIDE;
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
    if (!TrigLUT<W, H>::initialized) TrigLUT<W, H>::init();
    float cos_phi = TrigLUT<W, H>::cos_phi[y];
    float sin_phi = TrigLUT<W, H>::sin_phi[y];

    if (r_val < MIN_HORIZONTAL_PROJ)
      return false;

    float denom = r_val * sin_phi;
    if (std::abs(denom) < INTERVAL_DENOM_EPS)
      return false;

    float C_min = (cos_min_limit - ny * cos_phi) / denom;
    float C_max = (cos_max_limit - ny * cos_phi) / denom;
    float min_cos = std::max(-1.0f, C_min);
    float max_cos = std::min(1.0f, C_max);

    if (min_cos > max_cos)
      return true; // Empty row

    float angle_min = fast_acos(max_cos);
    float angle_max = fast_acos(min_cos);

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
    float d = dot(p, normal);
    // Early reject: outside bounding annulus
    if (d < cos_min_limit || d > cos_max_limit) {
      res = DistanceResult(100.0f, 0.0f, 100.0f, 0.0f, thickness);
      return;
    }
    float polar = fast_acos(hs::clamp(d, -1.0f, 1.0f));

    float dot_u = dot(p, u);
    float dot_w = dot(p, w);
    float azimuth = fast_atan2(dot_w, dot_u);
    if (azimuth < 0)
      azimuth += 2 * PI_F;

    float t_val = azimuth + phase;
    float t_norm = t_val / (2 * PI_F);
    float shift = shift_fn(t_norm);

    if constexpr (!ComputeUVs)
      t_norm = 0.0f;

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
    StaticCircularBuffer<std::pair<float, float>, 32> merged;

    bool hasA = a.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) { merged.push_back({start, end}); });

    bool hasB = b.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) { merged.push_back({start, end}); });

    // If neither shape provides intervals, fall back to full scan
    if (!hasA && !hasB)
      return false;

    // If only one shape returned intervals and the other fell back,
    // we must do full scan (the fallback shape needs all pixels)
    if (!hasA || !hasB)
      return false;

    if (merged.is_empty())
      return true; // Both reported intervals but none produced

    // Sort by start
    auto *data = &merged[0];
    size_t n = merged.size();
    for (size_t i = 1; i < n; ++i) {
      auto key = data[i];
      int j = static_cast<int>(i) - 1;
      while (j >= 0 && data[j].first > key.first) {
        data[j + 1] = data[j];
        --j;
      }
      data[j + 1] = key;
    }

    // Merge overlapping intervals
    float cur_start = data[0].first;
    float cur_end = data[0].second;
    for (size_t i = 1; i < n; ++i) {
      if (data[i].first <= cur_end) {
        cur_end = std::max(cur_end, data[i].second);
      } else {
        out(cur_start, cur_end);
        cur_start = data[i].first;
        cur_end = data[i].second;
      }
    }
    out(cur_start, cur_end);
    return true;
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
      return; // A is already closer, keep it
    res = resB;
  }
};

/**
 * @brief Smooth CSG Union operation using polynomial smooth minimum.
 * Shapes organically blend together within radius k (Inigo Quilez smin).
 */
template <typename A, typename B> struct SmoothUnion {
  const A &a;
  const B &b;
  float k; // Smoothing factor (e.g., 0.1)
  float thickness;
  static constexpr bool is_solid = true;

  SmoothUnion(const A &shapeA, const B &shapeB, float smoothness)
      : a(shapeA), b(shapeB), k(smoothness),
        thickness(std::max(shapeA.thickness, shapeB.thickness)) {}

  template <int H> Bounds get_vertical_bounds() const {
    auto b1 = a.template get_vertical_bounds<H>();
    auto b2 = b.template get_vertical_bounds<H>();
    // Expand bounds slightly to account for the blending radius (k)
    return {std::max(0, std::min(b1.y_min, b2.y_min) - 1),
            std::min(H - 1, std::max(b1.y_max, b2.y_max) + 1)};
  }

  // Conservative union of children's intervals, padded by k (in radians)
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    StaticCircularBuffer<std::pair<float, float>, 32> merged;
    float pad_px = k * W / (2 * PI_F);

    bool hasA = a.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) {
          merged.push_back({start - pad_px, end + pad_px});
        });

    bool hasB = b.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) {
          merged.push_back({start - pad_px, end + pad_px});
        });

    // If either child falls back to full-width, so must the blend
    if (!hasA || !hasB)
      return false;

    if (merged.is_empty())
      return true;

    // Sort by start
    auto *data = &merged[0];
    size_t n = merged.size();
    for (size_t i = 1; i < n; ++i) {
      auto key = data[i];
      int j = static_cast<int>(i) - 1;
      while (j >= 0 && data[j].first > key.first) {
        data[j + 1] = data[j];
        --j;
      }
      data[j + 1] = key;
    }

    // Merge overlapping intervals
    float cur_start = data[0].first;
    float cur_end = data[0].second;
    for (size_t i = 1; i < n; ++i) {
      if (data[i].first <= cur_end) {
        cur_end = std::max(cur_end, data[i].second);
      } else {
        out(cur_start, cur_end);
        cur_start = data[i].first;
        cur_end = data[i].second;
      }
    }
    out(cur_start, cur_end);
    return true;
  }

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

    // Polynomial smooth min (cubic)
    float h = std::max(k - std::abs(res.dist - resB.dist), 0.0f) / k;
    float m = h * h * h * k * (1.0f / 6.0f);

    if (res.dist < resB.dist) {
      res.dist -= m;
    } else {
      resB.dist -= m;
      res = resB;
    }
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
    return a.template get_vertical_bounds<H>();
  }

  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    StaticCircularBuffer<std::pair<float, float>, 32> intervalsA;
    StaticCircularBuffer<std::pair<float, float>, 32> intervalsB;

    bool hasA = a.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) { intervalsA.push_back({start, end}); });

    // If A falls back to full-width, we can't restrict
    if (!hasA)
      return false;

    // If A produced no intervals on this row, nothing to subtract from
    if (intervalsA.is_empty())
      return true;

    bool hasB = b.template get_horizontal_intervals<W, H>(
        y, [&](float start, float end) { intervalsB.push_back({start, end}); });

    // If B falls back (no interval info), B removes nothing we can
    // reason about — pass A's intervals through conservatively
    if (!hasB || intervalsB.is_empty()) {
      for (size_t i = 0; i < intervalsA.size(); ++i)
        out(intervalsA[i].first, intervalsA[i].second);
      return true;
    }

    // Set difference: A minus B
    // For each A interval, subtract all overlapping B intervals
    for (size_t ai = 0; ai < intervalsA.size(); ++ai) {
      float cur_start = intervalsA[ai].first;
      float cur_end = intervalsA[ai].second;

      for (size_t bi = 0; bi < intervalsB.size(); ++bi) {
        float bs = intervalsB[bi].first;
        float be = intervalsB[bi].second;

        // No overlap
        if (be <= cur_start || bs >= cur_end)
          continue;

        // B covers the left portion: emit nothing, shrink from left
        if (bs <= cur_start) {
          cur_start = be;
          if (cur_start >= cur_end)
            break;
          continue;
        }

        // B covers the right portion: shrink from right
        if (be >= cur_end) {
          cur_end = bs;
          if (cur_start >= cur_end)
            break;
          continue;
        }

        // B splits A: emit left piece, continue with right piece
        out(cur_start, bs);
        cur_start = be;
      }

      if (cur_start < cur_end)
        out(cur_start, cur_end);
    }
    return true;
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
 * @brief Angular domain repetition modifier.
 * Folds azimuthal angle to create N copies of a shape around an arbitrary axis
 * for constant cost (single distance evaluation).
 */
template <typename Shape> struct AngularRepeat {
  const Shape &shape;
  Vector axis, u, w; // axis = rotation axis, u/w = derived perpendicular plane
  int repetitions;
  float thickness;
  static constexpr bool is_solid = Shape::is_solid;

  /// Repeat around an arbitrary axis.
  AngularRepeat(const Shape &s, int reps, const Vector &ax)
      : shape(s), axis(ax), repetitions(reps), thickness(s.thickness) {
    // Derive perpendicular plane via Gram-Schmidt
    Vector ref = (fabsf(ax.y) < 0.9f) ? Vector(0, 1, 0) : Vector(1, 0, 0);
    float d = ref.x * ax.x + ref.y * ax.y + ref.z * ax.z;
    u = Vector(ref.x - d * ax.x, ref.y - d * ax.y, ref.z - d * ax.z);
    float len = u.length();
    if (len > TOLERANCE)
      u = u * (1.0f / len);
    w = Vector(ax.y * u.z - ax.z * u.y, ax.z * u.x - ax.x * u.z,
               ax.x * u.y - ax.y * u.x);
  }

  /// Repeat around the Y-axis (legacy behavior).
  AngularRepeat(const Shape &s, int reps)
      : AngularRepeat(s, reps, Vector(0, 1, 0)) {}

  template <int H> Bounds get_vertical_bounds() const {
    return shape.template get_vertical_bounds<H>();
  }

  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    // Full scan: repetitions cover the full azimuth
    return false;
  }

  DistanceResult distance(const Vector &p) const {
    DistanceResult res;
    distance<true>(p, res);
    return res;
  }

  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    // Project p into local coordinate system
    float local_u = p.x * u.x + p.y * u.y + p.z * u.z;
    float local_v = p.x * axis.x + p.y * axis.y + p.z * axis.z;
    float local_w = p.x * w.x + p.y * w.y + p.z * w.z;

    // Fold angle in the u-w plane
    float theta = atan2f(local_w, local_u);
    if (theta < 0)
      theta += 2 * PI_F;

    float sector = (2 * PI_F) / repetitions;
    float folded_theta = fmodf(theta + sector / 2.0f, sector) - sector / 2.0f;

    // Reconstruct folded local coordinates (preserving axis component)
    float r = sqrtf(local_u * local_u + local_w * local_w);
    float fu = r * cosf(folded_theta);
    float fw = r * sinf(folded_theta);

    // Project back to world space
    Vector folded_p(fu * u.x + local_v * axis.x + fw * w.x,
                    fu * u.y + local_v * axis.y + fw * w.y,
                    fu * u.z + local_v * axis.z + fw * w.z);

    shape.template distance<ComputeUVs>(folded_p, res);
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
  std::array<Vector, MAX_VERTS + 1> verts3D; // Original 3D vertices (+1 wrap)
  std::array<Vector, MAX_VERTS> edgeNormals; // Great circle plane normals
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
  std::span<Vector> verts3D;
  std::span<Vector> edgeNormals;

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
      float y = vertices[idx].y;
      if (y < min_y_val)
        min_y_val = y;
      if (y > max_y_val)
        max_y_val = y;
    }

    float min_phi_check = fast_acos(hs::clamp(max_y_val, -1.0f, 1.0f));
    float max_phi_check = fast_acos(hs::clamp(min_y_val, -1.0f, 1.0f));
    float margin_check = thickness + BOUNDS_MARGIN;

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
    if (count > FaceScratchBuffer::MAX_VERTS) {
      hs::log("Face: vertex count exceeds MAX_VERTS, truncating");
      count = FaceScratchBuffer::MAX_VERTS;
    }

    center = Vector(0, 0, 0);
    for (int i = 0; i < count; ++i)
      center = center + vertices[indices[i]];
    center.normalize();

    basisV = center;
    Vector ref = (std::abs(dot(center, X_AXIS)) > 0.9f) ? Y_AXIS : X_AXIS;
    basisU = cross(center, ref).normalized();
    basisW = cross(center, basisU).normalized();

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

    // Store 3D vertices and edge normals for angular distance computation
    for (int i = 0; i < count; ++i)
      scratch.verts3D[i] = vertices[indices[i]];
    scratch.verts3D[count] = scratch.verts3D[0];
    verts3D = std::span<Vector>(scratch.verts3D.data(), count + 1);

    for (int i = 0; i < count; ++i) {
      Vector n = cross(scratch.verts3D[i], scratch.verts3D[i + 1]);
      float len = n.magnitude();
      scratch.edgeNormals[i] =
          (len > 1e-9f) ? n * (1.0f / len) : Vector(0, 0, 0);
    }
    edgeNormals = std::span<Vector>(scratch.edgeNormals.data(), count);

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

    if (size < radius * MIN_SIZE_RADIUS_RATIO)
      size = radius * MIN_SIZE_RADIUS_RATIO;

    float vertex_phi_span = max_phi_check - min_phi_check;
    bool use_fast_bounds = (vertex_phi_span < FAST_BOUNDS_PHI_THRESHOLD);

    int planes_count = 0;

    if (use_fast_bounds) {
      compute_fast_bounds(scratch, vertices, indices, count, thickness,
                          min_phi_check, max_phi_check, h_virt, height,
                          y_min, y_max);
    } else {
      planes_count = compute_full_bounds(scratch, vertices, indices, count,
                                         center, thickness, h_virt, height,
                                         y_min, y_max);
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
        scratch.intervals[interval_count++] = {0.0f, endT};
        scratch.intervals[interval_count++] = {startT, 2 * PI_F};
      }
    } else {
      full_width = true;
    }
    intervals = std::span<std::pair<float, float>>(scratch.intervals.data(),
                                                   interval_count);

    // Pole containment: if a face encircles a pole, extend vertical bounds
    // to cover it.
    if (center.y > 0.01f) {
      float inv_c = 1.0f / center.y;
      float ppx = basisU.y * inv_c;
      float ppy = basisW.y * inv_c;
      bool pole_inside = false;
      for (int i = 0; i < count; ++i) {
        if ((poly2D[i].y > ppy) != (poly2D[i + 1].y > ppy)) {
          float ix = poly2D[i].x +
                     (ppy - poly2D[i].y) * edgeVectors[i].x * invEdgeJ[i];
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
    if (center.y < -0.01f) {
      float cos_sp = -center.y;
      float inv_c = 1.0f / cos_sp;
      float ppx = -basisU.y * inv_c;
      float ppy = -basisW.y * inv_c;
      bool pole_inside = false;
      for (int i = 0; i < count; ++i) {
        if ((poly2D[i].y > ppy) != (poly2D[i + 1].y > ppy)) {
          float ix = poly2D[i].x +
                     (ppy - poly2D[i].y) * edgeVectors[i].x * invEdgeJ[i];
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

  /// Fast-path bounds: uses vertex-based phi for small faces.
  static void compute_fast_bounds(
      FaceScratchBuffer &scratch, std::span<const Vector> vertices,
      std::span<const uint16_t> indices, int count, float thickness,
      float min_phi_check, float max_phi_check,
      int h_virt, int height, int &y_min_out, int &y_max_out) {
    for (int i = 0; i < count; ++i) {
      Vector edge = scratch.poly2D[(i + 1) % count] - scratch.poly2D[i];
      scratch.edgeVectors[i] = edge;
      float edge_len_sq = dot(edge, edge);
      scratch.edgeLengthsSq[i] = edge_len_sq;
      scratch.invEdgeLengthsSq[i] =
          (edge_len_sq > 1e-12f) ? (1.0f / edge_len_sq) : 0.0f;
      scratch.invEdgeJ[i] =
          (std::abs(edge.y) > 1e-12f) ? (1.0f / edge.y) : 0.0f;
      float theta = fast_atan2(vertices[indices[i]].z, vertices[indices[i]].x);
      if (theta < 0)
        theta += 2 * PI_F;
      scratch.thetas[i] = theta;
    }
    float margin = thickness + BOUNDS_MARGIN;
    y_min_out = std::max(
        0, static_cast<int>(floorf(
               (std::max(0.0f, min_phi_check - margin) * (h_virt - 1)) / PI_F)));
    y_max_out = std::min(
        height - 1,
        static_cast<int>(ceilf(
            (std::min(PI_F, max_phi_check + margin) * (h_virt - 1)) / PI_F)));
  }

  /// Full-path bounds: arc extrema + pole containment for large faces.
  /// Returns the number of planes computed.
  static int compute_full_bounds(
      FaceScratchBuffer &scratch, std::span<const Vector> vertices,
      std::span<const uint16_t> indices, int count,
      const Vector &center, float thickness,
      int h_virt, int height, int &y_min_out, int &y_max_out) {
    float min_phi = 100.0f;
    float max_phi = -100.0f;
    int planes_count = 0;
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
          (std::abs(edge.y) > 1e-12f) ? (1.0f / edge.y) : 0.0f;
      Vector normal = cross(v1, v2);
      float lenSq = dot(normal, normal);
      if (lenSq > 1e-12f)
        scratch.planes[planes_count++] = normal.normalized();
      float phi_val = fast_acos(hs::clamp(v1.y, -1.0f, 1.0f));
      if (phi_val < min_phi)
        min_phi = phi_val;
      if (phi_val > max_phi)
        max_phi = phi_val;
      // Arc Extrema Logic
      if (planes_count > 0) {
        const Vector &n = scratch.planes[planes_count - 1];
        float ny = n.y;
        if (std::abs(ny) < 0.99999f) {
          float nx = n.x;
          float nz = n.z;
          float tx = -nx * ny;
          float ty = 1.0f - ny * ny;
          float tz = -nz * ny;
          float tLenSq = tx * tx + ty * ty + tz * tz;
          if (tLenSq > 1e-12f) {
            float invLen = 1.0f / sqrtf(tLenSq);
            float ptx = tx * invLen;
            float pty = ty * invLen;
            float ptz = tz * invLen;
            float cx1 = (v1.y * ptz - v1.z * pty) * nx +
                        (v1.z * ptx - v1.x * ptz) * ny +
                        (v1.x * pty - v1.y * ptx) * nz;
            float cx2 = (pty * v2.z - ptz * v2.y) * nx +
                        (ptz * v2.x - ptx * v2.z) * ny +
                        (ptx * v2.y - pty * v2.x) * nz;
            if (cx1 > 0 && cx2 > 0) {
              float phiTop = fast_acos(hs::clamp(pty, -1.0f, 1.0f));
              if (phiTop < min_phi)
                min_phi = phiTop;
            }
            if (cx1 < 0 && cx2 < 0) {
              float phiBot = fast_acos(hs::clamp(-pty, -1.0f, 1.0f));
              if (phiBot > max_phi)
                max_phi = phiBot;
            }
          }
        }
      }
      float theta = fast_atan2(v1.z, v1.x);
      if (theta < 0)
        theta += 2 * PI_F;
      scratch.thetas[i] = theta;
    }
    bool npInside = (planes_count > 0);
    bool spInside = (planes_count > 0);
    for (int pi = 0; pi < planes_count; ++pi) {
      float center_side = dot(center, scratch.planes[pi]);
      if ((scratch.planes[pi].y > 0) != (center_side > 0))
        npInside = false;
      if ((-scratch.planes[pi].y > 0) != (center_side > 0))
        spInside = false;
    }
    if (npInside)
      min_phi = 0.0f;
    if (spInside)
      max_phi = PI_F;
    float margin = thickness + BOUNDS_MARGIN;
    y_min_out = std::max(
        0, static_cast<int>(floorf(
               (std::max(0.0f, min_phi - margin) * (h_virt - 1)) / PI_F)));
    y_max_out = std::min(
        height - 1,
        static_cast<int>(
            ceilf((std::min(PI_F, max_phi + margin) * (h_virt - 1)) / PI_F)));
    return planes_count;
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

      float wx = px - Vi.x;
      float wy = py - Vi.y;

      float t = (wx * edge.x + wy * edge.y) * invEdgeLengthsSq[i];
      float clampVal = hs::clamp(t, 0.0f, 1.0f);

      float bx = wx - edge.x * clampVal;
      float by = wy - edge.y * clampVal;
      float distSq = bx * bx + by * by;

      if (distSq < d)
        d = distSq;

      if ((Vi.y > py) != (Vnext.y > py)) {
        float intersect_x = Vi.x + (py - Vi.y) * edge.x * invEdgeJ[i];
        if (px < intersect_x) {
          inside = !inside;
        }
      }
    }

    float s = inside ? -1.0f : 1.0f;
    float plane_dist = s * sqrtf(d);
    // Convert gnomonic planar distance to angular for correct AA blending
    float angular_dist_raw = fast_atan2(plane_dist, 1.0f);
    float angular_dist = angular_dist_raw - thickness;

    res = DistanceResult(angular_dist, 0.0f, angular_dist_raw, 0.0f, size);
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
  static constexpr bool is_solid = true;

  Polygon(const Basis &b, float r, float th, int s, float ph, int h_virt,
          int height)
      : basis(b), thickness(th), sides(s), phase(ph) {
    apothem = thickness * cosf(PI_F / sides);
    nx = basis.v.x;
    ny = basis.v.y;
    nz = basis.v.z;
    R_val = sqrtf(nx * nx + nz * nz);
    alpha_angle = atan2f(nz, nx);

    float center_phi = acosf(std::max(-1.0f, std::min(1.0f, ny)));
    float margin = thickness + BOUNDS_MARGIN_WIDE;
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

    if (R_val < MIN_HORIZONTAL_PROJ)
      return false;

    float ang_high = thickness;
    float D_min = cosf(ang_high);

    float denom = R_val * sin_phi;
    if (std::abs(denom) < INTERVAL_DENOM_EPS)
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

    float dot_u = dot(p, basis.u);
    float dot_w = dot(p, basis.w);
    float azimuth = fast_atan2(dot_w, dot_u);
    if (azimuth < 0)
      azimuth += 2 * PI_F;
    azimuth += phase;

    float sector = 2 * PI_F / sides;
    float local = wrap(azimuth + sector / 2.0f, sector) - sector / 2.0f;

    float dist_edge = polar * cosf(local) - apothem;
    float t_val = ComputeUVs ? polar / thickness : 0.0f;

    res = DistanceResult(dist_edge, t_val, polar, 0.0f, apothem);
  }
};

/**
 * @brief Calculates signed distance to a spherical polygon (great circle
 * edges). Uses sector folding + precomputed great circle plane normal for
 * O(1) per-pixel distance, with exact angular distances for smooth AA.
 */
struct SphericalPolygon {
  const Basis &basis;
  int sides;
  float phase;
  float circumradius; // angular distance from center to vertex
  float edge_nv;      // edge normal · center
  float edge_nu;      // edge normal · u-axis
  int y_min, y_max;
  float nx, ny, nz, R_val, alpha_angle;
  static constexpr bool is_solid = true;

  SphericalPolygon(const Basis &b, float radius, int s, float ph, int h_virt,
                   int height)
      : basis(b), sides(s), phase(ph) {
    circumradius = radius * (PI_F / 2.0f);

    // Build canonical edge: between vertices at azimuth ±π/n from
    // the sector bisector (u-axis), at angular distance circumradius
    float half_step = PI_F / sides;
    float sinR = sinf(circumradius);
    float cosR = cosf(circumradius);
    float cos_hs = cosf(half_step);
    float sin_hs = sinf(half_step);

    Vector v1 = basis.v * cosR + (basis.u * cos_hs + basis.w * sin_hs) * sinR;
    Vector v2 = basis.v * cosR + (basis.u * cos_hs - basis.w * sin_hs) * sinR;

    // Normal pointing outward (away from polygon interior)
    Vector en = cross(v2, v1);
    float len = en.magnitude();
    if (len > 1e-9f)
      en = en * (1.0f / len);
    // Ensure outward: dot(center, n) should be negative
    if (dot(en, basis.v) > 0)
      en = -en;

    edge_nv = dot(en, basis.v);
    edge_nu = dot(en, basis.u);

    // Vertical/horizontal bounds (same as SDF::Polygon)
    nx = basis.v.x;
    ny = basis.v.y;
    nz = basis.v.z;
    R_val = sqrtf(nx * nx + nz * nz);
    alpha_angle = atan2f(nz, nx);

    float center_phi = acosf(std::max(-1.0f, std::min(1.0f, ny)));
    float margin = circumradius + BOUNDS_MARGIN_WIDE;
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

    if (R_val < MIN_HORIZONTAL_PROJ)
      return false;

    float ang_high = circumradius;
    float D_min = cosf(ang_high);
    float denom = R_val * sin_phi;
    if (std::abs(denom) < INTERVAL_DENOM_EPS)
      return false;

    float C_min = (D_min - ny * cos_phi) / denom;
    if (C_min > 1.0f)
      return true;
    if (C_min < -1.0f)
      return false;

    float d_alpha = acosf(C_min);
    float f_x1 = (alpha_angle - d_alpha) * W / (2 * PI_F);
    float f_x2 = (alpha_angle + d_alpha) * W / (2 * PI_F);
    out(floorf(f_x1), ceilf(f_x2));
    return true;
  }

  DistanceResult distance(const Vector &p) const {
    DistanceResult res;
    distance<true>(p, res);
    return res;
  }

  template <bool ComputeUVs = true>
  void distance(const Vector &p, DistanceResult &res) const {
    float polar = angle_between(p, basis.v);

    float dot_u = dot(p, basis.u);
    float dot_w = dot(p, basis.w);
    float azimuth = fast_atan2(dot_w, dot_u);
    if (azimuth < 0)
      azimuth += 2 * PI_F;
    azimuth += phase;

    float sector = 2 * PI_F / sides;
    float local = wrap(azimuth + sector / 2.0f, sector) - sector / 2.0f;

    // Angular distance to the nearest great circle edge via precomputed normal
    // cos(local) is even, so sector folding works automatically
    float sinP = sinf(polar);
    float cosP = cosf(polar);
    float dp = edge_nv * cosP + edge_nu * cosf(local) * sinP;
    float dist_edge = asinf(hs::clamp(dp, -1.0f, 1.0f));

    float t_val = ComputeUVs ? polar / circumradius : 0.0f;
    res = DistanceResult(dist_edge, t_val, polar, 0.0f, circumradius);
  }
};

/**
 * @brief Calculates signed distance to a star shape.
 * Returns:
 *  dist: Signed distance from edge (negative inside)
 *  t: Normalized polar distance (polar / thickness)
 *  raw_dist: Polar distance from center
 */
struct Star {
  const Basis &basis;
  int sides;
  float phase;
  static constexpr bool is_solid = true;

  float nx, ny, planeD;
  float thickness;

  // Scan
  float scanNy, scanNx, scanNz, scanR, scanAlpha;
  int yMin, yMax;

  Star(const Basis &b, float radius, int s, float ph, int h_virt, int height)
      : basis(b), sides(s), phase(ph) {
    float outerRadius = radius * (PI_F / 2.0f);
    float innerRadius = outerRadius * STAR_INNER_RATIO;
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

    scanNx = basis.v.x;
    scanNy = basis.v.y;
    scanNz = basis.v.z;
    scanR = sqrtf(scanNx * scanNx + scanNz * scanNz);
    scanAlpha = atan2f(scanNz, scanNx);

    float centerPhi = acosf(std::max(-1.0f, std::min(1.0f, basis.v.y)));
    float margin = outerRadius + BOUNDS_MARGIN_WIDE;
    yMin = std::max(
        0, static_cast<int>(floorf(
               (std::max(0.0f, centerPhi - margin) * (height - 1)) / PI_F)));
    yMax = std::min(
        height - 1,
        static_cast<int>(
            ceilf((std::min(PI_F, centerPhi + margin) * (height - 1)) / PI_F)));
  }

  template <int H> Bounds get_vertical_bounds() const { return {yMin, yMax}; }

  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    // Bounding circle
    float phi = y_to_phi<H>(static_cast<float>(y));
    float cosPhi = cosf(phi);
    float sinPhi = sinf(phi);

    if (scanR < MIN_HORIZONTAL_PROJ)
      return false;

    float pixelWidth = 2.0f * PI_F / W;
    float D_min = cosf(thickness + pixelWidth);
    float denom = scanR * sinPhi;

    if (std::abs(denom) < INTERVAL_DENOM_EPS)
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
    float azimuth = fast_atan2(dotW, dotU);
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
  float scan_nx, scan_ny, scan_nz, scan_R, scan_alpha;
  int y_min, y_max;
  static constexpr bool is_solid = true;

  Flower(const Basis &b, float radius, int s, float ph, int h_virt, int height)
      : basis(b), sides(s), phase(ph) {
    float outer = radius * (PI_F / 2.0f);
    apothem = PI_F - outer;
    thickness = outer;
    antipode = -basis.v;

    scan_nx = antipode.x;
    scan_ny = antipode.y;
    scan_nz = antipode.z;
    scan_R = sqrtf(scan_nx * scan_nx + scan_nz * scan_nz);
    scan_alpha = atan2f(scan_nz, scan_nx);

    float center_phi = acosf(std::max(-1.0f, std::min(1.0f, antipode.y)));
    float margin = thickness + BOUNDS_MARGIN_WIDE;
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

    if (scan_R < MIN_HORIZONTAL_PROJ)
      return false;

    float denom = scan_R * sin_phi;
    if (std::abs(denom) < INTERVAL_DENOM_EPS)
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

    float dot_u = dot(p, basis.u);
    float dot_w = dot(p, basis.w);
    float azimuth = fast_atan2(dot_w, dot_u);
    if (azimuth < 0)
      azimuth += 2 * PI_F;
    azimuth += phase;

    float sector = 2 * PI_F / sides;
    float local = wrap(azimuth + sector / 2.0f, sector) - sector / 2.0f;

    float dist_edge = polar * cosf(local) - apothem;
    float t_val = ComputeUVs ? scan_dist / thickness : 0.0f;

    res = DistanceResult(-dist_edge, t_val, scan_dist, 0.0f, thickness);
  }
};

/**
 * @brief Signed distance to a spherical harmonic "blob".
 *
 * Vertical and horizontal culling are intentionally disabled:
 * spherical harmonic lobes for arbitrary (l,m) can occupy any
 * region of the sphere, so no static bounding is possible
 * without per-(l,m) analysis.  Full-sphere scan is correct.
 */
struct HarmonicBlob {
  int l;
  int m;
  float amplitude;
  Quaternion inv_q;
  HarmonicWaveFn harmonic_fn;
  static constexpr bool is_solid = true;

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
    Vector v = rotate(p, inv_q);
    float phi = acosf(std::max(-1.0f, std::min(1.0f, v.y)));
    float theta = fast_atan2(v.z, v.x);
    float harmonic_val = harmonic_fn(l, m, theta, phi);

    float lobe_radius = 1.0f + std::abs(harmonic_val) * amplitude;
    float d = 1.0f - lobe_radius;

    float t_val = 0.0f;
    if constexpr (ComputeUVs) {
      t_val = tanhf(std::abs(harmonic_val) * amplitude);
    }

    res = DistanceResult(d, t_val, harmonic_val, 0.0f, 1.0f);
  }

  /// Full-sphere scan (see struct doc).
  template <int W, int H, typename OutputIt>
  bool get_horizontal_intervals(int y, OutputIt out) const {
    return false;
  }
};

struct Line {
  Vector a, b;
  float thickness;

  Vector n;
  float len;
  float phi_min, phi_max; // Precomputed vertical bounds (radians)
  static constexpr bool is_solid = false;

  Line(const Vector &start, const Vector &end, float th)
      : a(start), b(end), thickness(th) {
    len = angle_between(a, b);
    if (len < 1e-6f) {
      n = Vector(0, 0, 0);
    } else {
      n = cross(a, b).normalized();
    }

    // Compute vertical bounds from endpoint phi values + thickness
    float phi_a = acosf(std::max(-1.0f, std::min(1.0f, a.y)));
    float phi_b = acosf(std::max(-1.0f, std::min(1.0f, b.y)));
    float margin = thickness + BOUNDS_MARGIN;
    phi_min = std::max(0.0f, std::min(phi_a, phi_b) - margin);
    phi_max = std::min(PI_F, std::max(phi_a, phi_b) + margin);
  }

  template <int H> Bounds get_vertical_bounds() const {
    constexpr int H_VIRT = H + hs::H_OFFSET;
    int y_min = std::max(
        0, static_cast<int>(floorf((phi_min * (H_VIRT - 1)) / PI_F)));
    int y_max = std::min(
        H - 1, static_cast<int>(ceilf((phi_max * (H_VIRT - 1)) / PI_F)));
    return {y_min, y_max};
  }

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
    constexpr int H_VIRT = H + hs::H_OFFSET;
    float phi = y_to_phi<H_VIRT>(static_cast<float>(y));
    float cos_phi = cosf(phi);
    float sin_phi = sinf(phi);

    // Bounding cap centered on the midpoint of the line segment
    Vector mid = (a + b).normalized();
    float mid_ny = mid.y;
    float mid_r = sqrtf(mid.x * mid.x + mid.z * mid.z);
    float mid_alpha = atan2f(mid.z, mid.x);

    if (mid_r < MIN_HORIZONTAL_PROJ)
      return false;

    float cap_radius = len * 0.5f + thickness;
    float D_min = cosf(cap_radius);
    float denom = mid_r * sin_phi;
    if (std::abs(denom) < INTERVAL_DENOM_EPS)
      return false;

    float C_min = (D_min - mid_ny * cos_phi) / denom;
    if (C_min > 1.0f)
      return true; // Row is outside the bounding cap
    if (C_min < -1.0f)
      return false; // Cap covers full width

    float d_alpha = acosf(C_min);
    float f_x1 = (mid_alpha - d_alpha) * W / (2 * PI_F);
    float f_x2 = (mid_alpha + d_alpha) * W / (2 * PI_F);
    out(floorf(f_x1), ceilf(f_x2));
    return true;
  }
};

// ============================================================================
// 3D Volumetric SDF Shapes (for Scan::Volume raymarching)
// ============================================================================

/**
 * @brief 3D Torus signed distance field.
 *
 * The torus ring lies in the XZ plane with symmetry axis along Y.
 * Major radius R = ring centerline distance from origin.
 * Minor radius r = tube cross-section radius.
 *
 * Unlike the 2D spherical SDF shapes above (which use DistanceResult),
 * 3D volumetric shapes return plain float distances and operate in
 * Cartesian ray-space.
 */
struct Torus {
  float R; ///< Major radius (ring centerline)
  float r; ///< Minor radius (tube cross-section)

  /// Signed distance from point p to the torus surface.
  /// Negative = inside, positive = outside.
  float distance(const Vector &p) const {
    float q = sqrtf(p.x * p.x + p.z * p.z) - R;
    return sqrtf(q * q + p.y * p.y) - r;
  }

  /// Surface normal at a point near the torus surface.
  Vector normal(const Vector &p) const {
    float xz_len = sqrtf(p.x * p.x + p.z * p.z);
    float scale = (xz_len > TOLERANCE) ? R / xz_len : 0.0f;
    return (p - Vector(p.x * scale, 0.0f, p.z * scale)).normalized();
  }

  /// Populate a Fragment's registers for shading.
  /// v0 = ring angle (0-1, for palette lookup)
  /// v1,v2,v3 = surface normal (x,y,z)
  /// size = signed distance (for AA edge alpha)
  void populate(const Vector &p, Fragment &frag) const {
    Vector n = normal(p);
    frag.v0 = (fast_atan2(p.z, p.x) + PI_F) / (2.0f * PI_F);
    frag.v1 = n.x;
    frag.v2 = n.y;
    frag.v3 = n.z;
  }
};

/**
 * @brief Domain warp functions for composing with WarpedSDF.
 */
namespace Warp {

/**
 * @brief Oscillates the Y coordinate sinusoidally around the azimuthal
 * angle θ = atan2(z, x).
 *
 * Produces twisted/undulating geometry when composed with a torus.
 * Provides an analytical Lipschitz bound for safe sphere tracing.
 */
struct Twist {
  int twist;       ///< Number of oscillations around the ring
  float amplitude; ///< Vertical displacement magnitude
  float R;         ///< Major radius (needed for Lipschitz bound)

  /// Precomputed context: s = sqrtf(x² + z²), shared across apply/lipschitz.
  using Ctx = float;

  Ctx make_ctx(const Vector &p) const { return sqrtf(p.x * p.x + p.z * p.z); }

  /// Warp the domain: displace Y by amplitude * sin(twist * θ).
  Vector apply(const Vector &p, Ctx /*s*/) const {
    float theta = fast_atan2(p.z, p.x);
    return Vector(
        p.x, p.y - amplitude * sinf(static_cast<float>(twist) * theta), p.z);
  }

  /// Analytical Lipschitz constant at point p.
  float lipschitz(const Vector & /*p*/, Ctx s) const {
    if (twist == 0)
      return 1.0f;
    return 1.0f + static_cast<float>(twist) * amplitude /
                      (2.0f * std::max(s, R * 0.5f));
  }

  /// Maximum possible inflation of the bounding volume.
  float bounding_inflation() const { return amplitude; }

  /// Analytical normal correction via chain rule (called once per hit).
  Vector correct_normal(const Vector &p, const Vector &base_n, Ctx s) const {
    if (twist == 0 || amplitude < TOLERANCE)
      return base_n;
    float inv_s = (s > TOLERANCE) ? 1.0f / s : 0.0f;
    float theta = fast_atan2(p.z, p.x);
    float n_theta = static_cast<float>(twist) * theta;
    float dh_dtheta = -amplitude * static_cast<float>(twist) * cosf(n_theta);
    float inv_s2 = inv_s * inv_s;

    float dh_dx = dh_dtheta * (-p.z) * inv_s2;
    float dh_dz = dh_dtheta * p.x * inv_s2;

    return Vector(base_n.x - base_n.y * dh_dx, base_n.y,
                  base_n.z - base_n.y * dh_dz)
        .normalized();
  }
};

} // namespace Warp

/**
 * @brief Composable domain-warped volume SDF.
 *
 * Applies a Warp to the input domain of any 3D volume SDF, with:
 * - Bounding fast-path (skips warp trig when far from surface)
 * - Analytical Lipschitz correction for safe sphere tracing
 * - Normal correction via the warp's chain rule
 *
 * The Warp must provide:
 *   Ctx make_ctx(const Vector &p) const;
 *   Vector apply(const Vector &p, const Ctx &ctx) const;
 *   float lipschitz(const Vector &p, const Ctx &ctx) const;
 *   float bounding_inflation() const;
 * Optionally:
 *   Vector correct_normal(const Vector &p, const Vector &n, const Ctx&) const;
 */
template <typename SDF, typename Warp> struct WarpedVolume {
  SDF base;
  Warp warp;

  /// Cheap lower-bound: base distance minus warp's maximum displacement.
  float bounding_distance(const Vector &p) const {
    return base.distance(p) - warp.bounding_inflation();
  }

  /// Raw SDF distance (no Lipschitz correction). Use for surface projection.
  float raw_distance(const Vector &p) const {
    auto ctx = warp.make_ctx(p);
    return base.distance(warp.apply(p, ctx));
  }

  /// March-safe distance with bounding fast-path and Lipschitz correction.
  float distance(const Vector &p) const {
    float bd = bounding_distance(p);
    if (bd > warp.bounding_inflation())
      return bd;

    auto ctx = warp.make_ctx(p);
    float d = base.distance(warp.apply(p, ctx));

    if (d > 0.0f) {
      float lip = warp.lipschitz(p, ctx);
      if (lip > 1.0f)
        d /= lip;
    }
    return d;
  }

  /// Surface normal with warp chain-rule correction.
  Vector normal(const Vector &p) const {
    auto ctx = warp.make_ctx(p);
    Vector base_n = base.normal(warp.apply(p, ctx));
    if constexpr (requires { warp.correct_normal(p, base_n, ctx); }) {
      return warp.correct_normal(p, base_n, ctx);
    }
    return base_n;
  }

  /// Populate a Fragment's registers for shading.
  void populate(const Vector &p, Fragment &frag) const {
    Vector n = normal(p);
    frag.v0 = (fast_atan2(p.z, p.x) + PI_F) / (2.0f * PI_F);
    frag.v1 = n.x;
    frag.v2 = n.y;
    frag.v3 = n.z;
  }
};

} // namespace SDF
