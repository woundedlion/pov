/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#include <functional>
#include <utility>
#include <type_traits>
#include <algorithm>
#include <cmath>
#include <array>
#include <concepts>
#include "geometry.h"
#include "color.h"
#include "constants.h"
#include "canvas.h"
#include "animation.h"

namespace Plot {

/**
 * @brief Single point primitive.
 */
struct Point {
  /**
   * @brief Draws a single point.
   */
  static void draw(PipelineRef pipeline, Canvas &canvas, const Fragment &f,
                   FragmentShaderFn fragment_shader) {
    Fragment f_copy = f;
    f_copy.color = Color4(0, 0, 0, 0);
    f_copy.blend = 0;
    fragment_shader(f_copy.pos, f_copy);
    pipeline.plot(canvas, f_copy.pos, f_copy.color.color, f_copy.age,
                  f_copy.color.alpha, f_copy.blend);
  }
};

/**
 * @brief Core rasterization logic for 3D lines and curves.
 * Adapts step size based on screen-space density to avoid aliasing.
 */
// --- Strategy Helpers ---

// Planar Strategy
template <typename ProcessSegmentFn>
static void
rasterize_planar_strategy(const Fragment &curr, const Fragment &next,
                          const Basis &planar_basis, bool isLastSegment,
                          ProcessSegmentFn &&process_segment) {
  // match JS: PlanarStrategy(planarBasis)
  const Vector &u = planar_basis.u;
  const Vector &center = planar_basis.v;
  const Vector &w = planar_basis.w;

  auto project = [&](const Vector &p) -> std::pair<float, float> {
    float R = angle_between(p, center);
    if (R < 1e-5f)
      return {0.0f, 0.0f};
    float x = dot(p, u);
    float y = dot(p, w);
    float theta = atan2f(y, x);
    return {R * cosf(theta), R * sinf(theta)};
  };

  auto proj1 = project(curr.pos);
  auto proj2 = project(next.pos);
  float dx = proj2.first - proj1.first;
  float dy = proj2.second - proj1.second;
  float dist = sqrtf(dx * dx + dy * dy);

  auto map_planar = [=](float t) {
    float Px = proj1.first + dx * t;
    float Py = proj1.second + dy * t;

    // Unproject
    float R = sqrtf(Px * Px + Py * Py);
    if (R < 1e-5f)
      return center;

    float theta = atan2f(Py, Px);
    Vector axis = (u * cosf(theta)) + (w * sinf(theta));
    return (center * cosf(R)) + (axis * sinf(R));
  };

  process_segment(map_planar, curr, next, dist, isLastSegment);
}

// Geodesic Strategy
template <typename ProcessSegmentFn>
static void rasterize_geodesic_strategy(const Fragment &curr,
                                        const Fragment &next,
                                        bool isLastSegment,
                                        ProcessSegmentFn &&process_segment) {
  // GEODESIC STRATEGY (Standard)
  Vector v1 = curr.pos;
  Vector v2 = next.pos;
  float total_dist = angle_between(v1, v2);

  if (total_dist < 0.001f) {
    auto map_degenerate = [=](float t) { return v1; };
    process_segment(map_degenerate, curr, next, total_dist, isLastSegment);
  } else {
    Vector axis;
    if (std::abs(PI_F - total_dist) < TOLERANCE) {
      axis = (std::abs(dot(v1, X_AXIS)) > 0.999f) ? cross(v1, Y_AXIS)
                                                  : cross(v1, X_AXIS);
    } else {
      axis = cross(v1, v2).normalize();
    }

    Vector v_perp = cross(axis, v1);

    auto map_geodesic = [=](float t) {
      float ang = total_dist * t;
      float s = sinf(ang);
      float c = cosf(ang);
      return (v1 * c) + (v_perp * s);
    };
    process_segment(map_geodesic, curr, next, total_dist, isLastSegment);
  }
}

template <int W, int H>
static void rasterize(PipelineRef pipeline, Canvas &canvas,
                      const Fragments &points, FragmentShaderFn fragment_shader,
                      bool close_loop = false, float age = 0.0f,
                      const Basis *planar_basis = nullptr) {
  size_t len = points.size();
  if (len < 2)
    return;

  size_t count = close_loop ? len : len - 1;
  ScopedScratch _sc(scratch_arena_a);
  ArenaVector<float> _steps_cache;
  size_t max_cache = std::min(std::max(len * 4, (size_t)256), (size_t)2048);
  _steps_cache.bind(scratch_arena_a, max_cache);

  auto process_segment = [&](auto &&map, const Fragment &curr,
                             const Fragment &next, float total_dist,
                             bool isLastSegment) {
    // Handle Degenerate Segment
    if (total_dist < 1e-5f) {
      bool shouldOmit = (close_loop) ? true : !isLastSegment;
      if (!shouldOmit) {
        Fragment f_copy = curr;
        // Set temp values for shader
        f_copy.pos = curr.pos;
        f_copy.color = Color4(0, 0, 0, 0);
        f_copy.blend = 0;

        fragment_shader(curr.pos, f_copy);
        pipeline.plot(canvas, curr.pos, f_copy.color.color, f_copy.age,
                      f_copy.color.alpha, f_copy.blend);
      }
      return;
    }

    // 1. SIMULATION PHASE
    _steps_cache.clear();
    float sim_dist = 0.0f;
    const float base_step = (2.0f * PI_F) / W;

    // Calculate initial point for density
    Vector p_temp = map(0.0f);

    while (sim_dist < total_dist) {
      // Adaptive step size based on distortion (y-component)
      float scale_factor =
          std::max(0.05f, sqrtf(std::max(0.0f, 1.0f - p_temp.y * p_temp.y)));
      float step = base_step * scale_factor;

      _steps_cache.push_back(step);
      sim_dist += step;

      if (sim_dist < total_dist) {
        p_temp = map(sim_dist / total_dist);
      }
    }

    float scale = (sim_dist > 0.0f) ? (total_dist / sim_dist) : 0.0f;
    bool omitLast = (close_loop) ? true : !isLastSegment;

    if (omitLast && _steps_cache.is_empty())
      return;

    // 2. DRAWING PHASE
    {
      Vector start_pos = map(0.0f);
      Fragment f = Fragment::lerp(curr, next, 0.0f);
      f.pos = start_pos;
      f.color = Color4(0, 0, 0, 0);
      f.blend = 0;

      fragment_shader(start_pos, f);
      pipeline.plot(canvas, start_pos, f.color.color, f.age, f.color.alpha,
                    f.blend);
    }

    size_t loop_limit =
        omitLast ? _steps_cache.size() - 1 : _steps_cache.size();
    float current_dist = 0.0f;

    for (size_t j = 0; j < loop_limit; j++) {
      float step = _steps_cache[j] * scale;
      current_dist += step;

      float t = (total_dist > 0.0f) ? (current_dist / total_dist) : 1.0f;

      Vector p = map(t);
      Fragment f = Fragment::lerp(curr, next, t);
      f.pos = p;
      f.color = Color4(0, 0, 0, 0);
      f.blend = 0;

      fragment_shader(p, f);
      pipeline.plot(canvas, p, f.color.color, f.age, f.color.alpha, f.blend);
    }
  };

  for (size_t i = 0; i < count; i++) {
    const Fragment &curr = points[i];
    const Fragment &next = points[(i + 1) % len];
    bool isLastSegment = (i == count - 1);

    // --- Interpolation Strategy Selection ---
    if (planar_basis) {
      rasterize_planar_strategy(curr, next, *planar_basis, isLastSegment,
                                process_segment);
    } else {
      rasterize_geodesic_strategy(curr, next, isLastSegment, process_segment);
    }
  }
}

/**
 * @brief Draws a geodesic line between two points.
 * Registers:
 *  v0: line progress (0..1)
 *  v1: Arc Length (radians)
 */
struct Line {
  /**
   * @brief Samples a geodesic line between two points.
   * @param v1 Start point.
   * @param v2 End point.
   * @return Fragments Two fragments representing the line endpoints.
   */
  static void sample(Fragments &points, const Fragment &f1, const Fragment &f2,
                     int density = 1) {
    if (density < 1)
      density = 1;

    float angle = angle_between(f1.pos, f2.pos);
    if (std::abs(angle) < 0.0001f) {
      points.push_back(f1);
      points.push_back(f1); // Draw at least a dot
      return;
    }

    Vector axis = cross(f1.pos, f2.pos).normalize();

    for (int i = 0; i <= density; ++i) {
      float t = static_cast<float>(i) / density;

      Fragment f = Fragment::lerp(f1, f2, t);
      if (i == 0)
        f.pos = f1.pos;
      else if (i == density)
        f.pos = f2.pos;
      else {
        Quaternion q = make_rotation(axis, angle * t);
        f.pos = rotate(f1.pos, q);
      }

      f.v0 = t;
      f.v1 = angle * t;
      f.v2 = 0.0f; // Standard: Single lines have v2=0
      points.push_back(f);
    }
  }

  /**
   * @brief Draws a geodesic line.
   * @tparam W Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param f1 Start fragment.
   * @param f2 End fragment.
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Fragment &f1,
                   const Fragment &f2, FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader) {
    ScopedScratch _frag(scratch_arena_a);
    Fragments points;
    points.bind(scratch_arena_a, 4);
    sample(points, f1, f2);

    if (vertex_shader) {
      for (auto &p : points) {
        vertex_shader(p);
      }
    }
    rasterize<W, H>(pipeline, canvas, points, fragment_shader, false, 0.0f,
                    nullptr);
  }

  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Fragment &f1,
                   const Fragment &f2, FragmentShaderFn fragment_shader) {
    draw<W, H>(pipeline, canvas, f1, f2, fragment_shader, {});
  }
};

/**
 * @brief Draws vertices without connection.
 */
struct Vertices {
  /**
   * @brief Draws a list of vertices as points.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param points List of points.
   * @param fragment_shader Shader function.
   */
  static void draw(PipelineRef pipeline, Canvas &canvas, const auto &points,
                   FragmentShaderFn fragment_shader) {
    for (const Fragment &p : points) {
      Fragment f = p;
      f.color = Color4(0, 0, 0, 0);
      f.blend = 0;

      fragment_shader(f.pos, f);
      pipeline.plot(canvas, f.pos, f.color.color, f.age, f.color.alpha,
                    f.blend);
    }
  }
};

/**
 * @brief Multiline primitive (Polyline).
 * Registers:
 *  v0: Path Progress (0.0 -> 1.0)
 *  v1: Cumulative Arc Length (radians)
 *  v2: Vertex Index
 */
struct Multiline {
  /**
   * @brief Samples a multiline path from a container of vertices.
   * @param vertices Iterable container of Vector.
   * @param closed If true, connects the last point to the first.
   * @return Fragments List of fragments with arc-length parameterization.
   */
  static void sample(Fragments &points, const auto &vertices,
                     bool closed = false) {
    auto it = std::begin(vertices);
    auto end = std::end(vertices);

    if (it == end)
      return;

    // 1. Calculate total length to normalize v0
    float total_len = 0.0f;
    Fragment first = *it;
    Fragment prev = first;
    auto len_it = it;
    ++len_it;

    for (; len_it != end; ++len_it) {
      const Fragment &curr = *len_it;
      total_len += angle_between(prev.pos, curr.pos);
      prev = curr;
    }
    if (closed) {
      total_len += angle_between(prev.pos, first.pos);
    }

    if (total_len < 1e-6f)
      total_len = 1.0f; // Avoid divide by zero

    // 2. Generate fragments
    float current_len = 0.0f;
    it = std::begin(vertices); // Reset iterator
    prev = *it;

    // Add first point
    Fragment f = prev;
    f.v0 = 0.0f;
    f.v1 = 0.0f;
    f.v2 = 0.0f; // Index 0
    points.push_back(f);

    ++it;
    int idx = 1;
    for (; it != end; ++it) {
      const Fragment &curr = *it;
      float dist = angle_between(prev.pos, curr.pos);
      current_len += dist;

      f = curr;
      f.v0 = current_len / total_len;
      f.v1 = current_len;
      f.v2 = static_cast<float>(idx++);
      points.push_back(f);
      prev = curr;
    }

    if (closed) {
      float dist = angle_between(prev.pos, first.pos);
      current_len += dist;
      f = first;
      f.v0 = 1.0f; // Explicitly 1.0
      f.v1 = current_len;
      f.v2 = static_cast<float>(idx);
      points.push_back(f);
    }
  }

  /**
   * @brief Draws a multiline path.
   * @tparam W Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param vertices Iterable container of Vector.
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader.
   * @param closed If true, connects the last point to the first.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const auto &vertices,
                   FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader, bool closed = false) {
    ScopedScratch _frag(scratch_arena_a);
    Fragments points;
    points.bind(scratch_arena_a, vertices.size() + 2);
    sample(points, vertices, closed);

    if (vertex_shader) {
      for (auto &p : points) {
        vertex_shader(p);
      }
    }

    // We manually closed it if requested, so pass false to rasterize
    rasterize<W, H>(pipeline, canvas, points, fragment_shader, false, 0.0f,
                    nullptr);
  }

  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const auto &vertices,
                   FragmentShaderFn fragment_shader, bool closed = false) {
    draw<W, H>(pipeline, canvas, vertices, fragment_shader, {}, closed);
  }
};

/**
 * @brief Ring primitives.
 * Registers:
 *  v0: Angular progress (0.0 -> 1.0)
 *  v1: Arc Length (radians)
 *  v2: Index
 */
struct Ring {
  /**
   * @brief Calculates a 3D point on a ring given basis and progress.
   */
  static Vector calcPoint(float a, float radius, const Vector &u,
                          const Vector &v, const Vector &w) {
    auto d = sqrtf((1 - radius) * (1 - radius));
    return Vector(d * v.x + radius * u.x * cosf(a) + radius * w.x * sinf(a),
                  d * v.y + radius * u.y * cosf(a) + radius * w.y * sinf(a),
                  d * v.z + radius * u.z * cosf(a) + radius * w.z * sinf(a))
        .normalize();
  }

  static void sample(Fragments &points, const Basis &basis, float radius,
                     int num_samples, float phase = 0) {
    auto res = get_antipode(basis, radius);
    const Basis &work_basis = res.first;
    float work_radius = res.second;

    const Vector &v = work_basis.v;
    const Vector &u = work_basis.u;
    const Vector &w = work_basis.w;

    const float theta_eq = work_radius * (PI_F / 2.0f);
    const float r_val = sinf(theta_eq);
    const float d_val = cosf(theta_eq);
    const float arc_scale = sinf(work_radius);

    const float step = 2.0f * PI_F / num_samples;

    for (int i = 0; i < num_samples; i++) {
      float theta = i * step;
      float t = theta + phase;
      Vector u_temp = (u * cosf(t)) + (w * sinf(t));

      Fragment f;
      f.pos = ((v * d_val) + (u_temp * r_val)).normalize();
      f.v0 = static_cast<float>(i) / num_samples;
      f.v1 = theta * arc_scale;
      f.v2 = static_cast<float>(i);
      f.age = 0;

      points.push_back(f);
    }

    // Manual Close (Overlap)
    if (num_samples > 0) {
      Fragment f;
      float theta = 2.0f * PI_F;
      float t = theta + phase;
      Vector u_temp = (u * cosf(t)) + (w * sinf(t));
      f.pos = ((v * d_val) + (u_temp * r_val)).normalize();

      f.v0 = 1.0f;
      f.v1 = theta * arc_scale;
      f.v2 = static_cast<float>(num_samples);
      f.age = 0;
      points.push_back(f);
    }
  }

  /**
   * @brief Draws a ring.
   * @tparam W Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param basis Orientation basis.
   * @param radius Ring radius (radians).
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader.
   * @param phase Rotation phase.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader, float phase = 0) {
    ScopedScratch _frag(scratch_arena_a);
    Fragments points;
    points.bind(scratch_arena_a, W + 2);
    // Use W samples for smooth circles (fixes pinching at poles)
    sample(points, basis, radius, W, phase);

    if (vertex_shader) {
      for (auto &p : points) {
        vertex_shader(p);
      }
    }
    rasterize<W, H>(pipeline, canvas, points, fragment_shader, true, 0.0f,
                    nullptr);
  }

  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, FragmentShaderFn fragment_shader,
                   float phase = 0) {
    draw<W, H>(pipeline, canvas, basis, radius, fragment_shader, {}, phase);
  }
};

/**
 * @brief Planar Polygon.
 * Registers:
 *  v0: Perimeter progress (0.0 -> 1.0)
 *  v1: Arc Length (radians)
 *  v2: Vertex index
 */
struct PlanarPolygon {
  /**
   * @brief Samples a planar polygon.
   */
  static void sample(Fragments &points, const Basis &basis, float radius,
                     int num_sides, float phase = 0) {
    size_t start_idx = points.size();
    Ring::sample(points, basis, radius, num_sides, phase + PI_F / num_sides);
    float cumul = 0.0f;
    for (size_t i = start_idx; i < points.size(); i++) {
      points[i].v1 = cumul;
      if (i < points.size() - 1) {
        cumul += angle_between(points[i].pos, points[i + 1].pos);
      }
    }
  }

  /**
   * @brief Draws a planar polygon.
   * @tparam W Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param basis Orientation basis.
   * @param radius Polygon radius.
   * @param num_sides Number of sides.
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader.
   * @param phase Rotation phase.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, int num_sides,
                   FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader, float phase = 0) {
    ScopedScratch _frag(scratch_arena_a);
    Fragments points;
    points.bind(scratch_arena_a, num_sides + 2);
    sample(points, basis, radius, num_sides, phase);

    if (vertex_shader) {
      for (auto &p : points) {
        vertex_shader(p);
      }
    }
    rasterize<W, H>(pipeline, canvas, points, fragment_shader, true, 0.0f,
                    &basis);
  }

  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, int num_sides,
                   FragmentShaderFn fragment_shader, float phase = 0) {
    draw<W, H>(pipeline, canvas, basis, radius, num_sides, fragment_shader, {},
               phase);
  }
};

/**
 * @brief Spherical Polygon (Geodesic edges).
 * Registers:
 *  v0: Perimeter progress (0.0 -> 1.0)
 *  v1: Arc Length (radians)
 *  v2: Vertex index
 */
struct SphericalPolygon {
  /**
   * @brief Samples a spherical polygon.
   */
  static void sample(Fragments &points, const Basis &basis, float radius,
                     int num_sides, float phase = 0) {
    size_t start_idx = points.size();
    Ring::sample(points, basis, radius, num_sides, phase + PI_F / num_sides);

    // Re-calculate v1 to be true geodesic chord length
    float cumulative_length = 0.0f;
    for (size_t i = start_idx; i < points.size(); ++i) {
      points[i].v2 = static_cast<float>(i - start_idx);

      if (i > start_idx) {
        cumulative_length += angle_between(points[i - 1].pos, points[i].pos);
      }
      points[i].v1 = cumulative_length;
    }
  }

  /**
   * @brief Draws a spherical polygon.
   * @tparam W Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param basis Orientation basis.
   * @param radius Polygon radius.
   * @param num_sides Number of sides.
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader.
   * @param phase Rotation phase.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, int num_sides,
                   FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader, float phase = 0) {
    ScopedScratch _frag(scratch_arena_a);
    Fragments points;
    points.bind(scratch_arena_a, num_sides + 2);
    sample(points, basis, radius, num_sides, phase);

    if (vertex_shader) {
      for (auto &p : points) {
        vertex_shader(p);
      }
    }
    rasterize<W, H>(pipeline, canvas, points, fragment_shader, true, 0.0f,
                    nullptr);
  }

  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, int num_sides,
                   FragmentShaderFn fragment_shader, float phase = 0) {
    draw<W, H>(pipeline, canvas, basis, radius, num_sides, fragment_shader, {},
               phase);
  }
};

/**
 * @brief Distorted Ring.
 * Registers:
 *  v0: Angular progress (0.0 -> 1.0)
 *  v1: Arc Length (radians)
 *  v2: Index
 */
struct DistortedRing {
  /**
   * @brief Calculates a single point on a distorted ring.
   */
  static Vector fn_point(ScalarFn shift_fn, const Basis &basis, float radius,
                         float angle) {
    auto res = get_antipode(basis, radius);
    const Basis &work_basis = res.first;
    float work_radius = res.second;

    const Vector &v = work_basis.v;
    const Vector &u = work_basis.u;
    const Vector &w = work_basis.w;

    auto vi = Ring::calcPoint(angle, work_radius, u, v, w);
    auto vp = Ring::calcPoint(angle, 1, u, v, w);
    Vector axis = cross(v, vp).normalize();
    return rotate(vi, make_rotation(axis, shift_fn(angle * PI_F / 2)));
  }

  template <int W>
  /**
   * @brief Samples a distorted ring.
   */
  static void sample(Fragments &points, const Basis &basis, float radius,
                     ScalarFn shift_fn, float phase = 0) {
    auto res = get_antipode(basis, radius);
    const Basis &work_basis = res.first;
    float work_radius = res.second;

    const Vector &v = work_basis.v;
    const Vector &u = work_basis.u;
    const Vector &w = work_basis.w;

    const float theta_eq = work_radius * (PI_F / 2.0f);
    const float r_val = sinf(theta_eq);
    const float d_val = cosf(theta_eq);

    const int num_samples = W;
    const float step = 2.0f * PI_F / num_samples;

    float cumulative_len = 0.0f;
    Vector last_pos;

    for (int i = 0; i < num_samples; i++) {
      float theta = i * step;
      float t = theta + phase;
      Vector u_temp = (u * cosf(t)) + (w * sinf(t));

      float shift = shift_fn(theta / (2.0f * PI_F));
      float cos_shift = cosf(shift);
      float sin_shift = sinf(shift);

      float v_scale = d_val * cos_shift - r_val * sin_shift;
      float u_scale = r_val * cos_shift + d_val * sin_shift;

      Fragment f;
      f.pos = ((v * v_scale) + (u_temp * u_scale)).normalize();

      if (i > 0) {
        cumulative_len += angle_between(last_pos, f.pos);
      }
      last_pos = f.pos;

      f.v0 = static_cast<float>(i) / num_samples;
      f.v1 = cumulative_len;
      f.v2 = static_cast<float>(i);
      f.age = 0;
      points.push_back(f);
    }

    // Manual Close (Overlap)
    if (num_samples > 0) {
      Fragment f;
      const Fragment &first = points[0];
      f.pos = first.pos;
      cumulative_len += angle_between(last_pos, f.pos);
      f.v0 = 1.0f;
      f.v1 = cumulative_len;
      f.v2 = static_cast<float>(num_samples);
      f.age = 0;
      points.push_back(f);
    }
  }

  /**
   * @brief Draws a distorted ring.
   * @tparam W Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param basis Orientation basis.
   * @param radius Base radius.
   * @param shift_fn Distortion function.
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader.
   * @param phase Rotation phase.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, ScalarFn shift_fn,
                   FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader, float phase = 0) {
    ScopedScratch _frag(scratch_arena_a);
    Fragments points;
    points.bind(scratch_arena_a, W + 2);
    sample<W>(points, basis, radius, shift_fn, phase);

    if (vertex_shader) {
      for (auto &p : points) {
        vertex_shader(p);
      }
    }
    rasterize<W, H>(pipeline, canvas, points, fragment_shader, true, 0.0f,
                    nullptr);
  }

  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, ScalarFn shift_fn,
                   FragmentShaderFn fragment_shader, float phase = 0) {
    draw<W, H>(pipeline, canvas, basis, radius, shift_fn, fragment_shader, {},
               phase);
  }
};

/**
 * @brief Fibonacci Spiral primitive.
 */
struct Spiral {
  /**
   * @brief Samples a Fibonacci spiral.
   */
  static void sample(Fragments &fragments, int n, float eps) {
    float cumulative_len = 0.0f;
    Vector last_pos;

    for (int i = 0; i < n; i++) {
      Vector pos = fib_spiral(n, eps, i);

      if (i > 0) {
        cumulative_len += angle_between(last_pos, pos);
      }
      last_pos = pos;

      Fragment f;
      f.pos = pos;
      f.v0 = (n > 1) ? static_cast<float>(i) / (n - 1) : 0.0f;
      f.v1 = cumulative_len;
      f.age = 0;
      fragments.push_back(f);
    }
  }

  /**
   * @brief Draws a Fibonacci spiral.
   * @tparam W Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param n Number of points.
   * @param eps Parameter.
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, int n, float eps,
                   FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader) {
    ScopedScratch _frag(scratch_arena_a);
    Fragments frags;
    frags.bind(scratch_arena_a, n);
    sample(frags, n, eps);

    if (vertex_shader) {
      for (auto &f : frags) {
        vertex_shader(f);
      }
    }
    rasterize<W, H>(pipeline, canvas, frags, fragment_shader, false, 0.0f,
                    nullptr);
  }

  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, int n, float eps,
                   FragmentShaderFn fragment_shader) {
    draw<W, H>(pipeline, canvas, n, eps, fragment_shader, {});
  }
};

/**
 * @brief Star shape.
 * Registers:
 *  v0: Perimeter progress (0.0 -> 1.0)
 *  v1: Arc Length (radians)
 *  v2: Vertex index
 */
struct Star {
  /**
   * @brief Samples a star shape.
   */
  static void sample(Fragments &points, const Basis &basis, float radius,
                     int num_sides, float phase = 0) {
    auto res = get_antipode(basis, radius);
    const Basis &work_basis = res.first;
    float work_radius = res.second;

    const Vector &v = work_basis.v;
    const Vector &u = work_basis.u;
    const Vector &w = work_basis.w;

    float outer_radius = work_radius * (PI_F / 2.0f);
    float inner_radius = outer_radius * 0.382f;
    float angle_step = PI_F / num_sides;

    float cumulative_len = 0.0f;
    size_t start_idx = points.size();

    for (int i = 0; i < num_sides * 2; i++) {
      float theta = phase + i * angle_step;
      float r = (i % 2 == 0) ? outer_radius : inner_radius;

      float sin_r = sinf(r);
      float cos_r = cosf(r);
      float cos_t = cosf(theta);
      float sin_t = sinf(theta);

      Vector p = (v * cos_r) + (u * (cos_t * sin_r)) + (w * (sin_t * sin_r));
      p.normalize();

      Fragment f;
      f.pos = p;
      if (i > 0) {
        cumulative_len += angle_between(points.back().pos, p);
      }

      f.v0 = static_cast<float>(i) / (num_sides * 2);
      f.v1 = cumulative_len;
      f.v2 = static_cast<float>(i);
      f.age = 0;
      points.push_back(f);
    }

    // Manual Close
    if (points.size() > start_idx) {
      const Fragment &first = points[start_idx];
      Fragment last = first;
      cumulative_len += angle_between(points.back().pos, first.pos);
      last.v0 = 1.0f;
      last.v1 = cumulative_len;
      last.v2 = static_cast<float>(num_sides * 2);
      points.push_back(last);
    }
  }

  /**
   * @brief Draws a star.
   * @tparam W Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param basis Orientation basis.
   * @param radius Outer radius.
   * @param num_sides Number of points.
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader.
   * @param phase Rotation phase.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, int num_sides,
                   FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader, float phase = 0) {
    ScopedScratch _frag(scratch_arena_a);
    Fragments points;
    points.bind(scratch_arena_a, num_sides * 2 + 2);
    sample(points, basis, radius, num_sides, phase);

    if (vertex_shader) {
      for (auto &p : points) {
        vertex_shader(p);
      }
    }

    Vector center = basis.v;
    if (radius > 1.0f)
      center = -center;

    Vector v = center;
    Vector ref = (std::abs(dot(v, X_AXIS)) > 0.9f) ? Y_AXIS : X_AXIS;
    Vector u = cross(v, ref).normalize();
    Vector w = cross(v, u).normalize();
    Basis planar_basis = {u, v, w};
    rasterize<W, H>(pipeline, canvas, points, fragment_shader, true, 0.0f,
                    &planar_basis);
  }

  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, int num_sides,
                   FragmentShaderFn fragment_shader, float phase = 0) {
    draw<W, H>(pipeline, canvas, basis, radius, num_sides, fragment_shader, {},
               phase);
  }
};

/**
 * @brief Flower shape.
 * Registers:
 *  v0: Perimeter progress (0.0 -> 1.0)
 *  v1: Arc Length (radians)
 *  v2: Vertex index
 */
/**
 * @brief Flower primitive.
 */
struct Flower {
  /**
   * @brief Samples a flower shape.
   */
  static void sample(Fragments &points, const Basis &basis, float radius,
                     int num_sides, float phase = 0) {
    auto res = get_antipode(basis, radius);
    const Basis &work_basis = res.first;
    float work_radius = res.second;

    const Vector &v = work_basis.v;
    const Vector &u = work_basis.u;
    const Vector &w = work_basis.w;

    float desired_outer_radius = work_radius * (PI_F / 2.0f);
    float apothem = PI_F - desired_outer_radius;
    float safe_apothem = std::min(apothem, PI_F - 1e-4f);
    float angle_step = PI_F / num_sides;

    float cumulative_length = 0.0f;
    size_t start_idx = points.size();

    for (int i = 0; i < num_sides * 2; i++) {
      float theta = phase + i * angle_step;
      float R = safe_apothem;

      float sin_r = sinf(R);
      float cos_r = cosf(R);
      float cos_t = cosf(theta);
      float sin_t = sinf(theta);

      // Unproject Polar -> Sphere
      Vector p = (v * cos_r) + (u * (cos_t * sin_r)) + (w * (sin_t * sin_r));
      p.normalize(); // Ensure unit vector

      Fragment f;
      f.pos = p;

      if (i > 0) {
        cumulative_length += angle_between(points.back().pos, p);
      }

      f.v0 = static_cast<float>(i) / (num_sides * 2);
      f.v1 = cumulative_length;
      f.v2 = static_cast<float>(i);
      f.age = 0;
      points.push_back(f);
    }

    // Close Loop
    if (points.size() > start_idx) {
      const Fragment &first = points[start_idx];
      Fragment last = first; // copy
      cumulative_length += angle_between(points.back().pos, first.pos);
      last.v0 = 1.0f;
      last.v1 = cumulative_length;
      last.v2 = static_cast<float>(num_sides * 2);
      points.push_back(last);
    }
  }

  /**
   * @brief Draws a flower.
   * @tparam W Rasterization resolution.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param basis Orientation basis.
   * @param radius Outer radius.
   * @param num_sides Number of petals.
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader.
   * @param phase Rotation phase.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, int num_sides,
                   FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader, float phase = 0) {
    ScopedScratch _frag(scratch_arena_a);
    Fragments points;
    points.bind(scratch_arena_a, num_sides * 2 + 2);
    sample(points, basis, radius, num_sides, phase);

    if (vertex_shader) {
      for (auto &p : points) {
        vertex_shader(p);
      }
    }

    auto res = get_antipode(basis, radius);
    const Basis &work_basis = res.first;
    Vector center = work_basis.v;

    // Construct a planar basis aligned with the center
    Vector ref = (std::abs(dot(center, X_AXIS)) > 0.9f) ? Y_AXIS : X_AXIS;
    Vector u_p = cross(center, ref).normalize();
    Vector w_p = cross(center, u_p).normalize();
    Basis planar_basis = {u_p, center, w_p};
    rasterize<W, H>(pipeline, canvas, points, fragment_shader, true, 0.0f,
                    &planar_basis);
  }

  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const Basis &basis,
                   float radius, int num_sides,
                   FragmentShaderFn fragment_shader, float phase = 0) {
    draw<W, H>(pipeline, canvas, basis, radius, num_sides, fragment_shader, {},
               phase);
  }
};

/**
 * @brief Mesh drawing.
 * Registers:
 *  v0: Edge Progress t (0.0 -> 1.0) per edge
 *  v1: Cumulative Arc Length (radians) per edge
 *  v2: Vertex index
 */
struct Mesh {
  /**
   * @brief Draws a mesh (wireframe).
   * @tparam W Rasterization resolution.
   * @tparam MeshT Mesh type.
   * @param pipeline Render pipeline.
   * @param canvas Target canvas.
   * @param mesh The mesh to draw.
   * @param fragment_shader Shader function.
   * @param vertex_shader Optional vertex shader.
   */
  /**
   * @brief Samples edges of a mesh.
   * @tparam MeshT Mesh type.
   * @param mesh The mesh to sample.
   * @param density Sampling density per edge.
   * @return List of sampled edges (each edge is a list of Fragments).
   */
  template <int W, int H, typename MeshT>
  static void draw(PipelineRef pipeline, Canvas &canvas, const MeshT &mesh,
                   FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader) {
    int edge_index = 0;

    StaticCircularBuffer<uint32_t, 4096> visited;

    auto process_edge = [&](int u, int v) {
      int small = std::min(u, v);
      int large = std::max(u, v);
      uint32_t key =
          (static_cast<uint32_t>(small) << 16) | static_cast<uint32_t>(large);

      for (auto curr : visited) {
        if (curr == key)
          return;
      }
      visited.push_back(key);

      Fragment fu;
      fu.pos = mesh.vertices[u];
      Fragment fv;
      fv.pos = mesh.vertices[v];

      ScopedScratch _edge(scratch_arena_a);
      Fragments points;
      points.bind(scratch_arena_a, 16);
      Line::sample(points, fu, fv, 10);

      if (vertex_shader) {
        for (auto &p : points) {
          p.v2 = static_cast<float>(edge_index); // Edge Index
          vertex_shader(p);
        }
      } else {
        for (auto &p : points) {
          p.v2 = static_cast<float>(edge_index);
        }
      }
      rasterize<W, H>(pipeline, canvas, points, fragment_shader, false, 0.0f,
                      nullptr);

      edge_index++;
    };

    size_t offset = 0;
    const uint8_t *fc = mesh.get_face_counts_data();
    size_t num_f = mesh.get_face_counts_size();
    const uint16_t *fi = mesh.get_faces_data();

    for (size_t i = 0; i < num_f; ++i) {
      int count = fc[i];
      for (int k = 0; k < count; ++k) {
        int u = fi[offset + k];
        int v = fi[offset + (k + 1) % count];
        process_edge(u, v);
      }
      offset += count;
    }
  }

  template <int W, int H, typename MeshT>
  static void draw(PipelineRef pipeline, Canvas &canvas, const MeshT &mesh,
                   FragmentShaderFn fragment_shader) {
    draw<W, H>(pipeline, canvas, mesh, fragment_shader, {});
  }
};

/**
 * @brief Particle System trails.
 * Registers:
 *  v0: Trail Progress (0.0=Head -> 1.0=Tail)
 *  v1: Trail Arc Length
 *  v2: Particle ID
 *  v3: Normalized TTL
 */
struct ParticleSystem {
  /**
   * @brief Samples particle trails.
   */
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const auto &system,
                   FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader) {
    int count = system.active_count;

    for (int i = 0; i < count; ++i) {
      const auto &p = system.pool[i];
      ScopedScratch _trail(scratch_arena_a);
      Fragments trail;
      trail.bind(scratch_arena_a, 64);
      float cumulative_len = 0.0f;
      Vector last_pos;
      bool first = true;

      tween(p.history, [&](const Vector &v, float t) {
        if (!first) {
          cumulative_len += angle_between(last_pos, v);
        }
        last_pos = v;
        first = false;

        Fragment f;
        f.pos = v;
        f.v0 = t;
        f.v1 = cumulative_len;
        f.v2 = static_cast<float>(i);
        f.v3 = static_cast<float>(p.life) / static_cast<float>(system.max_life);
        f.age = 0;
        f.color = Color4(0, 0, 0, 0);
        trail.push_back(f);
      });

      if (!trail.is_empty()) {
        if (vertex_shader) {
          for (auto &f : trail) {
            vertex_shader(f);
          }
        }
        rasterize<W, H>(pipeline, canvas, trail, fragment_shader, false, 0.0f,
                        nullptr);
      }
    }
  }

  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas, const auto &system,
                   FragmentShaderFn fragment_shader) {
    draw<W, H>(pipeline, canvas, system, fragment_shader, {});
  }
};

/**
 * @brief Single cubic Bézier curve on the sphere.
 * Registers:
 *  v0: curve progress (0→1)
 *  v1: cumulative arc length (radians)
 *  v2: 0 (single segment)
 */
struct Bezier {
  static void sample(Fragments &points,
                     const Vector &p0, const Vector &p1,
                     const Vector &p2, const Vector &p3,
                     int num_samples,
                     SplineMode mode = SplineMode::Geodesic) {
    float cumulative_len = 0.0f;
    Vector last_pos = p0;

    for (int i = 0; i <= num_samples; ++i) {
      float t = static_cast<float>(i) / num_samples;
      Vector pos = Spline::cubic(p0, p1, p2, p3, t, mode);

      if (i > 0)
        cumulative_len += angle_between(last_pos, pos);
      last_pos = pos;

      Fragment f;
      f.pos = pos;
      f.v0 = t;
      f.v1 = cumulative_len;
      f.v2 = 0.0f;
      f.age = 0;
      points.push_back(f);
    }
  }

  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas,
                   const Vector &p0, const Vector &p1,
                   const Vector &p2, const Vector &p3,
                   FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader = {},
                   int num_samples = 32,
                   SplineMode mode = SplineMode::Geodesic) {
    ScopedScratch _frag(scratch_arena_a);
    Fragments points;
    points.bind(scratch_arena_a, num_samples + 1);
    sample(points, p0, p1, p2, p3, num_samples, mode);

    if (vertex_shader) {
      for (auto &f : points) vertex_shader(f);
    }
    rasterize<W, H>(pipeline, canvas, points, fragment_shader,
                     false, 0.0f, nullptr);
  }

  /// Convenience overload taking Fragment references (uses .pos).
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas,
                   const Fragment &f0, const Fragment &f1,
                   const Fragment &f2, const Fragment &f3,
                   FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader = {},
                   int num_samples = 32,
                   SplineMode mode = SplineMode::Geodesic) {
    draw<W, H>(pipeline, canvas, f0.pos, f1.pos, f2.pos, f3.pos,
               fragment_shader, vertex_shader, num_samples, mode);
  }
};

/**
 * @brief Catmull-Rom spline chain on the sphere.
 * Registers:
 *  v0: global progress (0→1 across full chain)
 *  v1: cumulative arc length (radians)
 *  v2: segment index
 */
struct SplineChain {
  template <int W, int H>
  static void draw(PipelineRef pipeline, Canvas &canvas,
                   const Fragments &control_points,
                   float tension,
                   FragmentShaderFn fragment_shader,
                   VertexShaderRef vertex_shader = {},
                   bool closed = false,
                   int samples_per_segment = 16,
                   SplineMode mode = SplineMode::Geodesic) {
    size_t n = control_points.size();
    if (n < 2) return;

    ScopedScratch _frag(scratch_arena_a);
    Fragments points;
    points.bind(scratch_arena_a, n * samples_per_segment + 1);

    float cumulative_len = 0.0f;
    Vector last_pos;
    bool first_point = true;

    size_t seg_count = closed ? n : n - 1;
    for (size_t i = 0; i < seg_count; ++i) {
      size_t i0 = (i == 0 && !closed) ? 0 : (i - 1 + n) % n;
      size_t i1 = i;
      size_t i2 = (i + 1) % n;
      size_t i3 = (i + 2 >= n && !closed) ? n - 1 : (i + 2) % n;

      Vector cp1, cp2;
      Spline::catmull_rom_tangents(
          control_points[i0].pos, control_points[i1].pos,
          control_points[i2].pos, control_points[i3].pos,
          tension, cp1, cp2);

      int start_j = (i == 0) ? 0 : 1;
      for (int j = start_j; j <= samples_per_segment; ++j) {
        float local_t = static_cast<float>(j) / samples_per_segment;
        Vector pos = Spline::cubic(
            control_points[i1].pos, cp1, cp2,
            control_points[i2].pos, local_t, mode);

        if (!first_point)
          cumulative_len += angle_between(last_pos, pos);
        last_pos = pos;
        first_point = false;

        Fragment f = Fragment::lerp(control_points[i1],
                                    control_points[i2], local_t);
        f.pos = pos;
        f.v0 = (static_cast<float>(i) + local_t) / seg_count;
        f.v1 = cumulative_len;
        f.v2 = static_cast<float>(i);
        points.push_back(f);
      }
    }

    if (vertex_shader) {
      for (auto &f : points) vertex_shader(f);
    }
    rasterize<W, H>(pipeline, canvas, points, fragment_shader,
                     false, 0.0f, nullptr);
  }
};

} // namespace Plot
