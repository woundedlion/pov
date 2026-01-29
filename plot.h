/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#include <vector>
#include <algorithm>
#include <cmath>
#include <array>
#include <concepts>
#include "geometry.h"
#include "color.h"
#include "led.h"

/**
 * @brief The Plot struct contains vector-based (thin) drawing primitives.
 */
template <int W>
struct Plot {

  /**
   * @brief Draws a single point.
   */
  struct Point {
    static void draw(auto& pipeline, Canvas& canvas, const Vector& v, const auto& color) {
      // Support both Pixel and Color4
      if constexpr (std::is_same_v<std::decay_t<decltype(color)>, Pixel>) {
        pipeline.plot(canvas, v, color, 0, 1.0f);
      }
      else {
        pipeline.plot(canvas, v, color.color, 0, color.alpha);
      }
    }
  };

  /**
   * @brief Draws a geodesic line between two points.
   */
  struct Line {
    static void draw(auto& pipeline, Canvas& canvas, const Vector& v1, const Vector& v2, ColorFn auto color_fn,
      float start = 0.0f, float end = 1.0f, bool long_way = false, bool omit_last = false)
    {
      Vector u(v1);
      Vector v(v2);
      float a = angle_between(u, v);
      Vector w;

      if (std::abs(a) < TOLERANCE) {
        if (!omit_last) {
          auto c = color_fn(u, 0.0f);
          pipeline.plot(canvas, u, c.color, 0, c.alpha);
        }
        return;
      }
      else if (std::abs(PI_F - a) < TOLERANCE) {
        if (std::abs(dot(v, X_AXIS)) > 0.9999f) {
          w = cross(u, Y_AXIS).normalize();
        }
        else {
          w = cross(u, X_AXIS).normalize();
        }
      }
      else {
        w = cross(u, v).normalize();
      }

      if (long_way) {
        a = 2 * PI_F - a;
        w = -w;
      }

      if (std::abs(start) > TOLERANCE) {
        Quaternion q = make_rotation(w, start * a);
        u = rotate(u, q).normalize();
      }
      a *= std::abs(end - start);

      // Simulation Phase
      Vector sim_u = u;
      float sim_angle = 0;
      std::array<float, W / 2> steps;
      size_t step_count = 0;
      const float base_step = 2 * PI_F / W;
      size_t max_steps = W / 2;

      while (sim_angle < a && step_count < max_steps) {
        float scale_factor = std::max(0.05f, sqrtf(std::max(0.0f, 1.0f - sim_u.j * sim_u.j)));
        float step = base_step * scale_factor;
        steps[step_count++] = step;
        sim_angle += step;
        Quaternion q = make_rotation(w, step);
        sim_u = rotate(sim_u, q).normalize();
      }

      // Drawing phase
      float scale = (sim_angle > TOLERANCE) ? (a / sim_angle) : 1.0f;
      if (step_count == 0) {
        if (!omit_last) {
          auto c = color_fn(u, 0.0f);
          pipeline.plot(canvas, u, c.color, 0, c.alpha);
        }
        return;
      }
      float current_angle = 0;

      auto c0 = color_fn(u, 0.0f);
      pipeline.plot(canvas, u, c0.color, 0, c0.alpha);

      size_t loop_limit = omit_last ? step_count - 1 : step_count;
      for (size_t i = 0; i < loop_limit; i++) {
        float step = steps[i] * scale;
        Quaternion q = make_rotation(w, step);
        u = rotate(u, q).normalize();
        current_angle += step;
        float t = (a > 0) ? (current_angle / a) : 1;

        auto c = color_fn(u, t);
        pipeline.plot(canvas, u, c.color, 0, c.alpha);
      }
    }
  };

  /**
   * @brief Helper to rasterize a list of points into segments.
   */
  static void rasterize(auto& pipeline, Canvas& canvas, const Points& points, ColorFn auto color_fn, bool close_loop = false) {
    size_t len = points.size();
    if (len == 0) return;

    size_t count = close_loop ? len : len - 1;
    for (size_t i = 0; i < count; i++) {
      const Vector& p1 = points[i];
      const Vector& p2 = points[(i + 1) % len];

      auto segment_color_fn = [&](const Vector& p, float sub_t) {
        float global_t = (i + sub_t) / count;
        return color_fn(p, global_t);
        };

      Line::draw(pipeline, canvas, p1, p2, segment_color_fn, 0.0f, 1.0f, false, true);
    }
  }

  /**
   * @brief Draws vertices without connection.
   */
  struct Vertices {
    static void draw(auto& pipeline, Canvas& canvas, const Points& points, ColorFn auto color_fn) {
      for (const Vector& v : points) {
        auto c = color_fn(v, 0);
        pipeline.plot(canvas, v, c.color, 0, c.alpha);
      }
    }
  };

  struct Polyhedron {
    static void sample(Points& out, const VertexList& vertices, const AdjacencyList& edges) {
      for (size_t i = 0; i < edges.size(); ++i) {
        for (auto j : edges[i]) {
          out.push_back(Vector(vertices[i]).normalize());
          out.push_back(Vector(vertices[j]).normalize());
        }
      }
    }

    static void draw(auto& pipeline, Canvas& canvas, const VertexList& vertices, const AdjacencyList& edges, ColorFn auto color_fn) {
      for (size_t i = 0; i < edges.size(); ++i) {
        Vector a(vertices[i]);
        for (auto j : edges[i]) {
          if (i < j) {
            Vector b(vertices[j]);
            Line::draw(pipeline, canvas, a, b, color_fn, 0, 1, false, true);
          }
        }
      }
    }
  };

  /**
   * @brief Ring primitives.
   */
  struct Ring {
    static Vector calcPoint(float a, float radius, const Vector& u, const Vector& v, const Vector& w) {
      auto d = sqrtf((1 - radius) * (1 - radius));
      return Vector(
        d * v.i + radius * u.i * cosf(a) + radius * w.i * sinf(a),
        d * v.j + radius * u.j * cosf(a) + radius * w.j * sinf(a),
        d * v.k + radius * u.k * cosf(a) + radius * w.k * sinf(a)
      ).normalize();
    }

    static void sample(Points& points, const Basis& basis, float radius, int num_samples, float phase = 0) {
      const Vector& v = basis.v;
      const Vector& u = basis.u;
      const Vector& w = basis.w;

      Vector v_dir = (radius > 1.0f) ? -v : v;
      float r_eff = (radius > 1.0f) ? (2.0f - radius) : radius;

      const float theta_eq = r_eff * (PI_F / 2.0f);
      const float r_val = sinf(theta_eq);
      const float d_val = cosf(theta_eq);

      const float step = 2.0f * PI_F / num_samples;
      for (int i = 0; i < num_samples; i++) {
        float theta = i * step;
        float t = theta + phase;
        Vector u_temp = (u * cosf(t)) + (w * sinf(t));
        points.push_back(((v_dir * d_val) + (u_temp * r_val)).normalize());
      }
    }

    static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, ColorFn auto color_fn, float phase = 0) {
      Points points;
      sample(points, basis, radius, W / 4, phase);
      rasterize(pipeline, canvas, points, color_fn, true);
    }
  };

  struct PlanarLine {
    static void draw(auto& pipeline, Canvas& canvas, const Vector& v1, const Vector& v2, const Vector& center, ColorFn auto color_fn) {
      // Basis for projection
      Vector ref_axis = (std::abs(dot(center, X_AXIS)) > 0.9999f) ? Y_AXIS : X_AXIS;
      Vector v = center; // The 'pole' for the projection
      Vector ref = (std::abs(dot(v, X_AXIS)) > 0.9f) ? Y_AXIS : X_AXIS;
      Vector u = cross(v, ref).normalize();
      Vector w = cross(v, u).normalize();

      auto project = [&](const Vector& p) -> std::pair<float, float> {
        float R = angle_between(p, v);        
        if (R < 0.0001f) return { 0.0f, 0.0f };
        float x = dot(p, u);
        float y = dot(p, w);
        float theta = atan2f(y, x);
        return { R * cosf(theta), R * sinf(theta) };
      };

      auto p1 = project(v1);
      auto p2 = project(v2);

      float dx = p1.first - p2.first;
      float dy = p1.second - p2.second;
      float dist = sqrtf(dx * dx + dy * dy);
      
      int num_steps = std::max(2, static_cast<int>(ceilf(dist * W / (2.0f * PI_F))));

      for (int i = 0; i < num_steps; i++) {
        float t = static_cast<float>(i) / (num_steps - 1);
        float Px = p1.first + (p2.first - p1.first) * t;
        float Py = p1.second + (p2.second - p1.second) * t;

        float R = sqrtf(Px * Px + Py * Py);
        float theta = atan2f(Py, Px);

        Vector point = v;
        if (R > 0.0001f) {
           float sinR = sinf(R);
           float cosR = cosf(R);
           float cosT = cosf(theta);
           float sinT = sinf(theta);

           // dir = u*cosT + w*sinT
           Vector dir = (u * cosT) + (w * sinT);
           // p = v*cosR + dir*sinR
           point = (v * cosR) + (dir * sinR);
        }
        point = point.normalize(); // Ensure unit vector

        auto c = color_fn(point, t); 
        pipeline.plot(canvas, point, c.color, 0, c.alpha);
      }
    }
  };

  struct Polygon {
    static void sample(Points& points, const Basis& basis, float radius, int num_sides, float phase = 0) {
       // Offset by half-sector to align with Scan.Polygon edges
       const float offset = PI_F / num_sides;
       Ring::sample(points, basis, radius, num_sides, phase + offset);
    }
    static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, int num_sides, ColorFn auto color_fn, float phase = 0) {
      Points points;
      sample(points, basis, radius, num_sides, phase);
      
      Vector center = basis.v;
      if (radius > 1.0f) center = -center;
      
      size_t len = points.size();
      if (len < 2) return;

      for (size_t i = 0; i < len; i++) {
        const Vector& p1 = points[i];
        const Vector& p2 = points[(i + 1) % len];
        PlanarLine::draw(pipeline, canvas, p1, p2, center, [&](const Vector& p, float t) {
            return color_fn(p, t);
        });
      }
    }
  };

  struct DistortedRing {
    static Vector fn_point(ScalarFn auto f, const Basis& basis, float radius, float angle) {
      Vector v = basis.v; Vector u = basis.u; Vector w = basis.w;
      if (radius > 1) { v = -v; u = -u; radius = 2 - radius; }

      auto vi = Ring::calcPoint(angle, radius, u, v, w);
      auto vp = Ring::calcPoint(angle, 1, u, v, w);
      Vector axis = cross(v, vp).normalize();
      return rotate(vi, make_rotation(axis, f(angle * PI_F / 2)));
    }

    static void sample(Points& points, const Basis& basis, float radius, ScalarFn auto shift_fn, float phase = 0) {
      Vector v = basis.v; Vector u = basis.u; Vector w = basis.w;

      float v_sign = 1.0f;
      if (radius > 1.0f) { v_sign = -1.0f; radius = 2.0f - radius; }

      const float theta_eq = radius * (PI_F / 2.0f);
      const float r_val = sinf(theta_eq);
      const float d_val = cosf(theta_eq);

      const int num_samples = W;
      const float step = 2.0f * PI_F / num_samples;
      for (int i = 0; i < num_samples; i++) {
        float theta = i * step;
        float t = theta + phase;
        Vector u_temp = (u * cosf(t)) + (w * sinf(t));

        float shift = shift_fn(theta / (2.0f * PI_F));
        float cos_shift = cosf(shift);
        float sin_shift = sinf(shift);

        float v_scale = (v_sign * d_val) * cos_shift - r_val * sin_shift;
        float u_scale = r_val * cos_shift + (v_sign * d_val) * sin_shift;

        points.push_back(((v * v_scale) + (u_temp * u_scale)).normalize());
      }
    }

    static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, ScalarFn auto shift_fn, ColorFn auto color_fn, float phase = 0) {
      Points points;
      sample(points, basis, radius, shift_fn, phase);
      rasterize(pipeline, canvas, points, color_fn, true);
    }
  };

  struct Spiral {
    static void draw(auto& pipeline, Canvas& canvas, int n, float eps, ColorFn auto color_fn) {
      for (int i = 0; i < n; ++i) {
        Vector v = fib_spiral(n, eps, i);
        auto c = color_fn(v, 0);
        pipeline.plot(canvas, v, c.color, 0, c.alpha);
      }
    }
  };

  struct Star {
    static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, int num_sides, ColorFn auto color_fn, float phase = 0) {
      Vector v = basis.v;
      Vector u = basis.u;
      Vector w = basis.w;

      // Handle backside
      if (radius > 1.0f) {
        v = -v;
        u = -u;
        radius = 2.0f - radius;
      }

      float outer_radius = radius * (PI_F / 2.0f);
      float inner_radius = outer_radius * 0.382f; // Pentagram ratio

      Points points;
      float angle_step = PI_F / num_sides;

      for (int i = 0; i < num_sides * 2; i++) {
        float theta = phase + i * angle_step;
        float r = (i % 2 == 0) ? outer_radius : inner_radius;

        float sin_r = sinf(r);
        float cos_r = cosf(r);
        float cos_t = cosf(theta);
        float sin_t = sinf(theta);

        // p = v*cosR + u*cosT*sinR + w*sinT*sinR
        Vector p = (v * cos_r) + (u * (cos_t * sin_r)) + (w * (sin_t * sin_r));
        points.push_back(p.normalize());
      }

      // Projection center is v (Front)
      Vector center = v;
      size_t len = points.size();
      
      for (size_t i = 0; i < len; i++) {
        const Vector& p1 = points[i];
        const Vector& p2 = points[(i + 1) % len];
        PlanarLine::draw(pipeline, canvas, p1, p2, center, [&](const Vector& p, float t) {
          return color_fn(p, t);
        });
      }
    }
  };

  struct Flower {
    static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, int num_sides, ColorFn auto color_fn, float phase = 0) {
      Vector v = basis.v; 
      Vector u = basis.u; 
      Vector w = basis.w;

      if (radius > 1.0f) {
        v = -v;
        u = -u;
        radius = 2.0f - radius;
      }

      float desired_outer_radius = radius * (PI_F / 2.0f);
      float apothem = PI_F - desired_outer_radius;
      float angle_step = PI_F / num_sides;

      Points points;
      int num_segments = std::max(2, W / num_sides);

      for (int i = 0; i < num_sides; i++) {
        float sector_center = phase + i * 2 * angle_step;

        for (int j = 0; j < num_segments; j++) {
          float t = static_cast<float>(j) / num_segments;
          float local_phi = -angle_step + t * (2 * angle_step);
          
          float R = apothem / cosf(local_phi);
          if (R > PI_F) R = PI_F;

          float theta = sector_center + local_phi;

          float sin_r = sinf(R);
          float cos_r = cosf(R);
          float cos_t = cosf(theta);
          float sin_t = sinf(theta);

          Vector p = (v * cos_r) + (u * (cos_t * sin_r)) + (w * (sin_t * sin_r));
          points.push_back(p.normalize());
        }
      }

      if (!points.is_empty()) points.push_back(points[0]);
      rasterize(pipeline, canvas, points, color_fn, false);
    }
  };

};
