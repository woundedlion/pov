/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#include <vector>
#include <algorithm>
#include <cmath>
#include <array>
#include <span>
#include "geometry.h"
#include "color.h"
#include "led.h" // For H, W, H_VIRT
#include <concepts>
#include "filter.h"

// Forward declaration to resolve circular dependency
float quintic_kernel(float t);

static void tween(Tweenable auto& history, auto draw_fn) {
  size_t s = history.length();
  size_t start = (s > 1) ? 1 : 0;
  for (size_t i = start; i < s; ++i) {
    draw_fn(history.get(i), static_cast<float>((s - 1 - i)) / s);
  }
}

static void deep_tween(auto& trail, TweenFn auto drawFn) {
  float dt = (trail.length() > 0) ? (1.0f / static_cast<float>(trail.length())) : 0.0f;
  tween(trail, [&](const auto& frame, float t) {
    tween(frame, [&](const auto& q, float subT) {
      float globalT = t + subT * dt;
      drawFn(q, globalT);
      });
    });
}

/**
 * @brief Creates a basis { u, v, w } from an orientation and normal.
 */
struct Basis {
  Vector u, v, w;
};

static Basis make_basis(const Quaternion& orientation, const Vector& normal) {
  Vector ref_axis = (std::abs(dot(normal, X_AXIS)) > 0.9999f) ? Y_AXIS : X_AXIS;
  Vector v = rotate(normal, orientation).normalize();
  Vector ref = rotate(ref_axis, orientation).normalize();
  Vector u = cross(v, ref).normalize();
  Vector w = cross(v, u).normalize();
  return { u, v, w };
}

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
   * @brief Represents a customizable path.
   * Retains internal buffer for state, but draws to pipeline.
   */
  class Path {
  public:
    Path() {}

    static void draw(auto& pipeline, Canvas& canvas, const Path& path, ColorFn auto color) {
      size_t samples = path.points.size();
      for (size_t i = 0; i < samples; ++i) {
        auto v = path.get_point(static_cast<float>(i) / samples);
        auto c = color(v, i / (samples - 1.0f));
        pipeline.plot(canvas, v, c.color, 0, c.alpha);
      }
    }

    Path& append_segment(PlotFn auto plot, float domain, float samples, ScalarFn auto easing) {
      if (!points.is_empty()) points.pop_back();
      for (float t = 0; t <= samples; t++) {
        points.push_back(plot(easing(t / samples) * domain));
      }
      return *this;
    }

    Vector get_point(float t) const {
      if (points.is_empty()) return Vector(0, 0, 0);
      float raw_index = t * (points.size() - 1);
      size_t i = static_cast<size_t>(raw_index);
      float f = raw_index - i;
      if (i >= points.size() - 1) return points.back();
      const Vector& p1 = points[i];
      const Vector& p2 = points[i + 1];
      return p1 * (1.0f - f) + p2 * f;
    }

    size_t num_points() const { return points.size(); }
    void collapse() { if (points.size() > 1) points = { points.back() }; }
    const Points& get_points() const {
      static Points temp;
      temp.clear();
      for (size_t i = 0; i < points.size(); ++i) temp.push_back(points[i]);
      return temp;
    }

  private:
    StaticCircularBuffer<Vector, 4096> points;
  };

  /**
   * @brief Represents a path defined by a single procedural function.
   * Matches interface of Path for use in Motion animations.
   */
  template <PlotFn F>
  struct ProceduralPath {
    F f;

    ProceduralPath(F path_fn) : f(path_fn) {}

    Vector get_point(float t) const {
      return f(t);
    }
  };

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

  struct SDF {
    struct Bounds {
      int y_min, y_max;
    };
    struct Interval {
      int start, end;
    };
    struct DistanceResult {
      float dist; // Signed distance (negative inside)
      float t;    // Normalized parameter (0-1) or angle
      float raw_dist; // Unsigned or supplementary distance
    };

    struct Ring {
      const Basis& basis;
      float radius;
      float thickness;
      float phase;
      
      Vector normal, u, w;
      float nx, ny, nz;
      float target_angle, center_phi;
      float cos_max, cos_min, cos_target, inv_sin_target;

       // Optional optimization for "Scan Full Row" check
      float r_val;
      float alpha_angle; // 'alpha' in JS

      Ring(const Basis& b, float r, float th, float ph = 0)
        : basis(b), radius(r), thickness(th), phase(ph)
      {
         normal = basis.v;
         u = basis.u;
         w = basis.w;
         nx = normal.i; ny = normal.j; nz = normal.k;

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
         r_val = sqrtf(nx*nx + nz*nz);
         alpha_angle = atan2f(nx, nz);
      }

      Bounds get_vertical_bounds() const {
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

        float f_phi_min = std::max(0.0f, phi_min - thickness);
        float f_phi_max = std::min(PI_F, phi_max + thickness);
        
        int y_min = std::max(0, static_cast<int>(floorf((f_phi_min * (H_VIRT - 1)) / PI_F)));
        int y_max = std::min(H - 1, static_cast<int>(ceilf((f_phi_max * (H_VIRT - 1)) / PI_F)));
        return { y_min, y_max };
      }

      // Output iterator or callback for intervals? Using a fixed size callback/vector might be complex.
      // Let's return a small static vector or use a callback. 
      // JS returns array of objects. C++: we'll use a callback or `std::vector` (alloc). 
      // Since this is for rasterization optimization, let's use a callback concept or simple optional check.
      // For this port, let's use a `std::vector` or similar. To avoid allocs, maybe pass a buffer.
      // JS logic: addWindow pushes to `intervals`.
      template<typename OutputIt>
      bool get_horizontal_intervals(int y, OutputIt out) const {
        float phi = y_to_phi(static_cast<float>(y));
        float cos_phi = cosf(phi);
        float sin_phi = sinf(phi);

        if (r_val < 0.01f) return false; 
        
        float denom = r_val * sin_phi;
        if (std::abs(denom) < 0.000001f) return false; 

        float C_min = (cos_min - ny * cos_phi) / denom;
        float C_max = (cos_max - ny * cos_phi) / denom;
        float min_cos = std::max(-1.0f, C_min);
        float max_cos = std::min(1.0f, C_max);

        if (min_cos > max_cos) return true; // Empty row

        float angle_min = acosf(max_cos);
        float angle_max = acosf(min_cos);

        // Windows: [alpha - angle_max, alpha - angle_min] and [alpha + angle_min, alpha + angle_max]
        out(alpha_angle - angle_max, alpha_angle - angle_min);
        out(alpha_angle + angle_min, alpha_angle + angle_max);
        
        return true;
      }
      
      // Let's implement distance.
      DistanceResult distance(const Vector& p) const {
         float d = dot(p, normal);
         if (d < cos_min || d > cos_max) return { 100.0f, 0.0f, 100.0f };
         
         float dist = 0;
         if (inv_sin_target != 0) {
            dist = std::abs(d - cos_target) * inv_sin_target;
         } else {
            float polar = acosf(std::max(-1.0f, std::min(1.0f, d)));
            dist = std::abs(polar - target_angle);
         }

         float dot_u = dot(p, u);
         float dot_w = dot(p, w);
         float azimuth = atan2f(dot_w, dot_u);
         if (azimuth < 0) azimuth += 2 * PI_F;
         azimuth += phase;
         float t = azimuth / (2 * PI_F);
         
         return { dist - thickness, t, dist };
      }
    };

    struct DistortedRing {
       const Basis& basis;
       float radius;
       float thickness;
       std::function<float(float)> shift_fn;
       float max_distortion;
       float phase;
       
       Vector normal, u, w;
       float nx, ny, nz;
       float target_angle, center_phi;
       float max_thickness;

       DistortedRing(const Basis& b, float r, float th, std::function<float(float)> sf, float md, float ph)
         : basis(b), radius(r), thickness(th), shift_fn(sf), max_distortion(md), phase(ph)
       {
          normal = basis.v; u = basis.u; w = basis.w;
          nx = normal.i; ny = normal.j; nz = normal.k;
          target_angle = radius * (PI_F / 2.0f);
          center_phi = acosf(ny);
          max_thickness = thickness + max_distortion;
       }

       Bounds get_vertical_bounds() const {
          // Similar logic to Ring but with max_thickness
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

          int y_min = std::max(0, static_cast<int>(floorf((f_phi_min * (H_VIRT - 1)) / PI_F)));
          int y_max = std::min(H - 1, static_cast<int>(ceilf((f_phi_max * (H_VIRT - 1)) / PI_F)));
           return { y_min, y_max };
       }

       template<typename OutputIt>
       bool get_horizontal_intervals(int y, OutputIt out) const { return false; }

      DistanceResult distance(const Vector& p) const {
         float polar = angle_between(p, normal);
         float dot_u = dot(p, u);
         float dot_w = dot(p, w);
         float azimuth = atan2f(dot_w, dot_u);
         if (azimuth < 0) azimuth += 2 * PI_F;
         
         float t_val = azimuth + phase;
         float norm_t = t_val / (2 * PI_F);
         
         float shift = shift_fn(norm_t);
         float local_target = target_angle + shift;
         float dist = std::abs(polar - local_target);
         
         return { dist - thickness, azimuth / (2 * PI_F), dist };
      }
    };
    
    template <typename A, typename B>
    struct Intersection {
       const A& a;
       const B& b;
       float thickness;

       Intersection(const A& shapeA, const B& shapeB)
          : a(shapeA), b(shapeB), thickness(std::min(shapeA.thickness, shapeB.thickness)) {}

       Bounds get_vertical_bounds() const {
          auto b1 = a.get_vertical_bounds();
          auto b2 = b.get_vertical_bounds();
          return { std::max(b1.y_min, b2.y_min), std::min(b1.y_max, b2.y_max) };
       }

       template<typename OutputIt>
       bool get_horizontal_intervals(int y, OutputIt out) const {
          std::vector<std::pair<float, float>> iA;
          bool handledA = a.get_horizontal_intervals(y, [&](float t1, float t2) {
              iA.push_back({t1, t2});
          });

          std::vector<std::pair<float, float>> iB;
          bool handledB = b.get_horizontal_intervals(y, [&](float t1, float t2) {
              iB.push_back({t1, t2});
          });

          if (!handledA && !handledB) return false;

          if (!handledA) {
             for (const auto& iv : iB) out(iv.first, iv.second);
             return true;
          }
          if (!handledB) {
             for (const auto& iv : iA) out(iv.first, iv.second);
             return true;
          }

          size_t idxA = 0, idxB = 0;
          while (idxA < iA.size() && idxB < iB.size()) {
             auto ivA = iA[idxA];
             auto ivB = iB[idxB];
             
             float start = std::max(ivA.first, ivB.first);
             float end = std::min(ivA.second, ivB.second);
             
             if (start < end) {
                 out(start, end);
             }
             
             if (ivA.second < ivB.second) idxA++;
             else idxB++;
          }
          return true;
       }

       DistanceResult distance(const Vector& p) const {
          auto resA = a.distance(p);
          auto resB = b.distance(p);
          if (resA.dist > resB.dist) return resA;
          return resB;
       }
    };

    struct Polygon {
       const Basis& basis;
       float thickness;
       int sides;
       float phase;
       float apothem;
       float nx, ny, nz, R_val, alpha_angle;
       int y_min, y_max;
       
       Polygon(const Basis& b, float r, float th, int s, float ph)
         : basis(b), thickness(th), sides(s), phase(ph)
       {
          apothem = thickness * cosf(PI_F / sides);
          nx = basis.v.i; ny = basis.v.j; nz = basis.v.k;
          R_val = sqrtf(nx*nx + nz*nz);
          alpha_angle = atan2f(nx, nz);
          
          float center_phi = acosf(std::max(-1.0f, std::min(1.0f, ny)));
          float margin = thickness + 0.1f;
          y_min = std::max(0, static_cast<int>(floorf((std::max(0.0f, center_phi - margin) * (H_VIRT - 1)) / PI_F)));
          y_max = std::min(H - 1, static_cast<int>(ceilf((std::min(PI_F, center_phi + margin) * (H_VIRT - 1)) / PI_F)));
       }
       
       Bounds get_vertical_bounds() const { return { y_min, y_max }; }
       
       template<typename OutputIt>
       bool get_horizontal_intervals(int y, OutputIt out) const {
           float phi = y_to_phi(static_cast<float>(y));
           float cos_phi = cosf(phi);
           float sin_phi = sinf(phi);

           if (R_val < 0.01f) return false;

           float ang_high = thickness;
           float D_min = cosf(ang_high); 

           float denom = R_val * sin_phi;
           if (std::abs(denom) < 0.000001f) return false;

           float C_min = (D_min - ny * cos_phi) / denom;
           
           if (C_min > 1.0f) return true; // Completely outside the cap
           if (C_min < -1.0f) return false; // Full scan fallback (simplification)

           float d_alpha = acosf(C_min);
           out(alpha_angle - d_alpha, alpha_angle + d_alpha);
           return true;
       }
       
       DistanceResult distance(const Vector& p) const {
          float polar = angle_between(p, basis.v);
          float dot_u = dot(p, basis.u);
          float dot_w = dot(p, basis.w);
          float azimuth = atan2f(dot_w, dot_u);
          if (azimuth < 0) azimuth += 2 * PI_F;
          azimuth += phase;
          
          float sector = 2 * PI_F / sides;
          float local = wrap(azimuth + sector/2.0f, sector) - sector/2.0f;
          
          float dist_edge = polar * cosf(local) - apothem;
          return { dist_edge, polar / thickness, polar };
       }
    };
    
    struct Star {
       const Basis& basis;
       int sides;
       float phase;
       float thickness;
       float nx, ny, plane_d;
       float scan_nx, scan_ny, scan_nz, scan_R, scan_alpha;
       int y_min, y_max;
       
       Star(const Basis& b, float radius, int s, float ph)
         : basis(b), sides(s), phase(ph)
       {
         float outer = radius * (PI_F / 2.0f);
         float inner = outer * 0.382f;
         float angle_step = PI_F / sides;
         
         float vT = outer;
         float vVx = inner * cosf(angle_step);
         float vVy = inner * sinf(angle_step);
         float dx = vVx - vT;
         float dy = vVy;
         float len = sqrtf(dx*dx + dy*dy);
         nx = -dy / len;
         ny = dx / len;
         plane_d = -(nx * vT);
         thickness = outer;
         
         scan_nx = basis.v.i; scan_ny = basis.v.j; scan_nz = basis.v.k;
         scan_R = sqrtf(scan_nx*scan_nx + scan_nz*scan_nz);
         scan_alpha = atan2f(scan_nx, scan_nz);
         
         float center_phi = acosf(std::max(-1.0f, std::min(1.0f, basis.v.j)));
         float margin = outer + 0.1f;
         y_min = std::max(0, static_cast<int>(floorf((std::max(0.0f, center_phi - margin) * (H_VIRT - 1)) / PI_F)));
         y_max = std::min(H - 1, static_cast<int>(ceilf((std::min(PI_F, center_phi + margin) * (H_VIRT - 1)) / PI_F)));
       }
       
       Bounds get_vertical_bounds() const { return { y_min, y_max }; }

       template<typename OutputIt>
       bool get_horizontal_intervals(int y, OutputIt out) const { return false; }
       
       DistanceResult distance(const Vector& p) const {
          float polar = angle_between(p, basis.v);
          float dot_u = dot(p, basis.u);
          float dot_w = dot(p, basis.w);
          float azimuth = atan2f(dot_w, dot_u);
          if (azimuth < 0) azimuth += 2 * PI_F;
          azimuth += phase;
          
          float sector = 2 * PI_F / sides;
          float local = wrap(azimuth + sector/2.0f, sector) - sector/2.0f;
          local = std::abs(local);
          
          float px = polar * cosf(local);
          float py = polar * sinf(local);
          float dist_edge = px * nx + py * ny + plane_d;
          
          return { -dist_edge, polar / thickness, polar };
       }
    };
    
    struct Flower {
       const Basis& basis;
       int sides;
       float phase;
       float thickness;
       float apothem;
       Vector antipode;
       int y_min, y_max;
       
       Flower(const Basis& b, float radius, int s, float ph)
         : basis(b), sides(s), phase(ph)
       {
          float outer = radius * (PI_F / 2.0f);
          apothem = PI_F - outer;
          thickness = outer;
          antipode = -basis.v;
          
          float center_phi = acosf(std::max(-1.0f, std::min(1.0f, antipode.j)));
          float margin = thickness + 0.1f;
          y_min = std::max(0, static_cast<int>(floorf((std::max(0.0f, center_phi - margin) * (H_VIRT - 1)) / PI_F)));
          y_max = std::min(H - 1, static_cast<int>(ceilf((std::min(PI_F, center_phi + margin) * (H_VIRT - 1)) / PI_F)));
       }
       
       Bounds get_vertical_bounds() const { return { y_min, y_max }; }

       template<typename OutputIt>
       bool get_horizontal_intervals(int y, OutputIt out) const { return false; }
       
       DistanceResult distance(const Vector& p) const {
          float scan_dist = angle_between(p, antipode);
          float polar = PI_F - scan_dist;
          
          float dot_u = dot(p, basis.u);
          float dot_w = dot(p, basis.w);
          float azimuth = atan2f(dot_w, dot_u);
          if (azimuth < 0) azimuth += 2 * PI_F;
          azimuth += phase;
          
          float sector = 2 * PI_F / sides;
          float local = wrap(azimuth + sector/2.0f, sector) - sector/2.0f;
          
          float dist_edge = polar * cosf(local) - apothem;
          return { -dist_edge, scan_dist / thickness, scan_dist };
       }
    };
  };

  /**
   * @brief The Scan struct contains volumetric (raster) drawing primitives.
   */
  template <int W>
  struct Scan {
    
    static void rasterize(auto& pipeline, Canvas& canvas, const auto& shape, ColorFn auto color_fn, bool debug_bb = false) {
       auto bounds = shape.get_vertical_bounds();
       
       for (int y = bounds.y_min; y <= bounds.y_max; ++y) {
          bool handled = shape.get_horizontal_intervals(y, [&](float t1, float t2) {
             // Convert angle interval to pixel interval and scan
             int x1 = static_cast<int>(floorf((t1 * W) / (2 * PI_F)));
             int x2 = static_cast<int>(ceilf((t2 * W) / (2 * PI_F)));
             
             for (int x = x1; x <= x2; ++x) {
                int wx = wrap(x, W);
                process_pixel(wx, y, pipeline, canvas, shape, color_fn, debug_bb);
             }
          });
          
          if (!handled) {
             for (int x = 0; x < W; ++x) {
                process_pixel(x, y, pipeline, canvas, shape, color_fn, debug_bb);
             }
          }
       }
    }
    
    static void process_pixel(int x, int y, auto& pipeline, Canvas& canvas, const auto& shape, ColorFn auto color_fn, bool debug_bb) {
       // Wrap handled by pixel_to_vector implicitly if x is out of bounds? No, usually x must be in range.
       // The generic scan loop passes 0..W-1.
       const Vector& p = pixel_to_vector<W>(x, y);
       
       auto result = shape.distance(p);
       float d = result.dist;
       
       float pixel_width = 2.0f * PI_F / W;
       float threshold = pixel_width;
       
       if (d < threshold) {
          float t_aa = 0.5f - d / (2.0f * pixel_width);
          float alpha = quintic_kernel(std::max(0.0f, std::min(1.0f, t_aa)));
          
          if (alpha <= 0.001f) return;
          
          auto c = color_fn(p, result.t); 
          // color_fn in C++ typically returns Color4 directly (or Pixel)
          // Adjust logic to handle 3-arg color_fn if needed? 
          // Existing code uses 2 args (p, t). JS uses (p, t, rawDist).
          // We'll stick to 2 args for compatibility or wrapper.
          
          if constexpr (std::is_same_v<std::decay_t<decltype(c)>, Pixel>) {
              pipeline.plot(canvas, x, y, c.color, 0.0f, alpha); // Pixel doesn't have alpha usually?
          } else {
             pipeline.plot(canvas, x, y, c.color, 0.0f, c.alpha * alpha);
          }
       }
    }

    struct DistortedRing {
      static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, float thickness, 
                       std::function<float(float)> shift_fn, float amplitude, ColorFn auto color_fn, float phase = 0, bool debug_bb = false) 
      {
        SDF::DistortedRing shape(basis, radius, thickness, shift_fn, amplitude, phase);
        Scan::rasterize(pipeline, canvas, shape, color_fn, debug_bb);
      }
    };

    struct Point {
       static void draw(auto& pipeline, Canvas& canvas, const Vector& p, float thickness, ColorFn auto color_fn) {
          // JS Scan.Point uses Scan.Ring with radius 0
          Basis basis = make_basis(Quaternion(), p);
          Ring::draw(pipeline, canvas, basis, 0.0f, thickness, color_fn);
       }
    };
    


    struct Polygon {
       static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, int sides, ColorFn auto color_fn, float phase = 0, bool debug_bb = false) {
         Basis eff_basis = basis;
         float effective_radius = radius;
         if (radius > 1.0f) {
            eff_basis.v = -eff_basis.v;
            eff_basis.u = -eff_basis.u;
            effective_radius = 2.0f - radius;
         }
         float thickness = effective_radius * (PI_F / 2.0f);
         
         SDF::Polygon shape(eff_basis, effective_radius, thickness, sides, phase);
         Scan::rasterize(pipeline, canvas, shape, color_fn, debug_bb);
       }
    };

    struct Star {
       static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, int sides, ColorFn auto color_fn, float phase = 0, bool debug_bb = false) {
         Basis eff_basis = basis;
         float effective_radius = radius;
         if (radius > 1.0f) {
            eff_basis.v = -eff_basis.v;
            eff_basis.u = -eff_basis.u;
            effective_radius = 2.0f - radius;
         }
         SDF::Star shape(eff_basis, effective_radius, sides, phase);
         Scan::rasterize(pipeline, canvas, shape, color_fn, debug_bb);
       }
    };

    struct Flower {
       static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, int sides, ColorFn auto color_fn, float phase = 0, bool debug_bb = false) {
         Basis eff_basis = basis;
         float effective_radius = radius;
         if (radius > 1.0f) {
            eff_basis.v = -eff_basis.v;
            eff_basis.u = -eff_basis.u;
            effective_radius = 2.0f - radius;
         }
         SDF::Flower shape(eff_basis, effective_radius, sides, phase);
         Scan::rasterize(pipeline, canvas, shape, color_fn, debug_bb);
       }
    };

    struct Ring {
       static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, float thickness, ColorFn auto color_fn, float phase = 0, bool debug_bb = false) {
          Basis eff_basis = basis;
          float effective_radius = radius;
          if (radius > 1.0f) {
             eff_basis.v = -eff_basis.v; 
             effective_radius = 2.0f - radius;
          }
          SDF::Ring shape(eff_basis, effective_radius, thickness, phase);
          Scan::rasterize(pipeline, canvas, shape, color_fn, debug_bb);
       }
         // Overload for Vector normal inputs (legacy support helper)
       static void draw(auto& pipeline, Canvas& canvas, const Vector& normal, float radius, float thickness, ColorFn auto color_fn, float phase = 0, bool debug_bb = false) {
          Basis basis = make_basis(Quaternion(), normal);
          draw(pipeline, canvas, basis, radius, thickness, color_fn, phase, debug_bb);
       }
    };
    
    struct Field {
        static void draw(auto& pipeline, Canvas& canvas, ColorFn auto color_fn) {
            for (int y = 0; y < H; ++y) {
               for (int x = 0; x < W; ++x) {
                  Vector p = pixel_to_vector<W>(x, y);
                  auto c = color_fn(p, 0.0f);
                  if constexpr (std::is_same_v<std::decay_t<decltype(c)>, Pixel>) {
                      pipeline.plot(canvas, x, y, c.color, 0.0f, 1.0f);
                  } else {
                      if (c.alpha > 0.001f) {
                          pipeline.plot(canvas, x, y, c.color, 0.0f, c.alpha);
                      }
                  }
               }
            }
        }
    };
  };
