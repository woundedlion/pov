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

    static void sample(Points& points, const Quaternion& orientation, const Vector& normal, float radius, int num_samples, float phase = 0) {
      Vector ref_axis = (std::abs(dot(normal, X_AXIS)) > 0.9999f) ? Y_AXIS : X_AXIS;
      Vector v = rotate(normal, orientation).normalize();
      Vector ref = rotate(ref_axis, orientation).normalize();
      Vector u = cross(v, ref).normalize();
      Vector w = cross(v, u).normalize();

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

    static void draw(auto& pipeline, Canvas& canvas, const Quaternion& orientation, const Vector& normal, float radius, ColorFn auto color_fn, float phase = 0) {
      Points points;
      sample(points, orientation, normal, radius, W / 4, phase);
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
    static void sample(Points& points, const Quaternion& orientation, const Vector& normal, float radius, int num_sides, float phase = 0) {
       // Offset by half-sector to align with Scan.Polygon edges
       const float offset = PI_F / num_sides;
       Ring::sample(points, orientation, normal, radius, num_sides, phase + offset);
    }
    static void draw(auto& pipeline, Canvas& canvas, const Quaternion& orientation, const Vector& normal, float radius, int num_sides, ColorFn auto color_fn, float phase = 0) {
      Points points;
      sample(points, orientation, normal, radius, num_sides, phase);
      
      Vector center = rotate(normal, orientation).normalize();
      
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
    static Vector fn_point(ScalarFn auto f, const Vector& normal, float radius, float angle) {
      Vector v(normal);
      if (radius > 1) { v = -v; radius = 2 - radius; }
      Vector u = (std::abs(dot(v, X_AXIS)) > 0.99995f) ? cross(v, Y_AXIS).normalize() : cross(v, X_AXIS).normalize();
      Vector w(cross(v, u));

      auto vi = Ring::calcPoint(angle, radius, u, v, w);
      auto vp = Ring::calcPoint(angle, 1, u, v, w);
      Vector axis = cross(v, vp).normalize();
      return rotate(vi, make_rotation(axis, f(angle * PI_F / 2)));
    }

    static void sample(Points& points, const Quaternion& orientation, const Vector& normal, float radius, ScalarFn auto shift_fn, float phase = 0) {
      Vector ref_axis = (std::abs(dot(normal, X_AXIS)) > 0.9999f) ? Y_AXIS : X_AXIS;
      Vector v = rotate(normal, orientation).normalize();
      Vector ref = rotate(ref_axis, orientation).normalize();
      Vector u = cross(v, ref).normalize();
      Vector w = cross(v, u).normalize();

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

    static void draw(auto& pipeline, Canvas& canvas, const Quaternion& orientation, const Vector& normal, float radius, ScalarFn auto shift_fn, ColorFn auto color_fn, float phase = 0) {
      Points points;
      sample(points, orientation, normal, radius, shift_fn, phase);
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
    static void draw(auto& pipeline, Canvas& canvas, const Quaternion& orientation, const Vector& normal, float radius, int num_sides, ColorFn auto color_fn, float phase = 0) {
      Basis basis = make_basis(orientation, normal);
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
    static void draw(auto& pipeline, Canvas& canvas, const Quaternion& orientation, const Vector& normal, float radius, int num_sides, ColorFn auto color_fn, float phase = 0) {
      Basis basis = make_basis(orientation, normal);
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


/**
 * @brief The Scan struct contains volumetric (raster) drawing primitives.
 */
template <int W>
struct Scan {

  struct DistortedRing {
    struct Context {
      Vector normal;
      float radius;
      float thickness;
      float amplitude;
      std::function<float(float)> shift_fn;
      float nx, ny, nz;
      float target_angle, r_val, alpha, center_phi;
      Vector u, w;
      float max_thickness;
      bool debug_bb;
    };

    static void draw(auto& pipeline, Canvas& canvas, const Quaternion& orientation, const Vector& normal,
      float radius, float thickness, ScalarFn auto shift_fn, float amplitude, ColorFn auto color_fn, bool debug_bb = false)
    {
      // Conservative bounds
      float max_thickness = thickness + amplitude;

      // Basis Construction matching Plot.Ring => Stable Twist
      Vector ref_axis = (std::abs(dot(normal, X_AXIS)) > 0.9999f) ? Y_AXIS : X_AXIS;

      // Calculate Basis
      Vector v = rotate(normal, orientation).normalize();
      Vector ref = rotate(ref_axis, orientation).normalize();
      Vector u = cross(v, ref).normalize();
      Vector w = cross(v, u).normalize();

      Context ctx;
      ctx.normal = v;
      ctx.radius = radius;
      ctx.thickness = thickness;
      ctx.amplitude = amplitude;
      ctx.shift_fn = shift_fn;
      ctx.nx = v.i; ctx.ny = v.j; ctx.nz = v.k;
      ctx.target_angle = radius * (PI_F / 2.0f);
      ctx.r_val = sqrtf(ctx.nx * ctx.nx + ctx.nz * ctx.nz);
      ctx.alpha = atan2f(ctx.nx, ctx.nz);
      ctx.center_phi = acosf(ctx.ny);
      ctx.u = u;
      ctx.w = w;
      ctx.max_thickness = max_thickness;
      ctx.debug_bb = debug_bb;

      float a1 = ctx.center_phi - ctx.target_angle;
      float a2 = ctx.center_phi + ctx.target_angle;
      float p1 = acosf(cosf(a1));
      float p2 = acosf(cosf(a2));
      float min_p = std::min(p1, p2);
      float max_p = std::max(p1, p2);

      float phi_min = std::max(0.0f, min_p - max_thickness);
      float phi_max = std::min(PI_F, max_p + max_thickness);

      if (phi_min > phi_max) return;

      int y_min = std::max(0, static_cast<int>(floorf((phi_min * (H_VIRT - 1)) / PI_F)));
      int y_max = std::min(H - 1, static_cast<int>(ceilf((phi_max * (H_VIRT - 1)) / PI_F)));

      for (int y = y_min; y <= y_max; y++) {
        scan_row(pipeline, canvas, y, ctx, color_fn);
      }
    }

  private:
    static void scan_row(auto& pipeline, Canvas& canvas, int y, const Context& ctx, ColorFn auto& color_fn) {
      float phi = y_to_phi(static_cast<float>(y));
      float cos_phi = cosf(phi);
      float sin_phi = sinf(phi);

      float ang_low = std::max(0.0f, ctx.target_angle - ctx.max_thickness);
      float ang_high = std::min(PI_F, ctx.target_angle + ctx.max_thickness);
      float D_max = cosf(ang_low);
      float D_min = cosf(ang_high);

      if (ctx.r_val < 0.01f) {
        scan_full_row(pipeline, canvas, y, ctx, color_fn);
        return;
      }

      float denom = ctx.r_val * sin_phi;
      if (std::abs(denom) < 0.000001f) {
        scan_full_row(pipeline, canvas, y, ctx, color_fn);
        return;
      }

      float C_min = (D_min - ctx.ny * cos_phi) / denom;
      float C_max = (D_max - ctx.ny * cos_phi) / denom;
      float min_cos = std::max(-1.0f, C_min);
      float max_cos = std::min(1.0f, C_max);

      if (min_cos > max_cos) return;

      float angle_min = acosf(max_cos);
      float angle_max = acosf(min_cos);

      const float pixel_width = 2.0f * PI_F / W;
      const float safe_threshold = pixel_width;

      if (angle_min <= safe_threshold) {
        scan_window(pipeline, canvas, y, ctx.alpha - angle_max, ctx.alpha + angle_max, ctx, color_fn);
      }
      else if (angle_max >= PI_F - safe_threshold) {
        scan_window(pipeline, canvas, y, ctx.alpha + angle_min, ctx.alpha + 2 * PI_F - angle_min, ctx, color_fn);
      }
      else {
        scan_window(pipeline, canvas, y, ctx.alpha - angle_max, ctx.alpha - angle_min, ctx, color_fn);
        scan_window(pipeline, canvas, y, ctx.alpha + angle_min, ctx.alpha + angle_max, ctx, color_fn);
      }
    }

    static void scan_full_row(auto& pipeline, Canvas& canvas, int y, const Context& ctx, ColorFn auto& color_fn) {
      for (int x = 0; x < W; x++) {
        process_pixel(pipeline, canvas, x, y, ctx, color_fn);
      }
    }

    static void scan_window(auto& pipeline, Canvas& canvas, int y, float t1, float t2, const Context& ctx, ColorFn auto& color_fn) {
      int x1 = static_cast<int>(floorf((t1 * W) / (2 * PI_F)));
      int x2 = static_cast<int>(ceilf((t2 * W) / (2 * PI_F)));
      for (int x = x1; x <= x2; x++) {
        int wx = wrap(x, W);
        process_pixel(pipeline, canvas, wx, y, ctx, color_fn);
      }
    }

    static void process_pixel(auto& pipeline, Canvas& canvas, int x, int y, const Context& ctx, ColorFn auto& color_fn) {
      const Vector& p = pixel_to_vector<W>(x, y);

      float polar_angle = angle_between(p, ctx.normal);
      float dot_u = dot(p, ctx.u);
      float dot_w = dot(p, ctx.w);
      float azimuth = atan2f(dot_w, dot_u);
      if (azimuth < 0) azimuth += 2 * PI_F;

      float norm_azimuth = azimuth / (2 * PI_F);
      float shift = ctx.shift_fn(norm_azimuth);
      float local_target = ctx.target_angle + shift;

      float dist = std::abs(polar_angle - local_target);

      if (dist < ctx.thickness) {
        float dist_t = dist / ctx.thickness;
        float aa_alpha = quintic_kernel(1.0f - dist_t);
        if (aa_alpha <= 0.001f) return;
        Color4 c = color_fn(p, norm_azimuth);        
        pipeline.plot(canvas, x, y, c.color, 0.0f, c.alpha * aa_alpha);
      }
    }
  };

  struct Ring {
    struct Context {
      Vector normal;
      float radius;
      float thickness;
      float nx, ny, nz;
      float target_angle, r_val, alpha, center_phi;
      Vector u, w;
      float start_angle, end_angle;
      bool check_sector;
      const std::span<const Vector> clip_planes;
      float cos_min, cos_max, inv_sin_target;
    };

    static void draw(auto& pipeline, Canvas& canvas, const Vector& normal, float radius, float thickness, ColorFn auto color_fn,
      float start_angle = 0, float end_angle = 2 * PI_F, std::span<const Vector> clip_planes = {})
    {
      Vector n = normal;
      if (radius > 1.0f) {
        n = -n;
        radius = 2.0f - radius;
      }

      Context ctx;
      ctx.normal = n;
      ctx.radius = radius;
      ctx.thickness = thickness;
      ctx.target_angle = radius * (PI_F / 2.0f);
      
      // Optimizations
      float min_angle = std::max(0.0f, ctx.target_angle - thickness);
      float max_angle = std::min(PI_F, ctx.target_angle + thickness);
      ctx.cos_min = cosf(min_angle); 
      ctx.cos_max = cosf(max_angle);
      ctx.inv_sin_target = (std::abs(sinf(ctx.target_angle)) > 0.001f) ? (1.0f / sinf(ctx.target_angle)) : 0.0f;

      // Bounding box setup...
      ctx.nx = n.i; ctx.ny = n.j; ctx.nz = n.k;
      ctx.r_val = sqrtf(ctx.nx * ctx.nx + ctx.nz * ctx.nz);
      ctx.alpha = atan2f(ctx.nx, ctx.nz);
      ctx.center_phi = acosf(ctx.ny);

      Vector ref = (std::abs(dot(n, X_AXIS)) > 0.9999f) ? Y_AXIS : X_AXIS;
      ctx.u = cross(n, ref).normalize();
      ctx.w = cross(n, ctx.u).normalize();
      
      ctx.start_angle = start_angle;
      ctx.end_angle = end_angle;
      ctx.check_sector = (std::abs(end_angle - start_angle) < (2 * PI_F - 0.001f));
      ctx.clip_planes = clip_planes;

      float a1 = ctx.center_phi - ctx.target_angle;
      float a2 = ctx.center_phi + ctx.target_angle;
      float p1 = acosf(cosf(a1));
      float p2 = acosf(cosf(a2));
      float phi_min = std::max(0.0f, std::min(p1, p2) - thickness);
      float phi_max = std::min(PI_F, std::max(p1, p2) + thickness);

      int y_min = std::max(0, static_cast<int>(floorf((phi_min * (H_VIRT - 1)) / PI_F)));
      int y_max = std::min(H - 1, static_cast<int>(ceilf((phi_max * (H_VIRT - 1)) / PI_F)));

      for (int y = y_min; y <= y_max; y++) {
        scan_full(pipeline, canvas, y, ctx, color_fn);
      }
    }

  private:

    static void scan_full(auto& pipeline, Canvas& canvas, int y, const Context& ctx, ColorFn auto& color_fn) {
      for (int x = 0; x < W; ++x) {
        process_pixel(pipeline, canvas, x, y, ctx, color_fn);
      }
    }

    static void process_pixel(auto& pipeline, Canvas& canvas, int x, int y, const Context& ctx, ColorFn auto& color_fn) {
      const Vector& p = pixel_to_vector<W>(x, y);

      for (const auto& cp : ctx.clip_planes) {
        if (dot(p, cp) < 0) return;
      }

      float cos_angle = dot(p, ctx.normal);
      if (cos_angle > ctx.cos_min || cos_angle < ctx.cos_max) return;

      float polar_angle = acosf(cos_angle);
      float dist = std::abs(polar_angle - ctx.target_angle);
      if (dist < ctx.thickness) {
        float t = 0;
        float dot_u = dot(p, ctx.u);
        float dot_w = dot(p, ctx.w);
        float azimuth = atan2f(dot_w, dot_u);
        if (azimuth < 0) azimuth += 2 * PI_F;
        t = azimuth / (2 * PI_F);

        if (ctx.check_sector) {
          bool inside = (ctx.start_angle <= ctx.end_angle)
            ? (azimuth >= ctx.start_angle && azimuth <= ctx.end_angle)
            : (azimuth >= ctx.start_angle || azimuth <= ctx.end_angle);
          if (!inside) return;
        }

        float thickness_t = dist / ctx.thickness;
        float alpha_factor = quintic_kernel(1.0f - thickness_t);
        if (alpha_factor <= 0.001f) return;

        Color4 c = color_fn(p, t);
        pipeline.plot(canvas, x, y, c.color, 0.0f, c.alpha * alpha_factor);
      }
    }
  };

  struct Line {
    static void draw(auto& pipeline, Canvas& canvas, const Vector& v1, const Vector& v2, float thickness, ColorFn auto color_fn) {
      Vector normal = cross(v1, v2).normalize();
      if (dot(normal, normal) < 1e-6) return;
      Vector c1 = cross(normal, v1);
      Vector c2 = cross(v2, normal);
      std::array<Vector, 2> clips = { c1, c2 };
      Ring::draw(pipeline, canvas, normal, 1.0f, thickness, color_fn, 0, 2 * PI_F, clips);
    }
  };

  struct Point {
    static void draw(auto& pipeline, Canvas& canvas, const Vector& p, float thickness, ColorFn auto color_fn) {
      Ring::draw(pipeline, canvas, p, 0.0f, thickness, color_fn);
    }
  };

  struct Polygon {
    struct Context {
      Vector normal;
      float radius;
      float thickness;
      int sides;
      float apothem;
      float nx, ny, nz;
      float R, alpha;
      Vector u, w;
      float pixel_width;
      bool debug_bb;
    };

    static void draw(auto& pipeline, Canvas& canvas, const Quaternion& orientation, const Vector& normal,
      float radius, int sides, ColorFn auto color_fn, float phase = 0, bool debug_bb = false)
    {
      float thickness = radius * (PI_F / 2.0f);
      if (thickness <= 0.0001f) return;

      // Basis Construction
      Basis basis = make_basis(orientation, normal);
      Vector v = basis.v; Vector u = basis.u; Vector w = basis.w;

      // Backside Logic
      if (radius > 1.0f) {
        v = -v;
        u = -u;
        radius = 2.0f - radius;
        thickness = radius * (PI_F / 2.0f);
      }

      float nx = v.i;
      float ny = v.j;
      float nz = v.k;

      float R = sqrtf(nx * nx + nz * nz);
      float alpha = atan2f(nx, nz);

      float angle = PI_F / sides;
      float apothem = thickness * cosf(angle);

      float pixel_width = 2.0f * PI_F / W;

      Context ctx;
      ctx.normal = v;
      ctx.radius = radius;
      ctx.thickness = thickness;
      ctx.sides = sides;
      ctx.apothem = apothem;
      ctx.nx = nx; ctx.ny = ny; ctx.nz = nz;
      ctx.R = R; ctx.alpha = alpha;
      ctx.u = u; ctx.w = w;
      ctx.pixel_width = pixel_width;
      ctx.debug_bb = debug_bb;
      
      // Update: Scan::Polygon needs to respect phase!
      // But Scan::Polygon uses global scan sweep (y->phi) and then checks azimuth within the sector.
      // The `azimuth` calculation is `atan2(dot_w, dot_u)`.
      // The basis {u, v, w} comes from `make_basis(orientation, normal)`.
      // If we add `phase`, we effectively rotate the polygon around the normal or shift the azimuth check.
      // In JS, Twist rotates the polygon. Here, we can just add phase to azimuth check in process_pixel.
      // But process_pixel is static inside struct. We need to store phase in Context.
      // ... WAIT, context definition is separate.
    }

  private:
    static void scan_row(auto& pipeline, Canvas& canvas, int y, const Context& ctx, ColorFn auto& color_fn) {
      float phi = y_to_phi(static_cast<float>(y));
      float cos_phi = cosf(phi);
      float sin_phi = sinf(phi);

      if (ctx.R < 0.01f) {
        scan_full_row(pipeline, canvas, y, ctx, color_fn);
        return;
      }

      float ang_high = ctx.thickness;
      float D_min = cosf(ang_high); // cos(thickness)

      float denom = ctx.R * sin_phi;
      if (std::abs(denom) < 0.000001f) {
        scan_full_row(pipeline, canvas, y, ctx, color_fn);
        return;
      }

      float C_min = (D_min - ctx.ny * cos_phi) / denom;
      
      if (C_min > 1.0f) return; // Completely outside the cap
      if (C_min < -1.0f) {
        scan_full_row(pipeline, canvas, y, ctx, color_fn);
        return;
      }

      float d_alpha = acosf(C_min);
      scan_window(pipeline, canvas, y, ctx.alpha - d_alpha, ctx.alpha + d_alpha, ctx, color_fn);
    }

    static void scan_full_row(auto& pipeline, Canvas& canvas, int y, const Context& ctx, ColorFn auto& color_fn) {
      for (int x = 0; x < W; x++) {
        process_pixel(pipeline, canvas, x, y, ctx, color_fn);
      }
    }

    static void scan_window(auto& pipeline, Canvas& canvas, int y, float t1, float t2, const Context& ctx, ColorFn auto& color_fn) {
      int x1 = static_cast<int>(floorf((t1 * W) / (2 * PI_F)));
      int x2 = static_cast<int>(ceilf((t2 * W) / (2 * PI_F)));

      if (x2 - x1 >= W) {
        scan_full_row(pipeline, canvas, y, ctx, color_fn);
        return;
      }

      for (int x = x1; x <= x2; x++) {
        int wx = wrap(x, W);
        process_pixel(pipeline, canvas, wx, y, ctx, color_fn);
      }
    }

    static void process_pixel(auto& pipeline, Canvas& canvas, int x, int y, const Context& ctx, ColorFn auto& color_fn) {
      const Vector& p = pixel_to_vector<W>(x, y);
      
      float polar_angle = angle_between(p, ctx.normal);
      if (polar_angle > ctx.thickness + ctx.pixel_width) return;

      float dot_u = dot(p, ctx.u);
      float dot_w = dot(p, ctx.w);
      float azimuth = atan2f(dot_w, dot_u);
      if (azimuth < 0) azimuth += 2 * PI_F;

      // SDF Logic
      float sector_angle = 2 * PI_F / ctx.sides;
      float local_azimuth = wrap(azimuth + sector_angle / 2.0f, sector_angle) - sector_angle / 2.0f;
      float dist_to_edge = polar_angle * cosf(local_azimuth) - ctx.apothem;

      if (dist_to_edge < ctx.pixel_width) {
        float alpha = 1.0f;
        if (dist_to_edge > -ctx.pixel_width) {
          float t = (dist_to_edge + ctx.pixel_width) / (2.0f * ctx.pixel_width);
          alpha = quintic_kernel(1.0f - t);
        }

        if (alpha <= 0.001f) return;

        float t = polar_angle / ctx.thickness; 
        
        Color4 c = color_fn(p, t);
        
        pipeline.plot(canvas, x, y, c.color, 0.0f, c.alpha * alpha);
      }
    }
  };

  struct Circle {
    static void draw(auto& pipeline, Canvas& canvas, const Quaternion& orientation, const Vector& normal, float radius, ColorFn auto color_fn, float phase = 0, bool debug_bb = false) {
      float thickness = radius * (PI_F / 2.0f);
      // Ring with radius 0 and thickness = radius
      Ring::draw(pipeline, canvas, rotate(normal, orientation), 0.0f, thickness, color_fn); 
    }
  };

  struct Star {
    struct Context {
      Vector scan_normal, u, w;
      float thickness;
      int sides;
      float nx, ny, plane_d; // Edge normal precalc
      float pixel_width;
      float phase;
      bool debug_bb;
    };

    static void draw(auto& pipeline, Canvas& canvas, const Quaternion& orientation, const Vector& normal, float radius, int sides, ColorFn auto color_fn, float phase = 0, bool debug_bb = false) {
      Basis basis = make_basis(orientation, normal);
      Vector v = basis.v; Vector u = basis.u; Vector w = basis.w;

      if (radius > 1.0f) {
        v = -v; u = -u; radius = 2.0f - radius;
      }

      float outer_radius = radius * (PI_F / 2.0f);
      float inner_radius = outer_radius * 0.382f;
      float angle_step = PI_F / sides;

      // Precompute Edge Normal for SDF
      float vT = outer_radius;
      float vVx = inner_radius * cosf(angle_step);
      float vVy = inner_radius * sinf(angle_step);
      float dx = vVx - vT;
      float dy = vVy;
      float len = sqrtf(dx*dx + dy*dy);
      float nx = -dy / len;
      float ny = dx / len;
      float plane_d = -(nx * vT); 

      float thickness = outer_radius;
      if (thickness <= 0.0001f) return;
      float pixel_width = 2.0f * PI_F / W;

      Context ctx = { v, u, w, thickness, sides, nx, ny, plane_d, pixel_width, phase, debug_bb };

      // Scan bounding box around normal
      float center_phi = acosf(v.j);
      float phi_min = std::max(0.0f, center_phi - thickness - pixel_width);
      float phi_max = std::min(PI_F, center_phi + thickness + pixel_width);

      int y_min = std::max(0, static_cast<int>(floorf((phi_min * (H_VIRT - 1)) / PI_F)));
      int y_max = std::min(H - 1, static_cast<int>(ceilf((phi_max * (H_VIRT - 1)) / PI_F)));

      for (int y = y_min; y <= y_max; y++) {
        for (int x = 0; x < W; x++) {
          process_pixel(pipeline, canvas, x, y, ctx, color_fn);
        }
      }
    }

    static void process_pixel(auto& pipeline, Canvas& canvas, int x, int y, const Context& ctx, ColorFn auto& color_fn) {
      const Vector& p = pixel_to_vector<W>(x, y);
      float scan_dist = angle_between(p, ctx.scan_normal);
      if (scan_dist > ctx.thickness + ctx.pixel_width) return;

      float dot_u = dot(p, ctx.u);
      float dot_w = dot(p, ctx.w);
      float azimuth = atan2f(dot_w, dot_u);
      if (azimuth < 0) azimuth += 2 * PI_F;
      azimuth += ctx.phase;

      float sector_angle = 2 * PI_F / ctx.sides;
      float local_azimuth = wrap(azimuth + sector_angle / 2.0f, sector_angle) - sector_angle / 2.0f;
      local_azimuth = std::abs(local_azimuth);

      float px = scan_dist * cosf(local_azimuth);
      float py = scan_dist * sinf(local_azimuth);

      // SDF distance
      float dist_to_edge = px * ctx.nx + py * ctx.ny + ctx.plane_d;

      if (dist_to_edge > -ctx.pixel_width) {
        float alpha = 1.0f;
        if (dist_to_edge < ctx.pixel_width) {
          float t = (dist_to_edge + ctx.pixel_width) / (2.0f * ctx.pixel_width);
          alpha = quintic_kernel(t);
        }
        if (alpha <= 0.001f) return;

        float t = scan_dist / ctx.thickness;
        Color4 c = color_fn(p, t);
        pipeline.plot(canvas, x, y, c.color, 0.0f, c.alpha * alpha);
      }
    }
  };

  struct Flower {
    struct Context {
      Vector scan_normal, u, w;
      float thickness;
      int sides;
      float apothem;
      float pixel_width;
      float phase;
      bool debug_bb;
    };

    static void draw(auto& pipeline, Canvas& canvas, const Quaternion& orientation, const Vector& normal, float radius, int sides, ColorFn auto color_fn, float phase = 0, bool debug_bb = false) {
      Basis basis = make_basis(orientation, normal);
      Vector v = basis.v; Vector u = basis.u; Vector w = basis.w;

      if (radius > 1.0f) {
        v = -v; u = -u; radius = 2.0f - radius;
      }

      float desired_outer_radius = radius * (PI_F / 2.0f);
      float apothem = PI_F - desired_outer_radius;
      float thickness = desired_outer_radius;
      
      if (thickness <= 0.0001f) return;

      // Invert normal to scan relative to antipode
      Vector scan_v = -v;
      float pixel_width = 2.0f * PI_F / W;

      Context ctx = { scan_v, u, w, thickness, sides, apothem, pixel_width, phase, debug_bb };

      float center_phi = acosf(scan_v.j);
      float phi_min = std::max(0.0f, center_phi - thickness - pixel_width);
      float phi_max = std::min(PI_F, center_phi + thickness + pixel_width);

      int y_min = std::max(0, static_cast<int>(floorf((phi_min * (H_VIRT - 1)) / PI_F)));
      int y_max = std::min(H - 1, static_cast<int>(ceilf((phi_max * (H_VIRT - 1)) / PI_F)));

      for (int y = y_min; y <= y_max; y++) {
        for (int x = 0; x < W; x++) {
          process_pixel(pipeline, canvas, x, y, ctx, color_fn);
        }
      }
    }

    static void process_pixel(auto& pipeline, Canvas& canvas, int x, int y, const Context& ctx, ColorFn auto& color_fn) {
      const Vector& p = pixel_to_vector<W>(x, y);
      float scan_dist = angle_between(p, ctx.scan_normal);
      if (scan_dist > ctx.thickness + ctx.pixel_width) return;

      // Distance from Antipode
      float polar_angle = PI_F - scan_dist; 

      float dot_u = dot(p, ctx.u);
      float dot_w = dot(p, ctx.w);
      float azimuth = atan2f(dot_w, dot_u);
      if (azimuth < 0) azimuth += 2 * PI_F;
      azimuth += ctx.phase;

      float sector_angle = 2 * PI_F / ctx.sides;
      float local_azimuth = wrap(azimuth + sector_angle / 2.0f, sector_angle) - sector_angle / 2.0f;

      // SDF vs Antipodal Polygon
      float dist_to_edge = polar_angle * cosf(local_azimuth) - ctx.apothem;

      // If distToEdge > 0, we are OUTSIDE the polygon at S (forming hole at N)
      if (dist_to_edge > -ctx.pixel_width) {
        float alpha = 1.0f;
        if (dist_to_edge < ctx.pixel_width) {
          float t = (dist_to_edge + ctx.pixel_width) / (2.0f * ctx.pixel_width);
          alpha = quintic_kernel(t);
        }
        if (alpha <= 0.001f) return;

        float t = scan_dist / ctx.thickness;
        Color4 c = color_fn(p, t);
        pipeline.plot(canvas, x, y, c.color, 0.0f, c.alpha * alpha);
      }
    }
  };

  struct Field {
    static void draw(auto& pipeline, Canvas& canvas, ColorFn auto color_fn) {
      for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
          Vector p = pixel_to_vector<W>(x, y);
          Color4 c = color_fn(p, 0.0f);
          if (c.alpha > 0.001f) {
            pipeline.plot(canvas, x, y, c.color, 0.0f, c.alpha);
          }
        }
      }
    }
  };

};