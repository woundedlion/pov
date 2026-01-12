/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#include <vector>
#include <algorithm>
#include <cmath>
#include <array>
#include "geometry.h"
#include "color.h"
#include "led.h" // For H, W, H_VIRT
#include <concepts>

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

  struct Polygon {
    static void sample(Points& points, const Quaternion& orientation, const Vector& normal, float radius, int num_sides, float phase = 0) {
      Ring::sample(points, orientation, normal, radius, num_sides, phase);
    }
    static void draw(auto& pipeline, Canvas& canvas, const Quaternion& orientation, const Vector& normal, float radius, int num_sides, ColorFn auto color_fn, float phase = 0) {
      Points points;
      sample(points, orientation, normal, radius, num_sides, phase);
      rasterize(pipeline, canvas, points, color_fn, true);
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


};


/**
 * @brief The Scan struct contains volumetric (raster) drawing primitives.
 */
template <int W>
struct Scan {

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
      const std::vector<Vector>* clip_planes;
    };

    static void draw(auto& pipeline, Canvas& canvas, const Vector& normal, float radius, float thickness, ColorFn auto color_fn,
      float start_angle = 0, float end_angle = 2 * PI_F, const std::vector<Vector>* clip_planes = nullptr)
    {
      Context ctx;
      ctx.normal = normal;
      ctx.radius = radius;
      ctx.thickness = thickness;
      ctx.nx = normal.i; ctx.ny = normal.j; ctx.nz = normal.k;
      ctx.target_angle = radius * (PI_F / 2.0f);
      ctx.r_val = sqrtf(ctx.nx * ctx.nx + ctx.nz * ctx.nz);
      ctx.alpha = atan2f(ctx.nx, ctx.nz);
      ctx.center_phi = acosf(ctx.ny);

      Vector ref = (std::abs(dot(normal, X_AXIS)) > 0.9999f) ? Y_AXIS : X_AXIS;
      ctx.u = cross(normal, ref).normalize();
      ctx.w = cross(normal, ctx.u).normalize();

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

      if (ctx.clip_planes) {
        for (const auto& cp : *ctx.clip_planes) if (dot(p, cp) < 0) return;
      }

      float polar_angle = angle_between(p, ctx.normal);
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
      std::vector<Vector> clips = { c1, c2 };
      Ring::draw(pipeline, canvas, normal, 1.0f, thickness, color_fn, 0, 2 * PI_F, &clips);
    }
  };

  struct Point {
    static void draw(auto& pipeline, Canvas& canvas, const Vector& p, float thickness, ColorFn auto color_fn) {
      Ring::draw(pipeline, canvas, p, 0.0f, thickness, color_fn);
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