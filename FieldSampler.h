/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <cmath>
#include <algorithm>
#include <vector>
#include <array>
#include "geometry.h"
#include "color.h"
#include "filter.h" 
#include "led.h" // For Canvas, H, W

// Context structure to avoid passing many arguments recursively
struct FSRingContext {
  Vector normal;
  float radius;
  Color4 color;
  float thickness;

  // Derived properties
  float nx, ny, nz;
  float target_angle;
  float r; // sqrt(nx^2 + nz^2)
  float alpha; // atan2(nx, nz)
  float center_phi;
  Vector u, w; // Basis vectors

  // Sector / Arc properties
  float start_angle;
  float end_angle;
  bool check_sector;
  const std::vector<Vector>* clip_planes;
};

/**
 * @brief Stateless helper class for rendering "thick" rings and arcs on the spherical display.
 */
template <int W>
class FSRing {
public:
  /**
   * @brief Draws a ring (or partial arc) onto the canvas.
   */
  static void draw(Canvas& canvas, const Vector& normal, float radius, const Color4& color,
    float thickness, float start_angle = 0, float end_angle = 2 * PI_F,
    const std::vector<Vector>* clip_planes = nullptr,
    float min_phi_limit = -1.0f, float max_phi_limit = -1.0f)
  {
    // Construct Basis
    Vector ref = X_AXIS;
    if (std::abs(dot(normal, ref)) > 0.9999f) {
      ref = Y_AXIS;
    }
    Vector u = cross(normal, ref).normalize();
    Vector w = cross(normal, u).normalize();

    // Setup Context
    FSRingContext ctx;
    ctx.normal = normal;
    ctx.radius = radius;
    ctx.color = color;
    ctx.thickness = thickness;
    ctx.nx = normal.i;
    ctx.ny = normal.j;
    ctx.nz = normal.k;
    ctx.target_angle = radius * (PI_F / 2.0f);
    ctx.r = sqrtf(ctx.nx * ctx.nx + ctx.nz * ctx.nz);
    ctx.alpha = atan2f(ctx.nx, ctx.nz);
    ctx.center_phi = acosf(ctx.ny);
    ctx.u = u;
    ctx.w = w;
    ctx.start_angle = start_angle;
    ctx.end_angle = end_angle;
    ctx.check_sector = (std::abs(end_angle - start_angle) < (2 * PI_F - 0.001f));
    ctx.clip_planes = clip_planes;

    // Calculate Vertical Bounds
    float a1 = ctx.center_phi - ctx.target_angle;
    float a2 = ctx.center_phi + ctx.target_angle;
    float p1 = acosf(cosf(a1));
    float p2 = acosf(cosf(a2));

    float min_p = std::min(p1, p2);
    float max_p = std::max(p1, p2);

    float phi_min = std::max(0.0f, min_p - thickness);
    float phi_max = std::min(PI_F, max_p + thickness);

    if (min_phi_limit >= 0) phi_min = std::max(phi_min, min_phi_limit);
    if (max_phi_limit >= 0) phi_max = std::min(phi_max, max_phi_limit);

    if (phi_min > phi_max) return;

    int y_min = std::max(0, static_cast<int>(floorf((phi_min * (H_VIRT - 1)) / PI_F)));
    int y_max = std::min(H - 1, static_cast<int>(ceilf((phi_max * (H_VIRT - 1)) / PI_F)));

    for (int y = y_min; y <= y_max; y++) {
      scanRow(canvas, y, ctx);
    }
  }

private:
  static void scanRow(Canvas& canvas, int y, const FSRingContext& ctx) {
    const Vector& p = pixel_to_vector<W>(0, y);
    float cos_phi = p.k;
    float sin_phi = p.i;

    // Singularity (Poles or Vertical Normal)
    if (ctx.r < 0.01f) {
      scanFullRow(canvas, y, ctx);
      return;
    }

    // General Intersection
    float ang_low = std::max(0.0f, ctx.target_angle - ctx.thickness);
    float ang_high = std::min(PI_F, ctx.target_angle + ctx.thickness);

    float d_max = cosf(ang_low);
    float d_min = cosf(ang_high);

    float denom = ctx.r * sin_phi;
    if (std::abs(denom) < 0.000001f) {
      scanFullRow(canvas, y, ctx);
      return;
    }

    float c_min = (d_min - ctx.ny * cos_phi) / denom;
    float c_max = (d_max - ctx.ny * cos_phi) / denom;

    float min_cos = std::max(-1.0f, c_min);
    float max_cos = std::min(1.0f, c_max);
    if (min_cos > max_cos) return;

    float angle_min = acosf(max_cos);
    float angle_max = acosf(min_cos);

    // Generate scan windows
    if (angle_min <= 0.0001f) {
      scanWindow(canvas, y, ctx.alpha - angle_max, ctx.alpha + angle_max, ctx);
    }
    else if (angle_max >= PI_F - 0.0001f) {
      scanWindow(canvas, y, ctx.alpha + angle_min, ctx.alpha + 2 * PI_F - angle_min, ctx);
    }
    else {
      scanWindow(canvas, y, ctx.alpha - angle_max, ctx.alpha - angle_min, ctx);
      scanWindow(canvas, y, ctx.alpha + angle_min, ctx.alpha + angle_max, ctx);
    }
  }

  static void scanFullRow(Canvas& canvas, int y, const FSRingContext& ctx) {
    for (int x = 0; x < W; x++) {
      processPixel(canvas, x, y, ctx);
    }
  }

  static void scanWindow(Canvas& canvas, int y, float t1, float t2, const FSRingContext& ctx) {
    int x1 = static_cast<int>(floorf((t1 * W) / (2 * PI_F)));
    int x2 = static_cast<int>(ceilf((t2 * W) / (2 * PI_F)));

    for (int x = x1; x <= x2; x++) {
      int wx = wrap(x, W);
      processPixel(canvas, wx, y, ctx);
    }
  }

  static void processPixel(Canvas& canvas, int x, int y, const FSRingContext& ctx) {
    const Vector& p = pixel_to_vector<W>(x, y);

    // Clipping Planes
    if (ctx.clip_planes) {
      for (const auto& cp : *ctx.clip_planes) {
        if (dot(p, cp) < 0) return;
      }
    }

    float polar_angle = angle_between(p, ctx.normal);
    float dist = std::abs(polar_angle - ctx.target_angle);

    if (dist < ctx.thickness) {
      // Sector Check
      if (ctx.check_sector) {
        float dot_u = dot(p, ctx.u);
        float dot_w = dot(p, ctx.w);
        float azimuth = atan2f(dot_w, dot_u);
        if (azimuth < 0) azimuth += 2 * PI_F;

        bool inside = false;
        if (ctx.start_angle <= ctx.end_angle) {
          inside = (azimuth >= ctx.start_angle && azimuth <= ctx.end_angle);
        }
        else {
          inside = (azimuth >= ctx.start_angle || azimuth <= ctx.end_angle);
        }
        if (!inside) return;
      }

      // Render
      float t = dist / ctx.thickness;
      float alpha_factor = quintic_kernel(1.0f - t);
      float final_alpha = ctx.color.alpha * alpha_factor;

      Pixel& out_color = canvas(XY(x, y));
      out_color = blend_alpha(final_alpha)(out_color, ctx.color.color);
    }
  }
};

/**
 * @brief Stateless helper for rendering anti-aliased dots (comets).
 */
template <int W>
class FSPoint {
public:
  static void draw(Canvas& canvas, const Vector& pos, const Color4& color, float thickness) {
    FSRing<W>::draw(canvas, pos, 0.0f, color, thickness, 0, 2 * PI_F);
  }
};

/**
 * @brief Stateless helper for rendering geodesic lines.
 */
template <int W>
class FSLine {
public:
  static void draw(Canvas& canvas, const Vector& v1, const Vector& v2, const Color4& color, float thickness) {
    Vector normal = cross(v1, v2).normalize();
    if (dot(normal, normal) < 0.000001f) return;

    // Clipping Planes
    Vector c1 = cross(normal, v1);
    Vector c2 = cross(v2, normal);
    std::vector<Vector> clips = { c1, c2 };

    float max_y = std::max(v1.j, v2.j);
    float min_y = std::min(v1.j, v2.j);

    Vector apex_plane_normal = cross(normal, Y_AXIS);
    if (dot(apex_plane_normal, apex_plane_normal) > 0.0001f) {
      float d1 = dot(v1, apex_plane_normal);
      float d2 = dot(v2, apex_plane_normal);

      // Segment crosses the apex plane
      if (d1 * d2 <= 0) {
        float global_max_y = sqrtf(1.0f - normal.j * normal.j);
        if (v1.j + v2.j > 0) {
          max_y = global_max_y;
        }
        else {
          min_y = -global_max_y;
        }
      }
    }

    float clamped_max_y = std::clamp(max_y, -1.0f, 1.0f);
    float clamped_min_y = std::clamp(min_y, -1.0f, 1.0f);

    float min_phi = acosf(clamped_max_y) - thickness;
    float max_phi = acosf(clamped_min_y) + thickness;

    FSRing<W>::draw(canvas, normal, 1.0f, color, thickness, 0, 2 * PI_F, &clips, min_phi, max_phi);
  }
};

/**
 * @brief Main FieldSampler class acting as a wrapper for drawing operations.
 */
template <int W>
class FieldSampler {
public:
  FieldSampler() {}

  void drawPoints(Canvas& canvas, const std::vector<Dot>& points, float thickness) {
    for (const auto& pt : points) {
      FSPoint<W>::draw(canvas, pt.position, pt.color, thickness);
    }
  }

  void drawLine(Canvas& canvas, const Vector& v1, const Vector& v2, const Color4& color, float thickness) {
    FSLine<W>::draw(canvas, v1, v2, color, thickness);
  }

  void drawRing(Canvas& canvas, const Vector& normal, float radius, const Color4& color, float thickness, float start_angle = 0, float end_angle = 2 * PI_F) {
    FSRing<W>::draw(canvas, normal, radius, color, thickness, start_angle, end_angle);
  }
};
