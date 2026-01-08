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

 // Helpers
inline float y_to_phi(int y) {
  return (y * PI_F) / (H_VIRT - 1);
}

// Context structure to avoid passing many arguments recursively
struct FSRingContext {
  Vector normal;
  float radius;
  Color4 color;
  float thickness;

  // Derived properties
  float nx, ny, nz;
  float targetAngle;
  float R; // sqrt(nx^2 + nz^2)
  float alpha; // atan2(nx, nz)
  float centerPhi;
  Vector u, w; // Basis vectors

  // Sector / Arc properties
  float startAngle;
  float endAngle;
  bool checkSector;
  const std::vector<Vector>* clipPlanes;
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
    float thickness, float startAngle = 0, float endAngle = 2 * PI_F,
    const std::vector<Vector>* clipPlanes = nullptr,
    float minPhiLimit = -1.0f, float maxPhiLimit = -1.0f)
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
    ctx.targetAngle = radius * (PI_F / 2.0f);
    ctx.R = sqrtf(ctx.nx * ctx.nx + ctx.nz * ctx.nz);
    ctx.alpha = atan2f(ctx.nx, ctx.nz);
    ctx.centerPhi = acosf(ctx.ny);
    ctx.u = u;
    ctx.w = w;
    ctx.startAngle = startAngle;
    ctx.endAngle = endAngle;
    ctx.checkSector = (std::abs(endAngle - startAngle) < (2 * PI_F - 0.001f));
    ctx.clipPlanes = clipPlanes;

    // Calculate Vertical Bounds
    float a1 = ctx.centerPhi - ctx.targetAngle;
    float a2 = ctx.centerPhi + ctx.targetAngle;
    float p1 = acosf(cosf(a1));
    float p2 = acosf(cosf(a2));

    float minP = std::min(p1, p2);
    float maxP = std::max(p1, p2);

    float phiMin = std::max(0.0f, minP - thickness);
    float phiMax = std::min(PI_F, maxP + thickness);

    if (minPhiLimit >= 0) phiMin = std::max(phiMin, minPhiLimit);
    if (maxPhiLimit >= 0) phiMax = std::min(phiMax, maxPhiLimit);

    if (phiMin > phiMax) return;

    int yMin = std::max(0, static_cast<int>(floorf((phiMin * (H_VIRT - 1)) / PI_F)));
    int yMax = std::min(H - 1, static_cast<int>(ceilf((phiMax * (H_VIRT - 1)) / PI_F)));

    for (int y = yMin; y <= yMax; y++) {
      scanRow(canvas, y, ctx);
    }
  }

private:
  static void scanRow(Canvas& canvas, int y, const FSRingContext& ctx) {
    float phi = y_to_phi(y);
    float cosPhi = cosf(phi);
    float sinPhi = sinf(phi);

    // Singularity (Poles or Vertical Normal)
    if (ctx.R < 0.01f) {
      scanFullRow(canvas, y, ctx);
      return;
    }

    // General Intersection
    float ang_low = std::max(0.0f, ctx.targetAngle - ctx.thickness);
    float ang_high = std::min(PI_F, ctx.targetAngle + ctx.thickness);

    float D_max = cosf(ang_low);
    float D_min = cosf(ang_high);

    float denom = ctx.R * sinPhi;
    if (std::abs(denom) < 0.000001f) {
      scanFullRow(canvas, y, ctx);
      return;
    }

    float C_min = (D_min - ctx.ny * cosPhi) / denom;
    float C_max = (D_max - ctx.ny * cosPhi) / denom;

    float minCos = std::max(-1.0f, C_min);
    float maxCos = std::min(1.0f, C_max);
    if (minCos > maxCos) return;

    float angleMin = acosf(maxCos);
    float angleMax = acosf(minCos);

    // Generate scan windows
    if (angleMin <= 0.0001f) {
      scanWindow(canvas, y, ctx.alpha - angleMax, ctx.alpha + angleMax, ctx);
    }
    else if (angleMax >= PI_F - 0.0001f) {
      scanWindow(canvas, y, ctx.alpha + angleMin, ctx.alpha + 2 * PI_F - angleMin, ctx);
    }
    else {
      scanWindow(canvas, y, ctx.alpha - angleMax, ctx.alpha - angleMin, ctx);
      scanWindow(canvas, y, ctx.alpha + angleMin, ctx.alpha + angleMax, ctx);
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
    Vector p = pixel_to_vector<W>(x, y);

    // Clipping Planes
    if (ctx.clipPlanes) {
      for (const auto& cp : *ctx.clipPlanes) {
        if (dot(p, cp) < 0) return;
      }
    }

    float polarAngle = angle_between(p, ctx.normal);
    float dist = std::abs(polarAngle - ctx.targetAngle);

    if (dist < ctx.thickness) {
      // Sector Check
      if (ctx.checkSector) {
        float dotU = dot(p, ctx.u);
        float dotW = dot(p, ctx.w);
        float azimuth = atan2f(dotW, dotU);
        if (azimuth < 0) azimuth += 2 * PI_F;

        bool inside = false;
        if (ctx.startAngle <= ctx.endAngle) {
          inside = (azimuth >= ctx.startAngle && azimuth <= ctx.endAngle);
        }
        else {
          inside = (azimuth >= ctx.startAngle || azimuth <= ctx.endAngle);
        }
        if (!inside) return;
      }

      // Render
      float t = dist / ctx.thickness;
      float alphaFactor = quintic_kernel(1.0f - t);
      float finalAlpha = ctx.color.alpha * alphaFactor;

      Pixel& outColor = canvas(XY(x, y));
      outColor = blend_alpha(finalAlpha)(outColor, ctx.color.color);
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

    float maxY = std::max(v1.j, v2.j);
    float minY = std::min(v1.j, v2.j);

    Vector apexPlaneNormal = cross(normal, Y_AXIS);
    if (dot(apexPlaneNormal, apexPlaneNormal) > 0.0001f) {
      float d1 = dot(v1, apexPlaneNormal);
      float d2 = dot(v2, apexPlaneNormal);

      // Segment crosses the apex plane
      if (d1 * d2 <= 0) {
        float globalMaxY = sqrtf(1.0f - normal.j * normal.j);
        if (v1.j + v2.j > 0) {
          maxY = globalMaxY;
        }
        else {
          minY = -globalMaxY;
        }
      }
    }

    float clampedMaxY = std::clamp(maxY, -1.0f, 1.0f);
    float clampedMinY = std::clamp(minY, -1.0f, 1.0f);

    float minPhi = acosf(clampedMaxY) - thickness;
    float maxPhi = acosf(clampedMinY) + thickness;

    FSRing<W>::draw(canvas, normal, 1.0f, color, thickness, 0, 2 * PI_F, &clips, minPhi, maxPhi);
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

  void drawRing(Canvas& canvas, const Vector& normal, float radius, const Color4& color, float thickness, float startAngle = 0, float endAngle = 2 * PI_F) {
    FSRing<W>::draw(canvas, normal, radius, color, thickness, startAngle, endAngle);
  }
};
