/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#include <vector>
#include <algorithm>
#include <cmath>
#include "geometry.h"
#include "color.h"
#include "led.h"
#include "filter.h"

// Forward declare quintic_kernel if not available (it IS in filter.h, but we might need it visible here)
// However, scan.h includes filter.h, so it should be fine.

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
      float alpha_angle; /**< The azimuth angle of the normal vector in the XZ plane. */

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

      /**
       * @brief Computes the horizontal scanline intervals for this shape at a given y-coordinate.
       * @param y The vertical pixel coordinate.
       * @param out Output iterator or callback accepting (float start, float end).
       * @return True if intervals were found and reported.
       */
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

       /**
        * @brief Computes Horizontal intervals.
        * @details Not implemented for DistortedRing (returns false to force full scan).
        */
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
        
        if constexpr (std::is_same_v<std::decay_t<decltype(c)>, Pixel>) {
            pipeline.plot(canvas, x, y, c.color, 0.0f, alpha);
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
