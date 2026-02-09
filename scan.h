/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include "geometry.h"
#include "color.h"
#include "constants.h"
#include "filter.h"

namespace SDF {
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
      float aux; // Auxiliary value (e.g. barycentric coordinate)
      
      struct Weights {
          float a = 0, b = 0, c = 0;
          int i0 = 0, i1 = 0, i2 = 0;
          bool valid = false;
      } weights;
    };

    /**
     * @brief Calculates signed distance to a ring.
     * Returns:
     *  dist: Signed distance (negative inside)
     *  t: Normalized parameter (0-1) corresponding to angle/2PI
     *  raw_dist: Unsigned distance to centerline
     */
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
      
      /**
       * @brief Computes signed distance to the ring.
       * @param p Point on sphere (normalized).
       * @return DistanceResult {dist, t, raw_dist}.
       */
      DistanceResult distance(const Vector& p) const {
         float d = dot(p, normal);
         if (d < cos_min || d > cos_max) return { 100.0f, 0.0f, 100.0f, 0.0f };
         
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
         
         return { dist - thickness, t, dist, 0.0f };
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
       
       // Optimization
       float r_val;
       float alpha_angle;
       float cos_max_limit, cos_min_limit;

       DistortedRing(const Basis& b, float r, float th, std::function<float(float)> sf, float md, float ph)
         : basis(b), radius(r), thickness(th), shift_fn(sf), max_distortion(md), phase(ph)
       {
          normal = basis.v; u = basis.u; w = basis.w;
          nx = normal.i; ny = normal.j; nz = normal.k;
          target_angle = radius * (PI_F / 2.0f);
          center_phi = acosf(ny);
          max_thickness = thickness + max_distortion;
          
          r_val = sqrtf(nx*nx + nz*nz);
          alpha_angle = atan2f(nx, nz);
          
          float ang_min = std::max(0.0f, target_angle - max_thickness);
          float ang_max = std::min(PI_F, target_angle + max_thickness);
          cos_max_limit = cosf(ang_min);
          cos_min_limit = cosf(ang_max);
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

          float margin = max_thickness + 0.1f;
          float f_phi_min = std::max(0.0f, phi_min - margin);
          float f_phi_max = std::min(PI_F, phi_max + margin);

          int y_min = std::max(0, static_cast<int>(floorf((f_phi_min * (H_VIRT - 1)) / PI_F)));
          int y_max = std::min(H - 1, static_cast<int>(ceilf((f_phi_max * (H_VIRT - 1)) / PI_F)));
           return { y_min, y_max };
       }

       template<typename OutputIt>
       bool get_horizontal_intervals(int y, OutputIt out) const {
         float phi = y_to_phi(static_cast<float>(y));
         float cos_phi = cosf(phi);
         float sin_phi = sinf(phi);

         if (r_val < 0.01f) return false; 
         
         float denom = r_val * sin_phi;
         if (std::abs(denom) < 0.000001f) return false; 

         float C_min = (cos_min_limit - ny * cos_phi) / denom;
         float C_max = (cos_max_limit - ny * cos_phi) / denom;
         float min_cos = std::max(-1.0f, C_min);
         float max_cos = std::min(1.0f, C_max);

         if (min_cos > max_cos) return true; // Empty row

         float angle_min = acosf(max_cos);
         float angle_max = acosf(min_cos);

         out(alpha_angle - angle_max, alpha_angle - angle_min);
         out(alpha_angle + angle_min, alpha_angle + angle_max);
         return true;
       }

      /**
       * @brief Computes signed distance to the distorted ring.
       * @param p Point on sphere (normalized).
       * @return DistanceResult {dist, t, raw_dist}.
       */
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
         
         return { dist - thickness, azimuth / (2 * PI_F), dist, 0.0f };
      }
    };

    template <typename A, typename B>
    struct Union {
       const A& a;
       const B& b;
       float thickness;

       Union(const A& shapeA, const B& shapeB)
          : a(shapeA), b(shapeB), thickness(std::max(shapeA.thickness, shapeB.thickness)) {}

       Bounds get_vertical_bounds() const {
          auto b1 = a.get_vertical_bounds();
          auto b2 = b.get_vertical_bounds();
          return { std::min(b1.y_min, b2.y_min), std::max(b1.y_max, b2.y_max) };
       }

       template<typename OutputIt>
       bool get_horizontal_intervals(int y, OutputIt out) const {
           // Union intervals are complex; fallback to full scan.
           return false; 
       }

       /**
        * @brief Signed distance to Union.
        */
       DistanceResult distance(const Vector& p) const {
          auto resA = a.distance(p);
          auto resB = b.distance(p);
          if (resA.dist < resB.dist) return resA; // Min distance
          return resB;
       }
    };

    template <typename A, typename B>
    struct Subtract {
       const A& a;
       const B& b;
       float thickness;

       Subtract(const A& shapeA, const B& shapeB)
          : a(shapeA), b(shapeB), thickness(shapeA.thickness) {}

       Bounds get_vertical_bounds() const {
          return a.get_vertical_bounds(); // Subtraction is bounded by A
       }

       template<typename OutputIt>
       bool get_horizontal_intervals(int y, OutputIt out) const {
           // Subtraction intervals delegate to A (conservative)
           return a.get_horizontal_intervals(y, out);
       }

       /**
        * @brief Signed distance to Subtraction (A - B).
        */
       DistanceResult distance(const Vector& p) const {
          auto resA = a.distance(p);
          auto resB = b.distance(p);
          // Max(A, -B)
          if (-resB.dist > resA.dist) {
             resB.dist = -resB.dist;
             return resB;
          }
          return resA;
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
           StaticCircularBuffer<std::pair<float, float>, 32> intervalsA;
           StaticCircularBuffer<std::pair<float, float>, 32> intervalsB;

           bool hasA = a.get_horizontal_intervals(y, [&](float start, float end) {
               intervalsA.push_back({start, end});
           });
           
           bool hasB = b.get_horizontal_intervals(y, [&](float start, float end) {
               intervalsB.push_back({start, end});
           });

           if (!hasA) {
               return b.get_horizontal_intervals(y, out);
           }
           if (!hasB) {
               return a.get_horizontal_intervals(y, out);
           }

           if (intervalsA.is_empty() || intervalsB.is_empty()) return true;

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
       DistanceResult distance(const Vector& p) const {
          auto resA = a.distance(p);
          auto resB = b.distance(p);
          // Max(A, B)
          if (resA.dist > resB.dist) return resA;
          return resB;
       }
    };

    struct FaceScratchBuffer {
        static constexpr int MAX_VERTS = 16;
        std::array<Vector, MAX_VERTS> poly2D;
        std::array<Vector, MAX_VERTS> edgeVectors;
        std::array<float, MAX_VERTS> edgeLengthsSq;
        std::array<Vector, MAX_VERTS> planes;
        std::array<std::pair<float, float>, 4> intervals; 
        std::array<float, MAX_VERTS> thetas; 
    };

    struct Face {
        Vector center;
        Vector basisV, basisU, basisW;
        int count;
        float thickness;
        
        std::span<Vector> poly2D; 
        std::span<Vector> edgeVectors; 
        std::span<float> edgeLengthsSq;
        std::span<Vector> planes; 
        
        int y_min, y_max;
        std::span<std::pair<float, float>> intervals; 
        bool full_width;
        
        Face(std::span<const Vector> vertices, std::span<const int> indices, float th, FaceScratchBuffer& scratch) 
          : thickness(th), full_width(true)
        {
           count = indices.size();
           if (count > FaceScratchBuffer::MAX_VERTS) count = FaceScratchBuffer::MAX_VERTS;

           center = Vector(0,0,0);
           for (size_t i=0; i<count; ++i) center = center + vertices[indices[i]];
           center.normalize();
           
           basisV = center;
           Vector ref = (std::abs(dot(center, X_AXIS)) > 0.9f) ? Y_AXIS : X_AXIS;
           basisU = cross(center, ref).normalize();
           basisW = cross(center, basisU).normalize();
           
           // Project 2D
           for(int i=0; i<count; ++i) {
              const Vector& v = vertices[indices[i]];
              float d = dot(v, basisV);
              scratch.poly2D[i].i = dot(v, basisU) / d; // x
              scratch.poly2D[i].j = dot(v, basisW) / d; // y
              scratch.poly2D[i].k = 0;
           }
           poly2D = std::span<Vector>(scratch.poly2D.data(), count);
           
           // Edges and Planes
           int planes_count = 0;
           
           float min_phi = 100.0f;
           float max_phi = -100.0f;
           
           for(int i=0; i<count; ++i) {
              int idx1 = indices[i];
              int idx2 = indices[(i+1)%count];
              const Vector& v1 = vertices[idx1];
              const Vector& v2 = vertices[idx2];
              
              // Edge 2D
              Vector edge = scratch.poly2D[(i+1)%count] - scratch.poly2D[i];
              scratch.edgeVectors[i] = edge;
              scratch.edgeLengthsSq[i] = dot(edge, edge);
              
              // Plane Normal
              Vector normal = cross(v1, v2);
              float lenSq = dot(normal, normal);
              if (lenSq > 1e-12f) {
                  scratch.planes[planes_count++] = normal.normalize();
              }
              
               // Vertex Bounds
               float phi_val = acosf(std::clamp(v1.j, -1.0f, 1.0f)); 
              if (phi_val < min_phi) min_phi = phi_val;
              if (phi_val > max_phi) max_phi = phi_val;
              
               // Arc Extrema Logic
               if (planes_count > 0) {
                  const Vector& n = scratch.planes[planes_count - 1];
                  float ny = n.j; 
                  if (std::abs(ny) < 0.99999f) {
                      float nx = n.i;
                      float nz = n.k;
                      float tx = -nx * ny;
                      float ty = 1.0f - ny * ny;
                      float tz = -nz * ny;
                      float tLenSq = tx*tx + ty*ty + tz*tz;
                      if (tLenSq > 1e-12f) {
                          float invLen = 1.0f / sqrtf(tLenSq);
                          float ptx = tx * invLen;
                          float pty = ty * invLen;
                          float ptz = tz * invLen;
                          
                          float cx1 = (v1.j * ptz - v1.k * pty) * nx + (v1.k * ptx - v1.i * ptz) * ny + (v1.i * pty - v1.j * ptx) * nz;
                          float cx2 = (pty * v2.k - ptz * v2.j) * nx + (ptz * v2.i - ptx * v2.k) * ny + (ptx * v2.j - pty * v2.i) * nz;
                          
                          if (cx1 > 0 && cx2 > 0) {
                              float phiTop = acosf(std::clamp(pty, -1.0f, 1.0f));
                              if (phiTop < min_phi) min_phi = phiTop;
                          }
                          if (cx1 < 0 && cx2 < 0) {
                              float phiBot = acosf(std::clamp(-pty, -1.0f, 1.0f));
                              if (phiBot > max_phi) max_phi = phiBot;
                          }
                      }
                  }
               }

              // Thetas
              float theta = atan2f(v1.i, v1.k); 
              if (theta < 0) theta += 2 * PI_F;
              scratch.thetas[i] = theta;
           }
           
           edgeVectors = std::span<Vector>(scratch.edgeVectors.data(), count);
           edgeLengthsSq = std::span<float>(scratch.edgeLengthsSq.data(), count);
           planes = std::span<Vector>(scratch.planes.data(), planes_count);
           
           // Pole Logic
           bool npInside = true;
           bool spInside = true;
           for(const auto& p : planes) {
               if (p.j < 0) npInside = false;
               if (p.j > 0) spInside = false;
           }
           if (npInside) min_phi = 0.0f;
           if (spInside) max_phi = PI_F;
           
           // Vertical Bounds
           float margin = thickness + 0.05f;
           y_min = std::max(0, static_cast<int>(floorf((std::max(0.0f, min_phi - margin) * (H_VIRT - 1)) / PI_F)));
           y_max = std::min(H - 1, static_cast<int>(ceilf((std::min(PI_F, max_phi + margin) * (H_VIRT - 1)) / PI_F)));
           
           // Horizontal Interval Logic
           std::sort(scratch.thetas.begin(), scratch.thetas.begin() + count);
           float maxGap = 0;
           float gapStart = 0;
           for(size_t i=0; i<count; ++i) {
               float next = (i+1 < count) ? scratch.thetas[i+1] : (scratch.thetas[0] + 2*PI_F);
               float diff = next - scratch.thetas[i];
               if (diff > maxGap) {
                   maxGap = diff;
                   gapStart = scratch.thetas[i];
               }
           }
           
           int interval_count = 0;
           if (maxGap > PI_F) {
               full_width = false;
               float startT = fmodf(gapStart + maxGap, 2*PI_F);
               float endT = gapStart;
               
               if (startT <= endT) {
                   scratch.intervals[interval_count++] = {startT, endT};
               } else {
                   // Split wrapped interval
                   scratch.intervals[interval_count++] = {startT, 2*PI_F};
                   scratch.intervals[interval_count++] = {0.0f, endT};
               }
           } else {
               full_width = true;
           }
           intervals = std::span<std::pair<float, float>>(scratch.intervals.data(), interval_count);
        }
        
        Bounds get_vertical_bounds() const { return { y_min, y_max }; }
        
        template<typename OutputIt>
        bool get_horizontal_intervals(int y, OutputIt out) const {
           if (full_width) return false;
           for(const auto& iv : intervals) {
               out(iv.first, iv.second);
           }
           return true; 
        }

        /**
         * @brief Signed distance to planar Face.
         */
        DistanceResult distance(const Vector& p) const {
            float cosAngle = dot(p, center);
            if (cosAngle <= 0.01f) return { 100.0f, 0.0f, 100.0f, 0.0f };

            float invCos = 1.0f / cosAngle;
            float px = dot(p, basisU) * invCos;
            float py = dot(p, basisW) * invCos;
            
            float d = FLT_MAX;
            int winding = 0;
            
            for(int i=0; i<count; ++i) {
               const Vector& V_curr = poly2D[i];
               const Vector& V_next = poly2D[(i+1)%count];
               const Vector& edge = edgeVectors[i]; 
               
               float ex = edge.i;
               float ey = edge.j;
               
               float wx = px - V_curr.i;
               float wy = py - V_curr.j;
               
               float dotWE = wx * ex + wy * ey;
               float dotEE = edgeLengthsSq[i];
               
               float clampVal = 0;
               if (dotEE > 1e-12f) {
                   clampVal = std::max(0.0f, std::min(1.0f, dotWE / dotEE));
               }
               
               float bx = wx - ex * clampVal;
               float by = wy - ey * clampVal;
               float distSq = bx*bx + by*by;
               
               if (distSq < d) d = distSq;
               
               bool isUpward = (V_curr.j <= py) && (V_next.j > py);
               bool isDownward = (V_curr.j > py) && (V_next.j <= py);
               
               if (isUpward || isDownward) {
                   float crossVal = ex * wy - ey * wx;
                   if (isUpward) {
                       if (crossVal > 0) winding++;
                   } else {
                       if (crossVal < 0) winding--;
                   }
               }
            }
            
            float s = (winding != 0) ? -1.0f : 1.0f;
            float planeDist = s * sqrtf(d);
            
            DistanceResult res = { planeDist - thickness, 0.0f, planeDist, 0.0f };
            
            // Barycentric Weights
            if (count == 3) {
                 const Vector& v0 = poly2D[0];
                 const Vector& v1 = poly2D[1];
                 const Vector& v2 = poly2D[2];
                 float denom = (v1.j - v2.j) * (v0.i - v2.i) + (v2.i - v1.i) * (v0.j - v2.j);
                 if (std::abs(denom) > 1e-12f) {
                     float invDenom = 1.0f / denom;
                     res.weights.a = ((v1.j - v2.j) * (px - v2.i) + (v2.i - v1.i) * (py - v2.j)) * invDenom;
                     res.weights.b = ((v2.j - v0.j) * (px - v2.i) + (v0.i - v2.i) * (py - v2.j)) * invDenom;
                     res.weights.c = 1.0f - res.weights.a - res.weights.b;
                     res.weights.i0 = 0;
                     res.weights.i1 = 1;
                     res.weights.i2 = 2;
                     res.weights.valid = true;
                 }
            } else {
                 // Triangle Fan
                 for (int i = 1; i < count - 1; i++) {
                     const Vector& v0 = poly2D[0];
                     const Vector& v1 = poly2D[i];
                     const Vector& v2 = poly2D[i + 1];
                     
                     float denom = (v1.j - v2.j) * (v0.i - v2.i) + (v2.i - v1.i) * (v0.j - v2.j);
                     if (std::abs(denom) > 1e-12f) {
                         float invDenom = 1.0f / denom;
                         float wA = ((v1.j - v2.j) * (px - v2.i) + (v2.i - v1.i) * (py - v2.j)) * invDenom;
                         float wB = ((v2.j - v0.j) * (px - v2.i) + (v0.i - v2.i) * (py - v2.j)) * invDenom;
                         float wC = 1.0f - wA - wB;
                         
                         if (wA >= -0.01f && wB >= -0.01f && wC >= -0.01f) {
                             res.weights.a = wA;
                             res.weights.b = wB;
                             res.weights.c = wC;
                             res.weights.i0 = 0;
                             res.weights.i1 = i;
                             res.weights.i2 = i + 1;
                             res.weights.valid = true;
                             break;
                         }
                     }
                 }
                 if (!res.weights.valid) {
                      res.weights.a = 1; res.weights.b = 0; res.weights.c = 0;
                      res.weights.i0 = 0; res.weights.i1 = 1; res.weights.i2 = 2;
                      res.weights.valid = true;
                 }
            }
            
            if (res.weights.valid) {
                res.t = res.weights.a;
                res.aux = res.weights.b;
            }
            return res;
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
       
       /**
        * @brief Signed distance to Polygon edge.
        */
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
          return { dist_edge, polar / thickness, polar, 0.0f };
       }
    };

    /**
     * @brief Calculates signed distance to a star shape.
     * Returns:
     *  dist: Signed distance from edge (negative inside)
     *  t: Normalized polar distance (polar / thickness)
     *  raw_dist: Polar distance from center
     */
    template <int W>
    struct Star {
        const Basis& basis;
        float radius;
        int sides;
        float phase;
        
        float nx, ny, plane_d; // Edge plane equation (2D)
        float thickness;
        
        // Scan optimization
        float scan_nx, scan_ny, scan_nz;
        float scan_r_val, scan_alpha_angle;
        int y_min, y_max;

        Star(const Basis& b, float r, int s, float ph)
            : basis(b), radius(r), sides(s), phase(ph)
        {
            float outer_radius = radius * (PI_F / 2.0f);
            float inner_radius = outer_radius * 0.382f;
            float angle_step = PI_F / sides;

            float v_t = outer_radius;
            float v_vx = inner_radius * cosf(angle_step);
            float v_vy = inner_radius * sinf(angle_step);

            // Vector from Tip to Valley
            float dx = v_vx - 0.0f; // Tip x is 0
            float dy = v_vy - v_t;  // Tip y is v_t
            float len = sqrtf(dx * dx + dy * dy);
            
            // Normal to edge
            nx = -dy / len;
            ny = dx / len;
            plane_d = -(nx * 0.0f + ny * v_t);
            thickness = outer_radius;

            // Scan Bounds setup
            scan_nx = basis.v.i;
            scan_ny = basis.v.j;
            scan_nz = basis.v.k;
            scan_r_val = sqrtf(scan_nx * scan_nx + scan_nz * scan_nz);
            scan_alpha_angle = atan2f(scan_nx, scan_nz);

            float center_phi = acosf(std::max(-1.0f, std::min(1.0f, basis.v.j)));
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

            if (scan_r_val < 0.01f) return false;

            float ang_high = thickness + (2.0f * PI_F / W); // Margin for AA
            float D_min = cosf(ang_high);
            float denom = scan_r_val * sin_phi;

            if (std::abs(denom) < 0.000001f) return false;

            float C_min = (D_min - scan_ny * cos_phi) / denom;

            if (C_min > 1.0f) return true;
            if (C_min < -1.0f) return false;

            float d_alpha = acosf(C_min);
            out(scan_alpha_angle - d_alpha, scan_alpha_angle + d_alpha);
            return true;
        }

        DistanceResult distance(const Vector& p) const {
            float scan_dist = angle_between(p, basis.v);
            float dot_u = dot(p, basis.u);
            float dot_w = dot(p, basis.w);
            float azimuth = atan2f(dot_w, dot_u);
            if (azimuth < 0) azimuth += 2 * PI_F;

            azimuth += phase;

            float sector_angle = 2 * PI_F / sides;
            float local_azimuth = wrap(azimuth + sector_angle / 2.0f, sector_angle) - sector_angle / 2.0f;
            local_azimuth = std::abs(local_azimuth);

            float px = scan_dist * cosf(local_azimuth);
            float py = scan_dist * sinf(local_azimuth);
            
            float dist_to_edge = px * nx + py * ny + plane_d;
            
            return { -dist_to_edge, azimuth / (2 * PI_F), scan_dist, 0.0f };
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
       const Basis& basis;
       int sides;
       float phase;
       float thickness;
       float apothem;
       Vector antipode;
       // Optimization
       float scan_nx, scan_ny, scan_nz, scan_R, scan_alpha;
       int y_min, y_max;
       
       Flower(const Basis& b, float radius, int s, float ph)
         : basis(b), sides(s), phase(ph)
       {
          float outer = radius * (PI_F / 2.0f);
          apothem = PI_F - outer;
          thickness = outer;
          antipode = -basis.v;
          
          scan_nx = antipode.i; scan_ny = antipode.j; scan_nz = antipode.k;
          scan_R = sqrtf(scan_nx*scan_nx + scan_nz*scan_nz);
          scan_alpha = atan2f(scan_nx, scan_nz);
          
          float center_phi = acosf(std::max(-1.0f, std::min(1.0f, antipode.j)));
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

         if (scan_R < 0.01f) return false; 
         
         float denom = scan_R * sin_phi;
         if (std::abs(denom) < 0.000001f) return false; 

         float ang_max = thickness; 
         float cos_limit = cosf(ang_max); 
         
         float C_min = (cos_limit - scan_ny * cos_phi) / denom;
         float C_max = (1.0f - scan_ny * cos_phi) / denom;
         
         float min_cos = std::max(-1.0f, C_min);
         float max_cos = std::min(1.0f, C_max);

         if (min_cos > max_cos) return true; 
         
         float angle_min = acosf(max_cos);
         float angle_max = acosf(min_cos);

         out(scan_alpha - angle_max, scan_alpha - angle_min);
         out(scan_alpha + angle_min, scan_alpha + angle_max);
         return true;
       }
       
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
          return { -dist_edge, scan_dist / thickness, scan_dist, 0.0f };
       }
    };

    struct HarmonicBlob {
        int l;
        int m;
        float amplitude;
        Quaternion inv_q;
        std::function<float(int, int, float, float)> harmonic_fn;
        bool is_solid = true;

        HarmonicBlob(int l, int m, float amplitude, const Quaternion& orientation, std::function<float(int, int, float, float)> harmonic_fn)
            : l(l), m(m), amplitude(amplitude), harmonic_fn(harmonic_fn)
        {
            inv_q = orientation.inverse();
        }

        Bounds get_vertical_bounds() const { return { 0, H - 1 }; }

        DistanceResult distance(const Vector& p) const {
             // Transform world pixel position to local harmonic space
             Vector v = rotate(p, inv_q);
             float phi = acosf(std::max(-1.0f, std::min(1.0f, v.j)));
             float theta = atan2f(v.k, v.i);
             float harmonic_val = harmonic_fn(l, m, theta, phi);
 
             float lobe_radius = 1.0f + std::abs(harmonic_val) * amplitude;
             float d = 1.0f - lobe_radius;
 
             return { d, tanhf(std::abs(harmonic_val) * amplitude), harmonic_val, 0.0f };
        }
        
        template<typename OutputIt>
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
        
        Line(const Vector& start, const Vector& end, float th)
           : a(start), b(end), thickness(th)
        {
            len = angle_between(a, b);
            if (len < 1e-6f) {
                n = Vector(0,0,0);
            } else {
                n = cross(a, b).normalize();
            }
        }
        
        Bounds get_vertical_bounds() const { return { 0, H - 1 }; }
        
        DistanceResult distance(const Vector& p) const {
             if (len < 1e-6f) {
                 float dist = angle_between(p, a);
                 return { dist - thickness, 0.0f, dist, 0.0f };
             }
             
             float d_plane = dot(p, n);
             Vector p_proj_plane = p - n * d_plane;
             float proj_mag = p_proj_plane.magnitude();

             if (proj_mag < 1e-6f) {
                 float dA = angle_between(p, a);
                 float dB = angle_between(p, b);
                 float dist = std::min(dA, dB);
                 return { dist - thickness, 0.0f, dist, 0.0f };
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
             
             return { dist_seg - thickness, 0.0f, dist_seg, 0.0f };
        }
        
        template<typename OutputIt>
        bool get_horizontal_intervals(int y, OutputIt out) const {
            return false;
        }
    };
}

/**
 * @brief The Scan struct contains volumetric (raster) drawing primitives.
 * 
 * General Register Mapping for Scan Primitives:
 *  v0: Normalized parameter t (0-1) or angle
 *  v1: Raw Distance or Supplementary Value
 */
namespace Scan {
  
  template <int W>
  static void rasterize(auto& pipeline, Canvas& canvas, const auto& shape, FragmentShaderFn auto fragment_shader, bool debug_bb = false) {
     auto bounds = shape.get_vertical_bounds();
     
     for (int y = bounds.y_min; y <= bounds.y_max; ++y) {
        bool handled = shape.get_horizontal_intervals(y, [&](float t1, float t2) {
              int x1 = static_cast<int>(floorf((t1 * W) / (2 * PI_F)));
           int x2 = static_cast<int>(ceilf((t2 * W) / (2 * PI_F)));
           
           for (int x = x1; x <= x2; ++x) {
              int wx = wrap(x, W);
              process_pixel<W>(wx, y, pipeline, canvas, shape, fragment_shader, debug_bb);
           }
        });
        
        if (!handled) {
           for (int x = 0; x < W; ++x) {
              process_pixel<W>(x, y, pipeline, canvas, shape, fragment_shader, debug_bb);
           }
        }
     }
  }
  
  template <int W>
  static void process_pixel(int x, int y, auto& pipeline, Canvas& canvas, const auto& shape, FragmentShaderFn auto fragment_shader, bool debug_bb) {
     const Vector& p = pixel_to_vector<W>(x, y);
     
     auto result = shape.distance(p);
     float d = result.dist;
     
     float pixel_width = 2.0f * PI_F / W;
     float threshold = pixel_width;
     
     if (d < threshold) {
        float t_aa = 0.5f - d / (2.0f * pixel_width);
        float alpha = quintic_kernel(std::max(0.0f, std::min(1.0f, t_aa)));
        
        if (alpha <= 0.001f) return;
        
        Fragment frag;
        frag.pos = p;
        frag.v0 = result.t;
        frag.v1 = result.raw_dist;
        frag.v2 = 0.0f;
        frag.v3 = result.aux;
        frag.age = 0;
        
        Fragment f_out = fragment_shader(p, frag);
        
        if (f_out.color.alpha > 0.001f) {
            pipeline.plot(canvas, x, y, f_out.color.color, f_out.age, f_out.color.alpha * alpha, f_out.blend);
        }
     }
  }

  struct DistortedRing {
    template <int W>
    static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, float thickness, 
                     std::function<float(float)> shift_fn, float amplitude, FragmentShaderFn auto fragment_shader, float phase = 0, bool debug_bb = false)  
    {
      SDF::DistortedRing shape(basis, radius, thickness, shift_fn, amplitude, phase);
      Scan::rasterize<W>(pipeline, canvas, shape, fragment_shader, debug_bb);
    }
  };

   struct PlanarPolygon {
      template <int W>
      static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, int sides, FragmentShaderFn auto fragment_shader, float phase = 0, bool debug_bb = false) {
        auto res = get_antipode(basis, radius);
        float thickness = res.second * (PI_F / 2.0f);
        
        SDF::Polygon shape(res.first, res.second, thickness, sides, phase);
        Scan::rasterize<W>(pipeline, canvas, shape, fragment_shader, debug_bb);
      }
   };

   struct Line {
      template <int W>
      static void draw(auto& pipeline, Canvas& canvas, const Vector& v1, const Vector& v2, float thickness, FragmentShaderFn auto fragment_shader, bool debug_bb = false) {
          SDF::Line shape(v1, v2, thickness);
          Scan::rasterize<W>(pipeline, canvas, shape, fragment_shader, debug_bb);
      }
   };



  struct Ring {
    template <int W>
    static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, float thickness, FragmentShaderFn auto fragment_shader, float phase = 0, bool debug_bb = false) {
        auto res = get_antipode(basis, radius);
        SDF::Ring shape(res.first, res.second, thickness, phase);
        Scan::rasterize<W>(pipeline, canvas, shape, fragment_shader, debug_bb);
     }

      // Overload for Vector normal inputs
     template <int W>
     static void draw(auto& pipeline, Canvas& canvas, const Vector& normal, float radius, float thickness, FragmentShaderFn auto fragment_shader, float phase = 0, bool debug_bb = false) {
        Basis basis = make_basis(Quaternion(), normal);
        draw<W>(pipeline, canvas, basis, radius, thickness, fragment_shader, phase, debug_bb);
     }
  };

   struct Circle {
      template <int W>
      static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, FragmentShaderFn auto fragment_shader, bool debug_bb = false) {
          // Circle is a Ring with inner radius 0
          float th = radius * (PI_F / 2.0f);
          Ring::draw<W>(pipeline, canvas, basis, 0.0f, th, fragment_shader, 0, debug_bb);
      }

      template <int W>
      static void draw(auto& pipeline, Canvas& canvas, const Vector& normal, float radius, FragmentShaderFn auto fragment_shader, bool debug_bb = false) {
          Basis basis = make_basis(Quaternion(), normal);
          draw<W>(pipeline, canvas, basis, radius, fragment_shader, debug_bb);
      }
   };
  
   struct Point {
      template <int W>
      static void draw(auto& pipeline, Canvas& canvas, const Vector& p, float thickness, FragmentShaderFn auto fragment_shader) {
         // Scan.Point uses Scan.Ring with radius 0
         Basis basis = make_basis(Quaternion(), p);
         Ring::draw<W>(pipeline, canvas, basis, 0.0f, thickness, fragment_shader);
      }
   };
   
    struct Star {
       template <int W>
       static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, int sides, FragmentShaderFn auto fragment_shader, float phase = 0, bool debug_bb = false) {
           auto res = get_antipode(basis, radius);
           SDF::Star<W> shape(res.first, res.second, sides, phase);
           Scan::rasterize<W>(pipeline, canvas, shape, fragment_shader, debug_bb);
       }
    };
    
    struct Flower {
       template <int W>
       static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, int sides, FragmentShaderFn auto fragment_shader, float phase = 0, bool debug_bb = false) {
           auto res = get_antipode(basis, radius);
           SDF::Flower shape(res.first, res.second, sides, phase);
           Scan::rasterize<W>(pipeline, canvas, shape, fragment_shader, debug_bb);
       }
    };
   
   struct SphericalPolygon {
        template <int W>
        static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, int sides, FragmentShaderFn auto fragment_shader, float phase = 0, bool debug_bb = false) {
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
             
             std::vector<Vector> vertices;
             vertices.reserve(sides);
             
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
                 vertices.push_back(pos);
             }
             
             // 3. Create Indices (0, 1, ..., sides-1)
             std::vector<int> indices(sides);
             std::iota(indices.begin(), indices.end(), 0);
             
             // 4. Create SDF::Face
             SDF::FaceScratchBuffer scratch;
             SDF::Face shape(vertices, indices, 0.0f, scratch);
             
             Scan::rasterize<W>(pipeline, canvas, shape, fragment_shader, debug_bb);
        }
   };

   struct HarmonicBlob {
       template <int W>
       static void draw(auto& pipeline, Canvas& canvas, int l, int m, float amplitude, const Quaternion& orientation, 
                        std::function<float(int, int, float, float)> harmonic_fn, FragmentShaderFn auto fragment_shader, bool debug_bb = false)
       {
           SDF::HarmonicBlob shape(l, m, amplitude, orientation, harmonic_fn);
           Scan::rasterize<W>(pipeline, canvas, shape, fragment_shader, debug_bb);
       }
   };

  struct Mesh {
      template <int W, typename MeshT, typename F>
      static void draw(auto& pipeline, Canvas& canvas, const MeshT& mesh, F fragment_shader, bool debug_bb = false) {
          SDF::FaceScratchBuffer scratch;
          
          size_t idx_offset = 0;
          for (size_t i = 0; i < mesh.num_faces; ++i) {
             size_t count = mesh.face_counts[i];
             std::span<const Vector> verts(mesh.vertices.data(), mesh.num_vertices);
             std::span<const int> indices(&mesh.faces[idx_offset], count);
             
             SDF::Face shape(verts, indices, 0.0f, scratch);
             idx_offset += count;
             
             auto wrapper = [&](const Vector& p, const Fragment& f_in) {
                 Fragment f = f_in;
                 f.v2 = static_cast<float>(i);
                 return fragment_shader(p, f);
             };
             
             Scan::rasterize<W>(pipeline, canvas, shape, wrapper, debug_bb);
          }
      }
      
      // Overload for MeshState
      template <int W>
      static void draw(auto& pipeline, Canvas& canvas, const MeshState& mesh, auto fragment_shader, bool debug_bb = false) {
          SDF::FaceScratchBuffer scratch;
          size_t face_offset = 0;
          
          for (size_t i = 0; i < mesh.num_faces; ++i) {
             size_t count = mesh.face_counts[i];
             
             // Create Spans from MeshState
             std::span<const Vector> verts(mesh.vertices.data(), mesh.num_vertices);
             std::span<const int> indices(mesh.faces + face_offset, count);
             
             SDF::Face shape(verts, indices, 0.0f, scratch);
             
             auto wrapper = [&](const Vector& p, const Fragment& f_in) {
                 Fragment f = f_in;
                 f.v2 = static_cast<float>(i);
                 return fragment_shader(p, f);
             };
             
             Scan::rasterize<W>(pipeline, canvas, shape, wrapper, debug_bb);
             face_offset += count;
          }
      }
  };
}
