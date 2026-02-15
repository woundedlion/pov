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
#include "constants.h"
#include "canvas.h"
#include "animation.h"

namespace Plot {

  // Reusable buffers to prevent allocation in tight loops
  static std::vector<float> _steps_cache;
  struct CachePoint { Vector p; float t; };
  static std::vector<CachePoint> _point_cache;

  struct Point {
    static void draw(auto& pipeline, Canvas& canvas, const Vector& v, FragmentShaderFn auto fragment_shader, float age = 0.0f) {
      Fragment f; f.pos = v; f.v0 = 0; f.v1 = 0; f.v2 = 0; f.v3 = 0; f.age = age; f.color = Color4(0,0,0,0); f.blend = 0;
      fragment_shader(v, f);
      pipeline.plot(canvas, v, f.color.color, f.age, f.color.alpha, f.blend);
    }
  };

  template <int W, int H>
  static void rasterize(auto& pipeline, Canvas& canvas, const Fragments& points, FragmentShaderFn auto fragment_shader, bool close_loop = false, float age = 0.0f, const Basis* planar_basis = nullptr) {
    size_t len = points.size();
    if (len < 2) return;

    size_t count = close_loop ? len : len - 1;
    _steps_cache.reserve(W);

    // Optimized segment processor
    auto process_segment = [&](auto&& map, const Fragment& curr, const Fragment& next, float total_dist, bool isLastSegment) {
        // Handle Degenerate Segment
        if (total_dist < 1e-5f) {
            bool shouldOmit = (close_loop) ? false : !isLastSegment; 
            if (!shouldOmit) {
                 Fragment f_copy = curr;
                 // Set temp values for shader
                 f_copy.pos = curr.pos;
                 f_copy.age = age;
                 f_copy.color = Color4(0,0,0,0);
                 f_copy.blend = 0;

                 fragment_shader(curr.pos, f_copy);
                 pipeline.plot(canvas, curr.pos, f_copy.color.color, age, f_copy.color.alpha, f_copy.blend);
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
            // match JS: Math.max(0.05, Math.sqrt(Math.max(0, 1.0 - pTemp.y * pTemp.y)));
            float scale_factor = std::max(0.05f, sqrtf(std::max(0.0f, 1.0f - p_temp.j * p_temp.j)));
            float step = base_step * scale_factor;
            
            _steps_cache.push_back(step);
            sim_dist += step;
            
            if (sim_dist < total_dist) {
                p_temp = map(sim_dist / total_dist);
            }
        }
        
        float scale = (sim_dist > 0.0f) ? (total_dist / sim_dist) : 0.0f;
        bool omitLast = (close_loop) ? false : !isLastSegment;
        
        if (omitLast && _steps_cache.empty()) return;

        // 2. DRAWING PHASE
        
        // Draw Start Point
        {
            Vector start_pos = map(0.0f);
            Fragment f = Fragment::lerp(curr, next, 0.0f);
            f.pos = start_pos;
            f.age = age;
            f.color = Color4(0,0,0,0);
            f.blend = 0;
            
            fragment_shader(start_pos, f);
            pipeline.plot(canvas, start_pos, f.color.color, age, f.color.alpha, f.blend);
        }

        size_t loop_limit = omitLast ? _steps_cache.size() - 1 : _steps_cache.size();
        float current_dist = 0.0f;
        
        for (size_t j = 0; j < loop_limit; j++) {
             float step = _steps_cache[j] * scale;
             current_dist += step;
             
             float t = (total_dist > 0.0f) ? (current_dist / total_dist) : 1.0f;
             
             Vector p = map(t);
             Fragment f = Fragment::lerp(curr, next, t);
             f.pos = p;
             f.age = age;
             f.color = Color4(0,0,0,0);
             f.blend = 0;

             fragment_shader(p, f);
             pipeline.plot(canvas, p, f.color.color, age, f.color.alpha, f.blend);
        }
    };

    for (size_t i = 0; i < count; i++) {
        const Fragment& curr = points[i];
        const Fragment& next = points[(i + 1) % len];
        bool isLastSegment = (i == count - 1);
        
        // --- Interpolation Strategy Selection ---
        if (planar_basis) {
             // PLANAR STRATEGY
             // match JS: PlanarStrategy(planarBasis)
             const Vector& u = planar_basis->u;
             const Vector& center = planar_basis->v;
             const Vector& w = planar_basis->w;
             
             auto project = [&](const Vector& p) -> std::pair<float, float> {
                 float R = angle_between(p, center);
                 if (R < 1e-5f) return {0.0f, 0.0f};
                 float x = dot(p, u);
                 float y = dot(p, w);
                 float theta = atan2f(y, x);
                 return {R * cosf(theta), R * sinf(theta)};
             };
             
             auto proj1 = project(curr.pos);
             auto proj2 = project(next.pos);
             float dx = proj2.first - proj1.first;
             float dy = proj2.second - proj1.second;
             float dist = sqrtf(dx*dx + dy*dy);
             
             auto map_planar = [=](float t) {
                  float Px = proj1.first + dx * t;
                  float Py = proj1.second + dy * t;
                  
                  // Unproject
                  float R = sqrtf(Px*Px + Py*Py);
                  if (R < 1e-5f) return center;
                  
                  float theta = atan2f(Py, Px);
                  Vector axis = (u * cosf(theta)) + (w * sinf(theta));
                  return (center * cosf(R)) + (axis * sinf(R));
             };
             
             process_segment(map_planar, curr, next, dist, isLastSegment);

        } else {
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
                     axis = (std::abs(dot(v1, X_AXIS)) > 0.999f) ? cross(v1, Y_AXIS) : cross(v1, X_AXIS);
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
    static Fragments sample(const Vector& v1, const Vector& v2, int density = 1) {
       if (density < 1) density = 1;
       Fragments points;
       points.reserve(density + 1);
       
       float angle = angle_between(v1, v2);
       if (std::abs(angle) < 0.0001f) {
           Fragment f; f.pos = v1; f.v0 = 0.0f; f.v1 = 0.0f; f.v2 = 0.0f;
           points.push_back(f);
           points.push_back(f); // Draw at least a dot
           return points;
       }
       
       Vector axis = cross(v1, v2).normalize();
       
       for(int i=0; i<=density; ++i) {
           float t = static_cast<float>(i) / density;
           
           Fragment f;
           if (i == 0) f.pos = v1;
           else if (i == density) f.pos = v2;
           else {
               Quaternion q = make_rotation(axis, angle * t);
               f.pos = rotate(v1, q);
           }
           
           f.v0 = t;
           f.v1 = angle * t;
           f.v2 = 0.0f; // Standard: Single lines have v2=0
           points.push_back(f);
       }
       return points;
    }

    /**
     * @brief Draws a geodesic line.
     * @tparam W Rasterization resolution.
     * @param pipeline Render pipeline.
     * @param canvas Target canvas.
     * @param v1 Start point.
     * @param v2 End point.
     * @param fragment_shader Shader function.
     * @param vertex_shader Optional vertex shader.
     */
    template <int W, int H>
    static void draw(auto& pipeline, Canvas& canvas, const Vector& v1, const Vector& v2, FragmentShaderFn auto fragment_shader, VertexShaderFn auto vertex_shader) {
       Fragments points = sample(v1, v2);
       
       if constexpr (!std::is_same_v<decltype(vertex_shader), NullVertexShader>) {
           for(auto& p : points) {
               vertex_shader(p);
           }
       }
       
       rasterize<W, H>(pipeline, canvas, points, fragment_shader, false, 0.0f, nullptr);
    }

    template <int W, int H>
    static void draw(auto& pipeline, Canvas& canvas, const Vector& v1, const Vector& v2, FragmentShaderFn auto fragment_shader) {
       draw<W, H>(pipeline, canvas, v1, v2, fragment_shader, NullVertexShader{});
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
    static void draw(auto& pipeline, Canvas& canvas, const Points& points, FragmentShaderFn auto fragment_shader) {
      for (const Vector& v : points) {
        Fragment f;
        f.pos = v;
        // Basic default registers for points
        f.v0 = 0.0f; f.v1 = 0.0f; f.v2 = 0.0f; f.v3 = 0.0f; f.age = 0.0f; 
        
        fragment_shader(v, f);
        pipeline.plot(canvas, v, f.color.color, f.age, f.color.alpha, f.blend);
      }
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
    static Vector calcPoint(float a, float radius, const Vector& u, const Vector& v, const Vector& w) {
      auto d = sqrtf((1 - radius) * (1 - radius));
      return Vector(
        d * v.i + radius * u.i * cosf(a) + radius * w.i * sinf(a),
        d * v.j + radius * u.j * cosf(a) + radius * w.j * sinf(a),
        d * v.k + radius * u.k * cosf(a) + radius * w.k * sinf(a)
      ).normalize();
    }

    static void sample(Fragments& points, const Basis& basis, float radius, int num_samples, float phase = 0) {
      auto res = get_antipode(basis, radius);
      const Basis& work_basis = res.first;
      float work_radius = res.second;

      const Vector& v = work_basis.v;
      const Vector& u = work_basis.u;
      const Vector& w = work_basis.w;

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
    static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, FragmentShaderFn auto fragment_shader, VertexShaderFn auto vertex_shader, float phase = 0) {
      Fragments points;
      // Use W samples for smooth circles (fixes pinching at poles)
      sample(points, basis, radius, W, phase);
      
      if constexpr (!std::is_same_v<decltype(vertex_shader), NullVertexShader>) {
          for(auto& p : points) {
              vertex_shader(p);
          }
      }
      rasterize<W, H>(pipeline, canvas, points, fragment_shader, true, 0.0f, nullptr);
    }

    template <int W, int H>
    static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, FragmentShaderFn auto fragment_shader, float phase = 0) {
      draw<W, H>(pipeline, canvas, basis, radius, fragment_shader, NullVertexShader{}, phase);
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
    static void sample(Fragments& points, const Basis& basis, float radius, int num_sides, float phase = 0) {
       size_t start_idx = points.size();       
       Ring::sample(points, basis, radius, num_sides, phase + PI_F / num_sides);
       float cumul = 0.0f;
       for (size_t i = start_idx; i < points.size(); i++) {
           points[i].v1 = cumul;
           if (i < points.size() - 1) {
               cumul += angle_between(points[i].pos, points[i+1].pos);
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
    static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, int num_sides, FragmentShaderFn auto fragment_shader, VertexShaderFn auto vertex_shader, float phase = 0) {
      Fragments points;
      sample(points, basis, radius, num_sides, phase);
      
      if constexpr (!std::is_same_v<decltype(vertex_shader), NullVertexShader>) {
          for(auto& p : points) {
              vertex_shader(p);
          }
      }
      
      rasterize<W, H>(pipeline, canvas, points, fragment_shader, true, 0.0f, &basis);
    }

    template <int W, int H>
    static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, int num_sides, FragmentShaderFn auto fragment_shader, float phase = 0) {
      draw<W, H>(pipeline, canvas, basis, radius, num_sides, fragment_shader, NullVertexShader{}, phase);
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
    static void sample(Fragments& points, const Basis& basis, float radius, int num_sides, float phase = 0) {
       size_t start_idx = points.size();
       Ring::sample(points, basis, radius, num_sides, phase + PI_F / num_sides);
       
       // Re-calculate v1 to be true geodesic chord length
       float cumulative_length = 0.0f;
       for (size_t i = start_idx; i < points.size(); ++i) {
           points[i].v2 = static_cast<float>(i - start_idx);
           
           if (i > start_idx) {
               cumulative_length += angle_between(points[i-1].pos, points[i].pos);
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
    static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, int num_sides, FragmentShaderFn auto fragment_shader, VertexShaderFn auto vertex_shader, float phase = 0) {
       Fragments points;
       sample(points, basis, radius, num_sides, phase);
       
       if constexpr (!std::is_same_v<decltype(vertex_shader), NullVertexShader>) {
           for(auto& p : points) {
               vertex_shader(p);
           }
       }
       rasterize<W, H>(pipeline, canvas, points, fragment_shader, true, 0.0f, nullptr);
    }

    template <int W, int H>
    static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, int num_sides, FragmentShaderFn auto fragment_shader, float phase = 0) {
       draw<W, H>(pipeline, canvas, basis, radius, num_sides, fragment_shader, NullVertexShader{}, phase);
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
    static Vector fn_point(std::function<float(float)> shift_fn, const Basis& basis, float radius, float angle) {
      auto res = get_antipode(basis, radius);
      const Basis& work_basis = res.first;
      float work_radius = res.second;
      
      const Vector& v = work_basis.v;
      const Vector& u = work_basis.u;
      const Vector& w = work_basis.w;

      auto vi = Ring::calcPoint(angle, work_radius, u, v, w);
      auto vp = Ring::calcPoint(angle, 1, u, v, w);
      Vector axis = cross(v, vp).normalize();
      return rotate(vi, make_rotation(axis, shift_fn(angle * PI_F / 2)));
    }

    template <int W>
    static void sample(Fragments& points, const Basis& basis, float radius, std::function<float(float)> shift_fn, float phase = 0) {
      auto res = get_antipode(basis, radius);
      const Basis& work_basis = res.first;
      float work_radius = res.second;

      const Vector& v = work_basis.v;
      const Vector& u = work_basis.u;
      const Vector& w = work_basis.w;

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
          const Fragment& first = points[0];
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
    static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, std::function<float(float)> shift_fn, FragmentShaderFn auto fragment_shader, VertexShaderFn auto vertex_shader, float phase = 0) {
      Fragments points;
      sample<W>(points, basis, radius, shift_fn, phase);
      
      if constexpr (!std::is_same_v<decltype(vertex_shader), NullVertexShader>) {
          for(auto& p : points) {
              vertex_shader(p);
          }
      }
      rasterize<W, H>(pipeline, canvas, points, fragment_shader, true, 0.0f, nullptr);
    }
    
    template <int W, int H>
    static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, std::function<float(float)> shift_fn, FragmentShaderFn auto fragment_shader, float phase = 0) {
      draw<W, H>(pipeline, canvas, basis, radius, shift_fn, fragment_shader, NullVertexShader{}, phase);
    }
  };

  struct Spiral {
    /**
     * @brief Samples a Fibonacci spiral.
     */
    static Fragments sample(int n, float eps) {
        Fragments fragments;
        fragments.reserve(n);
        
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
        return fragments;
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
    static void draw(auto& pipeline, Canvas& canvas, int n, float eps, FragmentShaderFn auto fragment_shader, VertexShaderFn auto vertex_shader) {
      Fragments frags = sample(n, eps);
      
      if constexpr (!std::is_same_v<decltype(vertex_shader), NullVertexShader>) {
          for(auto& f : frags) {
              vertex_shader(f);
          }
      }
      
      rasterize<W, H>(pipeline, canvas, frags, fragment_shader, false, 0.0f, nullptr);
    }

    template <int W, int H>
    static void draw(auto& pipeline, Canvas& canvas, int n, float eps, FragmentShaderFn auto fragment_shader) {
      draw<W, H>(pipeline, canvas, n, eps, fragment_shader, NullVertexShader{});
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
    static void sample(Fragments& points, const Basis& basis, float radius, int num_sides, float phase = 0) {
      auto res = get_antipode(basis, radius);
      const Basis& work_basis = res.first;
      float work_radius = res.second;

      const Vector& v = work_basis.v;
      const Vector& u = work_basis.u;
      const Vector& w = work_basis.w;

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
           const Fragment& first = points[start_idx];
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
    static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, int num_sides, FragmentShaderFn auto fragment_shader, VertexShaderFn auto vertex_shader, float phase = 0) {
      Fragments points;
      sample(points, basis, radius, num_sides, phase);
      
      if constexpr (!std::is_same_v<decltype(vertex_shader), NullVertexShader>) {
          for(auto& p : points) {
              vertex_shader(p);
          }
      }
      
      Vector center = basis.v;
      if (radius > 1.0f) center = -center;
      
      Vector v = center;
      Vector ref = (std::abs(dot(v, X_AXIS)) > 0.9f) ? Y_AXIS : X_AXIS;
      Vector u = cross(v, ref).normalize();
      Vector w = cross(v, u).normalize();
      Basis planar_basis = { u, v, w };

      rasterize<W, H>(pipeline, canvas, points, fragment_shader, true, 0.0f, &planar_basis);
    }

    template <int W, int H>
    static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, int num_sides, FragmentShaderFn auto fragment_shader, float phase = 0) {
      draw<W, H>(pipeline, canvas, basis, radius, num_sides, fragment_shader, NullVertexShader{}, phase);
    }
  };

  /**
   * @brief Flower shape.
   * Registers:
   *  v0: Perimeter progress (0.0 -> 1.0)
   *  v1: Arc Length (radians)
   *  v2: Vertex index
   */
  struct Flower {
    /**
     * @brief Samples a flower shape.
     */
    static void sample(Fragments& points, const Basis& basis, float radius, int num_sides, float phase = 0) {
        auto res = get_antipode(basis, radius);
        const Basis& work_basis = res.first;
        float work_radius = res.second;

        const Vector& v = work_basis.v;
        const Vector& u = work_basis.u;
        const Vector& w = work_basis.w;

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
             const Fragment& first = points[start_idx];
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
    static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, int num_sides, FragmentShaderFn auto fragment_shader, VertexShaderFn auto vertex_shader, float phase = 0) {
       Fragments points;
       sample(points, basis, radius, num_sides, phase);

       if constexpr (!std::is_same_v<decltype(vertex_shader), NullVertexShader>) {
           for(auto& p : points) {
               vertex_shader(p);
           }
       }
       
       auto res = get_antipode(basis, radius);
       const Basis& work_basis = res.first;
       Vector center = work_basis.v;
       
       // Construct a planar basis aligned with the center
       Vector ref = (std::abs(dot(center, X_AXIS)) > 0.9f) ? Y_AXIS : X_AXIS;
       Vector u_p = cross(center, ref).normalize();
       Vector w_p = cross(center, u_p).normalize();
       Basis planar_basis = { u_p, center, w_p };

       rasterize<W, H>(pipeline, canvas, points, fragment_shader, true, 0.0f, &planar_basis);
    }

    template <int W, int H>
    static void draw(auto& pipeline, Canvas& canvas, const Basis& basis, float radius, int num_sides, FragmentShaderFn auto fragment_shader, float phase = 0) {
      draw<W, H>(pipeline, canvas, basis, radius, num_sides, fragment_shader, NullVertexShader{}, phase);
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
    template <typename MeshT>
    static std::vector<Fragments> sample(const MeshT& mesh, int density = 10) {
       std::set<std::pair<int, int>> visited;
       std::vector<Fragments> result;
       
       auto process_edge = [&](int u, int v) {
           int small = std::min(u, v);
           int large = std::max(u, v);
           
           if (visited.find({small, large}) == visited.end()) {
               visited.insert({small, large});
               result.push_back(Line::sample(mesh.vertices[u], mesh.vertices[v], density));
           }
       };

       if constexpr (std::is_same_v<MeshT, MeshState>) {
           size_t offset = 0;
           for (size_t i = 0; i < mesh.face_counts.size(); ++i) {
               int count = mesh.face_counts[i];
               for (int k = 0; k < count; ++k) {
                   int u = mesh.faces[offset + k];
                   int v = mesh.faces[offset + (k + 1) % count];
                   process_edge(u, v); // Helper lambda
               }
               offset += count;
           }
       } else {
           for (const auto& face : mesh.faces) {
               size_t count = face.size();
               for (size_t i = 0; i < count; ++i) {
                   int u = face[i];
                   int v = face[(i + 1) % count];
                   process_edge(u, v);
               }
           }
       }
       
       return result;
    }

    template <int W, int H, typename MeshT>
    static void draw(auto& pipeline, Canvas& canvas, const MeshT& mesh, FragmentShaderFn auto fragment_shader, VertexShaderFn auto vertex_shader) {
       auto edges = sample(mesh, 10);
       
       for (size_t i = 0; i < edges.size(); ++i) {
           auto& points = edges[i];
           if constexpr (!std::is_same_v<decltype(vertex_shader), NullVertexShader>) {
               for(auto& p : points) {
                   p.v2 = static_cast<float>(i); // Edge Index
                   vertex_shader(p);
               }
           } else {
               for(auto& p : points) {
                   p.v2 = static_cast<float>(i);
               }
           }
           rasterize<W, H>(pipeline, canvas, points, fragment_shader, false, 0.0f, nullptr);
       }
    }

    template <int W, int H, typename MeshT>
    static void draw(auto& pipeline, Canvas& canvas, const MeshT& mesh, FragmentShaderFn auto fragment_shader) {
        draw<W, H>(pipeline, canvas, mesh, fragment_shader, NullVertexShader{});
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
      template <int W, int H>
      static std::vector<Fragments> sample(const auto& system) {
          std::vector<Fragments> trails;
          int count = system.active_count;
          trails.reserve(count);

          for (int i = 0; i < count; ++i) {
              const auto& p = system.pool[i];
              Fragments trail;
              float cumulative_len = 0.0f;
              Vector last_pos; 
              bool first = true;

              deep_tween(p.history, [&](const Quaternion& q, float t) {
                  Vector v = rotate(p.position, q);
                  
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
                  f.v3 = p.life / p.max_life;
                  f.age = 0;
                  f.color = Color4(0,0,0,0);
                  trail.push_back(f);
              });
              
              if (!trail.empty()) {
                  trails.push_back(trail);
              }
          }
          return trails;
      }

      template <int W, int H>
      static void draw(auto& pipeline, Canvas& canvas, const auto& system, FragmentShaderFn auto fragment_shader, VertexShaderFn auto vertex_shader) {
          auto trails = sample<W, H>(system);
          
          for (auto& trail : trails) {
              if constexpr (!std::is_same_v<decltype(vertex_shader), NullVertexShader>) {
                  for (auto& f : trail) {
                      vertex_shader(f);
                  }
              }
              rasterize<W, H>(pipeline, canvas, trail, fragment_shader, false, 0.0f, nullptr);
          }
      }
      
      template <int W, int H>
      static void draw(auto& pipeline, Canvas& canvas, const auto& system, FragmentShaderFn auto fragment_shader) {
          draw<W, H>(pipeline, canvas, system, fragment_shader, NullVertexShader{});
      }
   };

}
