/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <tuple>
#include <utility>
#include <type_traits>
#include <cmath>
#include <vector>
#include <algorithm>
#include "geometry.h" 
#include "color.h" 
#include "static_circular_buffer.h"
#include "canvas.h"

 // Filter Traits
struct Is2D {
  static constexpr bool is_2d = true;
  static constexpr bool has_history = false;
};
struct Is3D {
  static constexpr bool is_2d = false;
  static constexpr bool has_history = false;
};

struct Is2DWithHistory {
  static constexpr bool is_2d = true;
  static constexpr bool has_history = true;
};

struct Is3DWithHistory {
  static constexpr bool is_2d = false;
  static constexpr bool has_history = true;
};

template <int W, int H, typename... Filters>
struct Pipeline;

/**
 * @brief Plots a color to the canvas, ensuring the y-coordinate is within bounds (0 to H-1).
 */
inline void plot_virtual(Canvas& canvas, int x, int y, const Pixel& c) {
  if (y >= 0 && y < canvas.height()) {
    canvas(x, y) = c;
  }
}

// Base Case: Canvas Sink
template <int W, int H>
struct Pipeline<W, H> {
  static constexpr bool is_2d = true;

  // 2D Sink
  void plot(Canvas& cv, int x, int y, const Pixel& c, float age, float alpha, uint8_t tag = 0) {
    int xi = wrap(x, W);
    Pixel p;
    switch(tag) {
        case BLEND_ADD: 
            p = blend_add_alpha(alpha)(cv(xi, y), c); 
            break;
        case BLEND_MAX:
             p = blend_max_alpha(alpha)(cv(xi, y), c);
            break;
        default:
            p = blend_alpha(alpha)(cv(xi, y), c); // BLEND_OVER
            break;
    }
    plot_virtual(cv, xi, y, p);
  }

  void plot(Canvas& cv, float x, float y, const Pixel& c, float age, float alpha, uint8_t tag = 0) {
    int xi = static_cast<int>(std::round(x));
    int yi = static_cast<int>(std::round(y));
    xi = wrap(xi, W);
    
    Pixel p;
    switch(tag) {
        case BLEND_ADD: 
            p = blend_add_alpha(alpha)(cv(xi, yi), c); 
            break;
        case BLEND_MAX:
             p = blend_max_alpha(alpha)(cv(xi, yi), c);
            break;
        default:
            p = blend_alpha(alpha)(cv(xi, yi), c); // BLEND_OVER
            break;
    }
    plot_virtual(cv, xi, yi, p);
  }

  // 3D Sink
  void plot(Canvas& cv, const Vector& v, const Pixel& c, float age, float alpha, uint8_t tag = 0) {
    auto p = vector_to_pixel<W, H>(v);
    plot(cv, p.x, p.y, c, age, alpha, tag);
  }

  // History
  void flush(Canvas&, TrailFn auto, float) {}
};

// Recursive Case
template <int W, int H, typename Head, typename... Tail>
struct Pipeline<W, H, Head, Tail...> : public Head {
  using Next = Pipeline<W, H, Tail...>;
  Next next;

  Pipeline(Head h, Tail... t) : Head(std::move(h)), next(std::move(t)...) {}
  Pipeline() = default;

  // 2D Plot (Float)
  void plot(Canvas& cv, float x, float y, const Pixel& c, float age, float alpha, uint8_t tag = 0) {
    if constexpr (Head::is_2d) {
      Head::plot(x, y, c, age, alpha, tag,
        [&](float x, float y, const Pixel& c, float age, float alpha, uint8_t tag) {
          next.plot(cv, x, y, c, age, alpha, tag);
        });
    }
    else { // 2D -> 3D Mismatch
      Vector v = pixel_to_vector<W, H>(x, y);
      plot(cv, v, c, age, alpha, tag);
    }
  }

  // 2D Plot (Int) - useful for Scan
  void plot(Canvas& cv, int x, int y, const Pixel& c, float age, float alpha, uint8_t tag = 0) {
      // route to float
      plot(cv, static_cast<float>(x), static_cast<float>(y), c, age, alpha, tag);
  }

  // 3D Plot
  void plot(Canvas& cv, const Vector& v, const Pixel& c, float age, float alpha, uint8_t tag = 0) {
    if constexpr (!Head::is_2d) {
      Head::plot(v, c, age, alpha, tag,
        [&](const Vector& v, const Pixel& c, float age, float alpha, uint8_t tag) {
          next.plot(cv, v, c, age, alpha, tag);
        });
    }
    else { // 3D -> 2D Mismatch
      auto p = vector_to_pixel<W, H>(v);
      plot(cv, p.x, p.y, c, age, alpha, tag);
    }
  }

  void flush(Canvas& cv, TrailFn auto trailFn, float alpha) {
    if constexpr (Head::has_history) {
      Head::flush(trailFn, alpha,
        [&](auto... args) {
          next.plot(cv, args...);
        });
    }
    next.flush(cv, trailFn, alpha);
  }
};

///////////////////////////////////////////////////////////////////////////////
// Namespace: Filter
///////////////////////////////////////////////////////////////////////////////

namespace Filter {

// ----------------------------------------------------------------------------
// World Space Filters (3D)
// ----------------------------------------------------------------------------
namespace World {

template <int W>
class Orient : public Is3D {
public:
  Orient(Orientation<W>& orientation) : orientation(orientation) {}

  void plot(const Vector& v, const Pixel& color, float age, float alpha, uint8_t tag, auto pass) {
    tween(orientation, [&](const Quaternion& q, float t) {
      pass(rotate(v, q), color, age + (1.0f - t), alpha, tag);
    });
  }
private:
  Orientation<W>& orientation;
};

template <int W>
class OrientSlice : public Is3D {
public:
  OrientSlice(const std::vector<Orientation<W>>& orientations, const Vector& axis) 
      : enabled(true), axis(axis), orientations(orientations) {}

  void plot(const Vector& v, const Pixel& color, float age, float alpha, uint8_t tag, auto pass) {
    if (!enabled) {
        pass(v, color, age, alpha, tag);
        return;
    }

    float projection = v.i * axis.i + v.j * axis.j + v.k * axis.k; // dot product
    float dot_val = std::max(-1.0f, std::min(1.0f, projection));
    float t = 1.0f - acosf(dot_val) / PI_F;
    
    size_t count = orientations.size();
    if (count == 0) {
        pass(v, color, age, alpha, tag);
        return;
    }
    
    size_t idx = static_cast<size_t>(floorf(t * count));
    if (idx >= count) idx = count - 1;
    
    // Pass to selected orientation
    const Orientation<W>& q = orientations[idx];
    tween(q, [&](const Quaternion& rot, float tween_t) {
        pass(rotate(v, rot), color, age + (1.0f - tween_t), alpha, tag);
    });
  }

  bool enabled;
  Vector axis;

private:
  const std::vector<Orientation<W>>& orientations;
};

template <int W>
class Hole : public Is3D {
public:
  Hole(const Vector& origin, float radius) : origin(origin), radius(radius) {}
  void plot(const Vector& v, const Pixel& color, float age, float alpha, uint8_t tag, auto pass) {
    float d = angle_between(v, origin);
    if (d > radius) pass(v, color, age, alpha, tag);
    else {
      float t = d / radius;
      pass(v, color * quintic_kernel(t), age, alpha, tag);
    }
  }
private:
  Vector origin;
  float radius;
};

template <int W>
class HoleRef : public Is3D {
public:
  HoleRef(const Vector& origin, float radius) : origin(origin), radius(radius) {}
  void plot(const Vector& v, const Pixel& color, float age, float alpha, uint8_t tag, auto pass) {
    float d = angle_between(v, origin);
    if (d > radius) pass(v, color, age, alpha, tag);
    else {
      float t = d / radius;
      pass(v, color * quintic_kernel(t), age, alpha, tag);
    }
  }
private:
  const Vector& origin;
  float radius;
};

template <int W>
class Replicate : public Is3D {
public:
  Replicate(int count) : count(std::clamp(count, 1, W)), step(make_rotation(Y_AXIS, 2 * PI_F / count)) {}
  void plot(const Vector& v, const Pixel& color, float age, float alpha, uint8_t tag, auto pass) {
    Vector r = v;
    pass(r, color, age, alpha, tag);
    for (int i = 1; i < count; i++) {
      r = rotate(r, step);
      pass(r, color, age, alpha, tag);
    }
  }
private:
  int count;
  Quaternion step;
};

template <int W>
class Mobius : public Is3D {
public:
  Mobius(MobiusParams& params) : params(params) {}
  void plot(const Vector& v, const Pixel& color, float age, float alpha, uint8_t tag, auto pass) {
    pass(inv_stereo(mobius(stereo(v), params)), color, age, alpha, tag);
  }
private:
  MobiusParams& params;
};

/**
 * @brief Manages 3D world-space trails.
 */
template <int W, int Capacity>
class Trails : public Is3DWithHistory {
public:
  struct Item {
    Vector v;
    float ttl;
    uint8_t tag;
  };

  Trails(int lifetime) : lifetime(lifetime) {}

  void plot(const Vector& v, const Pixel& color, float age, float alpha, uint8_t tag, auto pass) {
      pass(v, color, age, alpha, tag); // Pass through current frame

      float ttl = lifetime - age;
      if (ttl > 0) {
          items.push_back({ v, ttl, tag }); 
      }
  }

  void flush(TrailFn auto trailFn, float alpha, auto pass) {
      // Age (In-place)
      size_t count = items.size();
      for (size_t i = 0; i < count; ++i) {
          items[i].ttl -= 1.0f;
      }

      // Cleanup (Head-only, as it's a circular buffer and earliest pushed are at head)
      // JS behavior: "Remove Dead" check from front
      while (!items.is_empty() && items.front().ttl <= 0) {
          items.pop();
      }

      // Draw
      // JS: Iterates entire buffer (which is efficiently packed now)
      // No sort needed if we accept draw order = insertion order (standard for trails)
      for (size_t i = 0; i < items.size(); ++i) {
          const auto& item = items[i];
          float t = 1.0f - (item.ttl / static_cast<float>(lifetime));
          Color4 c = trailFn(item.v, t); // Trail function returns color + alpha
          
          // Note: JS logic ignores original point alpha and uses trace alpha * global flush alpha
          if (c.alpha > 0.001f) {
             pass(item.v, c.color, lifetime - item.ttl, c.alpha * alpha, item.tag); 
          }
      }
  }

private:
  StaticCircularBuffer<Item, Capacity> items;
  int lifetime;
};

} // namespace World

// ----------------------------------------------------------------------------
// Screen Space Filters (2D)
// ----------------------------------------------------------------------------
namespace Screen {

template <int W, int H>
class AntiAlias : public Is2D {
public:
  AntiAlias() {}
  void plot(float x, float y, const Pixel& c, float age, float alpha, uint8_t tag, auto pass) {
    float x_i, y_i;
    float x_m = std::modf(x, &x_i);
    float y_m = std::modf(y, &y_i);
    float xs = quintic_kernel(x_m);
    float ys = quintic_kernel(y_m);

    float v00 = (1 - xs) * (1 - ys);
    float v10 = xs * (1 - ys);
    float v01 = (1 - xs) * ys;
    float v11 = xs * ys;

    if (v00 > 1e-4) pass(x_i, y_i, c, age, alpha * v00, tag);
    if (v10 > 1e-4) pass(wrap((x_i + 1), W), y_i, c, age, alpha * v10, tag);
    if (y_i < H - 1) {
      if (v01 > 1e-4) pass(x_i, y_i + 1, c, age, alpha * v01, tag);
      if (v11 > 1e-4) pass(wrap((x_i + 1), W), y_i + 1, c, age, alpha * v11, tag);
    }
  }
};

/**
 * @brief Manages 2D screen-space trails.
 */
template <int W, int MAX_PIXELS = 1024>
class Trails : public Is2DWithHistory {
public:
  Trails(int lifetime) : lifetime(lifetime), num_pixels(0) {}

  void plot(float x, float y, const Pixel& color, float age, float alpha, uint8_t tag, auto pass) {
    if (alpha <= 0.001f) return;

    if (age >= 0) {
      if (num_pixels < MAX_PIXELS) {
        ttls[num_pixels++] = { x, y, static_cast<float>(lifetime - age), tag };
      }
    }
    if (age <= 0) {
      pass(x, y, color, age, alpha, tag);
    }
  }

  void flush(TrailFn auto trailFn, float alpha, auto pass) {
    for (int i = 0; i < num_pixels; ++i) {
      Color4 color = trailFn(ttls[i].x, ttls[i].y, 1 - (ttls[i].ttl / lifetime));
      if (color.alpha > 0.001f) {
        pass(ttls[i].x, ttls[i].y, color.color, lifetime - ttls[i].ttl, alpha * color.alpha, ttls[i].tag);
      }
    }
    decay();
  }

  void decay() {
    for (int i = 0; i < num_pixels; ++i) {
      if (--ttls[i].ttl < TOLERANCE) {
        ttls[i] = ttls[--num_pixels];
        i--;
      }
    }
  }

private:
  struct DecayPixel { float x, y, ttl; uint8_t tag; };
  int lifetime;
  std::array<DecayPixel, MAX_PIXELS> ttls;
  int num_pixels;
};


template <int W, int Capacity>
class Temporal : public Is2DWithHistory {
public:
  using TTLFn = std::function<float(float x, float y)>;
  
  Temporal(TTLFn ttl_fn) : ttl_fn(ttl_fn) {}
  
  void plot(float x, float y, const Pixel& color, float age, float alpha, uint8_t tag, auto pass) {
      float delay = ttl_fn(x, y);
      if (delay <= 1e-4f) {
          pass(x, y, color, age, alpha, tag);
      } else {
          if (items.size() < Capacity) {
              items.push_back({x, y, color, delay, age, alpha, tag});
          } else {
              // Buffer full, bypass
              pass(x, y, color, age, alpha, tag);
          }
      }
  }
  
  void flush(TrailFn auto trailFn, float alpha, auto pass) {
      // Process pending
      for (size_t i = 0; i < items.size(); ++i) {
          auto& item = items[i];
          item.delay -= 1.0f;
          if (item.delay <= 0) {
              pass(item.x, item.y, item.c, item.age, item.alpha * alpha, item.tag);
              // Swap remove
              items[i] = items.back();
              items.pop_back();
              i--;
          }
      }
  }

private:
  struct Item { float x, y; Pixel c; float delay; float age; float alpha; uint8_t tag; };
  StaticCircularBuffer<Item, Capacity> items;
  TTLFn ttl_fn;
};

/**
 * @brief Applies a variable 3x3 Gaussian Blur.
 */
template <int W, int H>
class Blur : public Is2D {
public:
  Blur(float factor = 1.0f) {
    update(factor);
  }

  void update(float factor) {
    float f = std::clamp(factor, 0.0f, 1.0f);
    // Gaussian reference: Corner=1/16, Edge=2/16, Center=4/16
    float c = 1.0f - (0.75f * f); 
    float e = 0.125f * f;        
    float d = 0.0625f * f;       

    kernel = { d, e, d, e, c, e, d, e, d };
  }

  void plot(float x, float y, const Pixel& color, float age, float alpha, uint8_t tag, auto pass) {
    int cx = static_cast<int>(std::round(x));
    int cy = static_cast<int>(std::round(y));

    int k = 0;
    for (int dy = -1; dy <= 1; dy++) {
      int ny = cy + dy;

      if (ny >= 0 && ny < H) {
        for (int dx = -1; dx <= 1; dx++) {
          float weight = kernel[k++];
          if (weight > 1e-5f) {
            pass(static_cast<float>(wrap(cx + dx, W)), static_cast<float>(ny), color, age, alpha * weight, tag);
          }
        }
      }
      else {
        k += 3; // Skip row
      }
    }
  }

private:
  std::array<float, 9> kernel;
};

/**
 * @brief Screen Space Slew Limiter (Sparse/Circular Buffer).
 */
template <int W, int Capacity>
class Slew : public Is2DWithHistory {
public:
    Slew(float rise = 1.0f, float fall = 0.05f) : rise(rise), fall(fall) {}

    void plot(float x, float y, const Pixel& color, float age, float alpha, uint8_t tag, auto pass) {
        pass(x, y, color, age, alpha, tag);

        // Add to decay buffer
        if (items.size() < Capacity && alpha > 0.001f) {
            items.push_back({ x, y, color, alpha, tag });
        }
    }

    void flush(TrailFn auto, float globalAlpha, auto pass) {
        size_t count = items.size();
        for(size_t i = 0; i < count; ++i) {
            auto item = items.front();
            items.pop();

            // Decay
            item.alpha -= fall;

            if (item.alpha > 0.00001f) {
                // Still alive
                pass(item.x, item.y, item.c, 0, item.alpha * globalAlpha, item.tag);
                items.push_back(item);
            }
        }
    }

private:
    struct Item { float x, y; Pixel c; float alpha; uint8_t tag; };
    StaticCircularBuffer<Item, Capacity> items;
    float rise, fall;
};

} // namespace Screen

// ----------------------------------------------------------------------------
// Pixel Space Filters (1:1 with buffer)
// ----------------------------------------------------------------------------
namespace Pix {

template <int W>
class ChromaticShift : public Is2D {
public:
  ChromaticShift() {}
  
  void plot(float x, float y, const ::Pixel& c, float age, float alpha, uint8_t tag, auto pass) {      
      pass(x, y, c, age, alpha, tag);
      
      ::Pixel r_col = c; r_col.g = 0; r_col.b = 0;
      ::Pixel g_col = c; g_col.r = 0; g_col.b = 0;
      ::Pixel b_col = c; b_col.r = 0; b_col.g = 0;
      
      pass(static_cast<float>(wrap(static_cast<int>(x) + 1, W)), y, r_col, age, alpha, tag);
      pass(static_cast<float>(wrap(static_cast<int>(x) + 2, W)), y, g_col, age, alpha, tag);
      pass(static_cast<float>(wrap(static_cast<int>(x) + 3, W)), y, b_col, age, alpha, tag);
  }
};

} // namespace Pix

} // namespace Filter