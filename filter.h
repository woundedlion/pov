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
#include "static_circular_buffer.h"



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

template <int W, typename... Filters>
struct Pipeline;

/**
 * @brief Plots a color to the canvas, ensuring the y-coordinate is within bounds (0 to H-1).
 */
void plot_virtual(Canvas& canvas, int x, int y, const Pixel& c) {
  if (y >= 0 && y < H) {
    canvas(XY(x, y)) = c;
  }
}

// Base Case: Canvas Sink
template <int W>
struct Pipeline<W> {
  static constexpr bool is_2d = true;

  // 2D Sink
  void plot(Canvas& cv, int x, int y, const Pixel& c, float age, float alpha, uint8_t tag = 0) {
    int xi = wrap(x, W);
    auto p = blend_alpha(alpha)(cv(xi, y), c); 
    plot_virtual(cv, xi, y, p);
  }

  void plot(Canvas& cv, float x, float y, const Pixel& c, float age, float alpha, uint8_t tag = 0) {
    int xi = static_cast<int>(std::round(x));
    int yi = static_cast<int>(std::round(y));
    xi = wrap(xi, W);
    auto p = blend_alpha(alpha)(cv(xi, yi), c);
    plot_virtual(cv, xi, yi, p);
  }

  // 3D Sink
  void plot(Canvas& cv, const Vector& v, const Pixel& c, float age, float alpha, uint8_t tag = 0) {
    auto p = vector_to_pixel<W>(v);
    plot(cv, p.x, p.y, c, age, alpha, tag);
  }

  // History
  void flush(Canvas&, TrailFn auto, float) {}
};

// Recursive Case
template <int W, typename Head, typename... Tail>
struct Pipeline<W, Head, Tail...> : public Head {
  using Next = Pipeline<W, Tail...>;
  Next next;

  Pipeline(Head h, Tail... t) : Head(std::move(h)), next(std::move(t)...) {}
  Pipeline() = default;

  // 2D Plot (Float)
  void plot(Canvas& cv, float x, float y, const Pixel& c, float age, float alpha, uint8_t tag = 0) {
    if constexpr (Head::is_2d) {
      Head::plot(x, y, c, age, alpha,
        [&](float x, float y, const Pixel& c, float age, float alpha) {
          next.plot(cv, x, y, c, age, alpha, tag);
        });
    }
    else { // 2D -> 3D Mismatch
      Vector v = pixel_to_vector<W>(x, y);
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
      Head::plot(v, c, age, alpha,
        [&](const Vector& v, const Pixel& c, float age, float alpha) {
          next.plot(cv, v, c, age, alpha, tag);
        });
    }
    else { // 3D -> 2D Mismatch
      auto p = vector_to_pixel<W>(v);
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

float quintic_kernel(float t) {
  t = std::clamp(t, 0.0f, 1.0f);
  return t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f);
}

template <int W>
class FilterOrient : public Is3D {
public:
  FilterOrient(Orientation& orientation) : orientation(orientation) {}

  void plot(const Vector& v, const Pixel& color, float age, float alpha, auto pass) {
    tween(orientation, [&](const Quaternion& q, float t) {
      pass(rotate(v, q), color, age + (1.0f - t), alpha);
    });
  }
private:
  Orientation& orientation;
};

template <int W>
class FilterOrientSlice : public Is3D {
public:
  FilterOrientSlice(const std::vector<Orientation>& orientations, const Vector& axis) 
      : enabled(true), axis(axis), orientations(orientations) {}

  void plot(const Vector& v, const Pixel& color, float age, float alpha, auto pass) {
    if (!enabled) {
        pass(v, color, age, alpha);
        return;
    }

    float projection = v.i * axis.i + v.j * axis.j + v.k * axis.k; // dot product
    float dot_val = std::max(-1.0f, std::min(1.0f, projection));
    float t = 1.0f - acosf(dot_val) / PI_F;
    
    size_t count = orientations.size();
    if (count == 0) {
        pass(v, color, age, alpha);
        return;
    }
    
    size_t idx = static_cast<size_t>(floorf(t * count));
    if (idx >= count) idx = count - 1;
    
    // Pass to selected orientation
    const Orientation& q = orientations[idx];
    tween(q, [&](const Quaternion& rot, float tween_t) {
        pass(rotate(v, rot), color, age + (1.0f - tween_t), alpha);
    });
  }

  bool enabled;
  Vector axis;

private:
  const std::vector<Orientation>& orientations;
};

template <int W>
class FilterHole : public Is3D {
public:
  FilterHole(const Vector& origin, float radius) : origin(origin), radius(radius) {}
  void plot(const Vector& v, const Pixel& color, float age, float alpha, auto pass) {
    float d = angle_between(v, origin);
    if (d > radius) pass(v, color, age, alpha);
    else {
      float t = d / radius;
      pass(v, color * quintic_kernel(t), age, alpha);
    }
  }
private:
  Vector origin;
  float radius;
};

template <int W>
class FilterHoleRef : public Is3D {
public:
  FilterHoleRef(const Vector& origin, float radius) : origin(origin), radius(radius) {}
  void plot(const Vector& v, const Pixel& color, float age, float alpha, auto pass) {
    float d = angle_between(v, origin);
    if (d > radius) pass(v, color, age, alpha);
    else {
      float t = d / radius;
      pass(v, color * quintic_kernel(t), age, alpha);
    }
  }
private:
  const Vector& origin;
  float radius;
};

template <int W>
class FilterReplicate : public Is3D {
public:
  FilterReplicate(int count) : count(std::clamp(count, 1, W)), step(make_rotation(Y_AXIS, 2 * PI_F / count)) {}
  void plot(const Vector& v, const Pixel& color, float age, float alpha, auto pass) {
    Vector r = v;
    pass(r, color, age, alpha);
    for (int i = 1; i < count; i++) {
      r = rotate(r, step);
      pass(r, color, age, alpha);
    }
  }
private:
  int count;
  Quaternion step;
};

template <int W>
class FilterMobius : public Is3D {
public:
  FilterMobius(MobiusParams& params) : params(params) {}
  void plot(const Vector& v, const Pixel& color, float age, float alpha, auto pass) {
    pass(inv_stereo(mobius(stereo(v), params)), color, age, alpha);
  }
private:
  MobiusParams& params;
};

///////////////////////////////////////////////////////////////////////////////
// 2D Filters
///////////////////////////////////////////////////////////////////////////////

template <int W>
class FilterAntiAlias : public Is2D {
public:
  FilterAntiAlias() {}
  void plot(float x, float y, const Pixel& c, float age, float alpha, auto pass) {
    float x_i, y_i;
    float x_m = std::modf(x, &x_i);
    float y_m = std::modf(y, &y_i);
    float xs = quintic_kernel(x_m);
    float ys = quintic_kernel(y_m);

    float v00 = (1 - xs) * (1 - ys);
    float v10 = xs * (1 - ys);
    float v01 = (1 - xs) * ys;
    float v11 = xs * ys;

    if (v00 > 1e-4) pass(x_i, y_i, c, age, alpha * v00);
    if (v10 > 1e-4) pass(wrap((x_i + 1), W), y_i, c, age, alpha * v10);
    if (y_i < H - 1) {
      if (v01 > 1e-4) pass(x_i, y_i + 1, c, age, alpha * v01);
      if (v11 > 1e-4) pass(wrap((x_i + 1), W), y_i + 1, c, age, alpha * v11);
    }
  }
};

template <int W>
class FilterChromaticShift : public Is2D {
public:
  FilterChromaticShift() {}
  
  void plot(float x, float y, const Pixel& c, float age, float alpha, auto pass) {      
      pass(x, y, c, age, alpha);
      
      Pixel r_col = c; r_col.g = 0; r_col.b = 0;
      Pixel g_col = c; g_col.r = 0; g_col.b = 0;
      Pixel b_col = c; b_col.r = 0; b_col.g = 0;
      
      pass(static_cast<float>(wrap(static_cast<int>(x) + 1, W)), y, r_col, age, alpha);
      pass(static_cast<float>(wrap(static_cast<int>(x) + 2, W)), y, g_col, age, alpha);
      pass(static_cast<float>(wrap(static_cast<int>(x) + 3, W)), y, b_col, age, alpha);
  }
};

/**
 * @brief Renamed from FilterDecay. Manages 2D screen-space trails.
 */
template <int W, int MAX_PIXELS>
class FilterScreenTrails : public Is2DWithHistory {
public:
  FilterScreenTrails(int lifetime) : lifetime(lifetime), num_pixels(0) {}

  void plot(float x, float y, const Pixel& color, float age, float alpha, auto pass) {
    if (age >= 0) {
      if (num_pixels < MAX_PIXELS) {
        ttls[num_pixels++] = { x, y, static_cast<float>(lifetime - age) };
      }
    }
    if (age <= 0) {
      pass(x, y, color, age, alpha);
    }
  }

  void flush(TrailFn auto trailFn, float alpha, auto pass) {
    for (int i = 0; i < num_pixels; ++i) {
      Color4 color = trailFn(ttls[i].x, ttls[i].y, 1 - (ttls[i].ttl / lifetime));
      pass(ttls[i].x, ttls[i].y, color.color, lifetime - ttls[i].ttl, alpha * color.alpha);
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
  struct DecayPixel { float x, y, ttl; };
  int lifetime;
  std::array<DecayPixel, MAX_PIXELS> ttls;
  int num_pixels;
};


/**
 * @brief Applies a variable 3x3 Gaussian Blur.
 */
template <int W>
class FilterGaussianBlur : public Is2D {
public:
  FilterGaussianBlur(float factor = 1.0f) {
    update(factor);
  }

  void update(float factor) {
    float f = std::clamp(factor, 0.0f, 1.0f);

    // Interpolate weights between Identity (Center=1) and Gaussian (Center=0.25)
    // Gaussian reference: Corner=1/16, Edge=2/16, Center=4/16
    float c = 1.0f - (0.75f * f); // Center weight: 1.0 -> 0.25
    float e = 0.125f * f;        // Edge weight:   0.0 -> 0.125
    float d = 0.0625f * f;       // Diagonal weight: 0.0 -> 0.0625

    kernel = {
      d, e, d,
      e, c, e,
      d, e, d
    };
  }

  void plot(float x, float y, const Pixel& color, float age, float alpha, auto pass) {
    int cx = static_cast<int>(std::round(x));
    int cy = static_cast<int>(std::round(y));

    int k = 0;
    for (int dy = -1; dy <= 1; dy++) {
      int ny = cy + dy;

      if (ny >= 0 && ny < H) {
        for (int dx = -1; dx <= 1; dx++) {
          float weight = kernel[k++];
          if (weight > TOLERANCE) {
            pass(static_cast<float>(wrap(cx + dx, W)), static_cast<float>(ny), color, age, alpha * weight);
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
 * @brief New FilterWorldTrails. Manages 3D world-space trails.
 */
template <int W, int Capacity>
class FilterWorldTrails : public Is3DWithHistory {
public:
  struct Item {
    Vector v;
    Color4 color;
    float ttl;
    float alpha;
  };

  FilterWorldTrails(int lifetime) : lifetime(lifetime) {}

  void plot(const Vector& v, const Pixel& color, float age, float alpha, auto pass) {
      float ttl = lifetime - age;
      if (ttl > 0) {
          items.push_back({ v, Color4(color, 1.0f), ttl, alpha }); 
          // Note: Pixel is color. Color4 is different? 
          // Color4 has alpha. Pixel is just color?
          // Pixel color is usually CRGB, no alpha. Alpha is separate.
      }
      
      if (age <= 0) {
          pass(v, color, age, alpha);
      }
  }

  void flush(TrailFn auto trailFn, float alpha, auto pass) {
      // Sort
      std::sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
          return a.ttl < b.ttl;
      });

      for (auto& item : items) {
          float t = (lifetime - item.ttl) / static_cast<float>(lifetime);
          Color4 color = trailFn(item.v, t);
          pass(item.v, color.color, 0, item.alpha * color.alpha * alpha); // age 0 because its a 'fresh' plot from trail
          item.ttl -= 1.0f;
      }

      // Cleanup
      while (!items.is_empty() && items.front().ttl <= 0) {
          items.pop();
      }
  }

private:
  StaticCircularBuffer<Item, Capacity> items;
  int lifetime;
};

template <int W, int Capacity>
class FilterTemporal : public Is2DWithHistory {
public:
  using TTLFn = std::function<float(float x, float y)>;
  
  FilterTemporal(TTLFn ttl_fn) : ttl_fn(ttl_fn) {}
  
  void plot(float x, float y, const Pixel& color, float age, float alpha, auto pass) {
      float delay = ttl_fn(x, y);
      if (delay <= 1e-4f) {
          pass(x, y, color, age, alpha);
      } else {
          if (items.size() < Capacity) {
              items.push_back({x, y, color, delay, age, alpha});
          } else {
              // Buffer full, bypass
              pass(x, y, color, age, alpha);
          }
      }
  }
  
  void flush(TrailFn auto trailFn, float alpha, auto pass) {
      // Process pending
      for (size_t i = 0; i < items.size(); ++i) {
          auto& item = items[i];
          item.delay -= 1.0f;
          if (item.delay <= 0) {
              pass(item.x, item.y, item.c, item.age, item.alpha * alpha);
              // Swap remove
              items[i] = items.back();
              items.pop_back();
              i--;
          }
      }
  }

private:
  struct Item { float x, y; Pixel c; float delay; float age; float alpha; };
  StaticCircularBuffer<Item, Capacity> items;
  TTLFn ttl_fn;
};