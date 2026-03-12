/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <tuple>
#include <utility>
#include <type_traits>
#include <cmath>

#include <span>
#include <algorithm>
#include "geometry.h"
#include "color.h"
#include "static_circular_buffer.h"
#include "canvas.h"
#include "concepts.h"

using PassFn2D =
    FunctionRef<void(float, float, const Pixel &, float, float, uint8_t)>;
using PassFn3D =
    FunctionRef<void(const Vector &, const Pixel &, float, float, uint8_t)>;

// Filter Traits
/** @brief Trait indicating a filter operates in 2D screen space. */
struct Is2D {
  static constexpr bool is_2d = true;
  static constexpr bool has_history = false;
};
/** @brief Trait indicating a filter operates in 3D world space. */
struct Is3D {
  static constexpr bool is_2d = false;
  static constexpr bool has_history = false;
};

/** @brief Trait indicating a 2D filter that maintains state/history. */
struct Is2DWithHistory {
  static constexpr bool is_2d = true;
  static constexpr bool has_history = true;
};

/** @brief Trait indicating a 3D filter that maintains state/history. */
struct Is3DWithHistory {
  static constexpr bool is_2d = false;
  static constexpr bool has_history = true;
};

/**
 * @brief Recursive template pipeline for processing render commands.
 * Chains filters together, routing 2D and 3D plot commands.
 */
template <int W, int H, typename... Filters> struct Pipeline;

/**
 * @brief Plots a color to the canvas, ensuring the y-coordinate is within
 * bounds (0 to H-1).
 */
inline void plot_virtual(Canvas &canvas, int x, int y, const Pixel &c) {
  if (y >= 0 && y < canvas.height()) {
    canvas(x, y) = c;
  }
}

// Base Case: Canvas Sink
/**
 * @brief Terminal node of the pipeline. Writes final pixels to the Canvas.
 */
template <int W, int H> struct Pipeline<W, H> {
  static constexpr bool is_2d = true;

  // 2D Sink
  void plot(Canvas &cv, int x, int y, const Pixel &c, float age, float alpha,
            uint8_t tag = 0) {
    int xi = fast_wrap(x, W);
    Pixel p;
    switch (tag) {
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

  void plot(Canvas &cv, float x, float y, const Pixel &c, float age,
            float alpha, uint8_t tag = 0) {
    int xi = static_cast<int>(std::round(x));
    int yi = static_cast<int>(std::round(y));
    xi = fast_wrap(xi, W);

    Pixel p;
    switch (tag) {
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
  void plot(Canvas &cv, const Vector &v, const Pixel &c, float age, float alpha,
            uint8_t tag = 0) {
    auto p = vector_to_pixel<W, H>(v);
    plot(cv, p.x, p.y, c, age, alpha, tag);
  }

  // History
  void flush(Canvas &, const ScreenTrailFn &, float) {}
  void flush(Canvas &, const WorldTrailFn &, float) {}
};

// Recursive Case
template <int W, int H, typename Head, typename... Tail>
struct Pipeline<W, H, Head, Tail...> : public Head {
  using Next = Pipeline<W, H, Tail...>;
  Next next;

  // Forwarding Reference Constructor
  template <typename HArg, typename... TArgs>
  Pipeline(HArg &&h, TArgs &&...t)
      : Head(std::forward<HArg>(h)), next(std::forward<TArgs>(t)...) {}

  // Partial Constructor (Head only, Tail default constructed)
  template <typename HArg>
  explicit Pipeline(HArg &&h) : Head(std::forward<HArg>(h)) {}

  Pipeline() = default;

  // 2D Plot (Float)
  void plot(Canvas &cv, float x, float y, const Pixel &c, float age,
            float alpha, uint8_t tag = 0) {
    if constexpr (Head::is_2d) {
      Head::plot(x, y, c, age, alpha, tag,
                 [&](float x, float y, const Pixel &c, float age, float alpha,
                     uint8_t tag) { next.plot(cv, x, y, c, age, alpha, tag); });
    } else { // 2D -> 3D Mismatch
      Vector v = pixel_to_vector<W, H>(x, y);
      plot(cv, v, c, age, alpha, tag);
    }
  }

  // 2D Plot (Int) - useful for Scan
  void plot(Canvas &cv, int x, int y, const Pixel &c, float age, float alpha,
            uint8_t tag = 0) {
    // route to float
    plot(cv, static_cast<float>(x), static_cast<float>(y), c, age, alpha, tag);
  }

  // 3D Plot
  void plot(Canvas &cv, const Vector &v, const Pixel &c, float age, float alpha,
            uint8_t tag = 0) {
    if constexpr (!Head::is_2d) {
      Head::plot(v, c, age, alpha, tag,
                 [&](const Vector &v, const Pixel &c, float age, float alpha,
                     uint8_t tag) { next.plot(cv, v, c, age, alpha, tag); });
    } else { // 3D -> 2D Mismatch
      auto p = vector_to_pixel<W, H>(v);
      plot(cv, p.x, p.y, c, age, alpha, tag);
    }
  }

  void flush(Canvas &cv, const ScreenTrailFn &trailFn, float alpha) {
    if constexpr (Head::has_history) {
      if constexpr (Head::is_2d) {
        Head::flush(cv, trailFn, alpha,
                    [&](auto... args) { next.plot(cv, args...); });
      }
    }
    next.flush(cv, trailFn, alpha);
  }

  void flush(Canvas &cv, const WorldTrailFn &trailFn, float alpha) {
    if constexpr (Head::has_history) {
      if constexpr (!Head::is_2d) {
        Head::flush(trailFn, alpha,
                    [&](auto... args) { next.plot(cv, args...); });
      }
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

/**
 * @brief Rotates 3D points based on a dynamic Orientation.
 */
template <int W> class Orient : public Is3D {
public:
  Orient(Orientation<W> &orientation) : orientation(orientation) {}

  void plot(const Vector &v, const Pixel &color, float age, float alpha,
            uint8_t tag, PassFn3D pass) {
    tween(orientation, [&](const Quaternion &q, float t) {
      pass(rotate(v, q), color, age + (1.0f - t), alpha, tag);
    });
  }

private:
  Orientation<W> &orientation;
};

/**
 * @brief Selects an orientation from a list based on the point's projection
 * onto an axis. Useful for slicing objects with different rotations.
 */
template <int W> class OrientSlice : public Is3D {
public:
  OrientSlice(std::span<const Orientation<W>> orientations, const Vector &axis)
      : enabled(true), axis(axis), orientations(orientations) {}

  void plot(const Vector &v, const Pixel &color, float age, float alpha,
            uint8_t tag, PassFn3D pass) {
    if (!enabled) {
      pass(v, color, age, alpha, tag);
      return;
    }

    float projection =
        v.x * axis.x + v.y * axis.y + v.z * axis.z; // dot product
    float dot_val = std::max(-1.0f, std::min(1.0f, projection));
    float t = 1.0f - acosf(dot_val) / PI_F;

    size_t count = orientations.size();
    if (count == 0) {
      pass(v, color, age, alpha, tag);
      return;
    }

    size_t idx = static_cast<size_t>(floorf(t * count));
    if (idx >= count)
      idx = count - 1;

    // Pass to selected orientation
    const Orientation<W> &q = orientations[idx];
    tween(q, [&](const Quaternion &rot, float tween_t) {
      pass(rotate(v, rot), color, age + (1.0f - tween_t), alpha, tag);
    });
  }

  bool enabled;
  Vector axis;

private:
  std::span<const Orientation<W>> orientations;
};

/**
 * @brief Creates a spherical hole by masking points within a radius.
 * Parameterized on OriginT to allow storage by value (Vector) or reference
 * (std::reference_wrapper).
 */
template <int W, typename OriginT = Vector> class Hole : public Is3D {
public:
  Hole(OriginT origin, float radius) : origin(origin), radius(radius) {}
  void plot(const Vector &v, const Pixel &color, float age, float alpha,
            uint8_t tag, PassFn3D pass) {
    // Unwrap if reference_wrapper, or use directly if Vector
    const Vector &o = origin;
    float d = angle_between(v, o);
    if (d > radius)
      pass(v, color, age, alpha, tag);
    else {
      float t = d / radius;
      pass(v, color * quintic_kernel(t), age, alpha, tag);
    }
  }

private:
  OriginT origin;
  float radius;
};

/**
 * @brief Alias for Hole with reference storage.
 */
template <int W> using HoleRef = Hole<W, std::reference_wrapper<const Vector>>;

/**
 * @brief Replicates geometry by rotating it around the Y-axis.
 */
template <int W> class Replicate : public Is3D {
public:
  Replicate(int count)
      : count(hs::clamp(count, 1, W)),
        step(make_rotation(Y_AXIS, 2 * PI_F / count)) {}
  void plot(const Vector &v, const Pixel &color, float age, float alpha,
            uint8_t tag, PassFn3D pass) {
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

/**
 * @brief Applies a Mobius transformation to 3D points.
 */
template <int W> class Mobius : public Is3D {
public:
  Mobius(MobiusParams &params) : params(params) {}
  void plot(const Vector &v, const Pixel &color, float age, float alpha,
            uint8_t tag, PassFn3D pass) {
    pass(inv_stereo(mobius(stereo(v), params)), color, age, alpha, tag);
  }

private:
  MobiusParams &params;
};

/**
 * @brief Manages 3D world-space trails.
 */
template <int W, int Capacity> class Trails : public Is3DWithHistory {
public:
  struct Item {
    int16_t x, y, z; // Quantized unit vector (6 bytes)
    uint8_t ttl;     // Remaining lifetime in frames (1 byte)
    uint8_t pad_;    // Alignment (1 byte) — total 8 bytes
  };

  Trails(int lifetime) : lifetime(lifetime) {}

  /// Allocate ring buffer storage from persistent arena.
  /// Must be called from effect init(), not constructor (arenas aren't ready
  /// yet).
  void init_storage(Persistent arena) {
    items_ = static_cast<Item *>(
        arena.allocate(Capacity * sizeof(Item), alignof(Item)));
    head_ = tail_ = count_ = 0;
  }

  void plot(const Vector &v, const Pixel &color, float age, float alpha,
            uint8_t tag, PassFn3D pass) {
    pass(v, color, age, alpha, tag); // Pass through current frame

    int ttl = lifetime - static_cast<int>(age);
    if (ttl > 0 && items_) {
      push_back(encode(v, static_cast<uint8_t>(ttl)));
    }
  }

  void flush(const WorldTrailFn &trailFn, float alpha, PassFn3D pass) {
    // Age
    for (size_t i = 0; i < count_; ++i) {
      auto &item = at(i);
      if (item.ttl > 0)
        item.ttl--;
    }

    while (count_ > 0 && at(0).ttl == 0) {
      pop_front();
    }

    // Draw
    for (size_t i = 0; i < count_; ++i) {
      const auto &item = at(i);
      Vector v = decode(item);
      float t =
          1.0f - (static_cast<float>(item.ttl) / static_cast<float>(lifetime));
      Color4 c = trailFn(v, t);

      if (c.alpha > 0.001f) {
        pass(v, c.color, static_cast<float>(lifetime - item.ttl),
             c.alpha * alpha, 0);
      }
    }
  }

private:
  Item *items_ = nullptr;
  size_t head_ = 0, tail_ = 0, count_ = 0;
  int lifetime;

  static constexpr float Q = 32767.0f;
  static Item encode(const Vector &v, uint8_t ttl) {
    return {static_cast<int16_t>(v.x * Q), static_cast<int16_t>(v.y * Q),
            static_cast<int16_t>(v.z * Q), ttl, 0};
  }
  static Vector decode(const Item &item) {
    constexpr float INV_Q = 1.0f / Q;
    return Vector(item.x * INV_Q, item.y * INV_Q, item.z * INV_Q);
  }

  Item &at(size_t i) { return items_[(head_ + i) % Capacity]; }
  const Item &at(size_t i) const { return items_[(head_ + i) % Capacity]; }

  void push_back(const Item &item) {
    if (count_ == Capacity) {
      pop_front(); // Drop oldest on overflow
    }
    items_[tail_] = item;
    tail_ = (tail_ + 1) % Capacity;
    count_++;
  }

  void pop_front() {
    head_ = (head_ + 1) % Capacity;
    count_--;
  }
};

} // namespace World

// ----------------------------------------------------------------------------
// Screen Space Filters (2D)
// ----------------------------------------------------------------------------
namespace Screen {

/**
 * @brief Applies 2D anti-aliasing to sub-pixel coordinates.
 * Distributes intensity to 4 nearest neighbors using quintic kernel.
 */
template <int W, int H> class AntiAlias : public Is2D {
public:
  AntiAlias() {}
  void plot(float x, float y, const Pixel &c, float age, float alpha,
            uint8_t tag, PassFn2D pass) {
    float x_i, y_i;
    float x_m = std::modf(x, &x_i);
    float y_m = std::modf(y, &y_i);
    float xs = quintic_kernel(x_m);
    float ys = quintic_kernel(y_m);

    float v00 = (1 - xs) * (1 - ys);
    float v10 = xs * (1 - ys);
    float v01 = (1 - xs) * ys;
    float v11 = xs * ys;

    if (v00 > 1e-8)
      pass(x_i, y_i, c, age, alpha * v00, tag);
    if (v10 > 1e-8)
      pass(fast_wrap((x_i + 1), W), y_i, c, age, alpha * v10, tag);
    if (y_i < H - 1) {
      if (v01 > 1e-8)
        pass(x_i, y_i + 1, c, age, alpha * v01, tag);
      if (v11 > 1e-8)
        pass(fast_wrap((x_i + 1), W), y_i + 1, c, age, alpha * v11, tag);
    }
  }
};

/**
 * @brief Manages 2D screen-space trails.
 */
template <int W, int MAX_PIXELS = 1024> class Trails : public Is2DWithHistory {
public:
  Trails(int lifetime) : lifetime(lifetime) {}

  void init_storage(Persistent arena) {
    ttls_ = static_cast<DecayPixel *>(
        arena.allocate(MAX_PIXELS * sizeof(DecayPixel), alignof(DecayPixel)));
    num_pixels = 0;
  }

  void plot(float x, float y, const Pixel &color, float age, float alpha,
            uint8_t tag, PassFn2D pass) {
    if (alpha <= 0.001f)
      return;

    if (age >= 0) {
      if (ttls_ && num_pixels < MAX_PIXELS) {
        ttls_[num_pixels++] = {x, y, static_cast<float>(lifetime - age), tag};
      }
    }
    if (age <= 0) {
      pass(x, y, color, age, alpha, tag);
    }
  }

  void flush(Canvas &, const ScreenTrailFn &trailFn, float alpha, PassFn2D pass) {
    for (int i = 0; i < num_pixels; ++i) {
      Color4 color =
          trailFn(ttls_[i].x, ttls_[i].y, 1 - (ttls_[i].ttl / lifetime));
      if (color.alpha > 0.001f) {
        pass(ttls_[i].x, ttls_[i].y, color.color, lifetime - ttls_[i].ttl,
             alpha * color.alpha, ttls_[i].tag);
      }
    }
    decay();
  }

  void decay() {
    for (int i = 0; i < num_pixels; ++i) {
      if (--ttls_[i].ttl < TOLERANCE) {
        ttls_[i] = ttls_[--num_pixels];
        i--;
      }
    }
  }

private:
  struct DecayPixel {
    float x, y, ttl;
    uint8_t tag;
  };
  int lifetime;
  DecayPixel *ttls_ = nullptr;
  int num_pixels = 0;
};

/**
 * @brief Applies a variable 3x3 Gaussian Blur.
 */
template <int W, int H> class Blur : public Is2D {
public:
  Blur(float factor = 1.0f) { update(factor); }

  void update(float factor) {
    float f = hs::clamp(factor, 0.0f, 1.0f);
    // Gaussian reference: Corner=1/16, Edge=2/16, Center=4/16
    float c = 1.0f - (0.75f * f);
    float e = 0.125f * f;
    float d = 0.0625f * f;

    kernel = {d, e, d, e, c, e, d, e, d};
  }

  void plot(float x, float y, const ::Pixel &color, float age, float alpha,
            uint8_t tag, PassFn2D pass) {
    int cx = static_cast<int>(std::round(x));
    int cy = static_cast<int>(std::round(y));

    int k = 0;
    for (int dy = -1; dy <= 1; dy++) {
      int ny = cy + dy;

      if (ny >= 0 && ny < H) {
        for (int dx = -1; dx <= 1; dx++) {
          float weight = kernel[k++];
          if (weight > 1e-5f) {
            pass(static_cast<float>(wrap(cx + dx, W)), static_cast<float>(ny),
                 color, age, alpha * weight, tag);
          }
        }
      } else {
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
template <int W, int Capacity> class Slew : public Is2DWithHistory {
public:
  Slew(float rise = 1.0f, float fall = 0.05f) : rise(rise), fall(fall) {}

  void init_storage(Persistent arena) {
    items_ = static_cast<Item *>(
        arena.allocate(Capacity * sizeof(Item), alignof(Item)));
    head_ = tail_ = count_ = 0;
  }

  void plot(float x, float y, const ::Pixel &color, float age, float alpha,
            uint8_t tag, PassFn2D pass) {
    pass(x, y, color, age, alpha, tag);

    // Add to decay buffer
    if (items_ && count_ < Capacity && alpha > 0.001f) {
      items_[tail_] = {x, y, color, alpha, tag};
      tail_ = (tail_ + 1) % Capacity;
      count_++;
    }
  }

  void flush(Canvas &, const ScreenTrailFn &, float globalAlpha, PassFn2D pass) {
    size_t n = count_;
    for (size_t i = 0; i < n; ++i) {
      auto item = items_[head_];
      head_ = (head_ + 1) % Capacity;
      count_--;

      // Decay
      item.alpha -= fall;

      if (item.alpha > 0.00001f) {
        // Still alive — re-enqueue
        pass(item.x, item.y, item.c, 0, item.alpha * globalAlpha, item.tag);
        items_[tail_] = item;
        tail_ = (tail_ + 1) % Capacity;
        count_++;
      }
    }
  }

private:
  struct Item {
    float x, y;
    ::Pixel c;
    float alpha;
    uint8_t tag;
  };
  Item *items_ = nullptr;
  size_t head_ = 0, tail_ = 0, count_ = 0;
  float rise, fall;
};

} // namespace Screen

// ----------------------------------------------------------------------------
// Pixel Space Filters (1:1 with buffer)
// ----------------------------------------------------------------------------
namespace Pixel {

/**
 * @brief ::Pixel-space feedback filter (stateless).
 * Creates infinite trails and "watery" distortion by sampling the
 * previous frame from the Canvas front buffer with bilinear interpolation,
 * applying spatial distortion and fade, then blending into the back buffer.
 * Uses Canvas double-buffering — no internal frame storage needed.
 */
template <int W, int H, typename SpaceTransformFn>
class Feedback : public Is2DWithHistory {
public:
  Feedback(SpaceTransformFn transform_fn, float fade = 0.95f)
      : transform_fn(transform_fn), fade(fade) {}

  void set_fade(float f) { fade = f; }

  /// Pass-through: current-frame pixels go straight to the canvas sink.
  void plot(float x, float y, const ::Pixel &color, float age, float alpha,
            uint8_t tag, PassFn2D pass) {
    pass(x, y, color, age, alpha, tag);
  }

  /// Blend distorted previous frame (front buffer) into current frame (back buffer).
  void flush(Canvas &cv, const ScreenTrailFn &, float alpha, PassFn2D) {
    for (int y = 0; y < H; ++y) {
      for (int x = 0; x < W; ++x) {
        // Project 2D pixel to 3D sphere
        Vector v = pixel_to_vector<W, H>(x, y);

        // Apply 3D spatial transformation (e.g., NoiseTransformer)
        Vector v_dist = transform_fn(v);

        // Map back to continuous 2D coordinates
        Spherical s(v_dist);
        float bx = (s.theta * W) / (2.0f * PI_F);
        float by = phi_to_y<H>(s.phi);

        // Sample previous frame with bilinear interpolation and fade
        ::Pixel p = sample_bilinear_prev(cv, bx, by) * fade;

        // Blend into back buffer
        if (p.r | p.g | p.b) {
          int xi = fast_wrap(x, W);
          if (y >= 0 && y < cv.height()) {
            cv(xi, y) = blend_alpha(alpha)(cv(xi, y), p);
          }
        }
      }
    }
  }

private:
  /// Bilinear sample from the Canvas front buffer (previous frame).
  static ::Pixel sample_bilinear_prev(const Canvas &cv, float bx, float by) {
    float fx0 = std::floor(bx);
    float fy0 = std::floor(by);
    int x0 = static_cast<int>(fx0);
    int y0 = static_cast<int>(fy0);
    float fx = bx - fx0;
    float fy = by - fy0;

    int x1 = x0 + 1;
    int y1 = y0 + 1;

    // Safe wrapping for X (cylindrical)
    if (x0 < 0)
      x0 = W - 1 - ((-x0 - 1) % W);
    else if (x0 >= W)
      x0 %= W;

    if (x1 < 0)
      x1 = W - 1 - ((-x1 - 1) % W);
    else if (x1 >= W)
      x1 %= W;

    // Read from front buffer; black if out-of-bounds Y
    ::Pixel p00 = (y0 >= 0 && y0 < H) ? cv.prev(x0, y0) : ::Pixel(0, 0, 0);
    ::Pixel p10 = (y0 >= 0 && y0 < H) ? cv.prev(x1, y0) : ::Pixel(0, 0, 0);
    ::Pixel p01 = (y1 >= 0 && y1 < H) ? cv.prev(x0, y1) : ::Pixel(0, 0, 0);
    ::Pixel p11 = (y1 >= 0 && y1 < H) ? cv.prev(x1, y1) : ::Pixel(0, 0, 0);

    float w00 = (1.0f - fx) * (1.0f - fy);
    float w10 = fx * (1.0f - fy);
    float w01 = (1.0f - fx) * fy;
    float w11 = fx * fy;

    ::Pixel result = (p00 * w00);
    result += (p10 * w10);
    result += (p01 * w01);
    result += (p11 * w11);

    return result;
  }

  SpaceTransformFn transform_fn;
  float fade;
};

template <int W> class ChromaticShift : public Is2D {
public:
  ChromaticShift() {}

  void plot(float x, float y, const ::Pixel &c, float age, float alpha,
            uint8_t tag, PassFn2D pass) {
    pass(x, y, c, age, alpha, tag);

    ::Pixel r_col = c;
    r_col.g = 0;
    r_col.b = 0;
    ::Pixel g_col = c;
    g_col.r = 0;
    g_col.b = 0;
    ::Pixel b_col = c;
    b_col.r = 0;
    b_col.g = 0;

    pass(static_cast<float>(wrap(static_cast<int>(x) + 1, W)), y, r_col, age,
         alpha, tag);
    pass(static_cast<float>(wrap(static_cast<int>(x) + 2, W)), y, g_col, age,
         alpha, tag);
    pass(static_cast<float>(wrap(static_cast<int>(x) + 3, W)), y, b_col, age,
         alpha, tag);
  }
};

} // namespace Pixel

} // namespace Filter