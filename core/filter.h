/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the Polyform Noncommercial License 1.0.0
 */
#pragma once

#include <tuple>
#include <utility>
#include <type_traits>
#include <cmath>
#include <cassert>

#include <span>
#include <algorithm>
#include "geometry.h"
#include "color.h"
#include "static_circular_buffer.h"
#include "canvas.h"
#include "concepts.h"
#include "memory.h"
#include "styles.h"

using PassFn2D =
    FunctionRef<void(float, float, const Pixel &, float, float)>;
using PassFn3D =
    FunctionRef<void(const Vector &, const Pixel &, float, float)>;

// Filter Traits
//
// `is_terminal` marks a filter that writes the Canvas directly during flush()
// and ignores its downstream `pass` callback (e.g. Pixel::Feedback). Such a
// filter swallows the stream, so anything chained after it would silently never
// run — the Pipeline static-asserts a terminal filter is the last stage. Almost
// all filters forward via `pass` and leave this false.
/** @brief Trait indicating a filter operates in 2D screen space. */
struct Is2D {
  static constexpr bool is_2d = true;
  static constexpr bool has_history = false;
  static constexpr bool is_terminal = false;
};
/** @brief Trait indicating a filter operates in 3D world space. */
struct Is3D {
  static constexpr bool is_2d = false;
  static constexpr bool has_history = false;
  static constexpr bool is_terminal = false;
};

/** @brief Trait indicating a 2D filter that maintains state/history. */
struct Is2DWithHistory {
  static constexpr bool is_2d = true;
  static constexpr bool has_history = true;
  static constexpr bool is_terminal = false;
};

/** @brief Trait indicating a 3D filter that maintains state/history. */
struct Is3DWithHistory {
  static constexpr bool is_2d = false;
  static constexpr bool has_history = true;
  static constexpr bool is_terminal = false;
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

  // Type-safe filter accessor (base case: T not found). The guard is a
  // dependent-false: !sizeof(T*) is always false but cannot be proven so until
  // the template is instantiated, so it fires only when get<T>() is actually
  // named on a pipeline lacking T. (sizeof(T*) — not sizeof(T) — so the intended
  // "not found" diagnostic also wins for an incomplete T.)
  template <typename T> T &get() {
    static_assert(!sizeof(T *), "Filter type T not found in Pipeline");
  }
  template <typename T> const T &get() const {
    static_assert(!sizeof(T *), "Filter type T not found in Pipeline");
  }

  // 2D Sink
  void plot(Canvas &cv, int x, int y, const Pixel &c, float age, float alpha) {
    HS_PROFILE(filter_blend);
    // fast_wrap() only corrects a single +/-W offset, so the sink relies on the
    // producer keeping x within [-W, 2W). Every current producer does (3D path
    // full-wraps via vector_to_pixel; AntiAlias/Blur/ChromaticShift wrap before
    // forwarding). This precondition is enforced debug-only: it is stripped
    // under NDEBUG on device (zero hot-loop cost) but fires in the native test
    // suite / WASM-debug if a NaN coord or a future out-of-range filter slips in.
    assert(x >= -W && x < 2 * W);
    if (!cv.clip().contains_y(y)) return;
    int xi = fast_wrap(x, W);
    if (!cv.clip().contains_x(xi)) return;
    Pixel p = blend_alpha(alpha)(cv(xi, y), c);
    plot_virtual(cv, xi, y, p);
  }

  void plot(Canvas &cv, float x, float y, const Pixel &c, float age,
            float alpha) {
    // Non-finite coords would make the int casts below UB and bypass the wrap,
    // and round() must land within fast_wrap()'s [-W, 2W) window. Both hold for
    // every finite producer; assert debug-only (stripped on device, fires in
    // tests / WASM-debug) so a NaN-producing effect faults at the bench.
    assert(std::isfinite(x) && std::isfinite(y));
    int xi = static_cast<int>(std::round(x));
    int yi = static_cast<int>(std::round(y));
    assert(xi >= -W && xi < 2 * W);
    if (!cv.clip().contains_y(yi)) return;
    xi = fast_wrap(xi, W);
    if (!cv.clip().contains_x(xi)) return;
    Pixel p = blend_alpha(alpha)(cv(xi, yi), c);
    plot_virtual(cv, xi, yi, p);
  }

  // 3D Sink
  void plot(Canvas &cv, const Vector &v, const Pixel &c, float age, float alpha) {
    auto p = vector_to_pixel<W, H>(v);
    plot(cv, p.x, p.y, c, age, alpha);
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

  // Type-safe filter accessor
  template <typename T> T &get() {
    if constexpr (std::is_same_v<Head, T>) {
      return static_cast<T &>(*this);
    } else {
      return next.template get<T>();
    }
  }
  template <typename T> const T &get() const {
    if constexpr (std::is_same_v<Head, T>) {
      return static_cast<const T &>(*this);
    } else {
      return next.template get<T>();
    }
  }

  // 2D Plot (Float)
  void plot(Canvas &cv, float x, float y, const Pixel &c, float age,
            float alpha) {
    if constexpr (Head::is_2d) {
      Head::plot(x, y, c, age, alpha,
                 [&](float x, float y, const Pixel &c, float age, float alpha) {
                   next.plot(cv, x, y, c, age, alpha);
                 });
    } else { // 2D -> 3D Mismatch
      Vector v = pixel_to_vector<W, H>(x, y);
      plot(cv, v, c, age, alpha);
    }
  }

  // 2D Plot (Int)
  void plot(Canvas &cv, int x, int y, const Pixel &c, float age, float alpha) {
    plot(cv, static_cast<float>(x), static_cast<float>(y), c, age, alpha);
  }

  // 3D Plot
  void plot(Canvas &cv, const Vector &v, const Pixel &c, float age, float alpha) {
    if constexpr (!Head::is_2d) {
      Head::plot(v, c, age, alpha,
                 [&](const Vector &v, const Pixel &c, float age, float alpha) {
                   next.plot(cv, v, c, age, alpha);
                 });
    } else { // 3D -> 2D Mismatch
      auto p = vector_to_pixel<W, H>(v);
      plot(cv, p.x, p.y, c, age, alpha);
    }
  }

  // Flush contract: a history-bearing filter must define exactly the flush
  // signature matching its dimension trait — the if-constexpr dispatch below
  // only ever names the matching one (the other call sits in a discarded
  // branch, so a single-flush filter composes fine). These asserts make that
  // contract explicit: a trait/signature mismatch fails here with a readable
  // message instead of a deep overload-resolution error further down.
  static_assert(
      !Head::has_history || Head::is_2d ||
          requires(Head h, const WorldTrailFn &w, PassFn3D p) {
            h.flush(w, 1.0f, p);
          },
      "3D history filter must define "
      "flush(const WorldTrailFn&, float, PassFn3D)");
  static_assert(
      !Head::has_history || !Head::is_2d ||
          requires(Head h, Canvas &cv, const ScreenTrailFn &s, PassFn2D p) {
            h.flush(cv, s, 1.0f, p);
          },
      "2D history filter must define "
      "flush(Canvas&, const ScreenTrailFn&, float, PassFn2D)");

  // A terminal filter (one that writes the Canvas directly in flush() and
  // ignores its downstream pass, e.g. Pixel::Feedback) swallows the stream, so
  // any stage chained after it would silently never run. Enforce it is last.
  static_assert(
      !Head::is_terminal || sizeof...(Tail) == 0,
      "A terminal filter (e.g. Pixel::Feedback) writes the Canvas directly and "
      "ignores downstream filters — it must be the last stage in the Pipeline.");

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
            PassFn3D pass) {
    tween(orientation, [&](const Quaternion &q, float t) {
      pass(rotate(v, q), color, age + (1.0f - t), alpha);
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
            PassFn3D pass) {
    if (!enabled) {
      pass(v, color, age, alpha);
      return;
    }

    float projection =
        v.x * axis.x + v.y * axis.y + v.z * axis.z; // dot product
    float dot_val = std::max(-1.0f, std::min(1.0f, projection));
    float t = 1.0f - fast_acos(dot_val) / PI_F;

    size_t count = orientations.size();
    if (count == 0) {
      pass(v, color, age, alpha);
      return;
    }

    size_t idx = static_cast<size_t>(floorf(t * count));
    if (idx >= count)
      idx = count - 1;

    // Pass to selected orientation
    const Orientation<W> &q = orientations[idx];
    tween(q, [&](const Quaternion &rot, float tween_t) {
      pass(rotate(v, rot), color, age + (1.0f - tween_t), alpha);
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
            PassFn3D pass) {
    const Vector &o = origin;
    float d = angle_between(v, o);
    if (d > radius)
      pass(v, color, age, alpha);
    else {
      float t = d / radius;
      pass(v, color * quintic_kernel(t), age, alpha);
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
            PassFn3D pass) {
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

/**
 * @brief Replicates geometry onto the vertices of a solid.
 * Precomputes rotation quaternions from vertex[0] to each other vertex.
 * Each copy is emitted with age offset by its vertex index, so downstream
 * age-driven palettes can tint copies differently.
 */
template <int W, int N> class VertexReplicate : public Is3D {
public:
  /// Build from a vertex array. Computes rotations from vertices[0] → each.
  template <typename VertexArray>
  VertexReplicate(const VertexArray &vertices) {
    for (int i = 0; i < N; ++i)
      rotations[i] = make_rotation(vertices[0], vertices[i]);
  }

  void plot(const Vector &v, const Pixel &color, float age, float alpha,
            PassFn3D pass) {
    for (int i = 0; i < N; ++i) {
      pass(rotate(v, rotations[i]), color, age + static_cast<float>(i), alpha);
    }
  }

  std::array<Quaternion, N> rotations;
};

/**
 * @brief Applies a Mobius transformation to 3D points.
 */
template <int W> class Mobius : public Is3D {
public:
  Mobius(MobiusParams &params) : params(params) {}
  void plot(const Vector &v, const Pixel &color, float age, float alpha,
            PassFn3D pass) {
    pass(inv_stereo(mobius(stereo(v), params)), color, age, alpha);
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
  static_assert(sizeof(Item) == 8, "World::Trails::Item must be 8 bytes");

  // lifetime is a per-frame divisor (fade alpha); a zero/negative trail length
  // is a cold authoring error that would feed inf/NaN into blend_alpha. The
  // upper bound is structural: ttl is stored as uint8_t and encode() truncates,
  // so lifetime > 255 would silently wrap the trail length. Trap both at
  // construction (cold path) rather than producing garbage per pixel.
  Trails(int lifetime) : lifetime(lifetime) {
    HS_CHECK(lifetime > 0 && lifetime <= 255);
  }

  /// Allocate ring buffer storage from persistent arena.
  /// Must be called from effect init(), not constructor (arenas aren't ready
  /// yet).
  void init_storage(Arena &arena) {
    items_ = static_cast<Item *>(
        arena.allocate(Capacity * sizeof(Item), alignof(Item)));
    head_ = tail_ = count_ = 0;
  }

  void plot(const Vector &v, const Pixel &color, float age, float alpha,
            PassFn3D pass) {
    pass(v, color, age, alpha); // Pass through current frame

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
             c.alpha * alpha);
      }
    }
  }

  size_t size() const { return count_; }

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
  template <typename PassFnT>
  void plot(float x, float y, const Pixel &c, float age, float alpha,
            PassFnT &&pass) {
    float y_i;
    float y_m = std::modf(y, &y_i);

    // Floor-based integer/fraction split for X so a sub-pixel coordinate just
    // left of the theta=0 seam (x < 0, or x >= W) wraps to the correct columns
    // and keeps a fractional weight in [0, 1). std::modf truncates toward zero,
    // which collapses a negative fraction onto the wrong column and leaves the
    // left neighbor (x0) unwrapped, contrary to the sink's "wrap before
    // forwarding" contract.
    float x_floor = floorf(x);
    float x_m = x - x_floor; // always in [0, 1)

    // Spherical density compensation: scale X fractional by sin(phi).
    // At poles, all columns converge to the same point — snap to nearest.
    // At equator, sin(phi)=1 so behavior is unchanged.
    if (!TrigLUT<W, H>::initialized) TrigLUT<W, H>::init();
    int yi0 = hs::clamp(static_cast<int>(y_i), 0,
                        TrigLUT<W, H>::H_VIRT - 1);
    int yi1 = hs::clamp(static_cast<int>(y_i) + 1, 0,
                        TrigLUT<W, H>::H_VIRT - 1);
    float sin_phi = TrigLUT<W, H>::sin_phi[yi0]
        + (TrigLUT<W, H>::sin_phi[yi1] - TrigLUT<W, H>::sin_phi[yi0]) * y_m;
    float x_frac = hs::clamp(x_m * sin_phi, 0.0f, 1.0f);

    float xs = x_frac;  // sin(phi) already provides smooth pole collapse
    float ys = quintic_kernel(y_m);

    float v00 = (1 - xs) * (1 - ys);
    float v10 = xs * (1 - ys);
    float v01 = (1 - xs) * ys;
    float v11 = xs * ys;

    // Clamp Y neighbors to valid range (reflect at poles instead of dropping)
    int y0 = hs::clamp(static_cast<int>(y_i), 0, H - 1);
    int y1 = hs::clamp(static_cast<int>(y_i) + 1, 0, H - 1);
    int x0 = fast_wrap(static_cast<int>(x_floor), W);
    int x1 = fast_wrap(static_cast<int>(x_floor) + 1, W);

    if (v00 > 1e-8f)
      pass(static_cast<float>(x0), static_cast<float>(y0), c, age, alpha * v00);
    if (v10 > 1e-8f)
      pass(static_cast<float>(x1), static_cast<float>(y0), c, age, alpha * v10);
    if (v01 > 1e-8f)
      pass(static_cast<float>(x0), static_cast<float>(y1), c, age, alpha * v01);
    if (v11 > 1e-8f)
      pass(static_cast<float>(x1), static_cast<float>(y1), c, age, alpha * v11);
  }
};

/**
 * @brief Manages 2D screen-space trails.
 */
template <int W, int MAX_PIXELS = 1024> class Trails : public Is2DWithHistory {
public:
  // See World::Trails: lifetime divides per-frame; trap a zero/negative trail
  // length at construction (cold) rather than feeding inf/NaN into the fade.
  Trails(int lifetime) : lifetime(lifetime) { HS_CHECK(lifetime > 0); }

  void init_storage(Arena &arena) {
    ttls_ = static_cast<DecayPixel *>(
        arena.allocate(MAX_PIXELS * sizeof(DecayPixel), alignof(DecayPixel)));
    num_pixels = 0;
  }

  void plot(float x, float y, const Pixel &color, float age, float alpha,
            PassFn2D pass) {
    if (alpha <= 0.001f)
      return;

    if (age >= 0) {
      // Only seed live trail points. age >= lifetime means the point is already
      // dead; seeding it stores a negative ttl, which pushes t = 1-(ttl/lifetime)
      // above 1 (out of range for trailFn) and the pass-age above lifetime.
      // Mirror World::Trails' ttl>0 gate.
      float ttl = static_cast<float>(lifetime) - age;
      if (ttl > 0.0f && ttls_ && num_pixels < MAX_PIXELS) {
        ttls_[num_pixels++] = {x, y, ttl};
      }
    }
    if (age <= 0) {
      pass(x, y, color, age, alpha);
    }
  }

  void flush(Canvas &, const ScreenTrailFn &trailFn, float alpha,
             PassFn2D pass) {
    for (int i = 0; i < num_pixels; ++i) {
      Color4 color =
          trailFn(ttls_[i].x, ttls_[i].y, 1 - (ttls_[i].ttl / lifetime));
      if (color.alpha > 0.001f) {
        pass(ttls_[i].x, ttls_[i].y, color.color, lifetime - ttls_[i].ttl,
             alpha * color.alpha);
      }
    }
    decay();
  }

  void decay() {
    for (int i = 0; i < num_pixels; ++i) {
      if (--ttls_[i].ttl <= 0.0f) {
        ttls_[i] = ttls_[--num_pixels];
        i--;
      }
    }
  }

private:
  struct DecayPixel {
    float x, y, ttl;
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
            PassFn2D pass) {
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
                 color, age, alpha * weight);
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

} // namespace Screen

// ----------------------------------------------------------------------------
// Pixel Space Filters (1:1 with buffer)
// ----------------------------------------------------------------------------
namespace Pixel {

/**
 * @brief Style-aware feedback filter. Takes a `::Feedback::Style&` directly —
 * drop-in pipeline filter. The Style's spatial transform (noise, melt, etc.)
 * is smooth and expensive (noise + atan2 + acos per call), so the warp field
 * is computed on a coarse W/DS x H/DS grid (allocated from scratch_arena_a
 * per flush) and bilinearly upsampled. DS is read from style.downsample.
 *
 * Operates on the full canvas during flush() (integer coordinates, reads
 * cv.prev) — this is a pixel-space filter despite the 2D trait.
 *
 * TERMINAL: flush() composites directly into the Canvas and ignores its
 * `PassFn2D pass` callback, so it does NOT forward the stream downstream. It
 * must therefore be the last stage of the Pipeline — `is_terminal` makes the
 * Pipeline static-assert this (a non-terminal Feedback would silently drop
 * every filter chained after it).
 *
 * Usage:
 *   ::Feedback::Style style = ::Feedback::Style::Smoke();
 *   Pipeline<W, H, ..., Filter::Pixel::Feedback<W, H>> filters(
 *       ..., Filter::Pixel::Feedback<W, H>(style));
 */
template <int W, int H>
class Feedback : public Is2DWithHistory {
public:
  // Writes the Canvas directly in flush() and ignores `pass` — must be last.
  static constexpr bool is_terminal = true;

  explicit Feedback(::Feedback::Style &style) : style_(&style) {}

  /// Pass-through: current-frame pixels go straight to next filter.
  void plot(float x, float y, const ::Pixel &color, float age, float alpha,
            PassFn2D pass) {
    pass(x, y, color, age, alpha);
  }

  /// Blend distorted previous frame into current frame via the Style's transforms.
  void flush(Canvas &cv, const ScreenTrailFn &, float alpha, PassFn2D) {
    if (!enabled_) return;

    const int ds = style_->downsample;
    if (ds <= 0 || W % ds != 0 || H % ds != 0) return;
    const int hw = W / ds;
    const int hh = H / ds;

    scan_row_active(cv);
    if (is_frame_empty()) return;

    // Allocate coarse warp deltas (signed 8.8 fixed-point) from scratch.
    ScratchScope scope(scratch_arena_a);
    auto *dx = static_cast<int16_t *>(scope.get_arena().allocate(
        hh * hw * sizeof(int16_t), alignof(int16_t)));
    auto *dy = static_cast<int16_t *>(scope.get_arena().allocate(
        hh * hw * sizeof(int16_t), alignof(int16_t)));
    // No OOM guard here: Arena::allocate traps on over-allocation and never
    // returns null, so a "skip this frame" fallback would be dead code that
    // also implies a recovery path the fail-fast model deliberately rejects.

    // 1) Populate coarse warp field. Delta encoding keeps the bilerp safe
    //    across the longitude seam — neighbouring absolute bx values can
    //    straddle x=W ↔ x=0, but their deltas are continuous.
    constexpr float Q = 256.0f;
    for (int cy = 0; cy < hh; ++cy) {
      int y = cy * ds;
      for (int cx = 0; cx < hw; ++cx) {
        int x = cx * ds;
        Vector v_dist = style_->space_fn(pixel_to_vector<W, H>(x, y), *style_);
        Spherical s(v_dist);
        float bx = (s.theta * W) / (2.0f * PI_F);
        float by = phi_to_y<H>(s.phi);
        float ddx = bx - x;
        float ddy = by - y;
        if (ddx > W * 0.5f)       ddx -= W;
        else if (ddx < -W * 0.5f) ddx += W;
        dx[cy * hw + cx] = static_cast<int16_t>(
            hs::clamp(ddx * Q, -32767.0f, 32767.0f));
        dy[cy * hw + cx] = static_cast<int16_t>(
            hs::clamp(ddy * Q, -32767.0f, 32767.0f));
      }
    }

    // 2) Sample at full res, bilerping the coarse warp field per pixel.
    constexpr float INV_Q = 1.0f / Q;
    const float inv_ds = 1.0f / ds;
    const float fade = style_->fade;
    for (int y = 0; y < H; ++y) {
      int cy0 = y / ds;
      int cy1 = (cy0 + 1 < hh) ? cy0 + 1 : hh - 1;
      float fy = (y - cy0 * ds) * inv_ds;
      float wy0 = 1.0f - fy, wy1 = fy;

      for (int x = 0; x < W; ++x) {
        int cx0 = x / ds;
        int cx1 = (cx0 + 1) % hw;
        float fx = (x - cx0 * ds) * inv_ds;
        float wx0 = 1.0f - fx, wx1 = fx;

        int i00 = cy0 * hw + cx0;
        int i10 = cy0 * hw + cx1;
        int i01 = cy1 * hw + cx0;
        int i11 = cy1 * hw + cx1;

        float w00 = wx0 * wy0, w10 = wx1 * wy0;
        float w01 = wx0 * wy1, w11 = wx1 * wy1;

        float ddx = (dx[i00] * w00 + dx[i10] * w10
                   + dx[i01] * w01 + dx[i11] * w11) * INV_Q;
        float ddy = (dy[i00] * w00 + dy[i10] * w10
                   + dy[i01] * w01 + dy[i11] * w11) * INV_Q;

        ::Pixel sample = sample_bilinear_prev(cv, x + ddx, y + ddy);
        ::Pixel p = style_->color_fn(sample, fade, *style_);

        if (p.r | p.g | p.b) {
          int xi = fast_wrap(x, W);
          if (y >= 0 && y < cv.height()) {
            cv(xi, y) = blend_alpha(alpha)(cv(xi, y), p);
          }
        }
      }
    }
  }

  /// Enable/disable feedback (disabled = skip flush entirely).
  void set_enabled(bool e) { enabled_ = e; }

  /// Access the bound Style.
  ::Feedback::Style &style() { return *style_; }
  const ::Feedback::Style &style() const { return *style_; }

private:
  // --- Row-active bitmask for read-side optimizations ---
  static constexpr int ROW_BITMASK_SIZE = (H + 7) / 8;
  uint8_t row_active_[ROW_BITMASK_SIZE] = {};

  static bool is_row_active(const uint8_t *bitmask, int row) {
    return (bitmask[row >> 3] >> (row & 7)) & 1;
  }

  /// Scan the previous frame buffer and set one bit per row that has
  /// at least one non-black pixel.
  void scan_row_active(const Canvas &cv) {
    std::memset(row_active_, 0, ROW_BITMASK_SIZE);
    for (int y = 0; y < H; ++y) {
      for (int x = 0; x < W; ++x) {
        ::Pixel p = cv.prev(x, y);
        if (p.r | p.g | p.b) {
          row_active_[y >> 3] |= (1 << (y & 7));
          break;
        }
      }
    }
  }

  bool is_frame_empty() const {
    for (int i = 0; i < ROW_BITMASK_SIZE; ++i) {
      if (row_active_[i]) return false;
    }
    return true;
  }

  /// Bilinear sample from the Canvas front buffer (previous frame).
  /// Early-outs to black when both source rows are inactive.
  ::Pixel sample_bilinear_prev(const Canvas &cv, float bx, float by) const {
    float fy0 = std::floor(by);
    int y0 = static_cast<int>(fy0);
    int y1 = y0 + 1;

    bool r0_active = (y0 >= 0 && y0 < H) && is_row_active(row_active_, y0);
    bool r1_active = (y1 >= 0 && y1 < H) && is_row_active(row_active_, y1);
    if (!r0_active && !r1_active) return ::Pixel(0, 0, 0);

    float fx0 = std::floor(bx);
    int x0 = static_cast<int>(fx0);
    float fx = bx - fx0;
    float fy = by - fy0;
    int x1 = x0 + 1;

    if (x0 < 0)       x0 = W - 1 - ((-x0 - 1) % W);
    else if (x0 >= W) x0 %= W;
    if (x1 < 0)       x1 = W - 1 - ((-x1 - 1) % W);
    else if (x1 >= W) x1 %= W;

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

  ::Feedback::Style *style_;
  bool enabled_ = true;
};

/**
 * @brief Splits RGB into per-channel copies offset by 1/2/3 columns,
 * producing a chromatic-aberration fringe.
 */
template <int W> class ChromaticShift : public Is2D {
public:
  ChromaticShift() {}

  void plot(float x, float y, const ::Pixel &c, float age, float alpha,
            PassFn2D pass) {
    pass(x, y, c, age, alpha);

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
         alpha);
    pass(static_cast<float>(wrap(static_cast<int>(x) + 2, W)), y, g_col, age,
         alpha);
    pass(static_cast<float>(wrap(static_cast<int>(x) + 3, W)), y, b_col, age,
         alpha);
  }
};

} // namespace Pixel

} // namespace Filter

