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

/** @brief Callback that forwards a 2D plot (x, y, pixel, age, alpha) downstream. */
using PassFn2D =
    FunctionRef<void(float, float, const Pixel &, float, float)>;
/** @brief Callback that forwards a 3D plot (vector, pixel, age, alpha) downstream. */
using PassFn3D =
    FunctionRef<void(const Vector &, const Pixel &, float, float)>;

/**
 * @brief Trait indicating a filter operates in 2D screen space.
 * @details `is_terminal` marks a filter that writes the Canvas directly during
 * flush() and ignores its downstream `pass` callback (e.g. Pixel::Feedback).
 * Such a filter swallows the stream, so anything chained after it would
 * silently never run — the Pipeline static-asserts a terminal filter is the
 * last stage. Almost all filters forward via `pass` and leave this false.
 */
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
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @tparam Filters Ordered list of filter stages to chain.
 * @details Chains filters together, routing 2D and 3D plot commands.
 */
template <int W, int H, typename... Filters> struct Pipeline;

/**
 * @brief Plots a color to the canvas, clipping the y-coordinate to bounds.
 * @param canvas Target canvas to write into.
 * @param x Column index (assumed already wrapped into range).
 * @param y Row index; the write is skipped unless y is in [0, height).
 * @param c Color to store at (x, y).
 */
inline void plot_virtual(Canvas &canvas, int x, int y, const Pixel &c) {
  if (y >= 0 && y < canvas.height()) {
    canvas(x, y) = c;
  }
}

/**
 * @brief Terminal node of the pipeline (base case). Writes final pixels to the
 * Canvas.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H> struct Pipeline<W, H> {
  static constexpr bool is_2d = true;

  /**
   * @brief Type-safe filter accessor (base case: T not found).
   * @tparam T Filter type being looked up.
   * @return Never returns; instantiation is a hard error.
   * @details The guard is a dependent-false: !sizeof(T*) is always false but
   * cannot be proven so until the template is instantiated, so it fires only
   * when get<T>() is actually named on a pipeline lacking T. (sizeof(T*) — not
   * sizeof(T) — so the intended "not found" diagnostic also wins for an
   * incomplete T.)
   */
  template <typename T> T &get() {
    static_assert(!sizeof(T *), "Filter type T not found in Pipeline");
  }
  /**
   * @brief Type-safe const filter accessor (base case: T not found).
   * @tparam T Filter type being looked up.
   * @return Never returns; instantiation is a hard error.
   */
  template <typename T> const T &get() const {
    static_assert(!sizeof(T *), "Filter type T not found in Pipeline");
  }

  /**
   * @brief Writes an integer-coordinate 2D sample to the canvas (sink).
   * @param cv Target canvas.
   * @param x Column in [-W, 2W); wrapped into [0, W) before writing.
   * @param y Row index in pixels.
   * @param c Source color to blend in.
   * @param alpha Blend alpha in [0, 1].
   * @details The unnamed float parameter is the unused age channel.
   */
  void plot(Canvas &cv, int x, int y, const Pixel &c, float, float alpha) {
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

  /**
   * @brief Writes a float-coordinate 2D sample to the canvas (sink).
   * @param cv Target canvas.
   * @param x Column; rounded then required to land in [-W, 2W) and wrapped.
   * @param y Row; rounded to nearest pixel.
   * @param c Source color to blend in.
   * @param alpha Blend alpha in [0, 1].
   * @details The unnamed float parameter is the unused age channel.
   */
  void plot(Canvas &cv, float x, float y, const Pixel &c, float,
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

  /**
   * @brief Projects a 3D point to screen space and writes it (sink).
   * @param cv Target canvas.
   * @param v World-space point on the unit sphere.
   * @param c Source color to blend in.
   * @param age Temporal age channel (frames).
   * @param alpha Blend alpha in [0, 1].
   */
  void plot(Canvas &cv, const Vector &v, const Pixel &c, float age, float alpha) {
    auto p = vector_to_pixel<W, H>(v);
    plot(cv, p.x, p.y, c, age, alpha);
  }

  /**
   * @brief Screen-trail flush no-op (sink has no history).
   * @details Unused Canvas, ScreenTrailFn and alpha parameters.
   */
  void flush(Canvas &, const ScreenTrailFn &, float) {}
  /**
   * @brief World-trail flush no-op (sink has no history).
   * @details Unused Canvas, WorldTrailFn and alpha parameters.
   */
  void flush(Canvas &, const WorldTrailFn &, float) {}
};

/**
 * @brief Recursive pipeline case: applies Head, then forwards to the Tail.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @tparam Head Filter stage applied at this level.
 * @tparam Tail Remaining filter stages.
 */
template <int W, int H, typename Head, typename... Tail>
struct Pipeline<W, H, Head, Tail...> : public Head {
  using Next = Pipeline<W, H, Tail...>;
  Next next;

  /**
   * @brief Forwarding-reference constructor: builds Head and the Tail in place.
   * @tparam HArg Argument type forwarded to Head's constructor.
   * @tparam TArgs Argument types forwarded to the remaining stages.
   * @param h Argument forwarded to Head's constructor.
   * @param t Arguments forwarded to the Tail pipeline's constructors.
   * @details The requires guard stops it from hijacking copy construction:
   * copying from a non-const Pipeline lvalue would otherwise prefer this
   * template (an exact Pipeline& match) over the implicit copy ctor (which adds
   * const) and try to construct Head(Pipeline&). Excluding Pipeline lets the
   * copy ctor win.
   */
  template <typename HArg, typename... TArgs>
    requires(!std::is_same_v<std::remove_cvref_t<HArg>, Pipeline>)
  Pipeline(HArg &&h, TArgs &&...t)
      : Head(std::forward<HArg>(h)), next(std::forward<TArgs>(t)...) {}

  /**
   * @brief Partial constructor: builds Head only, default-constructing the Tail.
   * @tparam HArg Argument type forwarded to Head's constructor.
   * @param h Argument forwarded to Head's constructor.
   * @details Same guard as the variadic ctor: as an explicit single-arg
   * forwarding ctor it would hijack direct copy construction the same way.
   */
  template <typename HArg>
    requires(!std::is_same_v<std::remove_cvref_t<HArg>, Pipeline>)
  explicit Pipeline(HArg &&h) : Head(std::forward<HArg>(h)) {}

  /** @brief Default-constructs every stage in the pipeline. */
  Pipeline() = default;

  /**
   * @brief Type-safe filter accessor: finds the stage of type T.
   * @tparam T Filter type to retrieve.
   * @return Reference to the stage of type T (recurses into the Tail if Head is
   * not T).
   */
  template <typename T> T &get() {
    if constexpr (std::is_same_v<Head, T>) {
      return static_cast<T &>(*this);
    } else {
      return next.template get<T>();
    }
  }
  /**
   * @brief Type-safe const filter accessor: finds the stage of type T.
   * @tparam T Filter type to retrieve.
   * @return Const reference to the stage of type T.
   */
  template <typename T> const T &get() const {
    if constexpr (std::is_same_v<Head, T>) {
      return static_cast<const T &>(*this);
    } else {
      return next.template get<T>();
    }
  }

  /**
   * @brief Routes a float-coordinate 2D plot through Head, else converts to 3D.
   * @param cv Target canvas.
   * @param x Column coordinate in pixels.
   * @param y Row coordinate in pixels.
   * @param c Source color.
   * @param age Temporal age channel (frames).
   * @param alpha Blend alpha in [0, 1].
   * @details If Head is 2D it processes directly; otherwise the point is lifted
   * to a world vector and dispatched to the 3D path.
   */
  void plot(Canvas &cv, float x, float y, const Pixel &c, float age,
            float alpha) {
    if constexpr (Head::is_2d) {
      Head::plot(x, y, c, age, alpha,
                 [&](float x, float y, const Pixel &c, float age, float alpha) {
                   next.plot(cv, x, y, c, age, alpha);
                 });
    } else { // 2D -> 3D mismatch: lift to world space
      Vector v = pixel_to_vector<W, H>(x, y);
      plot(cv, v, c, age, alpha);
    }
  }

  /**
   * @brief Integer-coordinate 2D plot overload; forwards to the float path.
   * @param cv Target canvas.
   * @param x Column coordinate in pixels.
   * @param y Row coordinate in pixels.
   * @param c Source color.
   * @param age Temporal age channel (frames).
   * @param alpha Blend alpha in [0, 1].
   */
  void plot(Canvas &cv, int x, int y, const Pixel &c, float age, float alpha) {
    plot(cv, static_cast<float>(x), static_cast<float>(y), c, age, alpha);
  }

  /**
   * @brief Routes a 3D plot through Head, else projects to 2D.
   * @param cv Target canvas.
   * @param v World-space point on the unit sphere.
   * @param c Source color.
   * @param age Temporal age channel (frames).
   * @param alpha Blend alpha in [0, 1].
   * @details If Head is 3D it processes directly; otherwise the point is
   * projected to screen space and dispatched to the 2D path.
   */
  void plot(Canvas &cv, const Vector &v, const Pixel &c, float age, float alpha) {
    if constexpr (!Head::is_2d) {
      Head::plot(v, c, age, alpha,
                 [&](const Vector &v, const Pixel &c, float age, float alpha) {
                   next.plot(cv, v, c, age, alpha);
                 });
    } else { // 3D -> 2D mismatch: project to screen space
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

  /**
   * @brief Flushes 2D history for this stage, then recurses into the Tail.
   * @param cv Target canvas.
   * @param trailFn Callback producing trail color/alpha per screen point.
   * @param alpha Global blend alpha in [0, 1].
   * @details Only a 2D history-bearing Head emits; other Heads pass through.
   */
  void flush(Canvas &cv, const ScreenTrailFn &trailFn, float alpha) {
    if constexpr (Head::has_history) {
      if constexpr (Head::is_2d) {
        Head::flush(cv, trailFn, alpha,
                    [&](auto... args) { next.plot(cv, args...); });
      }
    }
    next.flush(cv, trailFn, alpha);
  }

  /**
   * @brief Flushes 3D history for this stage, then recurses into the Tail.
   * @param cv Target canvas.
   * @param trailFn Callback producing trail color/alpha per world point.
   * @param alpha Global blend alpha in [0, 1].
   * @details Only a 3D history-bearing Head emits; other Heads pass through.
   */
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

// ----------------------------------------------------------------------------
// Namespace: Filter
// ----------------------------------------------------------------------------

namespace Filter {

// ----------------------------------------------------------------------------
// World Space Filters (3D)
// ----------------------------------------------------------------------------
namespace World {

/**
 * @brief Rotates 3D points based on a dynamic Orientation.
 *
 * Age-channel contract: this filter sweeps the orientation's intra-frame
 * SLERP history (`tween` walks t: 0→1, oldest→newest sub-position) and offsets
 * `age` by the *fractional* `(1 - t)`, spreading the sweep across one frame's
 * worth of age — i.e. genuine temporal motion blur, where the trailing
 * sub-positions fade as if one frame older. This is the only filter that
 * adjusts age: the Replicate filters pass it through untouched, since their
 * copies are spatial, not temporal (see their notes).
 */
template <int W> class Orient : public Is3D {
public:
  /**
   * @brief Binds the filter to a live orientation source.
   * @param orientation Orientation whose SLERP history drives the rotation.
   */
  Orient(Orientation<> &orientation) : orientation(orientation) {}

  /**
   * @brief Rotates and re-emits the point across the orientation's tween sweep.
   * @param v World-space point to rotate.
   * @param color Source color, forwarded unchanged.
   * @param age Incoming age (frames); offset by the fractional (1 - t) per tween step.
   * @param alpha Blend alpha in [0, 1], forwarded unchanged.
   * @param pass Downstream 3D callback.
   */
  void plot(const Vector &v, const Pixel &color, float age, float alpha,
            PassFn3D pass) {
    tween(orientation, [&](const Quaternion &q, float t) {
      pass(rotate(v, q), color, age + (1.0f - t), alpha);
    });
  }

private:
  Orientation<> &orientation; /**< Live orientation source driving the rotation. */
};

/**
 * @brief Selects an orientation from a list based on the point's projection
 * onto an axis. Useful for slicing objects with different rotations.
 */
template <int W> class OrientSlice : public Is3D {
public:
  /**
   * @brief Binds the slice selector to an orientation list and a slicing axis.
   * @param orientations Candidate orientations, indexed by axis projection.
   * @param axis Unit axis the point is projected onto to pick an orientation.
   */
  OrientSlice(std::span<const Orientation<>> orientations, const Vector &axis)
      : enabled(true), axis(axis), orientations(orientations) {}

  /**
   * @brief Selects an orientation by axis projection and re-emits the point.
   * @param v World-space point to rotate.
   * @param color Source color, forwarded unchanged.
   * @param age Incoming age (frames); offset by fractional (1 - t) per tween step.
   * @param alpha Blend alpha in [0, 1], forwarded unchanged.
   * @param pass Downstream 3D callback.
   * @details Passes through untouched when disabled or the orientation list is empty.
   */
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
    const Orientation<> &q = orientations[idx];
    tween(q, [&](const Quaternion &rot, float tween_t) {
      pass(rotate(v, rot), color, age + (1.0f - tween_t), alpha);
    });
  }

  bool enabled;  /**< When false, the filter passes points through unrotated. */
  Vector axis;   /**< Unit axis points are projected onto to select a slice. */

private:
  std::span<const Orientation<>> orientations; /**< Candidate orientations indexed by projection. */
};

/**
 * @brief Creates a spherical hole by masking points within a radius.
 * @tparam W Canvas width in pixels.
 * @tparam OriginT Storage type for the hole center: by value (Vector) or by
 * reference (std::reference_wrapper).
 */
template <int W, typename OriginT = Vector> class Hole : public Is3D {
public:
  /**
   * @brief Constructs a hole centered at @p origin with angular @p radius.
   * @param origin Center of the hole (unit vector).
   * @param radius Angular radius of the hole in radians.
   */
  Hole(OriginT origin, float radius) : origin(origin), radius(radius) {}
  /**
   * @brief Attenuates points near the hole center, leaving others unchanged.
   * @param v World-space point to test.
   * @param color Source color; scaled by a quintic falloff inside the radius.
   * @param age Temporal age channel (frames), forwarded unchanged.
   * @param alpha Blend alpha in [0, 1], forwarded unchanged.
   * @param pass Downstream 3D callback.
   */
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
  OriginT origin; /**< Center of the hole (unit vector), stored by value or ref. */
  float radius;   /**< Angular radius of the hole in radians. */
};

/**
 * @brief Alias for Hole with reference (std::reference_wrapper) center storage.
 * @tparam W Canvas width in pixels.
 */
template <int W> using HoleRef = Hole<W, std::reference_wrapper<const Vector>>;

/**
 * @brief Replicates geometry by rotating it around the Y-axis.
 *
 * Age-channel contract: every copy is emitted with the *same* source `age`.
 * The replicas are simultaneous spatial mirrors of a single instant, so they
 * must share one age — there is deliberately no per-copy offset here, matching
 * `VertexReplicate`; only `Orient` adjusts age, via its motion-blur tween.
 */
template <int W> class Replicate : public Is3D {
public:
  /**
   * @brief Builds a replicator emitting @p count evenly-spaced Y-axis copies.
   * @param count Desired copy count; clamped to [1, W].
   * @details Uses the clamped member, not the parameter: unqualified `count`
   * here resolves to the ctor parameter, so a raw count > W would span only a
   * fraction of the circle and count == 0 would feed inf into make_rotation.
   * `this->count` is initialized first (declared first).
   */
  Replicate(int count)
      : count(hs::clamp(count, 1, W)),
        step(make_rotation(Y_AXIS, 2 * PI_F / this->count)) {}
  /**
   * @brief Emits the point plus count-1 rotated copies around the Y axis.
   * @param v World-space point to replicate.
   * @param color Source color, forwarded unchanged to every copy.
   * @param age Temporal age channel (frames), shared by every copy.
   * @param alpha Blend alpha in [0, 1], forwarded unchanged.
   * @param pass Downstream 3D callback.
   */
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
  int count;        /**< Number of copies emitted, in [1, W]. */
  Quaternion step;  /**< Per-copy Y-axis rotation (2*pi / count). */
};

/**
 * @brief Replicates geometry onto the vertices of a solid.
 * Precomputes rotation quaternions from vertex[0] to each other vertex.
 * Every copy carries the source point's age unchanged: age is the temporal
 * channel (motion blur, trail TTL), not a per-copy discriminator — replication
 * is spatial and must not steer age-driven palettes. This matches Replicate;
 * only Orient/OrientSlice adjust age, via their fractional motion-blur tween.
 */
template <int W, int N> class VertexReplicate : public Is3D {
public:
  /**
   * @brief Builds from a vertex array, precomputing rotations vertices[0] → each.
   * @tparam VertexArray Indexable container of N unit vectors.
   * @param vertices Vertex positions; rotations map vertices[0] onto each vertex.
   */
  template <typename VertexArray>
  VertexReplicate(const VertexArray &vertices) {
    for (int i = 0; i < N; ++i)
      rotations[i] = make_rotation(vertices[0], vertices[i]);
  }

  /**
   * @brief Emits one rotated copy of the point per stored vertex rotation.
   * @param v World-space point to replicate.
   * @param color Source color, forwarded unchanged to every copy.
   * @param age Temporal age channel (frames), shared by every copy.
   * @param alpha Blend alpha in [0, 1], forwarded unchanged.
   * @param pass Downstream 3D callback.
   */
  void plot(const Vector &v, const Pixel &color, float age, float alpha,
            PassFn3D pass) {
    for (int i = 0; i < N; ++i) {
      pass(rotate(v, rotations[i]), color, age, alpha);
    }
  }

  std::array<Quaternion, N> rotations; /**< Rotation from vertices[0] to each vertex. */
};

/**
 * @brief Applies a Mobius transformation to 3D points.
 */
template <int W> class Mobius : public Is3D {
public:
  /**
   * @brief Binds the filter to a live Mobius parameter set.
   * @param params Mobius transform parameters applied per point.
   */
  Mobius(MobiusParams &params) : params(params) {}
  /**
   * @brief Stereographically projects, applies the Mobius map, and re-emits.
   * @param v World-space point on the unit sphere.
   * @param color Source color, forwarded unchanged.
   * @param age Temporal age channel (frames), forwarded unchanged.
   * @param alpha Blend alpha in [0, 1], forwarded unchanged.
   * @param pass Downstream 3D callback.
   */
  void plot(const Vector &v, const Pixel &color, float age, float alpha,
            PassFn3D pass) {
    pass(inv_stereo(mobius(stereo(v), params)), color, age, alpha);
  }

private:
  MobiusParams &params; /**< Live Mobius transform parameters. */
};

/**
 * @brief Manages 3D world-space trails.
 */
template <int W, int Capacity> class Trails : public Is3DWithHistory {
public:
  /** @brief One quantized trail sample: unit vector plus remaining lifetime. */
  struct Item {
    int16_t x, y, z; /**< Quantized unit vector components (6 bytes). */
    uint8_t ttl;     /**< Remaining lifetime in frames (1 byte). */
    uint8_t pad_;    /**< Padding for 8-byte alignment (1 byte). */
  };
  static_assert(sizeof(Item) == 8, "World::Trails::Item must be 8 bytes");

  /**
   * @brief Constructs a world trail buffer with the given fade lifetime.
   * @param lifetime Per-frame fade divisor in frames; must be in [1, 255].
   * @details lifetime is a per-frame divisor (fade alpha); a zero/negative
   * trail length is a cold authoring error that would feed inf/NaN into
   * blend_alpha. The upper bound is structural: ttl is stored as uint8_t and
   * encode() truncates, so lifetime > 255 would silently wrap the trail length.
   * Trap both at construction (cold path) rather than producing garbage per pixel.
   */
  Trails(int lifetime) : lifetime(lifetime) {
    HS_CHECK(lifetime > 0 && lifetime <= 255);
  }

  /**
   * @brief Retunes the trail length at runtime (e.g. from a "Trail Len" slider).
   * @param new_lifetime New fade divisor in frames; must be in [1, 255].
   * @details Same bounds as the constructor; a per-frame caller is expected to
   * clamp into [1,255] first, so this trap fires only on a genuine authoring
   * error, not on slider motion. Already-buffered points keep their encoded ttl
   * and age out under the new length within a few frames.
   */
  void set_lifetime(int new_lifetime) {
    HS_CHECK(new_lifetime > 0 && new_lifetime <= 255);
    lifetime = new_lifetime;
  }

  /**
   * @brief Allocates ring-buffer storage from the persistent arena.
   * @param arena Persistent arena supplying Capacity Item slots.
   * @details Must be called from effect init(), not the constructor (arenas
   * aren't ready yet).
   */
  void init_storage(Arena &arena) {
    items_ = static_cast<Item *>(
        arena.allocate(Capacity * sizeof(Item), alignof(Item)));
    head_ = tail_ = count_ = 0;
  }

  /**
   * @brief Forwards the current sample and seeds a fading trail point.
   * @param v World-space point on the unit sphere.
   * @param color Source color, forwarded unchanged this frame.
   * @param age Incoming age (frames); ttl = lifetime - age, seeded only if positive.
   * @param alpha Blend alpha in [0, 1], forwarded unchanged.
   * @param pass Downstream 3D callback.
   */
  void plot(const Vector &v, const Pixel &color, float age, float alpha,
            PassFn3D pass) {
    pass(v, color, age, alpha); // Pass through current frame

    int ttl = lifetime - static_cast<int>(age);
    if (ttl > 0 && items_) {
      push_back(encode(v, static_cast<uint8_t>(ttl)));
    }
  }

  /**
   * @brief Ages every buffered point one frame and re-emits the survivors.
   * @param trailFn Callback producing trail color/alpha from (point, t).
   * @param alpha Global blend alpha in [0, 1].
   * @param pass Downstream 3D callback.
   */
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
      // Heterogeneous TTLs let an item expire behind a younger one; the pop loop
      // above only frees front-contiguous dead items, so an expired mid-buffer
      // item lingers until it reaches the head. Skip it rather than redraw a
      // frozen ghost at t = 1.0 every frame (it is still popped once it reaches
      // the front).
      if (item.ttl == 0)
        continue;
      Vector v = decode(item);
      float t =
          1.0f - (static_cast<float>(item.ttl) / static_cast<float>(lifetime));
      Color4 c = trailFn(v, t);

      if (c.alpha > 0.001f) {
        // lifetime can shrink mid-run via set_lifetime() while older points
        // still carry a ttl encoded under the previous (larger) lifetime, which
        // makes (lifetime - ttl) negative; clamp to honor the age >= 0 contract.
        int age = lifetime - static_cast<int>(item.ttl);
        if (age < 0)
          age = 0;
        pass(v, c.color, static_cast<float>(age), c.alpha * alpha);
      }
    }
  }

  /**
   * @brief Returns the number of live trail points currently buffered.
   * @return Count of buffered Item entries.
   */
  size_t size() const { return count_; }

private:
  Item *items_ = nullptr;                 /**< Ring-buffer storage (arena-owned). */
  size_t head_ = 0, tail_ = 0, count_ = 0; /**< Ring-buffer head, tail, and live count. */
  int lifetime;                           /**< Per-frame fade divisor in frames. */

  static constexpr float Q = 32767.0f;    /**< Quantization scale for unit-vector components. */
  /**
   * @brief Encodes a unit vector and ttl into a quantized Item.
   * @param v World-space point; each component is clamped to [-1, 1] before scaling.
   * @param ttl Remaining lifetime in frames.
   * @return Packed Item with quantized coordinates.
   */
  static Item encode(const Vector &v, uint8_t ttl) {
    // Saturate each component to the unit cube before quantizing. An upstream
    // World filter (Mobius warp, ripple) can transiently push a component past
    // 1; without the clamp v.x*Q overflows int16 and wraps to a garbage point
    // that lands on the opposite side of the sphere. Clamping pins the stray
    // sample to the surface instead.
    return {static_cast<int16_t>(hs::clamp(v.x, -1.0f, 1.0f) * Q),
            static_cast<int16_t>(hs::clamp(v.y, -1.0f, 1.0f) * Q),
            static_cast<int16_t>(hs::clamp(v.z, -1.0f, 1.0f) * Q), ttl, 0};
  }
  /**
   * @brief Decodes a quantized Item back into a unit vector.
   * @param item Packed trail sample.
   * @return Reconstructed world-space point.
   */
  static Vector decode(const Item &item) {
    constexpr float INV_Q = 1.0f / Q;
    return Vector(item.x * INV_Q, item.y * INV_Q, item.z * INV_Q);
  }

  /**
   * @brief Returns the i-th live item by age (0 = oldest).
   * @param i Index into the live range [0, count_).
   * @return Mutable reference to the buffered Item.
   */
  Item &at(size_t i) { return items_[(head_ + i) % Capacity]; }
  /**
   * @brief Returns the i-th live item by age (0 = oldest).
   * @param i Index into the live range [0, count_).
   * @return Const reference to the buffered Item.
   */
  const Item &at(size_t i) const { return items_[(head_ + i) % Capacity]; }

  /**
   * @brief Appends an item, evicting the oldest when at capacity.
   * @param item Encoded trail sample to push.
   */
  void push_back(const Item &item) {
    if (count_ == Capacity) {
      pop_front(); // Drop oldest on overflow
    }
    items_[tail_] = item;
    tail_ = (tail_ + 1) % Capacity;
    count_++;
  }

  /** @brief Drops the oldest buffered item. */
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
  // Cold one-time touch of the phi LUT so the per-sub-pixel plot() path can read
  // it unconditionally. In production init_geometry_luts() has already filled the
  // table; this guarded call is the lazy fallback for callers that skip engine
  // setup. Hoisting it here keeps the hot loop free of the init branch.
  AntiAlias() {
    if (!TrigLUT<W, H>::initialized) TrigLUT<W, H>::init();
  }
  /**
   * @brief Splats a sub-pixel sample across its four nearest pixel neighbors.
   * @tparam PassFnT Downstream 2D callback type.
   * @param x Sub-pixel column coordinate.
   * @param y Sub-pixel row coordinate.
   * @param c Source color, forwarded to each tap.
   * @param age Temporal age channel (frames), forwarded unchanged.
   * @param alpha Blend alpha in [0, 1]; scaled per tap by its bilinear weight.
   * @param pass Downstream 2D callback receiving each weighted tap.
   * @details X fractional weight is scaled by sin(phi) for spherical density
   * compensation; both axes are eased with a quintic kernel.
   */
  template <typename PassFnT>
  void plot(float x, float y, const Pixel &c, float age, float alpha,
            PassFnT &&pass) {
    // age/alpha are non-negative by contract; a negative alpha would silently
    // drop all four taps (alpha*weight < 1e-8f). Debug-only trap — stripped on
    // device, fires in the native tests / WASM-debug.
    assert(age >= 0.0f && alpha >= 0.0f);
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
    // At poles, all columns converge to the same point, so the compression
    // drives x_frac->0 and the sample collapses onto its left (floor) column
    // x0 — a floor-snap, NOT a round-to-nearest (x_m=0.9 still lands on x0).
    // At equator, sin(phi)=1 so behavior is unchanged.
    // (LUT init hoisted to the constructor — see above.)
    int yi0 = hs::clamp(static_cast<int>(y_i), 0,
                        TrigLUT<W, H>::H_VIRT - 1);
    int yi1 = hs::clamp(static_cast<int>(y_i) + 1, 0,
                        TrigLUT<W, H>::H_VIRT - 1);
    float sin_phi = TrigLUT<W, H>::sin_phi[yi0]
        + (TrigLUT<W, H>::sin_phi[yi1] - TrigLUT<W, H>::sin_phi[yi0]) * y_m;
    float x_frac = hs::clamp(x_m * sin_phi, 0.0f, 1.0f);

    // Ease both axes with the same quintic kernel so the bilinear weight has
    // vanishing 1st/2nd derivatives at the cell edges — without it, a point
    // drifting across a column boundary ramps linearly in X while Y stays
    // smooth (most visible at the equator, where sin(phi)=1 and x_frac==x_m).
    // sin(phi) is an orthogonal density compensation, not a substitute for the
    // easing; at the poles x_frac->0 so quintic_kernel(0)=0 collapses the
    // sample onto its left (floor) column x0.
    float xs = quintic_kernel(x_frac);
    float ys = quintic_kernel(y_m);

    float v00 = (1 - xs) * (1 - ys);
    float v10 = xs * (1 - ys);
    float v01 = (1 - xs) * ys;
    float v11 = xs * ys;

    // Clip Y neighbors to the physical row range rather than clamping them. With
    // H_OFFSET > 0 the LEDs stop short of the south pole, so a virtual sub-pole
    // row (y >= H) must be dropped — clamping it onto row H-1 would stretch the
    // bottom of the image onto the last LED row at full weight instead of
    // clipping it. The sink applies the same contains_y() clip downstream; doing
    // it here stops a clamped tap from defeating it. (Host builds set
    // H_OFFSET = 0, so no out-of-range row is produced and behavior is unchanged
    // — which is why this device-only divergence is invisible to the host suite.)
    int y0 = static_cast<int>(y_i);
    int y1 = y0 + 1;
    int x0 = fast_wrap(static_cast<int>(x_floor), W);
    int x1 = fast_wrap(static_cast<int>(x_floor) + 1, W);

    bool y0_ok = y0 >= 0 && y0 < H;
    bool y1_ok = y1 >= 0 && y1 < H;

    if (y0_ok && v00 > 1e-8f)
      pass(static_cast<float>(x0), static_cast<float>(y0), c, age, alpha * v00);
    if (y0_ok && v10 > 1e-8f)
      pass(static_cast<float>(x1), static_cast<float>(y0), c, age, alpha * v10);
    if (y1_ok && v01 > 1e-8f)
      pass(static_cast<float>(x0), static_cast<float>(y1), c, age, alpha * v01);
    if (y1_ok && v11 > 1e-8f)
      pass(static_cast<float>(x1), static_cast<float>(y1), c, age, alpha * v11);
  }
};

/**
 * @brief Manages 2D screen-space trails.
 */
template <int W, int MAX_PIXELS = 1024> class Trails : public Is2DWithHistory {
public:
  /**
   * @brief Constructs a screen trail buffer with the given fade lifetime.
   * @param lifetime Per-frame fade divisor in frames; must be positive.
   * @details lifetime divides per-frame; trap a zero/negative trail length at
   * construction (cold) rather than feeding inf/NaN into the fade.
   */
  Trails(int lifetime) : lifetime(lifetime) { HS_CHECK(lifetime > 0); }

  /**
   * @brief Allocates the decay-pixel storage from the persistent arena.
   * @param arena Persistent arena supplying MAX_PIXELS DecayPixel slots.
   */
  void init_storage(Arena &arena) {
    ttls_ = static_cast<DecayPixel *>(
        arena.allocate(MAX_PIXELS * sizeof(DecayPixel), alignof(DecayPixel)));
    num_pixels = 0;
  }

  /**
   * @brief Forwards the current sample and seeds a fading screen trail point.
   * @param x Column coordinate in pixels.
   * @param y Row coordinate in pixels.
   * @param color Source color, forwarded unchanged this frame.
   * @param age Incoming age (frames); ttl = lifetime - age, seeded only if positive.
   * @param alpha Blend alpha in [0, 1]; samples with alpha <= 0.001 are dropped.
   * @param pass Downstream 2D callback.
   */
  void plot(float x, float y, const Pixel &color, float age, float alpha,
            PassFn2D pass) {
    if (alpha <= 0.001f)
      return;

    pass(x, y, color, age, alpha); // Pass through current frame (mirror World::Trails)

    // Only seed live trail points. age >= lifetime gives a non-positive ttl: the
    // point is already dead, and seeding it would push t = 1-(ttl/lifetime) out
    // of trailFn's range and the re-emitted age past lifetime. The ttl>0 gate
    // mirrors World::Trails and gates only the seed, not the forward above, so an
    // already-aged emission still paints this frame.
    float ttl = static_cast<float>(lifetime) - age;
    // At capacity this drops the NEWEST point — the opposite of World::Trails,
    // whose ring buffer evicts the oldest so the freshest motion stays visible.
    // The asymmetry is deliberate: this buffer is an unordered array (decay()
    // swap-removes dead slots), so there is no front-of-queue "oldest" to evict
    // without an O(n) min-ttl scan on this per-point hot path. Saturation of the
    // MAX_PIXELS=1024 budget is a transient edge regime that clears within a few
    // frames as points decay; paying a per-point scan to reshape it is not worth
    // it. A saturated buffer briefly favors its existing trail over new seeds.
    if (ttl > 0.0f && ttls_ && num_pixels < MAX_PIXELS) {
      ttls_[num_pixels++] = {x, y, ttl};
    }
  }

  /**
   * @brief Re-emits each buffered trail point colored by @p trailFn.
   * @param trailFn Callback producing trail color/alpha from (x, y, t).
   * @param alpha Global blend alpha in [0, 1].
   * @param pass Downstream 2D callback.
   * @details The unused Canvas parameter satisfies the 2D flush signature; ages
   * all points one frame via decay() after emission.
   */
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

  /**
   * @brief Ages every point one frame and swap-removes dead slots.
   * @details Unordered compaction: a dead slot is overwritten by the last live
   * point and the count shrinks.
   */
  void decay() {
    for (int i = 0; i < num_pixels; ++i) {
      if (--ttls_[i].ttl <= 0.0f) {
        ttls_[i] = ttls_[--num_pixels];
        i--;
      }
    }
  }

private:
  /** @brief One screen trail point: position plus remaining lifetime. */
  struct DecayPixel {
    float x, y, ttl; /**< Pixel position and remaining lifetime in frames. */
  };
  int lifetime;               /**< Per-frame fade divisor in frames. */
  DecayPixel *ttls_ = nullptr; /**< Arena-owned array of live trail points. */
  int num_pixels = 0;          /**< Number of live points in ttls_. */
};

/**
 * @brief Applies a variable 3x3 Gaussian blur.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H> class Blur : public Is2D {
public:
  /**
   * @brief Constructs a blur with the given initial strength.
   * @param factor Blur strength in [0, 1] (0 = identity, 1 = full Gaussian).
   */
  Blur(float factor = 1.0f) { update(factor); }

  /**
   * @brief Rebuilds the 3x3 kernel for a new blur strength.
   * @param factor Blur strength; clamped to [0, 1].
   */
  void update(float factor) {
    float f = hs::clamp(factor, 0.0f, 1.0f);
    // Gaussian reference: Corner=1/16, Edge=2/16, Center=4/16
    float c = 1.0f - (0.75f * f);
    float e = 0.125f * f;
    float d = 0.0625f * f;

    kernel = {d, e, d, e, c, e, d, e, d};
  }

  /**
   * @brief Splats the sample across its 3x3 neighborhood weighted by the kernel.
   * @param x Column coordinate in pixels.
   * @param y Row coordinate in pixels.
   * @param color Source color, forwarded to each tap.
   * @param age Temporal age channel (frames), forwarded unchanged.
   * @param alpha Blend alpha in [0, 1]; scaled per tap by its kernel weight.
   * @param pass Downstream 2D callback.
   */
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
  std::array<float, 9> kernel; /**< Row-major 3x3 blur weights, summing to 1. */
};

} // namespace Screen

// ----------------------------------------------------------------------------
// Pixel Space Filters (1:1 with buffer)
// ----------------------------------------------------------------------------
namespace Pixel {

/**
 * @brief Style-aware terminal feedback filter that warps the previous frame.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details Takes a `::Feedback::Style&` directly — drop-in pipeline filter. The
 * Style's spatial transform (noise, melt, etc.) is smooth and expensive (noise
 * + atan2 + acos per call), so the warp field is computed on a coarse
 * W/DS x H/DS grid (allocated from scratch_arena_a per flush) and bilinearly
 * upsampled. DS is read from style.downsample.
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
  /** @brief Marks this as terminal: flush() writes the Canvas and ignores `pass`. */
  static constexpr bool is_terminal = true;

  /**
   * @brief Binds the filter to a live feedback Style.
   * @param style Style supplying the spatial warp and color transforms.
   */
  explicit Feedback(::Feedback::Style &style) : style_(&style) {}

  /**
   * @brief Pass-through: current-frame pixels go straight to the next filter.
   * @param x Column coordinate in pixels.
   * @param y Row coordinate in pixels.
   * @param color Source color, forwarded unchanged.
   * @param age Temporal age channel (frames), forwarded unchanged.
   * @param alpha Blend alpha in [0, 1], forwarded unchanged.
   * @param pass Downstream 2D callback.
   */
  void plot(float x, float y, const ::Pixel &color, float age, float alpha,
            PassFn2D pass) {
    pass(x, y, color, age, alpha);
  }

  /**
   * @brief Blends the distorted previous frame into the current frame.
   * @param cv Target canvas (reads cv.prev, writes the front buffer).
   * @param alpha Global blend alpha in [0, 1].
   * @details Computes a coarse warp field via the Style's space_fn, bilinearly
   * upsamples it, then composites the warped previous frame. Honors the
   * segment's cylindrical clip and ignores the unused ScreenTrailFn / PassFn2D
   * parameters (terminal stage). No-op when disabled.
   */
  void flush(Canvas &cv, const ScreenTrailFn &, float alpha, PassFn2D) {
    if (!enabled_) return;

    const int ds = style_->downsample;
    // Cold authoring/config guard: a downsample that doesn't divide the
    // resolution would silently no-op the whole effect. Trap so the typo
    // surfaces instead of feedback simply vanishing (use enabled_ to turn
    // feedback off). Runs once per flush(), never in the per-pixel loop.
    HS_CHECK(ds > 0 && W % ds == 0 && H % ds == 0,
             "feedback downsample %d must be > 0 and divide %dx%d", ds, W, H);
    const int hw = W / ds;
    const int hh = H / ds;

    if (!any_pixel_lit(cv)) return;

    // Allocate coarse warp deltas (signed 8.8 fixed-point) from scratch.
    //
    // SCRATCH ARENA CONTRACT (load-bearing): this flush and Plot::rasterize
    // both checkpoint scratch_arena_a TODAY — Plot::rasterize opens its own
    // reciprocal ScratchScope on it (see the matching contract note at
    // plot.h rasterize). They share it safely ONLY because plot() and flush()
    // never overlap in time — every effect runs them as separate phases (see
    // MeshFeedback::draw_frame: flush() fully returns, closing this
    // ScratchScope, before any Plot::Mesh::draw opens its own).
    // The dx/dy below stay live for the whole flush, so nothing that also draws
    // from scratch_arena_a may run between here and the end of this function.
    // If a future filter ever allocates from scratch_arena_a inside its own
    // plot() path AND interleaves with a feedback flush, move one of the two to
    // scratch_arena_b. (scratch_arena_b is not used here because it would then
    // share with the mesh-generation path, which is no safer.)
    ScratchScope scope(scratch_arena_a);
    auto *dx = static_cast<int16_t *>(scope.get_arena().allocate(
        hh * hw * sizeof(int16_t), alignof(int16_t)));
    auto *dy = static_cast<int16_t *>(scope.get_arena().allocate(
        hh * hw * sizeof(int16_t), alignof(int16_t)));
    // No OOM guard here: Arena::allocate traps on over-allocation and never
    // returns null, so a "skip this frame" fallback would be dead code that
    // also implies a recovery path the fail-fast model deliberately rejects.

    // Step 2 (sampling) honors the segment's cylindrical x-clip; precompute that
    // clip here so step 1 populates ONLY the coarse columns the sampling pass
    // will read. Without this, a feedback effect on segmented Phantasm recomputes
    // the whole sphere's warp field on every board (4x wasted space_fn) and then
    // uses only its own band.
    const auto &cr = cv.clip();
    const int y_lo = cr.render_y_start();
    const int y_hi = cr.render_y_end();
    // A full-width x range — or one whose margin expansion wraps to cover the
    // whole circumference (rs == re) — needs no x clipping; XClip folds both of
    // those cases into `active`, so this stays in lockstep with scan.h.
    const auto xc = cr.x_clip();

    // Mark the coarse columns the clipped sampling pass actually touches. Step 2
    // reads, for each in-band full-res x, coarse column cx0 = x/ds and its
    // seam-wrapping right neighbour cx1; mark exactly those (the clip test
    // mirrors step 2's below). Unmarked columns are never sampled, so leaving
    // their dx/dy uninitialized is safe. Sized to W (>= hw, since ds >= 1).
    bool col_used[W];
    if (xc.active) {
      for (int cx = 0; cx < hw; ++cx) col_used[cx] = false;
      for (int x = 0; x < W; ++x) {
        if (xc.clipped(x)) continue;
        int cx0 = x / ds;
        int cx1 = (cx0 + 1 < hw) ? cx0 + 1 : 0;
        col_used[cx0] = true;
        col_used[cx1] = true;
      }
    }

    // 1) Populate coarse warp field. Delta encoding keeps the bilerp safe
    //    across the longitude seam — neighbouring absolute bx values can
    //    straddle x=W ↔ x=0, but their deltas are continuous.
    constexpr float Q = 256.0f;
    // Step 2 samples only y in [y_lo, y_hi), reading coarse rows cy0 = y/ds and
    // cy1 = cy0+1 (clamped to hh-1). Populate exactly that coarse-row band so a
    // y-band-clipped segment (segmented Phantasm) skips space_fn for rows it
    // never composites — the y-axis counterpart to the col_used column pruning
    // above. Rows outside the band stay uninitialized but are never read,
    // mirroring that pruning. A full canvas spans all hh rows, so it does the
    // same work either way.
    const int cy_lo = y_lo / ds;
    int cy_hi = ((y_hi - 1) / ds) + 1;
    if (cy_hi > hh - 1) cy_hi = hh - 1;
    for (int cy = cy_lo; cy <= cy_hi; ++cy) {
      int y = cy * ds;
      for (int cx = 0; cx < hw; ++cx) {
        if (xc.active && !col_used[cx]) continue;
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
    //    Honor the segment clip exactly like every other rasterizer (scan.h):
    //    iterate y over the margin-expanded render band and skip x columns
    //    outside the cylindrical x clip (clip computed above for step 1).
    //    Without this, a feedback effect on segmented Phantasm would have every
    //    board composite the whole sphere instead of just its own band — wrong
    //    output and 4x wasted work.
    constexpr float INV_Q = 1.0f / Q;
    const float inv_ds = 1.0f / ds;
    const float fade = style_->fade;
    // Precompute hue_shift's rotation once per frame; the default hue_fade
    // color_fn reads style_->hue_ca/hue_sa instead of recomputing sin/cos of
    // this frame-constant angle for every one of the W*H pixels below.
    style_->sync_hue();
    for (int y = y_lo; y < y_hi; ++y) {
      int cy0 = y / ds;
      int cy1 = (cy0 + 1 < hh) ? cy0 + 1 : hh - 1;
      float fy = (y - cy0 * ds) * inv_ds;
      float wy0 = 1.0f - fy, wy1 = fy;
      const int row0 = cy0 * hw, row1 = cy1 * hw;

      // cx0 / cx1 / fx advance with x rather than a per-pixel divide+modulo
      // (ds divides W, guaranteed above): cx0 ticks every ds columns, fx ramps
      // 0..(ds-1)/ds and resets, cx1 wraps to 0 at the longitude seam. `sub`
      // tracks x - cx0*ds, so fx here is bit-identical to (x - cx0*ds)*inv_ds.
      // The bookkeeping advances every column even for clipped x so the ramp
      // stays in lockstep with the absolute coordinate.
      int cx0 = 0, sub = 0;
      for (int x = 0; x < W; ++x) {
        int cx1 = (cx0 + 1 < hw) ? cx0 + 1 : 0;
        float fx = sub * inv_ds;

        if (!xc.clipped(x)) {
          float wx0 = 1.0f - fx, wx1 = fx;

          int i00 = row0 + cx0;
          int i10 = row0 + cx1;
          int i01 = row1 + cx0;
          int i11 = row1 + cx1;

          float w00 = wx0 * wy0, w10 = wx1 * wy0;
          float w01 = wx0 * wy1, w11 = wx1 * wy1;

          float ddx = (dx[i00] * w00 + dx[i10] * w10
                     + dx[i01] * w01 + dx[i11] * w11) * INV_Q;
          float ddy = (dy[i00] * w00 + dy[i10] * w10
                     + dy[i01] * w01 + dy[i11] * w11) * INV_Q;

          ::Pixel sample = sample_bilinear_prev(cv, x + ddx, y + ddy);
          ::Pixel p = style_->color_fn(sample, fade, *style_);

          // x is already in [0,W) and y is within the render band, so no
          // x-wrap or y-bounds guard is needed before the write.
          if (p.r | p.g | p.b)
            cv(x, y) = blend_alpha(alpha)(cv(x, y), p);
        }

        if (++sub == ds) { sub = 0; ++cx0; }
      }
    }
  }

  /**
   * @brief Enables or disables feedback.
   * @param e When false, flush() is skipped entirely.
   */
  void set_enabled(bool e) { enabled_ = e; }

  /**
   * @brief Accesses the bound Style.
   * @return Mutable reference to the bound feedback Style.
   */
  ::Feedback::Style &style() { return *style_; }
  /**
   * @brief Accesses the bound Style.
   * @return Const reference to the bound feedback Style.
   */
  const ::Feedback::Style &style() const { return *style_; }

private:
  /**
   * @brief Tests whether the previous frame has any non-black pixel.
   * @param cv Canvas whose previous-frame buffer is scanned.
   * @return True on the first lit pixel found, false if the frame is all black.
   * @details Returns on the first lit pixel, so a non-empty frame (the common
   * warm case) costs ~one read; only a fully black frame — exactly when the
   * expensive flush is skipped anyway — pays a full scan. A whole-frame test
   * rather than per-row culling: a warped trail lights at least one pixel in
   * nearly every row, so a row mask would saturate without ever saving a sample.
   */
  static bool any_pixel_lit(const Canvas &cv) {
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x) {
        ::Pixel p = cv.prev(x, y);
        if (p.r | p.g | p.b) return true;
      }
    return false;
  }

  /**
   * @brief Bilinearly samples the Canvas front buffer (previous frame).
   * @param cv Canvas whose previous-frame buffer is sampled.
   * @param bx Fractional column; wrapped across the longitude seam.
   * @param by Fractional row; out-of-range rows contribute black.
   * @return Bilinearly interpolated pixel color.
   */
  ::Pixel sample_bilinear_prev(const Canvas &cv, float bx, float by) const {
    float fy0 = std::floor(by);
    int y0 = static_cast<int>(fy0);
    int y1 = y0 + 1;

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

  ::Feedback::Style *style_;  /**< Bound feedback Style (non-owning). */
  bool enabled_ = true;       /**< When false, flush() is skipped entirely. */
};

/**
 * @brief Splits RGB into per-channel copies offset by 1/2/3 columns,
 * producing a chromatic-aberration fringe.
 */
template <int W> class ChromaticShift : public Is2D {
public:
  /** @brief Constructs the chromatic-shift filter (stateless). */
  ChromaticShift() {}

  /**
   * @brief Emits the source pixel plus R/G/B copies offset by 1/2/3 columns.
   * @param x Column coordinate in pixels.
   * @param y Row coordinate in pixels.
   * @param c Source color; split into single-channel copies.
   * @param age Temporal age channel (frames), forwarded unchanged.
   * @param alpha Blend alpha in [0, 1], forwarded unchanged.
   * @param pass Downstream 2D callback.
   */
  void plot(float x, float y, const ::Pixel &c, float age, float alpha,
            PassFn2D pass) {
    // age/alpha are non-negative by contract. Debug-only trap — stripped on
    // device, fires in the native tests / WASM-debug.
    assert(age >= 0.0f && alpha >= 0.0f);
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

