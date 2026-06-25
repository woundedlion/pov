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
#include <bitset>
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
 * @details `is_terminal`: writes the Canvas directly in flush() and ignores its
 * `pass` callback (must be the last stage). `emits_nonunit_world` /
 * `requires_unit_world_input`: a non-unit-emitting world stage must not precede
 * a unit-assuming one. `crosses_segments`: per-frame state reads pixels outside
 * the worker's segment band, so the effect must render the full canvas; defaults
 * to `has_history` (fail-safe).
 */
struct Is2D {
  static constexpr bool is_2d = true;
  static constexpr bool has_history = false;
  static constexpr bool is_terminal = false;
  static constexpr bool emits_nonunit_world = false;
  static constexpr bool requires_unit_world_input = false;
  static constexpr bool crosses_segments = has_history;
};
/** @brief Trait indicating a filter operates in 3D world space. */
struct Is3D {
  static constexpr bool is_2d = false;
  static constexpr bool has_history = false;
  static constexpr bool is_terminal = false;
  static constexpr bool emits_nonunit_world = false;
  static constexpr bool requires_unit_world_input = false;
  static constexpr bool crosses_segments = has_history;
};

/** @brief Trait indicating a 2D filter that maintains state/history. */
struct Is2DWithHistory {
  static constexpr bool is_2d = true;
  static constexpr bool has_history = true;
  static constexpr bool is_terminal = false;
  static constexpr bool emits_nonunit_world = false;
  static constexpr bool requires_unit_world_input = false;
  static constexpr bool crosses_segments = has_history;
};

/** @brief Trait indicating a 3D filter that maintains state/history. */
struct Is3DWithHistory {
  static constexpr bool is_2d = false;
  static constexpr bool has_history = true;
  static constexpr bool is_terminal = false;
  static constexpr bool emits_nonunit_world = false;
  static constexpr bool requires_unit_world_input = false;
  static constexpr bool crosses_segments = has_history;
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
  static constexpr bool any_crosses_segments = false;

  /**
   * @brief Type-safe filter accessor (base case: T not found).
   * @tparam T Filter type being looked up.
   * @return Never returns; instantiation is a hard error.
   * @details Dependent-false guard: fires only when get<T>() is named on a
   * pipeline lacking T.
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
    // Producer must keep x in [-W, 2W); fast_wrap corrects only a single ±W offset.
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
    // Non-finite coords make the int casts below UB and bypass the wrap.
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

  static constexpr bool any_crosses_segments =
      Head::crosses_segments || Next::any_crosses_segments;

  /**
   * @brief Forwarding-reference constructor: builds Head and the Tail in place.
   * @tparam HArg Argument type forwarded to Head's constructor.
   * @tparam TArgs Argument types forwarded to the remaining stages.
   * @param h Argument forwarded to Head's constructor.
   * @param t Arguments forwarded to the Tail pipeline's constructors.
   * @details The requires guard excludes Pipeline so this template does not
   * hijack copy construction.
   */
  template <typename HArg, typename... TArgs>
    requires(!std::is_same_v<std::remove_cvref_t<HArg>, Pipeline>)
  Pipeline(HArg &&h, TArgs &&...t)
      : Head(std::forward<HArg>(h)), next(std::forward<TArgs>(t)...) {}

  /**
   * @brief Partial constructor: builds Head only, default-constructing the Tail.
   * @tparam HArg Argument type forwarded to Head's constructor.
   * @param h Argument forwarded to Head's constructor.
   * @details Same Pipeline-excluding guard as the variadic ctor.
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
                 [&](float nx, float ny, const Pixel &nc, float nage,
                     float nalpha) {
                   next.plot(cv, nx, ny, nc, nage, nalpha);
                 });
    } else {
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
                 [&](const Vector &nv, const Pixel &nc, float nage,
                     float nalpha) {
                   next.plot(cv, nv, nc, nage, nalpha);
                 });
    } else {
      auto p = vector_to_pixel<W, H>(v);
      plot(cv, p.x, p.y, c, age, alpha);
    }
  }

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

  static_assert(
      !Head::is_terminal || sizeof...(Tail) == 0,
      "A terminal filter (e.g. Pixel::Feedback) writes the Canvas directly and "
      "ignores downstream filters — it must be the last stage in the Pipeline.");

  static_assert(
      !Head::is_2d || (... && Tail::is_2d),
      "Filter ordering: a screen-space (2D) filter (Screen::* / Pixel::*) must "
      "not precede a world-space (3D) filter (World::*) — World filters operate "
      "before screen projection. Reorder so every World::* stage comes first.");

  static_assert(
      !Head::emits_nonunit_world || (... && !Tail::requires_unit_world_input),
      "Filter ordering: a World stage that emits non-unit-length points "
      "(World::Trails) must not precede a World stage that requires unit-length "
      "input (World::Mobius / World::Hole). Reorder so the unit-assuming stage "
      "runs first, or renormalize the trail re-emission.");

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

namespace Filter {

namespace World {

/**
 * @brief Rotates 3D points based on a dynamic Orientation.
 * @details Sweeps the orientation's intra-frame SLERP history and offsets `age`
 * by the fractional `(1 - t)`, producing temporal motion blur. The only filter
 * that adjusts age.
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
      : enabled(true), axis(axis.normalized()), orientations(orientations) {}

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

    float projection = v.x * axis.x + v.y * axis.y + v.z * axis.z;
    float dot_val = std::max(-1.0f, std::min(1.0f, projection));
    float t = hs::clamp(1.0f - fast_acos(dot_val) / PI_F, 0.0f, 1.0f);

    size_t count = orientations.size();
    if (count == 0) {
      pass(v, color, age, alpha);
      return;
    }

    size_t idx = static_cast<size_t>(floorf(t * count));
    if (idx >= count)
      idx = count - 1;

    const Orientation<> &q = orientations[idx];
    tween(q, [&](const Quaternion &rot, float tween_t) {
      pass(rotate(v, rot), color, age + (1.0f - tween_t), alpha);
    });
  }

  /**
   * @brief Sets the slicing axis, renormalizing to enforce the unit-length
   * contract that the projection bucket math assumes.
   * @param a New slicing axis (any non-zero length; renormalized internally).
   */
  void set_axis(const Vector &a) { axis = a.normalized(); }

  /**
   * @brief Accesses the current (unit-length) slicing axis.
   * @return The unit axis points are projected onto to select a slice.
   */
  const Vector &get_axis() const { return axis; }

  bool enabled;  /**< When false, the filter passes points through unrotated. */

private:
  Vector axis;   /**< Unit axis points are projected onto to select a slice. */
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
  static constexpr bool requires_unit_world_input = true;
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
 * @details Every copy shares the source `age` (replication is spatial, not
 * temporal).
 */
template <int W> class Replicate : public Is3D {
public:
  /**
   * @brief Builds a replicator emitting @p count evenly-spaced Y-axis copies.
   * @param count Desired copy count; clamped to [1, W].
   * @details `this->count` (the clamped member, declared/initialized first)
   * feeds make_rotation, so count == 0 cannot feed inf into it.
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
 * @details Precomputes rotation quaternions from vertex[0] to each other vertex.
 * Every copy carries the source age unchanged (replication is spatial).
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
  static constexpr bool requires_unit_world_input = true;
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
  static constexpr bool emits_nonunit_world = true;

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
   * @details Upper bound is structural: ttl is a uint8_t, so lifetime > 255
   * would wrap the trail length.
   */
  Trails(int lifetime) : lifetime(lifetime) {
    HS_CHECK(lifetime > 0 && lifetime <= 255);
  }

  /**
   * @brief Retunes the trail length at runtime (e.g. from a "Trail Len" slider).
   * @param new_lifetime New fade divisor in frames; must be in [1, 255].
   * @details Same bounds as the constructor; buffered points keep their ttl and
   * age out under the new length within a few frames.
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
    pass(v, color, age, alpha);

    // round, not truncate (ttl is an integer byte)
    int ttl = lifetime - static_cast<int>(age + 0.5f);
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
    for (size_t i = 0; i < count_; ++i) {
      auto &item = at(i);
      if (item.ttl > 0)
        item.ttl--;
    }

    while (count_ > 0 && at(0).ttl == 0) {
      pop_front();
    }

    for (size_t i = 0; i < count_; ++i) {
      const auto &item = at(i);
      if (item.ttl == 0)
        continue;
      Vector v = decode(item);
      float t = hs::clamp(
          1.0f - (static_cast<float>(item.ttl) / static_cast<float>(lifetime)),
          0.0f, 1.0f);
      Color4 c = trailFn(v, t);

      if (c.alpha > 0.001f) {
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
    // clamp before scaling: an unclamped component past 1 overflows int16
    return {static_cast<int16_t>(hs::clamp(v.x, -1.0f, 1.0f) * Q),
            static_cast<int16_t>(hs::clamp(v.y, -1.0f, 1.0f) * Q),
            static_cast<int16_t>(hs::clamp(v.z, -1.0f, 1.0f) * Q), ttl, 0};
  }
  /**
   * @brief Decodes a quantized Item back into a near-unit vector.
   * @param item Packed trail sample.
   * @return Reconstructed world-space point.
   * @note Only *near* unit length (int16 quantization), and not renormalized;
   * a World::Trails must not precede a unit-assuming World filter.
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
      pop_front();
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

namespace Screen {

/**
 * @brief Applies 2D anti-aliasing to sub-pixel coordinates.
 * @details Distributes intensity to the 4 nearest neighbors using a quintic kernel.
 */
template <int W, int H> class AntiAlias : public Is2D {
public:
  /**
   * @brief Splats a sub-pixel sample across its four nearest pixel neighbors.
   * @tparam PassFnT Downstream 2D callback type.
   * @param x Sub-pixel column coordinate.
   * @param y Sub-pixel row coordinate.
   * @param c Source color, forwarded to each tap.
   * @param age Temporal age channel (frames), forwarded unchanged.
   * @param alpha Blend alpha in [0, 1]; scaled per tap by its bilinear weight.
   * @param pass Downstream 2D callback receiving each weighted tap.
   * @details Both axes are eased with a quintic kernel; the splat is uniform in
   * framebuffer space at every latitude (no sin(phi) density compensation).
   * @p pass is a forwarding-reference template so the densest fan-out in the
   * family inlines its taps.
   */
  template <typename PassFnT>
  void plot(float x, float y, const Pixel &c, float age, float alpha,
            PassFnT &&pass) {
    assert(age >= 0.0f && alpha >= 0.0f);
    float y_i = floorf(y);
    float y_m = y - y_i;

    float x_floor = floorf(x);
    float x_m = x - x_floor;

    int yi = static_cast<int>(y_i);

    float xs = quintic_kernel(x_m);
    float ys = quintic_kernel(y_m);

    int y0 = yi;
    int y1 = y0 + 1;
    int x0 = fast_wrap(static_cast<int>(x_floor), W);
    int x1 = fast_wrap(static_cast<int>(x_floor) + 1, W);

    bool y0_ok = y0 >= 0 && y0 < H;
    bool y1_ok = y1 >= 0 && y1 < H;

    float wy0 = 1.0f - ys;
    float wy1 = ys;
    if (y0_ok && !y1_ok) {
      wy0 = 1.0f;
      wy1 = 0.0f;
    } else if (!y0_ok && y1_ok) {
      wy0 = 0.0f;
      wy1 = 1.0f;
    }

    float v00 = (1 - xs) * wy0;
    float v10 = xs * wy0;
    float v01 = (1 - xs) * wy1;
    float v11 = xs * wy1;

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
  static constexpr bool crosses_segments = false;

  /**
   * @brief Constructs a screen trail buffer with the given fade lifetime.
   * @param lifetime Per-frame fade divisor in frames; must be positive.
   */
  Trails(int lifetime) : lifetime(lifetime) { HS_CHECK(lifetime > 0); }

  /**
   * @brief Allocates the decay-pixel storage from the persistent arena.
   * @param arena Persistent arena supplying MAX_PIXELS DecayPixel slots.
   */
  void init_storage(Arena &arena) {
    points_ = static_cast<DecayPixel *>(
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

    pass(x, y, color, age, alpha);

    float ttl = static_cast<float>(lifetime) - age;
    if (ttl > 0.0f && points_ && num_pixels < MAX_PIXELS) {
      points_[num_pixels++] = {x, y, ttl};
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
      float t = hs::clamp(1.0f - (points_[i].ttl / lifetime), 0.0f, 1.0f);
      Color4 color = trailFn(points_[i].x, points_[i].y, t);
      if (color.alpha > 0.001f) {
        float age = lifetime - points_[i].ttl;
        if (age < 0.0f)
          age = 0.0f;
        pass(points_[i].x, points_[i].y, color.color, age, alpha * color.alpha);
      }
    }
    decay();
  }

  /**
   * @brief Ages every point one frame and swap-removes dead slots.
   * @details Unordered compaction: a dead slot is overwritten by the last live
   * point. ttl decrements by 1 per frame (whole-frame model), so a point
   * survives ceil(ttl) frames.
   */
  void decay() {
    for (int i = 0; i < num_pixels; ++i) {
      if (--points_[i].ttl <= 0.0f) {
        points_[i] = points_[--num_pixels];
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
  DecayPixel *points_ = nullptr; /**< Arena-owned array of live trail points. */
  int num_pixels = 0;          /**< Number of live points in points_. */
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

    float inv = 1.0f;
    if (cy - 1 < 0 || cy + 1 >= H) {
      float wsum = 0.0f;
      for (int dy = -1; dy <= 1; dy++) {
        int ny = cy + dy;
        if (ny >= 0 && ny < H) {
          int r = (dy + 1) * 3;
          wsum += kernel[r] + kernel[r + 1] + kernel[r + 2];
        }
      }
      if (wsum > 1e-5f)
        inv = 1.0f / wsum;
    }

    int k = 0;
    for (int dy = -1; dy <= 1; dy++) {
      int ny = cy + dy;

      if (ny >= 0 && ny < H) {
        for (int dx = -1; dx <= 1; dx++) {
          float weight = kernel[k++] * inv;
          if (weight > 1e-5f) {
            pass(static_cast<float>(fast_wrap(cx + dx, W)),
                 static_cast<float>(ny), color, age, alpha * weight);
          }
        }
      } else {
        k += 3;
      }
    }
  }

private:
  std::array<float, 9> kernel; /**< Row-major 3x3 blur weights, summing to 1. */
};

} // namespace Screen

namespace Pixel {

/**
 * @brief Style-aware terminal feedback filter that warps the previous frame.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 * @details The Style's spatial warp is computed on a coarse W/DS x H/DS grid
 * (DS = style.downsample, allocated from scratch_arena_a per flush) and
 * bilinearly upsampled. Operates on the full canvas during flush(). TERMINAL:
 * flush() composites directly into the Canvas and ignores its `pass` callback,
 * so it must be the last Pipeline stage.
 */
template <int W, int H>
class Feedback : public Is2DWithHistory {
public:
  /** @brief Marks this as terminal: flush() writes the Canvas and ignores `pass`. */
  static constexpr bool is_terminal = true;

  static_assert(::Feedback::Style{}.downsample > 0 &&
                    W % ::Feedback::Style{}.downsample == 0 &&
                    H % ::Feedback::Style{}.downsample == 0,
                "Feedback<W,H>: default style downsample must be > 0 and divide "
                "W and H");

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
   * upsamples it, then composites the warped previous frame, honoring the
   * segment's cylindrical clip. No-op when disabled.
   */
  void flush(Canvas &cv, const ScreenTrailFn &, float alpha, PassFn2D) {
    if (!enabled_) return;

    const int ds = style_->downsample;
    HS_CHECK(ds > 0 && W % ds == 0 && H % ds == 0,
             "feedback downsample %d must be > 0 and divide %dx%d", ds, W, H);
    const int hw = W / ds;
    const int hh = H / ds;

    if (!any_pixel_lit(cv)) return;

    // LIFO scope reclaims dx/dy on return; must not be read after that.
    ScratchScope scope(scratch_arena_a);
    auto *dx = static_cast<int16_t *>(scope.get_arena().allocate(
        hh * hw * sizeof(int16_t), alignof(int16_t)));
    auto *dy = static_cast<int16_t *>(scope.get_arena().allocate(
        hh * hw * sizeof(int16_t), alignof(int16_t)));

    const auto &cr = cv.clip();
    const int y_lo = cr.render_y_start();
    const int y_hi = cr.render_y_end();
    const auto xc = cr.x_clip();

    // Mark the coarse columns the clipped sampling pass touches.
    std::bitset<W> col_used;
    if (xc.active) {
      for (int x = 0; x < W; ++x) {
        if (xc.clipped(x)) continue;
        int cx0 = x / ds;
        int cx1 = (cx0 + 1 < hw) ? cx0 + 1 : 0;
        col_used[cx0] = true;
        col_used[cx1] = true;
      }
    }

    // 1) Populate coarse warp field as seam-continuous deltas (int16 1/128 px).
    constexpr float Q = 128.0f;
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
    constexpr float INV_Q = 1.0f / Q;
    const float inv_ds = 1.0f / ds;
    const float fade = style_->fade;
    // A non-finite fade would saturate every channel to white (NaN clamps to hi)
    // and persist; trap it here rather than flash the buffer.
    HS_CHECK(std::isfinite(fade), "feedback fade is non-finite");
    style_->sync_hue();
    for (int y = y_lo; y < y_hi; ++y) {
      int cy0 = y / ds;
      int cy1 = (cy0 + 1 < hh) ? cy0 + 1 : hh - 1;
      // bilerping uninitialized scratch is silent corruption, not a crash
      HS_CHECK(cy0 >= cy_lo && cy1 <= cy_hi,
               "feedback warp row %d outside populated band [%d,%d]", cy1, cy_lo,
               cy_hi);
      float fy = (y - cy0 * ds) * inv_ds;
      float wy0 = 1.0f - fy, wy1 = fy;
      const int row0 = cy0 * hw, row1 = cy1 * hw;

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

          // Re-center the three taps within W/2 of d00 for a seam-safe blend.
          constexpr float WQ = static_cast<float>(W) * Q;
          constexpr float HALF_WQ = WQ * 0.5f;
          float d00 = dx[i00], d10 = dx[i10], d01 = dx[i01], d11 = dx[i11];
          d10 += (d10 - d00 > HALF_WQ) ? -WQ : (d10 - d00 < -HALF_WQ ? WQ : 0.0f);
          d01 += (d01 - d00 > HALF_WQ) ? -WQ : (d01 - d00 < -HALF_WQ ? WQ : 0.0f);
          d11 += (d11 - d00 > HALF_WQ) ? -WQ : (d11 - d00 < -HALF_WQ ? WQ : 0.0f);

          float ddx = (d00 * w00 + d10 * w10
                     + d01 * w01 + d11 * w11) * INV_Q;
          float ddy = (dy[i00] * w00 + dy[i10] * w10
                     + dy[i01] * w01 + dy[i11] * w11) * INV_Q;

          ::Pixel sample = sample_bilinear_prev(cv, x + ddx, y + ddy);
          ::Pixel p = style_->color_fn(sample, fade, *style_);

          // write black too, to overwrite the stale double-buffer frame
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
   * @details Scans only this segment's clip band so another board's lit pixels
   * do not gate this board's flush.
   */
  static bool any_pixel_lit(const Canvas &cv) {
    const auto &cr = cv.clip();
    const auto xc = cr.x_clip();
    for (int y = cr.render_y_start(); y < cr.render_y_end(); ++y)
      for (int x = 0; x < W; ++x) {
        if (xc.active && xc.clipped(x)) continue;
        ::Pixel p = cv.prev(x, y);
        if (p.r | p.g | p.b) return true;
      }
    return false;
  }

  /**
   * @brief Bilinearly samples the Canvas front buffer (previous frame).
   * @param cv Canvas whose previous-frame buffer is sampled.
   * @param bx Fractional column in [-W, 2W) (producer contract); wrapped across
   *   the longitude seam by the family's single-step fast_wrap.
   * @param by Fractional row; out-of-range rows contribute black.
   * @return Bilinearly interpolated pixel color, round-to-nearest per channel.
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

    // Producer must keep bx in [-W, 2W); fast_wrap corrects only a single step.
    x0 = fast_wrap(x0, W);
    x1 = fast_wrap(x1, W);

    ::Pixel p00 = (y0 >= 0 && y0 < H) ? cv.prev(x0, y0) : ::Pixel(0, 0, 0);
    ::Pixel p10 = (y0 >= 0 && y0 < H) ? cv.prev(x1, y0) : ::Pixel(0, 0, 0);
    ::Pixel p01 = (y1 >= 0 && y1 < H) ? cv.prev(x0, y1) : ::Pixel(0, 0, 0);
    ::Pixel p11 = (y1 >= 0 && y1 < H) ? cv.prev(x1, y1) : ::Pixel(0, 0, 0);

    float w00 = (1.0f - fx) * (1.0f - fy);
    float w10 = fx * (1.0f - fy);
    float w01 = (1.0f - fx) * fy;
    float w11 = fx * fy;

    // clamp guards the cast against overshoot and maps NaN to the hi bound.
    float r = p00.r * w00 + p10.r * w10 + p01.r * w01 + p11.r * w11;
    float g = p00.g * w00 + p10.g * w10 + p01.g * w01 + p11.g * w11;
    float b = p00.b * w00 + p10.b * w10 + p01.b * w01 + p11.b * w11;
    return ::Pixel(static_cast<uint16_t>(hs::clamp(r, 0.0f, 65535.0f) + 0.5f),
                   static_cast<uint16_t>(hs::clamp(g, 0.0f, 65535.0f) + 0.5f),
                   static_cast<uint16_t>(hs::clamp(b, 0.0f, 65535.0f) + 0.5f));
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

    int xi = fast_wrap(static_cast<int>(x), W);
    pass(static_cast<float>(fast_wrap(xi + 1, W)), y, r_col, age, alpha);
    pass(static_cast<float>(fast_wrap(xi + 2, W)), y, g_col, age, alpha);
    pass(static_cast<float>(fast_wrap(xi + 3, W)), y, b_col, age, alpha);
  }
};

} // namespace Pixel

} // namespace Filter

