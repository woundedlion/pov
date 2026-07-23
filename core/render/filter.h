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
#include "math/geometry.h"
#include "color/color.h"
#include "engine/static_circular_buffer.h"
#include "render/canvas.h"
#include "engine/concepts.h"
#include "engine/memory.h"
#include "engine/styles.h"

/** @brief Callback that forwards a 2D plot (x, y, pixel, age, alpha) downstream. */
using PassFn2D =
    FunctionRef<void(float, float, const Pixel &, float, float)>;
/** @brief Callback that forwards a 3D plot (vector, pixel, age, alpha) downstream. */
using PassFn3D =
    FunctionRef<void(const Vector &, const Pixel &, float, float)>;

/**
 * @brief Trait indicating a filter operates in 2D screen space.
 * @details `is_terminal`: writes the Canvas directly in flush() and ignores its
 * `pass` callback (must be the last stage). `terminal_replaces`: a terminal that
 * overwrites the whole frame (Feedback's opaque store), so no history-bearing
 * stage may precede it — its flush emissions would be clobbered.
 * `emits_nonunit_world` /
 * `requires_unit_world_input`: a non-unit-emitting world stage must not precede
 * a unit-assuming one. `crosses_segments`: per-frame state reads pixels outside
 * the worker's segment band, so the effect must render the full canvas; defaults
 * to `has_history` (fail-safe).
 */
struct Is2D {
  static constexpr bool is_2d = true;
  static constexpr bool has_history = false;
  static constexpr bool is_terminal = false;
  static constexpr bool terminal_replaces = false;
  static constexpr bool emits_nonunit_world = false;
  static constexpr bool requires_unit_world_input = false;
  static constexpr bool crosses_segments = has_history;
};
/** @brief Trait indicating a filter operates in 3D world space. */
struct Is3D {
  static constexpr bool is_2d = false;
  static constexpr bool has_history = false;
  static constexpr bool is_terminal = false;
  static constexpr bool terminal_replaces = false;
  static constexpr bool emits_nonunit_world = false;
  static constexpr bool requires_unit_world_input = false;
  static constexpr bool crosses_segments = has_history;
};

/** @brief Trait indicating a 2D filter that maintains state/history. */
struct Is2DWithHistory {
  static constexpr bool is_2d = true;
  static constexpr bool has_history = true;
  static constexpr bool is_terminal = false;
  static constexpr bool terminal_replaces = false;
  static constexpr bool emits_nonunit_world = false;
  static constexpr bool requires_unit_world_input = false;
  static constexpr bool crosses_segments = has_history;
};

/** @brief Trait indicating a 3D filter that maintains state/history. */
struct Is3DWithHistory {
  static constexpr bool is_2d = false;
  static constexpr bool has_history = true;
  static constexpr bool is_terminal = false;
  static constexpr bool terminal_replaces = false;
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
 * @brief Probe callable for the has_world_cull detection below.
 */
struct PipelineCullEdgeProbe {
  bool operator()(const Vector &, const Vector &, const Basis *) const {
    return true;
  }
};

/**
 * @brief Terminal node of the pipeline (base case). Writes final pixels to the
 * Canvas.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
HS_O3_BEGIN
template <int W, int H> struct Pipeline<W, H> {
  static constexpr bool is_2d = true;
  static constexpr bool any_crosses_segments = false;
  /** @brief No stage re-emits clip-cull edges (see the recursive case). */
  static constexpr bool has_world_cull = false;

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
    Pixel &dst = cv(xi, y);
    dst = blend_alpha(alpha)(dst, c);
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
    // fast_wrap corrects only a single ±W offset, so xi must land in [-W, 2W).
    assert(xi >= -W && xi < 2 * W);
    if (!cv.clip().contains_y(yi)) return;
    xi = fast_wrap(xi, W);
    if (!cv.clip().contains_x(xi)) return;
    cv(xi, yi) = blend_alpha(alpha)(cv(xi, yi), c);
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

  /**
   * @brief Clip-cull terminal: the edge has cleared every world stage, so run
   *        the rasterizer's row-span vs clip-band test on it.
   * @tparam Pred Predicate `bool(const Vector&, const Vector&, const Basis*)`.
   * @return pred(a, b, planar_basis).
   */
  template <typename Pred>
  bool could_intersect_clip(const Vector &a, const Vector &b,
                            const Basis *planar_basis, Pred &&pred) const {
    return pred(a, b, planar_basis);
  }
};
HS_O3_END

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
   * @brief True when any stage overrides cull_edge (re-emits clip-cull edges
   *        under plot-time rotations), so a caller may not precompute per-point
   *        screen coordinates from the raw geometry.
   */
  static constexpr bool has_world_cull =
      requires(const Head &h, const Vector &v, const Basis *pb) {
        h.cull_edge(v, v, pb, PipelineCullEdgeProbe{});
      } || Next::has_world_cull;

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
   * @details Unlike the filter-less sink's int overload (which wraps directly),
   * a filtered pipeline promotes to float so the int sample takes the same path
   * as every filter stage. Both agree for in-range ints.
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

  /**
   * @brief Clip-cull: routes the edge through Head's world transform, then Tail.
   * @tparam Pred Predicate `bool(const Vector&, const Vector&, const Basis*)`.
   * @details A stage that moves world geometry (World::Orient) overrides
   *          cull_edge to re-emit the edge under each rotation it applies at
   *          plot() time, so the rasterizer culls by the RENDERED latitude, not
   *          the source geometry; every other stage forwards the edge unchanged.
   *          Returns true once any transformed copy could intersect the band.
   */
  template <typename Pred>
  bool could_intersect_clip(const Vector &a, const Vector &b,
                            const Basis *planar_basis, Pred &&pred) const {
    auto forward = [&](const Vector &fa, const Vector &fb, const Basis *fpb) {
      return next.could_intersect_clip(fa, fb, fpb, pred);
    };
    if constexpr (requires {
                    std::declval<const Head &>().cull_edge(a, b, planar_basis,
                                                           forward);
                  })
      return Head::cull_edge(a, b, planar_basis, forward);
    else
      return forward(a, b, planar_basis);
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
      !Head::has_history || !(... || Tail::terminal_replaces),
      "Filter ordering: a history-bearing stage's flush emissions would be "
      "overwritten by a frame-replacing terminal filter (Pixel::Feedback's "
      "opaque store owns the whole frame). Drop the history stage, or run the "
      "terminal in a compositing mode.");

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
   * @tparam PassFnT Downstream callback type; a forwarding reference so the
   * filter chain inlines with no per-point indirect call.
   * @param pass Downstream 3D callback.
   */
  template <typename PassFnT>
  void plot(const Vector &v, const Pixel &color, float age, float alpha,
            PassFnT &&pass) {
    tween(orientation, [&](const Quaternion &q, float t) {
      pass(rotate(v, q), color, age + (1.0f - t), alpha);
    });
  }

  /**
   * @brief Re-emits a clip-cull edge under the rotation(s) applied at plot time.
   * @tparam FwdFn Downstream cull continuation
   *         `bool(const Vector&, const Vector&, const Basis*)`.
   * @param a,b Edge endpoints in world space (pre-rotation).
   * @param pb Optional planar basis, rotated alongside the endpoints.
   * @param forward Tail-of-pipeline cull continuation.
   * @return True if any tweened copy of the edge could intersect the clip band.
   * @details Mirrors plot()'s tween so the cull spans the same motion-blur sweep
   *          the renderer draws. Without it the rasterizer would cull by the
   *          un-rotated latitude and drop geometry an off-axis orientation moves
   *          into a segment band (docs/segmented_stateful_effects_spec.md).
   */
  template <typename FwdFn>
  bool cull_edge(const Vector &a, const Vector &b, const Basis *pb,
                 FwdFn &&forward) const {
    bool hit = false;
    tween(orientation, [&](const Quaternion &q, float) {
      if (hit)
        return;
      if (pb) {
        Basis rb = rotate(*pb, q);
        hit = forward(rotate(a, q), rotate(b, q), &rb);
      } else {
        hit = forward(rotate(a, q), rotate(b, q), nullptr);
      }
    });
    return hit;
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
  static constexpr bool requires_unit_world_input = true;
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
   * @tparam PassFnT Downstream callback type; a forwarding reference so the
   * filter chain inlines with no per-point indirect call.
   * @param pass Downstream 3D callback.
   * @details Passes through untouched when disabled or the orientation list is empty.
   */
  template <typename PassFnT>
  void plot(const Vector &v, const Pixel &color, float age, float alpha,
            PassFnT &&pass) {
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
   * @brief Re-emits a clip-cull edge under every candidate slice's rotation.
   * @tparam FwdFn Downstream cull continuation
   *         `bool(const Vector&, const Vector&, const Basis*)`.
   * @param a,b Edge endpoints in world space (pre-rotation).
   * @param pb Optional planar basis, rotated alongside the endpoints.
   * @param forward Tail-of-pipeline cull continuation.
   * @return True if any candidate slice's tweened copy could intersect the band.
   * @details The endpoints may fall in different slices, so bound conservatively
   *          over all candidates rather than replicating the per-point selector.
   */
  template <typename FwdFn>
  bool cull_edge(const Vector &a, const Vector &b, const Basis *pb,
                 FwdFn &&forward) const {
    if (!enabled || orientations.empty())
      return forward(a, b, pb);
    for (const Orientation<> &o : orientations) {
      bool hit = false;
      tween(o, [&](const Quaternion &q, float) {
        if (hit)
          return;
        if (pb) {
          Basis rb = rotate(*pb, q);
          hit = forward(rotate(a, q), rotate(b, q), &rb);
        } else {
          hit = forward(rotate(a, q), rotate(b, q), nullptr);
        }
      });
      if (hit)
        return true;
    }
    return false;
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
   * @tparam PassFnT Downstream callback type; a forwarding reference so the
   * filter chain inlines with no per-point indirect call.
   * @param pass Downstream 3D callback.
   */
  template <typename PassFnT>
  void plot(const Vector &v, const Pixel &color, float age, float alpha,
            PassFnT &&pass) {
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
   * @tparam PassFnT Downstream callback type; a forwarding reference so the
   * filter chain inlines with no per-point indirect call.
   * @param pass Downstream 3D callback.
   */
  template <typename PassFnT>
  void plot(const Vector &v, const Pixel &color, float age, float alpha,
            PassFnT &&pass) {
    Vector r = v;
    pass(r, color, age, alpha);
    for (int i = 1; i < count; i++) {
      // renormalize so repeated rotation can't drift copies off the unit sphere
      r = rotate(r, step).normalized();
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
   * @tparam PassFnT Downstream callback type; a forwarding reference so the
   * filter chain inlines with no per-point indirect call.
   * @param pass Downstream 3D callback.
   */
  template <typename PassFnT>
  void plot(const Vector &v, const Pixel &color, float age, float alpha,
            PassFnT &&pass) {
    for (int i = 0; i < N; ++i) {
      pass(rotate(v, rotations[i]), color, age, alpha);
    }
  }

private:
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
   * @tparam PassFnT Downstream callback type; a forwarding reference so the
   * filter chain inlines with no per-point indirect call.
   * @param pass Downstream 3D callback.
   */
  template <typename PassFnT>
  void plot(const Vector &v, const Pixel &color, float age, float alpha,
            PassFnT &&pass) {
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
    uint8_t pad;     /**< Padding for 8-byte alignment (1 byte). */
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
    items = arena.allocate_n<Item>(Capacity);
    head = tail = count = 0;
  }

  /**
   * @brief Forwards the current sample and seeds a fading trail point.
   * @param v World-space point on the unit sphere.
   * @param color Source color, forwarded unchanged this frame.
   * @param age Incoming age (frames); ttl = lifetime - age, seeded only if positive.
   * @param alpha Blend alpha in [0, 1], forwarded unchanged.
   * @tparam PassFnT Downstream callback type; a forwarding reference so the
   * filter chain inlines with no per-point indirect call.
   * @param pass Downstream 3D callback.
   */
  template <typename PassFnT>
  void plot(const Vector &v, const Pixel &color, float age, float alpha,
            PassFnT &&pass) {
    pass(v, color, age, alpha);

    // round, not truncate (ttl is an integer byte)
    int ttl = lifetime - static_cast<int>(age + 0.5f);
    if (ttl > 0 && items) {
      push_back(encode(v, static_cast<uint8_t>(ttl)));
    }
  }

  /**
   * @brief Re-emits each buffered point colored by @p trailFn, then ages every
   * point one frame and culls the dead.
   * @param trailFn Callback producing trail color/alpha from (point, t).
   * @param alpha Global blend alpha in [0, 1].
   * @tparam PassFnT Downstream callback type; a forwarding reference so the
   * filter chain inlines with no per-point indirect call.
   * @param pass Downstream 3D callback.
   * @details Emits before aging (matching Screen::Trails::flush), so a point
   * still renders on the frame its ttl reaches 1 rather than being culled unseen.
   */
  template <typename PassFnT>
  void flush(const WorldTrailFn &trailFn, float alpha, PassFnT &&pass) {
    for (size_t i = 0; i < count; ++i) {
      const auto &item = at(i);
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

    for (size_t i = 0; i < count;) {
      Item &item = at(i);
      if (item.ttl > 0)
        item.ttl--;
      if (item.ttl == 0) {
        // swap-remove the logical-last live item (index count-1) into the dead
        // slot; only tail retreats, head stays put.
        item = at(count - 1);
        tail = (tail + Capacity - 1) % Capacity;
        count--;
      } else {
        ++i;
      }
    }
  }

  /**
   * @brief Returns the number of live trail points currently buffered.
   * @return Count of buffered Item entries.
   */
  size_t size() const { return count; }

private:
  Item *items = nullptr;                 /**< Ring-buffer storage (arena-owned). */
  size_t head = 0, tail = 0, count = 0; /**< Ring-buffer head, tail, and live count. */
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
   * @param i Index into the live range [0, count).
   * @return Mutable reference to the buffered Item.
   */
  Item &at(size_t i) { return items[(head + i) % Capacity]; }
  /**
   * @brief Returns the i-th live item by age (0 = oldest).
   * @param i Index into the live range [0, count).
   * @return Const reference to the buffered Item.
   */
  const Item &at(size_t i) const { return items[(head + i) % Capacity]; }

  /**
   * @brief Appends an item, evicting the oldest when at capacity.
   * @param item Encoded trail sample to push.
   */
  void push_back(const Item &item) {
    if (count == Capacity) {
      pop_front();
    }
    items[tail] = item;
    tail = (tail + 1) % Capacity;
    count++;
  }

  /** @brief Drops the oldest buffered item. */
  void pop_front() {
    head = (head + 1) % Capacity;
    count--;
  }
};

} // namespace World

namespace Screen {

/**
 * @brief Applies 2D anti-aliasing to sub-pixel coordinates.
 * @details Distributes intensity to the 4 nearest neighbors using a quintic kernel.
 */
HS_O3_BEGIN
template <int W, int H> class AntiAlias : public Is2D {
public:
  /**
   * @brief Splats a sub-pixel sample across its four nearest pixel neighbors.
   * @tparam PassFnT Downstream 2D callback type.
   * @param x Sub-pixel column coordinate.
   * @param y Sub-pixel row coordinate.
   * @param c Source color, forwarded to each tap.
   * @param age Temporal age channel (frames), forwarded unchanged.
   * @param alpha Blend alpha in [0, 1]; scaled per tap by its quintic-eased splat weight.
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
    int x1 = fast_wrap(x0 + 1, W);

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

    // Skip negligible splats. Cutoff is looser in Blur (1e-5): these are raw
    // bilinear coverage products, Blur's are normalized 3x3 kernel taps.
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
HS_O3_END

/**
 * @brief Terminal four-tap anti-alias sink with direct framebuffer writes.
 * @details Opt-in replacement for `Pipeline<W, H, AntiAlias<W, H>>` when no
 * downstream filter is required. It preserves AntiAlias tap ordering and the
 * base sink's q16 source-over blend while resolving rows, columns and clipping
 * once per sample.
 */
HS_O3_BEGIN
template <int W, int H> class DirectAntiAliasSink : public Is2D {
public:
  static constexpr bool any_crosses_segments = false;
  static constexpr bool has_world_cull = false;
  static constexpr bool direct_raster_path = true;

  /** @brief Splats one screen-space sample directly into the Canvas. */
  void plot(Canvas &cv, float x, float y, const Pixel &c, float age,
            float alpha) {
    assert(age >= 0.0f && alpha >= 0.0f);
    (void)age;

    const float y_floor = floorf(y);
    const float x_floor = floorf(x);
    const float xs = quintic_kernel(x - x_floor);
    const float ys = quintic_kernel(y - y_floor);

    const int y0 = static_cast<int>(y_floor);
    const int y1 = y0 + 1;
    const int x0 = fast_wrap(static_cast<int>(x_floor), W);
    const int x1 = fast_wrap(x0 + 1, W);

    const bool y0_physical = y0 >= 0 && y0 < H;
    const bool y1_physical = y1 >= 0 && y1 < H;
    float wy0 = 1.0f - ys;
    float wy1 = ys;
    if (y0_physical && !y1_physical) {
      wy0 = 1.0f;
      wy1 = 0.0f;
    } else if (!y0_physical && y1_physical) {
      wy0 = 0.0f;
      wy1 = 1.0f;
    }

    const float v00 = (1.0f - xs) * wy0;
    const float v10 = xs * wy0;
    const float v01 = (1.0f - xs) * wy1;
    const float v11 = xs * wy1;
    constexpr float TAP_CUTOFF = 1e-8f;

    const ClipRegion &cr = cv.clip();
    const ClipRegion::XClip xc = cr.x_clip();
    const bool x0_ok = !xc.clipped(x0);
    const bool x1_ok = !xc.clipped(x1);
    const bool y0_ok = y0_physical && cr.contains_y(y0);
    const bool y1_ok = y1_physical && cr.contains_y(y1);
    Pixel *const base = cv.data();
    Pixel *const row0 = y0_ok ? base + y0 * W : nullptr;
    Pixel *const row1 = y1_ok ? base + y1 * W : nullptr;
    Pixel *const dst00 = row0 && x0_ok && v00 > TAP_CUTOFF ? row0 + x0 : nullptr;
    Pixel *const dst10 = row0 && x1_ok && v10 > TAP_CUTOFF ? row0 + x1 : nullptr;
    Pixel *const dst01 = row1 && x0_ok && v01 > TAP_CUTOFF ? row1 + x0 : nullptr;
    Pixel *const dst11 = row1 && x1_ok && v11 > TAP_CUTOFF ? row1 + x1 : nullptr;
    const uint16_t a00 = dst00 ? tap_alpha_q16(alpha, v00) : 0;
    const uint16_t a10 = dst10 ? tap_alpha_q16(alpha, v10) : 0;
    const uint16_t a01 = dst01 ? tap_alpha_q16(alpha, v01) : 0;
    const uint16_t a11 = dst11 ? tap_alpha_q16(alpha, v11) : 0;

    if (dst00)
      blend_tap(dst00, c, a00);
    if (dst10)
      blend_tap(dst10, c, a10);
    if (dst01)
      blend_tap(dst01, c, a01);
    if (dst11)
      blend_tap(dst11, c, a11);
  }

  /** @brief Integer-coordinate overload matching a filtered Pipeline. */
  void plot(Canvas &cv, int x, int y, const Pixel &c, float age, float alpha) {
    plot(cv, static_cast<float>(x), static_cast<float>(y), c, age, alpha);
  }

  /** @brief Projects a world point, then applies the direct screen-space splat. */
  void plot(Canvas &cv, const Vector &v, const Pixel &c, float age,
            float alpha) {
    const PixelCoords p = vector_to_pixel<W, H>(v);
    plot(cv, p.x, p.y, c, age, alpha);
  }

  /** @brief Stateless screen flush no-op. */
  void flush(Canvas &, const ScreenTrailFn &, float) {}
  /** @brief Stateless world flush no-op. */
  void flush(Canvas &, const WorldTrailFn &, float) {}

  /** @brief Terminal clip-cull predicate forwarding. */
  template <typename Pred>
  bool could_intersect_clip(const Vector &a, const Vector &b,
                            const Basis *planar_basis, Pred &&pred) const {
    return pred(a, b, planar_basis);
  }

private:
  static uint16_t tap_alpha_q16(float alpha, float weight) {
    return static_cast<uint16_t>(
        hs::clamp(alpha * weight * 65535.0f + 0.5f, 0.0f, 65535.0f));
  }

  static __attribute__((always_inline)) void blend_tap(Pixel *dst,
                                                        const Pixel &src,
                                                        uint16_t alpha_q16) {
    *dst = dst->lerp16(src, alpha_q16);
  }
};
HS_O3_END

/**
 * @brief Manages 2D screen-space trails.
 */
template <int W, int MAX_PIXELS = 1024> class Trails : public Is2DWithHistory {
public:
  // Trail points are seeded from and re-emitted into the same band, so they
  // never sample a neighbor segment.
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
    points_ = arena.allocate_n<DecayPixel>(MAX_PIXELS);
    num_pixels = 0;
  }

  /**
   * @brief Forwards the current sample and seeds a fading screen trail point.
   * @param x Column coordinate in pixels.
   * @param y Row coordinate in pixels.
   * @param color Source color, forwarded unchanged this frame.
   * @param age Incoming age (frames); ttl = lifetime - age, seeded only if positive.
   * @param alpha Blend alpha in [0, 1]; samples with alpha <= 0.001 seed no
   * trail point but are still forwarded downstream.
   * @tparam PassFnT Downstream callback type; a forwarding reference so the
   * filter chain inlines with no per-point indirect call.
   * @param pass Downstream 2D callback.
   */
  template <typename PassFnT>
  void plot(float x, float y, const Pixel &color, float age, float alpha,
            PassFnT &&pass) {
    pass(x, y, color, age, alpha);

    if (alpha <= 0.001f)
      return;

    float ttl = static_cast<float>(lifetime) - age;
    if (ttl > 0.0f && points_) {
      if (num_pixels == MAX_PIXELS) {
        // At capacity: O(1) drop of slot 0. Saturation eviction order differs
        // by domain (World::Trails drops the FIFO-oldest, Screen drops slot 0);
        // per-point ttl fade absorbs the transient either way.
        num_pixels--;
        points_[0] = points_[num_pixels];
      }
      points_[num_pixels++] = {x, y, ttl};
    }
  }

  /**
   * @brief Re-emits each buffered trail point colored by @p trailFn.
   * @param trailFn Callback producing trail color/alpha from (x, y, t).
   * @param alpha Global blend alpha in [0, 1].
   * @tparam PassFnT Downstream callback type; a forwarding reference so the
   * filter chain inlines with no per-point indirect call.
   * @param pass Downstream 2D callback.
   * @details The unused Canvas parameter satisfies the 2D flush signature; ages
   * all points one frame via decay() after emission.
   */
  template <typename PassFnT>
  void flush(Canvas &, const ScreenTrailFn &trailFn, float alpha,
             PassFnT &&pass) {
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
   * @tparam PassFnT Downstream callback type; a forwarding reference so the
   * filter chain inlines with no per-point indirect call.
   * @param pass Downstream 2D callback.
   */
  template <typename PassFnT>
  void plot(float x, float y, const ::Pixel &color, float age, float alpha,
            PassFnT &&pass) {
    int cx = fast_wrap(static_cast<int>(std::round(x)), W);
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
          // Normalized kernel tap; cutoff is tighter in AntiAlias (1e-8) whose
          // weights are raw bilinear coverage products.
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
 * bilinearly upsampled. flush() iterates the full pixel grid within the active
 * clip band. TERMINAL:
 * flush() composites directly into the Canvas and ignores its `pass` callback,
 * so it must be the last Pipeline stage.
 */
template <int W, int H>
class Feedback : public Is2DWithHistory {
public:
  /** @brief Marks this as terminal: flush() writes the Canvas and ignores `pass`. */
  static constexpr bool is_terminal = true;
  /** @brief Opaque store owns the frame: no history stage may precede it. */
  static constexpr bool terminal_replaces = true;

  // Covers only the default-constructed Style; a runtime-swapped style's
  // downsample is validated in flush() (the HS_CHECK below).
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
   * @tparam PassFnT Downstream callback type; a forwarding reference so the
   * filter chain inlines with no per-point indirect call.
   * @param pass Downstream 2D callback.
   */
  template <typename PassFnT>
  void plot(float x, float y, const ::Pixel &color, float age, float alpha,
            PassFnT &&pass) {
    pass(x, y, color, age, alpha);
  }

  /**
   * @brief Blends the distorted previous frame into the current frame.
   * @param cv Target canvas (reads cv.prev, writes the back (current-draw) buffer).
   * @param alpha Global blend alpha in [0, 1].
   * @details Computes a coarse warp field via the Style's space_fn, bilinearly
   * upsamples it, then composites the warped previous frame, honoring the
   * segment's cylindrical clip. No-op when disabled.
   */
  HS_O3_BEGIN
  void flush(Canvas &cv, const ScreenTrailFn &, float alpha, PassFn2D) {
    if (!enabled_) return;

    const int ds = style_->downsample;
    HS_CHECK(ds > 0 && W % ds == 0 && H % ds == 0,
             "feedback downsample %d must be > 0 and divide %dx%d", ds, W, H);
    // The hoisted base pointers below are indexed with the template-constant W,
    // which only matches the canvas accessors when the strides agree.
    HS_CHECK(cv.width() == W,
             "feedback canvas width %d must equal template W %d", cv.width(), W);
    const int hw = W / ds;
    const int hh = H / ds;

    {
      HS_PROFILE(feedback_litscan);
      if (!any_pixel_lit(cv)) return;
    }

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

    constexpr float Q = 128.0f;
    const int cy_lo = y_lo / ds;
    int cy_hi = ((y_hi - 1) / ds) + 1;
    if (cy_hi > hh - 1) cy_hi = hh - 1;
    HS_CHECK(cy_hi >= cy_lo, "feedback coarse band inverted: [%d,%d]", cy_lo,
             cy_hi);

    // The warp field is a pure function of the key's fields for the stock
    // space transforms (their only time axis is noise time * speed), so equal
    // keys mean the cached field can be reused verbatim. Static presets
    // (speed 0) then skip the whole space_fn pass. An arbitrary user SpaceFn
    // may read state the key can't see; never cache it.
    const NoiseParams *np = style_->noise;
    const bool cacheable =
        warp_dx_ && !xc.active && ds == CACHE_DS &&
        (style_->space_fn == &::Feedback::noise_warp ||
         style_->space_fn == &::Feedback::melt_warp ||
         style_->space_fn == &::Feedback::identity_warp);
    const WarpKey key{style_->space_fn,
                      np,
                      style_->amplitude,
                      style_->frequency,
                      style_->speed,
                      style_->scale,
                      np ? np->time * np->speed : 0.0f,
                      cy_lo,
                      cy_hi};

    // LIFO scope reclaims a scratch dx/dy on return; must not be read after
    // that. Cacheable flushes use the init_storage() arrays instead.
    ScratchScope scope(scratch_arena_a);
    int16_t *dx, *dy;
    bool populate = true;
    if (cacheable) {
      dx = warp_dx_;
      dy = warp_dy_;
      populate = !(warp_cache_valid_ && key == warp_key_);
      warp_key_ = key;
      warp_cache_valid_ = true;
    } else {
      dx = scope.get_arena().allocate_n<int16_t>(hh * hw);
      dy = scope.get_arena().allocate_n<int16_t>(hh * hw);
    }

    // Populate coarse warp field as seam-continuous deltas (int16 1/128 px).
    {
      HS_PROFILE(feedback_populate);
      if (populate)
      for (int cy = cy_lo; cy <= cy_hi; ++cy) {
        int y = cy * ds;
        for (int cx = 0; cx < hw; ++cx) {
          if (xc.active && !col_used[cx]) continue;
          int x = cx * ds;
          Vector v_dist;
          {
            HS_PROFILE_DEEP(fb_pop_warp);
            v_dist = style_->space_fn(pixel_to_vector<W, H>(x, y), *style_);
          }
          HS_PROFILE_DEEP(fb_pop_project);
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
    }

    // Sample at full res, bilerping the coarse warp field per pixel.
    constexpr float INV_Q = 1.0f / Q;
    const float inv_ds = 1.0f / ds;
    const float fade = style_->fade;
    // A non-finite fade would saturate every channel to white (NaN clamps to hi)
    // and persist; trap it here rather than flash the buffer.
    HS_CHECK(std::isfinite(fade), "feedback fade is non-finite");
    style_->sync_hue();
    // The stock color transforms map black to black exactly, so near-black
    // samples can skip them and write black; an arbitrary user ColorFn may
    // not (e.g. a glow floor), so it sees every sample.
    const bool black_skips_color =
        style_->color_fn == &::Feedback::hue_fade ||
        style_->color_fn == &::Feedback::plain_fade;
    // Samples with every channel below this are invisible (~0.1% linear,
    // ~3/255 sRGB), yet round-to-nearest fade keeps them alive forever (at
    // FADE_MAX 0.99, v*fade re-rounds to v for v < 50): without the cut the
    // buffer saturates with immortal dim trails that defeat the skip and
    // never let a region go dark.
    constexpr float NEAR_BLACK = 64.0f;
    // lerp16 at frac 65535 returns the source exactly, so full alpha is a
    // plain store.
    const auto blend = blend_alpha(alpha);
    const bool opaque = alpha >= 1.0f;
    // Dispatch the color transform once per flush so the stock transforms
    // inline into the pixel loop instead of an indirect call per pixel.
    // The four warp taps and their seam-unify are constant across the ds pixels
    // of a coarse cell; only the horizontal fraction fx moves. Precompute per
    // cell the row-weighted, INV_Q-scaled left corner and left->right slope so
    // the per-pixel warp collapses to one fma in fx.
    constexpr float WQ = static_cast<float>(W) * Q;
    constexpr float HALF_WQ = WQ * 0.5f;
    // Hoisted once: the canvas accessors each do a relaxed atomic load plus a
    // runtime-width multiply, which no compiler can lift out of the pixel loop.
    const ::Pixel *prev = cv.prev_data();
    ::Pixel *cur = cv.data();
    // The clip band is row-independent, so the unclipped column runs resolve
    // once here and the pixel loop carries no per-pixel clip test.
    int runs[2][2];
    int nruns = 0;
    if (!xc.active) {
      runs[nruns][0] = 0;
      runs[nruns][1] = W;
      ++nruns;
    } else if (xc.wrap) {
      if (xc.re > 0) {
        runs[nruns][0] = 0;
        runs[nruns][1] = xc.re;
        ++nruns;
      }
      if (xc.rs < W) {
        runs[nruns][0] = xc.rs;
        runs[nruns][1] = W;
        ++nruns;
      }
    } else {
      runs[nruns][0] = xc.rs;
      runs[nruns][1] = xc.re;
      ++nruns;
    }
    auto composite = [&](auto &&color_px, auto &&color_px2, auto pair_on) {
      constexpr bool PAIR = decltype(pair_on)::value;
      for (int y = y_lo; y < y_hi; ++y) {
        const int row = y * W;
        int cy0 = y / ds;
        int cy1 = (cy0 + 1 < hh) ? cy0 + 1 : hh - 1;
        // bilerping uninitialized scratch is silent corruption, not a crash
        HS_CHECK(cy0 >= cy_lo && cy1 <= cy_hi,
                 "feedback warp row %d outside populated band [%d,%d]", cy1,
                 cy_lo, cy_hi);
        float fy = (y - cy0 * ds) * inv_ds;
        float wy0 = 1.0f - fy, wy1 = fy;
        const int row0 = cy0 * hw, row1 = cy1 * hw;

        for (int r = 0; r < nruns; ++r) {
          const int xs = runs[r][0], xe = runs[r][1];
          int cx0 = xs / ds, sub = xs - cx0 * ds;
          float leftx = 0.0f, slopex = 0.0f, lefty = 0.0f, slopey = 0.0f;
          auto cell = [&]() {
            HS_PROFILE_DEEP(fb_comp_cell);
            int cx1 = (cx0 + 1 < hw) ? cx0 + 1 : 0;
            int i00 = row0 + cx0, i10 = row0 + cx1;
            int i01 = row1 + cx0, i11 = row1 + cx1;
            // Unify the far taps onto d00's wrap branch (each within W/2 of the
            // anchor); the anchor must be a tap — a computed mid-value can unify
            // nothing and sweep the blend across the seam to the far hemisphere.
            float d00 = dx[i00], d10 = dx[i10], d01 = dx[i01], d11 = dx[i11];
            d10 += (d10 - d00 > HALF_WQ) ? -WQ : (d10 - d00 < -HALF_WQ ? WQ : 0.0f);
            d01 += (d01 - d00 > HALF_WQ) ? -WQ : (d01 - d00 < -HALF_WQ ? WQ : 0.0f);
            d11 += (d11 - d00 > HALF_WQ) ? -WQ : (d11 - d00 < -HALF_WQ ? WQ : 0.0f);
            leftx = (d00 * wy0 + d01 * wy1) * INV_Q;
            slopex = (d10 * wy0 + d11 * wy1) * INV_Q - leftx;
            lefty = (dy[i00] * wy0 + dy[i01] * wy1) * INV_Q;
            slopey = (dy[i10] * wy0 + dy[i11] * wy1) * INV_Q - lefty;
          };
          // A clipped run can open mid-cell; the loop only refreshes at sub 0.
          if (sub != 0) cell();

          for (int x = xs; x < xe;) {
            if (sub == 0) cell();

            // Two pixels in flight give the in-order FPU two independent
            // dependency chains to overlap. A pair stays inside one cell so the
            // coefficients are shared; a ragged tail (or ds 1) drops to scalar.
            if constexpr (PAIR) {
              if (ds - sub >= 2 && xe - x >= 2) {
                float fx0 = sub * inv_ds;
                float fx1 = (sub + 1) * inv_ds;
                float ddx0 = leftx + slopex * fx0;
                float ddy0 = lefty + slopey * fx0;
                float ddx1 = leftx + slopex * fx1;
                float ddy1 = lefty + slopey * fx1;

                float sr0, sg0, sb0, sr1, sg1, sb1;
                {
                  HS_PROFILE_DEEP(fb_comp_sample);
                  sample_bilinear_prev(prev, x + ddx0, y + ddy0, sr0, sg0, sb0);
                  sample_bilinear_prev(prev, x + 1 + ddx1, y + ddy1, sr1, sg1,
                                       sb1);
                }
                ::Pixel p0(0, 0, 0), p1(0, 0, 0);
                {
                  HS_PROFILE_DEEP(fb_comp_color);
                  color_px2(sr0, sg0, sb0, sr1, sg1, sb1, p0, p1);
                }
                // Both lanes transform unconditionally and the sub-threshold
                // ones are selected to black afterwards; branching on the skip
                // would split the lanes back onto separate code paths.
                const bool black0 = black_skips_color && sr0 < NEAR_BLACK &&
                                    sg0 < NEAR_BLACK && sb0 < NEAR_BLACK;
                const bool black1 = black_skips_color && sr1 < NEAR_BLACK &&
                                    sg1 < NEAR_BLACK && sb1 < NEAR_BLACK;
                p0 = black0 ? ::Pixel(0, 0, 0) : p0;
                p1 = black1 ? ::Pixel(0, 0, 0) : p1;

                HS_PROFILE_DEEP(fb_comp_write);
                ::Pixel &dst0 = cur[row + x];
                dst0 = opaque ? p0 : blend(dst0, p0);
                ::Pixel &dst1 = cur[row + x + 1];
                dst1 = opaque ? p1 : blend(dst1, p1);

                x += 2;
                sub += 2;
                if (sub == ds) { sub = 0; ++cx0; }
                continue;
              }
            }

            float fx = sub * inv_ds;
            float ddx = leftx + slopex * fx;
            float ddy = lefty + slopey * fx;

            float sr, sg, sb;
            {
              HS_PROFILE_DEEP(fb_comp_sample);
              sample_bilinear_prev(prev, x + ddx, y + ddy, sr, sg, sb);
            }
            ::Pixel p(0, 0, 0);
            if (!(black_skips_color && sr < NEAR_BLACK && sg < NEAR_BLACK &&
                  sb < NEAR_BLACK)) {
              HS_PROFILE_DEEP(fb_comp_color);
              p = color_px(sr, sg, sb);
            }

            // write black too, to overwrite the stale double-buffer frame
            HS_PROFILE_DEEP(fb_comp_write);
            ::Pixel &dst = cur[row + x];
            dst = opaque ? p : blend(dst, p);

            ++x;
            if (++sub == ds) { sub = 0; ++cx0; }
          }
        }
      }
    };
    // The second callback is unused when pairing is off, so the scalar entry
    // just repeats the first.
    auto composite_scalar = [&](auto &&color_px) {
      composite(color_px, color_px, std::false_type{});
    };
    // The stock transforms run on the sampler's float channels directly (one
    // quantization at the write); an arbitrary ColorFn keeps the Pixel
    // interface, so the sample is quantized first. An identity hue rotation
    // makes hue_fade a pure fade, so it takes the plain path.
    HS_PROFILE(feedback_composite);
    const bool hue_identity =
        style_->hue_ca == 1.0f && style_->hue_sa == 0.0f;
    if (style_->color_fn == &::Feedback::hue_fade && !hue_identity) {
      // Pre-scale the frame-constant rotation by cbrt(fade/65535), folding the
      // fade and the u16 normalization into the one matrix (see hue_fade).
      float k[9];
      const float sc = fast_cbrt(fade * (1.0f / 65535.0f));
      for (int i = 0; i < 9; ++i)
        k[i] = style_->hue_k[i] * sc;
      composite(
          [&](float r, float g, float b) {
            return ::Feedback::hue_fade_apply(k, r, g, b);
          },
          [&](float r0, float g0, float b0, float r1, float g1, float b1,
              ::Pixel &p0, ::Pixel &p1) {
            ::Feedback::hue_fade_apply2(k, r0, g0, b0, r1, g1, b1, p0, p1);
          },
          std::true_type{});
    } else if (style_->color_fn == &::Feedback::plain_fade ||
               style_->color_fn == &::Feedback::hue_fade) {
      auto plain = [&](float r, float g, float b) {
        return ::Pixel(quantize16(r * fade), quantize16(g * fade),
                       quantize16(b * fade));
      };
      // A plain fade is three multiplies, so pairing it buys only the sampler
      // overlap and does not earn its ITCM.
      composite_scalar(plain);
    } else {
      auto general = [&](float r, float g, float b) {
        return style_->color_fn(
            ::Pixel(quantize16(r), quantize16(g), quantize16(b)), fade,
            *style_);
      };
      // A general ColorFn is an indirect call per pixel, so it stays scalar.
      composite_scalar(general);
    }
  }
  HS_O3_END

  /**
   * @brief Enables or disables feedback.
   * @param e When false, flush() is skipped entirely.
   */
  void set_enabled(bool e) { enabled_ = e; }

  /**
   * @brief Allocates the warp-field cache from the persistent arena.
   * @param arena Persistent arena supplying 2 * CACHE_CELLS int16 slots.
   * @details Must be called from effect init(), not the constructor (arenas
   * aren't ready yet), and again after any compaction that resets the arena —
   * the cache is derived data, so it just re-populates on the next flush.
   * Without storage every flush renders uncached.
   */
  HS_COLD_MEMBER void init_storage(Arena &arena) {
    warp_dx_ = arena.allocate_n<int16_t>(CACHE_CELLS);
    warp_dy_ = arena.allocate_n<int16_t>(CACHE_CELLS);
    warp_cache_valid_ = false;
  }

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
   * @param prev Base of the previous-frame buffer, row-major with stride W.
   * @param bx Fractional column in [-W, 2W) (producer contract); wrapped across
   *   the longitude seam by the family's single-step fast_wrap.
   * @param by Fractional row; out-of-range rows contribute black.
   * @param r Out: interpolated red on the [0, 65535] scale, unquantized.
   * @param g Out: interpolated green.
   * @param b Out: interpolated blue.
   */
  HS_O3_FN
  void sample_bilinear_prev(const ::Pixel *prev, float bx, float by, float &r,
                            float &g, float &b) const {
    float fy0 = std::floor(by);
    int y0 = static_cast<int>(fy0);
    int y1 = y0 + 1;

    float fx0 = std::floor(bx);
    int x0 = static_cast<int>(fx0);
    float fx = bx - fx0;
    float fy = by - fy0;

    // Producer must keep bx in [-W, 2W); fast_wrap corrects only a single step.
    // Derive x1 from the wrapped x0: a pre-wrap x0 == 2W-1 would form x1 == 2W,
    // tripping fast_wrap's x < 2W assert (and returning column W in release).
    x0 = fast_wrap(x0, W);
    int x1 = fast_wrap(x0 + 1, W);

    ::Pixel p00 = (y0 >= 0 && y0 < H) ? prev[y0 * W + x0] : ::Pixel(0, 0, 0);
    ::Pixel p10 = (y0 >= 0 && y0 < H) ? prev[y0 * W + x1] : ::Pixel(0, 0, 0);
    ::Pixel p01 = (y1 >= 0 && y1 < H) ? prev[y1 * W + x0] : ::Pixel(0, 0, 0);
    ::Pixel p11 = (y1 >= 0 && y1 < H) ? prev[y1 * W + x1] : ::Pixel(0, 0, 0);

    float w00 = (1.0f - fx) * (1.0f - fy);
    float w10 = fx * (1.0f - fy);
    float w01 = (1.0f - fx) * fy;
    float w11 = fx * fy;

    r = p00.r * w00 + p10.r * w10 + p01.r * w01 + p11.r * w11;
    g = p00.g * w00 + p10.g * w10 + p01.g * w01 + p11.g * w11;
    b = p00.b * w00 + p10.b * w10 + p01.b * w01 + p11.b * w11;
  }

  /** @brief Quantizes an unclamped [0, 65535]-scale channel to a Pixel
   *  component, round-to-nearest; NaN maps to the hi bound. */
  static uint16_t quantize16(float v) {
    // clamp guards the cast against overshoot and maps NaN to the hi bound.
    return static_cast<uint16_t>(hs::clamp(v, 0.0f, 65535.0f) + 0.5f);
  }

  /** @brief Coarse grid downsample the warp cache is sized for (the default
   *  Style's; every preset keeps it). Other values render uncached. */
  static constexpr int CACHE_DS = ::Feedback::Style{}.downsample;
  /** @brief Cell count of the cached coarse grid. */
  static constexpr int CACHE_CELLS = (W / CACHE_DS) * (H / CACHE_DS);

public:
  /** @brief Persistent bytes init_storage() reserves (two int16 warp fields). */
  static constexpr size_t STORAGE_BYTES = 2 * CACHE_CELLS * sizeof(int16_t);

private:
  /** @brief Inputs the coarse warp field is a pure function of (stock
   *  transforms only); equal keys make the cached field reusable. */
  struct WarpKey {
    ::Feedback::SpaceFn space_fn;
    const NoiseParams *noise;
    float amplitude, frequency, speed, scale, time;
    int cy_lo, cy_hi;
    bool operator==(const WarpKey &) const = default;
  };

  ::Feedback::Style *style_;  /**< Bound feedback Style (non-owning). */
  bool enabled_ = true;       /**< When false, flush() is skipped entirely. */
  WarpKey warp_key_{};        /**< Key the cached warp field was built for. */
  bool warp_cache_valid_ = false; /**< True once warp_dx_/warp_dy_ hold a field. */
  int16_t *warp_dx_ = nullptr; /**< Arena-owned cached column deltas. */
  int16_t *warp_dy_ = nullptr; /**< Arena-owned cached row deltas. */
};

/**
 * @brief Splits RGB into per-channel copies offset by 1/2/3 columns,
 * producing a chromatic-aberration fringe.
 */
template <int W> class ChromaticShift : public Is2D {
  // fast_wrap corrects only one ±W step, so the +1/+2/+3 column offsets stay in
  // a single wrap of [0,W) only for W >= 4.
  static_assert(W >= 4, "ChromaticShift requires W >= 4 for fast_wrap offsets");

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
   * @tparam PassFnT Downstream callback type; a forwarding reference so the
   * filter chain inlines with no per-point indirect call.
   * @param pass Downstream 2D callback.
   */
  template <typename PassFnT>
  void plot(float x, float y, const ::Pixel &c, float age, float alpha,
            PassFnT &&pass) {
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

    int xi = fast_wrap(static_cast<int>(std::round(x)), W);
    pass(static_cast<float>(fast_wrap(xi + 1, W)), y, r_col, age, alpha);
    pass(static_cast<float>(fast_wrap(xi + 2, W)), y, g_col, age, alpha);
    pass(static_cast<float>(fast_wrap(xi + 3, W)), y, b_col, age, alpha);
  }
};

} // namespace Pixel

} // namespace Filter

