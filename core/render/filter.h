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
#include "math/spherical_field.h"
#include "color/color.h"
#include "engine/static_circular_buffer.h"
#include "render/canvas.h"
#include "engine/concepts.h"
#include "engine/memory.h"
#include "engine/styles.h"

/**
 * @file filter.h
 * @brief The filter pipeline: Pipeline composition plus the World, Screen and
 * Pixel filter families.
 * @details The stage roster is a library surface, deliberately wider than the
 * set of stages the shipping effects instantiate; a stage with no current user
 * is composable inventory, not dead code.
 */

/** @brief Callback that forwards a 2D plot (x, y, pixel, age, alpha) downstream. */
using PassFn2D = FunctionRef<void(float, float, const Pixel &, float, float)>;
/** @brief Callback that forwards a 3D plot (vector, pixel, age, alpha) downstream. */
using PassFn3D = FunctionRef<void(const Vector &, const Pixel &, float, float)>;

/**
 * @brief Trait base every filter stage derives from, tagging its domain and
 * whether it carries state across frames.
 * @tparam Is2d True for a screen-space stage, false for a world-space one.
 * @tparam HasHistory True when the stage keeps state between frames.
 * @details `is_terminal`: writes the Canvas directly, so it must be the last
 * stage and its flush takes neither a trail nor a `pass` callback
 * (`flush(Canvas&, float)`). `terminal_replaces`: a terminal that
 * overwrites the whole frame (Feedback's opaque store), so no history-bearing
 * stage may precede it — its flush emissions would be clobbered — and the
 * effect must flush BEFORE the frame's plot() calls, not after as a
 * non-replacing terminal allows: at alpha >= 1 the flush writes every
 * destination pixel, erasing anything already plotted.
 * `emits_nonunit_world` /
 * `requires_unit_world_input`: a non-unit-emitting world stage must not precede
 * a unit-assuming one. `crosses_segments`: output can move between worker
 * segment bands, so the effect must render the full canvas. `reads_outside_band`:
 * the stage samples framebuffer pixels outside the display band, so stale pixels
 * outside that band must also be cleared. Both default to `has_history`
 * (fail-safe). `segment_margin`: how many pixels the stage's output can land
 * away from the plotted position, i.e. how far the render bounds must be padded
 * past the display band for a segment worker to still write the stage's taps.
 * A stage overrides any of these by redeclaring it.
 */
template <bool Is2d, bool HasHistory> struct FilterTraits {
  static constexpr int domain_rank = Is2d ? 1 : 0;
  static constexpr bool is_2d = Is2d;
  static constexpr bool has_history = HasHistory;
  static constexpr bool is_terminal = false;
  static constexpr bool terminal_replaces = false;
  static constexpr bool emits_nonunit_world = false;
  static constexpr bool requires_unit_world_input = false;
  static constexpr bool crosses_segments = has_history;
  static constexpr bool reads_outside_band = has_history;
  static constexpr int segment_margin = 0;
};

/** @brief Trait indicating a filter operates in 2D screen space. */
using Is2D = FilterTraits<true, false>;
/** @brief Trait indicating a filter operates in 3D world space. */
using Is3D = FilterTraits<false, false>;
/** @brief Trait indicating a 2D filter that maintains state/history. */
using Is2DWithHistory = FilterTraits<true, true>;
/** @brief Trait indicating a 3D filter that maintains state/history. */
using Is3DWithHistory = FilterTraits<false, true>;

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

/** @brief Whether a filter stage transforms edges before clip culling. */
template <typename Stage>
inline constexpr bool has_cull_edge =
    requires(const Stage &stage, const Vector &v, const Basis *pb) {
      stage.cull_edge(v, v, pb, PipelineCullEdgeProbe{});
    };

/**
 * @brief Terminal node of the pipeline (base case). Writes final pixels to the
 * Canvas.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
HS_O3_BEGIN
template <int W, int H> struct Pipeline<W, H> {
  static constexpr int domain_rank = 2;
  static constexpr bool is_2d = true;
  static constexpr bool is_terminal = false;
  static constexpr bool any_crosses_segments = false;
  static constexpr bool any_reads_outside_band = false;
  static constexpr int max_segment_margin = 0;
  static constexpr bool any_2d_history = false;
  static constexpr bool any_3d_history = false;
  static constexpr bool any_2d_trail_history = false;
  static constexpr bool any_terminal_history = false;
  /** @brief No stage re-emits clip-cull edges (see the recursive case). */
  static constexpr bool has_world_cull = false;
  /** @brief No stage runs in world space (see the recursive case). */
  static constexpr bool has_world_stage = false;

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
    // Producer must keep x in [-W, 2W); fast_wrap corrects only a single ±W offset.
    assert(x >= -W && x < 2 * W);
    if (!cv.clip().contains_y(y))
      return;
    int xi = fast_wrap(x, W);
    if (!cv.clip().contains_x(xi))
      return;
    plot_in_bounds(cv, xi, y, c, 0.0f, alpha);
  }

  /** @brief Writes a sample whose integer coordinates are already clip-tested. */
  void plot_in_bounds(Canvas &cv, int x, int y, const Pixel &c, float,
                      float alpha) {
    HS_PROFILE(filter_blend);
    assert(x >= 0 && x < W && cv.clip().contains_x(x));
    assert(cv.clip().contains_y(y));
    Pixel &dst = cv(x, y);
    if (alpha >= 1.0f) {
      dst = c;
      return;
    }
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
  void plot(Canvas &cv, float x, float y, const Pixel &c, float, float alpha) {
    // Non-finite coords make the int casts below UB and bypass the wrap.
    assert(std::isfinite(x) && std::isfinite(y));
    int xi = static_cast<int>(std::round(x));
    int yi = static_cast<int>(std::round(y));
    // fast_wrap corrects only a single ±W offset, so xi must land in [-W, 2W).
    assert(xi >= -W && xi < 2 * W);
    if (!cv.clip().contains_y(yi))
      return;
    xi = fast_wrap(xi, W);
    if (!cv.clip().contains_x(xi))
      return;
    plot_in_bounds(cv, xi, yi, c, 0.0f, alpha);
  }

  /**
   * @brief Projects a 3D point to screen space and writes it (sink).
   * @param cv Target canvas.
   * @param v World-space point on the unit sphere.
   * @param c Source color to blend in.
   * @param age Temporal age channel (frames).
   * @param alpha Blend alpha in [0, 1].
   */
  void plot(Canvas &cv, const Vector &v, const Pixel &c, float age,
            float alpha) {
    auto p = vector_to_pixel<W, H>(v);
    plot(cv, p.x, p.y, c, age, alpha);
  }

  /**
   * @brief Trail flush (base case: no stage to flush).
   * @tparam TrailFn ScreenTrailFn or WorldTrailFn.
   * @return Never returns; instantiation is a hard error.
   * @details Dependent-false guard: fires only when flush() is named on a
   * filterless pipeline. The recursion runs through flush_stages(), so this
   * overload is reachable from outside the pipeline only.
   */
  template <typename TrailFn> void flush(Canvas &, const TrailFn &, float) {
    static_assert(
        !sizeof(TrailFn *),
        "Wrong flush() domain: this Pipeline has no filter stages, so "
        "it carries no history and this overload emits nothing. Drop "
        "the flush() call, or add the history stage it was meant for.");
  }
  /**
   * @brief Terminal-stage flush (base case: no stage to flush).
   * @tparam T Unused; carries the dependent-false guard.
   * @return Never returns; instantiation is a hard error.
   */
  template <typename T = void> void flush(Canvas &, float) {
    static_assert(
        !sizeof(T *),
        "Wrong flush() domain: this Pipeline has no filter stages, so "
        "it has no terminal stage and this overload emits nothing. "
        "Drop the flush() call, or add Pixel::Feedback.");
  }
  /** @brief Terminates the recursive screen-trail flush walk. */
  void flush_stages(Canvas &, const ScreenTrailFn &, float) {}
  /** @brief Terminates the recursive world-trail flush walk. */
  void flush_stages(Canvas &, const WorldTrailFn &, float) {}
  /** @brief Terminates the recursive terminal-stage flush walk. */
  void flush_stages(Canvas &, float) {}

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
struct Pipeline<W, H, Head, Tail...> : private Head {
  using Next = Pipeline<W, H, Tail...>;
  Next next;

  static constexpr int domain_rank = Head::domain_rank;
  static constexpr bool is_2d = Head::is_2d;
  static constexpr bool is_terminal = Head::is_terminal || Next::is_terminal;
  static constexpr bool any_crosses_segments =
      Head::crosses_segments || Next::any_crosses_segments;
  static constexpr bool any_reads_outside_band =
      Head::reads_outside_band || Next::any_reads_outside_band;
  static constexpr int max_segment_margin =
      Head::segment_margin > Next::max_segment_margin
          ? Head::segment_margin
          : Next::max_segment_margin;

  static constexpr bool crosses_segments = any_crosses_segments;
  static constexpr bool reads_outside_band = any_reads_outside_band;
  static constexpr int segment_margin = max_segment_margin;

  static constexpr bool any_2d_history =
      (Head::has_history && Head::is_2d) || Next::any_2d_history;
  static constexpr bool any_3d_history =
      (Head::has_history && !Head::is_2d) || Next::any_3d_history;
  // A terminal stage composites into the Canvas itself, so its flush takes no
  // trail callback; only these stages need one.
  static constexpr bool any_2d_trail_history =
      (Head::has_history && Head::is_2d && !Head::is_terminal) ||
      Next::any_2d_trail_history;
  static constexpr bool any_terminal_history =
      (Head::has_history && Head::is_terminal) || Next::any_terminal_history;

  /**
   * @brief True when any stage overrides cull_edge (re-emits clip-cull edges
   *        under plot-time rotations), so a caller may not precompute per-point
   *        screen coordinates from the raw geometry.
   */
  static constexpr bool has_world_cull =
      has_cull_edge<Head> || Next::has_world_cull;

  /**
   * @brief True when any stage runs in world space, so a screen-space
   *        coordinate handed to plot() is lifted back through pixel_to_vector
   *        (not an inverse of vector_to_pixel) before that stage sees it, and a
   *        caller may not substitute precomputed screen coordinates for the
   *        point's world position.
   */
  static constexpr bool has_world_stage = !Head::is_2d || Next::has_world_stage;

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
      Head::plot(
          x, y, c, age, alpha,
          [&](float nx, float ny, const Pixel &nc, float nage, float nalpha) {
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
  void plot(Canvas &cv, const Vector &v, const Pixel &c, float age,
            float alpha) {
    if constexpr (!Head::is_2d) {
      Head::plot(v, c, age, alpha,
                 [&](const Vector &nv, const Pixel &nc, float nage,
                     float nalpha) { next.plot(cv, nv, nc, nage, nalpha); });
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
    if constexpr (has_cull_edge<Head>)
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
      !Head::has_history || !Head::is_2d || Head::is_terminal ||
          requires(Head h, Canvas &cv, const ScreenTrailFn &s, PassFn2D p) {
            h.flush(cv, s, 1.0f, p);
          },
      "2D history filter must define "
      "flush(Canvas&, const ScreenTrailFn&, float, PassFn2D)");
  static_assert(
      !Head::has_history || !Head::is_terminal || requires(Head h, Canvas &cv) {
        h.flush(cv, 1.0f);
      }, "terminal history filter must define flush(Canvas&, float)");

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
      Head::domain_rank <= Next::domain_rank,
      "Filter ordering: filter domains must be non-decreasing (World, then "
      "Screen, then Pixel). Reorder Pixel::* stages after Screen::* stages.");

  static_assert(
      !Head::emits_nonunit_world || (... && !Tail::requires_unit_world_input),
      "Filter ordering: a World stage that emits non-unit-length points "
      "(World::Trails) must not precede a World stage that requires unit-length "
      "input (World::Mobius / World::Hole). Reorder so the unit-assuming stage "
      "runs first, or renormalize the trail re-emission.");

  /**
   * @brief Flushes every 2D history stage in the pipeline.
   * @param cv Target canvas.
   * @param trailFn Callback producing trail color/alpha per screen point.
   * @param alpha Global blend alpha in [0, 1].
   */
  void flush(Canvas &cv, const ScreenTrailFn &trailFn, float alpha) {
    static_assert(
        any_2d_history,
        "Wrong flush() domain: this Pipeline has no 2D history stage, so the "
        "ScreenTrailFn overload emits nothing. Aging happens inside flush() — "
        "a 3D history stage (World::Trails) left unflushed fills its ring "
        "buffer to capacity and never decays. Pass a WorldTrailFn instead.");
    static_assert(
        !any_3d_history,
        "Incomplete flush(): this Pipeline also carries a 3D history stage "
        "(World::Trails) that this overload leaves unflushed, so its ring "
        "buffer fills to capacity and never decays. Pass both callbacks: "
        "flush(cv, worldTrailFn, screenTrailFn, alpha).");
    flush_stages(cv, trailFn, alpha);
  }

  /**
   * @brief Flushes every 3D history stage in the pipeline.
   * @param cv Target canvas.
   * @param trailFn Callback producing trail color/alpha per world point.
   * @param alpha Global blend alpha in [0, 1].
   */
  void flush(Canvas &cv, const WorldTrailFn &trailFn, float alpha) {
    static_assert(
        any_3d_history,
        "Wrong flush() domain: this Pipeline has no 3D history stage, so the "
        "WorldTrailFn overload emits nothing. Aging happens inside flush() — "
        "a 2D history stage (Screen::Trails, Pixel::Feedback) left unflushed "
        "never decays. Pass a ScreenTrailFn instead.");
    static_assert(
        !any_2d_history,
        "Incomplete flush(): this Pipeline also carries a 2D history stage "
        "(Screen::Trails, Pixel::Feedback) that this overload leaves "
        "unflushed, so it never decays. Pass both callbacks: "
        "flush(cv, worldTrailFn, screenTrailFn, alpha).");
    flush_stages(cv, trailFn, alpha);
  }

  /**
   * @brief Flushes both history domains of a pipeline that carries each.
   * @param cv Target canvas.
   * @param worldFn Callback producing trail color/alpha per world point.
   * @param screenFn Callback producing trail color/alpha per screen point.
   * @param alpha Global blend alpha in [0, 1].
   * @details The 3D pass runs first: World stages precede Screen stages, so its
   * re-emissions reach the 2D history stage before that stage emits and ages.
   */
  void flush(Canvas &cv, const WorldTrailFn &worldFn,
             const ScreenTrailFn &screenFn, float alpha) {
    static_assert(
        any_3d_history && any_2d_history,
        "Wrong flush() domain: this Pipeline carries history in only one "
        "domain, so one of these callbacks emits nothing. Pass the single "
        "callback that domain needs.");
    flush_stages(cv, worldFn, alpha);
    flush_stages(cv, screenFn, alpha);
  }

  /**
   * @brief Flushes every terminal history stage in the pipeline.
   * @param cv Target canvas.
   * @param alpha Global blend alpha in [0, 1].
   * @details For pipelines whose only history is terminal (Pixel::Feedback):
   * a terminal composites into the Canvas itself and needs no trail callback.
   */
  void flush(Canvas &cv, float alpha) {
    static_assert(
        any_terminal_history,
        "Wrong flush() domain: this Pipeline has no terminal history stage, so "
        "this overload emits nothing. Pass a ScreenTrailFn or WorldTrailFn.");
    static_assert(
        !any_2d_trail_history,
        "This Pipeline carries a trail-bearing 2D history stage that this "
        "overload would leave unflushed (and therefore undecayed). Pass a "
        "ScreenTrailFn instead — it flushes the terminal stage too.");
    flush_stages(cv, alpha);
  }

  /**
   * @brief Flushes 2D history for this stage, then recurses into the Tail.
   * @param cv Target canvas.
   * @param trailFn Callback producing trail color/alpha per screen point.
   * @param alpha Global blend alpha in [0, 1].
   * @details Only a 2D history-bearing Head emits; other Heads pass through.
   * Recursion target of flush(), below the domain assert — a Tail is free to
   * carry no history of its own.
   */
  void flush_stages(Canvas &cv, const ScreenTrailFn &trailFn, float alpha) {
    if constexpr (Head::has_history) {
      if constexpr (Head::is_2d) {
        if constexpr (Head::is_terminal) {
          Head::flush(cv, alpha);
        } else {
          Head::flush(cv, trailFn, alpha,
                      [&](auto... args) { next.plot(cv, args...); });
        }
      }
    }
    next.flush_stages(cv, trailFn, alpha);
  }

  /**
   * @brief Flushes a terminal history stage, then recurses into the Tail.
   * @param cv Target canvas.
   * @param alpha Global blend alpha in [0, 1].
   * @details Recursion target of flush(), below the domain asserts.
   */
  void flush_stages(Canvas &cv, float alpha) {
    if constexpr (Head::has_history && Head::is_terminal) {
      Head::flush(cv, alpha);
    }
    next.flush_stages(cv, alpha);
  }

  /**
   * @brief Flushes 3D history for this stage, then recurses into the Tail.
   * @param cv Target canvas.
   * @param trailFn Callback producing trail color/alpha per world point.
   * @param alpha Global blend alpha in [0, 1].
   * @details Only a 3D history-bearing Head emits; other Heads pass through.
   * Recursion target of flush(), below the domain assert.
   */
  void flush_stages(Canvas &cv, const WorldTrailFn &trailFn, float alpha) {
    if constexpr (Head::has_history) {
      if constexpr (!Head::is_2d) {
        Head::flush(trailFn, alpha,
                    [&](auto... args) { next.plot(cv, args...); });
      }
    }
    next.flush_stages(cv, trailFn, alpha);
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
class Orient : public Is3D {
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
  Orientation<>
      &orientation; /**< Live orientation source driving the rotation. */
};

/**
 * @brief Selects an orientation from a list based on the point's projection
 * onto an axis. Useful for slicing objects with different rotations.
 */
class OrientSlice : public Is3D {
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

  bool enabled; /**< When false, the filter passes points through unrotated. */

private:
  Vector axis; /**< Unit axis points are projected onto to select a slice. */
  std::span<const Orientation<>>
      orientations; /**< Candidate orientations indexed by projection. */
};

/**
 * @brief Creates a spherical hole by masking points within a radius.
 */
class Hole : public Is3D {
public:
  static constexpr bool requires_unit_world_input = true;
  /**
   * @brief Constructs a hole centered at @p origin with angular @p radius.
   * @param origin Center of the hole (unit vector).
   * @param radius Angular radius of the hole in radians.
   */
  Hole(const Vector &origin, float radius) : origin(origin), radius(radius) {}

  /** @brief Moves the hole center to a new unit vector. */
  void set_origin(const Vector &new_origin) { origin = new_origin; }

  /** @brief Changes the hole's angular radius in radians. */
  void set_radius(float new_radius) { radius = new_radius; }
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
    float d = angle_between(v, origin);
    if (d > radius)
      pass(v, color, age, alpha);
    else {
      float t = d / radius;
      pass(v, color * quintic_kernel(t), age, alpha);
    }
  }

private:
  Vector origin; /**< Center of the hole (unit vector). */
  float radius;  /**< Angular radius of the hole in radians. */
};

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
   * feeds make_rotation, so count == 0 cannot feed inf into it. The ceiling is
   * W because W copies already sit one equatorial pixel column apart; beyond
   * that they land on the same column.
   */
  Replicate(int count) { set_count(count); }

  /**
   * @brief Changes the number of evenly spaced copies.
   * @param new_count Desired copy count; clamped to [1, W].
   */
  void set_count(int new_count) {
    count = hs::clamp(new_count, 1, W);
    step = make_rotation(Y_AXIS, 2 * PI_F / count);
  }
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

  /**
   * @brief Re-emits a clip-cull edge under every copy's Y-axis rotation.
   * @tparam FwdFn Downstream cull continuation
   *         `bool(const Vector&, const Vector&, const Basis*)`.
   * @param a,b Edge endpoints in world space (pre-rotation).
   * @param pb Optional planar basis, rotated alongside the endpoints.
   * @param forward Tail-of-pipeline cull continuation.
   * @return True if any copy of the edge could intersect the clip band.
   * @details Mirrors plot()'s rotation loop so the cull sees the longitudes the
   *          copies are drawn at, not the source geometry's.
   */
  template <typename FwdFn>
  bool cull_edge(const Vector &a, const Vector &b, const Basis *pb,
                 FwdFn &&forward) const {
    if (forward(a, b, pb))
      return true;
    Vector ra = a, rb = b;
    Basis rp;
    if (pb)
      rp = *pb;
    for (int i = 1; i < count; i++) {
      ra = rotate(ra, step).normalized();
      rb = rotate(rb, step).normalized();
      if (pb) {
        rp = rotate(rp, step);
        if (forward(ra, rb, &rp))
          return true;
      } else if (forward(ra, rb, nullptr)) {
        return true;
      }
    }
    return false;
  }

private:
  int count;       /**< Number of copies emitted, in [1, W]. */
  Quaternion step; /**< Per-copy Y-axis rotation (2*pi / count). */
};

/**
 * @brief Replicates geometry onto the vertices of a solid.
 * @details Precomputes rotation quaternions from vertex[0] to each other vertex.
 * Every copy carries the source age unchanged (replication is spatial).
 */
template <int N> class VertexReplicate : public Is3D {
public:
  /**
   * @brief Builds from a vertex array, precomputing rotations vertices[0] → each.
   * @tparam VertexArray Indexable container of N unit vectors.
   * @param vertices Vertex positions; rotations map vertices[0] onto each vertex.
   */
  template <typename VertexArray> VertexReplicate(const VertexArray &vertices) {
    set_vertices(vertices);
  }

  /**
   * @brief Rebuilds the rotations from a new vertex array.
   * @tparam VertexArray Indexable container of N unit vectors.
   * @param vertices Vertex positions; rotations map vertices[0] onto each vertex.
   */
  template <typename VertexArray>
  void set_vertices(const VertexArray &vertices) {
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

  /**
   * @brief Re-emits a clip-cull edge under every stored vertex rotation.
   * @tparam FwdFn Downstream cull continuation
   *         `bool(const Vector&, const Vector&, const Basis*)`.
   * @param a,b Edge endpoints in world space (pre-rotation).
   * @param pb Optional planar basis, rotated alongside the endpoints.
   * @param forward Tail-of-pipeline cull continuation.
   * @return True if any vertex copy of the edge could intersect the clip band.
   * @details Mirrors plot(). The rotations move latitude, so culling by the
   *          un-rotated endpoints would drop copies the fan-out places inside a
   *          segment band (docs/segmented_stateful_effects_spec.md).
   */
  template <typename FwdFn>
  bool cull_edge(const Vector &a, const Vector &b, const Basis *pb,
                 FwdFn &&forward) const {
    for (int i = 0; i < N; ++i) {
      if (pb) {
        Basis rb = rotate(*pb, rotations[i]);
        if (forward(rotate(a, rotations[i]), rotate(b, rotations[i]), &rb))
          return true;
      } else if (forward(rotate(a, rotations[i]), rotate(b, rotations[i]),
                         nullptr)) {
        return true;
      }
    }
    return false;
  }

private:
  std::array<Quaternion, N>
      rotations; /**< Rotation from vertices[0] to each vertex. */
};

/**
 * @brief Applies a Mobius transformation to 3D points.
 */
class Mobius : public Is3D {
public:
  static constexpr bool requires_unit_world_input = true;
  /**
   * @brief The map is non-rigid, so no rotation-mirroring cull_edge can bound
   *        the image of an edge; the effect must render the full canvas.
   */
  static constexpr bool crosses_segments = true;
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
 * @details Its Screen:: namesake is not the same filter in another domain — it
 * gates seeding on alpha where this one seeds regardless, evicts at capacity by
 * moving the last live point into slot 0 where this one drops the ring's
 * logical head, and overrides crosses_segments to false where this one keeps
 * the fail-safe has_history default.
 */
template <int Capacity> class Trails : public Is3DWithHistory {
public:
  static constexpr bool emits_nonunit_world = true;
  static constexpr bool reads_outside_band = false;

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
   * @param alpha Blend alpha in [0, 1], forwarded unchanged and NOT gated: a
   * transparent sample still consumes a ring slot, unlike Screen::Trails::plot.
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
    HS_CHECK(items, "World::Trails needs init_storage() from effect init()");
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
  Item *items = nullptr; /**< Ring-buffer storage (arena-owned). */
  size_t head = 0, tail = 0,
         count = 0; /**< Ring-buffer head, tail, and live count. */
  int lifetime;     /**< Per-frame fade divisor in frames. */

  static constexpr float Q =
      32767.0f; /**< Quantization scale for unit-vector components. */
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
   * @brief Returns the i-th logical live item.
   * @param i Index into the live range [0, count).
   * @return Mutable reference to the buffered Item.
   */
  Item &at(size_t i) { return items[(head + i) % Capacity]; }
  /**
   * @brief Returns the i-th logical live item.
   * @param i Index into the live range [0, count).
   * @return Const reference to the buffered Item.
   */
  const Item &at(size_t i) const { return items[(head + i) % Capacity]; }

  /**
   * @brief Appends an item, evicting a live item of arbitrary age at capacity.
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

  /** @brief Drops the logical head, whose age is arbitrary after compaction. */
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
  static constexpr int segment_margin = 1;

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
    // fast_wrap below corrects only a single ±W offset on floorf(x).
    assert(x >= -W && x < 2 * W);
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
  static constexpr bool any_reads_outside_band = false;
  static constexpr int segment_margin = 1;
  static constexpr int max_segment_margin = segment_margin;
  static constexpr bool has_world_cull = false;
  static constexpr bool has_world_stage = false;
  static constexpr bool direct_raster_path = true;

  /** @brief Caches the current frame's framebuffer and clip bounds. */
  void prepare(Canvas &cv) {
    base = cv.data();
    const ClipRegion &cr = cv.clip();
    clip_stamp = cr;
    const ClipRegion::XClip xc = cr.x_clip();
    for (int x = 0; x < W; ++x)
      x_visible[x] = !xc.clipped(x);
    const int y_start = cr.render_y_start();
    const int y_end = cr.render_y_end();
    for (int y = 0; y < H; ++y)
      y_visible[y] = y >= y_start && y < y_end;
  }

  /**
   * @brief True when the cached framebuffer base and clip state still match @p cv.
   * @param cv Canvas the caller is about to draw into.
   * @details The base pointer alone repeats every other frame when the canvas
   * double-buffers, so the clip bounds are stamped alongside it.
   */
  bool prepared_for(Canvas &cv) const {
    return base == cv.data() && clip_stamp == cv.clip();
  }

  /** @brief Splats one screen-space sample directly into the Canvas. */
  void plot(Canvas &cv, float x, float y, const Pixel &c, float age,
            float alpha) {
    assert(age >= 0.0f && alpha >= 0.0f);
    assert(prepared_for(cv));
    // fast_wrap below corrects only a single ±W offset on floorf(x).
    assert(x >= -W && x < 2 * W);
    (void)age;
    (void)cv;

    HS_MSP_STALL_START(aa_start);
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

    const bool x0_ok = x_visible[x0];
    const bool x1_ok = x_visible[x1];
    const bool y0_ok = y0_physical && y_visible[y0];
    const bool y1_ok = y1_physical && y_visible[y1];
    const uint8_t tap_mask =
        static_cast<uint8_t>((y0_ok && x0_ok && v00 > TAP_CUTOFF) |
                             ((y0_ok && x1_ok && v10 > TAP_CUTOFF) << 1) |
                             ((y1_ok && x0_ok && v01 > TAP_CUTOFF) << 2) |
                             ((y1_ok && x1_ok && v11 > TAP_CUTOFF) << 3));
    HS_MSP_STALL_STOP(aa_weights, aa_start);

#ifdef HS_PROFILE_MINDSPLATTER_COUNTS
    const unsigned tap_count = static_cast<unsigned>(tap_mask & 0x01) +
                               static_cast<unsigned>((tap_mask >> 1) & 0x01) +
                               static_cast<unsigned>((tap_mask >> 2) & 0x01) +
                               static_cast<unsigned>((tap_mask >> 3) & 0x01);
    ++hs::g_mindsplatter_counts.aa_tap_masks[tap_count];
#endif

    HS_MSP_STALL_START(blend_start);
    if (tap_mask == 0x0f) {
      HS_MSP_COUNT(interior_splats);
      blend_four(base + y0 * W, base + y1 * W, x0, x1, c, alpha, v00, v10, v01,
                 v11);
      HS_MSP_STALL_STOP(framebuffer_blend, blend_start);
      return;
    }
    HS_MSP_COUNT(clip_boundary_splats);
    blend_masked(base, x0, x1, y0, y1, c, alpha, v00, v10, v01, v11, tap_mask);
    HS_MSP_STALL_STOP(framebuffer_blend, blend_start);
  }

  /** @brief Integer-coordinate overload matching a filtered Pipeline. */
  void plot(Canvas &cv, int x, int y, const Pixel &c, float age, float alpha) {
    plot(cv, static_cast<float>(x), static_cast<float>(y), c, age, alpha);
  }

  /** @brief Projects a world point, then applies the direct screen-space splat. */
  void plot(Canvas &cv, const Vector &v, const Pixel &c, float age,
            float alpha) {
    HS_MSP_STALL_START(projection_start);
    const PixelCoords p = vector_to_pixel<W, H>(v);
    HS_MSP_STALL_STOP(projection, projection_start);
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
  Pixel *base = nullptr;
  ClipRegion clip_stamp{};
  // Zero-init: a plot() before prepare() is a masked no-op, not a read of
  // indeterminate bytes through a null base.
  std::array<uint8_t, W> x_visible{};
  std::array<uint8_t, H> y_visible{};

  static __attribute__((always_inline)) uint16_t tap_alpha_q16(float alpha,
                                                               float weight) {
    return static_cast<uint16_t>(
        hs::clamp(alpha * weight * 65535.0f + 0.5f, 0.0f, 65535.0f));
  }

  static __attribute__((always_inline)) void
  blend_four(Pixel *row0, Pixel *row1, int x0, int x1, const Pixel &src,
             float alpha, float v00, float v10, float v01, float v11) {
    const uint16_t a00 = tap_alpha_q16(alpha, v00);
    const uint16_t a10 = tap_alpha_q16(alpha, v10);
    const uint16_t a01 = tap_alpha_q16(alpha, v01);
    const uint16_t a11 = tap_alpha_q16(alpha, v11);
    blend_tap(row0 + x0, src, a00);
    blend_tap(row0 + x1, src, a10);
    blend_tap(row1 + x0, src, a01);
    blend_tap(row1 + x1, src, a11);
  }

  HS_NOINLINE_NOCLONE static void blend_masked(Pixel *base, int x0, int x1,
                                               int y0, int y1, const Pixel &src,
                                               float alpha, float v00,
                                               float v10, float v01, float v11,
                                               uint8_t tap_mask) {
    if (tap_mask & 0x01)
      blend_tap(base + y0 * W + x0, src, tap_alpha_q16(alpha, v00));
    if (tap_mask & 0x02)
      blend_tap(base + y0 * W + x1, src, tap_alpha_q16(alpha, v10));
    if (tap_mask & 0x04)
      blend_tap(base + y1 * W + x0, src, tap_alpha_q16(alpha, v01));
    if (tap_mask & 0x08)
      blend_tap(base + y1 * W + x1, src, tap_alpha_q16(alpha, v11));
  }

  static __attribute__((always_inline)) void
  blend_tap(Pixel *dst, const Pixel &src, uint16_t alpha_q16) {
    *dst = dst->lerp16(src, alpha_q16);
  }
};
HS_O3_END

/**
 * @brief Manages 2D screen-space trails.
 * @details Its World:: namesake is not the same filter in another domain — it
 * seeds regardless of alpha where this one gates, evicts at capacity by
 * dropping its ring's logical head where this one moves the last live point
 * into slot 0, and keeps the fail-safe has_history crosses_segments default
 * where this one overrides it to false.
 */
template <int MAX_PIXELS = 1024> class Trails : public Is2DWithHistory {
public:
  // Trail points are seeded from and re-emitted into the same band, so they
  // never sample a neighbor segment.
  static constexpr bool crosses_segments = false;
  static constexpr bool reads_outside_band = false;

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
    points = arena.allocate_n<DecayPixel>(MAX_PIXELS);
    num_pixels = 0;
  }

  /**
   * @brief Forwards the current sample and seeds a fading screen trail point.
   * @param x Column coordinate in pixels.
   * @param y Row coordinate in pixels.
   * @param color Source color, forwarded unchanged this frame.
   * @param age Incoming age (frames); ttl = lifetime - age, seeded only if positive.
   * @param alpha Blend alpha in [0, 1]; samples with alpha <= 0.001 seed no
   * trail point but are still forwarded downstream. World::Trails::plot has no
   * such gate and seeds regardless of alpha.
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
    if (ttl > 0.0f && points) {
      if (num_pixels == MAX_PIXELS) {
        // At capacity: O(1) drop of slot 0. Both trail domains may evict a point
        // of arbitrary age after swap-removal compacts an expired mid-buffer
        // slot; per-point ttl fade absorbs the transient.
        num_pixels--;
        points[0] = points[num_pixels];
      }
      points[num_pixels++] = {x, y, ttl};
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
    HS_CHECK(points, "Screen::Trails needs init_storage() from effect init()");
    for (int i = 0; i < num_pixels; ++i) {
      float t = hs::clamp(1.0f - (points[i].ttl / lifetime), 0.0f, 1.0f);
      Color4 color = trailFn(points[i].x, points[i].y, t);
      if (color.alpha > 0.001f) {
        float age = lifetime - points[i].ttl;
        if (age < 0.0f)
          age = 0.0f;
        pass(points[i].x, points[i].y, color.color, age, alpha * color.alpha);
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
      if (--points[i].ttl <= 0.0f) {
        points[i] = points[--num_pixels];
        i--;
      }
    }
  }

private:
  /** @brief One screen trail point: position plus remaining lifetime. */
  struct DecayPixel {
    float x, y, ttl; /**< Pixel position and remaining lifetime in frames. */
  };
  int lifetime;                 /**< Per-frame fade divisor in frames. */
  DecayPixel *points = nullptr; /**< Arena-owned array of live trail points. */
  int num_pixels = 0;           /**< Number of live points in points. */
};

/**
 * @brief Applies a variable 3x3 Gaussian blur.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H> class Blur : public Is2D {
public:
  static constexpr int segment_margin = 1;

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
 * @details The Style's spatial warp is computed on a spherical latitude-ring
 * field, then interpolated within and between rings.
 * flush() iterates the full pixel grid within the active clip band. TERMINAL:
 * flush() composites directly into the Canvas rather than re-emitting downstream,
 * so it must be the last Pipeline stage. The effect must call flush() BEFORE
 * the frame's plot() calls (see `terminal_replaces`); flushing last, as a
 * non-replacing terminal permits, blanks the frame at alpha >= 1.
 */
template <int W, int H> class Feedback : public Is2DWithHistory {
public:
  static constexpr int domain_rank = 2;
  /** @brief Marks this as terminal: flush() writes the Canvas directly. */
  static constexpr bool is_terminal = true;
  /** @brief Opaque store owns the frame: no history stage may precede it. */
  static constexpr bool terminal_replaces = true;

  // Covers only the default-constructed Style; a runtime-swapped style's
  // downsample is validated in flush() (the HS_CHECK below).
  static_assert(
      ::Feedback::Style{}.downsample > 0 &&
          W % ::Feedback::Style{}.downsample == 0,
      "Feedback<W,H>: default style downsample must be > 0 and divide "
      "W");

  /**
   * @brief Binds the filter to a live feedback Style.
   * @param style Style supplying the spatial warp and color transforms.
   */
  explicit Feedback(::Feedback::Style &style) : feedback_style(&style) {}

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
  void flush(Canvas &cv, float alpha) {
    if (!enabled)
      return;

    const Animation::NoiseParams *bound = feedback_style->noise;
    HS_CHECK(!bound || (bound->amplitude == feedback_style->amplitude &&
                        bound->frequency == feedback_style->frequency &&
                        bound->speed == feedback_style->speed &&
                        bound->scale == feedback_style->scale),
             "feedback style scalars never reached the bound NoiseParams; call "
             "Style::sync_noise() after changing them");

    const CoarseGrid grid = make_coarse_grid(cv);

    {
      HS_PROFILE(feedback_litscan);
      if (!any_pixel_lit(cv))
        return;
    }

    const RenderBand band = make_render_band(cv.clip(), grid);

    ScratchScope scope(scratch_arena_a);
    WarpField warp = select_warp_field(scope.get_arena(), grid, band);
    populate_warp_field(grid, band, warp);
    ::Pixel *filtered_row = scope.get_arena().allocate_n<::Pixel>(W);
    composite_previous_frame(cv, alpha, grid, band, warp, filtered_row);
  }

private:
  using SphereField = hs::SphericalFieldLayout<W, H>;

  struct CoarseGrid {
    int downsample;
    SphereField field;
    int columns;
    int field_rows;
  };

  struct RenderBand {
    int y_begin;
    int y_end;
    int field_y_begin;
    int field_y_end;
    ClipRegion::XClip x_clip;
    std::bitset<W> coarse_columns_used;
  };

  struct WarpControl {
    int16_t x;
    int16_t y;
  };

  struct WarpField {
    int16_t *x_offsets;
    int16_t *y_offsets;
    WarpControl *controls;
    bool needs_population;
  };

  struct PixelAccumulator {
    uint32_t r = 0;
    uint32_t g = 0;
    uint32_t b = 0;

    void add(const ::Pixel &pixel) {
      r += pixel.r;
      g += pixel.g;
      b += pixel.b;
    }

    void remove(const ::Pixel &pixel) {
      r -= pixel.r;
      g -= pixel.g;
      b -= pixel.b;
    }

    ::Pixel average(int width) const {
      const uint32_t round = static_cast<uint32_t>(width / 2);
      return ::Pixel(static_cast<uint16_t>((r + round) / width),
                     static_cast<uint16_t>((g + round) / width),
                     static_cast<uint16_t>((b + round) / width));
    }
  };

  struct ColumnRun {
    int begin;
    int end;
  };

  struct ColumnRuns {
    ColumnRun items[2];
    int count;
  };

  __attribute__((always_inline)) CoarseGrid
  make_coarse_grid(const Canvas &cv) const {
    const int downsample = feedback_style->downsample;
    HS_CHECK(downsample > 0 && W % downsample == 0,
             "feedback downsample %d must be > 0 and divide width %d",
             downsample, W);
    HS_CHECK(cv.width() == W,
             "feedback canvas width %d must equal template W %d", cv.width(),
             W);
    HS_CHECK(cv.height() == H,
             "feedback canvas height %d must equal template H %d", cv.height(),
             H);
    const int columns = W / downsample;
    const int south_infill = std::max(downsample - hs::H_OFFSET, 0);
    const SphereField field(downsample, downsample, south_infill, columns);
    return {downsample, field, columns, field.ring_count()};
  }

  static __attribute__((always_inline)) RenderBand
  make_render_band(const ClipRegion &clip, const CoarseGrid &grid) {
    RenderBand band{};
    band.y_begin = clip.render_y_start();
    band.y_end = clip.render_y_end();
    band.x_clip = clip.x_clip();

    band.field_y_begin = grid.field.ring_index_at_or_before(band.y_begin);
    band.field_y_end = grid.field.ring_index_at_or_after(band.y_end - 1);
    HS_CHECK(band.field_y_end >= band.field_y_begin,
             "feedback field band inverted: [%d,%d]", band.field_y_begin,
             band.field_y_end);

    if (band.x_clip.active) {
      for (int x = 0; x < W; ++x) {
        if (band.x_clip.clipped(x))
          continue;
        const int left = x / grid.downsample;
        const int right = (left + 1 < grid.columns) ? left + 1 : 0;
        band.coarse_columns_used[left] = true;
        band.coarse_columns_used[right] = true;
      }
    }
    return band;
  }

  static __attribute__((always_inline)) ColumnRuns
  make_column_runs(const ClipRegion::XClip &clip) {
    ColumnRuns runs{};
    if (!clip.active) {
      runs.items[runs.count++] = {0, W};
    } else if (clip.wrap) {
      if (clip.re > 0)
        runs.items[runs.count++] = {0, clip.re};
      if (clip.rs < W)
        runs.items[runs.count++] = {clip.rs, W};
    } else {
      runs.items[runs.count++] = {clip.rs, clip.re};
    }
    return runs;
  }

  __attribute__((always_inline)) WarpField select_warp_field(
      Arena &scratch, const CoarseGrid &grid, const RenderBand &band) {
    const bool stock_transform =
        feedback_style->space_fn == &::Feedback::noise_warp ||
        feedback_style->space_fn == &::Feedback::melt_warp;
    const bool cacheable = cached_warp_x && !band.x_clip.active &&
                           grid.downsample == CACHE_DOWNSAMPLE &&
                           stock_transform;

    const Animation::NoiseParams *noise = feedback_style->noise;
    const WarpKey key{feedback_style->space_fn,
                      noise,
                      feedback_style->amplitude,
                      feedback_style->frequency,
                      feedback_style->speed,
                      feedback_style->scale,
                      noise ? noise->time * noise->speed : 0.0f,
                      band.field_y_begin,
                      band.field_y_end};

    if (!cacheable) {
      const int cells = grid.field_rows * grid.columns;
      return {scratch.allocate_n<int16_t>(cells),
              scratch.allocate_n<int16_t>(cells),
              scratch.allocate_n<WarpControl>(grid.field.sample_count()), true};
    }

    const bool needs_population = !(warp_cache_valid && key == cached_warp_key);
    cached_warp_key = key;
    warp_cache_valid = true;
    return {cached_warp_x, cached_warp_y,
            needs_population
                ? scratch.allocate_n<WarpControl>(grid.field.sample_count())
                : nullptr,
            needs_population};
  }

  __attribute__((always_inline)) void
  populate_warp_field(const CoarseGrid &grid, const RenderBand &band,
                      const WarpField &warp) const {
    HS_PROFILE(feedback_populate);
    if (!warp.needs_population)
      return;

    hs::SphericalField<WarpControl, W, H> compact(warp.controls, grid.field);
    compact.populate(
        band.field_y_begin, band.field_y_end,
        [&](const Vector &position,
            const typename SphereField::Coordinates &point) {
          Vector distorted;
          {
            HS_PROFILE_DEEP(fb_pop_warp);
            distorted = feedback_style->space_fn(position, *feedback_style);
          }
          HS_PROFILE_DEEP(fb_pop_project);
          // Both ends go through project() so its fast-trig error cancels; a
          // pole row's single backing vector has no azimuth to round-trip.
          // Rings stop at y = H - 1, which is the south pole only when
          // H_OFFSET == 0; the device's sub-pole rows carry no field samples.
          // h_offset_renorm_check compiles this file with the device H_OFFSET.
          bool pole_row = point.y == 0.0f;
          if constexpr (hs::H_OFFSET == 0) {
            constexpr float SOUTH_POLE_ROW = H + hs::H_OFFSET - 1;
            pole_row = pole_row || point.y == SOUTH_POLE_ROW;
          }
          const auto projected = grid.field.project(distorted);
          const auto origin = pole_row ? point : grid.field.project(position);
          float x_offset = projected.x - origin.x;
          const float y_offset = projected.y - origin.y;
          if (x_offset > W * 0.5f)
            x_offset -= W;
          else if (x_offset < -W * 0.5f)
            x_offset += W;
          return WarpControl{static_cast<int16_t>(hs::clamp(
                                 x_offset * WARP_SCALE, -32767.0f, 32767.0f)),
                             static_cast<int16_t>(hs::clamp(
                                 y_offset * WARP_SCALE, -32767.0f, 32767.0f))};
        });

    constexpr float WRAP_PERIOD = static_cast<float>(W) * WARP_SCALE;
    constexpr float HALF_WRAP_PERIOD = WRAP_PERIOD * 0.5f;
    {
      HS_PROFILE_DEEP(fb_pop_expand);
      auto ring = grid.field.ring(band.field_y_begin);
      for (int field_y = band.field_y_begin; field_y <= band.field_y_end;
           ++field_y, ring = grid.field.next_ring(ring)) {
        for (int coarse_x = 0; coarse_x < grid.columns; ++coarse_x) {
          if (band.x_clip.active && !band.coarse_columns_used[coarse_x])
            continue;
          const int x = coarse_x * grid.downsample;
          const auto longitude = grid.field.longitude_bounded(ring, x);
          const WarpControl a = warp.controls[longitude.left];
          WarpControl b = warp.controls[longitude.right];
          float bx = b.x;
          bx += (bx - a.x > HALF_WRAP_PERIOD)
                    ? -WRAP_PERIOD
                    : (bx - a.x < -HALF_WRAP_PERIOD ? WRAP_PERIOD : 0.0f);
          const int index = field_y * grid.columns + coarse_x;
          // Keep the stored offset canonical: the seam correction can lift the
          // interpolant a full period out, past both the int16_t range and the
          // single-step wrap sample_bilinear contracts for.
          const float offset_x =
              hs::lerp(static_cast<float>(a.x), bx, longitude.mix);
          warp.x_offsets[index] = static_cast<int16_t>(
              offset_x > HALF_WRAP_PERIOD
                  ? offset_x - WRAP_PERIOD
                  : (offset_x < -HALF_WRAP_PERIOD ? offset_x + WRAP_PERIOD
                                                  : offset_x));
          warp.y_offsets[index] = static_cast<int16_t>(hs::lerp(
              static_cast<float>(a.y), static_cast<float>(b.y), longitude.mix));
        }
      }
    }
  }

  __attribute__((always_inline)) void
  composite_previous_frame(Canvas &cv, float alpha, const CoarseGrid &grid,
                           const RenderBand &band, const WarpField &warp,
                           ::Pixel *filtered_row) {
    const int downsample = grid.downsample;
    const int coarse_columns = grid.columns;
    const int row_begin = band.y_begin;
    const int row_end = band.y_end;
    const int16_t *x_offsets = warp.x_offsets;
    const int16_t *y_offsets = warp.y_offsets;
    constexpr float INVERSE_WARP_SCALE = 1.0f / WARP_SCALE;
    const float inverse_downsample = 1.0f / grid.downsample;
    const float fade = feedback_style->fade;
    HS_CHECK(std::isfinite(fade) && fade >= 0.0f,
             "feedback fade %f must be finite and non-negative", fade);
    feedback_style->sync_hue();
    const bool black_skips_color =
        feedback_style->color_fn == &::Feedback::hue_fade;
    // Round-to-nearest fade otherwise makes sub-50 channels immortal at
    // FADE_MAX.
    constexpr float NEAR_BLACK = 64.0f;
    const auto blend = blend_alpha(alpha);
    const bool opaque = alpha >= 1.0f;
    constexpr float WRAP_PERIOD = static_cast<float>(W) * WARP_SCALE;
    constexpr float HALF_WRAP_PERIOD = WRAP_PERIOD * 0.5f;
    const ::Pixel *previous = cv.prev_data();
    ::Pixel *current = cv.data();
    ::Pixel poles[SphereField::POLE_COUNT];
    poles[0] = select_pole_sample(previous);
    if constexpr (SphereField::POLE_COUNT == 2)
      poles[1] = select_pole_sample(previous + (H - 1) * W);
    const ColumnRuns runs = make_column_runs(band.x_clip);
    int field_y0 = band.field_y_begin;
    int field_y1 = field_y0 + (field_y0 < band.field_y_end ? 1 : 0);
    const auto control_ring0 = grid.field.ring(field_y0);
    auto control_ring1 = control_ring0;
    if (field_y1 > field_y0)
      control_ring1 = grid.field.next_ring(control_ring0);
    int control_y0 = control_ring0.y;
    int control_y1 = control_ring1.y;
    auto composite_pixels = [&](auto &&transform_pixel, auto &&transform_pair,
                                auto pair_pixels) {
      constexpr bool PAIR_PIXELS = decltype(pair_pixels)::value;
      for (int y = row_begin; y < row_end; ++y) {
        const int row = y * W;
        // The infill bands are dense on both sides, but the southern one backs
        // a pole only when H_OFFSET is 0; the device's last row is mid-latitude.
        bool infill_band = y > 0 && y < downsample;
        if constexpr (hs::H_OFFSET == 0)
          infill_band = infill_band || (y >= H - downsample && y < H - 1);
        const bool filter_output = !band.x_clip.active && infill_band &&
                                   grid.field.longitude_filter_width(y) > 1;
        const bool defer_filter = filter_output && !opaque;
        ::Pixel *output = defer_filter ? filtered_row : current + row;
        while (y > control_y1) {
          field_y0 = field_y1;
          control_y0 = control_y1;
          if (field_y1 < grid.field_rows - 1) {
            ++field_y1;
            control_ring1 = grid.field.next_ring(control_ring1);
            control_y1 = control_ring1.y;
          }
        }
        // Interpolating outside the populated band silently corrupts pixels.
        HS_CHECK(field_y0 >= band.field_y_begin && field_y1 <= band.field_y_end,
                 "feedback warp row %d outside populated band [%d,%d]",
                 field_y1, band.field_y_begin, band.field_y_end);
        const int control_height = control_y1 - control_y0;
        const float fy =
            control_height > 0
                ? static_cast<float>(y - control_y0) / control_height
                : 0.0f;
        const float wy0 = 1.0f - fy, wy1 = fy;
        const int row0 = field_y0 * coarse_columns;
        const int row1 = field_y1 * coarse_columns;

        for (int r = 0; r < runs.count; ++r) {
          const int xs = runs.items[r].begin;
          const int xe = runs.items[r].end;
          int cx0 = xs / downsample;
          int sub = xs - cx0 * downsample;
          float leftx = 0.0f, slopex = 0.0f;
          float lefty = 0.0f, slopey = 0.0f;
          auto cell = [&]() {
            HS_PROFILE_DEEP(fb_comp_cell);
            const int cx1 = (cx0 + 1 < coarse_columns) ? cx0 + 1 : 0;
            const int i00 = row0 + cx0, i10 = row0 + cx1;
            const int i01 = row1 + cx0, i11 = row1 + cx1;
            float d00 = x_offsets[i00], d10 = x_offsets[i10];
            float d01 = x_offsets[i01], d11 = x_offsets[i11];
            d10 += (d10 - d00 > HALF_WRAP_PERIOD)
                       ? -WRAP_PERIOD
                       : (d10 - d00 < -HALF_WRAP_PERIOD ? WRAP_PERIOD : 0.0f);
            d01 += (d01 - d00 > HALF_WRAP_PERIOD)
                       ? -WRAP_PERIOD
                       : (d01 - d00 < -HALF_WRAP_PERIOD ? WRAP_PERIOD : 0.0f);
            d11 += (d11 - d00 > HALF_WRAP_PERIOD)
                       ? -WRAP_PERIOD
                       : (d11 - d00 < -HALF_WRAP_PERIOD ? WRAP_PERIOD : 0.0f);
            leftx = (d00 * wy0 + d01 * wy1) * INVERSE_WARP_SCALE;
            slopex = (d10 * wy0 + d11 * wy1) * INVERSE_WARP_SCALE - leftx;
            lefty = (y_offsets[i00] * wy0 + y_offsets[i01] * wy1) *
                    INVERSE_WARP_SCALE;
            slopey = (y_offsets[i10] * wy0 + y_offsets[i11] * wy1) *
                         INVERSE_WARP_SCALE -
                     lefty;
          };
          if (sub != 0)
            cell();

          for (int x = xs; x < xe;) {
            if (sub == 0)
              cell();

            if constexpr (PAIR_PIXELS) {
              if (downsample - sub >= 2 && xe - x >= 2) {
                const float fx0 = sub * inverse_downsample;
                const float fx1 = (sub + 1) * inverse_downsample;
                const float ddx0 = leftx + slopex * fx0;
                const float ddy0 = lefty + slopey * fx0;
                const float ddx1 = leftx + slopex * fx1;
                const float ddy1 = lefty + slopey * fx1;

                float sr0, sg0, sb0, sr1, sg1, sb1;
                {
                  HS_PROFILE_DEEP(fb_comp_sample);
                  sample_bilinear_prev(grid.field, previous, poles, x + ddx0,
                                       y + ddy0, sr0, sg0, sb0);
                  sample_bilinear_prev(grid.field, previous, poles,
                                       x + 1 + ddx1, y + ddy1, sr1, sg1, sb1);
                }
                ::Pixel p0(0, 0, 0), p1(0, 0, 0);
                {
                  HS_PROFILE_DEEP(fb_comp_color);
                  transform_pair(sr0, sg0, sb0, sr1, sg1, sb1, p0, p1);
                }
                // Keep both lanes on the paired transform path.
                const bool black0 = black_skips_color && sr0 < NEAR_BLACK &&
                                    sg0 < NEAR_BLACK && sb0 < NEAR_BLACK;
                const bool black1 = black_skips_color && sr1 < NEAR_BLACK &&
                                    sg1 < NEAR_BLACK && sb1 < NEAR_BLACK;
                p0 = black0 ? ::Pixel(0, 0, 0) : p0;
                p1 = black1 ? ::Pixel(0, 0, 0) : p1;

                HS_PROFILE_DEEP(fb_comp_write);
                ::Pixel &dst0 = output[x];
                dst0 = (opaque || defer_filter) ? p0 : blend(dst0, p0);
                ::Pixel &dst1 = output[x + 1];
                dst1 = (opaque || defer_filter) ? p1 : blend(dst1, p1);

                x += 2;
                sub += 2;
                if (sub == downsample) {
                  sub = 0;
                  ++cx0;
                }
                continue;
              }
            }

            const float fx = sub * inverse_downsample;
            const float ddx = leftx + slopex * fx;
            const float ddy = lefty + slopey * fx;

            float sr, sg, sb;
            {
              HS_PROFILE_DEEP(fb_comp_sample);
              sample_bilinear_prev(grid.field, previous, poles, x + ddx,
                                   y + ddy, sr, sg, sb);
            }
            ::Pixel p(0, 0, 0);
            if (!(black_skips_color && sr < NEAR_BLACK && sg < NEAR_BLACK &&
                  sb < NEAR_BLACK)) {
              HS_PROFILE_DEEP(fb_comp_color);
              p = transform_pixel(sr, sg, sb);
            }

            // Black must overwrite the stale double-buffer frame.
            HS_PROFILE_DEEP(fb_comp_write);
            ::Pixel &dst = output[x];
            dst = (opaque || defer_filter) ? p : blend(dst, p);

            ++x;
            if (++sub == downsample) {
              sub = 0;
              ++cx0;
            }
          }
        }
        if (filter_output) {
          HS_PROFILE_DEEP(fb_comp_filter);
          if (!defer_filter)
            std::copy_n(current + row, W, filtered_row);
          grid.field.template reconstruct_longitude_row<PixelAccumulator>(
              filtered_row, y, [&](int x, const ::Pixel &pixel) {
                ::Pixel &dst = current[row + x];
                dst = opaque ? pixel : blend(dst, pixel);
              });
        }
      }
    };
    dispatch_color_transform(fade, composite_pixels);
  }

  template <typename CompositeFnT>
  __attribute__((always_inline)) void
  dispatch_color_transform(float fade, CompositeFnT &&composite_pixels) {
    auto composite_scalar = [&](auto &&transform_pixel) {
      composite_pixels(transform_pixel, transform_pixel, std::false_type{});
    };

    HS_PROFILE(feedback_composite);
    const bool hue_identity =
        feedback_style->hue_ca == 1.0f && feedback_style->hue_sa == 0.0f;
    if (feedback_style->color_fn == &::Feedback::hue_fade && !hue_identity) {
      float k[9];
      const float sc = fast_cbrt(fade * (1.0f / 65535.0f));
      for (int i = 0; i < 9; ++i)
        k[i] = feedback_style->hue_k[i] * sc;
      composite_pixels(
          [&](float r, float g, float b) {
            return ::Feedback::hue_fade_apply(k, r, g, b);
          },
          [&](float r0, float g0, float b0, float r1, float g1, float b1,
              ::Pixel &p0, ::Pixel &p1) {
            ::Feedback::hue_fade_apply2(k, r0, g0, b0, r1, g1, b1, p0, p1);
          },
          std::true_type{});
    } else if (feedback_style->color_fn == &::Feedback::hue_fade) {
      auto plain = [&](float r, float g, float b) {
        return ::Pixel(quantize16(r * fade), quantize16(g * fade),
                       quantize16(b * fade));
      };
      composite_scalar(plain);
    } else {
      auto general = [&](float r, float g, float b) {
        return feedback_style->color_fn(
            ::Pixel(quantize16(r), quantize16(g), quantize16(b)), fade,
            *feedback_style);
      };
      composite_scalar(general);
    }
  }
  HS_O3_END

public:
  /**
   * @brief Enables or disables feedback.
   * @param value When false, flush() is skipped entirely.
   */
  void set_enabled(bool value) { enabled = value; }

  /**
   * @brief Allocates the warp-field cache from the persistent arena.
   * @param arena Persistent arena supplying 2 * CACHE_CELLS int16 slots.
   * @details Must be called from effect init(), not the constructor (arenas
   * aren't ready yet), and again after any compaction that resets the arena —
   * the cache is derived data, so it just re-populates on the next flush.
   * Without storage every flush renders uncached.
   */
  HS_COLD_MEMBER void init_storage(Arena &arena) {
    cached_warp_x = arena.allocate_n<int16_t>(CACHE_CELLS);
    cached_warp_y = arena.allocate_n<int16_t>(CACHE_CELLS);
    warp_cache_valid = false;
  }

  /**
   * @brief Accesses the bound Style.
   * @return Mutable reference to the bound feedback Style.
   */
  ::Feedback::Style &style() { return *feedback_style; }
  /**
   * @brief Accesses the bound Style.
   * @return Const reference to the bound feedback Style.
   */
  const ::Feedback::Style &style() const { return *feedback_style; }

private:
  static constexpr float WARP_SCALE = 128.0f;

  /**
   * @brief Tests whether the previous frame has any non-black pixel.
   * @param cv Canvas whose previous-frame buffer is scanned.
   * @return True on the first lit pixel found, false if the frame is all black.
   * @details Scans only this segment's clip band so another board's lit pixels
   * do not gate this board's flush.
   */
  static bool any_pixel_lit(const Canvas &cv) {
    const auto &clip = cv.clip();
    const auto x_clip = clip.x_clip();
    const ::Pixel *previous = cv.prev_data();
    for (int y = clip.render_y_start(); y < clip.render_y_end(); ++y) {
      const int row = y * W;
      for (int x = 0; x < W; ++x) {
        if (x_clip.active && x_clip.clipped(x))
          continue;
        const ::Pixel pixel = previous[row + x];
        if (pixel.r | pixel.g | pixel.b)
          return true;
      }
    }
    return false;
  }

  /**
   * @brief Chooses one color for a longitude-aliased pole row.
   * @param pole_row Base of the pole row in the previous frame.
   */
  HS_O3_FN static ::Pixel select_pole_sample(const ::Pixel *pole_row) {
    ::Pixel selected = pole_row[0];
    uint32_t selected_energy =
        static_cast<uint32_t>(selected.r) + selected.g + selected.b;
    for (int x = 1; x < W; ++x) {
      const ::Pixel candidate = pole_row[x];
      const uint32_t energy =
          static_cast<uint32_t>(candidate.r) + candidate.g + candidate.b;
      if (energy > selected_energy) {
        selected = candidate;
        selected_energy = energy;
      }
    }
    return selected;
  }

  /**
   * @brief Bilinearly samples the Canvas front buffer (previous frame).
   * @param field Spherical topology and interpolation policy.
   * @param prev Base of the previous-frame buffer, row-major with stride W.
   * @param poles Shared values for every aliased column of the pole rows.
   * @param bx Fractional column in [-W, 2W).
   * @param by Fractional row; north crossings reflect with a half-turn.
   * @param r Out: interpolated red on the [0, 65535] scale, unquantized.
   * @param g Out: interpolated green.
   * @param b Out: interpolated blue.
   */
  HS_O3_FN
  void sample_bilinear_prev(const SphereField &field, const ::Pixel *prev,
                            const ::Pixel *poles, float bx, float by, float &r,
                            float &g, float &b) const {
    field.sample_bilinear_rgb(prev, poles, bx, by, r, g, b);
  }

  /** @brief Quantizes an unclamped [0, 65535]-scale channel to a Pixel
   *  component, round-to-nearest; NaN maps to the hi bound. */
  static uint16_t quantize16(float v) {
    // clamp guards the cast against overshoot and maps NaN to the hi bound.
    return static_cast<uint16_t>(hs::clamp(v, 0.0f, 65535.0f) + 0.5f);
  }

  /** @brief Coarse grid downsample the warp cache is sized for (the default
   *  Style's; every preset keeps it). Other values render uncached. */
  static constexpr int CACHE_DOWNSAMPLE = ::Feedback::Style{}.downsample;
  static constexpr int CACHE_SOUTH_INFILL =
      CACHE_DOWNSAMPLE > hs::H_OFFSET ? CACHE_DOWNSAMPLE - hs::H_OFFSET : 0;
  static constexpr int CACHE_COLUMNS = W / CACHE_DOWNSAMPLE;
  /** @brief Cell count of the cached spherical warp field. */
  static constexpr int CACHE_CELLS =
      CACHE_COLUMNS * SphereField(CACHE_DOWNSAMPLE, CACHE_DOWNSAMPLE,
                                  CACHE_SOUTH_INFILL, CACHE_COLUMNS)
                          .ring_count();

public:
  /** @brief Persistent bytes init_storage() reserves (two int16 warp fields). */
  static constexpr size_t STORAGE_BYTES = 2 * CACHE_CELLS * sizeof(int16_t);

private:
  /** @brief Inputs the coarse warp field is a pure function of (stock
   *  transforms only); equal keys make the cached field reusable. */
  struct WarpKey {
    ::Feedback::SpaceFn space_fn;
    const Animation::NoiseParams *noise;
    float amplitude;
    float frequency;
    float speed;
    float scale;
    float time;
    int field_y_begin;
    int field_y_end;
    bool operator==(const WarpKey &) const = default;
  };

  ::Feedback::Style *feedback_style; /**< Bound feedback Style (non-owning). */
  bool enabled = true;               /**< When false, flush() is skipped. */
  WarpKey cached_warp_key{};         /**< Key for the cached warp field. */
  bool warp_cache_valid = false;    /**< True when the cached field is valid. */
  int16_t *cached_warp_x = nullptr; /**< Arena-owned cached column deltas. */
  int16_t *cached_warp_y = nullptr; /**< Arena-owned cached row deltas. */
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
  static constexpr int domain_rank = 2;
  /**
   * @brief The +1/+2/+3 column taps land outside the plotted position, so a
   *        segment worker needs 3 columns of render margin to write them.
   */
  static constexpr int segment_margin = 3;
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
