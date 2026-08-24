/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <utility>
#include <type_traits>
#include <cmath>
#include <cassert>

#include <span>
#include <algorithm>
#include <bitset>
#include "math/geometry.h"
#include "math/stereographic.h"
#include "color/color.h"
#include "render/canvas.h"
#include "engine/concepts.h"
#include "engine/memory.h"

/**
 * @file pipeline.h
 * @brief Pipeline composition: the filter-domain traits every stage derives
 * from, the pipeline-level trait surface, and the recursive Pipeline template
 * that folds a stage list down onto the canvas sink.
 */

/** @brief Callback that forwards a 2D plot (x, y, pixel, age, alpha) downstream. */
using PassFn2D = FunctionRef<void(float, float, const ::Pixel &, float, float)>;
/** @brief Callback that forwards a 3D plot (vector, pixel, age, alpha) downstream. */
using PassFn3D =
    FunctionRef<void(const Vector &, const ::Pixel &, float, float)>;

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
 * a unit-assuming one. `emits_pixel_centers` / `requires_subpixel_input`: a
 * screen stage that rounds its taps to pixel centers must not precede one whose
 * whole job is the sub-pixel fraction, which the rounding would make an exact
 * identity. `crosses_segments`: output can move between worker
 * segment bands, so the effect must render the full canvas. `reads_outside_band`:
 * the stage samples framebuffer pixels outside the display band, so stale pixels
 * outside that band must also be cleared. Both default to `has_history`
 * (fail-safe). `segment_margin`: how many pixels the stage's output can land
 * away from the plotted position, i.e. how far the render bounds must be padded
 * past the display band for a segment worker to still write the stage's taps.
 * `is_pipeline`: the type is a whole pipeline rather than a stage, so it may not
 * be listed inside a `Pipeline<>`.
 * `world_transform_is_identity`: the stage forwards world points unmoved (it may
 * still mask or recolor), so the clip cull may run against the source geometry.
 * A world stage that moves points must instead define `cull_edge` or set
 * `crosses_segments`.
 * A stage overrides any of these by redeclaring it.
 */
template <bool Is2d, bool HasHistory> struct FilterTraits {
  static constexpr int domain_rank = Is2d ? 1 : 0;
  static constexpr bool is_2d = Is2d;
  static constexpr bool has_history = HasHistory;
  static constexpr bool is_pipeline = false;
  static constexpr bool is_terminal = false;
  static constexpr bool terminal_replaces = false;
  static constexpr bool emits_nonunit_world = false;
  static constexpr bool requires_unit_world_input = false;
  static constexpr bool emits_pixel_centers = false;
  static constexpr bool requires_subpixel_input = false;
  static constexpr bool crosses_segments = has_history;
  static constexpr bool reads_outside_band = has_history;
  static constexpr bool world_transform_is_identity = is_2d;
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
 * @brief Trait base for a screen-space sink that stands in for a whole Pipeline
 * instead of running as a stage inside one.
 * @details Keeps the stage vocabulary the Pipeline folds read, so a sink handed
 * to `Pipeline<>` still reaches the ordering asserts and is rejected there by
 * `is_pipeline` rather than deep inside plot() overload resolution.
 */
struct IsPipelineSink : FilterTraits<true, false> {
  static constexpr bool is_pipeline = true;
};

/**
 * @brief The pipeline-level trait surface a whole pipeline answers, as opposed
 *        to the stage vocabulary in FilterTraits.
 * @details The single authority for the list. A direct sink stands in for a
 * Pipeline<> and hand-mirrors every member, so a new fold member belongs here
 * too: the readers in plot_cull.h are `requires`-guarded and silently fall back
 * to their defaults on a type that never grew it, which reaches sink-based
 * effects as a wrong EffectConfig rather than a compile error.
 */
template <typename T>
concept PipelineFoldSurface = requires {
  requires std::is_same_v<decltype(T::is_pipeline), const bool>;
  requires std::is_same_v<decltype(T::any_crosses_segments), const bool>;
  requires std::is_same_v<decltype(T::any_reads_outside_band), const bool>;
  requires std::is_same_v<decltype(T::any_2d_history), const bool>;
  requires std::is_same_v<decltype(T::any_3d_history), const bool>;
  requires std::is_same_v<decltype(T::any_2d_trail_history), const bool>;
  requires std::is_same_v<decltype(T::any_terminal_history), const bool>;
  requires std::is_same_v<decltype(T::has_world_cull), const bool>;
  requires std::is_same_v<decltype(T::has_world_stage), const bool>;
  requires std::is_same_v<decltype(T::direct_raster_path), const bool>;
  requires std::is_same_v<decltype(T::segment_margin), const int>;
  requires std::is_same_v<decltype(T::total_segment_margin), const int>;
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
  static constexpr bool is_pipeline = true;
  static constexpr bool direct_raster_path = false;
  static constexpr bool is_terminal = false;
  // Stage vocabulary, so a Pipeline nested inside a Pipeline<> reaches the
  // is_pipeline diagnostic instead of failing in the trait folds first.
  static constexpr bool has_history = false;
  static constexpr bool terminal_replaces = false;
  static constexpr bool emits_nonunit_world = false;
  static constexpr bool requires_unit_world_input = false;
  static constexpr bool emits_pixel_centers = false;
  static constexpr bool requires_subpixel_input = false;
  static constexpr bool crosses_segments = false;
  static constexpr bool reads_outside_band = false;
  static constexpr bool world_transform_is_identity = true;
  static constexpr int segment_margin = 0;
  static constexpr bool any_crosses_segments = false;
  static constexpr bool any_reads_outside_band = false;
  static constexpr int total_segment_margin = 0;
  static constexpr bool any_2d_history = false;
  static constexpr bool any_3d_history = false;
  static constexpr bool any_2d_trail_history = false;
  static constexpr bool any_terminal_history = false;
  /** @brief No stage re-emits clip-cull edges (see the recursive case). */
  static constexpr bool has_world_cull = false;
  /** @brief No stage runs in world space (see the recursive case). */
  static constexpr bool has_world_stage = false;
  /** @brief Occurrences of stage type T in this pipeline (base case: none). */
  template <typename T> static constexpr int stage_count = 0;

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
  void plot(Canvas &cv, int x, int y, const ::Pixel &c, float, float alpha) {
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
  void plot_in_bounds(Canvas &cv, int x, int y, const ::Pixel &c, float,
                      float alpha) {
    HS_PROFILE(filter_blend);
    assert(x >= 0 && x < W && cv.clip().contains_x(x));
    assert(cv.clip().contains_y(y));
    ::Pixel &dst = cv(x, y);
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
  void plot(Canvas &cv, float x, float y, const ::Pixel &c, float,
            float alpha) {
    // Non-finite coords make the int casts below UB and bypass the wrap.
    assert(std::isfinite(x) && std::isfinite(y));
    // y never wraps; bounded only so the cast below stays in range.
    assert(y >= -H && y < 2 * H);
    const float xr = std::round(x);
    // fast_wrap corrects only a single ±W offset, so xr must land in [-W, 2W).
    assert(xr >= -W && xr < 2 * W);
    int xi = static_cast<int>(xr);
    int yi = static_cast<int>(std::round(y));
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
  void plot(Canvas &cv, const Vector &v, const ::Pixel &c, float age,
            float alpha) {
    auto p = vector_to_pixel<W, H>(v);
    plot(cv, p.x, p.y, c, age, alpha);
  }

  /**
   * @brief Trail flush (base case: no stage to flush).
   * @tparam TrailFn ScreenTrailFn or WorldTrailFn.
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
  static constexpr bool is_pipeline = true;
  static constexpr bool direct_raster_path = false;
  static constexpr bool is_terminal = Head::is_terminal || Next::is_terminal;
  static constexpr bool any_crosses_segments =
      Head::crosses_segments || Next::any_crosses_segments;
  static constexpr bool any_reads_outside_band =
      Head::reads_outside_band || Next::any_reads_outside_band;
  // Sum, not max: each stage displaces the taps the stages before it already
  // displaced, so chained spreading stages compose additively.
  static constexpr int total_segment_margin =
      Head::segment_margin + Next::total_segment_margin;

  static constexpr bool crosses_segments = any_crosses_segments;
  static constexpr bool reads_outside_band = any_reads_outside_band;
  static constexpr int segment_margin = total_segment_margin;

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

  // Stage vocabulary, so a Pipeline nested inside a Pipeline<> reaches the
  // is_pipeline diagnostic instead of failing in the trait folds first.
  static constexpr bool has_history = any_2d_history || any_3d_history;
  static constexpr bool terminal_replaces =
      Head::terminal_replaces || Next::terminal_replaces;
  static constexpr bool emits_nonunit_world =
      Head::emits_nonunit_world || Next::emits_nonunit_world;
  static constexpr bool requires_unit_world_input =
      Head::requires_unit_world_input || Next::requires_unit_world_input;
  static constexpr bool emits_pixel_centers =
      Head::emits_pixel_centers || Next::emits_pixel_centers;
  static constexpr bool requires_subpixel_input =
      Head::requires_subpixel_input || Next::requires_subpixel_input;
  static constexpr bool world_transform_is_identity =
      Head::world_transform_is_identity && Next::world_transform_is_identity;

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

  /** @brief Occurrences of stage type T across Head and the Tail. */
  template <typename T>
  static constexpr int stage_count =
      (std::is_same_v<Head, T> ? 1 : 0) + Next::template stage_count<T>;

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

  /** @brief Default-constructs every stage in the pipeline. */
  Pipeline() = default;

  /**
   * @brief Type-safe filter accessor: finds the stage of type T.
   * @tparam T Filter type to retrieve.
   * @return Reference to the stage of type T (recurses into the Tail if Head is
   * not T).
   */
  template <typename T> T &get() {
    static_assert(
        stage_count<T> <= 1,
        "Ambiguous get<T>(): this stage type is listed more than once, so "
        "get<T>() resolves to the first occurrence and the later ones are "
        "unreachable. Give each slot a distinct type (a trivial subclass "
        "inheriting the stage's constructors).");
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
    static_assert(
        stage_count<T> <= 1,
        "Ambiguous get<T>(): this stage type is listed more than once, so "
        "get<T>() resolves to the first occurrence and the later ones are "
        "unreachable. Give each slot a distinct type (a trivial subclass "
        "inheriting the stage's constructors).");
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
  void plot(Canvas &cv, float x, float y, const ::Pixel &c, float age,
            float alpha) {
    if constexpr (Head::is_2d) {
      Head::plot(
          x, y, c, age, alpha,
          [&](float nx, float ny, const ::Pixel &nc, float nage, float nalpha) {
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
  void plot(Canvas &cv, int x, int y, const ::Pixel &c, float age,
            float alpha) {
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
  void plot(Canvas &cv, const Vector &v, const ::Pixel &c, float age,
            float alpha) {
    if constexpr (!Head::is_2d) {
      Head::plot(v, c, age, alpha,
                 [&](const Vector &nv, const ::Pixel &nc, float nage,
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
      !Head::is_pipeline,
      "Not a filter stage: this type is a whole pipeline (a nested Pipeline, or "
      "a direct sink such as Screen::DirectAntiAliasSink) — it writes the "
      "Canvas itself and takes no downstream `pass` callback. Use it on its "
      "own, or list the stage it replaces (Screen::AntiAlias).");

  static_assert(
      Head::world_transform_is_identity || has_cull_edge<Head> ||
          Head::crosses_segments,
      "A World stage that moves geometry must let the clip cull follow it, or "
      "the rasterizer culls edges by their SOURCE latitude and drops the ones "
      "the transform moves into a segment band. Override cull_edge() to "
      "re-emit the edge under the transform (World::Orient), or set "
      "crosses_segments = true so the effect renders the full canvas "
      "(World::Mobius). A stage that only masks or recolors declares "
      "world_transform_is_identity = true.");

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
      !Head::emits_pixel_centers || (... && !Tail::requires_subpixel_input),
      "Filter ordering: a screen stage that rounds its taps to pixel centers "
      "(Screen::Blur) must not precede one that reads the sub-pixel fraction "
      "(Screen::AntiAlias) — the rounding makes that stage an exact identity. "
      "Reorder so the sub-pixel stage runs first.");

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
    static_assert(
        any_2d_trail_history,
        "Discarded flush() callback: this Pipeline's only 2D history is a "
        "terminal stage (Pixel::Feedback), which composites into the Canvas "
        "itself and takes no trail callback. Pass flush(cv, alpha) instead.");
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
    static_assert(
        !any_3d_history,
        "Incomplete flush(): this Pipeline also carries a 3D history stage "
        "(World::Trails) that this overload leaves unflushed, so its ring "
        "buffer fills to capacity and never decays. Pass both callbacks: "
        "flush(cv, worldTrailFn, screenTrailFn, alpha).");
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
