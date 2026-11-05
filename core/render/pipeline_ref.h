/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/**
 * @file pipeline_ref.h
 * @brief Borrowed pipeline render handle and plot-provider dispatch.
 */

#include "core/engine/concepts.h"

/**
 * @brief Whether P writes the framebuffer through a cached base.
 * @tparam P Pipeline type; one without the member is not a direct-raster path.
 * @return `P::direct_raster_path`, or false when P lacks it.
 */
template <typename P> consteval bool pipeline_direct_raster_path() {
  if constexpr (requires { P::direct_raster_path; })
    return P::direct_raster_path;
  else
    return false;
}

/**
 * @brief Whether P exposes the three plot() overloads PipelineRef erases.
 * @tparam P Candidate pipeline type.
 * @details Spelled with the same argument types the erasing thunks pass, so it
 * accepts exactly what they compile against.
 */
template <typename P>
concept Plottable =
    requires(P &p, Canvas &cv, const math::Vector &v, const Pixel &c) {
      p.plot(cv, 0.0f, 0.0f, c, 0.0f, 0.0f);
      p.plot(cv, 0, 0, c, 0.0f, 0.0f);
      p.plot(cv, v, c, 0.0f, 0.0f);
    };

/**
 * @brief Dispatches a clip query through world stages or a bare plot provider.
 * @param pipeline Pipeline queried.
 * @param a Edge start (unit sphere point, pre-world-transform).
 * @param b Edge end (unit sphere point, pre-world-transform).
 * @param pb Optional planar basis for the edge (null = geodesic).
 * @param pred Rasterizer's screen-row-span vs clip-band test.
 * @return Whether the edge could intersect the clip band.
 */
template <typename PipelineT, typename Pred>
inline bool pipeline_could_intersect_clip(PipelineT &pipeline,
                                          const math::Vector &a,
                                          const math::Vector &b,
                                          const math::Basis *pb, Pred &&pred) {
  if constexpr (requires {
                  pipeline.could_intersect_clip(a, b, pb,
                                                std::forward<Pred>(pred));
                }) {
    return pipeline.could_intersect_clip(a, b, pb, std::forward<Pred>(pred));
  } else {
    static_assert(
        !requires { PipelineT::any_crosses_segments; },
        "pipeline exposes any_crosses_segments but not "
        "could_intersect_clip (signature drift)");
    return pred(a, b, pb);
  }
}

/** @brief Visits the poses used to space a still raster's samples.
 * @param pipeline Rendering pipeline.
 * @param a,b Source position and tangent.
 * @param pb Optional source basis.
 * @param pred Sampling-pose visitor.
 * @return The visitor's combined verdict. */
template <typename PipelineT, typename Pred>
inline bool pipeline_sampling_edge(PipelineT &pipeline, const math::Vector &a,
                                   const math::Vector &b, const math::Basis *pb,
                                   Pred &&pred) {
  if constexpr (requires { pipeline.sampling_edge(a, b, pb, pred); })
    return pipeline.sampling_edge(a, b, pb, pred);
  else
    return pipeline_could_intersect_clip(pipeline, a, b, pb, pred);
}

/**
 * @brief Non-owning, type-erased handle to a rasterizer pipeline.
 * @details Forwards 2D screen-space and 3D world-space plot() calls to the
 * wrapped object; borrows the target and must not outlive it. A direct-raster
 * sink converts only through the Canvas-taking constructor, which checks
 * prepared_for().
 */
class PipelineRef {
  /** @brief Tag selecting the type-erasing constructor. */
  template <bool Timed> struct Erase {};

  void *ctx;
  void (*plot2d)(void *, Canvas &, float, float, const Pixel &, float, float);
  void (*plot2d_int)(void *, Canvas &, int, int, const Pixel &, float, float);
  void (*plot3d)(void *, Canvas &, const math::Vector &, const Pixel &, float,
                 float);
  bool (*cull)(void *, const math::Vector &, const math::Vector &,
               const math::Basis *, CullEdgePredRef);
  bool (*sampling)(void *, const math::Vector &, const math::Vector &,
                   const math::Basis *,
                   CullEdgePredRef); ///< Still sampling poses.
  float (*motion_bound)(void *);     ///< Angular world-stage motion bound.
  /** @brief Shared entries used only by geometry-callback draws. */
  struct TimedEntries {
    /** @brief Evaluates world stages without projection. */
    math::Vector (*world)(void *, Canvas &, const math::Vector &, float);
    /** @brief Projects an explicitly timed point. */
    math::PixelCoords (*probe)(void *, Canvas &, const math::Vector &, float);
    /** @brief Plots an explicitly timed point. */
    float (*plot)(void *, Canvas &, const math::Vector &, const Pixel &, float,
                  float, float);
    /** @brief Plots a sample after the world stages. */
    float (*transformed_plot)(void *, Canvas &, const math::Vector &,
                              const Pixel &, float, float, float);
    /** @brief Plots an already-projected timed sample. */
    float (*projected_plot)(void *, Canvas &, math::PixelCoords, const Pixel &,
                            float, float, float);
    /** @brief Plots a projected sample's anti-alias taps in one column; null
     *  unless the first screen stage anti-aliases. */
    void (*column_plot)(void *, Canvas &, math::PixelCoords, int, const Pixel &,
                        float, float, float);
    /** @brief Records an unprojected timed sample in source history. */
    void (*record)(void *, Canvas &, const math::Vector &, float, float);
    bool rigid;             ///< World stages can run at controls.
    bool dense_screen_safe; ///< Screen stages retain a one-pixel footprint.
    /** @brief Returns the spatial copy count. */
    int (*copies)(void *);
    /** @brief Selects one spatial copy. */
    void (*select)(void *, int);
  };
  /** @brief Shared timed entries, absent on ordinary refs. */
  const TimedEntries *timed_entries = nullptr;

  /** @brief Binds ordinary plot and cull entries. */
  template <typename T>
  PipelineRef(T &t, Erase<false>)
      : ctx(std::addressof(t)), world_transform_is_identity([]() HS_HOT_INLINE {
          if constexpr (requires { T::world_transform_is_identity; })
            return T::world_transform_is_identity;
          else
            return true;
        }()),
        has_temporal_motion([]() HS_HOT_INLINE {
          if constexpr (requires { T::has_temporal_motion; })
            return T::has_temporal_motion;
          else
            return true;
        }()) {
    plot2d = [](void *pipeline, Canvas &cv, float x, float y, const Pixel &c,
                float age, float alpha) HS_HOT_INLINE {
      static_cast<T *>(pipeline)->plot(cv, x, y, c, age, alpha);
    };
    plot2d_int = [](void *pipeline, Canvas &cv, int x, int y, const Pixel &c,
                    float age, float alpha) HS_HOT_INLINE {
      static_cast<T *>(pipeline)->plot(cv, x, y, c, age, alpha);
    };
    plot3d = [](void *pipeline, Canvas &cv, const math::Vector &v,
                const Pixel &c, float age, float alpha) HS_HOT_INLINE {
      static_cast<T *>(pipeline)->plot(cv, v, c, age, alpha);
    };
    cull = still_cull_entry<T>;
    sampling = sampling_entry<T>;
    motion_bound = [](void *pipeline) HS_HOT_INLINE {
      if constexpr (requires(T &p) { p.sweep_arc(); })
        return static_cast<T *>(pipeline)->sweep_arc();
      else
        return 0.0f;
    };
  }

  /** @brief Applies an erased pipeline's ordinary edge cull.
   * @tparam T Concrete pipeline type.
   * @param pipeline Borrowed pipeline.
   * @param a,b Edge endpoints.
   * @param pb Optional planar basis.
   * @param pred Downstream visibility predicate.
   * @return Whether the edge may intersect the clip. */
  template <typename T>
  HS_O3_FN HS_NOINLINE_NOCLONE static bool
  still_cull_entry(void *pipeline, const math::Vector &a, const math::Vector &b,
                   const math::Basis *pb, CullEdgePredRef pred) {
    return pipeline_could_intersect_clip(*static_cast<T *>(pipeline), a, b, pb,
                                         pred);
  }

  /** @brief Visits still-raster sampling poses without testing visibility.
   * @tparam T Concrete pipeline type.
   * @param pipeline Borrowed pipeline.
   * @param a,b Source position and tangent.
   * @param pb Optional source basis.
   * @param pred Sampling-pose visitor.
   * @return The visitor's combined verdict. */
  template <typename T>
  HS_HOT_FLASH_MEMBER static bool
  sampling_entry(void *pipeline, const math::Vector &a, const math::Vector &b,
                 const math::Basis *pb, CullEdgePredRef pred) {
    if constexpr (requires(T &p) { p.sampling_edge(a, b, pb, pred); })
      return static_cast<T *>(pipeline)->sampling_edge(a, b, pb, pred);
    else
      return still_cull_entry<T>(pipeline, a, b, pb, pred);
  }

  /** @brief Evaluates an erased pipeline's world stages.
   * @tparam T Concrete pipeline type.
   * @param pipeline Borrowed pipeline.
   * @param cv Target canvas.
   * @param v Source point.
   * @param time Frame fraction.
   * @return World-space point. */
  template <typename T>
  HS_O3_FN HS_NOINLINE_NOCLONE static math::Vector
  timed_world_entry(void *pipeline, Canvas &cv, const math::Vector &v,
                    float time) {
    if constexpr (requires(T &p) { p.world_at(cv, v, time); })
      return static_cast<T *>(pipeline)->world_at(cv, v, time);
    else
      return v;
  }

  /** @brief Projects an erased pipeline without writing.
   * @tparam T Concrete pipeline type.
   * @param pipeline Borrowed pipeline.
   * @param cv Target canvas.
   * @param v Source point.
   * @param time Frame fraction.
   * @return Projected column and row. */
  template <typename T>
  HS_O3_FN HS_NOINLINE_NOCLONE static math::PixelCoords
  timed_probe_entry(void *pipeline, Canvas &cv, const math::Vector &v,
                    float time) {
    return static_cast<T *>(pipeline)->probe(cv, v, time);
  }

  /** @brief Emits an erased timed point through its stages.
   * @tparam T Concrete pipeline type.
   * @param pipeline Borrowed pipeline.
   * @param cv Target canvas.
   * @param v Source point.
   * @param c Source color.
   * @param age Age register.
   * @param alpha Opacity.
   * @param time Frame fraction.
   * @return Projected column. */
  template <typename T>
  HS_O3_FN HS_NOINLINE_NOCLONE static float
  timed_plot_entry(void *pipeline, Canvas &cv, const math::Vector &v,
                   const Pixel &c, float age, float alpha, float time) {
    return static_cast<T *>(pipeline)->plot_at(cv, v, c, age, alpha, time);
  }

  /** @brief Dispatches an already-transformed sample.
   * @tparam T Concrete pipeline.
   * @param pipeline Borrowed pipeline.
   * @param cv Frame context.
   * @param v Transformed point.
   * @param c,age,alpha Shading registers.
   * @param t Frame fraction.
   * @return Projected column. */
  template <typename T>
  HS_NOINLINE_NOCLONE static float
  transformed_plot_entry(void *pipeline, Canvas &cv, const math::Vector &v,
                         const Pixel &c, float age, float alpha, float t) {
    if constexpr (requires(T &p) {
                    p.plot_transformed_at(cv, v, c, age, alpha, t);
                  })
      return static_cast<T *>(pipeline)->plot_transformed_at(cv, v, c, age,
                                                             alpha, t);
    else
      return timed_plot_entry<T>(pipeline, cv, v, c, age, alpha, t);
  }

  /** @brief Dispatches an already-projected timed sample.
   * @tparam T Concrete pipeline.
   * @param pipeline Borrowed pipeline.
   * @param cv Frame context.
   * @param p Projected column and row.
   * @param c,age,alpha Shading registers.
   * @param t Frame fraction.
   * @return Projected column. */
  template <typename T>
  HS_NOINLINE_NOCLONE static float
  projected_plot_entry(void *pipeline, Canvas &cv, math::PixelCoords p,
                       const Pixel &c, float age, float alpha, float t) {
    if constexpr (requires(T &v) {
                    v.plot_projected_at(cv, p, c, age, alpha, t);
                  })
      return static_cast<T *>(pipeline)->plot_projected_at(cv, p, c, age, alpha,
                                                           t);
    else {
      static_cast<T *>(pipeline)->plot(cv, p.x, p.y, c, age, alpha);
      return p.x;
    }
  }

  /** @brief Dispatches a projected sample's taps in one column.
   * @tparam T Concrete pipeline.
   * @param pipeline Borrowed pipeline.
   * @param cv Frame context.
   * @param p Projected column and row.
   * @param column Wrapped destination column.
   * @param c,age,alpha Shading registers.
   * @param t Frame fraction. */
  template <typename T>
  HS_HOT_FLASH_MEMBER static void
  column_plot_entry(void *pipeline, Canvas &cv, math::PixelCoords p, int column,
                    const Pixel &c, float age, float alpha, float t) {
    static_cast<T *>(pipeline)->plot_column_at(cv, p, column, c, age, alpha, t);
  }

  /** @brief Records a timed sample before its screen taps.
   * @tparam T Concrete pipeline.
   * @param pipeline Borrowed pipeline.
   * @param cv Frame context.
   * @param v Source point.
   * @param age Incoming sample age.
   * @param t Frame fraction. */
  template <typename T>
  HS_HOT_FLASH_MEMBER static void
  record_sample_entry(void *pipeline, Canvas &cv, const math::Vector &v,
                      float age, float t) {
    static_cast<T *>(pipeline)->record_world_sample(cv, v, age, t);
  }

  /** @brief Applies an erased pipeline's timed cull chain.
   * @tparam T Concrete pipeline type.
   * @param pipeline Borrowed pipeline.
   * @param a,b Source edge endpoints.
   * @param pb Optional planar basis.
   * @param pred Downstream visibility predicate.
   * @return Whether the edge can intersect the clip. */
  template <typename T>
  HS_HOT_FLASH_MEMBER static bool
  timed_cull_entry(void *pipeline, const math::Vector &a, const math::Vector &b,
                   const math::Basis *pb, CullEdgePredRef pred) {
    if constexpr (requires(T &p) {
                    p.could_intersect_timed_clip(a, b, pb, pred);
                  })
      return static_cast<T *>(pipeline)->could_intersect_timed_clip(a, b, pb,
                                                                    pred);
    return pipeline_could_intersect_clip(*static_cast<T *>(pipeline), a, b, pb,
                                         pred);
  }

  /** @brief Adds timed entries to the ordinary plot interface. */
  template <typename T>
  HS_HOT_FLASH_MEMBER PipelineRef(T &t, Erase<true>)
      : PipelineRef(t, Erase<false>{}) {
    cull = timed_cull_entry<T>;
    static constexpr TimedEntries entries{
        timed_world_entry<T>,
        timed_probe_entry<T>,
        timed_plot_entry<T>,
        transformed_plot_entry<T>,
        projected_plot_entry<T>,
        []() HS_HOT_INLINE {
          if constexpr (requires { requires T::column_taps_ready; })
            return column_plot_entry<T>;
          else
            return static_cast<decltype(&column_plot_entry<T>)>(nullptr);
        }(),
        []() HS_HOT_INLINE {
          if constexpr (requires { requires T::any_3d_history; })
            return record_sample_entry<T>;
          else
            return static_cast<decltype(&record_sample_entry<T>)>(nullptr);
        }(),
        []() HS_HOT_INLINE {
          if constexpr (requires { T::rigid_world_stages; })
            return T::rigid_world_stages;
          else
            return false;
        }(),
        []() HS_HOT_INLINE {
          if constexpr (requires {
                          T::total_segment_margin;
                          T::any_2d_history;
                          T::any_reads_outside_band;
                        })
            return T::total_segment_margin <= 1 && !T::any_2d_history &&
                   !T::any_reads_outside_band;
          else
            return false;
        }(),
        [](void *pipeline) HS_HOT_INLINE {
          if constexpr (requires(T &p) { p.copy_count(); })
            return static_cast<T *>(pipeline)->copy_count();
          else
            return 1;
        },
        [](void *pipeline, int copy) HS_HOT_INLINE {
          if constexpr (requires(T &p) { p.select_copy(copy); })
            static_cast<T *>(pipeline)->select_copy(copy);
        }};
    timed_entries = &entries;
  }

public:
  /** @brief Moving-stage angular bound of the pipeline.
   * @return Bound in radians. */
  float sweep_arc() const { return motion_bound(ctx); }
  /** @brief Number of independently walked spatial copies.
   * @return Copy count, or one for an untimed ref. */
  int copy_count() const {
    return timed_entries ? timed_entries->copies(ctx) : 1;
  }
  /** @brief Selects one fixed spatial copy.
   * @param copy Flattened copy index. */
  void select_copy(int copy) const {
    if (timed_entries)
      timed_entries->select(ctx, copy);
  }
  /** @brief Records one source sample before projected taps.
   * @param cv Frame context.
   * @param v Source position.
   * @param age Incoming sample age.
   * @param t Frame fraction.
   * @param column Owned tap column, or -1 for a complete splat.
   * @param x Projected sample column. */
  HS_HOT_INLINE void record_world_sample(Canvas &cv, const math::Vector &v,
                                         float age, float t, int column = -1,
                                         float x = 0.0f) const {
    if (timed_entries && timed_entries->record) {
      const int first = static_cast<int>(x);
      if (column < 0 ||
          column == (first < cv.width() ? first : first - cv.width()))
        timed_entries->record(ctx, cv, v, age, t);
    }
  }

  /** @brief Whether world transformations can run at controls.
   * @return True for a fully rigid, history-free world chain. */
  bool rigid_world_stages() const {
    return timed_entries && timed_entries->rigid;
  }

  /** @brief Whether screen stages permit a one-pixel dense-run cull pad.
   * @return True for bounded, history-free screen chains. */
  bool dense_screen_safe() const {
    return timed_entries && timed_entries->dense_screen_safe;
  }

  /** @brief Whether the first screen stage supports column-owned taps.
   * @return True when timed column plotting is available. */
  bool column_taps_ready() const {
    return timed_entries && timed_entries->column_plot;
  }

  /** @brief Plots after already-evaluated world transformations.
   * @param cv Frame context.
   * @param v Transformed point.
   * @param c,age,alpha Shading registers.
   * @param t Frame fraction.
   * @return Projected column. */
  HS_HOT_INLINE float plot_transformed_at(Canvas &cv, const math::Vector &v,
                                          const Pixel &c, float age,
                                          float alpha, float t) const {
    assert(timed_entries && timed_entries->transformed_plot);
    return timed_entries->transformed_plot(ctx, cv, v, c, age, alpha, t);
  }

  /** @brief Emits a projected sample through screen stages.
   * @param cv Frame context.
   * @param p Projected column and row.
   * @param c,age,alpha Shading registers.
   * @param t Frame fraction.
   * @return Projected column. */
  HS_HOT_INLINE float plot_projected_at(Canvas &cv, math::PixelCoords p,
                                        const Pixel &c, float age, float alpha,
                                        float t) const {
    assert(timed_entries && timed_entries->projected_plot);
    return timed_entries->projected_plot(ctx, cv, p, c, age, alpha, t);
  }

  /** @brief Emits a projected sample's anti-alias taps in one column.
   * @param cv Frame context.
   * @param p Projected column and row.
   * @param column Wrapped destination column.
   * @param c,age,alpha Shading registers.
   * @param t Frame fraction.
   * @details The erased pipeline's first screen stage anti-aliases. */
  HS_HOT_INLINE void plot_column_at(Canvas &cv, math::PixelCoords p, int column,
                                    const Pixel &c, float age, float alpha,
                                    float t) const {
    assert(timed_entries && timed_entries->column_plot);
    timed_entries->column_plot(ctx, cv, p, column, c, age, alpha, t);
  }

  /** @brief Whether the erased pipeline leaves world positions unchanged. */
  bool world_transform_is_identity;
  /** @brief Whether a stage evaluates world motion over a frame. */
  bool has_temporal_motion;

  /**
   * @brief Wraps any object exposing 2D and 3D plot() methods.
   * @tparam T Pipeline type; must satisfy Plottable and be neither PipelineRef
   *         nor a direct-raster sink.
   * @param t Pipeline object whose address is stored; must outlive this ref.
   */
  template <typename T>
    requires(Plottable<T> && !std::same_as<std::decay_t<T>, PipelineRef> &&
             !pipeline_direct_raster_path<std::decay_t<T>>())
  PipelineRef(T &t) : PipelineRef(t, Erase<false>{}) {}

  /** @brief Borrows a pipeline with explicit-time projection and plotting.
   * @param t Pipeline exposing probe and plot_at.
   * @details The pipeline must outlive this ref. */
  template <typename T>
    requires(!pipeline_direct_raster_path<std::decay_t<T>>())
  PipelineRef(T &t, Spacetime::Timed) : PipelineRef(t, Erase<true>{}) {}

  /**
   * @brief Erases a direct-raster sink for a draw into @p cv.
   * @param t Sink whose address is stored; must outlive this ref.
   * @param cv Canvas the draw writes into.
   * @details Checks prepared_for(@p cv): the sink writes through a framebuffer
   * base cached by prepare(), and a stale base is the buffer being scanned out.
   */
  template <typename T>
    requires(Plottable<T> && pipeline_direct_raster_path<std::decay_t<T>>())
  PipelineRef(T &t, Canvas &cv) : PipelineRef(t, Erase<false>{}) {
    HS_CHECK(t.prepared_for(cv),
             "direct raster pipeline not prepared for this canvas");
  }

  /** @brief Borrows a prepared direct sink with explicit-time entry points.
   * @param t Direct sink.
   * @param cv Prepared canvas. */
  template <typename T>
    requires(pipeline_direct_raster_path<std::decay_t<T>>())
  PipelineRef(T &t, Canvas &cv, Spacetime::Timed)
      : PipelineRef(t, Erase<true>{}) {
    HS_CHECK(t.prepared_for(cv),
             "timed direct raster pipeline not prepared for this canvas");
  }

  /** @brief Evaluates world stages without projection or history recording.
   * @param cv Frame context.
   * @param v Source point.
   * @param t Frame fraction.
   * @return Transformed world point. */
  math::Vector world_at(Canvas &cv, const math::Vector &v, float t) const {
    assert(timed_entries);
    return timed_entries->world(ctx, cv, v, t);
  }

  /** @brief Projects without plotting or recording.
   * @param cv Frame context.
   * @param v World point.
   * @param t Frame fraction.
   * @return Column and row through the world stages. */
  math::PixelCoords probe(Canvas &cv, const math::Vector &v, float t) const {
    assert(timed_entries);
    return timed_entries->probe(ctx, cv, v, t);
  }

  /** @brief Plots a world point at its supplied frame fraction.
   * @param cv Frame context.
   * @param v World point.
   * @param c Source color.
   * @param age Age register.
   * @param alpha Source opacity.
   * @param t Frame fraction.
   * @return Projected column. */
  float plot_at(Canvas &cv, const math::Vector &v, const Pixel &c, float age,
                float alpha, float t) const {
    assert(timed_entries);
    return timed_entries->plot(ctx, cv, v, c, age, alpha, t);
  }

  /**
   * @brief Plots a pixel at floating-point screen coordinates.
   * @param cv Target canvas to draw into.
   * @param x Column position in pixels (screen space).
   * @param y Row position in pixels (screen space).
   * @param c Source color to plot.
   * @param age Temporal age in frames.
   * @param alpha Coverage/opacity in [0, 1].
   */
  void plot(Canvas &cv, float x, float y, const Pixel &c, float age,
            float alpha) const {
    plot2d(ctx, cv, x, y, c, age, alpha);
  }
  /**
   * @brief Plots a pixel at integer screen coordinates.
   * @param cv Target canvas to draw into.
   * @param x Column position in pixels (screen space).
   * @param y Row position in pixels (screen space).
   * @param c Source color to plot.
   * @param age Temporal age in frames.
   * @param alpha Coverage/opacity in [0, 1].
   * @details Preserves integer coordinates when forwarding to the wrapped
   * pipeline.
   */
  void plot(Canvas &cv, int x, int y, const Pixel &c, float age,
            float alpha) const {
    plot2d_int(ctx, cv, x, y, c, age, alpha);
  }
  /**
   * @brief Plots a pixel at a 3D world-space position.
   * @param cv Target canvas to draw into.
   * @param v World-space position to project and plot.
   * @param c Source color to plot.
   * @param age Temporal age in frames.
   * @param alpha Coverage/opacity in [0, 1].
   */
  void plot(Canvas &cv, const math::Vector &v, const Pixel &c, float age,
            float alpha) const {
    plot3d(ctx, cv, v, c, age, alpha);
  }

  /**
   * @brief Clip-cull query forwarded to the wrapped pipeline.
   * @param a Edge start (unit sphere point, pre-world-transform).
   * @param b Edge end (unit sphere point, pre-world-transform).
   * @param pb Optional planar basis for the edge (null = geodesic).
   * @param pred Rasterizer's screen-row-span vs clip-band test.
   * @return Whether any world-transformed copy of the edge could intersect the
   *         clip band (see Pipeline::could_intersect_clip).
   */
  bool could_intersect_clip(const math::Vector &a, const math::Vector &b,
                            const math::Basis *pb, CullEdgePredRef pred) const {
    return cull(ctx, a, b, pb, pred);
  }

  /** @brief Visits the poses used to space a still raster's samples.
   * @param a,b Source position and tangent.
   * @param pb Optional source basis.
   * @param pred Sampling-pose visitor.
   * @return The visitor's combined verdict. */
  bool sampling_edge(const math::Vector &a, const math::Vector &b,
                     const math::Basis *pb, CullEdgePredRef pred) const {
    return sampling(ctx, a, b, pb, pred);
  }
};
