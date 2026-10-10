/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by core/render/plot/raster.h.

/**
 * @brief Starts a raster sample at @p src's position.
 * @tparam INTERP Copy @p src's registers; otherwise every register takes its
 *         Fragment default, including age = 0, and @p src's registers are not
 *         read.
 * @param f Output sample, positioned at src.pos with transparent black color.
 * @param src Source control point.
 */
template <bool INTERP>
__attribute__((always_inline)) inline void seed_fragment(Fragment &f,
                                                         const Fragment &src) {
  if constexpr (INTERP)
    f = src;
  else
    f = Fragment{};
  f.pos = src.pos;
  f.color = Color4(0, 0, 0, 0);
}

/**
 * @brief Starts a raster sample between two control points.
 * @tparam INTERP Interpolate the endpoints' registers at @p t; otherwise every
 *         register takes its Fragment default, including age = 0, and no
 *         endpoint register is read.
 * @param f Output sample, positioned at @p pos with transparent black color.
 * @param a Segment start.
 * @param b Segment end.
 * @param t Register interpolation factor.
 * @param pos Sample position.
 */
template <bool INTERP>
__attribute__((always_inline)) inline void
lerp_fragment(Fragment &f, const Fragment &a, const Fragment &b, float t,
              const math::Vector &pos) {
  if constexpr (INTERP)
    f = Fragment::lerp_registers(a, b, t);
  else
    f = Fragment{};
  f.pos = pos;
  f.color = Color4(0, 0, 0, 0);
}

/**
 * @brief Runs the fragment shader on one raster sample.
 * @param fragment_shader Fragment shader.
 * @param pos Sample position passed to the shader.
 * @param f Sample fragment, shaded in place.
 */
template <typename FragmentShaderT>
__attribute__((always_inline)) inline void
shade_fragment(FragmentShaderT &fragment_shader, const math::Vector &pos,
               Fragment &f) {
  HS_PLOT_COUNT(shader_calls);
  HS_PLOT_STALL_START(shade_start);
  HS_PLOT_RENDER_COUNT(fragment_shader_calls);
  fragment_shader(pos, f);
  HS_PLOT_STALL_STOP(shade_palette, shade_start);
}

/**
 * @brief Shades one raster sample and plots it at its world position.
 * @param pipeline Pipeline that plots the sample.
 * @param canvas Target canvas.
 * @param fragment_shader Fragment shader.
 * @param pos Sample position, shaded and plotted.
 * @param f Sample fragment, shaded in place.
 * @param alpha_scale Balanced-sampling step ratio applied to the shaded alpha
 *        through balanced_sample_alpha; none plots the shaded alpha.
 */
template <typename PipelineT, typename FragmentShaderT>
__attribute__((always_inline)) inline void
shade_and_plot(PipelineT &pipeline, Canvas &canvas,
               FragmentShaderT &fragment_shader, const math::Vector &pos,
               Fragment &f, std::optional<float> alpha_scale = std::nullopt) {
  shade_fragment(fragment_shader, pos, f);
  if (alpha_scale)
    f.color.alpha = balanced_sample_alpha(f.color.alpha, *alpha_scale);
  HS_PLOT_COUNT(plotted_samples);
  pipeline.plot(canvas, pos, f.color.color, f.age, f.color.alpha);
}

/**
 * @brief Rendered-arc v0/v1 registers for the in-flight segment of a
 *        planar-basis polyline.
 * @tparam DERIVE Derive the registers; false compiles the stamp out.
 * @details Issued by PlanarArcTable::bind_and_measure and positioned per
 * segment by PlanarArcTable::advance.
 */
template <bool DERIVE> class PlanarArcStamp {
public:
  /**
   * @brief Rewrites @p f's v0/v1 from the rendered arc; no-op when inactive.
   * @param f Sample fragment.
   * @param d Arc drawn so far within the segment.
   */
  __attribute__((always_inline)) void stamp(Fragment &f, float d) const {
    if constexpr (!DERIVE)
      return;
    if (!active)
      return;
    float arc = seg_base + d;
    f.v1 = arc;
    if (total_arc > math::EPS_GEOMETRIC)
      f.v0 = arc / total_arc;
  }

private:
  template <bool> friend class PlanarArcTable;
  /** Rendered arc at the segment's start. */
  float seg_base = 0.0f;
  /** Rendered arc of the whole polyline. */
  float total_arc = 0.0f;
  /** False for a geodesic polyline, which keeps its source v0/v1. */
  bool active = false;
};

/**
 * @brief Per-segment rendered arc lengths and antipode-seam flags of a
 *        planar-basis polyline.
 * @tparam DERIVE Derive v0/v1 from the rendered arc; false compiles the table
 *         out.
 * @details Under a planar basis the rendered edge bows longer than the
 * geodesic chord, so v1 is the rendered arc reached and v0 that arc over the
 * polyline's total. The table owns the caches and the arc reached so far; it
 * stays inactive for a geodesic polyline.
 */
template <bool DERIVE> class PlanarArcTable {
public:
  /**
   * @brief Binds the caches in scratch_arena_a and measures every segment.
   * @param points Polyline control points.
   * @param segment_next Returns the end point of segment i.
   * @param count Rasterized segment count.
   * @param basis Planar basis, or null for a geodesic polyline.
   * @return The v0/v1 stamp, positioned by advance() before each segment.
   */
  template <typename SegmentNextT>
  __attribute__((always_inline)) PlanarArcStamp<DERIVE>
  bind_and_measure(const Fragments &points, SegmentNextT &&segment_next,
                   size_t count, const math::Basis *basis) {
    PlanarArcStamp<DERIVE> out;
    if constexpr (DERIVE) {
      planar_basis = basis;
      if (basis == nullptr)
        return out;
      arc_cache.bind(scratch_arena_a, count);
      seam_cache.bind(scratch_arena_a, count);
      const math::Vector &pcenter = basis->v;
      for (size_t i = 0; i < count; i++) {
        const math::Vector &a = points[i].pos;
        const math::Vector &b = segment_next(i).pos;
        const bool seam = math::dot(a, pcenter) < -COS_PLANAR_ANTIPODE ||
                          math::dot(b, pcenter) < -COS_PLANAR_ANTIPODE;
        seam_cache.push_back(seam ? 1 : 0);
        float seg =
            seam ? unit_arc_length(a, b) : planar_arc_length(a, b, *basis);
        arc_cache.push_back(seg);
        total_arc += seg;
      }
      out.total_arc = total_arc;
      out.active = true;
    }
    return out;
  }

  /** @brief True when the caches are bound for a planar polyline. */
  __attribute__((always_inline)) bool is_active() const {
    if constexpr (DERIVE)
      return planar_basis != nullptr;
    else
      return false;
  }

  /** @brief Segment @p i has an endpoint at the basis antipode. */
  __attribute__((always_inline)) bool seam(size_t i) const {
    return seam_cache[i] != 0;
  }

  /**
   * @brief Moves @p stamp's arc origin to the start of segment @p i.
   * @details Called for every segment in order, drawn or culled, so v0/v1
   * span the full curve.
   * @param i Segment index.
   * @param stamp Stamp returned by bind_and_measure.
   */
  __attribute__((always_inline)) void advance(size_t i,
                                              PlanarArcStamp<DERIVE> &stamp) {
    if (!is_active())
      return;
    stamp.seg_base = cumul;
    cumul += arc_cache[i];
  }

private:
  ArenaVector<float> arc_cache;
  ArenaVector<uint8_t> seam_cache;
  float total_arc = 0.0f;
  float cumul = 0.0f;
  const math::Basis *planar_basis = nullptr;
};
