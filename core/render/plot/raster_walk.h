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
