/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_plot_scan.h.

// ---------------------------------------------------------------------------
// Headless rasterize() harness. A stub pipeline records plotted world positions;
// a no-op Effect supplies a Canvas with a full (unclipped) clip band.
// ---------------------------------------------------------------------------

/**
 * @brief Pipeline stub: records each plotted position; ignores color/age.
 * @details Carries both plot() overloads so it can also back a type-erased
 * PipelineRef (the Plot::ParticleSystem draw() entry points take one); only
 * the 3D world-space overload records — the 2D screen-space form is unused by
 * these paths.
 */
struct CapturePipeline {
  std::vector<math::Vector>
      plotted; /**< World positions handed to plot(), in order. */
  void plot(Canvas &, const math::Vector &v, const Pixel &, float, float) {
    plotted.push_back(v);
  }
  void plot(Canvas &, float, float, const Pixel &, float, float) {}
};

/** @brief Non-erased capture sink exposing the direct-path cull traits. */
struct DirectCapturePipeline : CapturePipeline {
  static constexpr bool has_world_cull = false;
  static constexpr bool direct_raster_path = true;
};

/** @brief Pipeline stub recording world positions and effective alpha. */
struct AlphaCapturePipeline {
  std::vector<math::Vector> plotted;
  std::vector<float> alphas;
  void plot(Canvas &, const math::Vector &v, const Pixel &, float,
            float alpha) {
    plotted.push_back(v);
    alphas.push_back(alpha);
  }
  void plot(Canvas &, float, float, const Pixel &, float, float) {}
};

/** @brief Identity fragment shader (leaves the fragment untouched). */
inline void noop_shader(const math::Vector &, Fragment &) {}

/**
 * @brief Largest angular gap (radians) between consecutive recorded positions.
 * @param pts Plotted positions in plot order.
 * @param wrap When true, also measures the gap from the last point back to the
 *             first (closed-loop seam continuity).
 */
inline float max_consecutive_gap(const std::vector<math::Vector> &pts,
                                 bool wrap) {
  float worst = 0.0f;
  for (size_t i = 1; i < pts.size(); ++i)
    worst = std::max(worst, math::angle_between(pts[i - 1], pts[i]));
  if (wrap && pts.size() >= 2)
    worst = std::max(worst, math::angle_between(pts.back(), pts.front()));
  return worst;
}

/**
 * @brief Largest seam-wrapped framebuffer-space gap between raster samples.
 * @details Longitude is ignored inside the two-row pole cap, where azimuth is
 * undefined and all columns meet.
 */
template <int W, int H>
inline float max_projected_gap(const std::vector<math::Vector> &points) {
  float worst = 0.0f;
  for (size_t i = 1; i < points.size(); ++i) {
    const math::PixelCoords a = math::vector_to_pixel<W, H>(points[i - 1]);
    const math::PixelCoords b = math::vector_to_pixel<W, H>(points[i]);
    float dx = std::abs(a.x - b.x);
    dx = std::min(dx, static_cast<float>(W) - dx);
    if (a.y < 2.0f || b.y < 2.0f || a.y > H - 3.0f || b.y > H - 3.0f)
      dx = 0.0f;
    worst = std::max(worst, std::hypot(dx, a.y - b.y));
  }
  return worst;
}
