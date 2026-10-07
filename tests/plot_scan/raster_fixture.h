/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Headless rasterize() harness.
// ---------------------------------------------------------------------------

/**
 * @brief Pipeline stub: records each plotted position; ignores color/age.
 * @details Carries both plot() overloads so it can also back a type-erased
 * PipelineRef; only the 3D world-space overload records.
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

/** @brief Pixel counts from a clipped render-band/reference comparison. */
struct RenderBandDiff {
  int lit = 0;            /**< Lit reference pixels in the render band. */
  int margin_lit = 0;     /**< Lit reference pixels outside the display band. */
  int row_margin_lit = 0; /**< Lit reference pixels above or below display. */
  int diff = 0;           /**< Pixels differing from the reference. */
  int margin_diff = 0;    /**< Differing pixels outside the display band. */
  int first_x = -1;       /**< Column of the first difference, or -1. */
  int first_y = -1;       /**< Row of the first difference, or -1. */
};

/** @brief Initializes every pixel in the pending frame to black. */
template <int W, int H> inline void initialize_parity_frame(Canvas &canvas) {
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x)
      canvas(x, y) = Pixel{};
}

/** @brief Compares a frame with its full-canvas reference over the render band. */
template <int W>
inline RenderBandDiff render_band_diff(const hs_test::StubEffect &fx,
                                       const std::vector<Pixel> &reference) {
  RenderBandDiff out;
  const ClipRegion &clip = fx.clip();
  for (int y = clip.render_y_start(); y < clip.render_y_end(); ++y) {
    for (int x = 0; x < W; ++x) {
      if (!clip.contains_x(x))
        continue;
      const Pixel &actual = fx.get_pixel(x, y);
      const Pixel &expected = reference[static_cast<size_t>(y) * W + x];
      const bool in_display_row = y >= clip.y_start && y < clip.y_end;
      const bool in_display =
          in_display_row && x >= clip.x_start && x < clip.x_end;
      const bool lit = (expected.r | expected.g | expected.b) != 0;
      out.lit += lit;
      out.margin_lit += !in_display && lit;
      out.row_margin_lit += !in_display_row && lit;
      if (actual == expected)
        continue;
      if (out.first_x < 0) {
        out.first_x = x;
        out.first_y = y;
      }
      ++out.diff;
      out.margin_diff += !in_display;
    }
  }
  return out;
}

/** @brief Reports one render-band parity result with its first mismatch. */
inline void expect_render_band_parity(const char *label,
                                      const RenderBandDiff &diff) {
  if (diff.diff != 0)
    std::printf("  [%s] diff=%d margin_diff=%d first=(%d,%d)\n", label,
                diff.diff, diff.margin_diff, diff.first_x, diff.first_y);
  HS_CONTEXT(label, diff.first_x, diff.first_y);
  HS_EXPECT_EQ(diff.diff, 0);
}
