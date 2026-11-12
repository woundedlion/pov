/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ---------------------------------------------------------------------------
// Headless rasterize() harness.
// ---------------------------------------------------------------------------

/** @brief Scoped override of the raster step budget. */
struct ScopedStepBudget {
  decltype(Plot::g_step_budget_override) saved = Plot::g_step_budget_override;
  explicit ScopedStepBudget(decltype(saved) budget) {
    Plot::g_step_budget_override = budget;
  }
  ~ScopedStepBudget() { Plot::g_step_budget_override = saved; }
};

/** @brief One planar star or flower contour to stroke. */
struct PlanarShapeCase {
  math::Quaternion orientation;
  float radius;
  int sides;
  float phase;
};

/** @brief Whether any bilinear screen tap of a world point reaches the clip. */
template <int W, int H>
inline bool taps_in_clip(const ClipRegion &clip, const math::Vector &point) {
  const auto screen = math::vector_to_pixel<W, H>(point);
  const int x = static_cast<int>(floorf(screen.x));
  const int y = static_cast<int>(floorf(screen.y));
  for (int dy = 0; dy <= 1; ++dy)
    for (int dx = 0; dx <= 1; ++dx)
      if (clip.contains_y(y + dy) && clip.contains_x((x + dx + W) % W))
        return true;
  return false;
}

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

/** @brief Binds a fragment trail from positions in the caller's plot arena.
 * @param out Destination trail.
 * @param points Positions copied into the fragments. */
inline void bind_trail(Fragments &out, std::span<const math::Vector> points) {
  out.bind(plot_arena(), points.size());
  for (const auto &point : points) {
    Fragment fragment;
    fragment.pos = point;
    out.push_back(fragment);
  }
}

/**
 * @brief Fragment shader painting full white at a fixed alpha.
 * @param alpha Alpha written to every fragment.
 */
inline auto white_shader(float alpha) {
  return [alpha](const math::Vector &, Fragment &f) {
    f.color = Color4(Pixel(65535, 65535, 65535), alpha);
  };
}

/** @brief The four W x H quadrants as {y0, y1, x0, x1}, row-major. */
template <int W, int H>
constexpr std::array<std::array<int, 4>, 4> quadrant_clips() {
  return {{{0, H / 2, 0, W / 2},
           {0, H / 2, W / 2, W},
           {H / 2, H, 0, W / 2},
           {H / 2, H, W / 2, W}}};
}

/**
 * @brief Two opposite quadrants, a mid-width wedge and a seam-edge wedge as
 *        {y0, y1, x0, x1}.
 */
template <int W, int H>
constexpr std::array<std::array<int, 4>, 4> quadrant_wedge_clips() {
  return {{{0, H / 2, 0, W / 2},
           {H / 2, H, W / 2, W},
           {0, H, 10, 34},
           {0, H, 0, 8}}};
}

/**
 * @brief Captures a single-pass planar raster of white alpha-0.4 fragments.
 * @tparam Policy Sampling policy; SELECTABLE also enables balanced sampling.
 * @param points Polyline to rasterize.
 * @param basis Planar projection basis.
 * @return The recorded positions and alphas.
 */
template <int W, int H, Plot::RasterSamplingPolicy Policy>
inline AlphaCapturePipeline
capture_planar_alpha_policy(const Fragments &points, const math::Basis &basis) {
  hs_test::StubEffect fx(W, H);
  AlphaCapturePipeline pipeline;
  Canvas canvas(fx);
  Plot::rasterize<
      W, H, Plot::RasterConfig{.single_pass = true, .sampling_policy = Policy}>(
      pipeline, canvas, points, white_shader(0.4f),
      {.projection = Plot::RasterProjection::planar(basis),
       .balanced_sampling = Policy == Plot::RasterSamplingPolicy::SELECTABLE});
  return pipeline;
}

/**
 * @brief Smallest angle (radians) from @p p to any of @p points.
 * @return PI_F when @p points is empty.
 */
inline float nearest_angle(const math::Vector &p,
                           const std::vector<math::Vector> &points) {
  float nearest = math::PI_F;
  for (const math::Vector &q : points)
    nearest = std::min(nearest, math::angle_between(p, q));
  return nearest;
}

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
    worst = fold_worst(worst, std::hypot(dx, a.y - b.y));
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
      const bool lit = !is_black(expected);
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

/** @brief Compares clipped render bands with a full-canvas reference.
 * @tparam W Canvas width.
 * @tparam H Canvas height.
 * @tparam Render Render callback type.
 * @tparam N Clip count.
 * @param label Failure context.
 * @param render Draws into a stub effect.
 * @param clips Display rectangles.
 * @param minimum_lit Per-clip exclusive lit floor, or -1 to omit it.
 * @return Summed comparison counts. */
template <int W, int H, typename Render, size_t N>
inline RenderBandDiff
expect_clip_band_parity(const char *label, Render &&render,
                        const std::array<std::array<int, 4>, N> &clips,
                        int minimum_lit = -1) {
  std::vector<Pixel> reference(static_cast<size_t>(W) * H);
  {
    hs_test::StubEffect fx(W, H);
    render(fx);
    fx.advance_display();
    hs_test::capture_frame<W, H>(fx, reference);
  }
  RenderBandDiff total;
  for (const auto &clip : clips) {
    hs_test::StubEffect fx(W, H);
    fx.set_clip(clip[0], clip[1], clip[2], clip[3]);
    render(fx);
    fx.advance_display();
    const RenderBandDiff diff = render_band_diff<W>(fx, reference);
    expect_render_band_parity(label, diff);
    if (minimum_lit >= 0)
      HS_EXPECT_GT(diff.lit, minimum_lit);
    total.lit += diff.lit;
    total.margin_lit += diff.margin_lit;
    total.row_margin_lit += diff.row_margin_lit;
    total.diff += diff.diff;
    total.margin_diff += diff.margin_diff;
  }
  return total;
}

/** @brief Counts qualifying pixels lacking a wrapped 3x3 neighbour.
 * @param from Source frame.
 * @param to Neighbour frame.
 * @param qualifies Source-pixel predicate.
 * @param lit Neighbour-pixel predicate.
 * @param y0,y1 Row interval.
 * @param x0,x1 Column interval.
 * @param pole_rows Rows excluded at each pole.
 * @return Number of uncovered source pixels. */
template <int W, int H, typename Element, typename Qualifies, typename Lit>
inline size_t pixels_without_lit_neighbor(const std::vector<Element> &from,
                                          const std::vector<Element> &to,
                                          Qualifies qualifies, Lit lit,
                                          int y0 = 0, int y1 = H, int x0 = 0,
                                          int x1 = W, int pole_rows = 0) {
  size_t count = 0;
  for (int y = std::max(y0, pole_rows); y < std::min(y1, H - pole_rows); ++y)
    for (int x = x0; x < x1; ++x) {
      if (!qualifies(from[static_cast<size_t>(y) * W + x]))
        continue;
      bool nearby = false;
      for (int dy = -1; dy <= 1 && !nearby; ++dy)
        for (int dx = -1; dx <= 1 && !nearby; ++dx) {
          const int sy = y + dy;
          if (sy >= 0 && sy < H)
            nearby = lit(to[static_cast<size_t>(sy) * W + (x + dx + W) % W]);
        }
      count += !nearby;
    }
  return count;
}
