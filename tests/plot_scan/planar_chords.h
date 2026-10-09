/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Plot::PlanarChords  — chord-walked planar polylines
// ============================================================================

/** @brief One planar star or flower contour to stroke. */
struct PlanarShapeCase {
  math::Quaternion orientation;
  float radius;
  int sides;
  float phase;
};

/**
 * @brief Strokes @p star into @p fx, through PlanarChords or, for the
 *        reference, a balanced Plot::rasterize of the same closed polyline.
 */
template <int W, int H>
inline std::vector<Pixel> render_planar_chord_star(hs_test::StubEffect &fx,
                                                   const PlanarShapeCase &star,
                                                   bool chords) {
  using Star = Plot::Star<Plot::PlanarProjection>;
  const math::Basis basis = math::make_basis(star.orientation, math::X_AXIS);
  const Color4 color(Pixel(65535, 65535, 65535), 0.3f);
  auto shader = [&](const math::Vector &, Fragment &f) { f.color = color; };
  {
    ScratchScope sc(plot_arena());
    Plot::PlanarChords<W, H> planar_chords;
    planar_chords.init_storage(plot_arena(), star.sides * 2);
    Fragments points;
    points.bind(plot_arena(), static_cast<size_t>(star.sides * 2 + 2));
    math::Basis projection_basis;
    const math::Basis &planar_basis = *Plot::PlanarProjection::edge_basis(
        basis, star.radius, projection_basis);
    Star::sample_chart_positions(
        points, planar_chords.chart_x(), planar_chords.chart_y(), basis,
        star.radius, star.sides, star.phase, Star::radius_trig(star.radius),
        Star::step_trig(star.sides), planar_basis);
    for (int i = 0; i < star.sides * 2; ++i) {
      const float x = math::dot(points[i].pos, planar_basis.u);
      const float y = math::dot(points[i].pos, planar_basis.w);
      const float radial = hypotf(x, y);
      const float angle =
          atan2f(radial, math::dot(points[i].pos, planar_basis.v));
      HS_EXPECT_NEAR(planar_chords.chart_x()[i], x * angle / radial, 2e-5f);
      HS_EXPECT_NEAR(planar_chords.chart_y()[i], y * angle / radial, 2e-5f);
    }
    Filter::Screen::DirectAntiAliasSink<W, H> sink;
    Canvas canvas(fx);
    initialize_parity_frame<W, H>(canvas);
    sink.prepare(canvas);
    if (chords) {
      planar_chords.prepare(canvas.clip());
      planar_chords.draw_closed(sink, canvas, points, star.sides * 2,
                                planar_basis, color, shader);
    } else {
      Plot::rasterize<W, H, Plot::PLANAR_CHORD_RASTER_CONFIG>(
          sink, canvas, points, shader,
          {.projection = Plot::RasterProjection::planar(planar_basis),
           .omit_end = true,
           .balanced_sampling = true});
    }
  }
  fx.advance_display();
  std::vector<Pixel> frame(static_cast<size_t>(W) * H);
  hs_test::capture_frame<W, H>(fx, frame);
  return frame;
}

inline uint64_t planar_chord_energy(const std::vector<Pixel> &frame) {
  uint64_t energy = 0;
  for (const Pixel &p : frame)
    energy += static_cast<uint64_t>(p.r) + p.g + p.b;
  return energy;
}

template <int W, int H>
inline void expect_covered_within_one_pixel(const std::vector<Pixel> &reference,
                                            const std::vector<Pixel> &candidate,
                                            int y0, int y1, int x0, int x1,
                                            int pole_rows, double max_drift) {
  uint64_t reference_energy = 0, candidate_energy = 0;
  size_t uncovered = 0;
  for (int y = y0; y < y1; ++y)
    for (int x = x0; x < x1; ++x) {
      const size_t i = static_cast<size_t>(y) * W + x;
      const Pixel &p = reference[i];
      const Pixel &q = candidate[i];
      reference_energy += static_cast<uint64_t>(p.r) + p.g + p.b;
      candidate_energy += static_cast<uint64_t>(q.r) + q.g + q.b;
      // Pole rows turn a sub-pixel sample-phase shift into whole columns.
      if (static_cast<uint32_t>(p.r) + p.g + p.b < 12288 || y < pole_rows ||
          y >= H - pole_rows)
        continue;
      bool covered = false;
      for (int dy = -1; dy <= 1 && !covered; ++dy)
        for (int dx = -1; dx <= 1 && !covered; ++dx) {
          const int sy = y + dy;
          if (sy < 0 || sy >= H)
            continue;
          covered = !is_black(
              candidate[static_cast<size_t>(sy) * W + (x + dx + W) % W]);
        }
      uncovered += !covered;
    }
  HS_EXPECT_EQ(uncovered, size_t{0});
  if (reference_energy == 0) {
    HS_EXPECT_EQ(candidate_energy, uint64_t{0});
    return;
  }
  const double drift = std::fabs(static_cast<double>(candidate_energy) -
                                 static_cast<double>(reference_energy)) /
                       static_cast<double>(reference_energy);
  HS_EXPECT_LT(drift, max_drift);
}

/** @brief Stars away from the poles, near a pole, and past the equator. */
inline std::array<PlanarShapeCase, 4> planar_chord_stars() {
  return {{
      {math::Quaternion(0.81f, 0.32f, -0.29f, 0.39f).normalized(), 0.45f, 7,
       0.3f},
      {math::Quaternion(0.72f, -0.41f, 0.18f, 0.53f).normalized(), 0.8f, 5,
       1.1f},
      {math::make_rotation(math::X_AXIS, math::Y_AXIS), 0.12f, 7, 0.6f},
      {math::Quaternion(0.93f, -0.11f, 0.24f, 0.25f).normalized(), 1.35f, 9,
       2.0f},
  }};
}

/**
 * @brief Verifies PlanarChords strokes a star as bright as the balanced
 *        adaptive walk and retains bright reference coverage within one pixel.
 */
inline void test_planar_chords_match_rasterize_brightness() {
  constexpr int W = 288, H = 144;
  hs_test::StubEffect fx(W, H);
  for (const PlanarShapeCase &star : planar_chord_stars()) {
    const auto reference = render_planar_chord_star<W, H>(fx, star, false);
    const auto chords = render_planar_chord_star<W, H>(fx, star, true);
    HS_EXPECT_GT(planar_chord_energy(reference), uint64_t{0});
    expect_covered_within_one_pixel<W, H>(reference, chords, 0, H, 0, W, 0,
                                          0.026);
  }
}

/** @brief What one band-split render of a flower produced. */
struct BandSplitFrame {
  std::vector<Pixel> pixels;
  size_t skipped_edges = 0; /**< Split edges flagged invisible. */
};

/**
 * @brief Rasterizes @p flower into @p fx, whole or through PlanarBandSplit
 *        against the clip @p clip.
 */
template <int W, int H>
inline BandSplitFrame
render_band_split_flower(hs_test::StubEffect &fx, const PlanarShapeCase &flower,
                         const ClipRegion &clip, bool split) {
  constexpr int PIECES = 8;
  const math::Basis basis = math::make_basis(flower.orientation, math::X_AXIS);
  const math::Basis planar_basis = Plot::planar_chart_basis(
      math::get_antipode(basis, flower.radius).first.v);
  const Color4 color(Pixel(65535, 65535, 65535), 0.3f);
  auto shader = [&](const math::Vector &, Fragment &f) { f.color = color; };
  fx.set_clip(clip.y_start, clip.y_end, clip.x_start, clip.x_end);
  BandSplitFrame frame;
  {
    ScratchScope sc(plot_arena());
    const int edges = flower.sides * 2;
    const int max_points =
        Plot::PlanarBandSplit<W, H>::max_points(edges, PIECES);
    Plot::PlanarBandSplit<W, H> band_split;
    band_split.init_storage(plot_arena(), max_points);
    Fragments ring;
    ring.bind(plot_arena(), static_cast<size_t>(edges + 1));
    Plot::Flower::sample(ring, basis, flower.radius, flower.sides,
                         flower.phase);
    Filter::Screen::DirectAntiAliasSink<W, H> sink;
    Canvas canvas(fx);
    initialize_parity_frame<W, H>(canvas);
    sink.prepare(canvas);
    constexpr Plot::RasterConfig CONFIG{
        .single_pass = true,
        .derive_planar_arc_registers = false,
        .interpolate_registers = false,
        .sampling_policy = Plot::RasterSamplingPolicy::SELECTABLE};
    if (split) {
      Fragments path;
      path.bind(plot_arena(), static_cast<size_t>(max_points));
      const auto flags =
          band_split.split(path, ring, edges, PIECES, planar_basis,
                           Plot::ClipBand<W, H>::of(canvas.clip()));
      for (uint8_t flag : flags)
        frame.skipped_edges += flag == 0;
      Plot::rasterize<W, H, CONFIG>(
          sink, canvas, path, shader,
          {.projection = Plot::RasterProjection::planar(planar_basis, flags),
           .omit_end = true});
    } else {
      Plot::rasterize<W, H, CONFIG>(
          sink, canvas, ring, shader,
          {.projection = Plot::RasterProjection::planar(planar_basis),
           .omit_end = true});
    }
  }
  fx.advance_display();
  frame.pixels.resize(static_cast<size_t>(W) * H);
  hs_test::capture_frame<W, H>(fx, frame.pixels);
  return frame;
}

/**
 * @brief Verifies a band-split flower draws each quadrant as the whole flower
 *        does: the same strokes to a fraction of a pixel, with pieces that
 *        cannot reach the band skipped.
 */
inline void test_planar_band_split_matches_whole_polyline() {
  constexpr int W = 288, H = 144;
  constexpr int POLE_ROWS = 3;
  hs_test::StubEffect fx(W, H);
  const ClipRegion full{0, H, 0, W, 1, W, H};
  const ClipRegion quadrants[] = {
      {0, H / 2, 0, W / 2, 1, W, H},
      {0, H / 2, W / 2, W, 1, W, H},
      {H / 2, H, 0, W / 2, 1, W, H},
      {H / 2, H, W / 2, W, 1, W, H},
  };
  const std::array<PlanarShapeCase, 4> flowers = {{
      {math::Quaternion(0.81f, 0.32f, -0.29f, 0.39f).normalized(), 0.6f, 3,
       0.4f},
      {math::Quaternion(0.72f, -0.41f, 0.18f, 0.53f).normalized(), 1.2f, 3,
       1.3f},
      {math::make_rotation(math::X_AXIS, math::Y_AXIS), 0.9f, 5, 0.2f},
      {math::Quaternion(0.93f, -0.11f, 0.24f, 0.25f).normalized(), 1.6f, 4,
       2.1f},
  }};
  size_t skipped = 0;
  for (const PlanarShapeCase &flower : flowers) {
    const auto whole =
        render_band_split_flower<W, H>(fx, flower, full, false).pixels;
    for (const ClipRegion &clip : quadrants) {
      const auto tile = render_band_split_flower<W, H>(fx, flower, clip, true);
      skipped += tile.skipped_edges;
      expect_covered_within_one_pixel<W, H>(whole, tile.pixels, clip.y_start,
                                            clip.y_end, clip.x_start,
                                            clip.x_end, POLE_ROWS, 0.035);
    }
  }
  HS_EXPECT_GT(skipped, size_t{0});
}

/**
 * @brief Verifies PlanarChords' band-split pole runs draw each quadrant as the
 *        whole stroke does, to a fraction of a pixel.
 */
inline void test_planar_chords_pole_split_matches_whole_star() {
  constexpr int W = 288, H = 144;
  constexpr int POLE_ROWS = 3;
  hs_test::StubEffect fx(W, H);
  const int quadrants[4][4] = {{0, H / 2, 0, W / 2},
                               {0, H / 2, W / 2, W},
                               {H / 2, H, 0, W / 2},
                               {H / 2, H, W / 2, W}};
  const std::array<PlanarShapeCase, 5> stars = {{
      {math::make_rotation(math::X_AXIS, math::Y_AXIS), 0.12f, 7, 0.6f},
      {math::make_rotation(math::X_AXIS, math::Y_AXIS), 0.35f, 7, 1.7f},
      {math::make_rotation(math::X_AXIS, -math::Y_AXIS), 0.22f, 5, 0.9f},
      {math::Quaternion(0.72f, -0.41f, 0.18f, 0.53f).normalized(), 1.8f, 7,
       2.4f},
      {math::Quaternion(0.93f, -0.11f, 0.24f, 0.25f).normalized(), 0.3f, 9,
       0.2f},
  }};
  size_t split_changed = 0;
  for (const PlanarShapeCase &star : stars) {
    fx.set_clip(0, H, 0, W);
    const auto whole = render_planar_chord_star<W, H>(fx, star, true);
    for (const auto &q : quadrants) {
      fx.set_clip(q[0], q[1], q[2], q[3]);
      Plot::g_planar_chords_split_pole_runs = false;
      const auto unsplit = render_planar_chord_star<W, H>(fx, star, true);
      Plot::g_planar_chords_split_pole_runs = true;
      const auto tile = render_planar_chord_star<W, H>(fx, star, true);
      for (int y = q[0]; y < q[1]; ++y)
        for (int x = q[2]; x < q[3]; ++x) {
          const size_t i = static_cast<size_t>(y) * W + x;
          split_changed += !(tile[i] == unsplit[i]);
        }
      expect_covered_within_one_pixel<W, H>(whole, tile, q[0], q[1], q[2], q[3],
                                            POLE_ROWS, 0.005);
    }
  }
  fx.set_clip(0, H, 0, W);
  HS_EXPECT_GT(split_changed, size_t{0});
}
