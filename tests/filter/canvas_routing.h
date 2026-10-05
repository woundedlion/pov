/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// Included by tests/test_filter.h.

// ============================================================================
// End-to-end pipeline routing through a live Canvas
//
// A Canvas spins in its ctor while !buffer_free(), so every drawn frame must
// advance_display() before the next Canvas is constructed.
// ============================================================================

/**
 * @brief Verifies the bare 2D sink: int + float overloads, exact (alpha=1)
 *        write, x-wrap, and clip.
 */
inline void test_pipeline_sink_2d_plot_blends_wraps_clips() {
  constexpr int W = 16, H = 8;
  Pipeline<W, H> pipe;
  // fx and fx2 alias the same static double buffer (Effect's single-live guard),
  // so scope the first effect closed before constructing the second.
  {
    hs_test::StubEffect fx(W, H);
    {
      Canvas c(fx);
      // Integer overload, full alpha over a black buffer -> exact source colour.
      pipe.plot(c, 4, 3, Pixel(100, 200, 300), 0.0f, 1.0f);
      // x = W+2 is in the sink's [-W, 2W) contract window and wraps to column 2.
      pipe.plot(c, W + 2, 5, Pixel(10, 20, 30), 0.0f, 1.0f);
      // Float overload rounds to the nearest pixel.
      pipe.plot(c, 7.4f, 1.6f, Pixel(1, 2, 3), 0.0f, 1.0f);
    }
    fx.advance_display();
    HS_EXPECT_PIXEL(fx.get_pixel(4, 3), 100, 200, 300);
    HS_EXPECT_PIXEL(fx.get_pixel(2, 5), 10, 20, 30); // wrapped
    HS_EXPECT_PIXEL(fx.get_pixel(7, 2), 1, 2, 3);    // rounded
  }

  // A second frame: a plot outside the clip band is dropped.
  hs_test::StubEffect fx2(W, H);
  fx2.set_clip(2, 5, 0, W); // rows [2,5)
  fx2.set_margin(0);
  {
    Canvas c(fx2);
    pipe.plot(c, 8, 6, Pixel(500, 0, 0), 0.0f, 1.0f); // row 6 outside band
    pipe.plot(c, 8, 3, Pixel(0, 500, 0), 0.0f, 1.0f); // row 3 inside band
  }
  fx2.advance_display();
  HS_EXPECT_TRUE(is_black(fx2.get_pixel(8, 6)));
  HS_EXPECT_PIXEL(fx2.get_pixel(8, 3), 0, 500, 0);
}

inline void test_pipeline_composition_alpha_and_draw_order() {
  constexpr int W = 16, H = 8;
  for (bool reverse : {false, true}) {
    hs_test::StubEffect fx(W, H);
    Pipeline<W, H, Filter::Screen::AntiAlias<W, H>,
             Filter::Pixel::ChromaticShift<W>>
        pipe;
    const Pixel first = reverse ? Pixel(0, 0, 40000) : Pixel(40000, 0, 0);
    const Pixel second = reverse ? Pixel(40000, 0, 0) : Pixel(0, 0, 40000);
    {
      Canvas canvas(fx);
      pipe.plot(canvas, 4.0f, 3.5f, first, 0.0f, 1.0f);
      pipe.plot(canvas, 4.0f, 3.5f, second, 0.0f, 1.0f);
    }
    fx.advance_display();
    for (int y : {3, 4}) {
      const Pixel &center = fx.get_pixel(4, y);
      // Each antialias tap covers half a pixel: C = first/4 + second/2.
      HS_EXPECT_NEAR(center.r, reverse ? 20000 : 10000, 2);
      HS_EXPECT_NEAR(center.b, reverse ? 10000 : 20000, 2);
      HS_EXPECT_EQ(center.g, 0);
      HS_EXPECT_NEAR(fx.get_pixel(5, y).r, reverse ? 5000 : 4375, 2);
      HS_EXPECT_NEAR(fx.get_pixel(7, y).b, reverse ? 4375 : 5000, 2);
      HS_EXPECT_TRUE(is_black(fx.get_pixel(6, y)));
    }
    HS_EXPECT_EQ(count_lit_canvas(fx), size_t{6});
  }
}

/**
 * @brief Verifies the 3D sink overload routes a unit vector via vector_to_pixel
 *        and writes it.
 */
inline void test_pipeline_sink_3d_plot_routes_to_canvas() {
  constexpr int W = 32, H = 16;
  hs_test::StubEffect fx(W, H);
  Pipeline<W, H> pipe;
  math::Vector v = math::Vector(0.6f, 0.4f, 0.69f).normalized();
  {
    Canvas c(fx);
    pipe.plot(c, v, Pixel(40000, 20000, 10000), 0.0f, 1.0f);
  }
  fx.advance_display();

  // Exactly one pixel lit, carrying the source colour, at the coordinate
  // vector_to_pixel predicts.
  HS_EXPECT_EQ(count_lit_canvas(fx), (size_t)1);
  math::PixelCoords pc = math::vector_to_pixel<W, H>(v);
  int ex = math::fast_wrap(static_cast<int>(std::round(pc.x)), W);
  int ey = static_cast<int>(std::round(pc.y));
  HS_EXPECT_PIXEL(fx.get_pixel(ex, ey), 40000, 20000, 10000);
}

/**
 * @brief Verifies a World filter (3D) head fans out through the 3D->2D sink
 *        conversion.
 */
inline void test_pipeline_world_replicate_fans_out() {
  constexpr int W = 32, H = 16;
  // Replicate(2): original + one copy rotated 180 deg about Y (same latitude,
  // longitude + W/2) -> two distinct columns -> two lit pixels.
  Pipeline<W, H, Filter::World::Replicate<W>> pipe(
      Filter::World::Replicate<W>(2));
  math::Vector v = math::Vector(0.6f, 0.4f, 0.69f).normalized();
  {
    hs_test::StubEffect fx(W, H);
    {
      Canvas c(fx);
      pipe.plot(c, v, Pixel(60000, 60000, 60000), 0.0f, 1.0f);
    }
    fx.advance_display();
    HS_EXPECT_EQ(count_lit_canvas(fx), (size_t)2);
  }

  pipe.get<Filter::World::Replicate<W>>().set_count(3);
  {
    hs_test::StubEffect fx(W, H);
    {
      Canvas c(fx);
      pipe.plot(c, v, Pixel(60000, 60000, 60000), 0.0f, 1.0f);
    }
    fx.advance_display();
    HS_EXPECT_EQ(count_lit_canvas(fx), (size_t)3);
  }

  pipe.get<Filter::World::Replicate<W>>().set_count(0);
  {
    hs_test::StubEffect fx(W, H);
    {
      Canvas c(fx);
      pipe.plot(c, v, Pixel(60000, 60000, 60000), 0.0f, 1.0f);
    }
    fx.advance_display();
    HS_EXPECT_EQ(count_lit_canvas(fx), (size_t)1);
  }
}

/**
 * @brief Verifies a 2D coordinate into a 3D-headed pipeline round-trips through
 *        the float-plot mismatch branch.
 * @details The mismatch branch maps pixel_to_vector, the identity filter passes
 *          it, and the sink maps it back.
 */
inline void test_pipeline_2d_into_3d_head_roundtrips() {
  constexpr int W = 32, H = 16;
  hs_test::StubEffect fx(W, H);
  // Replicate(1) clamps to a single emission -> identity pass-through.
  Pipeline<W, H, Filter::World::Replicate<W>> pipe(
      Filter::World::Replicate<W>(1));
  const int px = 16, py = 8; // interior, away from poles/seam
  {
    Canvas c(fx);
    pipe.plot(c, static_cast<float>(px), static_cast<float>(py),
              Pixel(0, 60000, 0), 0.0f, 1.0f);
  }
  fx.advance_display();

  // One pixel lit, within a pixel of the input after the trig round-trip.
  HS_EXPECT_EQ(count_lit_canvas(fx), (size_t)1);
  bool near_input = false;
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x)
      if (!is_black(fx.get_pixel(x, y)))
        near_input |= (std::abs(x - px) <= 1 && std::abs(y - py) <= 1);
  HS_EXPECT_TRUE(near_input);
}

/**
 * @brief Verifies a Screen filter (2D) head forwards to the sink.
 * @details An integer coordinate collapses AntiAlias to a single tap, so the
 *          routed result is one exact pixel.
 */
inline void test_pipeline_screen_antialias_routes_to_sink() {
  constexpr int W = 32, H = 16;
  hs_test::StubEffect fx(W, H);
  Pipeline<W, H, Filter::Screen::AntiAlias<W, H>> pipe;
  {
    Canvas c(fx);
    pipe.plot(c, 10.0f, 5.0f, Pixel(50000, 0, 25000), 0.0f, 1.0f);
  }
  fx.advance_display();
  HS_EXPECT_EQ(count_lit_canvas(fx), (size_t)1);
  HS_EXPECT_PIXEL(fx.get_pixel(10, 5), 50000, 0, 25000);
}

struct AASinkSample {
  float x;
  float y;
  Pixel color;
  float alpha;
};

template <typename Sink>
inline std::vector<Pixel>
render_aa_sink_case(int w, int h, int y0, int y1, int x0, int x1, int margin,
                    const std::vector<AASinkSample> &samples) {
  std::vector<Pixel> frame;
  {
    hs_test::StubEffect fx(w, h);
    fx.set_clip(y0, y1, x0, x1);
    fx.set_margin(margin);
    Sink sink;
    {
      Canvas cv(fx);
      if constexpr (requires { sink.prepare(cv); })
        sink.prepare(cv);
      for (int y = 0; y < h; ++y)
        for (int x = 0; x < w; ++x) {
          const uint32_t n = static_cast<uint32_t>(y * w + x);
          cv(x, y) = Pixel(static_cast<uint16_t>(n * 977u + 13u),
                           static_cast<uint16_t>(n * 313u + 101u),
                           static_cast<uint16_t>(n * 701u + 7u));
        }
      for (const AASinkSample &s : samples)
        sink.plot(cv, s.x, s.y, s.color, 0.0f, s.alpha);
    }
    fx.advance_display();
    frame.reserve(static_cast<size_t>(w * h));
    for (int y = 0; y < h; ++y)
      for (int x = 0; x < w; ++x)
        frame.push_back(fx.get_pixel(x, y));
  }
  return frame;
}

inline void test_direct_antialias_sink_stale_clip() {
  constexpr int W = 17, H = 9;
  hs_test::StubEffect fx(W, H);
  Filter::Screen::DirectAntiAliasSink<W, H> sink;
  const ::Pixel *original = nullptr;
  {
    Canvas cv(fx);
    sink.prepare(cv);
    original = cv.data();
    HS_EXPECT_TRUE(sink.prepared_for(cv));
  }
  fx.advance_display();
  {
    Canvas cv(fx);
    HS_EXPECT_NE(cv.data(), original);
    HS_EXPECT_FALSE(sink.prepared_for(cv));
  }
  fx.advance_display();
  fx.set_clip(2, H - 2, 2, W - 2);
  {
    Canvas cv(fx);
    HS_EXPECT_EQ(cv.data(), original);
    HS_EXPECT_FALSE(sink.prepared_for(cv));
    sink.prepare(cv);
    HS_EXPECT_TRUE(sink.prepared_for(cv));
  }
  fx.advance_display();
}

/**
 * @brief Proves the opt-in direct AA sink is framebuffer-identical to the
 * generic AntiAlias pipeline across poles, seams, clips and random splats, and
 * that every clip case deposits samples into the frame it compares.
 */
inline void test_direct_antialias_sink_framebuffer_parity() {
  constexpr int W = 17;
  constexpr int H = 9;
  using Reference = Pipeline<W, H, Filter::Screen::AntiAlias<W, H>>;
  using Direct = Filter::Screen::DirectAntiAliasSink<W, H>;
  static_assert(!Direct::has_world_cull);
  static_assert(!Direct::any_crosses_segments);

  std::vector<AASinkSample> samples;
  constexpr std::array<float, 12> xs{-17.0f,  -16.999f, -0.999f, -0.001f,
                                     0.0f,    0.001f,   7.5f,    15.999f,
                                     16.999f, 17.001f,  33.001f, 33.999f};
  constexpr std::array<float, 11> ys{-1.001f, -1.0f,  -0.999f, -0.001f,
                                     0.0f,    0.001f, 4.5f,    7.999f,
                                     8.0f,    8.999f, 9.001f};
  constexpr std::array<float, 5> alphas{0.0f, 0.00001f, 0.25f, 1.0f, 1.2f};
  size_t sequence = 0;
  for (float x : xs)
    for (float y : ys) {
      const uint16_t n = static_cast<uint16_t>(sequence++ * 4051u);
      samples.push_back({x, y,
                         Pixel(n, static_cast<uint16_t>(n * 3u),
                               static_cast<uint16_t>(65535u - n)),
                         alphas[sequence % alphas.size()]});
    }

  constexpr std::array<float, 5> interior_fractions{0.1f, 0.25f, 0.5f, 0.75f,
                                                    0.9f};
  for (int yi = 1; yi < H - 1; ++yi) {
    for (int xi = 0; xi < W; ++xi) {
      for (float fraction : interior_fractions) {
        const uint16_t n = static_cast<uint16_t>(sequence++ * 4051u);
        samples.push_back({static_cast<float>(xi) + fraction,
                           static_cast<float>(yi) + (1.0f - fraction),
                           Pixel(n, static_cast<uint16_t>(n * 3u),
                                 static_cast<uint16_t>(65535u - n)),
                           alphas[sequence % alphas.size()]});
      }
    }
  }

  hs::Pcg32 rng(0x6d2b79f5u);
  auto random_unit = [&rng]() { return rand_uniform(rng, 0.0f, 1.0f); };
  for (int i = 0; i < 512; ++i) {
    const float x = -W + random_unit() * (3.0f * W - 0.002f);
    const float y = -1.1f + random_unit() * (H + 1.2f);
    const uint16_t r = static_cast<uint16_t>(rng());
    const uint16_t g = static_cast<uint16_t>(rng());
    const uint16_t b = static_cast<uint16_t>(rng());
    samples.push_back({x, y, Pixel(r, g, b), random_unit() * 1.25f});
  }

  struct ClipCase {
    int y0, y1, x0, x1, margin;
  };
  constexpr std::array<ClipCase, 7> clips{{
      {0, H, 0, W, 0},
      {2, 7, 3, 13, 0},
      {2, 7, 3, 13, 1},
      {2, 7, 0, 2, 2},
      {2, 7, 15, W, 2},
      {4, 5, 8, 9, 4},
      {0, H, 4, 12, W - 1},
  }};

  for (const ClipCase &clip : clips) {
    const auto background = render_aa_sink_case<Reference>(
        W, H, clip.y0, clip.y1, clip.x0, clip.x1, clip.margin, {});
    const auto expected = render_aa_sink_case<Reference>(
        W, H, clip.y0, clip.y1, clip.x0, clip.x1, clip.margin, samples);
    const auto actual = render_aa_sink_case<Direct>(
        W, H, clip.y0, clip.y1, clip.x0, clip.x1, clip.margin, samples);
    HS_EXPECT_EQ(actual.size(), expected.size());
    HS_EXPECT_EQ(background.size(), expected.size());
    int changed = 0;
    for (size_t i = 0; i < expected.size(); ++i) {
      HS_EXPECT_EQ(actual[i], expected[i]);
      if (expected[i] != background[i])
        ++changed;
    }
    HS_EXPECT_GT(changed, 0);
  }
}

/**
 * @brief Verifies the Pixel::Feedback::flush warp-field path.
 * @details A Smoke style with no bound NoiseParams makes noise_warp an identity
 *          map, so the flush blends the previous frame back at the same location,
 *          faded by style.fade — deterministic.
 */
inline void test_feedback_flush_blends_prev_frame() {
  constexpr int W = 32, H = 16; // both divisible by Smoke's downsample (4)
  hs_test::StubEffect fx(W, H);
  ::Feedback::Style style = ::Feedback::Style::Smoke(); // noise stays nullptr
  Pipeline<W, H, Filter::Pixel::Feedback<W, H>> pipe{
      Filter::Pixel::Feedback<W, H>(style)};

  // Frame 1: a bright horizontal band becomes the "previous" frame.
  {
    Canvas c(fx);
    for (int y = 6; y <= 9; ++y)
      for (int x = 0; x < W; ++x)
        c(x, y) = Pixel(40000, 40000, 40000);
  }
  fx.advance_display();

  // Frame 2: empty buffer; flush warps+blends the prev band into it.
  {
    Canvas c(fx);
    (void)pipe.begin_frame(c, 1.0f);
  }
  fx.advance_display();

  // Interior band row carried over and was faded (fade=0.9 -> ~36000), not the
  // original brightness. A row far from the band stays black (identity warp,
  // no spill that far).
  const Pixel &band = fx.get_pixel(W / 2, 8);
  HS_EXPECT_FALSE(is_black(band));
  HS_EXPECT_NEAR(band.r, 36000, 400);
  HS_EXPECT_TRUE(is_black(fx.get_pixel(W / 2, 0)));

  // Disabled feedback short-circuits flush: a fresh frame stays black.
  pipe.get<Filter::Pixel::Feedback<W, H>>().set_enabled(false);
  {
    Canvas c(fx);
    (void)pipe.begin_frame(c, 1.0f);
  }
  fx.advance_display();
  HS_EXPECT_TRUE(is_black(fx.get_pixel(W / 2, 8)));
}

/**
 * @brief Asserts an aliased pole row carries one physical sample's color.
 * @param fx Effect holding the composited frame.
 * @param w Canvas width.
 * @param y Pole row index.
 * @param pole Color seeded into the pole row before the flush.
 * @details The callers use unity-fade PLAIN composition, whose quantize16
 * preserves the seeded u16 channels.
 */
inline void expect_pole_row_collapsed(hs_test::StubEffect &fx, int w, int y,
                                      const Pixel &pole) {
  for (int x = 0; x < w; ++x) {
    const Pixel px = fx.get_pixel(x, y);
    HS_EXPECT_EQ(px.r, pole.r);
    HS_EXPECT_EQ(px.g, pole.g);
    HS_EXPECT_EQ(px.b, pole.b);
  }
}

/**
 * @brief Verifies the aliased north-pole row retains a lone physical sample.
 */
inline void test_feedback_north_pole_uses_single_physical_sample() {
  constexpr int W = 32, H = 16;
  constexpr Pixel POLE(12000, 30000, 50000);
  hs_test::StubEffect fx(W, H);

  ::Feedback::Style style{};
  style.noise = nullptr; // unbound noise_warp is the identity
  style.fade = 1.0f;
  style.downsample = 4;
  Pipeline<W, H, Filter::Pixel::Feedback<W, H>> pipe{
      Filter::Pixel::Feedback<W, H>(style)};

  {
    Canvas c(fx);
    c(W / 3, 0) = POLE;
  }
  fx.advance_display();
  {
    Canvas c(fx);
    (void)pipe.begin_frame(c, 1.0f);
  }
  fx.advance_display();

  expect_pole_row_collapsed(fx, W, 0, POLE);
  for (int x = 0; x < W; ++x)
    HS_EXPECT_TRUE(is_black(fx.get_pixel(x, 1)));
}

/**
 * @brief Verifies the aliased south-pole row retains a lone physical sample.
 * @details Row H-1 collapses to the south pole in the ideal display profile
 *          with H_OFFSET == 0.
 */
inline void test_feedback_south_pole_uses_single_physical_sample() {
  constexpr int W = 32, H = 16;
  constexpr Pixel POLE(12000, 30000, 50000);
  hs_test::StubEffect fx(W, H);

  ::Feedback::Style style{};
  style.noise = nullptr; // unbound noise_warp is the identity
  style.fade = 1.0f;
  style.downsample = 4;
  Pipeline<W, H, Filter::Pixel::Feedback<W, H>> pipe{
      Filter::Pixel::Feedback<W, H>(style)};

  {
    Canvas c(fx);
    c(W / 3, H - 1) = POLE;
  }
  fx.advance_display();
  {
    Canvas c(fx);
    (void)pipe.begin_frame(c, 1.0f);
  }
  fx.advance_display();

  expect_pole_row_collapsed(fx, W, H - 1, POLE);
  for (int x = 0; x < W; ++x)
    HS_EXPECT_TRUE(is_black(fx.get_pixel(x, H - 2)));
}

/**
 * @brief Verifies polar reconstruction suppresses longitude aliasing.
 * @details Rows 1 and H-2 sit in the dense infill bands and reconstruct at
 *          their spherical footprint; rows 4 and H-5 lie under the half-
 *          resolution latitude outside that band, where each column pair
 *          composites as its box average, so the stripe flattens to its mean
 *          there too; rows 8 and H-9 have a one-pixel footprint and pass the
 *          stripe through.
 */
inline void test_feedback_polar_rows_use_spherical_footprint() {
  constexpr int W = 64, H = 34;
  hs_test::StubEffect fx(W, H);

  ::Feedback::Style style{};
  style.noise = nullptr; // unbound noise_warp is the identity
  style.fade = 1.0f;
  style.downsample = 4;
  Pipeline<W, H, Filter::Pixel::Feedback<W, H>> pipe{
      Filter::Pixel::Feedback<W, H>(style)};

  {
    Canvas c(fx);
    for (int x = 0; x < W; ++x) {
      const Pixel stripe((x & 1) ? 60000 : 0, 0, 0);
      c(x, 1) = stripe;
      c(x, 4) = stripe;
      c(x, 8) = stripe;
      c(x, H - 2) = stripe;
      c(x, H - 5) = stripe;
      c(x, H - 9) = stripe;
    }
  }
  fx.advance_display();
  {
    Canvas c(fx);
    (void)pipe.begin_frame(c, 1.0f);
  }
  fx.advance_display();

  auto expect_reconstructed = [&](int y) {
    uint16_t polar_min = 65535;
    uint16_t polar_max = 0;
    for (int x = 0; x < W; ++x) {
      const uint16_t value = fx.get_pixel(x, y).r;
      polar_min = std::min(polar_min, value);
      polar_max = std::max(polar_max, value);
    }
    HS_EXPECT_GT(polar_min, 25000);
    HS_EXPECT_LT(polar_max, 35000);
  };
  auto expect_stripe = [&](int y) {
    for (int x = 0; x < W; ++x) {
      if (x & 1)
        HS_EXPECT_GT(fx.get_pixel(x, y).r, 58000);
      else
        HS_EXPECT_LT(fx.get_pixel(x, y).r, 2000);
    }
  };

  expect_reconstructed(1);
  expect_reconstructed(4);
  expect_stripe(8);
  expect_reconstructed(H - 2);
  expect_reconstructed(H - 5);
  expect_stripe(H - 9);

  // A zero pole_half_res composites every column, so row 4 keeps the stripe.
  style.pole_half_res = 0.0f;
  {
    Canvas c(fx);
    for (int x = 0; x < W; ++x)
      c(x, 4) = Pixel((x & 1) ? 60000 : 0, 0, 0);
  }
  fx.advance_display();
  {
    Canvas c(fx);
    (void)pipe.begin_frame(c, 1.0f);
  }
  fx.advance_display();
  expect_stripe(4);
}

/**
 * @brief Verifies flush() honors the segment clip like every other rasterizer.
 * @details On segmented hardware each board owns a Y-band, and a feedback flush
 *          that iterated the full canvas would composite the whole sphere into
 *          every board's buffer (wrong output + wasted work). Here the prev frame
 *          is bright at every row but the clip restricts rendering to a sub-band;
 *          rows outside the margin-expanded render band must stay untouched.
 */
inline void test_feedback_flush_respects_clip() {
  constexpr int W = 32, H = 16; // both divisible by Smoke's downsample (4)
  hs_test::StubEffect fx(W, H);
  ::Feedback::Style style = ::Feedback::Style::Smoke(); // identity warp
  Pipeline<W, H, Filter::Pixel::Feedback<W, H>> pipe{
      Filter::Pixel::Feedback<W, H>(style)};

  // Frame 1: bright across the WHOLE canvas becomes the "previous" frame.
  {
    Canvas c(fx);
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x)
        c(x, y) = Pixel(40000, 40000, 40000);
  }
  fx.advance_display();

  // Restrict this board to the y-band [8,12). With the default margin of 1 the
  // render band is [7,13); rows outside it must not be written.
  fx.set_clip(8, 12, 0, W);

  // Frame 2: empty buffer; clipped flush warps+blends only within the band.
  {
    Canvas c(fx);
    (void)pipe.begin_frame(c, 1.0f);
  }
  fx.advance_display();

  // Inside the band: carried over and faded.
  HS_EXPECT_FALSE(is_black(fx.get_pixel(W / 2, 9)));
  HS_EXPECT_FALSE(is_black(fx.get_pixel(W / 2, 7)));
  HS_EXPECT_FALSE(is_black(fx.get_pixel(W / 2, 12)));
  HS_EXPECT_TRUE(is_black(fx.get_pixel(W / 2, 6)));
  HS_EXPECT_TRUE(is_black(fx.get_pixel(W / 2, 13)));
  // Outside the render band: untouched despite the prev frame being lit there.
  HS_EXPECT_TRUE(is_black(fx.get_pixel(W / 2, 0)));
  HS_EXPECT_TRUE(is_black(fx.get_pixel(W / 2, 5)));
  HS_EXPECT_TRUE(is_black(fx.get_pixel(W / 2, 15)));
}

/**
 * @brief Verifies the Pixel::Feedback::flush warp path under a NON-identity warp,
 *        end to end through the coarse-grid + bilinear-upsample pipeline.
 * @details The identity-warp flush tests above have no bound NoiseParams,
 *          so space_fn collapses to the identity map: the coarse warp
 *          field is all-zero and the bilinear upsample is exercised only on a
 *          degenerate (constant-zero) field — exactly the most bug-prone part left
 *          uncovered. This drives melt_warp instead, a deterministic-without-noise
 *          transform that slerps every sample direction toward the north pole by
 *          drip = speed * 0.04, so the previous frame "drips" south by a known
 *          amount. The spherical ring lattice has W/downsample samples on the
 *          equator ring and fewer toward the poles, bilinearly interpolated per
 *          pixel, so a correct displacement here
 *          proves the whole path (space_fn -> coarse deltas -> bilerp -> sample).
 *
 *          The displacement is predicted with the SAME production helpers the
 *          flush uses as its oracle: for output row y the warp samples source row
 *          by(y) = phi_to_y(Spherical(slerp(pixel_to_vector(0,y), +Y, drip)).phi),
 *          which is < y (north of y) and independent of x (slerp toward a pole
 *          preserves longitude). A bright source band at row R therefore re-appears
 *          at the output row y* where by(y*) == R, strictly south of R. The test
 *          asserts the output band's brightest row matches y* (within the 1px the
 *          coarse-grid bilerp can shift it) and that the band's original row went
 *          dark — i.e. the content actually moved, ruling out an identity warp.
 */
inline void test_feedback_flush_melt_warp_displaces_south() {
  constexpr int W = 64, H = 64; // both divisible by the downsample (4)
  hs_test::StubEffect fx(W, H);

  // The default zero hue_shift isolates the SPATIAL warp under test from any hue
  // rotation; noise stays nullptr so melt_warp is fully deterministic (the noise
  // wobble branch is gated on a bound NoiseParams). speed=6 -> drip=0.24 gives a
  // clearly multi-pixel southward shift.
  ::Feedback::Style style{};
  style.space_fn = &::Feedback::melt_warp;
  style.noise = nullptr;
  style.speed = 6.0f;
  style.fade = 0.9f;
  style.downsample = 4;
  const float drip = style.speed * 0.04f;

  Pipeline<W, H, Filter::Pixel::Feedback<W, H>> pipe{
      Filter::Pixel::Feedback<W, H>(style)};

  // Frame 1: a full-width bright band (3 rows thick) becomes the "previous" frame.
  constexpr int R = 40; // band center, southern hemisphere
  {
    Canvas c(fx);
    for (int y = R - 1; y <= R + 1; ++y)
      for (int x = 0; x < W; ++x)
        c(x, y) = Pixel(40000, 40000, 40000);
  }
  fx.advance_display();

  // Frame 2: empty buffer; the melt flush warps + blends the band into it.
  {
    Canvas c(fx);
    (void)pipe.begin_frame(c, 1.0f);
  }
  fx.advance_display();

  // Oracle: the output row that samples source row R. by(y) = warped source row
  // for output row y, computed with the production helpers (x-independent, so use
  // column 0). Pick the y whose source row is closest to the band center.
  const math::Vector NORTH(0.0f, 1.0f, 0.0f);
  auto by = [&](int y) {
    math::Vector v = math::pixel_to_vector<W, H>(0, y);
    return math::phi_to_y<H>(math::Spherical(math::slerp(v, NORTH, drip)).phi);
  };
  int oracle_y = R;
  float best = static_cast<float>(H);
  for (int y = 0; y < H; ++y) {
    float d = std::abs(by(y) - static_cast<float>(R));
    if (d < best) {
      best = d;
      oracle_y = y;
    }
  }

  // Find the output band by per-row brightness (band is full width, so summing
  // each row is robust to the bilerp's small per-pixel spread).
  auto row_sum = [&](int y) {
    uint64_t s = 0;
    for (int x = 0; x < W; ++x)
      s += fx.get_pixel(x, y).r;
    return s;
  };
  int peak_y = 0;
  uint64_t peak = 0;
  for (int y = 0; y < H; ++y) {
    uint64_t s = row_sum(y);
    if (s > peak) {
      peak = s;
      peak_y = y;
    }
  }

  // The band drifted south to the predicted row (within the coarse bilerp's 1px),
  // strictly past where an identity warp would have left it.
  HS_EXPECT_GT(peak, (uint64_t)0); // the warped band actually rendered
  HS_EXPECT_TRUE(std::abs(peak_y - oracle_y) <= 1);
  HS_EXPECT_GT(oracle_y, R);   // melt drifts south (y increases)
  HS_EXPECT_GT(peak_y, R + 2); // a real, multi-pixel displacement
  // The band's original location is now dark: its source row is north of the band,
  // so nothing bright maps back onto row R — confirming the content moved, not an
  // identity passthrough that would have left the band in place.
  HS_EXPECT_LT(row_sum(R), peak / 4);
}

inline math::Vector north_cap_rotation_warp(const math::Vector &v,
                                            const ::Feedback::Style &) {
  constexpr float ANGLE = 0.1f;
  const float c = std::cos(ANGLE);
  const float s = std::sin(ANGLE);
  return math::Vector(c * v.x - s * v.y, s * v.x + c * v.y, v.z);
}

inline float animated_cap_rotation_angle = 0.0f;

inline math::Vector animated_cap_rotation_warp(const math::Vector &v,
                                               const ::Feedback::Style &) {
  const float c = std::cos(animated_cap_rotation_angle);
  const float s = std::sin(animated_cap_rotation_angle);
  return math::Vector(c * v.x - s * v.y, s * v.x + c * v.y, v.z);
}

struct FloatRowAccumulator {
  float sum = 0.0f;

  void add(float value) { sum += value; }
  void remove(float value) { sum -= value; }
  float average(int width) const { return sum / width; }
};

template <int W, int H>
inline std::array<float, W>
expected_feedback_source_row(int y, int downsample,
                             const ::Feedback::Style &style) {
  using Layout = hs::SphericalFieldLayout<W, H>;
  const Layout layout(downsample, downsample, downsample, W / downsample);
  std::array<float, W> source{};
  std::array<float, W> reconstructed{};
  // Cap rows reconstruct each pixel's own target, so the reference is the
  // exact warp target row before the infill filter.
  for (int x = 0; x < W; ++x)
    source[x] =
        layout.project(style.space_fn(math::pixel_to_vector<W, H>(x, y), style))
            .y;
  // The render only longitude-filters the dense infill rows.
  const bool north_infill =
      y < downsample && (!Layout::HAS_NORTH_POLE || y > 0);
  const bool south_infill =
      y >= H - downsample && (!Layout::HAS_SOUTH_POLE || y < H - 1);
  if (!north_infill && !south_infill)
    return source;
  layout.template reconstruct_longitude_row<FloatRowAccumulator>(
      source.data(), y, [&](int x, float value) { reconstructed[x] = value; });
  return reconstructed;
}

/**
 * @brief Verifies north-cap rows use warp values evaluated at their latitude.
 */
inline void test_feedback_north_cap_uses_exact_control_rows() {
  constexpr int W = 64, H = 64;
  constexpr int DOWNSAMPLE = 4;
  constexpr float ROW_SCALE = 1000.0f;
  hs_test::StubEffect fx(W, H);

  ::Feedback::Style style{};
  style.space_fn = &north_cap_rotation_warp;
  style.fade = 1.0f;
  style.downsample = DOWNSAMPLE;

  Pipeline<W, H, Filter::Pixel::Feedback<W, H>> pipe{
      Filter::Pixel::Feedback<W, H>(style)};

  {
    Canvas c(fx);
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x)
        c(x, y) = Pixel(static_cast<uint16_t>(y * ROW_SCALE), 0, 0);
  }
  fx.advance_display();

  {
    Canvas c(fx);
    (void)pipe.begin_frame(c, 1.0f);
  }
  fx.advance_display();

  float first_offset = 0.0f;
  bool offsets_vary = false;
  for (int y = 1; y < DOWNSAMPLE; ++y) {
    const auto expected =
        expected_feedback_source_row<W, H>(y, DOWNSAMPLE, style);
    const float sampled_y = fx.get_pixel(0, y).r / ROW_SCALE;
    HS_EXPECT_NEAR(sampled_y, expected[0], 0.02f);
    const float offset = sampled_y - static_cast<float>(y);
    if (y == 1)
      first_offset = offset;
    else if (std::abs(offset - first_offset) > 1e-3f)
      offsets_vary = true;
  }
  // Independent of the reference row: the claim is that each cap row evaluates
  // the warp at its OWN latitude, so a collapse onto one shared control row —
  // which the reference would follow if it were derived the same way — shows up
  // as an identical displacement on every cap row.
  HS_EXPECT_TRUE(offsets_vary);
}

/**
 * @brief Verifies animated polar controls use the compositor's longitudes.
 */
inline void test_feedback_animated_cap_controls_match_compositor_lattice() {
  constexpr int W = 64, H = 34;
  constexpr int DOWNSAMPLE = 4;
  constexpr float ROW_SCALE = 1000.0f;
  hs_test::StubEffect fx(W, H);

  ::Feedback::Style style{};
  style.space_fn = &animated_cap_rotation_warp;
  style.fade = 1.0f;
  style.downsample = DOWNSAMPLE;

  Pipeline<W, H, Filter::Pixel::Feedback<W, H>> pipe{
      Filter::Pixel::Feedback<W, H>(style)};

  for (float angle : {0.15f, 0.35f, -0.25f}) {
    animated_cap_rotation_angle = angle;
    {
      Canvas c(fx);
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
          c(x, y) = Pixel(static_cast<uint16_t>(y * ROW_SCALE), 0, 0);
    }
    fx.advance_display();

    {
      Canvas c(fx);
      (void)pipe.begin_frame(c, 1.0f);
    }
    fx.advance_display();

    for (int y = 1; y < DOWNSAMPLE; ++y) {
      const auto expected =
          expected_feedback_source_row<W, H>(y, DOWNSAMPLE, style);
      for (int x = 0; x < W; x += DOWNSAMPLE) {
        const float sampled_y = fx.get_pixel(x, y).r / ROW_SCALE;
        HS_EXPECT_NEAR(sampled_y, expected[x], 0.02f);
      }
    }
  }
}

/**
 * @brief Verifies projected pole displacement retains its longitude branch.
 */
inline void test_feedback_poles_resolve_one_source_longitude() {
  constexpr int W = 64, H = 64;
  hs_test::StubEffect fx(W, H);
  ::Feedback::Style style{};
  style.space_fn = &north_cap_rotation_warp;
  style.fade = 1.0f;
  style.downsample = 4;
  Pipeline<W, H, Filter::Pixel::Feedback<W, H>> pipe{
      Filter::Pixel::Feedback<W, H>(style)};

  {
    Canvas c(fx);
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x)
        c(x, y) = Pixel(static_cast<uint16_t>(x * 900 + 1000), 0, 0);
  }
  fx.advance_display();
  {
    Canvas c(fx);
    (void)pipe.begin_frame(c, 1.0f);
  }
  fx.advance_display();

  for (int y : {0, H - 1}) {
    const uint16_t expected = fx.get_pixel(0, y).r;
    HS_EXPECT_NEAR(expected, y == 0 ? 32 * 900 + 1000 : 1000, 300.0f);
    for (int x = 4; x < W; x += 4)
      HS_EXPECT_NEAR(fx.get_pixel(x, y).r, expected, 2.0f);
  }
}

inline math::Vector midpoint_test_warp(const math::Vector &v,
                                       const ::Feedback::Style &) {
  const math::Spherical S(v);
  return math::Vector(math::Spherical(S.theta + 0.1f * sinf(S.theta),
                                      S.phi + 0.02f * sinf(S.theta)));
}

inline void test_feedback_half_res_warp_uses_pair_midpoint() {
  constexpr int W = 64, H = 34;
  std::array<Pixel, W> full{};
  for (bool half_res : {false, true}) {
    hs_test::StubEffect fx(W, H);
    ::Feedback::Style style{};
    style.space_fn = &midpoint_test_warp;
    style.fade = 1.0f;
    style.downsample = 4;
    style.pole_half_res = half_res ? 2.0f : 0.0f;
    Pipeline<W, H, Filter::Pixel::Feedback<W, H>> pipe{
        Filter::Pixel::Feedback<W, H>(style)};
    {
      Canvas c(fx);
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
          c(x, y) = Pixel(static_cast<uint16_t>(x * 900),
                          static_cast<uint16_t>(y * 1500), 0);
    }
    fx.advance_display();
    {
      Canvas c(fx);
      (void)pipe.begin_frame(c, 1.0f);
    }
    fx.advance_display();
    if (!half_res) {
      for (int x = 0; x < W; ++x)
        full[x] = fx.get_pixel(x, 8);
    } else {
      auto midpoint = [&](int x, auto channel) {
        return (full[x].*channel + full[x + 1].*channel) * 0.5f;
      };
      for (int x = 8; x < 24; x += 2) {
        for (auto channel : {&Pixel::r, &Pixel::g}) {
          const float OWN = midpoint(x, channel);
          const float LEFT = midpoint(x - 2, channel);
          const float RIGHT = midpoint(x + 2, channel);
          HS_EXPECT_NEAR(fx.get_pixel(x, 8).*channel,
                         0.75f * OWN + 0.25f * LEFT, 4.0f);
          HS_EXPECT_NEAR(fx.get_pixel(x + 1, 8).*channel,
                         0.75f * OWN + 0.25f * RIGHT, 4.0f);
        }
      }
    }
  }
}

inline math::Vector metric_row_test_warp(const math::Vector &v,
                                         const ::Feedback::Style &) {
  const math::Spherical s(v);
  const float delta = 0.08f * std::sin(s.phi) * std::sin(18.0f * s.phi);
  return math::Vector(math::Spherical(s.theta, s.phi + delta));
}

/**
 * @brief Verifies the feedback row lattice follows its spherical metric.
 */
inline void test_feedback_spherical_ring_control_rows() {
  constexpr int W = 64, H = 64;
  constexpr int DOWNSAMPLE = 4;
  constexpr float ROW_SCALE = 1000.0f;
  hs_test::StubEffect fx(W, H);

  ::Feedback::Style style{};
  style.space_fn = &metric_row_test_warp;
  style.fade = 1.0f;
  style.downsample = DOWNSAMPLE;

  Pipeline<W, H, Filter::Pixel::Feedback<W, H>> pipe{
      Filter::Pixel::Feedback<W, H>(style)};

  {
    Canvas c(fx);
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x)
        c(x, y) = Pixel(static_cast<uint16_t>(y * ROW_SCALE), 0, 0);
  }
  fx.advance_display();

  {
    Canvas c(fx);
    (void)pipe.begin_frame(c, 1.0f);
  }
  fx.advance_display();

  constexpr hs::SphericalFieldLayout<W, H> layout(DOWNSAMPLE, DOWNSAMPLE,
                                                  DOWNSAMPLE, W / DOWNSAMPLE);
  for (int ring_index = 0; ring_index < layout.ring_count(); ++ring_index) {
    const int y = layout.ring(ring_index).y;
    const math::Vector warped =
        metric_row_test_warp(math::pixel_to_vector<W, H>(0, y), style);
    const float expected_y = math::phi_to_y<H>(math::Spherical(warped).phi);
    const float sampled_y = fx.get_pixel(0, y).r / ROW_SCALE;
    HS_EXPECT_NEAR(sampled_y, expected_y, 0.02f);
  }

  // Warp fields and projected origins span all cells; cap offsets span polar cells.
  constexpr size_t BYTES_PER_ROW =
      (W / DOWNSAMPLE) *
      (2 * sizeof(int16_t) + sizeof(decltype(layout)::Coordinates));
  size_t cap_rows = 0;
  for (int i = 0; i < layout.ring_count(); ++i)
    if (decltype(layout)::latitude_sine(layout.ring(i).y) <
        Filter::Pixel::Feedback<W, H>::POLAR_TARGET_SINE)
      ++cap_rows;
  if constexpr (HS_RUNTIME_DISPLAY_GEOMETRY)
    cap_rows = static_cast<size_t>(layout.ring_count());
  constexpr size_t STORAGE = Filter::Pixel::Feedback<W, H>::STORAGE_BYTES;
  const size_t EXPECTED =
      BYTES_PER_ROW * static_cast<size_t>(layout.ring_count()) +
      cap_rows * (W / DOWNSAMPLE) * 2 * sizeof(int16_t);
  HS_EXPECT_EQ(STORAGE, EXPECTED);
  HS_EXPECT_LE(layout.sample_count(), (W / DOWNSAMPLE) * layout.ring_count());
}

/**
 * @brief Verifies the compact ring field has directionally balanced error.
 * @details The last bound anchors the compact field against metric_approximate,
 * a baseline stepping sin(phi)-scaled rows instead of reading the ring table.
 * The compact field measures 1.370x that baseline's mean polar error, so the
 * bound below leaves ~4% for libm drift and nothing for a regression.
 */
inline void test_feedback_spherical_field_angular_error() {
  constexpr int W = 288, H = 144;
  constexpr int DOWNSAMPLE = 4;
  constexpr hs::SphericalFieldLayout<W, H> layout(DOWNSAMPLE, DOWNSAMPLE,
                                                  DOWNSAMPLE, W / DOWNSAMPLE);
  struct Offset {
    float x;
    float y;
  };
  struct TangentError {
    float east;
    float down;
  };
  // Log map of `measured` at `reference`, resolved on the equirect east/down
  // basis at (x, y): a great-circle error split into its azimuthal and
  // meridional parts.
  auto tangent_error = [](const math::Vector &reference,
                          const math::Vector &measured, int x, int y) {
    const float theta = (2.0f * math::PI_F * x) / W;
    const float phi = math::y_to_phi<H>(y);
    const math::Vector east(-std::sin(theta), 0.0f, std::cos(theta));
    const math::Vector down(std::cos(phi) * std::cos(theta), -std::sin(phi),
                            std::cos(phi) * std::sin(theta));
    const float c = hs::clamp(math::dot(reference, measured), -1.0f, 1.0f);
    const math::Vector u = measured - reference * c;
    const float len_sq = math::dot(u, u);
    if (len_sq < math::EPS_NORMALIZE_SQ)
      return TangentError{0.0f, 0.0f};
    const math::Vector delta =
        u * (math::fast_acos(c) * math::fast_rsqrt(len_sq));
    return TangentError{math::dot(delta, east), math::dot(delta, down)};
  };
  Animation::NoiseParams noise;
  ::Feedback::Style style{};
  style.noise = &noise;
  style.amplitude = 3.99f;
  style.frequency = 0.37036f;
  style.speed = 2.625f;
  style.scale = 50.0f;
  std::vector<Offset> controls(layout.sample_count());
  std::vector<Offset> expanded(layout.ring_count() * (W / DOWNSAMPLE));
  hs::SphericalField<Offset, W, H> field(controls.data(), layout);
  double east_sq = 0.0;
  double down_sq = 0.0;
  double polar_error = 0.0;
  double equator_error = 0.0;
  double metric_polar_error = 0.0;
  double polar_weight = 0.0;
  double equator_weight = 0.0;

  for (int frame = 0; frame < 8; ++frame) {
    noise.time = 7.0f + frame * 3.25f;
    style.sync_noise();
    field.populate(
        0, layout.ring_count() - 1,
        [&](const math::Vector &position, const auto &point) {
          const math::Spherical warped(style.space_fn(position, style));
          float dx = warped.theta * W / (2.0f * math::PI_F) - point.x;
          if (dx > W * 0.5f)
            dx -= W;
          else if (dx < -W * 0.5f)
            dx += W;
          return Offset{dx, math::phi_to_y<H>(warped.phi) - point.y};
        });

    auto circular_lerp = [](const Offset &a, Offset b, float t) {
      if (b.x - a.x > W * 0.5f)
        b.x -= W;
      else if (b.x - a.x < -W * 0.5f)
        b.x += W;
      return Offset{hs::lerp(a.x, b.x, t), hs::lerp(a.y, b.y, t)};
    };
    for (int ring_index = 0; ring_index < layout.ring_count(); ++ring_index) {
      const auto ring = layout.ring(ring_index);
      for (int x = 0; x < W; x += DOWNSAMPLE) {
        const auto longitude = layout.longitude_bounded(ring, x);
        expanded[ring_index * (W / DOWNSAMPLE) + x / DOWNSAMPLE] =
            circular_lerp(controls[longitude.left], controls[longitude.right],
                          longitude.mix);
      }
    }

    auto approximate = [&](int x, int y) {
      const auto latitude = layout.row(y);
      const int lower_ring = layout.ring_index_at_or_before(y);
      const int upper_ring = layout.ring_index_at_or_after(y);
      const int coarse_x = x / DOWNSAMPLE;
      const int next_x = coarse_x + 1 < W / DOWNSAMPLE ? coarse_x + 1 : 0;
      const float fx = static_cast<float>(x % DOWNSAMPLE) / DOWNSAMPLE;
      const int row0 = lower_ring * (W / DOWNSAMPLE);
      const int row1 = upper_ring * (W / DOWNSAMPLE);
      const Offset lower =
          circular_lerp(expanded[row0 + coarse_x], expanded[row0 + next_x], fx);
      const Offset upper =
          circular_lerp(expanded[row1 + coarse_x], expanded[row1 + next_x], fx);
      const Offset offset = circular_lerp(lower, upper, latitude.mix);
      float source_x = std::fmod(x + offset.x, static_cast<float>(W));
      if (source_x < 0.0f)
        source_x += W;
      return math::pixel_to_vector<W, H>(source_x, y + offset.y);
    };

    auto exact_offset = [&](int x, int y) {
      const math::Spherical warped(
          style.space_fn(math::pixel_to_vector<W, H>(x, y), style));
      float dx = warped.theta * W / (2.0f * math::PI_F) - x;
      if (dx > W * 0.5f)
        dx -= W;
      else if (dx < -W * 0.5f)
        dx += W;
      return Offset{dx, math::phi_to_y<H>(warped.phi) - y};
    };
    auto metric_approximate = [&](int x, int y) {
      int y0 = 0;
      int y1 = 0;
      while (y1 < y) {
        y0 = y1;
        int step = static_cast<int>(
            DOWNSAMPLE * std::sin(math::y_to_phi<H>(y1)) + 0.5f);
        step = hs::clamp(step, 1, DOWNSAMPLE);
        y1 = std::min(y1 + step, H - 1);
      }
      const int x0 = (x / DOWNSAMPLE) * DOWNSAMPLE;
      const int x1 = (x0 + DOWNSAMPLE) % W;
      const float fx = static_cast<float>(x - x0) / DOWNSAMPLE;
      const float fy = y1 == y0 ? 0.0f : static_cast<float>(y - y0) / (y1 - y0);
      const Offset lower =
          circular_lerp(exact_offset(x0, y0), exact_offset(x1, y0), fx);
      const Offset upper =
          circular_lerp(exact_offset(x0, y1), exact_offset(x1, y1), fx);
      const Offset offset = circular_lerp(lower, upper, fy);
      float source_x = std::fmod(x + offset.x, static_cast<float>(W));
      if (source_x < 0.0f)
        source_x += W;
      return math::pixel_to_vector<W, H>(source_x, y + offset.y);
    };

    for (int y = 5; y < 88; y += 2)
      for (int x = 3; x < W; x += 8) {
        const math::Vector source = math::pixel_to_vector<W, H>(x, y);
        const math::Vector exact = style.space_fn(source, style);
        const math::Vector got = approximate(x, y);
        const TangentError error = tangent_error(exact, got, x, y);
        const double area_weight = std::sin(math::y_to_phi<H>(y));
        east_sq += area_weight * error.east * error.east;
        down_sq += area_weight * error.down * error.down;
        const float angular_error = math::angle_between(exact, got);
        if (y < 32) {
          polar_error += area_weight * angular_error;
          metric_polar_error +=
              area_weight *
              math::angle_between(exact, metric_approximate(x, y));
          polar_weight += area_weight;
        } else if (y >= 56) {
          equator_error += area_weight * angular_error;
          equator_weight += area_weight;
        }
      }
  }

  const double direction_ratio = east_sq / down_sq;
  polar_error /= polar_weight;
  equator_error /= equator_weight;
  metric_polar_error /= polar_weight;
  HS_EXPECT_GT(direction_ratio, 0.7);
  HS_EXPECT_LT(direction_ratio, 1.3);
  HS_EXPECT_LT(polar_error, equator_error * 1.35);
  HS_EXPECT_LT(equator_error, polar_error * 1.35);
  HS_EXPECT_LT(polar_error, metric_polar_error * 1.42);
}

/** @brief Displaces longitude alone, alternating a near-half-turn against a
 * small opposing step so adjacent controls straddle the seam correction. */
inline math::Vector opposed_seam_warp(const math::Vector &v,
                                      const ::Feedback::Style &) {
  const math::Spherical s(v);
  const int band =
      static_cast<int>(std::floor(s.theta * (256.0f / (2.0f * math::PI_F))));
  const float shift = (band & 1) ? (0.999f * math::PI_F) : (-0.1f * math::PI_F);
  return math::Vector(math::Spherical(s.theta + shift, s.phi));
}

/**
 * @brief Verifies a longitude-only warp samples its own latitude row.
 * @details populate_warp_field's seam correction can lift an interpolated
 * offset past the [-W, 2W) domain sample_bilinear contracts for; the resulting
 * out-of-range column index walks the flat framebuffer into a neighbouring row.
 * The band count and opposing step are tuned to place a control pair at that
 * extreme -- the excursion needs a near-half-turn against a small reverse step.
 */
inline void test_feedback_seam_warp_keeps_its_latitude_row() {
  constexpr int W = 288, H = 32;
  constexpr int DOWNSAMPLE = 4;
  constexpr float ROW_SCALE = 1000.0f;
  hs_test::StubEffect fx(W, H);

  ::Feedback::Style style{};
  style.space_fn = &opposed_seam_warp;
  style.fade = 1.0f;
  style.downsample = DOWNSAMPLE;

  Pipeline<W, H, Filter::Pixel::Feedback<W, H>> pipe{
      Filter::Pixel::Feedback<W, H>(style)};

  {
    Canvas c(fx);
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x)
        c(x, y) = Pixel(static_cast<uint16_t>(y * ROW_SCALE), 0, 0);
  }
  fx.advance_display();

  fx.set_margin(0);
  fx.set_clip(0, H, 0, W);
  {
    Canvas c(fx);
    (void)pipe.begin_frame(c, 1.0f);
  }
  fx.advance_display();

  // Polar rows interpolate cap-plane offsets rather than column offsets, and
  // this warp's alternating half-turns put their targets across the pole.
  for (int y = 1; y < H - 1; ++y) {
    if (hs::SphericalFieldLayout<W, H>::latitude_sine(y) <
        Filter::Pixel::Feedback<W, H>::POLAR_TARGET_SINE)
      continue;
    for (int x = 0; x < W; ++x) {
      const float sampled_row = fx.get_pixel(x, y).r / ROW_SCALE;
      HS_EXPECT_NEAR(sampled_row, static_cast<float>(y), 0.25f);
    }
  }
}

/**
 * @brief Verifies cached top-band clips share fully populated cap controls.
 */
inline void test_feedback_cached_north_cap_clips_share_control_rows() {
  constexpr int W = 64, H = 64;
  constexpr int DOWNSAMPLE = 4;
  constexpr float ROW_SCALE = 1000.0f;
  hs_test::StubEffect fx(W, H);
  Animation::NoiseParams noise;
  ScratchScope persistent_scope(persistent_arena);

  ::Feedback::Style style{};
  style.space_fn = &::Feedback::noise_warp;
  style.noise = &noise;
  style.amplitude = 10.0f;
  style.frequency = 0.2f;
  style.speed = 0.0f;
  style.scale = 8.0f;
  style.fade = 1.0f;
  style.downsample = DOWNSAMPLE;
  style.sync_noise();

  Pipeline<W, H, Filter::Pixel::Feedback<W, H>> pipe{
      Filter::Pixel::Feedback<W, H>(style)};
  pipe.get<Filter::Pixel::Feedback<W, H>>().init_storage(persistent_arena);
  auto seed_rows = [&]() {
    {
      Canvas c(fx);
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
          c(x, y) = Pixel(static_cast<uint16_t>(y * ROW_SCALE), 0, 0);
    }
    fx.advance_display();
  };
  // The reference is the exact target row; the lattice interpolates it.
  auto expect_rows = [&](int begin, int end) {
    for (int y = begin; y < end; ++y) {
      const auto expected =
          expected_feedback_source_row<W, H>(y, DOWNSAMPLE, style);
      const float sampled_y = fx.get_pixel(0, y).r / ROW_SCALE;
      HS_EXPECT_NEAR(sampled_y, expected[0], 0.03f);
    }
  };

  fx.set_margin(0);
  seed_rows();
  fx.set_clip(0, 2, 0, W);
  {
    Canvas c(fx);
    (void)pipe.begin_frame(c, 1.0f);
  }
  fx.advance_display();
  expect_rows(0, 2);

  seed_rows();
  fx.set_clip(2, 4, 0, W);
  {
    Canvas c(fx);
    (void)pipe.begin_frame(c, 1.0f);
  }
  fx.advance_display();
  expect_rows(2, 4);
}

/**
 * @brief Verifies the polar rows land every pixel on its warp target.
 * @details Encodes each pixel's own direction into the previous frame, flushes
 * once through the identity colour path, and decodes where each output pixel
 * sampled from. A strong static twist near the poles moves targets far in
 * longitude between lattice rings, where interpolating equirect offsets missed
 * by several degrees; the 3D target reconstruction stays within a fraction of
 * the 1.26 degree row pitch.
 */
inline void test_feedback_polar_rows_hit_their_targets() {
  constexpr int W = 288, H = 144;
  hs_test::StubEffect fx(W, H);
  ScratchScope persistent_scope(persistent_arena);
  Animation::NoiseParams noise;
  noise.noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
  noise.set_seed(12345);
  ::Feedback::Style style = ::Feedback::Style::LooseWormhole();
  style.noise = &noise;
  style.sync_noise();
  style.hue_shift = 0.0f;
  style.fade = 1.0f;
  style.pole_half_res = 0.0f;
  Pipeline<W, H, Filter::Pixel::Feedback<W, H>> pipe{
      Filter::Pixel::Feedback<W, H>(style)};
  pipe.get<Filter::Pixel::Feedback<W, H>>().init_storage(persistent_arena);

  auto encode = [](float c) {
    return static_cast<uint16_t>((c + 1.0f) * 0.5f * 65535.0f + 0.5f);
  };
  auto decode = [](uint16_t c) { return c / 65535.0f * 2.0f - 1.0f; };
  {
    Canvas c(fx);
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x) {
        const math::Vector v = math::pixel_to_vector<W, H>(x, y);
        c(x, y) = Pixel(encode(v.x), encode(v.y), encode(v.z));
      }
  }
  fx.advance_display();
  {
    Canvas c(fx);
    (void)pipe.begin_frame(c, 1.0f);
  }
  fx.advance_display();

  float worst = 0.0f;
  for (int y = 1; y < H - 1; ++y) {
    if (hs::SphericalFieldLayout<W, H>::latitude_sine(y) >= 0.45f)
      continue;
    for (int x = 0; x < W; ++x) {
      const Pixel p = fx.get_pixel(x, y);
      const math::Vector got =
          math::Vector(decode(p.r), decode(p.g), decode(p.b)).normalized();
      const math::Vector want =
          style.space_fn(math::pixel_to_vector<W, H>(x, y), style).normalized();
      worst = std::max(worst,
                       std::acos(hs::clamp(math::dot(got, want), -1.0f, 1.0f)));
    }
  }
  HS_EXPECT_LT(worst * 180.0f / math::PI_F, 0.5f);
}

/**
 * @brief Warp-cache parity: an init_storage'd Feedback filter must render
 *        exactly what an uncached one does, frame for frame.
 * @details Drives a cached and an uncached pipeline through identical frames:
 *          a static style (later frames hit the cache), a key-field mutation
 *          (amplitude), a generator seed change that no Style scalar mirrors,
 *          advancing noise time under nonzero speed (key changes every frame),
 *          and a mid-run init_storage() re-allocation (the effect's
 *          post-compaction path). Every frame must match the uncached reference
 *          pixel-exactly — reuse serves the same int16 deltas the populate pass
 *          would have written.
 */
inline void test_feedback_warp_cache_matches_uncached() {
  constexpr int W = 64, H = 64; // both divisible by the downsample (4)
  constexpr int FRAMES = 8;
  // Frames 3 and 4 share a band, so frame 4's seed change is the only key
  // difference between them and lands on a populated cache entry.
  constexpr int CLIP_BEGIN[FRAMES] = {20, 36, 24, 40, 40, 20, 36, 24};

  // Only one Effect may be alive at a time (shared static buffers), so the
  // two pipelines run sequentially over recorded frames.
  auto run = [&](bool cached) {
    std::vector<std::vector<Pixel>> frames;
    hs_test::StubEffect fx(W, H);
    ScratchScope persistent_scope(persistent_arena);
    Animation::NoiseParams np;
    ::Feedback::Style s{};
    s.noise = &np;
    s.amplitude = 3.0f;
    s.frequency = 0.3f;
    s.speed = 0.0f;
    s.scale = 5.0f;
    s.fade = 0.9f;

    Pipeline<W, H, Filter::Pixel::Feedback<W, H>> pipe{
        Filter::Pixel::Feedback<W, H>(s)};
    const size_t CACHE_OFFSET = persistent_arena.get_offset();
    if (cached)
      pipe.template get<Filter::Pixel::Feedback<W, H>>().init_storage(
          persistent_arena);

    {
      Canvas c(fx);
      for (int y = 20; y <= 44; ++y)
        for (int x = 8; x < W; x += 3)
          c(x, y) = Pixel(40000, 20000 + 400 * y, 60000);
    }
    fx.advance_display();

    for (int frame = 0; frame < FRAMES; ++frame) {
      fx.set_margin(0);
      fx.set_clip(CLIP_BEGIN[frame], CLIP_BEGIN[frame] + 8, 0, W);
      if (frame == 3)
        s.amplitude = 4.5f; // key change: repopulate
      if (frame == 4)
        np.set_seed(4242); // seed-only change: repopulate
      if (frame == 5)
        s.speed = 1.0f; // time-varying: miss every frame
      if (frame == 6 && cached) {
        persistent_arena.set_offset(CACHE_OFFSET);
        pipe.template get<Filter::Pixel::Feedback<W, H>>().init_storage(
            persistent_arena);
      }
      s.sync_noise();

      {
        Canvas c(fx);
        (void)pipe.begin_frame(c, 1.0f);
      }
      fx.advance_display();

      frames.emplace_back();
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
          frames.back().push_back(fx.get_pixel(x, y));

      np.time += np.speed; // mirrors Animation::Noise
    }
    return frames;
  };

  auto ref = run(false);
  auto got = run(true);
  for (int frame = 0; frame < FRAMES; ++frame) {
    HS_EXPECT_TRUE(
        std::any_of(ref[frame].begin(), ref[frame].end(),
                    [](const Pixel &p) { return p.r || p.g || p.b; }));
    int mismatches = 0;
    for (size_t i = 0; i < ref[frame].size(); ++i)
      if (!(ref[frame][i] == got[frame][i]))
        ++mismatches;
    HS_EXPECT_EQ(mismatches, 0);
  }
}

/**
 * @brief Space warp whose longitudinal displacement oscillates across the
 *        ±W/2 wrap: theta += PI + 0.3*sin(theta). Adjacent coarse warp-field
 *        columns land on opposite wrap branches while the true (unwrapped)
 *        field stays smooth.
 */
inline math::Vector antipodal_ripple_warp(const math::Vector &v,
                                          const ::Feedback::Style &) {
  math::Spherical s(v);
  return math::Vector(
      math::Spherical(s.theta + math::PI_F + 0.3f * std::sin(s.theta), s.phi));
}

/**
 * @brief Verifies the warp-field bilerp stays on one wrap branch when the
 *        coarse taps straddle the ±W/2 cut.
 * @details The antipodal ripple displaces every column by W/2 ± 1.5px, so the
 *          step-1 seam wrap flips sign between adjacent coarse columns
 *          (~ +16px vs ~ -15px at W=32) while the true field is smooth. The
 *          tap re-centering must unify the four taps onto one branch; a blend
 *          that sweeps across the cut instead samples near-zero displacements
 *          — pixels from the wrong side of the sphere, rendering as a
 *          longitudinal streak of foreign color. The previous frame encodes
 *          longitude seam-continuously (r=cos, g=sin), so each output pixel's
 *          actual sample source can be decoded and checked against the warp
 *          evaluated directly.
 */
inline void test_feedback_flush_straddled_taps_stay_on_branch() {
  constexpr int W = 32, H = 16; // both divisible by the downsample (4)
  constexpr float TWO_PI = 2.0f * math::PI_F;
  hs_test::StubEffect fx(W, H);

  ::Feedback::Style style{};
  style.space_fn = &antipodal_ripple_warp;
  style.fade = 1.0f;
  style.noise = nullptr;
  style.downsample = 4;

  Pipeline<W, H, Filter::Pixel::Feedback<W, H>> pipe{
      Filter::Pixel::Feedback<W, H>(style)};

  // Frame 1: seam-continuous longitude encoding becomes the "previous" frame.
  {
    Canvas c(fx);
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x) {
        float th = TWO_PI * x / W;
        c(x, y) = Pixel(
            static_cast<uint16_t>((std::cos(th) * 0.5f + 0.5f) * 65535.0f),
            static_cast<uint16_t>((std::sin(th) * 0.5f + 0.5f) * 65535.0f), 0);
      }
  }
  fx.advance_display();

  // Frame 2: empty buffer; flush pulls the warped prev frame into it.
  {
    Canvas c(fx);
    (void)pipe.begin_frame(c, 1.0f);
  }
  fx.advance_display();

  // Decode each pixel's sampled longitude and compare with the warp evaluated
  // directly. Interior rows dodge the pole rows' vertical clamp; the 2px
  // tolerance covers the coarse-grid bilerp and int16 quantization, while a
  // cross-branch sweep is off by ~W/4 at the straddle columns.
  float max_err = 0.0f;
  for (int y = 4; y < 12; ++y)
    for (int x = 0; x < W; ++x) {
      const Pixel &p = fx.get_pixel(x, y);
      float c = p.r / 65535.0f * 2.0f - 1.0f;
      float s = p.g / 65535.0f * 2.0f - 1.0f;
      float decoded_x = std::atan2(s, c) / TWO_PI * W;
      float th = TWO_PI * x / W;
      float expected_x = (th + math::PI_F + 0.3f * std::sin(th)) / TWO_PI * W;
      float d = std::fmod(decoded_x - expected_x, static_cast<float>(W));
      if (d > W * 0.5f)
        d -= W;
      if (d < -W * 0.5f)
        d += W;
      max_err = std::max(max_err, std::fabs(d));
    }
  HS_EXPECT_LT(max_err, 2.0f);
}
