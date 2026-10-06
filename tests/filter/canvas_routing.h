/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

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
