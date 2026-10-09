/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Segmented-mode rendering bound
// ============================================================================

/**
 * @brief The base Effect::needs_full_frame() defaults to false.
 * @details Effects whose filter pipeline folds any_crosses_segments true force
 *          full-frame rendering.
 */
inline void test_effect_needs_full_frame_default_false() {
  constexpr int W = 8, H = 8;
  hs_test::StubEffect fx(W, H);
  HS_EXPECT_FALSE(fx.needs_full_frame());
}

/**
 * @brief Screen::Trails parity fixture: a fresh trail buffer that plots one
 *        frame's seeds and flushes them onto an effect.
 * @details Seeds span every row and both x halves, plus one point that sweeps
 *          rows across frames. Trail brightness tracks remaining lifetime.
 */
struct ScreenTrailRig {
  static constexpr int W = 32, H = 16, MAXP = 512;
  static constexpr int LIFETIME = 4; // trail fade length (frames)
  using Trails = Filter::Screen::Trails<MAXP>;

  struct Seed {
    int x, y;
    Pixel c;
  };

  /** @brief Remaining-lifetime trail shade. */
  struct Shade {
    Color4 operator()(float, float, float t) const {
      const uint16_t v = static_cast<uint16_t>((1.0f - t) * 50000.0f);
      return Color4(Pixel(v, v, v), 1.0f);
    }
  };

  /** @brief Seeds plotted on frame @p f. */
  static std::array<Seed, 6> seeds(int f) {
    return {{
        {3, 1, Pixel(10000, 0, 0)},
        {12, 4, Pixel(0, 20000, 0)},
        {20, 7, Pixel(0, 0, 30000)},
        {7, 10, Pixel(15000, 15000, 0)},
        {25, 13, Pixel(0, 25000, 25000)},
        {17, (f * 3) % H, Pixel(40000, 40000, 40000)},
    }};
  }

  static inline uint8_t buf[MAXP * 32];
  Arena arena{buf, sizeof(buf)};
  Pipeline<W, H, Trails> pipe{Trails(LIFETIME)};

  ScreenTrailRig() { pipe.get<Trails>().init_storage(arena); }

  /** @brief Plots and flushes frame @p f's seeds onto @p fx. */
  void frame(hs_test::StubEffect &fx, int f) {
    Canvas c(fx);
    for (const Seed &s : seeds(f))
      pipe.plot(c, s.x, s.y, s.c, 0.0f, 1.0f);
    const Shade shade;
    pipe.flush(c, ScreenTrailFn(shade), 1.0f);
  }
};

/**
 * @brief Proves Screen::Trails has reach 0: under a FIXED clip a banded render
 *        matches the full-frame one byte-for-byte.
 */
inline void test_screen_trails_banded_matches_full() {
  constexpr int W = ScreenTrailRig::W, H = ScreenTrailRig::H;
  constexpr int K = 4; // frames driven
  constexpr int MID = H / 2;

  // Effect instances alias the same static double buffer (single-live guard),
  // so each run is scoped closed before the next.
  auto run = [&](int cy0, int cy1, Pixel out[H][W]) {
    ScreenTrailRig rig;
    hs_test::StubEffect fx(W, H);
    fx.set_clip(cy0, cy1, 0, W);
    for (int f = 0; f < K; ++f) {
      rig.frame(fx, f);
      fx.advance_display();
    }
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x)
        out[y][x] = fx.get_pixel(x, y);
  };

  static Pixel full[H][W], band_top[H][W], band_bot[H][W];
  run(0, H, full);       // single full-canvas instance
  run(0, MID, band_top); // worker A: top band
  run(MID, H, band_bot); // worker B: bottom band

  // Stitch each worker's DISPLAY band and require byte-identity with the full
  // instance: reach 0 => band clipping drops nothing.
  bool identical = true;
  int lit = 0;
  for (int y = 0; y < H; ++y)
    for (int x = 0; x < W; ++x) {
      const Pixel &want = full[y][x];
      const Pixel &got = (y < MID) ? band_top[y][x] : band_bot[y][x];
      if (!(got.r == want.r && got.g == want.g && got.b == want.b))
        identical = false;
      if (want.r | want.g | want.b)
        ++lit;
    }
  HS_EXPECT_TRUE(identical);
  HS_EXPECT_GT(lit, 0);
}

/**
 * @brief Proves Screen::Trails holds up under a clip that MOVES between frames,
 *        the segmented driver's per-frame arm-half alternation.
 * @details The clip alternates between [0, W/2) and [W/2, W) on successive
 *          frames; each frame's written half must equal that half of a
 *          full-canvas instance, so a point seeded while its half was clipped
 *          away must still re-emit once the clip swings back.
 */
inline void test_screen_trails_alternating_clip_matches_full() {
  constexpr int W = ScreenTrailRig::W, H = ScreenTrailRig::H;
  constexpr int K = 6; // frames driven
  constexpr int MID = W / 2;

  // flip alternates the x clip per frame; the reference run leaves the effect
  // at full canvas. Effect instances alias the same static double buffer, so
  // each run is scoped closed before the next.
  auto run = [&](bool flip, Pixel out[K][H][W]) {
    ScreenTrailRig rig;
    hs_test::StubEffect fx(W, H);
    for (int f = 0; f < K; ++f) {
      // Ahead of the Canvas: its stale-pixel clear honours the clip set here.
      if (flip)
        fx.set_clip(0, H, (f % 2) ? MID : 0, (f % 2) ? W : MID);
      rig.frame(fx, f);
      fx.advance_display();
      for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
          out[f][y][x] = fx.get_pixel(x, y);
    }
  };

  static Pixel full[K][H][W], flipped[K][H][W];
  run(false, full);
  run(true, flipped);

  // Per frame, the half the alternating instance owned must match the full
  // instance byte-for-byte.
  bool identical = true;
  int lit = 0;
  for (int f = 0; f < K; ++f) {
    const int x0 = (f % 2) ? MID : 0;
    const int x1 = (f % 2) ? W : MID;
    for (int y = 0; y < H; ++y)
      for (int x = x0; x < x1; ++x) {
        const Pixel &want = full[f][y][x];
        const Pixel &got = flipped[f][y][x];
        if (!(got.r == want.r && got.g == want.g && got.b == want.b))
          identical = false;
        if (want.r | want.g | want.b)
          ++lit;
      }
  }
  HS_EXPECT_TRUE(identical);
  HS_EXPECT_GT(lit, 0);
}

/**
 * @brief Proves a band-clipped feedback effect DIVERGES from the full-frame
 *        render — i.e. why crosses_segments forces full-frame for Pixel::Feedback.
 * @details Feedback reads cv.prev at unbounded warp offsets. A melt warp drips
 *          content south across the segment boundary from rows a worker clipped
 *          to the bottom band never seeded.
 */
inline void test_feedback_banded_diverges_from_full() {
  constexpr int W = 64, H = 64; // divisible by the downsample (4)
  using FB = Filter::Pixel::Feedback<W, H>;
  constexpr int K = 3;
  constexpr int MID = H / 2;

  auto run = [&](int cy0, int cy1, Pixel out[H][W]) {
    // melt_warp with noise disabled drips deterministically south.
    ::Feedback::Style style{};
    style.space_fn = &::Feedback::melt_warp;
    style.noise = nullptr;
    style.speed = 6.0f;
    style.fade = 0.9f;
    style.downsample = 4;
    Pipeline<W, H, FB> pipe{FB(style)};

    hs_test::StubEffect fx(W, H);
    fx.set_clip(cy0, cy1, 0, W);

    // Frame 0: seed a bright band straddling the boundary, THROUGH the clipped
    // pipeline, so a band-clipped worker only seeds rows inside its band.
    {
      Canvas c(fx);
      auto &frame = pipe.begin_frame(c, 1.0f);
      for (int y = MID - 4; y < MID + 4; ++y)
        for (int x = 0; x < W; ++x)
          frame.plot(c, x, y, Pixel(40000, 40000, 40000), 0.0f, 1.0f);
    }
    fx.advance_display();

    for (int f = 0; f < K; ++f) {
      {
        Canvas c(fx);
        (void)pipe.begin_frame(c, 1.0f);
      }
      fx.advance_display();
    }
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x)
        out[y][x] = fx.get_pixel(x, y);
  };

  static Pixel full[H][W], band_bot[H][W];
  run(0, H, full); // full-frame: what needs_full_frame() yields per worker
  run(MID, H, band_bot); // a band-clipped worker (the un-gated path)

  // The bottom band must differ, and the full frame must have lit it.
  bool differs = false;
  int full_bot_lit = 0;
  for (int y = MID; y < H; ++y)
    for (int x = 0; x < W; ++x) {
      const Pixel &a = full[y][x];
      const Pixel &b = band_bot[y][x];
      if (a.r != b.r || a.g != b.g || a.b != b.b)
        differs = true;
      if (a.r | a.g | a.b)
        ++full_bot_lit;
    }
  HS_EXPECT_TRUE(differs);
  HS_EXPECT_GT(full_bot_lit, 0);
}
