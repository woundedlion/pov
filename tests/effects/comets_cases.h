/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Comets white-box cases. */

/**
 * @brief White-box accessor for Comets' Lissajous-loop closing snap.
 * @details Every authored entry must close: path_fn(domain) == path_fn(0).
 */
struct CometsWhiteBox {
  /**
   * @brief Verifies every authored function table entry closes its loop.
   * @details The snapped endpoint must coincide with the t=0 start (0,1,0),
   *          and the snap must stay positive.
   */
  static void check_paths_close() {
    using C = Comets<DEFAULT_W, DEFAULT_H>;
    int idx = 0;
    for (const PresetEntry<math::LissajousParams> &entry : C::PRESETS) {
      const math::LissajousParams &cfg = entry.params;
      const float cd = C::closing_domain(cfg);
      HS_EXPECT_GT(cd, 0.0f); // floor-at-1 keeps the head moving
      // m2*domain >= PI keeps the floor from rewriting an authored domain.
      HS_EXPECT_GE(cfg.m2 * cfg.domain, math::PI_F);
      const math::Vector start = math::lissajous(cfg.m1, cfg.m2, cfg.a, 0.0f);
      const math::Vector end = math::lissajous(cfg.m1, cfg.m2, cfg.a, cd);
      const float gap = (end - start).magnitude();
      if (gap > 1e-3f)
        std::printf("  COMETS entry %d does not close: gap=%.6f\n", idx,
                    static_cast<double>(gap));
      HS_EXPECT_LT(gap, 1e-3f);
      ++idx;
    }
  }

  using C = Comets<SMALL_W, SMALL_H>;
  static constexpr int wipe_frames() { return C::WIPE_FRAMES; }
  static void roll_palette_over(C &c) { c.update_palette(); }
  static int wipe_frames_remaining(const C &c) {
    return c.wipe.frames_remaining;
  }
  static const GenerativePalette::Snapshot &palette_start(const C &c) {
    return c.wipe.start;
  }
  static const GenerativePalette::Snapshot &palette_target(const C &c) {
    return c.wipe.target;
  }
  static const math::Quaternion &node_orientation(const C &c) {
    return c.node->orientation.get();
  }
  static size_t trail_length(const C &c) { return c.node->trail.length(); }
  static int dwell_frames() { return C::PRESET_DWELL_FRAMES; }

  /**
   * @brief A disc recorded this frame leaves dark the columns swept before
   *        its sub-pose was reached.
   */
  static void check_births_skip_earlier_columns() {
    constexpr int W = 288, H = 144;
    using F = Comets<W, H>;
    reset_effect_globals();
    hs::random().seed(1337u);
    pin_frame_clock(0);
    F fx;
    fx.init();
    fx.params.cycle_duration = F::CYCLE_DURATION_MIN;
    constexpr float PIXEL = math::TWO_PI_F / W;
    int checked = 0;
    for (int f = 0; f < 48; ++f) {
      pin_frame_clock(f);
      fx.draw_frame();
      fx.advance_display();
      const auto &trail = fx.node->trail;
      const auto &newest = trail.get(trail.length() - 1);
      const int last = newest.length() - 1;
      const float reach = fx.params.thickness + 3.0f * PIXEL;
      const auto centre = [&](const math::Quaternion &q) {
        return fx.orientation.orient(math::rotate(fx.node->v, q));
      };
      for (int j = 1; j < last; ++j) {
        const math::Vector c = centre(newest.get(j));
        const math::PixelCoords px = math::vector_to_pixel<W, H>(c);
        const int x = static_cast<int>(lroundf(px.x)) % W;
        const int y = static_cast<int>(lroundf(px.y));
        const float fraction = Spacetime::SamplingSurface<W>::time(x);
        if (y < 0 || y >= H ||
            fraction >= static_cast<float>(j) / last - 2.0f / (W / 2))
          continue;
        bool clear = true;
        for (size_t i = 0; i < trail.length() && clear; ++i) {
          const auto &record = trail.get(i);
          for (int k = 0; k < record.length() && clear; ++k) {
            const bool same_record = i + 1 == trail.length();
            if (same_record && static_cast<float>(k) / last > fraction)
              continue;
            const float angle = math::fast_acos(
                hs::clamp(math::dot(centre(record.get(k)), c), -1.0f, 1.0f));
            clear = angle > reach;
          }
        }
        if (!clear)
          continue;
        ++checked;
        HS_EXPECT_TRUE(is_black(fx.get_pixel(x, y)));
      }
    }
    HS_EXPECT_TRUE(checked > 0);
  }

  /**
   * @brief Scans the trail at a camera held still.
   * @param c Effect whose trail is drawn.
   * @param palette_moves Whether the LUT moved since the last frame.
   * @return Displayed frame.
   */
  static std::vector<Pixel> draw_trail(C &c, bool palette_moves) {
    {
      Canvas canvas(c);
      c.draw_trail(canvas, palette_moves);
    }
    c.advance_display();
    std::vector<Pixel> out;
    capture_frame<SMALL_W, SMALL_H>(c, out);
    return out;
  }
  static void hold_camera(C &c) { c.orientation.collapse(); }
  static const BakedPalette &baked_palette(const C &c) {
    return c.baked_palette.view();
  }
  static const BakedPalette &drawn_palette(const C &c) {
    return c.drawn_palette.view();
  }
  static void rebake_crossfade(C &c, const BakedPalette &from,
                               const BakedPalette &to, float w) {
    c.baked_palette.rebake_crossfade(from, to, w);
  }
  static void restore_palette(C &c, const BakedPalette &saved) {
    c.baked_palette.rebake_copy(saved);
  }
  static bool wipe_pending(const C &c) { return c.wipe.pending; }
};

/** @brief Channel bound of a Comets column against its crossfaded-LUT
 *  oracle. */
inline constexpr int COMETS_WIPE_BOUND = 2;

/**
 * @brief Whether two baked palettes hold the same entries.
 * @param a First palette.
 * @param b Second palette.
 * @return True when every entry's colour and alpha match.
 */
inline bool same_palette_entries(const BakedPalette &a, const BakedPalette &b) {
  for (int i = 0; i < BakedPalette::LUT_SIZE; ++i) {
    const float t = static_cast<float>(i) / (BakedPalette::LUT_SIZE - 1);
    const Color4 x = a.get(t);
    const Color4 y = b.get(t);
    if (!(x.color == y.color) || x.alpha != y.alpha)
      return false;
  }
  return true;
}

/**
 * @brief Each column of a trail drawn during a palette wipe shades its discs
 *        from the LUT crossfaded to its sweep time, on both arms.
 */
inline void test_comets_wipe_shades_at_column_time() {
  using WB = CometsWhiteBox;
  reset_effect_globals();
  WB::C effect;
  effect.init();
  hs_test::run_frames(effect, 24);
  WB::roll_palette_over(effect);
  hs_test::run_frames(effect, 12);
  WB::hold_camera(effect);
  HS_EXPECT_FALSE(same_palette_entries(WB::drawn_palette(effect),
                                       WB::baked_palette(effect)));

  ScratchScope scope(scratch_arena_b);
  BakedPaletteStorage saved;
  saved.clone_from(WB::baked_palette(effect), scratch_arena_b);
  const std::vector<Pixel> drawn = WB::draw_trail(effect, true);
  const std::vector<Pixel> expected =
      column_time_oracle<SMALL_W, SMALL_H>([&](float u) {
        WB::rebake_crossfade(effect, WB::drawn_palette(effect), saved, u);
        const std::vector<Pixel> frame = WB::draw_trail(effect, false);
        WB::restore_palette(effect, saved);
        return frame;
      });
  HS_EXPECT_TRUE(expect_frames_within("comets wipe", drawn, expected,
                                      COMETS_WIPE_BOUND) > 30);
}

/**
 * @brief During a wipe the drawn LUT is the one the previous frame scanned;
 *        the arming frame and frames without a wipe hold the palette.
 */
inline void test_comets_wipe_history() {
  using WB = CometsWhiteBox;
  reset_effect_globals();
  WB::C effect;
  effect.init();
  hs_test::run_frames(effect, 4);
  WB::roll_palette_over(effect);
  HS_EXPECT_TRUE(WB::wipe_pending(effect));
  ScratchScope scope(scratch_arena_b);
  BakedPaletteStorage before;
  before.clone_from(WB::baked_palette(effect), scratch_arena_b);
  hs_test::run_frames(effect, 1);
  HS_EXPECT_TRUE(same_palette_entries(WB::baked_palette(effect), before));
  for (int f = 0; f < 3; ++f) {
    before.rebake_copy(WB::baked_palette(effect));
    hs_test::run_frames(effect, 1);
    HS_EXPECT_TRUE(same_palette_entries(WB::drawn_palette(effect), before));
    HS_EXPECT_FALSE(same_palette_entries(WB::baked_palette(effect), before));
  }
}

/**
 * @brief Field-wise equality for two generative-palette snapshots.
 * @param a First snapshot.
 * @param b Second snapshot.
 * @return true when every authored key and axis range matches.
 */
inline bool palette_snapshots_equal(const GenerativePalette::Snapshot &a,
                                    const GenerativePalette::Snapshot &b) {
  if (a.key_count != b.key_count || a.lightness_low != b.lightness_low ||
      a.lightness_high != b.lightness_high || a.chroma_low != b.chroma_low ||
      a.chroma_high != b.chroma_high ||
      a.lightness_curve != b.lightness_curve ||
      a.chroma_curve != b.chroma_curve)
    return false;
  for (size_t i = 0; i < a.keys.size(); ++i)
    if (a.keys[i].bytes != b.keys[i].bytes)
      return false;
  return true;
}

/**
 * @brief Pins Comets' mid-wipe rollover skip.
 * @details A rollover while a ColorWipe is animating is dropped, since the
 *          live ColorWipe holds references to palette_start/palette_target. A
 *          rollover arms again once the wipe drains.
 */
inline void test_comets_rollover_skipped_mid_wipe() {
  using WB = CometsWhiteBox;
  reset_effect_globals();

  WB::C effect;
  effect.init();

  WB::roll_palette_over(effect);
  HS_EXPECT_EQ(WB::wipe_frames_remaining(effect), WB::wipe_frames());
  const GenerativePalette::Snapshot start = WB::palette_start(effect);
  const GenerativePalette::Snapshot target = WB::palette_target(effect);

  WB::roll_palette_over(effect);
  HS_EXPECT_EQ(WB::wipe_frames_remaining(effect), WB::wipe_frames());
  HS_EXPECT_TRUE(palette_snapshots_equal(start, WB::palette_start(effect)));
  HS_EXPECT_TRUE(palette_snapshots_equal(target, WB::palette_target(effect)));

  // step_wipe_rebake burns one armed frame plus WIPE_FRAMES stepped ones.
  hs_test::run_frames(effect, WB::wipe_frames() + 1);
  HS_EXPECT_EQ(WB::wipe_frames_remaining(effect), 0);

  WB::roll_palette_over(effect);
  HS_EXPECT_EQ(WB::wipe_frames_remaining(effect), WB::wipe_frames());
}

/** @brief Manual Comets selection restarts the authored path deterministically. */
inline void test_comets_manual_preset_restarts_path() {
  using WB = CometsWhiteBox;
  reset_effect_globals();

  WB::C effect;
  effect.init();
  hs_test::run_frames(effect, 7);

  constexpr size_t PRESET = 4;
  HS_EXPECT_TRUE(effect.selectPreset(PRESET));
  HS_EXPECT_EQ(WB::node_orientation(effect), math::Quaternion());
  HS_EXPECT_EQ(WB::trail_length(effect), 0u);

  effect.draw_frame();
  effect.advance_display();
  const math::Quaternion first_step = WB::node_orientation(effect);
  HS_EXPECT_EQ(WB::trail_length(effect), 1u);

  hs_test::run_frames(effect, 7);
  HS_EXPECT_TRUE(effect.selectPreset(PRESET));
  HS_EXPECT_EQ(WB::node_orientation(effect), math::Quaternion());
  HS_EXPECT_EQ(WB::trail_length(effect), 0u);

  effect.draw_frame();
  effect.advance_display();
  HS_EXPECT_EQ(WB::node_orientation(effect), first_step);

  // Manual selection pauses the automatic dwell countdown as well as
  // selecting the requested path.
  hs_test::run_frames(effect, WB::dwell_frames());
  HS_EXPECT_EQ(effect.getPresetIndex(), PRESET);
}
