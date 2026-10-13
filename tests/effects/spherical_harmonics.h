/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 *
 * SphericalHarmonics white-box tests.
 */

inline constexpr std::array<int, 24> SH_PRESET_MODES{
    6,  1,  2,  3,  4,  5,  7,  8,  9,  10, 11, 12,
    13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};

/**
 * @brief White-box accessor for SphericalHarmonics' morph chain and shader inputs.
 * @details Befriended in effects/SphericalHarmonics.h.
 */
struct SphericalHarmonicsWhiteBox {
  using SH = SphericalHarmonics<SMALL_W, SMALL_H>;
  using Field = SH::HarmonicField;
  using Pipeline = SH::RenderPipeline;
  using PipelineFrame = Pipeline::Frame;

  static int current_idx(const SH &fx) { return fx.current_idx; }
  static int next_idx(const SH &fx) { return fx.next_idx; }
  static float morph_alpha(const SH &fx) { return fx.morph_alpha; }
  static math::Quaternion orientation(const SH &fx) {
    return fx.orientation.get();
  }
  static const BakedPalette &palette(const SH &fx) { return fx.baked_palette; }
  static float amplitude(const SH &fx) { return fx.params.amplitude; }
  static constexpr int max_mode_idx() { return SH::MAX_MODE_IDX; }

  static PipelineFrame prepare_pipeline(const SH &fx, int l1, int m1, int l2,
                                        int m2, float blend,
                                        const math::Quaternion &orientation,
                                        float amplitude) {
    return Pipeline::prepare({{l1, m1, l2, m2, blend, orientation},
                              {&fx.baked_palette.view(), amplitude}});
  }

  static Color4 shade_pipeline(const math::Vector &view,
                               const PipelineFrame &frame) {
    return Pipeline::evaluate(view, frame.ctx, frame.prepared);
  }

  static Color4 shade_direct(const SH &fx, float value, float amplitude) {
    return SH::colorize_harmonic(value, {&fx.baked_palette.view(), amplitude});
  }

  // Pinning both morph endpoints on one mode makes the blend an identity, so a
  // frame renders one pure harmonic whatever morph_alpha has reached.
  static void pin_mode(SH &fx, int idx) {
    fx.current_idx = idx;
    fx.next_idx = idx;
  }
  static void set_next_idx(SH &fx, int idx) { fx.next_idx = idx; }
  static void step_timeline(SH &fx, Canvas &canvas) {
    fx.timeline.step(canvas);
  }
};

/**
 * @brief Pins HarmonicField's blend endpoints and local frame.
 */
inline void test_sh_field_blend_endpoints_and_frame() {
  using Field = SphericalHarmonicsWhiteBox::Field;
  const math::Quaternion spin =
      math::make_rotation(math::Vector(0.3f, 0.8f, -0.5f).normalized(), 1.1f);
  const math::Quaternion identity;

  constexpr int LA = 1, MA = 0, LB = 3, MB = 2;
  Field mix_start(LA, MA, LB, MB, 0.0f, spin);
  Field mix_end(LA, MA, LB, MB, 1.0f, spin);
  Field pure_a(LA, MA, LA, MA, 0.0f, spin);
  Field pure_b(LB, MB, LB, MB, 0.0f, spin);
  Field pure_a_unrotated(LA, MA, LA, MA, 0.0f, identity);

  constexpr int PHI_STEPS = 24, THETA_STEPS = 32;
  int positives = 0, negatives = 0;
  double worst_start = 0.0, worst_end = 0.0, worst_frame = 0.0;
  for (int i = 0; i <= PHI_STEPS; ++i) {
    const float phi = math::PI_F * i / PHI_STEPS;
    for (int j = 0; j < THETA_STEPS; ++j) {
      const float theta = 2.0f * math::PI_F * j / THETA_STEPS;
      const math::Vector p(sinf(phi) * cosf(theta), cosf(phi),
                           sinf(phi) * sinf(theta));
      const float start = mix_start.sample(p);
      const float end = mix_end.sample(p);
      const float a = pure_a.sample(p);
      const float b = pure_b.sample(p);
      const float unspun = pure_a_unrotated.sample(p);
      // Rotating the sample by the same quaternion the shape carries must land
      // back on the unrotated shape's value.
      const float spun = pure_a.sample(math::rotate(p, spin));

      if (a > 0.02f)
        ++positives;
      if (a < -0.02f)
        ++negatives;
      worst_start =
          hs_test::fold_worst<double>(worst_start, std::fabs(start - a));
      worst_end = hs_test::fold_worst<double>(worst_end, std::fabs(end - b));
      worst_frame =
          hs_test::fold_worst<double>(worst_frame, std::fabs(spun - unspun));
    }
  }

  // A dipole must take both signs, or the endpoint and frame checks could pass
  // on a zero field.
  HS_EXPECT_GT(positives, 0);
  HS_EXPECT_GT(negatives, 0);

  if (worst_start > 0.0 || worst_end > 1e-5 || worst_frame > 1e-5)
    std::printf("  SH field: blend0 err=%g blend1 err=%g frame err=%g\n",
                worst_start, worst_end, worst_frame);
  // Blend 0 is exact; blend 1 reaches the second mode within lerp rounding.
  HS_EXPECT_EQ(worst_start, 0.0);
  HS_EXPECT_LT(worst_end, 1e-5);
  HS_EXPECT_LT(worst_frame, 1e-5);
}

/** @brief Keeps every shipped harmonic morph inside SampleSphere's signed range. */
inline void test_sh_field_stays_inside_unit_range() {
  using WB = SphericalHarmonicsWhiteBox;
  const math::Quaternion orientation =
      math::make_rotation(math::Vector(0.3f, 0.8f, -0.5f).normalized(), 1.1f);
  constexpr std::array<float, 5> BLENDS{0.0f, 0.25f, 0.5f, 0.75f, 1.0f};
  constexpr int PHI_STEPS = 48;
  constexpr int THETA_STEPS = 64;
  float peak = 0.0f;

  for (int first = 1; first <= WB::max_mode_idx(); ++first) {
    const int second = first % WB::max_mode_idx() + 1;
    auto [l1, m1] = SHMath::decode_lm(first);
    auto [l2, m2] = SHMath::decode_lm(second);
    for (float blend : BLENDS) {
      const WB::Field field(l1, m1, l2, m2, blend, orientation);
      for (int i = 0; i <= PHI_STEPS; ++i) {
        const float phi = math::PI_F * i / PHI_STEPS;
        for (int j = 0; j < THETA_STEPS; ++j) {
          const float theta = math::TWO_PI_F * j / THETA_STEPS;
          const math::Vector point(sinf(phi) * cosf(theta), cosf(phi),
                                   sinf(phi) * sinf(theta));
          peak = hs_test::fold_worst(peak, std::fabs(field.sample(point)));
        }
      }
    }
  }

  std::printf("  SH field bound: peak %.6f\n", static_cast<double>(peak));
  HS_EXPECT_LT(peak, 1.0f);
}

/** @brief Bounds color drift from the signed-to-unit carrier round trip. */
inline void test_sh_pullback_matches_direct_colorize() {
  using WB = SphericalHarmonicsWhiteBox;
  reset_effect_globals();
  WB::SH fx;
  fx.init();

  const math::Quaternion orientation =
      math::make_rotation(math::Vector(-0.4f, 0.2f, 0.9f).normalized(), 0.83f);
  constexpr std::array<float, 3> BLENDS{0.0f, 0.5f, 1.0f};
  constexpr std::array<float, 3> AMPLITUDES{0.2f, 3.2f, 7.0f};
  constexpr int PHI_STEPS = 24;
  constexpr int THETA_STEPS = 32;
  int max_channel_delta = 0;
  int differing_pixels = 0;

  for (int first = 1; first <= WB::max_mode_idx(); ++first) {
    const int second = first % WB::max_mode_idx() + 1;
    auto [l1, m1] = SHMath::decode_lm(first);
    auto [l2, m2] = SHMath::decode_lm(second);
    for (float blend : BLENDS) {
      const WB::Field direct_field(l1, m1, l2, m2, blend, orientation);
      for (float amplitude : AMPLITUDES) {
        const WB::PipelineFrame frame = WB::prepare_pipeline(
            fx, l1, m1, l2, m2, blend, orientation, amplitude);
        for (int i = 0; i <= PHI_STEPS; ++i) {
          const float phi = math::PI_F * i / PHI_STEPS;
          for (int j = 0; j < THETA_STEPS; ++j) {
            const float theta = math::TWO_PI_F * j / THETA_STEPS;
            const math::Vector point(sinf(phi) * cosf(theta), cosf(phi),
                                     sinf(phi) * sinf(theta));
            const Color4 direct =
                WB::shade_direct(fx, direct_field.sample(point), amplitude);
            const Color4 pullback = WB::shade_pipeline(point, frame);
            const int red = std::abs(static_cast<int>(direct.color.r) -
                                     static_cast<int>(pullback.color.r));
            const int green = std::abs(static_cast<int>(direct.color.g) -
                                       static_cast<int>(pullback.color.g));
            const int blue = std::abs(static_cast<int>(direct.color.b) -
                                      static_cast<int>(pullback.color.b));
            const int delta = std::max({red, green, blue});
            max_channel_delta = std::max(max_channel_delta, delta);
            differing_pixels += static_cast<int>(delta != 0);
            HS_EXPECT_EQ(direct.alpha, pullback.alpha);
          }
        }
      }
    }
  }

  std::printf("  SH pullback parity: max channel delta %d, differing %d\n",
              max_channel_delta, differing_pixels);
  HS_EXPECT_LE(max_channel_delta, 1);
}

/**
 * @brief Renders one frame of a mode-pinned SphericalHarmonics and inspects it.
 * @param idx Flat harmonic index pinned on both morph endpoints.
 * @param amplitude Value written to the Amplitude slider before the frame.
 * @param inspect Callback receiving (effect, field the shader saw, amplitude).
 */
template <typename FnT>
inline void sh_render_pinned_mode(int idx, float amplitude, FnT &&inspect) {
  using WB = SphericalHarmonicsWhiteBox;
  reset_effect_globals();
  pin_frame_clock(0);
  WB::SH fx;
  fx.init();
  WB::pin_mode(fx, idx);
  HS_EXPECT_EQ(fx.updateParameter("Amplitude", amplitude),
               ParamSetResult::APPLIED);
  fx.draw_frame();
  fx.advance_display();

  auto [l, m] = SHMath::decode_lm(idx);
  WB::Field field(l, m, l, m, WB::morph_alpha(fx), WB::orientation(fx));
  inspect(fx, field, WB::amplitude(fx));
  hs::clear_mock_time();
}

/**
 * @brief Pins the diverging palette split and the ambient-occlusion shaping.
 * @details Classifies each rendered pixel by the field's sign and magnitude: a
 *          saturated positive lobe wears the palette's top color unmodified, the
 *          negative lobe wears it with red and blue swapped and green dimmed,
 *          and below the AO falloff every pixel is dimmed under its palette
 *          entry.
 */
inline void test_sh_polarity_split_and_ao_shaping() {
  using WB = SphericalHarmonicsWhiteBox;
  // idx 2 is (l=1, m=0): one nodal circle, so both polarities cover a wide band.
  constexpr int DIPOLE_IDX = 2;

  // Amplitude 10 saturates the palette across most of both lobes.
  sh_render_pinned_mode(
      DIPOLE_IDX, 10.0f,
      [](const WB::SH &fx, const WB::Field &field, float amp) {
        constexpr float SATURATED =
            1.5f; // |val| * amp; palette saturates at 1.0
        Pixel pos(0, 0, 0), neg(0, 0, 0);
        int pos_n = 0, neg_n = 0, pos_split = 0, neg_split = 0;
        uint64_t pos_green = 0, neg_green = 0;
        for (int y = 0; y < SMALL_H; ++y)
          for (int x = 0; x < SMALL_W; ++x) {
            const float val =
                field.sample(math::pixel_to_vector<SMALL_W, SMALL_H>(x, y));
            const Pixel &px = fx.get_pixel(x, y);
            const bool positive = val > 0.0f;
            // A dipole's two lobes are antipodal mirrors, so their magnitude
            // distributions match and their green totals differ only by the
            // negative recolor's scale.
            (positive ? pos_green : neg_green) += px.g;

            const float mag = std::fabs(val) * amp;
            if (mag < SATURATED)
              continue;
            Pixel &slot = positive ? pos : neg;
            int &count = positive ? pos_n : neg_n;
            int &split = positive ? pos_split : neg_split;
            if (count++ == 0)
              slot = px;
            else if (px.r != slot.r || px.g != slot.g || px.b != slot.b)
              ++split;
          }

        std::printf("  SH polarity: pos=%d (%u,%u,%u) neg=%d (%u,%u,%u) "
                    "green %llu vs %llu\n",
                    pos_n, pos.r, pos.g, pos.b, neg_n, neg.r, neg.g, neg.b,
                    static_cast<unsigned long long>(pos_green),
                    static_cast<unsigned long long>(neg_green));
        HS_EXPECT_GT(pos_n, 0);
        HS_EXPECT_GT(neg_n, 0);
        // Past saturation the palette index, the seam weight and the AO factor are
        // all pinned, so every pixel of a lobe must carry one single color.
        HS_EXPECT_EQ(pos_split, 0);
        HS_EXPECT_EQ(neg_split, 0);

        // Positive lobe: the palette color, undimmed.
        const Pixel top = WB::palette(fx).get(1.0f).color;
        HS_EXPECT_EQ(pos.r, top.r);
        HS_EXPECT_EQ(pos.g, top.g);
        HS_EXPECT_EQ(pos.b, top.b);
        // Negative lobe: red and blue traded.
        HS_EXPECT_EQ(neg.r, pos.b);
        HS_EXPECT_EQ(neg.b, pos.r);
        // ...and the trade is visible, i.e. the two lobes are not the same color.
        HS_EXPECT_TRUE(pos.r != pos.b);
        // The green dimming is invisible at saturation (this palette's top color
        // has none), so pin it over the whole mirrored lobe instead.
        HS_EXPECT_GT(pos_green, 0u);
        HS_EXPECT_TRUE(neg_green * 10 < pos_green * 9);
      });

  // Amplitude 0.6 keeps every pixel below the AO falloff, so the whole frame
  // must read as a dimmed copy of the palette rather than the palette itself.
  sh_render_pinned_mode(
      DIPOLE_IDX, 0.6f,
      [](const WB::SH &fx, const WB::Field &field, float amp) {
        // Only sample where the palette entry is bright enough that a missing
        // occlusion factor would be unambiguous.
        constexpr uint32_t PALETTE_FLOOR = 4096;
        int probed = 0, undimmed = 0;
        for (int y = 0; y < SMALL_H; ++y)
          for (int x = 0; x < SMALL_W; ++x) {
            const float val =
                field.sample(math::pixel_to_vector<SMALL_W, SMALL_H>(x, y));
            const float mag = std::fabs(val) * amp;
            const Pixel want = WB::palette(fx).get(std::min(1.0f, mag)).color;
            const uint32_t want_energy =
                static_cast<uint32_t>(want.r) + want.g + want.b;
            if (want_energy < PALETTE_FLOOR)
              continue;
            const Pixel &px = fx.get_pixel(x, y);
            const uint32_t got_energy =
                static_cast<uint32_t>(px.r) + px.g + px.b;
            ++probed;
            // The occlusion factor cannot reach 0.8 at this amplitude.
            if (got_energy * 5 >= want_energy * 4)
              ++undimmed;
          }
        std::printf("  SH ambient occlusion: probed=%d undimmed=%d\n", probed,
                    undimmed);
        HS_EXPECT_GT(probed, 0);
        HS_EXPECT_EQ(undimmed, 0);
      });
}

/**
 * @brief Verifies the morph chain re-arms itself and keeps advancing modes.
 * @details start_morph() schedules a Transition whose then() callback commits
 *          the target and calls start_morph() again; three commits mean the
 *          callback re-armed twice.
 */
inline void test_sh_morph_chain_rearms() {
  using WB = SphericalHarmonicsWhiteBox;
  reset_effect_globals();
  pin_frame_clock(0);
  WB::SH fx;
  fx.init();

  const int seed = WB::current_idx(fx);
  HS_EXPECT_EQ(fx.getPresetIndex(), 0u);
  HS_EXPECT_GT(seed, 0); // never the constant harmonic
  HS_EXPECT_TRUE(WB::next_idx(fx) != seed);

  constexpr int FRAMES = 260;
  int commits = 0, held = seed, alpha_out_of_range = 0, self_blend = 0;
  int rearmed_at_zero = 0;
  float alpha_peak = 0.0f;
  std::vector<int> visited{seed};
  for (int f = 0; f < FRAMES; ++f) {
    pin_frame_clock(f);
    fx.draw_frame();
    fx.advance_display();

    const float alpha = WB::morph_alpha(fx);
    alpha_peak = std::max(alpha_peak, alpha);
    if (!(alpha >= 0.0f && alpha <= 1.0f))
      ++alpha_out_of_range;
    if (WB::next_idx(fx) == WB::current_idx(fx))
      ++self_blend;

    const int now = WB::current_idx(fx);
    if (now != held) {
      ++commits;
      held = now;
      const size_t expected_preset = static_cast<size_t>(
          std::find(SH_PRESET_MODES.begin(), SH_PRESET_MODES.end(), now) -
          SH_PRESET_MODES.begin());
      HS_EXPECT_EQ(fx.getPresetIndex(), expected_preset);
      // A committed leg rewinds the blend and schedules the next one.
      if (alpha == 0.0f)
        ++rearmed_at_zero;
      if (std::find(visited.begin(), visited.end(), now) == visited.end())
        visited.push_back(now);
    }
  }
  hs::clear_mock_time();

  std::printf("  SH morph: %d commits over %d frames, %zu distinct modes, "
              "alpha peak %.3f\n",
              commits, FRAMES, visited.size(), static_cast<double>(alpha_peak));
  HS_EXPECT_GE(commits, 3); // the then() callback re-armed at least twice
  HS_EXPECT_EQ(rearmed_at_zero, commits);
  HS_EXPECT_GE(visited.size(), 3u);
  HS_EXPECT_EQ(alpha_out_of_range, 0);
  HS_EXPECT_EQ(self_blend, 0); // blending a mode into itself would freeze it
  HS_EXPECT_GT(alpha_peak, 0.9f);
}

/** @brief Maps every runtime preset to its stable harmonic mode. */
inline void test_sh_preset_mode_mapping() {
  using WB = SphericalHarmonicsWhiteBox;
  reset_effect_globals();
  WB::SH fx;
  fx.init();

  HS_EXPECT_TRUE(fx.getParameters().find("Amplitude") != nullptr);
  HS_EXPECT_TRUE(fx.getParameters().find("Debug BB") == nullptr);
  HS_EXPECT_EQ(fx.updateParameter("Debug BB", 1.0f),
               ParamSetResult::UNKNOWN_PARAM);

  HS_EXPECT_EQ(fx.getPresetIndex(), 0u);
  HS_EXPECT_EQ(WB::current_idx(fx), SH_PRESET_MODES[0]);
  for (size_t preset = 0; preset < SH_PRESET_MODES.size(); ++preset) {
    HS_EXPECT_TRUE(fx.selectPreset(preset));
    HS_EXPECT_EQ(fx.getPresetIndex(), preset);
    HS_EXPECT_EQ(WB::current_idx(fx), SH_PRESET_MODES[preset]);
    HS_EXPECT_EQ(WB::next_idx(fx), SH_PRESET_MODES[preset]);
    HS_EXPECT_EQ(WB::morph_alpha(fx), 0.0f);
  }
}

/** @brief Keeps a manual harmonic selection past the automatic leg it replaces. */
inline void test_sh_manual_preset_replaces_inflight_morph() {
  using WB = SphericalHarmonicsWhiteBox;
  reset_effect_globals();
  pin_frame_clock(0);
  WB::SH fx;
  fx.init();

  HS_EXPECT_EQ(fx.getPresetCount(), SH_PRESET_MODES.size());
  HS_EXPECT_EQ(fx.getPresetIndex(), 0u);
  HS_EXPECT_EQ(WB::current_idx(fx), SH_PRESET_MODES[0]);
  for (int frame = 0; frame < 8; ++frame) {
    pin_frame_clock(frame);
    fx.draw_frame();
    fx.advance_display();
  }
  HS_EXPECT_GT(WB::morph_alpha(fx), 0.0f);
  const int replaced_target = WB::next_idx(fx);

  const size_t selected_preset = replaced_target == SH_PRESET_MODES.back()
                                     ? SH_PRESET_MODES.size() - 2
                                     : SH_PRESET_MODES.size() - 1;
  const int selected_mode = SH_PRESET_MODES[selected_preset];
  HS_EXPECT_NE(selected_mode, replaced_target);
  HS_EXPECT_TRUE(fx.selectPreset(selected_preset));
  HS_EXPECT_TRUE(fx.animations_paused());
  HS_EXPECT_EQ(fx.getPresetIndex(), selected_preset);
  HS_EXPECT_EQ(WB::current_idx(fx), selected_mode);
  HS_EXPECT_EQ(WB::next_idx(fx), selected_mode);
  HS_EXPECT_EQ(WB::morph_alpha(fx), 0.0f);

  fx.draw_frame();
  fx.advance_display();
  HS_EXPECT_EQ(WB::morph_alpha(fx), 0.0f);

  fx.setAnimationsPaused(false);
  for (int frame = 0; frame < 64; ++frame) {
    pin_frame_clock(frame + 9);
    fx.draw_frame();
    fx.advance_display();
  }
  hs::clear_mock_time();

  HS_EXPECT_EQ(WB::current_idx(fx), selected_mode);
  HS_EXPECT_EQ(fx.getPresetIndex(), selected_preset);
  HS_EXPECT_TRUE(WB::next_idx(fx) != selected_mode);
}
