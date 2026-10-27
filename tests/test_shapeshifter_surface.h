/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

/** @file test_shapeshifter_surface.h
 * @brief Column-time and edge-mask parity for the continuous star preset. */

struct SurfaceMotionCase {
  math::Quaternion start;               ///< Frame-start orientation.
  math::Quaternion end;                 ///< Frame-end orientation.
  float phase = 0.37f;                  ///< Frame-end waveform phase, turns.
  float amplitude = 1.0f;               ///< Radial twist amplitude.
  float step = 0.016f;                  ///< Waveform advance, turns.
  const char *name = "surface fixture"; ///< Diagnostic label.
};

/** @brief Dense star state used by the surface oracles.
 * @return P1 geometry and brightness. */
inline OracleState surface_state() {
  OracleState state;
  state.shape = OracleEffect::ShapeType::PLANAR_STAR;
  state.count = 288;
  state.sides = 7;
  state.spacing = OracleEffect::RadiusSpacing::SCREEN_BALANCED;
  state.alpha = 0.625f;
  return state;
}

/** @brief Replays the natural or reset trajectory to a frozen pose.
 * @param reset Uses the preset-reset sequence.
 * @param frame_number Last replay frame.
 * @return Captured frame endpoints and waveform phase. */
inline SurfaceMotionCase surface_replay_pose(bool reset, int frame_number) {
  reset_globals();
  hs::seed_effect_stream(hs::effect_seed_key<OracleEffect>("ShapeShifter"));
  hs::set_mock_time(0, 0);
  SurfaceMotionCase pose;
  {
    OracleEffect effect;
    effect.init();
    if (reset)
      effect.profile_select_preset(0);
    effect.set_clip(0, 72, 0, 144);
    for (int frame = 1; frame <= frame_number; ++frame) {
      hs::set_mock_time((frame - 1) * 62500ul / 1000, (frame - 1) * 62500ul);
      {
        Canvas canvas(effect);
        ShapeShifterWhiteBox::replay_step(effect, canvas);
      }
      effect.advance_display();
    }
    pose.start = ShapeShifterWhiteBox::orientation(effect).get(0.0f);
    pose.end = ShapeShifterWhiteBox::orientation(effect).get();
    pose.phase = ShapeShifterWhiteBox::phase(effect);
    pose.step = ShapeShifterWhiteBox::twist_step(effect);
  }
  Timeline().clear();
  hs::clear_mock_time();
  return pose;
}

inline std::vector<float> surface_kernel_motion(const SurfaceMotionCase &pose,
                                                bool timed, float fraction,
                                                float displacement,
                                                OracleClip clip = {}) {
  reset_globals();
  std::vector<float> out;
  {
    OracleEffect effect;
    effect.init();
    auto state = surface_state();
    state.phase = pose.phase;
    state.amplitude = pose.amplitude;
    ShapeShifterWhiteBox::configure(effect, state);
    auto &orientation = ShapeShifterWhiteBox::orientation(effect);
    orientation.set(pose.start);
    orientation.push(pose.end);
    ShapeShifterWhiteBox::set_twist(effect, pose.phase, pose.step);
    if (!timed) {
      const auto frozen = orientation.get(fraction);
      orientation.set(frozen);
      ShapeShifterWhiteBox::set_twist(
          effect, pose.phase - pose.step * (1 - fraction), 0);
    }
    if (timed || fraction == 1) {
      effect.set_clip(clip.y0, clip.y1, clip.x0, clip.x1);
      {
        Canvas canvas(effect);
        ShapeShifterWhiteBox::surface_kernel(effect, canvas, timed, nullptr,
                                             displacement);
      }
    } else {
      const int column = Spacetime::surface_column(fraction, 288);
      out.resize(288 * 144);
      for (int x : {column, column + 144}) {
        effect.set_clip(0, 144, x, x + 1);
        {
          Canvas canvas(effect);
          ShapeShifterWhiteBox::surface_kernel(effect, canvas, false, nullptr,
                                               displacement);
        }
        effect.advance_display();
        for (int y = 0; y < 144; ++y) {
          const auto p = effect.get_pixel(x, y);
          out[y * 288 + x] = (float(p.r) + p.g + p.b) / (3 * 65535);
        }
      }
    }
    if (timed || fraction == 1) {
      effect.advance_display();
      out = capture_luminance<288, 144>(effect);
    }
  }
  Timeline().clear();
  return out;
}

inline void test_surface_clips() {
  const auto pose = surface_replay_pose(true, 77);
  for (bool timed : {false, true}) {
    const auto full = surface_kernel_motion(pose, timed, 1, 0);
    std::vector<float> stitched(full.size());
    for (int y = 0; y < 144; y += 72)
      for (int x = 0; x < 288; x += 144) {
        const OracleClip clip{y, y + 72, x, x + 144};
        const auto tile = surface_kernel_motion(pose, timed, 1, 0, clip);
        for (int row = y; row < y + 72; ++row)
          for (int column = x; column < x + 144; ++column)
            stitched[row * 288 + column] = tile[row * 288 + column];
      }
    HS_EXPECT_EQ(stitched, full);
  }
}

inline void test_surface_seams() {
  SurfaceMotionCase fixture;
  fixture.start = math::make_rotation(math::X_AXIS, .27f);
  fixture.end =
      math::make_rotation(math::Vector(.2f, 1, .1f).normalized(), .03f) *
      fixture.start;
  auto zero = fixture;
  zero.step = 0;
  zero.amplitude = 0;
  zero.name = "surface zero twist";
  auto reset = surface_replay_pose(true, 77);
  reset.name = "surface reset";
  auto reversed = fixture;
  reversed.start = fixture.end;
  reversed.end = fixture.start;
  reversed.step = -fixture.step;
  reversed.name = "surface reverse";
  for (auto pose : {fixture, reversed, zero, reset})
    seam_compare_render<288, 144>(
        [&](bool timed, float t, float displacement) {
          return surface_kernel_motion(pose, timed, t, displacement);
        },
        pose.name, 1, true);
}

inline void surface_mask_case(OracleEffect &effect, const char *name,
                              bool write_images) {
  for (bool timed : {false, true})
    for (bool board : {false, true}) {
      effect.set_clip(0, board ? 72 : 144, 0, board ? 144 : 288);
      OracleFrame images[2];
      for (bool mask : {false, true}) {
        Scan::SurfaceAreaStats stats;
        {
          Canvas canvas(effect);
          ShapeShifterWhiteBox::surface_kernel(effect, canvas, timed, &stats, 0,
                                               mask);
        }
        effect.advance_display();
        for (int y = 0; y < 144; ++y)
          for (int x = 0; x < 288; ++x)
            images[mask].pixels.push_back(effect.get_pixel(x, y));
        if (mask) {
          std::printf(
              "SURFACE_MASK %s timed=%d board=%d pixels=%llu hits=%llu edges=%llu intervals=%llu max=%d prepares=%llu candidates=%llu pieces=%llu writes=%llu\n",
              name, int(timed), int(board),
              static_cast<unsigned long long>(stats.pixels),
              static_cast<unsigned long long>(stats.mask_hits),
              static_cast<unsigned long long>(stats.edge_pixels),
              static_cast<unsigned long long>(stats.intervals),
              stats.max_intervals,
              static_cast<unsigned long long>(stats.mask_prepares),
              static_cast<unsigned long long>(stats.mask_candidates),
              static_cast<unsigned long long>(stats.mask_pieces),
              static_cast<unsigned long long>(stats.mask_row_writes));
          HS_EXPECT_GT(stats.mask_hits, 0ull);
          HS_EXPECT_LT(stats.mask_hits, stats.pixels);
          HS_EXPECT_GT(stats.edge_pixels, 0ull);
          if (write_images)
            write_oracle_ppm(images[mask],
                             std::string(shapeshifter_image_prefix()) + "/" +
                                 name + (timed ? "-timed" : "-frozen") +
                                 (board ? "-board.ppm" : "-full.ppm"));
        }
      }
      const auto difference =
          compare_buffers(images[0].pixels, images[1].pixels);
      std::printf("SURFACE_MASK_DIFF %s timed=%d board=%d pixels=%zu\n", name,
                  int(timed), int(board), difference.different_pixels);
#if defined(HS_TEST_FAST_MATH)
      // Pixels with no candidate intervals take the point query, whose root
      // fast-math rounds apart from the integrator's by ulps of the fraction.
      HS_EXPECT_EQ(difference.coverage_mismatches, size_t{0});
      HS_EXPECT_LE(difference.max_absolute_error, 1u);
#else
      HS_EXPECT_EQ(difference.different_pixels, size_t{0});
#endif
    }
}

inline void surface_mask_cases(bool write_images) {
  reset_globals();
  {
    OracleEffect effect;
    effect.init();
    auto state = surface_state();
    state.phase = .37f;
    state.clip = {0, 144, 0, 288};
    ShapeShifterWhiteBox::configure(effect, state);
    auto &orientation = ShapeShifterWhiteBox::orientation(effect);
    const auto start = math::make_rotation(math::X_AXIS, .27f);
    orientation.set(start);
    orientation.push(
        math::make_rotation(math::Vector(.2f, 1, .1f).normalized(), .03f) *
        start);
    ShapeShifterWhiteBox::set_twist(effect, .37f, .016f);
    surface_mask_case(effect, "fixture", write_images);
  }
  Timeline().clear();
  for (bool reset : {false, true}) {
    reset_globals();
    hs::seed_effect_stream(hs::effect_seed_key<OracleEffect>("ShapeShifter"));
    hs::set_mock_time(0, 0);
    {
      OracleEffect effect;
      effect.init();
      if (reset)
        effect.profile_select_preset(0);
      const int last = reset ? 77 : 410;
      for (int frame = 1; frame <= last; ++frame) {
        hs::set_mock_time((frame - 1) * 62500ul / 1000, (frame - 1) * 62500ul);
        {
          Canvas canvas(effect);
          ShapeShifterWhiteBox::replay_step(effect, canvas);
        }
        effect.advance_display();
      }
      surface_mask_case(effect, reset ? "reset-77" : "shipping-410",
                        write_images);
    }
    Timeline().clear();
  }
  hs::clear_mock_time();
}

/** @brief Mask pruning preserves every rendered channel at all three poses. */
inline void test_surface_mask_parity() { surface_mask_cases(false); }

/** @brief Writes the surface comparison fixtures. */
inline void test_surface_images() { surface_mask_cases(true); }

/** @brief P1 dispatch reaches the continuous surface renderer. */
inline void test_surface_preset_dispatch() {
  auto state = surface_state();
  state.phase = 0.37f;
  state.orientation = math::make_rotation(math::X_AXIS, .27f);
  const auto comparison = compare_renders(
      state, candidate_renderer(), [](OracleEffect &effect, Canvas &canvas) {
        ShapeShifterWhiteBox::surface_kernel(effect, canvas, false);
      });
  HS_EXPECT_EQ(comparison.error.different_pixels, size_t{0});
  HS_EXPECT_GT(frame_energy(comparison.candidate), uint64_t{0});
}

/** @brief Samples board work throughout the complete natural P1 visit. */
inline void test_surface_natural_census() {
  reset_globals();
  hs::seed_effect_stream(hs::effect_seed_key<OracleEffect>("ShapeShifter"));
  hs::set_mock_time(0, 0);
  int samples = 0;
  {
    OracleEffect effect;
    effect.init();
    effect.set_clip(0, 72, 0, 144);
    for (int frame = 1; frame <= 471; ++frame) {
      hs::set_mock_time((frame - 1) * 62500ul / 1000, (frame - 1) * 62500ul);
      const bool selected =
          (frame - 1) % 16 == 0 || frame == 89 || frame == 410 || frame == 471;
      Scan::SurfaceAreaStats stats;
      {
        Canvas canvas(effect);
        ShapeShifterWhiteBox::replay_step(effect, canvas);
        if (selected) {
          HS_EXPECT_EQ(ShapeShifterWhiteBox::preset_index(effect), size_t{0});
          ShapeShifterWhiteBox::surface_kernel(effect, canvas, true, &stats);
        }
      }
      effect.advance_display();
      if (selected) {
        ++samples;
        HS_EXPECT_GE(stats.pixels, uint64_t{144 * 72});
        HS_EXPECT_GT(stats.mask_prepares, uint64_t{0});
        std::printf("SURFACE_CENSUS frame=%d pixels=%llu hits=%llu edges=%llu "
                    "point=%llu outside=%llu intervals=%llu max=%d maxpoint=%d "
                    "prepares=%llu candidates=%llu pieces=%llu writes=%llu\n",
                    frame, static_cast<unsigned long long>(stats.pixels),
                    static_cast<unsigned long long>(stats.mask_hits),
                    static_cast<unsigned long long>(stats.edge_pixels),
                    static_cast<unsigned long long>(stats.point_intervals),
                    static_cast<unsigned long long>(stats.outside_intervals),
                    static_cast<unsigned long long>(stats.intervals),
                    stats.max_intervals, stats.max_point_intervals,
                    static_cast<unsigned long long>(stats.mask_prepares),
                    static_cast<unsigned long long>(stats.mask_candidates),
                    static_cast<unsigned long long>(stats.mask_pieces),
                    static_cast<unsigned long long>(stats.mask_row_writes));
      }
    }
  }
  Timeline().clear();
  hs::clear_mock_time();
  HS_EXPECT_EQ(samples, 33);
}

/** @brief The maximum authored amplitude fits the production arena partition. */
inline void test_surface_max_amplitude() {
  auto state = surface_state();
  state.amplitude = 10;
  state.phase = .37f;
  state.clip = {0, 72, 0, 144};
  state.orientation = math::make_rotation(math::X_AXIS, .27f);
  for (bool timed : {false, true}) {
    const auto image =
        capture_frame(state, [&](OracleEffect &effect, Canvas &canvas) {
          const ScopedTemporalSampling sampling(timed);
          scratch_arena_a.reset_high_water_mark();
          scratch_arena_b.reset_high_water_mark();
          ShapeShifterWhiteBox::render_candidate(effect, canvas);
        });
    const auto a = scratch_arena_a.get_high_water_mark();
    const auto b = scratch_arena_b.get_high_water_mark();
    HS_EXPECT_GT(frame_energy(image), uint64_t{0});
    HS_EXPECT_GT(a, size_t{0});
    HS_EXPECT_GT(b, size_t{0});
    HS_EXPECT_LE(a, size_t{32 * 1024});
    HS_EXPECT_LE(b, size_t{32 * 1024});
    std::printf("  surface amplitude10 timed=%d scratch_a=%zu scratch_b=%zu\n",
                int(timed), a, b);
  }
}
