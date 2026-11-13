/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
/** @file
 * @brief Fishbowl white-box cases. */

/** @brief White-box accessor for Fishbowl's trail node. */
struct FishbowlWhiteBox {
  /** @brief Renders warped trail history or independently frozen endpoints.
   * @param timed Selects production temporal rendering.
   * @param fraction Frozen reference fraction.
   * @param displacement Reference row and column shift.
   * @param full Selects full-history rollover.
   * @return Frame luminance. */
  static std::vector<float> render_temporal_trail(bool timed, float fraction,
                                                  float displacement,
                                                  bool full) {
    using FX = Fishbowl<DEFAULT_W, DEFAULT_H>;
    reset_globals();
    FX fx;
    fx.init();
    fx.node->trail.clear();
    auto q_at = [](int i) {
      const auto v = math::pixel_to_vector<DEFAULT_W, DEFAULT_H>(
          119.0f + i * 3.7f, 65.0f + 5.0f * sinf(i * .3f));
      return math::make_rotation(math::cross(math::Y_AXIS, v).normalized(),
                                 acosf(v.y));
    };
    const int prior_records = full ? FX::TRAIL_LENGTH : 2;
    for (int record = 0; record < prior_records; ++record) {
      math::Orientation<FX::ORIENTATION_SUBSTEPS> pose;
      pose.set(q_at(record * 3));
      for (int j = 1; j < 4; ++j)
        pose.push(q_at(record * 3 + j));
      fx.node->trail.record(pose);
    }
    fx.orientation.set(math::Quaternion());
    fx.noise_xform.prepare_frame();
    auto &noise =
        const_cast<Animation::NoiseParams &>(fx.noise_xform.active_params(0));
    noise.time = .25f;
    const auto old_noise = noise;
    {
      Canvas canvas(fx);
      fx.render_trail(canvas, false, fx.cycle_phase);
    }
    fx.advance_display();
    math::Orientation<FX::ORIENTATION_SUBSTEPS> newest;
    newest.set(q_at(prior_records * 3));
    for (int j = 1; j < 4; ++j)
      newest.push(q_at(prior_records * 3 + j));
    fx.node->trail.record(newest);
    if (full)
      fx.history_head = 1;
    const auto end_pose = math::make_rotation(math::Y_AXIS, .02f);
    fx.orientation.set(end_pose);
    noise.time = .28f;
    fx.palette_speed = .1f;
    fx.cycle_phase = .1f;
    {
      Canvas canvas(fx);
      if (timed) {
        fx.render_trail(canvas, full, 0.0f);
      } else {
        ScratchScope guard(scratch_arena_a);
        ArenaVector<typename FX::TrailVertex> points(scratch_arena_a,
                                                     FX::MAX_FRAGMENTS);
        const float active = full ? FX::TRAIL_LENGTH : 3.0f;
        const float fill = active / FX::TRAIL_LENGTH;
        const auto displaced =
            displaced_reference_rotation<DEFAULT_W, DEFAULT_H>(
                fraction, displacement, displacement);
        math::Quaternion previous_q;
        math::Vector previous_pos;
        float previous_t = 0;
        bool have_previous = false;
        auto current = [&](const math::Quaternion &q) {
          return noise_transform(
              math::rotate(math::rotate(fx.node->v, q), end_pose), noise);
        };
        auto append = [&](const math::Quaternion &q, const math::Vector &end,
                          float t) {
          const auto old =
              math::Snorm3::encode(
                  noise_transform(math::rotate(fx.node->v, q), old_noise))
                  .decode();
          points.push_back(
              {math::rotate((old + (end - old) * fraction).normalized(),
                            displaced),
               t * fill});
        };
        deep_tween(fx.node->trail, [&](const math::Quaternion &q, float t) {
          const auto pos = current(q);
          if (have_previous) {
            const auto mid_q = math::slerp(previous_q, q, .5f);
            const auto mid = current(mid_q);
            if (math::angle_between(mid, math::slerp(previous_pos, pos, .5f)) >
                .25f * math::RADIANS_PER_COLUMN<DEFAULT_W>)
              append(mid_q, mid, .5f * (previous_t + t));
          }
          append(q, pos, t);
          previous_q = q;
          previous_pos = pos;
          previous_t = t;
          have_previous = true;
        });
        auto *ends = scratch_arena_a.allocate_n<math::Vector>(points.size());
        for (size_t i = 0; i < points.size(); ++i)
          ends[i] = points[i].pos;
        // The frozen trail is drawn through the timed route without motion,
        // shaded at the fraction throughout.
        auto shade = [&](float palette_t, float) {
          const float age = palette_t / fill;
          Color4 color = fx.shade_trail(
              palette_t -
                  (1.0f - fraction) * (.1f / fx.params.scale_factor -
                                       (full ? 1.0f / FX::TRAIL_LENGTH : 0.0f)),
              age);
          if (age > (active - 1.0f + fraction) / active)
            color.alpha = 0.0f;
          return color;
        };
        fx.draw_timed_trail(canvas, points, ends, 0.0f, shade);
      }
    }
    fx.advance_display();
    return capture_mean_channel<DEFAULT_W, DEFAULT_H>(fx);
  }

  using EffectType = Fishbowl<SMALL_W, SMALL_H>;
  template <int W, int H>
  static size_t tween_vertices(const Fishbowl<W, H> &fx) {
    size_t n = 0;
    deep_tween(fx.node->trail, [&](const math::Quaternion &, float) { ++n; });
    return n;
  }
  static Color4 sample_fire(const EffectType &fx, float t) {
    return fx.fire_palette.get(fx.duty_mod.modify(t));
  }
  static Color4 fire_gradient(const EffectType &fx, float u) {
    return fx.fire_palette.get(u);
  }
  static float palette_fill_scale(size_t trail_length) {
    return EffectType::palette_fill_scale(trail_length);
  }
  static float palette_phase_speed(float cycle_speed, float scale_factor,
                                   size_t trail_length) {
    return EffectType::palette_phase_speed(cycle_speed, scale_factor,
                                           trail_length);
  }
  static Color4 shade_trail(const EffectType &fx, float palette_t,
                            float age_t) {
    return fx.shade_trail(palette_t, age_t);
  }
  static bool needs_adaptive_midpoint(const math::Vector &a,
                                      const math::Vector &mid,
                                      const math::Vector &b) {
    return EffectType::needs_adaptive_midpoint(a, mid, b);
  }
};

/** @brief Snapshots enforce finite values and schema identity. */
inline void test_fishbowl_snapshot_admission() {
  reset_globals();
  using WB = FishbowlWhiteBox;
  WB::EffectType fx;
  fx.init();

  const auto snapshot = fx.serialize_parameters();
  auto edited = snapshot;
  edited.params.speed = 1.0f;
  HS_EXPECT_TRUE(fx.restore_parameters(edited));
  HS_EXPECT_EQ(fx.serialize_parameters().params.speed, 1.0f);
  edited.params.speed = std::numeric_limits<float>::quiet_NaN();
  HS_EXPECT_FALSE(fx.restore_parameters(edited));
  edited = snapshot;
  ++edited.schema_version;
  HS_EXPECT_FALSE(fx.restore_parameters(edited));
  HS_EXPECT_TRUE(fx.restore_parameters(snapshot));
}

/** @brief Preset selection restores authored controls and pauses animation. */
inline void test_fishbowl_preset_restores_params() {
  reset_globals();
  using WB = FishbowlWhiteBox;
  WB::EffectType fx;
  fx.init();

  const auto snapshot = fx.serialize_parameters();
  auto value = [&](const char *name) {
    const auto *def = required_param(fx, name);

    return def ? def->get() : -1.0f;
  };

  hs_test::expect_params_in_range(fx, "Fishbowl preset");
  HS_EXPECT_EQ(fx.getPresetCount(), 1u);
  HS_EXPECT_EQ(fx.getPresetIndex(), 0u);

  HS_EXPECT_EQ(fx.updateParameter("Cycle Dur", 111.0f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(fx.updateParameter("Duty Cycle", 0.25f),
               ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(fx.selectPreset(0));
  HS_EXPECT_TRUE(fx.animations_paused());

  HS_EXPECT_EQ(value("Cycle Dur"), snapshot.params.cycle_duration);
  HS_EXPECT_EQ(value("Duty Cycle"), snapshot.params.duty_cycle);
}

/** @brief The fire gradient fills the lit duty fraction and is black after it. */
inline void test_fishbowl_fire_duty_cycle() {
  reset_globals();
  using WB = FishbowlWhiteBox;
  WB::EffectType fx;
  fx.init();

  auto value = [&](const char *name) {
    const auto *def = required_param(fx, name);

    return def ? def->get() : -1.0f;
  };

  HS_EXPECT_EQ(fx.updateParameter("Duty Cycle", 0.5f), ParamSetResult::APPLIED);
  const Pixel BLACK(0, 0, 0);
  for (float t : {0.20f, 0.39f}) {
    HS_CONTEXT(t < 0.3f ? "early lit" : "late lit");
    const Color4 lit = WB::sample_fire(fx, t);
    HS_EXPECT_EQ(lit.color, WB::fire_gradient(fx, t / 0.5f).color);
    HS_EXPECT_TRUE(lit.color != BLACK);
  }
  HS_EXPECT_EQ(WB::sample_fire(fx, 0.50f).color, BLACK);
  HS_EXPECT_EQ(WB::sample_fire(fx, 0.75f).color, BLACK);
  const float dark_palette_t = 0.75f / value("Scale Factor");
  HS_EXPECT_EQ(WB::shade_trail(fx, dark_palette_t, 0.50f).alpha, 0.0f);

  HS_EXPECT_EQ(fx.updateParameter("Duty Cycle", 0.25f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(WB::sample_fire(fx, 0.30f).color, BLACK);
}

/** @brief Palette phase compensates trail fill consistently. */
inline void test_fishbowl_palette_phase_compensates_fill() {
  reset_globals();
  using WB = FishbowlWhiteBox;
  WB::EffectType fx;
  fx.init();

  auto value = [&](const char *name) {
    const auto *def = required_param(fx, name);

    return def ? def->get() : -1.0f;
  };

  HS_EXPECT_EQ(WB::palette_fill_scale(1), 1.0f / WB::EffectType::TRAIL_LENGTH);
  HS_EXPECT_EQ(WB::palette_fill_scale(WB::EffectType::TRAIL_LENGTH / 2),
               static_cast<float>(WB::EffectType::TRAIL_LENGTH / 2) /
                   WB::EffectType::TRAIL_LENGTH);
  HS_EXPECT_EQ(WB::palette_fill_scale(WB::EffectType::TRAIL_LENGTH), 1.0f);
  const float SAMPLE_INDEX = 12.0f;
  const float short_coord = (SAMPLE_INDEX / 32.0f) * WB::palette_fill_scale(32);
  const float long_coord = (SAMPLE_INDEX / 96.0f) * WB::palette_fill_scale(96);
  HS_EXPECT_NEAR(short_coord, long_coord, 1e-7f);

  const float cycle_speed = value("Cycle Speed");
  const float scale_factor = value("Scale Factor");
  const float coordinate_shift = scale_factor / WB::EffectType::TRAIL_LENGTH;
  const float filling_speed = WB::palette_phase_speed(
      cycle_speed, scale_factor, WB::EffectType::TRAIL_LENGTH - 1);
  const float saturated_speed = WB::palette_phase_speed(
      cycle_speed, scale_factor, WB::EffectType::TRAIL_LENGTH);
  HS_EXPECT_NEAR(filling_speed, cycle_speed - coordinate_shift, 1e-7f);
  HS_EXPECT_NEAR(saturated_speed - coordinate_shift, filling_speed, 1e-7f);
}

/** @brief Adaptive midpoint detection accepts the spherical midpoint. */
inline void test_fishbowl_adaptive_midpoint() {
  reset_globals();
  using WB = FishbowlWhiteBox;
  WB::EffectType fx;
  fx.init();

  const math::Vector A = math::X_AXIS;
  const math::Vector B = math::Y_AXIS;
  HS_EXPECT_FALSE(WB::needs_adaptive_midpoint(A, math::slerp(A, B, 0.5f), B));
  HS_EXPECT_TRUE(WB::needs_adaptive_midpoint(A, math::Z_AXIS, B));
}

/**
 * @brief Verifies the scratch-A split's predicted worst cases really cover a
 *        saturated Fishbowl frame, still and timed.
 * @details STILL_SCRATCH_A covers the vertex buffer, the Multiline fragment
 *          buffer it binds and rasterize's sub-step cache; TIMED_SCRATCH_A
 *          covers the vertices, their frame-end positions, the drawn vertices
 *          and one column-walk chunk.
 */
template <int W, int H>
inline void check_fishbowl_scratch_estimate_covers_peak() {
  using WB = FishbowlWhiteBox;
  using EffectType = Fishbowl<W, H>;
  const ScopedTemporalSampling temporal_sampling(!Spacetime::flat_surface());
  for (const bool timed : {false, true}) {
    reset_globals();
    Spacetime::set_flat_surface(!timed);
    EffectType fx;
    fx.init();
    scratch_arena_a.reset_high_water_mark();

    size_t worst_vertices = 0;
    for (int i = 0; i < EffectType::TRAIL_LENGTH + 4; ++i) {
      fx.draw_frame();
      fx.advance_display();
      worst_vertices = std::max(worst_vertices, WB::tween_vertices(fx));
    }

    HS_EXPECT_LE(scratch_arena_a.get_high_water_mark(),
                 timed ? EffectType::TIMED_SCRATCH_A
                       : EffectType::STILL_SCRATCH_A);
    // The trail fills, so the peak above is a saturated frame and not a
    // warm-up.
    HS_EXPECT_GT(worst_vertices, (size_t)EffectType::TRAIL_LENGTH);
    HS_EXPECT_LE(worst_vertices, (size_t)EffectType::MAX_FRAGMENTS);
  }
}

/** @brief Fishbowl scratch bounds cover both roster resolutions. */
inline void test_fishbowl_scratch_estimate_covers_peak() {
  check_fishbowl_scratch_estimate_covers_peak<SMALL_W, SMALL_H>();
  check_fishbowl_scratch_estimate_covers_peak<DEFAULT_W, DEFAULT_H>();
}
