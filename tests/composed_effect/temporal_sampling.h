/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @file temporal_sampling.h
 * @brief Column-time interpolation and retained-resource tests.
 */

/** @brief Samples spatial state and clocks between consecutive frame endpoints.
 * @tparam Frame Composed temporal frame context.
 * @param from Previous frame, with retained palette and hue tables.
 * @param to Current frame.
 * @param fraction Fraction of one frame, in [0, 1].
 * @param constant_params Interior parameter sample when both endpoints agree.
 * @return Column context with resources borrowed from the endpoints. */
template <typename Frame>
HS_FLASH_INLINE Frame sample_frame_reference(
    const Frame &from, const Frame &to, float fraction,
    const decltype(Frame::params) *constant_params = nullptr) {
  if (fraction <= 0.0f) {
    Frame frame = from;
    Pullback::hold_discrete_lenses(frame, from, to);
    return frame;
  }
  if (fraction >= 1.0f)
    return to;
  Frame frame = to;
  frame.projection_conjugate =
      math::slerp(from.projection_conjugate, to.projection_conjugate, fraction);
  frame.outer_conjugate =
      math::slerp(from.outer_conjugate, to.outer_conjugate, fraction);
  frame.params = constant_params
                     ? *constant_params
                     : Pullback::interpolate(from.params, to.params, fraction);
  frame.palette_mapping = Pullback::Color::PaletteMappingWeights::lerp(
      from.palette_mapping, to.palette_mapping, fraction);
  frame.palette_oscillation_phase =
      Spacetime::sample_clock(from.palette_oscillation_phase,
                              to.palette_oscillation_phase, fraction, 1.0f);
  frame.params.visit([&]<typename Resource>(const auto &family) {
    auto &clock = frame.resources.template get<Resource::KEY>();
    const auto &a = from.resources.template get<Resource::KEY>();
    const auto &b = to.resources.template get<Resource::KEY>();
    if constexpr (Resource::KIND == Pullback::ResourceKind::SOURCE) {
      clock.primary = Spacetime::sample_clock(a.primary, b.primary, fraction,
                                              math::TWO_PI_F);
      clock.secondary = Spacetime::sample_clock(a.secondary, b.secondary,
                                                fraction, math::TWO_PI_F);
      clock.angle =
          Spacetime::sample_clock(a.angle, b.angle, fraction, math::TWO_PI_F);
      clock.noise_time =
          Spacetime::sample_clock(a.noise_time, b.noise_time, fraction, 1.0f);
    } else if constexpr (Resource::KIND == Pullback::ResourceKind::WARP) {
      clock.phase = Spacetime::sample_clock(a.phase, b.phase, fraction, 1.0f);
      if constexpr (std::is_same_v<typename Resource::Family,
                                   Pullback::AffineParams>)
        clock.rotation = Spacetime::sample_clock(a.rotation, b.rotation,
                                                 fraction, math::TWO_PI_F);
    } else if constexpr (Resource::KIND == Pullback::ResourceKind::SURFACE) {
      if constexpr (std::is_same_v<typename Resource::Family,
                                   Pullback::PeriodicRippleParams>)
        clock.phase =
            Spacetime::sample_clock(a.phase, b.phase, fraction, family.period);
      else
        clock.phase = Spacetime::sample_clock(a.phase, b.phase, fraction, 1.0f);
    } else if constexpr (std::is_same_v<typename Resource::Family,
                                        Pullback::MobiusLensParams>) {
      const auto &start = from.params.template get<Resource::KEY>();
      const auto &end = to.params.template get<Resource::KEY>();
      frame.params.template get<Resource::KEY>() =
          Pullback::lens_phase_continuous(start, end, a, b)
              ? Pullback::sample_lens_phase(
                    end,
                    Spacetime::sample_clock(a.angle, b.angle, fraction,
                                            math::TWO_PI_F),
                    a.magnitude + (b.magnitude - a.magnitude) * fraction)
              : end;
    }
  });
  frame.temporal = {from.palette,
                    from.hue_rotation_lut,
                    from.hue_noise_lut,
                    fraction,
                    from.hue_rotation_linear,
                    from.palette_alpha,
                    from.hue_noise_linear};
  return frame;
}

/** @brief Custom shading remains authoritative with column factoring enabled. */
class TemporalShadeProbe
    : public Pullback::ComposedEffect<96, 20, TemporalShadeProbe,
                                      GridSpaceSpec> {
public:
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"custom-shade"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr Params initial_params() {
    return GridSpace<96, 20>::initial_params();
  }
  static consteval auto choreography_hooks() {
    return Pullback::ComposedEffect<96, 20, TemporalShadeProbe,
                                    GridSpaceSpec>::choreography_hooks()
        .template extend<
            TemporalShadeProbe,
            Control::HookBinding<Control::ChoreographyHook::INITIAL_PARAMS,
                                 &TemporalShadeProbe::initial_params>>();
  }
  static Color4 shade(const math::Vector &, const Frame &) {
    return Color4(Pixel(1000, 2000, 3000), 0.5f);
  }
};

/** @brief Both arms use absolute column timing through clock rollover. */
inline void test_temporal_clock_sampling() {
  HS_CONTEXT("temporal clock sampling");
  for (int x = 0; x < 144; ++x) {
    HS_EXPECT_NEAR(Spacetime::SamplingSurface<288>::time(x), x / 144.0f, 0.0f);
    HS_EXPECT_NEAR(Spacetime::surface_time(x + 144, 288), x / 144.0f, 0.0f);
  }
  HS_EXPECT_NEAR(Spacetime::sample_clock(0.99f, 0.01f, 0.5f, 1.0f), 1.0f,
                 1e-6f);
  HS_EXPECT_NEAR(Spacetime::sample_clock(0.01f, 0.99f, 0.5f, 1.0f), 0.0f,
                 1e-6f);
  const float LAST = Spacetime::sample_clock(0.98f, 0.0f, 143.0f / 144, 1.0f);
  HS_EXPECT_NEAR(1.0f - LAST, 0.02f / 144, 1e-6f);

  using FX = GridSpace<96, 20>;
  typename FX::FrameState from{}, to{};
  from.params = FX::initial_params();
  to.params = from.params;
  from.resources.template get<"outer_warp">().phase = 0.99f;
  to.resources.template get<"outer_warp">().phase = 0.01f;
  from.resources.template get<"outer_warp">().rotation = 6.25f;
  to.resources.template get<"outer_warp">().rotation = 0.03f;
  to.outer_conjugate = math::make_rotation(math::Y_AXIS, 0.4f);
  const auto middle = Spacetime::sample_frame(from, to, 0.5f);
  HS_EXPECT_NEAR(middle.resources.template get<"outer_warp">().phase, 1.0f,
                 1e-6f);
  HS_EXPECT_NEAR(middle.resources.template get<"outer_warp">().rotation,
                 6.281592653589793f, 1e-6f);
  const math::Quaternion EXPECTED = math::make_rotation(math::Y_AXIS, 0.2f);
  HS_EXPECT_NEAR(math::dot(middle.outer_conjugate, EXPECTED), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(Spacetime::sample_frame(from, to, 0.0f)
                     .resources.template get<"outer_warp">()
                     .phase,
                 0.99f, 0.0f);
  HS_EXPECT_NEAR(Spacetime::sample_frame(from, to, 1.0f)
                     .resources.template get<"outer_warp">()
                     .phase,
                 0.01f, 0.0f);

  from.resources.template get<"source">().primary = 6.2f;
  to.resources.template get<"source">().primary = 0.1f;
  from.resources.template get<"source">().secondary = 0.1f;
  to.resources.template get<"source">().secondary = 6.2f;
  from.resources.template get<"source">().angle = 6.25f;
  to.resources.template get<"source">().angle = 0.03f;
  from.resources.template get<"source">().noise_time = 0.99f;
  to.resources.template get<"source">().noise_time = 0.01f;
  from.palette_oscillation_phase = 0.02f;
  to.palette_oscillation_phase = 0.98f;
  for (int variant = 0; variant < 3; ++variant) {
    if (variant == 1) {
      to.params.template get<"source">().lattice_cell_scale = 2.7f;
      to.params.template get<"outer_warp">().rotation_rate = 0.3f;
      to.params.template get<"color">().palette_mapping =
          Pullback::Color::PaletteMapping::CUP;
    }
    if (variant == 2)
      from.params = to.params;
    for (float angle : {0.0f, 0.0001f, 0.02f, 0.4f, 3.2f, 6.2f}) {
      to.projection_conjugate = math::make_rotation(math::X_AXIS, angle);
      to.outer_conjugate = -math::make_rotation(math::Y_AXIS, angle);
      const Pullback::FrameInterval sampler(from, to);
      for (int column = 0; column <= 144; ++column) {
        const float FRACTION = column / 144.0f;
        const auto actual = sampler.at(FRACTION);
        const auto expected = sample_frame_reference(from, to, FRACTION);
        const auto check_rotation = [](const math::Quaternion &a,
                                       const math::Quaternion &b) {
          HS_EXPECT_NEAR(a.r, b.r, 1e-7f);
          HS_EXPECT_NEAR(a.v.x, b.v.x, 1e-7f);
          HS_EXPECT_NEAR(a.v.y, b.v.y, 1e-7f);
          HS_EXPECT_NEAR(a.v.z, b.v.z, 1e-7f);
        };
        check_rotation(actual.projection_conjugate,
                       expected.projection_conjugate);
        check_rotation(actual.outer_conjugate, expected.outer_conjugate);
        HS_EXPECT_NEAR(actual.palette_oscillation_phase,
                       expected.palette_oscillation_phase, 0.0f);
        actual.params.visit([&]<typename Resource>(const auto &family) {
          const auto &expected_family =
              expected.params.template get<Resource::KEY>();
          for (const auto &field : Resource::Family::FIELDS)
            HS_EXPECT_NEAR(family.*field.member, expected_family.*field.member,
                           0.0f);
        });
        HS_EXPECT_NEAR(actual.resources.template get<"outer_warp">().rotation,
                       expected.resources.template get<"outer_warp">().rotation,
                       0.0f);
        if (column == 72)
          HS_EXPECT_NEAR(actual.resources.template get<"outer_warp">().rotation,
                         6.281592653589793f, 1e-6f);
        const auto &a = actual.resources.template get<"source">();
        const auto &b = expected.resources.template get<"source">();
        HS_EXPECT_NEAR(a.primary, b.primary, 0.0f);
        HS_EXPECT_NEAR(a.secondary, b.secondary, 0.0f);
        HS_EXPECT_NEAR(a.angle, b.angle, 0.0f);
        HS_EXPECT_NEAR(a.noise_time, b.noise_time, 0.0f);
        HS_EXPECT_NEAR(actual.resources.template get<"outer_warp">().phase,
                       expected.resources.template get<"outer_warp">().phase,
                       0.0f);
      }
    }
  }

  using SurfaceParams = Pullback::ComposedDetail::ParameterSet<
      Pullback::ComposedDetail::ResourceList<
          Pullback::ParameterResource<"surface", Pullback::SurfaceNoiseParams,
                                      Pullback::ResourceKind::SURFACE>>>;
  Pullback::FrameIntervalState<SurfaceParams> surface_from{}, surface_to{};
  surface_from.resources.template get<"surface">().phase = 0.99f;
  surface_to.resources.template get<"surface">().phase = 0.01f;
  const Pullback::FrameInterval surface_sampler(surface_from, surface_to);
  for (int column = 0; column <= 144; ++column) {
    const float FRACTION = column / 144.0f;
    const auto actual = surface_sampler.at(FRACTION);
    const auto expected =
        Spacetime::sample_frame(surface_from, surface_to, FRACTION);
    HS_EXPECT_NEAR(actual.resources.template get<"surface">().phase,
                   expected.resources.template get<"surface">().phase, 0.0f);
    if (column == 72) {
      HS_EXPECT_NEAR(actual.resources.template get<"surface">().phase, 1.0f,
                     1e-6f);
      HS_EXPECT_NEAR(expected.resources.template get<"surface">().phase, 1.0f,
                     1e-6f);
    }
  }

  using RippleParams = Pullback::ComposedDetail::ParameterSet<
      Pullback::ComposedDetail::ResourceList<
          Pullback::ParameterResource<"surface", Pullback::PeriodicRippleParams,
                                      Pullback::ResourceKind::SURFACE>>>;
  Pullback::FrameIntervalState<RippleParams> ripple_from{}, ripple_to{};
  ripple_from.params.template get<"surface">().period = 2.0f;
  ripple_from.resources.template get<"surface">().phase = 1.99f;
  ripple_to.resources.template get<"surface">().phase = 0.01f;
  for (float period : {2.0f, 3.0f}) {
    ripple_to.params.template get<"surface">().period = period;
    const Pullback::FrameInterval sampler(ripple_from, ripple_to);
    for (int column = 0; column <= 144; ++column) {
      const float FRACTION = column / 144.0f;
      const auto actual = sampler.at(FRACTION);
      const auto expected =
          Spacetime::sample_frame(ripple_from, ripple_to, FRACTION);
      HS_EXPECT_NEAR(actual.resources.template get<"surface">().phase,
                     expected.resources.template get<"surface">().phase, 0.0f);
    }
  }
}

/** @brief Affine scale follows the periodic warp cycle. */
inline void test_temporal_affine_cycle() {
  HS_CONTEXT("temporal affine cycle");
  for (float phase : {0.0f, 0.25f, 0.5f, 0.99f, 1.01f}) {
    Pullback::AffineParams warp;
    warp.translation_x = 2.0f;
    warp.translation_y = -3.0f;
    for (float scale : {1.0f, 0.7f, 1.8f}) {
      warp.scale_x = scale;
      const auto actual = Pullback::Warp::prepare(warp, phase, 0.3f, 1.2f);
      const float CYCLE = cosf(math::TWO_PI_F * math::wrap_t(phase));
      HS_EXPECT_NEAR(actual.transform.affine.scale_x, powf(scale, CYCLE), 0.0f);
      HS_EXPECT_NEAR(actual.transform.affine.scale_y, powf(1.0f, CYCLE), 0.0f);
      HS_EXPECT_NEAR(actual.transform.affine.shear, warp.shear * CYCLE, 0.0f);
    }
  }
}

/** @brief Nested sampling guards each restore the mode they entered with. */
inline void test_scoped_temporal_sampling_restores_entering_mode() {
  HS_CONTEXT("scoped temporal sampling restores entering mode");
  Spacetime::set_flat_surface(true);
  {
    const ScopedTemporalSampling outer(false);
    HS_EXPECT_TRUE(Spacetime::flat_surface());
    {
      const ScopedTemporalSampling inner(true);
      HS_EXPECT_FALSE(Spacetime::flat_surface());
    }
    HS_EXPECT_TRUE(Spacetime::flat_surface());
  }
  HS_EXPECT_TRUE(Spacetime::flat_surface());
  {
    const ScopedTemporalSampling column_basis(true);
    HS_EXPECT_FALSE(Spacetime::flat_surface());
  }
  HS_EXPECT_TRUE(Spacetime::flat_surface());
  Spacetime::set_flat_surface(false);
}

/** @brief Custom shading is authoritative with temporal sampling enabled. */
inline void test_temporal_custom_shade_is_authoritative() {
  HS_CONTEXT("temporal custom shade is authoritative");
  reset_globals();
  TemporalShadeProbe custom;
  custom.init();
  for (bool enabled : {true, false, true}) {
    Spacetime::set_flat_surface(!enabled);
    custom.draw_frame();
    custom.advance_display();
    HS_EXPECT_TRUE(custom.get_pixel(0, 0) == Pixel(500, 1000, 1500));
    HS_EXPECT_TRUE(custom.get_pixel(47, 19) == Pixel(500, 1000, 1500));
  }
}

/** @brief Endpoints with equal palette tables colour as the current one, and
 *  any differing entry interpolates. */
inline void test_temporal_equal_palettes_color_as_current() {
  HS_CONTEXT("temporal equal palettes");
  reset_globals();
  using FX = GridSpace<96, 20>;
  FX effect;
  effect.init();
  effect.draw_frame();
  effect.advance_display();
  const auto current = effect.frame_for_test();
  Spacetime::FrameHistory<FX::FrameState, Pullback::HueMode::NOISE, true, true>
      history;
  history.capture(current, persistent_arena);
  using PathPalette = Pullback::Color::GeneratedPalette<
      Pullback::ColorProvider<FX::Binding, Pullback::HueMode::PATH_LENGTH,
                              Pullback::Color::BrightnessEnvelope::CUP>>;
  const Pullback::FrameInterval EQUAL(history.previous, current);
  auto sample = EQUAL.at(0.5f);
  HS_EXPECT_TRUE(Spacetime::sample_frame(history.previous, current, 0.5f)
                     .temporal.same_palette);
  HS_EXPECT_TRUE(sample.temporal.same_palette);
  sample.params.template get<"color">().hue_shift_amount = 0.3f;
  const auto prepared = PathPalette::prepare(sample);
  HS_EXPECT_NEAR(prepared.fraction, 1.0f, 0.0f);
  auto interpolated = sample;
  interpolated.temporal.same_palette = false;
  const auto blended = PathPalette::prepare(interpolated);
  HS_EXPECT_NEAR(blended.fraction, 0.5f, 0.0f);
  for (int index = 0; index < 257; ++index) {
    Pullback::FieldSample field{};
    field.value = static_cast<float>(index) / 256;
    field.coverage = 0.75f;
    field.path_length = static_cast<float>(index) * 0.37f;
    const Color4 A = PathPalette::apply(field, sample, prepared);
    const Color4 B = PathPalette::apply(field, interpolated, blended);
    HS_EXPECT_TRUE(A.color == B.color);
    HS_EXPECT_NEAR(A.alpha, B.alpha, 0.0f);
  }

  auto rotation = history.hue_rotation;
  rotation[7].r ^= 1;
  auto rotated = current;
  rotated.hue_rotation_lut = rotation.data();
  HS_EXPECT_TRUE(!Pullback::FrameInterval(history.previous, rotated)
                      .at(0.5f)
                      .temporal.same_palette);
  BakedPaletteStorage other;
  other.bake(persistent_arena, Palettes::RICH_SUNSET);
  auto repainted = current;
  repainted.palette = &other.view();
  HS_EXPECT_TRUE(!Pullback::FrameInterval(history.previous, repainted)
                      .at(0.5f)
                      .temporal.same_palette);
}

/** @brief Column samples interpolate parameters as the whole set does when
 *  only some families move. */
inline void test_temporal_moving_families_interpolate() {
  HS_CONTEXT("temporal moving families");
  reset_globals();
  using FX = GridSpace<96, 20>;
  FX effect;
  effect.init();
  effect.draw_frame();
  effect.advance_display();
  const auto current = effect.frame_for_test();
  Spacetime::FrameHistory<FX::FrameState, Pullback::HueMode::NOISE, true, true>
      history;
  history.capture(current, persistent_arena);
  auto moved = current;
  moved.params.template get<"color">().hue_shift_amount += 0.25f;
  const Pullback::FrameInterval SAMPLER(history.previous, moved);
  for (const float fraction : {0.1f, 0.37f, 0.5f, 0.9f}) {
    HS_CONTEXT("fraction", static_cast<double>(fraction));
    const auto sampled = SAMPLER.at(fraction).params;
    const auto expected =
        Pullback::interpolate(history.previous.params, moved.params, fraction);
    HS_EXPECT_NEAR(sampled.template get<"color">().hue_shift_amount,
                   expected.template get<"color">().hue_shift_amount, 0.0f);
    HS_EXPECT_TRUE(std::memcmp(&sampled, &expected, sizeof(sampled)) == 0);
  }
}

/** @brief Retained color tables survive subsequent bakes. */
inline void test_temporal_owned_resources() {
  HS_CONTEXT("temporal owned color resources");
  reset_globals();
  using FX = GridSpace<96, 20>;
  FX effect;
  effect.init();
  effect.draw_frame();
  effect.advance_display();
  const auto current = effect.frame_for_test();
  Spacetime::FrameHistory<FX::FrameState, Pullback::HueMode::NOISE, true, true>
      history;
  history.capture(current, persistent_arena);
  HS_EXPECT_TRUE(history.previous.palette != current.palette);
  HS_EXPECT_TRUE(history.previous.hue_noise_lut != current.hue_noise_lut);
  HS_EXPECT_TRUE(history.previous.hue_rotation_lut != current.hue_rotation_lut);
  HS_EXPECT_TRUE(current.hue_rotation_linear != nullptr);
  HS_EXPECT_TRUE(history.previous.hue_rotation_linear !=
                 current.hue_rotation_linear);
  const Pixel PALETTE = history.previous.palette->get(0.4f).color;
  const auto NOISE = history.hue_noise;
  const auto ROTATION = history.hue_rotation;
  const auto LINEAR = history.hue_rotation_linear;
  const auto LINEAR_NOISE = history.hue_noise_linear;
  auto equal_frame = current;
  equal_frame.hue_noise_lut = NOISE.data();
  equal_frame.hue_noise_linear = LINEAR_NOISE.data();
  const Pullback::FrameInterval EQUAL_SAMPLER(history.previous, equal_frame);
  auto equal_sample = EQUAL_SAMPLER.at(0.5f);
  HS_EXPECT_TRUE(equal_sample.temporal.same_noise);
  auto different_noise = NOISE;
  different_noise[0] ^= 1;
  auto different_frame = equal_frame;
  different_frame.hue_noise_lut = different_noise.data();
  const Pullback::FrameInterval DIFFERENT_SAMPLER(history.previous,
                                                  different_frame);
  HS_EXPECT_TRUE(!DIFFERENT_SAMPLER.at(0.5f).temporal.same_noise);
  different_frame.hue_noise_lut = nullptr;
  const Pullback::FrameInterval NO_NOISE_SAMPLER(history.previous,
                                                 different_frame);
  HS_EXPECT_TRUE(!NO_NOISE_SAMPLER.at(0.5f).temporal.same_noise);
  using LinearPalette = Pullback::Color::GeneratedPalette<
      Pullback::ColorProvider<FX::Binding, Pullback::HueMode::NOISE,
                              Pullback::Color::BrightnessEnvelope::NONE>>;
  using ShapedPalette = Pullback::Color::GeneratedPalette<
      Pullback::ColorProvider<FX::Binding, Pullback::HueMode::NOISE,
                              Pullback::Color::BrightnessEnvelope::CUP>>;
  const auto verify_reuse = [&]<typename Palette>() {
    auto prepared = Palette::prepare(equal_sample);
    HS_EXPECT_TRUE(prepared.same_noise);
    constexpr bool LINEAR = std::is_same_v<Palette, LinearPalette>;
    HS_EXPECT_TRUE(prepared.linear_noise == LINEAR);
    auto separate = prepared;
    separate.same_noise = false;
    for (int index = 0; index < 257; ++index) {
      Pullback::FieldSample field{};
      field.value = static_cast<float>(index) / 256;
      field.coverage = 0.75f;
      field.sphere = math::Vector(sinf(index * 0.731f), cosf(index * 0.471f),
                                  sinf(index * 0.927f));
      const Color4 A = Palette::apply(field, equal_sample, prepared);
      const Color4 B = Palette::apply(field, equal_sample, separate);
      HS_EXPECT_TRUE(A.color == B.color);
      HS_EXPECT_NEAR(A.alpha, B.alpha, 0.0f);
      if (prepared.bounded_noise) {
        auto rounded = prepared;
        rounded.bounded_noise = false;
        const Color4 REFERENCE = Palette::apply(field, equal_sample, rounded);
        HS_EXPECT_NEAR(A.color.r, REFERENCE.color.r, 2.0f);
        HS_EXPECT_NEAR(A.color.g, REFERENCE.color.g, 2.0f);
        HS_EXPECT_NEAR(A.color.b, REFERENCE.color.b, 2.0f);
        HS_EXPECT_NEAR(A.alpha, REFERENCE.alpha, 0.0f);
      }
    }
  };
  verify_reuse.template operator()<LinearPalette>();
  equal_sample.params.template get<"color">().brightness_bottom = 0.4f;
  equal_sample.params.template get<"color">().brightness_top = 1.3f;
  HS_EXPECT_TRUE(!ShapedPalette::prepare(equal_sample).bounded_noise);
  verify_reuse.template operator()<ShapedPalette>();
  equal_sample.params.template get<"color">().brightness_bottom = 0.345f;
  equal_sample.params.template get<"color">().brightness_top = 1.0f;
  HS_EXPECT_TRUE(ShapedPalette::prepare(equal_sample).bounded_noise);
  verify_reuse.template operator()<ShapedPalette>();
  equal_sample.hue_rotation_linear = nullptr;
  equal_sample.hue_noise_linear = nullptr;
  equal_sample.temporal.hue_rotation_linear = nullptr;
  equal_sample.temporal.hue_noise_linear = nullptr;
  verify_reuse.template operator()<LinearPalette>();
  verify_reuse.template operator()<ShapedPalette>();
  equal_sample.temporal.fraction = 0.0f;
  HS_EXPECT_TRUE(!ShapedPalette::prepare(equal_sample).bounded_noise);
  equal_sample.temporal.fraction = 1.0f;
  HS_EXPECT_TRUE(!ShapedPalette::prepare(equal_sample).bounded_noise);
  HS_EXPECT_TRUE(current.hue_noise_linear != nullptr);
  HS_EXPECT_TRUE(history.previous.hue_noise_linear != current.hue_noise_linear);
  for (size_t index = 0; index < LINEAR_NOISE.size(); ++index)
    HS_EXPECT_NEAR(LINEAR_NOISE[index], NOISE[index], 0.0f);
  HS_EXPECT_NEAR(history.previous.palette_alpha, current.palette_alpha, 0.0f);
  HS_EXPECT_NEAR(current.palette_alpha, current.palette->get_alpha(0.4f), 0.0f);
  for (size_t index = 0; index < LINEAR.size(); ++index) {
    HS_EXPECT_NEAR(LINEAR[index].x, ROTATION[index].r, 0.0f);
    HS_EXPECT_NEAR(LINEAR[index].y, ROTATION[index].g, 0.0f);
    HS_EXPECT_NEAR(LINEAR[index].z, ROTATION[index].b, 0.0f);
  }
  const Pullback::Color::HueRotationLutView UNCACHED{ROTATION.data(), true};
  const Pullback::Color::HueRotationLutView CACHED{ROTATION.data(), true,
                                                   LINEAR.data()};
  const HueNoiseLutView NOISE_UNCACHED{NOISE.data(), true};
  const HueNoiseLutView NOISE_CACHED{NOISE.data(), true, LINEAR_NOISE.data()};
#if defined(HS_TEST_FAST_MATH)
  // Fast-math reassociates each inlined bilinear lerp independently.
  constexpr float NOISE_TOLERANCE = 4.0f * FLT_EPSILON;
#else
  constexpr float NOISE_TOLERANCE = 0.0f;
#endif
  for (int index = 0; index < 4096; ++index) {
    const HueRotationValueSample VALUE(static_cast<float>(index % 257) / 256);
    const float HUE = static_cast<float>(index - 2048) / 511;
    const auto A = sample_hue_rotation_linear(UNCACHED, VALUE, HUE);
    const auto B = sample_hue_rotation_linear(CACHED, VALUE, HUE);
    HS_EXPECT_NEAR(A.x, B.x, 0.0f);
    HS_EXPECT_NEAR(A.y, B.y, 0.0f);
    HS_EXPECT_NEAR(A.z, B.z, 0.0f);
    const math::Vector DIRECTION(sinf(index * 0.731f), cosf(index * 0.471f),
                                 sinf(index * 0.927f));
    const auto ADDRESS = prepare_hue_noise_sample(DIRECTION);
    HS_EXPECT_NEAR(ADDRESS.sample(NOISE_UNCACHED), ADDRESS.sample(NOISE_CACHED),
                   NOISE_TOLERANCE);
  }
  for (int frame = 0; frame < 4; ++frame) {
    effect.draw_frame();
    effect.advance_display();
  }
  HS_EXPECT_TRUE(history.previous.palette->get(0.4f).color == PALETTE);
  HS_EXPECT_TRUE(history.hue_noise == NOISE);
  HS_EXPECT_TRUE(history.hue_rotation == ROTATION);
  HS_EXPECT_TRUE(std::memcmp(history.hue_rotation_linear.data(), LINEAR.data(),
                             sizeof(LINEAR)) == 0);
  HS_EXPECT_TRUE(std::memcmp(history.hue_noise_linear.data(),
                             LINEAR_NOISE.data(), sizeof(LINEAR_NOISE)) == 0);
}

/** @brief Factored projection stays within the RGB16 rounding budget. */
inline void test_gridspace_factored_columns() {
  HS_CONTEXT("GridSpace factored columns");
  reset_globals();
  using FX = GridSpace<288, 144>;
  using Pipeline = FX::RenderPipeline;
  static_assert(Pipeline::LINEAR_COLUMNS);
  using GroupedPrefix = Pullback::Detail::TypeList<
      Pipeline::StageAt<0>,
      Pullback::Stage::Placed<Pullback::CodeEmission::INLINE_ONLY,
                              Pipeline::StageAt<0>, Pipeline::StageAt<1>>,
      Pipeline::StageAt<2>>;
  static_assert(!Pullback::Detail::ColumnPrefix<GroupedPrefix>::SUPPORTED);
  FX effect;
  effect.init();
  effect.draw_frame();
  effect.advance_display();
  Spacetime::FrameHistory<FX::FrameState, Pullback::HueMode::NOISE, true>
      history;
  int maximum_error = 0;
  uint64_t error_sum = 0;
  uint64_t channels = 0;
  for (int step = 0; step < 8; ++step) {
    history.capture(effect.frame_for_test(), persistent_arena);
    effect.draw_frame();
    effect.advance_display();
    const auto current = effect.frame_for_test();
    for (int x = 0; x < 288; ++x) {
      const float TIME = step == 0   ? 0.0f
                         : step == 1 ? 1.0f
                                     : static_cast<float>(x % 144) / 144;
      const auto frame = Pipeline::prepare(
          Spacetime::sample_frame(history.previous, current, TIME));
      const float COSINE = math::TrigLUT<288, 144>::cos_theta(x);
      const float SINE = math::TrigLUT<288, 144>::sin_theta[x];
      const auto basis = Pipeline::column_basis(frame, COSINE, SINE);
      for (int y = 0; y < 144; ++y) {
        const float SIN_PHI = math::TrigLUT<288, 144>::sin_phi[y];
        const float COS_PHI = math::TrigLUT<288, 144>::cos_phi[y];
        const math::Vector VIEW(SIN_PHI * COSINE, COS_PHI, SIN_PHI * SINE);
        const Color4 EXPECTED = Pipeline::shade(VIEW, frame);
        const Color4 ACTUAL =
            Pipeline::shade_column(SIN_PHI, COS_PHI, frame, basis);
        const Pixel A = ACTUAL.color * ACTUAL.alpha;
        const Pixel B = EXPECTED.color * EXPECTED.alpha;
        for (int error : {std::abs(int(A.r) - B.r), std::abs(int(A.g) - B.g),
                          std::abs(int(A.b) - B.b)}) {
          maximum_error = std::max(maximum_error, error);
          error_sum += error;
          ++channels;
        }
      }
    }
  }
  HS_EXPECT_TRUE(maximum_error <= 16);
  HS_EXPECT_TRUE(static_cast<double>(error_sum) / channels < 0.05);
}

/** @brief GridSpace shades each clipped column at its own sweep time. */
inline void test_gridspace_temporal_pixels() {
  HS_CONTEXT("GridSpace column-time pixels");
  reset_globals();
  constexpr int W = 96;
  constexpr int H = 20;
  using FX = GridSpace<W, H>;
  FX effect;
  effect.init();
  effect.draw_frame();
  effect.advance_display();
  Spacetime::FrameHistory<FX::FrameState, Pullback::HueMode::NOISE, true>
      history;
  history.capture(effect.frame_for_test(), persistent_arena);
  effect.set_clip(3, 17, 39, 83);
  effect.draw_frame();
  effect.advance_display();
  const auto current = effect.frame_for_test();
  const auto frozen = FX::RenderPipeline::prepare(current);
  int changed = 0;
  const auto expect_pixel = [](const Pixel &actual, const Pixel &expected) {
    HS_EXPECT_NEAR(actual.r, expected.r, 16);
    HS_EXPECT_NEAR(actual.g, expected.g, 16);
    HS_EXPECT_NEAR(actual.b, expected.b, 16);
  };
  for (int x = 39; x < 83; ++x) {
    const float TIME = static_cast<float>(x % (W / 2)) / (W / 2);
    const auto frame = FX::RenderPipeline::prepare(
        Spacetime::sample_frame(history.previous, current, TIME));
    for (int y = 3; y < 17; ++y) {
      const float SIN_PHI = math::TrigLUT<W, H>::sin_phi[y];
      const math::Vector VIEW(SIN_PHI * math::TrigLUT<W, H>::cos_theta(x),
                              math::TrigLUT<W, H>::cos_phi[y],
                              SIN_PHI * math::TrigLUT<W, H>::sin_theta[x]);
      const Color4 EXPECTED = FX::shade(VIEW, frame);
      const Pixel PIXEL = EXPECTED.color * EXPECTED.alpha;
      expect_pixel(effect.get_pixel(x, y), PIXEL);
      const Color4 CURRENT = FX::shade(VIEW, frozen);
      if (PIXEL != CURRENT.color * CURRENT.alpha)
        ++changed;
    }
  }
  HS_EXPECT_TRUE(changed > 0);

  Spacetime::set_flat_surface(true);
  effect.draw_frame();
  effect.advance_display();
  const auto instantaneous =
      FX::RenderPipeline::prepare(effect.frame_for_test());
  for (int x = 39; x < 83; ++x) {
    for (int y = 3; y < 17; ++y) {
      const float SIN_PHI = math::TrigLUT<W, H>::sin_phi[y];
      const math::Vector VIEW(SIN_PHI * math::TrigLUT<W, H>::cos_theta(x),
                              math::TrigLUT<W, H>::cos_phi[y],
                              SIN_PHI * math::TrigLUT<W, H>::sin_theta[x]);
      const Color4 EXPECTED = FX::shade(VIEW, instantaneous);
      expect_pixel(effect.get_pixel(x, y), EXPECTED.color * EXPECTED.alpha);
    }
  }

  history.capture(effect.frame_for_test(), persistent_arena);
  Spacetime::set_flat_surface(false);
  effect.draw_frame();
  effect.advance_display();
  const auto resumed = effect.frame_for_test();
  for (int x = 39; x < 83; ++x) {
    const auto frame = FX::RenderPipeline::prepare(Spacetime::sample_frame(
        history.previous, resumed, Spacetime::SamplingSurface<W>::time(x)));
    for (int y = 3; y < 17; ++y) {
      const float SIN_PHI = math::TrigLUT<W, H>::sin_phi[y];
      const math::Vector VIEW(SIN_PHI * math::TrigLUT<W, H>::cos_theta(x),
                              math::TrigLUT<W, H>::cos_phi[y],
                              SIN_PHI * math::TrigLUT<W, H>::sin_theta[x]);
      const Color4 EXPECTED = FX::shade(VIEW, frame);
      expect_pixel(effect.get_pixel(x, y), EXPECTED.color * EXPECTED.alpha);
    }
  }
}

/**
 * @brief Authored effects retain endpoints across clipped draws and runtime toggles.
 * @tparam FX Composed effect specialization.
 * @param name Effect name for failure context.
 */
template <typename FX>
inline void check_temporal_roster_pixels(const char *name) {
  HS_CONTEXT(name);
  HS_EXPECT_TRUE(FX::TEMPORAL_SAMPLING);
  reset_globals();
  constexpr int W = 96;
  constexpr int H = 20;
  constexpr int TOLERANCE = FX::Spec::FACTORED_COLUMNS ? 16 : 2;
  FX effect;
  effect.init();
  Spacetime::FrameHistory<typename FX::FrameState, FX::Metadata::HUE, true,
                          FX::FLOAT_HUE_TABLES>
      history;
  const auto view = [](int x, int y) {
    const float SIN_PHI = math::TrigLUT<W, H>::sin_phi[y];
    return math::Vector(SIN_PHI * math::TrigLUT<W, H>::cos_theta(x),
                        math::TrigLUT<W, H>::cos_phi[y],
                        SIN_PHI * math::TrigLUT<W, H>::sin_theta[x]);
  };
  const auto expect_pixel = [](const Pixel &actual, const Color4 &expected) {
    const Pixel PIXEL = expected.color * expected.alpha;
    HS_EXPECT_NEAR(actual.r, PIXEL.r, TOLERANCE);
    HS_EXPECT_NEAR(actual.g, PIXEL.g, TOLERANCE);
    HS_EXPECT_NEAR(actual.b, PIXEL.b, TOLERANCE);
  };

  effect.draw_frame();
  effect.advance_display();
  const auto initial = effect.frame_for_test();
  const auto first = FX::RenderPipeline::prepare(initial);
  for (int x : {0, 1, 47, 48, 49, 95})
    for (int y : {1, 7, 18})
      expect_pixel(effect.get_pixel(x, y), FX::shade(view(x, y), first));
  history.capture(initial, persistent_arena);
  HS_EXPECT_TRUE(history.previous.palette != initial.palette);
  if constexpr (FX::Metadata::HUE == Pullback::HueMode::NOISE) {
    HS_EXPECT_TRUE(history.previous.hue_noise_lut != initial.hue_noise_lut);
    HS_EXPECT_TRUE(history.previous.hue_rotation_lut !=
                   initial.hue_rotation_lut);
  }

  effect.set_clip(3, 17, 39, 83);
  for (bool enabled : {true, false, false, true}) {
    Spacetime::set_flat_surface(!enabled);
    effect.draw_frame();
    effect.advance_display();
    const auto current = effect.frame_for_test();
    const auto instantaneous = FX::RenderPipeline::prepare(current);
    const auto previous = FX::RenderPipeline::prepare(history.previous);
    const auto start = FX::RenderPipeline::prepare(
        Spacetime::sample_frame(history.previous, current, 0.0f));
    const auto end = FX::RenderPipeline::prepare(
        Spacetime::sample_frame(history.previous, current, 1.0f));
    for (int x : {0, 47, 48, 95}) {
      const auto direction = view(x, 7);
      const Color4 BEFORE = FX::shade(direction, previous);
      const Color4 AFTER = FX::shade(direction, instantaneous);
      expect_pixel(BEFORE.color * BEFORE.alpha, FX::shade(direction, start));
      expect_pixel(AFTER.color * AFTER.alpha, FX::shade(direction, end));
    }
    int changed = 0;
    for (int x = 39; x < 83; ++x) {
      const float FRACTION = static_cast<float>(x % (W / 2)) / (W / 2);
      const auto frame = FX::RenderPipeline::prepare(
          Spacetime::sample_frame(history.previous, current, FRACTION));
      for (int y = 3; y < 17; ++y) {
        const auto direction = view(x, y);
        const Color4 TEMPORAL = FX::shade(direction, frame);
        const Color4 CURRENT = FX::shade(direction, instantaneous);
        expect_pixel(effect.get_pixel(x, y), enabled ? TEMPORAL : CURRENT);
        if (TEMPORAL.color * TEMPORAL.alpha != CURRENT.color * CURRENT.alpha)
          ++changed;
      }
    }
    HS_EXPECT_TRUE(changed > 0);
    history.capture(current, persistent_arena);
  }
  Spacetime::set_flat_surface(false);
}

/** @brief The compatible composed roster uses column-time rendering. */
inline void test_composed_temporal_roster() {
#define HS_COMPOSED_TEMPORAL(name, seconds)                                    \
  check_temporal_roster_pixels<name<96, 20>>(#name);
  HS_SHADER_PRODUCT_GROUP(HS_COMPOSED_TEMPORAL)
#undef HS_COMPOSED_TEMPORAL
}

/** @brief MobiusGrid with temporal sampling disabled matches its committed
 *  render through a Lerp departure, clips and a manual preset. */
inline void test_mobius_grid_instantaneous_pixel_golden() {
  using FX = MobiusGrid<SMALL_W, SMALL_H>;
  const ScopedTemporalSampling instantaneous(false);
  reset_effect_globals();
  FX effect;
  effect.init();
  const int departure = FX::PRESET_DWELL_FRAMES +
                        Segue::Preset::frames(FX::preset_departure(0)) + 20;
  uint64_t hash = HASH_SEED64;
  for (int frame = 0; frame < departure + 40; ++frame) {
    pin_frame_clock(frame);
    if (frame == 30)
      effect.set_clip(2, SMALL_H - 2, 40, 60);
    if (frame == 40) {
      effect.set_clip(0, SMALL_H, 0, 8);
      effect.set_margin(4);
    }
    if (frame == 50) {
      effect.set_margin(0);
      effect.set_clip(0, SMALL_H, 0, SMALL_W);
    }
    if (frame == departure)
      HS_EXPECT_TRUE(effect.selectPreset(0));
    if (frame == departure + 20)
      effect.setAnimationsPaused(false);
    effect.draw_frame();
    effect.advance_display();
    hash = frame_hash<SMALL_W, SMALL_H>(effect, hash);
  }
  hs::clear_mock_time();
  expect_pixel_golden("MobiusGrid 96x20", hash, 0x7aba7e31fa280f26ull);
}

/** @brief Composed history type of an effect.
 * @tparam FX Composed effect specialization. */
template <typename FX>
using ComposedHistory =
    Spacetime::FrameHistory<typename FX::FrameState, FX::Metadata::HUE, true,
                            FX::FLOAT_HUE_TABLES>;

/** @brief MobiusGrid's lens b at a continuous frame time.
 * @param time Frames since init; draw k spans [k, k + 1].
 * @return The circular warp's b at the eased phase of @p time. */
inline math::Complex mobius_grid_b_at(double time) {
  constexpr double PERIOD = 160.0;
  const double cycle = std::fmod(time, PERIOD) / PERIOD;
  const float angle =
      static_cast<float>(2.0 * 3.14159265358979323846 *
                         math::ease_in_out_sin(static_cast<float>(cycle)));
  return {std::cos(angle), -std::sin(angle)};
}

/** @brief MobiusGrid's lens follows its eased animation phase at every
 *  column, across the frame seam and the phase cycle wrap. */
inline void test_mobius_grid_lens_follows_phase() {
  using FX = MobiusGrid<SMALL_W, SMALL_H>;
  reset_globals();
  FX effect;
  effect.init();
  ComposedHistory<FX> history;
  double worst = 0.0;
  double worst_endpoint = 0.0;
  double held = 0.0;
  for (int draw = 0; draw < 165; ++draw) {
    pin_frame_clock(draw);
    effect.draw_frame();
    effect.advance_display();
    const auto current = effect.frame_for_test();
    const math::Complex endpoint =
        current.params.template get<"lens">().mobius.b;
    const math::Complex expected_endpoint = mobius_grid_b_at(draw + 1.0);
    worst_endpoint = std::max<double>(
        worst_endpoint, std::hypot(endpoint.re - expected_endpoint.re,
                                   endpoint.im - expected_endpoint.im));
    if (history.initialized &&
        (draw < 4 || (draw > 77 && draw < 82) || draw > 155)) {
      const Pullback::FrameInterval sampler(history.previous, current);
      for (int x = 0; x < SMALL_W; ++x) {
        const float fraction = Spacetime::SamplingSurface<SMALL_W>::time(x);
        const auto sampled = sampler.at(fraction);
        const auto reference =
            Spacetime::sample_frame(history.previous, current, fraction);
        const math::MobiusParams &lens =
            sampled.params.template get<"lens">().mobius;
        verify_mobius_equal(lens,
                            reference.params.template get<"lens">().mobius);
        const math::Complex expected = mobius_grid_b_at(draw + fraction);
        worst = std::max<double>(worst, std::hypot(lens.b.re - expected.re,
                                                   lens.b.im - expected.im));
        held = std::max<double>(held, std::hypot(endpoint.re - expected.re,
                                                 endpoint.im - expected.im));
        HS_EXPECT_EQ(lens.a.re,
                     current.params.template get<"lens">().mobius.a.re);
        HS_EXPECT_EQ(lens.d.re,
                     current.params.template get<"lens">().mobius.d.re);
      }
    }
    history.capture(current, persistent_arena);
  }
  hs::clear_mock_time();
  std::printf(
      "  MobiusGrid lens phase: worst %.3g, endpoints %.3g, held %.3g\n", worst,
      worst_endpoint, held);
  HS_EXPECT_LT(worst_endpoint, 1e-5);
  HS_EXPECT_LT(worst, 5e-4);
  HS_EXPECT_GT(held, 0.02);
}

/** @brief Unsupported lens transitions take the current lens in every
 *  column: a snapped b, changed a, c or d, and a degenerate path sample. */
inline void test_mobius_grid_lens_unsupported_transitions() {
  using FX = MobiusGrid<SMALL_W, SMALL_H>;
  reset_globals();
  FX effect;
  effect.init();
  for (int draw = 0; draw < 10; ++draw) {
    pin_frame_clock(draw);
    effect.draw_frame();
    effect.advance_display();
  }
  const auto before = effect.frame_for_test();
  HS_EXPECT_TRUE(before.resources.template get<"lens">().generated);
  const auto lens_at = [](const auto &from, const auto &to, float fraction) {
    return Pullback::FrameInterval(from, to)
        .at(fraction)
        .params.template get<"lens">()
        .mobius;
  };

  HS_EXPECT_TRUE(effect.selectPreset(1));
  const auto snapped = effect.frame_for_test();
  HS_EXPECT_FALSE(snapped.resources.template get<"lens">().generated);
  for (const float fraction : {0.0f, 0.5f}) {
    verify_mobius_equal(lens_at(before, snapped, fraction),
                        snapped.params.template get<"lens">().mobius);
    verify_mobius_equal(Spacetime::sample_frame(before, snapped, fraction)
                            .params.template get<"lens">()
                            .mobius,
                        snapped.params.template get<"lens">().mobius);
  }

  auto moved = before;
  auto &moved_phase = moved.resources.template get<"lens">();
  moved_phase.angle += 0.2f;
  Animation::MobiusWarpCircular::write_b(
      moved.params.template get<"lens">().mobius, moved_phase.magnitude,
      moved_phase.angle);
  HS_EXPECT_TRUE(Pullback::lens_generated(moved.params.template get<"lens">(),
                                          moved_phase));
  const math::MobiusParams halfway = lens_at(before, moved, 0.5f);
  HS_EXPECT_NEAR(
      std::atan2(-halfway.b.im, halfway.b.re),
      std::remainder(before.resources.template get<"lens">().angle + 0.1f,
                     math::TWO_PI_F),
      1e-5f);
  for (auto field : {&math::MobiusParams::a, &math::MobiusParams::c,
                     &math::MobiusParams::d}) {
    auto edited = moved;
    (edited.params.template get<"lens">().mobius.*field).im += 0.01f;
    verify_mobius_equal(lens_at(before, edited, 0.5f),
                        edited.params.template get<"lens">().mobius);
  }

  auto degenerate_from = before;
  auto degenerate_to = before;
  for (auto *frame : {&degenerate_from, &degenerate_to}) {
    math::MobiusParams &mobius = frame->params.template get<"lens">().mobius;
    mobius.a = {1.0f, 0.0f};
    mobius.c = {1.0f, 0.0f};
    mobius.d = {1.0f, 0.0f};
  }
  degenerate_from.resources.template get<"lens">().angle = -0.1f;
  degenerate_to.resources.template get<"lens">().angle = 0.1f;
  for (auto *frame : {&degenerate_from, &degenerate_to}) {
    auto &phase = frame->resources.template get<"lens">();
    phase.magnitude = 1.0f;
    Animation::MobiusWarpCircular::write_b(
        frame->params.template get<"lens">().mobius, 1.0f, phase.angle);
    HS_EXPECT_TRUE(Pullback::valid(frame->params.template get<"lens">()));
  }
  verify_mobius_equal(lens_at(degenerate_from, degenerate_to, 0.5f),
                      degenerate_to.params.template get<"lens">().mobius);
  const math::MobiusParams quarter =
      lens_at(degenerate_from, degenerate_to, 0.25f);
  HS_EXPECT_NEAR(quarter.b.re, std::cos(-0.05f), 1e-5f);
  HS_EXPECT_NEAR(quarter.b.im, -std::sin(-0.05f), 1e-5f);
}

/** @brief Segmented MobiusGrid engines, including clips crossing both arm
 *  boundaries, draw the same pixels as one full-frame engine. */
inline void test_mobius_grid_segmented_matches_single_engine() {
  using FX = MobiusGrid<SMALL_W, SMALL_H>;
  struct Segment {
    int x0, x1, margin;
  };
  constexpr std::array<Segment, 4> SEGMENTS{
      {{0, 40, 0}, {40, 60, 0}, {60, SMALL_W, 0}, {0, 8, 4}}};
  constexpr int DRAWS = 24;
  const auto render = [](const Segment &segment, auto &&inspect) {
    reset_globals();
    FX engine;
    engine.init();
    engine.set_clip(0, SMALL_H, segment.x0, segment.x1);
    engine.set_margin(segment.margin);
    for (int draw = 0; draw < DRAWS; ++draw) {
      pin_frame_clock(draw);
      engine.draw_frame();
      engine.advance_display();
      inspect(engine, draw);
    }
    hs::clear_mock_time();
  };
  std::vector<std::vector<Pixel>> single(DRAWS);
  render({0, SMALL_W, 0}, [&](const FX &engine, int draw) {
    capture_frame<SMALL_W, SMALL_H>(engine, single[draw]);
  });
  int mismatches = 0;
  int compared = 0;
  for (const Segment &segment : SEGMENTS)
    render(segment, [&](const FX &engine, int draw) {
      for (int column = segment.x0 - segment.margin;
           column < segment.x1 + segment.margin; ++column) {
        const int x = (column + SMALL_W) % SMALL_W;
        for (int y = 0; y < SMALL_H; ++y) {
          ++compared;
          mismatches += engine.get_pixel(x, y) != single[draw][y * SMALL_W + x];
        }
      }
    });
  HS_EXPECT_EQ(compared, DRAWS * (SMALL_W + 16) * SMALL_H);
  HS_EXPECT_EQ(mismatches, 0);
}
