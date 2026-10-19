/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

/** @file temporal_sampling.h
 * @brief Column-time interpolation and retained-resource tests.
 */

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
    HS_EXPECT_NEAR(Temporal::sweep_fraction(x, 288), x / 144.0f, 0.0f);
    HS_EXPECT_NEAR(Temporal::sweep_fraction(x + 144, 288), x / 144.0f, 0.0f);
  }
  HS_EXPECT_NEAR(Pullback::sample_clock(0.99f, 0.01f, 0.5f, 1.0f), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(Pullback::sample_clock(0.01f, 0.99f, 0.5f, 1.0f), 0.0f, 1e-6f);
  const float LAST = Pullback::sample_clock(0.98f, 0.0f, 143.0f / 144, 1.0f);
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
  const auto middle = Pullback::sample_frame(from, to, 0.5f);
  HS_EXPECT_NEAR(middle.resources.template get<"outer_warp">().phase, 1.0f,
                 1e-6f);
  HS_EXPECT_NEAR(middle.resources.template get<"outer_warp">().rotation,
                 6.281592653589793f, 1e-6f);
  const math::Quaternion EXPECTED = math::make_rotation(math::Y_AXIS, 0.2f);
  HS_EXPECT_NEAR(math::dot(middle.outer_conjugate, EXPECTED), 1.0f, 1e-6f);
  HS_EXPECT_NEAR(Pullback::sample_frame(from, to, 0.0f)
                     .resources.template get<"outer_warp">()
                     .phase,
                 0.99f, 0.0f);
  HS_EXPECT_NEAR(Pullback::sample_frame(from, to, 1.0f)
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
      const Pullback::TemporalFrameSampler sampler(from, to);
      for (int column = 0; column <= 144; ++column) {
        const float FRACTION = column / 144.0f;
        const auto actual = sampler.sample(FRACTION);
        const auto expected = Pullback::sample_frame(from, to, FRACTION);
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
  Pullback::TemporalFrameState<SurfaceParams> surface_from{}, surface_to{};
  surface_from.resources.template get<"surface">().phase = 0.99f;
  surface_to.resources.template get<"surface">().phase = 0.01f;
  const Pullback::TemporalFrameSampler surface_sampler(surface_from,
                                                       surface_to);
  for (int column = 0; column <= 144; ++column) {
    const float FRACTION = column / 144.0f;
    const auto actual = surface_sampler.sample(FRACTION);
    const auto expected =
        Pullback::sample_frame(surface_from, surface_to, FRACTION);
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
  Pullback::TemporalFrameState<RippleParams> ripple_from{}, ripple_to{};
  ripple_from.params.template get<"surface">().period = 2.0f;
  ripple_from.resources.template get<"surface">().phase = 1.99f;
  ripple_to.resources.template get<"surface">().phase = 0.01f;
  for (float period : {2.0f, 3.0f}) {
    ripple_to.params.template get<"surface">().period = period;
    const Pullback::TemporalFrameSampler sampler(ripple_from, ripple_to);
    for (int column = 0; column <= 144; ++column) {
      const float FRACTION = column / 144.0f;
      const auto actual = sampler.sample(FRACTION);
      const auto expected =
          Pullback::sample_frame(ripple_from, ripple_to, FRACTION);
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

/** @brief Custom shading is authoritative with temporal sampling enabled. */
inline void test_temporal_custom_shade_is_authoritative() {
  HS_CONTEXT("temporal custom shade is authoritative");
  reset_globals();
  TemporalShadeProbe custom;
  custom.init();
  for (bool enabled : {true, false, true}) {
    Temporal::set_temporal_sampling_enabled(enabled);
    custom.draw_frame();
    custom.advance_display();
    HS_EXPECT_TRUE(custom.get_pixel(0, 0) == Pixel(500, 1000, 1500));
    HS_EXPECT_TRUE(custom.get_pixel(47, 19) == Pixel(500, 1000, 1500));
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
  Pullback::TemporalHistory<FX::FrameState, Pullback::HueMode::NOISE, true,
                            true>
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
  const Pullback::TemporalFrameSampler EQUAL_SAMPLER(history.previous,
                                                     equal_frame);
  auto equal_sample = EQUAL_SAMPLER.sample(0.5f);
  HS_EXPECT_TRUE(equal_sample.temporal.same_noise);
  auto different_noise = NOISE;
  different_noise[0] ^= 1;
  auto different_frame = equal_frame;
  different_frame.hue_noise_lut = different_noise.data();
  const Pullback::TemporalFrameSampler DIFFERENT_SAMPLER(history.previous,
                                                         different_frame);
  HS_EXPECT_TRUE(!DIFFERENT_SAMPLER.sample(0.5f).temporal.same_noise);
  different_frame.hue_noise_lut = nullptr;
  const Pullback::TemporalFrameSampler NO_NOISE_SAMPLER(history.previous,
                                                        different_frame);
  HS_EXPECT_TRUE(!NO_NOISE_SAMPLER.sample(0.5f).temporal.same_noise);
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
  Pullback::TemporalHistory<FX::FrameState, Pullback::HueMode::NOISE, true>
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
          Pullback::sample_frame(history.previous, current, TIME));
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
  Pullback::TemporalHistory<FX::FrameState, Pullback::HueMode::NOISE, true>
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
        Pullback::sample_frame(history.previous, current, TIME));
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

  Temporal::set_temporal_sampling_enabled(false);
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
  Temporal::set_temporal_sampling_enabled(true);
  effect.draw_frame();
  effect.advance_display();
  const auto resumed = effect.frame_for_test();
  for (int x = 39; x < 83; ++x) {
    const auto frame = FX::RenderPipeline::prepare(Pullback::sample_frame(
        history.previous, resumed, Temporal::sweep_fraction(x, W)));
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
  if constexpr (std::is_same_v<FX, MobiusGrid<96, 20>>) {
    HS_EXPECT_FALSE(FX::TEMPORAL_SAMPLING);
  } else {
    HS_EXPECT_TRUE(FX::TEMPORAL_SAMPLING);
    reset_globals();
    constexpr int W = 96;
    constexpr int H = 20;
    constexpr int TOLERANCE = FX::Spec::FACTORED_COLUMNS ? 16 : 2;
    FX effect;
    effect.init();
    Pullback::TemporalHistory<typename FX::FrameState, FX::Metadata::HUE, true,
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
      Temporal::set_temporal_sampling_enabled(enabled);
      effect.draw_frame();
      effect.advance_display();
      const auto current = effect.frame_for_test();
      const auto instantaneous = FX::RenderPipeline::prepare(current);
      const auto previous = FX::RenderPipeline::prepare(history.previous);
      const auto start = FX::RenderPipeline::prepare(
          Pullback::sample_frame(history.previous, current, 0.0f));
      const auto end = FX::RenderPipeline::prepare(
          Pullback::sample_frame(history.previous, current, 1.0f));
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
            Pullback::sample_frame(history.previous, current, FRACTION));
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
    Temporal::set_temporal_sampling_enabled(true);
  }
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
