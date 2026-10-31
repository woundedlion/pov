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
 * @return Column context with resources borrowed from the endpoints. */
template <typename Frame>
HS_FLASH_INLINE Frame sample_frame_reference(const Frame &from, const Frame &to,
                                             float fraction) {
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
  frame.params = Pullback::interpolate(from.params, to.params, fraction);
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
      if constexpr (std::is_same_v<
                        typename Resource::Family,
                        Pullback::SphereAnimation::BurstRippleParams>) {
        const uint32_t DELTA = b.frames - a.frames;
        for (auto &wave : clock.waves)
          wave.age -= DELTA * (1.0f - fraction);
        clock.radial_fraction = DELTA == 1 ? fraction : 1.0f;
        if (DELTA > 1 ||
            std::memcmp(&family, &from.params.template get<Resource::KEY>(),
                        sizeof(family)) != 0 ||
            std::memcmp(&family, &to.params.template get<Resource::KEY>(),
                        sizeof(family)) != 0)
          clock.radial_tables = nullptr;
      } else if constexpr (std::is_same_v<typename Resource::Family,
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
  frame.temporal = {from.palette, from.hue_rotation_lut, fraction};
  return frame;
}

template <bool Persistent, bool Temporal>
struct PreparedFrameProbeSpec : GridSpaceSpec {
  static constexpr bool PERSISTENT_PREPARED_FRAME = Persistent;
  static constexpr bool TEMPORAL_SAMPLING = Temporal;
};

template <bool Persistent, bool Temporal>
class PreparedFrameProbe
    : public Pullback::ComposedEffect<
          SMALL_W, SMALL_H, PreparedFrameProbe<Persistent, Temporal>,
          PreparedFrameProbeSpec<Persistent, Temporal>> {
  using Base =
      Pullback::ComposedEffect<SMALL_W, SMALL_H,
                               PreparedFrameProbe<Persistent, Temporal>,
                               PreparedFrameProbeSpec<Persistent, Temporal>>;

public:
  static constexpr std::string_view EFFECT_ID = "prepared-frame-probe";
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"probe"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr auto initial_params() {
    return GridSpace<SMALL_W, SMALL_H>::initial_params();
  }
  static consteval auto choreography_hooks() {
    return Base::choreography_hooks()
        .template extend<
            PreparedFrameProbe,
            Control::HookBinding<Control::ChoreographyHook::INITIAL_PARAMS,
                                 &PreparedFrameProbe::initial_params>>();
  }
};

template <bool Persistent, bool Temporal>
inline std::vector<Pixel> prepared_frame_pixels(bool sampling) {
  reset_effect_globals();
  const ScopedTemporalSampling mode(sampling);
  using FX = PreparedFrameProbe<Persistent, Temporal>;
  FX effect;
  effect.init();
  const typename FX::Frame *retained = nullptr;
  if constexpr (Persistent)
    HS_EXPECT_TRUE(ComposedFrameWhiteBox::prepared_frame(effect) == nullptr);
  std::vector<Pixel> pixels;
  for (int frame = 0; frame < 3; ++frame) {
    effect.draw_frame();
    if constexpr (Persistent) {
      const auto *prepared = ComposedFrameWhiteBox::prepared_frame(effect);
      HS_EXPECT_TRUE(prepared != nullptr);
      if (retained)
        HS_EXPECT_TRUE(prepared == retained);
      retained = prepared;
      const typename FX::Frame owned(prepared->ctx);
      for (const math::Vector &direction :
           {math::UP, math::Vector(1, 0, 0), math::Vector(0, 0, 1)}) {
        const Color4 actual = FX::shade(direction, *prepared);
        const Color4 expected = FX::shade(direction, owned);
        HS_EXPECT_EQ(actual.color, expected.color);
        HS_EXPECT_EQ(actual.alpha, expected.alpha);
      }
    }
    for (int y = 0; y < SMALL_H; ++y)
      for (int x = 0; x < SMALL_W; ++x)
        pixels.push_back(effect.get_pixel(x, y));
    effect.advance_display();
  }
  return pixels;
}

inline void test_composed_persistent_prepared_frame() {
  static_assert(!Pullback::Spec::PERSISTENT_PREPARED_FRAME);
  static_assert(MillefioriSpec::PERSISTENT_PREPARED_FRAME);
  for (bool sampling : {false, true}) {
    HS_EXPECT_TRUE((prepared_frame_pixels<true, true>(sampling) ==
                    prepared_frame_pixels<false, true>(sampling)));
    HS_EXPECT_TRUE((prepared_frame_pixels<true, false>(sampling) ==
                    prepared_frame_pixels<false, false>(sampling)));
  }
}

inline void test_composed_document_modulation_members() {
  namespace Document = hs_wasm::document_parameters;
  reset_effect_globals();
  Millefiori<SMALL_W, SMALL_H> effect;
  effect.init();
  const auto read = [&] {
    ColorModulation::Config config;
    size_t matches = 0;
    effect.visit_document_modulations(
        [&](std::string_view resource, const auto &current) {
          if (resource == "field") {
            config = current;
            ++matches;
          }
        });
    HS_EXPECT_EQ(matches, size_t{1});
    return config;
  };
  const auto initial = read();
  HS_EXPECT_EQ(initial.source, Pattern::SourceKind::SPHERE_SMOOTH_GRID);
  const auto write = [&](std::string_view resource, std::string_view id,
                         float value) {
    return effect.edit_document_modulation(resource, [&](auto &candidate) {
      return Document::write(candidate, id, value);
    });
  };
  HS_EXPECT_EQ(write("field", "grid-speed", 0.03f), ParamSetResult::APPLIED);
  HS_EXPECT_EQ(read().grid.speed, 0.03f);
  HS_EXPECT_EQ(read().noise.speed, initial.noise.speed);
  HS_EXPECT_TRUE(effect.animations_paused());
  HS_EXPECT_EQ(write("field", "pattern-speed", 0.04f),
               ParamSetResult::READONLY);
  HS_EXPECT_EQ(write("field", "modulation-source", 0.0f),
               ParamSetResult::READONLY);
  HS_EXPECT_EQ(write("modulate", "grid-speed", 0.04f),
               ParamSetResult::READONLY);
  HS_EXPECT_EQ(write("missing", "grid-speed", 0.04f),
               ParamSetResult::UNKNOWN_PARAM);
  HS_EXPECT_EQ(write("field", "missing", 0.04f), ParamSetResult::UNKNOWN_PARAM);
  const auto before = read();
  HS_EXPECT_EQ(
      write("field", "grid-speed", std::numeric_limits<float>::quiet_NaN()),
      ParamSetResult::NON_FINITE);
  HS_EXPECT_EQ(write("field", "grid-frequency", -1.0f),
               ParamSetResult::INADMISSIBLE);
  HS_EXPECT_TRUE(read() == before);
  HS_EXPECT_EQ(effect.edit_document_modulation("field",
                                               [](auto &candidate) {
                                                 candidate.source =
                                                     Pattern::SourceKind::NOISE;
                                                 return ParamSetResult::APPLIED;
                                               }),
               ParamSetResult::INADMISSIBLE);
  HS_EXPECT_TRUE(read() == before);

  auto config = initial;
  config.grid.axis = math::Vector(1, 0, 0);
  const auto valid = config;
  HS_EXPECT_EQ(Document::write(config, "grid-axis-x", 0.0f),
               ParamSetResult::INADMISSIBLE);
  HS_EXPECT_TRUE(config == valid);
  for (const auto source :
       {Pattern::SourceKind::NOISE, Pattern::SourceKind::DIRECTIONAL_SINE,
        Pattern::SourceKind::SPHERE_INTERFERENCE}) {
    config = ColorModulation::Config{};
    config.source = source;
    config.evaluation = Pattern::Evaluation::DIRECT;
    const char *id = source == Pattern::SourceKind::NOISE ? "pattern-speed"
                     : source == Pattern::SourceKind::DIRECTIONAL_SINE
                         ? "sine-speed"
                         : "interference-speed";
    HS_EXPECT_EQ(Document::write(config, id, 0.0005f), ParamSetResult::APPLIED);
    if (source == Pattern::SourceKind::NOISE)
      HS_EXPECT_EQ(Document::write(config, id, 0.01f),
                   ParamSetResult::INADMISSIBLE);
    HS_EXPECT_EQ(Document::write(config, "grid-speed", 0.01f),
                 ParamSetResult::READONLY);
  }
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
        sample_frame_reference(surface_from, surface_to, FRACTION);
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
          sample_frame_reference(ripple_from, ripple_to, FRACTION);
      HS_EXPECT_NEAR(actual.resources.template get<"surface">().phase,
                     expected.resources.template get<"surface">().phase, 0.0f);
      if (column == 72 && period == 3.0f)
        HS_EXPECT_NEAR(actual.resources.template get<"surface">().phase, 2.25f,
                       1e-6f);
    }
  }
}

/** @brief Affine scales and shear follow the periodic warp cycle. */
inline void test_temporal_affine_cycle() {
  HS_CONTEXT("temporal affine cycle");
  for (float phase : {0.0f, 0.25f, 0.5f, 0.99f, 1.01f}) {
    Pullback::AffineParams warp;
    for (const auto &params :
         {std::array{1.0f, 1.0f, 0.0f}, std::array{0.7f, 1.8f, 0.4f},
          std::array{1.8f, 0.7f, -0.6f}}) {
      warp.scale_x = params[0];
      warp.scale_y = params[1];
      warp.shear = params[2];
      const auto actual = Pullback::Warp::prepare(warp, phase, 0.3f, 1.2f);
      const float CYCLE = cosf(math::TWO_PI_F * math::wrap_t(phase));
      HS_EXPECT_NEAR(actual.transform.affine.scale_x, powf(warp.scale_x, CYCLE),
                     0.0f);
      HS_EXPECT_NEAR(actual.transform.affine.scale_y, powf(warp.scale_y, CYCLE),
                     0.0f);
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
  Spacetime::FrameHistory<FX::FrameState, true, true> history;
  history.capture(current, persistent_arena);
  using PathPalette = Pullback::Color::GeneratedPalette<Pullback::ColorProvider<
      FX::Binding, Pullback::Color::BrightnessEnvelope::CUP>>;
  const Pullback::FrameInterval EQUAL(history.previous, current);
  auto sample = EQUAL.at(0.5f);
  HS_EXPECT_TRUE(Spacetime::sample_frame(history.previous, current, 0.5f)
                     .temporal.same_palette);
  HS_EXPECT_TRUE(sample.temporal.same_palette);
  sample.params.template get<"modulation">().modulation.mechanic.amount = 0.3f;
  sample.params.template get<"modulation">().modulation.enabled = false;
  const auto prepared = PathPalette::prepare(sample);
  HS_EXPECT_NEAR(prepared.fraction, 1.0f, 0.0f);
  auto interpolated = sample;
  interpolated.temporal.same_palette = false;
  interpolated.temporal.fraction = 0.5f;
  const auto blended = PathPalette::prepare(interpolated);
  HS_EXPECT_NEAR(blended.fraction, 0.5f, 0.0f);
  for (int index = 0; index < 257; ++index) {
    Pullback::FieldSample field{};
    field.value = static_cast<float>(index) / 256;
    field.coverage = 0.75f;
    field.path_length = static_cast<float>(index) * 0.37f;
    const Color4 A =
        Pullback::Color::finish(PathPalette::apply(field, sample, prepared));
    const Color4 B = Pullback::Color::finish(
        PathPalette::apply(field, interpolated, blended));
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
  Spacetime::FrameHistory<FX::FrameState, true, true> history;
  history.capture(current, persistent_arena);
  auto moved = current;
  moved.params.template get<"modulation">().modulation.mechanic.amount += 0.25f;
  auto &start_config =
      history.previous.params.template get<"modulation">().modulation;
  auto &end_config = moved.params.template get<"modulation">().modulation;
  start_config.noise.scale = 1.0f;
  start_config.noise.speed = -0.001f;
  end_config.noise.scale = 4.0f;
  end_config.noise.speed = 0.001f;
  const auto MIDPOINT =
      Pullback::interpolate(history.previous.params, moved.params, 0.5f);
  HS_EXPECT_NEAR(MIDPOINT.template get<"modulation">().modulation.noise.scale,
                 2.0f, 1e-6f);
  HS_EXPECT_NEAR(MIDPOINT.template get<"modulation">().modulation.noise.speed,
                 0.0f, 1e-8f);

  for (bool variant_mismatch : {false, true}) {
    Pullback::Modulation::Config start, end;
    start.noise.basis = Pattern::NoiseBasis::FBM3;
    start.noise.variant = Pattern::NoiseVariant::OCTAVE_CONTOUR;
    start.noise.contrast = 2.0f;
    end = start;
    if (variant_mismatch) {
      start.noise.basis = end.noise.basis = Pattern::NoiseBasis::SIMPLEX;
      end.noise.variant = Pattern::NoiseVariant::CONFIGURED_SINGLE;
      end.noise.contrast = 0.0f;
    } else {
      end.noise.basis = Pattern::NoiseBasis::RIDGED3;
      end.noise.contrast = 6.0f;
    }
    HS_EXPECT_TRUE(Pullback::Modulation::valid(start));
    HS_EXPECT_TRUE(Pullback::Modulation::valid(end));
    for (float fraction : {.25f, .5f, .75f}) {
      const auto sampled =
          Pullback::Modulation::interpolate(start, end, fraction);
      HS_EXPECT_TRUE(Pullback::Modulation::valid(sampled));
      HS_EXPECT_EQ(sampled.noise.basis, end.noise.basis);
      HS_EXPECT_EQ(sampled.noise.variant, end.noise.variant);
      HS_EXPECT_EQ(sampled.noise.contrast, end.noise.contrast);
    }
  }
  const Pullback::FrameInterval SAMPLER(history.previous, moved);
  for (const float fraction : {0.1f, 0.37f, 0.5f, 0.9f}) {
    HS_CONTEXT("fraction", static_cast<double>(fraction));
    const auto sampled = SAMPLER.at(fraction).params;
    const auto expected =
        Pullback::interpolate(history.previous.params, moved.params, fraction);
    HS_EXPECT_NEAR(
        sampled.template get<"modulation">().modulation.mechanic.amount,
        expected.template get<"modulation">().modulation.mechanic.amount, 0.0f);
    HS_EXPECT_TRUE(std::memcmp(&sampled, &expected, sizeof(sampled)) == 0);
  }
}

/** @brief Retained color tables survive subsequent bakes. */
inline void test_temporal_owned_resources() {
  HS_CONTEXT("temporal owned modulation resources");
  reset_globals();
  using FX = GridSpace<96, 20>;
  FX effect;
  effect.init();
  auto params = FX::initial_params();
  params.template get<"modulation">().modulation.noise.speed = .001f;
  ComposedFrameWhiteBox::set_params(effect, params);
  effect.draw_frame();
  effect.advance_display();
  const auto current = effect.frame_for_test();
  Spacetime::FrameHistory<FX::FrameState, true, true> history;
  history.capture(current, persistent_arena);
  HS_EXPECT_TRUE(history.previous.palette != current.palette);
  HS_EXPECT_TRUE(
      history.previous.resources.template get<"modulation">().modulation !=
      current.resources.template get<"modulation">().modulation);
  HS_EXPECT_TRUE(
      history.previous.resources.template get<"modulation">()
          .modulation->cube.data !=
      current.resources.template get<"modulation">().modulation->cube.data);
  std::array<float, 8> retained_samples{}, live_samples{};
  const auto direction = [](size_t index) {
    return Pattern::SphereDirection::from(
        math::Vector(1.0f, float(index) * .3f - 1.0f, .37f));
  };
  for (size_t index = 0; index < retained_samples.size(); ++index) {
    retained_samples[index] = Pullback::Modulation::sample(
        *history.previous.resources.template get<"modulation">().modulation,
        direction(index));
    live_samples[index] = Pullback::Modulation::sample(
        *current.resources.template get<"modulation">().modulation,
        direction(index));
    HS_EXPECT_EQ(retained_samples[index], live_samples[index]);
  }
  const auto TABLE = history.resources.template get<"modulation">().table;
  const auto RECIPE =
      history.resources.template get<"modulation">().modulation.config;
  const Pixel PALETTE = history.previous.palette->get(0.4f).color;
  for (int frame = 0; frame < 4; ++frame) {
    effect.draw_frame();
    effect.advance_display();
  }
  bool live_changed = false;
  for (size_t index = 0; index < retained_samples.size(); ++index) {
    HS_EXPECT_EQ(
        Pullback::Modulation::sample(
            *history.previous.resources.template get<"modulation">().modulation,
            direction(index)),
        retained_samples[index]);
    live_changed |=
        Pullback::Modulation::sample(
            *current.resources.template get<"modulation">().modulation,
            direction(index)) != live_samples[index];
  }
  HS_EXPECT_TRUE(live_changed);
  HS_EXPECT_TRUE(history.resources.template get<"modulation">().table == TABLE);
  HS_EXPECT_TRUE(
      history.resources.template get<"modulation">().modulation.config ==
      RECIPE);
  HS_EXPECT_TRUE(history.previous.palette->get(0.4f).color == PALETTE);

  using namespace Pullback;
  Modulation::Config config;
  config.enabled = true;
  config.source = Pattern::SourceKind::DIRECTIONAL_SINE;
  config.evaluation = Pattern::Evaluation::DIRECT;
  config.sine.frequency = 0.0f;
  config.mechanic = {ColorModulation::Mechanic::OPACITY, 1.0f};
  Modulation::Resources<false> a, b;
  auto clock = Modulation::init(config);
  clock.phase.phase_turns = 0.75f;
  a.prepare(config, clock);
  clock.phase.phase_turns = 0.25f;
  b.prepare(config, clock);
  auto from = current;
  auto to = current;
  from.params.template get<"modulation">().modulation = config;
  to.params.template get<"modulation">().modulation = config;
  to.params.template get<"modulation">().modulation.source =
      Pattern::SourceKind::SPHERE_SMOOTH_GRID;
  from.params.template get<"modulation">().modulation.mechanic = {
      ColorModulation::Mechanic::OPACITY, 0.2f};
  to.params.template get<"modulation">().modulation.mechanic = {
      ColorModulation::Mechanic::OPACITY, 0.8f};
  from.resources.template get<"modulation">().modulation = &a.prepared;
  to.resources.template get<"modulation">().modulation = &b.prepared;
  for (float fraction : {0.0f, 0.25f, 0.5f, 1.0f}) {
    const auto COLUMN = FrameInterval(from, to).at(fraction);
    HS_EXPECT_TRUE(COLUMN.resources.template get<"modulation">().modulation ==
                   to.resources.template get<"modulation">().modulation);
    HS_EXPECT_NEAR(
        COLUMN.params.template get<"modulation">().modulation.mechanic.amount,
        0.8f, 0.0f);
  }
  config = {};
  config.enabled = false;
  config.source = Pattern::SourceKind::DIRECTIONAL_SINE;
  config.evaluation = Pattern::Evaluation::DIRECT;
  config.sine.speed = -0.125f;
  clock = Modulation::init(config);
  Modulation::advance(clock, config);
  HS_EXPECT_NEAR(clock.phase.phase_turns, 0.875f, 0.0f);
  config.mechanic.kind = ColorModulation::Mechanic::OPACITY;
  Modulation::synchronize(clock, config);
  HS_EXPECT_NEAR(clock.phase.phase_turns, 0.875f, 0.0f);
  config.source = Pattern::SourceKind::SPHERE_SMOOTH_GRID;
  Modulation::synchronize(clock, config);
  HS_EXPECT_NEAR(clock.phase.phase_turns, 0.0f, 0.0f);
  config.grid.axis_b = config.grid.axis;
  HS_EXPECT_TRUE(Modulation::valid(config));
  config.grid.axis_b = math::Vector(0.0f, 0.0f, 0.0f);
  HS_EXPECT_TRUE(!Modulation::valid(config));
  config.grid.axis_b = math::X_AXIS;
  config.version = 2;
  HS_EXPECT_TRUE(!Modulation::valid(config));

  config = {};
  config.evaluation = Pattern::Evaluation::DIRECT;
  config.noise.variant = Pattern::NoiseVariant::OCTAVE_CONTOUR;
  for (auto mechanic : {ColorModulation::Mechanic::HUE_ROTATION,
                        ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT,
                        ColorModulation::Mechanic::OPACITY}) {
    config.mechanic.kind = mechanic;
    HS_EXPECT_TRUE(Modulation::valid(config));
    HS_EXPECT_TRUE(Modulation::admitted(config));
    HS_EXPECT_FALSE(Modulation::device_qualified(config));
    config.enabled = false;
    HS_EXPECT_TRUE(Modulation::valid(config));
    HS_EXPECT_FALSE(Modulation::device_qualified(config));
    config.enabled = true;
  }
  config.evaluation = Pattern::Evaluation::BAKED_CUBE;
  HS_EXPECT_TRUE(Modulation::valid(config));
  config.mechanic.kind = ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT;
  HS_EXPECT_TRUE(Modulation::valid(config));
  HS_EXPECT_TRUE(Modulation::admitted(config));
  HS_EXPECT_FALSE(Modulation::device_qualified(config));
  config.noise.variant = Pattern::NoiseVariant::CONFIGURED_SINGLE;
  config.mechanic.kind = ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT;
  config.evaluation = Pattern::Evaluation::DIRECT;
  HS_EXPECT_TRUE(Modulation::valid(config));
  config.enabled = false;
  HS_EXPECT_TRUE(Modulation::valid(config));
  config.evaluation = Pattern::Evaluation::BAKED_CUBE;
  HS_EXPECT_TRUE(Modulation::valid(config));
  config = {};
  config.enabled = true;
  config.mechanic.amount = 0.4f;
  Modulation::Resources<> cube;
  clock = Modulation::init(config);
  cube.prepare(config, clock);
  HS_EXPECT_EQ(cube.bake_generation, uint32_t{1});
  const auto SCALAR_TABLE = cube.table;
  config.mechanic.kind = ColorModulation::Mechanic::OPACITY;
  config.mechanic.amount = 0.7f;
  cube.prepare(config, clock);
  HS_EXPECT_EQ(cube.bake_generation, uint32_t{1});
  HS_EXPECT_TRUE(cube.table == SCALAR_TABLE);
  config.sine.offset += 0.25f;
  cube.prepare(config, clock);
  HS_EXPECT_EQ(cube.bake_generation, uint32_t{1});
  ++config.seed;
  cube.prepare(config, clock);
  HS_EXPECT_EQ(cube.bake_generation, uint32_t{2});
  HS_EXPECT_TRUE(cube.table != SCALAR_TABLE);
  cube.reset();
  cube.prepare(config, clock);
  HS_EXPECT_EQ(cube.bake_generation, uint32_t{3});
  config.evaluation = Pattern::Evaluation::DIRECT;
  clock.phase.phase_turns = 0.375f;
  Modulation::synchronize(clock, config);
  HS_EXPECT_NEAR(clock.phase.phase_turns, 0.375f, 0.0f);
  const auto DIRECT = cube.prepare(config, clock);
  HS_EXPECT_TRUE(!DIRECT.cube.active);
  cube.prepare(config, clock);
  HS_EXPECT_NEAR(clock.phase.phase_turns, 0.375f, 0.0f);
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
  Spacetime::FrameHistory<FX::FrameState, true, true> history;
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
  Spacetime::FrameHistory<FX::FrameState, true, true> history;
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
  Spacetime::FrameHistory<typename FX::FrameState, true,
                          FX::Metadata::MODULATION>
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
  const auto initial = ComposedFrameWhiteBox::current_frame(effect);
  for (int x : {0, 1, 47, 48, 49, 95}) {
    const auto FIRST = FX::RenderPipeline::prepare(Spacetime::sample_frame(
        initial, initial, static_cast<float>(x % (W / 2)) / (W / 2)));
    for (int y : {1, 7, 18})
      expect_pixel(effect.get_pixel(x, y), FX::shade(view(x, y), FIRST));
  }
  history.capture(initial, persistent_arena);
  HS_EXPECT_TRUE(history.previous.palette != initial.palette);
  if constexpr (FX::Params::template HAS<"modulation">) {
    HS_EXPECT_TRUE(
        history.previous.resources.template get<"modulation">().modulation !=
        initial.resources.template get<"modulation">().modulation);
    HS_EXPECT_TRUE(history.previous.hue_rotation_lut !=
                   initial.hue_rotation_lut);
  }

  effect.set_clip(3, 17, 39, 83);
  for (bool enabled : {true, false, false, true}) {
    Spacetime::set_flat_surface(!enabled);
    effect.draw_frame();
    effect.advance_display();
    const auto current = ComposedFrameWhiteBox::current_frame(effect);
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
  reset_effect_globals();
  const ScopedTemporalSampling instantaneous(false);
  HS_EXPECT_TRUE(Spacetime::flat_surface());
  FX effect;
  effect.init();
  const int departure = FX::PRESET_DWELL_FRAMES +
                        Segue::Preset::frames(FX::preset_departure(0)) + 20;
  PixelGolden golden;
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
    fold_pixel_golden<SMALL_W, SMALL_H>(golden, effect);
  }
  hs::clear_mock_time();
  expect_pixel_golden(
      "MobiusGrid 96x20", golden,
      RecordedGolden{{0x6d66ce5d4876ef89ull,
                      {{6805330220ull, 6844874216ull, 7112697845ull},
                       201312982882ull,
                       983633955317ull}},
                     {0x69e38ef852f070daull,
                      {{6568081584ull, 6978293201ull, 6966434110ull},
                       195726651229ull,
                       957665657686ull}},
                     1e-3,
                     0.01,
                     {0x54b286aedc50e5b7ull,
                      {{7025810653ull, 6816434558ull, 7032280424ull},
                       208885953563ull,
                       984308556083ull}},
                     {0x8c65f65a113972deull,
                      {{6401113354ull, 6273991381ull, 6152479994ull},
                       179644043625ull,
                       913767027589ull}}});
}

/** @brief Composed history type of an effect.
 * @tparam FX Composed effect specialization. */
template <typename FX>
using ComposedHistory = Spacetime::FrameHistory<typename FX::FrameState, true,
                                                FX::Metadata::MODULATION>;

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
            sample_frame_reference(history.previous, current, fraction);
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

/** @brief Each baked modulation retains its own immutable cube and source clock. */
inline void test_temporal_stacked_modulation_resources() {
  reset_globals();
  StackedModulationProbe effect;
  effect.init();
  HS_EXPECT_TRUE(effect.getParameters().find("modulation.Modulation Amount") !=
                 nullptr);
  HS_EXPECT_TRUE(effect.getParameters().find("second.Modulation Amount") !=
                 nullptr);
  effect.draw_frame();
  effect.advance_display();
  const auto current = effect.frame_for_test();
  const auto *a = current.resources.template get<"modulation">().modulation;
  const auto *b = current.resources.template get<"second">().modulation;
  HS_EXPECT_TRUE(a != b);
  HS_EXPECT_TRUE(a->cube.data != b->cube.data);
  HS_EXPECT_NEAR(a->state.phase.phase_turns, 0.001f, 1e-7f);
  HS_EXPECT_NEAR(b->state.phase.phase_turns, 0.9995f, 1e-7f);
  Spacetime::FrameHistory<StackedModulationProbe::FrameState, true, true>
      history;
  history.capture(current, persistent_arena);
  const auto *old_a =
      history.previous.resources.template get<"modulation">().modulation;
  const auto *old_b =
      history.previous.resources.template get<"second">().modulation;
  HS_EXPECT_TRUE(old_a->cube.data != old_b->cube.data);
  HS_EXPECT_TRUE(old_a->cube.data != a->cube.data);
  HS_EXPECT_TRUE(old_b->cube.data != b->cube.data);
  const auto DIRECTION =
      Pattern::SphereDirection::from(math::Vector(0.3f, 0.8f, 0.2f));
  const float A = Pullback::Modulation::sample(*a, DIRECTION);
  const float B = Pullback::Modulation::sample(*b, DIRECTION);
  HS_EXPECT_TRUE(A != B);
  const auto BEFORE =
      StackedModulationProbe::RenderPipeline::prepare(history.previous);
  const auto PIXEL_BEFORE = StackedModulationProbe::shade(
      math::Vector(0.3f, 0.8f, 0.2f).normalized(), BEFORE);
  for (int i = 0; i < 4; ++i) {
    effect.draw_frame();
    effect.advance_display();
  }
  HS_EXPECT_EQ(Pullback::Modulation::sample(*old_a, DIRECTION), A);
  HS_EXPECT_EQ(Pullback::Modulation::sample(*old_b, DIRECTION), B);
  const auto AFTER =
      StackedModulationProbe::RenderPipeline::prepare(history.previous);
  const auto PIXEL_AFTER = StackedModulationProbe::shade(
      math::Vector(0.3f, 0.8f, 0.2f).normalized(), AFTER);
  HS_EXPECT_EQ(PIXEL_BEFORE.color, PIXEL_AFTER.color);
  HS_EXPECT_NEAR(PIXEL_BEFORE.alpha, PIXEL_AFTER.alpha, 0.0f);
}

/** @brief Millefiori samples stable analytic grid phases along the shortest arc. */
inline void test_millefiori_column_color_phase() {
  using namespace Pullback;
  using FX = Millefiori<SMALL_W, SMALL_H>;
  using Provider = ModulationProvider<typename FX::Binding, "field",
                                      Pattern::SourceKind::SPHERE_SMOOTH_GRID,
                                      Pattern::Evaluation::DIRECT, true>;
  reset_effect_globals();
  FX effect;
  effect.init();
  auto from = ComposedFrameWhiteBox::frame(effect);
  auto to = from;
  auto config = from.params.template get<"field">().modulation;
  config.enabled = true;
  config.mechanic = {ColorModulation::Mechanic::POSITIVE_LIGHTNESS_LIFT, 0.35f};
  config.grid.axis = math::Y_AXIS;
  config.grid.axis_b = math::X_AXIS;
  config.grid.frequency = 3.0f;
  config.grid.frequency_b = 2.0f;
  config.grid.offset = 0.0f;
  config.grid.offset_b = 0.25f;
  Modulation::Resources<false> previous, current;
  auto clock = Modulation::init(config);
  clock.phase.phase_turns = 0.99f;
  previous.prepare(config, clock);
  clock.phase.phase_turns = 0.01f;
  current.prepare(config, clock);
  from.params.template get<"field">().modulation = config;
  to.params.template get<"field">().modulation = config;
  from.resources.template get<"field">().modulation = &previous.prepared;
  to.resources.template get<"field">().modulation = &current.prepared;
  for (float fraction : {0.0f, 0.25f, 0.5f, 0.75f, 1.0f}) {
    const auto COLUMN = FrameInterval(from, to).at(fraction);
    const auto PREPARED = Provider::prepare(COLUMN);
    HS_EXPECT_EQ(PREPARED.ready, fraction > 0.0f && fraction < 1.0f);
    if (!PREPARED.ready)
      continue;
    HS_EXPECT_NEAR(PREPARED.grid.a.phase_turns, 0.99f + 0.02f * fraction,
                   1e-7f);
    HS_EXPECT_NEAR(PREPARED.grid.b.phase_turns, 0.24f + 0.02f * fraction,
                   1e-7f);
    for (const auto &direction :
         {math::X_AXIS, math::Y_AXIS, math::Z_AXIS,
          math::Vector(0.3f, 0.4f, 0.5f).normalized()}) {
      const float A = math::sinf_0_2pi(
          math::TWO_PI_F *
          math::wrap_t(3.0f * direction.y + 0.99f + 0.02f * fraction));
      const float B = math::sinf_0_2pi(
          math::TWO_PI_F *
          math::wrap_t(2.0f * direction.x + 0.24f + 0.02f * fraction));
      ColorSample sample{Color4(Pixel(2000, 4000, 6000), 0.37f),
                         {0.5f, 0.7f, direction, 2.0f}};
      Color::Modulate<Provider>::apply_phase_in_place(sample, COLUMN, PREPARED);
      HS_EXPECT_NEAR(sample.pending_lift, 0.35f * fmaxf(0.0f, A * B), 1e-6f);
      HS_EXPECT_EQ(sample.color.alpha, 0.37f);
    }
  }
  for (int change = 0; change < 8; ++change) {
    auto changed = to;
    auto &recipe = changed.params.template get<"field">().modulation;
    if (change == 0)
      recipe.grid.frequency += 0.1f;
    if (change == 1)
      recipe.grid.axis = math::Z_AXIS;
    if (change == 2)
      recipe.grid.offset_b += 0.1f;
    if (change == 3)
      recipe.grid.speed += 0.1f;
    if (change == 4)
      recipe.enabled = !recipe.enabled;
    if (change == 5)
      recipe.mechanic.kind = ColorModulation::Mechanic::HUE_ROTATION;
    if (change == 6)
      recipe.source = Pattern::SourceKind::DISPLACEMENT;
    if (change == 7)
      recipe.evaluation = Pattern::Evaluation::BAKED_CUBE;
    const auto COLUMN = FrameInterval(from, changed).at(0.5f);
    HS_EXPECT_FALSE(COLUMN.temporal.phase_compatible);
    HS_EXPECT_FALSE(Provider::prepare(COLUMN).ready);
  }
}

/** @brief Palette-base slots preserve the other endpoint and skip changed controls. */
inline void test_millefiori_palette_lab_lifetime() {
  using FX = Millefiori<SMALL_W, SMALL_H>;
  reset_effect_globals();
  Spacetime::set_flat_surface(false);
  pin_frame_clock(0);
  FX effect;
  effect.init();
  pin_frame_clock(1);
  effect.draw_frame();
  const auto FIRST = ComposedFrameWhiteBox::current_frame(effect);
  HS_EXPECT_TRUE(FIRST.palette_lab_lut && FIRST.palette_lab_lut->active);
  const auto FIRST_NODES = FIRST.palette_lab_lut->nodes;
  const auto FIRST_VIEW = FIRST.palette_lab_lut->view(0.37f);
  const auto FIRST_LAB = PaletteLabLut::sample(FIRST_VIEW, 0.43f);
  Pullback::ColorSample carrier{Color4(Pixel(0, 0, 0), 0.37f),
                                {0.5f, 0.7f, math::Y_AXIS, 2.0f}};
  carrier.cached_lab = FIRST_LAB;
  carrier.has_cached_lab = true;
  const auto FIRST_COLOR = Pullback::finish_color_sample(carrier);
  effect.advance_display();
  ComposedFrameWhiteBox::advance_palette(effect,
                                         GeneratedPaletteBank::FADE_FRAMES / 2);
  pin_frame_clock(2);
  effect.draw_frame();
  const auto SECOND = ComposedFrameWhiteBox::current_frame(effect);
  HS_EXPECT_TRUE(SECOND.palette_lab_lut != FIRST.palette_lab_lut);
  HS_EXPECT_EQ(std::memcmp(FIRST_NODES.data(),
                           FIRST.palette_lab_lut->nodes.data(),
                           sizeof(FIRST_NODES)),
               0);
  const auto SECOND_NODES = SECOND.palette_lab_lut->nodes;
  HS_EXPECT_TRUE(std::memcmp(FIRST_NODES.data(), SECOND_NODES.data(),
                             sizeof(FIRST_NODES)) != 0);
  effect.advance_display();
  ComposedFrameWhiteBox::advance_palette(effect,
                                         GeneratedPaletteBank::FADE_FRAMES / 2);
  pin_frame_clock(3);
  effect.draw_frame();
  const auto THIRD = ComposedFrameWhiteBox::current_frame(effect);
  HS_EXPECT_TRUE(THIRD.palette_lab_lut == FIRST.palette_lab_lut);
  HS_EXPECT_EQ(std::memcmp(SECOND_NODES.data(),
                           SECOND.palette_lab_lut->nodes.data(),
                           sizeof(SECOND_NODES)),
               0);
  HS_EXPECT_EQ(Pullback::finish_color_sample(carrier).color, FIRST_COLOR.color);
  HS_EXPECT_EQ(Pullback::finish_color_sample(carrier).alpha, FIRST_COLOR.alpha);
  const auto THIRD_NODES = THIRD.palette_lab_lut->nodes;
  auto params = THIRD.params;
  params.template get<"color">().palette_chroma += 0.01f;
  ComposedFrameWhiteBox::set_params(effect, params);
  effect.advance_display();
  pin_frame_clock(4);
  effect.draw_frame();
  HS_EXPECT_TRUE(ComposedFrameWhiteBox::current_frame(effect).palette_lab_lut ==
                 nullptr);
  HS_EXPECT_EQ(std::memcmp(THIRD_NODES.data(),
                           THIRD.palette_lab_lut->nodes.data(),
                           sizeof(THIRD_NODES)),
               0);
  HS_EXPECT_TRUE(ComposedFrameWhiteBox::frame(effect).palette_lab_lut ==
                 nullptr);
}

/** @brief Millefiori retains ripple endpoints and samples wave ages continuously. */
inline void test_millefiori_ripple_endpoints() {
  using FX = Millefiori<SMALL_W, SMALL_H>;
  reset_effect_globals();
  FX effect;
  effect.init();
  ComposedFrameWhiteBox::advance(effect);
  const auto from = ComposedFrameWhiteBox::frame(effect);
  const auto &first = from.resources.template get<"displace">();
  HS_EXPECT_EQ(first.frames, 1u);
  HS_EXPECT_EQ(first.burst_count, 1u);
  HS_EXPECT_EQ(first.seed, Pullback::Interp::instance_hash(
                               "displace", "sphere.displace.ripple.v2"));
  HS_EXPECT_NEAR(first.waves[0].age, 1.0f, 0.0f);
  HS_EXPECT_NEAR(first.waves[1].age, -15.0f, 0.0f);
  ComposedFrameWhiteBox::advance(effect);
  const auto to = ComposedFrameWhiteBox::frame(effect);
  const Pullback::FrameInterval interval(from, to);
  const auto quarter = interval.at(0.25f);
  HS_EXPECT_NEAR(quarter.resources.template get<"displace">().waves[0].age,
                 1.25f, 0.0f);
  for (int frame = 0; frame < 31; ++frame)
    ComposedFrameWhiteBox::advance(effect);
  const auto later = ComposedFrameWhiteBox::frame(effect);
  HS_EXPECT_NEAR(first.waves[0].age, 1.0f, 0.0f);
  HS_EXPECT_NEAR(later.resources.template get<"displace">().waves[2].age, 1.0f,
                 0.0f);
  const auto snapshot = effect.serialize_parameters();
  auto invalid = snapshot;
  invalid.params.template get<"displace">().burst_size = 4;
  HS_EXPECT_FALSE(effect.restore_parameters(invalid));
  invalid = snapshot;
  invalid.params.template get<"displace">().random_center = 2;
  HS_EXPECT_FALSE(effect.restore_parameters(invalid));
  HS_EXPECT_TRUE(effect.restore_parameters(snapshot));
}

/** @brief Ripple curves borrow current storage with interval-specific weights. */
inline void test_millefiori_radial_table_interval() {
  using FX = Millefiori<96, 20>;
  using namespace Pullback::SphereAnimation;
  static_assert(FX::Spec::RIPPLE_RADIAL_TABLE);
  reset_globals();
  FX effect;
  effect.init();
  auto previous = effect.frame_for_test();
  ComposedFrameWhiteBox::advance(effect);
  const auto CURRENT = effect.frame_for_test();
  previous.resources.template get<"displace">().radial_tables = nullptr;
  const auto *TABLE =
      CURRENT.resources.template get<"displace">().radial_tables;
  HS_EXPECT_TRUE(TABLE != nullptr);
  HS_EXPECT_TRUE(TABLE->enabled);
  HS_EXPECT_LT(CURRENT.resources.template get<"displace">().radial_fraction,
               0.0f);
  using Provider = Pullback::SurfaceProvider<FX::Binding, BurstRippleParams,
                                             true, "displace">;
  HS_EXPECT_TRUE(Provider::radial_tables(CURRENT) == nullptr);
  const Pullback::FrameInterval interval(previous, CURRENT);
  for (float fraction : {0.0f, 0.25f, 0.75f, 1.0f}) {
    const auto FRAME = interval.at(fraction);
    const auto &resource = FRAME.resources.template get<"displace">();
    HS_EXPECT_EQ(resource.radial_tables, TABLE);
    HS_EXPECT_EQ(resource.radial_fraction,
                 fraction > 0.0f && fraction < 1.0f ? fraction : -1.0f);
    HS_EXPECT_EQ(Provider::radial_tables(FRAME),
                 fraction > 0.0f && fraction < 1.0f ? TABLE : nullptr);
    HS_EXPECT_EQ(Provider::radial_fraction(FRAME), resource.radial_fraction);
  }
  const auto SAME = Pullback::FrameInterval(CURRENT, CURRENT).at(0.5f);
  HS_EXPECT_EQ(SAME.resources.template get<"displace">().radial_fraction, 1.0f);
  auto skipped = CURRENT;
  skipped.resources.template get<"displace">().frames += 2;
  for (float fraction : {0.0f, 0.5f}) {
    const auto FRAME = Pullback::FrameInterval(previous, skipped).at(fraction);
    HS_EXPECT_TRUE(FRAME.resources.template get<"displace">().radial_tables ==
                   nullptr);
  }
  previous.params.template get<"displace">().period += 1.0f;
  const auto CHANGED = Pullback::FrameInterval(previous, CURRENT).at(0.5f);
  HS_EXPECT_TRUE(CHANGED.resources.template get<"displace">().radial_tables ==
                 nullptr);
}
