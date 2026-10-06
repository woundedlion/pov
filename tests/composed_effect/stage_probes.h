/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */

// ============================================================================
// Stage and choreography probes on hand-built composed specs
// ============================================================================

template <bool Animated> struct RippleProbeSpec : Pullback::Spec {
  static constexpr bool ANIMATED_PROJECTION = Animated;
  static constexpr Pullback::HueMode HUE = Pullback::HueMode::PATH_LENGTH;
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Placed<
          Pullback::CodeEmission::INLINE_ONLY,
          Pullback::Stage::Displace<
              Pullback::Surface::PeriodicRipple<Pullback::SurfaceProvider<
                  B, Pullback::PeriodicRippleParams,
                  HUE == Pullback::HueMode::PATH_LENGTH>>>,
          Pullback::Stage::Project<
              typename Pullback::ProjectionPolicyFor<PROJECTION, B>::Type>>,
      Pullback::Stage::Sample<typename Pullback::SourcePolicyFor<
                                  Pullback::GridSourceParams, B>::Type,
                              Pullback::Weight::Projection,
                              Pullback::ProjectionCoverage::Weight>,
      Pullback::Stage::Colorize<Pullback::Color::GeneratedPalette<
          Pullback::ColorProvider<B, HUE, BRIGHTNESS>>>>;
};

template <int W, int H, bool AnimatedProjection = false>
class RippleProbe
    : public Pullback::ComposedEffect<W, H,
                                      RippleProbe<W, H, AnimatedProjection>,
                                      RippleProbeSpec<AnimatedProjection>> {
public:
  using Params = Pullback::ParamsFor<RippleProbeSpec<AnimatedProjection>>;
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"ripple"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;

  static constexpr Params initial_params() {
    Params value;
    value.template get<"surface">().period = 80.0f;
    value.template get<"surface">().strength = 0.15f;
    value.template get<"surface">().decay = 0.0f;
    value.template get<"surface">().thickness = 0.7f;
    return value;
  }
};

inline void test_composed_keyed_parameters() {
  using Params = Pullback::ComposedDetail::ParameterSet<
      Pullback::ComposedDetail::ResourceList<
          Pullback::ParameterResource<"value", Pullback::CutoutValueParams,
                                      Pullback::ResourceKind::VALUE>,
          Pullback::ParameterResource<"alternate", Pullback::CutoutValueParams,
                                      Pullback::ResourceKind::VALUE>>>;
  Params from;
  Params to;
  from.template get<"value">().cutout_threshold = 0.2f;
  from.template get<"alternate">().cutout_threshold = 0.4f;
  to.template get<"value">().cutout_threshold = 0.6f;
  to.template get<"alternate">().cutout_threshold = 0.8f;
  const Params blended = Pullback::interpolate(from, to, 0.5f);
  HS_EXPECT_NEAR(blended.template get<"value">().cutout_threshold, 0.4f, 1e-6f);
  HS_EXPECT_NEAR(blended.template get<"alternate">().cutout_threshold, 0.6f,
                 1e-6f);
  HS_EXPECT_TRUE(Pullback::valid(blended));
  size_t visited = 0;
  blended.visit([&]<typename Resource>(const auto &family) {
    HS_EXPECT_TRUE(&family == &blended.template get<Resource::KEY>());
    ++visited;
  });
  HS_EXPECT_EQ(visited, 2u);
}

struct AffineNamedSourceSpec : Pullback::Spec {
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B,
      Pullback::Stage::Warp<Pullback::Warp::AffineFrame<Pullback::WarpProvider<
          B, "affine", Pullback::AffineParams, false, "cells">>>,
      Pullback::Stage::Sample<Pullback::Source::PrimitiveLattice<
          Pullback::SourceProvider<B, Pullback::LatticeSourceParams, "cells">>>,
      Pullback::Stage::Colorize<Pullback::Color::GeneratedPalette<
          Pullback::ColorProvider<B, HUE, BRIGHTNESS>>>>;
};

inline void test_composed_affine_named_source() {
  using Params = Pullback::ParamsFor<AffineNamedSourceSpec>;
  using Frame = Pullback::FrameState<Params>;
  using Provider =
      Pullback::WarpProvider<Pullback::Binding<Frame>, "affine",
                             Pullback::AffineParams, false, "cells">;
  static_assert(Params::template HAS<"cells"> &&
                !Params::template HAS<"source">);
  Frame frame{};
  frame.params.template get<"cells">().lattice_cell_scale = 2.0f;
  auto &warp = frame.params.template get<"affine">();
  warp.translation_x = 4.0f;
  warp.translation_y = -8.0f;
  frame.resources.template get<"affine">().phase = 0.25f;
  const auto prepared = Provider::prepare(frame);
  HS_EXPECT_NEAR(prepared.transform.affine.translation_x, 0.5f, 1e-6f);
  HS_EXPECT_NEAR(prepared.transform.affine.translation_y, -1.0f, 1e-6f);
}

struct RepeatedStagesSpec : Pullback::Spec {
  static constexpr bool ANIMATED_PROJECTION = false;
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Placed<
          Pullback::CodeEmission::INLINE_ONLY,
          Pullback::Stage::Lens<
              Pullback::Lens::Mobius<Pullback::LensProvider<B, "lens_a">>>,
          Pullback::Stage::Lens<
              Pullback::Lens::Mobius<Pullback::LensProvider<B, "lens_b">>>,
          Pullback::Stage::Project<
              typename Pullback::ProjectionPolicyFor<PROJECTION, B>::Type>>,
      Pullback::Stage::Warp<typename Pullback::WarpPolicyFor<
          Pullback::VectorNoiseParams, B, "warp_a", false>::Type>,
      Pullback::Stage::Warp<typename Pullback::WarpPolicyFor<
          Pullback::VectorNoiseParams, B, "warp_b", false>::Type>,
      Pullback::Stage::Warp<typename Pullback::WarpPolicyFor<
          Pullback::WaveShearParams, B, "warp_c", false>::Type>,
      Pullback::Stage::Sample<typename Pullback::SourcePolicyFor<
          Pullback::GridSourceParams, B>::Type>,
      Pullback::Stage::Colorize<Pullback::Color::GeneratedPalette<
          Pullback::ColorProvider<B, HUE, BRIGHTNESS>>>>;
};
template <int W, int H>
class RepeatedStagesProbe
    : public Pullback::ComposedEffect<W, H, RepeatedStagesProbe<W, H>,
                                      RepeatedStagesSpec> {
public:
  using Params = Pullback::ParamsFor<RepeatedStagesSpec>;
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"repeated"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr Params initial_params() {
    Params value;
    value.template get<"warp_a">().speed = 0.003f;
    value.template get<"warp_b">().speed = 0.009f;
    value.template get<"warp_c">().speed = -0.012f;
    value.template get<"warp_a">().strength = 0.2f;
    value.template get<"warp_b">().strength = 0.4f;
    value.template get<"warp_c">().strength = 0.3f;
    return value;
  }
};

inline void test_composed_repeated_instances() {
  using FX = RepeatedStagesProbe<SMALL_W, SMALL_H>;
  static_assert(FX::Params::template HAS<"warp_a"> &&
                FX::Params::template HAS<"warp_b"> &&
                FX::Params::template HAS<"warp_c">);
  static_assert(!FX::Params::template HAS<"outer_warp"> &&
                !FX::Params::template HAS<"lens"> &&
                !FX::Params::template HAS<"surface">);
  static_assert(FX::RenderPipeline::STAGE_COUNT == 9 &&
                FX::RenderPipeline::NODE_COUNT == 7);
  static_assert(sizeof(FX::Params) ==
                sizeof(Pullback::GridSourceParams) +
                    sizeof(Pullback::ProjectionParams) +
                    2 * sizeof(Pullback::VectorNoiseParams) +
                    sizeof(Pullback::WaveShearParams) +
                    2 * sizeof(Pullback::MobiusLensParams) +
                    sizeof(Pullback::ColorParams));
  reset_effect_globals();
  FX effect;
  effect.init();
  const auto captured = effect.serialize_parameters();
  captured.params.visit([&]<typename Resource>(const auto &) {
    if constexpr (Pullback::HasFields<typename Resource::Family>)
      verify_family_rejection<FX, Resource>(effect, captured);
  });
  for (int step = 0; step < 7; ++step)
    ComposedFrameWhiteBox::advance(effect);
  const auto frame = ComposedFrameWhiteBox::frame(effect);
  HS_EXPECT_NEAR(frame.resources.template get<"warp_a">().phase, 0.021f, 1e-6f);
  HS_EXPECT_NEAR(frame.resources.template get<"warp_b">().phase, 0.063f, 1e-6f);
  HS_EXPECT_NEAR(frame.resources.template get<"warp_c">().phase,
                 math::wrap_t(-0.084f), 1e-6f);
  HS_EXPECT_TRUE(frame.resources.template get<"warp_a">().noise !=
                 frame.resources.template get<"warp_b">().noise);
  HS_EXPECT_NE(frame.resources.template get<"warp_a">().noise->GetNoise(
                   0.3f, 0.4f, 0.5f),
               frame.resources.template get<"warp_b">().noise->GetNoise(
                   0.3f, 0.4f, 0.5f));
  HS_EXPECT_EQ(effect.getParameters().size(),
               effect.getParameters().capacity());
  for (const char *name :
       {"warp_a.Planar Warp Speed", "warp_b.Planar Warp Speed",
        "warp_c.Planar Warp Speed", "lens_a.Mobius B Re", "lens_b.Mobius B Re"})
    HS_EXPECT_TRUE(effect.getParameters().find(name) != nullptr);
  HS_EXPECT_EQ(effect.updateParameter("lens_a.Mobius B Re", 0.2f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(
      effect.serialize_parameters().params.template get<"lens_a">().mobius.b.re,
      0.2f);
  HS_EXPECT_EQ(
      effect.serialize_parameters().params.template get<"lens_b">().mobius.b.re,
      0.0f);
  HS_EXPECT_EQ(effect.updateParameter("lens_b.Mobius B Re", -0.3f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(
      effect.serialize_parameters().params.template get<"lens_a">().mobius.b.re,
      0.2f);
  HS_EXPECT_EQ(
      effect.serialize_parameters().params.template get<"lens_b">().mobius.b.re,
      -0.3f);
  auto changed = captured;
  changed.params.template get<"lens_a">().mobius.b.re = 0.2f;
  changed.params.template get<"lens_b">().mobius.b.re = -0.3f;
  const auto middle =
      Pullback::interpolate(captured.params, changed.params, 0.5f);
  HS_EXPECT_NEAR(middle.template get<"lens_a">().mobius.b.re, 0.0f, 1e-6f);
  HS_EXPECT_NEAR(middle.template get<"lens_b">().mobius.b.re, 0.0f, 1e-6f);
  HS_EXPECT_TRUE(FX::valid_params(middle));
  verify_params_equal(
      Pullback::interpolate(captured.params, changed.params, 1.0f),
      changed.params);
  const auto render = [&](const typename FX::Params &params) {
    ComposedFrameWhiteBox::set_params(effect, params);
    const auto prepared =
        FX::RenderPipeline::prepare(ComposedFrameWhiteBox::frame(effect));
    std::array<Color4, 16> pixels{};
    for (size_t i = 0; i < pixels.size(); ++i)
      pixels[i] = FX::RenderPipeline::shade(
          math::Vector(0.2f + 0.03f * i, 0.4f, -0.8f).normalized(), prepared);
    return pixels;
  };
  const auto initial = render(captured.params);
  size_t differences = 0;
  const auto rendered = render(changed.params);
  for (size_t i = 0; i < initial.size(); ++i)
    differences += initial[i].color.r != rendered[i].color.r ||
                   initial[i].color.g != rendered[i].color.g ||
                   initial[i].color.b != rendered[i].color.b;
  HS_EXPECT_GT(differences, size_t{0});
  auto invalid = changed;
  invalid.params.template get<"lens_b">().mobius.c =
      invalid.params.template get<"lens_b">().mobius.a;
  invalid.params.template get<"lens_b">().mobius.d =
      invalid.params.template get<"lens_b">().mobius.b;
  HS_EXPECT_FALSE(effect.restore_parameters(invalid));
  verify_params_equal(effect.serialize_parameters().params, changed.params);
}

struct ProjectionWalkFootprint {
  size_t arena_bytes;
  int timeline_events;
};

template <bool AnimatedProjection>
ProjectionWalkFootprint measure_projection_walk_footprint() {
  reset_effect_globals();
  RippleProbe<SMALL_W, SMALL_H, AnimatedProjection> effect;
  effect.init();
  return {persistent_arena.get_offset(), Timeline::event_count()};
}

inline void test_composed_projection_walk_storage() {
  using Static = RippleProbe<SMALL_W, SMALL_H, false>;
  using Animated = RippleProbe<SMALL_W, SMALL_H, true>;
  static_assert(sizeof(Static) < sizeof(Animated));

  const ProjectionWalkFootprint static_footprint =
      measure_projection_walk_footprint<false>();
  const ProjectionWalkFootprint animated_footprint =
      measure_projection_walk_footprint<true>();
  HS_EXPECT_LT(static_footprint.arena_bytes, animated_footprint.arena_bytes);
  HS_EXPECT_EQ(static_footprint.timeline_events + 1,
               animated_footprint.timeline_events);
}

inline void test_composed_periodic_ripple_surface() {
  using FX = RippleProbe<SMALL_W, SMALL_H>;
  static_assert(!FX::HAS_SURFACE_NOISE);

  const auto RENDER = [](float strength) {
    reset_effect_globals();
    pin_frame_clock(0);
    FX effect;
    effect.init();
    HS_EXPECT_TRUE(effect.getParameters().find("Ripple Strength") != nullptr);
    HS_EXPECT_TRUE(effect.getParameters().find("Ripple Period") != nullptr);
    auto snapshot = effect.serialize_parameters();
    snapshot.params.template get<"surface">().strength = strength;
    HS_EXPECT_TRUE(effect.restore_parameters(snapshot));
    for (int frame = 0; frame < 20; ++frame) {
      pin_frame_clock(frame);
      effect.draw_frame();
      effect.advance_display();
    }
    std::vector<Pixel> pixels;
    capture_frame<SMALL_W, SMALL_H>(effect, pixels);
    return pixels;
  };
  const auto FLAT = RENDER(0.0f);
  const auto RIPPLE = RENDER(0.15f);
  size_t changed = 0;
  for (size_t i = 0; i < FLAT.size(); ++i)
    changed += FLAT[i].r != RIPPLE[i].r || FLAT[i].g != RIPPLE[i].g ||
               FLAT[i].b != RIPPLE[i].b;
  HS_EXPECT_GT(changed, size_t{0});
  hs::clear_mock_time();
}

template <typename SourceT> struct NoiseSourceProbeSpec : Pullback::Spec {
  static constexpr bool ANIMATED_PROJECTION = false;
  template <typename B>
  using Pipeline = Pullback::Pipeline<
      B, Pullback::Stage::Rotate<Pullback::OuterCameraProvider<B>>,
      Pullback::Stage::Placed<
          Pullback::CodeEmission::INLINE_ONLY,
          Pullback::Stage::Project<
              typename Pullback::ProjectionPolicyFor<PROJECTION, B>::Type>>,
      Pullback::Stage::Sample<
          typename Pullback::SourcePolicyFor<SourceT, B>::Type,
          Pullback::Weight::Projection, Pullback::ProjectionCoverage::Weight>,
      Pullback::Stage::Colorize<Pullback::Color::GeneratedPalette<
          Pullback::ColorProvider<B, HUE, BRIGHTNESS>>>>;
};

template <int W, int H, typename SourceT>
class NoiseSourceProbe
    : public Pullback::ComposedEffect<W, H, NoiseSourceProbe<W, H, SourceT>,
                                      NoiseSourceProbeSpec<SourceT>> {
public:
  using Params = Pullback::ParamsFor<NoiseSourceProbeSpec<SourceT>>;
  static constexpr std::array<std::string_view, 1> PRESET_IDS{"noise"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr int SOURCE_NOISE_SEED = 1337;

  static constexpr Params initial_params() {
    Params value;
    value.template get<"source">().noise_scale = 8.0f;
    value.template get<"source">().noise_contrast = 1.0f;
    value.template get<"source">().noise_time_rate = 1.0f / 128.0f;
    return value;
  }
};

template <int W, int H>
using NoHueProbe = NoiseSourceProbe<W, H, Pullback::ProjectedNoiseSourceParams>;

/**
 * @brief Both noise source families reach the derivation path.
 * @details Each source family implies its noise field; the plane-domain
 * and sphere-domain contours light the
 * frame and differ from each other on the same seed and parameters.
 */
inline void test_composed_noise_sources() {
  using Projected =
      NoiseSourceProbe<SMALL_W, SMALL_H, Pullback::ProjectedNoiseSourceParams>;
  using Spherical =
      NoiseSourceProbe<SMALL_W, SMALL_H, Pullback::SphericalNoiseSourceParams>;
  static_assert(Projected::HAS_SOURCE_NOISE);
  static_assert(Spherical::HAS_SOURCE_NOISE);

  std::array<Pixel, SMALL_W * SMALL_H> projected{};
  size_t lit = 0;
  {
    reset_effect_globals();
    Projected effect;
    effect.init();
    HS_EXPECT_TRUE(effect.getParameters().find("Source Noise Scale") !=
                   nullptr);
    HS_EXPECT_TRUE(effect.getParameters().find("Source Noise Speed") !=
                   nullptr);
    effect.draw_frame();
    effect.advance_display();
    for (int i = 0; i < SMALL_W * SMALL_H; ++i) {
      projected[i] = effect.display_buffer()[i];
      lit += projected[i].r != 0 || projected[i].g != 0 || projected[i].b != 0;
    }
  }
  HS_EXPECT_GT(lit, size_t(0));

  size_t spherical_lit = 0;
  size_t differing = 0;
  {
    reset_effect_globals();
    Spherical effect;
    effect.init();
    effect.draw_frame();
    effect.advance_display();
    for (int i = 0; i < SMALL_W * SMALL_H; ++i) {
      const Pixel &pixel = effect.display_buffer()[i];
      spherical_lit += pixel.r != 0 || pixel.g != 0 || pixel.b != 0;
      differing += pixel.r != projected[i].r || pixel.g != projected[i].g ||
                   pixel.b != projected[i].b;
    }
  }
  HS_EXPECT_GT(spherical_lit, size_t(0));
  HS_EXPECT_GT(differing, size_t(0));
}
/** @brief Parameter set of the ChoreographedEffect probes: one animated float. */
struct ChoreoProbeParams {
  float level = 0.0f;
};

/**
 * @brief Probe pinning a Segue::Preset::Fade departure.
 * @details Records every set_preset_opacity sample the base feeds it and the
 * parameter level while the opacity falls and rises, so the test can see the
 * departing preset hold until dark and the next one brighten in.
 */
template <int W, int H>
class FadeChoreoProbe
    : public ChoreographedEffect<FadeChoreoProbe<W, H>, ChoreoProbeParams> {
  using Choreography =
      ChoreographedEffect<FadeChoreoProbe<W, H>, ChoreoProbeParams>;
  friend Choreography;

public:
  using Params = ChoreoProbeParams;
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr Segue::Preset::Fade DEPARTURE{8};
  static constexpr uint16_t PRESET_DWELL_FRAMES = 12;
  static constexpr std::array<PresetEntry<Params>, 3> PRESETS = {
      {{{0.25f}, DEPARTURE}, {{0.5f}, DEPARTURE}, {{0.75f}, DEPARTURE}}};

  static bool valid_params(const Params &p) {
    return p.level >= 0.0f && p.level <= 1.0f;
  }

  FadeChoreoProbe() : Choreography(W, H) {}

  void init() override { this->begin_choreography(); }

  void draw_frame() override {
    Canvas canvas(*this);
    this->timeline.step(canvas);
    this->step_choreography();
  }

  /** @brief Live level, which a preset change snaps. */
  float level() const { return this->params.level; }

  int opacity_samples = 0;     /**< Opacity samples the base fed. */
  float min_opacity = 2.0f;    /**< Lowest opacity sample seen. */
  float last_opacity = -1.0f;  /**< Most recent opacity sample. */
  float dimming_level = -1.0f; /**< Level while the opacity last fell. */
  float rising_level = -1.0f;  /**< Level while the opacity last rose. */

private:
  void set_preset_opacity(float value) {
    ++opacity_samples;
    if (value < last_opacity && value > 0.0f)
      dimming_level = this->params.level;
    else if (value > last_opacity)
      rising_level = this->params.level;
    last_opacity = value;
    min_opacity = std::min(min_opacity, value);
  }
};

/**
 * @brief Probe pinning the Segue::Preset::Lerp transition hooks.
 * @details Counts transition_armed and blend_params calls and writes the blend
 * itself, so a cancelled crossfade shows up as a blend that stops writing.
 */
template <int W, int H, bool Pausable = false, bool Animated = true>
class LerpChoreoProbe
    : public ChoreographedEffect<LerpChoreoProbe<W, H, Pausable, Animated>,
                                 ChoreoProbeParams> {
  using Choreography =
      ChoreographedEffect<LerpChoreoProbe<W, H, Pausable, Animated>,
                          ChoreoProbeParams>;
  friend Choreography;

public:
  using Params = ChoreoProbeParams;
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr Segue::Preset::Lerp DEPARTURE{8, math::ease_linear,
                                                 Pausable};
  static constexpr uint16_t PRESET_DWELL_FRAMES = 12;
  static constexpr std::array<PresetEntry<Params>, 2> PRESETS = {
      {{{0.0f}, DEPARTURE}, {{1.0f}, DEPARTURE}}};

  static bool valid_params(const Params &p) {
    return p.level >= 0.0f && p.level <= 1.0f;
  }

  LerpChoreoProbe() : Choreography(W, H) {}

  void init() override {
    if constexpr (Animated)
      this->register_animated_param("Level", &this->params.level, 0.0f, 1.0f);
    else
      this->register_param("Level", &this->params.level, 0.0f, 1.0f);
    this->begin_choreography();
  }

  void draw_frame() override {
    Canvas canvas(*this);
    this->timeline.step(canvas);
    this->step_choreography();
  }

  /** @brief Live level, which the blend rewrites every frame it runs. */
  float level() const { return this->params.level; }

  int armed_count = 0;        /**< transition_armed calls seen. */
  float armed_target = -1.0f; /**< Level the last arming named. */
  int blend_calls = 0;        /**< blend_params calls seen. */

private:
  void transition_armed(const Params &target) {
    ++armed_count;
    armed_target = target.level;
  }

  void blend_params(float progress) {
    ++blend_calls;
    this->params.level = hs::lerp(this->transition.from.level,
                                  this->transition.to.level, progress);
  }
};

/** @brief Runs @p frames whole frames of @p effect, buffer swap included. */
template <typename FX> void run_probe_frames(FX &effect, int frames) {
  for (int frame = 0; frame < frames; ++frame) {
    effect.draw_frame();
    effect.advance_display();
  }
}

/**
 * @brief Pins a Segue::Preset::Fade departure.
 * @details After the dwell the departing preset holds while the opacity falls,
 * the next preset is adopted in the dark and holds while it rises back to
 * full, the cadence wraps the table, and a pause shows full opacity.
 */
inline void test_choreography_fade_departure() {
  using FX = FadeChoreoProbe<SMALL_W, SMALL_H>;
  constexpr int DWELL = FX::PRESET_DWELL_FRAMES;
  constexpr int FADE = FX::DEPARTURE.frames;

  reset_effect_globals();
  FX effect;
  effect.init();
  HS_EXPECT_EQ(effect.getPresetCount(), size_t{3});
  run_probe_frames(effect, DWELL);
  HS_EXPECT_EQ(effect.opacity_samples, 0);
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{0});
  HS_EXPECT_EQ(effect.level(), 0.25f);

  run_probe_frames(effect, FADE + 1);
  HS_EXPECT_EQ(effect.opacity_samples, FADE);
  HS_EXPECT_EQ(effect.dimming_level, 0.25f);
  HS_EXPECT_EQ(effect.rising_level, 0.5f);
  HS_EXPECT_LE(effect.min_opacity, 1e-3f);
  HS_EXPECT_EQ(effect.last_opacity, 1.0f);
  HS_EXPECT_EQ(effect.level(), 0.5f);

  // The cadence wraps the table back to the first preset.
  int frames = 0;
  while (effect.getPresetIndex() != 0 && frames < 4 * (DWELL + FADE)) {
    run_probe_frames(effect, 1);
    ++frames;
  }
  HS_EXPECT_EQ(effect.getPresetIndex(), size_t{0});
  run_probe_frames(effect, FADE + 1);
  HS_EXPECT_EQ(effect.level(), 0.25f);

  // Paused mid-fade: full opacity, and no advance while paused.
  HS_EXPECT_TRUE(effect.selectPreset(0));
  effect.setAnimationsPaused(false);
  run_probe_frames(effect, DWELL + 2);
  HS_EXPECT_LT(effect.last_opacity, 1.0f);
  const size_t HELD = effect.getPresetIndex();
  effect.setAnimationsPaused(true);
  run_probe_frames(effect, 2 * (DWELL + FADE));
  HS_EXPECT_EQ(effect.last_opacity, 1.0f);
  HS_EXPECT_EQ(effect.getPresetIndex(), HELD);
  HS_EXPECT_TRUE(effect.selectPreset(0));
  run_probe_frames(effect, 1);
  HS_EXPECT_EQ(effect.last_opacity, 1.0f);
  HS_EXPECT_EQ(effect.level(), 0.25f);
}

/**
 * @brief Pins the two Segue::Preset::Lerp transition hooks.
 * @details transition_armed fires once per automatic crossfade, with the
 * incoming preset; parameter_written ends the crossfade in flight, so
 * the blend stops rewriting the value the write just landed.
 */
inline void test_choreography_lerp_transition_hooks() {
  const auto check = []<bool Animated>() {
    using FX = LerpChoreoProbe<SMALL_W, SMALL_H, false, Animated>;
    reset_effect_globals();
    FX effect;
    effect.init();
    HS_EXPECT_EQ(effect.armed_count, 0);
    HS_EXPECT_EQ(effect.level(), 0.0f);

    // Retire the dwell; the next frame arms the crossfade to PRESETS[1].
    run_probe_frames(effect, FX::PRESET_DWELL_FRAMES);
    HS_EXPECT_EQ(effect.armed_count, 1);
    HS_EXPECT_EQ(effect.armed_target, 1.0f);
    HS_EXPECT_EQ(effect.getPresetIndex(), size_t{1});

    // Mid-crossfade the level is strictly between the endpoints.
    run_probe_frames(effect, FX::DEPARTURE.frames / 2);
    HS_EXPECT_GT(effect.blend_calls, 0);
    HS_EXPECT_GT(effect.level(), 0.0f);
    HS_EXPECT_LT(effect.level(), 1.0f);

    // The manual write cancels the crossfade: the lerp keeps stepping (the policy
    // is unpausable) but blend_params is never called again.
    const int blends = effect.blend_calls;
    HS_EXPECT_EQ(effect.updateParameter("Level", 0.3f),
                 ParamSetResult::APPLIED);
    run_probe_frames(effect, Animated ? 2 * FX::DEPARTURE.frames
                                      : FX::DEPARTURE.frames / 2);
    HS_EXPECT_EQ(effect.blend_calls, blends);
    HS_EXPECT_EQ(effect.level(), 0.3f);
    HS_EXPECT_EQ(effect.armed_count, 1);
    HS_EXPECT_EQ(effect.animations_paused(), Animated);
  };
  check.template operator()<true>();
  check.template operator()<false>();
}

inline void test_choreography_lerp_pause_policy() {
  const auto check = []<bool Pausable>() {
    using FX = LerpChoreoProbe<SMALL_W, SMALL_H, Pausable>;
    reset_effect_globals();
    FX effect;
    effect.init();
    run_probe_frames(effect, FX::PRESET_DWELL_FRAMES + 3);
    const float held = effect.level();
    HS_EXPECT_GT(held, 0.0f);
    HS_EXPECT_LT(held, 1.0f);
    effect.setAnimationsPaused(true);
    run_probe_frames(effect, FX::DEPARTURE.frames);
    HS_EXPECT_EQ(effect.level(), Pausable ? held : 1.0f);
    effect.setAnimationsPaused(false);
    run_probe_frames(effect, FX::DEPARTURE.frames);
    HS_EXPECT_EQ(effect.level(), 1.0f);
  };
  check.template operator()<true>();
  check.template operator()<false>();
}
