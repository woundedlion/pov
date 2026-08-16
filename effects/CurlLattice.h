/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <span>
#include <string_view>

#include "core/color/effect_palette_recipes.h"
#include "core/engine/engine.h"
#include "core/engine/fixed_pipeline.h"
#include "core/math/noise_field.h"
#include "core/render/pullback.h"

namespace hs_test {
namespace curl_lattice_tests {
struct CurlLatticeWhiteBox;
} // namespace curl_lattice_tests
} // namespace hs_test

/**
 * @brief Fixed folded-sinusoidal lattice displaced by sphere-space curl noise.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H> class CurlLattice : public Effect {
public:
  static constexpr std::string_view EFFECT_ID = "curl-lattice";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "487dc32c5a96a54ab9a01daa6a6d9b95bcb96fa6a79676f3a4f8fc28a303fcd3";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "9b4b0152cd159b90bd3663505f53c8b6bdbcce846a37594135eff1e46aa8eaa3";
  static constexpr std::array<std::string_view, 2> PRESET_IDS{"open-curl",
                                                              "dense-curl"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 3;
  static constexpr uint16_t TRANSITION_DURATION = 480;

  struct LatticeParams {
    float lattice_cell_scale = 0.710265636f;
    float lattice_shape_blend = 1.0f;
    float lattice_softness = 0.455532223f;
    float lattice_radius = 0.290762514f;
  };

  struct ProjectionParams {
    float pole_fade = 20.0f;
    float central_meridian = 0.0f;
    float spin_speed = 0.0f;
    float wander = 1.0f;
    float camera_wander = 1.0f;
  };

  struct SurfaceNoiseParams {
    float scale = 1.78815627f;
    float strength = 0.0759999976f;
    float speed = 0.0f;
  };

  struct ColorParams {
    float palette_chroma = 1.0f;
    float mapping_frequency = 1.0f;
    float mapping_phase = -0.165999994f;
    float phase_oscillation_depth = 0.0f;
    float phase_oscillation_speed = 0.0f;
    float brightness_depth = 1.0f;
    float opacity_low = 1.0f;
    float opacity_high = 1.0f;
    float hue_shift_amount = 0.268000007f;
    float hue_noise_scale = 2.0f;
    float hue_noise_speed = 0.0f;
    Pullback::Color::PaletteMapping palette_mapping =
        Pullback::Color::PaletteMapping::CUP;
  };

  struct Params {
    LatticeParams lattice;
    ProjectionParams projection;
    SurfaceNoiseParams surface_noise;
    ColorParams color;
  };

  struct ParameterSnapshot {
    uint32_t schema_version;
    Params params;
  };

  HS_COLD_MEMBER CurlLattice() : Effect(W, H, {.strobe = true}) {}

  HS_COLD_MEMBER void init() override {
    configure_presets(PRESETS.size());
    state = persistent_arena.make<State>();
    use_parameter_storage(persistent_arena.allocate_n<ParamDef>(PARAM_CAPACITY),
                          PARAM_CAPACITY);

    configure_noise(state->surface_noise, SURFACE_NOISE_SEED);
    configure_noise(state->color_noise, COLOR_NOISE_SEED);
    prepare_hue_noise_lut();

    init_gamut_lut(persistent_arena, GAMUT_ANGLE_STEPS, GAMUT_L_STEPS);
    palette_cycler.init_generated(persistent_arena, next_palette, this, 0, 600,
                                  ease_in_out_sin);

    timeline.add(0, Animation::RandomWalk<W>(projection_walk, UP,
                                             state->projection_walk_noise));
    timeline.add(
        0, Animation::RandomWalk<W>(outer_walk, UP, state->outer_walk_noise));

    register_animated_param("Lattice Cell Scale",
                            &params.lattice.lattice_cell_scale, CELL_MIN,
                            CELL_MAX);
    register_animated_param("Lattice Shape",
                            &params.lattice.lattice_shape_blend, 0.0f, 1.0f);
    register_animated_param("Lattice Softness",
                            &params.lattice.lattice_softness, SOFTNESS_MIN,
                            1.0f);
    register_animated_param("Lattice Radius", &params.lattice.lattice_radius,
                            RADIUS_MIN, RADIUS_MAX);
    register_animated_param("Pole Fade", &params.projection.pole_fade,
                            POLE_FADE_MIN, POLE_FADE_MAX);
    register_animated_param("Central Meridian",
                            &params.projection.central_meridian, 0.0f,
                            TWO_PI_F);
    register_animated_param("Projection Spin Speed",
                            &params.projection.spin_speed, SPIN_SPEED_MIN,
                            SPIN_SPEED_MAX);
    register_animated_param("Projection Wander", &params.projection.wander,
                            0.0f, 1.0f);
    register_animated_param("Camera Wander", &params.projection.camera_wander,
                            0.0f, 1.0f);
    register_animated_param("Surface Noise Scale", &params.surface_noise.scale,
                            SURFACE_SCALE_MIN, SURFACE_SCALE_MAX);
    register_animated_param("Surface Noise Strength",
                            &params.surface_noise.strength,
                            SURFACE_STRENGTH_MIN, SURFACE_STRENGTH_MAX);
    register_animated_param("Surface Noise Speed", &params.surface_noise.speed,
                            -NOISE_SPEED_MAX, NOISE_SPEED_MAX);
    register_animated_param("Palette Chroma", &params.color.palette_chroma,
                            0.0f, 1.0f);
    register_animated_param("Palette Mapping", &params.color.palette_mapping,
                            PALETTE_MAPPING_OPTIONS,
                            PALETTE_MAPPING_EXPORT_OPTIONS,
                            std::size(PALETTE_MAPPING_OPTIONS));
    register_animated_param("Mapping Frequency",
                            &params.color.mapping_frequency,
                            MAPPING_FREQUENCY_MIN, MAPPING_FREQUENCY_MAX);
    register_animated_param("Mapping Phase", &params.color.mapping_phase, -1.0f,
                            1.0f);
    register_animated_param("Phase Oscillation Depth",
                            &params.color.phase_oscillation_depth, 0.0f, 1.0f);
    register_animated_param(
        "Phase Oscillation Speed", &params.color.phase_oscillation_speed,
        -PHASE_OSCILLATION_SPEED_MAX, PHASE_OSCILLATION_SPEED_MAX);
    register_animated_param("Brightness Depth", &params.color.brightness_depth,
                            0.0f, 1.0f);
    register_animated_param("Value Opacity Low", &params.color.opacity_low,
                            0.0f, 1.0f);
    register_animated_param("Value Opacity High", &params.color.opacity_high,
                            0.0f, 1.0f);
    register_animated_param("Hue Shift Amount", &params.color.hue_shift_amount,
                            -HUE_SHIFT_AMOUNT_MAX, HUE_SHIFT_AMOUNT_MAX);
    register_animated_param("Hue Noise Scale", &params.color.hue_noise_scale,
                            HUE_NOISE_SCALE_MIN, HUE_NOISE_SCALE_MAX);
    register_animated_param("Hue Noise Speed", &params.color.hue_noise_speed,
                            -HUE_NOISE_SPEED_MAX, HUE_NOISE_SPEED_MAX);
  }

  HS_FLASH_MEMBER void draw_frame() override {
    Canvas canvas(*this);
    {
      HS_PROFILE(cl_timeline_step);
      timeline.step(canvas);
    }
    {
      HS_PROFILE(cl_advance);
      begin_automatic_transition();
      prepare_transition_value();
      advance_runtime();
      update_spatial_frames();
      update_palette_chroma();
      palette_cycler.step();
    }
    const FrameState frame = prepare_frame();
    {
      HS_PROFILE(cl_shader_draw);
      Scan::Shader::draw_cached<W, H, 1>(canvas, [&frame](const Vector &view) {
        return RenderPipeline::shade(view, frame);
      });
    }
    finish_transition_evaluation();
  }

  ParameterSnapshot serialize_parameters() const {
    return {PARAMETER_SCHEMA_VERSION, params};
  }

  bool restore_parameters(const ParameterSnapshot &snapshot) {
    if (snapshot.schema_version != PARAMETER_SCHEMA_VERSION ||
        !params_in_range(snapshot.params))
      return false;
    transition.active = false;
    params = snapshot.params;
    palette_mapping = Pullback::Color::PaletteMappingWeights::single(
        params.color.palette_mapping);
    return true;
  }

private:
  friend struct ::hs_test::curl_lattice_tests::CurlLatticeWhiteBox;

  static constexpr int GAMUT_ANGLE_STEPS = GAMUT_LUT_ANGLE_STEPS;
  static constexpr int GAMUT_L_STEPS = GAMUT_LUT_L_STEPS;
  static constexpr int32_t SURFACE_NOISE_SEED = 1337;
  static constexpr int32_t COLOR_NOISE_SEED = 6047;
  static constexpr float CELL_MIN = 1.0f / 64.0f;
  static constexpr float CELL_MAX = 8.0f;
  static constexpr float SOFTNESS_MIN = 1.0f / 1024.0f;
  static constexpr float RADIUS_MIN = 1.0f / 64.0f;
  static constexpr float RADIUS_MAX = 0.49f;
  static constexpr float POLE_FADE_MIN = 1.0f;
  static constexpr float POLE_FADE_MAX = 20.0f;
  static constexpr float SPIN_SPEED_MIN = 0.0f;
  static constexpr float SPIN_SPEED_MAX = 0.05f;
  static constexpr float SURFACE_SCALE_MIN = 1.0f / 64.0f;
  static constexpr float SURFACE_SCALE_MAX = 64.0f;
  static constexpr float SURFACE_STRENGTH_MIN = -0.5f;
  static constexpr float SURFACE_STRENGTH_MAX = 0.5f;
  static constexpr float NOISE_SPEED_MAX = 1.0f / 64.0f;
  static constexpr float MAPPING_FREQUENCY_MIN = 1.0f;
  static constexpr float MAPPING_FREQUENCY_MAX = 32.0f;
  static constexpr float PHASE_OSCILLATION_SPEED_MAX = 0.01f;
  static constexpr float HUE_SHIFT_AMOUNT_MAX = 1.0f;
  static constexpr float HUE_NOISE_SCALE_MIN = 1.0f / 64.0f;
  static constexpr float HUE_NOISE_SCALE_MAX = 8.0f;
  static constexpr float HUE_NOISE_SPEED_MAX = 0.001f;
  static constexpr uint32_t HUE_STEP = 159;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;
  static constexpr size_t PARAM_CAPACITY = 48;
  static constexpr const char *PALETTE_MAPPING_OPTIONS[] = {
      "Cup", "Bell", "Linear", "Reverse"};
  static constexpr const char *PALETTE_MAPPING_EXPORT_OPTIONS[] = {
      "Pullback::Color::PaletteMapping::CUP",
      "Pullback::Color::PaletteMapping::BELL",
      "Pullback::Color::PaletteMapping::LINEAR",
      "Pullback::Color::PaletteMapping::REVERSE"};

  struct Preset {
    std::string_view id;
    std::string_view display_name;
    Params params;
  };

  static constexpr std::array<Preset, 2> PRESETS{{
      {PRESET_IDS[0], "Open Curl", {{}, {}, {1.78815627f}, {}}},
      {PRESET_IDS[1], "Dense Curl", {{}, {}, {3.29720306f}, {}}},
  }};

  struct Transition {
    Params from{};
    Params to{};
    Pullback::Color::PaletteMappingWeights mapping_from{};
    Pullback::Color::PaletteMappingWeights mapping_to{};
    uint16_t evaluation = 0;
    uint16_t duration = TRANSITION_DURATION;
    bool active = false;
  };

  struct FrameState {
    Quaternion projection_conjugate;
    Quaternion outer_conjugate;
    Vector surface_loop_offset;
    const FastNoiseLite *surface_noise;
    const BakedPalette *palette;
    const Pixel *hue_rotation_lut;
    const int8_t *hue_noise_lut;
    Params params;
    Pullback::Color::PaletteMappingWeights palette_mapping;
    float palette_oscillation_phase;
  };

  struct Binding {
    using FrameState = CurlLattice::FrameState;
    using Instrumentation = Pullback::NoInstrumentation;
  };

  struct OuterCameraStateProvider {
    using Binding = CurlLattice::Binding;
    using FrameState = typename Binding::FrameState;

    static const Quaternion &conjugate(const FrameState &frame) {
      return frame.outer_conjugate;
    }
  };

  struct ProjectionStateProvider {
    using Binding = CurlLattice::Binding;
    using FrameState = typename Binding::FrameState;

    static const Quaternion &conjugate(const FrameState &frame) {
      return frame.projection_conjugate;
    }
    static float central_meridian(const FrameState &frame) {
      return frame.params.projection.central_meridian;
    }
    static float pole_fade(const FrameState &frame) {
      return frame.params.projection.pole_fade;
    }
  };

  struct SurfaceStateProvider {
    using Binding = CurlLattice::Binding;
    using FrameState = typename Binding::FrameState;

    static const FastNoiseLite &noise(const FrameState &frame) {
      return *frame.surface_noise;
    }
    static float scale(const FrameState &frame) {
      return frame.params.surface_noise.scale;
    }
    static const Vector &loop_offset(const FrameState &frame) {
      return frame.surface_loop_offset;
    }
    static float strength(const FrameState &frame) {
      return frame.params.surface_noise.strength;
    }
    static bool path_length_required(const FrameState &) { return false; }
  };

  struct SourceStateProvider {
    using Binding = CurlLattice::Binding;
    using FrameState = typename Binding::FrameState;

    static const LatticeParams &params(const FrameState &frame) {
      return frame.params.lattice;
    }
  };

  struct ColorStateProvider {
    using Binding = CurlLattice::Binding;
    using FrameState = typename Binding::FrameState;

    static Pullback::Color::PaletteMappingWeights
    mapping_weights(const FrameState &frame) {
      return frame.palette_mapping;
    }
    static float mapping_frequency(const FrameState &frame) {
      return frame.params.color.mapping_frequency;
    }
    static float mapping_phase(const FrameState &frame) {
      return frame.params.color.mapping_phase;
    }
    static float oscillation_depth(const FrameState &frame) {
      return frame.params.color.phase_oscillation_depth;
    }
    static float oscillation_phase(const FrameState &frame) {
      return frame.palette_oscillation_phase;
    }
    static const BakedPalette &palette(const FrameState &frame) {
      return *frame.palette;
    }
    static Pullback::Color::HueMode hue_mode(const FrameState &) {
      return Pullback::Color::HueMode::NOISE;
    }
    static float hue_shift_amount(const FrameState &frame) {
      return frame.params.color.hue_shift_amount;
    }
    static Pullback::Color::HueRotationLutView
    hue_rotation(const FrameState &frame) {
      return {frame.hue_rotation_lut, true};
    }
    static Pullback::Color::HueNoiseLutView hue_noise(const FrameState &frame) {
      return {frame.hue_noise_lut, true};
    }
    static Pullback::Color::BrightnessEnvelope
    brightness_envelope(const FrameState &) {
      return Pullback::Color::BrightnessEnvelope::CUP;
    }
    static float brightness_depth(const FrameState &frame) {
      return frame.params.color.brightness_depth;
    }
    static float opacity_low(const FrameState &frame) {
      return frame.params.color.opacity_low;
    }
    static float opacity_high(const FrameState &frame) {
      return frame.params.color.opacity_high;
    }
  };

  using OuterCameraStage =
      Pullback::Stage::OuterCamera<Binding, OuterCameraStateProvider>;
  using SurfaceImplementation = Pullback::Stage::SurfaceProject<
      Binding,
      Pullback::Surface::CurlNoise<SurfaceStateProvider, NoiseBasis::SIMPLEX,
                                   Pullback::Surface::Euler>,
      Pullback::Lens::Identity, Pullback::Surface::Identity,
      Pullback::Projection::FoldedSinusoidal<ProjectionStateProvider>>;
  using SurfaceStage =
      Pullback::Stage::Placed<Pullback::CodeEmission::OUT_OF_LINE_FLASH,
                              SurfaceImplementation>;
  using PlanarWarpStage =
      Pullback::Stage::PlanarWarp<Binding, Pullback::Warp::Identity,
                                  Pullback::Warp::Identity>;
  using SourceStage = Pullback::Stage::Source<
      Binding, Pullback::Source::PrimitiveLattice<SourceStateProvider>>;
  using MaterialStage =
      Pullback::Stage::Material<Binding, Pullback::Weight::Projection,
                                Pullback::Transfer::Linear,
                                Pullback::Coverage::Projection>;
  using ColorStage = Pullback::Stage::Color<
      Binding, Pullback::Color::GeneratedPalette<ColorStateProvider>>;
  using RenderPipeline =
      Pullback::Pipeline<Binding, OuterCameraStage, SurfaceStage,
                         PlanarWarpStage, SourceStage, MaterialStage,
                         ColorStage>;

  struct State {
    std::array<Pixel, Pullback::Color::HueRotationLutView::SIZE>
        hue_rotation_lut;
    std::array<int8_t, Pullback::Color::HueNoiseLutView::SIZE> hue_noise_lut;
    FastNoiseLite surface_noise;
    FastNoiseLite color_noise;
    FastNoiseLite projection_walk_noise;
    FastNoiseLite outer_walk_noise;
    float hue_noise_lut_scale = -1.0f;
    float hue_noise_lut_phase = -1.0f;
  };

  static bool in_range(float value, float minimum, float maximum) {
    return std::isfinite(value) && value >= minimum && value <= maximum;
  }

  static bool params_in_range(const Params &value) {
    static_assert(sizeof(Params) == 24 * sizeof(float));
    return in_range(value.lattice.lattice_cell_scale, CELL_MIN, CELL_MAX) &&
           in_range(value.lattice.lattice_shape_blend, 0.0f, 1.0f) &&
           in_range(value.lattice.lattice_softness, SOFTNESS_MIN, 1.0f) &&
           in_range(value.lattice.lattice_radius, RADIUS_MIN, RADIUS_MAX) &&
           in_range(value.projection.pole_fade, POLE_FADE_MIN, POLE_FADE_MAX) &&
           in_range(value.projection.central_meridian, 0.0f, TWO_PI_F) &&
           in_range(value.projection.spin_speed, SPIN_SPEED_MIN,
                    SPIN_SPEED_MAX) &&
           in_range(value.projection.wander, 0.0f, 1.0f) &&
           in_range(value.projection.camera_wander, 0.0f, 1.0f) &&
           in_range(value.surface_noise.scale, SURFACE_SCALE_MIN,
                    SURFACE_SCALE_MAX) &&
           in_range(value.surface_noise.strength, SURFACE_STRENGTH_MIN,
                    SURFACE_STRENGTH_MAX) &&
           in_range(value.surface_noise.speed, -NOISE_SPEED_MAX,
                    NOISE_SPEED_MAX) &&
           in_range(value.color.palette_chroma, 0.0f, 1.0f) &&
           static_cast<uint8_t>(value.color.palette_mapping) <=
               static_cast<uint8_t>(Pullback::Color::PaletteMapping::REVERSE) &&
           in_range(value.color.mapping_frequency, MAPPING_FREQUENCY_MIN,
                    MAPPING_FREQUENCY_MAX) &&
           in_range(value.color.mapping_phase, -1.0f, 1.0f) &&
           in_range(value.color.phase_oscillation_depth, 0.0f, 1.0f) &&
           in_range(value.color.phase_oscillation_speed,
                    -PHASE_OSCILLATION_SPEED_MAX,
                    PHASE_OSCILLATION_SPEED_MAX) &&
           in_range(value.color.brightness_depth, 0.0f, 1.0f) &&
           in_range(value.color.opacity_low, 0.0f, 1.0f) &&
           in_range(value.color.opacity_high, 0.0f, 1.0f) &&
           in_range(value.color.hue_shift_amount, -HUE_SHIFT_AMOUNT_MAX,
                    HUE_SHIFT_AMOUNT_MAX) &&
           in_range(value.color.hue_noise_scale, HUE_NOISE_SCALE_MIN,
                    HUE_NOISE_SCALE_MAX) &&
           in_range(value.color.hue_noise_speed, -HUE_NOISE_SPEED_MAX,
                    HUE_NOISE_SPEED_MAX);
  }

  static Params interpolate_params(const Params &from, const Params &to,
                                   float progress) {
    Params value;
#define HS_CURL_LINEAR(path)                                                   \
  value.path = FixedPipeline::linear(from.path, to.path, progress)
    HS_CURL_LINEAR(lattice.lattice_cell_scale);
    HS_CURL_LINEAR(lattice.lattice_shape_blend);
    HS_CURL_LINEAR(lattice.lattice_softness);
    HS_CURL_LINEAR(lattice.lattice_radius);
    HS_CURL_LINEAR(projection.pole_fade);
    value.projection.central_meridian = FixedPipeline::shortest_periodic(
        from.projection.central_meridian, to.projection.central_meridian,
        progress, TWO_PI_F);
    HS_CURL_LINEAR(projection.spin_speed);
    HS_CURL_LINEAR(projection.wander);
    HS_CURL_LINEAR(projection.camera_wander);
    HS_CURL_LINEAR(surface_noise.scale);
    HS_CURL_LINEAR(surface_noise.strength);
    HS_CURL_LINEAR(surface_noise.speed);
    HS_CURL_LINEAR(color.palette_chroma);
    value.color.palette_mapping =
        progress < 1.0f ? from.color.palette_mapping : to.color.palette_mapping;
    HS_CURL_LINEAR(color.mapping_frequency);
    HS_CURL_LINEAR(color.mapping_phase);
    HS_CURL_LINEAR(color.phase_oscillation_depth);
    HS_CURL_LINEAR(color.phase_oscillation_speed);
    HS_CURL_LINEAR(color.brightness_depth);
    HS_CURL_LINEAR(color.opacity_low);
    HS_CURL_LINEAR(color.opacity_high);
    HS_CURL_LINEAR(color.hue_shift_amount);
    HS_CURL_LINEAR(color.hue_noise_scale);
    HS_CURL_LINEAR(color.hue_noise_speed);
#undef HS_CURL_LINEAR
    return value;
  }

  HS_COLD_MEMBER bool apply_preset(const PresetChange &change) override {
    const Params target = PRESETS[change.to].params;
    if (change.origin != PresetChangeOrigin::AUTOMATIC) {
      transition.active = false;
      params = target;
      palette_mapping = Pullback::Color::PaletteMappingWeights::single(
          target.color.palette_mapping);
      preset_dwell_remaining = PRESET_DWELL_FRAMES;
      return true;
    }
    transition = {params,
                  target,
                  palette_mapping,
                  Pullback::Color::PaletteMappingWeights::single(
                      target.color.palette_mapping),
                  0,
                  TRANSITION_DURATION,
                  true};
    return true;
  }

  static void configure_noise(FastNoiseLite &noise, int32_t seed) {
    noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    noise.SetSeed(seed);
    noise.SetFrequency(1.0f);
  }

  HS_COLD_MEMBER void prepare_hue_noise_lut() {
    Pullback::Color::prepare_hue_noise_lut(
        std::span<int8_t, Pullback::Color::HueNoiseLutView::SIZE>(
            state->hue_noise_lut),
        state->color_noise, params.color.hue_noise_scale, hue_noise_phase);
    state->hue_noise_lut_scale = params.color.hue_noise_scale;
    state->hue_noise_lut_phase = hue_noise_phase;
  }

  HS_COLD_MEMBER void prepare_hue_rotation_lut() const {
    Pullback::Color::prepare_hue_rotation_lut(
        std::span<Pixel, Pullback::Color::HueRotationLutView::SIZE>(
            state->hue_rotation_lut),
        palette_cycler.palette());
  }

  HS_COLD_MEMBER void begin_automatic_transition() {
    if (anims_paused || transition.active)
      return;
    if (preset_dwell_remaining > 0 && --preset_dwell_remaining > 0)
      return;
    HS_CHECK(advancePreset());
#ifdef HS_PROFILE_ENABLE
    hs::log("Preset: %u/%u", static_cast<unsigned>(getPresetIndex() + 1),
            static_cast<unsigned>(getPresetCount()));
#endif
  }

  HS_COLD_MEMBER void prepare_transition_value() {
    if (!transition.active)
      return;
    const FixedPipeline::EdgeProgress progress =
        FixedPipeline::edge_progress(transition.evaluation, transition.duration,
                                     FixedPipeline::Easing::EASE_IN_OUT_SIN);
    params = interpolate_params(transition.from, transition.to, progress.eased);
    palette_mapping = Pullback::Color::PaletteMappingWeights::lerp(
        transition.mapping_from, transition.mapping_to, progress.eased);
  }

  HS_COLD_MEMBER void finish_transition_evaluation() {
    if (!transition.active || anims_paused)
      return;
    if (transition.evaluation == transition.duration) {
      params = transition.to;
      palette_mapping = transition.mapping_to;
      transition.active = false;
      preset_dwell_remaining = PRESET_DWELL_FRAMES;
      return;
    }
    ++transition.evaluation;
  }

  HS_COLD_MEMBER void update_spatial_frames() {
    const Quaternion projection = projection_walk.get();
    const Quaternion projection_delta =
        projection * projection_walk_previous.conjugate();
    projection_walk_previous = projection;
    projection_wander =
        (FixedPipeline::scaled_rotation_delta(projection_delta.normalized(),
                                              params.projection.wander) *
         projection_wander)
            .normalized();

    const Quaternion outer = outer_walk.get();
    const Quaternion outer_delta = outer * outer_walk_previous.conjugate();
    outer_walk_previous = outer;
    outer_wander =
        (FixedPipeline::scaled_rotation_delta(outer_delta.normalized(),
                                              params.projection.camera_wander) *
         outer_wander)
            .normalized();

    projection_conjugate = (make_rotation(Y_AXIS, projection_spin) *
                            base_orientation * projection_wander)
                               .conjugate();
    outer_conjugate = outer_wander.conjugate();
  }

  HS_COLD_MEMBER void advance_runtime() {
    projection_spin =
        fmodf(projection_spin + params.projection.spin_speed, TWO_PI_F);
    surface_noise_phase =
        wrap_t(surface_noise_phase + params.surface_noise.speed);
    hue_noise_phase = wrap_t(hue_noise_phase + params.color.hue_noise_speed);
    palette_oscillation_phase = wrap_t(palette_oscillation_phase +
                                       params.color.phase_oscillation_speed);
  }

  HS_COLD_MEMBER void update_palette_chroma() {
    if (palette_chroma == params.color.palette_chroma)
      return;
    palette_chroma = params.color.palette_chroma;
    palette_cycler.set_generated_chroma(palette_chroma);
  }

  HS_COLD_MEMBER FrameState prepare_frame() {
    HS_PROFILE(cl_prepare_frame);
    if (state->hue_noise_lut_scale != params.color.hue_noise_scale ||
        state->hue_noise_lut_phase != hue_noise_phase)
      prepare_hue_noise_lut();
    prepare_hue_rotation_lut();
    const float surface_angle = TWO_PI_F * wrap_t(surface_noise_phase);
    return {projection_conjugate,
            outer_conjugate,
            Vector(NOISE_LOOP_RADIUS * cosf(surface_angle),
                   NOISE_LOOP_RADIUS * sinf(surface_angle), 0.0f),
            &state->surface_noise,
            &palette_cycler.palette(),
            state->hue_rotation_lut.data(),
            state->hue_noise_lut.data(),
            params,
            palette_mapping,
            palette_oscillation_phase};
  }

  static void next_palette(void *context, uint32_t sequence,
                           GenerativePalette &out) {
    CurlLattice &effect = *static_cast<CurlLattice *>(context);
    if (sequence > 0)
      effect.palette_hue += HUE_STEP;
    out = GenerativePalette{PaletteRecipes::profile(
        PaletteDomain::STRAIGHT, PaletteHarmony::TRIADIC, AxisCurve::ASCENDING,
        PaletteRecipes::hue_turns(effect.palette_hue),
        effect.params.color.palette_chroma)};
  }

  static constexpr size_t FOOTPRINT_BYTES =
      gamut_lut_bytes(GAMUT_ANGLE_STEPS, GAMUT_L_STEPS) +
      PaletteCycler::generated_arena_bytes() +
      PARAM_CAPACITY * sizeof(ParamDef) + sizeof(State) + alignof(State);
  static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
                "CurlLattice persistent footprint exceeds the default "
                "partition");

  State *state = nullptr;
  Params params = PRESETS[0].params;
  Pullback::Color::PaletteMappingWeights palette_mapping =
      Pullback::Color::PaletteMappingWeights::single(
          params.color.palette_mapping);
  Transition transition;
  uint16_t preset_dwell_remaining = PRESET_DWELL_FRAMES;

  Timeline timeline;
  Orientation<> projection_walk;
  Orientation<> outer_walk;
  Quaternion projection_walk_previous;
  Quaternion outer_walk_previous;
  Quaternion projection_wander;
  Quaternion outer_wander;
  Quaternion projection_conjugate;
  Quaternion outer_conjugate;
  Quaternion base_orientation =
      make_rotation(Vector(0, 0, -1), Vector(0, -1, 0));

  float projection_spin = 0.0f;
  float surface_noise_phase = 0.0f;
  float hue_noise_phase = 0.0f;
  float palette_oscillation_phase = 0.0f;
  float palette_chroma = -1.0f;
  uint32_t palette_hue = 0;
  PaletteCycler palette_cycler;
};

REGISTER_EFFECT(CurlLattice)
