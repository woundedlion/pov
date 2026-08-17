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
namespace facet_grid_tests {
struct FacetGridWhiteBox;
} // namespace facet_grid_tests
} // namespace hs_test

/**
 * @brief Mirrored grids folded through a dodecahedral stereographic lens.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H> class FacetGrid : public Effect {
public:
  static constexpr std::string_view EFFECT_ID = "facet-grid";
  static constexpr std::string_view DESCRIPTOR_DIGEST =
      "966c4aad90ff062c5c6ab25b2308fec6acbe406cbe5885096d71bde83284c804";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "75a78672e4b2ad32eee8bd54cd7d7b9408db1d7bff79e084a1b58ccb33b5232d";
  static constexpr std::array<std::string_view, 4> PRESET_IDS{
      "coupled-grid", "direct-grid", "double-map", "stretched-grid"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 2;
  static constexpr uint16_t TRANSITION_DURATION = 480;

  struct SourceParams {
    float pattern_freq = 2.82629991f;
    float speed = 0.0f;
    float complexity = 0.513f;
    float pattern_mix = 0.0f;
    float secondary_rate = 0.8f;
    float angle_speed = 0.0269999988f;
  };

  struct ProjectionParams {
    float pole_fade = 3.432f;
    float spin_speed = 0.0f;
    float wander = 1.0f;
    float camera_wander = 1.0f;
  };

  struct MirrorParams {
    float speed = 0.00013f;
    float rotation = 0.0f;
    float cell_x = 1.0f;
    float cell_y = 0.997703135f;
    float offset_x = 0.0f;
    float offset_y = 0.0f;
  };

  struct ColorParams {
    float palette_chroma = 1.0f;
    float mapping_frequency = 1.0f;
    float mapping_phase = 0.0f;
    float phase_oscillation_depth = 0.0f;
    float phase_oscillation_speed = 0.0f;
    float opacity_low = 1.0f;
    float opacity_high = 1.0f;
    float hue_shift_amount = 0.366f;
    float hue_noise_scale = 1.47215629f;
    float hue_noise_speed = 0.0f;
    Pullback::Color::PaletteMapping palette_mapping =
        Pullback::Color::PaletteMapping::CUP;
  };

  struct Params {
    SourceParams source;
    ProjectionParams projection;
    MirrorParams mirror;
    ColorParams color;
  };

  struct ParameterSnapshot {
    uint32_t schema_version;
    Params params;
  };

  HS_COLD_MEMBER FacetGrid() : Effect(W, H, {.strobe = true}) {}

  HS_COLD_MEMBER void init() override {
    configure_presets(PRESETS.size());
    state = persistent_arena.make<State>();
    use_parameter_storage(persistent_arena.allocate_n<ParamDef>(PARAM_CAPACITY),
                          PARAM_CAPACITY);

    configure_noise(state->color_noise, COLOR_NOISE_SEED);
    prepare_hue_noise_lut();

    init_gamut_lut(persistent_arena, GAMUT_ANGLE_STEPS, GAMUT_L_STEPS);
    palette_cycler.init_generated(persistent_arena, next_palette, this, 0, 600,
                                  ease_in_out_sin);

    timeline.add(0, Animation::RandomWalk<W>(projection_walk, UP,
                                             state->projection_walk_noise));
    timeline.add(
        0, Animation::RandomWalk<W>(outer_walk, UP, state->outer_walk_noise));

    register_animated_param("Pattern Freq", &params.source.pattern_freq,
                            PATTERN_FREQ_MIN, PATTERN_FREQ_MAX);
    register_animated_param("Speed", &params.source.speed, SPEED_MIN,
                            SPEED_MAX);
    register_animated_param("Source Angle Speed", &params.source.angle_speed,
                            ANGLE_SPEED_MIN, ANGLE_SPEED_MAX);
    register_animated_param("Complexity", &params.source.complexity, 0.0f,
                            3.0f);
    register_animated_param("Pattern Mix", &params.source.pattern_mix, 0.0f,
                            1.0f);
    register_animated_param("Drift", &params.source.secondary_rate, 0.0f,
                            DRIFT_MAX);
    register_animated_param("Pole Fade", &params.projection.pole_fade,
                            POLE_FADE_MIN, POLE_FADE_MAX);
    register_animated_param("Projection Spin Speed",
                            &params.projection.spin_speed, 0.0f,
                            PROJECTION_SPIN_MAX);
    register_animated_param("Projection Wander", &params.projection.wander,
                            0.0f, 1.0f);
    register_animated_param("Camera Wander", &params.projection.camera_wander,
                            0.0f, 1.0f);
    register_animated_param("Planar Warp 2 Speed", &params.mirror.speed,
                            -MIRROR_SPEED_MAX, MIRROR_SPEED_MAX);
    register_animated_param("Planar Warp 2 Rotation", &params.mirror.rotation,
                            0.0f, TWO_PI_F);
    register_animated_param("Planar Warp 2 Cell X", &params.mirror.cell_x,
                            CELL_MIN, CELL_MAX);
    register_animated_param("Planar Warp 2 Cell Y", &params.mirror.cell_y,
                            CELL_MIN, CELL_MAX);
    register_animated_param("Planar Warp 2 Offset X", &params.mirror.offset_x,
                            -8.0f, 8.0f);
    register_animated_param("Planar Warp 2 Offset Y", &params.mirror.offset_y,
                            -8.0f, 8.0f);
    register_animated_param("Palette Chroma", &params.color.palette_chroma,
                            0.0f, 1.0f);
    register_animated_param("Palette Mapping", &params.color.palette_mapping,
                            PALETTE_MAPPING_OPTIONS,
                            PALETTE_MAPPING_EXPORT_OPTIONS,
                            std::size(PALETTE_MAPPING_OPTIONS));
    register_animated_param("Mapping Frequency",
                            &params.color.mapping_frequency, 1.0f, 32.0f);
    register_animated_param("Mapping Phase", &params.color.mapping_phase, -1.0f,
                            1.0f);
    register_animated_param("Phase Oscillation Depth",
                            &params.color.phase_oscillation_depth, 0.0f, 1.0f);
    register_animated_param("Phase Oscillation Speed",
                            &params.color.phase_oscillation_speed, -0.01f,
                            0.01f);
    register_animated_param("Value Opacity Low", &params.color.opacity_low,
                            0.0f, 1.0f);
    register_animated_param("Value Opacity High", &params.color.opacity_high,
                            0.0f, 1.0f);
    register_animated_param("Hue Shift Amount", &params.color.hue_shift_amount,
                            -1.0f, 1.0f);
    register_animated_param("Hue Noise Scale", &params.color.hue_noise_scale,
                            HUE_NOISE_SCALE_MIN, HUE_NOISE_SCALE_MAX);
    register_animated_param("Hue Noise Speed", &params.color.hue_noise_speed,
                            -HUE_NOISE_SPEED_MAX, HUE_NOISE_SPEED_MAX);
  }

  HS_FLASH_MEMBER void draw_frame() override {
    Canvas canvas(*this);
    {
      HS_PROFILE(fg_timeline_step);
      timeline.step(canvas);
    }
    {
      HS_PROFILE(fg_advance);
      begin_automatic_transition();
      prepare_transition_value();
      advance_runtime();
      update_spatial_frames();
      update_palette_chroma();
      palette_cycler.step();
    }
    const FrameState frame = prepare_frame();
    {
      HS_PROFILE(fg_shader_draw);
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
    preset_dwell_remaining = PRESET_DWELL_FRAMES;
    return true;
  }

private:
  friend struct ::hs_test::facet_grid_tests::FacetGridWhiteBox;

  static constexpr int GAMUT_ANGLE_STEPS = GAMUT_LUT_ANGLE_STEPS;
  static constexpr int GAMUT_L_STEPS = GAMUT_LUT_L_STEPS;
  static constexpr int32_t COLOR_NOISE_SEED = 6047;
  static constexpr float PATTERN_FREQ_MIN = 0.1f;
  static constexpr float PATTERN_FREQ_MAX = 20.0f;
  static constexpr float SPEED_MIN = 0.0f;
  static constexpr float SPEED_MAX = 0.5f;
  static constexpr float ANGLE_SPEED_MIN = 0.0f;
  static constexpr float ANGLE_SPEED_MAX = 0.03f;
  static constexpr float DRIFT_MAX = 1.25f;
  static constexpr float POLE_FADE_MIN = 1.0f;
  static constexpr float POLE_FADE_MAX = 20.0f;
  static constexpr float PROJECTION_SPIN_MAX = 0.04f;
  static constexpr float MIRROR_SPEED_MAX = 0.005f;
  static constexpr float CELL_MIN = 1.0f / 64.0f;
  static constexpr float CELL_MAX = 8.0f;
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

  static constexpr Params direct_grid_params() {
    Params value;
    value.source.complexity = 3.0f;
    value.source.pattern_mix = 1.0f;
    return value;
  }

  static constexpr Params double_map_params() {
    Params value = direct_grid_params();
    value.source.pattern_freq = 3.9407f;
    value.projection.wander = 0.165f;
    value.color.mapping_frequency = 2.0f;
    return value;
  }

  static constexpr Params stretched_grid_params() {
    Params value = direct_grid_params();
    value.source.pattern_freq = 2.9059f;
    value.projection.wander = 0.165f;
    value.mirror.speed = 0.0027299998f;
    value.mirror.rotation = 3.455752f;
    value.mirror.cell_x = 0.22321875f;
    value.mirror.cell_y = 5.085703f;
    value.color.mapping_frequency = 1.558f;
    return value;
  }

  static constexpr std::array<Preset, 4> PRESETS{{
      {PRESET_IDS[0], "Coupled Grid", {}},
      {PRESET_IDS[1], "Direct Grid", direct_grid_params()},
      {PRESET_IDS[2], "Double Map", double_map_params()},
      {PRESET_IDS[3], "Stretched Grid", stretched_grid_params()},
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

  struct PreparedSource {
    float primary;
    float secondary;
    float angle_cos;
    float angle_sin;
  };

  struct PreparedMirror {
    struct Transform {
      struct Mirror {
        float offset_x;
        float offset_y;
      } mirror;
    } transform;
    float rotation_cos;
    float rotation_sin;
  };

  struct FrameState {
    Quaternion projection_conjugate;
    Quaternion outer_conjugate;
    PreparedSource prepared_source;
    PreparedMirror prepared_mirror;
    const BakedPalette *palette;
    const Pixel *hue_rotation_lut;
    const int8_t *hue_noise_lut;
    Params params;
    Pullback::Color::PaletteMappingWeights palette_mapping;
    float palette_oscillation_phase;
  };

  struct Binding {
    using FrameState = FacetGrid::FrameState;
    using Instrumentation = Pullback::NoInstrumentation;
  };

  struct OuterCameraStateProvider {
    using Binding = FacetGrid::Binding;
    using FrameState = typename Binding::FrameState;
    static const Quaternion &conjugate(const FrameState &frame) {
      return frame.outer_conjugate;
    }
  };

  struct ProjectionStateProvider {
    using Binding = FacetGrid::Binding;
    using FrameState = typename Binding::FrameState;
    static const Quaternion &conjugate(const FrameState &frame) {
      return frame.projection_conjugate;
    }
    static float pole_fade(const FrameState &frame) {
      return frame.params.projection.pole_fade;
    }
  };

  struct MirrorStateProvider {
    using Binding = FacetGrid::Binding;
    using FrameState = typename Binding::FrameState;
    static const MirrorParams &params(const FrameState &frame) {
      return frame.params.mirror;
    }
    static const PreparedMirror &prepared(const FrameState &frame) {
      return frame.prepared_mirror;
    }
    static bool path_length_required(const FrameState &) { return false; }
  };

  struct SourceStateProvider {
    using Binding = FacetGrid::Binding;
    using FrameState = typename Binding::FrameState;
    static const SourceParams &params(const FrameState &frame) {
      return frame.params.source;
    }
    static const PreparedSource &prepared(const FrameState &frame) {
      return frame.prepared_source;
    }
  };

  struct ColorStateProvider {
    using Binding = FacetGrid::Binding;
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
      return {frame.hue_rotation_lut,
              frame.params.color.hue_shift_amount != 0.0f};
    }
    static Pullback::Color::HueNoiseLutView hue_noise(const FrameState &frame) {
      return {frame.hue_noise_lut, frame.params.color.hue_shift_amount != 0.0f};
    }
    static Pullback::Color::BrightnessEnvelope
    brightness_envelope(const FrameState &) {
      return Pullback::Color::BrightnessEnvelope::NONE;
    }
    static float brightness_depth(const FrameState &) { return 1.0f; }
    static float opacity_low(const FrameState &frame) {
      return frame.params.color.opacity_low;
    }
    static float opacity_high(const FrameState &frame) {
      return frame.params.color.opacity_high;
    }
  };

  using OuterCameraStage =
      Pullback::Stage::OuterCamera<Binding, OuterCameraStateProvider>;
  using SurfaceStage = Pullback::Stage::SurfaceProject<
      Binding, Pullback::Surface::Identity,
      Pullback::Lens::DodecahedralKaleidoscope, Pullback::Surface::Identity,
      Pullback::Projection::Stereographic<ProjectionStateProvider>>;
  using PlanarWarpStage = Pullback::Stage::PlanarWarp<
      Binding, Pullback::Warp::Identity,
      Pullback::Warp::MirrorTile<MirrorStateProvider>>;
  using SourceStage =
      Pullback::Stage::Source<Binding,
                              Pullback::Source::Grid<SourceStateProvider>>;
  using MaterialStage =
      Pullback::Stage::Material<Binding, Pullback::Weight::Projection,
                                Pullback::Transfer::Linear,
                                Pullback::Coverage::ProjectionSquared>;
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
    static_assert(sizeof(Params) == 27 * sizeof(float));
    return in_range(value.source.pattern_freq, PATTERN_FREQ_MIN,
                    PATTERN_FREQ_MAX) &&
           in_range(value.source.speed, SPEED_MIN, SPEED_MAX) &&
           in_range(value.source.complexity, 0.0f, 3.0f) &&
           in_range(value.source.pattern_mix, 0.0f, 1.0f) &&
           in_range(value.source.secondary_rate, 0.0f, DRIFT_MAX) &&
           in_range(value.source.angle_speed, ANGLE_SPEED_MIN,
                    ANGLE_SPEED_MAX) &&
           in_range(value.projection.pole_fade, POLE_FADE_MIN, POLE_FADE_MAX) &&
           in_range(value.projection.spin_speed, 0.0f, PROJECTION_SPIN_MAX) &&
           in_range(value.projection.wander, 0.0f, 1.0f) &&
           in_range(value.projection.camera_wander, 0.0f, 1.0f) &&
           in_range(value.mirror.speed, -MIRROR_SPEED_MAX, MIRROR_SPEED_MAX) &&
           in_range(value.mirror.rotation, 0.0f, TWO_PI_F) &&
           in_range(value.mirror.cell_x, CELL_MIN, CELL_MAX) &&
           in_range(value.mirror.cell_y, CELL_MIN, CELL_MAX) &&
           in_range(value.mirror.offset_x, -8.0f, 8.0f) &&
           in_range(value.mirror.offset_y, -8.0f, 8.0f) &&
           in_range(value.color.palette_chroma, 0.0f, 1.0f) &&
           static_cast<uint8_t>(value.color.palette_mapping) <=
               static_cast<uint8_t>(Pullback::Color::PaletteMapping::REVERSE) &&
           in_range(value.color.mapping_frequency, 1.0f, 32.0f) &&
           in_range(value.color.mapping_phase, -1.0f, 1.0f) &&
           in_range(value.color.phase_oscillation_depth, 0.0f, 1.0f) &&
           in_range(value.color.phase_oscillation_speed, -0.01f, 0.01f) &&
           in_range(value.color.opacity_low, 0.0f, 1.0f) &&
           in_range(value.color.opacity_high, 0.0f, 1.0f) &&
           in_range(value.color.hue_shift_amount, -1.0f, 1.0f) &&
           in_range(value.color.hue_noise_scale, HUE_NOISE_SCALE_MIN,
                    HUE_NOISE_SCALE_MAX) &&
           in_range(value.color.hue_noise_speed, -HUE_NOISE_SPEED_MAX,
                    HUE_NOISE_SPEED_MAX);
  }

  static Params interpolate_params(const Params &from, const Params &to,
                                   float progress) {
    Params value;
#define HS_FACET_LINEAR(path)                                                  \
  value.path = FixedPipeline::linear(from.path, to.path, progress)
    HS_FACET_LINEAR(source.pattern_freq);
    HS_FACET_LINEAR(source.speed);
    HS_FACET_LINEAR(source.complexity);
    HS_FACET_LINEAR(source.pattern_mix);
    HS_FACET_LINEAR(source.secondary_rate);
    HS_FACET_LINEAR(source.angle_speed);
    HS_FACET_LINEAR(projection.pole_fade);
    HS_FACET_LINEAR(projection.spin_speed);
    HS_FACET_LINEAR(projection.wander);
    HS_FACET_LINEAR(projection.camera_wander);
    HS_FACET_LINEAR(mirror.speed);
    value.mirror.rotation = FixedPipeline::shortest_periodic(
        from.mirror.rotation, to.mirror.rotation, progress, TWO_PI_F);
    HS_FACET_LINEAR(mirror.cell_x);
    HS_FACET_LINEAR(mirror.cell_y);
    HS_FACET_LINEAR(mirror.offset_x);
    HS_FACET_LINEAR(mirror.offset_y);
    HS_FACET_LINEAR(color.palette_chroma);
    value.color.palette_mapping =
        progress < 1.0f ? from.color.palette_mapping : to.color.palette_mapping;
    HS_FACET_LINEAR(color.mapping_frequency);
    HS_FACET_LINEAR(color.mapping_phase);
    HS_FACET_LINEAR(color.phase_oscillation_depth);
    HS_FACET_LINEAR(color.phase_oscillation_speed);
    HS_FACET_LINEAR(color.opacity_low);
    HS_FACET_LINEAR(color.opacity_high);
    HS_FACET_LINEAR(color.hue_shift_amount);
    HS_FACET_LINEAR(color.hue_noise_scale);
    HS_FACET_LINEAR(color.hue_noise_speed);
#undef HS_FACET_LINEAR
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
    source_primary = fmodf(source_primary + params.source.speed, TWO_PI_F);
    source_secondary = fmodf(
        source_secondary + params.source.speed * params.source.secondary_rate,
        TWO_PI_F);
    source_angle = fmodf(source_angle + params.source.angle_speed, TWO_PI_F);
    projection_spin =
        fmodf(projection_spin + params.projection.spin_speed, TWO_PI_F);
    mirror_phase = wrap_t(mirror_phase + params.mirror.speed);
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
    HS_PROFILE(fg_prepare_frame);
    if (state->hue_noise_lut_scale != params.color.hue_noise_scale ||
        state->hue_noise_lut_phase != hue_noise_phase)
      prepare_hue_noise_lut();
    prepare_hue_rotation_lut();
    const PreparedMirror mirror{
        {
            wrap_t(params.mirror.offset_x / params.mirror.cell_x +
                   mirror_phase) *
                params.mirror.cell_x,
            wrap_t(params.mirror.offset_y / params.mirror.cell_y) *
                params.mirror.cell_y,
        },
        cosf(params.mirror.rotation),
        sinf(params.mirror.rotation)};
    return {projection_conjugate,
            outer_conjugate,
            {source_primary, source_secondary, fast_cosf(source_angle),
             fast_sinf(source_angle)},
            mirror,
            &palette_cycler.palette(),
            state->hue_rotation_lut.data(),
            state->hue_noise_lut.data(),
            params,
            palette_mapping,
            palette_oscillation_phase};
  }

  static void next_palette(void *context, uint32_t sequence,
                           GenerativePalette &out) {
    FacetGrid &effect = *static_cast<FacetGrid *>(context);
    if (sequence > 0)
      effect.palette_hue += HUE_STEP;
    out = GenerativePalette{PaletteRecipes::profile(
        PaletteDomain::STRAIGHT, PaletteHarmony::ANALOGOUS,
        AxisCurve::ASCENDING, PaletteRecipes::hue_turns(effect.palette_hue),
        effect.params.color.palette_chroma)};
  }

  static constexpr size_t FOOTPRINT_BYTES =
      gamut_lut_bytes(GAMUT_ANGLE_STEPS, GAMUT_L_STEPS) +
      PaletteCycler::generated_arena_bytes() +
      PARAM_CAPACITY * sizeof(ParamDef) + sizeof(State) + alignof(State);
  static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
                "FacetGrid persistent footprint exceeds the default "
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

  float source_primary = 0.0f;
  float source_secondary = 0.0f;
  float source_angle = 0.0f;
  float projection_spin = 0.0f;
  float mirror_phase = 0.0f;
  float hue_noise_phase = 0.0f;
  float palette_oscillation_phase = 0.0f;
  float palette_chroma = -1.0f;
  uint32_t palette_hue = 0;
  PaletteCycler palette_cycler;
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(FacetGrid)
