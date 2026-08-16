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
      "54468343cd539059ed564bb192969114e01b54bd17fec194653b300a6f674f92";
  static constexpr std::string_view PRESET_BANK_DIGEST =
      "3f15f1ba74f7d58d9933d05daa3ddcb0dcef93f5d85d53afdbafbcc7497bef51";
  static constexpr std::array<std::string_view, 2> PRESET_IDS{"open-curl",
                                                              "dense-curl"};
  static constexpr uint32_t PARAMETER_SCHEMA_VERSION = 1;
  static constexpr uint16_t TRANSITION_DURATION = 480;

  struct ParameterSnapshot {
    uint32_t schema_version;
    float surface_noise_scale;
  };

  HS_COLD_MEMBER CurlLattice() : Effect(W, H, {.strobe = true}) {}

  HS_COLD_MEMBER void init() override {
    configure_presets(PRESETS.size());
    state = persistent_arena.make<State>();

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

    register_animated_param("Surface Noise Scale", &params.surface_noise_scale,
                            SURFACE_SCALE_MIN, SURFACE_SCALE_MAX);
  }

  HS_FLASH_MEMBER void draw_frame() override {
    Canvas canvas(*this);
    timeline.step(canvas);
    begin_automatic_transition();
    prepare_transition_value();
    update_spatial_frames();
    palette_cycler.step();

    const FrameState frame = prepare_frame();
    Scan::Shader::draw<W, H, 1>(canvas, [&frame](const Vector &view) {
      return RenderPipeline::shade(view, frame);
    });
    finish_transition_evaluation();
  }

  ParameterSnapshot serialize_parameters() const {
    return {PARAMETER_SCHEMA_VERSION, params.surface_noise_scale};
  }

  bool restore_parameters(const ParameterSnapshot &snapshot) {
    if (snapshot.schema_version != PARAMETER_SCHEMA_VERSION ||
        !std::isfinite(snapshot.surface_noise_scale) ||
        snapshot.surface_noise_scale < SURFACE_SCALE_MIN ||
        snapshot.surface_noise_scale > SURFACE_SCALE_MAX)
      return false;
    transition.active = false;
    params.surface_noise_scale = snapshot.surface_noise_scale;
    return true;
  }

private:
  friend struct ::hs_test::curl_lattice_tests::CurlLatticeWhiteBox;

  static constexpr int GAMUT_ANGLE_STEPS = GAMUT_LUT_ANGLE_STEPS;
  static constexpr int GAMUT_L_STEPS = GAMUT_LUT_L_STEPS;
  static constexpr int32_t SURFACE_NOISE_SEED = 1337;
  static constexpr int32_t COLOR_NOISE_SEED = 6047;
  static constexpr float SURFACE_SCALE_MIN = 1.0f / 64.0f;
  static constexpr float SURFACE_SCALE_MAX = 64.0f;
  static constexpr float SURFACE_NOISE_STRENGTH = 0.0759999976f;
  static constexpr float HUE_SHIFT_AMOUNT = 0.268000007f;
  static constexpr float HUE_NOISE_SCALE = 2.0f;
  static constexpr float MAPPING_PHASE = -0.165999994f;
  static constexpr float POLE_FADE = 20.0f;
  static constexpr uint32_t HUE_STEP = 159;
  static constexpr uint16_t PRESET_DWELL_FRAMES = 600;

  struct Params {
    float surface_noise_scale;
  };

  struct Preset {
    std::string_view id;
    std::string_view display_name;
    Params params;
  };

  static constexpr std::array<Preset, 2> PRESETS{{
      {PRESET_IDS[0], "Open Curl", {1.78815627f}},
      {PRESET_IDS[1], "Dense Curl", {3.29720306f}},
  }};

  struct Transition {
    Params from{};
    Params to{};
    uint16_t evaluation = 0;
    uint16_t duration = TRANSITION_DURATION;
    bool active = false;
  };

  struct LatticeParams {
    float lattice_cell_scale = 0.710265636f;
    float lattice_shape_blend = 1.0f;
    float lattice_softness = 0.455532223f;
    float lattice_radius = 0.290762514f;
  };

  struct FrameState {
    Quaternion projection_conjugate;
    Quaternion outer_conjugate;
    Vector surface_loop_offset;
    const FastNoiseLite *surface_noise;
    const BakedPalette *palette;
    const Pixel *hue_rotation_lut;
    const int8_t *hue_noise_lut;
    LatticeParams lattice;
    float surface_noise_scale;
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
    static float central_meridian(const FrameState &) { return 0.0f; }
    static float pole_fade(const FrameState &) { return POLE_FADE; }
  };

  struct SurfaceStateProvider {
    using Binding = CurlLattice::Binding;
    using FrameState = typename Binding::FrameState;

    static const FastNoiseLite &noise(const FrameState &frame) {
      return *frame.surface_noise;
    }
    static float scale(const FrameState &frame) {
      return frame.surface_noise_scale;
    }
    static const Vector &loop_offset(const FrameState &frame) {
      return frame.surface_loop_offset;
    }
    static float strength(const FrameState &) { return SURFACE_NOISE_STRENGTH; }
    static bool path_length_required(const FrameState &) { return false; }
  };

  struct SourceStateProvider {
    using Binding = CurlLattice::Binding;
    using FrameState = typename Binding::FrameState;

    static const LatticeParams &params(const FrameState &frame) {
      return frame.lattice;
    }
  };

  struct ColorStateProvider {
    using Binding = CurlLattice::Binding;
    using FrameState = typename Binding::FrameState;

    static Pullback::Color::PaletteMapping mapping(const FrameState &) {
      return Pullback::Color::PaletteMapping::CUP;
    }
    static float mapping_frequency(const FrameState &) { return 1.0f; }
    static float mapping_phase(const FrameState &) { return MAPPING_PHASE; }
    static float oscillation_depth(const FrameState &) { return 0.0f; }
    static float oscillation_phase(const FrameState &) { return 0.0f; }
    static Color4 palette(const FrameState &frame, float value) {
      return frame.palette->get(value);
    }
    static Pullback::Color::HueMode hue_mode(const FrameState &) {
      return Pullback::Color::HueMode::NOISE;
    }
    static float hue_shift_amount(const FrameState &) {
      return HUE_SHIFT_AMOUNT;
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
    static float brightness_depth(const FrameState &) { return 1.0f; }
    static float opacity_low(const FrameState &) { return 1.0f; }
    static float opacity_high(const FrameState &) { return 1.0f; }
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
  };

  HS_COLD_MEMBER bool apply_preset(const PresetChange &change) override {
    const Params target = PRESETS[change.to].params;
    if (change.origin != PresetChangeOrigin::AUTOMATIC) {
      transition.active = false;
      params = target;
      preset_dwell_remaining = PRESET_DWELL_FRAMES;
      return true;
    }
    transition = {params, target, 0, TRANSITION_DURATION, true};
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
        state->color_noise, HUE_NOISE_SCALE, 0.0f);
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
    if (preset_dwell_remaining > 0) {
      --preset_dwell_remaining;
      return;
    }
    HS_CHECK(advancePreset());
  }

  HS_COLD_MEMBER void prepare_transition_value() {
    if (!transition.active)
      return;
    const FixedPipeline::EdgeProgress progress =
        FixedPipeline::edge_progress(transition.evaluation, transition.duration,
                                     FixedPipeline::Easing::EASE_IN_OUT_SIN);
    params.surface_noise_scale = FixedPipeline::linear(
        transition.from.surface_noise_scale, transition.to.surface_noise_scale,
        progress.eased);
  }

  HS_COLD_MEMBER void finish_transition_evaluation() {
    if (!transition.active || anims_paused)
      return;
    if (transition.evaluation == transition.duration) {
      params = transition.to;
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
        (projection_delta.normalized() * projection_wander).normalized();

    const Quaternion outer = outer_walk.get();
    const Quaternion outer_delta = outer * outer_walk_previous.conjugate();
    outer_walk_previous = outer;
    outer_wander = (outer_delta.normalized() * outer_wander).normalized();

    projection_conjugate = (base_orientation * projection_wander).conjugate();
    outer_conjugate = outer_wander.conjugate();
  }

  HS_COLD_MEMBER FrameState prepare_frame() const {
    prepare_hue_rotation_lut();
    return {projection_conjugate,
            outer_conjugate,
            Vector(NOISE_LOOP_RADIUS, 0.0f, 0.0f),
            &state->surface_noise,
            &palette_cycler.palette(),
            state->hue_rotation_lut.data(),
            state->hue_noise_lut.data(),
            {},
            params.surface_noise_scale};
  }

  static void next_palette(void *context, uint32_t sequence,
                           GenerativePalette &out) {
    CurlLattice &effect = *static_cast<CurlLattice *>(context);
    if (sequence > 0)
      effect.palette_hue += HUE_STEP;
    out = GenerativePalette{PaletteRecipes::profile(
        PaletteDomain::STRAIGHT, PaletteHarmony::TRIADIC, AxisCurve::ASCENDING,
        PaletteRecipes::hue_turns(effect.palette_hue), 1.0f)};
  }

  static constexpr size_t FOOTPRINT_BYTES =
      gamut_lut_bytes(GAMUT_ANGLE_STEPS, GAMUT_L_STEPS) +
      PaletteCycler::generated_arena_bytes() + sizeof(State) + alignof(State);
  static_assert(FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
                "CurlLattice persistent footprint exceeds the default "
                "partition");

  State *state = nullptr;
  Params params = PRESETS[0].params;
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

  uint32_t palette_hue = 0;
  PaletteCycler palette_cycler;
};

REGISTER_EFFECT(CurlLattice)
