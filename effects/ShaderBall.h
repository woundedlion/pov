/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <cstdarg>
#include <cstdio>

/**
 * @file ShaderBall.h
 * @brief Typed pullback sphere shader with composable projection and material stages.
 */

#include "core/color/effect_palette_recipes.h"
#include "core/engine/engine.h"
#include "core/math/projections.h"

namespace hs_test {
namespace shaderball_tests {
struct ShaderBallWhiteBox;
} // namespace shaderball_tests
} // namespace hs_test

/**
 * @brief Slot-based sphere shader with an immutable per-frame pullback state.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H> class ShaderBall : public Effect {
public:
  static constexpr int GAMUT_ANGLE_STEPS = GAMUT_LUT_ANGLE_STEPS;
  static constexpr int GAMUT_L_STEPS = GAMUT_LUT_L_STEPS;

  HS_COLD_MEMBER ShaderBall() : Effect(W, H, {.strobe = true}) {}

  /** @brief Initializes slots, clocks, palette resources, and choreography. */
  HS_COLD_MEMBER void init() override {
    configure_presets(PRESETS.size());
#if HS_ENABLE_PARAM_GUI_BRIDGE
    set_parameter_updated_hook(&ShaderBall::dispatch_parameter_updated);
#endif
    state = persistent_arena.make<StateBundle>();
    use_parameter_storage(persistent_arena.allocate_n<ParamDef>(PARAM_CAPACITY),
                          PARAM_CAPACITY);
    requested_config = PRESETS[0];
    published_config = PRESETS[0];
    accepted_config = PRESETS[0];
    prepare_resource_union(PRESETS[0], PRESETS[0]);

    rebind_parameters();

    timeline.add(0, Animation::RandomWalk<W>(projection_walk, UP,
                                             state->projection_walk_noise));
    timeline.add(
        0, Animation::RandomWalk<W>(outer_walk, UP, state->outer_walk_noise));

    liquid_palette_cycler.init_generated(persistent_arena, next_liquid_palette,
                                         &liquid_rotation, PALETTE_DWELL_FRAMES,
                                         PALETTE_FADE_FRAMES, ease_in_out_sin);
    init_gamut_lut(persistent_arena, GAMUT_ANGLE_STEPS, GAMUT_L_STEPS);
    generated_palette_cycler.init_generated(
        persistent_arena, next_triadic_palette, &palette_hue,
        PALETTE_DWELL_FRAMES, PALETTE_FADE_FRAMES, ease_in_out_sin);

    enter_preset();
  }

  /** @brief Advances mutable state, snapshots it, and renders one frame. */
  void draw_frame() override {
    Canvas canvas(*this);
    {
      HS_PROFILE(sb_timeline_step);
      timeline.step(canvas);
    }
    advance_preset_choreography();

    apply_requested_config();
    prepare_param_morph();
    const WalkDeltas walk_deltas = sample_walk_deltas();
    if (state->transition.active) {
      advance_runtime(state->transition.from_runtime,
                      state->transition.from_config, walk_deltas);
      advance_runtime(state->transition.to_runtime, state->transition.to_config,
                      walk_deltas);
    } else {
      advance_runtime(runtime, {active_slots, blend.params}, walk_deltas);
    }
    liquid_palette_cycler.step();
    generated_palette_cycler.step();
#if HS_ENABLE_TEST_HOOKS
    ++liquid_palette_step_count;
    ++generated_palette_step_count;
#endif

    if (state->transition.active) {
      draw_through_clear_transition(canvas);
    } else {
      const FrameState frame = prepare_frame();
      FrameShader shader{&frame, 1.0f};
      HS_PROFILE(sb_shader_draw);
      Scan::Shader::draw<W, H, 1>(canvas, shader);
    }
    finish_transitions();
    publish_live_config();
  }

#if HS_ENABLE_EFFECT_CONTROL_API
  void profile_select_preset(size_t index) {
    HS_CHECK(index < PRESETS.size(),
             "ShaderBall profile preset index out of range");
    HS_CHECK(selectPreset(index), "ShaderBall profile preset selection failed");
    hs::log("Profile preset: %u/%u", static_cast<unsigned>(index),
            static_cast<unsigned>(PRESETS.size()));
  }
#endif

private:
  friend struct ::hs_test::shaderball_tests::ShaderBallWhiteBox;

  HS_COLD_MEMBER bool apply_preset(const PresetChange &change) override {
    const size_t index = change.to;
    if (change.origin == PresetChangeOrigin::AUTOMATIC) {
      const Choreo choreo = preset_choreo(change.from);
      const Preset &to = PRESETS[index];
      if (!try_apply_config(to, choreo.blend_frames, choreo.staggered, true))
        return false;
      requested_config = to;
      published_config = to;
      accepted_config = to;
#if HS_ENABLE_PARAM_GUI_BRIDGE
      pending_edit_count = 0;
#endif
      rebind_parameters();
      return true;
    }

    state->param_morph.active = false;
    state->transition.active = false;
    active_slots = PRESETS[index].slots;
    blend.params = PRESETS[index].params;
#if HS_ENABLE_PARAM_GUI_BRIDGE
    display_config = PRESETS[index];
#endif
    requested_config = PRESETS[index];
    published_config = PRESETS[index];
    accepted_config = PRESETS[index];
#if HS_ENABLE_PARAM_GUI_BRIDGE
    pending_edit_count = 0;
#endif
    runtime = {};
    HS_CHECK(prepare_resource_union(PRESETS[index], PRESETS[index]),
             "ShaderBall preset resources exceed capacity");
    rebind_parameters();
    return true;
  }

  HS_COLD_MEMBER void preset_changed(const PresetChange &) override {
    if (!state->param_morph.active && !state->transition.active)
      enter_preset();
  }

  enum class Function : uint8_t {
    TWIN_WAVE,
    RINGS,
    SPIRAL,
    GRID,
    NOISE_CONTOUR,
    PRIMITIVE_LATTICE
  };
  enum class Projection : uint8_t {
    SINUSOIDAL,
    STEREOGRAPHIC,
    GNOMONIC,
    BONNE,
    PEIRCE_QUINCUNCIAL,
    AIROCEAN,
    EQUIRECTANGULAR
  };
  enum class PeirceLayout : uint8_t { DIAMOND, SQUARE, HORIZONTAL, VERTICAL };
  enum class AiroceanLayout : uint8_t { VERTICAL, HORIZONTAL };
  enum class BonneHemisphere : uint8_t { NORTH, SOUTH };
  enum class GnomonicHemispherePolicy : uint8_t {
    FOLDED,
    FRONT_HEMISPHERE,
    BACK_HEMISPHERE
  };
  enum class SurfaceLens : uint8_t {
    NONE,
    GLITCH,
    TWIST,
    KALEIDOSCOPE,
    MOBIUS,
    TANGENT_NOISE,
    KALEIDOSCOPE_TETRAHEDRAL,
    KALEIDOSCOPE_OCTAHEDRAL,
    KALEIDOSCOPE_DODECAHEDRAL,
    KALEIDOSCOPE_TRIANGULAR_PRISM,
    KALEIDOSCOPE_SQUARE_PRISM,
    KALEIDOSCOPE_PENTAGONAL_PRISM,
    KALEIDOSCOPE_HEXAGONAL_PRISM,
    KALEIDOSCOPE_OCTAGONAL_PRISM
  };
  enum class NoiseBasis : uint8_t { SIMPLEX, FBM3, RIDGED3 };
  enum class WarpEnvelope : uint8_t { FLAT, PROJECTION_WEIGHT, EDGE_FADE };
  enum class PolarMode : uint8_t { LINEAR, LOGARITHMIC };
  enum class CurlIntegrator : uint8_t { EULER_1, MIDPOINT_2, MIDPOINT_4 };
  enum class WarpStageKind : uint8_t {
    NONE,
    LEGACY_STEREO_NOISE,
    AFFINE_FRAME,
    WAVE_SHEAR,
    VORTEX,
    VECTOR_NOISE,
    CURL_FLOW,
    MIRROR_TILE,
    POLAR_CHART
  };
  struct WarpStageSpec {
    WarpStageKind kind;
    NoiseBasis basis = NoiseBasis::SIMPLEX;
    WarpEnvelope envelope = WarpEnvelope::FLAT;
    PolarMode polar_mode = PolarMode::LINEAR;
    CurlIntegrator curl_integrator = CurlIntegrator::EULER_1;
    uint8_t polar_harmonic = 1;
    int32_t seed = 1337;
    uint8_t resource_id = 0;

    HS_COLD_MEMBER bool operator==(const WarpStageSpec &) const = default;
  };
  struct WarpProgram {
    WarpStageSpec outer;
    WarpStageSpec inner;

    HS_COLD_MEMBER bool operator==(const WarpProgram &) const = default;
  };
  enum class ProjectionFramePolicy : uint8_t { IDENTITY, SPIN_WANDER };
  enum class SignalWeight : uint8_t { NONE, PROJECTION };
  enum class ValueTransfer : uint8_t {
    LINEAR,
    RIDGE,
    ISO_CONTOUR,
    SMOOTH_BANDS
  };
  enum class CoveragePolicy : uint8_t {
    OPAQUE,
    PROJECTION_WEIGHT_SQUARED,
    VALUE_CUTOUT,
    EDGE_FADE,
    PROJECTION_WEIGHT
  };
  enum class Colorizer : uint8_t { GENERATED_TRIADIC, LIQUID, DEFORMATION_INK };

  struct Slots {
    Function function;
    Projection projection;
    ProjectionFramePolicy projection_frame;
    SurfaceLens surface_lens;
    WarpProgram warp_program;
    SignalWeight signal_weight;
    ValueTransfer value_transfer;
    CoveragePolicy coverage;
    Colorizer colorizer;
    PeirceLayout peirce_layout = PeirceLayout::SQUARE;
    AiroceanLayout airocean_layout = AiroceanLayout::VERTICAL;
    BonneHemisphere bonne_hemisphere = BonneHemisphere::NORTH;
    GnomonicHemispherePolicy gnomonic_hemisphere =
        GnomonicHemispherePolicy::FOLDED;

    HS_COLD_MEMBER bool operator==(const Slots &) const = default;
  };

  struct SourceParams {
    float pattern_freq = 1.0f;
    float speed = 0.0f;
    float complexity = 0.0f;
    float pattern_mix = 0.0f;
    float secondary_rate = 0.0f;
    float angle_rate = 0.0f;
    float noise_scale = 1.0f;
    float noise_contrast = 0.0f;
    float noise_time_rate = 0.0f;
    float lattice_cell_scale = 1.0f;
    float lattice_shape_blend = 0.0f;
    float lattice_softness = 0.05f;
    float lattice_radius = 0.25f;
    NoiseBasis noise_basis = NoiseBasis::SIMPLEX;
    int32_t noise_seed = 2927;
    uint8_t noise_resource_id = 2;

    constexpr SourceParams() = default;

    constexpr SourceParams(float pattern_freq, float speed, float complexity,
                           float pattern_mix, float secondary_rate,
                           float angle_rate = 0.0f)
        : pattern_freq(pattern_freq), speed(speed), complexity(complexity),
          pattern_mix(pattern_mix), secondary_rate(secondary_rate),
          angle_rate(angle_rate) {}

    HS_COLD_MEMBER bool operator==(const SourceParams &) const = default;

    HS_COLD_MEMBER void lerp(const SourceParams &a, const SourceParams &b,
                             float t) {
      pattern_freq = hs::lerp(a.pattern_freq, b.pattern_freq, t);
      speed = hs::lerp(a.speed, b.speed, t);
      complexity = hs::lerp(a.complexity, b.complexity, t);
      pattern_mix = hs::lerp(a.pattern_mix, b.pattern_mix, t);
      secondary_rate = hs::lerp(a.secondary_rate, b.secondary_rate, t);
      angle_rate = hs::lerp(a.angle_rate, b.angle_rate, t);
      noise_scale = hs::lerp(a.noise_scale, b.noise_scale, t);
      noise_contrast = hs::lerp(a.noise_contrast, b.noise_contrast, t);
      noise_time_rate = hs::lerp(a.noise_time_rate, b.noise_time_rate, t);
      lattice_cell_scale =
          hs::lerp(a.lattice_cell_scale, b.lattice_cell_scale, t);
      lattice_shape_blend =
          hs::lerp(a.lattice_shape_blend, b.lattice_shape_blend, t);
      lattice_softness = hs::lerp(a.lattice_softness, b.lattice_softness, t);
      lattice_radius = hs::lerp(a.lattice_radius, b.lattice_radius, t);
      noise_basis = t < 1.0f ? a.noise_basis : b.noise_basis;
    }
  };

  struct WarpStageParams {
    float scale = 1.0f;
    float strength = 0.0f;
    float time_scale = 0.1f;
    float translation_x = 0.0f;
    float translation_y = 0.0f;
    float rotation = 0.0f;
    float scale_x = 1.0f;
    float scale_y = 1.0f;
    float shear = 0.0f;
    float frequency = 1.0f;
    float field_angle = 0.0f;
    float center_x = 0.0f;
    float center_y = 0.0f;
    float radius = 1.0f;
    float turns = 0.0f;
    float center_orbit_radius = 0.0f;
    float vector_angle = 0.0f;
    float cell_x = 1.0f;
    float cell_y = 1.0f;
    float offset_x = 0.0f;
    float offset_y = 0.0f;
    float radial_scale = 1.0f;
    float radial_phase = 0.0f;
    float angular_phase = 0.0f;
    float edge_width = 0.1f;

    constexpr WarpStageParams() = default;

    constexpr WarpStageParams(float scale, float strength, float time_scale)
        : scale(scale), strength(strength), time_scale(time_scale) {}

    HS_COLD_MEMBER bool operator==(const WarpStageParams &) const = default;

    HS_COLD_MEMBER void lerp(const WarpStageParams &a, const WarpStageParams &b,
                             float t) {
      scale = hs::lerp(a.scale, b.scale, t);
      strength = hs::lerp(a.strength, b.strength, t);
      time_scale = hs::lerp(a.time_scale, b.time_scale, t);
      translation_x = hs::lerp(a.translation_x, b.translation_x, t);
      translation_y = hs::lerp(a.translation_y, b.translation_y, t);
      rotation = lerp_angle(a.rotation, b.rotation, t);
      scale_x = expf(hs::lerp(logf(a.scale_x), logf(b.scale_x), t));
      scale_y = expf(hs::lerp(logf(a.scale_y), logf(b.scale_y), t));
      shear = hs::lerp(a.shear, b.shear, t);
      frequency = hs::lerp(a.frequency, b.frequency, t);
      field_angle = lerp_angle(a.field_angle, b.field_angle, t);
      center_x = hs::lerp(a.center_x, b.center_x, t);
      center_y = hs::lerp(a.center_y, b.center_y, t);
      radius = hs::lerp(a.radius, b.radius, t);
      turns = hs::lerp(a.turns, b.turns, t);
      center_orbit_radius =
          hs::lerp(a.center_orbit_radius, b.center_orbit_radius, t);
      vector_angle = lerp_angle(a.vector_angle, b.vector_angle, t);
      cell_x = hs::lerp(a.cell_x, b.cell_x, t);
      cell_y = hs::lerp(a.cell_y, b.cell_y, t);
      offset_x = hs::lerp(a.offset_x, b.offset_x, t);
      offset_y = hs::lerp(a.offset_y, b.offset_y, t);
      radial_scale = hs::lerp(a.radial_scale, b.radial_scale, t);
      radial_phase = lerp_angle(a.radial_phase, b.radial_phase, t);
      angular_phase = lerp_angle(a.angular_phase, b.angular_phase, t);
      edge_width = hs::lerp(a.edge_width, b.edge_width, t);
    }

    HS_COLD_MEMBER static float lerp_angle(float a, float b, float t) {
      if (t == 0.0f)
        return a;
      if (t == 1.0f)
        return b;
      float delta = fmodf(b - a + PI_F, TWO_PI_F);
      if (delta < 0.0f)
        delta += TWO_PI_F;
      return fmodf(a + (delta - PI_F) * t + TWO_PI_F, TWO_PI_F);
    }
  };

  struct WarpParams {
    WarpStageParams outer;
    WarpStageParams inner;

    HS_COLD_MEMBER bool operator==(const WarpParams &) const = default;

    HS_COLD_MEMBER void lerp(const WarpParams &a, const WarpParams &b,
                             float t) {
      outer.lerp(a.outer, b.outer, t);
      inner.lerp(a.inner, b.inner, t);
    }
  };

  struct ProjectionParams {
    float pole_fade = 1.0f;
    float spin_rate = 0.0f;
    float wander = 0.0f;
    float central_meridian = 0.0f;
    float coordinate_scale = 1.0f;
    float bonne_standard_parallel = PI_F * 0.25f;
    float layout_scroll = 0.0f;

    constexpr ProjectionParams() = default;

    constexpr ProjectionParams(float pole_fade, float spin_rate)
        : pole_fade(pole_fade), spin_rate(spin_rate) {}
    constexpr ProjectionParams(float pole_fade, float spin_rate, float wander)
        : pole_fade(pole_fade), spin_rate(spin_rate), wander(wander) {}

    HS_COLD_MEMBER bool operator==(const ProjectionParams &) const = default;

    HS_COLD_MEMBER static float lerp_periodic(float a, float b, float t) {
      if (t == 0.0f)
        return a;
      if (t == 1.0f)
        return b;
      float delta = fmodf(b - a + 0.5f, 1.0f);
      if (delta < 0.0f)
        delta += 1.0f;
      return a + (delta - 0.5f) * t;
    }

    HS_COLD_MEMBER void lerp(const ProjectionParams &a,
                             const ProjectionParams &b, float t) {
      pole_fade = hs::lerp(a.pole_fade, b.pole_fade, t);
      spin_rate = hs::lerp(a.spin_rate, b.spin_rate, t);
      wander = hs::lerp(a.wander, b.wander, t);
      central_meridian = WarpStageParams::lerp_angle(a.central_meridian,
                                                     b.central_meridian, t);
      coordinate_scale = hs::lerp(a.coordinate_scale, b.coordinate_scale, t);
      bonne_standard_parallel =
          hs::lerp(a.bonne_standard_parallel, b.bonne_standard_parallel, t);
      layout_scroll = lerp_periodic(a.layout_scroll, b.layout_scroll, t);
    }
  };

  struct SurfaceLensParams {
    float mix = 0.0f;
    float amount = 1.0f;
    float noise_scale = 1.0f;
    float noise_rate = 0.0f;
    NoiseBasis noise_basis = NoiseBasis::SIMPLEX;
    int32_t noise_seed = 4253;
    uint8_t noise_resource_id = 3;
    MobiusParams mobius{0.7071067811865475f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                        0.7071067811865475f, 0.0f};

    constexpr SurfaceLensParams() = default;

    constexpr SurfaceLensParams(float mix) : mix(mix) {}

    HS_COLD_MEMBER bool operator==(const SurfaceLensParams &other) const {
      return mix == other.mix && amount == other.amount &&
             noise_scale == other.noise_scale &&
             noise_rate == other.noise_rate &&
             noise_basis == other.noise_basis &&
             noise_seed == other.noise_seed &&
             noise_resource_id == other.noise_resource_id &&
             mobius.a.re == other.mobius.a.re &&
             mobius.a.im == other.mobius.a.im &&
             mobius.b.re == other.mobius.b.re &&
             mobius.b.im == other.mobius.b.im &&
             mobius.c.re == other.mobius.c.re &&
             mobius.c.im == other.mobius.c.im &&
             mobius.d.re == other.mobius.d.re &&
             mobius.d.im == other.mobius.d.im;
    }

    HS_COLD_MEMBER void lerp(const SurfaceLensParams &a,
                             const SurfaceLensParams &b, float t) {
      mix = hs::lerp(a.mix, b.mix, t);
      amount = hs::lerp(a.amount, b.amount, t);
      noise_scale = hs::lerp(a.noise_scale, b.noise_scale, t);
      noise_rate = hs::lerp(a.noise_rate, b.noise_rate, t);
      noise_basis = t < 1.0f ? a.noise_basis : b.noise_basis;
      noise_seed = t < 1.0f ? a.noise_seed : b.noise_seed;
      noise_resource_id = t < 1.0f ? a.noise_resource_id : b.noise_resource_id;
      mobius = t < 1.0f ? a.mobius : b.mobius;
    }
  };

  struct ValueParams {
    float iso_level = 0.5f;
    float iso_width = 0.05f;
    uint8_t band_count = 4;
    float band_phase = 0.0f;
    float cutout_threshold = 0.5f;
    float cutout_softness = 0.05f;
    float edge_width = 0.1f;

    HS_COLD_MEMBER bool operator==(const ValueParams &) const = default;
    HS_COLD_MEMBER void lerp(const ValueParams &a, const ValueParams &b,
                             float t) {
      iso_level = hs::lerp(a.iso_level, b.iso_level, t);
      iso_width = hs::lerp(a.iso_width, b.iso_width, t);
      band_count = t < 1.0f ? a.band_count : b.band_count;
      band_phase = WarpStageParams::lerp_angle(a.band_phase, b.band_phase, t);
      cutout_threshold = hs::lerp(a.cutout_threshold, b.cutout_threshold, t);
      cutout_softness = hs::lerp(a.cutout_softness, b.cutout_softness, t);
      edge_width = hs::lerp(a.edge_width, b.edge_width, t);
    }
  };

  struct ColorizerParams {
    float breathe_depth = 0.0f;
    float cycle_speed = 0.0f;
    float hue_shift = 0.0f;
    float value_fade = 0.0f;
    float displacement_gain = 1.0f;
    float path_gain = 0.0f;
    float direction_gain = 0.0f;
    float displacement_norm = 1.0f;
    float path_norm = 1.0f;
    float direction_phase = 0.0f;

    constexpr ColorizerParams() = default;

    constexpr ColorizerParams(float breathe_depth, float cycle_speed,
                              float hue_shift, float value_fade)
        : breathe_depth(breathe_depth), cycle_speed(cycle_speed),
          hue_shift(hue_shift), value_fade(value_fade) {}

    HS_COLD_MEMBER bool operator==(const ColorizerParams &) const = default;

    HS_COLD_MEMBER void lerp(const ColorizerParams &a, const ColorizerParams &b,
                             float t) {
      breathe_depth = hs::lerp(a.breathe_depth, b.breathe_depth, t);
      cycle_speed = hs::lerp(a.cycle_speed, b.cycle_speed, t);
      hue_shift = hs::lerp(a.hue_shift, b.hue_shift, t);
      value_fade = hs::lerp(a.value_fade, b.value_fade, t);
      displacement_gain = hs::lerp(a.displacement_gain, b.displacement_gain, t);
      path_gain = hs::lerp(a.path_gain, b.path_gain, t);
      direction_gain = hs::lerp(a.direction_gain, b.direction_gain, t);
      displacement_norm = hs::lerp(a.displacement_norm, b.displacement_norm, t);
      path_norm = hs::lerp(a.path_norm, b.path_norm, t);
      direction_phase =
          WarpStageParams::lerp_angle(a.direction_phase, b.direction_phase, t);
    }
  };

  struct OuterCameraParams {
    float wander;

    HS_COLD_MEMBER bool operator==(const OuterCameraParams &) const = default;

    HS_COLD_MEMBER void lerp(const OuterCameraParams &a,
                             const OuterCameraParams &b, float t) {
      wander = hs::lerp(a.wander, b.wander, t);
    }
  };

  struct Params {
    SourceParams source;
    WarpParams warp;
    ProjectionParams projection;
    SurfaceLensParams surface_lens;
    ValueParams value;
    ColorizerParams colorizer;
    OuterCameraParams outer_camera;

    HS_COLD_MEMBER bool operator==(const Params &) const = default;

    HS_COLD_MEMBER void lerp(const Params &a, const Params &b, float t) {
      source.lerp(a.source, b.source, t);
      warp.lerp(a.warp, b.warp, t);
      projection.lerp(a.projection, b.projection, t);
      surface_lens.lerp(a.surface_lens, b.surface_lens, t);
      value.lerp(a.value, b.value, t);
      colorizer.lerp(a.colorizer, b.colorizer, t);
      outer_camera.lerp(a.outer_camera, b.outer_camera, t);
    }

    HS_COLD_MEMBER void lerp_staggered(const Params &a, const Params &b,
                                       float t) {
      const int phase_count =
          (a.source != b.source) + (a.warp != b.warp) +
          (a.projection != b.projection) + (a.surface_lens != b.surface_lens) +
          (a.value != b.value) + (a.colorizer != b.colorizer) +
          (a.outer_camera != b.outer_camera);
      int phase = 0;
      source = a.source;
      warp = a.warp;
      projection = a.projection;
      surface_lens = a.surface_lens;
      value = a.value;
      colorizer = a.colorizer;
      outer_camera = a.outer_camera;
      if (a.source != b.source)
        source.lerp(a.source, b.source, phase_t(t, phase++, phase_count));
      if (a.warp != b.warp)
        warp.lerp(a.warp, b.warp, phase_t(t, phase++, phase_count));
      if (a.projection != b.projection)
        projection.lerp(a.projection, b.projection,
                        phase_t(t, phase++, phase_count));
      if (a.surface_lens != b.surface_lens)
        surface_lens.lerp(a.surface_lens, b.surface_lens,
                          phase_t(t, phase++, phase_count));
      if (a.value != b.value)
        value.lerp(a.value, b.value, phase_t(t, phase++, phase_count));
      if (a.colorizer != b.colorizer)
        colorizer.lerp(a.colorizer, b.colorizer,
                       phase_t(t, phase++, phase_count));
      if (a.outer_camera != b.outer_camera)
        outer_camera.lerp(a.outer_camera, b.outer_camera,
                          phase_t(t, phase, phase_count));
    }

    HS_COLD_MEMBER static float phase_t(float t, int phase, int phase_count) {
      return ease_in_out_sin(hs::clamp(t * phase_count - phase, 0.0f, 1.0f));
    }
  };

  struct Config {
    Slots slots;
    Params params;

    HS_COLD_MEMBER bool operator==(const Config &) const = default;
  };
  using RequestedConfig = Config;
  using Preset = Config;

  HS_COLD_MEMBER void rebind_parameters() {
    reset_parameters();
    Slots &slots = requested_config.slots;
    register_animated_param("Function", &slots.function, FUNCTION_OPTIONS,
                            FUNCTION_EXPORT_OPTIONS, NUM_FUNCTIONS);
    const bool polar_topology =
        slots.warp_program.outer.kind == WarpStageKind::POLAR_CHART ||
        slots.warp_program.inner.kind == WarpStageKind::POLAR_CHART;
    register_source_controls(slots.function, requested_config.params.source,
                             !polar_topology);
    register_animated_param("Projection", &slots.projection, PROJECTION_OPTIONS,
                            PROJECTION_EXPORT_OPTIONS, NUM_PROJECTIONS);
    register_projection_controls(slots, requested_config.params);
    register_animated_param(
        "Projection Frame", &slots.projection_frame, PROJECTION_FRAME_OPTIONS,
        PROJECTION_FRAME_EXPORT_OPTIONS, NUM_PROJECTION_FRAMES);
    register_projection_frame_controls(slots.projection_frame,
                                       requested_config.params);
    register_animated_param("Camera Wander",
                            &requested_config.params.outer_camera.wander,
                            WANDER_MIN, WANDER_MAX);
    register_animated_param("Lens", &slots.surface_lens, LENS_OPTIONS,
                            LENS_EXPORT_OPTIONS, NUM_LENSES);
    register_lens_controls(slots.surface_lens,
                           requested_config.params.surface_lens);
    register_animated_param("Planar Warp 1", &slots.warp_program.outer.kind,
                            WARP_OPTIONS, WARP_EXPORT_OPTIONS, NUM_WARPS);
    register_stage_slot_controls(true, slots.warp_program.outer);
    register_active_warp_controls(true, slots.warp_program.outer,
                                  requested_config.params.warp.outer);
    register_animated_param("Planar Warp 2", &slots.warp_program.inner.kind,
                            WARP_OPTIONS, WARP_EXPORT_OPTIONS, NUM_WARPS);
    register_stage_slot_controls(false, slots.warp_program.inner);
    register_active_warp_controls(false, slots.warp_program.inner,
                                  requested_config.params.warp.inner);
    register_animated_param("Signal Weight", &slots.signal_weight,
                            SIGNAL_OPTIONS, SIGNAL_EXPORT_OPTIONS, NUM_SIGNALS);
    register_animated_param("Value Transfer", &slots.value_transfer,
                            VALUE_TRANSFER_OPTIONS,
                            VALUE_TRANSFER_EXPORT_OPTIONS, NUM_VALUE_TRANSFERS);
    register_value_transfer_controls(slots.value_transfer,
                                     requested_config.params.value);
    register_animated_param("Coverage", &slots.coverage, COVERAGE_OPTIONS,
                            COVERAGE_EXPORT_OPTIONS, NUM_COVERAGE_POLICIES);
    register_coverage_controls(slots.coverage, requested_config.params.value);
    register_animated_param("Colorizer", &slots.colorizer, COLORIZER_OPTIONS,
                            COLORIZER_EXPORT_OPTIONS, NUM_COLORIZERS);
    register_colorizer_controls(slots.colorizer,
                                requested_config.params.colorizer);
#if HS_ENABLE_PARAM_GUI_BRIDGE
    mirror_parameter_display_state(requested_config, display_config);
    for (size_t index = 0; index < pending_edit_count; ++index)
      show_requested_parameter_value(pending_edits[index].name);
#endif
    requested_schema_bound = true;
  }

#if HS_ENABLE_PARAM_GUI_BRIDGE
  struct PendingEdit {
    const char *name = nullptr;
    size_t offset = 0;
    size_t size = 0;
  };

  static void dispatch_parameter_updated(Effect *effect, const char *name,
                                         bool is_enum) {
    static_cast<ShaderBall *>(effect)->parameter_updated(name, is_enum);
  }

  void parameter_updated(const char *name, bool is_enum) {
    const ParamDef *parameter = getParameters().find(name);
    HS_CHECK(parameter != nullptr, "updated ShaderBall parameter disappeared");
    const uintptr_t target = reinterpret_cast<uintptr_t>(parameter->target);
    const uintptr_t requested = reinterpret_cast<uintptr_t>(&requested_config);
    const size_t size = parameter_target_size(*parameter);
    HS_CHECK(target >= requested &&
                 target + size <= requested + sizeof(requested_config),
             "ShaderBall parameter target lies outside requested config");
    const size_t offset = target - requested;
    const size_t before_count = pending_edit_count;
    const bool was_pending = pending_edit_at(offset) < pending_edit_count;
    remember_pending_edit(name, offset, size);
    refresh_accepted_config();
    const bool is_pending = pending_edit_at(offset) < pending_edit_count;
    if (before_count != pending_edit_count || was_pending != is_pending ||
        (is_enum && schema_selector(name)))
      rebind_parameters();
  }

  static size_t parameter_target_size(const ParamDef &parameter) {
    switch (parameter.target_type) {
    case ParamDef::TargetType::FLOAT:
      return sizeof(float);
    case ParamDef::TargetType::INT_I32:
    case ParamDef::TargetType::INT_U32:
      return sizeof(uint32_t);
    case ParamDef::TargetType::INT_I16:
    case ParamDef::TargetType::INT_U16:
      return sizeof(uint16_t);
    case ParamDef::TargetType::BOOL:
    case ParamDef::TargetType::INT_I8:
    case ParamDef::TargetType::INT_U8:
      return sizeof(uint8_t);
    }
    __builtin_unreachable();
  }

  size_t pending_edit_at(size_t offset) const {
    for (size_t index = 0; index < pending_edit_count; ++index)
      if (pending_edits[index].offset == offset)
        return index;
    return pending_edit_count;
  }

  void remember_pending_edit(const char *name, size_t offset, size_t size) {
    const size_t existing = pending_edit_at(offset);
    if (existing < pending_edit_count) {
      pending_edits[existing].name = name;
      pending_edits[existing].size = size;
      return;
    }
    HS_CHECK(pending_edit_count < pending_edits.size(),
             "ShaderBall pending edit capacity exceeded");
    pending_edits[pending_edit_count++] = {name, offset, size};
  }

  void copy_pending_value(Config &to, const Config &from,
                          const PendingEdit &edit) const {
    std::memcpy(reinterpret_cast<uint8_t *>(&to) + edit.offset,
                reinterpret_cast<const uint8_t *>(&from) + edit.offset,
                edit.size);
  }

  void erase_pending_edit(size_t index) {
    for (size_t next = index + 1; next < pending_edit_count; ++next)
      pending_edits[next - 1] = pending_edits[next];
    --pending_edit_count;
  }

  void refresh_accepted_config() {
    if (valid_config(requested_config)) {
      accepted_config = requested_config;
      pending_edit_count = 0;
      return;
    }

    Config candidate = requested_config;
    for (size_t index = 0; index < pending_edit_count; ++index)
      copy_pending_value(candidate, accepted_config, pending_edits[index]);
    if (valid_config(candidate))
      accepted_config = candidate;

    bool admitted;
    do {
      admitted = false;
      for (size_t index = 0; index < pending_edit_count;) {
        candidate = accepted_config;
        copy_pending_value(candidate, requested_config, pending_edits[index]);
        if (!valid_config(candidate)) {
          ++index;
          continue;
        }
        accepted_config = candidate;
        erase_pending_edit(index);
        admitted = true;
      }
    } while (admitted);
  }

  const char *parameter_warning(const char *name) const override {
    for (size_t index = 0; index < pending_edit_count; ++index) {
      const PendingEdit &edit = pending_edits[index];
      if (std::strcmp(edit.name, name) != 0)
        continue;
      return admission_warning(requested_config, edit.name);
    }
    return nullptr;
  }

  float accepted_parameter_value(const ParamDef &parameter) const override {
    const uintptr_t target = reinterpret_cast<uintptr_t>(parameter.target);
    const uintptr_t requested = reinterpret_cast<uintptr_t>(&requested_config);
    const size_t size = parameter_target_size(parameter);
    if (target < requested ||
        target + size > requested + sizeof(requested_config))
      return parameter.get_requested();
    const size_t offset = target - requested;
    return parameter.get_from(
        reinterpret_cast<const uint8_t *>(&accepted_config) + offset);
  }

  static bool schema_selector(const char *name) {
    return std::strcmp(name, "Function") == 0 ||
           std::strcmp(name, "Projection") == 0 ||
           std::strcmp(name, "Peirce Layout") == 0 ||
           std::strcmp(name, "Airocean Layout") == 0 ||
           std::strcmp(name, "Bonne Hemisphere") == 0 ||
           std::strcmp(name, "Gnomonic Hemisphere") == 0 ||
           std::strcmp(name, "Projection Frame") == 0 ||
           std::strcmp(name, "Lens") == 0 ||
           std::strcmp(name, "Planar Warp 1") == 0 ||
           std::strcmp(name, "Planar Warp 2") == 0 ||
           std::strcmp(name, "Value Transfer") == 0 ||
           std::strcmp(name, "Coverage") == 0 ||
           std::strcmp(name, "Colorizer") == 0;
  }
#endif

  HS_COLD_MEMBER void register_value_transfer_controls(ValueTransfer transfer,
                                                       ValueParams &params) {
    if (transfer == ValueTransfer::ISO_CONTOUR) {
      register_animated_param("Iso Level", &params.iso_level, 0.0f, 1.0f);
      register_animated_param("Iso Width", &params.iso_width, SOFTNESS_MIN,
                              0.5f);
    } else if (transfer == ValueTransfer::SMOOTH_BANDS) {
      register_animated_int_param("Band Count", &params.band_count, 1,
                                  BAND_COUNT_MAX);
      register_animated_param("Band Phase", &params.band_phase, 0.0f, TWO_PI_F);
    }
  }

  HS_COLD_MEMBER void register_coverage_controls(CoveragePolicy coverage,
                                                 ValueParams &params) {
    if (coverage == CoveragePolicy::VALUE_CUTOUT) {
      register_animated_param("Cutout Threshold", &params.cutout_threshold,
                              0.0f, 1.0f);
      register_animated_param("Cutout Softness", &params.cutout_softness,
                              SOFTNESS_MIN, 0.5f);
    } else if (coverage == CoveragePolicy::EDGE_FADE) {
      register_animated_param("Edge Fade Width", &params.edge_width, 0.0f,
                              0.5f);
    }
  }

  HS_COLD_MEMBER void register_stage_slot_controls(bool first,
                                                   WarpStageSpec &spec) {
    if (spec.kind == WarpStageKind::VECTOR_NOISE ||
        spec.kind == WarpStageKind::CURL_FLOW) {
      register_animated_param(first ? "Planar Warp 1 Noise Basis"
                                    : "Planar Warp 2 Noise Basis",
                              &spec.basis, NOISE_BASIS_OPTIONS,
                              NOISE_BASIS_EXPORT_OPTIONS, NUM_NOISE_BASES);
      register_animated_param(first ? "Planar Warp 1 Envelope"
                                    : "Planar Warp 2 Envelope",
                              &spec.envelope, WARP_ENVELOPE_OPTIONS,
                              WARP_ENVELOPE_EXPORT_OPTIONS, NUM_WARP_ENVELOPES);
    }
    if (spec.kind == WarpStageKind::CURL_FLOW)
      register_animated_param(first ? "Planar Warp 1 Curl Integrator"
                                    : "Planar Warp 2 Curl Integrator",
                              &spec.curl_integrator, CURL_INTEGRATOR_OPTIONS,
                              CURL_INTEGRATOR_EXPORT_OPTIONS,
                              NUM_CURL_INTEGRATORS);
    if (spec.kind == WarpStageKind::POLAR_CHART) {
      register_animated_param(first ? "Planar Warp 1 Polar Mode"
                                    : "Planar Warp 2 Polar Mode",
                              &spec.polar_mode, POLAR_MODE_OPTIONS,
                              POLAR_MODE_EXPORT_OPTIONS, NUM_POLAR_MODES);
      register_animated_int_param(first ? "Planar Warp 1 Polar Harmonic"
                                        : "Planar Warp 2 Polar Harmonic",
                                  &spec.polar_harmonic, 1, POLAR_HARMONIC_MAX);
    }
  }

  HS_COLD_MEMBER void register_source_controls(Function function,
                                               SourceParams &params,
                                               bool expose_pattern_frequency) {
    if (function == Function::NOISE_CONTOUR) {
      register_animated_param("Source Noise Scale", &params.noise_scale,
                              SOURCE_NOISE_SCALE_MIN, SOURCE_NOISE_SCALE_MAX);
      register_animated_param("Source Noise Contrast", &params.noise_contrast,
                              0.0f, 8.0f);
      register_animated_param("Source Noise Rate", &params.noise_time_rate,
                              SOURCE_NOISE_RATE_MIN, SOURCE_NOISE_RATE_MAX);
      register_animated_param("Source Noise Basis", &params.noise_basis,
                              NOISE_BASIS_OPTIONS, NOISE_BASIS_EXPORT_OPTIONS,
                              NUM_NOISE_BASES);
      return;
    }
    if (function == Function::PRIMITIVE_LATTICE) {
      register_animated_param("Lattice Cell Scale", &params.lattice_cell_scale,
                              CELL_MIN, CELL_MAX);
      register_animated_param("Lattice Shape", &params.lattice_shape_blend,
                              0.0f, 1.0f);
      register_animated_param("Lattice Softness", &params.lattice_softness,
                              SOFTNESS_MIN, 1.0f);
      register_animated_param("Lattice Radius", &params.lattice_radius,
                              1.0f / 64.0f, 0.49f);
      return;
    }
    if (expose_pattern_frequency)
      register_animated_param("Pattern Freq", &params.pattern_freq,
                              PATTERN_FREQ_MIN, PATTERN_FREQ_MAX);
    register_animated_param("Speed", &params.speed, SPEED_MIN, SPEED_MAX);
    register_animated_param("Source Angle Rate", &params.angle_rate,
                            WAVE_SPIN_MIN, WAVE_SPIN_MAX);
    if (function == Function::GRID) {
      register_animated_param("Complexity", &params.complexity, COMPLEXITY_MIN,
                              COMPLEXITY_MAX);
      register_animated_param("Pattern Mix", &params.pattern_mix,
                              PATTERN_MIX_MIN, PATTERN_MIX_MAX);
      register_animated_param("Drift", &params.secondary_rate, PHASE2_RATE_MIN,
                              PHASE2_RATE_MAX);
    }
  }

  HS_COLD_MEMBER void register_projection_controls(Slots &slots,
                                                   Params &params) {
    if (slots.projection == Projection::PEIRCE_QUINCUNCIAL)
      register_animated_param("Peirce Layout", &slots.peirce_layout,
                              PEIRCE_LAYOUT_OPTIONS,
                              PEIRCE_LAYOUT_EXPORT_OPTIONS, NUM_PEIRCE_LAYOUTS);
    if (slots.projection == Projection::AIROCEAN)
      register_animated_param(
          "Airocean Layout", &slots.airocean_layout, AIROCEAN_LAYOUT_OPTIONS,
          AIROCEAN_LAYOUT_EXPORT_OPTIONS, NUM_AIROCEAN_LAYOUTS);
    if (slots.projection == Projection::SINUSOIDAL ||
        slots.projection == Projection::EQUIRECTANGULAR ||
        slots.projection == Projection::STEREOGRAPHIC ||
        slots.projection == Projection::GNOMONIC)
      register_animated_param("Pole Fade", &params.projection.pole_fade,
                              POLE_FADE_MIN, POLE_FADE_MAX);
    if (slots.projection == Projection::SINUSOIDAL ||
        slots.projection == Projection::EQUIRECTANGULAR ||
        strict_projection(slots.projection)) {
      register_animated_param("Central Meridian",
                              &params.projection.central_meridian, 0.0f,
                              TWO_PI_F);
    }
    if (strict_projection(slots.projection)) {
      register_animated_param("Projection Scale",
                              &params.projection.coordinate_scale, 0.25f, 4.0f);
    }
    if (slots.projection == Projection::BONNE)
      register_animated_param(
          "Bonne Hemisphere", &slots.bonne_hemisphere, BONNE_HEMISPHERE_OPTIONS,
          BONNE_HEMISPHERE_EXPORT_OPTIONS, NUM_BONNE_HEMISPHERES);
    if (slots.projection == Projection::GNOMONIC)
      register_animated_param("Gnomonic Hemisphere", &slots.gnomonic_hemisphere,
                              GNOMONIC_HEMISPHERE_OPTIONS,
                              GNOMONIC_HEMISPHERE_EXPORT_OPTIONS,
                              NUM_GNOMONIC_HEMISPHERES);
    if (slots.projection == Projection::BONNE)
      register_animated_param("Bonne Standard Parallel",
                              &params.projection.bonne_standard_parallel, 1e-3f,
                              0.5f * PI_F);
    if (slots.projection == Projection::PEIRCE_QUINCUNCIAL &&
        (slots.peirce_layout == PeirceLayout::HORIZONTAL ||
         slots.peirce_layout == PeirceLayout::VERTICAL))
      register_animated_param("Projection Layout Scroll",
                              &params.projection.layout_scroll, -1.0f, 1.0f);
  }

  HS_COLD_MEMBER void
  register_projection_frame_controls(ProjectionFramePolicy frame,
                                     Params &params) {
    if (frame == ProjectionFramePolicy::SPIN_WANDER) {
      register_animated_param("Spin Rate", &params.projection.spin_rate,
                              SPIN_RATE_MIN, SPIN_RATE_MAX);
      register_animated_param("Projection Wander", &params.projection.wander,
                              WANDER_MIN, WANDER_MAX);
    }
  }

  HS_COLD_MEMBER void register_lens_controls(SurfaceLens lens,
                                             SurfaceLensParams &params) {
    if (lens == SurfaceLens::NONE)
      return;
    register_animated_param("Lens Mix", &params.mix, LENS_MIX_MIN,
                            LENS_MIX_MAX);
    if (lens == SurfaceLens::TANGENT_NOISE) {
      register_animated_param("Lens Amount", &params.amount, 0.0f, 4.0f);
      register_animated_param("Lens Noise Scale", &params.noise_scale,
                              LENS_NOISE_SCALE_MIN, LENS_NOISE_SCALE_MAX);
      register_animated_param("Lens Noise Rate", &params.noise_rate,
                              NOISE_RATE_MIN, NOISE_RATE_MAX);
      register_animated_param("Lens Noise Basis", &params.noise_basis,
                              NOISE_BASIS_OPTIONS, NOISE_BASIS_EXPORT_OPTIONS,
                              NUM_NOISE_BASES);
    } else if (lens == SurfaceLens::MOBIUS) {
      register_animated_param("Mobius A Real", &params.mobius.a.re, -8.0f,
                              8.0f);
      register_animated_param("Mobius A Imag", &params.mobius.a.im, -8.0f,
                              8.0f);
      register_animated_param("Mobius B Real", &params.mobius.b.re, -8.0f,
                              8.0f);
      register_animated_param("Mobius B Imag", &params.mobius.b.im, -8.0f,
                              8.0f);
      register_animated_param("Mobius C Real", &params.mobius.c.re, -8.0f,
                              8.0f);
      register_animated_param("Mobius C Imag", &params.mobius.c.im, -8.0f,
                              8.0f);
      register_animated_param("Mobius D Real", &params.mobius.d.re, -8.0f,
                              8.0f);
      register_animated_param("Mobius D Imag", &params.mobius.d.im, -8.0f,
                              8.0f);
    }
  }

  HS_COLD_MEMBER void register_active_warp_controls(bool first,
                                                    const WarpStageSpec &spec,
                                                    WarpStageParams &params) {
    if (spec.kind == WarpStageKind::NONE)
      return;
    const char *const *names =
        first ? FIRST_WARP_PARAM_NAMES : SECOND_WARP_PARAM_NAMES;
    const char *time_name = first ? "Planar Warp 1 Time" : "Planar Warp 2 Time";
    auto register_current = [&](const char *name, float *target, float minimum,
                                float maximum) {
      register_animated_param(name, target,
                              *target < minimum ? *target : minimum,
                              *target > maximum ? *target : maximum);
    };
    if (spec.kind == WarpStageKind::LEGACY_STEREO_NOISE ||
        spec.kind == WarpStageKind::WAVE_SHEAR ||
        spec.kind == WarpStageKind::VECTOR_NOISE ||
        spec.kind == WarpStageKind::CURL_FLOW) {
      const char *strength_name =
          first ? "Planar Warp 1 Strength" : "Planar Warp 2 Strength";
      const bool signed_strength = spec.kind == WarpStageKind::WAVE_SHEAR ||
                                   spec.kind == WarpStageKind::CURL_FLOW;
      float strength_max = 4.0f;
      if (spec.kind == WarpStageKind::LEGACY_STEREO_NOISE) {
        strength_max = 30.0f;
      } else if (spec.kind == WarpStageKind::CURL_FLOW) {
        strength_max = curl_strength_limit(spec, params);
      }
      register_current(strength_name, &params.strength,
                       signed_strength ? -strength_max : 0.0f, strength_max);
    }
    if (spec.kind == WarpStageKind::LEGACY_STEREO_NOISE) {
      register_current(first ? "Planar Warp 1 Scale" : "Planar Warp 2 Scale",
                       &params.scale, 0.1f, 100.0f);
      register_current(time_name, &params.time_scale, 0.05f, 1.0f);
      return;
    }
    if (spec.kind == WarpStageKind::WAVE_SHEAR ||
        spec.kind == WarpStageKind::VORTEX ||
        spec.kind == WarpStageKind::VECTOR_NOISE ||
        spec.kind == WarpStageKind::CURL_FLOW)
      register_current(time_name, &params.time_scale, NOISE_RATE_MIN,
                       NOISE_RATE_MAX);
    switch (spec.kind) {
    case WarpStageKind::AFFINE_FRAME: {
      for (int index = 0; index < 6; ++index) {
        float *targets[] = {&params.translation_x, &params.translation_y,
                            &params.rotation,      &params.scale_x,
                            &params.scale_y,       &params.shear};
        const float minimum[] = {-4.0f, -4.0f, 0.0f, 0.25f, 0.25f, -0.75f};
        const float maximum[] = {4.0f, 4.0f, TWO_PI_F, 4.0f, 4.0f, 0.75f};
        register_current(names[WARP_NAME_TRANSLATION_X + index], targets[index],
                         minimum[index], maximum[index]);
      }
      break;
    }
    case WarpStageKind::WAVE_SHEAR:
      register_current(names[WARP_NAME_FREQUENCY], &params.frequency, 0.0f,
                       64.0f);
      register_current(names[WARP_NAME_FIELD_ANGLE], &params.field_angle, 0.0f,
                       TWO_PI_F);
      break;
    case WarpStageKind::VORTEX:
      register_current(names[WARP_NAME_CENTER_X], &params.center_x, -4.0f,
                       4.0f);
      register_current(names[WARP_NAME_CENTER_Y], &params.center_y, -4.0f,
                       4.0f);
      register_current(names[WARP_NAME_RADIUS], &params.radius, 1.0f / 64.0f,
                       8.0f);
      register_current(names[WARP_NAME_TURNS], &params.turns, -4.0f, 4.0f);
      register_current(names[WARP_NAME_CENTER_ORBIT],
                       &params.center_orbit_radius, 0.0f, 4.0f);
      break;
    case WarpStageKind::VECTOR_NOISE:
    case WarpStageKind::CURL_FLOW:
      register_current(first ? "Planar Warp 1 Scale" : "Planar Warp 2 Scale",
                       &params.scale, 1.0f / 64.0f,
                       spec.kind == WarpStageKind::CURL_FLOW ? 16.0f : 64.0f);
      register_current(names[WARP_NAME_VECTOR_ANGLE], &params.vector_angle,
                       0.0f, TWO_PI_F);
      register_current(names[WARP_NAME_EDGE_WIDTH], &params.edge_width,
                       SOFTNESS_MIN, 0.5f);
      break;
    case WarpStageKind::MIRROR_TILE:
      register_current(names[WARP_NAME_ROTATION], &params.rotation, 0.0f,
                       TWO_PI_F);
      register_current(names[WARP_NAME_CELL_X], &params.cell_x, CELL_MIN,
                       CELL_MAX);
      register_current(names[WARP_NAME_CELL_Y], &params.cell_y, CELL_MIN,
                       CELL_MAX);
      register_current(names[WARP_NAME_OFFSET_X], &params.offset_x, -8.0f,
                       8.0f);
      register_current(names[WARP_NAME_OFFSET_Y], &params.offset_y, -8.0f,
                       8.0f);
      break;
    case WarpStageKind::POLAR_CHART:
      register_current(names[WARP_NAME_RADIAL_SCALE], &params.radial_scale,
                       1.0f / 64.0f, 16.0f);
      register_current(names[WARP_NAME_RADIAL_PHASE], &params.radial_phase,
                       0.0f, TWO_PI_F);
      register_current(names[WARP_NAME_ANGULAR_PHASE], &params.angular_phase,
                       0.0f, TWO_PI_F);
      break;
    case WarpStageKind::NONE:
    case WarpStageKind::LEGACY_STEREO_NOISE:
      break;
    }
  }

  HS_COLD_MEMBER void register_colorizer_controls(Colorizer colorizer,
                                                  ColorizerParams &params) {
    if (colorizer == Colorizer::LIQUID) {
      register_animated_param("Breathe Depth", &params.breathe_depth,
                              BREATHE_MIN, BREATHE_MAX);
      register_animated_param("Cycle Speed", &params.cycle_speed,
                              CYCLE_SPEED_MIN, CYCLE_SPEED_MAX);
      register_animated_param("Hue Shift", &params.hue_shift, HUE_SHIFT_MIN,
                              HUE_SHIFT_MAX);
      register_animated_param("Value Fade", &params.value_fade, VALUE_FADE_MIN,
                              VALUE_FADE_MAX);
    } else if (colorizer == Colorizer::DEFORMATION_INK) {
      register_animated_param("Ink Displacement", &params.displacement_gain,
                              -4.0f, 4.0f);
      register_animated_param("Ink Path", &params.path_gain, -4.0f, 4.0f);
      register_animated_param("Ink Direction", &params.direction_gain, -4.0f,
                              4.0f);
      register_animated_param("Ink Displacement Norm",
                              &params.displacement_norm, SOFTNESS_MIN, 32.0f);
      register_animated_param("Ink Path Norm", &params.path_norm, SOFTNESS_MIN,
                              32.0f);
      register_animated_param("Ink Direction Phase", &params.direction_phase,
                              0.0f, TWO_PI_F);
    }
  }

  struct Blend {
    Params params;
  };

  struct Choreo {
    uint16_t dwell_min;
    uint16_t dwell_max;
    uint16_t blend_frames;
    bool staggered;
  };

  struct SourceState {
    float primary;
    float secondary;
    float angle;
    float angle_cos;
    float angle_sin;
  };

  struct ProjectedLookup {
    Complex coords;
    uint8_t region_id;
    uint8_t component_id;
    uint8_t boundary_flags;
    float fade_edge_distance;
    float value_weight;
    uint8_t flags;
    uint8_t traits;
    uint8_t edge_class;
    float domain_coverage;

    constexpr ProjectedLookup(Complex coords, uint8_t region_id,
                              uint8_t component_id, uint8_t boundary_flags,
                              float fade_edge_distance, float value_weight,
                              uint8_t flags, uint8_t traits = 0,
                              uint8_t edge_class = 0,
                              float domain_coverage = 1.0f)
        : coords(coords), region_id(region_id), component_id(component_id),
          boundary_flags(boundary_flags),
          fade_edge_distance(fade_edge_distance), value_weight(value_weight),
          flags(flags), traits(traits), edge_class(edge_class),
          domain_coverage(domain_coverage) {}
  };

  struct SourceTraits {
    bool cartesian_x;
    bool cartesian_y;
    bool x_periodic;
    bool y_periodic;
    bool rotation_equivariant;
    bool polar_angle_compatible;
  };

  struct PlanarWarpResult {
    Complex coords;
    Complex net_delta;
    float deformation;
    float path_length;
  };

  struct PlanarWarpStageResult {
    Complex coords;
    Complex delta;
    float deformation;
    float path_length;
  };

  struct MaterialSample {
    float value;
    float coverage;
    Complex net_delta;
    float deformation;
    float path_length;
  };

  struct ClockState {
    float source_primary = 0.0f;
    float source_secondary = 0.0f;
    float source_angle = 0.0f;
    float warp_time = 0.0f;
    float projection_spin = 0.0f;
    float breathe_phase = 0.0f;
    float source_noise_time = 0.0f;
    float lens_noise_time = 0.0f;
    float warp_outer_phase = 0.0f;
    float warp_inner_phase = 0.0f;

    constexpr ClockState() = default;
    constexpr ClockState(float source_primary, float source_secondary,
                         float source_angle, float warp_time,
                         float projection_spin, float breathe_phase)
        : source_primary(source_primary), source_secondary(source_secondary),
          source_angle(source_angle), warp_time(warp_time),
          projection_spin(projection_spin), breathe_phase(breathe_phase) {}
  };

  struct PreparedTransforms {
    Quaternion projection_conj;
    Quaternion outer_conj;
  };

  struct WrappedNoisePhase {
    float current_time;
    float previous_time;
    float mix;
    bool blends;
  };

  struct PreparedWarpStage {
    float rotation_cos;
    float rotation_sin;
    float mirror_offset_x;
    float mirror_offset_y;
    float vortex_center_x;
    float vortex_center_y;
    float vortex_radius_sq;
    float vortex_angle_numerator;
    WrappedNoisePhase noise_phase;
  };

  struct PreparedWarpProgram {
    PreparedWarpStage outer;
    PreparedWarpStage inner;
  };

  struct ResourceBindings {
    const FastNoiseLite *outer_warp_noise;
    const FastNoiseLite *inner_warp_noise;
    const FastNoiseLite *source_noise;
    const FastNoiseLite *lens_noise;
    const BakedPalette *generated_palette;
    const BakedPalette *liquid_palette;
  };

  static constexpr size_t MAX_NOISE_RESOURCES = 8;

  struct NoiseResourceKey {
    NoiseBasis basis;
    int32_t seed;
    float generator_frequency;
    uint8_t resource_id;
    uint8_t channel_version;
    uint8_t octave_version;
    uint8_t stencil_version;

    HS_COLD_MEMBER bool operator==(const NoiseResourceKey &) const = default;
  };

  struct FrameState {
    Slots slots;
    Params params;
    ClockState clocks;
    SourceState prepared_source;
    float breathe_offset;
    PreparedTransforms transforms;
    PreparedWarpProgram prepared_warp;
    WrappedNoisePhase source_noise_phase;
    WrappedNoisePhase lens_noise_phase;
    WrappedNoisePhase legacy_noise_phase;
    ResourceBindings resources;
  };

  struct FrameShader {
    const FrameState *frame;
    float alpha;

    HS_FLASH_MEMBER Color4 operator()(const Vector &view) const {
      Color4 color = shade(view, *frame);
      color.alpha *= alpha;
      return color;
    }
  };

  struct LookRuntime {
    ClockState clocks{};
    Quaternion projection_wander;
    Quaternion outer_wander;
    PreparedTransforms transforms;
  };

  struct WalkDeltas {
    Quaternion projection;
    Quaternion outer;
  };

  struct ParamMorphRuntime {
    Params from;
    Params to;
    uint16_t elapsed = 0;
    uint16_t duration = 0;
    bool staggered = false;
    bool continue_choreo = false;
    bool active = false;
  };

  struct TransitionRuntime {
    Config from_config;
    Config to_config;
    LookRuntime from_runtime;
    LookRuntime to_runtime;
    uint16_t elapsed = 0;
    uint16_t duration = 0;
    bool continue_choreo = false;
    bool active = false;
  };

  struct StateBundle {
    std::array<FastNoiseLite, MAX_NOISE_RESOURCES> noise_resources;
    std::array<NoiseResourceKey, MAX_NOISE_RESOURCES> prepared_noise_keys{};
    FastNoiseLite projection_walk_noise;
    FastNoiseLite outer_walk_noise;
    ParamMorphRuntime param_morph;
    TransitionRuntime transition;
  };

  struct ThroughClearPhase {
    float alpha;
    bool from_endpoint;
    bool clear;
  };

  HS_COLD_MEMBER static constexpr bool warp_uses_noise(WarpStageKind kind) {
    return kind == WarpStageKind::LEGACY_STEREO_NOISE ||
           kind == WarpStageKind::VECTOR_NOISE ||
           kind == WarpStageKind::CURL_FLOW;
  }

  HS_COLD_MEMBER static constexpr bool seam_sensitive_warp(WarpStageKind kind) {
    return kind == WarpStageKind::VECTOR_NOISE ||
           kind == WarpStageKind::CURL_FLOW;
  }

  HS_COLD_MEMBER static constexpr NoiseResourceKey
  warp_resource_key(const WarpStageSpec &spec) {
    return {
        spec.basis,
        spec.seed,
        spec.kind == WarpStageKind::LEGACY_STEREO_NOISE ? WARP_NOISE_FREQUENCY
                                                        : 1.0f,
        spec.resource_id,
        static_cast<uint8_t>(spec.kind),
        static_cast<uint8_t>(spec.basis == NoiseBasis::SIMPLEX ? 1 : 3),
        static_cast<uint8_t>(spec.kind == WarpStageKind::CURL_FLOW ? 1 : 0)};
  }

  HS_COLD_MEMBER static constexpr NoiseResourceKey
  source_resource_key(const Config &config) {
    return {
        config.params.source.noise_basis,
        config.params.source.noise_seed,
        1.0f,
        config.params.source.noise_resource_id,
        32,
        static_cast<uint8_t>(
            config.params.source.noise_basis == NoiseBasis::SIMPLEX ? 1 : 3),
        0};
  }

  HS_COLD_MEMBER static constexpr NoiseResourceKey
  lens_resource_key(const Config &config) {
    return {config.params.surface_lens.noise_basis,
            config.params.surface_lens.noise_seed,
            1.0f,
            config.params.surface_lens.noise_resource_id,
            48,
            static_cast<uint8_t>(config.params.surface_lens.noise_basis ==
                                         NoiseBasis::SIMPLEX
                                     ? 1
                                     : 3),
            0};
  }

  HS_COLD_MEMBER static constexpr bool
  append_resource_key(const NoiseResourceKey &key,
                      std::array<NoiseResourceKey, MAX_NOISE_RESOURCES> &keys,
                      size_t &count) {
    for (size_t index = 0; index < count; ++index)
      if (keys[index] == key)
        return true;
      else if (keys[index].resource_id == key.resource_id)
        return false;
    if (count == keys.size())
      return false;
    keys[count++] = key;
    return true;
  }

  HS_COLD_MEMBER static constexpr bool append_config_resource_keys(
      const Config &config,
      std::array<NoiseResourceKey, MAX_NOISE_RESOURCES> &keys, size_t &count) {
    if (warp_uses_noise(config.slots.warp_program.outer.kind) &&
        !append_resource_key(warp_resource_key(config.slots.warp_program.outer),
                             keys, count))
      return false;
    if (warp_uses_noise(config.slots.warp_program.inner.kind) &&
        !append_resource_key(warp_resource_key(config.slots.warp_program.inner),
                             keys, count))
      return false;
    if (config.slots.function == Function::NOISE_CONTOUR &&
        !append_resource_key(source_resource_key(config), keys, count))
      return false;
    if (config.slots.surface_lens == SurfaceLens::TANGENT_NOISE &&
        !append_resource_key(lens_resource_key(config), keys, count))
      return false;
    return true;
  }

  HS_COLD_MEMBER static constexpr bool resource_union_fits(const Config &from,
                                                           const Config &to) {
    std::array<NoiseResourceKey, MAX_NOISE_RESOURCES> keys{};
    size_t count = 0;
    return append_config_resource_keys(from, keys, count) &&
           append_config_resource_keys(to, keys, count);
  }

  HS_COLD_MEMBER bool prepare_resource_union(const Config &from,
                                             const Config &to) {
    std::array<NoiseResourceKey, MAX_NOISE_RESOURCES> keys{};
    size_t count = 0;
    if (!append_config_resource_keys(from, keys, count) ||
        !append_config_resource_keys(to, keys, count))
      return false;
    prepared_noise_count = count;
    for (size_t index = 0; index < count; ++index) {
      state->prepared_noise_keys[index] = keys[index];
      state->noise_resources[index].SetNoiseType(
          FastNoiseLite::NoiseType_OpenSimplex2);
      state->noise_resources[index].SetSeed(keys[index].seed);
      state->noise_resources[index].SetFrequency(
          keys[index].generator_frequency);
    }
    return true;
  }

  HS_COLD_MEMBER const FastNoiseLite *
  resolve_resource(const NoiseResourceKey &key) const {
    for (size_t index = 0; index < prepared_noise_count; ++index)
      if (state->prepared_noise_keys[index] == key)
        return &state->noise_resources[index];
    return nullptr;
  }

  HS_COLD_MEMBER const FastNoiseLite *
  resolve_warp_resource(const WarpStageSpec &spec) const {
    return warp_uses_noise(spec.kind)
               ? resolve_resource(warp_resource_key(spec))
               : nullptr;
  }

  HS_COLD_MEMBER const FastNoiseLite *
  resolve_source_resource(const Config &config) const {
    return config.slots.function == Function::NOISE_CONTOUR
               ? resolve_resource(source_resource_key(config))
               : nullptr;
  }

  HS_COLD_MEMBER const FastNoiseLite *
  resolve_lens_resource(const Config &config) const {
    return config.slots.surface_lens == SurfaceLens::TANGENT_NOISE
               ? resolve_resource(lens_resource_key(config))
               : nullptr;
  }

  HS_FLASH_MEMBER static WrappedNoisePhase prepare_noise_phase(float turns) {
    const float wrapped_turns = wrap_t(turns);
    const float current_time = wrapped_turns * NOISE_NATIVE_PERIOD;
    if (wrapped_turns <= NOISE_WRAP_START)
      return {current_time, current_time, 0.0f, false};
    return {current_time, current_time - NOISE_NATIVE_PERIOD,
            ease_in_out_sin((wrapped_turns - NOISE_WRAP_START) /
                            (1.0f - NOISE_WRAP_START)),
            true};
  }

  HS_FLASH_MEMBER static WrappedNoisePhase
  prepare_legacy_noise_phase(float time) {
    const float blend_start = STEREO_NOISE_TIME_PERIOD - NOISE_WRAP_BLEND;
    if (time <= blend_start)
      return {time, time, 0.0f, false};
    return {time, time - STEREO_NOISE_TIME_PERIOD,
            ease_in_out_sin((time - blend_start) / NOISE_WRAP_BLEND), true};
  }

  HS_FLASH_MEMBER static PreparedWarpStage
  prepare_warp_stage(const WarpStageParams &params, float stage_phase) {
    const float rotation_cos = cosf(params.rotation);
    const float rotation_sin = sinf(params.rotation);
    const float orbit_phase = TWO_PI_F * stage_phase;
    return {
        rotation_cos,
        rotation_sin,
        wrap_t(params.offset_x / params.cell_x) * params.cell_x,
        wrap_t(params.offset_y / params.cell_y) * params.cell_y,
        params.center_x + params.center_orbit_radius * cosf(orbit_phase),
        params.center_y + params.center_orbit_radius * sinf(orbit_phase),
        params.radius * params.radius,
        TWO_PI_F * params.turns,
        prepare_noise_phase(stage_phase),
    };
  }

  HS_COLD_MEMBER FrameState prepare_frame() const {
    return prepare_frame({active_slots, blend.params}, runtime);
  }

  HS_COLD_MEMBER FrameState prepare_frame(const Config &config,
                                          const LookRuntime &look) const {
    const bool animated_projection =
        config.slots.projection_frame == ProjectionFramePolicy::SPIN_WANDER;
    const PreparedWarpProgram prepared_warp{
        prepare_warp_stage(config.params.warp.outer,
                           look.clocks.warp_outer_phase),
        prepare_warp_stage(config.params.warp.inner,
                           look.clocks.warp_inner_phase)};
    return {
        config.slots,
        config.params,
        look.clocks,
        {look.clocks.source_primary, look.clocks.source_secondary,
         look.clocks.source_angle, fast_cosf(look.clocks.source_angle),
         fast_sinf(look.clocks.source_angle)},
        fast_sinf(look.clocks.breathe_phase) *
            config.params.colorizer.breathe_depth,
        {animated_projection ? look.transforms.projection_conj : Quaternion(),
         look.transforms.outer_conj},
        prepared_warp,
        prepare_noise_phase(look.clocks.source_noise_time),
        prepare_noise_phase(look.clocks.lens_noise_time),
        prepare_legacy_noise_phase(look.clocks.warp_time),
        {resolve_warp_resource(config.slots.warp_program.outer),
         resolve_warp_resource(config.slots.warp_program.inner),
         resolve_source_resource(config), resolve_lens_resource(config),
         &generated_palette_cycler.palette(),
         &liquid_palette_cycler.palette()}};
  }

  static ThroughClearPhase through_clear_phase(uint16_t elapsed,
                                               uint16_t duration) {
    const uint16_t center = duration / 2;
    if (elapsed == center)
      return {0.0f, false, true};
    const bool from_endpoint = elapsed < center;
    const float phase = from_endpoint ? static_cast<float>(elapsed) / center
                                      : static_cast<float>(elapsed - center) /
                                            (duration - center);
    return {from_endpoint ? 1.0f - ease_in_out_sin(phase)
                          : ease_in_out_sin(phase),
            from_endpoint, false};
  }

  HS_FLASH_MEMBER void draw_through_clear_transition(Canvas &canvas) const {
    const ThroughClearPhase phase = through_clear_phase(
        state->transition.elapsed, state->transition.duration);
    if (phase.clear)
      return;
    const FrameState visible =
        phase.from_endpoint ? prepare_frame(state->transition.from_config,
                                            state->transition.from_runtime)
                            : prepare_frame(state->transition.to_config,
                                            state->transition.to_runtime);
    FrameShader shader{&visible, phase.alpha};
    HS_PROFILE(sb_shader_draw);
    Scan::Shader::draw<W, H, 1>(canvas, shader);
  }

  /**
   * @brief Shades one sphere sample by pulling it back to a source coordinate.
   * @param view Unit direction of the visible sphere point.
   * @param frame Immutable snapshot of slots, parameters, clocks, transforms,
   *        and palette resources for this frame.
   * @return Premultiplied-alpha colour for the sample.
   * @details Walks outer camera, surface lens, and projection backward. A
   * strict projection whose two lens branches land in different regions cannot
   * be joined in the plane, so the branches are shaded separately and their
   * outputs blended instead.
   */
  static Color4 shade(const Vector &view, const FrameState &frame) {
    const Vector outer_local = outer_camera_lookup(view, frame);
    if (strict_projection(frame.slots.projection) &&
        frame.slots.surface_lens != SurfaceLens::NONE &&
        frame.params.surface_lens.mix > 0.0f &&
        frame.params.surface_lens.mix < 1.0f)
      return shade_strict_lens_mix(outer_local, frame);
    const ProjectedLookup projected =
        surface_lens_project_lookup(outer_local, frame);
    return shade_projected(projected, frame);
  }

  HS_FLASH_MEMBER static Color4 shade_strict_lens_mix(const Vector &outer_local,
                                                      const FrameState &frame) {
    const ProjectedLookup direct = profiled_project_branch(outer_local, frame);
    const Vector lensed_direction = profiled_apply_lens(outer_local, frame);
    const ProjectedLookup lensed =
        profiled_project_branch(lensed_direction, frame);
    if (!projection_join_compatible(direct, lensed, frame.slots.projection,
                                    frame.params.projection.coordinate_scale))
      return blend_outputs(shade_projected(direct, frame),
                           shade_projected(lensed, frame),
                           frame.params.surface_lens.mix);
    return shade_projected(join_projected(direct, lensed,
                                          frame.params.surface_lens.mix,
                                          frame.slots.projection,
                                          frame.params.projection.pole_fade),
                           frame);
  }

  /**
   * @brief Runs the planar half of the pullback: warps, samples, shapes, and
   *        colorizes.
   * @param projected Plane coordinates and seam metadata for the sample.
   * @param frame Frame snapshot.
   * @return Premultiplied-alpha colour for the sample.
   */
  HS_FLASH_MEMBER static Color4
  shade_projected(const ProjectedLookup &projected, const FrameState &frame) {
    HS_SB_STAGE_MARK(stage_start);
    const PlanarWarpResult warped = planar_warp_lookup(projected, frame);
    HS_SB_STAGE_SPAN(planar_warp, stage_start);
    const Complex source_coords = condition_source_coords(warped.coords, frame);
    const float field = sample_source(source_coords, frame);
    HS_SB_STAGE_SPAN(source, stage_start);
    const MaterialSample material =
        shape_material(field, projected, warped, frame);
    HS_SB_STAGE_SPAN(material, stage_start);
    const Color4 color = colorize(material, frame);
    HS_SB_STAGE_SPAN(color, stage_start);
    return color;
  }

  static constexpr bool strict_projection(Projection projection) {
    return projection == Projection::BONNE ||
           projection == Projection::PEIRCE_QUINCUNCIAL ||
           projection == Projection::AIROCEAN;
  }

  static bool projection_edge_distance_required(const FrameState &frame) {
    const WarpProgram &program = frame.slots.warp_program;
    return frame.slots.coverage == CoveragePolicy::EDGE_FADE ||
           (program.outer.kind != WarpStageKind::NONE &&
            program.outer.envelope == WarpEnvelope::EDGE_FADE) ||
           (program.inner.kind != WarpStageKind::NONE &&
            program.inner.envelope == WarpEnvelope::EDGE_FADE);
  }

  static bool projection_join_compatible(const ProjectedLookup &a,
                                         const ProjectedLookup &b,
                                         Projection projection,
                                         float coordinate_scale = 1.0f) {
    if (strict_projection(projection)) {
      (void)a;
      (void)b;
      (void)coordinate_scale;
      return false;
    }
    return a.component_id == b.component_id && a.flags == b.flags &&
           ((a.boundary_flags | b.boundary_flags) &
            (BOUNDARY_CUT | BOUNDARY_SINGULAR)) == 0;
  }

  HS_FLASH_MEMBER static Color4 blend_outputs(const Color4 &from,
                                              const Color4 &to, float mix) {
    if (mix == 0.0f)
      return from;
    if (mix == 1.0f)
      return to;
    const float alpha = hs::lerp(from.alpha, to.alpha, mix);
    if (alpha == 0.0f)
      return Color4();
    const float from_weight = from.alpha * (1.0f - mix);
    const float to_weight = to.alpha * mix;
    const float inv_alpha = 1.0f / alpha;
    return Color4(
        Pixel(static_cast<uint16_t>(hs::clamp(
                  (from.color.r * from_weight + to.color.r * to_weight) *
                          inv_alpha +
                      0.5f,
                  0.0f, 65535.0f)),
              static_cast<uint16_t>(hs::clamp(
                  (from.color.g * from_weight + to.color.g * to_weight) *
                          inv_alpha +
                      0.5f,
                  0.0f, 65535.0f)),
              static_cast<uint16_t>(hs::clamp(
                  (from.color.b * from_weight + to.color.b * to_weight) *
                          inv_alpha +
                      0.5f,
                  0.0f, 65535.0f))),
        alpha);
  }

  static Vector outer_camera_lookup(const Vector &view,
                                    const FrameState &frame) {
    return rotate(view, frame.transforms.outer_conj);
  }

  static ProjectedLookup surface_lens_project_lookup(const Vector &v,
                                                     const FrameState &frame) {
    const Slots &slots = frame.slots;
    const float mix = frame.params.surface_lens.mix;
    if (slots.surface_lens == SurfaceLens::NONE || mix == 0.0f)
      return profiled_project_branch(v, frame);
    const Vector lensed = profiled_apply_lens(v, frame);
    if (mix == 1.0f)
      return profiled_project_branch(lensed, frame);
    return surface_lens_mixed_lookup(v, lensed, frame, mix);
  }

  HS_FLASH_MEMBER static ProjectedLookup
  surface_lens_mixed_lookup(const Vector &v, const Vector &lensed,
                            const FrameState &frame, float mix) {
    const ProjectedLookup direct = profiled_project_branch(v, frame);
    const ProjectedLookup distorted = profiled_project_branch(lensed, frame);
    return join_projected(direct, distorted, mix, frame.slots.projection,
                          frame.params.projection.pole_fade);
  }

  /**
   * @brief Rescales a projection kernel's plane coordinates into a lookup.
   * @param result Kernel output in the projection's native plane units.
   * @param coordinate_scale Plane-unit scale; applied to the coordinates and,
   *        as a magnitude, to the edge distance.
   * @return The kernel result carried over with scaled coordinates and a
   *         saturated value weight.
   */
  __attribute__((always_inline)) static ProjectedLookup
  scaled_kernel_lookup(const projections::ProjectionKernelResult &result,
                       float coordinate_scale) {
    return {{result.coords.re * coordinate_scale,
             result.coords.im * coordinate_scale},
            result.region_id,
            result.component_id,
            result.boundary_flags,
            result.fade_edge_distance * fabsf(coordinate_scale),
            1.0f,
            result.flags,
            result.traits,
            result.edge_class};
  }

  HS_FLASH_MEMBER static ProjectedLookup
  project_bonne(const Vector &local, const FrameState &frame) {
    return scaled_kernel_lookup(
        projections::bonne_projection(
            local, frame.params.projection.central_meridian,
            (frame.slots.bonne_hemisphere == BonneHemisphere::NORTH ? 1.0f
                                                                    : -1.0f) *
                frame.params.projection.bonne_standard_parallel),
        frame.params.projection.coordinate_scale);
  }

  HS_FLASH_MEMBER static ProjectedLookup
  project_peirce(const Vector &local, const FrameState &frame) {
    if (frame.slots.peirce_layout == PeirceLayout::SQUARE &&
        frame.params.projection.central_meridian == 0.0f &&
        projection_edge_distance_required(frame))
      return scaled_kernel_lookup(
          projections::peirce_projection_fast_square(local),
          frame.params.projection.coordinate_scale);
    return scaled_kernel_lookup(
        projections::peirce_projection(
            local, frame.params.projection.central_meridian,
            static_cast<uint8_t>(frame.slots.peirce_layout),
            frame.params.projection.layout_scroll,
            projection_edge_distance_required(frame)),
        frame.params.projection.coordinate_scale);
  }

  HS_FLASH_MEMBER static ProjectedLookup
  project_airocean(const Vector &local, const FrameState &frame) {
    return scaled_kernel_lookup(
        projections::airocean_projection(
            local, frame.params.projection.central_meridian,
            frame.slots.airocean_layout == AiroceanLayout::HORIZONTAL,
            projection_edge_distance_required(frame)),
        frame.params.projection.coordinate_scale);
  }

  HS_FLASH_MEMBER static ProjectedLookup
  project_sinusoidal(const Vector &local, const FrameState &frame) {
    return finalize_projection(
        local,
        folded_sinusoidal(local, frame.params.projection.central_meridian),
        Projection::SINUSOIDAL, frame.params.projection.pole_fade);
  }

  HS_FLASH_MEMBER static ProjectedLookup
  project_equirectangular(const Vector &local, const FrameState &frame) {
    return finalize_projection(
        local, equirectangular(local, frame.params.projection.central_meridian),
        Projection::EQUIRECTANGULAR, frame.params.projection.pole_fade);
  }

  HS_FLASH_MEMBER static ProjectedLookup
  project_gnomonic(const Vector &local, const FrameState &frame) {
    return finalize_projection(local, gnomonic(local), Projection::GNOMONIC,
                               frame.params.projection.pole_fade,
                               frame.slots.gnomonic_hemisphere);
  }

  HS_FLASH_MEMBER static ProjectedLookup
  project_nonstereographic(const Vector &local, const FrameState &frame) {
    if (frame.slots.projection == Projection::BONNE)
      return project_bonne(local, frame);
    if (frame.slots.projection == Projection::PEIRCE_QUINCUNCIAL)
      return project_peirce(local, frame);
    if (frame.slots.projection == Projection::AIROCEAN)
      return project_airocean(local, frame);
    if (frame.slots.projection == Projection::SINUSOIDAL)
      return project_sinusoidal(local, frame);
    if (frame.slots.projection == Projection::EQUIRECTANGULAR)
      return project_equirectangular(local, frame);
    if (frame.slots.projection == Projection::GNOMONIC)
      return project_gnomonic(local, frame);

    __builtin_unreachable();
  }

  HS_FLASH_MEMBER static ProjectedLookup
  project_branch(const Vector &v, const FrameState &frame) {
    const Vector local = rotate(v, frame.transforms.projection_conj);
    if (frame.slots.projection != Projection::STEREOGRAPHIC)
      return project_nonstereographic(local, frame);
    const Complex coords = stereo(local);
    const float r_sq = coords.re * coords.re + coords.im * coords.im;
    return {coords,
            0,
            0,
            BOUNDARY_SINGULAR,
            std::max(0.0f, 1.0f - local.y),
            pole_attenuation(r_sq, frame.params.projection.pole_fade),
            0};
  }

  __attribute__((always_inline)) static ProjectedLookup
  profiled_project_branch(const Vector &v, const FrameState &frame) {
    HS_SB_STAGE_MARK(stage_start);
    const ProjectedLookup projected = project_branch(v, frame);
    HS_SB_STAGE_SPAN(projection, stage_start);
    return projected;
  }

  HS_FLASH_MEMBER static ProjectedLookup
  finalize_projection(const Vector &local, const Complex &coords,
                      Projection projection, float pole_fade,
                      GnomonicHemispherePolicy gnomonic_hemisphere =
                          GnomonicHemispherePolicy::FOLDED) {
    const float r_sq = coords.re * coords.re + coords.im * coords.im;
    switch (projection) {
    case Projection::SINUSOIDAL:
      return {coords,
              static_cast<uint8_t>(local.z < 0.0f),
              0,
              0,
              1.0f,
              pole_attenuation(r_sq, pole_fade),
              PROJECTION_FLAG_FOLDED};
    case Projection::EQUIRECTANGULAR:
      return {coords,
              0,
              0,
              BOUNDARY_CUT,
              PI_F - std::fabs(coords.re),
              pole_attenuation(r_sq, pole_fade),
              0};
    case Projection::STEREOGRAPHIC:
      return {coords,
              0,
              0,
              BOUNDARY_SINGULAR,
              std::max(0.0f, 1.0f - local.y),
              pole_attenuation(r_sq, pole_fade),
              0};
    case Projection::GNOMONIC: {
      const bool in_domain =
          gnomonic_hemisphere == GnomonicHemispherePolicy::FOLDED ||
          (gnomonic_hemisphere == GnomonicHemispherePolicy::FRONT_HEMISPHERE
               ? local.y >= 0.0f
               : local.y < 0.0f);
      return {coords,
              static_cast<uint8_t>(local.y < 0.0f),
              static_cast<uint8_t>(local.y < 0.0f),
              static_cast<uint8_t>(BOUNDARY_CUT | BOUNDARY_SINGULAR),
              std::fabs(local.y),
              pole_attenuation(r_sq, pole_fade),
              0,
              0,
              0,
              in_domain ? 1.0f : 0.0f};
    }
    case Projection::BONNE:
    case Projection::PEIRCE_QUINCUNCIAL:
    case Projection::AIROCEAN:
      break;
    }
    __builtin_unreachable();
  }

  HS_FLASH_MEMBER static ProjectedLookup
  join_projected(const ProjectedLookup &direct, const ProjectedLookup &lensed,
                 float mix, Projection projection, float pole_fade) {
    if (mix == 0.0f)
      return direct;
    if (mix == 1.0f)
      return lensed;
    const Complex coords(hs::lerp(direct.coords.re, lensed.coords.re, mix),
                         hs::lerp(direct.coords.im, lensed.coords.im, mix));
    const ProjectedLookup *selected = nullptr;
    switch (projection) {
    case Projection::SINUSOIDAL:
      selected = mix < 0.5f ? &direct : &lensed;
      break;
    case Projection::EQUIRECTANGULAR:
      selected = mix < 0.5f ? &direct : &lensed;
      break;
    case Projection::STEREOGRAPHIC:
      selected = mix < 0.5f ? &direct : &lensed;
      break;
    case Projection::GNOMONIC:
      selected = mix < 0.5f ? &direct : &lensed;
      break;
    case Projection::BONNE:
    case Projection::PEIRCE_QUINCUNCIAL:
    case Projection::AIROCEAN:
      HS_CHECK(false, "strict projection joins require complete-output blend");
      __builtin_unreachable();
    }
    const float r_sq = coords.re * coords.re + coords.im * coords.im;
    return {coords,
            selected->region_id,
            selected->component_id,
            selected->boundary_flags,
            selected->fade_edge_distance,
            pole_attenuation(r_sq, pole_fade),
            selected->flags,
            selected->traits,
            selected->edge_class,
            hs::lerp(direct.domain_coverage, lensed.domain_coverage, mix)};
  }

  /**
   * @brief Pulls plane coordinates back through both warp stages.
   * @param projected Projection output; supplies the stage input coordinates
   *        and the weight and edge distance the stage envelopes read.
   * @param frame Frame snapshot.
   * @return Source-side coordinates plus the summed delta, deformation, and
   *         path length the colorizers consume.
   * @details Pullback order is Planar Warp 1 then Planar Warp 2, the reverse
   * of the authored order `source -> Warp 2 -> Warp 1 -> projection`.
   * Deformation is the
   * magnitude of the summed delta, except when a lone legacy stereo stage is
   * programmed and its own displacement is reported.
   */
  HS_FLASH_MEMBER static PlanarWarpResult
  planar_warp_lookup(const ProjectedLookup &projected,
                     const FrameState &frame) {
    const PlanarWarpStageResult outer = warp_stage_lookup(
        projected.coords, projected, frame.slots.warp_program.outer,
        frame.params.warp.outer, frame.clocks.warp_outer_phase,
        frame.resources.outer_warp_noise, frame.prepared_warp.outer, frame);
    const PlanarWarpStageResult inner = warp_stage_lookup(
        outer.coords, projected, frame.slots.warp_program.inner,
        frame.params.warp.inner, frame.clocks.warp_inner_phase,
        frame.resources.inner_warp_noise, frame.prepared_warp.inner, frame);
    const Complex net_delta(outer.delta.re + inner.delta.re,
                            outer.delta.im + inner.delta.im);
    const bool sole_legacy =
        (frame.slots.warp_program.outer.kind ==
             WarpStageKind::LEGACY_STEREO_NOISE &&
         frame.slots.warp_program.inner.kind == WarpStageKind::NONE) ||
        (frame.slots.warp_program.inner.kind ==
             WarpStageKind::LEGACY_STEREO_NOISE &&
         frame.slots.warp_program.outer.kind == WarpStageKind::NONE);
    const float deformation =
        sole_legacy
            ? outer.deformation + inner.deformation
            : sqrtf(net_delta.re * net_delta.re + net_delta.im * net_delta.im);
    return {inner.coords, net_delta, deformation,
            outer.path_length + inner.path_length};
  }

  HS_FLASH_MEMBER static PlanarWarpStageResult
  finish_closed_form_warp(const Complex &input, const Complex &output) {
    const Complex delta(output.re - input.re, output.im - input.im);
    const float deformation = sqrtf(delta.re * delta.re + delta.im * delta.im);
    return {output, delta, deformation, deformation};
  }

  HS_FLASH_MEMBER static PlanarWarpStageResult
  warp_affine_frame(const Complex &input, const WarpStageParams &params,
                    const PreparedWarpStage &prepared) {
    const float c = prepared.rotation_cos;
    const float s = prepared.rotation_sin;
    const float x = input.re - params.translation_x;
    const float y = input.im - params.translation_y;
    const float rx = c * x + s * y;
    const float ry = -s * x + c * y;
    const Complex output(rx / params.scale_x -
                             params.shear * ry / params.scale_y,
                         ry / params.scale_y);
    return finish_closed_form_warp(input, output);
  }

  HS_FLASH_MEMBER static PlanarWarpStageResult
  warp_wave_shear(const Complex &input, const WarpStageParams &params,
                  float stage_phase, float amplitude) {
    if (params.strength == 0.0f)
      return {input, Complex(), 0.0f, 0.0f};
    const float c = cosf(params.field_angle);
    const float s = sinf(params.field_angle);
    const float phase = params.frequency * (c * input.re + s * input.im) +
                        TWO_PI_F * stage_phase;
    const float offset = amplitude * sinf(phase);
    const Complex delta(-s * offset, c * offset);
    return {{input.re + delta.re, input.im + delta.im},
            delta,
            fabsf(offset),
            fabsf(offset)};
  }

  HS_FLASH_MEMBER static PlanarWarpStageResult
  warp_vortex(const Complex &input, const PreparedWarpStage &prepared) {
    const float x = input.re - prepared.vortex_center_x;
    const float y = input.im - prepared.vortex_center_y;
    const float r_sq = x * x + y * y;
    const float angle = prepared.vortex_angle_numerator /
                        (1.0f + r_sq / prepared.vortex_radius_sq);
    const float c = cosf(angle);
    const float s = sinf(angle);
    const Complex output(prepared.vortex_center_x + c * x - s * y,
                         prepared.vortex_center_y + s * x + c * y);
    return finish_closed_form_warp(input, output);
  }

  HS_FLASH_MEMBER static PlanarWarpStageResult
  warp_vector_noise(const Complex &input, const WarpStageSpec &spec,
                    const WarpStageParams &params, float amplitude,
                    const FastNoiseLite &noise,
                    const PreparedWarpStage &prepared) {
    if (params.strength == 0.0f)
      return {input, Complex(), 0.0f, 0.0f};
    const float nx = sample_wrapped_noise_basis(
        noise, spec.basis, params.scale * input.re, params.scale * input.im,
        prepared.noise_phase);
    float ny = sample_wrapped_noise_basis(
        noise, spec.basis, params.scale * input.re + 31.416f,
        params.scale * input.im - 47.853f, prepared.noise_phase, 11.0f);
    float zero_dc_x = nx;
    if (spec.basis == NoiseBasis::RIDGED3) {
      zero_dc_x -= sample_wrapped_noise_basis(
          noise, spec.basis, params.scale * input.re - 73.271f,
          params.scale * input.im + 19.119f, prepared.noise_phase, 5.0f);
      ny -= sample_wrapped_noise_basis(
          noise, spec.basis, params.scale * input.re + 61.731f,
          params.scale * input.im + 89.417f, prepared.noise_phase, -7.0f);
    }
    const float c = cosf(params.vector_angle);
    const float s = sinf(params.vector_angle);
    const Complex delta(amplitude * (c * zero_dc_x - s * ny),
                        amplitude * (s * zero_dc_x + c * ny));
    const float deformation = sqrtf(delta.re * delta.re + delta.im * delta.im);
    return {{input.re + delta.re, input.im + delta.im},
            delta,
            deformation,
            deformation};
  }

  HS_FLASH_MEMBER static PlanarWarpStageResult
  warp_curl_flow(const Complex &input, const WarpStageSpec &spec,
                 const WarpStageParams &params, float amplitude,
                 const FastNoiseLite &noise,
                 const PreparedWarpStage &prepared) {
    if (params.strength == 0.0f)
      return {input, Complex(), 0.0f, 0.0f};
    Complex delta;
    float path_length = 0.0f;
    const Complex output = curl_flow(input, noise, spec, params, amplitude,
                                     prepared.noise_phase, delta, path_length);
    const float deformation =
        spec.curl_integrator == CurlIntegrator::EULER_1
            ? path_length
            : sqrtf(delta.re * delta.re + delta.im * delta.im);
    return {output, delta, deformation, path_length};
  }

  HS_FLASH_MEMBER static PlanarWarpStageResult
  warp_polar_chart(const Complex &input, const WarpStageSpec &spec,
                   const WarpStageParams &params) {
    const float radius = sqrtf(input.re * input.re + input.im * input.im);
    const float radial = spec.polar_mode == PolarMode::LOGARITHMIC
                             ? logf(std::max(radius, 1.0f / 4096.0f))
                             : radius;
    const Complex output(params.radial_scale * radial + params.radial_phase,
                         static_cast<float>(spec.polar_harmonic) *
                                 fast_atan2(input.im, input.re) +
                             params.angular_phase);
    return finish_closed_form_warp(input, output);
  }

  /**
   * @brief Pulls plane coordinates back through one warp stage.
   * @param input Coordinates entering the stage.
   * @param projected Projection output, read only for the envelope weight and
   *        edge distance.
   * @param spec Stage kind and its discrete options.
   * @param params Stage parameters, already canonicalized.
   * @param stage_phase Wrapped noise phase for this stage's clock.
   * @param stage_noise Noise resource bound to this stage; may be null for
   *        kinds that sample no noise.
   * @param prepared Per-frame precomputation for this stage.
   * @param frame Frame snapshot.
   * @return Stage output coordinates, the delta it applied, its deformation,
   *         and the path length travelled.
   * @details Path length equals the deformation for the closed-form kinds and
   * the integrated arc length for curl flow.
   */
  HS_FLASH_MEMBER static PlanarWarpStageResult
  warp_stage_lookup(const Complex &input, const ProjectedLookup &projected,
                    const WarpStageSpec &spec, const WarpStageParams &params,
                    float stage_phase, const FastNoiseLite *stage_noise,
                    const PreparedWarpStage &prepared,
                    const FrameState &frame) {
    if (spec.kind == WarpStageKind::NONE)
      return {input, Complex(), 0.0f, 0.0f};
    if (spec.kind == WarpStageKind::LEGACY_STEREO_NOISE) {
      if (params.strength == 0.0f)
        return {input, Complex(), 0.0f, 0.0f};
      const float r_sq = projected.coords.re * projected.coords.re +
                         projected.coords.im * projected.coords.im;
      HS_SB_STAGE_MARK(legacy_start);
      const StereoWarpResult result = sample_wrapped_warp(
          input, r_sq, *stage_noise, params.scale, params.strength,
          frame.params.projection.pole_fade, frame.legacy_noise_phase);
      HS_SB_STAGE_SPAN(legacy_noise, legacy_start);
      return {result.coords, result.delta, result.displacement,
              result.displacement};
    }
    return warp_nonlegacy_stage(input, projected, spec, params, stage_phase,
                                *stage_noise, prepared);
  }

  HS_FLASH_MEMBER static PlanarWarpStageResult
  warp_nonlegacy_stage(const Complex &input, const ProjectedLookup &projected,
                       const WarpStageSpec &spec, const WarpStageParams &params,
                       float stage_phase, const FastNoiseLite &stage_noise,
                       const PreparedWarpStage &prepared) {
    const float envelope =
        warp_envelope(projected, spec.envelope, params.edge_width);
    const float amplitude = params.strength * envelope;
    switch (spec.kind) {
    case WarpStageKind::NONE:
    case WarpStageKind::LEGACY_STEREO_NOISE:
      break;
    case WarpStageKind::AFFINE_FRAME:
      return warp_affine_frame(input, params, prepared);
    case WarpStageKind::WAVE_SHEAR:
      return warp_wave_shear(input, params, stage_phase, amplitude);
    case WarpStageKind::VORTEX:
      return warp_vortex(input, prepared);
    case WarpStageKind::VECTOR_NOISE:
      return warp_vector_noise(input, spec, params, amplitude, stage_noise,
                               prepared);
    case WarpStageKind::CURL_FLOW:
      return warp_curl_flow(input, spec, params, amplitude, stage_noise,
                            prepared);
    case WarpStageKind::MIRROR_TILE: {
      HS_SB_STAGE_MARK(mirror_start);
      const PlanarWarpStageResult result =
          finish_closed_form_warp(input, mirror_tile(input, params, prepared));
      HS_SB_STAGE_SPAN(mirror_tile, mirror_start);
      return result;
    }
    case WarpStageKind::POLAR_CHART:
      return warp_polar_chart(input, spec, params);
    }
    __builtin_unreachable();
  }

  static float warp_envelope(const ProjectedLookup &projected,
                             WarpEnvelope envelope, float edge_width) {
    if (envelope == WarpEnvelope::PROJECTION_WEIGHT)
      return projected.value_weight;
    if (envelope == WarpEnvelope::EDGE_FADE)
      return cubic_kernel(projected.fade_edge_distance / edge_width);
    return 1.0f;
  }

  static float sample_noise_basis(const FastNoiseLite &noise, NoiseBasis basis,
                                  float x, float y, float time) {
    if (basis == NoiseBasis::SIMPLEX)
      return noise.GetNoiseSingle(x, y, time);
    float value = 0.0f;
    float normalization = 0.0f;
    float weight = 1.0f;
    for (int octave = 0; octave < 3; ++octave) {
      const float n = noise.GetNoiseSingle(x, y, time);
      value += weight * (basis == NoiseBasis::RIDGED3 ? 1.0f - fabsf(n) : n);
      normalization += weight;
      x *= 2.0f;
      y *= 2.0f;
      time *= 2.0f;
      weight *= 0.5f;
    }
    return value / normalization;
  }

  static float sample_wrapped_noise_basis(const FastNoiseLite &noise,
                                          NoiseBasis basis, float x, float y,
                                          const WrappedNoisePhase &phase,
                                          float time_offset = 0.0f) {
    const float current = sample_noise_basis(noise, basis, x, y,
                                             phase.current_time + time_offset);
    if (!phase.blends)
      return current;
    const float previous = sample_noise_basis(
        noise, basis, x, y, phase.previous_time + time_offset);
    return hs::lerp(current, previous, phase.mix);
  }

  HS_FLASH_MEMBER static Complex
  curl_flow(const Complex &input, const FastNoiseLite &noise,
            const WarpStageSpec &spec, const WarpStageParams &params,
            float distance, const WrappedNoisePhase &noise_phase,
            Complex &net_delta, float &path_length) {
    if (spec.curl_integrator == CurlIntegrator::EULER_1) {
      const Complex direction =
          curl_vector(input, noise, spec.basis, params.scale, noise_phase);
      net_delta = {distance * direction.re, distance * direction.im};
      path_length =
          sqrtf(net_delta.re * net_delta.re + net_delta.im * net_delta.im);
      return {input.re + net_delta.re, input.im + net_delta.im};
    }
    const int intervals =
        spec.curl_integrator == CurlIntegrator::MIDPOINT_2 ? 2 : 4;
    Complex q = input;
    const float step = distance / intervals;
    path_length = 0.0f;
    net_delta = {};
    for (int index = 0; index < intervals; ++index) {
      const Complex first =
          curl_vector(q, noise, spec.basis, params.scale, noise_phase);
      const Complex midpoint(q.re + 0.5f * step * first.re,
                             q.im + 0.5f * step * first.im);
      const Complex direction =
          curl_vector(midpoint, noise, spec.basis, params.scale, noise_phase);
      const Complex delta(step * direction.re, step * direction.im);
      q = {q.re + delta.re, q.im + delta.im};
      net_delta = {net_delta.re + delta.re, net_delta.im + delta.im};
      path_length += sqrtf(delta.re * delta.re + delta.im * delta.im);
    }
    return q;
  }

  HS_FLASH_MEMBER static Complex
  curl_vector(const Complex &p, const FastNoiseLite &noise, NoiseBasis basis,
              float scale, const WrappedNoisePhase &noise_phase) {
    if (basis == NoiseBasis::SIMPLEX && !noise_phase.blends)
      return simplex_curl_vector(p, noise, scale, noise_phase.current_time);
    constexpr float STENCIL = 1.0f / 64.0f;
    const float x = scale * p.re;
    const float y = scale * p.im;
    const float dx =
        (sample_wrapped_noise_basis(noise, basis, x + STENCIL, y, noise_phase) -
         sample_wrapped_noise_basis(noise, basis, x - STENCIL, y,
                                    noise_phase)) /
        (2.0f * STENCIL);
    const float dy =
        (sample_wrapped_noise_basis(noise, basis, x, y + STENCIL, noise_phase) -
         sample_wrapped_noise_basis(noise, basis, x, y - STENCIL,
                                    noise_phase)) /
        (2.0f * STENCIL);
    return {-scale * dy, scale * dx};
  }

  HS_FLASH_MEMBER static Complex simplex_curl_vector(const Complex &p,
                                                     const FastNoiseLite &noise,
                                                     float scale, float time) {
    constexpr float STENCIL = 1.0f / 64.0f;
    constexpr float ONE_THIRD_STENCIL = STENCIL / 3.0f;
    constexpr float TWO_THIRDS_STENCIL = 2.0f * STENCIL / 3.0f;
    const float x = scale * p.re;
    const float y = scale * p.im;
    const float rotation = (x + y + time) * (2.0f / 3.0f);
    const float tx = rotation - x;
    const float ty = rotation - y;
    const float tz = rotation - time;
    const float x_plus = noise.GetNoiseSingleTransformed(
        tx - ONE_THIRD_STENCIL, ty + TWO_THIRDS_STENCIL,
        tz + TWO_THIRDS_STENCIL);
    const float x_minus = noise.GetNoiseSingleTransformed(
        tx + ONE_THIRD_STENCIL, ty - TWO_THIRDS_STENCIL,
        tz - TWO_THIRDS_STENCIL);
    const float y_plus = noise.GetNoiseSingleTransformed(
        tx + TWO_THIRDS_STENCIL, ty - ONE_THIRD_STENCIL,
        tz + TWO_THIRDS_STENCIL);
    const float y_minus = noise.GetNoiseSingleTransformed(
        tx - TWO_THIRDS_STENCIL, ty + ONE_THIRD_STENCIL,
        tz - TWO_THIRDS_STENCIL);
    const float dx = (x_plus - x_minus) / (2.0f * STENCIL);
    const float dy = (y_plus - y_minus) / (2.0f * STENCIL);
    return {-scale * dy, scale * dx};
  }

  HS_FLASH_MEMBER static Complex
  mirror_tile(const Complex &input, const WarpStageParams &params,
              const PreparedWarpStage &prepared) {
    const float c = prepared.rotation_cos;
    const float s = prepared.rotation_sin;
    const float offset_x = prepared.mirror_offset_x;
    const float offset_y = prepared.mirror_offset_y;
    const float x = c * input.re + s * input.im + offset_x;
    const float y = -s * input.re + c * input.im + offset_y;
    const float folded_x =
        params.cell_x * (1.0f - 2.0f * fabsf(wrap_t(x / params.cell_x) - 0.5f));
    const float folded_y =
        params.cell_y * (1.0f - 2.0f * fabsf(wrap_t(y / params.cell_y) - 0.5f));
    return {c * folded_x - s * folded_y, s * folded_x + c * folded_y};
  }

  static StereoWarpResult sample_wrapped_warp(const Complex &z, float r_sq,
                                              const FastNoiseLite &noise,
                                              float scale, float strength,
                                              float pole_fade,
                                              const WrappedNoisePhase &phase) {
    const StereoWarpResult current = stereo_noise_warp(
        z, r_sq, noise, scale, strength, pole_fade, phase.current_time);
    if (!phase.blends) {
      HS_SB_STAGE_COUNT(++hs::g_shaderball_stage_cycles.legacy_single_samples);
      return current;
    }
    HS_SB_STAGE_COUNT(++hs::g_shaderball_stage_cycles.legacy_blended_samples);
    const StereoWarpResult next = stereo_noise_warp(
        z, r_sq, noise, scale, strength, pole_fade, phase.previous_time);
    return {{hs::lerp(current.coords.re, next.coords.re, phase.mix),
             hs::lerp(current.coords.im, next.coords.im, phase.mix)},
            {hs::lerp(current.delta.re, next.delta.re, phase.mix),
             hs::lerp(current.delta.im, next.delta.im, phase.mix)},
            hs::lerp(current.displacement, next.displacement, phase.mix)};
  }

  static Complex condition_source_coords(const Complex &coords,
                                         const FrameState &frame) {
    if (frame.slots.function == Function::NOISE_CONTOUR ||
        frame.slots.function == Function::PRIMITIVE_LATTICE)
      return coords;
    return stereo_pattern_args(coords, frame.params.source.pattern_freq);
  }

  /**
   * @brief Turns the signed source field into a shaped value and a coverage.
   * @param field Signed source sample, nominally in [-1, 1].
   * @param projected Projection output; supplies the signal weight, edge
   *        distance, and out-of-domain coverage.
   * @param warped Warp output, carried through for the deformation colorizer.
   * @param frame Frame snapshot.
   * @return Value in [0, 1], coverage in [0, 1], and the warp metadata.
   * @details The signal weight scales the signed field before the remap to
   * [0, 1], so it changes value rather than alpha; coverage is the separate
   * alpha channel.
   */
  HS_FLASH_MEMBER static MaterialSample
  shape_nontrivial_material(float value, const ProjectedLookup &projected,
                            const PlanarWarpResult &warped,
                            const FrameState &frame) {
    switch (frame.slots.value_transfer) {
    case ValueTransfer::LINEAR:
      break;
    case ValueTransfer::RIDGE:
      value = 1.0f - fabsf(2.0f * value - 1.0f);
      break;
    case ValueTransfer::ISO_CONTOUR: {
      const float distance = fabsf(value - frame.params.value.iso_level);
      value = 1.0f - smooth_ramp(frame.params.value.iso_width,
                                 2.0f * frame.params.value.iso_width, distance);
      break;
    }
    case ValueTransfer::SMOOTH_BANDS:
      value =
          0.5f -
          0.5f * cosf(TWO_PI_F *
                          static_cast<float>(frame.params.value.band_count) *
                          value +
                      frame.params.value.band_phase);
      break;
    }
    float coverage = 1.0f;
    switch (frame.slots.coverage) {
    case CoveragePolicy::OPAQUE:
      break;
    case CoveragePolicy::PROJECTION_WEIGHT_SQUARED:
      coverage = projected.value_weight * projected.value_weight;
      break;
    case CoveragePolicy::VALUE_CUTOUT:
      coverage = smooth_ramp(frame.params.value.cutout_threshold -
                                 frame.params.value.cutout_softness,
                             frame.params.value.cutout_threshold +
                                 frame.params.value.cutout_softness,
                             value);
      break;
    case CoveragePolicy::EDGE_FADE:
      coverage = frame.params.value.edge_width == 0.0f
                     ? static_cast<float>(projected.fade_edge_distance > 0.0f)
                     : smooth_ramp(0.0f, frame.params.value.edge_width,
                                   projected.fade_edge_distance);
      break;
    case CoveragePolicy::PROJECTION_WEIGHT:
      coverage = projected.value_weight;
      break;
    }
    coverage *= projected.domain_coverage;
    return {value, coverage, warped.net_delta, warped.deformation,
            warped.path_length};
  }

  static MaterialSample shape_material(float field,
                                       const ProjectedLookup &projected,
                                       const PlanarWarpResult &warped,
                                       const FrameState &frame) {
    const float weight = frame.slots.signal_weight == SignalWeight::PROJECTION
                             ? projected.value_weight
                             : 1.0f;
    const float value = hs::clamp((field * weight + 1.0f) * 0.5f, 0.0f, 1.0f);
    if (frame.slots.value_transfer == ValueTransfer::LINEAR &&
        frame.slots.coverage == CoveragePolicy::OPAQUE)
      return {value, projected.domain_coverage, warped.net_delta,
              warped.deformation, warped.path_length};
    return shape_nontrivial_material(value, projected, warped, frame);
  }

  /**
   * @brief Samples the selected source function at a planar coordinate.
   * @param p Conditioned source coordinates.
   * @param frame Frame snapshot.
   * @return Signed field value in [-1, 1].
   */
  HS_FLASH_MEMBER static float sample_noise_contour(const Complex &p,
                                                    const FrameState &frame) {
    const float n = hs::clamp(
        sample_wrapped_noise_basis(
            *frame.resources.source_noise, frame.params.source.noise_basis,
            frame.params.source.noise_scale * p.re,
            frame.params.source.noise_scale * p.im, frame.source_noise_phase),
        -1.0f, 1.0f);
    const float contrast = frame.params.source.noise_contrast;
    return n * (1.0f + contrast) / (1.0f + contrast * fabsf(n));
  }

  HS_FLASH_MEMBER static float sample_source(const Complex &p,
                                             const FrameState &frame) {
    if (frame.slots.function == Function::GRID)
      return grid(p, frame.params.source, frame.prepared_source);
    if (frame.slots.function == Function::NOISE_CONTOUR)
      return sample_noise_contour(p, frame);
    if (frame.slots.function == Function::PRIMITIVE_LATTICE)
      return primitive_lattice(p, frame.params.source);
    return sample_function(frame.slots.function, p, frame.prepared_source);
  }

  HS_FLASH_MEMBER static float primitive_lattice(const Complex &p,
                                                 const SourceParams &params) {
    const float x = wrap_t(params.lattice_cell_scale * p.re + 0.5f) - 0.5f;
    const float y = wrap_t(params.lattice_cell_scale * p.im + 0.5f) - 0.5f;
    const float circle = sqrtf(x * x + y * y) - params.lattice_radius;
    const float bx = fabsf(x) - params.lattice_radius;
    const float by = fabsf(y) - params.lattice_radius;
    const float square = sqrtf(std::max(bx, 0.0f) * std::max(bx, 0.0f) +
                               std::max(by, 0.0f) * std::max(by, 0.0f)) +
                         std::min(std::max(bx, by), 0.0f);
    const float distance = hs::lerp(circle, square, params.lattice_shape_blend);
    return 1.0f - 2.0f * smooth_ramp(-params.lattice_softness,
                                     params.lattice_softness, distance);
  }

  HS_FLASH_MEMBER static float smooth_ramp(float edge0, float edge1,
                                           float value) {
    return cubic_kernel((value - edge0) / (edge1 - edge0));
  }

  /**
   * @brief Maps a shaped material sample to a palette colour.
   * @param sample Shaped value, coverage, and warp metadata.
   * @param frame Frame snapshot.
   * @return Colour whose alpha is the palette alpha scaled by the coverage.
   * @details The deformation colorizer offsets the palette coordinate by the
   * normalized displacement, path length, and delta direction; the liquid
  * colorizer offsets it by the breathe phase.
  */
  HS_FLASH_MEMBER static Color4 colorize_generated(const MaterialSample &sample,
                                                   const FrameState &frame) {
    Color4 color = frame.resources.generated_palette->get(sample.value);
    color.alpha *= sample.coverage;
    return color;
  }

  HS_FLASH_MEMBER static Color4
  colorize_deformation(const MaterialSample &sample, const FrameState &frame) {
    const ColorizerParams &params = frame.params.colorizer;
    const float displacement =
        hs::clamp(sample.deformation / params.displacement_norm, 0.0f, 1.0f);
    const float path =
        hs::clamp(sample.path_length / params.path_norm, 0.0f, 1.0f);
    float direction = 0.0f;
    constexpr float DIRECTION_EPSILON_SQ = 1e-12f;
    if (sample.net_delta.re * sample.net_delta.re +
            sample.net_delta.im * sample.net_delta.im >
        DIRECTION_EPSILON_SQ)
      direction = params.direction_gain *
                  cosf(fast_atan2(sample.net_delta.im, sample.net_delta.re) -
                       params.direction_phase);
    const float u =
        wrap_t(sample.value + params.displacement_gain * displacement +
               params.path_gain * path + direction);
    Color4 color = frame.resources.generated_palette->get(u);
    color.alpha *= sample.coverage;
    return color;
  }

  HS_O3_FN static Color4 colorize(const MaterialSample &sample,
                                  const FrameState &frame) {
    if (frame.slots.colorizer == Colorizer::GENERATED_TRIADIC)
      return colorize_generated(sample, frame);
    if (frame.slots.colorizer == Colorizer::DEFORMATION_INK)
      return colorize_deformation(sample, frame);
    const float value = std::min(sample.value, ONE_BELOW_UNIT);
    const float u = wrap_t(value + frame.breathe_offset);
    Color4 color = frame.resources.liquid_palette->get(u);
    color.alpha *=
        sample.coverage * (1.0f - value * frame.params.colorizer.value_fade);
    if (frame.params.colorizer.hue_shift != 0.0f)
      color = hue_rotate_lut_gamut(color, -sample.deformation *
                                              frame.params.colorizer.hue_shift);
    return color;
  }

  HS_FLASH_MEMBER static OKLab gamut_clip_lut(OKLab lab) {
    lab.L = hs::clamp(lab.L, 0.0f, 1.0f);
    const float chroma_sq = lab.a * lab.a + lab.b * lab.b;
    if (!(chroma_sq > 1e-12f))
      return {lab.L, 0.0f, 0.0f};
    uint32_t inverse_bits;
    std::memcpy(&inverse_bits, &chroma_sq, sizeof(inverse_bits));
    inverse_bits = 0x5f3759dfu - (inverse_bits >> 1);
    float inverse_chroma;
    std::memcpy(&inverse_chroma, &inverse_bits, sizeof(inverse_chroma));
    inverse_chroma *= 1.5f - 0.5f * chroma_sq * inverse_chroma * inverse_chroma;
    const GamutLut &lut = g_gamut_lut;
    int angle_index =
        static_cast<int>(diamond_angle(lab.b, lab.a) * lut.angle_scale);
    int lightness_index = static_cast<int>(lab.L * lut.l_scale);
    angle_index = hs::clamp(angle_index, 0, lut.angle_steps - 1);
    lightness_index = hs::clamp(lightness_index, 0, lut.l_steps - 1);
    const float max_chroma = std::max(
        0.0f,
        static_cast<float>(
            lut.table[(lightness_index * lut.angle_steps + angle_index) * 2]) *
                GAMUT_LUT_INV_SCALE -
            GAMUT_CLIP_MARGIN);
    const float scale = max_chroma * inverse_chroma;
    return {lab.L, lab.a * scale, lab.b * scale};
  }

  static void hue_sincos(float turns, float &cosine, float &sine) {
    const float x = 2.0f * (turns - floorf(turns + 0.5f));
    const auto sine_turn = [](float value) {
      const float parabolic = 4.0f * value * (1.0f - fabsf(value));
      return parabolic * (0.775f + 0.225f * fabsf(parabolic));
    };
    sine = sine_turn(x);
    cosine = sine_turn(0.5f - fabsf(x));
  }

  static Color4 hue_rotate_lut_gamut(const Color4 &color, float amount) {
    float cosine, sine;
    hue_sincos(amount, cosine, sine);
    const LinRGB input = pixel_to_linrgb(color.color);
    OKLab lab = linear_rgb_to_oklab_fast(input.r, input.g, input.b);
    const float rotated_a = lab.a * cosine - lab.b * sine;
    const float rotated_b = lab.a * sine + lab.b * cosine;
    lab = {lab.L, rotated_a, rotated_b};
    LinRGB output = oklab_to_linear_rgb(lab);
    if (!linear_rgb_in_gamut(output.r, output.g, output.b))
      output = oklab_to_linear_rgb(gamut_clip_lut(lab));
    return Color4(Pixel(float_to_pixel16(output.r), float_to_pixel16(output.g),
                        float_to_pixel16(output.b)),
                  color.alpha);
  }

  HS_FLASH_MEMBER static Vector apply_lens(const Vector &v,
                                           const FrameState &frame) {
    switch (frame.slots.surface_lens) {
    case SurfaceLens::NONE:
    case SurfaceLens::GLITCH:
    case SurfaceLens::TWIST:
    case SurfaceLens::KALEIDOSCOPE:
    case SurfaceLens::KALEIDOSCOPE_TETRAHEDRAL:
    case SurfaceLens::KALEIDOSCOPE_OCTAHEDRAL:
    case SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL:
    case SurfaceLens::KALEIDOSCOPE_TRIANGULAR_PRISM:
    case SurfaceLens::KALEIDOSCOPE_SQUARE_PRISM:
    case SurfaceLens::KALEIDOSCOPE_PENTAGONAL_PRISM:
    case SurfaceLens::KALEIDOSCOPE_HEXAGONAL_PRISM:
    case SurfaceLens::KALEIDOSCOPE_OCTAGONAL_PRISM:
      return apply_frame_free_lens(v, frame.slots.surface_lens);
    case SurfaceLens::MOBIUS:
      return mobius_lens(v, frame.params.surface_lens.mobius);
    case SurfaceLens::TANGENT_NOISE:
      return tangent_noise_lens(v, frame);
    }
    __builtin_unreachable();
  }

  __attribute__((always_inline)) static Vector
  profiled_apply_lens(const Vector &v, const FrameState &frame) {
    HS_SB_STAGE_MARK(stage_start);
    const Vector lensed = apply_lens(v, frame);
    HS_SB_STAGE_SPAN(lens, stage_start);
    return lensed;
  }

  /**
   * @brief Applies a lens whose image depends on the direction alone.
   * @param v Unit direction on the sphere.
   * @param lens Lens to apply; MOBIUS and TANGENT_NOISE read FrameState
   *        parameters and are rejected here.
   * @return The lensed direction.
   */
  __attribute__((always_inline)) static Vector
  apply_frame_free_lens(const Vector &v, SurfaceLens lens) {
    switch (lens) {
    case SurfaceLens::NONE:
      return v;
    case SurfaceLens::GLITCH:
      return glitch_lens(v);
    case SurfaceLens::TWIST:
      return twist_lens(v);
    case SurfaceLens::KALEIDOSCOPE:
      return kaleidoscope_lens(v);
    case SurfaceLens::KALEIDOSCOPE_TETRAHEDRAL:
      return polyhedral_kaleidoscope_lens(v, TETRAHEDRAL_MIRRORS);
    case SurfaceLens::KALEIDOSCOPE_OCTAHEDRAL:
      return polyhedral_kaleidoscope_lens(v, OCTAHEDRAL_MIRRORS);
    case SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL:
      return polyhedral_kaleidoscope_lens(v, DODECAHEDRAL_MIRRORS);
    case SurfaceLens::KALEIDOSCOPE_TRIANGULAR_PRISM:
      return polyhedral_kaleidoscope_lens(v, TRIANGULAR_PRISM_MIRRORS);
    case SurfaceLens::KALEIDOSCOPE_SQUARE_PRISM:
      return polyhedral_kaleidoscope_lens(v, SQUARE_PRISM_MIRRORS);
    case SurfaceLens::KALEIDOSCOPE_PENTAGONAL_PRISM:
      return polyhedral_kaleidoscope_lens(v, PENTAGONAL_PRISM_MIRRORS);
    case SurfaceLens::KALEIDOSCOPE_HEXAGONAL_PRISM:
      return polyhedral_kaleidoscope_lens(v, HEXAGONAL_PRISM_MIRRORS);
    case SurfaceLens::KALEIDOSCOPE_OCTAGONAL_PRISM:
      return polyhedral_kaleidoscope_lens(v, OCTAGONAL_PRISM_MIRRORS);
    case SurfaceLens::MOBIUS:
    case SurfaceLens::TANGENT_NOISE:
      HS_CHECK(false, "frame-parameterized lens needs the FrameState overload");
      __builtin_unreachable();
    }
    __builtin_unreachable();
  }

  static Vector nlerp_unit(const Vector &a, const Vector &b, float t) {
    const Vector mixed(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t,
                       a.z + (b.z - a.z) * t);
    const float length_sq =
        mixed.x * mixed.x + mixed.y * mixed.y + mixed.z * mixed.z;
    if (length_sq < 1e-8f)
      return a;
    const float inv_length = 1.0f / sqrtf(length_sq);
    return Vector(mixed.x * inv_length, mixed.y * inv_length,
                  mixed.z * inv_length);
  }

  HS_FLASH_MEMBER static Vector mobius_lens(const Vector &v,
                                            const MobiusParams &params) {
    return mobius_transform(v, params);
  }

  HS_FLASH_MEMBER static Vector tangent_noise_lens(const Vector &v,
                                                   const FrameState &frame) {
    const float scale = frame.params.surface_lens.noise_scale;
    const WrappedNoisePhase &phase = frame.lens_noise_phase;
    const float x = scale * v.x;
    const float y = scale * v.y;
    const float z = scale * v.z;
    const float nx = sample_wrapped_noise_basis(
        *frame.resources.lens_noise, frame.params.surface_lens.noise_basis, x,
        y, phase, z);
    const float ny = sample_wrapped_noise_basis(
        *frame.resources.lens_noise, frame.params.surface_lens.noise_basis,
        x + 17.0f, y - 29.0f, phase, z + 43.0f);
    const float nz = sample_wrapped_noise_basis(
        *frame.resources.lens_noise, frame.params.surface_lens.noise_basis,
        x - 47.0f, y + 11.0f, phase, z - 23.0f);
    const float radial = nx * v.x + ny * v.y + nz * v.z;
    const Vector tangent(nx - radial * v.x, ny - radial * v.y,
                         nz - radial * v.z);
    const float amount = frame.params.surface_lens.amount;
    return nlerp_unit(v,
                      Vector(v.x + amount * tangent.x, v.y + amount * tangent.y,
                             v.z + amount * tangent.z),
                      1.0f);
  }

  static Vector twist_lens(const Vector &v) {
    const float angle = TWIST_RATE * v.y;
    const float c = fast_cosf(angle);
    const float s = fast_sinf(angle);
    return Vector(v.x * c - v.z * s, v.y, v.x * s + v.z * c);
  }

  HS_FLASH_MEMBER static Vector kaleidoscope_lens(const Vector &v) {
    constexpr float COS_60 = 0.5f;
    constexpr float SIN_60 = 0.8660254037844386f;
    constexpr float SQRT_3 = 1.7320508075688772f;
    float x = fabsf(v.x);
    float z = fabsf(v.z);
    if (z > SQRT_3 * x) {
      const float reflected_x = -COS_60 * x + SIN_60 * z;
      z = SIN_60 * x + COS_60 * z;
      x = reflected_x;
    }
    if (SQRT_3 * z > x) {
      const float reflected_x = COS_60 * x + SIN_60 * z;
      z = SIN_60 * x - COS_60 * z;
      x = reflected_x;
    }
    return Vector(x, v.y, z);
  }

  /**
   * @brief Folds a direction into a spherical reflection-group chamber.
   * @param v Unit direction on the sphere.
   * @param mirrors Inward unit normals of the chamber's simple mirrors.
   * @return A symmetry-equivalent direction inside the chamber.
   */
  template <size_t N>
  HS_O3_FN static Vector
  polyhedral_kaleidoscope_lens(Vector v, const std::array<Vector, N> &mirrors) {
    [[maybe_unused]] uint32_t reflections = 0;
    for (int reflection = 0; reflection < POLYHEDRAL_REFLECTION_LIMIT;
         ++reflection) {
      bool inside = true;
      for (const Vector &normal : mirrors) {
        const float distance = dot(v, normal);
        if (distance >= -POLYHEDRAL_MIRROR_EPS)
          continue;
        v.x -= 2.0f * distance * normal.x;
        v.y -= 2.0f * distance * normal.y;
        v.z -= 2.0f * distance * normal.z;
        HS_SB_STAGE_COUNT(++reflections);
        inside = false;
        break;
      }
      if (inside) {
        HS_SB_STAGE_COUNT(++hs::g_shaderball_stage_cycles.polyhedral_pixels);
        HS_SB_STAGE_COUNT(
            hs::g_shaderball_stage_cycles.polyhedral_reflections +=
            reflections);
        HS_SB_STAGE_COUNT(
            hs::g_shaderball_stage_cycles.polyhedral_max_reflections = std::max(
                hs::g_shaderball_stage_cycles.polyhedral_max_reflections,
                reflections));
        return v;
      }
    }
    HS_CHECK(false, "polyhedral kaleidoscope fold did not converge");
    return v;
  }

  /**
   * @brief Projects a sphere direction with the projections that take no
   *        `ProjectionParams`.
   * @param v Direction in the projection frame.
   * @param projection Must be equirectangular, stereographic, or gnomonic; the
   *        cartographic kernels read live parameters and are reached through
   *        `project_branch` instead.
   * @return Plane coordinates in the projection's native units.
   */
  HS_FLASH_MEMBER static Complex project_point(const Vector &v,
                                               Projection projection) {
    switch (projection) {
    case Projection::SINUSOIDAL:
      return folded_sinusoidal(v);
    case Projection::EQUIRECTANGULAR:
      return equirectangular(v);
    case Projection::STEREOGRAPHIC:
      return stereo(v);
    case Projection::GNOMONIC:
      return gnomonic(v);
    case Projection::BONNE:
    case Projection::PEIRCE_QUINCUNCIAL:
    case Projection::AIROCEAN:
      break;
    }
    __builtin_unreachable();
  }

  /**
   * @brief Folded sinusoidal (Sanson-Flamsteed) pseudocylindrical projection.
   * @param v Unit direction on the sphere.
   * @param central_meridian Longitude placed at the image's axis, in radians.
   * @return Plane coordinates in radians: absolute azimuth tapered by
   *         cos(latitude), against latitude.
   * @details Folding the azimuth about the central meridian maps both
   * hemispheres onto one image, so the antimeridian carries no seam; the
   * cos(latitude) taper collapses each pole to a point.
   */
  HS_FLASH_MEMBER static Complex
  folded_sinusoidal(const Vector &v, float central_meridian = 0.0f) {
    const float radius = sqrtf(v.x * v.x + v.z * v.z);
    return {std::fabs(projections::wrap_longitude(fast_atan2(v.z, v.x) -
                                                  central_meridian)) *
                radius,
            0.5f * PI_F - fast_acos(v.y)};
  }

  /**
   * @brief Equirectangular (plate carree) cylindrical projection.
   * @param v Unit direction on the sphere.
   * @param central_meridian Longitude placed at the image's axis, in radians.
   * @return Plane coordinates in radians: wrapped longitude against latitude.
   * @details Longitude is periodic, so the image is cut at the antimeridian and
   * each pole spreads across a full image row.
   */
  HS_FLASH_MEMBER static Complex
  equirectangular(const Vector &v, float central_meridian = 0.0f) {
    return {
        projections::wrap_longitude(fast_atan2(v.z, v.x) - central_meridian),
        0.5f * PI_F - fast_acos(v.y)};
  }

  HS_FLASH_MEMBER static Complex gnomonic(const Vector &v) {
    float y = v.y;
    if (std::fabs(y) < GNOMONIC_AXIS_EPS)
      y = y < 0.0f ? -GNOMONIC_AXIS_EPS : GNOMONIC_AXIS_EPS;
    return {v.x / y, v.z / y};
  }

  HS_FLASH_MEMBER static float sample_function(Function function,
                                               const Complex &p,
                                               const SourceState &source) {
    switch (function) {
    case Function::TWIN_WAVE:
      return twin_wave(p, source);
    case Function::RINGS:
      return rings(p, source);
    case Function::SPIRAL:
      return spiral(p, source);
    case Function::GRID:
    case Function::NOISE_CONTOUR:
    case Function::PRIMITIVE_LATTICE:
      break;
    }
    __builtin_unreachable();
  }

  HS_FLASH_MEMBER static float twin_wave(const Complex &p,
                                         const SourceState &source) {
    const float rotated = p.re * source.angle_cos + p.im * source.angle_sin;
    return 0.5f * (fast_sinf(p.re + source.primary) +
                   fast_sinf(rotated + source.primary));
  }

  HS_FLASH_MEMBER static float rings(const Complex &p,
                                     const SourceState &source) {
    return fast_sinf(sqrtf(p.re * p.re + p.im * p.im) - source.primary);
  }

  HS_FLASH_MEMBER static float spiral(const Complex &p,
                                      const SourceState &source) {
    const float radius = sqrtf(p.re * p.re + p.im * p.im);
    const float azimuth = fast_atan2(p.im, p.re);
    return fast_sinf(radius - SPIRAL_ARMS * (azimuth + source.angle) -
                     source.primary);
  }

  HS_O3_FN static float grid(const Complex &p, const SourceParams &params,
                             const SourceState &source) {
    const float x = p.re * source.angle_cos + p.im * source.angle_sin;
    const float y = -p.re * source.angle_sin + p.im * source.angle_cos;
    if (params.pattern_mix == 1.0f)
      return fast_sinf(x + source.primary) * fast_cosf(y - source.secondary);
    float re = x;
    float im = y;
    if (params.complexity != 0.0f) {
      re += params.complexity * fast_sinf(y + source.primary);
      im += params.complexity * fast_cosf(x - source.secondary);
    }
    const float coupled = fast_sinf(re) * fast_cosf(im);
    if (params.pattern_mix == 0.0f)
      return coupled;
    const float direct =
        fast_sinf(x + source.primary) * fast_cosf(y - source.secondary);
    return hs::lerp(coupled, direct, params.pattern_mix);
  }

  HS_COLD_MEMBER WalkDeltas sample_walk_deltas() {
#if HS_ENABLE_TEST_HOOKS
    ++walk_step_count;
#endif
    const Quaternion projection = projection_walk.get();
    const Quaternion projection_delta =
        projection * projection_walk_prev.conjugate();
    projection_walk_prev = projection;
    const Quaternion outer = outer_walk.get();
    const Quaternion outer_delta = outer * outer_walk_prev.conjugate();
    outer_walk_prev = outer;
    return {projection_delta.normalized(), outer_delta.normalized()};
  }

  HS_COLD_MEMBER void update_spatial_frames(LookRuntime &look,
                                            const Config &config,
                                            const WalkDeltas &deltas) const {
    look.projection_wander = (slerp(Quaternion(), deltas.projection,
                                    config.params.projection.wander) *
                              look.projection_wander)
                                 .normalized();
    look.outer_wander =
        (slerp(Quaternion(), deltas.outer, config.params.outer_camera.wander) *
         look.outer_wander)
            .normalized();
    look.transforms.projection_conj =
        (make_rotation(Y_AXIS, look.clocks.projection_spin) * base_orientation *
         look.projection_wander)
            .conjugate();
    look.transforms.outer_conj = look.outer_wander.conjugate();
  }

  HS_COLD_MEMBER void advance_runtime(LookRuntime &look, const Config &config,
                                      const WalkDeltas &deltas) const {
    const Params &params = config.params;
    const bool legacy_inner = config.slots.warp_program.inner.kind ==
                              WarpStageKind::LEGACY_STEREO_NOISE;
    const float warp_time_scale = legacy_inner ? params.warp.inner.time_scale
                                               : params.warp.outer.time_scale;
    look.clocks.warp_time =
        fmodf(look.clocks.warp_time + params.source.speed * warp_time_scale,
              STEREO_NOISE_TIME_PERIOD);
    look.clocks.source_primary =
        fmodf(look.clocks.source_primary + params.source.speed, TWO_PI_F);
    look.clocks.source_secondary =
        fmodf(look.clocks.source_secondary +
                  params.source.speed * params.source.secondary_rate,
              TWO_PI_F);
    look.clocks.source_angle =
        fmodf(look.clocks.source_angle + params.source.angle_rate, TWO_PI_F);
    look.clocks.projection_spin = fmodf(
        look.clocks.projection_spin + params.projection.spin_rate, TWO_PI_F);
    look.clocks.breathe_phase = fmodf(
        look.clocks.breathe_phase + params.colorizer.cycle_speed, TWO_PI_F);
    look.clocks.source_noise_time =
        wrap_t(look.clocks.source_noise_time + params.source.noise_time_rate);
    look.clocks.lens_noise_time =
        wrap_t(look.clocks.lens_noise_time + params.surface_lens.noise_rate);
    look.clocks.warp_outer_phase =
        wrap_t(look.clocks.warp_outer_phase + params.warp.outer.time_scale);
    look.clocks.warp_inner_phase =
        wrap_t(look.clocks.warp_inner_phase + params.warp.inner.time_scale);
    update_spatial_frames(look, config, deltas);
  }

  HS_COLD_MEMBER void prepare_param_morph() {
    if (!state->param_morph.active)
      return;
    const float mix =
        transition_mix(state->param_morph.elapsed, state->param_morph.duration);
    if (mix == 0.0f)
      blend.params = state->param_morph.from;
    else if (mix == 1.0f)
      blend.params = state->param_morph.to;
    else if (state->param_morph.staggered)
      blend.params.lerp_staggered(state->param_morph.from,
                                  state->param_morph.to, mix);
    else
      blend.params.lerp(state->param_morph.from, state->param_morph.to, mix);
  }

  static float transition_mix(uint16_t elapsed, uint16_t duration) {
    if (elapsed == 0)
      return 0.0f;
    if (elapsed >= duration)
      return 1.0f;
    return ease_in_out_sin(static_cast<float>(elapsed) / duration);
  }

  HS_COLD_MEMBER void apply_requested_config() {
#if HS_ENABLE_PARAM_GUI_BRIDGE
    if (!requested_schema_bound) {
      if (!valid_config(requested_config)) {
        reject_requested_config();
        return;
      }
      accepted_config = requested_config;
    } else {
      const size_t before_count = pending_edit_count;
      refresh_accepted_config();
      if (before_count != pending_edit_count)
        rebind_parameters();
    }
    const Config &next_config = accepted_config;
#else
    const Config &next_config = requested_config;
    if (!valid_config(next_config)) {
      reject_requested_config();
      return;
    }
#endif
    if (next_config == published_config)
      return;
    if (!prepare_resource_union(next_config, next_config)) {
      reject_requested_config();
      return;
    }
    if (state->transition.active)
      runtime = state->transition.elapsed * 2 < state->transition.duration
                    ? state->transition.from_runtime
                    : state->transition.to_runtime;
    state->transition.active = false;
    state->param_morph.active = false;
    active_slots = next_config.slots;
    blend.params = next_config.params;
#if HS_ENABLE_PARAM_GUI_BRIDGE
    display_config = next_config;
#endif
    published_config = next_config;
    accepted_config = next_config;
    if (!requested_schema_bound)
      rebind_parameters();
  }

  HS_COLD_MEMBER void reject_requested_config() {
    requested_config = published_config;
    accepted_config = published_config;
#if HS_ENABLE_PARAM_GUI_BRIDGE
    pending_edit_count = 0;
    display_config = published_config;
#endif
    rebind_parameters();
  }

  static void canonicalize_mobius(MobiusParams &params) {
    Complex *coefficients[] = {&params.a, &params.b, &params.c, &params.d};
    Complex pivot = params.a;
    float pivot_magnitude_sq = pivot.re * pivot.re + pivot.im * pivot.im;
    float norm_sq = pivot_magnitude_sq;
    for (size_t index = 1; index < std::size(coefficients); ++index) {
      const Complex candidate = *coefficients[index];
      const float magnitude_sq =
          candidate.re * candidate.re + candidate.im * candidate.im;
      norm_sq += magnitude_sq;
      if (magnitude_sq > pivot_magnitude_sq) {
        pivot = candidate;
        pivot_magnitude_sq = magnitude_sq;
      }
    }
    if (norm_sq == 0.0f || pivot_magnitude_sq == 0.0f)
      return;
    const float factor = 1.0f / sqrtf(norm_sq * pivot_magnitude_sq);
    const Complex inverse(pivot.re * factor, -pivot.im * factor);
    for (Complex *coefficient : coefficients)
      *coefficient = *coefficient * inverse;
  }

  HS_COLD_MEMBER bool try_apply_config(const Config &candidate,
                                       uint16_t duration, bool staggered,
                                       bool continue_choreo) {
    if (!valid_config(candidate) || duration == 0)
      return false;
    if (state->transition.active)
      return false;
    const Config current{active_slots, blend.params};
    if (!transition_admitted(current, candidate))
      return false;
    if (candidate == current) {
      state->param_morph.active = false;
      return true;
    }
    if (stable_topology(current, candidate)) {
      state->param_morph = {current.params, candidate.params, 0,   duration,
                            staggered,      continue_choreo,  true};
      return true;
    }
    if (!prepare_resource_union(current, current))
      return false;
    const uint16_t planned_duration =
        (duration & 1U) != 0 ? duration + 1 : duration;
    state->param_morph.active = false;
    state->transition = {current, candidate,        runtime,         runtime,
                         0,       planned_duration, continue_choreo, true};
    return true;
  }

  HS_COLD_MEMBER void finish_transitions() {
    if (state->transition.active) {
      if (state->transition.elapsed == state->transition.duration / 2)
        HS_CHECK(prepare_resource_union(state->transition.to_config,
                                        state->transition.to_config),
                 "through-clear destination resources exceed capacity");
      if (state->transition.elapsed < state->transition.duration) {
        ++state->transition.elapsed;
        return;
      }
      const bool continue_choreo = state->transition.continue_choreo;
      runtime = state->transition.to_runtime;
      active_slots = state->transition.to_config.slots;
      blend.params = state->transition.to_config.params;
      state->transition.active = false;
      if (continue_choreo)
        enter_preset();
      return;
    }
    if (!state->param_morph.active)
      return;
    if (state->param_morph.elapsed < state->param_morph.duration) {
      ++state->param_morph.elapsed;
      return;
    }
    const bool continue_choreo = state->param_morph.continue_choreo;
    blend.params = state->param_morph.to;
    state->param_morph.active = false;
    if (continue_choreo)
      enter_preset();
  }

  HS_COLD_MEMBER void publish_live_config() {
    if (anims_paused || state->transition.active || state->param_morph.active ||
        accepted_config != published_config)
      return;
    published_config = {active_slots, blend.params};
#if HS_ENABLE_PARAM_GUI_BRIDGE
    Config next_requested = published_config;
    for (size_t index = 0; index < pending_edit_count; ++index)
      copy_pending_value(next_requested, requested_config,
                         pending_edits[index]);
    requested_config = next_requested;
#else
    requested_config = published_config;
#endif
    accepted_config = published_config;
  }

#if HS_ENABLE_PARAM_GUI_BRIDGE
  HS_COLD_MEMBER void refresh_parameter_display() override {
    if (state->transition.active) {
      const float mix =
          transition_mix(state->transition.elapsed, state->transition.duration);
      display_config.slots = mix < 0.5f ? state->transition.from_config.slots
                                        : state->transition.to_config.slots;
      display_config.params.lerp(state->transition.from_config.params,
                                 state->transition.to_config.params, mix);
      return;
    }
    display_config = {active_slots, blend.params};
  }
#endif

  template <typename Enum>
  HS_COLD_MEMBER static constexpr bool enum_at_most(Enum value, Enum last) {
    return static_cast<uint8_t>(value) <= static_cast<uint8_t>(last);
  }

  HS_COLD_MEMBER static constexpr SourceTraits
  source_traits(Function function) {
    switch (function) {
    case Function::GRID:
      return {true, true, true, true, false, true};
    case Function::PRIMITIVE_LATTICE:
      return {true, true, true, true, false, false};
    case Function::TWIN_WAVE:
    case Function::RINGS:
    case Function::SPIRAL:
    case Function::NOISE_CONTOUR:
      return {true, true, false, false, false, false};
    }
    return {true, true, false, false, false, false};
  }

  HS_COLD_MEMBER static constexpr bool
  polar_source_compatible(const RequestedConfig &config,
                          const WarpStageSpec &polar) {
    const SourceTraits traits = source_traits(config.slots.function);
    if (!traits.y_periodic || !traits.polar_angle_compatible)
      return false;
    const float periods = config.params.source.pattern_freq *
                          static_cast<float>(polar.polar_harmonic);
    return periods == static_cast<float>(static_cast<int>(periods));
  }

  HS_COLD_MEMBER static constexpr bool
  strict_seam_compatible(const Config &config) {
    if (!strict_projection(config.slots.projection))
      return true;
    return config.slots.function != Function::NOISE_CONTOUR &&
           !seam_sensitive_warp(config.slots.warp_program.outer.kind) &&
           !seam_sensitive_warp(config.slots.warp_program.inner.kind);
  }

  HS_COLD_MEMBER static constexpr bool
  valid_config(const RequestedConfig &candidate) {
    const Slots &slots = candidate.slots;
    if (!enum_at_most(slots.function, Function::PRIMITIVE_LATTICE) ||
        !enum_at_most(slots.projection, Projection::EQUIRECTANGULAR) ||
        !enum_at_most(slots.peirce_layout, PeirceLayout::VERTICAL) ||
        !enum_at_most(slots.airocean_layout, AiroceanLayout::HORIZONTAL) ||
        !enum_at_most(slots.bonne_hemisphere, BonneHemisphere::SOUTH) ||
        !enum_at_most(slots.gnomonic_hemisphere,
                      GnomonicHemispherePolicy::BACK_HEMISPHERE) ||
        !enum_at_most(slots.projection_frame,
                      ProjectionFramePolicy::SPIN_WANDER) ||
        !enum_at_most(slots.surface_lens,
                      SurfaceLens::KALEIDOSCOPE_OCTAGONAL_PRISM) ||
        !enum_at_most(slots.warp_program.outer.kind,
                      WarpStageKind::POLAR_CHART) ||
        !enum_at_most(slots.warp_program.inner.kind,
                      WarpStageKind::POLAR_CHART) ||
        !valid_warp_spec(slots.warp_program.outer) ||
        !valid_warp_spec(slots.warp_program.inner) ||
        !enum_at_most(slots.signal_weight, SignalWeight::PROJECTION) ||
        !enum_at_most(slots.value_transfer, ValueTransfer::SMOOTH_BANDS) ||
        !enum_at_most(slots.coverage, CoveragePolicy::PROJECTION_WEIGHT) ||
        !enum_at_most(slots.colorizer, Colorizer::DEFORMATION_INK))
      return false;
    const int legacy_stages =
        (slots.warp_program.outer.kind == WarpStageKind::LEGACY_STEREO_NOISE) +
        (slots.warp_program.inner.kind == WarpStageKind::LEGACY_STEREO_NOISE);
    if (legacy_stages > 1)
      return false;
    if (legacy_stages == 1 && slots.projection != Projection::STEREOGRAPHIC)
      return false;
    const bool outer_polar =
        slots.warp_program.outer.kind == WarpStageKind::POLAR_CHART;
    const bool inner_polar =
        slots.warp_program.inner.kind == WarpStageKind::POLAR_CHART;
    if (inner_polar &&
        !polar_source_compatible(candidate, slots.warp_program.inner))
      return false;
    if (outer_polar &&
        (slots.warp_program.inner.kind != WarpStageKind::NONE ||
         !polar_source_compatible(candidate, slots.warp_program.outer)))
      return false;
    if (!strict_seam_compatible(candidate) ||
        !preset_in_ranges(candidate.params))
      return false;
    if (!valid_stage_tuple(slots.warp_program.outer,
                           candidate.params.warp.outer) ||
        !valid_stage_tuple(slots.warp_program.inner,
                           candidate.params.warp.inner) ||
        !safe_program_bounds(candidate))
      return false;
    if (slots.surface_lens == SurfaceLens::MOBIUS &&
        !valid_mobius(candidate.params.surface_lens.mobius))
      return false;
    return resource_union_fits(candidate, candidate);
  }

#if HS_ENABLE_PARAM_GUI_BRIDGE
  const char *begin_warning(const char *format, ...) const {
    va_list args;
    va_start(args, format);
    std::vsnprintf(warning_text.data(), warning_text.size(), format, args);
    va_end(args);
    return warning_text.data();
  }

  void append_warning(const char *format, ...) const {
    const size_t length = std::strlen(warning_text.data());
    if (length >= warning_text.size() - 1)
      return;
    va_list args;
    va_start(args, format);
    std::vsnprintf(warning_text.data() + length, warning_text.size() - length,
                   format, args);
    va_end(args);
  }

  bool append_range_warning(const char *label, float value, float minimum,
                            float maximum) const {
    if (value >= minimum && value <= maximum)
      return false;
    append_warning(" %s %.7g is outside [%.7g, %.7g].", label,
                   static_cast<double>(value), static_cast<double>(minimum),
                   static_cast<double>(maximum));
    return true;
  }

  const char *stage_tuple_warning(const char *position,
                                  const WarpStageSpec &spec,
                                  const WarpStageParams &params) const {
    begin_warning("%s %s rejected.", position,
                  WARP_OPTIONS[static_cast<uint8_t>(spec.kind)]);
    switch (spec.kind) {
    case WarpStageKind::NONE:
      break;
    case WarpStageKind::LEGACY_STEREO_NOISE:
      append_range_warning("Warp Scale", params.scale, 0.1f, 100.0f);
      append_range_warning("Warp Strength", params.strength, 0.0f, 30.0f);
      append_range_warning("Warp Time", params.time_scale, 0.05f, 1.0f);
      break;
    case WarpStageKind::AFFINE_FRAME:
      append_range_warning("Translate X", params.translation_x, -4.0f, 4.0f);
      append_range_warning("Translate Y", params.translation_y, -4.0f, 4.0f);
      append_range_warning("Scale X", params.scale_x, 0.25f, 4.0f);
      append_range_warning("Scale Y", params.scale_y, 0.25f, 4.0f);
      append_range_warning("Shear", params.shear, -0.75f, 0.75f);
      break;
    case WarpStageKind::WAVE_SHEAR:
      append_range_warning("Warp Strength", params.strength, -4.0f, 4.0f);
      append_range_warning("Frequency", params.frequency, 0.0f, 64.0f);
      append_range_warning("Warp Time", params.time_scale, NOISE_RATE_MIN,
                           NOISE_RATE_MAX);
      break;
    case WarpStageKind::VORTEX:
      append_range_warning("Radius", params.radius, 1.0f / 64.0f, 8.0f);
      append_range_warning("Turns", params.turns, -4.0f, 4.0f);
      append_range_warning("Orbit Radius", params.center_orbit_radius, 0.0f,
                           4.0f);
      append_range_warning("Warp Time", params.time_scale, NOISE_RATE_MIN,
                           NOISE_RATE_MAX);
      break;
    case WarpStageKind::VECTOR_NOISE:
      append_range_warning("Warp Strength", params.strength, 0.0f, 4.0f);
      append_range_warning("Warp Scale", params.scale, 1.0f / 64.0f, 64.0f);
      append_range_warning("Warp Time", params.time_scale, NOISE_RATE_MIN,
                           NOISE_RATE_MAX);
      break;
    case WarpStageKind::CURL_FLOW: {
      append_range_warning("Warp Strength", params.strength, -4.0f, 4.0f);
      append_range_warning("Warp Scale", params.scale, 1.0f / 64.0f, 16.0f);
      append_range_warning("Warp Time", params.time_scale, NOISE_RATE_MIN,
                           NOISE_RATE_MAX);
      const float strength_limit = curl_strength_limit(spec, params);
      if (abs_value(params.strength) > strength_limit)
        append_warning(
            " %s at Warp Scale %.7g requires |Warp Strength| <= %.9f; "
            "current value is %.7g.",
            CURL_INTEGRATOR_OPTIONS[static_cast<uint8_t>(spec.curl_integrator)],
            static_cast<double>(params.scale),
            static_cast<double>(strength_limit),
            static_cast<double>(params.strength));
      break;
    }
    case WarpStageKind::MIRROR_TILE:
      append_range_warning("Cell X", params.cell_x, CELL_MIN, CELL_MAX);
      append_range_warning("Cell Y", params.cell_y, CELL_MIN, CELL_MAX);
      break;
    case WarpStageKind::POLAR_CHART:
      append_range_warning("Radial Scale", params.radial_scale, 1.0f / 64.0f,
                           16.0f);
      break;
    }
    append_warning(" Set every listed control within its stated limit.");
    return warning_text.data();
  }

  const char *program_bounds_warning(const Config &candidate) const {
    float bound = projection_coordinate_bound(candidate);
    const WarpStageSpec stages[] = {candidate.slots.warp_program.outer,
                                    candidate.slots.warp_program.inner};
    const WarpStageParams params[] = {candidate.params.warp.outer,
                                      candidate.params.warp.inner};
    const char *positions[] = {"Planar Warp 1", "Planar Warp 2"};
    for (size_t index = 0; index < 2; ++index) {
      if (stages[index].kind == WarpStageKind::VECTOR_NOISE ||
          stages[index].kind == WarpStageKind::CURL_FLOW) {
        const float lattice_bound = params[index].scale * (bound + 100.0f);
        if (lattice_bound > NOISE_LATTICE_LIMIT) {
          const float scale_limit = NOISE_LATTICE_LIMIT / (bound + 100.0f);
          return begin_warning(
              "%s %s rejected: Warp Scale %.7g produces noise coordinate "
              "bound %.7g above %.7g. Set Warp Scale <= %.7g or choose a "
              "projection/lens with a smaller coordinate extent.",
              positions[index],
              WARP_OPTIONS[static_cast<uint8_t>(stages[index].kind)],
              static_cast<double>(params[index].scale),
              static_cast<double>(lattice_bound),
              static_cast<double>(NOISE_LATTICE_LIMIT),
              static_cast<double>(scale_limit));
        }
      }
      bound = stage_coordinate_bound(stages[index], params[index], bound);
      if (bound > WARP_COORD_LIMIT)
        return begin_warning(
            "%s %s rejected: its predicted coordinate bound %.7g exceeds "
            "%.7g. Reduce this warp's displacement/translation controls or "
            "choose a projection/lens with a smaller coordinate extent.",
            positions[index],
            WARP_OPTIONS[static_cast<uint8_t>(stages[index].kind)],
            static_cast<double>(bound), static_cast<double>(WARP_COORD_LIMIT));
    }
    const float source_bound = candidate.params.source.noise_scale * bound;
    return begin_warning(
        "Noise Contour rejected: Source Noise Scale %.7g produces noise "
        "coordinate bound %.7g above %.7g. Set Source Noise Scale <= %.7g "
        "or reduce the preceding warp extent.",
        static_cast<double>(candidate.params.source.noise_scale),
        static_cast<double>(source_bound),
        static_cast<double>(NOISE_LATTICE_LIMIT),
        static_cast<double>(NOISE_LATTICE_LIMIT / bound));
  }

  const char *resource_warning(const Config &candidate) const {
    std::array<NoiseResourceKey, 4> keys{};
    const char *owners[4]{};
    size_t count = 0;
    auto add = [&](const char *owner, const NoiseResourceKey &key) {
      for (size_t index = 0; index < count; ++index) {
        if (keys[index].resource_id != key.resource_id)
          continue;
        if (keys[index] == key)
          return false;
        begin_warning(
            "%s conflicts with %s: both use noise resource %u with different "
            "basis, seed, or sampling requirements. Disable one noise "
            "consumer or choose a non-noise Function, Lens, or Warp for one "
            "of them.",
            owner, owners[index], static_cast<unsigned>(key.resource_id));
        return true;
      }
      keys[count] = key;
      owners[count++] = owner;
      return false;
    };
    if (warp_uses_noise(candidate.slots.warp_program.outer.kind) &&
        add("Planar Warp 1",
            warp_resource_key(candidate.slots.warp_program.outer)))
      return warning_text.data();
    if (warp_uses_noise(candidate.slots.warp_program.inner.kind) &&
        add("Planar Warp 2",
            warp_resource_key(candidate.slots.warp_program.inner)))
      return warning_text.data();
    if (candidate.slots.function == Function::NOISE_CONTOUR &&
        add("Function Noise Contour", source_resource_key(candidate)))
      return warning_text.data();
    if (candidate.slots.surface_lens == SurfaceLens::TANGENT_NOISE &&
        add("Lens Tangent Noise", lens_resource_key(candidate)))
      return warning_text.data();
    return begin_warning(
        "The active noise consumers exceed the resource limit of %u. Disable "
        "one noise Function, Lens, or Warp.",
        static_cast<unsigned>(MAX_NOISE_RESOURCES));
  }

  const char *admission_warning(const Config &candidate,
                                const char *edited_name) const {
    const WarpStageSpec &outer = candidate.slots.warp_program.outer;
    const WarpStageSpec &inner = candidate.slots.warp_program.inner;
    const int legacy_stages =
        (outer.kind == WarpStageKind::LEGACY_STEREO_NOISE) +
        (inner.kind == WarpStageKind::LEGACY_STEREO_NOISE);
    if (legacy_stages > 1)
      return begin_warning(
          "Planar Warp 1 and Planar Warp 2 are both Stereo Noise, but only one "
          "Stereo Noise stage is allowed. Set either Planar Warp 1 or Planar "
          "Warp 2 "
          "to None or another warp.");
    if (legacy_stages == 1 &&
        candidate.slots.projection != Projection::STEREOGRAPHIC) {
      const char *position = outer.kind == WarpStageKind::LEGACY_STEREO_NOISE
                                 ? "Planar Warp 1"
                                 : "Planar Warp 2";
      return begin_warning(
          "%s Stereo Noise requires Projection = Stereographic; current "
          "Projection is %s. Select Stereographic or choose a different %s.",
          position,
          PROJECTION_OPTIONS[static_cast<uint8_t>(candidate.slots.projection)],
          position);
    }
    if (outer.kind == WarpStageKind::POLAR_CHART &&
        inner.kind != WarpStageKind::NONE)
      return begin_warning(
          "Planar Warp 1 Polar Chart cannot run while Planar Warp 2 is %s. Set "
          "Planar Warp 2 to None or choose a different Planar Warp 1.",
          WARP_OPTIONS[static_cast<uint8_t>(inner.kind)]);
    const WarpStageSpec *polar =
        outer.kind == WarpStageKind::POLAR_CHART   ? &outer
        : inner.kind == WarpStageKind::POLAR_CHART ? &inner
                                                   : nullptr;
    if (polar != nullptr && !polar_source_compatible(candidate, *polar)) {
      const char *position =
          polar == &outer ? "Planar Warp 1" : "Planar Warp 2";
      const SourceTraits traits = source_traits(candidate.slots.function);
      if (!traits.y_periodic || !traits.polar_angle_compatible)
        return begin_warning(
            "%s Polar Chart requires a polar-periodic Function; %s is not "
            "compatible. Select Grid or Primitive Lattice, or "
            "choose another %s.",
            position,
            FUNCTION_OPTIONS[static_cast<uint8_t>(candidate.slots.function)],
            position);
      const float periods = candidate.params.source.pattern_freq *
                            static_cast<float>(polar->polar_harmonic);
      const float nearest_periods = floorf(periods + 0.5f);
      const float suggested_frequency =
          nearest_periods / static_cast<float>(polar->polar_harmonic);
      return begin_warning(
          "%s Polar Chart requires Pattern Freq x Polar Harmonic to be a "
          "whole number. %.7g x %u = %.7g. Set Pattern Freq to %.7g or change "
          "%s Polar Harmonic.",
          position, static_cast<double>(candidate.params.source.pattern_freq),
          static_cast<unsigned>(polar->polar_harmonic),
          static_cast<double>(periods),
          static_cast<double>(suggested_frequency), position);
    }
    if (!strict_seam_compatible(candidate)) {
      begin_warning(
          "Projection %s requires seam-safe stages.",
          PROJECTION_OPTIONS[static_cast<uint8_t>(candidate.slots.projection)]);
      if (candidate.slots.function == Function::NOISE_CONTOUR)
        append_warning(" Function Noise Contour is not seam-safe.");
      if (seam_sensitive_warp(outer.kind))
        append_warning(" Planar Warp 1 %s is not seam-safe.",
                       WARP_OPTIONS[static_cast<uint8_t>(outer.kind)]);
      if (seam_sensitive_warp(inner.kind))
        append_warning(" Planar Warp 2 %s is not seam-safe.",
                       WARP_OPTIONS[static_cast<uint8_t>(inner.kind)]);
      append_warning(" Replace the named stage or select Folded Sinusoidal, "
                     "Stereographic, Gnomonic, or Equirectangular.");
      return warning_text.data();
    }
    if (!preset_in_ranges(candidate.params)) {
      const ParamDef *parameter = getParameters().find(edited_name);
      if (parameter != nullptr)
        return begin_warning(
            "%s %.7g is outside its registered range [%.7g, %.7g]. Set %s "
            "within that range.",
            edited_name, static_cast<double>(parameter->get_requested()),
            static_cast<double>(parameter->min),
            static_cast<double>(parameter->max), edited_name);
    }
    if (!valid_stage_tuple(outer, candidate.params.warp.outer))
      return stage_tuple_warning("Planar Warp 1", outer,
                                 candidate.params.warp.outer);
    if (!valid_stage_tuple(inner, candidate.params.warp.inner))
      return stage_tuple_warning("Planar Warp 2", inner,
                                 candidate.params.warp.inner);
    if (!safe_program_bounds(candidate))
      return program_bounds_warning(candidate);
    if (candidate.slots.surface_lens == SurfaceLens::MOBIUS &&
        !valid_mobius(candidate.params.surface_lens.mobius)) {
      const MobiusParams &m = candidate.params.surface_lens.mobius;
      const float det_re =
          m.a.re * m.d.re - m.a.im * m.d.im - m.b.re * m.c.re + m.b.im * m.c.im;
      const float det_im =
          m.a.re * m.d.im + m.a.im * m.d.re - m.b.re * m.c.im - m.b.im * m.c.re;
      return begin_warning(
          "Mobius Lens rejected: |A*D - B*C| is %.7g; it must be at least "
          "0.001. Adjust the requested Mobius coefficient until the "
          "determinant reaches 0.001 or more.",
          static_cast<double>(sqrtf(det_re * det_re + det_im * det_im)));
    }
    if (!resource_union_fits(candidate, candidate))
      return resource_warning(candidate);
    return begin_warning(
        "%s was rejected by an unclassified ShaderBall admission rule. Keep "
        "the requested value and report this exact configuration as a bug.",
        edited_name);
  }
#endif

  HS_COLD_MEMBER static constexpr bool
  valid_warp_spec(const WarpStageSpec &spec) {
    return enum_at_most(spec.basis, NoiseBasis::RIDGED3) &&
           enum_at_most(spec.envelope, WarpEnvelope::EDGE_FADE) &&
           enum_at_most(spec.polar_mode, PolarMode::LOGARITHMIC) &&
           enum_at_most(spec.curl_integrator, CurlIntegrator::MIDPOINT_4) &&
           spec.polar_harmonic >= 1 &&
           spec.polar_harmonic <= POLAR_HARMONIC_MAX;
  }

  HS_COLD_MEMBER static constexpr float abs_value(float value) {
    return value < 0.0f ? -value : value;
  }

  HS_COLD_MEMBER static constexpr int
  curl_intervals(CurlIntegrator integrator) {
    return integrator == CurlIntegrator::EULER_1      ? 1
           : integrator == CurlIntegrator::MIDPOINT_2 ? 2
                                                      : 4;
  }

  /**
   * @brief Conservative bound on the sampled gradient magnitude of a basis.
   * @return 64, per unit of scaled noise coordinate.
   * @details The bound is a property of the stencil, not of the basis, so the
   * parameter is unnamed. Every basis returns values in [-1, 1] and the curl
   * gradient is a central difference over the fixed stencil h = 1/64, so no
   * component can exceed 2 / (2h) = 1/h = 64. Widening the stencil or letting
   * a basis leave [-1, 1] invalidates this and every curl stability check
   * built on it.
   */
  HS_COLD_MEMBER static constexpr float noise_gradient_bound(NoiseBasis) {
    return 64.0f;
  }

  /**
   * @brief Largest curl-flow strength the stage stability inequality admits at
   *        a stage's live scale, basis, and integrator.
   * @details Solves `scale * |strength| * G / n <= 1/2` — the same inequality
   * `valid_stage_tuple` enforces — for `|strength|`, so the registered slider
   * spans exactly the admissible set instead of a range whose bulk is rejected.
   */
  HS_COLD_MEMBER static constexpr float
  curl_strength_limit(const WarpStageSpec &spec,
                      const WarpStageParams &params) {
    const float scale =
        params.scale > WARP_SCALE_MIN ? params.scale : WARP_SCALE_MIN;
    return 0.5f * static_cast<float>(curl_intervals(spec.curl_integrator)) /
           (scale * noise_gradient_bound(spec.basis));
  }

  HS_COLD_MEMBER static constexpr bool
  valid_stage_tuple(const WarpStageSpec &spec, const WarpStageParams &params) {
    switch (spec.kind) {
    case WarpStageKind::NONE:
      return true;
    case WarpStageKind::LEGACY_STEREO_NOISE:
      return params.scale >= 0.1f && params.scale <= 100.0f &&
             params.strength >= 0.0f && params.strength <= 30.0f &&
             params.time_scale >= 0.05f && params.time_scale <= 1.0f;
    case WarpStageKind::AFFINE_FRAME:
      return params.translation_x >= -4.0f && params.translation_x <= 4.0f &&
             params.translation_y >= -4.0f && params.translation_y <= 4.0f &&
             params.scale_x >= 0.25f && params.scale_x <= 4.0f &&
             params.scale_y >= 0.25f && params.scale_y <= 4.0f &&
             params.shear >= -0.75f && params.shear <= 0.75f;
    case WarpStageKind::WAVE_SHEAR:
      return params.strength >= -4.0f && params.strength <= 4.0f &&
             params.frequency >= 0.0f && params.frequency <= 64.0f &&
             params.time_scale >= NOISE_RATE_MIN &&
             params.time_scale <= NOISE_RATE_MAX;
    case WarpStageKind::VORTEX:
      return params.radius >= 1.0f / 64.0f && params.radius <= 8.0f &&
             params.turns >= -4.0f && params.turns <= 4.0f &&
             params.center_orbit_radius >= 0.0f &&
             params.center_orbit_radius <= 4.0f &&
             params.time_scale >= NOISE_RATE_MIN &&
             params.time_scale <= NOISE_RATE_MAX;
    case WarpStageKind::VECTOR_NOISE:
      return params.strength >= 0.0f && params.strength <= 4.0f &&
             params.scale >= 1.0f / 64.0f && params.scale <= 64.0f &&
             params.time_scale >= NOISE_RATE_MIN &&
             params.time_scale <= NOISE_RATE_MAX;
    case WarpStageKind::CURL_FLOW:
      return params.strength >= -4.0f && params.strength <= 4.0f &&
             params.scale >= 1.0f / 64.0f && params.scale <= 16.0f &&
             params.time_scale >= NOISE_RATE_MIN &&
             params.time_scale <= NOISE_RATE_MAX &&
             params.scale * abs_value(params.strength) *
                     noise_gradient_bound(spec.basis) /
                     curl_intervals(spec.curl_integrator) <=
                 0.5f;
    case WarpStageKind::MIRROR_TILE:
      return params.cell_x >= CELL_MIN && params.cell_x <= CELL_MAX &&
             params.cell_y >= CELL_MIN && params.cell_y <= CELL_MAX;
    case WarpStageKind::POLAR_CHART:
      return params.radial_scale >= 1.0f / 64.0f &&
             params.radial_scale <= 16.0f;
    }
    return false;
  }

  HS_COLD_MEMBER static constexpr float
  projection_coordinate_bound(const Config &config) {
    if (config.slots.surface_lens == SurfaceLens::MOBIUS)
      return STEREO_INF;
    if (config.slots.projection == Projection::STEREOGRAPHIC)
      return STEREO_INF;
    if (config.slots.projection == Projection::GNOMONIC)
      return 1.0f / GNOMONIC_AXIS_EPS;
    if (strict_projection(config.slots.projection))
      return 16.0f * config.params.projection.coordinate_scale;
    return 4.0f;
  }

  HS_COLD_MEMBER static constexpr float
  stage_coordinate_bound(const WarpStageSpec &spec,
                         const WarpStageParams &params, float input_bound) {
    switch (spec.kind) {
    case WarpStageKind::NONE:
    case WarpStageKind::LEGACY_STEREO_NOISE:
      return input_bound;
    case WarpStageKind::AFFINE_FRAME: {
      const float rotated =
          1.414214f * (input_bound + (abs_value(params.translation_x) >
                                              abs_value(params.translation_y)
                                          ? abs_value(params.translation_x)
                                          : abs_value(params.translation_y)));
      const float x_bound = rotated / params.scale_x +
                            abs_value(params.shear) * rotated / params.scale_y;
      const float y_bound = rotated / params.scale_y;
      return x_bound > y_bound ? x_bound : y_bound;
    }
    case WarpStageKind::WAVE_SHEAR:
      return input_bound + abs_value(params.strength);
    case WarpStageKind::VORTEX: {
      const float center_x =
          abs_value(params.center_x) + params.center_orbit_radius;
      const float center_y =
          abs_value(params.center_y) + params.center_orbit_radius;
      const float center_bound = center_x > center_y ? center_x : center_y;
      return 1.414214f * (input_bound + center_bound) + center_bound;
    }
    case WarpStageKind::VECTOR_NOISE:
      return input_bound + 1.414214f * params.strength;
    case WarpStageKind::CURL_FLOW:
      return input_bound + 1.414214f * abs_value(params.strength) *
                               params.scale * noise_gradient_bound(spec.basis);
    case WarpStageKind::MIRROR_TILE:
      return 1.414214f * (params.cell_x + params.cell_y);
    case WarpStageKind::POLAR_CHART: {
      const float radial =
          params.radial_scale * (spec.polar_mode == PolarMode::LOGARITHMIC
                                     ? 12.0f
                                     : input_bound) +
          TWO_PI_F;
      return radial > 17.0f * PI_F ? radial : 17.0f * PI_F;
    }
    }
    return WARP_COORD_LIMIT + 1.0f;
  }

  HS_COLD_MEMBER static constexpr bool
  safe_program_bounds(const Config &config) {
    float bound = projection_coordinate_bound(config);
    const WarpStageSpec stages[] = {config.slots.warp_program.outer,
                                    config.slots.warp_program.inner};
    const WarpStageParams params[] = {config.params.warp.outer,
                                      config.params.warp.inner};
    for (size_t index = 0; index < 2; ++index) {
      if ((stages[index].kind == WarpStageKind::VECTOR_NOISE ||
           stages[index].kind == WarpStageKind::CURL_FLOW) &&
          params[index].scale * (bound + 100.0f) > NOISE_LATTICE_LIMIT)
        return false;
      bound = stage_coordinate_bound(stages[index], params[index], bound);
      if (bound > WARP_COORD_LIMIT)
        return false;
    }
    if (config.slots.function == Function::NOISE_CONTOUR &&
        config.params.source.noise_scale * bound > NOISE_LATTICE_LIMIT)
      return false;
    return true;
  }

  HS_COLD_MEMBER static constexpr bool
  valid_mobius(const MobiusParams &params) {
    const float ad_re = params.a.re * params.d.re - params.a.im * params.d.im;
    const float ad_im = params.a.re * params.d.im + params.a.im * params.d.re;
    const float bc_re = params.b.re * params.c.re - params.b.im * params.c.im;
    const float bc_im = params.b.re * params.c.im + params.b.im * params.c.re;
    const float det_re = ad_re - bc_re;
    const float det_im = ad_im - bc_im;
    return coefficient_in_range(params.a) && coefficient_in_range(params.b) &&
           coefficient_in_range(params.c) && coefficient_in_range(params.d) &&
           det_re * det_re + det_im * det_im >= 1e-6f;
  }

  HS_COLD_MEMBER static constexpr bool
  coefficient_in_range(const Complex &coefficient) {
    return coefficient.re >= -8.0f && coefficient.re <= 8.0f &&
           coefficient.im >= -8.0f && coefficient.im <= 8.0f;
  }

  HS_COLD_MEMBER static constexpr float max_value(float a, float b) {
    return a > b ? a : b;
  }

  HS_COLD_MEMBER static constexpr float min_value(float a, float b) {
    return a < b ? a : b;
  }

  HS_COLD_MEMBER static constexpr float max_abs_value(float a, float b) {
    return max_value(abs_value(a), abs_value(b));
  }

  HS_COLD_MEMBER static constexpr void
  maximize_stage_path(WarpStageParams &out, const WarpStageParams &a,
                      const WarpStageParams &b) {
    out.translation_x = max_abs_value(a.translation_x, b.translation_x);
    out.translation_y = max_abs_value(a.translation_y, b.translation_y);
    out.scale_x = min_value(a.scale_x, b.scale_x);
    out.scale_y = min_value(a.scale_y, b.scale_y);
    out.shear = max_abs_value(a.shear, b.shear);
    out.strength = max_abs_value(a.strength, b.strength);
    out.scale = max_value(a.scale, b.scale);
    out.center_x = max_abs_value(a.center_x, b.center_x);
    out.center_y = max_abs_value(a.center_y, b.center_y);
    out.center_orbit_radius =
        max_value(a.center_orbit_radius, b.center_orbit_radius);
    out.cell_x = max_value(a.cell_x, b.cell_x);
    out.cell_y = max_value(a.cell_y, b.cell_y);
    out.radial_scale = max_value(a.radial_scale, b.radial_scale);
  }

  HS_COLD_MEMBER static constexpr bool safe_program_path(const Config &from,
                                                         const Config &to) {
    Config worst = from;
    worst.params.projection.coordinate_scale =
        max_value(from.params.projection.coordinate_scale,
                  to.params.projection.coordinate_scale);
    worst.params.source.noise_scale =
        max_value(from.params.source.noise_scale, to.params.source.noise_scale);
    maximize_stage_path(worst.params.warp.outer, from.params.warp.outer,
                        to.params.warp.outer);
    maximize_stage_path(worst.params.warp.inner, from.params.warp.inner,
                        to.params.warp.inner);
    return safe_program_bounds(worst);
  }

  HS_COLD_MEMBER static constexpr bool polar_pair_stable(const Config &from,
                                                         const Config &to) {
    const bool has_polar =
        from.slots.warp_program.outer.kind == WarpStageKind::POLAR_CHART ||
        from.slots.warp_program.inner.kind == WarpStageKind::POLAR_CHART;
    return !has_polar ||
           from.params.source.pattern_freq == to.params.source.pattern_freq;
  }

  HS_COLD_MEMBER static constexpr bool
  stable_parameter_path_admitted(const Config &from, const Config &to) {
    return curl_pair_stable(from.slots.warp_program.outer,
                            from.params.warp.outer, to.params.warp.outer) &&
           curl_pair_stable(from.slots.warp_program.inner,
                            from.params.warp.inner, to.params.warp.inner) &&
           polar_pair_stable(from, to) && safe_program_path(from, to);
  }

  HS_COLD_MEMBER static constexpr bool
  same_parameter_topology(const Config &from, const Config &to) {
    return from.slots == to.slots &&
           from.params.source.noise_basis == to.params.source.noise_basis &&
           from.params.source.noise_seed == to.params.source.noise_seed &&
           from.params.source.noise_resource_id ==
               to.params.source.noise_resource_id &&
           from.params.value.band_count == to.params.value.band_count &&
           from.params.surface_lens.noise_basis ==
               to.params.surface_lens.noise_basis &&
           from.params.surface_lens.noise_seed ==
               to.params.surface_lens.noise_seed &&
           from.params.surface_lens.noise_resource_id ==
               to.params.surface_lens.noise_resource_id &&
           from.params.surface_lens.mobius.a.re ==
               to.params.surface_lens.mobius.a.re &&
           from.params.surface_lens.mobius.a.im ==
               to.params.surface_lens.mobius.a.im &&
           from.params.surface_lens.mobius.b.re ==
               to.params.surface_lens.mobius.b.re &&
           from.params.surface_lens.mobius.b.im ==
               to.params.surface_lens.mobius.b.im &&
           from.params.surface_lens.mobius.c.re ==
               to.params.surface_lens.mobius.c.re &&
           from.params.surface_lens.mobius.c.im ==
               to.params.surface_lens.mobius.c.im &&
           from.params.surface_lens.mobius.d.re ==
               to.params.surface_lens.mobius.d.re &&
           from.params.surface_lens.mobius.d.im ==
               to.params.surface_lens.mobius.d.im;
  }

  HS_COLD_MEMBER static constexpr bool stable_topology(const Config &from,
                                                       const Config &to) {
    return valid_config(from) && valid_config(to) &&
           same_parameter_topology(from, to) &&
           stable_parameter_path_admitted(from, to);
  }

  HS_COLD_MEMBER static constexpr bool
  curl_pair_stable(const WarpStageSpec &spec, const WarpStageParams &a,
                   const WarpStageParams &b) {
    if (spec.kind != WarpStageKind::CURL_FLOW)
      return true;
    const float scale = a.scale > b.scale ? a.scale : b.scale;
    const float a_distance = abs_value(a.strength);
    const float b_distance = abs_value(b.strength);
    const float distance = a_distance > b_distance ? a_distance : b_distance;
    return scale * distance * noise_gradient_bound(spec.basis) /
               curl_intervals(spec.curl_integrator) <=
           0.5f;
  }

  /** @brief Reports whether both transition endpoints are admissible holds. */
  HS_COLD_MEMBER static constexpr bool transition_admitted(const Config &from,
                                                           const Config &to) {
    return valid_config(from) && valid_config(to);
  }

  HS_COLD_MEMBER static constexpr Choreo preset_choreo(size_t index) {
#ifdef HS_PROFILE_SHADERBALL_FAST_CYCLE
    (void)index;
    return {32, 32, 2, false};
#else
    return CHOREO[index];
#endif
  }

  HS_COLD_MEMBER void enter_preset() {
    const Choreo choreo = preset_choreo(getPresetIndex());
    preset_dwell_remaining = static_cast<uint16_t>(
        hs::rand_int(choreo.dwell_min, choreo.dwell_max + 1));
    preset_dwell_armed = true;
  }

  HS_COLD_MEMBER void advance_preset_choreography() {
    if (anims_paused || !preset_dwell_armed)
      return;
    if (preset_dwell_remaining > 0 && --preset_dwell_remaining > 0)
      return;
    preset_dwell_armed = false;
    begin_blend();
  }

  HS_COLD_MEMBER void begin_blend() {
    if (advancePreset()) {
#if defined(HS_PROFILE_ENABLE)
      hs::log("Preset: %u/%u", static_cast<unsigned>(getPresetIndex()),
              static_cast<unsigned>(PRESETS.size()));
#endif
    } else {
      preset_dwell_remaining = 1;
      preset_dwell_armed = true;
    }
  }

  static float golden_rotation(void *context, uint32_t sequence) {
    float &rotation = *static_cast<float *>(context);
    if (sequence > 0)
      rotation = wrap_t(rotation + GOLDEN_HUE_STEP);
    return rotation;
  }

  static void next_liquid_palette(void *context, uint32_t sequence,
                                  GenerativePalette &out) {
    out = GenerativePalette{EffectPaletteRecipes::shader_ball_liquid_at(
        golden_rotation(context, sequence))};
  }

  static void next_triadic_palette(void *context, uint32_t sequence,
                                   GenerativePalette &out) {
    uint32_t &hue = *static_cast<uint32_t *>(context);
    if (sequence > 0)
      hue += HUE_STEP;
    out = GenerativePalette{PaletteRecipes::profile(
        PaletteDomain::STRAIGHT, PaletteHarmony::TRIADIC, AxisCurve::ASCENDING,
        PaletteRecipes::hue_turns(hue))};
  }

  static constexpr uint8_t BOUNDARY_CUT = 1U << 0;
  static constexpr uint8_t BOUNDARY_SINGULAR = 1U << 1;
  static constexpr uint8_t PROJECTION_FLAG_FOLDED = 1U << 0;
  static constexpr float GNOMONIC_AXIS_EPS = 1e-3f;
  static constexpr float WARP_COORD_LIMIT = 65536.0f;
  static constexpr float NOISE_LATTICE_LIMIT = 1048576.0f;
  static constexpr float SPIRAL_ARMS = 3.0f;
  static constexpr float TWIST_RATE = 3.0f;
  static constexpr float KALEIDOSCOPE_SECTORS = 6.0f;
  static constexpr float POLYHEDRAL_MIRROR_EPS = 1e-6f;
  static constexpr int POLYHEDRAL_REFLECTION_LIMIT = 16;
  static constexpr std::array<Vector, 3> TETRAHEDRAL_MIRRORS = {
      Vector(1.0f, 0.0f, 0.0f), Vector(-0.5f, 0.8660254038f, 0.0f),
      Vector(0.0f, -0.5773502692f, 0.8164965809f)};
  static constexpr std::array<Vector, 3> OCTAHEDRAL_MIRRORS = {
      Vector(1.0f, 0.0f, 0.0f), Vector(-0.7071067812f, 0.7071067812f, 0.0f),
      Vector(0.0f, -0.7071067812f, 0.7071067812f)};
  static constexpr std::array<Vector, 3> DODECAHEDRAL_MIRRORS = {
      Vector(1.0f, 0.0f, 0.0f), Vector(-0.8090169944f, 0.3090169944f, -0.5f),
      Vector(0.0f, 0.0f, 1.0f)};
  static constexpr std::array<Vector, 3> TRIANGULAR_PRISM_MIRRORS = {
      Vector(0.0f, 1.0f, 0.0f), Vector(0.0f, 0.0f, 1.0f),
      Vector(0.8660254038f, 0.0f, -0.5f)};
  static constexpr std::array<Vector, 3> SQUARE_PRISM_MIRRORS = {
      Vector(0.0f, 1.0f, 0.0f), Vector(0.0f, 0.0f, 1.0f),
      Vector(0.7071067812f, 0.0f, -0.7071067812f)};
  static constexpr std::array<Vector, 3> PENTAGONAL_PRISM_MIRRORS = {
      Vector(0.0f, 1.0f, 0.0f), Vector(0.0f, 0.0f, 1.0f),
      Vector(0.5877852523f, 0.0f, -0.8090169944f)};
  static constexpr std::array<Vector, 3> HEXAGONAL_PRISM_MIRRORS = {
      Vector(0.0f, 1.0f, 0.0f), Vector(0.0f, 0.0f, 1.0f),
      Vector(0.5f, 0.0f, -0.8660254038f)};
  static constexpr std::array<Vector, 3> OCTAGONAL_PRISM_MIRRORS = {
      Vector(0.0f, 1.0f, 0.0f), Vector(0.0f, 0.0f, 1.0f),
      Vector(0.3826834324f, 0.0f, -0.9238795325f)};
  static constexpr float WARP_NOISE_FREQUENCY = 0.01f;
  static constexpr float STEREO_NOISE_TIME_PERIOD = 65536.0f;
  static constexpr float NOISE_NATIVE_PERIOD = 256.0f;
  static constexpr float NOISE_WRAP_START = 63.0f / 64.0f;
  static constexpr float ONE_BELOW_UNIT = 0x1.fffffep-1f;
  static constexpr float NOISE_WRAP_BLEND = 1024.0f;
  static constexpr float GOLDEN_HUE_STEP = 0.618034f;
  static constexpr uint32_t HUE_STEP = 159;
  static constexpr int PALETTE_DWELL_FRAMES = 0;
  static constexpr int PALETTE_FADE_FRAMES = 600;
  static constexpr size_t PARAM_CAPACITY = 64;

  static constexpr const char *FUNCTION_OPTIONS[] = {
      "Twin Wave", "Rings",         "Spiral",
      "Grid",      "Noise Contour", "Primitive Lattice"};
  static constexpr const char *FUNCTION_EXPORT_OPTIONS[] = {
      "Function::TWIN_WAVE",     "Function::RINGS",
      "Function::SPIRAL",        "Function::GRID",
      "Function::NOISE_CONTOUR", "Function::PRIMITIVE_LATTICE"};
  static constexpr int NUM_FUNCTIONS = std::size(FUNCTION_OPTIONS);
  static constexpr const char *PROJECTION_OPTIONS[] = {
      "Folded Sinusoidal",  "Stereographic",       "Gnomonic",       "Bonne",
      "Peirce Quincuncial", "Dymaxion / Airocean", "Equirectangular"};
  static constexpr const char *PROJECTION_EXPORT_OPTIONS[] = {
      "Projection::SINUSOIDAL",         "Projection::STEREOGRAPHIC",
      "Projection::GNOMONIC",           "Projection::BONNE",
      "Projection::PEIRCE_QUINCUNCIAL", "Projection::AIROCEAN",
      "Projection::EQUIRECTANGULAR"};
  static constexpr int NUM_PROJECTIONS = std::size(PROJECTION_OPTIONS);
  static constexpr const char *PEIRCE_LAYOUT_OPTIONS[] = {
      "Diamond", "Square", "Horizontal", "Vertical"};
  static constexpr const char *PEIRCE_LAYOUT_EXPORT_OPTIONS[] = {
      "PeirceLayout::DIAMOND", "PeirceLayout::SQUARE",
      "PeirceLayout::HORIZONTAL", "PeirceLayout::VERTICAL"};
  static constexpr int NUM_PEIRCE_LAYOUTS = std::size(PEIRCE_LAYOUT_OPTIONS);
  static constexpr const char *AIROCEAN_LAYOUT_OPTIONS[] = {"Vertical",
                                                            "Horizontal"};
  static constexpr const char *AIROCEAN_LAYOUT_EXPORT_OPTIONS[] = {
      "AiroceanLayout::VERTICAL", "AiroceanLayout::HORIZONTAL"};
  static constexpr int NUM_AIROCEAN_LAYOUTS =
      std::size(AIROCEAN_LAYOUT_OPTIONS);
  static constexpr const char *BONNE_HEMISPHERE_OPTIONS[] = {"North", "South"};
  static constexpr const char *BONNE_HEMISPHERE_EXPORT_OPTIONS[] = {
      "BonneHemisphere::NORTH", "BonneHemisphere::SOUTH"};
  static constexpr int NUM_BONNE_HEMISPHERES =
      std::size(BONNE_HEMISPHERE_OPTIONS);
  static constexpr const char *GNOMONIC_HEMISPHERE_OPTIONS[] = {
      "Folded", "Front Hemisphere", "Back Hemisphere"};
  static constexpr const char *GNOMONIC_HEMISPHERE_EXPORT_OPTIONS[] = {
      "GnomonicHemispherePolicy::FOLDED",
      "GnomonicHemispherePolicy::FRONT_HEMISPHERE",
      "GnomonicHemispherePolicy::BACK_HEMISPHERE"};
  static constexpr int NUM_GNOMONIC_HEMISPHERES =
      std::size(GNOMONIC_HEMISPHERE_OPTIONS);
  static constexpr const char *PROJECTION_FRAME_OPTIONS[] = {"Identity",
                                                             "Spin + Wander"};
  static constexpr const char *PROJECTION_FRAME_EXPORT_OPTIONS[] = {
      "ProjectionFramePolicy::IDENTITY", "ProjectionFramePolicy::SPIN_WANDER"};
  static constexpr int NUM_PROJECTION_FRAMES =
      std::size(PROJECTION_FRAME_OPTIONS);
  static constexpr const char *LENS_OPTIONS[] = {
      "None",
      "Glitch",
      "Twist",
      "Kaleidoscope (Azimuthal 6-fold)",
      "Mobius",
      "Tangent Noise",
      "Kaleidoscope (Tetrahedral)",
      "Kaleidoscope (Octahedral / Cubic)",
      "Kaleidoscope (Dodecahedral / Icosahedral)",
      "Kaleidoscope (Triangular Prism)",
      "Kaleidoscope (Square Prism)",
      "Kaleidoscope (Pentagonal Prism)",
      "Kaleidoscope (Hexagonal Prism)",
      "Kaleidoscope (Octagonal Prism)"};
  static constexpr const char *LENS_EXPORT_OPTIONS[] = {
      "SurfaceLens::NONE",
      "SurfaceLens::GLITCH",
      "SurfaceLens::TWIST",
      "SurfaceLens::KALEIDOSCOPE",
      "SurfaceLens::MOBIUS",
      "SurfaceLens::TANGENT_NOISE",
      "SurfaceLens::KALEIDOSCOPE_TETRAHEDRAL",
      "SurfaceLens::KALEIDOSCOPE_OCTAHEDRAL",
      "SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL",
      "SurfaceLens::KALEIDOSCOPE_TRIANGULAR_PRISM",
      "SurfaceLens::KALEIDOSCOPE_SQUARE_PRISM",
      "SurfaceLens::KALEIDOSCOPE_PENTAGONAL_PRISM",
      "SurfaceLens::KALEIDOSCOPE_HEXAGONAL_PRISM",
      "SurfaceLens::KALEIDOSCOPE_OCTAGONAL_PRISM"};
  static constexpr int NUM_LENSES = std::size(LENS_OPTIONS);
  static constexpr const char *WARP_OPTIONS[] = {
      "None",         "Stereo Noise", "Affine Frame", "Wave Shear", "Vortex",
      "Vector Noise", "Curl Flow",    "Mirror Tile",  "Polar Chart"};
  static constexpr const char *WARP_EXPORT_OPTIONS[] = {
      "WarpStageKind::NONE",         "WarpStageKind::LEGACY_STEREO_NOISE",
      "WarpStageKind::AFFINE_FRAME", "WarpStageKind::WAVE_SHEAR",
      "WarpStageKind::VORTEX",       "WarpStageKind::VECTOR_NOISE",
      "WarpStageKind::CURL_FLOW",    "WarpStageKind::MIRROR_TILE",
      "WarpStageKind::POLAR_CHART"};
  static constexpr int NUM_WARPS = std::size(WARP_OPTIONS);
  static constexpr const char *NOISE_BASIS_OPTIONS[] = {"Simplex", "FBM 3",
                                                        "Ridged 3"};
  static constexpr const char *NOISE_BASIS_EXPORT_OPTIONS[] = {
      "NoiseBasis::SIMPLEX", "NoiseBasis::FBM3", "NoiseBasis::RIDGED3"};
  static constexpr int NUM_NOISE_BASES = std::size(NOISE_BASIS_OPTIONS);
  static constexpr const char *POLAR_MODE_OPTIONS[] = {"Linear", "Logarithmic"};
  static constexpr const char *POLAR_MODE_EXPORT_OPTIONS[] = {
      "PolarMode::LINEAR", "PolarMode::LOGARITHMIC"};
  static constexpr int NUM_POLAR_MODES = std::size(POLAR_MODE_OPTIONS);
  static constexpr const char *CURL_INTEGRATOR_OPTIONS[] = {
      "Euler 1", "Midpoint 2", "Midpoint 4"};
  static constexpr const char *CURL_INTEGRATOR_EXPORT_OPTIONS[] = {
      "CurlIntegrator::EULER_1", "CurlIntegrator::MIDPOINT_2",
      "CurlIntegrator::MIDPOINT_4"};
  static constexpr int NUM_CURL_INTEGRATORS =
      std::size(CURL_INTEGRATOR_OPTIONS);
  static constexpr int POLAR_HARMONIC_MAX = 16;
  static constexpr int BAND_COUNT_MAX = 32;
  static constexpr const char *WARP_ENVELOPE_OPTIONS[] = {
      "Flat", "Projection Weight", "Edge Fade"};
  static constexpr const char *WARP_ENVELOPE_EXPORT_OPTIONS[] = {
      "WarpEnvelope::FLAT", "WarpEnvelope::PROJECTION_WEIGHT",
      "WarpEnvelope::EDGE_FADE"};
  static constexpr int NUM_WARP_ENVELOPES = std::size(WARP_ENVELOPE_OPTIONS);
  static constexpr const char *SIGNAL_OPTIONS[] = {"None", "Projection"};
  static constexpr const char *SIGNAL_EXPORT_OPTIONS[] = {
      "SignalWeight::NONE", "SignalWeight::PROJECTION"};
  static constexpr int NUM_SIGNALS = std::size(SIGNAL_OPTIONS);
  static constexpr const char *VALUE_TRANSFER_OPTIONS[] = {
      "Linear", "Ridge", "Iso Contour", "Smooth Bands"};
  static constexpr const char *VALUE_TRANSFER_EXPORT_OPTIONS[] = {
      "ValueTransfer::LINEAR", "ValueTransfer::RIDGE",
      "ValueTransfer::ISO_CONTOUR", "ValueTransfer::SMOOTH_BANDS"};
  static constexpr int NUM_VALUE_TRANSFERS = std::size(VALUE_TRANSFER_OPTIONS);
  static constexpr const char *COVERAGE_OPTIONS[] = {
      "Opaque", "Projection Weight Squared", "Value Cutout", "Edge Fade",
      "Projection Weight"};
  static constexpr const char *COVERAGE_EXPORT_OPTIONS[] = {
      "CoveragePolicy::OPAQUE", "CoveragePolicy::PROJECTION_WEIGHT_SQUARED",
      "CoveragePolicy::VALUE_CUTOUT", "CoveragePolicy::EDGE_FADE",
      "CoveragePolicy::PROJECTION_WEIGHT"};
  static constexpr int NUM_COVERAGE_POLICIES = std::size(COVERAGE_OPTIONS);
  static constexpr const char *COLORIZER_OPTIONS[] = {
      "Generated Triadic", "ShaderBall Liquid", "Deformation Ink"};
  static constexpr const char *COLORIZER_EXPORT_OPTIONS[] = {
      "Colorizer::GENERATED_TRIADIC", "Colorizer::LIQUID",
      "Colorizer::DEFORMATION_INK"};
  static constexpr int NUM_COLORIZERS = std::size(COLORIZER_OPTIONS);

  static constexpr float WARP_SCALE_MIN = 1.0f / 64.0f;
  static constexpr float WARP_SCALE_MAX = 100.0f;
  static constexpr float WARP_STRENGTH_MIN = -4.0f;
  static constexpr float WARP_STRENGTH_MAX = 30.0f;
  static constexpr float WARP_TIME_MIN = -1.0f / 64.0f;
  static constexpr float WARP_TIME_MAX = 1.0f;
  static constexpr float PATTERN_FREQ_MIN = 0.1f, PATTERN_FREQ_MAX = 20.0f;
  static constexpr float SPEED_MIN = 0.0f, SPEED_MAX = 5.0f;
  static constexpr float COMPLEXITY_MIN = 0.0f, COMPLEXITY_MAX = 3.0f;
  static constexpr float PATTERN_MIX_MIN = 0.0f, PATTERN_MIX_MAX = 1.0f;
  static constexpr float PHASE2_RATE_MIN = 0.0f, PHASE2_RATE_MAX = 2.0f;
  static constexpr float POLE_FADE_MIN = 1.0f, POLE_FADE_MAX = 20.0f;
  static constexpr float SPIN_RATE_MIN = 0.0f, SPIN_RATE_MAX = 0.05f;
  static constexpr float WANDER_MIN = 0.0f, WANDER_MAX = 1.0f;
  static constexpr float LENS_MIX_MIN = 0.0f, LENS_MIX_MAX = 1.0f;
  static constexpr float BREATHE_MIN = 0.0f, BREATHE_MAX = 0.3f;
  static constexpr float CYCLE_SPEED_MIN = 0.0f, CYCLE_SPEED_MAX = 1.0f;
  static constexpr float HUE_SHIFT_MIN = 0.0f, HUE_SHIFT_MAX = 1.0f;
  static constexpr float VALUE_FADE_MIN = 0.0f, VALUE_FADE_MAX = 1.0f;
  static constexpr float WAVE_SPIN_MIN = 0.0f, WAVE_SPIN_MAX = 0.05f;
  static constexpr float SOURCE_NOISE_SCALE_MIN = 0.0f;
  static constexpr float SOURCE_NOISE_SCALE_MAX = 2.0f;
  static constexpr float SOURCE_NOISE_RATE_MIN = -1.0f / 1024.0f;
  static constexpr float SOURCE_NOISE_RATE_MAX = 1.0f / 1024.0f;
  static constexpr float LENS_NOISE_SCALE_MIN = 1.0f / 64.0f;
  static constexpr float LENS_NOISE_SCALE_MAX = 8.0f;
  static constexpr float NOISE_RATE_MIN = -1.0f / 64.0f;
  static constexpr float NOISE_RATE_MAX = 1.0f / 64.0f;
  static constexpr float CELL_MIN = 1.0f / 64.0f;
  static constexpr float CELL_MAX = 8.0f;
  static constexpr float SOFTNESS_MIN = 1.0f / 1024.0f;
  static constexpr float ORBIT_SPIN_RATE = TWO_PI_F / 300.0f;

  HS_COLD_MEMBER static constexpr bool preset_in_ranges(const Params &p) {
    return warp_stage_params_in_ranges(p.warp.outer) &&
           warp_stage_params_in_ranges(p.warp.inner) &&
           p.source.pattern_freq >= PATTERN_FREQ_MIN &&
           p.source.pattern_freq <= PATTERN_FREQ_MAX &&
           p.source.speed >= SPEED_MIN && p.source.speed <= SPEED_MAX &&
           p.source.complexity >= COMPLEXITY_MIN &&
           p.source.complexity <= COMPLEXITY_MAX &&
           p.source.pattern_mix >= PATTERN_MIX_MIN &&
           p.source.pattern_mix <= PATTERN_MIX_MAX &&
           p.source.secondary_rate >= PHASE2_RATE_MIN &&
           p.source.secondary_rate <= PHASE2_RATE_MAX &&
           p.source.angle_rate >= WAVE_SPIN_MIN &&
           p.source.angle_rate <= WAVE_SPIN_MAX &&
           p.source.noise_scale >= SOURCE_NOISE_SCALE_MIN &&
           p.source.noise_scale <= SOURCE_NOISE_SCALE_MAX &&
           p.source.noise_contrast >= 0.0f && p.source.noise_contrast <= 8.0f &&
           p.source.noise_time_rate >= SOURCE_NOISE_RATE_MIN &&
           p.source.noise_time_rate <= SOURCE_NOISE_RATE_MAX &&
           p.source.lattice_cell_scale >= CELL_MIN &&
           p.source.lattice_cell_scale <= CELL_MAX &&
           p.source.lattice_shape_blend >= 0.0f &&
           p.source.lattice_shape_blend <= 1.0f &&
           p.source.lattice_softness >= SOFTNESS_MIN &&
           p.source.lattice_softness <= 1.0f &&
           p.source.lattice_radius >= 1.0f / 64.0f &&
           p.source.lattice_radius <= 0.49f &&
           enum_at_most(p.source.noise_basis, NoiseBasis::RIDGED3) &&
           p.projection.pole_fade >= POLE_FADE_MIN &&
           p.projection.pole_fade <= POLE_FADE_MAX &&
           p.projection.spin_rate >= SPIN_RATE_MIN &&
           p.projection.spin_rate <= SPIN_RATE_MAX &&
           p.projection.wander >= WANDER_MIN &&
           p.projection.wander <= WANDER_MAX &&
           p.projection.central_meridian >= 0.0f &&
           p.projection.central_meridian <= TWO_PI_F &&
           p.projection.coordinate_scale >= 0.25f &&
           p.projection.coordinate_scale <= 4.0f &&
           p.projection.bonne_standard_parallel >= 1e-3f &&
           p.projection.bonne_standard_parallel <= 0.5f * PI_F &&
           p.projection.layout_scroll >= -1.0f &&
           p.projection.layout_scroll <= 1.0f &&
           p.outer_camera.wander >= WANDER_MIN &&
           p.outer_camera.wander <= WANDER_MAX &&
           p.surface_lens.mix >= LENS_MIX_MIN &&
           p.surface_lens.mix <= LENS_MIX_MAX &&
           p.surface_lens.amount >= 0.0f && p.surface_lens.amount <= 4.0f &&
           p.surface_lens.noise_scale >= LENS_NOISE_SCALE_MIN &&
           p.surface_lens.noise_scale <= LENS_NOISE_SCALE_MAX &&
           p.surface_lens.noise_rate >= NOISE_RATE_MIN &&
           p.surface_lens.noise_rate <= NOISE_RATE_MAX &&
           p.value.iso_level >= 0.0f && p.value.iso_level <= 1.0f &&
           p.value.iso_width >= SOFTNESS_MIN && p.value.iso_width <= 0.5f &&
           p.value.band_count >= 1 && p.value.band_count <= BAND_COUNT_MAX &&
           p.value.band_phase >= 0.0f && p.value.band_phase <= TWO_PI_F &&
           p.value.cutout_threshold >= 0.0f &&
           p.value.cutout_threshold <= 1.0f &&
           p.value.cutout_softness >= SOFTNESS_MIN &&
           p.value.cutout_softness <= 0.5f && p.value.edge_width >= 0.0f &&
           p.value.edge_width <= 0.5f &&
           p.colorizer.breathe_depth >= BREATHE_MIN &&
           p.colorizer.breathe_depth <= BREATHE_MAX &&
           p.colorizer.cycle_speed >= CYCLE_SPEED_MIN &&
           p.colorizer.cycle_speed <= CYCLE_SPEED_MAX &&
           p.colorizer.hue_shift >= HUE_SHIFT_MIN &&
           p.colorizer.hue_shift <= HUE_SHIFT_MAX &&
           p.colorizer.value_fade >= VALUE_FADE_MIN &&
           p.colorizer.value_fade <= VALUE_FADE_MAX &&
           p.colorizer.displacement_gain >= -4.0f &&
           p.colorizer.displacement_gain <= 4.0f &&
           p.colorizer.path_gain >= -4.0f && p.colorizer.path_gain <= 4.0f &&
           p.colorizer.direction_gain >= -4.0f &&
           p.colorizer.direction_gain <= 4.0f &&
           p.colorizer.displacement_norm >= SOFTNESS_MIN &&
           p.colorizer.displacement_norm <= 32.0f &&
           p.colorizer.path_norm >= SOFTNESS_MIN &&
           p.colorizer.path_norm <= 32.0f &&
           p.colorizer.direction_phase >= 0.0f &&
           p.colorizer.direction_phase <= TWO_PI_F;
  }

  HS_COLD_MEMBER static constexpr bool
  warp_stage_params_in_ranges(const WarpStageParams &params) {
    return params.scale >= WARP_SCALE_MIN && params.scale <= WARP_SCALE_MAX &&
           params.strength >= WARP_STRENGTH_MIN &&
           params.strength <= WARP_STRENGTH_MAX &&
           params.time_scale >= WARP_TIME_MIN &&
           params.time_scale <= WARP_TIME_MAX &&
           params.translation_x >= -4.0f && params.translation_x <= 4.0f &&
           params.translation_y >= -4.0f && params.translation_y <= 4.0f &&
           params.rotation >= 0.0f && params.rotation <= TWO_PI_F &&
           params.scale_x >= 0.25f && params.scale_x <= 4.0f &&
           params.scale_y >= 0.25f && params.scale_y <= 4.0f &&
           params.shear >= -0.75f && params.shear <= 0.75f &&
           params.frequency >= 0.0f && params.frequency <= 64.0f &&
           params.field_angle >= 0.0f && params.field_angle <= TWO_PI_F &&
           params.center_x >= -4.0f && params.center_x <= 4.0f &&
           params.center_y >= -4.0f && params.center_y <= 4.0f &&
           params.radius >= 1.0f / 64.0f && params.radius <= 8.0f &&
           params.turns >= -4.0f && params.turns <= 4.0f &&
           params.vector_angle >= 0.0f && params.vector_angle <= TWO_PI_F &&
           params.cell_x >= CELL_MIN && params.cell_x <= CELL_MAX &&
           params.cell_y >= CELL_MIN && params.cell_y <= CELL_MAX &&
           params.offset_x >= -8.0f && params.offset_x <= 8.0f &&
           params.offset_y >= -8.0f && params.offset_y <= 8.0f &&
           params.radial_scale >= 1.0f / 64.0f &&
           params.radial_scale <= 16.0f && params.radial_phase >= 0.0f &&
           params.radial_phase <= TWO_PI_F && params.angular_phase >= 0.0f &&
           params.angular_phase <= TWO_PI_F &&
           params.edge_width >= SOFTNESS_MIN && params.edge_width <= 0.5f;
  }

  static constexpr Slots LIQUID_STEREO_SLOTS{
      Function::GRID,
      Projection::STEREOGRAPHIC,
      ProjectionFramePolicy::SPIN_WANDER,
      SurfaceLens::GLITCH,
      {{WarpStageKind::LEGACY_STEREO_NOISE}, {WarpStageKind::NONE}},
      SignalWeight::PROJECTION,
      ValueTransfer::LINEAR,
      CoveragePolicy::OPAQUE,
      Colorizer::LIQUID,
      PeirceLayout::SQUARE,
      AiroceanLayout::VERTICAL};
  static constexpr Slots KALEIDOSCOPE_LIQUID_STEREO_SLOTS = [] {
    Slots slots = LIQUID_STEREO_SLOTS;
    slots.surface_lens = SurfaceLens::KALEIDOSCOPE;
    return slots;
  }();

  /** @brief Index of a warp parameter name in the per-position name tables. */
  enum WarpParamName : uint8_t {
    WARP_NAME_TRANSLATION_X,
    WARP_NAME_TRANSLATION_Y,
    WARP_NAME_ROTATION,
    WARP_NAME_SCALE_X,
    WARP_NAME_SCALE_Y,
    WARP_NAME_SHEAR,
    WARP_NAME_FREQUENCY,
    WARP_NAME_FIELD_ANGLE,
    WARP_NAME_CENTER_X,
    WARP_NAME_CENTER_Y,
    WARP_NAME_RADIUS,
    WARP_NAME_TURNS,
    WARP_NAME_VECTOR_ANGLE,
    WARP_NAME_CELL_X,
    WARP_NAME_CELL_Y,
    WARP_NAME_OFFSET_X,
    WARP_NAME_OFFSET_Y,
    WARP_NAME_RADIAL_SCALE,
    WARP_NAME_RADIAL_PHASE,
    WARP_NAME_ANGULAR_PHASE,
    WARP_NAME_EDGE_WIDTH,
    WARP_NAME_CENTER_ORBIT,
    WARP_NAME_COUNT,
  };

  static constexpr const char *FIRST_WARP_PARAM_NAMES[] = {
      "Planar Warp 1 Translation X", "Planar Warp 1 Translation Y",
      "Planar Warp 1 Rotation",      "Planar Warp 1 Scale X",
      "Planar Warp 1 Scale Y",       "Planar Warp 1 Shear",
      "Planar Warp 1 Frequency",     "Planar Warp 1 Field Angle",
      "Planar Warp 1 Center X",      "Planar Warp 1 Center Y",
      "Planar Warp 1 Radius",        "Planar Warp 1 Turns",
      "Planar Warp 1 Vector Angle",  "Planar Warp 1 Cell X",
      "Planar Warp 1 Cell Y",        "Planar Warp 1 Offset X",
      "Planar Warp 1 Offset Y",      "Planar Warp 1 Radial Scale",
      "Planar Warp 1 Radial Phase",  "Planar Warp 1 Angular Phase",
      "Planar Warp 1 Edge Width",    "Planar Warp 1 Center Orbit"};
  static constexpr const char *SECOND_WARP_PARAM_NAMES[] = {
      "Planar Warp 2 Translation X", "Planar Warp 2 Translation Y",
      "Planar Warp 2 Rotation",      "Planar Warp 2 Scale X",
      "Planar Warp 2 Scale Y",       "Planar Warp 2 Shear",
      "Planar Warp 2 Frequency",     "Planar Warp 2 Field Angle",
      "Planar Warp 2 Center X",      "Planar Warp 2 Center Y",
      "Planar Warp 2 Radius",        "Planar Warp 2 Turns",
      "Planar Warp 2 Vector Angle",  "Planar Warp 2 Cell X",
      "Planar Warp 2 Cell Y",        "Planar Warp 2 Offset X",
      "Planar Warp 2 Offset Y",      "Planar Warp 2 Radial Scale",
      "Planar Warp 2 Radial Phase",  "Planar Warp 2 Angular Phase",
      "Planar Warp 2 Edge Width",    "Planar Warp 2 Center Orbit"};
  static_assert(sizeof(FIRST_WARP_PARAM_NAMES) / sizeof(const char *) ==
                    WARP_NAME_COUNT,
                "first warp name table must match WarpParamName");
  static_assert(sizeof(SECOND_WARP_PARAM_NAMES) / sizeof(const char *) ==
                    WARP_NAME_COUNT,
                "second warp name table must match WarpParamName");
  static constexpr Params
  authored_params(SourceParams source, WarpStageParams outer_warp,
                  ProjectionParams projection, SurfaceLensParams surface_lens,
                  ColorizerParams colorizer, OuterCameraParams outer_camera) {
    const WarpStageParams inner_warp{0.1f, 0.0f, outer_warp.time_scale};
    projection.wander = outer_camera.wander;
    return {source,      {outer_warp, inner_warp},
            projection,  surface_lens,
            {},          colorizer,
            outer_camera};
  }

  static constexpr Preset wave_shear_liquid_preset() {
    Slots slots = LIQUID_STEREO_SLOTS;
    slots.warp_program.outer.kind = WarpStageKind::WAVE_SHEAR;
    slots.coverage = CoveragePolicy::PROJECTION_WEIGHT_SQUARED;
    Params params = authored_params(
        {4.439f, 0.245f, 0.5f, 0.0f, 0.0f, 0.0f}, {1.0f, 0.5f, 1.0f / 64.0f},
        {1.0f, 0.0f, 0.0f}, {1.0f}, {0.133f, 0.05f, 0.0f, 0.0f}, {1.0f});
    params.projection.wander = 0.0f;
    return {slots, params};
  }

  static constexpr Preset kaleidoscope_mirror_preset() {
    const Slots slots{Function::TWIN_WAVE,
                      Projection::STEREOGRAPHIC,
                      ProjectionFramePolicy::SPIN_WANDER,
                      SurfaceLens::KALEIDOSCOPE,
                      {{WarpStageKind::NONE}, {WarpStageKind::MIRROR_TILE}},
                      SignalWeight::PROJECTION,
                      ValueTransfer::LINEAR,
                      CoveragePolicy::PROJECTION_WEIGHT_SQUARED,
                      Colorizer::LIQUID};
    const Params params = authored_params(
        {10.158f, 0.245f, 0.513f, 0.0f, 0.8f, 0.027f}, {0.1f, 0.0f, 0.5f},
        {4.971f, 0.0f, 1.0f}, {1.0f}, {0.15f, 0.0f, 0.0f, 0.0f}, {1.0f});
    return {slots, params};
  }

  static constexpr Preset gnomonic_grid_mirror_preset(SurfaceLens lens) {
    Slots slots{Function::GRID,
                Projection::GNOMONIC,
                ProjectionFramePolicy::IDENTITY,
                lens,
                {{WarpStageKind::MIRROR_TILE}, {WarpStageKind::NONE}},
                SignalWeight::PROJECTION,
                ValueTransfer::LINEAR,
                CoveragePolicy::EDGE_FADE,
                Colorizer::GENERATED_TRIADIC};
    slots.gnomonic_hemisphere = GnomonicHemispherePolicy::FOLDED;
    WarpStageParams outer_warp;
    outer_warp.rotation = 0.29530972f;
    outer_warp.cell_x = 5.381125f;
    outer_warp.cell_y = 1.0f;
    outer_warp.offset_x = 1.344f;
    outer_warp.offset_y = -1.456f;
    Params params =
        authored_params({3.565f, 0.235f, 0.0f, 1.0f, 1.0f, 0.0f}, outer_warp,
                        {1.4f, 0.0f}, {1.0f}, {}, {1.0f});
    params.value.edge_width = 0.5f;
    return {slots, params};
  }

  static constexpr Preset bonne_lattice_mirror_preset() {
    Slots slots{Function::PRIMITIVE_LATTICE,
                Projection::BONNE,
                ProjectionFramePolicy::SPIN_WANDER,
                SurfaceLens::KALEIDOSCOPE,
                {{WarpStageKind::MIRROR_TILE}, {WarpStageKind::NONE}},
                SignalWeight::PROJECTION,
                ValueTransfer::LINEAR,
                CoveragePolicy::EDGE_FADE,
                Colorizer::GENERATED_TRIADIC};
    slots.bonne_hemisphere = BonneHemisphere::NORTH;
    WarpStageParams outer_warp;
    outer_warp.rotation = 1.7215928f;
    outer_warp.cell_x = 5.381125f;
    outer_warp.cell_y = 1.0f;
    outer_warp.offset_x = 1.344f;
    outer_warp.offset_y = -1.456f;
    Params params =
        authored_params({}, outer_warp, {1.0f, 0.0f}, {1.0f}, {}, {1.0f});
    params.source.lattice_cell_scale = 1.1494062f;
    params.source.lattice_shape_blend = 1.0f;
    params.source.lattice_softness = 0.26372f;
    params.source.lattice_radius = 0.31164f;
    params.projection.central_meridian = 0.0f;
    params.projection.coordinate_scale = 1.0f;
    params.projection.bonne_standard_parallel = 0.001f;
    params.value.edge_width = 0.5f;
    return {slots, params};
  }

  static constexpr Preset peirce_lattice_preset() {
    Slots slots{Function::PRIMITIVE_LATTICE,
                Projection::PEIRCE_QUINCUNCIAL,
                ProjectionFramePolicy::SPIN_WANDER,
                SurfaceLens::KALEIDOSCOPE,
                {{WarpStageKind::NONE}, {WarpStageKind::NONE}},
                SignalWeight::PROJECTION,
                ValueTransfer::LINEAR,
                CoveragePolicy::EDGE_FADE,
                Colorizer::LIQUID};
    slots.peirce_layout = PeirceLayout::SQUARE;
    Params params = authored_params({}, {}, {1.0f, 0.0f}, {1.0f},
                                    {0.15f, 0.05f, 0.0f, 0.0f}, {1.0f});
    params.source.lattice_cell_scale = 2.2911718f;
    params.source.lattice_shape_blend = 1.0f;
    params.source.lattice_softness = 0.5804102f;
    params.source.lattice_radius = 0.25f;
    params.projection.central_meridian = 0.0f;
    params.projection.coordinate_scale = 1.0f;
    params.value.edge_width = 0.1f;
    return {slots, params};
  }

  static constexpr Preset kaleidoscope_edge_fade_liquid_preset() {
    Slots slots = KALEIDOSCOPE_LIQUID_STEREO_SLOTS;
    slots.coverage = CoveragePolicy::EDGE_FADE;
    Params params = authored_params(
        {4.116f, 0.1f, 0.5f, 0.0f, 0.8f}, {19.7803f, 16.74f, 0.5f},
        {1.4f, 0.0f}, {1.0f}, {0.15f, 0.05f, 0.657f, 0.0f}, {1.0f});
    params.value.edge_width = 0.2575f;
    return {slots, params};
  }

  static constexpr Preset dodecahedral_grid_preset() {
    const Slots slots{
        Function::GRID,
        Projection::STEREOGRAPHIC,
        ProjectionFramePolicy::SPIN_WANDER,
        SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL,
        {{WarpStageKind::MIRROR_TILE}, {WarpStageKind::LEGACY_STEREO_NOISE}},
        SignalWeight::PROJECTION,
        ValueTransfer::LINEAR,
        CoveragePolicy::EDGE_FADE,
        Colorizer::LIQUID};
    WarpStageParams outer_warp;
    outer_warp.cell_x = 1.8041f;
    outer_warp.cell_y = 1.7083f;
    Params params =
        authored_params({1.532f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f}, outer_warp,
                        {3.907f, 0.0387f, 0.0f}, {1.0f},
                        {0.25410002f, 0.00015458837f, 0.339f, 0.847f}, {1.0f});
    params.projection.wander = 0.0f;
    params.warp.inner = {24.8752f, 10.5f, 0.05f};
    params.value.edge_width = 0.0f;
    return {slots, params};
  }

  static constexpr Preset peirce_dodecahedral_liquid_preset() {
    Slots slots{Function::GRID,
                Projection::PEIRCE_QUINCUNCIAL,
                ProjectionFramePolicy::SPIN_WANDER,
                SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL,
                {{WarpStageKind::NONE}, {WarpStageKind::NONE}},
                SignalWeight::PROJECTION,
                ValueTransfer::LINEAR,
                CoveragePolicy::EDGE_FADE,
                Colorizer::LIQUID};
    slots.peirce_layout = PeirceLayout::SQUARE;
    Params params =
        authored_params({5.0f, 0.1f, 0.5f, 0.0f, 0.8f, 0.0f}, {}, {1.0f, 0.0f},
                        {1.0f}, {0.15f, 0.05f, 0.319f, 0.2f}, {1.0f});
    params.projection.central_meridian = 0.0f;
    params.projection.coordinate_scale = 1.0f;
    params.value.edge_width = 0.1f;
    return {slots, params};
  }

  static constexpr Preset dodecahedral_noise_liquid_preset() {
    Slots slots = KALEIDOSCOPE_LIQUID_STEREO_SLOTS;
    slots.surface_lens = SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL;
    const Params params = authored_params(
        {1.0f, 0.075f, 0.009122372f, 1.0f, 1.146f, 0.0f},
        {50.749298f, 30.0f, 0.4699f}, {1.5482996f, 0.020879198f, 0.0030917525f},
        {1.0f}, {0.25410002f, 0.00015458837f, 0.201f, 0.847f}, {0.0030917525f});
    return {slots, params};
  }

  static constexpr Preset dodecahedral_lattice_noise_preset() {
    const Slots slots{
        Function::PRIMITIVE_LATTICE,
        Projection::STEREOGRAPHIC,
        ProjectionFramePolicy::SPIN_WANDER,
        SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL,
        {{WarpStageKind::MIRROR_TILE}, {WarpStageKind::LEGACY_STEREO_NOISE}},
        SignalWeight::PROJECTION,
        ValueTransfer::LINEAR,
        CoveragePolicy::EDGE_FADE,
        Colorizer::LIQUID};
    SourceParams source;
    source.lattice_cell_scale = 0.57453126f;
    source.lattice_shape_blend = 1.0f;
    source.lattice_softness = 0.5454443f;
    source.lattice_radius = 0.25f;
    WarpStageParams outer_warp;
    outer_warp.rotation = 0.0f;
    outer_warp.cell_x = 1.0f;
    outer_warp.cell_y = 1.0f;
    outer_warp.offset_x = 0.0f;
    outer_warp.offset_y = 0.0f;
    Params params = authored_params(
        source, outer_warp, {1.5482996f, 0.020879198f, 0.421f}, {1.0f},
        {0.25410002f, 0.00015458837f, 0.562f, 0.847f}, {0.591f});
    params.projection.wander = 0.421f;
    params.warp.inner = {16.0f, 8.73f, 0.20772f};
    params.value.edge_width = 0.0775f;
    return {slots, params};
  }

  static constexpr Preset gnomonic_wave_shear_grid_preset() {
    Slots slots{Function::GRID,
                Projection::GNOMONIC,
                ProjectionFramePolicy::IDENTITY,
                SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL,
                {{WarpStageKind::WAVE_SHEAR}, {WarpStageKind::MIRROR_TILE}},
                SignalWeight::PROJECTION,
                ValueTransfer::LINEAR,
                CoveragePolicy::PROJECTION_WEIGHT_SQUARED,
                Colorizer::LIQUID};
    slots.gnomonic_hemisphere = GnomonicHemispherePolicy::FOLDED;
    WarpStageParams outer_warp;
    outer_warp.strength = -0.176f;
    outer_warp.time_scale = -0.00325f;
    outer_warp.frequency = 1.408f;
    outer_warp.field_angle = 2.2305307f;
    WarpStageParams inner_warp;
    inner_warp.rotation = 0.0f;
    inner_warp.cell_x = 1.0f;
    inner_warp.cell_y = 1.0f;
    inner_warp.offset_x = 0.0f;
    inner_warp.offset_y = 0.0f;
    Params params = authored_params(
        {6.3287f, 0.04f, 1.704f, 0.0f, 0.8f, 0.027f}, outer_warp,
        {2.311f, 0.0f}, {1.0f}, {0.15f, 0.0f, 0.721f, 0.854f}, {1.0f});
    params.warp.inner = inner_warp;
    return {slots, params};
  }

  static constexpr std::array<Preset, 23> PRESETS = {{
      {KALEIDOSCOPE_LIQUID_STEREO_SLOTS,
       authored_params({1.0f, 0.075f, 0.009122372f, 1.0f, 1.146f},
                       {50.749298f, 30.0f, 0.4699f}, {1.5482996f, 0.020879198f},
                       {1.0f}, {0.25410002f, 0.00015458837f, 0.201f, 0.847f},
                       {0.0030917525f})},
      {LIQUID_STEREO_SLOTS,
       authored_params({5.0f, 0.1f, 0.5f, 0.0f, 0.8f}, {3.0f, 0.5f, 0.5f},
                       {1.4f, 0.0f}, {1.0f}, {0.15f, 0.05f, 0.0f, 0.0f},
                       {1.0f})},
      {LIQUID_STEREO_SLOTS,
       authored_params({1.2f, 0.05f, 3.0f, 0.0f, 0.8f}, {3.0f, 0.5f, 0.5f},
                       {1.4f, 0.0f}, {1.0f}, {0.15f, 0.05f, 0.0f, 0.0f},
                       {1.0f})},
      {LIQUID_STEREO_SLOTS,
       authored_params({14.528f, 0.1f, 0.5f, 0.0f, 0.8f}, {3.0f, 1.479f, 0.5f},
                       {1.0f, 0.0f}, {1.0f}, {0.15f, 0.05f, 0.0f, 0.0f},
                       {1.0f})},
      {LIQUID_STEREO_SLOTS,
       authored_params({15.763f, 0.1f, 2.950552f, 0.0f, 0.8f},
                       {3.0f, 0.0f, 0.5f}, {1.0f, 0.0f}, {0.0f},
                       {0.15f, 0.05f, 0.0f, 0.0f}, {1.0f})},
      {LIQUID_STEREO_SLOTS,
       authored_params({3.28f, 0.1f, 2.463f, 0.0f, 0.8f}, {0.1f, 13.47f, 0.5f},
                       {1.209f, 0.03725f}, {0.0f},
                       {0.19710001f, 0.02f, 0.011f, 0.0f}, {0.252f})},
      {LIQUID_STEREO_SLOTS,
       authored_params(
           {15.972f, 0.1f, 2.7558982f, 0.536f, 0.8f},
           {1.8421826f, 5.377862f, 0.5f}, {1.0834427f, 0.014871964f}, {0.0f},
           {0.16880456f, 0.038022578f, 0.004391721f, 0.0f}, {0.70136297f})},
      {LIQUID_STEREO_SLOTS,
       authored_params({2.7f, 0.586f, 0.0f, 1.0f, 0.7f},
                       {47.752f, 11.55f, 0.3f}, {1.55f, ORBIT_SPIN_RATE},
                       {0.0f}, {0.0f, 0.0f, 0.097f, 1.0f}, {0.0f})},
      {LIQUID_STEREO_SLOTS,
       authored_params({8.0f, 0.30f, 0.0f, 1.0f, 0.7f}, {1.5f, 0.5f, 0.3f},
                       {2.0f, ORBIT_SPIN_RATE}, {0.0f},
                       {0.0f, 0.0f, 0.15f, 1.0f}, {0.0f})},
      {LIQUID_STEREO_SLOTS,
       authored_params({13.430163f, 0.385f, 0.0f, 1.0f, 0.7f},
                       {0.1f, 0.0f, 0.3f}, {1.627f, ORBIT_SPIN_RATE}, {1.0f},
                       {0.0f, 0.0f, 0.10404046f, 1.0f}, {0.0f})},
      {LIQUID_STEREO_SLOTS,
       authored_params({1.0f, 0.075f, 0.009122372f, 1.0f, 1.146f},
                       {38.761299f, 30.0f, 0.4699f}, {1.5482996f, 0.020879198f},
                       {0.0f}, {0.25410002f, 0.00015458837f, 0.201f, 0.847f},
                       {0.0030917525f})},
      wave_shear_liquid_preset(),
      kaleidoscope_mirror_preset(),
      gnomonic_grid_mirror_preset(SurfaceLens::KALEIDOSCOPE),
      gnomonic_grid_mirror_preset(SurfaceLens::GLITCH),
      bonne_lattice_mirror_preset(),
      peirce_lattice_preset(),
      kaleidoscope_edge_fade_liquid_preset(),
      dodecahedral_grid_preset(),
      peirce_dodecahedral_liquid_preset(),
      dodecahedral_noise_liquid_preset(),
      dodecahedral_lattice_noise_preset(),
      gnomonic_wave_shear_grid_preset(),
  }};
  static_assert(
      [] {
        for (const Preset &preset : PRESETS)
          if (!valid_config(preset))
            return false;
        return true;
      }(),
      "a ShaderBall preset lies outside its registered range");
  static_assert(
      [] {
        for (size_t index = 0; index < PRESETS.size(); ++index)
          if (!transition_admitted(PRESETS[index],
                                   PRESETS[(index + 1) % PRESETS.size()]))
            return false;
        return true;
      }(),
      "a ShaderBall preset edge lacks continuous transition admission");

  static constexpr std::array<Choreo, PRESETS.size()> CHOREO = [] {
    std::array<Choreo, PRESETS.size()> choreo;
    for (Choreo &entry : choreo)
      entry = {0, 0, 480, false};
    choreo[0] = {30, 90, 60, true};
    choreo[1] = {30, 90, 60, true};
    choreo[2] = {30, 90, 60, true};
    choreo[3] = {30, 90, 60, true};
    choreo[4] = {30, 90, 480, false};
    return choreo;
  }();
  static_assert(CHOREO.size() == PRESETS.size());

  Orientation<> projection_walk;
  Orientation<> outer_walk;
  Timeline timeline;
  size_t prepared_noise_count = 0;
  StateBundle *state = nullptr;

  Quaternion base_orientation =
      make_rotation(Vector(0, 0, -1), Vector(0, -1, 0));
  Quaternion projection_walk_prev;
  Quaternion outer_walk_prev;

  float liquid_rotation = 0.0f;
  uint32_t palette_hue = 0;
  PaletteCycler liquid_palette_cycler;
  PaletteCycler generated_palette_cycler;

  Slots active_slots = PRESETS[0].slots;
#if HS_ENABLE_PARAM_GUI_BRIDGE
  Config display_config = PRESETS[0];
  std::array<PendingEdit, PARAM_CAPACITY> pending_edits{};
  size_t pending_edit_count = 0;
  mutable std::array<char, 1024> warning_text{};
#endif
  RequestedConfig requested_config = PRESETS[0];
  Config published_config = PRESETS[0];
  Config accepted_config = PRESETS[0];
  bool requested_schema_bound = false;
  uint16_t preset_dwell_remaining = 0;
  bool preset_dwell_armed = false;
  Blend blend{PRESETS[0].params};
  LookRuntime runtime;
#if HS_ENABLE_TEST_HOOKS
  uint32_t walk_step_count = 0;
  uint32_t liquid_palette_step_count = 0;
  uint32_t generated_palette_step_count = 0;
#endif

  static constexpr size_t FOOTPRINT_BYTES =
      gamut_lut_bytes(GAMUT_ANGLE_STEPS, GAMUT_L_STEPS) +
      2 * PaletteCycler::generated_arena_bytes() +
      PARAM_CAPACITY * sizeof(ParamDef) + sizeof(StateBundle) +
      alignof(StateBundle);
  static_assert(
      FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
      "ShaderBall persistent footprint exceeds the default partition");
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(ShaderBall)
