/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

/**
 * @file ShadierBall.h
 * @brief Typed pullback sphere shader spanning ShadierBall and ShaderBall.
 */

#include "core/color/effect_palette_recipes.h"
#include "core/engine/engine.h"
#include "effects/shadierball_projections.h"

namespace hs_test {
namespace shadierball_tests {
struct ShadierBallWhiteBox;
} // namespace shadierball_tests
} // namespace hs_test

/**
 * @brief Slot-based sphere shader with an immutable per-frame pullback state.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H> class ShadierBall : public Effect {
public:
  static constexpr int GAMUT_ANGLE_STEPS = GAMUT_LUT_ANGLE_STEPS;
  static constexpr int GAMUT_L_STEPS = GAMUT_LUT_L_STEPS;

  HS_COLD_MEMBER ShadierBall() : Effect(W, H, {.strobe = true}) {}

  /** @brief Initializes slots, clocks, palette resources, and choreography. */
  HS_COLD_MEMBER void init() override {
#if defined(__EMSCRIPTEN__) || defined(HS_TEST_BUILD)
    set_parameter_updated_hook(&ShadierBall::dispatch_parameter_updated);
#endif
    use_parameter_storage(persistent_arena.allocate_n<ParamDef>(PARAM_CAPACITY),
                          PARAM_CAPACITY);
    requested_config = PRESETS[0];
    published_config = PRESETS[0];
    prepare_resource_union(PRESETS[0], PRESETS[0]);

    rebind_parameters();

    timeline.add(0, Animation::RandomWalk<W>(projection_walk, UP,
                                             projection_walk_noise));
    timeline.add(0, Animation::RandomWalk<W>(outer_walk, UP, outer_walk_noise));

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
      HS_PROFILE(sdb_timeline_step);
      timeline.step(canvas);
    }

    apply_requested_config();
    prepare_param_morph();
    const WalkDeltas walk_deltas = sample_walk_deltas();
    if (transition.active) {
      advance_runtime(transition.from_runtime, transition.from_config,
                      walk_deltas);
      advance_runtime(transition.to_runtime, transition.to_config, walk_deltas);
    } else {
      advance_runtime(runtime, {active_slots, blend.params}, walk_deltas);
    }
    liquid_palette_cycler.step();
    generated_palette_cycler.step();
#if defined(HS_TEST_BUILD)
    ++liquid_palette_step_count;
    ++generated_palette_step_count;
#endif

    if (transition.active) {
      if (transition.mode == TransitionMode::THROUGH_CLEAR)
        draw_through_clear_transition(canvas);
      else
        draw_dual_output_transition(canvas);
    } else {
      const FrameState frame = prepare_frame();
      auto shader = [&](const Vector &view)
                        HS_O3_FN -> Color4 { return shade(view, frame); };
      HS_PROFILE(sdb_shader_draw);
      Scan::Shader::draw<W, H, 1>(canvas, shader);
    }
    finish_transitions();
    publish_live_config();
  }

#if defined(HS_PROFILE_ENABLE) || defined(HS_TEST_BUILD)
  void profile_select_preset(size_t index) {
    HS_CHECK(index < PRESETS.size(),
             "ShadierBall profile preset index out of range");
    setAnimationsPaused(true);
    param_morph.active = false;
    transition.active = false;
    pending_transition.active = false;
    preset_index = index;
    active_slots = PRESETS[index].slots;
    blend.params = PRESETS[index].params;
    requested_config = PRESETS[index];
    published_config = PRESETS[index];
    runtime = {};
    HS_CHECK(prepare_resource_union(PRESETS[index], PRESETS[index]),
             "ShadierBall profile preset resources exceed capacity");
    rebind_parameters();
    hs::log("Profile preset: %u/%u", static_cast<unsigned>(index),
            static_cast<unsigned>(PRESETS.size()));
  }
#endif

private:
  friend struct ::hs_test::shadierball_tests::ShadierBallWhiteBox;

  enum class Function : uint8_t {
    TWIN_WAVE,
    RINGS,
    SPIRAL,
    GRID,
    COUPLED_DIRECT,
    NOISE_CONTOUR,
    PRIMITIVE_LATTICE
  };
  enum class Projection : uint8_t {
    EQUIRECTANGULAR,
    STEREOGRAPHIC,
    GNOMONIC,
    BONNE,
    PEIRCE_QUINCUNCIAL,
    AIROCEAN
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
    TANGENT_NOISE
  };
  enum class NoiseBasis : uint8_t { SIMPLEX, FBM3, RIDGED3 };
  enum class WarpEnvelope : uint8_t { FLAT, PROJECTION_WEIGHT, EDGE_FADE };
  enum class PolarMode : uint8_t { LINEAR, LOGARITHMIC };
  enum class CurlIntegrator : uint8_t { EULER_1, MIDPOINT_2, MIDPOINT_4 };
  enum class PolarHarmonic : uint8_t {
    H1,
    H2,
    H3,
    H4,
    H5,
    H6,
    H7,
    H8,
    H9,
    H10,
    H11,
    H12,
    H13,
    H14,
    H15,
    H16
  };
  enum class BandCount : uint8_t {
    C1,
    C2,
    C3,
    C4,
    C5,
    C6,
    C7,
    C8,
    C9,
    C10,
    C11,
    C12,
    C13,
    C14,
    C15,
    C16,
    C17,
    C18,
    C19,
    C20,
    C21,
    C22,
    C23,
    C24,
    C25,
    C26,
    C27,
    C28,
    C29,
    C30,
    C31,
    C32
  };
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
    PolarHarmonic polar_harmonic = PolarHarmonic::H1;
    int32_t seed = 1337;
    uint8_t resource_id = 0;

    bool operator==(const WarpStageSpec &) const = default;
  };
  struct WarpProgram {
    WarpStageSpec outer;
    WarpStageSpec inner;

    bool operator==(const WarpProgram &) const = default;
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
    EDGE_FADE
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

    bool operator==(const Slots &) const = default;
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

    bool operator==(const SourceParams &) const = default;

    void lerp(const SourceParams &a, const SourceParams &b, float t) {
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

    bool operator==(const WarpStageParams &) const = default;

    void lerp(const WarpStageParams &a, const WarpStageParams &b, float t) {
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

    static float lerp_angle(float a, float b, float t) {
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

    bool operator==(const WarpParams &) const = default;

    void lerp(const WarpParams &a, const WarpParams &b, float t) {
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

    bool operator==(const ProjectionParams &) const = default;

    static float lerp_periodic(float a, float b, float t) {
      if (t == 0.0f)
        return a;
      if (t == 1.0f)
        return b;
      float delta = fmodf(b - a + 0.5f, 1.0f);
      if (delta < 0.0f)
        delta += 1.0f;
      return a + (delta - 0.5f) * t;
    }

    void lerp(const ProjectionParams &a, const ProjectionParams &b, float t) {
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

    bool operator==(const SurfaceLensParams &other) const {
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

    void lerp(const SurfaceLensParams &a, const SurfaceLensParams &b, float t) {
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
    BandCount band_count = BandCount::C4;
    float band_phase = 0.0f;
    float cutout_threshold = 0.5f;
    float cutout_softness = 0.05f;
    float edge_width = 0.1f;

    bool operator==(const ValueParams &) const = default;
    void lerp(const ValueParams &a, const ValueParams &b, float t) {
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

    bool operator==(const ColorizerParams &) const = default;

    void lerp(const ColorizerParams &a, const ColorizerParams &b, float t) {
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

    bool operator==(const OuterCameraParams &) const = default;

    void lerp(const OuterCameraParams &a, const OuterCameraParams &b, float t) {
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

    bool operator==(const Params &) const = default;

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

    static float phase_t(float t, int phase, int phase_count) {
      return ease_in_out_sin(hs::clamp(t * phase_count - phase, 0.0f, 1.0f));
    }
  };

  struct Config {
    Slots slots;
    Params params;

    bool operator==(const Config &) const = default;
  };
  using RequestedConfig = Config;
  using Preset = Config;

  HS_COLD_MEMBER void rebind_parameters() {
    reset_parameters();
    Slots &slots = requested_config.slots;
    register_animated_param("Function", &slots.function, FUNCTION_OPTIONS,
                            FUNCTION_EXPORT_OPTIONS, NUM_FUNCTIONS);
    register_animated_param("Projection", &slots.projection, PROJECTION_OPTIONS,
                            PROJECTION_EXPORT_OPTIONS, NUM_PROJECTIONS);
    if (slots.projection == Projection::PEIRCE_QUINCUNCIAL)
      register_animated_param("Peirce Layout", &slots.peirce_layout,
                              PEIRCE_LAYOUT_OPTIONS,
                              PEIRCE_LAYOUT_EXPORT_OPTIONS, NUM_PEIRCE_LAYOUTS);
    if (slots.projection == Projection::AIROCEAN)
      register_animated_param(
          "Airocean Layout", &slots.airocean_layout, AIROCEAN_LAYOUT_OPTIONS,
          AIROCEAN_LAYOUT_EXPORT_OPTIONS, NUM_AIROCEAN_LAYOUTS);
    register_animated_param(
        "Projection Frame", &slots.projection_frame, PROJECTION_FRAME_OPTIONS,
        PROJECTION_FRAME_EXPORT_OPTIONS, NUM_PROJECTION_FRAMES);
    register_animated_param("Lens", &slots.surface_lens, LENS_OPTIONS,
                            LENS_EXPORT_OPTIONS, NUM_LENSES);
    register_animated_param("Outer Warp", &slots.warp_program.outer.kind,
                            WARP_OPTIONS, WARP_EXPORT_OPTIONS, NUM_WARPS);
    register_animated_param("Inner Warp", &slots.warp_program.inner.kind,
                            WARP_OPTIONS, WARP_EXPORT_OPTIONS, NUM_WARPS);
    register_stage_slot_controls("Outer", slots.warp_program.outer);
    register_stage_slot_controls("Inner", slots.warp_program.inner);
    register_animated_param("Signal Weight", &slots.signal_weight,
                            SIGNAL_OPTIONS, SIGNAL_EXPORT_OPTIONS, NUM_SIGNALS);
    register_animated_param("Value Transfer", &slots.value_transfer,
                            VALUE_TRANSFER_OPTIONS,
                            VALUE_TRANSFER_EXPORT_OPTIONS, NUM_VALUE_TRANSFERS);
    register_animated_param("Coverage", &slots.coverage, COVERAGE_OPTIONS,
                            COVERAGE_EXPORT_OPTIONS, NUM_COVERAGE_POLICIES);
    register_animated_param("Colorizer", &slots.colorizer, COLORIZER_OPTIONS,
                            COLORIZER_EXPORT_OPTIONS, NUM_COLORIZERS);
    register_active_parameter_controls(slots, requested_config.params);
    requested_schema_bound = true;
  }

#if defined(__EMSCRIPTEN__) || defined(HS_TEST_BUILD)
  enum class EditIntent : uint8_t {
    NONE,
    FUNCTION,
    PROJECTION,
    LENS,
    OUTER_WARP,
    INNER_WARP,
    OUTER_POLAR_HARMONIC,
    INNER_POLAR_HARMONIC,
  };

  static void dispatch_parameter_updated(Effect *effect, const char *name,
                                         bool is_enum) {
    static_cast<ShadierBall *>(effect)->parameter_updated(name, is_enum);
  }

  void parameter_updated(const char *name, bool is_enum) {
    if (!is_enum) {
      if (std::strncmp(name, "Mobius ", 7) == 0)
        canonicalize_mobius(requested_config.params.surface_lens.mobius);
      return;
    }
    const EditIntent intent = edit_intent(name);
    if (intent != EditIntent::NONE)
      canonicalize_selector_edit(requested_config, intent);
    if (schema_selector(name))
      rebind_parameters();
  }

  static EditIntent edit_intent(const char *name) {
    if (std::strcmp(name, "Function") == 0)
      return EditIntent::FUNCTION;
    if (std::strcmp(name, "Projection") == 0)
      return EditIntent::PROJECTION;
    if (std::strcmp(name, "Lens") == 0)
      return EditIntent::LENS;
    if (std::strcmp(name, "Outer Warp") == 0)
      return EditIntent::OUTER_WARP;
    if (std::strcmp(name, "Inner Warp") == 0)
      return EditIntent::INNER_WARP;
    if (std::strcmp(name, "Outer Polar Harmonic") == 0)
      return EditIntent::OUTER_POLAR_HARMONIC;
    if (std::strcmp(name, "Inner Polar Harmonic") == 0)
      return EditIntent::INNER_POLAR_HARMONIC;
    return EditIntent::NONE;
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
           std::strcmp(name, "Outer Warp") == 0 ||
           std::strcmp(name, "Inner Warp") == 0 ||
           std::strcmp(name, "Value Transfer") == 0 ||
           std::strcmp(name, "Coverage") == 0 ||
           std::strcmp(name, "Colorizer") == 0;
  }

  static void clear_strict_seam_incompatible(Config &config) {
    if (config.slots.function == Function::NOISE_CONTOUR)
      config.slots.function = Function::COUPLED_DIRECT;
    if (config.slots.surface_lens == SurfaceLens::TANGENT_NOISE)
      config.slots.surface_lens = SurfaceLens::NONE;
    if (seam_sensitive_warp(config.slots.warp_program.outer.kind))
      config.slots.warp_program.outer.kind = WarpStageKind::NONE;
    if (seam_sensitive_warp(config.slots.warp_program.inner.kind))
      config.slots.warp_program.inner.kind = WarpStageKind::NONE;
  }

  static void clear_incompatible_polar_stages(Config &config) {
    if (config.slots.warp_program.outer.kind == WarpStageKind::POLAR_CHART &&
        !polar_source_compatible(config, config.slots.warp_program.outer))
      config.slots.warp_program.outer.kind = WarpStageKind::NONE;
    if (config.slots.warp_program.inner.kind == WarpStageKind::POLAR_CHART &&
        !polar_source_compatible(config, config.slots.warp_program.inner))
      config.slots.warp_program.inner.kind = WarpStageKind::NONE;
  }

  static void fit_selected_warp_to_device_budget(Config &config,
                                                 bool outer_selected) {
    if (within_hold_device_budget(device_cost(config)))
      return;
    WarpStageSpec &sibling = outer_selected ? config.slots.warp_program.inner
                                            : config.slots.warp_program.outer;
    sibling.kind = WarpStageKind::NONE;
  }

  static void canonicalize_selector_edit(Config &config, EditIntent intent) {
    WarpStageSpec &outer = config.slots.warp_program.outer;
    WarpStageSpec &inner = config.slots.warp_program.inner;
    switch (intent) {
    case EditIntent::FUNCTION:
      if (strict_projection(config.slots.projection) &&
          config.slots.function == Function::NOISE_CONTOUR)
        config.slots.projection = Projection::EQUIRECTANGULAR;
      clear_incompatible_polar_stages(config);
      break;
    case EditIntent::PROJECTION:
      if (config.slots.projection != Projection::STEREOGRAPHIC) {
        if (outer.kind == WarpStageKind::LEGACY_STEREO_NOISE)
          outer.kind = WarpStageKind::NONE;
        if (inner.kind == WarpStageKind::LEGACY_STEREO_NOISE)
          inner.kind = WarpStageKind::NONE;
      }
      if (strict_projection(config.slots.projection))
        clear_strict_seam_incompatible(config);
      break;
    case EditIntent::LENS:
      if (strict_projection(config.slots.projection) &&
          config.slots.surface_lens == SurfaceLens::TANGENT_NOISE)
        config.slots.projection = Projection::EQUIRECTANGULAR;
      break;
    case EditIntent::OUTER_WARP:
      if (outer.kind == WarpStageKind::LEGACY_STEREO_NOISE) {
        config.slots.projection = Projection::STEREOGRAPHIC;
        if (inner.kind == WarpStageKind::LEGACY_STEREO_NOISE)
          inner.kind = WarpStageKind::NONE;
      } else if (seam_sensitive_warp(outer.kind) &&
                 strict_projection(config.slots.projection)) {
        config.slots.projection = Projection::EQUIRECTANGULAR;
      }
      if (outer.kind == WarpStageKind::POLAR_CHART) {
        inner.kind = WarpStageKind::NONE;
        if (!polar_source_compatible(config, outer)) {
          config.slots.function = Function::COUPLED_DIRECT;
          config.params.source.pattern_freq = 1.0f;
        }
      }
      if (warp_uses_noise(outer.kind))
        outer.resource_id = 0;
      canonicalize_warp_params(outer, config.params.warp.outer);
      fit_selected_warp_to_device_budget(config, true);
      break;
    case EditIntent::INNER_WARP:
      if (inner.kind == WarpStageKind::LEGACY_STEREO_NOISE) {
        config.slots.projection = Projection::STEREOGRAPHIC;
        if (outer.kind == WarpStageKind::LEGACY_STEREO_NOISE)
          outer.kind = WarpStageKind::NONE;
      } else if (seam_sensitive_warp(inner.kind) &&
                 strict_projection(config.slots.projection)) {
        config.slots.projection = Projection::EQUIRECTANGULAR;
      }
      if (inner.kind == WarpStageKind::POLAR_CHART) {
        if (outer.kind == WarpStageKind::POLAR_CHART)
          outer.kind = WarpStageKind::NONE;
        if (!polar_source_compatible(config, inner)) {
          config.slots.function = Function::COUPLED_DIRECT;
          config.params.source.pattern_freq = 1.0f;
        }
      }
      if (warp_uses_noise(inner.kind))
        inner.resource_id = 1;
      canonicalize_warp_params(inner, config.params.warp.inner);
      fit_selected_warp_to_device_budget(config, false);
      break;
    case EditIntent::OUTER_POLAR_HARMONIC:
      if (outer.kind == WarpStageKind::POLAR_CHART &&
          !polar_source_compatible(config, outer))
        config.params.source.pattern_freq = 1.0f;
      break;
    case EditIntent::INNER_POLAR_HARMONIC:
      if (inner.kind == WarpStageKind::POLAR_CHART &&
          !polar_source_compatible(config, inner))
        config.params.source.pattern_freq = 1.0f;
      break;
    case EditIntent::NONE:
      break;
    }
  }
#endif

  HS_COLD_MEMBER void register_active_parameter_controls(Slots &slots,
                                                         Params &params) {
    const bool polar_topology =
        slots.warp_program.outer.kind == WarpStageKind::POLAR_CHART ||
        slots.warp_program.inner.kind == WarpStageKind::POLAR_CHART;
    register_source_controls(slots.function, params.source, !polar_topology);
    register_projection_controls(slots, params);
    register_lens_controls(slots.surface_lens, params.surface_lens);
    register_active_warp_controls("Outer", slots.warp_program.outer,
                                  params.warp.outer);
    register_active_warp_controls("Inner", slots.warp_program.inner,
                                  params.warp.inner);
    if (slots.value_transfer == ValueTransfer::ISO_CONTOUR) {
      register_animated_param("Iso Level", &params.value.iso_level, 0.0f, 1.0f);
      register_animated_param("Iso Width", &params.value.iso_width,
                              SOFTNESS_MIN, 0.5f);
    } else if (slots.value_transfer == ValueTransfer::SMOOTH_BANDS) {
      register_animated_param("Band Count", &params.value.band_count,
                              BAND_COUNT_OPTIONS, BAND_COUNT_EXPORT_OPTIONS,
                              NUM_BAND_COUNTS);
      register_animated_param("Band Phase", &params.value.band_phase, 0.0f,
                              TWO_PI_F);
    }
    if (slots.coverage == CoveragePolicy::VALUE_CUTOUT) {
      register_animated_param("Cutout Threshold",
                              &params.value.cutout_threshold, 0.0f, 1.0f);
      register_animated_param("Cutout Softness", &params.value.cutout_softness,
                              SOFTNESS_MIN, 0.5f);
    } else if (slots.coverage == CoveragePolicy::EDGE_FADE) {
      register_animated_param("Edge Fade Width", &params.value.edge_width,
                              SOFTNESS_MIN, 0.5f);
    }
    register_colorizer_controls(slots.colorizer, params.colorizer);
    register_animated_param("Outer Wander", &params.outer_camera.wander,
                            WANDER_MIN, WANDER_MAX);
  }

  HS_COLD_MEMBER void register_stage_slot_controls(const char *position,
                                                   WarpStageSpec &spec) {
    const bool outer = position[0] == 'O';
    if (spec.kind == WarpStageKind::VECTOR_NOISE ||
        spec.kind == WarpStageKind::CURL_FLOW) {
      register_animated_param(outer ? "Outer Noise Basis" : "Inner Noise Basis",
                              &spec.basis, NOISE_BASIS_OPTIONS,
                              NOISE_BASIS_EXPORT_OPTIONS, NUM_NOISE_BASES);
      register_animated_param(outer ? "Outer Warp Envelope"
                                    : "Inner Warp Envelope",
                              &spec.envelope, WARP_ENVELOPE_OPTIONS,
                              WARP_ENVELOPE_EXPORT_OPTIONS, NUM_WARP_ENVELOPES);
    }
    if (spec.kind == WarpStageKind::CURL_FLOW)
      register_animated_param(
          outer ? "Outer Curl Integrator" : "Inner Curl Integrator",
          &spec.curl_integrator, CURL_INTEGRATOR_OPTIONS,
          CURL_INTEGRATOR_EXPORT_OPTIONS, NUM_CURL_INTEGRATORS);
    if (spec.kind == WarpStageKind::POLAR_CHART) {
      register_animated_param(outer ? "Outer Polar Mode" : "Inner Polar Mode",
                              &spec.polar_mode, POLAR_MODE_OPTIONS,
                              POLAR_MODE_EXPORT_OPTIONS, NUM_POLAR_MODES);
      register_animated_param(
          outer ? "Outer Polar Harmonic" : "Inner Polar Harmonic",
          &spec.polar_harmonic, POLAR_HARMONIC_OPTIONS,
          POLAR_HARMONIC_EXPORT_OPTIONS, NUM_POLAR_HARMONICS);
    }
  }

  HS_COLD_MEMBER void register_source_controls(Function function,
                                               SourceParams &params,
                                               bool expose_pattern_frequency) {
    if (function == Function::NOISE_CONTOUR) {
      register_animated_param("Source Noise Scale", &params.noise_scale,
                              NOISE_SCALE_MIN, NOISE_SCALE_MAX);
      register_animated_param("Source Noise Contrast", &params.noise_contrast,
                              0.0f, 8.0f);
      register_animated_param("Source Noise Rate", &params.noise_time_rate,
                              NOISE_RATE_MIN, NOISE_RATE_MAX);
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
                              SOFTNESS_MIN, 0.25f);
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
    if (function == Function::TWIN_WAVE ||
        function == Function::COUPLED_DIRECT) {
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
    if (slots.projection == Projection::EQUIRECTANGULAR ||
        slots.projection == Projection::STEREOGRAPHIC ||
        slots.projection == Projection::GNOMONIC)
      register_animated_param("Pole Fade", &params.projection.pole_fade,
                              POLE_FADE_MIN, POLE_FADE_MAX);
    if (slots.projection == Projection::EQUIRECTANGULAR ||
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
    if (slots.projection_frame == ProjectionFramePolicy::SPIN_WANDER) {
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
                              NOISE_SCALE_MIN, NOISE_SCALE_MAX);
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

  HS_COLD_MEMBER void register_active_warp_controls(const char *position,
                                                    const WarpStageSpec &spec,
                                                    WarpStageParams &params) {
    if (spec.kind == WarpStageKind::NONE)
      return;
    const bool outer = position[0] == 'O';
    const char *const *names =
        outer ? OUTER_WARP_PARAM_NAMES : INNER_WARP_PARAM_NAMES;
    const char *time_name = outer ? "Outer Warp Time" : "Inner Warp Time";
    if (spec.kind == WarpStageKind::LEGACY_STEREO_NOISE ||
        spec.kind == WarpStageKind::WAVE_SHEAR ||
        spec.kind == WarpStageKind::VECTOR_NOISE ||
        spec.kind == WarpStageKind::CURL_FLOW) {
      const char *strength_name =
          outer ? "Outer Warp Strength" : "Inner Warp Strength";
      const bool signed_strength = spec.kind == WarpStageKind::WAVE_SHEAR ||
                                   spec.kind == WarpStageKind::CURL_FLOW;
      const float strength_max =
          spec.kind == WarpStageKind::LEGACY_STEREO_NOISE ? 30.0f : 4.0f;
      register_animated_param(strength_name, &params.strength,
                              signed_strength ? -4.0f : 0.0f, strength_max);
    }
    if (spec.kind == WarpStageKind::LEGACY_STEREO_NOISE) {
      register_animated_param(outer ? "Outer Warp Scale" : "Inner Warp Scale",
                              &params.scale, 0.1f, 100.0f);
      register_animated_param(time_name, &params.time_scale, 0.05f, 1.0f);
      return;
    }
    if (spec.kind == WarpStageKind::WAVE_SHEAR ||
        spec.kind == WarpStageKind::VORTEX ||
        spec.kind == WarpStageKind::VECTOR_NOISE ||
        spec.kind == WarpStageKind::CURL_FLOW)
      register_animated_param(time_name, &params.time_scale, NOISE_RATE_MIN,
                              NOISE_RATE_MAX);
    switch (spec.kind) {
    case WarpStageKind::AFFINE_FRAME: {
      for (int index = 0; index < 6; ++index) {
        float *targets[] = {&params.translation_x, &params.translation_y,
                            &params.rotation,      &params.scale_x,
                            &params.scale_y,       &params.shear};
        const float minimum[] = {-4.0f, -4.0f, 0.0f, 0.25f, 0.25f, -0.75f};
        const float maximum[] = {4.0f, 4.0f, TWO_PI_F, 4.0f, 4.0f, 0.75f};
        register_animated_param(names[index], targets[index], minimum[index],
                                maximum[index]);
      }
      break;
    }
    case WarpStageKind::WAVE_SHEAR:
      register_animated_param(names[6], &params.frequency, 0.0f, 64.0f);
      register_animated_param(names[7], &params.field_angle, 0.0f, TWO_PI_F);
      break;
    case WarpStageKind::VORTEX:
      register_animated_param(names[8], &params.center_x, -4.0f, 4.0f);
      register_animated_param(names[9], &params.center_y, -4.0f, 4.0f);
      register_animated_param(names[10], &params.radius, 1.0f / 64.0f, 8.0f);
      register_animated_param(names[11], &params.turns, -4.0f, 4.0f);
      register_animated_param(names[21], &params.center_orbit_radius, 0.0f,
                              4.0f);
      break;
    case WarpStageKind::VECTOR_NOISE:
    case WarpStageKind::CURL_FLOW:
      register_animated_param(
          outer ? "Outer Warp Scale" : "Inner Warp Scale", &params.scale,
          1.0f / 64.0f, spec.kind == WarpStageKind::CURL_FLOW ? 16.0f : 64.0f);
      register_animated_param(names[12], &params.vector_angle, 0.0f, TWO_PI_F);
      register_animated_param(names[20], &params.edge_width, SOFTNESS_MIN,
                              0.5f);
      break;
    case WarpStageKind::MIRROR_TILE:
      register_animated_param(names[2], &params.rotation, 0.0f, TWO_PI_F);
      register_animated_param(names[13], &params.cell_x, CELL_MIN, CELL_MAX);
      register_animated_param(names[14], &params.cell_y, CELL_MIN, CELL_MAX);
      register_animated_param(names[15], &params.offset_x, -8.0f, 8.0f);
      register_animated_param(names[16], &params.offset_y, -8.0f, 8.0f);
      break;
    case WarpStageKind::POLAR_CHART:
      register_animated_param(names[17], &params.radial_scale, 1.0f / 64.0f,
                              16.0f);
      register_animated_param(names[18], &params.radial_phase, 0.0f, TWO_PI_F);
      register_animated_param(names[19], &params.angular_phase, 0.0f, TWO_PI_F);
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
    uint8_t flags;
  };

  struct PlanarWarpStageResult {
    Complex coords;
    Complex delta;
    float deformation;
    float path_length;
    uint8_t flags;
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

    bool operator==(const NoiseResourceKey &) const = default;
  };

  struct FrameState {
    Slots slots;
    Params params;
    ClockState clocks;
    SourceState prepared_source;
    float breathe_offset;
    PreparedTransforms transforms;
    PreparedWarpProgram prepared_warp;
    ResourceBindings resources;
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

  enum class TransitionMode : uint8_t { DUAL_OUTPUT, THROUGH_CLEAR };

  static constexpr uint16_t NOISE_CALL_POINTS = 4;
  static constexpr uint16_t HOLD_DEVICE_POINT_BUDGET = 40;

  struct DeviceCost {
    uint16_t fixed_points = 0;
    uint16_t noise_calls = 0;
    uint8_t shader_lookups = 1;

    constexpr uint16_t worst_case_noise_calls() const {
      return noise_calls * shader_lookups;
    }

    constexpr uint16_t worst_case_points() const {
      return (fixed_points + NOISE_CALL_POINTS * noise_calls) * shader_lookups;
    }
  };

  enum class CostTier : uint8_t { T0, T1, T2, T3 };

  struct TransitionRuntime {
    Config from_config;
    Config to_config;
    LookRuntime from_runtime;
    LookRuntime to_runtime;
    uint16_t elapsed = 0;
    uint16_t duration = 0;
    TransitionMode mode = TransitionMode::DUAL_OUTPUT;
    bool continue_choreo = false;
    bool active = false;

    TransitionRuntime() = default;
    TransitionRuntime(const Config &from_config, const Config &to_config,
                      const LookRuntime &from_runtime,
                      const LookRuntime &to_runtime, uint16_t elapsed,
                      uint16_t duration, bool continue_choreo, bool active,
                      TransitionMode mode = TransitionMode::DUAL_OUTPUT)
        : from_config(from_config), to_config(to_config),
          from_runtime(from_runtime), to_runtime(to_runtime), elapsed(elapsed),
          duration(duration), mode(mode), continue_choreo(continue_choreo),
          active(active) {}
  };

  struct PendingTransition {
    bool active = false;
  };

  struct DualOutputFrame {
    FrameState from;
    FrameState to;
    float mix;
  };

  struct ThroughClearPhase {
    float alpha;
    bool from_endpoint;
    bool clear;
  };

  static constexpr bool warp_uses_noise(WarpStageKind kind) {
    return kind == WarpStageKind::LEGACY_STEREO_NOISE ||
           kind == WarpStageKind::VECTOR_NOISE ||
           kind == WarpStageKind::CURL_FLOW;
  }

  static constexpr bool seam_sensitive_warp(WarpStageKind kind) {
    return kind == WarpStageKind::VECTOR_NOISE ||
           kind == WarpStageKind::CURL_FLOW;
  }

  static constexpr NoiseResourceKey
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

  static constexpr NoiseResourceKey source_resource_key(const Config &config) {
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

  static constexpr NoiseResourceKey lens_resource_key(const Config &config) {
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

  static constexpr bool
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

  static constexpr bool append_config_resource_keys(
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

  static constexpr bool resource_union_fits(const Config &from,
                                            const Config &to) {
    std::array<NoiseResourceKey, MAX_NOISE_RESOURCES> keys{};
    size_t count = 0;
    return append_config_resource_keys(from, keys, count) &&
           append_config_resource_keys(to, keys, count);
  }

  static constexpr uint16_t scalar_noise_calls(NoiseBasis basis) {
    return basis == NoiseBasis::SIMPLEX ? 1 : 3;
  }

  static constexpr uint16_t warp_noise_calls(const WarpStageSpec &spec) {
    if (spec.kind == WarpStageKind::LEGACY_STEREO_NOISE)
      return 6;
    if (spec.kind == WarpStageKind::VECTOR_NOISE)
      return 2 * scalar_noise_calls(spec.basis) *
             (spec.basis == NoiseBasis::RIDGED3 ? 4 : 2);
    if (spec.kind != WarpStageKind::CURL_FLOW)
      return 0;
    const uint16_t gradients =
        spec.curl_integrator == CurlIntegrator::EULER_1      ? 1
        : spec.curl_integrator == CurlIntegrator::MIDPOINT_2 ? 4
                                                             : 8;
    return 2 * gradients * 4 * scalar_noise_calls(spec.basis);
  }

  static constexpr uint16_t projection_fixed_points(Projection projection) {
    return projection == Projection::BONNE                ? 12
           : projection == Projection::PEIRCE_QUINCUNCIAL ? 36
           : projection == Projection::AIROCEAN           ? 40
                                                          : 0;
  }

  static constexpr uint16_t warp_fixed_points(WarpStageKind kind) {
    return kind == WarpStageKind::AFFINE_FRAME  ? 2
           : kind == WarpStageKind::WAVE_SHEAR  ? 4
           : kind == WarpStageKind::VORTEX      ? 6
           : kind == WarpStageKind::MIRROR_TILE ? 6
           : kind == WarpStageKind::POLAR_CHART ? 8
                                                : 0;
  }

  static constexpr uint16_t lens_fixed_points(SurfaceLens lens) {
    return lens == SurfaceLens::MOBIUS ? 12 : lens == SurfaceLens::NONE ? 0 : 4;
  }

  static constexpr DeviceCost device_cost(const Config &config) {
    DeviceCost result;
    result.fixed_points =
        projection_fixed_points(config.slots.projection) +
        warp_fixed_points(config.slots.warp_program.outer.kind) +
        warp_fixed_points(config.slots.warp_program.inner.kind) +
        lens_fixed_points(config.slots.surface_lens);
    result.noise_calls = warp_noise_calls(config.slots.warp_program.outer) +
                         warp_noise_calls(config.slots.warp_program.inner);
    if (config.slots.function == Function::NOISE_CONTOUR)
      result.noise_calls +=
          2 * scalar_noise_calls(config.params.source.noise_basis);
    if (config.slots.surface_lens == SurfaceLens::TANGENT_NOISE)
      result.noise_calls +=
          4 * scalar_noise_calls(config.params.surface_lens.noise_basis);
    if (strict_projection(config.slots.projection) &&
        config.slots.surface_lens != SurfaceLens::NONE &&
        config.params.surface_lens.mix > 0.0f &&
        config.params.surface_lens.mix < 1.0f)
      result.shader_lookups = 2;
    return result;
  }

  static constexpr CostTier cost_tier(const DeviceCost &cost) {
    const uint16_t points = cost.worst_case_points();
    return points <= 12   ? CostTier::T0
           : points <= 24 ? CostTier::T1
           : points <= 32 ? CostTier::T2
                          : CostTier::T3;
  }

  static constexpr bool within_hold_device_budget(const DeviceCost &cost) {
    return cost.worst_case_points() <= HOLD_DEVICE_POINT_BUDGET;
  }

  static constexpr bool hold_admitted(const Config &config) {
    return valid_config(config) &&
           within_hold_device_budget(device_cost(config));
  }

  static constexpr TransitionMode transition_mode(const Config &from,
                                                  const Config &to) {
    (void)from;
    (void)to;
    return TransitionMode::THROUGH_CLEAR;
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
      prepared_noise_keys[index] = keys[index];
      noise_resources[index].SetNoiseType(
          FastNoiseLite::NoiseType_OpenSimplex2);
      noise_resources[index].SetSeed(keys[index].seed);
      noise_resources[index].SetFrequency(keys[index].generator_frequency);
    }
    return true;
  }

  const FastNoiseLite *resolve_resource(const NoiseResourceKey &key) const {
    for (size_t index = 0; index < prepared_noise_count; ++index)
      if (prepared_noise_keys[index] == key)
        return &noise_resources[index];
    return nullptr;
  }

  const FastNoiseLite *resolve_warp_resource(const WarpStageSpec &spec) const {
    return warp_uses_noise(spec.kind)
               ? resolve_resource(warp_resource_key(spec))
               : nullptr;
  }

  const FastNoiseLite *resolve_source_resource(const Config &config) const {
    return config.slots.function == Function::NOISE_CONTOUR
               ? resolve_resource(source_resource_key(config))
               : nullptr;
  }

  const FastNoiseLite *resolve_lens_resource(const Config &config) const {
    return config.slots.surface_lens == SurfaceLens::TANGENT_NOISE
               ? resolve_resource(lens_resource_key(config))
               : nullptr;
  }

  static WrappedNoisePhase prepare_noise_phase(float turns) {
    const float wrapped_turns = wrap_t(turns);
    const float current_time = wrapped_turns * NOISE_NATIVE_PERIOD;
    if (wrapped_turns <= NOISE_WRAP_START)
      return {current_time, current_time, 0.0f, false};
    return {current_time, current_time - NOISE_NATIVE_PERIOD,
            ease_in_out_sin((wrapped_turns - NOISE_WRAP_START) /
                            (1.0f - NOISE_WRAP_START)),
            true};
  }

  static PreparedWarpStage prepare_warp_stage(const WarpStageParams &params,
                                              float stage_phase) {
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
        {resolve_warp_resource(config.slots.warp_program.outer),
         resolve_warp_resource(config.slots.warp_program.inner),
         resolve_source_resource(config), resolve_lens_resource(config),
         &generated_palette_cycler.palette(),
         &liquid_palette_cycler.palette()}};
  }

  HS_COLD_MEMBER DualOutputFrame prepare_dual_output_frame() const {
    return {prepare_frame(transition.from_config, transition.from_runtime),
            prepare_frame(transition.to_config, transition.to_runtime),
            transition_mix(transition.elapsed, transition.duration)};
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

  HS_NOINLINE_NOCLONE void draw_through_clear_transition(Canvas &canvas) const {
    const ThroughClearPhase phase =
        through_clear_phase(transition.elapsed, transition.duration);
    if (phase.clear)
      return;
    const FrameState visible =
        phase.from_endpoint
            ? prepare_frame(transition.from_config, transition.from_runtime)
            : prepare_frame(transition.to_config, transition.to_runtime);
    auto shader = [&](const Vector &view) HS_O3_FN -> Color4 {
      return shade_through_clear(view, &visible, phase);
    };
    HS_PROFILE(sdb_shader_draw);
    Scan::Shader::draw<W, H, 1>(canvas, shader);
  }

  HS_NOINLINE_NOCLONE void draw_dual_output_transition(Canvas &canvas) const {
    const DualOutputFrame frame = prepare_dual_output_frame();
    auto shader = [&](const Vector &view) HS_O3_FN -> Color4 {
      return shade_dual_output(view, frame);
    };
    HS_PROFILE(sdb_shader_draw);
    Scan::Shader::draw<W, H, 1>(canvas, shader);
  }

  static Color4 shade(const Vector &view, const FrameState &frame) {
    const Vector outer_local = outer_camera_lookup(view, frame);
    if (strict_projection(frame.slots.projection) &&
        frame.slots.surface_lens != SurfaceLens::NONE &&
        frame.params.surface_lens.mix > 0.0f &&
        frame.params.surface_lens.mix < 1.0f) {
      const ProjectedLookup direct = project_branch(outer_local, frame);
      const ProjectedLookup lensed =
          project_branch(apply_lens(outer_local, frame), frame);
      if (!projection_join_compatible(direct, lensed, frame.slots.projection,
                                      frame.params.projection.coordinate_scale))
        return blend_outputs(shade_projected(direct, frame),
                             shade_projected(lensed, frame),
                             frame.params.surface_lens.mix);
    }
    const ProjectedLookup projected =
        surface_lens_project_lookup(outer_local, frame);
    return shade_projected(projected, frame);
  }

  static Color4 shade_projected(const ProjectedLookup &projected,
                                const FrameState &frame) {
    const PlanarWarpResult warped = planar_warp_lookup(projected, frame);
    const Complex source_coords = condition_source_coords(warped.coords, frame);
    const float field = sample_source(source_coords, frame);
    const MaterialSample material =
        shape_material(field, projected, warped, frame);
    return colorize(material, frame);
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

  static Color4 shade_dual_output(const Vector &view,
                                  const DualOutputFrame &frame) {
    if (frame.mix == 0.0f)
      return shade(view, frame.from);
    if (frame.mix == 1.0f)
      return shade(view, frame.to);
    return blend_outputs(shade(view, frame.from), shade(view, frame.to),
                         frame.mix);
  }

  static Color4 shade_through_clear(const Vector &view,
                                    const FrameState *visible,
                                    const ThroughClearPhase &phase) {
    if (phase.clear)
      return Color4();
    HS_CHECK(visible != nullptr,
             "through-clear visible phase requires an endpoint frame");
    Color4 color = shade(view, *visible);
    color.alpha *= phase.alpha;
    return color;
  }

  static Color4 blend_outputs(const Color4 &from, const Color4 &to, float mix) {
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
      return project_branch(v, frame);
    const Vector lensed = apply_lens(v, frame);
    if (mix == 1.0f)
      return project_branch(lensed, frame);
    const ProjectedLookup direct = project_branch(v, frame);
    const ProjectedLookup distorted = project_branch(lensed, frame);
    return join_projected(direct, distorted, mix, frame.slots.projection,
                          frame.params.projection.pole_fade);
  }

  static ProjectedLookup project_branch(const Vector &v,
                                        const FrameState &frame) {
    const Vector local = rotate(v, frame.transforms.projection_conj);
    if (frame.slots.projection == Projection::BONNE) {
      const shadierball::ProjectionKernelResult result =
          shadierball::bonne_projection(
              local, frame.params.projection.central_meridian,
              (frame.slots.bonne_hemisphere == BonneHemisphere::NORTH ? 1.0f
                                                                      : -1.0f) *
                  frame.params.projection.bonne_standard_parallel);
      return {{result.coords.re * frame.params.projection.coordinate_scale,
               result.coords.im * frame.params.projection.coordinate_scale},
              result.region_id,
              result.component_id,
              result.boundary_flags,
              result.fade_edge_distance *
                  fabsf(frame.params.projection.coordinate_scale),
              1.0f,
              result.flags,
              result.traits,
              result.edge_class};
    }
    if (frame.slots.projection == Projection::PEIRCE_QUINCUNCIAL) {
      const shadierball::ProjectionKernelResult result =
          shadierball::peirce_projection(
              local, frame.params.projection.central_meridian,
              static_cast<uint8_t>(frame.slots.peirce_layout),
              frame.params.projection.layout_scroll,
              projection_edge_distance_required(frame));
      return {{result.coords.re * frame.params.projection.coordinate_scale,
               result.coords.im * frame.params.projection.coordinate_scale},
              result.region_id,
              result.component_id,
              result.boundary_flags,
              result.fade_edge_distance *
                  fabsf(frame.params.projection.coordinate_scale),
              1.0f,
              result.flags,
              result.traits,
              result.edge_class};
    }
    if (frame.slots.projection == Projection::AIROCEAN) {
      const shadierball::ProjectionKernelResult result =
          shadierball::airocean_projection(
              local, frame.params.projection.central_meridian,
              frame.slots.airocean_layout == AiroceanLayout::HORIZONTAL,
              projection_edge_distance_required(frame));
      return {{result.coords.re * frame.params.projection.coordinate_scale,
               result.coords.im * frame.params.projection.coordinate_scale},
              result.region_id,
              result.component_id,
              result.boundary_flags,
              result.fade_edge_distance *
                  fabsf(frame.params.projection.coordinate_scale),
              1.0f,
              result.flags,
              result.traits,
              result.edge_class};
    }
    const Complex coords =
        frame.slots.projection == Projection::EQUIRECTANGULAR
            ? equirectangular(local, frame.params.projection.central_meridian)
            : project_point(local, frame.slots.projection);
    return finalize_projection(local, coords, frame.slots.projection,
                               frame.params.projection.pole_fade,
                               frame.slots.gnomonic_hemisphere);
  }

  static ProjectedLookup
  finalize_projection(const Vector &local, const Complex &coords,
                      Projection projection, float pole_fade,
                      GnomonicHemispherePolicy gnomonic_hemisphere =
                          GnomonicHemispherePolicy::FOLDED) {
    const float r_sq = coords.re * coords.re + coords.im * coords.im;
    switch (projection) {
    case Projection::EQUIRECTANGULAR:
      return {coords,
              static_cast<uint8_t>(local.z < 0.0f),
              0,
              0,
              1.0f,
              pole_attenuation(r_sq, pole_fade),
              PROJECTION_FLAG_FOLDED};
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

  static ProjectedLookup join_projected(const ProjectedLookup &direct,
                                        const ProjectedLookup &lensed,
                                        float mix, Projection projection,
                                        float pole_fade) {
    if (mix == 0.0f)
      return direct;
    if (mix == 1.0f)
      return lensed;
    const Complex coords(hs::lerp(direct.coords.re, lensed.coords.re, mix),
                         hs::lerp(direct.coords.im, lensed.coords.im, mix));
    const ProjectedLookup *selected = nullptr;
    switch (projection) {
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

  static PlanarWarpResult planar_warp_lookup(const ProjectedLookup &projected,
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
            outer.path_length + inner.path_length,
            static_cast<uint8_t>(outer.flags | inner.flags)};
  }

  static PlanarWarpStageResult
  warp_stage_lookup(const Complex &input, const ProjectedLookup &projected,
                    const WarpStageSpec &spec, const WarpStageParams &params,
                    float stage_phase, const FastNoiseLite *stage_noise,
                    const PreparedWarpStage &prepared,
                    const FrameState &frame) {
    if (spec.kind == WarpStageKind::NONE)
      return {input, Complex(), 0.0f, 0.0f, 0};
    const float envelope =
        warp_envelope(projected, spec.envelope, params.edge_width);
    const float amplitude = params.strength * envelope;
    Complex output = input;
    float path_length = 0.0f;
    uint8_t flags = 0;
    switch (spec.kind) {
    case WarpStageKind::NONE:
      break;
    case WarpStageKind::LEGACY_STEREO_NOISE: {
      if (params.strength == 0.0f)
        return {input, Complex(), 0.0f, 0.0f, 0};
      const float r_sq = projected.coords.re * projected.coords.re +
                         projected.coords.im * projected.coords.im;
      const StereoWarpResult result = sample_wrapped_warp(
          input, r_sq, *stage_noise, params.scale, params.strength,
          frame.params.projection.pole_fade, frame.clocks.warp_time);
      return {result.coords, result.delta, result.displacement,
              result.displacement, 0};
    }
    case WarpStageKind::AFFINE_FRAME: {
      const float c = prepared.rotation_cos;
      const float s = prepared.rotation_sin;
      const float x = input.re - params.translation_x;
      const float y = input.im - params.translation_y;
      const float rx = c * x + s * y;
      const float ry = -s * x + c * y;
      const Complex transformed(rx / params.scale_x -
                                    params.shear * ry / params.scale_y,
                                ry / params.scale_y);
      output = transformed;
      break;
    }
    case WarpStageKind::WAVE_SHEAR: {
      if (params.strength == 0.0f)
        return {input, Complex(), 0.0f, 0.0f, 0};
      const float c = cosf(params.field_angle);
      const float s = sinf(params.field_angle);
      const float phase = params.frequency * (c * input.re + s * input.im) +
                          TWO_PI_F * stage_phase;
      const float offset = amplitude * sinf(phase);
      const Complex delta(-s * offset, c * offset);
      return {{input.re + delta.re, input.im + delta.im},
              delta,
              fabsf(offset),
              fabsf(offset),
              0};
    }
    case WarpStageKind::VORTEX: {
      const float center_x = prepared.vortex_center_x;
      const float center_y = prepared.vortex_center_y;
      const float x = input.re - center_x;
      const float y = input.im - center_y;
      const float r_sq = x * x + y * y;
      const float angle = prepared.vortex_angle_numerator /
                          (1.0f + r_sq / prepared.vortex_radius_sq);
      const float c = cosf(angle);
      const float s = sinf(angle);
      output = {center_x + c * x - s * y, center_y + s * x + c * y};
      break;
    }
    case WarpStageKind::VECTOR_NOISE: {
      if (params.strength == 0.0f)
        return {input, Complex(), 0.0f, 0.0f, 0};
      const float nx = sample_wrapped_noise_basis(
          *stage_noise, spec.basis, params.scale * input.re,
          params.scale * input.im, prepared.noise_phase);
      float ny = sample_wrapped_noise_basis(
          *stage_noise, spec.basis, params.scale * input.re + 31.416f,
          params.scale * input.im - 47.853f, prepared.noise_phase, 11.0f);
      float zero_dc_x = nx;
      if (spec.basis == NoiseBasis::RIDGED3) {
        zero_dc_x -= sample_wrapped_noise_basis(
            *stage_noise, spec.basis, params.scale * input.re - 73.271f,
            params.scale * input.im + 19.119f, prepared.noise_phase, 5.0f);
        ny -= sample_wrapped_noise_basis(
            *stage_noise, spec.basis, params.scale * input.re + 61.731f,
            params.scale * input.im + 89.417f, prepared.noise_phase, -7.0f);
      }
      const float c = cosf(params.vector_angle);
      const float s = sinf(params.vector_angle);
      const Complex delta(amplitude * (c * zero_dc_x - s * ny),
                          amplitude * (s * zero_dc_x + c * ny));
      const float deformation =
          sqrtf(delta.re * delta.re + delta.im * delta.im);
      return {{input.re + delta.re, input.im + delta.im},
              delta,
              deformation,
              deformation,
              0};
    }
    case WarpStageKind::CURL_FLOW: {
      if (params.strength == 0.0f)
        return {input, Complex(), 0.0f, 0.0f, 0};
      Complex delta;
      output = curl_flow(input, *stage_noise, spec, params, amplitude,
                         prepared.noise_phase, delta, path_length);
      const float deformation =
          spec.curl_integrator == CurlIntegrator::EULER_1
              ? path_length
              : sqrtf(delta.re * delta.re + delta.im * delta.im);
      return {output, delta, deformation, path_length, 0};
    }
    case WarpStageKind::MIRROR_TILE: {
      const Complex transformed = mirror_tile(input, params, prepared);
      output = transformed;
      flags = WARP_FLAG_REFLECTED;
      break;
    }
    case WarpStageKind::POLAR_CHART: {
      const float radius = sqrtf(input.re * input.re + input.im * input.im);
      const float radial = spec.polar_mode == PolarMode::LOGARITHMIC
                               ? logf(std::max(radius, 1.0f / 4096.0f))
                               : radius;
      const Complex transformed(
          params.radial_scale * radial + params.radial_phase,
          (static_cast<float>(spec.polar_harmonic) + 1.0f) *
                  fast_atan2(input.im, input.re) +
              params.angular_phase);
      output = transformed;
      if (radius < 1.0f / 4096.0f)
        flags = WARP_FLAG_SINGULAR;
      break;
    }
    }
    const Complex delta(output.re - input.re, output.im - input.im);
    const float deformation = sqrtf(delta.re * delta.re + delta.im * delta.im);
    if (path_length == 0.0f)
      path_length = deformation;
    return {output, delta, deformation, path_length, flags};
  }

  static float warp_envelope(const ProjectedLookup &projected,
                             WarpEnvelope envelope, float edge_width) {
    if (envelope == WarpEnvelope::PROJECTION_WEIGHT)
      return projected.value_weight;
    if (envelope == WarpEnvelope::EDGE_FADE)
      return cubic_kernel(projected.fade_edge_distance / edge_width);
    return 1.0f;
  }

  static PlanarWarpStageResult
  warp_stage_lookup(const Complex &input, const ProjectedLookup &projected,
                    const WarpStageSpec &spec, const WarpStageParams &params,
                    const FrameState &frame) {
    const PreparedWarpStage prepared =
        prepare_warp_stage(params, frame.clocks.warp_outer_phase);
    return warp_stage_lookup(input, projected, spec, params,
                             frame.clocks.warp_outer_phase,
                             frame.resources.outer_warp_noise, prepared, frame);
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

  static float sample_wrapped_noise_basis(const FastNoiseLite &noise,
                                          NoiseBasis basis, float x, float y,
                                          float turns,
                                          float time_offset = 0.0f) {
    return sample_wrapped_noise_basis(noise, basis, x, y,
                                      prepare_noise_phase(turns), time_offset);
  }

  static Complex curl_flow(const Complex &input, const FastNoiseLite &noise,
                           const WarpStageSpec &spec,
                           const WarpStageParams &params, float distance,
                           const WrappedNoisePhase &noise_phase,
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

  static Complex curl_vector(const Complex &p, const FastNoiseLite &noise,
                             NoiseBasis basis, float scale,
                             const WrappedNoisePhase &noise_phase) {
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

  static Complex simplex_curl_vector(const Complex &p,
                                     const FastNoiseLite &noise, float scale,
                                     float time) {
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

  static Complex curl_vector(const Complex &p, const FastNoiseLite &noise,
                             NoiseBasis basis, float scale, float time) {
    return curl_vector(p, noise, basis, scale, prepare_noise_phase(time));
  }

  static Complex mirror_tile(const Complex &input,
                             const WarpStageParams &params,
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

  static Complex mirror_tile(const Complex &input,
                             const WarpStageParams &params) {
    return mirror_tile(input, params, prepare_warp_stage(params, 0.0f));
  }

  static StereoWarpResult sample_wrapped_warp(const Complex &z, float r_sq,
                                              const FastNoiseLite &noise,
                                              float scale, float strength,
                                              float pole_fade, float time) {
    const StereoWarpResult current =
        stereo_noise_warp(z, r_sq, noise, scale, strength, pole_fade, time);
    const float blend_start = STEREO_NOISE_TIME_PERIOD - NOISE_WRAP_BLEND;
    if (time <= blend_start)
      return current;
    const float mix = ease_in_out_sin((time - blend_start) / NOISE_WRAP_BLEND);
    const StereoWarpResult next =
        stereo_noise_warp(z, r_sq, noise, scale, strength, pole_fade,
                          time - STEREO_NOISE_TIME_PERIOD);
    return {{hs::lerp(current.coords.re, next.coords.re, mix),
             hs::lerp(current.coords.im, next.coords.im, mix)},
            {hs::lerp(current.delta.re, next.delta.re, mix),
             hs::lerp(current.delta.im, next.delta.im, mix)},
            hs::lerp(current.displacement, next.displacement, mix)};
  }

  static Complex condition_source_coords(const Complex &coords,
                                         const FrameState &frame) {
    if (frame.slots.function == Function::NOISE_CONTOUR ||
        frame.slots.function == Function::PRIMITIVE_LATTICE)
      return coords;
    return stereo_pattern_args(coords, frame.params.source.pattern_freq);
  }

  static MaterialSample shape_material(float field,
                                       const ProjectedLookup &projected,
                                       const PlanarWarpResult &warped,
                                       const FrameState &frame) {
    const float weight = frame.slots.signal_weight == SignalWeight::PROJECTION
                             ? projected.value_weight
                             : 1.0f;
    float value = hs::clamp((field * weight + 1.0f) * 0.5f, 0.0f, 1.0f);
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
                          (static_cast<float>(frame.params.value.band_count) +
                           1.0f) *
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
      coverage = smooth_ramp(0.0f, frame.params.value.edge_width,
                             projected.fade_edge_distance);
      break;
    }
    coverage *= projected.domain_coverage;
    return {value, coverage, warped.net_delta, warped.deformation,
            warped.path_length};
  }

  static float sample_source(const Complex &p, const FrameState &frame) {
    if (frame.slots.function == Function::COUPLED_DIRECT)
      return sample_pattern(
          p, frame.params.source.complexity, frame.params.source.pattern_mix,
          frame.prepared_source.primary, frame.prepared_source.secondary);
    if (frame.slots.function == Function::NOISE_CONTOUR) {
      const float n = hs::clamp(
          sample_wrapped_noise_basis(*frame.resources.source_noise,
                                     frame.params.source.noise_basis,
                                     frame.params.source.noise_scale * p.re,
                                     frame.params.source.noise_scale * p.im,
                                     frame.clocks.source_noise_time),
          -1.0f, 1.0f);
      const float contrast = frame.params.source.noise_contrast;
      return n * (1.0f + contrast) / (1.0f + contrast * fabsf(n));
    }
    if (frame.slots.function == Function::PRIMITIVE_LATTICE)
      return primitive_lattice(p, frame.params.source);
    return sample_function(frame.slots.function, p, frame.prepared_source);
  }

  static float primitive_lattice(const Complex &p, const SourceParams &params) {
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

  static float smooth_ramp(float edge0, float edge1, float value) {
    return cubic_kernel((value - edge0) / (edge1 - edge0));
  }

  static Color4 colorize(const MaterialSample &sample,
                         const FrameState &frame) {
    if (frame.slots.colorizer == Colorizer::GENERATED_TRIADIC) {
      Color4 color = frame.resources.generated_palette->get(sample.value);
      color.alpha *= sample.coverage;
      return color;
    }
    if (frame.slots.colorizer == Colorizer::DEFORMATION_INK) {
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
    const float value = std::min(sample.value, ONE_BELOW_UNIT);
    const float u = wrap_t(value + frame.breathe_offset);
    Color4 color = frame.resources.liquid_palette->get(u);
    color.alpha *=
        sample.coverage * (1.0f - value * frame.params.colorizer.value_fade);
    if (frame.params.colorizer.hue_shift != 0.0f)
      color = hue_rotate(color, -sample.deformation *
                                    frame.params.colorizer.hue_shift);
    return color;
  }

  static Vector apply_lens(const Vector &v, const FrameState &frame) {
    const SurfaceLens lens = frame.slots.surface_lens;
    switch (lens) {
    case SurfaceLens::NONE:
      return v;
    case SurfaceLens::GLITCH:
      return glitch_lens(v);
    case SurfaceLens::TWIST:
      return twist_lens(v);
    case SurfaceLens::KALEIDOSCOPE:
      return kaleidoscope_lens(v);
    case SurfaceLens::MOBIUS: {
      return mobius_transform(v, frame.params.surface_lens.mobius);
    }
    case SurfaceLens::TANGENT_NOISE:
      return tangent_noise_lens(v, frame);
    }
    __builtin_unreachable();
  }

  static Vector apply_lens(const Vector &v, SurfaceLens lens) {
    switch (lens) {
    case SurfaceLens::NONE:
      return v;
    case SurfaceLens::GLITCH:
      return glitch_lens(v);
    case SurfaceLens::TWIST:
      return twist_lens(v);
    case SurfaceLens::KALEIDOSCOPE:
      return kaleidoscope_lens(v);
    case SurfaceLens::MOBIUS:
    case SurfaceLens::TANGENT_NOISE:
      return v;
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

  static Vector tangent_noise_lens(const Vector &v, const FrameState &frame) {
    const float scale = frame.params.surface_lens.noise_scale;
    const float time = frame.clocks.lens_noise_time;
    const float nx = sample_wrapped_noise_basis(
        *frame.resources.lens_noise, frame.params.surface_lens.noise_basis,
        scale * v.x, scale * v.y, time, scale * v.z);
    const float ny = sample_wrapped_noise_basis(
        *frame.resources.lens_noise, frame.params.surface_lens.noise_basis,
        scale * v.x + 17.0f, scale * v.y - 29.0f, time, scale * v.z + 43.0f);
    Vector tangent_a =
        fabsf(v.y) < 0.9f ? Vector(-v.z, 0.0f, v.x) : Vector(0.0f, v.z, -v.y);
    tangent_a = tangent_a.normalized();
    const Vector tangent_b(v.y * tangent_a.z - v.z * tangent_a.y,
                           v.z * tangent_a.x - v.x * tangent_a.z,
                           v.x * tangent_a.y - v.y * tangent_a.x);
    const float amount = frame.params.surface_lens.amount;
    return nlerp_unit(
        v,
        Vector(v.x + amount * (nx * tangent_a.x + ny * tangent_b.x),
               v.y + amount * (nx * tangent_a.y + ny * tangent_b.y),
               v.z + amount * (nx * tangent_a.z + ny * tangent_b.z)),
        1.0f);
  }

  static Vector twist_lens(const Vector &v) {
    const float angle = TWIST_RATE * v.y;
    const float c = fast_cosf(angle);
    const float s = fast_sinf(angle);
    return Vector(v.x * c - v.z * s, v.y, v.x * s + v.z * c);
  }

  static Vector kaleidoscope_lens(const Vector &v) {
    constexpr float SECTOR = TWO_PI_F / KALEIDOSCOPE_SECTORS;
    const float radius = sqrtf(v.x * v.x + v.z * v.z);
    float azimuth = fmodf(fast_atan2(v.z, v.x) + PI_F, SECTOR);
    if (azimuth > 0.5f * SECTOR)
      azimuth = SECTOR - azimuth;
    return Vector(radius * fast_cosf(azimuth), v.y,
                  radius * fast_sinf(azimuth));
  }

  static Complex project_point(const Vector &v, Projection projection) {
    switch (projection) {
    case Projection::EQUIRECTANGULAR:
      return equirectangular(v);
    case Projection::STEREOGRAPHIC:
      return stereo(v);
    case Projection::GNOMONIC:
      return gnomonic(v);
    case Projection::BONNE:
      return shadierball::bonne_projection(v, 0.0f, 0.25f * PI_F).coords;
    case Projection::PEIRCE_QUINCUNCIAL:
      return shadierball::peirce_projection(v, 0.0f, 1, 0.0f).coords;
    case Projection::AIROCEAN:
      return shadierball::airocean_projection(v, 0.0f, false).coords;
    }
    __builtin_unreachable();
  }

  static Complex equirectangular(const Vector &v,
                                 float central_meridian = 0.0f) {
    const float radius = sqrtf(v.x * v.x + v.z * v.z);
    return {std::fabs(shadierball::wrap_longitude(fast_atan2(v.z, v.x) -
                                                  central_meridian)) *
                radius,
            0.5f * PI_F - fast_acos(v.y)};
  }

  static Complex gnomonic(const Vector &v) {
    float y = v.y;
    if (std::fabs(y) < GNOMONIC_AXIS_EPS)
      y = y < 0.0f ? -GNOMONIC_AXIS_EPS : GNOMONIC_AXIS_EPS;
    return {v.x / y, v.z / y};
  }

  static float sample_function(Function function, const Complex &p,
                               const SourceState &source) {
    switch (function) {
    case Function::TWIN_WAVE:
      return twin_wave(p, source);
    case Function::RINGS:
      return rings(p, source);
    case Function::SPIRAL:
      return spiral(p, source);
    case Function::GRID:
      return grid(p, source);
    case Function::COUPLED_DIRECT:
    case Function::NOISE_CONTOUR:
    case Function::PRIMITIVE_LATTICE:
      break;
    }
    __builtin_unreachable();
  }

  static float twin_wave(const Complex &p, const SourceState &source) {
    const float rotated = p.re * source.angle_cos + p.im * source.angle_sin;
    return 0.5f * (fast_sinf(p.re + source.primary) +
                   fast_sinf(rotated + source.primary));
  }

  static float rings(const Complex &p, const SourceState &source) {
    return fast_sinf(sqrtf(p.re * p.re + p.im * p.im) - source.primary);
  }

  static float spiral(const Complex &p, const SourceState &source) {
    const float radius = sqrtf(p.re * p.re + p.im * p.im);
    const float azimuth = fast_atan2(p.im, p.re);
    return fast_sinf(radius - SPIRAL_ARMS * (azimuth + source.angle) -
                     source.primary);
  }

  static float grid(const Complex &p, const SourceState &source) {
    const float a = p.re * source.angle_cos + p.im * source.angle_sin;
    const float b = -p.re * source.angle_sin + p.im * source.angle_cos;
    return fast_sinf(a + source.primary) * fast_cosf(b - source.primary);
  }

  HS_O3_FN static float sample_pattern(const Complex &p, float complexity,
                                       float pattern_mix, float primary,
                                       float secondary) {
    if (pattern_mix == 1.0f)
      return fast_sinf(p.re + primary) * fast_cosf(p.im - secondary);
    float re = p.re;
    float im = p.im;
    if (complexity != 0.0f) {
      re += complexity * fast_sinf(p.im + primary);
      im += complexity * fast_cosf(p.re - secondary);
    }
    const float coupled = fast_sinf(re) * fast_cosf(im);
    if (pattern_mix == 0.0f)
      return coupled;
    const float direct =
        fast_sinf(p.re + primary) * fast_cosf(p.im - secondary);
    return hs::lerp(coupled, direct, pattern_mix);
  }

  HS_COLD_MEMBER WalkDeltas sample_walk_deltas() {
#if defined(HS_TEST_BUILD)
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
    if (!param_morph.active)
      return;
    const float mix = transition_mix(param_morph.elapsed, param_morph.duration);
    if (mix == 0.0f)
      blend.params = param_morph.from;
    else if (mix == 1.0f)
      blend.params = param_morph.to;
    else if (param_morph.staggered)
      blend.params.lerp_staggered(param_morph.from, param_morph.to, mix);
    else
      blend.params.lerp(param_morph.from, param_morph.to, mix);
  }

  static float transition_mix(uint16_t elapsed, uint16_t duration) {
    if (elapsed == 0)
      return 0.0f;
    if (elapsed >= duration)
      return 1.0f;
    return ease_in_out_sin(static_cast<float>(elapsed) / duration);
  }

  HS_COLD_MEMBER void apply_requested_config() {
    if (requested_config == published_config)
      return;
    if (!hold_admitted(requested_config)) {
      reject_requested_config();
      return;
    }
    const Config current{active_slots, blend.params};
    if (!transition.active &&
        same_parameter_topology(current, requested_config)) {
      if (!stable_parameter_path_admitted(current, requested_config)) {
        reject_requested_config();
        return;
      }
      param_morph.active = false;
      blend.params = requested_config.params;
      published_config = requested_config;
      return;
    }
    if (try_apply_config(requested_config, MANUAL_TRANSITION_FRAMES, false,
                         false)) {
      published_config = requested_config;
      if (!requested_schema_bound)
        rebind_parameters();
    } else if (!transition.active) {
      requested_config = published_config;
      rebind_parameters();
    }
  }

  HS_COLD_MEMBER void reject_requested_config() {
    requested_config = published_config;
    rebind_parameters();
    if (transition.active && transition.continue_choreo)
      anims_paused = false;
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

  static void canonicalize_warp_params(const WarpStageSpec &requested,
                                       WarpStageParams &params) {
    if (requested.kind == WarpStageKind::NONE)
      return;
    if (requested.kind == WarpStageKind::LEGACY_STEREO_NOISE) {
      params.scale = hs::clamp(params.scale, 0.1f, 100.0f);
      params.strength = hs::clamp(params.strength, 0.0f, 30.0f);
      params.time_scale = hs::clamp(params.time_scale, 0.05f, 1.0f);
      return;
    }
    if (requested.kind == WarpStageKind::CURL_FLOW)
      params.strength = 0.0f;
    params.time_scale =
        hs::clamp(params.time_scale, NOISE_RATE_MIN, NOISE_RATE_MAX);
    if (requested.kind == WarpStageKind::WAVE_SHEAR ||
        requested.kind == WarpStageKind::CURL_FLOW)
      params.strength = hs::clamp(params.strength, -4.0f, 4.0f);
    else if (requested.kind == WarpStageKind::VECTOR_NOISE)
      params.strength = hs::clamp(params.strength, 0.0f, 4.0f);
    else
      params.strength = hs::clamp(params.strength, 0.0f, 1.0f);
    params.scale =
        hs::clamp(params.scale, 1.0f / 64.0f,
                  requested.kind == WarpStageKind::CURL_FLOW ? 16.0f : 64.0f);
  }

  HS_COLD_MEMBER bool try_apply_config(const Config &candidate,
                                       uint16_t duration, bool staggered,
                                       bool continue_choreo) {
    if (!hold_admitted(candidate) || duration == 0)
      return false;
    if (transition.active) {
      if (continue_choreo)
        return false;
      if (retarget_transition_destination(candidate))
        return true;
      return false;
    }
    const Config current{active_slots, blend.params};
    if (!transition_admitted(current, candidate))
      return false;
    if (candidate == current) {
      param_morph.active = false;
      return true;
    }
    if (stable_topology(current, candidate)) {
      param_morph = {current.params, candidate.params, 0,   duration,
                     staggered,      continue_choreo,  true};
      return true;
    }
    const TransitionMode mode = transition_mode(current, candidate);
    if (mode == TransitionMode::DUAL_OUTPUT) {
      if (!prepare_resource_union(current, candidate))
        return false;
    } else if (!prepare_resource_union(current, current))
      return false;
    const uint16_t planned_duration =
        mode == TransitionMode::THROUGH_CLEAR && (duration & 1U) != 0
            ? duration + 1
            : duration;
    param_morph.active = false;
    transition = {current,          candidate,       runtime, runtime, 0,
                  planned_duration, continue_choreo, true,    mode};
    return true;
  }

  HS_COLD_MEMBER bool retarget_transition_destination(const Config &candidate) {
    if (!hold_admitted(candidate))
      return false;
    transition.continue_choreo = false;
    if (transition.elapsed == 0) {
      if (!transition_admitted(transition.from_config, candidate))
        return false;
      transition.mode = transition_mode(transition.from_config, candidate);
      if (transition.mode == TransitionMode::DUAL_OUTPUT) {
        if (!prepare_resource_union(transition.from_config, candidate))
          return false;
      } else if (!prepare_resource_union(transition.from_config,
                                         transition.from_config))
        return false;
      transition.to_config = candidate;
      return true;
    }
    if (!transition_admitted(transition.to_config, candidate))
      return false;
    pending_transition.active = true;
    return true;
  }

  HS_COLD_MEMBER void finish_transitions() {
    if (transition.active) {
      if (anims_paused && transition.continue_choreo) {
        const bool valid_manual_pending =
            requested_config != published_config &&
            valid_config(requested_config);
        if (!valid_manual_pending)
          return;
      }
      if (transition.mode == TransitionMode::THROUGH_CLEAR &&
          transition.elapsed == transition.duration / 2)
        HS_CHECK(
            prepare_resource_union(transition.to_config, transition.to_config),
            "through-clear destination resources exceed capacity");
      if (transition.elapsed < transition.duration) {
        ++transition.elapsed;
        return;
      }
      const bool continue_choreo = transition.continue_choreo;
      runtime = transition.to_runtime;
      active_slots = transition.to_config.slots;
      blend.params = transition.to_config.params;
      transition.active = false;
      if (pending_transition.active) {
        pending_transition.active = false;
        try_apply_config(requested_config, MANUAL_TRANSITION_FRAMES, false,
                         false);
        return;
      }
      if (continue_choreo)
        enter_preset();
      return;
    }
    if (!param_morph.active)
      return;
    if (anims_paused && param_morph.continue_choreo)
      return;
    if (param_morph.elapsed < param_morph.duration) {
      ++param_morph.elapsed;
      return;
    }
    const bool continue_choreo = param_morph.continue_choreo;
    blend.params = param_morph.to;
    param_morph.active = false;
    if (continue_choreo)
      enter_preset();
  }

  HS_COLD_MEMBER void publish_live_config() {
    if (anims_paused || transition.active || param_morph.active ||
        requested_config != published_config)
      return;
    published_config = {active_slots, blend.params};
    requested_config = published_config;
  }

  template <typename Enum>
  static constexpr bool enum_at_most(Enum value, Enum last) {
    return static_cast<uint8_t>(value) <= static_cast<uint8_t>(last);
  }

  static constexpr SourceTraits source_traits(Function function) {
    switch (function) {
    case Function::COUPLED_DIRECT:
      return {true, true, true, true, false, true};
    case Function::PRIMITIVE_LATTICE:
      return {true, true, true, true, false, false};
    case Function::TWIN_WAVE:
    case Function::RINGS:
    case Function::SPIRAL:
    case Function::GRID:
    case Function::NOISE_CONTOUR:
      return {true, true, false, false, false, false};
    }
    return {true, true, false, false, false, false};
  }

  static constexpr bool polar_source_compatible(const RequestedConfig &config,
                                                const WarpStageSpec &polar) {
    const SourceTraits traits = source_traits(config.slots.function);
    if (!traits.y_periodic || !traits.polar_angle_compatible)
      return false;
    const float periods = config.params.source.pattern_freq *
                          (static_cast<float>(polar.polar_harmonic) + 1.0f);
    return periods == static_cast<float>(static_cast<int>(periods));
  }

  static constexpr bool strict_seam_compatible(const Config &config) {
    if (!strict_projection(config.slots.projection))
      return true;
    return config.slots.function != Function::NOISE_CONTOUR &&
           config.slots.surface_lens != SurfaceLens::TANGENT_NOISE &&
           !seam_sensitive_warp(config.slots.warp_program.outer.kind) &&
           !seam_sensitive_warp(config.slots.warp_program.inner.kind);
  }

  static constexpr bool valid_config(const RequestedConfig &candidate) {
    const Slots &slots = candidate.slots;
    if (!enum_at_most(slots.function, Function::PRIMITIVE_LATTICE) ||
        !enum_at_most(slots.projection, Projection::AIROCEAN) ||
        !enum_at_most(slots.peirce_layout, PeirceLayout::VERTICAL) ||
        !enum_at_most(slots.airocean_layout, AiroceanLayout::HORIZONTAL) ||
        !enum_at_most(slots.bonne_hemisphere, BonneHemisphere::SOUTH) ||
        !enum_at_most(slots.gnomonic_hemisphere,
                      GnomonicHemispherePolicy::BACK_HEMISPHERE) ||
        !enum_at_most(slots.projection_frame,
                      ProjectionFramePolicy::SPIN_WANDER) ||
        !enum_at_most(slots.surface_lens, SurfaceLens::TANGENT_NOISE) ||
        !enum_at_most(slots.warp_program.outer.kind,
                      WarpStageKind::POLAR_CHART) ||
        !enum_at_most(slots.warp_program.inner.kind,
                      WarpStageKind::POLAR_CHART) ||
        !valid_warp_spec(slots.warp_program.outer) ||
        !valid_warp_spec(slots.warp_program.inner) ||
        !enum_at_most(slots.signal_weight, SignalWeight::PROJECTION) ||
        !enum_at_most(slots.value_transfer, ValueTransfer::SMOOTH_BANDS) ||
        !enum_at_most(slots.coverage, CoveragePolicy::EDGE_FADE) ||
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

  static constexpr bool valid_warp_spec(const WarpStageSpec &spec) {
    return enum_at_most(spec.basis, NoiseBasis::RIDGED3) &&
           enum_at_most(spec.envelope, WarpEnvelope::EDGE_FADE) &&
           enum_at_most(spec.polar_mode, PolarMode::LOGARITHMIC) &&
           enum_at_most(spec.curl_integrator, CurlIntegrator::MIDPOINT_4) &&
           enum_at_most(spec.polar_harmonic, PolarHarmonic::H16);
  }

  static constexpr float abs_value(float value) {
    return value < 0.0f ? -value : value;
  }

  static constexpr int curl_intervals(CurlIntegrator integrator) {
    return integrator == CurlIntegrator::EULER_1      ? 1
           : integrator == CurlIntegrator::MIDPOINT_2 ? 2
                                                      : 4;
  }

  static constexpr float noise_gradient_bound(NoiseBasis) { return 64.0f; }

  static constexpr bool valid_stage_tuple(const WarpStageSpec &spec,
                                          const WarpStageParams &params) {
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

  static constexpr float projection_coordinate_bound(const Config &config) {
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

  static constexpr float stage_coordinate_bound(const WarpStageSpec &spec,
                                                const WarpStageParams &params,
                                                float input_bound) {
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

  static constexpr bool safe_program_bounds(const Config &config) {
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

  static constexpr bool valid_mobius(const MobiusParams &params) {
    const float ad_re = params.a.re * params.d.re - params.a.im * params.d.im;
    const float ad_im = params.a.re * params.d.im + params.a.im * params.d.re;
    const float bc_re = params.b.re * params.c.re - params.b.im * params.c.im;
    const float bc_im = params.b.re * params.c.im + params.b.im * params.c.re;
    const float det_re = ad_re - bc_re;
    const float det_im = ad_im - bc_im;
    const Complex coefficients[] = {params.a, params.b, params.c, params.d};
    size_t pivot = 0;
    float pivot_magnitude_sq = coefficients[0].re * coefficients[0].re +
                               coefficients[0].im * coefficients[0].im;
    for (size_t index = 1; index < std::size(coefficients); ++index) {
      const float magnitude_sq =
          coefficients[index].re * coefficients[index].re +
          coefficients[index].im * coefficients[index].im;
      if (magnitude_sq > pivot_magnitude_sq) {
        pivot = index;
        pivot_magnitude_sq = magnitude_sq;
      }
    }
    const float norm_sq =
        params.a.re * params.a.re + params.a.im * params.a.im +
        params.b.re * params.b.re + params.b.im * params.b.im +
        params.c.re * params.c.re + params.c.im * params.c.im +
        params.d.re * params.d.re + params.d.im * params.d.im;
    const float norm_error = norm_sq > 1.0f ? norm_sq - 1.0f : 1.0f - norm_sq;
    const float pivot_imag = coefficients[pivot].im < 0.0f
                                 ? -coefficients[pivot].im
                                 : coefficients[pivot].im;
    const bool canonical = norm_error <= 2e-5f && pivot_imag <= 2e-5f &&
                           coefficients[pivot].re > 0.0f;
    return canonical && coefficient_in_range(params.a) &&
           coefficient_in_range(params.b) && coefficient_in_range(params.c) &&
           coefficient_in_range(params.d) &&
           det_re * det_re + det_im * det_im >= 1e-6f;
  }

  static constexpr bool coefficient_in_range(const Complex &coefficient) {
    return coefficient.re >= -8.0f && coefficient.re <= 8.0f &&
           coefficient.im >= -8.0f && coefficient.im <= 8.0f;
  }

  static constexpr float max_value(float a, float b) { return a > b ? a : b; }

  static constexpr float min_value(float a, float b) { return a < b ? a : b; }

  static constexpr float max_abs_value(float a, float b) {
    return max_value(abs_value(a), abs_value(b));
  }

  static constexpr void maximize_stage_path(WarpStageParams &out,
                                            const WarpStageParams &a,
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

  static constexpr bool safe_program_path(const Config &from,
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

  static constexpr bool hold_path_admitted(const Config &from,
                                           const Config &to) {
    Config worst = from;
    if (strict_projection(from.slots.projection) &&
        from.slots.surface_lens != SurfaceLens::NONE &&
        from.params.surface_lens.mix != to.params.surface_lens.mix)
      worst.params.surface_lens.mix = 0.5f;
    return within_hold_device_budget(device_cost(worst));
  }

  static constexpr bool polar_pair_stable(const Config &from,
                                          const Config &to) {
    const bool has_polar =
        from.slots.warp_program.outer.kind == WarpStageKind::POLAR_CHART ||
        from.slots.warp_program.inner.kind == WarpStageKind::POLAR_CHART;
    return !has_polar ||
           from.params.source.pattern_freq == to.params.source.pattern_freq;
  }

  static constexpr bool stable_parameter_path_admitted(const Config &from,
                                                       const Config &to) {
    return curl_pair_stable(from.slots.warp_program.outer,
                            from.params.warp.outer, to.params.warp.outer) &&
           curl_pair_stable(from.slots.warp_program.inner,
                            from.params.warp.inner, to.params.warp.inner) &&
           polar_pair_stable(from, to) && safe_program_path(from, to) &&
           hold_path_admitted(from, to);
  }

  static constexpr bool same_parameter_topology(const Config &from,
                                                const Config &to) {
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

  static constexpr bool stable_topology(const Config &from, const Config &to) {
    return hold_admitted(from) && hold_admitted(to) &&
           same_parameter_topology(from, to) &&
           stable_parameter_path_admitted(from, to);
  }

  static constexpr bool curl_pair_stable(const WarpStageSpec &spec,
                                         const WarpStageParams &a,
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

  static constexpr bool transition_admitted(const Config &from,
                                            const Config &to) {
    return hold_admitted(from) && hold_admitted(to);
  }

  HS_COLD_MEMBER void enter_preset() {
    const Choreo &choreo = CHOREO[preset_index];
    if (choreo.dwell_max == 0) {
      begin_blend();
      return;
    }
    timeline.add_pausable(
        0,
        Animation::RandomTimer(choreo.dwell_min, choreo.dwell_max,
                               [this](Canvas &) { begin_blend(); }),
        &anims_paused);
  }

  HS_COLD_MEMBER void begin_blend() {
    const Choreo &choreo = CHOREO[preset_index];
    const size_t next_index = (preset_index + 1) % PRESETS.size();
    const Preset &to = PRESETS[next_index];
    if (try_apply_config(to, choreo.blend_frames, choreo.staggered, true)) {
      preset_index = next_index;
#if defined(HS_PROFILE_ENABLE)
      hs::log("Preset: %u/%u", static_cast<unsigned>(preset_index),
              static_cast<unsigned>(PRESETS.size()));
#endif
      requested_config = to;
      published_config = to;
      rebind_parameters();
      if (!param_morph.active && !transition.active)
        enter_preset();
    } else {
      timeline.add_pausable(
          0, Animation::RandomTimer(1, 1, [this](Canvas &) { begin_blend(); }),
          &anims_paused);
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
  static constexpr uint8_t WARP_FLAG_REFLECTED = 1U << 0;
  static constexpr uint8_t WARP_FLAG_SINGULAR = 1U << 1;
  static constexpr float GNOMONIC_AXIS_EPS = 1e-3f;
  static constexpr float WARP_COORD_LIMIT = 65536.0f;
  static constexpr float NOISE_LATTICE_LIMIT = 1048576.0f;
  static constexpr float SPIRAL_ARMS = 3.0f;
  static constexpr float TWIST_RATE = 3.0f;
  static constexpr float KALEIDOSCOPE_SECTORS = 6.0f;
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
  static constexpr uint16_t MANUAL_TRANSITION_FRAMES = 60;
  static constexpr size_t PARAM_CAPACITY = 64;

  static constexpr const char *FUNCTION_OPTIONS[] = {
      "Twin Wave",        "Rings",         "Spiral",           "Grid",
      "Coupled / Direct", "Noise Contour", "Primitive Lattice"};
  static constexpr const char *FUNCTION_EXPORT_OPTIONS[] = {
      "Function::TWIN_WAVE",        "Function::RINGS",
      "Function::SPIRAL",           "Function::GRID",
      "Function::COUPLED_DIRECT",   "Function::NOISE_CONTOUR",
      "Function::PRIMITIVE_LATTICE"};
  static constexpr int NUM_FUNCTIONS = std::size(FUNCTION_OPTIONS);
  static constexpr const char *PROJECTION_OPTIONS[] = {
      "Equirectangular", "Stereographic",      "Gnomonic",
      "Bonne",           "Peirce Quincuncial", "Dymaxion / Airocean"};
  static constexpr const char *PROJECTION_EXPORT_OPTIONS[] = {
      "Projection::EQUIRECTANGULAR",    "Projection::STEREOGRAPHIC",
      "Projection::GNOMONIC",           "Projection::BONNE",
      "Projection::PEIRCE_QUINCUNCIAL", "Projection::AIROCEAN"};
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
      "None", "Glitch", "Twist", "Kaleidoscope", "Mobius", "Tangent Noise"};
  static constexpr const char *LENS_EXPORT_OPTIONS[] = {
      "SurfaceLens::NONE",   "SurfaceLens::GLITCH",
      "SurfaceLens::TWIST",  "SurfaceLens::KALEIDOSCOPE",
      "SurfaceLens::MOBIUS", "SurfaceLens::TANGENT_NOISE"};
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
  static constexpr const char *POLAR_HARMONIC_OPTIONS[] = {
      "1", "2",  "3",  "4",  "5",  "6",  "7",  "8",
      "9", "10", "11", "12", "13", "14", "15", "16"};
  static constexpr const char *POLAR_HARMONIC_EXPORT_OPTIONS[] = {
      "PolarHarmonic::H1",  "PolarHarmonic::H2",  "PolarHarmonic::H3",
      "PolarHarmonic::H4",  "PolarHarmonic::H5",  "PolarHarmonic::H6",
      "PolarHarmonic::H7",  "PolarHarmonic::H8",  "PolarHarmonic::H9",
      "PolarHarmonic::H10", "PolarHarmonic::H11", "PolarHarmonic::H12",
      "PolarHarmonic::H13", "PolarHarmonic::H14", "PolarHarmonic::H15",
      "PolarHarmonic::H16"};
  static constexpr int NUM_POLAR_HARMONICS = std::size(POLAR_HARMONIC_OPTIONS);
  static constexpr const char *BAND_COUNT_OPTIONS[] = {
      "1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9",  "10", "11",
      "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22",
      "23", "24", "25", "26", "27", "28", "29", "30", "31", "32"};
  static constexpr const char *BAND_COUNT_EXPORT_OPTIONS[] = {
      "BandCount::C1",  "BandCount::C2",  "BandCount::C3",  "BandCount::C4",
      "BandCount::C5",  "BandCount::C6",  "BandCount::C7",  "BandCount::C8",
      "BandCount::C9",  "BandCount::C10", "BandCount::C11", "BandCount::C12",
      "BandCount::C13", "BandCount::C14", "BandCount::C15", "BandCount::C16",
      "BandCount::C17", "BandCount::C18", "BandCount::C19", "BandCount::C20",
      "BandCount::C21", "BandCount::C22", "BandCount::C23", "BandCount::C24",
      "BandCount::C25", "BandCount::C26", "BandCount::C27", "BandCount::C28",
      "BandCount::C29", "BandCount::C30", "BandCount::C31", "BandCount::C32"};
  static constexpr int NUM_BAND_COUNTS = std::size(BAND_COUNT_OPTIONS);
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
      "Opaque", "Projection Weight Squared", "Value Cutout", "Edge Fade"};
  static constexpr const char *COVERAGE_EXPORT_OPTIONS[] = {
      "CoveragePolicy::OPAQUE", "CoveragePolicy::PROJECTION_WEIGHT_SQUARED",
      "CoveragePolicy::VALUE_CUTOUT", "CoveragePolicy::EDGE_FADE"};
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
  static constexpr float PATTERN_FREQ_MIN = 1.0f, PATTERN_FREQ_MAX = 20.0f;
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
  static constexpr float NOISE_SCALE_MIN = 1.0f / 64.0f;
  static constexpr float NOISE_SCALE_MAX = 64.0f;
  static constexpr float NOISE_RATE_MIN = -1.0f / 64.0f;
  static constexpr float NOISE_RATE_MAX = 1.0f / 64.0f;
  static constexpr float CELL_MIN = 1.0f / 64.0f;
  static constexpr float CELL_MAX = 8.0f;
  static constexpr float SOFTNESS_MIN = 1.0f / 1024.0f;
  static constexpr float ORBIT_SPIN_RATE = TWO_PI_F / 300.0f;

  static constexpr bool preset_in_ranges(const Params &p) {
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
           p.source.noise_scale >= NOISE_SCALE_MIN &&
           p.source.noise_scale <= NOISE_SCALE_MAX &&
           p.source.noise_contrast >= 0.0f && p.source.noise_contrast <= 8.0f &&
           p.source.noise_time_rate >= NOISE_RATE_MIN &&
           p.source.noise_time_rate <= NOISE_RATE_MAX &&
           p.source.lattice_cell_scale >= CELL_MIN &&
           p.source.lattice_cell_scale <= CELL_MAX &&
           p.source.lattice_shape_blend >= 0.0f &&
           p.source.lattice_shape_blend <= 1.0f &&
           p.source.lattice_softness >= SOFTNESS_MIN &&
           p.source.lattice_softness <= 0.25f &&
           p.source.lattice_radius >= 1.0f / 64.0f &&
           p.source.lattice_radius <= 0.5f - p.source.lattice_softness &&
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
           p.surface_lens.noise_scale >= NOISE_SCALE_MIN &&
           p.surface_lens.noise_scale <= NOISE_SCALE_MAX &&
           p.surface_lens.noise_rate >= NOISE_RATE_MIN &&
           p.surface_lens.noise_rate <= NOISE_RATE_MAX &&
           p.value.iso_level >= 0.0f && p.value.iso_level <= 1.0f &&
           p.value.iso_width >= SOFTNESS_MIN && p.value.iso_width <= 0.5f &&
           enum_at_most(p.value.band_count, BandCount::C32) &&
           p.value.band_phase >= 0.0f && p.value.band_phase <= TWO_PI_F &&
           p.value.cutout_threshold >= 0.0f &&
           p.value.cutout_threshold <= 1.0f &&
           p.value.cutout_softness >= SOFTNESS_MIN &&
           p.value.cutout_softness <= 0.5f &&
           p.value.edge_width >= SOFTNESS_MIN && p.value.edge_width <= 0.5f &&
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

  static constexpr bool
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
      Function::COUPLED_DIRECT,
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

  static constexpr const char *OUTER_WARP_PARAM_NAMES[] = {
      "Outer Translation X", "Outer Translation Y", "Outer Rotation",
      "Outer Scale X",       "Outer Scale Y",       "Outer Shear",
      "Outer Frequency",     "Outer Field Angle",   "Outer Center X",
      "Outer Center Y",      "Outer Radius",        "Outer Turns",
      "Outer Vector Angle",  "Outer Cell X",        "Outer Cell Y",
      "Outer Offset X",      "Outer Offset Y",      "Outer Radial Scale",
      "Outer Radial Phase",  "Outer Angular Phase", "Outer Edge Width",
      "Outer Center Orbit"};
  static constexpr const char *INNER_WARP_PARAM_NAMES[] = {
      "Inner Translation X", "Inner Translation Y", "Inner Rotation",
      "Inner Scale X",       "Inner Scale Y",       "Inner Shear",
      "Inner Frequency",     "Inner Field Angle",   "Inner Center X",
      "Inner Center Y",      "Inner Radius",        "Inner Turns",
      "Inner Vector Angle",  "Inner Cell X",        "Inner Cell Y",
      "Inner Offset X",      "Inner Offset Y",      "Inner Radial Scale",
      "Inner Radial Phase",  "Inner Angular Phase", "Inner Edge Width",
      "Inner Center Orbit"};
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

  static constexpr Preset diagnostic_preset(uint8_t index) {
    Slots slots{Function::COUPLED_DIRECT,
                Projection::EQUIRECTANGULAR,
                ProjectionFramePolicy::IDENTITY,
                SurfaceLens::NONE,
                {{WarpStageKind::NONE}, {WarpStageKind::NONE}},
                SignalWeight::NONE,
                ValueTransfer::LINEAR,
                CoveragePolicy::OPAQUE,
                Colorizer::DEFORMATION_INK};
    Params params{};
    if (index == 0) {
      slots.projection = Projection::BONNE;
      slots.warp_program.outer.kind = WarpStageKind::MIRROR_TILE;
      slots.warp_program.inner.kind = WarpStageKind::VORTEX;
      params.warp.outer.cell_x = 1.25f;
      params.warp.outer.cell_y = 0.75f;
      params.warp.inner.turns = 0.75f;
      params.warp.inner.center_orbit_radius = 0.2f;
      params.warp.inner.time_scale = 1.0f / 256.0f;
    } else if (index == 1) {
      slots.function = Function::PRIMITIVE_LATTICE;
      slots.projection = Projection::STEREOGRAPHIC;
      slots.warp_program.outer.kind = WarpStageKind::VECTOR_NOISE;
      params.warp.outer.strength = 0.6f;
      params.warp.outer.scale = 1.5f;
      params.warp.outer.time_scale = 0.0f;
      params.source.lattice_cell_scale = 2.0f;
    } else if (index == 2) {
      slots.function = Function::GRID;
      slots.projection = Projection::STEREOGRAPHIC;
      slots.warp_program.outer.kind = WarpStageKind::CURL_FLOW;
      slots.warp_program.outer.basis = NoiseBasis::SIMPLEX;
      slots.warp_program.outer.curl_integrator = CurlIntegrator::EULER_1;
      params.warp.outer.strength = 0.0078125f;
      params.warp.outer.scale = 1.0f;
      params.warp.outer.time_scale = 0.0f;
    } else if (index == 3) {
      slots.projection = Projection::PEIRCE_QUINCUNCIAL;
    } else {
      slots.projection = Projection::AIROCEAN;
    }
    if (slots.projection == Projection::BONNE ||
        slots.projection == Projection::PEIRCE_QUINCUNCIAL ||
        slots.projection == Projection::AIROCEAN) {
      slots.coverage = CoveragePolicy::EDGE_FADE;
      params.value.edge_width = 0.1f;
    }
    return {slots, params};
  }

  static constexpr std::array<Preset, 20> PRESETS = {{
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
       authored_params({14.262f, 0.586f, 0.0f, 1.0f, 0.7f}, {0.1f, 0.87f, 0.3f},
                       {3.527f, ORBIT_SPIN_RATE}, {0.0f},
                       {0.0f, 0.0f, 0.097f, 1.0f}, {0.0f})},
      {LIQUID_STEREO_SLOTS,
       authored_params({8.0f, 0.30f, 0.0f, 1.0f, 0.7f}, {1.5f, 0.5f, 0.3f},
                       {2.0f, ORBIT_SPIN_RATE}, {0.0f},
                       {0.0f, 0.0f, 0.15f, 1.0f}, {0.0f})},
      {LIQUID_STEREO_SLOTS,
       authored_params({7.878f, 0.562f, 0.0f, 1.0f, 0.7f},
                       {47.752f, 2.55f, 0.3f}, {2.843f, ORBIT_SPIN_RATE},
                       {0.0f}, {0.0f, 0.0f, 0.0f, 1.0f}, {0.0f})},
      {LIQUID_STEREO_SLOTS,
       authored_params({1.0f, 0.586f, 0.0f, 1.0f, 0.7f}, {100.0f, 8.67f, 0.3f},
                       {3.432f, ORBIT_SPIN_RATE}, {0.0f},
                       {0.0f, 0.0f, 0.636f, 1.0f}, {0.0f})},
      {LIQUID_STEREO_SLOTS,
       authored_params({13.430163f, 0.385f, 0.0f, 1.0f, 0.7f},
                       {0.1f, 0.0f, 0.3f}, {1.627f, ORBIT_SPIN_RATE}, {0.0f},
                       {0.0f, 0.0f, 0.10404046f, 1.0f}, {0.0f})},
      {LIQUID_STEREO_SLOTS,
       authored_params({13.430163f, 0.385f, 0.0f, 1.0f, 0.7f},
                       {0.1f, 0.0f, 0.3f}, {1.627f, ORBIT_SPIN_RATE}, {1.0f},
                       {0.0f, 0.0f, 0.10404046f, 1.0f}, {0.0f})},
      {LIQUID_STEREO_SLOTS,
       authored_params({1.0f, 0.075f, 0.009122372f, 1.0f, 1.146f},
                       {50.749298f, 30.0f, 0.4699f}, {1.5482996f, 0.020879198f},
                       {0.0f}, {0.25410002f, 0.00015458837f, 0.201f, 0.847f},
                       {0.0030917525f})},
      {LIQUID_STEREO_SLOTS,
       authored_params({1.0f, 0.075f, 0.009122372f, 1.0f, 1.146f},
                       {38.761299f, 30.0f, 0.4699f}, {1.5482996f, 0.020879198f},
                       {0.0f}, {0.25410002f, 0.00015458837f, 0.201f, 0.847f},
                       {0.0030917525f})},
      diagnostic_preset(0),
      diagnostic_preset(1),
      diagnostic_preset(2),
      diagnostic_preset(3),
      diagnostic_preset(4),
  }};
  static_assert(
      [] {
        for (const Preset &preset : PRESETS)
          if (!valid_config(preset))
            return false;
        return true;
      }(),
      "a ShadierBall preset lies outside its registered range");
  static_assert(
      [] {
        for (const Preset &preset : PRESETS)
          if (!hold_admitted(preset))
            return false;
        return true;
      }(),
      "a ShadierBall preset exceeds the calibrated device hold budget");
  static_assert(
      [] {
        for (size_t index = 0; index < PRESETS.size(); ++index)
          if (!transition_admitted(PRESETS[index],
                                   PRESETS[(index + 1) % PRESETS.size()]))
            return false;
        return true;
      }(),
      "a ShadierBall preset edge lacks continuous transition admission");

  static constexpr std::array<Choreo, 20> CHOREO = {{
      {30, 90, 60, true},   {30, 90, 60, true}, {30, 90, 60, true},
      {30, 90, 480, false}, {0, 0, 480, false}, {0, 0, 480, false},
      {0, 0, 480, false},   {0, 0, 480, false}, {0, 0, 480, false},
      {0, 0, 480, false},   {0, 0, 480, false}, {0, 0, 480, false},
      {0, 0, 480, false},   {0, 0, 480, false}, {0, 0, 480, false},
      {0, 0, 480, false},   {0, 0, 480, false}, {0, 0, 480, false},
      {0, 0, 480, false},   {0, 0, 480, false},
  }};
  static_assert(CHOREO.size() == PRESETS.size());

  Orientation<> projection_walk;
  Orientation<> outer_walk;
  Timeline timeline;
  std::array<FastNoiseLite, MAX_NOISE_RESOURCES> noise_resources;
  std::array<NoiseResourceKey, MAX_NOISE_RESOURCES> prepared_noise_keys{};
  size_t prepared_noise_count = 0;
  FastNoiseLite projection_walk_noise;
  FastNoiseLite outer_walk_noise;

  Quaternion base_orientation =
      make_rotation(Vector(0, 0, -1), Vector(0, -1, 0));
  Quaternion projection_walk_prev;
  Quaternion outer_walk_prev;

  float liquid_rotation = 0.0f;
  uint32_t palette_hue = 0;
  PaletteCycler liquid_palette_cycler;
  PaletteCycler generated_palette_cycler;

  Slots active_slots = LIQUID_STEREO_SLOTS;
  RequestedConfig requested_config{LIQUID_STEREO_SLOTS, PRESETS[0].params};
  Config published_config{LIQUID_STEREO_SLOTS, PRESETS[0].params};
  bool requested_schema_bound = false;
  size_t preset_index = 0;
  Blend blend{PRESETS[0].params};
  LookRuntime runtime;
  ParamMorphRuntime param_morph;
  TransitionRuntime transition;
  PendingTransition pending_transition;
#if defined(HS_TEST_BUILD)
  uint32_t walk_step_count = 0;
  uint32_t liquid_palette_step_count = 0;
  uint32_t generated_palette_step_count = 0;
#endif

  static constexpr size_t FOOTPRINT_BYTES =
      gamut_lut_bytes(GAMUT_ANGLE_STEPS, GAMUT_L_STEPS) +
      2 * PaletteCycler::generated_arena_bytes() +
      PARAM_CAPACITY * sizeof(ParamDef);
  static_assert(
      FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
      "ShadierBall persistent footprint exceeds the default partition");
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(ShadierBall)
