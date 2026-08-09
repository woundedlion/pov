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
    warp_noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
    warp_noise.SetFrequency(WARP_NOISE_FREQUENCY);

    requested_config = PRESETS[0];
    published_config = PRESETS[0];

    register_animated_param("Function", &requested_config.slots.function,
                            FUNCTION_OPTIONS, FUNCTION_EXPORT_OPTIONS,
                            NUM_FUNCTIONS);
    register_animated_param("Projection", &requested_config.slots.projection,
                            PROJECTION_OPTIONS, PROJECTION_EXPORT_OPTIONS,
                            NUM_PROJECTIONS);
    register_animated_param(
        "Projection Frame", &requested_config.slots.projection_frame,
        PROJECTION_FRAME_OPTIONS, PROJECTION_FRAME_EXPORT_OPTIONS,
        NUM_PROJECTION_FRAMES);
    register_animated_param("Lens", &requested_config.slots.surface_lens,
                            LENS_OPTIONS, LENS_EXPORT_OPTIONS, NUM_LENSES);
    register_animated_param("Outer Warp",
                            &requested_config.slots.warp_program.outer.kind,
                            WARP_OPTIONS, WARP_EXPORT_OPTIONS, NUM_WARPS);
    register_animated_param("Inner Warp",
                            &requested_config.slots.warp_program.inner.kind,
                            WARP_OPTIONS, WARP_EXPORT_OPTIONS, NUM_WARPS);
    register_animated_param("Signal Weight",
                            &requested_config.slots.signal_weight,
                            SIGNAL_OPTIONS, SIGNAL_EXPORT_OPTIONS, NUM_SIGNALS);
    register_animated_param("Value Transfer",
                            &requested_config.slots.value_transfer,
                            VALUE_TRANSFER_OPTIONS,
                            VALUE_TRANSFER_EXPORT_OPTIONS, NUM_VALUE_TRANSFERS);
    register_animated_param("Coverage", &requested_config.slots.coverage,
                            COVERAGE_OPTIONS, COVERAGE_EXPORT_OPTIONS,
                            NUM_COVERAGE_POLICIES);
    register_animated_param("Colorizer", &requested_config.slots.colorizer,
                            COLORIZER_OPTIONS, COLORIZER_EXPORT_OPTIONS,
                            NUM_COLORIZERS);

    Params &params = requested_config.params;
    register_animated_param("Outer Warp Scale", &params.warp.outer.scale,
                            WARP_SCALE_MIN, WARP_SCALE_MAX);
    register_animated_param("Outer Warp Strength", &params.warp.outer.strength,
                            WARP_STRENGTH_MIN, WARP_STRENGTH_MAX);
    register_animated_param("Outer Warp Time", &params.warp.outer.time_scale,
                            WARP_TIME_MIN, WARP_TIME_MAX);
    register_animated_param("Inner Warp Scale", &params.warp.inner.scale,
                            WARP_SCALE_MIN, WARP_SCALE_MAX);
    register_animated_param("Inner Warp Strength", &params.warp.inner.strength,
                            WARP_STRENGTH_MIN, WARP_STRENGTH_MAX);
    register_animated_param("Inner Warp Time", &params.warp.inner.time_scale,
                            WARP_TIME_MIN, WARP_TIME_MAX);
    register_animated_param("Pattern Freq", &params.source.pattern_freq,
                            PATTERN_FREQ_MIN, PATTERN_FREQ_MAX);
    register_animated_param("Speed", &params.source.speed, SPEED_MIN,
                            SPEED_MAX);
    register_animated_param("Complexity", &params.source.complexity,
                            COMPLEXITY_MIN, COMPLEXITY_MAX);
    register_animated_param("Pattern Mix", &params.source.pattern_mix,
                            PATTERN_MIX_MIN, PATTERN_MIX_MAX);
    register_animated_param("Drift", &params.source.secondary_rate,
                            PHASE2_RATE_MIN, PHASE2_RATE_MAX);
    register_animated_param("Pole Fade", &params.projection.pole_fade,
                            POLE_FADE_MIN, POLE_FADE_MAX);
    register_animated_param("Spin Rate", &params.projection.spin_rate,
                            SPIN_RATE_MIN, SPIN_RATE_MAX);
    register_animated_param("Source Angle Rate", &params.source.angle_rate,
                            WAVE_SPIN_MIN, WAVE_SPIN_MAX);
    register_animated_param("Projection Wander", &params.projection.wander,
                            WANDER_MIN, WANDER_MAX);
    register_animated_param("Outer Wander", &params.outer_camera.wander,
                            WANDER_MIN, WANDER_MAX);
    register_animated_param("Lens Mix", &params.surface_lens.mix, LENS_MIX_MIN,
                            LENS_MIX_MAX);
    register_animated_param("Breathe Depth", &params.colorizer.breathe_depth,
                            BREATHE_MIN, BREATHE_MAX);
    register_animated_param("Cycle Speed", &params.colorizer.cycle_speed,
                            CYCLE_SPEED_MIN, CYCLE_SPEED_MAX);
    register_animated_param("Hue Shift", &params.colorizer.hue_shift,
                            HUE_SHIFT_MIN, HUE_SHIFT_MAX);
    register_animated_param("Value Fade", &params.colorizer.value_fade,
                            VALUE_FADE_MIN, VALUE_FADE_MAX);

    timeline.add(0, Animation::RandomWalk<W>(projection_walk, UP,
                                             projection_walk_noise));
    timeline.add(0, Animation::RandomWalk<W>(outer_walk, UP, outer_walk_noise));

    liquid_palette_cycler.init_generated(persistent_arena, next_liquid_palette,
                                         &liquid_rotation, PALETTE_DWELL_FRAMES,
                                         PALETTE_FADE_FRAMES, ease_in_out_sin,
                                         &anims_paused);
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
      const TransitionFrame frame = prepare_transition_frame();
      auto shader = [&](const Vector &view) HS_O3_FN -> Color4 {
        return shade_transition(view, frame);
      };
      HS_PROFILE(sdb_shader_draw);
      Scan::Shader::draw<W, H, 1>(canvas, shader);
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

private:
  friend struct ::hs_test::shadierball_tests::ShadierBallWhiteBox;

  enum class Function : uint8_t {
    TWIN_WAVE,
    RINGS,
    SPIRAL,
    GRID,
    COUPLED_DIRECT
  };
  enum class Projection : uint8_t { EQUIRECTANGULAR, STEREOGRAPHIC, GNOMONIC };
  enum class SurfaceLens : uint8_t { NONE, GLITCH, TWIST, KALEIDOSCOPE };
  enum class WarpStageKind : uint8_t { NONE, LEGACY_STEREO_NOISE };
  struct WarpStageSpec {
    WarpStageKind kind;

    bool operator==(const WarpStageSpec &) const = default;
  };
  struct WarpProgram {
    WarpStageSpec outer;
    WarpStageSpec inner;

    bool operator==(const WarpProgram &) const = default;
  };
  enum class ProjectionFramePolicy : uint8_t { IDENTITY, SPIN_WANDER };
  enum class SignalWeight : uint8_t { NONE, PROJECTION };
  enum class ValueTransfer : uint8_t { LINEAR };
  enum class CoveragePolicy : uint8_t { OPAQUE, PROJECTION_WEIGHT_SQUARED };
  enum class Colorizer : uint8_t { GENERATED_TRIADIC, LIQUID };

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

    bool operator==(const Slots &) const = default;
  };

  struct SourceParams {
    float pattern_freq;
    float speed;
    float complexity;
    float pattern_mix;
    float secondary_rate;
    float angle_rate = 0.0f;

    bool operator==(const SourceParams &) const = default;

    void lerp(const SourceParams &a, const SourceParams &b, float t) {
      pattern_freq = hs::lerp(a.pattern_freq, b.pattern_freq, t);
      speed = hs::lerp(a.speed, b.speed, t);
      complexity = hs::lerp(a.complexity, b.complexity, t);
      pattern_mix = hs::lerp(a.pattern_mix, b.pattern_mix, t);
      secondary_rate = hs::lerp(a.secondary_rate, b.secondary_rate, t);
      angle_rate = hs::lerp(a.angle_rate, b.angle_rate, t);
    }
  };

  struct WarpStageParams {
    float scale;
    float strength;
    float time_scale;

    bool operator==(const WarpStageParams &) const = default;

    void lerp(const WarpStageParams &a, const WarpStageParams &b, float t) {
      scale = hs::lerp(a.scale, b.scale, t);
      strength = hs::lerp(a.strength, b.strength, t);
      time_scale = hs::lerp(a.time_scale, b.time_scale, t);
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
    float pole_fade;
    float spin_rate;
    float wander = 0.0f;

    bool operator==(const ProjectionParams &) const = default;

    void lerp(const ProjectionParams &a, const ProjectionParams &b, float t) {
      pole_fade = hs::lerp(a.pole_fade, b.pole_fade, t);
      spin_rate = hs::lerp(a.spin_rate, b.spin_rate, t);
      wander = hs::lerp(a.wander, b.wander, t);
    }
  };

  struct SurfaceLensParams {
    float mix;

    bool operator==(const SurfaceLensParams &) const = default;

    void lerp(const SurfaceLensParams &a, const SurfaceLensParams &b, float t) {
      mix = hs::lerp(a.mix, b.mix, t);
    }
  };

  struct ValueParams {
    bool operator==(const ValueParams &) const = default;
    void lerp(const ValueParams &, const ValueParams &, float) {}
  };

  struct ColorizerParams {
    float breathe_depth;
    float cycle_speed;
    float hue_shift;
    float value_fade;

    bool operator==(const ColorizerParams &) const = default;

    void lerp(const ColorizerParams &a, const ColorizerParams &b, float t) {
      breathe_depth = hs::lerp(a.breathe_depth, b.breathe_depth, t);
      cycle_speed = hs::lerp(a.cycle_speed, b.cycle_speed, t);
      hue_shift = hs::lerp(a.hue_shift, b.hue_shift, t);
      value_fade = hs::lerp(a.value_fade, b.value_fade, t);
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
    float source_primary;
    float source_secondary;
    float source_angle;
    float warp_time;
    float projection_spin;
    float breathe_phase;
  };

  struct PreparedTransforms {
    Quaternion projection_conj;
    Quaternion outer_conj;
  };

  struct ResourceBindings {
    const FastNoiseLite *warp_noise;
    const BakedPalette *generated_palette;
    const BakedPalette *liquid_palette;
  };

  struct FrameState {
    Slots slots;
    Params params;
    ClockState clocks;
    SourceState prepared_source;
    float breathe_offset;
    PreparedTransforms transforms;
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

  struct TransitionFrame {
    FrameState from;
    FrameState to;
    float mix;
  };

  HS_COLD_MEMBER FrameState prepare_frame() const {
    return prepare_frame({active_slots, blend.params}, runtime);
  }

  HS_COLD_MEMBER FrameState prepare_frame(const Config &config,
                                          const LookRuntime &look) const {
    const bool animated_projection =
        config.slots.projection_frame == ProjectionFramePolicy::SPIN_WANDER;
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
        {&warp_noise, &generated_palette_cycler.palette(),
         &liquid_palette_cycler.palette()}};
  }

  HS_COLD_MEMBER TransitionFrame prepare_transition_frame() const {
    return {prepare_frame(transition.from_config, transition.from_runtime),
            prepare_frame(transition.to_config, transition.to_runtime),
            transition_mix(transition.elapsed, transition.duration)};
  }

  static Color4 shade(const Vector &view, const FrameState &frame) {
    const Vector outer_local = outer_camera_lookup(view, frame);
    const ProjectedLookup projected =
        surface_lens_project_lookup(outer_local, frame);
    const PlanarWarpResult warped = planar_warp_lookup(projected, frame);
    const Complex source_coords = condition_source_coords(warped.coords, frame);
    const float field = sample_source(source_coords, frame);
    const MaterialSample material =
        shape_material(field, projected, warped, frame);
    return colorize(material, frame);
  }

  static Color4 shade_transition(const Vector &view,
                                 const TransitionFrame &frame) {
    if (frame.mix == 0.0f)
      return shade(view, frame.from);
    if (frame.mix == 1.0f)
      return shade(view, frame.to);
    return blend_outputs(shade(view, frame.from), shade(view, frame.to),
                         frame.mix);
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
    const Vector lensed = apply_lens(v, slots.surface_lens);
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
    return finalize_projection(
        local, project_point(local, frame.slots.projection),
        frame.slots.projection, frame.params.projection.pole_fade);
  }

  static ProjectedLookup finalize_projection(const Vector &local,
                                             const Complex &coords,
                                             Projection projection,
                                             float pole_fade) {
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
    case Projection::GNOMONIC:
      return {coords,
              static_cast<uint8_t>(local.y < 0.0f),
              static_cast<uint8_t>(local.y < 0.0f),
              static_cast<uint8_t>(BOUNDARY_CUT | BOUNDARY_SINGULAR),
              std::fabs(local.y),
              pole_attenuation(r_sq, pole_fade),
              0};
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
    }
    const float r_sq = coords.re * coords.re + coords.im * coords.im;
    return {coords,
            selected->region_id,
            selected->component_id,
            selected->boundary_flags,
            selected->fade_edge_distance,
            pole_attenuation(r_sq, pole_fade),
            selected->flags};
  }

  static PlanarWarpResult planar_warp_lookup(const ProjectedLookup &projected,
                                             const FrameState &frame) {
    const PlanarWarpStageResult outer = warp_stage_lookup(
        projected.coords, projected, frame.slots.warp_program.outer,
        frame.params.warp.outer, frame);
    const PlanarWarpStageResult inner = warp_stage_lookup(
        outer.coords, projected, frame.slots.warp_program.inner,
        frame.params.warp.inner, frame);
    const Complex net_delta(outer.delta.re + inner.delta.re,
                            outer.delta.im + inner.delta.im);
    const bool sole_legacy = (frame.slots.warp_program.outer.kind ==
                              WarpStageKind::LEGACY_STEREO_NOISE) !=
                             (frame.slots.warp_program.inner.kind ==
                              WarpStageKind::LEGACY_STEREO_NOISE);
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
                    const FrameState &frame) {
    if (spec.kind == WarpStageKind::NONE || params.strength == 0.0f)
      return {input, Complex(), 0.0f, 0.0f, 0};
    const float r_sq = projected.coords.re * projected.coords.re +
                       projected.coords.im * projected.coords.im;
    const StereoWarpResult result = sample_wrapped_warp(
        input, r_sq, *frame.resources.warp_noise, params.scale, params.strength,
        frame.params.projection.pole_fade, frame.clocks.warp_time);
    return {result.coords, result.delta, result.displacement,
            result.displacement, 0};
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
    }
    const float coverage =
        frame.slots.coverage == CoveragePolicy::PROJECTION_WEIGHT_SQUARED
            ? projected.value_weight * projected.value_weight
            : 1.0f;
    return {value, coverage, warped.net_delta, warped.deformation,
            warped.path_length};
  }

  static float sample_source(const Complex &p, const FrameState &frame) {
    if (frame.slots.function == Function::COUPLED_DIRECT)
      return sample_pattern(
          p, frame.params.source.complexity, frame.params.source.pattern_mix,
          frame.prepared_source.primary, frame.prepared_source.secondary);
    return sample_function(frame.slots.function, p, frame.prepared_source);
  }

  static Color4 colorize(const MaterialSample &sample,
                         const FrameState &frame) {
    if (frame.slots.colorizer == Colorizer::GENERATED_TRIADIC) {
      Color4 color = frame.resources.generated_palette->get(sample.value);
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
    }
    __builtin_unreachable();
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
    }
    __builtin_unreachable();
  }

  static Complex equirectangular(const Vector &v) {
    const float radius = sqrtf(v.x * v.x + v.z * v.z);
    return {std::fabs(fast_atan2(v.z, v.x)) * radius,
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
    if (try_apply_config(requested_config, MANUAL_TRANSITION_FRAMES, false,
                         false))
      published_config = requested_config;
    else if (!transition.active && valid_config(requested_config))
      requested_config = published_config;
  }

  HS_COLD_MEMBER bool try_apply_config(const Config &candidate,
                                       uint16_t duration, bool staggered,
                                       bool continue_choreo) {
    if (!valid_config(candidate) || duration == 0)
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
    if (candidate.slots == current.slots) {
      param_morph = {current.params, candidate.params, 0,   duration,
                     staggered,      continue_choreo,  true};
      return true;
    }
    param_morph.active = false;
    transition = {current, candidate, runtime,         runtime,
                  0,       duration,  continue_choreo, true};
    return true;
  }

  HS_COLD_MEMBER bool retarget_transition_destination(const Config &candidate) {
    if (!transition_admitted(transition.from_config, candidate))
      return false;
    transition.to_config = candidate;
    transition.continue_choreo = false;
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
      if (transition.elapsed < transition.duration) {
        ++transition.elapsed;
        return;
      }
      const bool continue_choreo = transition.continue_choreo;
      runtime = transition.to_runtime;
      active_slots = transition.to_config.slots;
      blend.params = transition.to_config.params;
      transition.active = false;
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
    if (anims_paused || transition.active ||
        requested_config != published_config)
      return;
    published_config = {active_slots, blend.params};
    requested_config = published_config;
  }

  template <typename Enum>
  static constexpr bool enum_at_most(Enum value, Enum last) {
    return static_cast<uint8_t>(value) <= static_cast<uint8_t>(last);
  }

  static constexpr bool valid_config(const RequestedConfig &candidate) {
    const Slots &slots = candidate.slots;
    if (!enum_at_most(slots.function, Function::COUPLED_DIRECT) ||
        !enum_at_most(slots.projection, Projection::GNOMONIC) ||
        !enum_at_most(slots.projection_frame,
                      ProjectionFramePolicy::SPIN_WANDER) ||
        !enum_at_most(slots.surface_lens, SurfaceLens::KALEIDOSCOPE) ||
        !enum_at_most(slots.warp_program.outer.kind,
                      WarpStageKind::LEGACY_STEREO_NOISE) ||
        !enum_at_most(slots.warp_program.inner.kind,
                      WarpStageKind::LEGACY_STEREO_NOISE) ||
        !enum_at_most(slots.signal_weight, SignalWeight::PROJECTION) ||
        !enum_at_most(slots.value_transfer, ValueTransfer::LINEAR) ||
        !enum_at_most(slots.coverage,
                      CoveragePolicy::PROJECTION_WEIGHT_SQUARED) ||
        !enum_at_most(slots.colorizer, Colorizer::LIQUID))
      return false;
    const int legacy_stages =
        (slots.warp_program.outer.kind == WarpStageKind::LEGACY_STEREO_NOISE) +
        (slots.warp_program.inner.kind == WarpStageKind::LEGACY_STEREO_NOISE);
    if (legacy_stages > 1)
      return false;
    if (legacy_stages == 1 && slots.projection != Projection::STEREOGRAPHIC)
      return false;
    if (!preset_in_ranges(candidate.params))
      return false;
    return true;
  }

  static constexpr bool transition_admitted(const Config &from,
                                            const Config &to) {
    return valid_config(from) && valid_config(to) && from.slots == to.slots;
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
      requested_config = to;
      published_config = to;
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
  static constexpr float GNOMONIC_AXIS_EPS = 1e-3f;
  static constexpr float SPIRAL_ARMS = 3.0f;
  static constexpr float TWIST_RATE = 3.0f;
  static constexpr float KALEIDOSCOPE_SECTORS = 6.0f;
  static constexpr float WARP_NOISE_FREQUENCY = 0.01f;
  static constexpr float STEREO_NOISE_TIME_PERIOD = 65536.0f;
  static constexpr float ONE_BELOW_UNIT = 0x1.fffffep-1f;
  static constexpr float NOISE_WRAP_BLEND = 1024.0f;
  static constexpr float GOLDEN_HUE_STEP = 0.618034f;
  static constexpr uint32_t HUE_STEP = 159;
  static constexpr int PALETTE_DWELL_FRAMES = 0;
  static constexpr int PALETTE_FADE_FRAMES = 600;
  static constexpr uint16_t MANUAL_TRANSITION_FRAMES = 60;

  static constexpr const char *FUNCTION_OPTIONS[] = {
      "Twin Wave", "Rings", "Spiral", "Grid", "Coupled / Direct"};
  static constexpr const char *FUNCTION_EXPORT_OPTIONS[] = {
      "Function::TWIN_WAVE", "Function::RINGS", "Function::SPIRAL",
      "Function::GRID", "Function::COUPLED_DIRECT"};
  static constexpr int NUM_FUNCTIONS = std::size(FUNCTION_OPTIONS);
  static constexpr const char *PROJECTION_OPTIONS[] = {
      "Equirectangular", "Stereographic", "Gnomonic"};
  static constexpr const char *PROJECTION_EXPORT_OPTIONS[] = {
      "Projection::EQUIRECTANGULAR", "Projection::STEREOGRAPHIC",
      "Projection::GNOMONIC"};
  static constexpr int NUM_PROJECTIONS = std::size(PROJECTION_OPTIONS);
  static constexpr const char *PROJECTION_FRAME_OPTIONS[] = {"Identity",
                                                             "Spin + Wander"};
  static constexpr const char *PROJECTION_FRAME_EXPORT_OPTIONS[] = {
      "ProjectionFramePolicy::IDENTITY", "ProjectionFramePolicy::SPIN_WANDER"};
  static constexpr int NUM_PROJECTION_FRAMES =
      std::size(PROJECTION_FRAME_OPTIONS);
  static constexpr const char *LENS_OPTIONS[] = {"None", "Glitch", "Twist",
                                                 "Kaleidoscope"};
  static constexpr const char *LENS_EXPORT_OPTIONS[] = {
      "SurfaceLens::NONE", "SurfaceLens::GLITCH", "SurfaceLens::TWIST",
      "SurfaceLens::KALEIDOSCOPE"};
  static constexpr int NUM_LENSES = std::size(LENS_OPTIONS);
  static constexpr const char *WARP_OPTIONS[] = {"None", "Stereo Noise"};
  static constexpr const char *WARP_EXPORT_OPTIONS[] = {
      "WarpStageKind::NONE", "WarpStageKind::LEGACY_STEREO_NOISE"};
  static constexpr int NUM_WARPS = std::size(WARP_OPTIONS);
  static constexpr const char *SIGNAL_OPTIONS[] = {"None", "Projection"};
  static constexpr const char *SIGNAL_EXPORT_OPTIONS[] = {
      "SignalWeight::NONE", "SignalWeight::PROJECTION"};
  static constexpr int NUM_SIGNALS = std::size(SIGNAL_OPTIONS);
  static constexpr const char *VALUE_TRANSFER_OPTIONS[] = {"Linear"};
  static constexpr const char *VALUE_TRANSFER_EXPORT_OPTIONS[] = {
      "ValueTransfer::LINEAR"};
  static constexpr int NUM_VALUE_TRANSFERS = std::size(VALUE_TRANSFER_OPTIONS);
  static constexpr const char *COVERAGE_OPTIONS[] = {
      "Opaque", "Projection Weight Squared"};
  static constexpr const char *COVERAGE_EXPORT_OPTIONS[] = {
      "CoveragePolicy::OPAQUE", "CoveragePolicy::PROJECTION_WEIGHT_SQUARED"};
  static constexpr int NUM_COVERAGE_POLICIES = std::size(COVERAGE_OPTIONS);
  static constexpr const char *COLORIZER_OPTIONS[] = {"Generated Triadic",
                                                      "ShaderBall Liquid"};
  static constexpr const char *COLORIZER_EXPORT_OPTIONS[] = {
      "Colorizer::GENERATED_TRIADIC", "Colorizer::LIQUID"};
  static constexpr int NUM_COLORIZERS = std::size(COLORIZER_OPTIONS);

  static constexpr float WARP_SCALE_MIN = 0.1f, WARP_SCALE_MAX = 100.0f;
  static constexpr float WARP_STRENGTH_MIN = 0.0f, WARP_STRENGTH_MAX = 30.0f;
  static constexpr float WARP_TIME_MIN = 0.05f, WARP_TIME_MAX = 1.0f;
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
           p.projection.pole_fade >= POLE_FADE_MIN &&
           p.projection.pole_fade <= POLE_FADE_MAX &&
           p.projection.spin_rate >= SPIN_RATE_MIN &&
           p.projection.spin_rate <= SPIN_RATE_MAX &&
           p.projection.wander >= WANDER_MIN &&
           p.projection.wander <= WANDER_MAX &&
           p.outer_camera.wander >= WANDER_MIN &&
           p.outer_camera.wander <= WANDER_MAX &&
           p.surface_lens.mix >= LENS_MIX_MIN &&
           p.surface_lens.mix <= LENS_MIX_MAX &&
           p.colorizer.breathe_depth >= BREATHE_MIN &&
           p.colorizer.breathe_depth <= BREATHE_MAX &&
           p.colorizer.cycle_speed >= CYCLE_SPEED_MIN &&
           p.colorizer.cycle_speed <= CYCLE_SPEED_MAX &&
           p.colorizer.hue_shift >= HUE_SHIFT_MIN &&
           p.colorizer.hue_shift <= HUE_SHIFT_MAX &&
           p.colorizer.value_fade >= VALUE_FADE_MIN &&
           p.colorizer.value_fade <= VALUE_FADE_MAX;
  }

  static constexpr bool
  warp_stage_params_in_ranges(const WarpStageParams &params) {
    return params.scale >= WARP_SCALE_MIN && params.scale <= WARP_SCALE_MAX &&
           params.strength >= WARP_STRENGTH_MIN &&
           params.strength <= WARP_STRENGTH_MAX &&
           params.time_scale >= WARP_TIME_MIN &&
           params.time_scale <= WARP_TIME_MAX;
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
      Colorizer::LIQUID};
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

  static constexpr std::array<Preset, 15> PRESETS = {{
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
        for (size_t index = 0; index < PRESETS.size(); ++index)
          if (!transition_admitted(PRESETS[index],
                                   PRESETS[(index + 1) % PRESETS.size()]))
            return false;
        return true;
      }(),
      "a ShadierBall preset edge lacks continuous transition admission");

  static constexpr std::array<Choreo, 15> CHOREO = {{
      {30, 90, 60, true},
      {30, 90, 60, true},
      {30, 90, 60, true},
      {30, 90, 480, false},
      {0, 0, 480, false},
      {0, 0, 480, false},
      {0, 0, 480, false},
      {0, 0, 480, false},
      {0, 0, 480, false},
      {0, 0, 480, false},
      {0, 0, 480, false},
      {0, 0, 480, false},
      {0, 0, 480, false},
      {0, 0, 480, false},
      {0, 0, 480, false},
  }};
  static_assert(CHOREO.size() == PRESETS.size());

  Orientation<> projection_walk;
  Orientation<> outer_walk;
  Timeline timeline;
  FastNoiseLite warp_noise;
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
  size_t preset_index = 0;
  Blend blend{PRESETS[0].params};
  LookRuntime runtime;
  ParamMorphRuntime param_morph;
  TransitionRuntime transition;
#if defined(HS_TEST_BUILD)
  uint32_t walk_step_count = 0;
  uint32_t liquid_palette_step_count = 0;
  uint32_t generated_palette_step_count = 0;
#endif

  static constexpr size_t FOOTPRINT_BYTES =
      gamut_lut_bytes(GAMUT_ANGLE_STEPS, GAMUT_L_STEPS) +
      2 * PaletteCycler::generated_arena_bytes();
  static_assert(
      FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
      "ShadierBall persistent footprint exceeds the default partition");
};

#include "core/engine/effect_registry.h"
REGISTER_EFFECT(ShadierBall)
