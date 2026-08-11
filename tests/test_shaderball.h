/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <limits>
#include "tests/test_effects.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace shaderball_tests {

using effects_tests::FRAME_MS;
using effects_tests::FRAME_US;
using effects_tests::reset_effect_globals;
using effects_tests::SMALL_H;
using effects_tests::SMALL_W;

/** @brief White-box access to ShaderBall's typed pipeline. */
struct ShaderBallWhiteBox {
  using SB = ShaderBall<SMALL_W, SMALL_H>;
  using Function = SB::Function;
  using Projection = SB::Projection;
  using BonneHemisphere = SB::BonneHemisphere;
  using GnomonicHemispherePolicy = SB::GnomonicHemispherePolicy;
  using ProjectionFramePolicy = SB::ProjectionFramePolicy;
  using SurfaceLens = SB::SurfaceLens;
  using NoiseBasis = SB::NoiseBasis;
  using WarpEnvelope = SB::WarpEnvelope;
  using PolarMode = SB::PolarMode;
  using CurlIntegrator = SB::CurlIntegrator;
  using WarpStageKind = SB::WarpStageKind;
  using WarpStageSpec = SB::WarpStageSpec;
  using WarpStageParams = SB::WarpStageParams;
  using ProjectionParams = SB::ProjectionParams;
  using SignalWeight = SB::SignalWeight;
  using ValueTransfer = SB::ValueTransfer;
  using CoveragePolicy = SB::CoveragePolicy;
  using Colorizer = SB::Colorizer;
  using Slots = SB::Slots;
  using Params = SB::Params;
  using RequestedConfig = SB::RequestedConfig;
  using SourceState = SB::SourceState;
  using FrameState = SB::FrameState;
  using ProjectedLookup = SB::ProjectedLookup;
  using PlanarWarpStageResult = SB::PlanarWarpStageResult;
  using PlanarWarpResult = SB::PlanarWarpResult;
  using MaterialSample = SB::MaterialSample;
  using ClockState = SB::ClockState;
  using LookRuntime = SB::LookRuntime;
  using WalkDeltas = SB::WalkDeltas;
  using ThroughClearPhase = SB::ThroughClearPhase;
  using CostTier = SB::CostTier;
  using DeviceCost = SB::DeviceCost;

  static constexpr float AXIS_EPS = SB::GNOMONIC_AXIS_EPS;
  static constexpr uint32_t HUE_STEP = SB::HUE_STEP;

  static ClockState clocks(const SB &sb) { return sb.runtime.clocks; }
  static void seed_clocks(SB &sb, float value) {
    sb.runtime.clocks = {value, value, value, value, value, value};
  }
  static FrameState frame(const SB &sb) { return sb.prepare_frame(); }
  static Slots active_slots(const SB &sb) { return sb.active_slots; }
  static RequestedConfig active_config(const SB &sb) {
    return {sb.active_slots, sb.blend.params};
  }
  static constexpr Slots liquid_stereo_slots() {
    return SB::LIQUID_STEREO_SLOTS;
  }
  static constexpr Slots legacy_slots() {
    return {Function::TWIN_WAVE,
            Projection::STEREOGRAPHIC,
            ProjectionFramePolicy::IDENTITY,
            SurfaceLens::GLITCH,
            {{WarpStageKind::NONE}, {WarpStageKind::NONE}},
            SignalWeight::PROJECTION,
            ValueTransfer::LINEAR,
            CoveragePolicy::PROJECTION_WEIGHT_SQUARED,
            Colorizer::GENERATED_TRIADIC};
  }
  static constexpr RequestedConfig legacy_config() {
    return {legacy_slots(),
            {{6.0f, 0.05f, 0.0f, 0.0f, 0.8f, 0.006f},
             {{0.1f, 0.0f, 0.5f}, {0.1f, 0.0f, 0.5f}},
             {2.0f, 0.0f, 0.0f},
             {1.0f},
             {},
             {0.0f, 0.0f, 0.0f, 0.0f},
             {0.25f}}};
  }
  static void request_slots(SB &sb, const Slots &slots) {
    sb.requested_config.slots = slots;
    sb.requested_schema_bound = false;
    sb.apply_requested_config();
  }
  static Slots requested_slots(const SB &sb) {
    return sb.requested_config.slots;
  }
  static const RequestedConfig &requested_config(const SB &sb) {
    return sb.requested_config;
  }
  static const RequestedConfig &published_config(const SB &sb) {
    return sb.published_config;
  }
  static void request_config(SB &sb, const RequestedConfig &config) {
    sb.requested_config = config;
    sb.requested_schema_bound = false;
    sb.apply_requested_config();
  }
  static bool slots_equal(const Slots &a, const Slots &b) { return a == b; }
  static constexpr bool valid_config(const RequestedConfig &config) {
    return SB::valid_config(config);
  }
  static constexpr bool seam_compatible(const RequestedConfig &config) {
    return SB::strict_seam_compatible(config);
  }
  static constexpr bool transition_admitted(const RequestedConfig &from,
                                            const RequestedConfig &to) {
    return SB::transition_admitted(from, to);
  }
  static constexpr bool hold_admitted(const RequestedConfig &config) {
    return SB::hold_admitted(config);
  }
  static constexpr bool stable_topology(const RequestedConfig &from,
                                        const RequestedConfig &to) {
    return SB::stable_topology(from, to);
  }
  static constexpr bool
  stable_parameter_path_admitted(const RequestedConfig &from,
                                 const RequestedConfig &to) {
    return SB::stable_parameter_path_admitted(from, to);
  }
  static constexpr uint16_t hold_device_point_budget() {
    return SB::HOLD_DEVICE_POINT_BUDGET;
  }
  static constexpr uint16_t noise_call_points() {
    return SB::NOISE_CALL_POINTS;
  }
  static constexpr float noise_time_period() {
    return SB::STEREO_NOISE_TIME_PERIOD;
  }
  static bool transition_active(const SB &sb) { return sb.transition.active; }
  static bool param_morph_active(const SB &sb) { return sb.param_morph.active; }
  static uint16_t param_morph_elapsed(const SB &sb) {
    return sb.param_morph.elapsed;
  }
  static const Params &live_params(const SB &sb) { return sb.blend.params; }
  static size_t preset_index(const SB &sb) { return sb.preset_index; }
  static float transition_mix(const SB &sb) {
    return SB::transition_mix(sb.transition.elapsed, sb.transition.duration);
  }
  static const LookRuntime &transition_from_runtime(const SB &sb) {
    return sb.transition.from_runtime;
  }
  static const LookRuntime &transition_to_runtime(const SB &sb) {
    return sb.transition.to_runtime;
  }
  static const RequestedConfig &transition_from_config(const SB &sb) {
    return sb.transition.from_config;
  }
  static const RequestedConfig &transition_to_config(const SB &sb) {
    return sb.transition.to_config;
  }
  static bool transition_continues_choreo(const SB &sb) {
    return sb.transition.continue_choreo;
  }
  static uint16_t transition_elapsed(const SB &sb) {
    return sb.transition.elapsed;
  }
  static NoiseBasis prepared_noise_basis(const SB &sb, uint8_t resource_id) {
    for (size_t index = 0; index < sb.prepared_noise_count; ++index)
      if (sb.prepared_noise_keys[index].resource_id == resource_id)
        return sb.prepared_noise_keys[index].basis;
    return static_cast<NoiseBasis>(0xff);
  }
  static constexpr DeviceCost device_cost(const RequestedConfig &config) {
    return SB::device_cost(config);
  }
  static constexpr CostTier cost_tier(const RequestedConfig &config) {
    return SB::cost_tier(SB::device_cost(config));
  }
  static void force_transition(SB &sb, const RequestedConfig &to,
                               uint16_t duration, bool continue_choreo) {
    const RequestedConfig from = active_config(sb);
    sb.param_morph.active = false;
    sb.transition = {from, to,       sb.runtime,      sb.runtime,
                     0,    duration, continue_choreo, true};
  }
  static const LookRuntime &runtime(const SB &sb) { return sb.runtime; }
  static Quaternion projection_walk(const SB &sb) {
    return sb.projection_walk.get();
  }
  static Quaternion outer_walk(const SB &sb) { return sb.outer_walk.get(); }
  static void advance_runtime(SB &sb, LookRuntime &runtime,
                              const RequestedConfig &config,
                              const WalkDeltas &deltas) {
    sb.advance_runtime(runtime, config, deltas);
  }
  static ThroughClearPhase through_clear_phase(uint16_t elapsed,
                                               uint16_t duration) {
    return SB::through_clear_phase(elapsed, duration);
  }
  static Color4 shade_through_clear(const Vector &view,
                                    const FrameState *visible,
                                    const ThroughClearPhase &phase) {
    return SB::shade_through_clear(view, visible, phase);
  }
  static void begin_blend(SB &sb) { sb.begin_blend(); }
  static void step_param_morph(SB &sb) {
    sb.prepare_param_morph();
    sb.advance_runtime(sb.runtime, {sb.active_slots, sb.blend.params},
                       {Quaternion(), Quaternion()});
    sb.finish_transitions();
  }
  static void settle_transition(SB &sb) {
    for (int frame = 0;
         frame < 1024 && (sb.transition.active || sb.param_morph.active);
         ++frame) {
      sb.draw_frame();
      sb.advance_display();
    }
  }
  static uint32_t walk_steps(const SB &sb) { return sb.walk_step_count; }
  static uint32_t liquid_palette_steps(const SB &sb) {
    return sb.liquid_palette_step_count;
  }
  static uint32_t generated_palette_steps(const SB &sb) {
    return sb.generated_palette_step_count;
  }
  static Color4 blend_outputs(const Color4 &from, const Color4 &to, float mix) {
    return SB::blend_outputs(from, to, mix);
  }
  static ProjectedLookup join(const ProjectedLookup &direct,
                              const ProjectedLookup &lensed, float mix,
                              Projection projection, float pole_fade) {
    return SB::join_projected(direct, lensed, mix, projection, pole_fade);
  }
  static bool join_compatible(const ProjectedLookup &direct,
                              const ProjectedLookup &lensed,
                              Projection projection,
                              float coordinate_scale = 1.0f) {
    return SB::projection_join_compatible(direct, lensed, projection,
                                          coordinate_scale);
  }
  static constexpr uint8_t boundary_cut() { return SB::BOUNDARY_CUT; }
  static constexpr uint8_t boundary_singular() { return SB::BOUNDARY_SINGULAR; }
  static constexpr uint8_t projection_folded() {
    return SB::PROJECTION_FLAG_FOLDED;
  }
  static Vector outer_lookup(const Vector &v, const FrameState &frame) {
    return SB::outer_camera_lookup(v, frame);
  }
  static ProjectedLookup surface_project(const Vector &v,
                                         const FrameState &frame) {
    return SB::surface_lens_project_lookup(v, frame);
  }
  static PlanarWarpResult warp(const ProjectedLookup &projected,
                               const FrameState &frame) {
    return SB::planar_warp_lookup(projected, frame);
  }
  /**
   * @brief Runs one warp stage against a frame.
   * @param inner Selects the inner stage's clock and noise resource; the outer
   *        stage's otherwise. A stage reads whichever pair it was programmed
   *        into, so the caller states which one `spec` describes.
   */
  static PlanarWarpStageResult
  warp_stage(const Complex &input, const ProjectedLookup &projected,
             const WarpStageSpec &spec, const WarpStageParams &params,
             const FrameState &frame, bool inner = false) {
    const float phase =
        inner ? frame.clocks.warp_inner_phase : frame.clocks.warp_outer_phase;
    const FastNoiseLite *noise = inner ? frame.resources.inner_warp_noise
                                       : frame.resources.outer_warp_noise;
    return SB::warp_stage_lookup(input, projected, spec, params, phase, noise,
                                 SB::prepare_warp_stage(params, phase), frame);
  }
  static MaterialSample material(const ProjectedLookup &projected,
                                 const PlanarWarpResult &warped,
                                 const FrameState &frame) {
    const Complex source_coords =
        SB::condition_source_coords(warped.coords, frame);
    const float field = SB::sample_source(source_coords, frame);
    return SB::shape_material(field, projected, warped, frame);
  }
  static MaterialSample shape(float field, const ProjectedLookup &projected,
                              const PlanarWarpResult &warped,
                              const FrameState &frame) {
    return SB::shape_material(field, projected, warped, frame);
  }
  static float edge_fade_width(const ProjectedLookup &projected,
                               const FrameState &frame) {
    return SB::edge_fade_width(projected, frame);
  }
  static Color4 colorize(const MaterialSample &sample,
                         const FrameState &frame) {
    return SB::colorize(sample, frame);
  }
  static Color4 shade(const Vector &v, const FrameState &frame) {
    return SB::shade(v, frame);
  }
  static Complex project_point(const Vector &v, Projection projection) {
    return SB::project_point(v, projection);
  }
  static ProjectedLookup
  finalize_projection(const Vector &v, const Complex &coords,
                      Projection projection, float pole_fade,
                      GnomonicHemispherePolicy hemisphere) {
    return SB::finalize_projection(v, coords, projection, pole_fade,
                                   hemisphere);
  }
  static void canonicalize_mobius(MobiusParams &params) {
    SB::canonicalize_mobius(params);
  }
  static Complex curl_vector(const Complex &p, const FastNoiseLite &noise,
                             NoiseBasis basis, float scale, float time) {
    return SB::curl_vector(p, noise, basis, scale,
                           SB::prepare_noise_phase(time));
  }
  static float wrapped_noise(const FastNoiseLite &noise, NoiseBasis basis,
                             float x, float y, float turns) {
    return SB::sample_wrapped_noise_basis(noise, basis, x, y, turns);
  }
  static ProjectionParams lerp_projection(const ProjectionParams &a,
                                          const ProjectionParams &b, float t) {
    ProjectionParams result;
    result.lerp(a, b, t);
    return result;
  }
  static Vector apply_lens(const Vector &v, SurfaceLens lens) {
    return SB::apply_frame_free_lens(v, lens);
  }
  static float sample_function(Function function, const Complex &p,
                               const SourceState &source) {
    return SB::sample_function(function, p, source);
  }
  static float sample_pattern(const Complex &p, float complexity, float mix,
                              float primary, float secondary) {
    return SB::sample_pattern(p, complexity, mix, primary, secondary);
  }
  static const auto &presets() { return SB::PRESETS; }
  static const auto &choreo() { return SB::CHOREO; }
  static void make_triadic(uint32_t &hue, uint32_t sequence,
                           GenerativePalette &out) {
    SB::next_triadic_palette(&hue, sequence, out);
  }
  static Pixel generated_color(const SB &sb, float value) {
    return sb.generated_palette_cycler.palette().get(value).color;
  }
  static Pixel liquid_color(const SB &sb, float value) {
    return sb.liquid_palette_cycler.palette().get(value).color;
  }
};

/** @brief Every named clock wraps in its native domain. */
inline void test_shaderball_clocks_wrapped() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  hs::set_mock_time(0, 0);
  WB::SB sb;
  sb.init();
  WB::seed_clocks(sb, WB::noise_time_period() * 4.0f);
  for (int frame = 0; frame < 32; ++frame) {
    hs::set_mock_time(frame * FRAME_MS, frame * FRAME_US);
    sb.draw_frame();
    sb.advance_display();
    const WB::ClockState clocks = WB::clocks(sb);
    HS_EXPECT_GE(clocks.warp_time, 0.0f);
    HS_EXPECT_LT(clocks.warp_time, WB::noise_time_period());
    for (float phase :
         {clocks.source_primary, clocks.source_secondary, clocks.source_angle,
          clocks.projection_spin, clocks.breathe_phase,
          clocks.source_noise_time, clocks.lens_noise_time,
          clocks.warp_outer_phase, clocks.warp_inner_phase}) {
      HS_EXPECT_GE(phase, 0.0f);
      HS_EXPECT_LT(phase, TWO_PI_F);
    }
  }
  hs::clear_mock_time();
}

/** @brief Pause freezes choreography while ambient motion keeps advancing. */
inline void test_shaderball_pause_semantics() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();
  WB::RequestedConfig ambient = WB::presets()[0];
  ambient.params.source.speed = 0.019f;
  ambient.params.source.secondary_rate = 0.37f;
  ambient.params.source.angle_rate = 0.007f;
  ambient.params.source.noise_time_rate = 0.004f;
  ambient.params.surface_lens.noise_rate = 0.003f;
  ambient.params.warp.outer.time_scale = 0.071f;
  ambient.params.warp.inner.time_scale = 0.003f;
  ambient.params.projection.spin_rate = 0.009f;
  ambient.params.colorizer.cycle_speed = 0.013f;
  WB::request_config(sb, ambient);
  sb.setAnimationsPaused(true);

  const WB::ClockState paused_clocks = WB::clocks(sb);
  const Quaternion paused_projection_walk = WB::projection_walk(sb);
  const Quaternion paused_outer_walk = WB::outer_walk(sb);
  const uint32_t paused_walk_steps = WB::walk_steps(sb);
  const uint32_t paused_liquid_steps = WB::liquid_palette_steps(sb);
  const uint32_t paused_generated_steps = WB::generated_palette_steps(sb);
  const Pixel paused_liquid_color = WB::liquid_color(sb, 0.25f);
  const Pixel paused_generated_color = WB::generated_color(sb, 0.25f);
  const size_t paused_preset = WB::preset_index(sb);
  for (int frame = 0; frame < 120; ++frame) {
    sb.draw_frame();
    sb.advance_display();
  }
  HS_EXPECT_NE(WB::clocks(sb).source_primary, paused_clocks.source_primary);
  HS_EXPECT_NE(WB::clocks(sb).source_secondary, paused_clocks.source_secondary);
  HS_EXPECT_NE(WB::clocks(sb).source_angle, paused_clocks.source_angle);
  HS_EXPECT_NE(WB::clocks(sb).warp_time, paused_clocks.warp_time);
  HS_EXPECT_NE(WB::clocks(sb).projection_spin, paused_clocks.projection_spin);
  HS_EXPECT_NE(WB::clocks(sb).breathe_phase, paused_clocks.breathe_phase);
  HS_EXPECT_NE(WB::clocks(sb).source_noise_time,
               paused_clocks.source_noise_time);
  HS_EXPECT_NE(WB::clocks(sb).lens_noise_time, paused_clocks.lens_noise_time);
  HS_EXPECT_NE(WB::clocks(sb).warp_outer_phase, paused_clocks.warp_outer_phase);
  HS_EXPECT_NE(WB::clocks(sb).warp_inner_phase, paused_clocks.warp_inner_phase);
  HS_EXPECT_TRUE(WB::projection_walk(sb) != paused_projection_walk);
  HS_EXPECT_TRUE(WB::outer_walk(sb) != paused_outer_walk);
  HS_EXPECT_EQ(WB::walk_steps(sb), paused_walk_steps + 120);
  HS_EXPECT_EQ(WB::liquid_palette_steps(sb), paused_liquid_steps + 120);
  HS_EXPECT_EQ(WB::generated_palette_steps(sb), paused_generated_steps + 120);
  const Pixel active_liquid_color = WB::liquid_color(sb, 0.25f);
  const Pixel active_generated_color = WB::generated_color(sb, 0.25f);
  HS_EXPECT_TRUE(active_liquid_color.r != paused_liquid_color.r ||
                 active_liquid_color.g != paused_liquid_color.g ||
                 active_liquid_color.b != paused_liquid_color.b);
  HS_EXPECT_TRUE(active_generated_color.r != paused_generated_color.r ||
                 active_generated_color.g != paused_generated_color.g ||
                 active_generated_color.b != paused_generated_color.b);
  HS_EXPECT_EQ(WB::preset_index(sb), paused_preset);
  HS_EXPECT_FALSE(WB::transition_active(sb));
  HS_EXPECT_FALSE(WB::param_morph_active(sb));

  sb.setAnimationsPaused(false);
  sb.draw_frame();
  sb.advance_display();
  HS_EXPECT_NE(WB::clocks(sb).source_primary, paused_clocks.source_primary);
  HS_EXPECT_TRUE(WB::projection_walk(sb) != paused_projection_walk);
  HS_EXPECT_TRUE(WB::outer_walk(sb) != paused_outer_walk);
  HS_EXPECT_EQ(WB::walk_steps(sb), paused_walk_steps + 121);
  HS_EXPECT_EQ(WB::liquid_palette_steps(sb), paused_liquid_steps + 121);
  HS_EXPECT_EQ(WB::generated_palette_steps(sb), paused_generated_steps + 121);

  for (int frame = 0; frame < 120 && WB::preset_index(sb) == paused_preset;
       ++frame) {
    sb.draw_frame();
    sb.advance_display();
  }
  HS_EXPECT_NE(WB::preset_index(sb), paused_preset);
  HS_EXPECT_TRUE(WB::param_morph_active(sb));
  sb.setAnimationsPaused(true);
  const uint16_t paused_morph_elapsed = WB::param_morph_elapsed(sb);
  for (int frame = 0; frame < 8; ++frame) {
    sb.draw_frame();
    sb.advance_display();
  }
  HS_EXPECT_EQ(WB::param_morph_elapsed(sb), paused_morph_elapsed);
  sb.setAnimationsPaused(false);
  sb.draw_frame();
  sb.advance_display();
  HS_EXPECT_EQ(WB::param_morph_elapsed(sb), paused_morph_elapsed + 1);
}

/** @brief A paused topology edit commits on the next frame. */
inline void test_shaderball_paused_selector_commit() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();
  sb.setAnimationsPaused(true);

  HS_EXPECT_TRUE(
      sb.updateParameter("Inner Warp",
                         static_cast<float>(WB::WarpStageKind::CURL_FLOW)) ==
      ParamSetResult::APPLIED);
  HS_EXPECT_EQ(WB::requested_slots(sb).warp_program.inner.kind,
               WB::WarpStageKind::CURL_FLOW);
  HS_EXPECT_EQ(WB::requested_slots(sb).warp_program.outer.kind,
               WB::WarpStageKind::NONE);
  HS_EXPECT_NE(WB::active_slots(sb).warp_program.inner.kind,
               WB::WarpStageKind::CURL_FLOW);
  sb.draw_frame();
  sb.advance_display();
  HS_EXPECT_FALSE(WB::transition_active(sb));
  HS_EXPECT_FALSE(WB::param_morph_active(sb));
  HS_EXPECT_EQ(WB::active_slots(sb).warp_program.inner.kind,
               WB::WarpStageKind::CURL_FLOW);
  HS_EXPECT_EQ(sb.getParameters().find("Inner Warp")->get(),
               static_cast<float>(WB::WarpStageKind::CURL_FLOW));
}

/** @brief Manual parameter edits commit on the next frame. */
inline void test_shaderball_manual_edit_timing() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();
  sb.setAnimationsPaused(true);
  sb.draw_frame();
  sb.advance_display();
  const float before_speed = WB::active_config(sb).params.source.speed;

  sb.setAnimationsPaused(false);
  HS_EXPECT_TRUE(sb.updateParameter("Speed", before_speed + 0.5f) ==
                 ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(sb.animations_paused());
  HS_EXPECT_EQ(WB::active_config(sb).params.source.speed, before_speed);
  sb.draw_frame();
  sb.advance_display();
  HS_EXPECT_FALSE(WB::param_morph_active(sb));
  HS_EXPECT_EQ(WB::active_config(sb).params.source.speed, before_speed + 0.5f);
  HS_EXPECT_EQ(WB::frame(sb).params.source.speed, before_speed + 0.5f);
  HS_EXPECT_TRUE(WB::requested_config(sb) == WB::active_config(sb));
  HS_EXPECT_TRUE(WB::published_config(sb) == WB::active_config(sb));

  for (WB::Projection projection :
       {WB::Projection::BONNE, WB::Projection::PEIRCE_QUINCUNCIAL,
        WB::Projection::AIROCEAN}) {
    HS_EXPECT_TRUE(
        sb.updateParameter("Projection", static_cast<float>(projection)) ==
        ParamSetResult::APPLIED);
    HS_EXPECT_NE(WB::active_slots(sb).projection, projection);
    sb.draw_frame();
    sb.advance_display();
    HS_EXPECT_FALSE(WB::transition_active(sb));
    HS_EXPECT_FALSE(WB::param_morph_active(sb));
    HS_EXPECT_EQ(WB::active_slots(sb).projection, projection);
    HS_EXPECT_EQ(WB::active_slots(sb).coverage, WB::CoveragePolicy::EDGE_FADE);
    HS_EXPECT_TRUE(WB::hold_admitted(WB::active_config(sb)));
    HS_EXPECT_TRUE(WB::requested_config(sb) == WB::active_config(sb));
    HS_EXPECT_TRUE(WB::published_config(sb) == WB::active_config(sb));
  }
  HS_EXPECT_EQ(WB::active_slots(sb).surface_lens,
               WB::liquid_stereo_slots().surface_lens);
  HS_EXPECT_EQ(WB::active_slots(sb).warp_program.outer.kind,
               WB::WarpStageKind::NONE);
  HS_EXPECT_EQ(WB::active_slots(sb).warp_program.inner.kind,
               WB::WarpStageKind::NONE);
}

/** @brief Pullback stages preserve their typed order and metadata. */
inline void test_shaderball_pipeline_contract() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();
  WB::FrameState frame = WB::frame(sb);
  frame.slots = WB::legacy_slots();
  frame.slots.surface_lens = WB::SurfaceLens::NONE;
  frame.slots.warp_program.outer.kind = WB::WarpStageKind::NONE;
  frame.slots.warp_program.inner.kind = WB::WarpStageKind::NONE;
  frame.transforms.outer_conj = make_rotation(X_AXIS, 0.3f);
  frame.transforms.projection_conj = make_rotation(Y_AXIS, -0.2f);
  const Vector view = Vector(0.4f, -0.8f, 0.3f).normalized();
  const Vector outer = WB::outer_lookup(view, frame);
  const Vector expected_outer = rotate(view, frame.transforms.outer_conj);
  HS_EXPECT_EQ(outer.x, expected_outer.x);
  HS_EXPECT_EQ(outer.y, expected_outer.y);
  HS_EXPECT_EQ(outer.z, expected_outer.z);

  const WB::ProjectedLookup projected = WB::surface_project(outer, frame);
  const Complex expected =
      stereo(rotate(expected_outer, frame.transforms.projection_conj));
  HS_EXPECT_EQ(projected.coords.re, expected.re);
  HS_EXPECT_EQ(projected.coords.im, expected.im);
  const float radius_sq = expected.re * expected.re + expected.im * expected.im;
  HS_EXPECT_EQ(projected.value_weight,
               pole_attenuation(radius_sq, frame.params.projection.pole_fade));
  HS_EXPECT_TRUE(projected.boundary_flags != 0);
  HS_EXPECT_TRUE(std::isfinite(projected.fade_edge_distance));

  const WB::PlanarWarpResult warped = WB::warp(projected, frame);
  HS_EXPECT_EQ(warped.coords.re, projected.coords.re);
  HS_EXPECT_EQ(warped.coords.im, projected.coords.im);
  HS_EXPECT_EQ(warped.deformation, 0.0f);

  frame.slots.warp_program.outer.kind = WB::WarpStageKind::LEGACY_STEREO_NOISE;
  frame.slots.warp_program.inner.kind = WB::WarpStageKind::NONE;
  frame.params.warp.outer = {1.0f, 0.7f, 0.5f};
  const WB::PlanarWarpStageResult outer_stage = WB::warp_stage(
      projected.coords, projected, frame.slots.warp_program.outer,
      frame.params.warp.outer, frame);
  const WB::PlanarWarpResult outer_only = WB::warp(projected, frame);
  HS_EXPECT_EQ(outer_only.coords.re, outer_stage.coords.re);
  HS_EXPECT_EQ(outer_only.coords.im, outer_stage.coords.im);
  HS_EXPECT_EQ(outer_only.net_delta.re, outer_stage.delta.re);
  HS_EXPECT_EQ(outer_only.net_delta.im, outer_stage.delta.im);
  HS_EXPECT_EQ(outer_only.path_length, outer_stage.path_length);
  HS_EXPECT_EQ(outer_stage.coords.re,
               projected.coords.re + outer_stage.delta.re);
  HS_EXPECT_EQ(outer_stage.coords.im,
               projected.coords.im + outer_stage.delta.im);

  frame.slots.warp_program.outer.kind = WB::WarpStageKind::NONE;
  frame.slots.warp_program.inner.kind = WB::WarpStageKind::LEGACY_STEREO_NOISE;
  frame.resources.inner_warp_noise = frame.resources.outer_warp_noise;
  frame.params.warp.inner = {2.0f, 1.1f, 0.5f};
  const WB::PlanarWarpStageResult inner_stage = WB::warp_stage(
      projected.coords, projected, frame.slots.warp_program.inner,
      frame.params.warp.inner, frame, /*inner=*/true);
  const WB::PlanarWarpResult inner_only = WB::warp(projected, frame);
  HS_EXPECT_EQ(inner_only.coords.re, inner_stage.coords.re);
  HS_EXPECT_EQ(inner_only.coords.im, inner_stage.coords.im);
  HS_EXPECT_EQ(inner_only.net_delta.re, inner_stage.delta.re);
  HS_EXPECT_EQ(inner_only.net_delta.im, inner_stage.delta.im);
  HS_EXPECT_EQ(inner_only.path_length, inner_stage.path_length);

  frame.resources.outer_warp_noise = nullptr;
  frame.params.warp.inner.strength = 0.0f;
  const WB::PlanarWarpStageResult zero_strength = WB::warp_stage(
      projected.coords, projected, frame.slots.warp_program.inner,
      frame.params.warp.inner, frame, /*inner=*/true);
  HS_EXPECT_EQ(zero_strength.coords.re, projected.coords.re);
  HS_EXPECT_EQ(zero_strength.coords.im, projected.coords.im);
  HS_EXPECT_EQ(zero_strength.delta.re, 0.0f);
  HS_EXPECT_EQ(zero_strength.delta.im, 0.0f);
  HS_EXPECT_EQ(zero_strength.deformation, 0.0f);

  frame.slots.warp_program.outer.kind = WB::WarpStageKind::NONE;
  frame.slots.warp_program.inner.kind = WB::WarpStageKind::NONE;
  const WB::MaterialSample material = WB::material(projected, warped, frame);
  HS_EXPECT_GE(material.value, 0.0f);
  HS_EXPECT_LE(material.value, 1.0f);
  HS_EXPECT_EQ(material.coverage,
               projected.value_weight * projected.value_weight);
  frame.slots.coverage = WB::CoveragePolicy::PROJECTION_WEIGHT;
  HS_EXPECT_EQ(WB::material(projected, warped, frame).coverage,
               projected.value_weight);
  const Color4 color = WB::colorize(material, frame);
  HS_EXPECT_TRUE(color.alpha >= 0.0f);

  frame.slots.signal_weight = WB::SignalWeight::NONE;
  HS_EXPECT_EQ(WB::shape(3.0f, projected, warped, frame).value, 1.0f);
  HS_EXPECT_EQ(WB::shape(-3.0f, projected, warped, frame).value, 0.0f);

  const WB::ProjectedLookup direct_meta{
      Complex(1.0f, 2.0f), 3, 4, WB::boundary_singular(), 0.25f, 0.1f, 0x10};
  const WB::ProjectedLookup lensed_meta{
      Complex(3.0f, 4.0f), 7, 8, WB::boundary_singular(), 0.75f, 0.9f, 0x20};
  const WB::ProjectedLookup joined = WB::join(
      direct_meta, lensed_meta, 0.75f, WB::Projection::STEREOGRAPHIC, 2.0f);
  const WB::ProjectedLookup direct_endpoint = WB::join(
      direct_meta, lensed_meta, 0.0f, WB::Projection::STEREOGRAPHIC, 2.0f);
  const WB::ProjectedLookup lensed_endpoint = WB::join(
      direct_meta, lensed_meta, 1.0f, WB::Projection::STEREOGRAPHIC, 2.0f);
  HS_EXPECT_EQ(direct_endpoint.coords.re, direct_meta.coords.re);
  HS_EXPECT_EQ(direct_endpoint.coords.im, direct_meta.coords.im);
  HS_EXPECT_EQ(direct_endpoint.region_id, direct_meta.region_id);
  HS_EXPECT_EQ(direct_endpoint.component_id, direct_meta.component_id);
  HS_EXPECT_EQ(direct_endpoint.boundary_flags, direct_meta.boundary_flags);
  HS_EXPECT_EQ(direct_endpoint.fade_edge_distance,
               direct_meta.fade_edge_distance);
  HS_EXPECT_EQ(direct_endpoint.value_weight, direct_meta.value_weight);
  HS_EXPECT_EQ(direct_endpoint.flags, direct_meta.flags);
  HS_EXPECT_EQ(lensed_endpoint.coords.re, lensed_meta.coords.re);
  HS_EXPECT_EQ(lensed_endpoint.coords.im, lensed_meta.coords.im);
  HS_EXPECT_EQ(lensed_endpoint.region_id, lensed_meta.region_id);
  HS_EXPECT_EQ(lensed_endpoint.component_id, lensed_meta.component_id);
  HS_EXPECT_EQ(lensed_endpoint.boundary_flags, lensed_meta.boundary_flags);
  HS_EXPECT_EQ(lensed_endpoint.fade_edge_distance,
               lensed_meta.fade_edge_distance);
  HS_EXPECT_EQ(lensed_endpoint.value_weight, lensed_meta.value_weight);
  HS_EXPECT_EQ(lensed_endpoint.flags, lensed_meta.flags);
  HS_EXPECT_EQ(joined.coords.re, 2.5f);
  HS_EXPECT_EQ(joined.coords.im, 3.5f);
  HS_EXPECT_EQ(joined.region_id, uint8_t(7));
  HS_EXPECT_EQ(joined.component_id, uint8_t(8));
  HS_EXPECT_EQ(joined.boundary_flags, WB::boundary_singular());
  HS_EXPECT_EQ(joined.fade_edge_distance, 0.75f);
  HS_EXPECT_EQ(joined.flags, uint8_t(0x20));
  HS_EXPECT_EQ(joined.value_weight,
               pole_attenuation(2.5f * 2.5f + 3.5f * 3.5f, 2.0f));
}

/** @brief Legacy projection and lens slots retain their shipped kernels. */
inline void test_shaderball_legacy_spatial_slots() {
  using WB = ShaderBallWhiteBox;
  const Vector directions[] = {Vector(1, 0, 0), Vector(0, 1, 0),
                               Vector(0, -1, 0), Vector(1, 1, 1).normalized(),
                               Vector(-1, 2, -3).normalized()};
  for (const Vector &v : directions) {
    const Complex stereo_actual =
        WB::project_point(v, WB::Projection::STEREOGRAPHIC);
    const Complex stereo_expected = stereo(v);
    HS_EXPECT_EQ(stereo_actual.re, stereo_expected.re);
    HS_EXPECT_EQ(stereo_actual.im, stereo_expected.im);

    const Complex equirect =
        WB::project_point(v, WB::Projection::EQUIRECTANGULAR);
    const float radius = sqrtf(v.x * v.x + v.z * v.z);
    HS_EXPECT_NEAR(equirect.re, std::fabs(fast_atan2(v.z, v.x)) * radius,
                   1e-6f);
    HS_EXPECT_EQ(equirect.im, 0.5f * PI_F - fast_acos(v.y));

    const Complex gnomonic = WB::project_point(v, WB::Projection::GNOMONIC);
    HS_EXPECT_TRUE(std::isfinite(gnomonic.re));
    HS_EXPECT_TRUE(std::isfinite(gnomonic.im));
    for (WB::SurfaceLens lens :
         {WB::SurfaceLens::GLITCH, WB::SurfaceLens::TWIST,
          WB::SurfaceLens::KALEIDOSCOPE})
      HS_EXPECT_NEAR(WB::apply_lens(v, lens).length(), 1.0f, 4e-3f);
  }

  reset_effect_globals();
  WB::SB landmark_sb;
  landmark_sb.init();
  WB::FrameState landmark_frame = WB::frame(landmark_sb);
  landmark_frame.slots.projection = WB::Projection::STEREOGRAPHIC;
  landmark_frame.slots.projection_frame = WB::ProjectionFramePolicy::IDENTITY;
  landmark_frame.slots.surface_lens = WB::SurfaceLens::NONE;
  landmark_frame.transforms.projection_conj = Quaternion();
  HS_EXPECT_EQ(WB::surface_project(Vector(0.0f, 1.0f, 0.0f), landmark_frame)
                   .fade_edge_distance,
               0.0f);
  HS_EXPECT_EQ(WB::surface_project(Vector(0.0f, -1.0f, 0.0f), landmark_frame)
                   .fade_edge_distance,
               2.0f);
  HS_EXPECT_EQ(WB::surface_project(Vector(1.0f, 0.0f, 0.0f), landmark_frame)
                   .boundary_flags,
               WB::boundary_singular());
  landmark_frame.slots.projection = WB::Projection::EQUIRECTANGULAR;
  const WB::ProjectedLookup equirectangular =
      WB::surface_project(Vector(1.0f, 0.0f, 0.0f), landmark_frame);
  HS_EXPECT_EQ(equirectangular.boundary_flags, uint8_t(0));
  HS_EXPECT_EQ(equirectangular.flags, WB::projection_folded());
  HS_EXPECT_EQ(
      WB::surface_project(Vector(0.0f, 0.0f, -1.0f), landmark_frame).region_id,
      uint8_t(1));
  landmark_frame.slots.projection = WB::Projection::GNOMONIC;
  HS_EXPECT_EQ(
      WB::surface_project(Vector(1.0f, 0.0f, 0.0f), landmark_frame)
          .boundary_flags,
      static_cast<uint8_t>(WB::boundary_cut() | WB::boundary_singular()));

  const Vector v(0.6f, 0.48f, 0.64f);
  const Complex direct = stereo(v);
  const Complex lensed = stereo(glitch_lens(v));
  WB::FrameState frame = WB::frame(landmark_sb);
  frame.slots.projection = WB::Projection::STEREOGRAPHIC;
  frame.slots.projection_frame = WB::ProjectionFramePolicy::IDENTITY;
  frame.slots.surface_lens = WB::SurfaceLens::GLITCH;
  frame.transforms.projection_conj = Quaternion();
  frame.params.surface_lens.mix = 0.0f;
  const WB::ProjectedLookup start = WB::surface_project(v, frame);
  frame.params.surface_lens.mix = 1.0f;
  const WB::ProjectedLookup end = WB::surface_project(v, frame);
  HS_EXPECT_EQ(start.coords.re, direct.re);
  HS_EXPECT_EQ(start.coords.im, direct.im);
  // glitch_lens' polar terms are FMA-contractable, so the reference above and
  // the pipeline's own call can round 2 ULP apart under -O2; a wrong lens
  // branch would miss by ~1.
  HS_EXPECT_NEAR(end.coords.re, lensed.re, 1e-6f);
  HS_EXPECT_NEAR(end.coords.im, lensed.im, 1e-6f);
}

/** @brief Paired projection seams use a broad under-fade and pixel over-fade. */
inline void test_shaderball_subduction_edge_fade() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();
  WB::FrameState frame = WB::frame(sb);
  frame.params.value.edge_width = 0.1f;
  frame.params.projection.coordinate_scale = 1.0f;
  WB::ProjectedLookup under{Complex(), 0,    0, WB::boundary_cut(),
                            0.01f,     1.0f, 0};
  WB::ProjectedLookup over = under;

  frame.slots.projection = WB::Projection::BONNE;
  over.region_id = 1;
  HS_EXPECT_EQ(WB::edge_fade_width(under, frame), 0.1f);
  HS_EXPECT_EQ(WB::edge_fade_width(over, frame), 1.0f / 20.0f);

  frame.slots.projection = WB::Projection::PEIRCE_QUINCUNCIAL;
  under.region_id = 2;
  over.region_id = 3;
  HS_EXPECT_EQ(WB::edge_fade_width(under, frame), 0.1f);
  HS_EXPECT_EQ(WB::edge_fade_width(over, frame), 1.0f / 20.0f);

  frame.slots.projection = WB::Projection::AIROCEAN;
  under.edge_class = 9;
  over.edge_class = 14;
  HS_EXPECT_EQ(WB::edge_fade_width(under, frame), 0.1f);
  HS_EXPECT_EQ(WB::edge_fade_width(over, frame), 1.0f / 20.0f);

  const uint8_t airocean_edge_pairs[][2] = {
      {9, 14},  {13, 37}, {16, 18}, {17, 41}, {20, 56}, {24, 29}, {38, 40},
      {42, 59}, {44, 45}, {47, 48}, {55, 58}, {62, 66}, {64, 67}};
  for (const auto &pair : airocean_edge_pairs) {
    HS_EXPECT_TRUE(shaderball::airocean_edge_is_under(pair[0]));
    HS_EXPECT_FALSE(shaderball::airocean_edge_is_under(pair[1]));
  }

  frame.slots.projection = WB::Projection::STEREOGRAPHIC;
  HS_EXPECT_EQ(WB::edge_fade_width(over, frame), 0.1f);
}

// A reference expression written out here and the implementation's own
// evaluation of it are two independent orderings, which -ffast-math is free to
// reassociate apart. sin/cos have unit derivative, so the result moves by at
// most the rounding of their argument; the widest argument these sweeps reach
// is SPIRAL's radius + 3*(azimuth + angle) + phase, under 24 radians.
constexpr float SOURCE_ARGUMENT_BOUND = 24.0f;
constexpr float SOURCE_DRIFT_BOUND =
    SOURCE_ARGUMENT_BOUND * std::numeric_limits<float>::epsilon();

/** @brief Legacy source functions retain their closed forms. */
inline void test_shaderball_legacy_sources() {
  using WB = ShaderBallWhiteBox;
  const float values[] = {-6.0f, -2.5f, -0.7f, 0.0f, 0.9f, 3.1f, 5.8f};
  const WB::SourceState source{0.9f, 1.1f, 0.7f, fast_cosf(0.7f),
                               fast_sinf(0.7f)};
  for (float re : values) {
    for (float im : values) {
      const Complex p(re, im);
      const float rotated = re * source.angle_cos + im * source.angle_sin;
      HS_EXPECT_EQ(WB::sample_function(WB::Function::TWIN_WAVE, p, source),
                   0.5f * (fast_sinf(re + source.primary) +
                           fast_sinf(rotated + source.primary)));
      HS_EXPECT_EQ(WB::sample_function(WB::Function::RINGS, p, source),
                   fast_sinf(sqrtf(re * re + im * im) - source.primary));
      const float radius = sqrtf(re * re + im * im);
      const float azimuth = fast_atan2(im, re);
      HS_EXPECT_NEAR(
          WB::sample_function(WB::Function::SPIRAL, p, source),
          fast_sinf(radius - 3.0f * (azimuth + source.angle) - source.primary),
          SOURCE_DRIFT_BOUND);
      const float b = -re * source.angle_sin + im * source.angle_cos;
      HS_EXPECT_NEAR(WB::sample_function(WB::Function::GRID, p, source),
                     fast_sinf(rotated + source.primary) *
                         fast_cosf(b - source.primary),
                     SOURCE_DRIFT_BOUND);
    }
  }
}

/** @brief Coupled/direct sampling reduces to both authored formulas. */
inline void test_shaderball_coupled_source() {
  using WB = ShaderBallWhiteBox;
  const float values[] = {-6.0f, -2.5f, -0.7f, 0.0f, 0.9f, 3.1f, 5.8f};
  const float phases[] = {0.0f, 0.8f, 2.4f, 4.9f};
  for (float re : values) {
    for (float im : values) {
      const Complex p(re, im);
      for (float primary : phases) {
        for (float secondary : phases) {
          for (float complexity : {0.5f, 1.7f, 3.0f}) {
            const float coupled =
                fast_sinf(re + complexity * fast_sinf(im + primary)) *
                fast_cosf(im + complexity * fast_cosf(re - secondary));
            HS_EXPECT_NEAR(
                WB::sample_pattern(p, complexity, 0.0f, primary, secondary),
                coupled, SOURCE_DRIFT_BOUND);
          }
          const float direct =
              fast_sinf(re + primary) * fast_cosf(im - secondary);
          HS_EXPECT_EQ(WB::sample_pattern(p, 0.0f, 1.0f, primary, secondary),
                       direct);
        }
      }
    }
  }
}

/** @brief Presets retain the legacy bank and add authored diagnostics. */
inline void test_shaderball_preset_bank() {
  using WB = ShaderBallWhiteBox;
  const auto &presets = WB::presets();
  const auto &choreo = WB::choreo();
  HS_EXPECT_EQ(presets.size(), size_t(21));
  HS_EXPECT_EQ(choreo.size(), presets.size());
  HS_EXPECT_EQ(presets[0].params.source.pattern_freq, 5.0f);
  HS_EXPECT_EQ(presets[0].params.warp.outer.scale, 3.0f);
  HS_EXPECT_EQ(presets[0].params.projection.pole_fade, 1.4f);
  HS_EXPECT_EQ(presets[0].params.surface_lens.mix, 1.0f);
  HS_EXPECT_EQ(presets[0].params.colorizer.breathe_depth, 0.15f);
  HS_EXPECT_EQ(presets[0].params.outer_camera.wander, 1.0f);
  HS_EXPECT_TRUE(choreo[0].staggered);
  HS_EXPECT_FALSE(choreo[6].staggered);
  HS_EXPECT_EQ(presets[15].slots.projection, WB::Projection::BONNE);
  HS_EXPECT_EQ(presets[16].slots.projection, WB::Projection::STEREOGRAPHIC);
  HS_EXPECT_EQ(presets[17].slots.projection, WB::Projection::STEREOGRAPHIC);
  HS_EXPECT_EQ(presets[18].slots.projection,
               WB::Projection::PEIRCE_QUINCUNCIAL);
  HS_EXPECT_EQ(presets[19].slots.projection, WB::Projection::AIROCEAN);
  HS_EXPECT_EQ(presets[17].slots.warp_program.outer.kind,
               WB::WarpStageKind::CURL_FLOW);
  HS_EXPECT_EQ(presets[17].slots.warp_program.outer.basis,
               WB::NoiseBasis::SIMPLEX);
  HS_EXPECT_EQ(presets[17].slots.warp_program.outer.curl_integrator,
               WB::CurlIntegrator::EULER_1);
  HS_EXPECT_EQ(presets[18].slots.warp_program.outer.kind,
               WB::WarpStageKind::NONE);
  HS_EXPECT_EQ(presets[18].slots.warp_program.inner.kind,
               WB::WarpStageKind::NONE);
  HS_EXPECT_EQ(presets[19].slots.surface_lens, WB::SurfaceLens::NONE);
  const auto &wave_shear = presets[20];
  HS_EXPECT_EQ(wave_shear.slots.function, WB::Function::COUPLED_DIRECT);
  HS_EXPECT_EQ(wave_shear.slots.projection, WB::Projection::STEREOGRAPHIC);
  HS_EXPECT_EQ(wave_shear.slots.projection_frame,
               WB::ProjectionFramePolicy::SPIN_WANDER);
  HS_EXPECT_EQ(wave_shear.slots.surface_lens, WB::SurfaceLens::GLITCH);
  HS_EXPECT_EQ(wave_shear.slots.warp_program.outer.kind,
               WB::WarpStageKind::WAVE_SHEAR);
  HS_EXPECT_EQ(wave_shear.slots.warp_program.inner.kind,
               WB::WarpStageKind::NONE);
  HS_EXPECT_EQ(wave_shear.slots.signal_weight, WB::SignalWeight::PROJECTION);
  HS_EXPECT_EQ(wave_shear.slots.value_transfer, WB::ValueTransfer::LINEAR);
  HS_EXPECT_EQ(wave_shear.slots.coverage,
               WB::CoveragePolicy::PROJECTION_WEIGHT_SQUARED);
  HS_EXPECT_EQ(wave_shear.slots.colorizer, WB::Colorizer::LIQUID);
  HS_EXPECT_EQ(wave_shear.params.source.pattern_freq, 4.439f);
  HS_EXPECT_EQ(wave_shear.params.source.speed, 0.245f);
  HS_EXPECT_EQ(wave_shear.params.source.angle_rate, 0.0f);
  HS_EXPECT_EQ(wave_shear.params.source.complexity, 0.5f);
  HS_EXPECT_EQ(wave_shear.params.source.pattern_mix, 0.0f);
  HS_EXPECT_EQ(wave_shear.params.source.secondary_rate, 0.0f);
  HS_EXPECT_EQ(wave_shear.params.projection.pole_fade, 1.0f);
  HS_EXPECT_EQ(wave_shear.params.projection.spin_rate, 0.0f);
  HS_EXPECT_EQ(wave_shear.params.projection.wander, 0.0f);
  HS_EXPECT_EQ(wave_shear.params.surface_lens.mix, 1.0f);
  HS_EXPECT_EQ(wave_shear.params.warp.outer.strength, 0.5f);
  HS_EXPECT_EQ(wave_shear.params.warp.outer.time_scale, 1.0f / 64.0f);
  HS_EXPECT_EQ(wave_shear.params.warp.outer.frequency, 1.0f);
  HS_EXPECT_EQ(wave_shear.params.warp.outer.field_angle, 0.0f);
  HS_EXPECT_EQ(wave_shear.params.colorizer.breathe_depth, 0.133f);
  HS_EXPECT_EQ(wave_shear.params.colorizer.cycle_speed, 0.05f);
  HS_EXPECT_EQ(wave_shear.params.colorizer.hue_shift, 0.0f);
  HS_EXPECT_EQ(wave_shear.params.colorizer.value_fade, 0.0f);
  HS_EXPECT_EQ(wave_shear.params.outer_camera.wander, 1.0f);
  for (size_t index = 0; index < presets.size(); ++index)
    HS_EXPECT_TRUE(WB::seam_compatible(presets[index]));
  for (size_t index = 0; index < 15; ++index)
    HS_EXPECT_EQ(presets[index].slots.projection,
                 WB::Projection::STEREOGRAPHIC);
  for (size_t index : {size_t(15), size_t(18), size_t(19)}) {
    HS_EXPECT_EQ(presets[index].slots.coverage, WB::CoveragePolicy::EDGE_FADE);
    HS_EXPECT_EQ(presets[index].params.value.edge_width, 0.1f);
  }
  HS_EXPECT_EQ(WB::device_cost(presets[15]).worst_case_points(), uint16_t(45));
  HS_EXPECT_EQ(WB::device_cost(presets[17]).worst_case_points(), uint16_t(43));
  HS_EXPECT_EQ(WB::device_cost(presets[18]).worst_case_points(), uint16_t(57));
  HS_EXPECT_EQ(WB::device_cost(presets[19]).worst_case_points(),
               WB::hold_device_point_budget());
  for (size_t index = 0; index < presets.size(); ++index) {
    const auto &preset = presets[index];
    if (index < 15)
      HS_EXPECT_TRUE(WB::slots_equal(preset.slots, WB::liquid_stereo_slots()));
    HS_EXPECT_TRUE(WB::valid_config(preset));
    if (index < 14)
      HS_EXPECT_TRUE(WB::slots_equal(preset.slots, presets[index + 1].slots));
  }

  reset_effect_globals();
  WB::SB sb;
  sb.init();
  HS_EXPECT_TRUE(
      WB::slots_equal(WB::active_slots(sb), WB::liquid_stereo_slots()));
}

/** @brief Whole-schema validation applies valid configs and rejects invalid. */
inline void test_shaderball_config_admission() {
  using WB = ShaderBallWhiteBox;
  {
    reset_effect_globals();
    WB::SB sb;
    sb.init();
    const WB::Slots original = WB::active_slots(sb);
    WB::RequestedConfig invalid_params = WB::active_config(sb);
    invalid_params.params.source.pattern_freq = 21.0f;
    const WB::RequestedConfig before_invalid_params = WB::active_config(sb);
    WB::request_config(sb, invalid_params);
    HS_EXPECT_TRUE(WB::active_config(sb) == before_invalid_params);
    HS_EXPECT_TRUE(WB::requested_slots(sb) == original);
    HS_EXPECT_FALSE(WB::transition_active(sb));
    HS_EXPECT_FALSE(WB::param_morph_active(sb));

    WB::RequestedConfig candidate = WB::legacy_config();
    const auto invalid_tag = static_cast<uint8_t>(0xff);
    candidate.slots.function = static_cast<WB::Function>(invalid_tag);
    HS_EXPECT_FALSE(WB::valid_config(candidate));
    candidate = WB::legacy_config();
    candidate.slots.projection = static_cast<WB::Projection>(invalid_tag);
    HS_EXPECT_FALSE(WB::valid_config(candidate));
    candidate = WB::legacy_config();
    candidate.slots.projection_frame =
        static_cast<WB::ProjectionFramePolicy>(invalid_tag);
    HS_EXPECT_FALSE(WB::valid_config(candidate));
    candidate = WB::legacy_config();
    candidate.slots.surface_lens = static_cast<WB::SurfaceLens>(invalid_tag);
    HS_EXPECT_FALSE(WB::valid_config(candidate));
    candidate = WB::legacy_config();
    candidate.slots.warp_program.outer.kind =
        static_cast<WB::WarpStageKind>(invalid_tag);
    HS_EXPECT_FALSE(WB::valid_config(candidate));
    candidate = WB::legacy_config();
    candidate.slots.warp_program.inner.kind =
        static_cast<WB::WarpStageKind>(invalid_tag);
    HS_EXPECT_FALSE(WB::valid_config(candidate));
    candidate = WB::legacy_config();
    candidate.slots.signal_weight = static_cast<WB::SignalWeight>(invalid_tag);
    HS_EXPECT_FALSE(WB::valid_config(candidate));
    candidate = WB::legacy_config();
    candidate.slots.value_transfer =
        static_cast<WB::ValueTransfer>(invalid_tag);
    HS_EXPECT_FALSE(WB::valid_config(candidate));
    candidate = WB::legacy_config();
    candidate.slots.coverage = static_cast<WB::CoveragePolicy>(invalid_tag);
    HS_EXPECT_FALSE(WB::valid_config(candidate));
    candidate = WB::legacy_config();
    candidate.slots.colorizer = static_cast<WB::Colorizer>(invalid_tag);
    HS_EXPECT_FALSE(WB::valid_config(candidate));

    const WB::RequestedConfig legacy_config = WB::legacy_config();
    HS_EXPECT_TRUE(WB::valid_config(legacy_config));
    HS_EXPECT_TRUE(
        WB::transition_admitted(WB::active_config(sb), legacy_config));
    WB::request_config(sb, legacy_config);
    HS_EXPECT_TRUE(WB::active_config(sb) == legacy_config);
    HS_EXPECT_FALSE(WB::transition_active(sb));
    HS_EXPECT_TRUE(WB::requested_slots(sb) == legacy_config.slots);

    WB::RequestedConfig resource_from = WB::legacy_config();
    resource_from.slots.function = WB::Function::NOISE_CONTOUR;
    resource_from.params.source.noise_basis = WB::NoiseBasis::SIMPLEX;
    WB::RequestedConfig resource_to = resource_from;
    resource_to.params.source.noise_basis = WB::NoiseBasis::FBM3;
    HS_EXPECT_TRUE(WB::valid_config(resource_from));
    HS_EXPECT_TRUE(WB::valid_config(resource_to));
    HS_EXPECT_TRUE(WB::transition_admitted(resource_from, resource_to));
    resource_to.params.source.noise_resource_id ^= 4U;
    HS_EXPECT_TRUE(WB::transition_admitted(resource_from, resource_to));

    WB::RequestedConfig shared_owner = WB::legacy_config();
    shared_owner.slots.warp_program.outer.kind =
        WB::WarpStageKind::VECTOR_NOISE;
    shared_owner.slots.warp_program.inner.kind =
        WB::WarpStageKind::VECTOR_NOISE;
    shared_owner.params.warp.outer.time_scale = 0.0f;
    shared_owner.params.warp.inner.time_scale = 0.0f;
    shared_owner.slots.warp_program.outer.resource_id = 6;
    shared_owner.slots.warp_program.inner.resource_id = 6;
    HS_EXPECT_TRUE(WB::valid_config(shared_owner));
    shared_owner.slots.warp_program.inner.basis = WB::NoiseBasis::FBM3;
    HS_EXPECT_FALSE(WB::valid_config(shared_owner));
    HS_EXPECT_FALSE(WB::transition_admitted(shared_owner, shared_owner));
  }

  {
    reset_effect_globals();
    WB::SB projection_change;
    projection_change.init();
    HS_EXPECT_TRUE(
        projection_change.updateParameter(
            "Projection", static_cast<float>(WB::Projection::BONNE)) ==
        ParamSetResult::APPLIED);
    projection_change.draw_frame();
    projection_change.advance_display();
    HS_EXPECT_FALSE(WB::transition_active(projection_change));
    HS_EXPECT_EQ(WB::active_slots(projection_change).projection,
                 WB::Projection::BONNE);
    HS_EXPECT_EQ(WB::active_slots(projection_change).coverage,
                 WB::CoveragePolicy::EDGE_FADE);
    HS_EXPECT_EQ(WB::active_slots(projection_change).warp_program.outer.kind,
                 WB::WarpStageKind::NONE);
  }

  {
    reset_effect_globals();
    WB::SB legacy_warp_change;
    legacy_warp_change.init();
    HS_EXPECT_TRUE(
        legacy_warp_change.updateParameter(
            "Projection", static_cast<float>(WB::Projection::BONNE)) ==
        ParamSetResult::APPLIED);
    legacy_warp_change.draw_frame();
    legacy_warp_change.advance_display();
    HS_EXPECT_FALSE(WB::transition_active(legacy_warp_change));
    HS_EXPECT_TRUE(
        legacy_warp_change.updateParameter(
            "Outer Warp",
            static_cast<float>(WB::WarpStageKind::LEGACY_STEREO_NOISE)) ==
        ParamSetResult::APPLIED);
    legacy_warp_change.draw_frame();
    legacy_warp_change.advance_display();
    HS_EXPECT_FALSE(WB::transition_active(legacy_warp_change));
    HS_EXPECT_EQ(WB::active_slots(legacy_warp_change).projection,
                 WB::Projection::STEREOGRAPHIC);
  }
}

/** @brief Selector intent is independent of the live worker destination. */
inline void test_shaderball_deterministic_gui_edits() {
  using WB = ShaderBallWhiteBox;
  struct Result {
    std::array<WB::RequestedConfig, 3> configs;
    std::array<uint64_t, 3> schema_hashes{};
    std::array<uint32_t, 3> generation_deltas{};
  };
  auto schema_hash = [](const WB::SB &sb) {
    uint64_t hash = 1469598103934665603ULL;
    auto append = [&](const char *text) {
      for (; *text != '\0'; ++text) {
        hash ^= static_cast<uint8_t>(*text);
        hash *= 1099511628211ULL;
      }
    };
    for (const Effect::ParamDef &def : sb.getParameters()) {
      append(def.name);
      hash ^= static_cast<uint32_t>(def.option_count);
      hash *= 1099511628211ULL;
      for (int option = 0; option < def.option_count; ++option)
        append(def.options[option]);
    }
    return hash;
  };
  auto run = [&](bool advance_worker) {
    reset_effect_globals();
    WB::SB sb;
    sb.init();
    Result result;
    const char *names[] = {"Projection", "Projection", "Outer Warp"};
    const float values[] = {
        static_cast<float>(WB::Projection::BONNE),
        static_cast<float>(WB::Projection::AIROCEAN),
        static_cast<float>(WB::WarpStageKind::LEGACY_STEREO_NOISE)};
    for (size_t index = 0; index < result.configs.size(); ++index) {
      const uint32_t before = sb.getParameterSchemaGeneration();
      HS_EXPECT_TRUE(sb.updateParameter(names[index], values[index]) ==
                     ParamSetResult::APPLIED);
      result.generation_deltas[index] =
          sb.getParameterSchemaGeneration() - before;
      result.configs[index] = WB::requested_config(sb);
      result.schema_hashes[index] = schema_hash(sb);
      if (advance_worker && index + 1 < result.configs.size()) {
        sb.draw_frame();
        sb.advance_display();
      }
    }
    return result;
  };

  const Result idle = run(false);
  const Result worker = run(true);
  for (size_t index = 0; index < idle.configs.size(); ++index) {
    HS_EXPECT_TRUE(idle.configs[index] == worker.configs[index]);
    HS_EXPECT_EQ(idle.schema_hashes[index], worker.schema_hashes[index]);
    HS_EXPECT_EQ(idle.generation_deltas[index],
                 worker.generation_deltas[index]);
  }
  HS_EXPECT_EQ(idle.configs[2].slots.projection, WB::Projection::STEREOGRAPHIC);
  HS_EXPECT_EQ(idle.configs[2].slots.warp_program.outer.kind,
               WB::WarpStageKind::LEGACY_STEREO_NOISE);
}

/** @brief Calibrated device cost admits measured-safe operator frontiers. */
inline void test_shaderball_work_admission() {
  using WB = ShaderBallWhiteBox;
  const auto &presets = WB::presets();
  for (const auto &preset : presets)
    HS_EXPECT_TRUE(WB::hold_admitted(preset));

  HS_EXPECT_EQ(WB::noise_call_points(), uint16_t(4));
  HS_EXPECT_EQ(WB::cost_tier(presets[15]), WB::CostTier::T2);
  HS_EXPECT_EQ(WB::cost_tier(presets[17]), WB::CostTier::T2);
  HS_EXPECT_EQ(WB::cost_tier(presets[18]), WB::CostTier::T3);
  HS_EXPECT_EQ(WB::cost_tier(presets[19]), WB::CostTier::T3);

  WB::RequestedConfig integrated_ridged = presets[17];
  integrated_ridged.slots.warp_program.outer.basis = WB::NoiseBasis::RIDGED3;
  integrated_ridged.slots.warp_program.outer.curl_integrator =
      WB::CurlIntegrator::MIDPOINT_2;
  HS_EXPECT_TRUE(WB::valid_config(integrated_ridged));
  HS_EXPECT_EQ(WB::device_cost(integrated_ridged).worst_case_noise_calls(),
               uint16_t(96));
  HS_EXPECT_EQ(WB::device_cost(integrated_ridged).worst_case_points(),
               uint16_t(395));
  HS_EXPECT_FALSE(WB::hold_admitted(integrated_ridged));

  WB::RequestedConfig peirce_polar = presets[18];
  peirce_polar.slots.warp_program.inner.kind = WB::WarpStageKind::POLAR_CHART;
  peirce_polar.slots.warp_program.inner.polar_harmonic = 2;
  peirce_polar.params.source.pattern_freq = 1.0f;
  HS_EXPECT_TRUE(WB::valid_config(peirce_polar));
  HS_EXPECT_EQ(WB::device_cost(peirce_polar).worst_case_points(), uint16_t(65));
  HS_EXPECT_FALSE(WB::hold_admitted(peirce_polar));

  WB::RequestedConfig airocean_mobius = presets[19];
  airocean_mobius.slots.surface_lens = WB::SurfaceLens::MOBIUS;
  airocean_mobius.params.surface_lens.mix = 1.0f;
  HS_EXPECT_TRUE(WB::valid_config(airocean_mobius));
  HS_EXPECT_EQ(WB::device_cost(airocean_mobius).worst_case_points(),
               uint16_t(73));
  HS_EXPECT_FALSE(WB::hold_admitted(airocean_mobius));

  for (WB::Projection projection :
       {WB::Projection::BONNE, WB::Projection::PEIRCE_QUINCUNCIAL,
        WB::Projection::AIROCEAN}) {
    WB::RequestedConfig strict = WB::legacy_config();
    strict.slots.projection = projection;
    strict.slots.warp_program.outer.kind = WB::WarpStageKind::NONE;
    strict.slots.surface_lens = WB::SurfaceLens::NONE;
    strict.params.surface_lens.mix = 0.0f;
    HS_EXPECT_TRUE(WB::hold_admitted(strict));
  }
  WB::RequestedConfig airo_mobius = WB::legacy_config();
  airo_mobius.slots.projection = WB::Projection::AIROCEAN;
  airo_mobius.slots.warp_program.outer.kind = WB::WarpStageKind::NONE;
  airo_mobius.slots.surface_lens = WB::SurfaceLens::MOBIUS;
  airo_mobius.params.surface_lens.mix = 0.5f;
  HS_EXPECT_TRUE(WB::valid_config(airo_mobius));
  HS_EXPECT_EQ(WB::device_cost(airo_mobius).noise_calls, uint16_t(0));
  HS_EXPECT_GT(WB::device_cost(airo_mobius).worst_case_points(),
               WB::hold_device_point_budget());
  HS_EXPECT_FALSE(WB::hold_admitted(airo_mobius));

  WB::RequestedConfig from = WB::legacy_config();
  from.slots.surface_lens = WB::SurfaceLens::NONE;
  from.params.surface_lens.mix = 0.0f;
  from.slots.warp_program.outer.kind = WB::WarpStageKind::CURL_FLOW;
  from.slots.warp_program.outer.basis = WB::NoiseBasis::SIMPLEX;
  from.slots.warp_program.outer.curl_integrator = WB::CurlIntegrator::EULER_1;
  from.params.warp.outer.strength = 0.0078125f;
  from.params.warp.outer.scale = 1.0f;
  from.params.warp.outer.time_scale = 0.0f;
  WB::RequestedConfig to = from;
  to.params.warp.outer.strength = 0.5f;
  to.params.warp.outer.scale = 1.0f / 64.0f;
  HS_EXPECT_TRUE(WB::hold_admitted(from));
  HS_EXPECT_TRUE(WB::hold_admitted(to));
  HS_EXPECT_TRUE(WB::transition_admitted(from, to));
  HS_EXPECT_FALSE(WB::stable_topology(from, to));
  HS_EXPECT_FALSE(WB::stable_parameter_path_admitted(from, to));
  for (size_t index = 0; index < presets.size(); ++index) {
    const auto &next = presets[(index + 1) % presets.size()];
    if (WB::stable_topology(presets[index], next))
      HS_EXPECT_TRUE(WB::stable_parameter_path_admitted(presets[index], next));
  }
}

/** @brief Folded and interrupted projections reject unproved noise seams. */
inline void test_shaderball_strict_seam_admission() {
  using WB = ShaderBallWhiteBox;
  for (WB::Projection projection :
       {WB::Projection::BONNE, WB::Projection::PEIRCE_QUINCUNCIAL,
        WB::Projection::AIROCEAN}) {
    WB::RequestedConfig config = WB::legacy_config();
    config.slots.projection = projection;
    config.slots.surface_lens = WB::SurfaceLens::NONE;
    config.params.surface_lens.mix = 0.0f;
    config.slots.warp_program.outer.kind = WB::WarpStageKind::NONE;
    HS_EXPECT_TRUE(WB::valid_config(config));

    config.slots.function = WB::Function::NOISE_CONTOUR;
    HS_EXPECT_FALSE(WB::valid_config(config));
    config.slots.function = WB::Function::COUPLED_DIRECT;

    config.slots.warp_program.outer.kind = WB::WarpStageKind::VECTOR_NOISE;
    config.params.warp.outer.scale = 1.0f;
    config.params.warp.outer.strength = 0.1f;
    config.params.warp.outer.time_scale = 0.0f;
    HS_EXPECT_FALSE(WB::valid_config(config));
    config.slots.warp_program.outer.kind = WB::WarpStageKind::CURL_FLOW;
    config.params.warp.outer.scale = 1.0f;
    config.params.warp.outer.strength = 1e-4f;
    HS_EXPECT_FALSE(WB::valid_config(config));
    config.slots.warp_program.outer.kind = WB::WarpStageKind::NONE;

    config.slots.surface_lens = WB::SurfaceLens::TANGENT_NOISE;
    config.params.surface_lens.mix = 0.5f;
    HS_EXPECT_FALSE(WB::valid_config(config));

    config.slots.projection = WB::Projection::EQUIRECTANGULAR;
    HS_EXPECT_TRUE(WB::valid_config(config));
  }
}

/** @brief Additive warp metadata retains sub-ULP displacement at large coordinates. */
inline void test_shaderball_additive_delta_precision() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();
  WB::FrameState frame = WB::frame(sb);
  frame.clocks.warp_outer_phase = 0.25f;
  WB::WarpStageSpec spec{WB::WarpStageKind::WAVE_SHEAR};
  WB::WarpStageParams params;
  params.strength = 0.001f;
  params.frequency = 0.0f;
  params.field_angle = 0.0f;
  const Complex input(32768.0f, 32768.0f);
  const WB::ProjectedLookup projected{input, 0, 0, 0, 1.0f, 1.0f, 0};
  const auto result = WB::warp_stage(input, projected, spec, params, frame);
  HS_EXPECT_NEAR(result.delta.re, 0.0f, 1e-8f);
  HS_EXPECT_NEAR(result.delta.im, 0.001f, 1e-7f);
  HS_EXPECT_NEAR(result.deformation, 0.001f, 1e-7f);
  HS_EXPECT_EQ(result.coords.im, input.im);
}

/** @brief Profiling can land every curated hold without choreography. */
inline void test_shaderball_profile_presets() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();
  const auto &presets = WB::presets();
  for (size_t index = 0; index < presets.size(); ++index) {
    sb.profile_select_preset(index);
    HS_EXPECT_TRUE(WB::active_config(sb) == presets[index]);
    HS_EXPECT_TRUE(WB::requested_config(sb) == presets[index]);
    HS_EXPECT_TRUE(WB::hold_admitted(WB::active_config(sb)));
    HS_EXPECT_FALSE(WB::transition_active(sb));
    HS_EXPECT_FALSE(WB::param_morph_active(sb));
    if (index == 15 || index == 18 || index == 19) {
      const auto projected = WB::surface_project(
          Vector(0.808122f, -0.303046f, 0.505076f), WB::frame(sb));
      HS_EXPECT_TRUE(std::isfinite(projected.fade_edge_distance));
    }
  }
}

/** @brief Every GUI enum option is writable and survives its handoff. */
inline void test_shaderball_gui_catalog() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();
  HS_EXPECT_LE(sb.getParameters().size(), size_t(64));
  const auto *projection = sb.getParameters().find("Projection");
  HS_EXPECT_TRUE(projection != nullptr);
  HS_EXPECT_EQ(projection->option_count, 6);
  HS_EXPECT_TRUE(std::strcmp(projection->options[3], "Bonne") == 0);
  HS_EXPECT_TRUE(std::strcmp(projection->options[4], "Peirce Quincuncial") ==
                 0);
  HS_EXPECT_TRUE(std::strcmp(projection->options[5], "Dymaxion / Airocean") ==
                 0);
  const auto *coverage = sb.getParameters().find("Coverage");
  HS_EXPECT_TRUE(coverage != nullptr);
  HS_EXPECT_EQ(coverage->option_count, 5);
  HS_EXPECT_TRUE(std::strcmp(coverage->options[4], "Projection Weight") == 0);
  HS_EXPECT_TRUE(
      sb.updateParameter("Function",
                         static_cast<float>(WB::Function::PRIMITIVE_LATTICE)) ==
      ParamSetResult::APPLIED);
  const auto *lattice_softness = sb.getParameters().find("Lattice Softness");
  HS_EXPECT_TRUE(lattice_softness != nullptr);
  HS_EXPECT_EQ(lattice_softness->max, 1.0f);
  HS_EXPECT_TRUE(sb.updateParameter("Lattice Softness", 1.0f) ==
                 ParamSetResult::APPLIED);
  sb.draw_frame();
  sb.advance_display();
  WB::settle_transition(sb);
  HS_EXPECT_EQ(WB::active_config(sb).params.source.lattice_softness, 1.0f);
  const uint32_t schema_before = sb.getParameterSchemaGeneration();
  HS_EXPECT_TRUE(sb.updateParameter(
                     "Projection", static_cast<float>(WB::Projection::BONNE)) ==
                 ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(sb.getParameterSchemaGeneration() > schema_before);
  HS_EXPECT_TRUE(sb.getParameters().find("Bonne Hemisphere") != nullptr);
  sb.draw_frame();
  sb.advance_display();
  WB::settle_transition(sb);

  WB::RequestedConfig gui_base = WB::legacy_config();
  gui_base.slots.function = WB::Function::COUPLED_DIRECT;
  gui_base.slots.surface_lens = WB::SurfaceLens::NONE;
  gui_base.params.surface_lens.mix = 0.0f;
  gui_base.slots.warp_program.outer.kind = WB::WarpStageKind::NONE;
  gui_base.slots.warp_program.inner.kind = WB::WarpStageKind::NONE;
  gui_base.params.source.pattern_freq = 1.0f;
  auto reset_gui = [&] {
    WB::request_config(sb, gui_base);
    WB::settle_transition(sb);
  };

  constexpr const char *ROOT_ENUMS[] = {
      "Function",   "Projection", "Projection Frame", "Lens",
      "Outer Warp", "Inner Warp", "Signal Weight",    "Value Transfer",
      "Coverage",   "Colorizer"};
  for (const char *name : ROOT_ENUMS) {
    const int option_count = sb.getParameters().find(name)->option_count;
    for (int option = 0; option < option_count; ++option) {
      reset_gui();
      HS_EXPECT_TRUE(sb.updateParameter(name, static_cast<float>(option)) ==
                     ParamSetResult::APPLIED);
      sb.draw_frame();
      sb.advance_display();
      WB::settle_transition(sb);
      HS_EXPECT_EQ(sb.getParameters().find(name)->get(),
                   static_cast<float>(option));
      HS_EXPECT_FALSE(WB::transition_active(sb));
      HS_EXPECT_FALSE(WB::param_morph_active(sb));
    }
  }

  auto select_and_set_all = [&](const char *root, int selection,
                                const char *subordinate) {
    reset_gui();
    HS_EXPECT_TRUE(sb.updateParameter(root, static_cast<float>(selection)) ==
                   ParamSetResult::APPLIED);
    sb.draw_frame();
    sb.advance_display();
    WB::settle_transition(sb);
    const auto *def = sb.getParameters().find(subordinate);
    HS_EXPECT(def != nullptr, subordinate);
    if (def == nullptr)
      return;
    // An enum's admissible values are its option indices; an integer param
    // carries a range instead, so sweep whichever the target declares.
    const int first = def->option_count > 0 ? 0 : static_cast<int>(def->min);
    const int last = def->option_count > 0 ? def->option_count - 1
                                           : static_cast<int>(def->max);
    for (int value = first; value <= last; ++value) {
      HS_EXPECT_TRUE(
          sb.updateParameter(subordinate, static_cast<float>(value)) ==
          ParamSetResult::APPLIED);
      sb.draw_frame();
      sb.advance_display();
      WB::settle_transition(sb);
      HS_EXPECT_EQ(sb.getParameters().find(subordinate)->get(),
                   static_cast<float>(value));
    }
  };
  auto select_and_reject_expensive_options = [&](const char *root,
                                                 int selection,
                                                 const char *subordinate) {
    reset_gui();
    HS_EXPECT_TRUE(sb.updateParameter(root, static_cast<float>(selection)) ==
                   ParamSetResult::APPLIED);
    sb.draw_frame();
    sb.advance_display();
    WB::settle_transition(sb);
    const auto *def = sb.getParameters().find(subordinate);
    HS_EXPECT(def != nullptr, subordinate);
    if (def == nullptr)
      return;
    HS_EXPECT_EQ(def->option_count, 3);
    for (int option = 0; option < def->option_count; ++option) {
      HS_EXPECT_TRUE(
          sb.updateParameter(subordinate, static_cast<float>(option)) ==
          ParamSetResult::APPLIED);
      sb.draw_frame();
      sb.advance_display();
      WB::settle_transition(sb);
      HS_EXPECT_EQ(sb.getParameters().find(subordinate)->get(), 0.0f);
    }
  };
  select_and_set_all("Function", 5, "Source Noise Basis");
  select_and_set_all("Projection", 2, "Gnomonic Hemisphere");
  select_and_set_all("Projection", 3, "Bonne Hemisphere");
  HS_EXPECT_TRUE(sb.updateParameter("Bonne Standard Parallel", 0.9f) ==
                 ParamSetResult::APPLIED);
  sb.draw_frame();
  sb.advance_display();
  WB::settle_transition(sb);
  HS_EXPECT_NEAR(
      WB::active_config(sb).params.projection.bonne_standard_parallel, 0.9f,
      1e-6f);
  select_and_set_all("Projection", 4, "Peirce Layout");
  for (int layout = 0; layout < 4; ++layout) {
    HS_EXPECT_TRUE(
        sb.updateParameter("Peirce Layout", static_cast<float>(layout)) ==
        ParamSetResult::APPLIED);
    const bool has_scroll =
        sb.getParameters().find("Projection Layout Scroll") != nullptr;
    HS_EXPECT_EQ(has_scroll, layout >= 2);
    sb.draw_frame();
    sb.advance_display();
    WB::settle_transition(sb);
  }
  select_and_set_all("Projection", 5, "Airocean Layout");
  select_and_reject_expensive_options("Outer Warp", 5, "Outer Noise Basis");
  select_and_set_all("Outer Warp", 5, "Outer Warp Envelope");
  select_and_reject_expensive_options("Outer Warp", 6, "Outer Curl Integrator");
  {
    reset_gui();
    HS_EXPECT_TRUE(
        sb.updateParameter("Outer Warp",
                           static_cast<float>(WB::WarpStageKind::CURL_FLOW)) ==
        ParamSetResult::APPLIED);
    sb.draw_frame();
    sb.advance_display();
    WB::settle_transition(sb);
    const auto *strength = sb.getParameters().find("Outer Warp Strength");
    HS_EXPECT_TRUE(strength != nullptr);
    // 0.5 / (scale 0.1 * gradient bound 64), one Euler interval.
    const float limit = strength->max;
    HS_EXPECT_NEAR(limit, 0.078125f, 1e-6f);
    HS_EXPECT_EQ(strength->min, -limit);
    HS_EXPECT_TRUE(sb.updateParameter("Outer Warp Strength", 4.0f) ==
                   ParamSetResult::APPLIED);
    sb.draw_frame();
    sb.advance_display();
    WB::settle_transition(sb);
    HS_EXPECT_EQ(WB::active_config(sb).params.warp.outer.strength, limit);
  }
  select_and_set_all("Outer Warp", 8, "Outer Polar Mode");
  select_and_set_all("Outer Warp", 8, "Outer Polar Harmonic");
  HS_EXPECT_TRUE(sb.getParameters().find("Pattern Freq") == nullptr);
  HS_EXPECT_TRUE(
      sb.updateParameter("Function", static_cast<float>(WB::Function::RINGS)) ==
      ParamSetResult::APPLIED);
  HS_EXPECT_EQ(WB::requested_config(sb).slots.warp_program.outer.kind,
               WB::WarpStageKind::NONE);
  HS_EXPECT_TRUE(sb.getParameters().find("Pattern Freq") != nullptr);
  select_and_set_all("Value Transfer", 3, "Band Count");
}

/** @brief New cartographic kernels preserve landmarks and stay finite. */
inline void test_shaderball_projection_catalog() {
  using WB = ShaderBallWhiteBox;
  const float standard_parallel = PI_F * 0.25f;
  const Vector bonne_origin(cosf(standard_parallel), sinf(standard_parallel),
                            0.0f);
  const auto bonne =
      shaderball::bonne_projection(bonne_origin, 0.0f, standard_parallel);
  HS_EXPECT_NEAR(bonne.coords.re, 0.0f, 2e-5f);
  HS_EXPECT_NEAR(bonne.coords.im, 0.0f, 2e-5f);

  const auto peirce = shaderball::peirce_projection(UP, 0.0f, 1, 0.0f);
  HS_EXPECT_NEAR(peirce.coords.re, 0.0f, 2e-5f);
  HS_EXPECT_NEAR(peirce.coords.im, 0.0f, 2e-5f);

  const auto &center = shaderball::AIROCEAN_CENTERS[0];
  const auto airocean = shaderball::airocean_projection(
      Vector(center.x, center.z, center.y), 0.0f, false);
  const auto &triangle = shaderball::AIROCEAN_PLANAR_FACES[0];
  HS_EXPECT_NEAR(airocean.coords.re,
                 (triangle[0].x + triangle[1].x + triangle[2].x) / 3.0f, 2e-4f);
  HS_EXPECT_NEAR(airocean.coords.im,
                 (triangle[0].y + triangle[1].y + triangle[2].y) / 3.0f, 2e-4f);
  HS_EXPECT_EQ(airocean.region_id, uint8_t(0));

  const Vector equator_zero(1.0f, 0.0f, 0.0f);
  const Vector equator_east(0.0f, 0.0f, 1.0f);
  const Vector antimeridian(-1.0f, 0.0f, 0.0f);
  const auto bonne_equator =
      shaderball::bonne_projection(equator_zero, 0.0f, standard_parallel);
  HS_EXPECT_NEAR(bonne_equator.coords.re, 0.0f, 2e-5f);
  HS_EXPECT_NEAR(bonne_equator.coords.im, -0.7853981634f, 2e-5f);
  const auto bonne_east =
      shaderball::bonne_projection(equator_east, 0.0f, standard_parallel);
  HS_EXPECT_NEAR(bonne_east.coords.re, 1.3758501640f, 3e-5f);
  HS_EXPECT_NEAR(bonne_east.coords.im, -0.1378413458f, 3e-5f);
  const auto bonne_shifted = shaderball::bonne_projection(
      equator_east, 0.5f * PI_F, standard_parallel);
  HS_EXPECT_NEAR(bonne_shifted.coords.re, 0.0f, 2e-5f);
  HS_EXPECT_NEAR(bonne_shifted.coords.im, -0.7853981634f, 2e-5f);
  const auto bonne_cut =
      shaderball::bonne_projection(antimeridian, 0.0f, standard_parallel);
  HS_EXPECT_NEAR(bonne_cut.fade_edge_distance, 0.0f, 2e-5f);
  const auto werner =
      shaderball::bonne_projection(equator_east, 0.0f, 0.5f * PI_F);
  HS_EXPECT_NEAR(werner.coords.re, 1.3217795320f, 3e-5f);
  HS_EXPECT_NEAR(werner.coords.im, -0.8487048774f, 3e-5f);

  constexpr float PEIRCE_K = 1.8540746773013719f;
  const Vector south_pole(0.0f, -1.0f, 0.0f);
  const auto peirce_diamond =
      shaderball::peirce_projection(south_pole, 0.0f, 0, 0.0f);
  const auto peirce_square =
      shaderball::peirce_projection(south_pole, 0.0f, 1, 0.0f);
  const auto peirce_horizontal =
      shaderball::peirce_projection(south_pole, 0.0f, 2, 0.0f);
  const auto peirce_vertical =
      shaderball::peirce_projection(south_pole, 0.0f, 3, 0.0f);
  HS_EXPECT_NEAR(peirce_diamond.coords.re, 2.0f * PEIRCE_K, 3e-5f);
  HS_EXPECT_NEAR(peirce_diamond.coords.im, 0.0f, 3e-5f);
  HS_EXPECT_NEAR(peirce_square.coords.re, PEIRCE_K * 1.4142135624f, 3e-5f);
  HS_EXPECT_NEAR(peirce_square.coords.im, PEIRCE_K * 1.4142135624f, 3e-5f);
  HS_EXPECT_NEAR(peirce_horizontal.coords.re, PEIRCE_K, 3e-5f);
  HS_EXPECT_NEAR(peirce_vertical.coords.im, PEIRCE_K, 3e-5f);
  const auto peirce_scroll0 =
      shaderball::peirce_projection(equator_east, 0.0f, 2, 0.0f);
  const auto peirce_scroll1 =
      shaderball::peirce_projection(equator_east, 0.0f, 2, 1.0f);
  HS_EXPECT_NEAR(peirce_scroll0.coords.re, peirce_scroll1.coords.re, 3e-5f);
  HS_EXPECT_NEAR(peirce_scroll0.coords.im, peirce_scroll1.coords.im, 3e-5f);
  WB::ProjectionParams scroll_from;
  WB::ProjectionParams scroll_to;
  scroll_from.layout_scroll = 0.9f;
  scroll_to.layout_scroll = -0.9f;
  const auto scroll_start = WB::lerp_projection(scroll_from, scroll_to, 0.0f);
  const auto scroll_mid = WB::lerp_projection(scroll_from, scroll_to, 0.5f);
  const auto scroll_end = WB::lerp_projection(scroll_from, scroll_to, 1.0f);
  HS_EXPECT_EQ(scroll_start.layout_scroll, 0.9f);
  HS_EXPECT_NEAR(scroll_mid.layout_scroll, 1.0f, 1e-6f);
  HS_EXPECT_EQ(scroll_end.layout_scroll, -0.9f);
  const auto peirce_scroll_mid = shaderball::peirce_projection(
      equator_east, 0.0f, 2, scroll_mid.layout_scroll);
  HS_EXPECT_NEAR(peirce_scroll_mid.coords.re, peirce_scroll0.coords.re, 3e-5f);
  HS_EXPECT_NEAR(peirce_scroll_mid.coords.im, peirce_scroll0.coords.im, 3e-5f);
  const auto peirce_zero_fade =
      shaderball::peirce_projection(equator_zero, 0.0f, 1, 0.0f);
  HS_EXPECT_NEAR(peirce_zero_fade.fade_edge_distance, 0.0f, 2e-5f);

  auto lon_lat = [](float longitude, float latitude) {
    const float cp = cosf(latitude);
    return Vector(cp * cosf(longitude), sinf(latitude), cp * sinf(longitude));
  };
  const Vector peirce_oracle_point =
      lon_lat(23.0f * PI_F / 180.0f, 28.0f * PI_F / 180.0f);
  const Complex peirce_oracles[] = {
      {0.4550257621f, -1.1120261272f},
      {1.1080730174f, -0.4645694134f},
      {-0.4349300830f, -1.1120261272f},
      {0.4550257621f, -2.0019819724f},
  };
  for (uint8_t layout = 0; layout < 4; ++layout) {
    const auto mapped =
        shaderball::peirce_projection(peirce_oracle_point, 0.0f, layout, 0.13f);
    HS_EXPECT_NEAR(mapped.coords.re, peirce_oracles[layout].re, 3e-5f);
    HS_EXPECT_NEAR(mapped.coords.im, peirce_oracles[layout].im, 3e-5f);
  }
  constexpr float TIE_EPS = 1e-5f;
  constexpr float SOUTH_LATITUDE = -0.4f;
  struct SectorTie {
    float longitude;
    uint8_t before;
    uint8_t exact;
  };
  const SectorTie sector_ties[] = {
      {-0.75f * PI_F, 1, 2},
      {-0.25f * PI_F, 2, 3},
      {0.25f * PI_F, 3, 4},
      {0.75f * PI_F, 4, 1},
  };
  for (const auto &tie : sector_ties) {
    const auto before = shaderball::peirce_projection(
        lon_lat(tie.longitude - TIE_EPS, SOUTH_LATITUDE), 0.0f, 0, 0.0f);
    const auto exact = shaderball::peirce_projection(
        lon_lat(tie.longitude, SOUTH_LATITUDE), 0.0f, 0, 0.0f);
    HS_EXPECT_EQ(before.region_id, tie.before);
    HS_EXPECT_EQ(exact.region_id, tie.exact);
  }

  constexpr float MERIDIAN = 0.37f;
  constexpr float MERIDIAN_LATITUDE = 0.28f;
  const Vector on_meridian = lon_lat(MERIDIAN, MERIDIAN_LATITUDE);
  const Vector zero_meridian = lon_lat(0.0f, MERIDIAN_LATITUDE);
  for (uint8_t layout = 0; layout < 4; ++layout) {
    const auto shifted =
        shaderball::peirce_projection(on_meridian, MERIDIAN, layout, 0.13f);
    const auto reference =
        shaderball::peirce_projection(zero_meridian, 0.0f, layout, 0.13f);
    HS_EXPECT_NEAR(shifted.coords.re, reference.coords.re, 2e-5f);
    HS_EXPECT_NEAR(shifted.coords.im, reference.coords.im, 2e-5f);
    HS_EXPECT_EQ(shifted.region_id, reference.region_id);
  }
  const auto airocean_shifted =
      shaderball::airocean_projection(on_meridian, MERIDIAN, false);
  const auto airocean_reference =
      shaderball::airocean_projection(zero_meridian, 0.0f, false);
  HS_EXPECT_NEAR(airocean_shifted.coords.re, airocean_reference.coords.re,
                 2e-5f);
  HS_EXPECT_NEAR(airocean_shifted.coords.im, airocean_reference.coords.im,
                 2e-5f);
  HS_EXPECT_EQ(airocean_shifted.region_id, airocean_reference.region_id);

  const float oracle_longitude = 23.0f * PI_F / 180.0f;
  const float oracle_latitude = 28.0f * PI_F / 180.0f;
  const Vector oracle_point(cosf(oracle_latitude) * cosf(oracle_longitude),
                            sinf(oracle_latitude),
                            cosf(oracle_latitude) * sinf(oracle_longitude));
  const auto airocean_oracle =
      shaderball::airocean_projection(oracle_point, 0.0f, false);
  HS_EXPECT_NEAR(airocean_oracle.coords.re, 2.1265288136f, 4e-5f);
  HS_EXPECT_NEAR(airocean_oracle.coords.im, 3.6817439808f, 4e-5f);
  const auto airocean_horizontal =
      shaderball::airocean_projection(oracle_point, 0.0f, true);
  HS_EXPECT_NEAR(airocean_horizontal.coords.re,
                 5.7830422333f - airocean_oracle.coords.im, 4e-5f);
  HS_EXPECT_NEAR(airocean_horizontal.coords.im, airocean_oracle.coords.re,
                 4e-5f);

  struct AiroceanOracle {
    Vector point;
    uint8_t face;
    Complex coords;
  };
  const AiroceanOracle face_oracles[] = {
      {Vector(-0.0913057694774583f, -0.7547095802227720f, 0.6496743076188996f),
       18, Complex(0.751563669073068f, 5.176900182371453f)},
      {Vector(-0.125688534945440f, -0.951056516295154f, 0.282301071545602f), 19,
       Complex(1.209325099973815f, 0.174221877187437f)},
      {Vector(-0.788802981658962f, -0.156434465040231f, 0.594405681562271f), 20,
       Complex(0.504542659515686f, 4.382713723716426f)},
      {Vector(-0.540579378237808f, 0.121869343405147f, 0.832419244709073f), 21,
       Complex(0.821525205008616f, 4.204640141604954f)},
      {Vector(-0.806181892476771f, 0.275637355816999f, 0.523540642472755f), 22,
       Complex(0.414678953083584f, 3.507358320698008f)},
  };
  for (const auto &oracle : face_oracles) {
    const auto mapped =
        shaderball::airocean_projection(oracle.point, 0.0f, false);
    HS_EXPECT_EQ(mapped.region_id, oracle.face);
    HS_EXPECT_NEAR(mapped.coords.re, oracle.coords.re, 2e-5f);
    HS_EXPECT_NEAR(mapped.coords.im, oracle.coords.im, 2e-5f);
  }

  const Vector glued_points[] = {
      {0.00321530370315651f, 0.902108328962722f, 0.431497653108545f},
      {0.00321964224570043f, 0.902105871397194f, 0.431502758617508f},
      {0.00321096516044884f, 0.902110786482306f, 0.431492547577607f},
  };
  const uint8_t glued_faces[] = {1, 1, 2};
  const uint8_t glued_edges[] = {0, 0, 2};
  const Complex glued_coords[] = {
      {1.365889495965044f, 3.417252228774369f},
      {1.365892745153149f, 3.417257856533251f},
      {1.365886246776939f, 3.417246601015487f},
  };
  const float glued_cut_distances[] = {65536.0f, 65536.0f, 65536.0f};
  for (size_t index = 0; index < std::size(glued_points); ++index) {
    const auto mapped =
        shaderball::airocean_projection(glued_points[index], 0.0f, false);
    HS_EXPECT_EQ(mapped.region_id, glued_faces[index]);
    HS_EXPECT_EQ(mapped.edge_class,
                 shaderball::airocean_edge_identity(glued_faces[index],
                                                    glued_edges[index]));
    HS_EXPECT_TRUE((mapped.traits & shaderball::projection_traits(
                                        shaderball::ProjectionTrait::GLUED)) !=
                   0);
    HS_EXPECT_NEAR(mapped.coords.re, glued_coords[index].re, 5e-5f);
    HS_EXPECT_NEAR(mapped.coords.im, glued_coords[index].im, 5e-5f);
    HS_EXPECT_NEAR(mapped.fade_edge_distance, glued_cut_distances[index],
                   5e-5f);
  }

  const Vector cut_points[] = {
      {0.456082461583071f, 0.767833736470940f, -0.449911259442795f},
      {0.456076125629434f, 0.767836786220573f, -0.449912477441232f},
      {0.456088797513479f, 0.767830686682203f, -0.449910041421443f},
  };
  const uint8_t cut_faces[] = {3, 3, 4};
  const uint8_t cut_edges[] = {0, 0, 2};
  const Complex cut_coords[] = {
      {1.821185994620058f, 2.628655560595667f},
      {1.821179496243847f, 2.628655560595667f},
      {2.276485742463179f, 2.891526744414116f},
  };
  const float cut_distances[] = {0.0f, 6.498376211e-6f, 6.498376211e-6f};
  for (size_t index = 0; index < std::size(cut_points); ++index) {
    const auto mapped =
        shaderball::airocean_projection(cut_points[index], 0.0f, false);
    HS_EXPECT_EQ(mapped.region_id, cut_faces[index]);
    HS_EXPECT_EQ(mapped.edge_class, shaderball::airocean_edge_identity(
                                        cut_faces[index], cut_edges[index]));
    HS_EXPECT_TRUE(
        (mapped.traits &
         shaderball::projection_traits(shaderball::ProjectionTrait::CUT)) != 0);
    HS_EXPECT_NEAR(mapped.coords.re, cut_coords[index].re, 5e-5f);
    HS_EXPECT_NEAR(mapped.coords.im, cut_coords[index].im, 5e-5f);
    HS_EXPECT_NEAR(mapped.fade_edge_distance, cut_distances[index], 5e-5f);
  }
  for (uint8_t face = 0; face < 23; ++face) {
    const auto &face_center = shaderball::AIROCEAN_CENTERS[face];
    const auto mapped = shaderball::airocean_projection(
        Vector(face_center.x, face_center.z, face_center.y), 0.0f, false);
    HS_EXPECT_EQ(mapped.region_id, face);
  }
  HS_EXPECT_TRUE(shaderball::airocean_edge_is_cut(14, 0));
  HS_EXPECT_EQ(shaderball::airocean_edge_identity(14, 0), uint8_t(42));
  HS_EXPECT_EQ(shaderball::airocean_edge_identity(18, 0), uint8_t(54));

  const auto japan_edge_point = [](float from_weight) {
    const auto &a = shaderball::AIROCEAN_FACES[14][0];
    const auto &b = shaderball::AIROCEAN_FACES[14][1];
    const auto &center = shaderball::AIROCEAN_CENTERS[14];
    const float to_weight = 1.0f - from_weight;
    return Vector(from_weight * a.x + to_weight * b.x + 1e-4f * center.x,
                  from_weight * a.z + to_weight * b.z + 1e-4f * center.z,
                  from_weight * a.y + to_weight * b.y + 1e-4f * center.y)
        .normalized();
  };
  const auto japan_cut =
      shaderball::airocean_projection(japan_edge_point(0.75f), 0.0f, false);
  const auto japan_glued =
      shaderball::airocean_projection(japan_edge_point(0.25f), 0.0f, false);
  HS_EXPECT_EQ(japan_cut.region_id, uint8_t(14));
  HS_EXPECT_EQ(japan_cut.edge_class, uint8_t(42));
  HS_EXPECT_TRUE(
      (japan_cut.traits &
       shaderball::projection_traits(shaderball::ProjectionTrait::CUT)) != 0);
  HS_EXPECT_LT(japan_cut.fade_edge_distance, 1e-3f);
  HS_EXPECT_EQ(japan_glued.region_id, uint8_t(14));
  HS_EXPECT_EQ(japan_glued.edge_class, uint8_t(54));
  HS_EXPECT_TRUE(
      (japan_glued.traits &
       shaderball::projection_traits(shaderball::ProjectionTrait::GLUED)) != 0);
  HS_EXPECT_GT(japan_glued.fade_edge_distance, 0.1f);
  auto same_point = [](const shaderball::AiroceanPoint &a,
                       const shaderball::AiroceanPoint &b) {
    return fabsf(a.x - b.x) <= 1e-6f && fabsf(a.y - b.y) <= 1e-6f;
  };
  for (uint8_t face = 0; face < 23; ++face) {
    for (uint8_t edge_index = 0; edge_index < 3; ++edge_index) {
      bool expected_cut = false;
      for (size_t index = 0; index < std::size(shaderball::AIROCEAN_CUT_FACES);
           ++index)
        expected_cut |= shaderball::AIROCEAN_CUT_FACES[index] == face &&
                        shaderball::AIROCEAN_CUT_EDGES[index] == edge_index;
      HS_EXPECT_EQ(shaderball::airocean_edge_is_cut(face, edge_index),
                   expected_cut);

      const auto &a = shaderball::AIROCEAN_PLANAR_FACES[face][edge_index];
      const auto &b =
          shaderball::AIROCEAN_PLANAR_FACES[face][(edge_index + 1) % 3];
      uint8_t expected_identity = face * 3 + edge_index;
      bool found = false;
      for (uint8_t candidate_face = 0; candidate_face < 23 && !found;
           ++candidate_face) {
        for (uint8_t candidate_edge = 0; candidate_edge < 3; ++candidate_edge) {
          const auto &c =
              shaderball::AIROCEAN_PLANAR_FACES[candidate_face][candidate_edge];
          const auto &d =
              shaderball::AIROCEAN_PLANAR_FACES[candidate_face]
                                               [(candidate_edge + 1) % 3];
          if ((same_point(a, c) && same_point(b, d)) ||
              (same_point(a, d) && same_point(b, c))) {
            expected_identity = candidate_face * 3 + candidate_edge;
            found = true;
            break;
          }
        }
      }
      HS_EXPECT_EQ(shaderball::airocean_edge_identity(face, edge_index),
                   expected_identity);
    }
  }

  const auto peirce_seam_a = shaderball::peirce_projection(
      lon_lat(0.25f * PI_F - TIE_EPS, SOUTH_LATITUDE), 0.0f, 1, 0.0f);
  const auto peirce_seam_b = shaderball::peirce_projection(
      lon_lat(0.25f * PI_F + TIE_EPS, SOUTH_LATITUDE), 0.0f, 1, 0.0f);
  const Complex peirce_seam_delta = peirce_seam_a.coords - peirce_seam_b.coords;
  const Complex airocean_seam_delta = cut_coords[1] - cut_coords[2];
  HS_EXPECT_GT(sqrtf(peirce_seam_delta.re * peirce_seam_delta.re +
                     peirce_seam_delta.im * peirce_seam_delta.im),
               2.0f);
  HS_EXPECT_LT(peirce_seam_a.fade_edge_distance, 2e-5f);
  HS_EXPECT_LT(peirce_seam_b.fade_edge_distance, 2e-5f);
  HS_EXPECT_GT(sqrtf(airocean_seam_delta.re * airocean_seam_delta.re +
                     airocean_seam_delta.im * airocean_seam_delta.im),
               0.5f);

  for (int latitude_index = -8; latitude_index <= 8; ++latitude_index) {
    const float latitude = latitude_index * (PI_F / 18.0f);
    for (int longitude_index = -18; longitude_index < 18; ++longitude_index) {
      const float longitude = longitude_index * (PI_F / 18.0f);
      const float cp = cosf(latitude);
      const Vector v(cp * cosf(longitude), sinf(latitude),
                     cp * sinf(longitude));
      const auto b = shaderball::bonne_projection(v, 0.37f, standard_parallel);
      HS_EXPECT_TRUE(std::isfinite(b.coords.re));
      HS_EXPECT_TRUE(std::isfinite(b.coords.im));
      HS_EXPECT_TRUE(std::isfinite(b.fade_edge_distance));
      for (uint8_t layout = 0; layout < 4; ++layout) {
        const auto p = shaderball::peirce_projection(v, -0.21f, layout, 0.13f);
        const auto without_edge =
            shaderball::peirce_projection(v, -0.21f, layout, 0.13f, false);
        HS_EXPECT_TRUE(std::isfinite(p.coords.re));
        HS_EXPECT_TRUE(std::isfinite(p.coords.im));
        HS_EXPECT_TRUE(std::isfinite(p.fade_edge_distance));
        HS_EXPECT_EQ(without_edge.coords.re, p.coords.re);
        HS_EXPECT_EQ(without_edge.coords.im, p.coords.im);
        HS_EXPECT_EQ(without_edge.region_id, p.region_id);
        HS_EXPECT_EQ(without_edge.edge_class, p.edge_class);
        HS_EXPECT_EQ(without_edge.fade_edge_distance, 65536.0f);
      }
      for (bool horizontal : {false, true}) {
        const auto a = shaderball::airocean_projection(v, 0.19f, horizontal);
        const auto without_edge =
            shaderball::airocean_projection(v, 0.19f, horizontal, false);
        HS_EXPECT_TRUE(std::isfinite(a.coords.re));
        HS_EXPECT_TRUE(std::isfinite(a.coords.im));
        HS_EXPECT_TRUE(std::isfinite(a.fade_edge_distance));
        HS_EXPECT_TRUE(a.region_id < 23);
        HS_EXPECT_EQ(without_edge.coords.re, a.coords.re);
        HS_EXPECT_EQ(without_edge.coords.im, a.coords.im);
        HS_EXPECT_EQ(without_edge.region_id, a.region_id);
        HS_EXPECT_EQ(without_edge.edge_class, a.edge_class);
        HS_EXPECT_EQ(without_edge.traits, a.traits);
        HS_EXPECT_EQ(without_edge.fade_edge_distance, 65536.0f);
      }
    }
  }
}

/** @brief Domain policies, gauges, and analytic admission reject unsafe tuples. */
inline void test_shaderball_projection_and_admission_contracts() {
  using WB = ShaderBallWhiteBox;
  const uint8_t periodic_traits = shaderball::projection_traits(
      shaderball::ProjectionTrait::GLUED, shaderball::ProjectionTrait::FOLDED,
      shaderball::ProjectionTrait::PERIODIC);
  const WB::ProjectedLookup horizontal_left{
      Complex(-3.7080f, 0.2f), 2, 0,   2, 0.1f, 1.0f, 1,
      periodic_traits,         4, 1.0f};
  const WB::ProjectedLookup horizontal_right{
      Complex(3.7080f, 0.2f), 2, 0, 2, 0.1f, 1.0f, 1, periodic_traits, 4, 1.0f};
  HS_EXPECT_FALSE(WB::join_compatible(horizontal_left, horizontal_right,
                                      WB::Projection::PEIRCE_QUINCUNCIAL));
  WB::ProjectedLookup horizontal_neighbor = horizontal_left;
  horizontal_neighbor.coords.re = -3.6f;
  HS_EXPECT_FALSE(WB::join_compatible(horizontal_left, horizontal_neighbor,
                                      WB::Projection::PEIRCE_QUINCUNCIAL));
  WB::ProjectedLookup vertical_bottom = horizontal_left;
  WB::ProjectedLookup vertical_top = horizontal_right;
  vertical_bottom.edge_class = 5;
  vertical_top.edge_class = 5;
  vertical_bottom.coords = Complex(0.2f, -3.7080f);
  vertical_top.coords = Complex(0.2f, 3.7080f);
  HS_EXPECT_FALSE(WB::join_compatible(vertical_bottom, vertical_top,
                                      WB::Projection::PEIRCE_QUINCUNCIAL));

  auto kernel_lookup = [](const shaderball::ProjectionKernelResult &result) {
    return WB::ProjectedLookup{result.coords,
                               result.region_id,
                               result.component_id,
                               result.boundary_flags,
                               result.fade_edge_distance,
                               1.0f,
                               result.flags,
                               result.traits,
                               result.edge_class,
                               1.0f};
  };
  auto lon_lat = [](float longitude, float latitude) {
    const float cp = cosf(latitude);
    return Vector(cp * cosf(longitude), sinf(latitude), cp * sinf(longitude));
  };
  const Vector sector_before = lon_lat(0.25f * PI_F - 1e-4f, -0.4f);
  const Vector sector_after = lon_lat(0.25f * PI_F + 1e-4f, -0.4f);
  const auto diamond_before = kernel_lookup(
      shaderball::peirce_projection(sector_before, 0.0f, 0, 0.0f));
  const auto diamond_after =
      kernel_lookup(shaderball::peirce_projection(sector_after, 0.0f, 0, 0.0f));
  HS_EXPECT_FALSE(WB::join_compatible(diamond_before, diamond_after,
                                      WB::Projection::PEIRCE_QUINCUNCIAL));
  const auto horizontal_before = kernel_lookup(
      shaderball::peirce_projection(sector_before, 0.0f, 2, 0.0f));
  const auto horizontal_after =
      kernel_lookup(shaderball::peirce_projection(sector_after, 0.0f, 2, 0.0f));
  HS_EXPECT_FALSE(WB::join_compatible(horizontal_before, horizontal_after,
                                      WB::Projection::PEIRCE_QUINCUNCIAL));

  const auto glued_side_a = kernel_lookup(shaderball::airocean_projection(
      Vector(0.00321964224570043f, 0.902105871397194f, 0.431502758617508f),
      0.0f, false));
  const auto glued_side_b = kernel_lookup(shaderball::airocean_projection(
      Vector(0.00321096516044884f, 0.902110786482306f, 0.431492547577607f),
      0.0f, false));
  HS_EXPECT_FALSE(WB::join_compatible(glued_side_a, glued_side_b,
                                      WB::Projection::AIROCEAN));
  const auto cut_side_a = kernel_lookup(shaderball::airocean_projection(
      Vector(0.456076125629434f, 0.767836786220573f, -0.449912477441232f), 0.0f,
      false));
  const auto cut_side_b = kernel_lookup(shaderball::airocean_projection(
      Vector(0.456088797513479f, 0.767830686682203f, -0.449910041421443f), 0.0f,
      false));
  HS_EXPECT_FALSE(
      WB::join_compatible(cut_side_a, cut_side_b, WB::Projection::AIROCEAN));

  const Vector front_neighbor(1.0f, 1e-5f, 0.0f);
  const Vector back_neighbor(1.0f, -1e-5f, 0.0f);
  const Vector axis(1.0f, 0.0f, 0.0f);
  const auto front = WB::finalize_projection(
      front_neighbor,
      WB::project_point(front_neighbor, WB::Projection::GNOMONIC),
      WB::Projection::GNOMONIC, 1.0f,
      WB::GnomonicHemispherePolicy::FRONT_HEMISPHERE);
  const auto back_clipped = WB::finalize_projection(
      back_neighbor, WB::project_point(back_neighbor, WB::Projection::GNOMONIC),
      WB::Projection::GNOMONIC, 1.0f,
      WB::GnomonicHemispherePolicy::FRONT_HEMISPHERE);
  const auto axis_front = WB::finalize_projection(
      axis, WB::project_point(axis, WB::Projection::GNOMONIC),
      WB::Projection::GNOMONIC, 1.0f,
      WB::GnomonicHemispherePolicy::FRONT_HEMISPHERE);
  const auto axis_back = WB::finalize_projection(
      axis, WB::project_point(axis, WB::Projection::GNOMONIC),
      WB::Projection::GNOMONIC, 1.0f,
      WB::GnomonicHemispherePolicy::BACK_HEMISPHERE);
  HS_EXPECT_EQ(front.domain_coverage, 1.0f);
  HS_EXPECT_EQ(back_clipped.domain_coverage, 0.0f);
  HS_EXPECT_EQ(axis_front.domain_coverage, 1.0f);
  HS_EXPECT_EQ(axis_back.domain_coverage, 0.0f);

  WB::RequestedConfig bonne = WB::legacy_config();
  bonne.slots.projection = WB::Projection::BONNE;
  bonne.params.projection.bonne_standard_parallel = 0.0f;
  HS_EXPECT_FALSE(WB::valid_config(bonne));
  bonne.params.projection.bonne_standard_parallel = 1e-3f;
  HS_EXPECT_TRUE(WB::valid_config(bonne));

  MobiusParams mobius(1.0f, 0.2f, 0.3f, -0.1f, -0.2f, 0.1f, 0.9f, -0.15f);
  WB::canonicalize_mobius(mobius);
  const Complex coefficients[] = {mobius.a, mobius.b, mobius.c, mobius.d};
  float norm_sq = 0.0f;
  size_t pivot = 0;
  float pivot_magnitude = 0.0f;
  for (size_t index = 0; index < std::size(coefficients); ++index) {
    const float magnitude = coefficients[index].re * coefficients[index].re +
                            coefficients[index].im * coefficients[index].im;
    norm_sq += magnitude;
    if (magnitude > pivot_magnitude) {
      pivot_magnitude = magnitude;
      pivot = index;
    }
  }
  HS_EXPECT_NEAR(norm_sq, 1.0f, 2e-5f);
  HS_EXPECT_NEAR(coefficients[pivot].im, 0.0f, 2e-5f);
  HS_EXPECT_GT(coefficients[pivot].re, 0.0f);
  WB::RequestedConfig mobius_config = WB::legacy_config();
  mobius_config.slots.surface_lens = WB::SurfaceLens::MOBIUS;
  mobius_config.params.surface_lens.mobius = mobius;
  HS_EXPECT_TRUE(WB::valid_config(mobius_config));

  WB::RequestedConfig curl = WB::legacy_config();
  curl.slots.warp_program.outer.kind = WB::WarpStageKind::CURL_FLOW;
  curl.slots.warp_program.outer.curl_integrator =
      WB::CurlIntegrator::MIDPOINT_4;
  curl.params.warp.outer.scale = 1.0f;
  curl.params.warp.outer.time_scale = 0.0f;
  curl.params.warp.outer.strength = 0.03f;
  HS_EXPECT_TRUE(WB::valid_config(curl));
  curl.params.warp.outer.strength = 0.04f;
  HS_EXPECT_FALSE(WB::valid_config(curl));

  WB::RequestedConfig affine = WB::legacy_config();
  affine.slots.warp_program.outer.kind = WB::WarpStageKind::AFFINE_FRAME;
  affine.params.warp.outer.scale_x = 0.25f;
  affine.params.warp.outer.scale_y = 0.25f;
  affine.params.warp.outer.shear = 0.75f;
  HS_EXPECT_FALSE(WB::valid_config(affine));

  FastNoiseLite noise;
  const Complex curl_one =
      WB::curl_vector(Complex(), noise, WB::NoiseBasis::SIMPLEX, 1.0f, 0.2f);
  const Complex curl_two =
      WB::curl_vector(Complex(), noise, WB::NoiseBasis::SIMPLEX, 2.0f, 0.2f);
  HS_EXPECT_NEAR(curl_two.re, 2.0f * curl_one.re, 1e-5f);
  HS_EXPECT_NEAR(curl_two.im, 2.0f * curl_one.im, 1e-5f);
  for (WB::NoiseBasis basis : {WB::NoiseBasis::SIMPLEX, WB::NoiseBasis::FBM3,
                               WB::NoiseBasis::RIDGED3}) {
    HS_EXPECT_NEAR(WB::wrapped_noise(noise, basis, 0.37f, -0.91f, 0.0f),
                   WB::wrapped_noise(noise, basis, 0.37f, -0.91f, 1.0f), 1e-6f);
    HS_EXPECT_NEAR(WB::wrapped_noise(noise, basis, 0.37f, -0.91f, 1.0f - 1e-6f),
                   WB::wrapped_noise(noise, basis, 0.37f, -0.91f, 1e-6f),
                   2e-3f);
  }

  reset_effect_globals();
  WB::SB sb;
  sb.init();
  WB::FrameState frame = WB::frame(sb);
  frame.resources.outer_warp_noise = nullptr;
  const Complex point(0.31f, -0.27f);
  const WB::ProjectedLookup projected(point, 0, 0, 0, 1.0f, 1.0f, 0);
  WB::WarpStageParams periodic_params;
  periodic_params.strength = 0.8f;
  periodic_params.frequency = 2.5f;
  periodic_params.turns = 0.75f;
  periodic_params.center_orbit_radius = 0.4f;
  for (WB::WarpStageKind kind :
       {WB::WarpStageKind::WAVE_SHEAR, WB::WarpStageKind::VORTEX}) {
    WB::WarpStageSpec spec{kind};
    frame.clocks.warp_outer_phase = 0.0f;
    const auto at_zero =
        WB::warp_stage(point, projected, spec, periodic_params, frame);
    frame.clocks.warp_outer_phase = 1.0f;
    const auto at_wrap =
        WB::warp_stage(point, projected, spec, periodic_params, frame);
    HS_EXPECT_NEAR(at_zero.coords.re, at_wrap.coords.re, 2e-5f);
    HS_EXPECT_NEAR(at_zero.coords.im, at_wrap.coords.im, 2e-5f);
  }
}

/** @brief Every shader catalog family produces finite bounded output. */
inline void test_shaderball_kernel_catalog() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();
  const Vector view = Vector(0.31f, 0.87f, -0.38f).normalized();
  WB::RequestedConfig config = WB::legacy_config();
  config.slots.surface_lens = WB::SurfaceLens::NONE;
  config.slots.warp_program.outer.kind = WB::WarpStageKind::NONE;
  config.slots.warp_program.inner.kind = WB::WarpStageKind::NONE;
  config.slots.coverage = WB::CoveragePolicy::OPAQUE;
  config.params.source.pattern_freq = 1.0f;

  auto check = [&](const WB::RequestedConfig &candidate) {
    HS_EXPECT_TRUE(WB::valid_config(candidate));
    WB::request_config(sb, candidate);
    const WB::RequestedConfig canonical = WB::requested_config(sb);
    HS_EXPECT_TRUE(canonical.slots == candidate.slots);
    WB::settle_transition(sb);
    HS_EXPECT_TRUE(WB::active_config(sb) == canonical);
    const Color4 color = WB::shade(view, WB::frame(sb));
    HS_EXPECT_TRUE(std::isfinite(color.alpha));
    HS_EXPECT_GE(color.alpha, 0.0f);
    HS_EXPECT_LE(color.alpha, 1.0f);
  };

  for (uint8_t value = 0; value <= 6; ++value) {
    config.slots.function = static_cast<WB::Function>(value);
    config.params.source.noise_basis = WB::NoiseBasis::FBM3;
    check(config);
  }
  config.slots.function = WB::Function::TWIN_WAVE;
  for (uint8_t value = 0; value <= 8; ++value) {
    config.slots.warp_program.outer.kind =
        static_cast<WB::WarpStageKind>(value);
    config.slots.projection = value == 1 ? WB::Projection::STEREOGRAPHIC
                                         : WB::Projection::EQUIRECTANGULAR;
    config.params.warp.outer.strength = value == 0   ? 0.0f
                                        : value == 6 ? 0.005f
                                                     : 0.35f;
    config.params.warp.outer.turns = 0.4f;
    config.params.warp.outer.translation_x = 0.2f;
    config.params.warp.outer.translation_y = -0.1f;
    if (value >= 3 && value <= 6)
      config.params.warp.outer.time_scale = 0.0f;
    config.slots.warp_program.outer.basis = WB::NoiseBasis::SIMPLEX;
    if (value == 8)
      config.slots.function = WB::Function::COUPLED_DIRECT;
    check(config);
  }
  config.slots.warp_program.outer.kind = WB::WarpStageKind::NONE;
  config.params.warp.outer.strength = 0.0f;
  config.slots.projection = WB::Projection::EQUIRECTANGULAR;
  for (uint8_t value = 0; value <= 5; ++value) {
    config.slots.surface_lens = static_cast<WB::SurfaceLens>(value);
    config.params.surface_lens.mix = value == 0 ? 0.0f : 0.6f;
    check(config);
  }
  config.slots.surface_lens = WB::SurfaceLens::NONE;
  config.params.surface_lens.mix = 0.0f;
  for (uint8_t transfer = 0; transfer <= 3; ++transfer) {
    config.slots.value_transfer = static_cast<WB::ValueTransfer>(transfer);
    for (uint8_t coverage = 0; coverage <= 4; ++coverage) {
      config.slots.coverage = static_cast<WB::CoveragePolicy>(coverage);
      for (uint8_t colorizer = 0; colorizer <= 2; ++colorizer) {
        config.slots.colorizer = static_cast<WB::Colorizer>(colorizer);
        check(config);
      }
    }
  }

  WB::FrameState frame = WB::frame(sb);
  frame.resources.outer_warp_noise = nullptr;
  WB::WarpStageParams zero_params;
  zero_params.strength = 0.0f;
  const Complex input(0.27f, -0.41f);
  for (uint8_t value = 0; value <= 6; ++value) {
    WB::WarpStageSpec spec{static_cast<WB::WarpStageKind>(value)};
    const auto identity = WB::warp_stage(input, {input, 0, 0, 0, 1.0f, 1.0f, 0},
                                         spec, zero_params, frame);
    HS_EXPECT_EQ(identity.coords.re, input.re);
    HS_EXPECT_EQ(identity.coords.im, input.im);
    HS_EXPECT_EQ(identity.deformation, 0.0f);
    HS_EXPECT_EQ(identity.path_length, 0.0f);
  }
  for (uint8_t value = 7; value <= 8; ++value) {
    WB::WarpStageSpec spec{static_cast<WB::WarpStageKind>(value)};
    const auto mapped = WB::warp_stage(input, {input, 0, 0, 0, 1.0f, 1.0f, 0},
                                       spec, zero_params, frame);
    HS_EXPECT_TRUE(std::isfinite(mapped.coords.re));
    HS_EXPECT_TRUE(std::isfinite(mapped.coords.im));
  }
}

/** @brief Adjacent preset edges use one live topology and authoritative clocks. */
inline void test_shaderball_stable_preset_transition() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();
  const WB::Params from = WB::presets()[0].params;
  const WB::Params to = WB::presets()[1].params;
  const size_t initial_index = WB::preset_index(sb);
  WB::begin_blend(sb);
  HS_EXPECT_EQ(WB::preset_index(sb), initial_index + 1);
  HS_EXPECT_TRUE(WB::param_morph_active(sb));
  HS_EXPECT_FALSE(WB::transition_active(sb));
  HS_EXPECT_TRUE(WB::live_params(sb) == from);

  float previous_phase = WB::clocks(sb).source_primary;
  for (int step = 0; step <= 60; ++step) {
    WB::step_param_morph(sb);
    const float live_speed = WB::live_params(sb).source.speed;
    const float phase = WB::clocks(sb).source_primary;
    HS_EXPECT_NEAR(phase, fmodf(previous_phase + live_speed, TWO_PI_F), 1e-6f);
    previous_phase = phase;
    if (WB::param_morph_elapsed(sb) == 6) {
      HS_EXPECT_LT(WB::live_params(sb).source.pattern_freq,
                   from.source.pattern_freq);
      HS_EXPECT_GT(WB::live_params(sb).source.pattern_freq,
                   to.source.pattern_freq);
      HS_EXPECT_EQ(WB::live_params(sb).warp, from.warp);
    }
    if (WB::param_morph_elapsed(sb) == 31) {
      HS_EXPECT_LT(WB::live_params(sb).source.pattern_freq,
                   from.source.pattern_freq);
      HS_EXPECT_GT(WB::live_params(sb).source.pattern_freq,
                   to.source.pattern_freq);
    }
  }
  HS_EXPECT_FALSE(WB::param_morph_active(sb));
  HS_EXPECT_TRUE(WB::live_params(sb) == to);

  const auto &presets = WB::presets();
  for (size_t index = 0; index < presets.size(); ++index) {
    HS_EXPECT_TRUE(WB::valid_config(presets[index]));
    const auto &next = presets[(index + 1) % presets.size()];
    if (index < 14)
      HS_EXPECT_TRUE(WB::slots_equal(presets[index].slots, next.slots));
    HS_EXPECT_TRUE(WB::transition_admitted(presets[index], next));
  }
}

/** @brief Discrete look changes preserve continuous dual-runtime handoff. */
inline void test_shaderball_discrete_transition() {
  using WB = ShaderBallWhiteBox;

  const Color4 from_color(Pixel(10000, 20000, 30000), 0.25f);
  const Color4 to_color(Pixel(40000, 50000, 60000), 0.75f);
  HS_EXPECT_TRUE(WB::blend_outputs(from_color, to_color, 0.0f).color ==
                 from_color.color);
  HS_EXPECT_TRUE(WB::blend_outputs(from_color, to_color, 1.0f).color ==
                 to_color.color);
  const Color4 middle = WB::blend_outputs(from_color, to_color, 0.5f);
  HS_EXPECT_EQ(middle.alpha, 0.5f);
  HS_EXPECT_EQ(middle.color.r, uint16_t(32500));
  HS_EXPECT_EQ(middle.color.g, uint16_t(42500));
  HS_EXPECT_EQ(middle.color.b, uint16_t(52500));
  const Color4 transparent = WB::blend_outputs(
      Color4(Pixel(65535, 0, 0), 0.0f), Color4(Pixel(0, 0, 65535), 0.0f), 0.5f);
  HS_EXPECT_EQ(transparent.alpha, 0.0f);
  HS_EXPECT_TRUE(transparent.color == Pixel());

  {
    reset_effect_globals();
    WB::SB sb;
    sb.init();
    const Vector view(0.2f, 0.9f, -0.3f);
    const WB::FrameState valid = WB::frame(sb);
    const Color4 expected = WB::shade(view, valid);
    const auto clear_phase = WB::through_clear_phase(30, 60);
    const Color4 clear = WB::shade_through_clear(view, nullptr, clear_phase);
    HS_EXPECT_EQ(clear.alpha, 0.0f);
    HS_EXPECT_TRUE(clear.color == Pixel());
    const auto from_phase = WB::through_clear_phase(15, 60);
    const Color4 from_only = WB::shade_through_clear(view, &valid, from_phase);
    HS_EXPECT_TRUE(std::isfinite(from_only.alpha));
    HS_EXPECT_GT(from_only.alpha, 0.0f);
    const auto to_phase = WB::through_clear_phase(45, 60);
    const Color4 to_only = WB::shade_through_clear(view, &valid, to_phase);
    HS_EXPECT_TRUE(std::isfinite(to_only.alpha));
    HS_EXPECT_GT(to_only.alpha, 0.0f);
    const auto through_start = WB::through_clear_phase(0, 60);
    const Color4 exact_start =
        WB::shade_through_clear(view, &valid, through_start);
    HS_EXPECT_TRUE(exact_start.color == expected.color);
    HS_EXPECT_EQ(exact_start.alpha, expected.alpha);
    const auto through_end = WB::through_clear_phase(60, 60);
    const Color4 exact_end = WB::shade_through_clear(view, &valid, through_end);
    HS_EXPECT_TRUE(exact_end.color == expected.color);
    HS_EXPECT_EQ(exact_end.alpha, expected.alpha);
  }

  {
    reset_effect_globals();
    WB::SB sb;
    sb.init();
    sb.setAnimationsPaused(true);
    WB::RequestedConfig from = WB::legacy_config();
    from.slots.function = WB::Function::NOISE_CONTOUR;
    from.params.source.noise_basis = WB::NoiseBasis::SIMPLEX;
    WB::request_config(sb, from);
    WB::settle_transition(sb);
    HS_EXPECT_EQ(WB::prepared_noise_basis(sb, 2), WB::NoiseBasis::SIMPLEX);

    WB::RequestedConfig to = from;
    to.params.source.noise_basis = WB::NoiseBasis::FBM3;
    WB::force_transition(sb, to, 60, false);
    HS_EXPECT_TRUE(WB::transition_active(sb));
    for (int frame = 0; frame < 30; ++frame) {
      sb.draw_frame();
      sb.advance_display();
    }
    HS_EXPECT_EQ(WB::transition_elapsed(sb), uint16_t(30));
    HS_EXPECT_EQ(WB::prepared_noise_basis(sb, 2), WB::NoiseBasis::SIMPLEX);
    sb.draw_frame();
    sb.advance_display();
    HS_EXPECT_EQ(WB::transition_elapsed(sb), uint16_t(31));
    HS_EXPECT_EQ(WB::prepared_noise_basis(sb, 2), WB::NoiseBasis::FBM3);
  }

  {
    reset_effect_globals();
    WB::SB sb;
    sb.init();
    WB::LookRuntime liquid_runtime;
    WB::LookRuntime generated_runtime;
    const WB::WalkDeltas deltas{make_rotation(Y_AXIS, 0.2f),
                                make_rotation(X_AXIS, 0.3f)};
    WB::advance_runtime(sb, liquid_runtime, WB::presets()[0], deltas);
    WB::advance_runtime(sb, generated_runtime, WB::legacy_config(), deltas);
    HS_EXPECT_EQ(liquid_runtime.clocks.source_primary, 0.1f);
    HS_EXPECT_EQ(generated_runtime.clocks.source_primary, 0.05f);
    HS_EXPECT_TRUE(liquid_runtime.projection_wander !=
                   generated_runtime.projection_wander);
    HS_EXPECT_TRUE(liquid_runtime.outer_wander !=
                   generated_runtime.outer_wander);

    const size_t original_index = WB::preset_index(sb);
    WB::force_transition(sb, WB::legacy_config(), 60, true);
    const WB::RequestedConfig captured_source = WB::transition_from_config(sb);
    HS_EXPECT_TRUE(WB::transition_active(sb));
    HS_EXPECT_EQ(WB::transition_mix(sb), 0.0f);
    WB::begin_blend(sb);
    HS_EXPECT_EQ(WB::preset_index(sb), original_index);
    const uint32_t walk_steps = WB::walk_steps(sb);
    const uint32_t liquid_steps = WB::liquid_palette_steps(sb);
    const uint32_t generated_steps = WB::generated_palette_steps(sb);
    sb.draw_frame();
    sb.advance_display();
    HS_EXPECT_EQ(WB::walk_steps(sb), walk_steps + 1);
    HS_EXPECT_EQ(WB::liquid_palette_steps(sb), liquid_steps + 1);
    HS_EXPECT_EQ(WB::generated_palette_steps(sb), generated_steps + 1);
    HS_EXPECT_EQ(WB::transition_from_runtime(sb).clocks.source_primary, 0.1f);
    HS_EXPECT_EQ(WB::transition_to_runtime(sb).clocks.source_primary, 0.05f);

    for (int frame_index = 1; frame_index < 20; ++frame_index) {
      sb.draw_frame();
      sb.advance_display();
    }
    const uint16_t elapsed_before_takeover = WB::transition_elapsed(sb);
    sb.setAnimationsPaused(true);
    const float visible_phase =
        WB::transition_from_runtime(sb).clocks.source_primary;
    WB::RequestedConfig queued = WB::legacy_config();
    queued.slots.function = WB::Function::RINGS;
    HS_EXPECT_TRUE(WB::valid_config(queued));
    HS_EXPECT_TRUE(WB::transition_admitted(captured_source, queued));
    WB::request_config(sb, queued);
    HS_EXPECT_GT(elapsed_before_takeover, uint16_t(0));
    HS_EXPECT_FALSE(WB::transition_active(sb));
    HS_EXPECT_TRUE(WB::active_config(sb) == queued);
    HS_EXPECT_TRUE(WB::requested_slots(sb) == queued.slots);
    HS_EXPECT_EQ(WB::clocks(sb).source_primary, visible_phase);
    sb.draw_frame();
    sb.advance_display();
    HS_EXPECT_NEAR(WB::clocks(sb).source_primary,
                   fmodf(visible_phase + queued.params.source.speed, TWO_PI_F),
                   1e-6f);

    const WB::RequestedConfig manual = WB::presets()[6];
    HS_EXPECT_TRUE(WB::transition_admitted(captured_source, manual));
    WB::request_config(sb, manual);
    HS_EXPECT_FALSE(WB::transition_active(sb));
    HS_EXPECT_TRUE(WB::active_config(sb) == manual);
    const float committed_phase = WB::clocks(sb).source_primary;
    sb.draw_frame();
    sb.advance_display();
    HS_EXPECT_NEAR(
        WB::clocks(sb).source_primary,
        fmodf(committed_phase + manual.params.source.speed, TWO_PI_F), 1e-6f);
  }
}

/**
 * @brief A pause taken inside a through-clear transition holds a drawn frame.
 * @details The through-clear midpoint writes no pixels and a paused
 * choreography transition stops advancing, so the hold has to land on the
 * nearer endpoint or the sphere stays black for as long as the pause lasts.
 */
inline void test_shaderball_paused_through_clear_holds_endpoint() {
  using WB = ShaderBallWhiteBox;
  auto lit_pixels = [](const WB::SB &sb) {
    size_t count = 0;
    for (int y = 0; y < SMALL_H; ++y)
      for (int x = 0; x < SMALL_W; ++x) {
        const Pixel &pixel = sb.get_pixel(x, y);
        if ((pixel.r | pixel.g | pixel.b) != 0)
          ++count;
      }
    return count;
  };
  auto run_frames = [](WB::SB &sb, int frames) {
    for (int frame = 0; frame < frames; ++frame) {
      sb.draw_frame();
      sb.advance_display();
    }
  };

  {
    reset_effect_globals();
    WB::SB sb;
    sb.init();
    WB::force_transition(sb, WB::legacy_config(), 60, true);
    run_frames(sb, 30);
    HS_EXPECT_EQ(WB::transition_elapsed(sb), uint16_t(30));

    sb.setAnimationsPaused(true);
    run_frames(sb, 1);
    HS_EXPECT_EQ(lit_pixels(sb), size_t(0)); // the cleared midpoint itself
    run_frames(sb, 4);
    HS_EXPECT_TRUE(WB::transition_active(sb));
    HS_EXPECT_EQ(WB::transition_elapsed(sb), uint16_t(0));
    HS_EXPECT_GT(lit_pixels(sb), size_t(0));

    sb.setAnimationsPaused(false);
    run_frames(sb, 1);
    HS_EXPECT_EQ(WB::transition_elapsed(sb), uint16_t(1));
  }

  {
    reset_effect_globals();
    WB::SB sb;
    sb.init();
    const size_t entry_preset = WB::preset_index(sb);
    WB::force_transition(sb, WB::legacy_config(), 60, true);
    run_frames(sb, 40);
    HS_EXPECT_EQ(WB::transition_elapsed(sb), uint16_t(40));

    sb.setAnimationsPaused(true);
    run_frames(sb, 3);
    HS_EXPECT_TRUE(WB::transition_active(sb));
    HS_EXPECT_EQ(WB::transition_elapsed(sb), uint16_t(60));
    HS_EXPECT_GT(lit_pixels(sb), size_t(0));

    sb.setAnimationsPaused(false);
    run_frames(sb, 1);
    HS_EXPECT_FALSE(WB::transition_active(sb));
    HS_EXPECT_TRUE(WB::active_config(sb) == WB::legacy_config());
    HS_EXPECT_EQ(WB::preset_index(sb), entry_preset);
  }
}

/** @brief Generated and liquid color resources remain independent owners. */
inline void test_shaderball_palette_resources() {
  using WB = ShaderBallWhiteBox;
  uint32_t hue = 0;
  GenerativePalette previous;
  WB::make_triadic(hue, 0, previous);
  for (uint32_t sequence = 1; sequence <= 8; ++sequence) {
    GenerativePalette next;
    WB::make_triadic(hue, sequence, next);
    HS_EXPECT_TRUE(previous.morph_compatible(next));
    previous = next;
  }
  HS_EXPECT_EQ(hue, uint32_t(8) * WB::HUE_STEP);

  reset_effect_globals();
  WB::SB sb;
  sb.init();
  const Pixel generated = WB::generated_color(sb, 0.25f);
  const Pixel liquid = WB::liquid_color(sb, 0.25f);
  HS_EXPECT_TRUE(generated.r != liquid.r || generated.g != liquid.g ||
                 generated.b != liquid.b);
}

/** @brief Module entry point for ShaderBall contract tests. */
inline int run_shaderball_tests() {
  ModuleFixture fixture("shaderball");
  test_shaderball_clocks_wrapped();
  test_shaderball_pause_semantics();
  test_shaderball_paused_selector_commit();
  test_shaderball_manual_edit_timing();
  test_shaderball_pipeline_contract();
  test_shaderball_legacy_spatial_slots();
  test_shaderball_subduction_edge_fade();
  test_shaderball_legacy_sources();
  test_shaderball_coupled_source();
  test_shaderball_preset_bank();
  test_shaderball_config_admission();
  test_shaderball_deterministic_gui_edits();
  test_shaderball_work_admission();
  test_shaderball_strict_seam_admission();
  test_shaderball_additive_delta_precision();
  test_shaderball_profile_presets();
  test_shaderball_gui_catalog();
  test_shaderball_projection_catalog();
  test_shaderball_projection_and_admission_contracts();
  test_shaderball_kernel_catalog();
  test_shaderball_stable_preset_transition();
  test_shaderball_discrete_transition();
  test_shaderball_paused_through_clear_holds_endpoint();
  test_shaderball_palette_resources();
  return fixture.result();
}

} // namespace shaderball_tests
} // namespace hs_test
