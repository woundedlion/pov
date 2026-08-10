/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include "tests/test_effects.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace shadierball_tests {

using effects_tests::FRAME_MS;
using effects_tests::FRAME_US;
using effects_tests::reset_effect_globals;
using effects_tests::SMALL_H;
using effects_tests::SMALL_W;

/** @brief White-box access to ShadierBall's typed pipeline. */
struct ShadierBallWhiteBox {
  using SDB = ShadierBall<SMALL_W, SMALL_H>;
  using Function = SDB::Function;
  using Projection = SDB::Projection;
  using BonneHemisphere = SDB::BonneHemisphere;
  using GnomonicHemispherePolicy = SDB::GnomonicHemispherePolicy;
  using ProjectionFramePolicy = SDB::ProjectionFramePolicy;
  using SurfaceLens = SDB::SurfaceLens;
  using NoiseBasis = SDB::NoiseBasis;
  using WarpEnvelope = SDB::WarpEnvelope;
  using PolarMode = SDB::PolarMode;
  using CurlIntegrator = SDB::CurlIntegrator;
  using PolarHarmonic = SDB::PolarHarmonic;
  using WarpStageKind = SDB::WarpStageKind;
  using WarpStageSpec = SDB::WarpStageSpec;
  using WarpStageParams = SDB::WarpStageParams;
  using ProjectionParams = SDB::ProjectionParams;
  using SignalWeight = SDB::SignalWeight;
  using ValueTransfer = SDB::ValueTransfer;
  using CoveragePolicy = SDB::CoveragePolicy;
  using Colorizer = SDB::Colorizer;
  using Slots = SDB::Slots;
  using Params = SDB::Params;
  using RequestedConfig = SDB::RequestedConfig;
  using SourceState = SDB::SourceState;
  using FrameState = SDB::FrameState;
  using ProjectedLookup = SDB::ProjectedLookup;
  using PlanarWarpStageResult = SDB::PlanarWarpStageResult;
  using PlanarWarpResult = SDB::PlanarWarpResult;
  using MaterialSample = SDB::MaterialSample;
  using ClockState = SDB::ClockState;
  using LookRuntime = SDB::LookRuntime;
  using WalkDeltas = SDB::WalkDeltas;
  using DualOutputFrame = SDB::DualOutputFrame;
  using ThroughClearPhase = SDB::ThroughClearPhase;
  using TransitionMode = SDB::TransitionMode;
  using WorkSignature = SDB::WorkSignature;
  using WorkTier = SDB::WorkTier;

  static constexpr float AXIS_EPS = SDB::GNOMONIC_AXIS_EPS;
  static constexpr uint32_t HUE_STEP = SDB::HUE_STEP;

  static ClockState clocks(const SDB &sdb) { return sdb.runtime.clocks; }
  static void seed_clocks(SDB &sdb, float value) {
    sdb.runtime.clocks = {value, value, value, value, value, value};
  }
  static FrameState frame(const SDB &sdb) { return sdb.prepare_frame(); }
  static Slots active_slots(const SDB &sdb) { return sdb.active_slots; }
  static RequestedConfig active_config(const SDB &sdb) {
    return {sdb.active_slots, sdb.blend.params};
  }
  static constexpr Slots liquid_stereo_slots() {
    return SDB::LIQUID_STEREO_SLOTS;
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
  static void request_slots(SDB &sdb, const Slots &slots) {
    sdb.requested_config.slots = slots;
    sdb.requested_schema_bound = false;
    sdb.apply_requested_config();
  }
  static Slots requested_slots(const SDB &sdb) {
    return sdb.requested_config.slots;
  }
  static const RequestedConfig &requested_config(const SDB &sdb) {
    return sdb.requested_config;
  }
  static void request_config(SDB &sdb, const RequestedConfig &config) {
    sdb.requested_config = config;
    sdb.requested_schema_bound = false;
    sdb.apply_requested_config();
  }
  static bool slots_equal(const Slots &a, const Slots &b) { return a == b; }
  static constexpr bool valid_config(const RequestedConfig &config) {
    return SDB::valid_config(config);
  }
  static constexpr bool transition_admitted(const RequestedConfig &from,
                                            const RequestedConfig &to) {
    return SDB::transition_admitted(from, to);
  }
  static constexpr bool hold_admitted(const RequestedConfig &config) {
    return SDB::hold_admitted(config);
  }
  static constexpr bool stable_topology(const RequestedConfig &from,
                                        const RequestedConfig &to) {
    return SDB::stable_topology(from, to);
  }
  static constexpr bool
  stable_parameter_path_admitted(const RequestedConfig &from,
                                 const RequestedConfig &to) {
    return SDB::stable_parameter_path_admitted(from, to);
  }
  static constexpr uint16_t hold_work_unit_ceiling() {
    return SDB::HOLD_WORK_UNIT_CEILING;
  }
  static constexpr float noise_time_period() {
    return SDB::STEREO_NOISE_TIME_PERIOD;
  }
  static bool transition_active(const SDB &sdb) {
    return sdb.transition.active;
  }
  static bool param_morph_active(const SDB &sdb) {
    return sdb.param_morph.active;
  }
  static uint16_t param_morph_elapsed(const SDB &sdb) {
    return sdb.param_morph.elapsed;
  }
  static const Params &live_params(const SDB &sdb) { return sdb.blend.params; }
  static size_t preset_index(const SDB &sdb) { return sdb.preset_index; }
  static float transition_mix(const SDB &sdb) {
    return SDB::transition_mix(sdb.transition.elapsed, sdb.transition.duration);
  }
  static const LookRuntime &transition_from_runtime(const SDB &sdb) {
    return sdb.transition.from_runtime;
  }
  static const LookRuntime &transition_to_runtime(const SDB &sdb) {
    return sdb.transition.to_runtime;
  }
  static const RequestedConfig &transition_from_config(const SDB &sdb) {
    return sdb.transition.from_config;
  }
  static const RequestedConfig &transition_to_config(const SDB &sdb) {
    return sdb.transition.to_config;
  }
  static bool transition_continues_choreo(const SDB &sdb) {
    return sdb.transition.continue_choreo;
  }
  static uint16_t transition_elapsed(const SDB &sdb) {
    return sdb.transition.elapsed;
  }
  static NoiseBasis prepared_noise_basis(const SDB &sdb, uint8_t resource_id) {
    for (size_t index = 0; index < sdb.prepared_noise_count; ++index)
      if (sdb.prepared_noise_keys[index].resource_id == resource_id)
        return sdb.prepared_noise_keys[index].basis;
    return static_cast<NoiseBasis>(0xff);
  }
  static constexpr TransitionMode transition_mode(const RequestedConfig &from,
                                                  const RequestedConfig &to) {
    return SDB::transition_mode(from, to);
  }
  static constexpr WorkSignature work_signature(const RequestedConfig &config) {
    return SDB::work_signature(config);
  }
  static bool pending_transition_active(const SDB &sdb) {
    return sdb.pending_transition.active;
  }
  static const RequestedConfig &pending_transition_config(const SDB &sdb) {
    return sdb.requested_config;
  }
  static void force_transition(SDB &sdb, const RequestedConfig &to,
                               uint16_t duration, bool continue_choreo) {
    const RequestedConfig from = active_config(sdb);
    sdb.param_morph.active = false;
    sdb.transition = {from, to,       sdb.runtime,     sdb.runtime,
                      0,    duration, continue_choreo, true};
  }
  static const LookRuntime &runtime(const SDB &sdb) { return sdb.runtime; }
  static void advance_runtime(SDB &sdb, LookRuntime &runtime,
                              const RequestedConfig &config,
                              const WalkDeltas &deltas) {
    sdb.advance_runtime(runtime, config, deltas);
  }
  static Color4 shade_dual_output(const Vector &view,
                                  const DualOutputFrame &frame) {
    return SDB::shade_dual_output(view, frame);
  }
  static ThroughClearPhase through_clear_phase(uint16_t elapsed,
                                               uint16_t duration) {
    return SDB::through_clear_phase(elapsed, duration);
  }
  static Color4 shade_through_clear(const Vector &view,
                                    const FrameState *visible,
                                    const ThroughClearPhase &phase) {
    return SDB::shade_through_clear(view, visible, phase);
  }
  static void begin_blend(SDB &sdb) { sdb.begin_blend(); }
  static void step_param_morph(SDB &sdb) {
    sdb.prepare_param_morph();
    sdb.advance_runtime(sdb.runtime, {sdb.active_slots, sdb.blend.params},
                        {Quaternion(), Quaternion()});
    sdb.finish_transitions();
  }
  static void settle_transition(SDB &sdb) {
    for (int frame = 0;
         frame < 1024 && (sdb.transition.active || sdb.param_morph.active);
         ++frame) {
      sdb.draw_frame();
      sdb.advance_display();
    }
  }
  static uint32_t walk_steps(const SDB &sdb) { return sdb.walk_step_count; }
  static uint32_t liquid_palette_steps(const SDB &sdb) {
    return sdb.liquid_palette_step_count;
  }
  static uint32_t generated_palette_steps(const SDB &sdb) {
    return sdb.generated_palette_step_count;
  }
  static Color4 blend_outputs(const Color4 &from, const Color4 &to, float mix) {
    return SDB::blend_outputs(from, to, mix);
  }
  static ProjectedLookup join(const ProjectedLookup &direct,
                              const ProjectedLookup &lensed, float mix,
                              Projection projection, float pole_fade) {
    return SDB::join_projected(direct, lensed, mix, projection, pole_fade);
  }
  static bool join_compatible(const ProjectedLookup &direct,
                              const ProjectedLookup &lensed,
                              Projection projection,
                              float coordinate_scale = 1.0f) {
    return SDB::projection_join_compatible(direct, lensed, projection,
                                           coordinate_scale);
  }
  static constexpr uint8_t boundary_cut() { return SDB::BOUNDARY_CUT; }
  static constexpr uint8_t boundary_singular() {
    return SDB::BOUNDARY_SINGULAR;
  }
  static constexpr uint8_t projection_folded() {
    return SDB::PROJECTION_FLAG_FOLDED;
  }
  static Vector outer_lookup(const Vector &v, const FrameState &frame) {
    return SDB::outer_camera_lookup(v, frame);
  }
  static ProjectedLookup surface_project(const Vector &v,
                                         const FrameState &frame) {
    return SDB::surface_lens_project_lookup(v, frame);
  }
  static PlanarWarpResult warp(const ProjectedLookup &projected,
                               const FrameState &frame) {
    return SDB::planar_warp_lookup(projected, frame);
  }
  static PlanarWarpStageResult warp_stage(const Complex &input,
                                          const ProjectedLookup &projected,
                                          const WarpStageSpec &spec,
                                          const WarpStageParams &params,
                                          const FrameState &frame) {
    return SDB::warp_stage_lookup(input, projected, spec, params, frame);
  }
  static MaterialSample material(const ProjectedLookup &projected,
                                 const PlanarWarpResult &warped,
                                 const FrameState &frame) {
    const Complex source_coords =
        SDB::condition_source_coords(warped.coords, frame);
    const float field = SDB::sample_source(source_coords, frame);
    return SDB::shape_material(field, projected, warped, frame);
  }
  static MaterialSample shape(float field, const ProjectedLookup &projected,
                              const PlanarWarpResult &warped,
                              const FrameState &frame) {
    return SDB::shape_material(field, projected, warped, frame);
  }
  static Color4 colorize(const MaterialSample &sample,
                         const FrameState &frame) {
    return SDB::colorize(sample, frame);
  }
  static Color4 shade(const Vector &v, const FrameState &frame) {
    return SDB::shade(v, frame);
  }
  static Complex project_point(const Vector &v, Projection projection) {
    return SDB::project_point(v, projection);
  }
  static ProjectedLookup
  finalize_projection(const Vector &v, const Complex &coords,
                      Projection projection, float pole_fade,
                      GnomonicHemispherePolicy hemisphere) {
    return SDB::finalize_projection(v, coords, projection, pole_fade,
                                    hemisphere);
  }
  static void canonicalize_mobius(MobiusParams &params) {
    SDB::canonicalize_mobius(params);
  }
  static Complex curl_vector(const Complex &p, const FastNoiseLite &noise,
                             NoiseBasis basis, float scale, float time) {
    return SDB::curl_vector(p, noise, basis, scale, time);
  }
  static float wrapped_noise(const FastNoiseLite &noise, NoiseBasis basis,
                             float x, float y, float turns) {
    return SDB::sample_wrapped_noise_basis(noise, basis, x, y, turns);
  }
  static ProjectionParams lerp_projection(const ProjectionParams &a,
                                          const ProjectionParams &b, float t) {
    ProjectionParams result;
    result.lerp(a, b, t);
    return result;
  }
  static Vector apply_lens(const Vector &v, SurfaceLens lens) {
    return SDB::apply_lens(v, lens);
  }
  static float sample_function(Function function, const Complex &p,
                               const SourceState &source) {
    return SDB::sample_function(function, p, source);
  }
  static float sample_pattern(const Complex &p, float complexity, float mix,
                              float primary, float secondary) {
    return SDB::sample_pattern(p, complexity, mix, primary, secondary);
  }
  static const auto &presets() { return SDB::PRESETS; }
  static const auto &choreo() { return SDB::CHOREO; }
  static void make_triadic(uint32_t &hue, uint32_t sequence,
                           GenerativePalette &out) {
    SDB::next_triadic_palette(&hue, sequence, out);
  }
  static Pixel generated_color(const SDB &sdb, float value) {
    return sdb.generated_palette_cycler.palette().get(value).color;
  }
  static Pixel liquid_color(const SDB &sdb, float value) {
    return sdb.liquid_palette_cycler.palette().get(value).color;
  }
};

/** @brief Every named clock wraps in its native domain. */
inline void test_shadierball_clocks_wrapped() {
  using WB = ShadierBallWhiteBox;
  reset_effect_globals();
  hs::set_mock_time(0, 0);
  WB::SDB sdb;
  sdb.init();
  WB::seed_clocks(sdb, WB::noise_time_period() * 4.0f);
  for (int frame = 0; frame < 32; ++frame) {
    hs::set_mock_time(frame * FRAME_MS, frame * FRAME_US);
    sdb.draw_frame();
    sdb.advance_display();
    const WB::ClockState clocks = WB::clocks(sdb);
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

/** @brief Pullback stages preserve their typed order and metadata. */
inline void test_shadierball_pipeline_contract() {
  using WB = ShadierBallWhiteBox;
  reset_effect_globals();
  WB::SDB sdb;
  sdb.init();
  WB::FrameState frame = WB::frame(sdb);
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
      frame.params.warp.inner, frame);
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
      frame.params.warp.inner, frame);
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
inline void test_shadierball_legacy_spatial_slots() {
  using WB = ShadierBallWhiteBox;
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
  WB::SDB landmark_sdb;
  landmark_sdb.init();
  WB::FrameState landmark_frame = WB::frame(landmark_sdb);
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
  WB::FrameState frame = WB::frame(landmark_sdb);
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
  HS_EXPECT_EQ(end.coords.re, lensed.re);
  HS_EXPECT_EQ(end.coords.im, lensed.im);
}

/** @brief Legacy source functions retain their closed forms. */
inline void test_shadierball_legacy_sources() {
  using WB = ShadierBallWhiteBox;
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
      HS_EXPECT_EQ(
          WB::sample_function(WB::Function::SPIRAL, p, source),
          fast_sinf(radius - 3.0f * (azimuth + source.angle) - source.primary));
      const float b = -re * source.angle_sin + im * source.angle_cos;
      HS_EXPECT_EQ(WB::sample_function(WB::Function::GRID, p, source),
                   fast_sinf(rotated + source.primary) *
                       fast_cosf(b - source.primary));
    }
  }
}

/** @brief Coupled/direct sampling reduces to both authored formulas. */
inline void test_shadierball_coupled_source() {
  using WB = ShadierBallWhiteBox;
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
            HS_EXPECT_EQ(
                WB::sample_pattern(p, complexity, 0.0f, primary, secondary),
                coupled);
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
inline void test_shadierball_preset_bank() {
  using WB = ShadierBallWhiteBox;
  const auto &presets = WB::presets();
  const auto &choreo = WB::choreo();
  HS_EXPECT_EQ(presets.size(), size_t(20));
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
  HS_EXPECT_EQ(WB::work_signature(presets[17]).worst_case_work_units(),
               WB::hold_work_unit_ceiling());
  for (size_t index = 0; index < presets.size(); ++index) {
    const auto &preset = presets[index];
    if (index < 15)
      HS_EXPECT_TRUE(WB::slots_equal(preset.slots, WB::liquid_stereo_slots()));
    HS_EXPECT_TRUE(WB::valid_config(preset));
    if (index < 14)
      HS_EXPECT_TRUE(WB::slots_equal(preset.slots, presets[index + 1].slots));
  }

  reset_effect_globals();
  WB::SDB sdb;
  sdb.init();
  HS_EXPECT_TRUE(
      WB::slots_equal(WB::active_slots(sdb), WB::liquid_stereo_slots()));
}

/** @brief Whole-schema validation applies valid configs and rejects invalid. */
inline void test_shadierball_config_admission() {
  using WB = ShadierBallWhiteBox;
  {
    reset_effect_globals();
    WB::SDB sdb;
    sdb.init();
    const WB::Slots original = WB::active_slots(sdb);
    WB::RequestedConfig invalid_params = WB::active_config(sdb);
    invalid_params.params.source.pattern_freq = 21.0f;
    const WB::RequestedConfig before_invalid_params = WB::active_config(sdb);
    WB::request_config(sdb, invalid_params);
    HS_EXPECT_TRUE(WB::active_config(sdb) == before_invalid_params);
    HS_EXPECT_TRUE(WB::requested_slots(sdb) == original);
    HS_EXPECT_FALSE(WB::transition_active(sdb));
    HS_EXPECT_FALSE(WB::param_morph_active(sdb));

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
        WB::transition_admitted(WB::active_config(sdb), legacy_config));
    WB::request_config(sdb, legacy_config);
    HS_EXPECT_TRUE(WB::slots_equal(WB::active_slots(sdb), original));
    HS_EXPECT_TRUE(WB::transition_active(sdb));
    HS_EXPECT_TRUE(WB::transition_to_config(sdb) == legacy_config);
    HS_EXPECT_TRUE(WB::requested_slots(sdb) == legacy_config.slots);

    WB::RequestedConfig resource_from = WB::legacy_config();
    resource_from.slots.warp_program.outer.kind =
        WB::WarpStageKind::VECTOR_NOISE;
    resource_from.slots.warp_program.outer.basis = WB::NoiseBasis::SIMPLEX;
    resource_from.params.warp.outer.time_scale = 0.0f;
    WB::RequestedConfig resource_to = resource_from;
    resource_to.slots.warp_program.outer.basis = WB::NoiseBasis::FBM3;
    HS_EXPECT_TRUE(WB::valid_config(resource_from));
    HS_EXPECT_TRUE(WB::valid_config(resource_to));
    HS_EXPECT_TRUE(WB::transition_admitted(resource_from, resource_to));
    HS_EXPECT_EQ(WB::transition_mode(resource_from, resource_to),
                 WB::TransitionMode::THROUGH_CLEAR);
    resource_to.slots.warp_program.outer.resource_id ^= 4U;
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
    WB::SDB projection_change;
    projection_change.init();
    HS_EXPECT_TRUE(
        projection_change.updateParameter(
            "Projection", static_cast<float>(WB::Projection::BONNE)) ==
        ParamSetResult::APPLIED);
    projection_change.draw_frame();
    projection_change.advance_display();
    HS_EXPECT_TRUE(WB::transition_active(projection_change));
    HS_EXPECT_EQ(WB::transition_to_config(projection_change).slots.projection,
                 WB::Projection::BONNE);
    HS_EXPECT_EQ(WB::transition_to_config(projection_change)
                     .slots.warp_program.outer.kind,
                 WB::WarpStageKind::NONE);
  }

  {
    reset_effect_globals();
    WB::SDB legacy_warp_change;
    legacy_warp_change.init();
    HS_EXPECT_TRUE(
        legacy_warp_change.updateParameter(
            "Projection", static_cast<float>(WB::Projection::BONNE)) ==
        ParamSetResult::APPLIED);
    legacy_warp_change.draw_frame();
    legacy_warp_change.advance_display();
    while (WB::transition_active(legacy_warp_change)) {
      legacy_warp_change.draw_frame();
      legacy_warp_change.advance_display();
    }
    HS_EXPECT_TRUE(
        legacy_warp_change.updateParameter(
            "Outer Warp",
            static_cast<float>(WB::WarpStageKind::LEGACY_STEREO_NOISE)) ==
        ParamSetResult::APPLIED);
    legacy_warp_change.draw_frame();
    legacy_warp_change.advance_display();
    HS_EXPECT_EQ(WB::transition_to_config(legacy_warp_change).slots.projection,
                 WB::Projection::STEREOGRAPHIC);
  }
}

/** @brief Selector intent is independent of the live worker destination. */
inline void test_shadierball_deterministic_gui_edits() {
  using WB = ShadierBallWhiteBox;
  struct Result {
    std::array<WB::RequestedConfig, 3> configs;
    std::array<uint64_t, 3> schema_hashes{};
    std::array<uint32_t, 3> generation_deltas{};
  };
  auto schema_hash = [](const WB::SDB &sdb) {
    uint64_t hash = 1469598103934665603ULL;
    auto append = [&](const char *text) {
      for (; *text != '\0'; ++text) {
        hash ^= static_cast<uint8_t>(*text);
        hash *= 1099511628211ULL;
      }
    };
    for (const Effect::ParamDef &def : sdb.getParameters()) {
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
    WB::SDB sdb;
    sdb.init();
    Result result;
    const char *names[] = {"Projection", "Projection", "Outer Warp"};
    const float values[] = {
        static_cast<float>(WB::Projection::BONNE),
        static_cast<float>(WB::Projection::AIROCEAN),
        static_cast<float>(WB::WarpStageKind::LEGACY_STEREO_NOISE)};
    for (size_t index = 0; index < result.configs.size(); ++index) {
      const uint32_t before = sdb.getParameterSchemaGeneration();
      HS_EXPECT_TRUE(sdb.updateParameter(names[index], values[index]) ==
                     ParamSetResult::APPLIED);
      result.generation_deltas[index] =
          sdb.getParameterSchemaGeneration() - before;
      result.configs[index] = WB::requested_config(sdb);
      result.schema_hashes[index] = schema_hash(sdb);
      if (advance_worker && index + 1 < result.configs.size()) {
        sdb.draw_frame();
        sdb.advance_display();
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

/** @brief Hold cost and stable morph admission are explicit conservative proofs. */
inline void test_shadierball_work_admission() {
  using WB = ShadierBallWhiteBox;
  for (const auto &preset : WB::presets())
    HS_EXPECT_TRUE(WB::hold_admitted(preset));

  WB::RequestedConfig one = WB::legacy_config();
  one.slots.surface_lens = WB::SurfaceLens::NONE;
  one.params.surface_lens.mix = 0.0f;
  one.slots.warp_program.outer.kind = WB::WarpStageKind::CURL_FLOW;
  one.slots.warp_program.outer.basis = WB::NoiseBasis::FBM3;
  one.slots.warp_program.outer.curl_integrator = WB::CurlIntegrator::MIDPOINT_2;
  one.params.warp.outer.strength = 1e-4f;
  one.params.warp.outer.scale = 1.0f;
  one.params.warp.outer.time_scale = 0.0f;
  HS_EXPECT_TRUE(WB::valid_config(one));
  HS_EXPECT_EQ(WB::work_signature(one).worst_case_noise_calls(),
               WB::hold_work_unit_ceiling());
  HS_EXPECT_EQ(WB::work_signature(one).worst_case_work_units(),
               WB::hold_work_unit_ceiling());
  HS_EXPECT_TRUE(WB::hold_admitted(one));

  WB::RequestedConfig dual = one;
  dual.slots.warp_program.outer.basis = WB::NoiseBasis::RIDGED3;
  dual.slots.warp_program.outer.curl_integrator =
      WB::CurlIntegrator::MIDPOINT_4;
  dual.slots.warp_program.inner = dual.slots.warp_program.outer;
  dual.slots.warp_program.inner.resource_id = 1;
  dual.params.warp.inner = dual.params.warp.outer;
  HS_EXPECT_TRUE(WB::valid_config(dual));
  HS_EXPECT_EQ(WB::work_signature(dual).worst_case_noise_calls(),
               uint16_t(384));
  HS_EXPECT_FALSE(WB::hold_admitted(dual));
  dual.slots.projection = WB::Projection::BONNE;
  dual.slots.surface_lens = WB::SurfaceLens::GLITCH;
  dual.params.surface_lens.mix = 0.5f;
  HS_EXPECT_EQ(WB::work_signature(dual).worst_case_noise_calls(),
               uint16_t(768));
  HS_EXPECT_GT(WB::work_signature(dual).worst_case_work_units(), uint16_t(768));
  HS_EXPECT_FALSE(WB::hold_admitted(dual));

  for (WB::Projection projection :
       {WB::Projection::BONNE, WB::Projection::PEIRCE_QUINCUNCIAL,
        WB::Projection::AIROCEAN}) {
    WB::RequestedConfig strict = WB::legacy_config();
    strict.slots.projection = projection;
    strict.slots.warp_program.outer.kind = WB::WarpStageKind::NONE;
    HS_EXPECT_TRUE(WB::hold_admitted(strict));
  }
  WB::RequestedConfig airo_mobius = WB::legacy_config();
  airo_mobius.slots.projection = WB::Projection::AIROCEAN;
  airo_mobius.slots.warp_program.outer.kind = WB::WarpStageKind::NONE;
  airo_mobius.slots.surface_lens = WB::SurfaceLens::MOBIUS;
  airo_mobius.params.surface_lens.mix = 0.5f;
  HS_EXPECT_TRUE(WB::valid_config(airo_mobius));
  HS_EXPECT_EQ(WB::work_signature(airo_mobius).noise_calls, uint16_t(0));
  HS_EXPECT_GT(WB::work_signature(airo_mobius).worst_case_work_units(),
               WB::hold_work_unit_ceiling());
  HS_EXPECT_FALSE(WB::hold_admitted(airo_mobius));

  WB::RequestedConfig from = WB::legacy_config();
  from.slots.surface_lens = WB::SurfaceLens::NONE;
  from.params.surface_lens.mix = 0.0f;
  from.slots.warp_program.outer.kind = WB::WarpStageKind::CURL_FLOW;
  from.slots.warp_program.outer.basis = WB::NoiseBasis::SIMPLEX;
  from.slots.warp_program.outer.curl_integrator =
      WB::CurlIntegrator::MIDPOINT_4;
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
  const auto &presets = WB::presets();
  for (size_t index = 0; index < presets.size(); ++index) {
    const auto &next = presets[(index + 1) % presets.size()];
    if (WB::stable_topology(presets[index], next))
      HS_EXPECT_TRUE(WB::stable_parameter_path_admitted(presets[index], next));
  }
}

/** @brief Folded and interrupted projections reject unproved noise seams. */
inline void test_shadierball_strict_seam_admission() {
  using WB = ShadierBallWhiteBox;
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
inline void test_shadierball_additive_delta_precision() {
  using WB = ShadierBallWhiteBox;
  reset_effect_globals();
  WB::SDB sdb;
  sdb.init();
  WB::FrameState frame = WB::frame(sdb);
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
inline void test_shadierball_profile_presets() {
  using WB = ShadierBallWhiteBox;
  reset_effect_globals();
  WB::SDB sdb;
  sdb.init();
  const auto &presets = WB::presets();
  for (size_t index = 0; index < presets.size(); ++index) {
    sdb.profile_select_preset(index);
    HS_EXPECT_TRUE(WB::active_config(sdb) == presets[index]);
    HS_EXPECT_TRUE(WB::requested_config(sdb) == presets[index]);
    HS_EXPECT_TRUE(WB::hold_admitted(WB::active_config(sdb)));
    HS_EXPECT_FALSE(WB::transition_active(sdb));
    HS_EXPECT_FALSE(WB::param_morph_active(sdb));
  }
}

/** @brief Every GUI enum option is writable and survives its handoff. */
inline void test_shadierball_gui_catalog() {
  using WB = ShadierBallWhiteBox;
  reset_effect_globals();
  WB::SDB sdb;
  sdb.init();
  HS_EXPECT_LE(sdb.getParameters().size(), size_t(64));
  const auto *projection = sdb.getParameters().find("Projection");
  HS_EXPECT_TRUE(projection != nullptr);
  HS_EXPECT_EQ(projection->option_count, 6);
  HS_EXPECT_TRUE(std::strcmp(projection->options[3], "Bonne") == 0);
  HS_EXPECT_TRUE(std::strcmp(projection->options[4], "Peirce Quincuncial") ==
                 0);
  HS_EXPECT_TRUE(std::strcmp(projection->options[5], "Dymaxion / Airocean") ==
                 0);
  const uint32_t schema_before = sdb.getParameterSchemaGeneration();
  HS_EXPECT_TRUE(sdb.updateParameter(
                     "Projection", static_cast<float>(WB::Projection::BONNE)) ==
                 ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(sdb.getParameterSchemaGeneration() > schema_before);
  HS_EXPECT_TRUE(sdb.getParameters().find("Bonne Hemisphere") != nullptr);
  sdb.draw_frame();
  sdb.advance_display();
  WB::settle_transition(sdb);

  constexpr const char *ROOT_ENUMS[] = {
      "Function",   "Projection", "Projection Frame", "Lens",
      "Outer Warp", "Inner Warp", "Signal Weight",    "Value Transfer",
      "Coverage",   "Colorizer"};
  for (const char *name : ROOT_ENUMS) {
    const int option_count = sdb.getParameters().find(name)->option_count;
    for (int option = 0; option < option_count; ++option) {
      if (std::strcmp(name, "Inner Warp") == 0 && option == 1) {
        HS_EXPECT_TRUE(sdb.updateParameter("Outer Warp", 0.0f) ==
                       ParamSetResult::APPLIED);
        sdb.draw_frame();
        sdb.advance_display();
        WB::settle_transition(sdb);
      } else if (std::strcmp(name, "Outer Warp") == 0 && option == 1) {
        HS_EXPECT_TRUE(sdb.updateParameter("Inner Warp", 0.0f) ==
                       ParamSetResult::APPLIED);
        sdb.draw_frame();
        sdb.advance_display();
        WB::settle_transition(sdb);
      }
      HS_EXPECT_TRUE(sdb.updateParameter(name, static_cast<float>(option)) ==
                     ParamSetResult::APPLIED);
      sdb.draw_frame();
      sdb.advance_display();
      WB::settle_transition(sdb);
      HS_EXPECT_EQ(sdb.getParameters().find(name)->get(),
                   static_cast<float>(option));
      HS_EXPECT_FALSE(WB::transition_active(sdb));
      HS_EXPECT_FALSE(WB::param_morph_active(sdb));
    }
  }

  auto select_and_set_all = [&](const char *root, int selection,
                                const char *subordinate) {
    HS_EXPECT_TRUE(sdb.updateParameter(root, static_cast<float>(selection)) ==
                   ParamSetResult::APPLIED);
    sdb.draw_frame();
    sdb.advance_display();
    WB::settle_transition(sdb);
    const auto *def = sdb.getParameters().find(subordinate);
    HS_EXPECT(def != nullptr, subordinate);
    if (def == nullptr)
      return;
    const int count = def->option_count;
    for (int option = 0; option < count; ++option) {
      HS_EXPECT_TRUE(
          sdb.updateParameter(subordinate, static_cast<float>(option)) ==
          ParamSetResult::APPLIED);
      sdb.draw_frame();
      sdb.advance_display();
      WB::settle_transition(sdb);
      HS_EXPECT_EQ(sdb.getParameters().find(subordinate)->get(),
                   static_cast<float>(option));
    }
  };
  select_and_set_all("Projection", 2, "Gnomonic Hemisphere");
  select_and_set_all("Projection", 3, "Bonne Hemisphere");
  HS_EXPECT_TRUE(sdb.updateParameter("Bonne Standard Parallel", 0.9f) ==
                 ParamSetResult::APPLIED);
  sdb.draw_frame();
  sdb.advance_display();
  WB::settle_transition(sdb);
  HS_EXPECT_NEAR(
      WB::active_config(sdb).params.projection.bonne_standard_parallel, 0.9f,
      1e-6f);
  select_and_set_all("Projection", 4, "Peirce Layout");
  for (int layout = 0; layout < 4; ++layout) {
    HS_EXPECT_TRUE(
        sdb.updateParameter("Peirce Layout", static_cast<float>(layout)) ==
        ParamSetResult::APPLIED);
    const bool has_scroll =
        sdb.getParameters().find("Projection Layout Scroll") != nullptr;
    HS_EXPECT_EQ(has_scroll, layout >= 2);
    sdb.draw_frame();
    sdb.advance_display();
    WB::settle_transition(sdb);
  }
  select_and_set_all("Projection", 5, "Airocean Layout");
  select_and_set_all("Outer Warp", 5, "Outer Noise Basis");
  select_and_set_all("Outer Warp", 5, "Outer Warp Envelope");
  HS_EXPECT_TRUE(sdb.updateParameter("Inner Warp", 0.0f) ==
                 ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(
      sdb.updateParameter("Outer Warp",
                          static_cast<float>(WB::WarpStageKind::CURL_FLOW)) ==
      ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(
      sdb.updateParameter("Outer Noise Basis",
                          static_cast<float>(WB::NoiseBasis::SIMPLEX)) ==
      ParamSetResult::APPLIED);
  sdb.draw_frame();
  sdb.advance_display();
  WB::settle_transition(sdb);
  select_and_set_all("Outer Warp", 6, "Outer Curl Integrator");
  select_and_set_all("Outer Warp", 8, "Outer Polar Mode");
  select_and_set_all("Outer Warp", 8, "Outer Polar Harmonic");
  HS_EXPECT_TRUE(sdb.getParameters().find("Pattern Freq") == nullptr);
  HS_EXPECT_TRUE(sdb.updateParameter("Function",
                                     static_cast<float>(WB::Function::RINGS)) ==
                 ParamSetResult::APPLIED);
  HS_EXPECT_EQ(WB::requested_config(sdb).slots.warp_program.outer.kind,
               WB::WarpStageKind::NONE);
  HS_EXPECT_TRUE(sdb.getParameters().find("Pattern Freq") != nullptr);
  select_and_set_all("Value Transfer", 3, "Band Count");
}

/** @brief New cartographic kernels preserve landmarks and stay finite. */
inline void test_shadierball_projection_catalog() {
  using WB = ShadierBallWhiteBox;
  const float standard_parallel = PI_F * 0.25f;
  const Vector bonne_origin(cosf(standard_parallel), sinf(standard_parallel),
                            0.0f);
  const auto bonne =
      shadierball::bonne_projection(bonne_origin, 0.0f, standard_parallel);
  HS_EXPECT_NEAR(bonne.coords.re, 0.0f, 2e-5f);
  HS_EXPECT_NEAR(bonne.coords.im, 0.0f, 2e-5f);

  const auto peirce = shadierball::peirce_projection(UP, 0.0f, 1, 0.0f);
  HS_EXPECT_NEAR(peirce.coords.re, 0.0f, 2e-5f);
  HS_EXPECT_NEAR(peirce.coords.im, 0.0f, 2e-5f);

  const auto &center = shadierball::AIROCEAN_CENTERS[0];
  const auto airocean = shadierball::airocean_projection(
      Vector(center.x, center.z, center.y), 0.0f, false);
  const auto &triangle = shadierball::AIROCEAN_PLANAR_FACES[0];
  HS_EXPECT_NEAR(airocean.coords.re,
                 (triangle[0].x + triangle[1].x + triangle[2].x) / 3.0f, 2e-4f);
  HS_EXPECT_NEAR(airocean.coords.im,
                 (triangle[0].y + triangle[1].y + triangle[2].y) / 3.0f, 2e-4f);
  HS_EXPECT_EQ(airocean.region_id, uint8_t(0));

  const Vector equator_zero(1.0f, 0.0f, 0.0f);
  const Vector equator_east(0.0f, 0.0f, 1.0f);
  const Vector antimeridian(-1.0f, 0.0f, 0.0f);
  const auto bonne_equator =
      shadierball::bonne_projection(equator_zero, 0.0f, standard_parallel);
  HS_EXPECT_NEAR(bonne_equator.coords.re, 0.0f, 2e-5f);
  HS_EXPECT_NEAR(bonne_equator.coords.im, -0.7853981634f, 2e-5f);
  const auto bonne_east =
      shadierball::bonne_projection(equator_east, 0.0f, standard_parallel);
  HS_EXPECT_NEAR(bonne_east.coords.re, 1.3758501640f, 3e-5f);
  HS_EXPECT_NEAR(bonne_east.coords.im, -0.1378413458f, 3e-5f);
  const auto bonne_shifted = shadierball::bonne_projection(
      equator_east, 0.5f * PI_F, standard_parallel);
  HS_EXPECT_NEAR(bonne_shifted.coords.re, 0.0f, 2e-5f);
  HS_EXPECT_NEAR(bonne_shifted.coords.im, -0.7853981634f, 2e-5f);
  const auto bonne_cut =
      shadierball::bonne_projection(antimeridian, 0.0f, standard_parallel);
  HS_EXPECT_NEAR(bonne_cut.fade_edge_distance, 0.0f, 2e-5f);
  const auto werner =
      shadierball::bonne_projection(equator_east, 0.0f, 0.5f * PI_F);
  HS_EXPECT_NEAR(werner.coords.re, 1.3217795320f, 3e-5f);
  HS_EXPECT_NEAR(werner.coords.im, -0.8487048774f, 3e-5f);

  constexpr float PEIRCE_K = 1.8540746773013719f;
  const Vector south_pole(0.0f, -1.0f, 0.0f);
  const auto peirce_diamond =
      shadierball::peirce_projection(south_pole, 0.0f, 0, 0.0f);
  const auto peirce_square =
      shadierball::peirce_projection(south_pole, 0.0f, 1, 0.0f);
  const auto peirce_horizontal =
      shadierball::peirce_projection(south_pole, 0.0f, 2, 0.0f);
  const auto peirce_vertical =
      shadierball::peirce_projection(south_pole, 0.0f, 3, 0.0f);
  HS_EXPECT_NEAR(peirce_diamond.coords.re, 2.0f * PEIRCE_K, 3e-5f);
  HS_EXPECT_NEAR(peirce_diamond.coords.im, 0.0f, 3e-5f);
  HS_EXPECT_NEAR(peirce_square.coords.re, PEIRCE_K * 1.4142135624f, 3e-5f);
  HS_EXPECT_NEAR(peirce_square.coords.im, PEIRCE_K * 1.4142135624f, 3e-5f);
  HS_EXPECT_NEAR(peirce_horizontal.coords.re, PEIRCE_K, 3e-5f);
  HS_EXPECT_NEAR(peirce_vertical.coords.im, PEIRCE_K, 3e-5f);
  const auto peirce_scroll0 =
      shadierball::peirce_projection(equator_east, 0.0f, 2, 0.0f);
  const auto peirce_scroll1 =
      shadierball::peirce_projection(equator_east, 0.0f, 2, 1.0f);
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
  const auto peirce_scroll_mid = shadierball::peirce_projection(
      equator_east, 0.0f, 2, scroll_mid.layout_scroll);
  HS_EXPECT_NEAR(peirce_scroll_mid.coords.re, peirce_scroll0.coords.re, 3e-5f);
  HS_EXPECT_NEAR(peirce_scroll_mid.coords.im, peirce_scroll0.coords.im, 3e-5f);
  const auto peirce_zero_fade =
      shadierball::peirce_projection(equator_zero, 0.0f, 1, 0.0f);
  HS_EXPECT_NEAR(peirce_zero_fade.fade_edge_distance, 0.0f, 2e-5f);

  auto lon_lat = [](float longitude, float latitude) {
    const float cp = cosf(latitude);
    return Vector(cp * cosf(longitude), sinf(latitude), cp * sinf(longitude));
  };
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
    const auto before = shadierball::peirce_projection(
        lon_lat(tie.longitude - TIE_EPS, SOUTH_LATITUDE), 0.0f, 0, 0.0f);
    const auto exact = shadierball::peirce_projection(
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
        shadierball::peirce_projection(on_meridian, MERIDIAN, layout, 0.13f);
    const auto reference =
        shadierball::peirce_projection(zero_meridian, 0.0f, layout, 0.13f);
    HS_EXPECT_NEAR(shifted.coords.re, reference.coords.re, 2e-5f);
    HS_EXPECT_NEAR(shifted.coords.im, reference.coords.im, 2e-5f);
    HS_EXPECT_EQ(shifted.region_id, reference.region_id);
  }
  const auto airocean_shifted =
      shadierball::airocean_projection(on_meridian, MERIDIAN, false);
  const auto airocean_reference =
      shadierball::airocean_projection(zero_meridian, 0.0f, false);
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
      shadierball::airocean_projection(oracle_point, 0.0f, false);
  HS_EXPECT_NEAR(airocean_oracle.coords.re, 2.1265288136f, 4e-5f);
  HS_EXPECT_NEAR(airocean_oracle.coords.im, 3.6817439808f, 4e-5f);
  const auto airocean_horizontal =
      shadierball::airocean_projection(oracle_point, 0.0f, true);
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
        shadierball::airocean_projection(oracle.point, 0.0f, false);
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
  const float glued_cut_distances[] = {0.525731112119133f, 0.525731112159295f,
                                       0.525731112159296f};
  for (size_t index = 0; index < std::size(glued_points); ++index) {
    const auto mapped =
        shadierball::airocean_projection(glued_points[index], 0.0f, false);
    HS_EXPECT_EQ(mapped.region_id, glued_faces[index]);
    HS_EXPECT_EQ(mapped.edge_class,
                 shadierball::airocean_edge_identity(glued_faces[index],
                                                     glued_edges[index]));
    HS_EXPECT_TRUE((mapped.traits & shadierball::projection_traits(
                                        shadierball::ProjectionTrait::GLUED)) !=
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
        shadierball::airocean_projection(cut_points[index], 0.0f, false);
    HS_EXPECT_EQ(mapped.region_id, cut_faces[index]);
    HS_EXPECT_EQ(mapped.edge_class, shadierball::airocean_edge_identity(
                                        cut_faces[index], cut_edges[index]));
    HS_EXPECT_TRUE((mapped.traits & shadierball::projection_traits(
                                        shadierball::ProjectionTrait::CUT)) !=
                   0);
    HS_EXPECT_NEAR(mapped.coords.re, cut_coords[index].re, 5e-5f);
    HS_EXPECT_NEAR(mapped.coords.im, cut_coords[index].im, 5e-5f);
    HS_EXPECT_NEAR(mapped.fade_edge_distance, cut_distances[index], 5e-5f);
  }
  for (uint8_t face = 0; face < 23; ++face) {
    const auto &face_center = shadierball::AIROCEAN_CENTERS[face];
    const auto mapped = shadierball::airocean_projection(
        Vector(face_center.x, face_center.z, face_center.y), 0.0f, false);
    HS_EXPECT_EQ(mapped.region_id, face);
  }
  HS_EXPECT_TRUE(shadierball::airocean_edge_is_cut(14, 0));
  HS_EXPECT_EQ(shadierball::airocean_edge_identity(14, 0), uint8_t(42));
  HS_EXPECT_EQ(shadierball::airocean_edge_identity(18, 0), uint8_t(54));

  for (int latitude_index = -8; latitude_index <= 8; ++latitude_index) {
    const float latitude = latitude_index * (PI_F / 18.0f);
    for (int longitude_index = -18; longitude_index < 18; ++longitude_index) {
      const float longitude = longitude_index * (PI_F / 18.0f);
      const float cp = cosf(latitude);
      const Vector v(cp * cosf(longitude), sinf(latitude),
                     cp * sinf(longitude));
      const auto b = shadierball::bonne_projection(v, 0.37f, standard_parallel);
      HS_EXPECT_TRUE(std::isfinite(b.coords.re));
      HS_EXPECT_TRUE(std::isfinite(b.coords.im));
      HS_EXPECT_TRUE(std::isfinite(b.fade_edge_distance));
      for (uint8_t layout = 0; layout < 4; ++layout) {
        const auto p = shadierball::peirce_projection(v, -0.21f, layout, 0.13f);
        HS_EXPECT_TRUE(std::isfinite(p.coords.re));
        HS_EXPECT_TRUE(std::isfinite(p.coords.im));
        HS_EXPECT_TRUE(std::isfinite(p.fade_edge_distance));
      }
      for (bool horizontal : {false, true}) {
        const auto a = shadierball::airocean_projection(v, 0.19f, horizontal);
        HS_EXPECT_TRUE(std::isfinite(a.coords.re));
        HS_EXPECT_TRUE(std::isfinite(a.coords.im));
        HS_EXPECT_TRUE(std::isfinite(a.fade_edge_distance));
        HS_EXPECT_TRUE(a.region_id < 23);
      }
    }
  }
}

/** @brief Domain policies, gauges, and analytic admission reject unsafe tuples. */
inline void test_shadierball_projection_and_admission_contracts() {
  using WB = ShadierBallWhiteBox;
  const uint8_t periodic_traits = shadierball::projection_traits(
      shadierball::ProjectionTrait::GLUED, shadierball::ProjectionTrait::FOLDED,
      shadierball::ProjectionTrait::PERIODIC);
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

  auto kernel_lookup = [](const shadierball::ProjectionKernelResult &result) {
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
      shadierball::peirce_projection(sector_before, 0.0f, 0, 0.0f));
  const auto diamond_after = kernel_lookup(
      shadierball::peirce_projection(sector_after, 0.0f, 0, 0.0f));
  HS_EXPECT_FALSE(WB::join_compatible(diamond_before, diamond_after,
                                      WB::Projection::PEIRCE_QUINCUNCIAL));
  const auto horizontal_before = kernel_lookup(
      shadierball::peirce_projection(sector_before, 0.0f, 2, 0.0f));
  const auto horizontal_after = kernel_lookup(
      shadierball::peirce_projection(sector_after, 0.0f, 2, 0.0f));
  HS_EXPECT_FALSE(WB::join_compatible(horizontal_before, horizontal_after,
                                      WB::Projection::PEIRCE_QUINCUNCIAL));

  const auto glued_side_a = kernel_lookup(shadierball::airocean_projection(
      Vector(0.00321964224570043f, 0.902105871397194f, 0.431502758617508f),
      0.0f, false));
  const auto glued_side_b = kernel_lookup(shadierball::airocean_projection(
      Vector(0.00321096516044884f, 0.902110786482306f, 0.431492547577607f),
      0.0f, false));
  HS_EXPECT_FALSE(WB::join_compatible(glued_side_a, glued_side_b,
                                      WB::Projection::AIROCEAN));
  const auto cut_side_a = kernel_lookup(shadierball::airocean_projection(
      Vector(0.456076125629434f, 0.767836786220573f, -0.449912477441232f), 0.0f,
      false));
  const auto cut_side_b = kernel_lookup(shadierball::airocean_projection(
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
  WB::SDB sdb;
  sdb.init();
  WB::FrameState frame = WB::frame(sdb);
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
inline void test_shadierball_kernel_catalog() {
  using WB = ShadierBallWhiteBox;
  reset_effect_globals();
  WB::SDB sdb;
  sdb.init();
  const Vector view = Vector(0.31f, 0.87f, -0.38f).normalized();
  WB::RequestedConfig config = WB::legacy_config();
  config.slots.surface_lens = WB::SurfaceLens::NONE;
  config.slots.warp_program.outer.kind = WB::WarpStageKind::NONE;
  config.slots.warp_program.inner.kind = WB::WarpStageKind::NONE;
  config.slots.coverage = WB::CoveragePolicy::OPAQUE;

  auto check = [&](const WB::RequestedConfig &candidate) {
    HS_EXPECT_TRUE(WB::valid_config(candidate));
    WB::request_config(sdb, candidate);
    const WB::RequestedConfig canonical = WB::requested_config(sdb);
    HS_EXPECT_TRUE(canonical.slots == candidate.slots);
    WB::settle_transition(sdb);
    HS_EXPECT_TRUE(WB::active_config(sdb) == canonical);
    const Color4 color = WB::shade(view, WB::frame(sdb));
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
    config.slots.warp_program.outer.basis = WB::NoiseBasis::RIDGED3;
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
    for (uint8_t coverage = 0; coverage <= 3; ++coverage) {
      config.slots.coverage = static_cast<WB::CoveragePolicy>(coverage);
      for (uint8_t colorizer = 0; colorizer <= 2; ++colorizer) {
        config.slots.colorizer = static_cast<WB::Colorizer>(colorizer);
        check(config);
      }
    }
  }

  WB::FrameState frame = WB::frame(sdb);
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
inline void test_shadierball_stable_preset_transition() {
  using WB = ShadierBallWhiteBox;
  reset_effect_globals();
  WB::SDB sdb;
  sdb.init();
  const WB::Params from = WB::presets()[0].params;
  const WB::Params to = WB::presets()[1].params;
  const size_t initial_index = WB::preset_index(sdb);
  WB::begin_blend(sdb);
  HS_EXPECT_EQ(WB::preset_index(sdb), initial_index + 1);
  HS_EXPECT_TRUE(WB::param_morph_active(sdb));
  HS_EXPECT_FALSE(WB::transition_active(sdb));
  HS_EXPECT_TRUE(WB::live_params(sdb) == from);

  float previous_phase = WB::clocks(sdb).source_primary;
  for (int step = 0; step <= 60; ++step) {
    WB::step_param_morph(sdb);
    const float live_speed = WB::live_params(sdb).source.speed;
    const float phase = WB::clocks(sdb).source_primary;
    HS_EXPECT_NEAR(phase, fmodf(previous_phase + live_speed, TWO_PI_F), 1e-6f);
    previous_phase = phase;
    if (WB::param_morph_elapsed(sdb) == 6) {
      HS_EXPECT_LT(WB::live_params(sdb).source.pattern_freq,
                   from.source.pattern_freq);
      HS_EXPECT_GT(WB::live_params(sdb).source.pattern_freq,
                   to.source.pattern_freq);
      HS_EXPECT_EQ(WB::live_params(sdb).warp, from.warp);
    }
    if (WB::param_morph_elapsed(sdb) == 31) {
      HS_EXPECT_LT(WB::live_params(sdb).source.pattern_freq,
                   from.source.pattern_freq);
      HS_EXPECT_GT(WB::live_params(sdb).source.pattern_freq,
                   to.source.pattern_freq);
    }
  }
  HS_EXPECT_FALSE(WB::param_morph_active(sdb));
  HS_EXPECT_TRUE(WB::live_params(sdb) == to);

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
inline void test_shadierball_discrete_transition() {
  using WB = ShadierBallWhiteBox;

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
    WB::SDB sdb;
    sdb.init();
    const Vector view(0.2f, 0.9f, -0.3f);
    const WB::FrameState valid = WB::frame(sdb);
    WB::FrameState inactive = valid;
    inactive.resources.liquid_palette = nullptr;
    const Color4 expected = WB::shade(view, valid);
    const Color4 at_start =
        WB::shade_dual_output(view, {valid, inactive, 0.0f});
    HS_EXPECT_TRUE(at_start.color == expected.color);
    HS_EXPECT_EQ(at_start.alpha, expected.alpha);
    const Color4 at_end = WB::shade_dual_output(view, {inactive, valid, 1.0f});
    HS_EXPECT_TRUE(at_end.color == expected.color);
    HS_EXPECT_EQ(at_end.alpha, expected.alpha);
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
    WB::SDB sdb;
    sdb.init();
    sdb.setAnimationsPaused(true);
    WB::RequestedConfig from = WB::legacy_config();
    from.slots.warp_program.outer.kind = WB::WarpStageKind::VECTOR_NOISE;
    from.slots.warp_program.outer.basis = WB::NoiseBasis::SIMPLEX;
    from.params.warp.outer.time_scale = 0.0f;
    WB::request_config(sdb, from);
    WB::settle_transition(sdb);
    HS_EXPECT_EQ(WB::prepared_noise_basis(sdb, 0), WB::NoiseBasis::SIMPLEX);

    WB::RequestedConfig to = from;
    to.slots.warp_program.outer.basis = WB::NoiseBasis::FBM3;
    WB::request_config(sdb, to);
    HS_EXPECT_TRUE(WB::transition_active(sdb));
    HS_EXPECT_EQ(WB::transition_mode(WB::transition_from_config(sdb),
                                     WB::transition_to_config(sdb)),
                 WB::TransitionMode::THROUGH_CLEAR);
    for (int frame = 0; frame < 30; ++frame) {
      sdb.draw_frame();
      sdb.advance_display();
    }
    HS_EXPECT_EQ(WB::transition_elapsed(sdb), uint16_t(30));
    HS_EXPECT_EQ(WB::prepared_noise_basis(sdb, 0), WB::NoiseBasis::SIMPLEX);
    sdb.draw_frame();
    sdb.advance_display();
    HS_EXPECT_EQ(WB::transition_elapsed(sdb), uint16_t(31));
    HS_EXPECT_EQ(WB::prepared_noise_basis(sdb, 0), WB::NoiseBasis::FBM3);
  }

  {
    reset_effect_globals();
    WB::SDB sdb;
    sdb.init();
    WB::LookRuntime liquid_runtime;
    WB::LookRuntime generated_runtime;
    const WB::WalkDeltas deltas{make_rotation(Y_AXIS, 0.2f),
                                make_rotation(X_AXIS, 0.3f)};
    WB::advance_runtime(sdb, liquid_runtime, WB::presets()[0], deltas);
    WB::advance_runtime(sdb, generated_runtime, WB::legacy_config(), deltas);
    HS_EXPECT_EQ(liquid_runtime.clocks.source_primary, 0.1f);
    HS_EXPECT_EQ(generated_runtime.clocks.source_primary, 0.05f);
    HS_EXPECT_TRUE(liquid_runtime.projection_wander !=
                   generated_runtime.projection_wander);
    HS_EXPECT_TRUE(liquid_runtime.outer_wander !=
                   generated_runtime.outer_wander);

    const size_t original_index = WB::preset_index(sdb);
    WB::force_transition(sdb, WB::legacy_config(), 60, true);
    const WB::RequestedConfig captured_source = WB::transition_from_config(sdb);
    HS_EXPECT_TRUE(WB::transition_active(sdb));
    HS_EXPECT_EQ(WB::transition_mix(sdb), 0.0f);
    WB::begin_blend(sdb);
    HS_EXPECT_EQ(WB::preset_index(sdb), original_index);
    const uint32_t walk_steps = WB::walk_steps(sdb);
    const uint32_t liquid_steps = WB::liquid_palette_steps(sdb);
    const uint32_t generated_steps = WB::generated_palette_steps(sdb);
    sdb.draw_frame();
    sdb.advance_display();
    HS_EXPECT_EQ(WB::walk_steps(sdb), walk_steps + 1);
    HS_EXPECT_EQ(WB::liquid_palette_steps(sdb), liquid_steps + 1);
    HS_EXPECT_EQ(WB::generated_palette_steps(sdb), generated_steps + 1);
    HS_EXPECT_EQ(WB::transition_from_runtime(sdb).clocks.source_primary, 0.1f);
    HS_EXPECT_EQ(WB::transition_to_runtime(sdb).clocks.source_primary, 0.05f);

    for (int frame_index = 1; frame_index < 20; ++frame_index) {
      sdb.draw_frame();
      sdb.advance_display();
    }
    const uint16_t elapsed_before_takeover = WB::transition_elapsed(sdb);
    sdb.setAnimationsPaused(true);
    const WB::RequestedConfig original_destination =
        WB::transition_to_config(sdb);
    WB::RequestedConfig queued = WB::legacy_config();
    queued.slots.function = WB::Function::RINGS;
    HS_EXPECT_TRUE(WB::valid_config(queued));
    HS_EXPECT_TRUE(WB::transition_admitted(captured_source, queued));
    WB::request_config(sdb, queued);
    HS_EXPECT_TRUE(WB::transition_from_config(sdb) == captured_source);
    HS_EXPECT_TRUE(WB::transition_to_config(sdb) == original_destination);
    HS_EXPECT_FALSE(WB::transition_continues_choreo(sdb));
    HS_EXPECT_EQ(WB::transition_elapsed(sdb), elapsed_before_takeover);
    HS_EXPECT_TRUE(WB::requested_slots(sdb) == queued.slots);
    HS_EXPECT_TRUE(WB::pending_transition_active(sdb));
    HS_EXPECT_TRUE(WB::pending_transition_config(sdb) == queued);
    sdb.draw_frame();
    sdb.advance_display();
    HS_EXPECT_EQ(WB::transition_elapsed(sdb), elapsed_before_takeover + 1);

    const uint16_t elapsed_before_retarget = WB::transition_elapsed(sdb);
    const WB::RequestedConfig manual = WB::presets()[6];
    HS_EXPECT_TRUE(WB::transition_admitted(captured_source, manual));
    WB::request_config(sdb, manual);
    HS_EXPECT_TRUE(WB::transition_active(sdb));
    HS_EXPECT_TRUE(WB::transition_from_config(sdb) == captured_source);
    HS_EXPECT_TRUE(WB::transition_to_config(sdb) == original_destination);
    HS_EXPECT_FALSE(WB::transition_continues_choreo(sdb));
    HS_EXPECT_EQ(WB::transition_elapsed(sdb), elapsed_before_retarget);
    HS_EXPECT_TRUE(WB::requested_slots(sdb) == manual.slots);
    HS_EXPECT_TRUE(WB::pending_transition_active(sdb));
    HS_EXPECT_TRUE(WB::pending_transition_config(sdb) == manual);

    while (WB::transition_active(sdb) &&
           WB::transition_to_config(sdb) == original_destination) {
      sdb.draw_frame();
      sdb.advance_display();
    }
    HS_EXPECT_TRUE(WB::transition_active(sdb));
    HS_EXPECT_TRUE(WB::transition_from_config(sdb) == original_destination);
    HS_EXPECT_TRUE(WB::transition_to_config(sdb) == manual);
    HS_EXPECT_EQ(WB::transition_mix(sdb), 0.0f);

    while (WB::transition_active(sdb) && WB::transition_mix(sdb) < 1.0f) {
      sdb.draw_frame();
      sdb.advance_display();
    }
    HS_EXPECT_TRUE(WB::transition_active(sdb));
    HS_EXPECT_EQ(WB::transition_mix(sdb), 1.0f);
    const float destination_before_endpoint =
        WB::transition_to_runtime(sdb).clocks.source_primary;
    sdb.draw_frame();
    sdb.advance_display();
    HS_EXPECT_FALSE(WB::transition_active(sdb));
    HS_EXPECT_TRUE(WB::active_config(sdb) == manual);
    const float committed_phase = WB::clocks(sdb).source_primary;
    HS_EXPECT_NEAR(
        committed_phase,
        fmodf(destination_before_endpoint + manual.params.source.speed,
              TWO_PI_F),
        1e-6f);
    sdb.draw_frame();
    sdb.advance_display();
    HS_EXPECT_NEAR(
        WB::clocks(sdb).source_primary,
        fmodf(committed_phase + manual.params.source.speed, TWO_PI_F), 1e-6f);
  }
}

/** @brief Generated and liquid color resources remain independent owners. */
inline void test_shadierball_palette_resources() {
  using WB = ShadierBallWhiteBox;
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
  WB::SDB sdb;
  sdb.init();
  const Pixel generated = WB::generated_color(sdb, 0.25f);
  const Pixel liquid = WB::liquid_color(sdb, 0.25f);
  HS_EXPECT_TRUE(generated.r != liquid.r || generated.g != liquid.g ||
                 generated.b != liquid.b);
}

/** @brief Module entry point for ShadierBall contract tests. */
inline int run_shadierball_tests() {
  ModuleFixture fixture("shadierball");
  test_shadierball_clocks_wrapped();
  test_shadierball_pipeline_contract();
  test_shadierball_legacy_spatial_slots();
  test_shadierball_legacy_sources();
  test_shadierball_coupled_source();
  test_shadierball_preset_bank();
  test_shadierball_config_admission();
  test_shadierball_deterministic_gui_edits();
  test_shadierball_work_admission();
  test_shadierball_strict_seam_admission();
  test_shadierball_additive_delta_precision();
  test_shadierball_profile_presets();
  test_shadierball_gui_catalog();
  test_shadierball_projection_catalog();
  test_shadierball_projection_and_admission_contracts();
  test_shadierball_kernel_catalog();
  test_shadierball_stable_preset_transition();
  test_shadierball_discrete_transition();
  test_shadierball_palette_resources();
  return fixture.result();
}

} // namespace shadierball_tests
} // namespace hs_test
