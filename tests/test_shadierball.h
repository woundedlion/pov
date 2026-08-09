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
  using ProjectionFramePolicy = SDB::ProjectionFramePolicy;
  using SurfaceLens = SDB::SurfaceLens;
  using WarpStageKind = SDB::WarpStageKind;
  using WarpStageSpec = SDB::WarpStageSpec;
  using WarpStageParams = SDB::WarpStageParams;
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
  using TransitionFrame = SDB::TransitionFrame;

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
    sdb.apply_requested_config();
  }
  static Slots requested_slots(const SDB &sdb) {
    return sdb.requested_config.slots;
  }
  static void request_config(SDB &sdb, const RequestedConfig &config) {
    sdb.requested_config = config;
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
  static TransitionFrame transition_frame(const SDB &sdb) {
    return sdb.prepare_transition_frame();
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
  static Color4 shade_transition(const Vector &view,
                                 const TransitionFrame &frame) {
    return SDB::shade_transition(view, frame);
  }
  static void begin_blend(SDB &sdb) { sdb.begin_blend(); }
  static void step_param_morph(SDB &sdb) {
    sdb.prepare_param_morph();
    sdb.advance_runtime(sdb.runtime, {sdb.active_slots, sdb.blend.params},
                        {Quaternion(), Quaternion()});
    sdb.finish_transitions();
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
          clocks.projection_spin, clocks.breathe_phase}) {
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

  frame.resources.warp_noise = nullptr;
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
    HS_EXPECT_EQ(equirect.re, std::fabs(fast_atan2(v.z, v.x)) * radius);
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

/** @brief Presets are phase-grouped and retain the 15-look choreography. */
inline void test_shadierball_preset_bank() {
  using WB = ShadierBallWhiteBox;
  const auto &presets = WB::presets();
  const auto &choreo = WB::choreo();
  HS_EXPECT_EQ(presets.size(), size_t(15));
  HS_EXPECT_EQ(choreo.size(), presets.size());
  HS_EXPECT_EQ(presets[0].params.source.pattern_freq, 5.0f);
  HS_EXPECT_EQ(presets[0].params.warp.outer.scale, 3.0f);
  HS_EXPECT_EQ(presets[0].params.projection.pole_fade, 1.4f);
  HS_EXPECT_EQ(presets[0].params.surface_lens.mix, 1.0f);
  HS_EXPECT_EQ(presets[0].params.colorizer.breathe_depth, 0.15f);
  HS_EXPECT_EQ(presets[0].params.outer_camera.wander, 1.0f);
  HS_EXPECT_TRUE(choreo[0].staggered);
  HS_EXPECT_FALSE(choreo[6].staggered);
  for (size_t index = 0; index < presets.size(); ++index) {
    const auto &preset = presets[index];
    HS_EXPECT_TRUE(WB::slots_equal(preset.slots, WB::liquid_stereo_slots()));
    HS_EXPECT_TRUE(WB::valid_config(preset));
    HS_EXPECT_TRUE(WB::slots_equal(
        preset.slots, presets[(index + 1) % presets.size()].slots));
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
  reset_effect_globals();
  WB::SDB sdb;
  sdb.init();
  const WB::Slots original = WB::active_slots(sdb);
  WB::Slots invalid = original;
  invalid.projection = WB::Projection::EQUIRECTANGULAR;
  invalid.warp_program.outer.kind = WB::WarpStageKind::LEGACY_STEREO_NOISE;
  WB::request_slots(sdb, invalid);
  HS_EXPECT_TRUE(WB::slots_equal(WB::active_slots(sdb), original));
  HS_EXPECT_TRUE(WB::slots_equal(WB::requested_slots(sdb), invalid));

  invalid = original;
  invalid.warp_program.inner.kind = WB::WarpStageKind::LEGACY_STEREO_NOISE;
  WB::request_slots(sdb, invalid);
  HS_EXPECT_TRUE(WB::slots_equal(WB::active_slots(sdb), original));
  HS_EXPECT_TRUE(WB::slots_equal(WB::requested_slots(sdb), invalid));

  WB::RequestedConfig invalid_params = WB::active_config(sdb);
  invalid_params.params.source.pattern_freq = 21.0f;
  const WB::RequestedConfig before_invalid_params = WB::active_config(sdb);
  WB::request_config(sdb, invalid_params);
  HS_EXPECT_TRUE(WB::active_config(sdb) == before_invalid_params);
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
  candidate.slots.value_transfer = static_cast<WB::ValueTransfer>(invalid_tag);
  HS_EXPECT_FALSE(WB::valid_config(candidate));
  candidate = WB::legacy_config();
  candidate.slots.coverage = static_cast<WB::CoveragePolicy>(invalid_tag);
  HS_EXPECT_FALSE(WB::valid_config(candidate));
  candidate = WB::legacy_config();
  candidate.slots.colorizer = static_cast<WB::Colorizer>(invalid_tag);
  HS_EXPECT_FALSE(WB::valid_config(candidate));

  const WB::RequestedConfig legacy_config = WB::legacy_config();
  HS_EXPECT_TRUE(WB::valid_config(legacy_config));
  HS_EXPECT_FALSE(
      WB::transition_admitted(WB::active_config(sdb), legacy_config));
  WB::request_config(sdb, legacy_config);
  HS_EXPECT_TRUE(WB::slots_equal(WB::active_slots(sdb), original));
  HS_EXPECT_FALSE(WB::transition_active(sdb));
  HS_EXPECT_TRUE(WB::active_config(sdb) == before_invalid_params);
  HS_EXPECT_TRUE(WB::requested_slots(sdb) == original);
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
    const Color4 at_start = WB::shade_transition(view, {valid, inactive, 0.0f});
    HS_EXPECT_TRUE(at_start.color == expected.color);
    HS_EXPECT_EQ(at_start.alpha, expected.alpha);
    const Color4 at_end = WB::shade_transition(view, {inactive, valid, 1.0f});
    HS_EXPECT_TRUE(at_end.color == expected.color);
    HS_EXPECT_EQ(at_end.alpha, expected.alpha);
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
    HS_EXPECT_EQ(WB::transition_frame(sdb).mix, 0.0f);
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
    const bool original_continuation = WB::transition_continues_choreo(sdb);
    WB::RequestedConfig queued = WB::legacy_config();
    queued.slots.function = WB::Function::RINGS;
    HS_EXPECT_TRUE(WB::valid_config(queued));
    HS_EXPECT_FALSE(WB::transition_admitted(captured_source, queued));
    WB::request_config(sdb, queued);
    HS_EXPECT_TRUE(WB::transition_from_config(sdb) == captured_source);
    HS_EXPECT_TRUE(WB::transition_to_config(sdb) == original_destination);
    HS_EXPECT_EQ(WB::transition_continues_choreo(sdb), original_continuation);
    HS_EXPECT_EQ(WB::transition_elapsed(sdb), elapsed_before_takeover);
    HS_EXPECT_TRUE(WB::requested_slots(sdb) == queued.slots);
    sdb.draw_frame();
    sdb.advance_display();
    HS_EXPECT_EQ(WB::transition_elapsed(sdb), elapsed_before_takeover + 1);

    const uint16_t elapsed_before_retarget = WB::transition_elapsed(sdb);
    const WB::RequestedConfig manual = WB::presets()[6];
    HS_EXPECT_TRUE(WB::transition_admitted(captured_source, manual));
    WB::request_config(sdb, manual);
    HS_EXPECT_TRUE(WB::transition_active(sdb));
    HS_EXPECT_TRUE(WB::transition_from_config(sdb) == captured_source);
    HS_EXPECT_TRUE(WB::transition_to_config(sdb) == manual);
    HS_EXPECT_FALSE(WB::transition_continues_choreo(sdb));
    HS_EXPECT_EQ(WB::transition_elapsed(sdb), elapsed_before_retarget);
    HS_EXPECT_TRUE(WB::requested_slots(sdb) == manual.slots);

    while (WB::transition_active(sdb) && WB::transition_frame(sdb).mix < 1.0f) {
      sdb.draw_frame();
      sdb.advance_display();
    }
    HS_EXPECT_TRUE(WB::transition_active(sdb));
    HS_EXPECT_EQ(WB::transition_frame(sdb).mix, 1.0f);
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
  test_shadierball_stable_preset_transition();
  test_shadierball_discrete_transition();
  test_shadierball_palette_resources();
  return fixture.result();
}

} // namespace shadierball_tests
} // namespace hs_test
