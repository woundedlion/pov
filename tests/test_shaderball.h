/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * Licensed under the PolyForm Noncommercial License 1.0.0
 */
#pragma once

#include <limits>
#include "pullback_manifest.generated.h"
#include "tests/test_effects.h"
#include "tests/test_fixture.h"
#include "tests/test_harness.h"

namespace hs_test {
namespace shaderball_tests {

using effects_tests::reset_effect_globals;
using effects_tests::SMALL_H;
using effects_tests::SMALL_W;

/** @brief White-box access to ShaderBall's typed pipeline. */
struct ShaderBallWhiteBox {
  using SB = ShaderBall<SMALL_W, SMALL_H>;
  using Function = SB::Function;
  using Projection = SB::Projection;
  using PeirceLayout = SB::PeirceLayout;
  using AiroceanLayout = SB::AiroceanLayout;
  using BonneHemisphere = SB::BonneHemisphere;
  using GnomonicHemispherePolicy = SB::GnomonicHemispherePolicy;
  using ProjectionFramePolicy = SB::ProjectionFramePolicy;
  using SurfaceLens = SB::SurfaceLens;
  using SurfaceNoise = SB::SurfaceNoise;
  using SurfaceNoisePlacement = SB::SurfaceNoisePlacement;
  using SurfaceCurlIntegrator = SB::SurfaceCurlIntegrator;
  using NoiseBasis = SB::NoiseBasis;
  using WarpEnvelope = SB::WarpEnvelope;
  using PolarMode = SB::PolarMode;
  using CurlIntegrator = SB::CurlIntegrator;
  using WarpStageKind = SB::WarpStageKind;
  using WarpStageSpec = SB::WarpStageSpec;
  using WarpStageParams = SB::WarpStageParams;
  using ProjectionParams = SB::ProjectionParams;
  using SurfaceNoiseParams = SB::SurfaceNoiseParams;
  using SignalWeight = SB::SignalWeight;
  using ValueTransfer = SB::ValueTransfer;
  using CoveragePolicy = SB::CoveragePolicy;
  using PaletteMode = SB::PaletteMode;
  using PaletteMapping = SB::PaletteMapping;
  using BrightnessEnvelope = SB::BrightnessEnvelope;
  using HueShiftMode = SB::HueShiftMode;
  using Slots = SB::Slots;
  using SourceParams = SB::SourceParams;
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
  using ConfigFieldId = SB::ConfigFieldId;
  using ConfigRestoreResult = SB::ConfigRestoreResult;
  using FullConfigSnapshot = SB::FullConfigSnapshot;
  using TopologyKey = SB::TopologyKey;
  using InversePipelineId = SB::InversePipelineId;
  using InverseStageKind = SB::InverseStageKind;
  using CodeEmission = SB::CodeEmission;
  using ApproximationOracleId = SB::ApproximationOracleId;
  using ProjectionStateProvider = SB::ProjectionStateProvider;
  using SurfaceStateProvider = SB::SurfaceStateProvider;
  using LensStateProvider = SB::LensStateProvider;
  using SourceStateProvider = SB::SourceStateProvider;
  using ValueStateProvider = SB::ValueStateProvider;
  using ColorStateProvider = SB::ColorStateProvider;
  template <bool Outer>
  using WarpStateProvider = SB::template WarpStateProvider<Outer>;

  static constexpr float AXIS_EPS = SB::GNOMONIC_AXIS_EPS;
  static constexpr uint32_t HUE_STEP = SB::HUE_STEP;

  static constexpr int NUM_FUNCTIONS = SB::NUM_FUNCTIONS;
  static constexpr int NUM_PROJECTIONS = SB::NUM_PROJECTIONS;
  static constexpr int NUM_PEIRCE_LAYOUTS = SB::NUM_PEIRCE_LAYOUTS;
  static constexpr int NUM_AIROCEAN_LAYOUTS = SB::NUM_AIROCEAN_LAYOUTS;
  static constexpr int NUM_BONNE_HEMISPHERES = SB::NUM_BONNE_HEMISPHERES;
  static constexpr int NUM_GNOMONIC_HEMISPHERES = SB::NUM_GNOMONIC_HEMISPHERES;
  static constexpr int NUM_PROJECTION_FRAMES = SB::NUM_PROJECTION_FRAMES;
  static constexpr int NUM_LENSES = SB::NUM_LENSES;
  static constexpr int NUM_SURFACE_NOISE = SB::NUM_SURFACE_NOISE;
  static constexpr int NUM_SURFACE_NOISE_PLACEMENTS =
      SB::NUM_SURFACE_NOISE_PLACEMENTS;
  static constexpr int NUM_SURFACE_CURL_INTEGRATORS =
      SB::NUM_SURFACE_CURL_INTEGRATORS;
  static constexpr int NUM_WARPS = SB::NUM_WARPS;
  static constexpr int NUM_NOISE_BASES = SB::NUM_NOISE_BASES;
  static constexpr int NUM_POLAR_MODES = SB::NUM_POLAR_MODES;
  static constexpr int NUM_CURL_INTEGRATORS = SB::NUM_CURL_INTEGRATORS;
  static constexpr int NUM_WARP_ENVELOPES = SB::NUM_WARP_ENVELOPES;
  static constexpr int NUM_SIGNALS = SB::NUM_SIGNALS;
  static constexpr int NUM_VALUE_TRANSFERS = SB::NUM_VALUE_TRANSFERS;
  static constexpr int NUM_COVERAGE_POLICIES = SB::NUM_COVERAGE_POLICIES;
  static constexpr int NUM_PALETTES = SB::NUM_PALETTES;
  static constexpr int NUM_PALETTE_MAPPINGS = SB::NUM_PALETTE_MAPPINGS;
  static constexpr int NUM_BRIGHTNESS_ENVELOPES = SB::NUM_BRIGHTNESS_ENVELOPES;
  static constexpr int NUM_HUE_SHIFT_MODES = SB::NUM_HUE_SHIFT_MODES;

  static ClockState clocks(const SB &sb) { return sb.runtime.clocks; }
  static void set_clocks(SB &sb, const ClockState &clocks) {
    sb.runtime.clocks = clocks;
  }
  static void seed_clocks(SB &sb, float value) {
    sb.runtime.clocks = {value, value, value, value, value};
  }
  static FrameState frame(const SB &sb) { return sb.prepare_frame(); }
  static FrameState preset_frame(SB &sb, size_t index) {
    const auto &preset = SB::PRESETS[index].config;
    HS_CHECK(sb.prepare_resource_union(preset, preset),
             "ShaderBall test preset resources must fit");
    return sb.prepare_frame(preset, sb.runtime);
  }
  /**
   * @brief Frame snapshot for an arbitrary config, bypassing GUI admission.
   * @details Kernel coverage sweeps combinations the authored roster does not
   * compile an inverse pipeline for, so they cannot be reached through
   * `request_config`; the reference shader and the individual stages still
   * serve them.
   */
  static FrameState config_frame(SB &sb, const RequestedConfig &config) {
    HS_CHECK(sb.prepare_resource_union(config, config),
             "ShaderBall test config resources must fit");
    return sb.prepare_frame(config, sb.runtime);
  }
  static Slots active_slots(const SB &sb) { return sb.active_slots; }
  static InversePipelineId active_pipeline(const SB &sb) {
    return sb.active_pipeline;
  }
  static RequestedConfig active_config(const SB &sb) {
    return {sb.active_slots, sb.blend.params};
  }
  static constexpr Slots generated_surface_noise_slots() {
    return SB::GENERATED_SURFACE_NOISE_SLOTS;
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
            PaletteMode::TRIADIC};
  }
  static constexpr RequestedConfig legacy_config() {
    return {legacy_slots(),
            {{6.0f, 0.05f, 0.0f, 0.0f, 0.8f, 0.006f},
             {{0.1f, 0.0f, 0.5f}, {0.1f, 0.0f, 0.5f}},
             {2.0f, 0.0f, 0.0f},
             {1.0f},
             {},
             {0.0f, 1.0f, 0.0f},
             {0.25f},
             {}}};
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
  static const RequestedConfig &display_config(const SB &sb) {
    return sb.display_config;
  }
  static const char *parameter_warning(const SB &sb, const char *name) {
    return sb.parameter_warning(name);
  }
  static void request_config(SB &sb, const RequestedConfig &config) {
    sb.requested_config = config;
    sb.requested_schema_bound = false;
    sb.apply_requested_config();
  }
  static constexpr size_t PARAM_CAPACITY = SB::PARAM_CAPACITY;
  /// Parameters @p config registers, independent of whether it is admissible.
  static size_t registered_parameter_count(SB &sb,
                                           const RequestedConfig &config) {
    sb.requested_config = config;
    sb.rebind_parameters();
    return sb.getParameters().size();
  }
  static bool slots_equal(const Slots &a, const Slots &b) { return a == b; }
  static constexpr bool valid_config(const RequestedConfig &config) {
    return SB::valid_config(config);
  }
  static bool admissible_config(const RequestedConfig &config) {
    return SB::admissible_config(config);
  }
  static constexpr bool seam_compatible(const RequestedConfig &config) {
    return SB::strict_seam_compatible(config);
  }
  static constexpr bool transition_admitted(const RequestedConfig &from,
                                            const RequestedConfig &to) {
    return SB::transition_admitted(from, to);
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
  static bool transition_active(const SB &sb) {
    return sb.state->transition.active;
  }
  static bool param_morph_active(const SB &sb) {
    return sb.state->param_morph.active;
  }
  static bool try_apply_config(SB &sb, const RequestedConfig &config,
                               uint16_t duration) {
    return sb.try_apply_config(config, duration, false, false);
  }
  static uint16_t param_morph_elapsed(const SB &sb) {
    return sb.state->param_morph.elapsed;
  }
  static const Params &live_params(const SB &sb) { return sb.blend.params; }
  static size_t preset_index(const SB &sb) { return sb.getPresetIndex(); }
  static float transition_mix(const SB &sb) {
    return SB::transition_mix(sb.state->transition.elapsed,
                              sb.state->transition.duration);
  }
  static const LookRuntime &transition_from_runtime(const SB &sb) {
    return sb.state->transition.from_runtime;
  }
  static const LookRuntime &transition_to_runtime(const SB &sb) {
    return sb.state->transition.to_runtime;
  }
  static const RequestedConfig &transition_from_config(const SB &sb) {
    return sb.state->transition.from_config;
  }
  static const RequestedConfig &transition_to_config(const SB &sb) {
    return sb.state->transition.to_config;
  }
  static bool transition_continues_choreo(const SB &sb) {
    return sb.state->transition.continue_choreo;
  }
  static uint16_t transition_elapsed(const SB &sb) {
    return sb.state->transition.elapsed;
  }
  static void set_transition_elapsed(SB &sb, uint16_t elapsed) {
    sb.state->transition.elapsed = elapsed;
  }
  static int32_t prepared_surface_noise_seed(const SB &sb) {
    for (size_t index = 0; index < sb.prepared_noise_count; ++index)
      if (sb.state->prepared_noise_keys[index].domain == NoiseDomain::SPHERE_3D)
        return sb.state->prepared_noise_keys[index].seed;
    return 0;
  }
  static void force_transition(SB &sb, const RequestedConfig &to,
                               uint16_t duration, bool continue_choreo) {
    const RequestedConfig from = active_config(sb);
    sb.state->param_morph.active = false;
    sb.state->transition = {from, to,       sb.runtime,      sb.runtime,
                            0,    duration, continue_choreo, true};
    const auto *from_program = SB::find_inverse_program(from);
    const auto *to_program = SB::find_inverse_program(to);
    sb.state->transition.from_pipeline =
        from_program == nullptr ? InversePipelineId::NONE : from_program->id;
    sb.state->transition.to_pipeline =
        to_program == nullptr ? InversePipelineId::NONE : to_program->id;
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
    if (phase.clear)
      return Color4();
    HS_CHECK(visible != nullptr,
             "through-clear visible phase requires an endpoint frame");
    const auto function = SB::resolve_shade_function(*visible);
    HS_CHECK(function != nullptr,
             "through-clear test topology has no compiled inverse pipeline");
    typename SB::FrameShader shader{visible, phase.alpha, function};
    return shader(view);
  }
  static void begin_blend(SB &sb) {
    sb.preset_dwell_armed = false;
    sb.begin_blend();
  }
  static void step_param_morph(SB &sb) {
    sb.prepare_param_morph();
    sb.advance_runtime(sb.runtime, {sb.active_slots, sb.blend.params},
                       {Quaternion(), Quaternion()});
    sb.finish_transitions();
  }
  static void finish_transition(SB &sb) { sb.finish_transitions(); }
  static void settle_transition(SB &sb) {
    for (int frame = 0; frame < 1024 && (sb.state->transition.active ||
                                         sb.state->param_morph.active);
         ++frame) {
      sb.draw_frame();
      sb.advance_display();
    }
  }
  static void refresh_display(SB &sb) { sb.refresh_parameter_display(); }
  static uint32_t walk_steps(const SB &sb) { return sb.walk_step_count; }
  static uint32_t generated_palette_steps(const SB &sb) {
    return sb.generated_palette_step_count;
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
             const FrameState &frame, bool inner = false,
             const Complex &source_period = Complex()) {
    const float phase =
        inner ? frame.clocks.warp_inner_phase : frame.clocks.warp_outer_phase;
    const FastNoiseLite *noise = inner ? frame.resources.inner_warp_noise
                                       : frame.resources.outer_warp_noise;
    return SB::warp_stage_lookup(
        input, projected, spec, params, phase, noise,
        SB::prepare_warp_stage(spec, params, phase, source_period),
        SB::tracks_displacement(frame));
  }
  static MaterialSample material(const ProjectedLookup &projected,
                                 const PlanarWarpResult &warped,
                                 const FrameState &frame) {
    const Complex source_coords =
        SB::condition_source_coords(warped.coords, frame);
    const float field = SB::sample_source(source_coords, projected, frame);
    return SB::shape_material(field, projected, warped, frame);
  }
  static float source(const ProjectedLookup &projected,
                      const FrameState &frame) {
    const Complex source_coords =
        SB::condition_source_coords(projected.coords, frame);
    return SB::sample_source(source_coords, projected, frame);
  }
  static MaterialSample shape(float field, const ProjectedLookup &projected,
                              const PlanarWarpResult &warped,
                              const FrameState &frame) {
    return SB::shape_material(field, projected, warped, frame);
  }
  static Color4 colorize(const MaterialSample &sample,
                         const FrameState &frame) {
    return SB::colorize(sample, frame);
  }
  static float palette_mapping(float value, PaletteMapping mapping,
                               float frequency, float phase = 0.0f) {
    return SB::palette_mapping_coordinate(value, mapping, frequency, phase);
  }
  static float brightness_envelope(float value, BrightnessEnvelope envelope,
                                   float depth) {
    return SB::brightness_envelope_gain(value, envelope, depth);
  }
  static Pixel prepared_hue_rotation(const FrameState &frame, float value,
                                     float amount) {
    return SB::sample_hue_rotation_lut(frame.prepared_hue_rotation, value,
                                       amount);
  }
  static Pixel direct_hue_rotation(const FrameState &frame, float value,
                                   float amount) {
    const Color4 base = frame.resources.generated_palette->get(value);
    return ::hue_rotate_lut_gamut(base, amount).color;
  }
  static float prepared_hue_noise(const FrameState &frame, const Vector &v) {
    return SB::sample_hue_noise_lut(frame.prepared_hue_noise, v);
  }
  static float direct_hue_noise(const FrameState &frame, const Vector &v) {
    const Vector q = noise_sphere_coordinate(
        v, frame.params.color.hue_noise_scale, frame.clocks.hue_noise_phase);
    return frame.resources.color_noise->GetNoiseSingle(q.x, q.y, q.z);
  }
  static Pixel prepared_hue_noise_color(const FrameState &frame,
                                        const Vector &v, float value) {
    const float palette_value = SB::palette_mapping_coordinate(
        value, frame.slots.palette_mapping,
        frame.params.color.mapping_frequency, frame.params.color.mapping_phase);
    const float amount =
        frame.params.color.hue_shift_amount * prepared_hue_noise(frame, v);
    return SB::sample_hue_rotation_lut(frame.prepared_hue_rotation,
                                       palette_value, amount);
  }
  static Pixel direct_hue_noise_color(const FrameState &frame, const Vector &v,
                                      float value) {
    const float palette_value = SB::palette_mapping_coordinate(
        value, frame.slots.palette_mapping,
        frame.params.color.mapping_frequency, frame.params.color.mapping_phase);
    const HueRotateBase base = make_hue_rotate_base(
        frame.resources.generated_palette->get(palette_value));
    const float amount =
        frame.params.color.hue_shift_amount * direct_hue_noise(frame, v);
    return hue_rotate_lut_gamut(base, amount).color;
  }
  static Color4 shade(const Vector &v, const FrameState &frame) {
    return SB::shade_dynamic(v, frame);
  }
  static Color4 pipeline_shade(const Vector &v, const FrameState &frame) {
    const auto function = SB::resolve_shade_function(frame);
    HS_CHECK(function != nullptr,
             "ShaderBall test topology has no compiled inverse pipeline");
    return function(v, frame);
  }
  static TopologyKey topology_key(const RequestedConfig &config) {
    return SB::make_topology_key(config);
  }
  template <WarpStageKind Outer, WarpStageKind Inner>
  static constexpr bool planar_warp_implements(const TopologyKey &key) {
    return SB::template SelectedPlanarWarpStage<Outer, Inner>::implements(key);
  }
  static size_t inverse_program_count() {
    return SB::inverse_programs().size();
  }
  static std::array<uint8_t, 28> topology_values(const TopologyKey &key) {
    const auto &[function, projection, projection_frame, surface_lens,
                 signal_weight, value_transfer, coverage, peirce_layout,
                 airocean_layout, bonne_hemisphere, gnomonic_hemisphere,
                 surface_noise, surface_noise_placement, surface_noise_basis,
                 surface_curl_integrator, source_noise_basis, outer_warp,
                 outer_warp_basis, outer_warp_envelope, outer_polar_mode,
                 outer_curl_integrator, outer_polar_harmonic, inner_warp,
                 inner_warp_basis, inner_warp_envelope, inner_polar_mode,
                 inner_curl_integrator, inner_polar_harmonic] = key;
    return {{
        static_cast<uint8_t>(function),
        static_cast<uint8_t>(projection),
        static_cast<uint8_t>(projection_frame),
        static_cast<uint8_t>(surface_lens),
        static_cast<uint8_t>(signal_weight),
        static_cast<uint8_t>(value_transfer),
        static_cast<uint8_t>(coverage),
        static_cast<uint8_t>(peirce_layout),
        static_cast<uint8_t>(airocean_layout),
        static_cast<uint8_t>(bonne_hemisphere),
        static_cast<uint8_t>(gnomonic_hemisphere),
        static_cast<uint8_t>(surface_noise),
        static_cast<uint8_t>(surface_noise_placement),
        static_cast<uint8_t>(surface_noise_basis),
        static_cast<uint8_t>(surface_curl_integrator),
        static_cast<uint8_t>(source_noise_basis),
        static_cast<uint8_t>(outer_warp),
        static_cast<uint8_t>(outer_warp_basis),
        static_cast<uint8_t>(outer_warp_envelope),
        static_cast<uint8_t>(outer_polar_mode),
        static_cast<uint8_t>(outer_curl_integrator),
        outer_polar_harmonic,
        static_cast<uint8_t>(inner_warp),
        static_cast<uint8_t>(inner_warp_basis),
        static_cast<uint8_t>(inner_warp_envelope),
        static_cast<uint8_t>(inner_polar_mode),
        static_cast<uint8_t>(inner_curl_integrator),
        inner_polar_harmonic,
    }};
  }
  static const TopologyKey &inverse_program_key(size_t index) {
    return SB::inverse_programs()[index].key;
  }
  static InversePipelineId inverse_program_id(size_t index) {
    return SB::inverse_programs()[index].id;
  }
  static const char *inverse_program_name(size_t index) {
    return SB::pipeline_name(inverse_program_id(index));
  }
  static InversePipelineId preset_program_id(size_t index) {
    return SB::PRESETS[index].pipeline;
  }
  static InversePipelineId inverse_program_id(const FrameState &frame) {
    const auto *program = SB::resolve_inverse_program(frame);
    return program == nullptr ? InversePipelineId::NONE : program->id;
  }
  static bool has_inverse_program(const RequestedConfig &config) {
    return SB::find_inverse_program(config) != nullptr;
  }
  static bool inverse_programs_well_formed() {
    const auto &programs = SB::inverse_programs();
    for (size_t index = 0; index < programs.size(); ++index) {
      if (programs[index].shade == nullptr)
        return false;
      for (size_t prior = 0; prior < index; ++prior)
        if (programs[index].id == programs[prior].id ||
            programs[index].key == programs[prior].key ||
            programs[index].shade == programs[prior].shade)
          return false;
    }
    return true;
  }
  static constexpr bool inverse_stage_contracts() {
    using ArityGate = typename SB::template InversePipelineValidation<
        false, typename SB::OuterCameraStage>;
    using Surface =
        typename SB::PeirceDodecahedralGridPipeline::Validation::SurfaceStage;
    using Color =
        typename SB::PeirceDodecahedralGridPipeline::Validation::ColorStage;
    return !ArityGate::ORDER && !ArityGate::CARRIERS &&
           SB::PeirceDodecahedralGridPipeline::Validation::ORDER &&
           SB::PeirceDodecahedralGridPipeline::Validation::RUN_RETURNS &&
           SB::PeirceDodecahedralGridPipeline::Validation::CARRIERS &&
           SB::PeirceDodecahedralGridPipeline::Validation::TERMINALS &&
           SB::PeirceDodecahedralGridPipeline::Validation::APPROXIMATIONS &&
           Surface::KIND == InverseStageKind::SURFACE_PROJECT &&
           Surface::EMISSION == CodeEmission::INLINE_ONLY &&
           Surface::ORACLE == ApproximationOracleId::PEIRCE_FAST_SQUARE &&
           Surface::NON_FLOATING_FIELDS_EXACT &&
           Surface::METRICS[0].limit == 1.2e-3f &&
           Surface::METRICS[1].limit == 2e-4f &&
           Surface::METRICS[2].limit == 128.0f &&
           Color::ORACLE == ApproximationOracleId::HUE_ROTATION_AND_NOISE_LUTS;
  }
  static float peirce_metric_limit(size_t index) {
    return SB::template SelectedSurfaceProjectStage<
               Projection::PEIRCE_QUINCUNCIAL,
               SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL>::METRICS[index]
        .limit;
  }
  static float color_metric_limit(size_t index) {
    return SB::ColorStage::METRICS[index].limit;
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
    return SB::curl_vector(p, noise, basis, scale, time);
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
  static Vector surface_noise(const Vector &v, const FrameState &frame) {
    return SB::apply_surface_noise(v, frame);
  }
  static Vector surface_curl_field(const Vector &v, const FrameState &frame) {
    return SB::surface_curl_field(v, frame);
  }
  static Vector dodecahedral_reference(const Vector &v) {
    return lenses::polyhedral_kaleidoscope_lens(v,
                                                lenses::DODECAHEDRAL_MIRRORS);
  }
  static float sample_function(Function function, const Complex &p,
                               const SourceState &source) {
    return SB::sample_function(function, p, source);
  }
  static float grid(const Complex &p, const SourceParams &params,
                    const SourceState &source) {
    return SB::grid(p, params, source);
  }
  static constexpr auto presets() {
    std::array<RequestedConfig, SB::PRESETS.size()> configs{};
    for (size_t index = 0; index < configs.size(); ++index)
      configs[index] = SB::PRESETS[index].config;
    return configs;
  }
  static constexpr InversePipelineId preset_pipeline(size_t index) {
    return SB::PRESETS[index].pipeline;
  }
  static constexpr NoiseChannelLayout
  warp_channel_layout(const WarpStageSpec &spec) {
    return SB::warp_resource_key(spec).channel_layout;
  }
  static const auto &choreo() { return SB::CHOREO; }
  static void make_triadic(uint32_t &hue, uint32_t sequence,
                           GenerativePalette &out) {
    SB::next_generated_palette(hue, sequence, PaletteHarmony::TRIADIC, 0.62f,
                               out);
  }
  static void make_complementary(uint32_t &hue, uint32_t sequence,
                                 GenerativePalette &out) {
    SB::next_generated_palette(hue, sequence, PaletteHarmony::COMPLEMENTARY,
                               0.62f, out);
  }
  static void make_analogous(uint32_t &hue, uint32_t sequence,
                             GenerativePalette &out) {
    SB::next_generated_palette(hue, sequence, PaletteHarmony::ANALOGOUS, 0.62f,
                               out);
  }
  static Pixel palette_color(const SB &sb, PaletteMode mode, float value) {
    return sb.palette_for(mode).get(value).color;
  }
  static float palette_alpha(const SB &sb, PaletteMode mode, float value) {
    return sb.palette_for(mode).get(value).alpha;
  }
  static Pixel generated_color(const SB &sb, float value) {
    return palette_color(sb, PaletteMode::TRIADIC, value);
  }
  static void update_palette_chroma(SB &sb, float chroma) {
    sb.update_palette_chroma(chroma);
  }
  static float visible_palette_chroma(const SB &sb) {
    return sb.visible_palette_chroma();
  }
};

inline uint32_t shaderball_float_payload(float value) {
  uint32_t payload;
  std::memcpy(&payload, &value, sizeof(payload));
  return payload;
}

inline bool
shaderball_snapshots_equal(const ShaderBallWhiteBox::FullConfigSnapshot &a,
                           const ShaderBallWhiteBox::FullConfigSnapshot &b) {
  return a.schema_version == b.schema_version && a.accepted == b.accepted &&
         a.requested == b.requested && a.pending == b.pending &&
         a.has_runtime == b.has_runtime && a.runtime == b.runtime;
}

/** @brief Full snapshots preserve every field and reject invalid input atomically. */
inline void test_shaderball_full_config_snapshot() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();

  WB::FullConfigSnapshot snapshot = sb.capture_full_config_snapshot();
  const size_t source_seed =
      static_cast<size_t>(WB::ConfigFieldId::SOURCE_NOISE_SEED);
  const size_t palette = static_cast<size_t>(WB::ConfigFieldId::SLOTS_PALETTE);
  const size_t inactive_phase =
      static_cast<size_t>(WB::ConfigFieldId::WARP_INNER_RADIAL_PHASE);
  const size_t hue_noise_speed =
      static_cast<size_t>(WB::ConfigFieldId::COLOR_HUE_NOISE_SPEED);
  const size_t palette_mapping =
      static_cast<size_t>(WB::ConfigFieldId::SLOTS_PALETTE_MAPPING);
  const size_t mapping_frequency =
      static_cast<size_t>(WB::ConfigFieldId::COLOR_MAPPING_FREQUENCY);

  snapshot.accepted[source_seed] = 0x80000000u;
  snapshot.requested[source_seed] = 0x80000000u;
  snapshot.accepted[inactive_phase] = shaderball_float_payload(5.75f);
  snapshot.requested[inactive_phase] = shaderball_float_payload(5.75f);
  snapshot.accepted[palette_mapping] =
      static_cast<uint32_t>(WB::PaletteMapping::CUP);
  snapshot.requested[palette_mapping] =
      static_cast<uint32_t>(WB::PaletteMapping::CUP);
  snapshot.accepted[mapping_frequency] = shaderball_float_payload(2.5f);
  snapshot.requested[mapping_frequency] = shaderball_float_payload(2.5f);
  snapshot.runtime = {0.25f, 0.5f, 0.75f, 0.0f, 1.25f, 1.5f,
                      1.75f, 0.0f, 2.25f, 2.5f, 2.75f, 3.0f};
  HS_EXPECT_EQ(sb.restore_full_config_snapshot(snapshot),
               WB::ConfigRestoreResult::APPLIED);
  HS_EXPECT_TRUE(
      shaderball_snapshots_equal(sb.capture_full_config_snapshot(), snapshot));

  snapshot.accepted[source_seed] = 0x7fffffffu;
  snapshot.requested[source_seed] = 0x7fffffffu;
  HS_EXPECT_EQ(sb.restore_full_config_snapshot(snapshot),
               WB::ConfigRestoreResult::APPLIED);
  HS_EXPECT_EQ(sb.capture_full_config_snapshot().accepted[source_seed],
               0x7fffffffu);

  snapshot.requested[palette] =
      static_cast<uint32_t>(WB::PaletteMode::COMPLEMENTARY);
  snapshot.pending[palette] = 1;
  HS_EXPECT_EQ(sb.restore_full_config_snapshot(snapshot),
               WB::ConfigRestoreResult::APPLIED);
  HS_EXPECT_TRUE(
      shaderball_snapshots_equal(sb.capture_full_config_snapshot(), snapshot));

  snapshot = sb.capture_full_config_snapshot();
  snapshot.accepted[hue_noise_speed] = shaderball_float_payload(0.01f);
  snapshot.requested[hue_noise_speed] = shaderball_float_payload(0.01f);
  HS_EXPECT_EQ(sb.restore_full_config_snapshot(snapshot),
               WB::ConfigRestoreResult::APPLIED);
  const WB::FullConfigSnapshot speed_clamped =
      sb.capture_full_config_snapshot();
  HS_EXPECT_EQ(speed_clamped.accepted[hue_noise_speed],
               shaderball_float_payload(0.001f));
  HS_EXPECT_EQ(speed_clamped.requested[hue_noise_speed],
               shaderball_float_payload(0.001f));

  const WB::FullConfigSnapshot before_failure =
      sb.capture_full_config_snapshot();
  WB::FullConfigSnapshot invalid = before_failure;
  invalid.schema_version = WB::SB::CONFIG_SCHEMA_VERSION + 1;
  HS_EXPECT_EQ(sb.restore_full_config_snapshot(invalid),
               WB::ConfigRestoreResult::UNSUPPORTED_VERSION);
  HS_EXPECT_TRUE(shaderball_snapshots_equal(sb.capture_full_config_snapshot(),
                                            before_failure));

  invalid = before_failure;
  invalid.accepted[static_cast<size_t>(WB::ConfigFieldId::SLOTS_FUNCTION)] =
      0xffffffffu;
  HS_EXPECT_EQ(sb.restore_full_config_snapshot(invalid),
               WB::ConfigRestoreResult::INVALID_VALUE);
  HS_EXPECT_TRUE(shaderball_snapshots_equal(sb.capture_full_config_snapshot(),
                                            before_failure));
}

/** @brief A mode edit clamps stale subordinate values to its new range. */
inline void test_shaderball_surface_noise_range_rebind() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();
  sb.setAnimationsPaused(true);
  const size_t noise_slot =
      static_cast<size_t>(WB::ConfigFieldId::SLOTS_SURFACE_NOISE);
  const size_t strength =
      static_cast<size_t>(WB::ConfigFieldId::SURFACE_NOISE_STRENGTH);
  const WB::FullConfigSnapshot settled = sb.capture_full_config_snapshot();

  WB::FullConfigSnapshot pending_curl = settled;
  pending_curl.accepted[noise_slot] =
      static_cast<uint32_t>(WB::SurfaceNoise::CURL);
  pending_curl.requested[noise_slot] =
      static_cast<uint32_t>(WB::SurfaceNoise::DIRECT);
  pending_curl.accepted[strength] = shaderball_float_payload(-0.25f);
  pending_curl.requested[strength] = shaderball_float_payload(-0.25f);
  pending_curl.pending[noise_slot] = 1;
  HS_EXPECT_EQ(sb.restore_full_config_snapshot(pending_curl),
               WB::ConfigRestoreResult::APPLIED);
  const WB::FullConfigSnapshot clamped = sb.capture_full_config_snapshot();
  HS_EXPECT_EQ(clamped.accepted[noise_slot],
               static_cast<uint32_t>(WB::SurfaceNoise::DIRECT));
  HS_EXPECT_EQ(clamped.requested[noise_slot],
               static_cast<uint32_t>(WB::SurfaceNoise::DIRECT));
  HS_EXPECT_EQ(clamped.accepted[strength], shaderball_float_payload(0.0f));
  HS_EXPECT_EQ(clamped.requested[strength], shaderball_float_payload(0.0f));
  HS_EXPECT_EQ(clamped.pending[noise_slot], uint8_t{0});
  sb.draw_frame();
  sb.advance_display();
  WB::refresh_display(sb);
  HS_EXPECT_EQ(WB::display_config(sb).slots.surface_noise,
               WB::SurfaceNoise::DIRECT);
  HS_EXPECT_NEAR(WB::display_config(sb).params.surface_noise.strength, 0.0f,
                 1e-6f);
  HS_EXPECT_TRUE(WB::parameter_warning(sb, "Surface Noise Strength") ==
                 nullptr);

  HS_EXPECT_EQ(sb.restore_full_config_snapshot(settled),
               WB::ConfigRestoreResult::APPLIED);
  HS_EXPECT_TRUE(
      shaderball_snapshots_equal(sb.capture_full_config_snapshot(), settled));
  HS_EXPECT_TRUE(WB::parameter_warning(sb, "Surface Noise Strength") ==
                 nullptr);
}

/** @brief Incompatible snapshot layouts are rejected atomically. */
inline void test_shaderball_incompatible_config_snapshot() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();

  const WB::FullConfigSnapshot current = sb.capture_full_config_snapshot();
  WB::FullConfigSnapshot incompatible = current;
  incompatible.schema_version = WB::SB::CONFIG_SCHEMA_VERSION - 1;
  HS_EXPECT_EQ(sb.restore_full_config_snapshot(incompatible),
               WB::ConfigRestoreResult::UNSUPPORTED_VERSION);
  HS_EXPECT_TRUE(
      shaderball_snapshots_equal(sb.capture_full_config_snapshot(), current));
}

/** @brief Dense live selectors preserve sparse snapshot storage IDs. */
inline void test_shaderball_selector_storage() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();
  const size_t lens =
      static_cast<size_t>(WB::ConfigFieldId::SLOTS_SURFACE_LENS);
  const size_t outer_warp =
      static_cast<size_t>(WB::ConfigFieldId::SLOTS_WARP_OUTER_KIND);
  const size_t pattern_freq =
      static_cast<size_t>(WB::ConfigFieldId::SOURCE_PATTERN_FREQ);
  for (uint32_t storage_id :
       {0u, 1u, 2u, 3u, 4u, 6u, 7u, 8u, 9u, 10u, 11u, 12u, 13u}) {
    WB::FullConfigSnapshot snapshot = sb.capture_full_config_snapshot();
    snapshot.accepted[lens] = storage_id;
    snapshot.requested[lens] = storage_id;
    HS_EXPECT_EQ(sb.restore_full_config_snapshot(snapshot),
                 WB::ConfigRestoreResult::APPLIED);
    HS_EXPECT_EQ(sb.capture_full_config_snapshot().accepted[lens], storage_id);
  }
  for (uint32_t storage_id : {0u, 2u, 3u, 4u, 5u, 6u, 7u, 8u}) {
    WB::FullConfigSnapshot snapshot = sb.capture_full_config_snapshot();
    snapshot.accepted[outer_warp] = storage_id;
    snapshot.requested[outer_warp] = storage_id;
    snapshot.accepted[pattern_freq] = shaderball_float_payload(1.0f);
    snapshot.requested[pattern_freq] = shaderball_float_payload(1.0f);
    const auto set_outer = [&](WB::ConfigFieldId id, float value) {
      const size_t index = static_cast<size_t>(id);
      snapshot.accepted[index] = shaderball_float_payload(value);
      snapshot.requested[index] = shaderball_float_payload(value);
    };
    set_outer(WB::ConfigFieldId::WARP_OUTER_SCALE, 1.0f);
    set_outer(WB::ConfigFieldId::WARP_OUTER_STRENGTH, 0.0f);
    set_outer(WB::ConfigFieldId::WARP_OUTER_SPEED, 0.0f);
    set_outer(WB::ConfigFieldId::WARP_OUTER_SCALE_X, 1.0f);
    set_outer(WB::ConfigFieldId::WARP_OUTER_SCALE_Y, 1.0f);
    set_outer(WB::ConfigFieldId::WARP_OUTER_FREQUENCY, 1.0f);
    set_outer(WB::ConfigFieldId::WARP_OUTER_RADIUS, 1.0f);
    set_outer(WB::ConfigFieldId::WARP_OUTER_CELL_X, 1.0f);
    set_outer(WB::ConfigFieldId::WARP_OUTER_CELL_Y, 1.0f);
    set_outer(WB::ConfigFieldId::WARP_OUTER_RADIAL_SCALE, 1.0f);
    set_outer(WB::ConfigFieldId::WARP_OUTER_EDGE_WIDTH, 0.1f);
    HS_EXPECT_EQ(sb.restore_full_config_snapshot(snapshot),
                 WB::ConfigRestoreResult::APPLIED);
    HS_EXPECT_EQ(sb.capture_full_config_snapshot().accepted[outer_warp],
                 storage_id);
  }
  for (const WB::RequestedConfig &preset : WB::presets()) {
    HS_EXPECT_NE(preset.slots.surface_lens, WB::SurfaceLens::TANGENT_NOISE);
    HS_EXPECT_NE(preset.slots.warp_program.outer.kind,
                 WB::WarpStageKind::LEGACY_STEREO_NOISE);
    HS_EXPECT_NE(preset.slots.warp_program.inner.kind,
                 WB::WarpStageKind::LEGACY_STEREO_NOISE);
    HS_EXPECT_TRUE(WB::valid_config(preset));
  }
}

/** @brief Every named clock wraps in its native domain. */
inline void test_shaderball_clocks_wrapped() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  hs::set_mock_time(0, 0);
  WB::SB sb;
  sb.init();
  sb.setAnimationsPaused(true);
  WB::seed_clocks(sb, TWO_PI_F * 4.0f);
  for (int frame = 0; frame < 32; ++frame) {
    hs::set_mock_time(frame * FRAME_MS, frame * FRAME_US);
    sb.draw_frame();
    sb.advance_display();
    const WB::ClockState clocks = WB::clocks(sb);
    for (float phase :
         {clocks.source_primary, clocks.source_secondary, clocks.source_angle,
          clocks.projection_spin, clocks.hue_noise_phase,
          clocks.source_noise_time, clocks.surface_noise_time,
          clocks.warp_outer_phase, clocks.warp_inner_phase}) {
      HS_EXPECT_GE(phase, 0.0f);
      HS_EXPECT_LT(phase, TWO_PI_F);
    }
  }
  hs::clear_mock_time();
}

/** @brief Pause gates preset selection while all live motion keeps advancing. */
inline void test_shaderball_pause_semantics() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();
  WB::RequestedConfig ambient = WB::presets()[6];
  ambient.params.source.speed = 0.019f;
  ambient.params.source.secondary_rate = 0.37f;
  ambient.params.source.angle_rate = 0.007f;
  ambient.params.source.noise_time_rate = 0.0004f;
  ambient.params.surface_noise.rate = 0.001f;
  ambient.params.warp.outer.speed = 0.01f;
  ambient.params.warp.inner.speed = 0.001f;
  ambient.params.projection.spin_rate = 0.009f;
  ambient.params.projection.wander = 1.0f;
  ambient.params.outer_camera.wander = 1.0f;
  ambient.params.color.hue_noise_speed = 0.001f;
  HS_EXPECT_TRUE(WB::valid_config(ambient));
  WB::request_config(sb, ambient);
  sb.setAnimationsPaused(true);

  const WB::ClockState paused_clocks = WB::clocks(sb);
  const Quaternion paused_projection_walk = WB::projection_walk(sb);
  const Quaternion paused_outer_walk = WB::outer_walk(sb);
  const uint32_t paused_walk_steps = WB::walk_steps(sb);
  const uint32_t paused_generated_steps = WB::generated_palette_steps(sb);
  const Pixel paused_generated_color = WB::generated_color(sb, 0.25f);
  const size_t paused_preset = WB::preset_index(sb);
  for (int frame = 0; frame < 120; ++frame) {
    sb.draw_frame();
    sb.advance_display();
  }
  HS_EXPECT_NE(WB::clocks(sb).source_primary, paused_clocks.source_primary);
  HS_EXPECT_NE(WB::clocks(sb).source_secondary, paused_clocks.source_secondary);
  HS_EXPECT_NE(WB::clocks(sb).source_angle, paused_clocks.source_angle);
  HS_EXPECT_NE(WB::clocks(sb).projection_spin, paused_clocks.projection_spin);
  HS_EXPECT_NE(WB::clocks(sb).hue_noise_phase, paused_clocks.hue_noise_phase);
  HS_EXPECT_NE(WB::clocks(sb).source_noise_time,
               paused_clocks.source_noise_time);
  HS_EXPECT_NE(WB::clocks(sb).surface_noise_time,
               paused_clocks.surface_noise_time);
  HS_EXPECT_NE(WB::clocks(sb).warp_outer_phase, paused_clocks.warp_outer_phase);
  HS_EXPECT_NE(WB::clocks(sb).warp_inner_phase, paused_clocks.warp_inner_phase);
  HS_EXPECT_TRUE(WB::projection_walk(sb) != paused_projection_walk);
  HS_EXPECT_TRUE(WB::outer_walk(sb) != paused_outer_walk);
  HS_EXPECT_EQ(WB::walk_steps(sb), paused_walk_steps + 120);
  HS_EXPECT_EQ(WB::generated_palette_steps(sb), paused_generated_steps + 120);
  const Pixel active_generated_color = WB::generated_color(sb, 0.25f);
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
  HS_EXPECT_EQ(WB::generated_palette_steps(sb), paused_generated_steps + 121);

  for (int frame = 0; frame < 120 && WB::preset_index(sb) == paused_preset;
       ++frame) {
    sb.draw_frame();
    sb.advance_display();
  }
  HS_EXPECT_NE(WB::preset_index(sb), paused_preset);
  HS_EXPECT_TRUE(WB::transition_active(sb) || WB::param_morph_active(sb));
  const size_t active_preset = WB::preset_index(sb);
  sb.setAnimationsPaused(true);
  for (int frame = 0; frame < 512; ++frame)
    WB::finish_transition(sb);
  HS_EXPECT_EQ(WB::preset_index(sb), active_preset);
  HS_EXPECT_FALSE(WB::param_morph_active(sb));
  HS_EXPECT_FALSE(WB::transition_active(sb));
  HS_EXPECT_TRUE(WB::active_config(sb) == WB::presets()[active_preset]);
  HS_EXPECT_TRUE(sb.animations_paused());
  sb.setAnimationsPaused(false);
  for (int frame = 0; frame < 120 && WB::preset_index(sb) == active_preset;
       ++frame) {
    sb.draw_frame();
    sb.advance_display();
  }
  HS_EXPECT_NE(WB::preset_index(sb), active_preset);
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

  const auto initial_coverage = WB::active_slots(sb).coverage;
  HS_EXPECT_TRUE(
      sb.updateParameter("Planar Warp 1",
                         static_cast<float>(WB::WarpStageKind::NONE)) ==
      ParamSetResult::APPLIED);
  sb.draw_frame();
  sb.advance_display();

  // Every lens below keeps the config valid, so each edit lands on the frame
  // after it is written.
  for (WB::SurfaceLens lens :
       {WB::SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL,
        WB::SurfaceLens::KALEIDOSCOPE, WB::SurfaceLens::GLITCH}) {
    HS_EXPECT_TRUE(sb.updateParameter("Lens", static_cast<float>(lens)) ==
                   ParamSetResult::APPLIED);
    HS_EXPECT_NE(WB::active_slots(sb).surface_lens, lens);
    sb.draw_frame();
    sb.advance_display();
    HS_EXPECT_FALSE(WB::transition_active(sb));
    HS_EXPECT_FALSE(WB::param_morph_active(sb));
    HS_EXPECT_EQ(WB::active_slots(sb).surface_lens, lens);
    HS_EXPECT_EQ(WB::active_slots(sb).coverage, initial_coverage);
    HS_EXPECT_TRUE(WB::valid_config(WB::active_config(sb)));
    HS_EXPECT_TRUE(WB::requested_config(sb) == WB::active_config(sb));
    HS_EXPECT_TRUE(WB::published_config(sb) == WB::active_config(sb));
  }
  HS_EXPECT_EQ(WB::active_slots(sb).surface_lens,
               WB::presets()[0].slots.surface_lens);
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

  frame.slots.warp_program.outer.kind = WB::WarpStageKind::NONE;
  frame.slots.warp_program.inner.kind = WB::WarpStageKind::NONE;
  const WB::MaterialSample material = WB::material(projected, warped, frame);
  HS_EXPECT_GE(material.value, 0.0f);
  HS_EXPECT_LE(material.value, 1.0f);
  HS_EXPECT_EQ(material.coverage,
               projected.value_weight * projected.value_weight);
  HS_EXPECT_EQ(material.path_length,
               projected.surface_path_length + warped.path_length);
  WB::PlanarWarpResult accumulated = warped;
  accumulated.path_length = 0.75f;
  HS_EXPECT_EQ(WB::shape(0.0f, projected, accumulated, frame).path_length,
               0.75f);
  frame.slots.coverage = WB::CoveragePolicy::PROJECTION_WEIGHT;
  HS_EXPECT_EQ(WB::material(projected, warped, frame).coverage,
               projected.value_weight);
  const Color4 color = WB::colorize(material, frame);
  HS_EXPECT_TRUE(color.alpha >= 0.0f);

  frame.slots.signal_weight = WB::SignalWeight::NONE;
  HS_EXPECT_EQ(WB::shape(3.0f, projected, warped, frame).value, 1.0f);
  HS_EXPECT_EQ(WB::shape(-3.0f, projected, warped, frame).value, 0.0f);

  WB::ProjectedLookup direct_meta{
      Complex(1.0f, 2.0f), 3, 4, WB::boundary_singular(), 0.25f, 0.1f, 0x10};
  WB::ProjectedLookup lensed_meta{
      Complex(3.0f, 4.0f), 7, 8, WB::boundary_singular(), 0.75f, 0.9f, 0x20};
  direct_meta.surface_path_length = 0.25f;
  lensed_meta.surface_path_length = 0.75f;
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
  HS_EXPECT_EQ(direct_endpoint.surface_path_length,
               direct_meta.surface_path_length);
  HS_EXPECT_EQ(lensed_endpoint.coords.re, lensed_meta.coords.re);
  HS_EXPECT_EQ(lensed_endpoint.coords.im, lensed_meta.coords.im);
  HS_EXPECT_EQ(lensed_endpoint.region_id, lensed_meta.region_id);
  HS_EXPECT_EQ(lensed_endpoint.component_id, lensed_meta.component_id);
  HS_EXPECT_EQ(lensed_endpoint.boundary_flags, lensed_meta.boundary_flags);
  HS_EXPECT_EQ(lensed_endpoint.fade_edge_distance,
               lensed_meta.fade_edge_distance);
  HS_EXPECT_EQ(lensed_endpoint.value_weight, lensed_meta.value_weight);
  HS_EXPECT_EQ(lensed_endpoint.flags, lensed_meta.flags);
  HS_EXPECT_EQ(lensed_endpoint.surface_path_length,
               lensed_meta.surface_path_length);
  HS_EXPECT_EQ(joined.coords.re, 2.5f);
  HS_EXPECT_EQ(joined.coords.im, 3.5f);
  HS_EXPECT_EQ(joined.region_id, uint8_t(7));
  HS_EXPECT_EQ(joined.component_id, uint8_t(8));
  HS_EXPECT_EQ(joined.boundary_flags, WB::boundary_singular());
  HS_EXPECT_EQ(joined.fade_edge_distance, 0.75f);
  HS_EXPECT_EQ(joined.flags, uint8_t(0x20));
  HS_EXPECT_EQ(joined.value_weight,
               pole_attenuation(2.5f * 2.5f + 3.5f * 3.5f, 2.0f));
  HS_EXPECT_EQ(joined.surface_path_length, 0.625f);
}

/** Unit-length allowance on a lens output. The lenses renormalize through a
 * fast reciprocal square root, and glitch/twist additionally displace before
 * it, so the result is unit to a few parts in a thousand rather than to ULP. */
constexpr float LENS_UNIT_TOL = 4e-3f;

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

    const Complex sinusoidal = WB::project_point(v, WB::Projection::SINUSOIDAL);
    const float radius = sqrtf(v.x * v.x + v.z * v.z);
    HS_EXPECT_NEAR(sinusoidal.re, std::fabs(fast_atan2(v.z, v.x)) * radius,
                   1e-6f);
    HS_EXPECT_EQ(sinusoidal.im, 0.5f * PI_F - fast_acos(v.y));

    const Complex gnomonic = WB::project_point(v, WB::Projection::GNOMONIC);
    HS_EXPECT_TRUE(std::isfinite(gnomonic.re));
    HS_EXPECT_TRUE(std::isfinite(gnomonic.im));
    for (WB::SurfaceLens lens :
         {WB::SurfaceLens::GLITCH, WB::SurfaceLens::TWIST,
          WB::SurfaceLens::KALEIDOSCOPE})
      HS_EXPECT_NEAR(WB::apply_lens(v, lens).length(), 1.0f, LENS_UNIT_TOL);
  }

  reset_effect_globals();
  WB::SB landmark_sb;
  landmark_sb.init();
  WB::FrameState landmark_frame = WB::frame(landmark_sb);
  landmark_frame.slots.projection = WB::Projection::STEREOGRAPHIC;
  landmark_frame.slots.projection_frame = WB::ProjectionFramePolicy::IDENTITY;
  landmark_frame.slots.surface_lens = WB::SurfaceLens::NONE;
  landmark_frame.slots.surface_noise = WB::SurfaceNoise::NONE;
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
  landmark_frame.slots.projection = WB::Projection::SINUSOIDAL;
  const WB::ProjectedLookup sinusoidal =
      WB::surface_project(Vector(1.0f, 0.0f, 0.0f), landmark_frame);
  HS_EXPECT_EQ(sinusoidal.boundary_flags, uint8_t(0));
  HS_EXPECT_EQ(sinusoidal.flags, WB::projection_folded());
  HS_EXPECT_EQ(
      WB::surface_project(Vector(0.0f, 0.0f, -1.0f), landmark_frame).region_id,
      uint8_t(1));
  landmark_frame.slots.projection = WB::Projection::GNOMONIC;
  HS_EXPECT_EQ(
      WB::surface_project(Vector(1.0f, 0.0f, 0.0f), landmark_frame)
          .boundary_flags,
      static_cast<uint8_t>(WB::boundary_cut() | WB::boundary_singular()));

  const Vector v(0.6f, 0.48f, 0.64f);
  const Complex lensed = stereo(lenses::glitch_lens(v));
  WB::FrameState frame = WB::frame(landmark_sb);
  frame.slots.projection = WB::Projection::STEREOGRAPHIC;
  frame.slots.projection_frame = WB::ProjectionFramePolicy::IDENTITY;
  frame.slots.surface_lens = WB::SurfaceLens::GLITCH;
  frame.slots.surface_noise = WB::SurfaceNoise::NONE;
  frame.transforms.projection_conj = Quaternion();
  frame.params.surface_lens.mix = 0.0f;
  const WB::ProjectedLookup start = WB::surface_project(v, frame);
  frame.params.surface_lens.mix = 1.0f;
  const WB::ProjectedLookup end = WB::surface_project(v, frame);
  HS_EXPECT_NEAR(start.coords.re, lensed.re, 1e-6f);
  HS_EXPECT_NEAR(start.coords.im, lensed.im, 1e-6f);
  // glitch_lens' polar terms are FMA-contractable, so the reference above and
  // the pipeline's own call can round 2 ULP apart under -O2; a wrong lens
  // branch would miss by ~1.
  HS_EXPECT_NEAR(end.coords.re, lensed.re, 1e-6f);
  HS_EXPECT_NEAR(end.coords.im, lensed.im, 1e-6f);
}

/** @brief The reflection-only six-sector fold matches the legacy polar map. */
inline void test_shaderball_kaleidoscope_reflection_fold() {
  using WB = ShaderBallWhiteBox;
  constexpr float SECTOR = TWO_PI_F / 6.0f;
  auto reference = [&](const Vector &v) {
    const float radius = sqrtf(v.x * v.x + v.z * v.z);
    float azimuth = fmodf(fast_atan2(v.z, v.x) + PI_F, SECTOR);
    if (azimuth > 0.5f * SECTOR)
      azimuth = SECTOR - azimuth;
    return Vector(radius * fast_cosf(azimuth), v.y,
                  radius * fast_sinf(azimuth));
  };

  float max_coordinate_error = 0.0f;
  float max_length_error = 0.0f;
  for (int latitude_step = -64; latitude_step <= 64; ++latitude_step) {
    const float latitude = latitude_step * (0.5f * PI_F / 64.0f);
    const float radius = cosf(latitude);
    const float y = sinf(latitude);
    for (int longitude_step = 0; longitude_step < 1440; ++longitude_step) {
      const float longitude = longitude_step * (TWO_PI_F / 1440.0f);
      const Vector input(radius * cosf(longitude), y, radius * sinf(longitude));
      const Vector expected = reference(input);
      const Vector actual =
          WB::apply_lens(input, WB::SurfaceLens::KALEIDOSCOPE);
      max_coordinate_error = std::max(
          max_coordinate_error,
          std::max(fabsf(actual.x - expected.x), fabsf(actual.z - expected.z)));
      max_length_error =
          std::max(max_length_error, fabsf(actual.length() - input.length()));
    }
  }
  HS_EXPECT_LT(max_coordinate_error, 6e-3f);
  HS_EXPECT_LT(max_length_error, 2e-6f);

  constexpr float BOUNDARY_EPSILON = 1e-5f;
  for (int boundary = -6; boundary <= 6; ++boundary) {
    const float angle = boundary * (PI_F / 6.0f);
    Vector folded[3];
    for (int side = -1; side <= 1; ++side) {
      const float sample_angle = angle + side * BOUNDARY_EPSILON;
      folded[side + 1] =
          WB::apply_lens(Vector(cosf(sample_angle), 0.0f, sinf(sample_angle)),
                         WB::SurfaceLens::KALEIDOSCOPE);
      HS_EXPECT_GE(folded[side + 1].z, -2e-6f);
      HS_EXPECT_GE(folded[side + 1].x,
                   1.7320508075688772f * folded[side + 1].z - 2e-6f);
    }
    HS_EXPECT_NEAR(folded[0].x, folded[2].x, 2e-5f);
    HS_EXPECT_NEAR(folded[0].z, folded[2].z, 2e-5f);
  }
}

/** @brief Polyhedral lenses fold every simple-mirror orbit into one chamber. */
inline void test_shaderball_polyhedral_kaleidoscopes() {
  using WB = ShaderBallWhiteBox;
  struct Symmetry {
    WB::SurfaceLens lens;
    std::array<Vector, 3> mirrors;
  };
  static constexpr Symmetry SYMMETRIES[] = {
      {WB::SurfaceLens::KALEIDOSCOPE_TETRAHEDRAL,
       {Vector(1.0f, 0.0f, 0.0f), Vector(-0.5f, 0.8660254038f, 0.0f),
        Vector(0.0f, -0.5773502692f, 0.8164965809f)}},
      {WB::SurfaceLens::KALEIDOSCOPE_OCTAHEDRAL,
       {Vector(1.0f, 0.0f, 0.0f), Vector(-0.7071067812f, 0.7071067812f, 0.0f),
        Vector(0.0f, -0.7071067812f, 0.7071067812f)}},
      {WB::SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL,
       {Vector(1.0f, 0.0f, 0.0f), Vector(-0.8090169944f, 0.3090169944f, -0.5f),
        Vector(0.0f, 0.0f, 1.0f)}},
      {WB::SurfaceLens::KALEIDOSCOPE_TRIANGULAR_PRISM,
       {Vector(0.0f, 1.0f, 0.0f), Vector(0.0f, 0.0f, 1.0f),
        Vector(0.8660254038f, 0.0f, -0.5f)}},
      {WB::SurfaceLens::KALEIDOSCOPE_SQUARE_PRISM,
       {Vector(0.0f, 1.0f, 0.0f), Vector(0.0f, 0.0f, 1.0f),
        Vector(0.7071067812f, 0.0f, -0.7071067812f)}},
      {WB::SurfaceLens::KALEIDOSCOPE_PENTAGONAL_PRISM,
       {Vector(0.0f, 1.0f, 0.0f), Vector(0.0f, 0.0f, 1.0f),
        Vector(0.5877852523f, 0.0f, -0.8090169944f)}},
      {WB::SurfaceLens::KALEIDOSCOPE_HEXAGONAL_PRISM,
       {Vector(0.0f, 1.0f, 0.0f), Vector(0.0f, 0.0f, 1.0f),
        Vector(0.5f, 0.0f, -0.8660254038f)}},
      {WB::SurfaceLens::KALEIDOSCOPE_OCTAGONAL_PRISM,
       {Vector(0.0f, 1.0f, 0.0f), Vector(0.0f, 0.0f, 1.0f),
        Vector(0.3826834324f, 0.0f, -0.9238795325f)}}};
  const Vector direction = Vector(-0.371f, 0.557f, -0.743f).normalized();
  Vector folded[std::size(SYMMETRIES)];

  for (size_t index = 0; index < std::size(SYMMETRIES); ++index) {
    const Symmetry &symmetry = SYMMETRIES[index];
    folded[index] = WB::apply_lens(direction, symmetry.lens);
    HS_EXPECT_NEAR(folded[index].length(), 1.0f, 1e-5f);
    const Vector idempotent = WB::apply_lens(folded[index], symmetry.lens);
    HS_EXPECT_NEAR(idempotent.x, folded[index].x, 1e-6f);
    HS_EXPECT_NEAR(idempotent.y, folded[index].y, 1e-6f);
    HS_EXPECT_NEAR(idempotent.z, folded[index].z, 1e-6f);

    for (const Vector &normal : symmetry.mirrors) {
      const float distance = dot(direction, normal);
      const Vector reflected(direction.x - 2.0f * distance * normal.x,
                             direction.y - 2.0f * distance * normal.y,
                             direction.z - 2.0f * distance * normal.z);
      const Vector equivalent = WB::apply_lens(reflected, symmetry.lens);
      HS_EXPECT_NEAR(equivalent.x, folded[index].x, 2e-5f);
      HS_EXPECT_NEAR(equivalent.y, folded[index].y, 2e-5f);
      HS_EXPECT_NEAR(equivalent.z, folded[index].z, 2e-5f);
    }
  }

  for (size_t index = 1; index < std::size(folded); ++index)
    HS_EXPECT_TRUE(folded[index - 1] != folded[index]);

  for (int latitude_step = -64; latitude_step <= 64; ++latitude_step) {
    const float latitude = latitude_step * (0.5f * PI_F / 64.0f);
    const float radius = cosf(latitude);
    const float y = sinf(latitude);
    for (int longitude_step = 0; longitude_step < 1440; ++longitude_step) {
      const float longitude = longitude_step * (TWO_PI_F / 1440.0f);
      const Vector input(radius * cosf(longitude), y, radius * sinf(longitude));
      const Vector expected = WB::dodecahedral_reference(input);
      const Vector actual =
          WB::apply_lens(input, WB::SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL);
      HS_EXPECT_EQ(actual.x, expected.x);
      HS_EXPECT_EQ(actual.y, expected.y);
      HS_EXPECT_EQ(actual.z, expected.z);
    }
  }

  reset_effect_globals();
  WB::SB sb;
  sb.init();
  const auto *lens = sb.getParameters().find("Lens");
  HS_EXPECT_TRUE(lens != nullptr);
  HS_EXPECT_EQ(lens->option_count, 13);
  HS_EXPECT_EQ(std::string_view(lens->options[3]),
               std::string_view("Kaleidoscope (Azimuthal 6-fold)"));
  HS_EXPECT_EQ(std::string_view(lens->options[5]),
               std::string_view("Kaleidoscope (Tetrahedral)"));
  HS_EXPECT_EQ(std::string_view(lens->options[6]),
               std::string_view("Kaleidoscope (Octahedral / Cubic)"));
  HS_EXPECT_EQ(std::string_view(lens->options[7]),
               std::string_view("Kaleidoscope (Dodecahedral / Icosahedral)"));
  HS_EXPECT_EQ(std::string_view(lens->options[8]),
               std::string_view("Kaleidoscope (Triangular Prism)"));
  HS_EXPECT_EQ(std::string_view(lens->options[9]),
               std::string_view("Kaleidoscope (Square Prism)"));
  HS_EXPECT_EQ(std::string_view(lens->options[10]),
               std::string_view("Kaleidoscope (Pentagonal Prism)"));
  HS_EXPECT_EQ(std::string_view(lens->options[11]),
               std::string_view("Kaleidoscope (Hexagonal Prism)"));
  HS_EXPECT_EQ(std::string_view(lens->options[12]),
               std::string_view("Kaleidoscope (Octagonal Prism)"));
}

/** Longitude roundtrip allowance. Equirectangular is deliberately folded and
 * tapered near the poles, so the recovered longitude carries the taper's
 * residual; latitude comes straight out of an asin and is far tighter. */
constexpr float EQUIRECT_LON_TOL = 5e-3f;
/** Latitude roundtrip allowance: one float asin of a normalized component. */
constexpr float EQUIRECT_LAT_TOL = 3e-4f;

/** @brief Equirectangular is unfolded, periodic, and cut at the antimeridian. */
inline void test_shaderball_equirectangular_projection() {
  using WB = ShaderBallWhiteBox;
  auto lon_lat = [](float longitude, float latitude) {
    const float cp = cosf(latitude);
    return Vector(cp * cosf(longitude), sinf(latitude), cp * sinf(longitude));
  };
  for (float longitude : {-2.5f, -0.75f, 0.0f, 0.75f, 2.5f})
    for (float latitude : {-1.2f, -0.3f, 0.0f, 0.3f, 1.2f}) {
      const Complex coords = WB::project_point(lon_lat(longitude, latitude),
                                               WB::Projection::EQUIRECTANGULAR);
      HS_EXPECT_NEAR(coords.re, longitude, EQUIRECT_LON_TOL);
      HS_EXPECT_NEAR(coords.im, latitude, EQUIRECT_LAT_TOL);
    }

  reset_effect_globals();
  WB::SB sb;
  sb.init();
  WB::FrameState frame = WB::frame(sb);
  frame.slots.projection = WB::Projection::EQUIRECTANGULAR;
  frame.slots.projection_frame = WB::ProjectionFramePolicy::IDENTITY;
  frame.slots.surface_lens = WB::SurfaceLens::NONE;
  frame.slots.surface_noise = WB::SurfaceNoise::NONE;
  frame.transforms.projection_conj = Quaternion();
  frame.params.projection.central_meridian = 0.0f;

  const WB::ProjectedLookup prime =
      WB::surface_project(Vector(1.0f, 0.0f, 0.0f), frame);
  HS_EXPECT_EQ(prime.flags, uint8_t(0));
  HS_EXPECT_EQ(prime.boundary_flags, WB::boundary_cut());
  HS_EXPECT_EQ(prime.region_id, uint8_t(0));
  HS_EXPECT_NEAR(prime.coords.re, 0.0f, 1e-5f);
  HS_EXPECT_NEAR(prime.fade_edge_distance, PI_F, 1e-5f);

  const WB::ProjectedLookup seam =
      WB::surface_project(Vector(-1.0f, 0.0f, 0.0f), frame);
  HS_EXPECT_NEAR(seam.fade_edge_distance, 0.0f, 1e-5f);

  const WB::ProjectedLookup east =
      WB::surface_project(Vector(0.0f, 0.0f, 1.0f), frame);
  const WB::ProjectedLookup west =
      WB::surface_project(Vector(0.0f, 0.0f, -1.0f), frame);
  HS_EXPECT_NEAR(east.coords.re, 0.5f * PI_F, 5e-3f);
  HS_EXPECT_NEAR(west.coords.re, -0.5f * PI_F, 5e-3f);
  HS_EXPECT_EQ(east.region_id, west.region_id);

  frame.params.projection.central_meridian = 0.5f * PI_F;
  const WB::ProjectedLookup recentred =
      WB::surface_project(Vector(0.0f, 0.0f, 1.0f), frame);
  HS_EXPECT_NEAR(recentred.coords.re, 0.0f, 5e-3f);
  frame.params.projection.central_meridian = 0.0f;

  WB::RequestedConfig config = WB::legacy_config();
  config.slots.projection = WB::Projection::EQUIRECTANGULAR;
  config.slots.warp_program.outer.kind = WB::WarpStageKind::NONE;
  config.slots.surface_lens = WB::SurfaceLens::NONE;
  config.params.surface_lens.mix = 0.0f;
  HS_EXPECT_TRUE(WB::valid_config(config));
  HS_EXPECT_FALSE(WB::has_inverse_program(config));
  HS_EXPECT_EQ(
      sb.updateParameter("Projection",
                         static_cast<float>(WB::Projection::EQUIRECTANGULAR)),
      ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(sb.getParameters().find("Pole Fade") != nullptr);
  HS_EXPECT_TRUE(sb.getParameters().find("Central Meridian") != nullptr);
  HS_EXPECT_TRUE(WB::parameter_warning(sb, "Projection") == nullptr);
  sb.draw_frame();
  sb.advance_display();
  WB::settle_transition(sb);
  HS_EXPECT_EQ(WB::active_slots(sb).projection,
               WB::Projection::EQUIRECTANGULAR);
  HS_EXPECT_EQ(WB::active_pipeline(sb), WB::InversePipelineId::NONE);
  const Color4 color = WB::shade(Vector(0.31f, 0.87f, -0.38f).normalized(),
                                 WB::config_frame(sb, config));
  HS_EXPECT_TRUE(std::isfinite(color.alpha));
  HS_EXPECT_GE(color.alpha, 0.0f);
  HS_EXPECT_LE(color.alpha, 1.0f);
}

/** @brief Both sides of every projection seam use the authored fade width. */
inline void test_shaderball_flush_edge_fade() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();
  WB::FrameState frame = WB::frame(sb);
  frame.slots.coverage = WB::CoveragePolicy::EDGE_FADE;
  frame.params.value.edge_width = 0.1f;
  const WB::PlanarWarpResult warped{};
  WB::ProjectedLookup under{Complex(), 0,    0, WB::boundary_cut(),
                            0.01f,     1.0f, 0};
  WB::ProjectedLookup over = under;
  const auto coverage = [&](const WB::ProjectedLookup &projected) {
    return WB::shape(0.0f, projected, warped, frame).coverage;
  };

  frame.slots.projection = WB::Projection::BONNE;
  over.region_id = 1;
  HS_EXPECT_NEAR(coverage(under), 0.028f, 1e-6f);
  HS_EXPECT_NEAR(coverage(over), 0.028f, 1e-6f);

  frame.slots.projection = WB::Projection::PEIRCE_QUINCUNCIAL;
  under.region_id = 2;
  over.region_id = 3;
  HS_EXPECT_NEAR(coverage(under), 0.028f, 1e-6f);
  HS_EXPECT_NEAR(coverage(over), 0.028f, 1e-6f);

  frame.slots.projection = WB::Projection::AIROCEAN;
  under.edge_class = 9;
  over.edge_class = 14;
  HS_EXPECT_NEAR(coverage(under), 0.028f, 1e-6f);
  HS_EXPECT_NEAR(coverage(over), 0.028f, 1e-6f);

  frame.params.value.edge_width = 0.0f;
  HS_EXPECT_EQ(coverage(under), 1.0f);
  under.fade_edge_distance = 0.0f;
  HS_EXPECT_EQ(coverage(under), 0.0f);
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
    }
  }
}

/** @brief Grid sampling reduces to its coupled and direct endpoints. */
inline void test_shaderball_coupled_source() {
  using WB = ShaderBallWhiteBox;
  const float values[] = {-6.0f, -2.5f, -0.7f, 0.0f, 0.9f, 3.1f, 5.8f};
  const float phases[] = {0.0f, 0.8f, 2.4f, 4.9f};
  for (float re : values) {
    for (float im : values) {
      const Complex p(re, im);
      for (float primary : phases) {
        for (float secondary : phases) {
          const WB::SourceState source{primary, secondary, 0.0f, 1.0f, 0.0f};
          for (float complexity : {0.5f, 1.7f, 3.0f}) {
            WB::SourceParams params;
            params.complexity = complexity;
            const float coupled =
                fast_sinf(re + complexity * fast_sinf(im + primary)) *
                fast_cosf(im + complexity * fast_cosf(re - secondary));
            HS_EXPECT_NEAR(WB::grid(p, params, source), coupled,
                           SOURCE_DRIFT_BOUND);
          }
          WB::SourceParams params;
          params.pattern_mix = 1.0f;
          const float direct =
              fast_sinf(re + primary) * fast_cosf(im - secondary);
          HS_EXPECT_EQ(WB::grid(p, params, source), direct);
        }
      }
    }
  }
}
/** @brief Presets retain the curated topology roster and generated palette path. */
inline void test_shaderball_preset_bank() {
  using WB = ShaderBallWhiteBox;
  const auto &presets = WB::presets();
  const auto &choreo = WB::choreo();
  HS_EXPECT_EQ(presets.size(), size_t(18));
  HS_EXPECT_EQ(choreo.size(), presets.size());
  for (const auto &entry : choreo) {
    HS_EXPECT_FALSE(entry.staggered);
    HS_EXPECT_EQ(entry.dwell_min, uint16_t(0));
    HS_EXPECT_EQ(entry.dwell_max, uint16_t(0));
    HS_EXPECT_EQ(entry.blend_frames, uint16_t(480));
  }

  bool has_hue_shift = false;
  for (size_t index = 0; index < presets.size(); ++index) {
    const auto &preset = presets[index];
    const WB::PaletteMode expected_palette = index == 10 || index >= 12
                                                 ? WB::PaletteMode::ANALOGOUS
                                                 : WB::PaletteMode::TRIADIC;
    HS_EXPECT_EQ(preset.slots.palette, expected_palette);
    HS_EXPECT_TRUE(WB::seam_compatible(preset));
    HS_EXPECT_TRUE(WB::valid_config(preset));
    HS_EXPECT_TRUE(WB::has_inverse_program(preset));
    HS_EXPECT_LE(fabsf(preset.params.color.hue_noise_speed), 0.001f);
    has_hue_shift |= preset.params.color.hue_shift_amount != 0.0f;
  }
  HS_EXPECT_TRUE(has_hue_shift);

  const auto &wave_shear = presets[0];
  HS_EXPECT_EQ(wave_shear.slots.warp_program.outer.kind,
               WB::WarpStageKind::WAVE_SHEAR);
  HS_EXPECT_EQ(wave_shear.params.warp.outer.speed, 1.0f / 64.0f);
  const auto &mirror = presets[1];
  HS_EXPECT_EQ(mirror.slots.warp_program.inner.kind,
               WB::WarpStageKind::MIRROR_TILE);
  HS_EXPECT_EQ(mirror.params.warp.inner.speed, 0.0f);
  const auto &affine_lattice = presets[7];
  HS_EXPECT_EQ(affine_lattice.slots.function, WB::Function::PRIMITIVE_LATTICE);
  HS_EXPECT_EQ(affine_lattice.slots.projection, WB::Projection::GNOMONIC);
  HS_EXPECT_EQ(affine_lattice.slots.warp_program.outer.kind,
               WB::WarpStageKind::AFFINE_FRAME);
  HS_EXPECT_EQ(affine_lattice.slots.value_transfer,
               WB::ValueTransfer::ISO_CONTOUR);
  HS_EXPECT_EQ(affine_lattice.slots.coverage,
               WB::CoveragePolicy::PROJECTION_WEIGHT);
  HS_EXPECT_EQ(affine_lattice.params.source.lattice_cell_scale, 1.22925f);
  HS_EXPECT_EQ(affine_lattice.params.source.lattice_softness, 0.1608203f);
  HS_EXPECT_EQ(affine_lattice.params.warp.outer.translation_x, 4.0f);
  HS_EXPECT_EQ(affine_lattice.params.warp.outer.translation_y, 4.0f);
  HS_EXPECT_EQ(affine_lattice.params.warp.outer.scale_x, 1.0f);
  HS_EXPECT_EQ(affine_lattice.params.warp.outer.scale_y, 1.0f);
  HS_EXPECT_EQ(affine_lattice.params.outer_camera.wander, 1.0f);
  const auto &fine_curl = presets[8];
  const auto &coarse_curl = presets[9];
  HS_EXPECT_EQ(fine_curl.slots.function, WB::Function::PRIMITIVE_LATTICE);
  HS_EXPECT_EQ(fine_curl.slots.projection, WB::Projection::SINUSOIDAL);
  HS_EXPECT_EQ(fine_curl.slots.surface_noise, WB::SurfaceNoise::CURL);
  HS_EXPECT_EQ(fine_curl.slots.surface_noise_placement,
               WB::SurfaceNoisePlacement::BEFORE_LENS);
  HS_EXPECT_TRUE(WB::slots_equal(fine_curl.slots, coarse_curl.slots));
  HS_EXPECT_EQ(fine_curl.params.surface_noise.scale, 1.78815627f);
  HS_EXPECT_EQ(coarse_curl.params.surface_noise.scale, 3.29720306f);
  const auto &polar_wave = presets[10];
  HS_EXPECT_EQ(polar_wave.slots.surface_lens,
               WB::SurfaceLens::KALEIDOSCOPE_TRIANGULAR_PRISM);
  HS_EXPECT_EQ(polar_wave.slots.warp_program.outer.kind,
               WB::WarpStageKind::POLAR_CHART);
  HS_EXPECT_EQ(polar_wave.slots.warp_program.inner.kind,
               WB::WarpStageKind::WAVE_SHEAR);
  HS_EXPECT_EQ(polar_wave.slots.palette, WB::PaletteMode::ANALOGOUS);
  const auto &vector_mirror = presets[11];
  HS_EXPECT_EQ(vector_mirror.slots.function, WB::Function::GRID);
  HS_EXPECT_EQ(vector_mirror.slots.projection, WB::Projection::GNOMONIC);
  HS_EXPECT_EQ(vector_mirror.slots.surface_lens,
               WB::SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL);
  HS_EXPECT_EQ(vector_mirror.slots.warp_program.outer.kind,
               WB::WarpStageKind::VECTOR_NOISE);
  HS_EXPECT_EQ(vector_mirror.slots.warp_program.inner.kind,
               WB::WarpStageKind::MIRROR_TILE);
  HS_EXPECT_EQ(WB::warp_channel_layout(vector_mirror.slots.warp_program.outer),
               NoiseChannelLayout::DIRECT_VECTOR_V2);
  HS_EXPECT_EQ(vector_mirror.params.source.pattern_freq, 4.9755f);
  HS_EXPECT_EQ(vector_mirror.params.warp.outer.strength, 0.138f);
  HS_EXPECT_EQ(vector_mirror.params.warp.outer.field_angle, 2.23053074f);
  HS_EXPECT_EQ(vector_mirror.params.warp.inner.speed, 0.00327999983f);
  HS_EXPECT_EQ(vector_mirror.params.color.hue_shift_amount, 0.721f);
  const auto &inner_mirror = presets[12];
  HS_EXPECT_EQ(inner_mirror.slots.function, WB::Function::GRID);
  HS_EXPECT_EQ(inner_mirror.slots.projection, WB::Projection::STEREOGRAPHIC);
  HS_EXPECT_EQ(inner_mirror.slots.projection_frame,
               WB::ProjectionFramePolicy::SPIN_WANDER);
  HS_EXPECT_EQ(inner_mirror.slots.surface_lens,
               WB::SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL);
  HS_EXPECT_EQ(inner_mirror.slots.warp_program.outer.kind,
               WB::WarpStageKind::NONE);
  HS_EXPECT_EQ(inner_mirror.slots.warp_program.inner.kind,
               WB::WarpStageKind::MIRROR_TILE);
  HS_EXPECT_EQ(inner_mirror.slots.coverage,
               WB::CoveragePolicy::PROJECTION_WEIGHT_SQUARED);
  HS_EXPECT_EQ(inner_mirror.slots.palette_mapping, WB::PaletteMapping::CUP);
  HS_EXPECT_EQ(inner_mirror.params.source.pattern_freq, 2.82629991f);
  HS_EXPECT_EQ(inner_mirror.params.source.complexity, 0.513f);
  HS_EXPECT_EQ(inner_mirror.params.source.secondary_rate, 0.8f);
  HS_EXPECT_EQ(inner_mirror.params.source.angle_rate, 0.0269999988f);
  const auto &hexagonal_mirror = presets[13];
  HS_EXPECT_EQ(hexagonal_mirror.slots.function, WB::Function::TWIN_WAVE);
  HS_EXPECT_EQ(hexagonal_mirror.slots.projection,
               WB::Projection::STEREOGRAPHIC);
  HS_EXPECT_EQ(hexagonal_mirror.slots.surface_lens,
               WB::SurfaceLens::KALEIDOSCOPE_HEXAGONAL_PRISM);
  HS_EXPECT_EQ(hexagonal_mirror.slots.warp_program.outer.kind,
               WB::WarpStageKind::NONE);
  HS_EXPECT_EQ(hexagonal_mirror.slots.warp_program.inner.kind,
               WB::WarpStageKind::MIRROR_TILE);
  HS_EXPECT_EQ(hexagonal_mirror.slots.palette, WB::PaletteMode::ANALOGOUS);
  HS_EXPECT_EQ(hexagonal_mirror.slots.palette_mapping,
               WB::PaletteMapping::BELL);
  HS_EXPECT_EQ(hexagonal_mirror.params.source.pattern_freq, 3.881f);
  HS_EXPECT_EQ(hexagonal_mirror.params.source.speed, 0.128598228f);
  HS_EXPECT_EQ(hexagonal_mirror.params.color.hue_shift_amount, 0.226f);
  HS_EXPECT_EQ(hexagonal_mirror.params.color.hue_noise_scale, 1.47215629f);
  HS_EXPECT_EQ(hexagonal_mirror.params.color.hue_noise_speed, 0.000138f);
  HS_EXPECT_EQ(hexagonal_mirror.params.color.palette_chroma, 1.0f);
  HS_EXPECT_EQ(hexagonal_mirror.params.color.mapping_frequency, 1.341f);
  HS_EXPECT_EQ(hexagonal_mirror.params.color.mapping_phase, -1.0f);
  const auto &complex_grid = presets[14];
  HS_EXPECT_EQ(complex_grid.slots.function, WB::Function::GRID);
  HS_EXPECT_EQ(complex_grid.slots.projection, WB::Projection::STEREOGRAPHIC);
  HS_EXPECT_EQ(complex_grid.slots.surface_lens,
               WB::SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL);
  HS_EXPECT_EQ(complex_grid.slots.warp_program.outer.kind,
               WB::WarpStageKind::NONE);
  HS_EXPECT_EQ(complex_grid.slots.warp_program.inner.kind,
               WB::WarpStageKind::MIRROR_TILE);
  HS_EXPECT_EQ(complex_grid.params.source.pattern_freq, 2.82629991f);
  HS_EXPECT_EQ(complex_grid.params.source.complexity, 3.0f);
  HS_EXPECT_EQ(complex_grid.params.source.pattern_mix, 1.0f);
  HS_EXPECT_EQ(complex_grid.params.warp.inner.speed, 0.00013f);
  HS_EXPECT_EQ(complex_grid.params.warp.inner.cell_y, 0.997703135f);
  HS_EXPECT_EQ(complex_grid.params.color.hue_shift_amount, 0.366f);
  HS_EXPECT_EQ(complex_grid.params.color.hue_noise_scale, 1.47215629f);
  HS_EXPECT_EQ(complex_grid.params.color.palette_chroma, 1.0f);
  const auto &double_mapping_grid = presets[15];
  HS_EXPECT_EQ(double_mapping_grid.slots.function, WB::Function::GRID);
  HS_EXPECT_EQ(double_mapping_grid.slots.projection,
               WB::Projection::STEREOGRAPHIC);
  HS_EXPECT_EQ(double_mapping_grid.slots.surface_lens,
               WB::SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL);
  HS_EXPECT_EQ(double_mapping_grid.slots.warp_program.inner.kind,
               WB::WarpStageKind::MIRROR_TILE);
  HS_EXPECT_EQ(double_mapping_grid.params.source.pattern_freq, 3.9407f);
  HS_EXPECT_EQ(double_mapping_grid.params.source.complexity, 3.0f);
  HS_EXPECT_EQ(double_mapping_grid.params.source.pattern_mix, 1.0f);
  HS_EXPECT_EQ(double_mapping_grid.params.projection.wander, 0.165f);
  HS_EXPECT_EQ(double_mapping_grid.params.color.mapping_frequency, 2.0f);
  HS_EXPECT_EQ(double_mapping_grid.params.color.palette_chroma, 1.0f);
  const auto &equirectangular_grid = presets[16];
  HS_EXPECT_EQ(equirectangular_grid.slots.function, WB::Function::GRID);
  HS_EXPECT_EQ(equirectangular_grid.slots.projection,
               WB::Projection::EQUIRECTANGULAR);
  HS_EXPECT_EQ(equirectangular_grid.slots.surface_lens,
               WB::SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL);
  HS_EXPECT_EQ(equirectangular_grid.slots.warp_program.inner.kind,
               WB::WarpStageKind::MIRROR_TILE);
  HS_EXPECT_EQ(equirectangular_grid.params.source.pattern_freq, 3.9407f);
  HS_EXPECT_EQ(equirectangular_grid.params.source.complexity, 3.0f);
  HS_EXPECT_EQ(equirectangular_grid.params.projection.pole_fade, 2.14f);
  HS_EXPECT_EQ(equirectangular_grid.params.projection.wander, 0.165f);
  HS_EXPECT_EQ(equirectangular_grid.params.color.mapping_frequency, 2.0f);
  HS_EXPECT_EQ(equirectangular_grid.params.color.palette_chroma, 1.0f);
  const auto &single_mapping_equirectangular_grid = presets[17];
  HS_EXPECT_EQ(single_mapping_equirectangular_grid.slots.function,
               WB::Function::GRID);
  HS_EXPECT_EQ(single_mapping_equirectangular_grid.slots.projection,
               WB::Projection::EQUIRECTANGULAR);
  HS_EXPECT_EQ(single_mapping_equirectangular_grid.slots.surface_lens,
               WB::SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL);
  HS_EXPECT_EQ(single_mapping_equirectangular_grid.params.source.pattern_freq,
               3.9407f);
  HS_EXPECT_EQ(single_mapping_equirectangular_grid.params.projection.pole_fade,
               2.14f);
  HS_EXPECT_EQ(single_mapping_equirectangular_grid.params.projection.wander,
               0.165f);
  HS_EXPECT_EQ(
      single_mapping_equirectangular_grid.params.color.mapping_frequency, 1.0f);
  HS_EXPECT_EQ(single_mapping_equirectangular_grid.params.color.palette_chroma,
               1.0f);
  HS_EXPECT_EQ(inner_mirror.params.warp.outer.scale, 0.1f);
  HS_EXPECT_EQ(inner_mirror.params.warp.outer.speed, 0.5f);
  HS_EXPECT_EQ(inner_mirror.params.warp.inner.scale, 0.1f);
  HS_EXPECT_EQ(inner_mirror.params.warp.inner.speed, 0.00013f);
  HS_EXPECT_EQ(inner_mirror.params.warp.inner.cell_y, 0.997703135f);
  HS_EXPECT_EQ(inner_mirror.params.projection.pole_fade, 3.432f);
  HS_EXPECT_EQ(inner_mirror.params.projection.wander, 1.0f);
  HS_EXPECT_EQ(inner_mirror.params.surface_lens.mix, 1.0f);
  HS_EXPECT_EQ(inner_mirror.params.color.hue_shift_amount, 0.366f);
  HS_EXPECT_EQ(inner_mirror.params.color.hue_noise_scale, 1.47215629f);
  HS_EXPECT_EQ(inner_mirror.params.color.palette_chroma, 1.0f);
  HS_EXPECT_EQ(inner_mirror.params.outer_camera.wander, 1.0f);

  reset_effect_globals();
  WB::SB sb;
  sb.init();
  HS_EXPECT_TRUE(WB::slots_equal(WB::active_slots(sb), presets[0].slots));
}
/** @brief A staggered morph walks each changed parameter group in its own
 *         time slice. */
inline void test_shaderball_staggered_param_morph() {
  using WB = ShaderBallWhiteBox;
  WB::Params from;
  from.color.palette_chroma = 0.5f;
  WB::Params to = from;
  to.source.speed = 0.5f;
  to.color.palette_chroma = 0.75f;

  auto staggered = [&](float t) {
    WB::Params result;
    result.lerp_staggered(from, to, t);
    return result;
  };

  const WB::Params quarter = staggered(0.25f);
  HS_EXPECT_GT(quarter.source.speed, from.source.speed);
  HS_EXPECT_LT(quarter.source.speed, to.source.speed);
  HS_EXPECT_EQ(quarter.color.palette_chroma, from.color.palette_chroma);

  const WB::Params half = staggered(0.5f);
  HS_EXPECT_EQ(half.source.speed, to.source.speed);
  HS_EXPECT_EQ(half.color.palette_chroma, from.color.palette_chroma);
  HS_EXPECT_EQ(half.projection.pole_fade, from.projection.pole_fade);

  const WB::Params three_quarters = staggered(0.75f);
  HS_EXPECT_EQ(three_quarters.source.speed, to.source.speed);
  HS_EXPECT_GT(three_quarters.color.palette_chroma, from.color.palette_chroma);
  HS_EXPECT_LT(three_quarters.color.palette_chroma, to.color.palette_chroma);

  const WB::Params end = staggered(1.0f);
  HS_EXPECT_EQ(end.source.speed, to.source.speed);
  HS_EXPECT_EQ(end.color.palette_chroma, to.color.palette_chroma);

  WB::Params parallel;
  parallel.lerp(from, to, 0.5f);
  HS_EXPECT_LT(parallel.source.speed, half.source.speed);
  HS_EXPECT_GT(parallel.color.palette_chroma, half.color.palette_chroma);
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

    WB::RequestedConfig low_frequency = WB::active_config(sb);
    low_frequency.params.source.pattern_freq = 0.1f;
    HS_EXPECT_TRUE(WB::valid_config(low_frequency));
    low_frequency.params.source.pattern_freq = 0.099f;
    HS_EXPECT_FALSE(WB::valid_config(low_frequency));

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
    candidate.slots.palette = static_cast<WB::PaletteMode>(invalid_tag);
    HS_EXPECT_FALSE(WB::valid_config(candidate));
    candidate = WB::legacy_config();
    candidate.slots.palette_mapping =
        static_cast<WB::PaletteMapping>(invalid_tag);
    HS_EXPECT_FALSE(WB::valid_config(candidate));
    candidate = WB::legacy_config();
    candidate.slots.brightness_envelope =
        static_cast<WB::BrightnessEnvelope>(invalid_tag);
    HS_EXPECT_FALSE(WB::valid_config(candidate));
    candidate = WB::legacy_config();
    candidate.params.color.mapping_frequency = 0.99f;
    HS_EXPECT_FALSE(WB::valid_config(candidate));
    candidate.params.color.mapping_frequency = 32.01f;
    HS_EXPECT_FALSE(WB::valid_config(candidate));
    candidate = WB::legacy_config();
    candidate.params.color.palette_chroma = 1.01f;
    HS_EXPECT_FALSE(WB::valid_config(candidate));
    candidate = WB::legacy_config();
    candidate.params.color.mapping_phase = 1.01f;
    HS_EXPECT_FALSE(WB::valid_config(candidate));
    candidate = WB::legacy_config();
    candidate.params.color.phase_oscillation_depth = 1.01f;
    HS_EXPECT_FALSE(WB::valid_config(candidate));
    candidate = WB::legacy_config();
    candidate.params.color.phase_oscillation_speed = 0.0101f;
    HS_EXPECT_FALSE(WB::valid_config(candidate));
    candidate = WB::legacy_config();
    candidate.params.color.brightness_depth = 1.01f;
    HS_EXPECT_FALSE(WB::valid_config(candidate));
    candidate = WB::legacy_config();
    candidate.params.color.value_opacity_low = -0.01f;
    HS_EXPECT_FALSE(WB::valid_config(candidate));
    candidate = WB::legacy_config();
    candidate.params.color.value_opacity_high = 1.01f;
    HS_EXPECT_FALSE(WB::valid_config(candidate));

    const WB::RequestedConfig legacy_config = WB::legacy_config();
    HS_EXPECT_TRUE(WB::valid_config(legacy_config));
    HS_EXPECT_TRUE(
        WB::transition_admitted(WB::active_config(sb), legacy_config));
    HS_EXPECT_FALSE(WB::has_inverse_program(legacy_config));
    WB::request_config(sb, legacy_config);
    HS_EXPECT_TRUE(WB::active_config(sb) == legacy_config);
    HS_EXPECT_EQ(WB::active_pipeline(sb), WB::InversePipelineId::NONE);

    const WB::RequestedConfig compiled = WB::presets()[0];
    HS_EXPECT_TRUE(WB::valid_config(compiled));
    HS_EXPECT_TRUE(WB::has_inverse_program(compiled));
    HS_EXPECT_TRUE(WB::transition_admitted(WB::active_config(sb), compiled));
    WB::request_config(sb, compiled);
    HS_EXPECT_TRUE(WB::active_config(sb) == compiled);
    HS_EXPECT_FALSE(WB::transition_active(sb));
    HS_EXPECT_TRUE(WB::requested_slots(sb) == compiled.slots);

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
    shared_owner.params.warp.outer.speed = 0.0f;
    shared_owner.params.warp.inner.speed = 0.0f;
    shared_owner.slots.warp_program.outer.resource_id = 6;
    shared_owner.slots.warp_program.inner.resource_id = 6;
    HS_EXPECT_TRUE(WB::valid_config(shared_owner));
    shared_owner.slots.warp_program.inner.basis = WB::NoiseBasis::FBM3;
    HS_EXPECT_TRUE(WB::valid_config(shared_owner));
    HS_EXPECT_TRUE(WB::transition_admitted(shared_owner, shared_owner));
  }

  {
    reset_effect_globals();
    WB::SB lens_change;
    lens_change.init();
    const auto dodecahedral = WB::SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL;
    HS_EXPECT_TRUE(
        lens_change.updateParameter("Lens", static_cast<float>(dodecahedral)) ==
        ParamSetResult::APPLIED);
    lens_change.draw_frame();
    lens_change.advance_display();
    HS_EXPECT_FALSE(WB::transition_active(lens_change));
    HS_EXPECT_EQ(WB::active_slots(lens_change).surface_lens, dodecahedral);
    HS_EXPECT_EQ(WB::active_slots(lens_change).warp_program.outer.kind,
                 WB::WarpStageKind::WAVE_SHEAR);
    HS_EXPECT_TRUE(WB::parameter_warning(lens_change, "Lens") == nullptr);
    HS_EXPECT_TRUE(
        lens_change.updateParameter(
            "Planar Warp 1", static_cast<float>(WB::WarpStageKind::NONE)) ==
        ParamSetResult::APPLIED);
    lens_change.draw_frame();
    lens_change.advance_display();
    HS_EXPECT_EQ(WB::active_slots(lens_change).surface_lens, dodecahedral);
    HS_EXPECT_EQ(WB::active_slots(lens_change).warp_program.outer.kind,
                 WB::WarpStageKind::NONE);
    HS_EXPECT_TRUE(WB::parameter_warning(lens_change, "Lens") == nullptr);
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
    uint64_t hash = hs_test::FNV1A64_BASIS;
    auto append = [&](const char *text) {
      for (; *text != '\0'; ++text)
        hash = hs_test::fnv1a64_byte(hash, static_cast<uint8_t>(*text));
    };
    for (const Effect::ParamDef &def : sb.getParameters()) {
      append(def.name);
      hash =
          hs_test::fnv1a64_byte(hash, static_cast<uint8_t>(def.option_count));
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
    const char *names[] = {"Projection", "Projection", "Planar Warp 1"};
    const float values[] = {static_cast<float>(WB::Projection::BONNE),
                            static_cast<float>(WB::Projection::AIROCEAN),
                            static_cast<float>(WB::WarpStageKind::POLAR_CHART)};
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
  HS_EXPECT_EQ(idle.configs[2].slots.projection, WB::Projection::AIROCEAN);
  HS_EXPECT_EQ(idle.configs[2].slots.warp_program.outer.kind,
               WB::WarpStageKind::POLAR_CHART);
}

/** @brief Mode-specific ranges clamp stale subordinate values. */
inline void test_shaderball_mode_specific_parameter_warnings() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();

  HS_EXPECT_EQ(
      sb.updateParameter("Planar Warp 1",
                         static_cast<float>(WB::WarpStageKind::VECTOR_NOISE)),
      ParamSetResult::APPLIED);
  HS_EXPECT_EQ(sb.updateParameter("Planar Warp 1 Scale", 1.0f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(sb.updateParameter("Planar Warp 1 Speed", 0.0f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(sb.updateParameter("Planar Warp 1 Strength", 1.0f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(
      sb.updateParameter("Planar Warp 1",
                         static_cast<float>(WB::WarpStageKind::CURL_FLOW)),
      ParamSetResult::APPLIED);

  const auto *strength = sb.getParameters().find("Planar Warp 1 Strength");
  HS_EXPECT_TRUE(strength != nullptr);
  const float limit = 0.5f / 4.0f;
  HS_EXPECT_NEAR(strength->min, -limit, 1e-7f);
  HS_EXPECT_NEAR(strength->max, limit, 1e-7f);
  HS_EXPECT_EQ(strength->get_requested(), limit);
  HS_EXPECT_TRUE(WB::parameter_warning(sb, "Planar Warp 1 Strength") ==
                 nullptr);
  HS_EXPECT_TRUE(WB::parameter_warning(sb, "Planar Warp 1") == nullptr);
  sb.draw_frame();
  sb.advance_display();
  HS_EXPECT_EQ(WB::active_slots(sb).warp_program.outer.kind,
               WB::WarpStageKind::CURL_FLOW);
  HS_EXPECT_EQ(WB::active_config(sb).params.warp.outer.strength, limit);
  HS_EXPECT_EQ(WB::active_pipeline(sb), WB::InversePipelineId::NONE);
}

/** @brief A function edit preserves both warp stages in the dodecahedral hold. */
inline void test_shaderball_dodecahedral_lattice_edit() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();
  HS_EXPECT_TRUE(sb.selectPreset(11));

  const WB::RequestedConfig before = WB::requested_config(sb);
  HS_EXPECT_EQ(before.slots.function, WB::Function::GRID);
  HS_EXPECT_EQ(before.slots.surface_lens,
               WB::SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL);
  HS_EXPECT_EQ(before.slots.warp_program.outer.kind,
               WB::WarpStageKind::VECTOR_NOISE);
  HS_EXPECT_EQ(before.slots.warp_program.inner.kind,
               WB::WarpStageKind::MIRROR_TILE);

  HS_EXPECT_TRUE(
      sb.updateParameter("Function",
                         static_cast<float>(WB::Function::PRIMITIVE_LATTICE)) ==
      ParamSetResult::APPLIED);
  WB::RequestedConfig expected = before;
  expected.slots.function = WB::Function::PRIMITIVE_LATTICE;
  HS_EXPECT_TRUE(WB::requested_config(sb) == expected);
  HS_EXPECT_TRUE(WB::parameter_warning(sb, "Function") == nullptr);
  HS_EXPECT_EQ(sb.getParameters().find("Planar Warp 1")->get(),
               static_cast<float>(WB::WarpStageKind::VECTOR_NOISE));
  HS_EXPECT_EQ(sb.getParameters().find("Planar Warp 2")->get(),
               static_cast<float>(WB::WarpStageKind::MIRROR_TILE));

  sb.draw_frame();
  sb.advance_display();
  WB::refresh_display(sb);
  HS_EXPECT_TRUE(WB::active_config(sb) == expected);
  HS_EXPECT_EQ(sb.getParameters().find("Planar Warp 1")->get(),
               static_cast<float>(WB::WarpStageKind::VECTOR_NOISE));
  HS_EXPECT_EQ(sb.getParameters().find("Planar Warp 2")->get(),
               static_cast<float>(WB::WarpStageKind::MIRROR_TILE));
}

/** @brief Rejected Polar Chart selectors expose controls needed for admission. */
inline void test_shaderball_polar_gui_repair() {
  using WB = ShaderBallWhiteBox;
  auto repair = [](bool first, bool lattice) {
    reset_effect_globals();
    WB::SB sb;
    sb.init();
    WB::RequestedConfig base = WB::legacy_config();
    base.slots.function =
        lattice ? WB::Function::PRIMITIVE_LATTICE : WB::Function::GRID;
    base.slots.warp_program.outer.kind = WB::WarpStageKind::NONE;
    base.slots.warp_program.inner.kind = WB::WarpStageKind::NONE;
    base.params.source.pattern_freq = 1.3f;
    base.params.source.lattice_cell_scale = 1.0f;
    WB::request_config(sb, base);
    const WB::RequestedConfig rendered = WB::active_config(sb);
    const char *root = first ? "Planar Warp 1" : "Planar Warp 2";
    const char *mode =
        first ? "Planar Warp 1 Polar Mode" : "Planar Warp 2 Polar Mode";
    const char *harmonic =
        first ? "Planar Warp 1 Polar Harmonic" : "Planar Warp 2 Polar Harmonic";
    const char *radial_scale =
        first ? "Planar Warp 1 Radial Scale" : "Planar Warp 2 Radial Scale";
    const char *density = lattice ? "Lattice Cell Scale" : "Pattern Freq";
    const float repaired = lattice ? 8.0f / (TWO_PI_F * 2.0f) : 1.5f;

    HS_EXPECT_TRUE(
        sb.updateParameter(
            root, static_cast<float>(WB::WarpStageKind::POLAR_CHART)) ==
        ParamSetResult::APPLIED);
    HS_EXPECT_TRUE(WB::active_config(sb) == rendered);
    HS_EXPECT_TRUE(sb.getParameters().find(density) != nullptr);
    HS_EXPECT_TRUE(sb.getParameters().find(mode) != nullptr);
    HS_EXPECT_TRUE(sb.getParameters().find(harmonic) != nullptr);
    HS_EXPECT_TRUE(sb.getParameters().find(radial_scale) != nullptr);
    HS_EXPECT_TRUE(WB::parameter_warning(sb, root) != nullptr);
    HS_EXPECT_EQ(sb.getParameters().find(root)->get_requested(),
                 static_cast<float>(WB::WarpStageKind::POLAR_CHART));
    HS_EXPECT_EQ(static_cast<const Effect &>(sb).accepted_parameter_value(
                     *sb.getParameters().find(root)),
                 static_cast<float>(WB::WarpStageKind::NONE));

    HS_EXPECT_TRUE(sb.updateParameter(harmonic, 2.0f) ==
                   ParamSetResult::APPLIED);
    HS_EXPECT_TRUE(WB::active_config(sb) == rendered);
    HS_EXPECT_TRUE(WB::parameter_warning(sb, root) != nullptr);
    HS_EXPECT_TRUE(sb.updateParameter(density, repaired) ==
                   ParamSetResult::APPLIED);
    HS_EXPECT_TRUE(WB::parameter_warning(sb, root) == nullptr);
    sb.draw_frame();
    sb.advance_display();
    const WB::RequestedConfig active = WB::active_config(sb);
    const WB::WarpStageSpec &polar = first ? active.slots.warp_program.outer
                                           : active.slots.warp_program.inner;
    HS_EXPECT_EQ(polar.kind, WB::WarpStageKind::POLAR_CHART);
    HS_EXPECT_EQ(WB::active_pipeline(sb), WB::InversePipelineId::NONE);
    HS_EXPECT_EQ(lattice ? active.params.source.lattice_cell_scale
                         : active.params.source.pattern_freq,
                 repaired);
  };

  repair(true, false);
  repair(false, true);

  reset_effect_globals();
  WB::SB shear;
  shear.init();
  HS_EXPECT_TRUE(shear.selectPreset(10));
  shear.draw_frame();
  shear.advance_display();
  WB::settle_transition(shear);
  HS_EXPECT_EQ(WB::active_slots(shear).warp_program.outer.kind,
               WB::WarpStageKind::POLAR_CHART);
  HS_EXPECT_EQ(WB::active_slots(shear).warp_program.inner.kind,
               WB::WarpStageKind::WAVE_SHEAR);
  HS_EXPECT_TRUE(WB::parameter_warning(shear, "Planar Warp 1") == nullptr);
  HS_EXPECT_EQ(shear.updateParameter("Function",
                                     static_cast<float>(WB::Function::RINGS)),
               ParamSetResult::APPLIED);
  const char *shear_warning = WB::parameter_warning(shear, "Function");
  HS_EXPECT_TRUE(shear_warning != nullptr);
  HS_EXPECT_TRUE(std::strstr(shear_warning, "polar-periodic Function") !=
                 nullptr);
  HS_EXPECT_TRUE(std::strstr(shear_warning, "Planar Warp 2") == nullptr);
}

/** @brief Structural admission accepts curated holds and heavy stage tuples. */
inline void test_shaderball_structural_admission() {
  using WB = ShaderBallWhiteBox;
  const auto &presets = WB::presets();
  for (const auto &preset : presets)
    HS_EXPECT_TRUE(WB::valid_config(preset));

  WB::RequestedConfig integrated_ridged = WB::legacy_config();
  integrated_ridged.slots.warp_program.outer.kind =
      WB::WarpStageKind::CURL_FLOW;
  integrated_ridged.slots.warp_program.outer.basis = WB::NoiseBasis::RIDGED3;
  integrated_ridged.slots.warp_program.outer.curl_integrator =
      WB::CurlIntegrator::MIDPOINT_2;
  integrated_ridged.params.warp.outer.scale = 1.0f;
  integrated_ridged.params.warp.outer.strength = 0.0078125f;
  integrated_ridged.params.warp.outer.speed = 0.0f;
  HS_EXPECT_TRUE(WB::valid_config(integrated_ridged));

  WB::RequestedConfig peirce_polar = WB::legacy_config();
  peirce_polar.slots.function = WB::Function::GRID;
  peirce_polar.slots.projection = WB::Projection::PEIRCE_QUINCUNCIAL;
  peirce_polar.slots.coverage = WB::CoveragePolicy::EDGE_FADE;
  peirce_polar.params.value.edge_width = 0.1f;
  peirce_polar.slots.warp_program.inner.kind = WB::WarpStageKind::POLAR_CHART;
  peirce_polar.slots.warp_program.inner.polar_harmonic = 2;
  peirce_polar.params.warp.inner.speed = 0.0f;
  peirce_polar.params.source.pattern_freq = 1.0f;
  HS_EXPECT_TRUE(WB::valid_config(peirce_polar));
  peirce_polar.slots.warp_program.outer.kind = WB::WarpStageKind::MIRROR_TILE;
  HS_EXPECT_FALSE(WB::valid_config(peirce_polar));
  peirce_polar.slots.warp_program.outer.kind = WB::WarpStageKind::NONE;

  WB::RequestedConfig airocean_mobius = WB::legacy_config();
  airocean_mobius.slots.projection = WB::Projection::AIROCEAN;
  airocean_mobius.slots.coverage = WB::CoveragePolicy::EDGE_FADE;
  airocean_mobius.params.value.edge_width = 0.1f;
  airocean_mobius.slots.surface_lens = WB::SurfaceLens::MOBIUS;
  airocean_mobius.params.surface_lens.mix = 1.0f;
  HS_EXPECT_TRUE(WB::valid_config(airocean_mobius));

  for (WB::Projection projection :
       {WB::Projection::BONNE, WB::Projection::PEIRCE_QUINCUNCIAL,
        WB::Projection::AIROCEAN}) {
    WB::RequestedConfig strict = WB::legacy_config();
    strict.slots.projection = projection;
    strict.slots.warp_program.outer.kind = WB::WarpStageKind::NONE;
    strict.slots.surface_lens = WB::SurfaceLens::NONE;
    strict.params.surface_lens.mix = 0.0f;
    HS_EXPECT_TRUE(WB::valid_config(strict));
  }
  WB::RequestedConfig airo_mobius = WB::legacy_config();
  airo_mobius.slots.projection = WB::Projection::AIROCEAN;
  airo_mobius.slots.warp_program.outer.kind = WB::WarpStageKind::NONE;
  airo_mobius.slots.surface_lens = WB::SurfaceLens::MOBIUS;
  airo_mobius.params.surface_lens.mix = 0.5f;
  HS_EXPECT_TRUE(WB::valid_config(airo_mobius));

  WB::RequestedConfig from = WB::legacy_config();
  from.slots.surface_lens = WB::SurfaceLens::NONE;
  from.params.surface_lens.mix = 0.0f;
  from.slots.warp_program.outer.kind = WB::WarpStageKind::CURL_FLOW;
  from.slots.warp_program.outer.basis = WB::NoiseBasis::SIMPLEX;
  from.slots.warp_program.outer.curl_integrator = WB::CurlIntegrator::EULER_1;
  from.params.warp.outer.strength = 0.0078125f;
  from.params.warp.outer.scale = 1.0f;
  from.params.warp.outer.speed = 0.0f;
  WB::RequestedConfig to = from;
  to.params.warp.outer.strength = 0.5f;
  to.params.warp.outer.scale = 1.0f / 64.0f;
  HS_EXPECT_TRUE(WB::valid_config(from));
  HS_EXPECT_TRUE(WB::valid_config(to));
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
    config.slots.function = WB::Function::GRID;

    config.slots.warp_program.outer.kind = WB::WarpStageKind::VECTOR_NOISE;
    config.params.warp.outer.scale = 1.0f;
    config.params.warp.outer.strength = 0.1f;
    config.params.warp.outer.speed = 0.0f;
    HS_EXPECT_FALSE(WB::valid_config(config));
    config.slots.warp_program.outer.kind = WB::WarpStageKind::CURL_FLOW;
    config.params.warp.outer.scale = 1.0f;
    config.params.warp.outer.strength = 1e-4f;
    HS_EXPECT_FALSE(WB::valid_config(config));
    config.slots.warp_program.outer.kind = WB::WarpStageKind::NONE;

    config.slots.surface_lens = WB::SurfaceLens::TANGENT_NOISE;
    config.params.surface_lens.mix = 0.5f;
    HS_EXPECT_FALSE(WB::valid_config(config));

    config.slots.surface_lens = WB::SurfaceLens::NONE;
    config.slots.projection = WB::Projection::SINUSOIDAL;
    HS_EXPECT_TRUE(WB::valid_config(config));
  }
}

/** @brief Additive warp metadata retains sub-ULP displacement at large
 *         coordinates, and reports zero when no colorizer reads it. */
inline void test_shaderball_additive_delta_precision() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();
  WB::FrameState frame = WB::frame(sb);
  frame.clocks.warp_outer_phase = 0.25f;
  frame.slots.hue_shift = WB::HueShiftMode::WARP_DISPLACEMENT;
  frame.prepared_hue_rotation.active = true;
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
  HS_EXPECT_NEAR(result.path_length, 0.001f, 1e-7f);
  HS_EXPECT_EQ(result.coords.im, input.im);

  frame.prepared_hue_rotation.active = false;
  const auto untracked = WB::warp_stage(input, projected, spec, params, frame);
  HS_EXPECT_EQ(untracked.coords.re, result.coords.re);
  HS_EXPECT_EQ(untracked.coords.im, result.coords.im);
  HS_EXPECT_EQ(untracked.delta.re, result.delta.re);
  HS_EXPECT_EQ(untracked.delta.im, result.delta.im);
  HS_EXPECT_EQ(untracked.path_length, 0.0f);
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
    HS_EXPECT_TRUE(WB::valid_config(WB::active_config(sb)));
    HS_EXPECT_FALSE(WB::transition_active(sb));
    HS_EXPECT_FALSE(WB::param_morph_active(sb));
    if (index == 5 || index == 11) {
      const auto projected = WB::surface_project(
          Vector(0.808122f, -0.303046f, 0.505076f), WB::frame(sb));
      HS_EXPECT_TRUE(std::isfinite(projected.fade_edge_distance));
    }
  }
}

/** @brief Manual navigation wraps and resumes automatic preset selection. */
inline void test_shaderball_manual_preset_navigation() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();
  HS_EXPECT_EQ(sb.getPresetCount(), size_t(18));
  HS_EXPECT_EQ(sb.getPresetIndex(), size_t(0));
  HS_EXPECT_TRUE(sb.previousPreset());
  HS_EXPECT_EQ(sb.getPresetIndex(), size_t(17));
  HS_EXPECT_TRUE(WB::active_config(sb) == WB::presets()[17]);
  HS_EXPECT_TRUE(sb.nextPreset());
  HS_EXPECT_EQ(sb.getPresetIndex(), size_t(0));
  HS_EXPECT_TRUE(WB::active_config(sb) == WB::presets()[0]);
  HS_EXPECT_TRUE(sb.animations_paused());

  sb.setAnimationsPaused(false);
  for (int frame = 0;
       frame < 200 && !WB::transition_active(sb) && !WB::param_morph_active(sb);
       ++frame) {
    sb.draw_frame();
    sb.advance_display();
  }
  HS_EXPECT_TRUE(WB::transition_active(sb) || WB::param_morph_active(sb));

  HS_EXPECT_TRUE(sb.synchronizePreset(5));
  HS_EXPECT_EQ(sb.getPresetIndex(), size_t(5));
  HS_EXPECT_TRUE(WB::active_config(sb) == WB::presets()[5]);
  HS_EXPECT_FALSE(sb.animations_paused());
  HS_EXPECT_FALSE(WB::transition_active(sb));
  HS_EXPECT_FALSE(WB::param_morph_active(sb));

  HS_EXPECT_TRUE(sb.selectPreset(6));
  HS_EXPECT_EQ(sb.getPresetIndex(), size_t(6));
  for (int frame = 0; frame < 8; ++frame) {
    sb.draw_frame();
    sb.advance_display();
  }
  HS_EXPECT_EQ(sb.getPresetIndex(), size_t(6));

  sb.setAnimationsPaused(false);
  sb.draw_frame();
  sb.advance_display();
  HS_EXPECT_EQ(sb.getPresetIndex(), size_t(7));
  HS_EXPECT_TRUE(WB::transition_active(sb) || WB::param_morph_active(sb));
}

/** @brief GUI values follow the rendered preset transition. */
inline void test_shaderball_preset_gui_transition() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();

  HS_EXPECT_TRUE(sb.selectPreset(8));
  sb.setAnimationsPaused(false);
  WB::begin_blend(sb);
  sb.setAnimationsPaused(true);
  HS_EXPECT_EQ(sb.getPresetIndex(), size_t(9));
  HS_EXPECT_TRUE(WB::param_morph_active(sb));
  const auto *noise_scale = sb.getParameters().find("Surface Noise Scale");
  HS_EXPECT_TRUE(noise_scale != nullptr);
  HS_EXPECT_NEAR(noise_scale->get(), 1.78815627f, 1e-6f);

  bool saw_intermediate = false;
  for (int frame = 0; frame < 512 && WB::param_morph_active(sb); ++frame) {
    sb.draw_frame();
    sb.advance_display();
    WB::refresh_display(sb);
    const float displayed = noise_scale->get();
    saw_intermediate |= displayed > 1.78815627f && displayed < 3.29720306f;
  }
  HS_EXPECT_TRUE(saw_intermediate);
  HS_EXPECT_NEAR(noise_scale->get(), 3.29720306f, 1e-5f);

  HS_EXPECT_TRUE(sb.selectPreset(6));
  sb.setAnimationsPaused(false);
  WB::begin_blend(sb);
  sb.setAnimationsPaused(true);
  HS_EXPECT_EQ(sb.getPresetIndex(), size_t(7));
  HS_EXPECT_TRUE(WB::transition_active(sb));
  const auto *function = sb.getParameters().find("Function");
  const auto *projection_wander = sb.getParameters().find("Projection Wander");
  HS_EXPECT_TRUE(function != nullptr);
  HS_EXPECT_TRUE(projection_wander != nullptr);
  HS_EXPECT_EQ(function->get(), static_cast<float>(WB::Function::GRID));
  HS_EXPECT_NEAR(projection_wander->get(), 1.0f, 1e-6f);
  saw_intermediate = false;
  for (int frame = 0; frame < 1024 && WB::transition_active(sb); ++frame) {
    sb.draw_frame();
    sb.advance_display();
    WB::refresh_display(sb);
    const float displayed = projection_wander->get();
    saw_intermediate |= displayed > 0.0030917525f && displayed < 1.0f;
  }
  HS_EXPECT_TRUE(saw_intermediate);
  HS_EXPECT_NEAR(projection_wander->get(), 0.0030917525f, 1e-5f);
  HS_EXPECT_EQ(function->get(),
               static_cast<float>(WB::Function::PRIMITIVE_LATTICE));
}

/** @brief No reachable slot combination registers past the parameter arena. */
inline void test_shaderball_parameter_capacity() {
  using WB = ShaderBallWhiteBox;
  using Slots = WB::Slots;
  reset_effect_globals();
  WB::SB sb;
  sb.init();
  WB::RequestedConfig densest = WB::legacy_config();

  // Registration is per slot, so sweeping each slot in turn and keeping its
  // densest option converges on the densest reachable configuration.
  auto maximize = [&](auto apply, int option_count) {
    size_t best_count = WB::registered_parameter_count(sb, densest);
    int best_option = -1;
    for (int option = 0; option < option_count; ++option) {
      WB::RequestedConfig candidate = densest;
      apply(candidate.slots, option);
      const size_t count = WB::registered_parameter_count(sb, candidate);
      HS_EXPECT_LE(count, WB::PARAM_CAPACITY);
      if (count > best_count) {
        best_count = count;
        best_option = option;
      }
    }
    if (best_option >= 0)
      apply(densest.slots, best_option);
  };

  for (int pass = 0; pass < 2; ++pass) {
    maximize([](Slots &s, int o) { s.function = static_cast<WB::Function>(o); },
             WB::NUM_FUNCTIONS);
    maximize(
        [](Slots &s, int o) { s.projection = static_cast<WB::Projection>(o); },
        WB::NUM_PROJECTIONS);
    maximize([](Slots &s,
                int o) { s.peirce_layout = static_cast<WB::PeirceLayout>(o); },
             WB::NUM_PEIRCE_LAYOUTS);
    maximize(
        [](Slots &s, int o) {
          s.airocean_layout = static_cast<WB::AiroceanLayout>(o);
        },
        WB::NUM_AIROCEAN_LAYOUTS);
    maximize(
        [](Slots &s, int o) {
          s.bonne_hemisphere = static_cast<WB::BonneHemisphere>(o);
        },
        WB::NUM_BONNE_HEMISPHERES);
    maximize(
        [](Slots &s, int o) {
          s.gnomonic_hemisphere = static_cast<WB::GnomonicHemispherePolicy>(o);
        },
        WB::NUM_GNOMONIC_HEMISPHERES);
    maximize(
        [](Slots &s, int o) {
          s.projection_frame = static_cast<WB::ProjectionFramePolicy>(o);
        },
        WB::NUM_PROJECTION_FRAMES);
    maximize([](Slots &s,
                int o) { s.surface_lens = static_cast<WB::SurfaceLens>(o); },
             WB::NUM_LENSES);
    maximize([](Slots &s,
                int o) { s.surface_noise = static_cast<WB::SurfaceNoise>(o); },
             WB::NUM_SURFACE_NOISE);
    maximize(
        [](Slots &s, int o) {
          s.surface_noise_placement = static_cast<WB::SurfaceNoisePlacement>(o);
        },
        WB::NUM_SURFACE_NOISE_PLACEMENTS);
    for (bool outer : {true, false}) {
      auto stage = [outer](Slots &s) -> WB::WarpStageSpec & {
        return outer ? s.warp_program.outer : s.warp_program.inner;
      };
      maximize(
          [stage](Slots &s, int o) {
            stage(s).kind = static_cast<WB::WarpStageKind>(o);
          },
          WB::NUM_WARPS);
      maximize(
          [stage](Slots &s, int o) {
            stage(s).basis = static_cast<WB::NoiseBasis>(o);
          },
          WB::NUM_NOISE_BASES);
      maximize(
          [stage](Slots &s, int o) {
            stage(s).envelope = static_cast<WB::WarpEnvelope>(o);
          },
          WB::NUM_WARP_ENVELOPES);
      maximize(
          [stage](Slots &s, int o) {
            stage(s).polar_mode = static_cast<WB::PolarMode>(o);
          },
          WB::NUM_POLAR_MODES);
      maximize(
          [stage](Slots &s, int o) {
            stage(s).curl_integrator = static_cast<WB::CurlIntegrator>(o);
          },
          WB::NUM_CURL_INTEGRATORS);
    }
    maximize([](Slots &s,
                int o) { s.signal_weight = static_cast<WB::SignalWeight>(o); },
             WB::NUM_SIGNALS);
    maximize(
        [](Slots &s, int o) {
          s.value_transfer = static_cast<WB::ValueTransfer>(o);
        },
        WB::NUM_VALUE_TRANSFERS);
    maximize([](Slots &s,
                int o) { s.coverage = static_cast<WB::CoveragePolicy>(o); },
             WB::NUM_COVERAGE_POLICIES);
    maximize(
        [](Slots &s, int o) { s.palette = static_cast<WB::PaletteMode>(o); },
        WB::NUM_PALETTES);
    maximize(
        [](Slots &s, int o) {
          s.palette_mapping = static_cast<WB::PaletteMapping>(o);
        },
        WB::NUM_PALETTE_MAPPINGS);
    maximize(
        [](Slots &s, int o) {
          s.brightness_envelope = static_cast<WB::BrightnessEnvelope>(o);
        },
        WB::NUM_BRIGHTNESS_ENVELOPES);
    maximize(
        [](Slots &s, int o) { s.hue_shift = static_cast<WB::HueShiftMode>(o); },
        WB::NUM_HUE_SHIFT_MODES);
  }

  HS_EXPECT_LE(WB::registered_parameter_count(sb, densest), WB::PARAM_CAPACITY);
}

/** @brief Every GUI enum option is writable and survives its handoff. */
inline void test_shaderball_gui_catalog() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();
  HS_EXPECT_LE(sb.getParameters().size(), size_t(80));
  const auto *function = sb.getParameters().find("Function");
  HS_EXPECT_EQ(function->option_count, 7);
  HS_EXPECT_TRUE(
      std::strcmp(function->options[4], "Noise Contour (Projected)") == 0);
  HS_EXPECT_TRUE(std::strcmp(function->options[6], "Noise Contour (Sphere)") ==
                 0);
  HS_EXPECT_TRUE(sb.getParameters().find("Complexity") != nullptr);
  HS_EXPECT_TRUE(sb.getParameters().find("Pattern Mix") != nullptr);
  HS_EXPECT_TRUE(sb.getParameters().find("Drift") != nullptr);
  auto parameter_index = [&](const char *name) {
    size_t index = 0;
    for (const auto &parameter : sb.getParameters()) {
      if (std::strcmp(parameter.name, name) == 0)
        return index;
      ++index;
    }
    return sb.getParameters().size();
  };
  HS_EXPECT_LT(parameter_index("Function"), parameter_index("Pattern Freq"));
  HS_EXPECT_LT(parameter_index("Pattern Freq"), parameter_index("Projection"));
  HS_EXPECT_LT(parameter_index("Projection"), parameter_index("Pole Fade"));
  HS_EXPECT_LT(parameter_index("Pole Fade"),
               parameter_index("Projection Frame"));
  HS_EXPECT_LT(parameter_index("Projection Frame"),
               parameter_index("Projection Spin Speed"));
  HS_EXPECT_LT(parameter_index("Camera Wander"),
               parameter_index("Surface Noise"));
  HS_EXPECT_LT(parameter_index("Surface Noise"), parameter_index("Lens"));
  HS_EXPECT_TRUE(sb.getParameters().find("Lens Mix") == nullptr);
  HS_EXPECT_LT(parameter_index("Lens"), parameter_index("Planar Warp 1"));
  HS_EXPECT_LT(parameter_index("Planar Warp 1"),
               parameter_index("Planar Warp 2"));
  HS_EXPECT_LT(parameter_index("Palette"), parameter_index("Palette Chroma"));
  HS_EXPECT_LT(parameter_index("Palette Chroma"),
               parameter_index("Palette Mapping"));
  HS_EXPECT_LT(parameter_index("Palette Mapping"),
               parameter_index("Mapping Frequency"));
  HS_EXPECT_LT(parameter_index("Mapping Frequency"),
               parameter_index("Brightness Envelope"));
  HS_EXPECT_LT(parameter_index("Brightness Envelope"),
               parameter_index("Hue Shift Mode"));
  HS_EXPECT_LT(parameter_index("Hue Shift Mode"),
               parameter_index("Hue Shift Amount"));
  HS_EXPECT_EQ(sb.getParameters().find("Palette")->option_count, 3);
  const auto *palette_mapping = sb.getParameters().find("Palette Mapping");
  HS_EXPECT_TRUE(palette_mapping != nullptr);
  HS_EXPECT_EQ(palette_mapping->option_count, 4);
  HS_EXPECT_TRUE(std::strcmp(palette_mapping->options[2], "Linear") == 0);
  const auto *brightness_envelope =
      sb.getParameters().find("Brightness Envelope");
  HS_EXPECT_TRUE(brightness_envelope != nullptr);
  HS_EXPECT_EQ(brightness_envelope->option_count, 5);
  HS_EXPECT_TRUE(std::strcmp(brightness_envelope->options[0], "None") == 0);
  HS_EXPECT_TRUE(std::strcmp(brightness_envelope->options[3], "Ascending") ==
                 0);
  const auto *mapping_frequency = sb.getParameters().find("Mapping Frequency");
  HS_EXPECT_TRUE(mapping_frequency != nullptr);
  HS_EXPECT_EQ(mapping_frequency->min, 1.0f);
  HS_EXPECT_EQ(mapping_frequency->max, 32.0f);
  const auto *hue_shift = sb.getParameters().find("Hue Shift Mode");
  HS_EXPECT_TRUE(hue_shift != nullptr);
  HS_EXPECT_EQ(hue_shift->option_count, 3);
  HS_EXPECT_TRUE(
      std::strcmp(hue_shift->options[2], "Total Warp Displacement") == 0);
  const auto *projection = sb.getParameters().find("Projection");
  HS_EXPECT_TRUE(projection != nullptr);
  HS_EXPECT_EQ(projection->option_count, 7);
  HS_EXPECT_TRUE(std::strcmp(projection->options[0], "Folded Sinusoidal") == 0);
  HS_EXPECT_TRUE(std::strcmp(projection->options[3], "Bonne") == 0);
  HS_EXPECT_TRUE(std::strcmp(projection->options[4], "Peirce Quincuncial") ==
                 0);
  HS_EXPECT_TRUE(std::strcmp(projection->options[5], "Dymaxion / Airocean") ==
                 0);
  HS_EXPECT_TRUE(std::strcmp(projection->options[6], "Equirectangular") == 0);
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
  HS_EXPECT_LT(parameter_index("Projection"),
               parameter_index("Bonne Hemisphere"));
  HS_EXPECT_LT(parameter_index("Bonne Standard Parallel"),
               parameter_index("Projection Frame"));
  sb.draw_frame();
  sb.advance_display();
  WB::settle_transition(sb);

  WB::RequestedConfig gui_base = WB::presets()[0];
  gui_base.params.source.pattern_freq = 5.0f;
  auto reset_gui = [&] {
    WB::request_config(sb, gui_base);
    WB::settle_transition(sb);
  };

  reset_gui();
  HS_EXPECT_TRUE(
      sb.updateParameter("Coverage",
                         static_cast<float>(WB::CoveragePolicy::EDGE_FADE)) ==
      ParamSetResult::APPLIED);
  const auto *edge_fade_width = sb.getParameters().find("Edge Fade Width");
  HS_EXPECT_TRUE(edge_fade_width != nullptr);
  HS_EXPECT_EQ(edge_fade_width->min, 0.0f);
  HS_EXPECT_TRUE(sb.updateParameter("Edge Fade Width", 0.0f) ==
                 ParamSetResult::APPLIED);
  sb.draw_frame();
  sb.advance_display();
  HS_EXPECT_EQ(WB::active_config(sb).params.value.edge_width, 0.0f);

  constexpr const char *ROOT_ENUMS[] = {
      "Function",      "Projection",          "Projection Frame",
      "Lens",          "Planar Warp 1",       "Planar Warp 2",
      "Signal Weight", "Value Transfer",      "Coverage",
      "Palette",       "Brightness Envelope", "Hue Shift Mode"};
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
  select_and_set_all("Function", 4, "Source Noise Basis");
  const auto *source_noise_scale =
      sb.getParameters().find("Source Noise Scale");
  const auto *source_noise_rate = sb.getParameters().find("Source Noise Speed");
  HS_EXPECT_TRUE(source_noise_scale != nullptr);
  HS_EXPECT_TRUE(source_noise_rate != nullptr);
  HS_EXPECT_EQ(source_noise_scale->min, 0.0f);
  HS_EXPECT_EQ(source_noise_scale->max, 2.0f);
  HS_EXPECT_EQ(source_noise_rate->min, -1.0f / 1024.0f);
  HS_EXPECT_EQ(source_noise_rate->max, 1.0f / 1024.0f);
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
  select_and_set_all("Planar Warp 1", 4, "Planar Warp 1 Noise Basis");
  select_and_set_all("Planar Warp 1", 4, "Planar Warp 1 Envelope");
  select_and_set_all("Planar Warp 1", 5, "Planar Warp 1 Curl Integrator");
  {
    reset_gui();
    HS_EXPECT_TRUE(sb.updateParameter(
                       "Planar Warp 1",
                       static_cast<float>(WB::WarpStageKind::VECTOR_NOISE)) ==
                   ParamSetResult::APPLIED);
    HS_EXPECT_TRUE(sb.updateParameter("Planar Warp 1 Scale", 0.1f) ==
                   ParamSetResult::APPLIED);
    HS_EXPECT_TRUE(
        sb.updateParameter("Planar Warp 1",
                           static_cast<float>(WB::WarpStageKind::CURL_FLOW)) ==
        ParamSetResult::APPLIED);
    sb.draw_frame();
    sb.advance_display();
    WB::settle_transition(sb);
    const auto *strength = sb.getParameters().find("Planar Warp 1 Strength");
    HS_EXPECT_TRUE(strength != nullptr);
    const float limit = strength->max;
    HS_EXPECT_NEAR(limit, 1.0f, 1e-6f);
    HS_EXPECT_EQ(strength->min, -limit);
    HS_EXPECT_TRUE(sb.updateParameter("Planar Warp 1 Strength", 4.0f) ==
                   ParamSetResult::APPLIED);
    sb.draw_frame();
    sb.advance_display();
    WB::settle_transition(sb);
    HS_EXPECT_EQ(WB::active_config(sb).params.warp.outer.strength, limit);
  }
  select_and_set_all("Planar Warp 1", 7, "Planar Warp 1 Polar Mode");
  select_and_set_all("Planar Warp 1", 7, "Planar Warp 1 Polar Harmonic");
  HS_EXPECT_TRUE(sb.getParameters().find("Pattern Freq") != nullptr);
  HS_EXPECT_TRUE(
      sb.updateParameter("Function", static_cast<float>(WB::Function::RINGS)) ==
      ParamSetResult::APPLIED);
  HS_EXPECT_EQ(WB::requested_config(sb).slots.warp_program.outer.kind,
               WB::WarpStageKind::POLAR_CHART);
  HS_EXPECT_TRUE(WB::parameter_warning(sb, "Function") != nullptr);
  HS_EXPECT_TRUE(
      sb.updateParameter("Planar Warp 1",
                         static_cast<float>(WB::WarpStageKind::NONE)) ==
      ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(WB::parameter_warning(sb, "Function") == nullptr);
  sb.draw_frame();
  sb.advance_display();
  HS_EXPECT_EQ(WB::active_slots(sb).function, WB::Function::RINGS);
  HS_EXPECT_EQ(WB::active_pipeline(sb), WB::InversePipelineId::NONE);
  HS_EXPECT_TRUE(sb.getParameters().find("Pattern Freq") != nullptr);
  select_and_set_all("Value Transfer", 3, "Band Count");
  HS_EXPECT_LT(parameter_index("Value Transfer"),
               parameter_index("Band Count"));
  HS_EXPECT_LT(parameter_index("Band Phase"), parameter_index("Coverage"));
}

/** @brief Polyhedral lenses expose controls at their chamber scale. */
inline void test_shaderball_lens_domain_ranges() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();
  WB::RequestedConfig config = WB::legacy_config();
  config.slots.function = WB::Function::GRID;
  config.slots.projection_frame = WB::ProjectionFramePolicy::SPIN_WANDER;
  config.slots.surface_lens = WB::SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL;
  config.slots.surface_noise = WB::SurfaceNoise::DIRECT;
  config.params.surface_noise.strength = 0.1f;
  config.params.warp.outer.scale = 4.0f;
  WB::request_config(sb, config);

  auto parameter = [&](const char *name) {
    const auto *definition = sb.getParameters().find(name);
    HS_EXPECT_TRUE(definition != nullptr);
    return definition;
  };

  HS_EXPECT_TRUE(
      sb.updateParameter("Lens", static_cast<float>(WB::SurfaceLens::NONE)) ==
      ParamSetResult::APPLIED);
  HS_EXPECT_EQ(parameter("Pattern Freq")->max, 20.0f);
  HS_EXPECT_EQ(parameter("Speed")->max, 5.0f);
  HS_EXPECT_EQ(parameter("Hue Noise Speed")->max, 0.001f);
  HS_EXPECT_TRUE(sb.updateParameter("Speed", 5.0f) == ParamSetResult::APPLIED);

  HS_EXPECT_TRUE(sb.updateParameter(
                     "Lens", static_cast<float>(
                                 WB::SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL)) ==
                 ParamSetResult::APPLIED);
  HS_EXPECT_EQ(parameter("Pattern Freq")->max, 20.0f);
  HS_EXPECT_EQ(parameter("Speed")->max, 0.5f);
  HS_EXPECT_EQ(parameter("Speed")->get_requested(), 0.5f);
  HS_EXPECT_EQ(parameter("Source Angle Speed")->max, 0.03f);
  HS_EXPECT_EQ(parameter("Drift")->max, 1.25f);
  HS_EXPECT_EQ(parameter("Projection Spin Speed")->max, 0.04f);
  HS_EXPECT_EQ(parameter("Surface Noise Scale")->max, 8.0f);
  HS_EXPECT_EQ(parameter("Surface Noise Speed")->max, 0.002f);
  HS_EXPECT_EQ(parameter("Hue Noise Scale")->max, 2.0f);
  HS_EXPECT_EQ(parameter("Hue Noise Speed")->max, 0.001f);

  HS_EXPECT_TRUE(
      sb.updateParameter("Planar Warp 1",
                         static_cast<float>(WB::WarpStageKind::WAVE_SHEAR)) ==
      ParamSetResult::APPLIED);
  HS_EXPECT_EQ(parameter("Planar Warp 1 Speed")->max, 0.005f);
  HS_EXPECT_EQ(parameter("Planar Warp 1 Frequency")->max, 8.0f);

  HS_EXPECT_TRUE(
      sb.updateParameter("Planar Warp 1",
                         static_cast<float>(WB::WarpStageKind::VECTOR_NOISE)) ==
      ParamSetResult::APPLIED);
  HS_EXPECT_EQ(parameter("Planar Warp 1 Scale")->max, 1.0f);
  HS_EXPECT_EQ(parameter("Planar Warp 1 Strength")->max, 1.0f);

  HS_EXPECT_TRUE(
      sb.updateParameter("Lens", static_cast<float>(WB::SurfaceLens::NONE)) ==
      ParamSetResult::APPLIED);
  HS_EXPECT_EQ(parameter("Planar Warp 1 Scale")->max, 4.0f);

  HS_EXPECT_TRUE(
      sb.updateParameter("Planar Warp 1",
                         static_cast<float>(WB::WarpStageKind::CURL_FLOW)) ==
      ParamSetResult::APPLIED);
  HS_EXPECT_EQ(parameter("Planar Warp 1 Scale")->max, 2.0f);
  HS_EXPECT_EQ(parameter("Planar Warp 1 Strength")->max, 0.125f);

  HS_EXPECT_TRUE(
      sb.updateParameter("Function",
                         static_cast<float>(WB::Function::PRIMITIVE_LATTICE)) ==
      ParamSetResult::APPLIED);
  HS_EXPECT_EQ(parameter("Lattice Cell Scale")->max, 8.0f);
}

/** @brief New cartographic kernels preserve landmarks and stay finite. */
inline void test_shaderball_projection_catalog() {
  using WB = ShaderBallWhiteBox;
  const float standard_parallel = PI_F * 0.25f;
  const Vector bonne_origin(cosf(standard_parallel), sinf(standard_parallel),
                            0.0f);
  const auto bonne =
      projections::bonne_projection(bonne_origin, 0.0f, standard_parallel);
  HS_EXPECT_NEAR(bonne.coords.re, 0.0f, 2e-5f);
  HS_EXPECT_NEAR(bonne.coords.im, 0.0f, 2e-5f);

  const auto peirce = projections::peirce_projection(UP, 0.0f, 1, 0.0f);
  HS_EXPECT_NEAR(peirce.coords.re, 0.0f, 2e-5f);
  HS_EXPECT_NEAR(peirce.coords.im, 0.0f, 2e-5f);

  const auto &center = projections::AIROCEAN_CENTERS[0];
  const auto airocean = projections::airocean_projection(
      Vector(center.x, center.z, center.y), 0.0f, false);
  const auto &triangle = projections::AIROCEAN_PLANAR_FACES[0];
  HS_EXPECT_NEAR(airocean.coords.re,
                 (triangle[0].x + triangle[1].x + triangle[2].x) / 3.0f, 2e-4f);
  HS_EXPECT_NEAR(airocean.coords.im,
                 (triangle[0].y + triangle[1].y + triangle[2].y) / 3.0f, 2e-4f);
  HS_EXPECT_EQ(airocean.region_id, uint8_t(0));

  const Vector equator_zero(1.0f, 0.0f, 0.0f);
  const Vector equator_east(0.0f, 0.0f, 1.0f);
  const Vector antimeridian(-1.0f, 0.0f, 0.0f);
  const auto bonne_equator =
      projections::bonne_projection(equator_zero, 0.0f, standard_parallel);
  HS_EXPECT_NEAR(bonne_equator.coords.re, 0.0f, 2e-5f);
  HS_EXPECT_NEAR(bonne_equator.coords.im, -0.7853981634f, 2e-5f);
  const auto bonne_east =
      projections::bonne_projection(equator_east, 0.0f, standard_parallel);
  HS_EXPECT_NEAR(bonne_east.coords.re, 1.3758501640f, 3e-5f);
  HS_EXPECT_NEAR(bonne_east.coords.im, -0.1378413458f, 3e-5f);
  const auto bonne_shifted = projections::bonne_projection(
      equator_east, 0.5f * PI_F, standard_parallel);
  HS_EXPECT_NEAR(bonne_shifted.coords.re, 0.0f, 2e-5f);
  HS_EXPECT_NEAR(bonne_shifted.coords.im, -0.7853981634f, 2e-5f);
  const auto bonne_cut =
      projections::bonne_projection(antimeridian, 0.0f, standard_parallel);
  HS_EXPECT_NEAR(bonne_cut.fade_edge_distance, 0.0f, 2e-5f);
  const auto werner =
      projections::bonne_projection(equator_east, 0.0f, 0.5f * PI_F);
  HS_EXPECT_NEAR(werner.coords.re, 1.3217795320f, 3e-5f);
  HS_EXPECT_NEAR(werner.coords.im, -0.8487048774f, 3e-5f);

  constexpr float PEIRCE_K = 1.8540746773013719f;
  const Vector south_pole(0.0f, -1.0f, 0.0f);
  const auto peirce_diamond =
      projections::peirce_projection(south_pole, 0.0f, 0, 0.0f);
  const auto peirce_square =
      projections::peirce_projection(south_pole, 0.0f, 1, 0.0f);
  const auto peirce_horizontal =
      projections::peirce_projection(south_pole, 0.0f, 2, 0.0f);
  const auto peirce_vertical =
      projections::peirce_projection(south_pole, 0.0f, 3, 0.0f);
  HS_EXPECT_NEAR(peirce_diamond.coords.re, 2.0f * PEIRCE_K, 3e-5f);
  HS_EXPECT_NEAR(peirce_diamond.coords.im, 0.0f, 3e-5f);
  HS_EXPECT_NEAR(peirce_square.coords.re, PEIRCE_K * 1.4142135624f, 3e-5f);
  HS_EXPECT_NEAR(peirce_square.coords.im, PEIRCE_K * 1.4142135624f, 3e-5f);
  HS_EXPECT_NEAR(peirce_horizontal.coords.re, PEIRCE_K, 3e-5f);
  HS_EXPECT_NEAR(peirce_vertical.coords.im, PEIRCE_K, 3e-5f);
  const auto peirce_scroll0 =
      projections::peirce_projection(equator_east, 0.0f, 2, 0.0f);
  const auto peirce_scroll1 =
      projections::peirce_projection(equator_east, 0.0f, 2, 1.0f);
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
  const auto peirce_scroll_mid = projections::peirce_projection(
      equator_east, 0.0f, 2, scroll_mid.layout_scroll);
  HS_EXPECT_NEAR(peirce_scroll_mid.coords.re, peirce_scroll0.coords.re, 3e-5f);
  HS_EXPECT_NEAR(peirce_scroll_mid.coords.im, peirce_scroll0.coords.im, 3e-5f);
  const auto peirce_zero_fade =
      projections::peirce_projection(equator_zero, 0.0f, 1, 0.0f);
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
    const auto mapped = projections::peirce_projection(peirce_oracle_point,
                                                       0.0f, layout, 0.13f);
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
    const auto before = projections::peirce_projection(
        lon_lat(tie.longitude - TIE_EPS, SOUTH_LATITUDE), 0.0f, 0, 0.0f);
    const auto exact = projections::peirce_projection(
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
        projections::peirce_projection(on_meridian, MERIDIAN, layout, 0.13f);
    const auto reference =
        projections::peirce_projection(zero_meridian, 0.0f, layout, 0.13f);
    HS_EXPECT_NEAR(shifted.coords.re, reference.coords.re, 2e-5f);
    HS_EXPECT_NEAR(shifted.coords.im, reference.coords.im, 2e-5f);
    HS_EXPECT_EQ(shifted.region_id, reference.region_id);
  }
  const auto airocean_shifted =
      projections::airocean_projection(on_meridian, MERIDIAN, false);
  const auto airocean_reference =
      projections::airocean_projection(zero_meridian, 0.0f, false);
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
      projections::airocean_projection(oracle_point, 0.0f, false);
  HS_EXPECT_NEAR(airocean_oracle.coords.re, 2.1265288136f, 4e-5f);
  HS_EXPECT_NEAR(airocean_oracle.coords.im, 3.6817439808f, 4e-5f);
  const auto airocean_horizontal =
      projections::airocean_projection(oracle_point, 0.0f, true);
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
        projections::airocean_projection(oracle.point, 0.0f, false);
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
        projections::airocean_projection(glued_points[index], 0.0f, false);
    HS_EXPECT_EQ(mapped.region_id, glued_faces[index]);
    HS_EXPECT_EQ(mapped.edge_class,
                 projections::airocean_edge_identity(glued_faces[index],
                                                     glued_edges[index]));
    HS_EXPECT_TRUE((mapped.traits & projections::projection_traits(
                                        projections::ProjectionTrait::GLUED)) !=
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
        projections::airocean_projection(cut_points[index], 0.0f, false);
    HS_EXPECT_EQ(mapped.region_id, cut_faces[index]);
    HS_EXPECT_EQ(mapped.edge_class, projections::airocean_edge_identity(
                                        cut_faces[index], cut_edges[index]));
    HS_EXPECT_TRUE((mapped.traits & projections::projection_traits(
                                        projections::ProjectionTrait::CUT)) !=
                   0);
    HS_EXPECT_NEAR(mapped.coords.re, cut_coords[index].re, 5e-5f);
    HS_EXPECT_NEAR(mapped.coords.im, cut_coords[index].im, 5e-5f);
    HS_EXPECT_NEAR(mapped.fade_edge_distance, cut_distances[index], 5e-5f);
  }
  for (uint8_t face = 0; face < 23; ++face) {
    const auto &face_center = projections::AIROCEAN_CENTERS[face];
    const auto mapped = projections::airocean_projection(
        Vector(face_center.x, face_center.z, face_center.y), 0.0f, false);
    HS_EXPECT_EQ(mapped.region_id, face);
  }
  HS_EXPECT_TRUE(projections::airocean_edge_is_cut(14, 0));
  HS_EXPECT_EQ(projections::airocean_edge_identity(14, 0), uint8_t(42));
  HS_EXPECT_EQ(projections::airocean_edge_identity(18, 0), uint8_t(54));

  const auto japan_edge_point = [](float from_weight) {
    const auto &a = projections::AIROCEAN_FACES[14][0];
    const auto &b = projections::AIROCEAN_FACES[14][1];
    const auto &center = projections::AIROCEAN_CENTERS[14];
    const float to_weight = 1.0f - from_weight;
    return Vector(from_weight * a.x + to_weight * b.x + 1e-4f * center.x,
                  from_weight * a.z + to_weight * b.z + 1e-4f * center.z,
                  from_weight * a.y + to_weight * b.y + 1e-4f * center.y)
        .normalized();
  };
  const auto japan_cut =
      projections::airocean_projection(japan_edge_point(0.75f), 0.0f, false);
  const auto japan_glued =
      projections::airocean_projection(japan_edge_point(0.25f), 0.0f, false);
  HS_EXPECT_EQ(japan_cut.region_id, uint8_t(14));
  HS_EXPECT_EQ(japan_cut.edge_class, uint8_t(42));
  HS_EXPECT_TRUE(
      (japan_cut.traits &
       projections::projection_traits(projections::ProjectionTrait::CUT)) != 0);
  HS_EXPECT_LT(japan_cut.fade_edge_distance, 1e-3f);
  HS_EXPECT_EQ(japan_glued.region_id, uint8_t(14));
  HS_EXPECT_EQ(japan_glued.edge_class, uint8_t(54));
  HS_EXPECT_TRUE(
      (japan_glued.traits & projections::projection_traits(
                                projections::ProjectionTrait::GLUED)) != 0);
  HS_EXPECT_GT(japan_glued.fade_edge_distance, 0.1f);
  auto same_point = [](const projections::AiroceanPoint &a,
                       const projections::AiroceanPoint &b) {
    return fabsf(a.x - b.x) <= 1e-6f && fabsf(a.y - b.y) <= 1e-6f;
  };
  for (uint8_t face = 0; face < 23; ++face) {
    for (uint8_t edge_index = 0; edge_index < 3; ++edge_index) {
      bool expected_cut = false;
      for (size_t index = 0; index < std::size(projections::AIROCEAN_CUT_FACES);
           ++index)
        expected_cut |= projections::AIROCEAN_CUT_FACES[index] == face &&
                        projections::AIROCEAN_CUT_EDGES[index] == edge_index;
      HS_EXPECT_EQ(projections::airocean_edge_is_cut(face, edge_index),
                   expected_cut);

      const auto &a = projections::AIROCEAN_PLANAR_FACES[face][edge_index];
      const auto &b =
          projections::AIROCEAN_PLANAR_FACES[face][(edge_index + 1) % 3];
      uint8_t expected_identity = face * 3 + edge_index;
      bool found = false;
      for (uint8_t candidate_face = 0; candidate_face < 23 && !found;
           ++candidate_face) {
        for (uint8_t candidate_edge = 0; candidate_edge < 3; ++candidate_edge) {
          const auto &c = projections::AIROCEAN_PLANAR_FACES[candidate_face]
                                                            [candidate_edge];
          const auto &d =
              projections::AIROCEAN_PLANAR_FACES[candidate_face]
                                                [(candidate_edge + 1) % 3];
          if ((same_point(a, c) && same_point(b, d)) ||
              (same_point(a, d) && same_point(b, c))) {
            expected_identity = candidate_face * 3 + candidate_edge;
            found = true;
            break;
          }
        }
      }
      HS_EXPECT_EQ(projections::airocean_edge_identity(face, edge_index),
                   expected_identity);
    }
  }

  const auto peirce_seam_a = projections::peirce_projection(
      lon_lat(0.25f * PI_F - TIE_EPS, SOUTH_LATITUDE), 0.0f, 1, 0.0f);
  const auto peirce_seam_b = projections::peirce_projection(
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
      const auto b = projections::bonne_projection(v, 0.37f, standard_parallel);
      HS_EXPECT_TRUE(std::isfinite(b.coords.re));
      HS_EXPECT_TRUE(std::isfinite(b.coords.im));
      HS_EXPECT_TRUE(std::isfinite(b.fade_edge_distance));
      for (uint8_t layout = 0; layout < 4; ++layout) {
        const auto p = projections::peirce_projection(v, -0.21f, layout, 0.13f);
        const auto without_edge =
            projections::peirce_projection(v, -0.21f, layout, 0.13f, false);
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
        const auto a = projections::airocean_projection(v, 0.19f, horizontal);
        const auto without_edge =
            projections::airocean_projection(v, 0.19f, horizontal, false);
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

/** Worst single-channel gap measured between the LUT mapper and the exact
 * gamut refinement, 16-bit scale. */
constexpr uint16_t MEASURED_HUE_LUT_MAX_ERROR = 5529;
/** Mean per-channel gap measured over the same sweep. */
constexpr uint64_t MEASURED_HUE_LUT_MEAN_ERROR = 203;
/** Headroom over the measured single-pixel worst case. */
constexpr uint16_t HUE_LUT_PEAK_HEADROOM = 2;
/** Headroom over the measured whole-sweep mean. */
constexpr uint64_t HUE_LUT_MEAN_HEADROOM = 2;
constexpr uint16_t HUE_LUT_MAX_CHANNEL_ERROR =
    MEASURED_HUE_LUT_MAX_ERROR * HUE_LUT_PEAK_HEADROOM;
constexpr uint64_t HUE_LUT_MEAN_CHANNEL_ERROR =
    MEASURED_HUE_LUT_MEAN_ERROR * HUE_LUT_MEAN_HEADROOM;

/** @brief The LUT-only hue mapper stays close to the exact gamut refinement. */
inline void test_shaderball_hue_rotate_lut_gamut() {
  uint16_t max_channel_error = 0;
  uint64_t total_error = 0;
  uint64_t channels = 0;
  for (uint32_t red = 0; red <= 65535; red += 8191) {
    for (uint32_t green = 0; green <= 65535; green += 8191) {
      for (uint32_t blue = 0; blue <= 65535; blue += 8191) {
        const Color4 input(Pixel(red, green, blue), 0.75f);
        for (int amount_step = -32; amount_step <= 32; ++amount_step) {
          const float amount = amount_step * (1.0f / 16.0f);
          const Color4 exact = hue_rotate(input, amount);
          const Color4 fast = ::hue_rotate_lut_gamut(input, amount);
          const uint16_t exact_channels[] = {exact.color.r, exact.color.g,
                                             exact.color.b};
          const uint16_t fast_channels[] = {fast.color.r, fast.color.g,
                                            fast.color.b};
          for (int channel = 0; channel < 3; ++channel) {
            const uint16_t a = exact_channels[channel];
            const uint16_t b = fast_channels[channel];
            const uint16_t error = a > b ? a - b : b - a;
            max_channel_error = std::max(max_channel_error, error);
            total_error += error;
            ++channels;
          }
          HS_EXPECT_EQ(fast.alpha, exact.alpha);
        }
      }
    }
  }
  std::printf("  [hue-lut] max=%u mean=%llu\n", max_channel_error,
              static_cast<unsigned long long>(total_error / channels));
  HS_EXPECT_LE(max_channel_error, HUE_LUT_MAX_CHANNEL_ERROR);
  HS_EXPECT_LE(total_error / channels, HUE_LUT_MEAN_CHANNEL_ERROR);
}

/** @brief The prepared hue field tracks direct palette conversion. */
inline void test_shaderball_prepared_hue_rotation() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();
  const WB::FrameState frame = WB::preset_frame(sb, 11);
  uint16_t max_channel_error = 0;
  uint64_t total_error = 0;
  uint64_t channels = 0;
  for (int value_step = 0; value_step <= 64; ++value_step) {
    const float value = value_step / 64.0f;
    for (int amount_step = -128; amount_step <= 128; ++amount_step) {
      const float amount = amount_step / 64.0f;
      const Pixel direct = WB::direct_hue_rotation(frame, value, amount);
      const Pixel prepared = WB::prepared_hue_rotation(frame, value, amount);
      const uint16_t direct_channels[] = {direct.r, direct.g, direct.b};
      const uint16_t prepared_channels[] = {prepared.r, prepared.g, prepared.b};
      for (int channel = 0; channel < 3; ++channel) {
        const uint16_t a = direct_channels[channel];
        const uint16_t b = prepared_channels[channel];
        const uint16_t error = a > b ? a - b : b - a;
        max_channel_error = std::max(max_channel_error, error);
        total_error += error;
        ++channels;
      }
    }
  }
  HS_EXPECT_LE(max_channel_error, uint16_t(5400));
  HS_EXPECT_LE(total_error / channels, uint64_t(220));
}

/** @brief The prepared spherical hue field tracks its simplex source. */
inline void test_shaderball_prepared_hue_noise() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();
  float max_error = 0.0f;
  double total_error = 0.0;
  uint64_t samples = 0;
  static constexpr std::array PRESETS_AND_PHASES = {
      std::pair{size_t(7), 13.0f}, std::pair{size_t(8), 14.0f},
      std::pair{size_t(11), 17.0f}};
  for (const auto &[preset, phase] : PRESETS_AND_PHASES) {
    WB::ClockState clocks = WB::clocks(sb);
    clocks.hue_noise_phase = 0.137f * phase;
    WB::set_clocks(sb, clocks);
    const WB::FrameState frame = WB::preset_frame(sb, preset);
    HS_EXPECT_TRUE(frame.prepared_hue_noise.active);
    for (int latitude_step = -48; latitude_step <= 48; ++latitude_step) {
      const float latitude = latitude_step * (0.5f * PI_F / 48.0f);
      const float radius = cosf(latitude);
      for (int longitude_step = -192; longitude_step < 192; ++longitude_step) {
        const float longitude = longitude_step * (PI_F / 192.0f);
        const Vector direction(radius * cosf(longitude), sinf(latitude),
                               radius * sinf(longitude));
        const float error = fabsf(WB::prepared_hue_noise(frame, direction) -
                                  WB::direct_hue_noise(frame, direction));
        max_error = std::max(max_error, error);
        total_error += error;
        ++samples;
      }
    }
  }
  HS_EXPECT_LE(max_error, 0.15f);
  HS_EXPECT_LE(total_error / samples, 0.013);
}

inline const PullbackManifest::OracleMetric &
pullback_oracle_metric(std::string_view oracle, std::string_view domain,
                       std::string_view aggregation) {
  for (const auto &metric : PullbackManifest::ORACLE_METRICS)
    if (metric.oracle_id == oracle && metric.domain == domain &&
        metric.aggregation == aggregation)
      return metric;
  HS_CHECK(false, "pullback oracle metric is missing from generated manifest");
  return PullbackManifest::ORACLE_METRICS[0];
}

/** @brief Prepared hue noise and rotation stay within the color-stage budget. */
inline void test_shaderball_prepared_hue_noise_color() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();
  uint16_t max_channel_error = 0;
  uint64_t total_error = 0;
  uint64_t channels = 0;
  static constexpr std::array PRESETS_AND_PHASES = {
      std::pair{size_t(7), 13.0f}, std::pair{size_t(8), 14.0f},
      std::pair{size_t(11), 17.0f}};
  for (const auto &[preset, phase] : PRESETS_AND_PHASES) {
    WB::ClockState clocks = WB::clocks(sb);
    clocks.hue_noise_phase = 0.137f * phase;
    WB::set_clocks(sb, clocks);
    const WB::FrameState frame = WB::preset_frame(sb, preset);
    for (int latitude_step = -24; latitude_step <= 24; ++latitude_step) {
      const float latitude = latitude_step * (0.5f * PI_F / 24.0f);
      const float radius = cosf(latitude);
      for (int longitude_step = -96; longitude_step < 96; ++longitude_step) {
        const float longitude = longitude_step * (PI_F / 96.0f);
        const Vector direction(radius * cosf(longitude), sinf(latitude),
                               radius * sinf(longitude));
        for (int value_step = 0; value_step <= 16; ++value_step) {
          const float value = value_step / 16.0f;
          const Pixel exact =
              WB::direct_hue_noise_color(frame, direction, value);
          const Pixel prepared =
              WB::prepared_hue_noise_color(frame, direction, value);
          const uint16_t exact_channels[] = {exact.r, exact.g, exact.b};
          const uint16_t prepared_channels[] = {prepared.r, prepared.g,
                                                prepared.b};
          for (int channel = 0; channel < 3; ++channel) {
            const uint16_t a = exact_channels[channel];
            const uint16_t b = prepared_channels[channel];
            const uint16_t error = a > b ? a - b : b - a;
            max_channel_error = std::max(max_channel_error, error);
            total_error += error;
            ++channels;
          }
        }
      }
    }
  }
  std::printf("  [hue-noise-color] max=%u mean=%llu\n", max_channel_error,
              static_cast<unsigned long long>(total_error / channels));
  const auto &maximum = pullback_oracle_metric("HUE_ROTATION_AND_NOISE_LUTS",
                                               "COLOR_CHANNEL", "MAXIMUM");
  const auto &mean = pullback_oracle_metric("HUE_ROTATION_AND_NOISE_LUTS",
                                            "COLOR_CHANNEL", "MEAN");
  HS_EXPECT_TRUE(maximum.measured);
  HS_EXPECT_TRUE(mean.measured);
  HS_EXPECT_EQ(max_channel_error,
               static_cast<uint16_t>(maximum.measured_baseline));
  HS_EXPECT_EQ(total_error / channels,
               static_cast<uint64_t>(mean.measured_baseline));
  HS_EXPECT_EQ(WB::color_metric_limit(0), maximum.accepted_limit);
  HS_EXPECT_EQ(WB::color_metric_limit(1), mean.accepted_limit);
}

/** @brief Fast square Peirce stays within renderer error and seam budgets. */
inline void test_shaderball_fast_peirce_square() {
  float max_coordinate_error = 0.0f;
  float max_edge_error = 0.0f;
  bool metadata_matches = true;
  for (int latitude_step = -64; latitude_step <= 64; ++latitude_step) {
    const float latitude = latitude_step * (0.5f * PI_F / 64.0f);
    const float radius = cosf(latitude);
    for (int longitude_step = -256; longitude_step < 256; ++longitude_step) {
      const float longitude = longitude_step * (PI_F / 256.0f);
      const Vector input(radius * cosf(longitude), sinf(latitude),
                         radius * sinf(longitude));
      const auto exact =
          projections::peirce_projection(input, 0.0f, 1, 0.0f, true);
      const auto fast = projections::peirce_projection_fast_square(input);
      max_coordinate_error =
          std::max(max_coordinate_error,
                   std::max(fabsf(fast.coords.re - exact.coords.re),
                            fabsf(fast.coords.im - exact.coords.im)));
      max_edge_error =
          std::max(max_edge_error,
                   fabsf(fast.fade_edge_distance - exact.fade_edge_distance));
      metadata_matches &= fast.region_id == exact.region_id &&
                          fast.component_id == exact.component_id &&
                          fast.boundary_flags == exact.boundary_flags &&
                          fast.flags == exact.flags &&
                          fast.traits == exact.traits &&
                          fast.edge_class == exact.edge_class;
    }
  }
  std::printf("  [peirce-fast-square] coord=%.9g edge=%.9g\n",
              static_cast<double>(max_coordinate_error),
              static_cast<double>(max_edge_error));
  const auto &coordinate = pullback_oracle_metric(
      "PEIRCE_FAST_SQUARE", "PROJECTED_COORDINATE", "MAXIMUM");
  const auto &edge = pullback_oracle_metric(
      "PEIRCE_FAST_SQUARE", "PROJECTED_EDGE_DISTANCE", "MAXIMUM");
  HS_EXPECT_TRUE(coordinate.measured);
  HS_EXPECT_TRUE(edge.measured);
  HS_EXPECT_EQ(max_coordinate_error, coordinate.measured_baseline);
  HS_EXPECT_EQ(max_edge_error, edge.measured_baseline);
  HS_EXPECT_EQ(ShaderBallWhiteBox::peirce_metric_limit(0),
               coordinate.accepted_limit);
  HS_EXPECT_EQ(ShaderBallWhiteBox::peirce_metric_limit(1), edge.accepted_limit);
  HS_EXPECT_TRUE(metadata_matches);

  constexpr float TIE_EPSILON = 2e-6f;
  constexpr float LATITUDE = -0.41f;
  const float radius = cosf(LATITUDE);
  for (float boundary :
       {-0.75f * PI_F, -0.25f * PI_F, 0.25f * PI_F, 0.75f * PI_F}) {
    // Half the snap width: at the band edge itself the exact kernel's snap
    // turns on atan2f round-trip error and goes different ways per quadrant,
    // so the two kernels' tie predicates are not comparable there.
    for (float offset : {-0.5f * TIE_EPSILON, 0.0f, 0.5f * TIE_EPSILON}) {
      const float longitude = boundary + offset;
      const Vector input(radius * cosf(longitude), sinf(LATITUDE),
                         radius * sinf(longitude));
      const auto exact =
          projections::peirce_projection(input, 0.0f, 1, 0.0f, true);
      const auto fast = projections::peirce_projection_fast_square(input);
      HS_EXPECT_EQ(fast.region_id, exact.region_id);
      HS_EXPECT_EQ(fast.edge_class, exact.edge_class);
      HS_EXPECT_NEAR(fast.coords.re, exact.coords.re, 1.2e-3f);
      HS_EXPECT_NEAR(fast.coords.im, exact.coords.im, 1.2e-3f);
      HS_EXPECT_NEAR(fast.fade_edge_distance, exact.fade_edge_distance, 2e-4f);
    }
  }
}

/** @brief The inverse manifest has unique, canonical, selectable programs. */
inline void test_shaderball_operator_catalog_census() {
  using WB = ShaderBallWhiteBox;
  static_assert(WB::NUM_FUNCTIONS == 7);
  static_assert(WB::NUM_PROJECTIONS == 7);
  static_assert(WB::NUM_LENSES == 13);
  static_assert(WB::NUM_SURFACE_NOISE == 3);
  static_assert(WB::NUM_SURFACE_CURL_INTEGRATORS == 3);
  static_assert(WB::NUM_WARPS == 8);
  static_assert(WB::NUM_NOISE_BASES == 3);
  static_assert(WB::NUM_POLAR_MODES == 2);
  static_assert(WB::NUM_CURL_INTEGRATORS == 3);
  static_assert(WB::NUM_WARP_ENVELOPES == 3);
  static_assert(WB::NUM_SIGNALS == 2);
  static_assert(WB::NUM_VALUE_TRANSFERS == 4);
  static_assert(WB::NUM_COVERAGE_POLICIES == 5);
  static_assert(WB::NUM_PALETTE_MAPPINGS == 4);
  static_assert(WB::NUM_BRIGHTNESS_ENVELOPES == 5);
  static_assert(WB::NUM_HUE_SHIFT_MODES == 3);
  static_assert(static_cast<uint8_t>(WB::SurfaceLens::TANGENT_NOISE) == 255);
  static_assert(static_cast<uint8_t>(WB::WarpStageKind::LEGACY_STEREO_NOISE) ==
                255);

  using SurfaceDirect = Pullback::Surface::DirectNoise<WB::SurfaceStateProvider,
                                                       NoiseBasis::RIDGED3>;
  using SurfaceCurl =
      Pullback::Surface::CurlNoise<WB::SurfaceStateProvider, NoiseBasis::FBM3,
                                   Pullback::Surface::Midpoint2x>;
  using Mobius = Pullback::Lens::Mobius<WB::LensStateProvider>;
  using Equirectangular =
      Pullback::Projection::Equirectangular<WB::ProjectionStateProvider>;
  using Vortex = Pullback::Warp::Vortex<WB::WarpStateProvider<true>>;
  using Curl =
      Pullback::Warp::CurlFlow<WB::WarpStateProvider<true>, NoiseBasis::RIDGED3,
                               Pullback::Warp::Midpoint4,
                               Pullback::Warp::EdgeFadeEnvelope>;
  using ProjectedNoise =
      Pullback::Source::ProjectedNoise<WB::SourceStateProvider,
                                       NoiseBasis::FBM3>;
  using SphericalNoise =
      Pullback::Source::SphericalNoise<WB::SourceStateProvider,
                                       NoiseBasis::SIMPLEX>;
  using Generated = Pullback::Color::GeneratedPalette<WB::ColorStateProvider>;
  static_assert(std::is_empty_v<SurfaceDirect> &&
                std::is_empty_v<SurfaceCurl> && std::is_empty_v<Mobius> &&
                std::is_empty_v<Equirectangular> && std::is_empty_v<Vortex> &&
                std::is_empty_v<Curl> && std::is_empty_v<ProjectedNoise> &&
                std::is_empty_v<SphericalNoise> && std::is_empty_v<Generated>);
  HS_EXPECT_TRUE(true);
}

inline void test_shaderball_inverse_pipeline_manifest() {
  using WB = ShaderBallWhiteBox;
  static_assert(WB::inverse_stage_contracts());
  reset_effect_globals();
  WB::SB sb;
  sb.init();

  HS_EXPECT_EQ(WB::inverse_program_count(), PullbackManifest::PROGRAMS.size());
  HS_EXPECT_TRUE(WB::inverse_programs_well_formed());
  uint32_t compiled_preset_mask = 0;
  for (size_t index = 0; index < PullbackManifest::PROGRAMS.size(); ++index) {
    const auto &expected = PullbackManifest::PROGRAMS[index];
    HS_EXPECT_EQ(static_cast<size_t>(WB::inverse_program_id(index)), index);
    HS_EXPECT_TRUE(expected.id == WB::inverse_program_name(index));
    HS_EXPECT_TRUE(expected.topology_key ==
                   WB::topology_values(WB::inverse_program_key(index)));
    uint32_t preset_mask = 0;
    for (size_t preset = 0; preset < WB::presets().size(); ++preset)
      if (WB::preset_program_id(preset) == WB::inverse_program_id(index))
        preset_mask |= 1U << preset;
    HS_EXPECT_EQ(preset_mask, expected.preset_mask);
    compiled_preset_mask |= preset_mask;
  }
  HS_EXPECT_EQ(compiled_preset_mask, uint32_t(0x3ffff));
  const auto &peirce_framebuffer =
      pullback_oracle_metric("PEIRCE_FAST_SQUARE", "FRAMEBUFFER", "MAXIMUM");
  const auto &hue_framebuffer = pullback_oracle_metric(
      "HUE_ROTATION_AND_NOISE_LUTS", "FRAMEBUFFER", "MAXIMUM");
  HS_EXPECT_TRUE(peirce_framebuffer.measured);
  HS_EXPECT_TRUE(hue_framebuffer.measured);
  HS_EXPECT_EQ(peirce_framebuffer.measured_baseline, 100.0f);
  HS_EXPECT_EQ(hue_framebuffer.measured_baseline, 2189.0f);
  HS_EXPECT_LE(peirce_framebuffer.measured_baseline,
               peirce_framebuffer.accepted_limit);
  HS_EXPECT_LE(hue_framebuffer.measured_baseline,
               hue_framebuffer.accepted_limit);
  HS_EXPECT_EQ(WB::peirce_metric_limit(2), peirce_framebuffer.accepted_limit);
  HS_EXPECT_EQ(WB::color_metric_limit(2), hue_framebuffer.accepted_limit);
  for (const WB::RequestedConfig &preset : WB::presets())
    HS_EXPECT_TRUE(WB::has_inverse_program(preset));

  WB::RequestedConfig canonical = WB::presets()[1];
  const WB::TopologyKey expected_key = WB::topology_key(canonical);
  canonical.slots.peirce_layout = WB::PeirceLayout::DIAMOND;
  canonical.slots.bonne_hemisphere = WB::BonneHemisphere::SOUTH;
  canonical.slots.gnomonic_hemisphere =
      WB::GnomonicHemispherePolicy::BACK_HEMISPHERE;
  canonical.params.source.noise_basis = WB::NoiseBasis::RIDGED3;
  canonical.slots.warp_program.outer.basis = WB::NoiseBasis::FBM3;
  canonical.slots.warp_program.outer.envelope = WB::WarpEnvelope::EDGE_FADE;
  canonical.slots.warp_program.outer.polar_mode = WB::PolarMode::LOGARITHMIC;
  canonical.slots.warp_program.outer.curl_integrator =
      WB::CurlIntegrator::MIDPOINT_4;
  canonical.slots.warp_program.outer.polar_harmonic = 13;
  HS_EXPECT_TRUE(WB::topology_key(canonical) == expected_key);

  WB::RequestedConfig no_surface_noise = WB::presets()[4];
  const WB::TopologyKey no_surface_noise_key =
      WB::topology_key(no_surface_noise);
  no_surface_noise.slots.surface_noise_placement =
      WB::SurfaceNoisePlacement::BEFORE_LENS;
  no_surface_noise.params.surface_noise.basis = WB::NoiseBasis::RIDGED3;
  no_surface_noise.params.surface_noise.integrator =
      WB::SurfaceCurlIntegrator::MIDPOINT_2X;
  HS_EXPECT_TRUE(WB::topology_key(no_surface_noise) == no_surface_noise_key);

  WB::TopologyKey shear_key = WB::topology_key(WB::presets()[0]);
  HS_EXPECT_EQ(shear_key.outer_warp, WB::WarpStageKind::WAVE_SHEAR);
  HS_EXPECT_TRUE(
      (WB::planar_warp_implements<WB::WarpStageKind::WAVE_SHEAR,
                                  WB::WarpStageKind::NONE>(shear_key)));
  for (WB::WarpEnvelope envelope :
       {WB::WarpEnvelope::PROJECTION_WEIGHT, WB::WarpEnvelope::EDGE_FADE}) {
    shear_key.outer_warp_envelope = envelope;
    HS_EXPECT_FALSE(
        (WB::planar_warp_implements<WB::WarpStageKind::WAVE_SHEAR,
                                    WB::WarpStageKind::NONE>(shear_key)));
  }

  WB::RequestedConfig unsupported = WB::presets()[0];
  unsupported.slots.surface_lens = WB::SurfaceLens::TWIST;
  HS_EXPECT_FALSE(WB::has_inverse_program(unsupported));
  WB::request_config(sb, unsupported);
  HS_EXPECT_TRUE(WB::active_config(sb) == unsupported);
  HS_EXPECT_TRUE(WB::requested_config(sb) == unsupported);
  HS_EXPECT_EQ(WB::active_pipeline(sb), WB::InversePipelineId::NONE);
}

/** @brief Compiled inverse programs match the host reference within rounding. */
inline void test_shaderball_inverse_program_equivalence() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();
  for (size_t preset_index = 0; preset_index < WB::presets().size();
       ++preset_index) {
    const WB::FrameState frame = WB::preset_frame(sb, preset_index);
    HS_CONTEXT("preset", static_cast<long long>(preset_index));
    for (int latitude_step = -32; latitude_step <= 32; ++latitude_step) {
      const float latitude = latitude_step * (0.5f * PI_F / 32.0f);
      const float radius = cosf(latitude);
      // Per-ring, not per-sample: a frame pushed inside the 256-step longitude
      // loop costs more than the shading it labels.
      HS_CONTEXT("latitude step", latitude_step);
      for (int longitude_step = 0; longitude_step < 256; ++longitude_step) {
        const float longitude = longitude_step * (TWO_PI_F / 256.0f);
        const Vector view(radius * cosf(longitude), sinf(latitude),
                          radius * sinf(longitude));
        const Color4 expected = WB::shade(view, frame);
        const Color4 actual = WB::pipeline_shade(view, frame);
        const uint16_t red_error = actual.color.r > expected.color.r
                                       ? actual.color.r - expected.color.r
                                       : expected.color.r - actual.color.r;
        const uint16_t green_error = actual.color.g > expected.color.g
                                         ? actual.color.g - expected.color.g
                                         : expected.color.g - actual.color.g;
        const uint16_t blue_error = actual.color.b > expected.color.b
                                        ? actual.color.b - expected.color.b
                                        : expected.color.b - actual.color.b;
        HS_EXPECT_LE(std::max({red_error, green_error, blue_error}),
                     uint16_t(1));
        HS_EXPECT_NEAR(actual.alpha, expected.alpha, 1e-6f);
      }
    }
  }

  WB::FrameState peirce = WB::preset_frame(sb, 5);
  peirce.params.projection.central_meridian = 0.375f;
  for (int step = 0; step < 64; ++step) {
    const float longitude = step * (TWO_PI_F / 64.0f);
    const Vector view(cosf(longitude), 0.25f, sinf(longitude));
    const Vector normalized = view.normalized();
    const Color4 expected = WB::shade(normalized, peirce);
    const Color4 actual = WB::pipeline_shade(normalized, peirce);
    const uint16_t red_error = actual.color.r > expected.color.r
                                   ? actual.color.r - expected.color.r
                                   : expected.color.r - actual.color.r;
    const uint16_t green_error = actual.color.g > expected.color.g
                                     ? actual.color.g - expected.color.g
                                     : expected.color.g - actual.color.g;
    const uint16_t blue_error = actual.color.b > expected.color.b
                                    ? actual.color.b - expected.color.b
                                    : expected.color.b - actual.color.b;
    HS_EXPECT_LE(std::max({red_error, green_error, blue_error}), uint16_t(1));
    HS_EXPECT_NEAR(actual.alpha, expected.alpha, 1e-6f);
  }
}

/** @brief Domain policies, gauges, and analytic admission reject unsafe tuples. */
inline void test_shaderball_projection_and_admission_contracts() {
  using WB = ShaderBallWhiteBox;
  const uint8_t periodic_traits = projections::projection_traits(
      projections::ProjectionTrait::GLUED, projections::ProjectionTrait::FOLDED,
      projections::ProjectionTrait::PERIODIC);
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

  auto kernel_lookup = [](const projections::ProjectionKernelResult &result) {
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
      projections::peirce_projection(sector_before, 0.0f, 0, 0.0f));
  const auto diamond_after = kernel_lookup(
      projections::peirce_projection(sector_after, 0.0f, 0, 0.0f));
  HS_EXPECT_FALSE(WB::join_compatible(diamond_before, diamond_after,
                                      WB::Projection::PEIRCE_QUINCUNCIAL));
  const auto horizontal_before = kernel_lookup(
      projections::peirce_projection(sector_before, 0.0f, 2, 0.0f));
  const auto horizontal_after = kernel_lookup(
      projections::peirce_projection(sector_after, 0.0f, 2, 0.0f));
  HS_EXPECT_FALSE(WB::join_compatible(horizontal_before, horizontal_after,
                                      WB::Projection::PEIRCE_QUINCUNCIAL));

  const auto glued_side_a = kernel_lookup(projections::airocean_projection(
      Vector(0.00321964224570043f, 0.902105871397194f, 0.431502758617508f),
      0.0f, false));
  const auto glued_side_b = kernel_lookup(projections::airocean_projection(
      Vector(0.00321096516044884f, 0.902110786482306f, 0.431492547577607f),
      0.0f, false));
  HS_EXPECT_FALSE(WB::join_compatible(glued_side_a, glued_side_b,
                                      WB::Projection::AIROCEAN));
  const auto cut_side_a = kernel_lookup(projections::airocean_projection(
      Vector(0.456076125629434f, 0.767836786220573f, -0.449912477441232f), 0.0f,
      false));
  const auto cut_side_b = kernel_lookup(projections::airocean_projection(
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
  mobius_config.params.surface_lens.mobius.a.re *= 2.0f;
  HS_EXPECT_TRUE(WB::valid_config(mobius_config));

  WB::RequestedConfig curl = WB::legacy_config();
  curl.slots.warp_program.outer.kind = WB::WarpStageKind::CURL_FLOW;
  curl.slots.warp_program.outer.curl_integrator =
      WB::CurlIntegrator::MIDPOINT_4;
  curl.params.warp.outer.scale = 1.0f;
  curl.params.warp.outer.speed = 0.0f;
  curl.params.warp.outer.strength = 0.5f;
  HS_EXPECT_TRUE(WB::valid_config(curl));
  curl.params.warp.outer.strength = 0.5001f;
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
  HS_EXPECT_NEAR(curl_two.re, curl_one.re, 1e-5f);
  HS_EXPECT_NEAR(curl_two.im, curl_one.im, 1e-5f);
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

/** @brief Frame-like planar warps complete one documented cycle per phase. */
inline void test_shaderball_planar_warp_animation() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();
  WB::FrameState frame = WB::frame(sb);
  frame.resources.outer_warp_noise = nullptr;
  const Complex input(0.31f, -0.27f);
  const WB::ProjectedLookup projected(input, 0, 0, 0, 1.0f, 1.0f, 0);

  auto sample = [&](WB::WarpStageKind kind, WB::WarpStageParams params,
                    float phase) {
    frame.clocks.warp_outer_phase = phase;
    const Complex source_period = kind == WB::WarpStageKind::AFFINE_FRAME
                                      ? Complex(1.0f, 1.0f)
                                      : Complex();
    return WB::warp_stage(input, projected, WB::WarpStageSpec{kind}, params,
                          frame, false, source_period);
  };
  auto expect_cycle = [&](WB::WarpStageKind kind, WB::WarpStageParams params) {
    const auto start = sample(kind, params, 0.0f);
    const auto quarter = sample(kind, params, 0.25f);
    const auto wrapped = sample(kind, params, 1.0f);
    HS_EXPECT_TRUE(fabsf(start.coords.re - quarter.coords.re) > 1e-4f ||
                   fabsf(start.coords.im - quarter.coords.im) > 1e-4f);
    HS_EXPECT_NEAR(start.coords.re, wrapped.coords.re, 2e-5f);
    HS_EXPECT_NEAR(start.coords.im, wrapped.coords.im, 2e-5f);
  };

  WB::WarpStageParams affine;
  affine.translation_x = 0.4f;
  affine.translation_y = -0.2f;
  affine.rotation = 0.6f;
  affine.scale_x = 1.5f;
  affine.scale_y = 0.75f;
  affine.shear = 0.25f;
  expect_cycle(WB::WarpStageKind::AFFINE_FRAME, affine);
  const auto affine_quarter =
      sample(WB::WarpStageKind::AFFINE_FRAME, affine, 0.25f);
  const auto affine_three_quarters =
      sample(WB::WarpStageKind::AFFINE_FRAME, affine, 0.75f);
  HS_EXPECT_TRUE(fabsf(affine_quarter.coords.re -
                       affine_three_quarters.coords.re) > 1e-4f ||
                 fabsf(affine_quarter.coords.im -
                       affine_three_quarters.coords.im) > 1e-4f);
  float WB::WarpStageParams::*const affine_fields[] = {
      &WB::WarpStageParams::translation_x, &WB::WarpStageParams::translation_y,
      &WB::WarpStageParams::rotation,      &WB::WarpStageParams::scale_x,
      &WB::WarpStageParams::scale_y,       &WB::WarpStageParams::shear,
  };
  for (size_t index = 0; index < std::size(affine_fields); ++index) {
    WB::WarpStageParams one_value;
    one_value.*affine_fields[index] = index == 3 || index == 4 ? 1.5f : 0.3f;
    const auto first = sample(WB::WarpStageKind::AFFINE_FRAME, one_value, 0.0f);
    const auto second =
        sample(WB::WarpStageKind::AFFINE_FRAME, one_value, 0.25f);
    HS_EXPECT_TRUE(fabsf(first.coords.re - second.coords.re) > 1e-4f ||
                   fabsf(first.coords.im - second.coords.im) > 1e-4f);
  }
  const auto affine_identity_start =
      sample(WB::WarpStageKind::AFFINE_FRAME, WB::WarpStageParams{}, 0.0f);
  const auto affine_identity_quarter =
      sample(WB::WarpStageKind::AFFINE_FRAME, WB::WarpStageParams{}, 0.25f);
  HS_EXPECT_NEAR(affine_identity_start.coords.re, input.re, 1e-6f);
  HS_EXPECT_NEAR(affine_identity_start.coords.im, input.im, 1e-6f);
  HS_EXPECT_NEAR(affine_identity_quarter.coords.re, input.re, 1e-6f);
  HS_EXPECT_NEAR(affine_identity_quarter.coords.im, input.im, 1e-6f);

  WB::WarpStageParams affine_scale;
  affine_scale.scale_x = 2.0f;
  affine_scale.scale_y = 2.0f;
  const auto scale_start =
      sample(WB::WarpStageKind::AFFINE_FRAME, affine_scale, 0.0f);
  const auto scale_quarter =
      sample(WB::WarpStageKind::AFFINE_FRAME, affine_scale, 0.25f);
  HS_EXPECT_NEAR(scale_start.coords.re, input.re * 0.5f, 1e-6f);
  HS_EXPECT_NEAR(scale_start.coords.im, input.im * 0.5f, 1e-6f);
  HS_EXPECT_NEAR(scale_quarter.coords.re, input.re, 1e-6f);
  HS_EXPECT_NEAR(scale_quarter.coords.im, input.im, 1e-6f);

  WB::WarpStageParams affine_scroll;
  affine_scroll.translation_x = 1.0f;
  affine_scroll.translation_y = -2.0f;
  const auto scroll_quarter =
      sample(WB::WarpStageKind::AFFINE_FRAME, affine_scroll, 0.25f);
  const auto scroll_half =
      sample(WB::WarpStageKind::AFFINE_FRAME, affine_scroll, 0.5f);
  const auto scroll_three_quarters =
      sample(WB::WarpStageKind::AFFINE_FRAME, affine_scroll, 0.75f);
  HS_EXPECT_NEAR(scroll_half.coords.re - scroll_quarter.coords.re,
                 scroll_three_quarters.coords.re - scroll_half.coords.re,
                 1e-6f);
  HS_EXPECT_NEAR(scroll_half.coords.im - scroll_quarter.coords.im,
                 scroll_three_quarters.coords.im - scroll_half.coords.im,
                 1e-6f);

  WB::RequestedConfig lattice_scroll = WB::presets()[7];
  lattice_scroll.params.source.lattice_cell_scale = 2.5f;
  lattice_scroll.params.warp.outer.translation_x = 2.0f;
  lattice_scroll.params.warp.outer.translation_y = -1.0f;
  lattice_scroll.params.warp.outer.rotation = 0.0f;
  lattice_scroll.params.warp.outer.scale_x = 1.0f;
  lattice_scroll.params.warp.outer.scale_y = 1.0f;
  lattice_scroll.params.warp.outer.shear = 0.0f;
  auto lattice_field = [&](float x, float y, float phase) {
    WB::ClockState clocks = WB::clocks(sb);
    clocks.warp_outer_phase = phase;
    WB::set_clocks(sb, clocks);
    const WB::FrameState lattice_frame = WB::config_frame(sb, lattice_scroll);
    const WB::ProjectedLookup lattice_projected(Complex(x, y), 0, 0, 0, 1.0f,
                                                1.0f, 0);
    const WB::PlanarWarpResult warped =
        WB::warp(lattice_projected, lattice_frame);
    return std::pair{
        warped.coords,
        WB::source(WB::ProjectedLookup(warped.coords, 0, 0, 0, 1.0f, 1.0f, 0),
                   lattice_frame)};
  };
  for (int index = 0; index < 16; ++index) {
    const float x = -0.7f + 0.093f * index;
    const float y = 0.51f - 0.071f * index;
    const auto start = lattice_field(x, y, 0.0f);
    const auto seam = lattice_field(x, y, 1.0f - 1e-6f);
    HS_EXPECT_TRUE(fabsf(start.first.re - seam.first.re) > 0.5f ||
                   fabsf(start.first.im - seam.first.im) > 0.2f);
    HS_EXPECT_NEAR(start.second, seam.second, 1e-3f);
  }
  HS_EXPECT_TRUE(WB::valid_config(lattice_scroll));
  WB::RequestedConfig changed_winding = lattice_scroll;
  changed_winding.params.warp.outer.translation_x = 1.0f;
  HS_EXPECT_TRUE(WB::valid_config(changed_winding));
  HS_EXPECT_TRUE(WB::transition_admitted(lattice_scroll, changed_winding));
  HS_EXPECT_FALSE(
      WB::stable_parameter_path_admitted(lattice_scroll, changed_winding));
  WB::RequestedConfig incompatible_scroll = lattice_scroll;
  incompatible_scroll.params.warp.outer.translation_x = 0.5f;
  HS_EXPECT_FALSE(WB::valid_config(incompatible_scroll));
  incompatible_scroll = lattice_scroll;
  incompatible_scroll.slots.function = WB::Function::GRID;
  HS_EXPECT_FALSE(WB::valid_config(incompatible_scroll));
  incompatible_scroll = lattice_scroll;
  incompatible_scroll.slots.warp_program.inner.kind =
      WB::WarpStageKind::MIRROR_TILE;
  HS_EXPECT_FALSE(WB::valid_config(incompatible_scroll));
  incompatible_scroll = lattice_scroll;
  incompatible_scroll.slots.hue_shift = WB::HueShiftMode::WARP_DISPLACEMENT;
  incompatible_scroll.params.color.hue_shift_amount = 0.5f;
  HS_EXPECT_FALSE(WB::valid_config(incompatible_scroll));

  HS_EXPECT_TRUE(sb.selectPreset(7));
  HS_EXPECT_EQ(sb.updateParameter("Planar Warp 1 Translation X", 0.6f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(WB::requested_config(sb).params.warp.outer.translation_x, 1.0f);
  HS_EXPECT_TRUE(WB::parameter_warning(sb, "Planar Warp 1 Translation X") ==
                 nullptr);

  WB::WarpStageParams mirror;
  mirror.cell_x = 1.0f;
  mirror.cell_y = 1.0f;
  expect_cycle(WB::WarpStageKind::MIRROR_TILE, mirror);

  WB::WarpStageParams polar;
  polar.radial_scale = 1.0f;
  const auto polar_start = sample(WB::WarpStageKind::POLAR_CHART, polar, 0.0f);
  const auto polar_quarter =
      sample(WB::WarpStageKind::POLAR_CHART, polar, 0.25f);
  HS_EXPECT_NEAR(polar_start.coords.re, polar_quarter.coords.re, 1e-6f);
  HS_EXPECT_NEAR(polar_quarter.coords.im - polar_start.coords.im, 0.5f * PI_F,
                 2e-5f);

  for (WB::WarpStageKind kind :
       {WB::WarpStageKind::AFFINE_FRAME, WB::WarpStageKind::MIRROR_TILE,
        WB::WarpStageKind::POLAR_CHART}) {
    HS_EXPECT_EQ(sb.updateParameter("Planar Warp 1", static_cast<float>(kind)),
                 ParamSetResult::APPLIED);
    HS_EXPECT_TRUE(sb.getParameters().find("Planar Warp 1 Speed") != nullptr);
  }
}

/** @brief Valid uncompiled GUI edits render through the dynamic backend. */
inline void test_shaderball_dynamic_backend_admission() {
  using WB = ShaderBallWhiteBox;
  HS_EXPECT_TRUE(WB::has_inverse_program(WB::presets()[0]));
  {
    reset_effect_globals();
    WB::SB sb;
    sb.init();
    HS_EXPECT_EQ(
        sb.updateParameter("Lens", static_cast<float>(WB::SurfaceLens::TWIST)),
        ParamSetResult::APPLIED);
    sb.draw_frame();
    sb.advance_display();
    HS_EXPECT_EQ(WB::active_slots(sb).surface_lens, WB::SurfaceLens::TWIST);
    HS_EXPECT_EQ(WB::requested_config(sb).slots.surface_lens,
                 WB::SurfaceLens::TWIST);
    HS_EXPECT_TRUE(WB::parameter_warning(sb, "Lens") == nullptr);
    HS_EXPECT_EQ(WB::active_pipeline(sb), WB::InversePipelineId::NONE);
    const WB::FullConfigSnapshot dynamic = sb.capture_full_config_snapshot();
    HS_EXPECT_TRUE(sb.selectPreset(1));
    HS_EXPECT_EQ(sb.restore_full_config_snapshot(dynamic),
                 WB::ConfigRestoreResult::APPLIED);
    HS_EXPECT_EQ(WB::active_pipeline(sb), WB::InversePipelineId::NONE);
    sb.draw_frame();
    sb.advance_display();
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

  // The sweep covers kernel families the authored roster compiles no inverse
  // pipeline for, so each candidate is shaded through the reference pullback
  // on a frame prepared straight from the config.
  auto check = [&](const WB::RequestedConfig &candidate) {
    HS_EXPECT_TRUE(WB::valid_config(candidate));
    const Color4 color = WB::shade(view, WB::config_frame(sb, candidate));
    HS_EXPECT_TRUE(std::isfinite(color.alpha));
    HS_EXPECT_GE(color.alpha, 0.0f);
    HS_EXPECT_LE(color.alpha, 1.0f);
  };

  for (uint8_t value = 0; value <= 5; ++value) {
    config.slots.function = static_cast<WB::Function>(value);
    config.params.source.noise_basis = WB::NoiseBasis::FBM3;
    check(config);
  }
  for (uint8_t value = 0; value <= 7; ++value) {
    config.slots.function =
        value == static_cast<uint8_t>(WB::WarpStageKind::AFFINE_FRAME)
            ? WB::Function::PRIMITIVE_LATTICE
        : value == static_cast<uint8_t>(WB::WarpStageKind::POLAR_CHART)
            ? WB::Function::GRID
            : WB::Function::TWIN_WAVE;
    config.slots.warp_program.outer.kind =
        static_cast<WB::WarpStageKind>(value);
    config.slots.projection = WB::Projection::SINUSOIDAL;
    config.params.warp.outer.strength = value == 0   ? 0.0f
                                        : value == 5 ? 0.005f
                                                     : 0.35f;
    config.params.warp.outer.turns = 0.4f;
    config.params.warp.outer.translation_x = value == 1 ? 1.0f : 0.0f;
    config.params.warp.outer.translation_y = value == 1 ? -1.0f : 0.0f;
    config.params.warp.outer.speed = 0.0f;
    config.slots.warp_program.outer.basis = WB::NoiseBasis::SIMPLEX;
    check(config);
  }
  config.slots.warp_program.outer.kind = WB::WarpStageKind::NONE;
  config.params.warp.outer.strength = 0.0f;
  config.slots.projection = WB::Projection::SINUSOIDAL;
  for (uint8_t value = 0; value <= 12; ++value) {
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
      for (uint8_t palette = 0; palette <= 2; ++palette) {
        config.slots.palette = static_cast<WB::PaletteMode>(palette);
        check(config);
      }
    }
  }
  for (uint8_t envelope = 0; envelope <= 3; ++envelope) {
    config.slots.brightness_envelope =
        static_cast<WB::BrightnessEnvelope>(envelope);
    check(config);
  }

  WB::FrameState frame = WB::frame(sb);
  frame.resources.outer_warp_noise = nullptr;
  WB::WarpStageParams zero_params;
  zero_params.strength = 0.0f;
  const Complex input(0.27f, -0.41f);
  for (uint8_t value = 0; value <= 5; ++value) {
    WB::WarpStageSpec spec{static_cast<WB::WarpStageKind>(value)};
    const auto identity = WB::warp_stage(input, {input, 0, 0, 0, 1.0f, 1.0f, 0},
                                         spec, zero_params, frame);
    HS_EXPECT_EQ(identity.coords.re, input.re);
    HS_EXPECT_EQ(identity.coords.im, input.im);
    HS_EXPECT_EQ(identity.path_length, 0.0f);
  }
  for (uint8_t value = 6; value <= 7; ++value) {
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
  HS_EXPECT_TRUE(sb.selectPreset(8));
  const WB::Params from = WB::presets()[8].params;
  const WB::Params to = WB::presets()[9].params;
  const size_t initial_index = WB::preset_index(sb);
  WB::begin_blend(sb);
  HS_EXPECT_EQ(WB::preset_index(sb), initial_index + 1);
  HS_EXPECT_TRUE(WB::param_morph_active(sb));
  HS_EXPECT_FALSE(WB::transition_active(sb));
  HS_EXPECT_TRUE(WB::live_params(sb) == from);

  float previous_phase = WB::clocks(sb).source_primary;
  for (int step = 0; step <= 480; ++step) {
    WB::step_param_morph(sb);
    const float live_speed = WB::live_params(sb).source.speed;
    const float phase = WB::clocks(sb).source_primary;
    HS_EXPECT_NEAR(phase, fmodf(previous_phase + live_speed, TWO_PI_F), 1e-6f);
    previous_phase = phase;
    if (WB::param_morph_elapsed(sb) == 6) {
      HS_EXPECT_GT(WB::live_params(sb).surface_noise.scale,
                   from.surface_noise.scale);
      HS_EXPECT_LT(WB::live_params(sb).surface_noise.scale,
                   to.surface_noise.scale);
      HS_EXPECT_EQ(WB::live_params(sb).warp, from.warp);
    }
  }
  HS_EXPECT_FALSE(WB::param_morph_active(sb));
  HS_EXPECT_TRUE(WB::live_params(sb) == to);

  const auto &presets = WB::presets();
  for (size_t index = 0; index < presets.size(); ++index) {
    HS_EXPECT_TRUE(WB::valid_config(presets[index]));
    const auto &next = presets[(index + 1) % presets.size()];
    HS_EXPECT_TRUE(WB::transition_admitted(presets[index], next));
  }
}

/** @brief Discrete look changes preserve continuous dual-runtime handoff. */
inline void test_shaderball_discrete_transition() {
  using WB = ShaderBallWhiteBox;

  const Color4 from_color(Pixel(10000, 20000, 30000), 0.25f);
  const Color4 to_color(Pixel(40000, 50000, 60000), 0.75f);
  HS_EXPECT_TRUE(::blend_outputs(from_color, to_color, 0.0f).color ==
                 from_color.color);
  HS_EXPECT_TRUE(::blend_outputs(from_color, to_color, 1.0f).color ==
                 to_color.color);
  const Color4 middle = ::blend_outputs(from_color, to_color, 0.5f);
  HS_EXPECT_EQ(middle.alpha, 0.5f);
  HS_EXPECT_EQ(middle.color.r, uint16_t(32500));
  HS_EXPECT_EQ(middle.color.g, uint16_t(42500));
  HS_EXPECT_EQ(middle.color.b, uint16_t(52500));
  const Color4 transparent = ::blend_outputs(
      Color4(Pixel(65535, 0, 0), 0.0f), Color4(Pixel(0, 0, 65535), 0.0f), 0.5f);
  HS_EXPECT_EQ(transparent.alpha, 0.0f);
  HS_EXPECT_EQ(transparent.color, Pixel());

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
    HS_EXPECT_EQ(clear.color, Pixel());
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
    HS_EXPECT_EQ(exact_start.color, expected.color);
    HS_EXPECT_EQ(exact_start.alpha, expected.alpha);
    const auto through_end = WB::through_clear_phase(60, 60);
    const Color4 exact_end = WB::shade_through_clear(view, &valid, through_end);
    HS_EXPECT_EQ(exact_end.color, expected.color);
    HS_EXPECT_EQ(exact_end.alpha, expected.alpha);
  }

  {
    reset_effect_globals();
    WB::SB sb;
    sb.init();
    sb.setAnimationsPaused(true);
    // The noise seed is outside the topology key, so both endpoints compile
    // the same pipeline over different prepared fields.
    WB::RequestedConfig from = WB::presets()[8];
    WB::request_config(sb, from);
    WB::settle_transition(sb);
    const int32_t from_seed = from.params.surface_noise.seed;
    HS_EXPECT_EQ(WB::prepared_surface_noise_seed(sb), from_seed);

    WB::RequestedConfig to = from;
    to.params.surface_noise.seed = from_seed + 91;
    WB::force_transition(sb, to, 60, false);
    HS_EXPECT_TRUE(WB::transition_active(sb));
    for (int frame = 0; frame < 30; ++frame) {
      sb.draw_frame();
      sb.advance_display();
    }
    HS_EXPECT_EQ(WB::transition_elapsed(sb), uint16_t(30));
    HS_EXPECT_EQ(WB::prepared_surface_noise_seed(sb), from_seed);
    const float source_phase_at_clear =
        WB::transition_from_runtime(sb).clocks.source_primary;
    const float destination_phase_at_clear =
        WB::transition_to_runtime(sb).clocks.source_primary;
    sb.draw_frame();
    sb.advance_display();
    HS_EXPECT_EQ(WB::transition_elapsed(sb), uint16_t(31));
    HS_EXPECT_EQ(WB::prepared_surface_noise_seed(sb), from_seed + 91);
    HS_EXPECT_EQ(WB::transition_from_runtime(sb).clocks.source_primary,
                 source_phase_at_clear);
    HS_EXPECT_NEAR(
        WB::transition_to_runtime(sb).clocks.source_primary,
        fmodf(destination_phase_at_clear + to.params.source.speed, TWO_PI_F),
        1e-6f);
  }

  {
    reset_effect_globals();
    WB::SB sb;
    sb.init();
    WB::LookRuntime authored_runtime;
    WB::LookRuntime generated_runtime;
    const WB::WalkDeltas deltas{make_rotation(Y_AXIS, 0.2f),
                                make_rotation(X_AXIS, 0.3f)};
    const WB::RequestedConfig authored = WB::presets()[11];
    WB::advance_runtime(sb, authored_runtime, authored, deltas);
    WB::advance_runtime(sb, generated_runtime, WB::legacy_config(), deltas);
    HS_EXPECT_EQ(authored_runtime.clocks.source_primary,
                 authored.params.source.speed);
    HS_EXPECT_EQ(generated_runtime.clocks.source_primary, 0.05f);
    HS_EXPECT_TRUE(authored_runtime.projection_wander !=
                   generated_runtime.projection_wander);
    HS_EXPECT_TRUE(authored_runtime.outer_wander !=
                   generated_runtime.outer_wander);

    const size_t original_index = WB::preset_index(sb);
    // Both transition endpoints must compile a pipeline, so the destination
    // is an authored topology carrying the generated-look source speed.
    WB::RequestedConfig discrete = WB::presets()[1];
    discrete.params.source.speed = 0.05f;
    WB::force_transition(sb, discrete, 60, true);
    const WB::RequestedConfig captured_source = WB::transition_from_config(sb);
    HS_EXPECT_TRUE(WB::transition_active(sb));
    HS_EXPECT_EQ(WB::transition_mix(sb), 0.0f);
    WB::begin_blend(sb);
    HS_EXPECT_EQ(WB::preset_index(sb), original_index);
    const uint32_t walk_steps = WB::walk_steps(sb);
    const uint32_t generated_steps = WB::generated_palette_steps(sb);
    sb.draw_frame();
    sb.advance_display();
    HS_EXPECT_EQ(WB::walk_steps(sb), walk_steps + 1);
    HS_EXPECT_EQ(WB::generated_palette_steps(sb), generated_steps + 1);
    HS_EXPECT_EQ(WB::transition_from_runtime(sb).clocks.source_primary, 0.245f);
    HS_EXPECT_EQ(WB::transition_to_runtime(sb).clocks.source_primary, 0.05f);

    for (int frame_index = 1; frame_index < 20; ++frame_index) {
      sb.draw_frame();
      sb.advance_display();
    }
    const uint16_t elapsed_before_takeover = WB::transition_elapsed(sb);
    sb.setAnimationsPaused(true);
    const float visible_phase =
        WB::transition_from_runtime(sb).clocks.source_primary;
    const WB::RequestedConfig queued = WB::presets()[2];
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

/** @brief Pause does not stretch an in-flight through-clear transition. */
inline void test_shaderball_pause_does_not_hold_through_clear() {
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

  reset_effect_globals();
  WB::SB sb;
  sb.init();
  const size_t entry_preset = WB::preset_index(sb);
  // Both through-clear endpoints must compile a pipeline.
  const WB::RequestedConfig destination = WB::presets()[1];
  WB::force_transition(sb, destination, 60, true);
  run_frames(sb, 30);
  HS_EXPECT_EQ(WB::transition_elapsed(sb), uint16_t(30));

  sb.setAnimationsPaused(true);
  run_frames(sb, 1);
  HS_EXPECT_EQ(lit_pixels(sb), size_t(0));
  HS_EXPECT_EQ(WB::transition_elapsed(sb), uint16_t(31));
  run_frames(sb, 4);
  HS_EXPECT_TRUE(WB::transition_active(sb));
  HS_EXPECT_EQ(WB::transition_elapsed(sb), uint16_t(35));
  HS_EXPECT_GT(lit_pixels(sb), size_t(0));
  run_frames(sb, 26);
  HS_EXPECT_FALSE(WB::transition_active(sb));
  HS_EXPECT_TRUE(WB::active_config(sb) == destination);
  HS_EXPECT_EQ(WB::preset_index(sb), entry_preset);
}
/** @brief Generated harmonies remain independent and morph-compatible. */
inline void test_shaderball_palette_resources() {
  using WB = ShaderBallWhiteBox;
  auto check_sequence = [](auto make) {
    uint32_t hue = 0;
    GenerativePalette previous;
    make(hue, 0, previous);
    for (uint32_t sequence = 1; sequence <= 8; ++sequence) {
      GenerativePalette next;
      make(hue, sequence, next);
      HS_EXPECT_TRUE(previous.morph_compatible(next));
      previous = next;
    }
    HS_EXPECT_EQ(hue, uint32_t(8) * WB::HUE_STEP);
  };
  check_sequence(WB::make_triadic);
  check_sequence(WB::make_complementary);
  check_sequence(WB::make_analogous);

  reset_effect_globals();
  WB::SB sb;
  sb.init();
  const Pixel triadic = WB::palette_color(sb, WB::PaletteMode::TRIADIC, 0.25f);
  const Pixel complementary =
      WB::palette_color(sb, WB::PaletteMode::COMPLEMENTARY, 0.25f);
  const Pixel analogous =
      WB::palette_color(sb, WB::PaletteMode::ANALOGOUS, 0.25f);
  for (WB::PaletteMode mode :
       {WB::PaletteMode::TRIADIC, WB::PaletteMode::COMPLEMENTARY,
        WB::PaletteMode::ANALOGOUS})
    for (int step = 0; step <= 32; ++step)
      HS_EXPECT_EQ(WB::palette_alpha(sb, mode, step / 32.0f), 1.0f);
  HS_EXPECT_TRUE(triadic.r != complementary.r || triadic.g != complementary.g ||
                 triadic.b != complementary.b);
  HS_EXPECT_TRUE(triadic.r != analogous.r || triadic.g != analogous.g ||
                 triadic.b != analogous.b);
  WB::update_palette_chroma(sb, 1.0f);
  const Pixel saturated =
      WB::palette_color(sb, WB::PaletteMode::TRIADIC, 0.25f);
  HS_EXPECT_TRUE(triadic != saturated);

  WB::RequestedConfig destination = WB::presets()[1];
  destination.params.color.palette_chroma = 0.91f;
  WB::force_transition(sb, destination, 60, false);
  HS_EXPECT_EQ(WB::visible_palette_chroma(sb),
               WB::presets()[0].params.color.palette_chroma);
  WB::set_transition_elapsed(sb, 30);
  HS_EXPECT_EQ(WB::visible_palette_chroma(sb),
               destination.params.color.palette_chroma);
}

/** @brief Palette mapping and brightness operate on separate color stages. */
inline void test_shaderball_brightness_envelopes() {
  using WB = ShaderBallWhiteBox;
  using Mapping = WB::PaletteMapping;
  using Envelope = WB::BrightnessEnvelope;

  HS_EXPECT_EQ(WB::palette_mapping(0.0f, Mapping::CUP, 1.0f), 1.0f);
  HS_EXPECT_EQ(WB::palette_mapping(0.5f, Mapping::CUP, 1.0f), 0.0f);
  HS_EXPECT_EQ(WB::palette_mapping(0.0f, Mapping::BELL, 1.0f), 0.0f);
  HS_EXPECT_EQ(WB::palette_mapping(0.5f, Mapping::BELL, 1.0f), 1.0f);
  HS_EXPECT_EQ(WB::palette_mapping(0.25f, Mapping::LINEAR, 1.0f), 0.25f);
  HS_EXPECT_EQ(WB::palette_mapping(0.25f, Mapping::REVERSE, 1.0f), 0.75f);
  HS_EXPECT_NEAR(WB::palette_mapping(0.25f, Mapping::CUP, 2.0f, 0.1f), 0.2f,
                 1e-6f);

  HS_EXPECT_EQ(WB::brightness_envelope(0.0f, Envelope::NONE, 1.0f), 1.0f);
  HS_EXPECT_EQ(WB::brightness_envelope(0.0f, Envelope::CUP, 1.0f), 1.0f);
  HS_EXPECT_EQ(WB::brightness_envelope(0.5f, Envelope::CUP, 1.0f), 0.0f);
  HS_EXPECT_EQ(WB::brightness_envelope(0.0f, Envelope::BELL, 1.0f), 0.0f);
  HS_EXPECT_EQ(WB::brightness_envelope(0.5f, Envelope::BELL, 1.0f), 1.0f);
  HS_EXPECT_EQ(WB::brightness_envelope(0.25f, Envelope::ASCENDING, 1.0f),
               0.25f);
  HS_EXPECT_EQ(WB::brightness_envelope(1.0f, Envelope::ASCENDING, 1.0f), 1.0f);
  HS_EXPECT_EQ(WB::brightness_envelope(0.25f, Envelope::DESCENDING, 1.0f),
               0.75f);
  HS_EXPECT_EQ(WB::brightness_envelope(0.5f, Envelope::CUP, 0.25f), 0.75f);

  reset_effect_globals();
  WB::SB sb;
  sb.init();
  WB::RequestedConfig config = WB::legacy_config();
  config.slots.hue_shift = WB::HueShiftMode::NONE;
  const WB::MaterialSample sample{
      0.37f, 0.6f, Vector(0.31f, 0.87f, -0.38f).normalized(), 0.0f};
  const Color4 default_color =
      WB::colorize(sample, WB::config_frame(sb, config));
  HS_EXPECT_EQ(default_color.color,
               WB::palette_color(sb, config.slots.palette, sample.value));
  HS_EXPECT_NEAR(default_color.alpha, 0.6f, 1e-6f);

  config.slots.palette_mapping = Mapping::CUP;
  config.slots.hue_shift = WB::HueShiftMode::WARP_DISPLACEMENT;
  config.params.color.hue_shift_amount = 0.25f;
  WB::MaterialSample shifted_sample = sample;
  shifted_sample.value = 0.5f;
  shifted_sample.path_length = 1.0f;
  WB::FrameState frame = WB::config_frame(sb, config);
  const Color4 shifted_without_brightness = WB::colorize(shifted_sample, frame);
  HS_EXPECT_EQ(shifted_without_brightness.color,
               WB::prepared_hue_rotation(frame, 0.0f, 0.25f));

  config.slots.brightness_envelope = Envelope::CUP;
  config.params.color.brightness_depth = 0.5f;
  frame = WB::config_frame(sb, config);
  const Color4 shifted = WB::colorize(shifted_sample, frame);
  const Pixel expected = shifted_without_brightness.color * 0.5f;
  HS_EXPECT_EQ(shifted.color, expected);
  const auto normalized_ratio = [](uint16_t a, uint16_t b) {
    return static_cast<float>(a) / static_cast<float>(b);
  };
  HS_EXPECT_NEAR(normalized_ratio(shifted.color.r, shifted.color.g),
                 normalized_ratio(shifted_without_brightness.color.r,
                                  shifted_without_brightness.color.g),
                 0.01f);
  HS_EXPECT_NEAR(normalized_ratio(shifted.color.b, shifted.color.g),
                 normalized_ratio(shifted_without_brightness.color.b,
                                  shifted_without_brightness.color.g),
                 0.01f);
  config.params.color.value_opacity_low = 0.2f;
  config.params.color.value_opacity_high = 0.8f;
  const Color4 faded =
      WB::colorize(shifted_sample, WB::config_frame(sb, config));
  HS_EXPECT_NEAR(faded.alpha, sample.coverage * 0.5f, 1e-6f);
}

/** @brief Colorize hue modes consume their selected pipeline carrier. */
inline void test_shaderball_hue_shift_modes() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();
  WB::RequestedConfig base = WB::legacy_config();
  base.slots.warp_program.outer.kind = WB::WarpStageKind::NONE;
  base.slots.warp_program.inner.kind = WB::WarpStageKind::NONE;
  base.slots.hue_shift = WB::HueShiftMode::NONE;
  base.params.color = {1.0f, 2.0f, 0.0f};
  HS_EXPECT_TRUE(WB::valid_config(base));

  const WB::MaterialSample sample{
      0.37f, 1.0f, Vector(0.31f, 0.87f, -0.38f).normalized(), 0.0f};
  const Color4 plain = WB::colorize(sample, WB::config_frame(sb, base));
  WB::RequestedConfig noisy = base;
  noisy.slots.hue_shift = WB::HueShiftMode::NOISE;
  noisy.params.color.hue_shift_amount = -1.0f;
  HS_EXPECT_TRUE(WB::valid_config(noisy));
  const WB::FrameState noise_frame = WB::config_frame(sb, noisy);
  HS_EXPECT_TRUE(noise_frame.resources.color_noise != nullptr);
  const Color4 noise_shifted = WB::colorize(sample, noise_frame);
  HS_EXPECT_TRUE(plain.color.r != noise_shifted.color.r ||
                 plain.color.g != noise_shifted.color.g ||
                 plain.color.b != noise_shifted.color.b);

  WB::RequestedConfig displaced = base;
  displaced.slots.hue_shift = WB::HueShiftMode::WARP_DISPLACEMENT;
  displaced.params.color.hue_shift_amount = -1.0f;
  HS_EXPECT_TRUE(WB::valid_config(displaced));
  WB::RequestedConfig max_displaced = displaced;
  max_displaced.params.color.hue_shift_amount = 4.0f;
  HS_EXPECT_TRUE(WB::valid_config(max_displaced));
  const WB::FrameState displacement_frame = WB::config_frame(sb, displaced);
  HS_EXPECT_TRUE(displacement_frame.resources.color_noise == nullptr);
  const Color4 undisplaced = WB::colorize(sample, displacement_frame);
  HS_EXPECT_EQ(undisplaced.color, plain.color);
  WB::MaterialSample warped_sample = sample;
  warped_sample.path_length = 0.25f;
  const Color4 displacement_shifted =
      WB::colorize(warped_sample, displacement_frame);
  HS_EXPECT_TRUE(plain.color.r != displacement_shifted.color.r ||
                 plain.color.g != displacement_shifted.color.g ||
                 plain.color.b != displacement_shifted.color.b);
  noisy.params.color.hue_shift_amount = 4.0f;
  HS_EXPECT_FALSE(WB::valid_config(noisy));
  noisy.params.color.hue_shift_amount = 1.0f;

  const WB::TopologyKey topology = WB::topology_key(base);
  for (WB::PaletteMode palette :
       {WB::PaletteMode::TRIADIC, WB::PaletteMode::COMPLEMENTARY,
        WB::PaletteMode::ANALOGOUS}) {
    noisy.slots.palette = palette;
    HS_EXPECT_TRUE(WB::valid_config(noisy));
    HS_EXPECT_TRUE(WB::topology_key(noisy) == topology);
    const Color4 color = WB::colorize(sample, WB::config_frame(sb, noisy));
    HS_EXPECT_TRUE(std::isfinite(color.alpha));
  }

  for (const char *name : {"Hue Shift Mode", "Hue Shift Amount",
                           "Hue Noise Scale", "Hue Noise Speed"})
    HS_EXPECT_TRUE(sb.getParameters().find(name) != nullptr);
  HS_EXPECT_EQ(sb.getParameters().find("Hue Shift Amount")->max, 1.0f);
  HS_EXPECT_EQ(sb.updateParameter(
                   "Hue Shift Mode",
                   static_cast<float>(WB::HueShiftMode::WARP_DISPLACEMENT)),
               ParamSetResult::APPLIED);
  HS_EXPECT_TRUE(sb.getParameters().find("Hue Shift Amount") != nullptr);
  HS_EXPECT_EQ(sb.getParameters().find("Hue Shift Amount")->max, 4.0f);
  HS_EXPECT_TRUE(sb.getParameters().find("Hue Noise Scale") == nullptr);
  HS_EXPECT_TRUE(sb.getParameters().find("Hue Noise Speed") == nullptr);
  HS_EXPECT_EQ(sb.updateParameter("Hue Shift Amount", 4.0f),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(sb.updateParameter("Hue Shift Mode",
                                  static_cast<float>(WB::HueShiftMode::NOISE)),
               ParamSetResult::APPLIED);
  HS_EXPECT_EQ(sb.getParameters().find("Hue Shift Amount")->max, 1.0f);
  HS_EXPECT_EQ(sb.getParameters().find("Hue Shift Amount")->get_requested(),
               1.0f);

  sb.setAnimationsPaused(true);
  HS_EXPECT_EQ(sb.updateParameter("Hue Noise Speed", 0.0f),
               ParamSetResult::APPLIED);
  WB::settle_transition(sb);
  const float stopped_phase = WB::clocks(sb).hue_noise_phase;
  for (int frame = 0; frame < 120; ++frame) {
    sb.draw_frame();
    sb.advance_display();
  }
  HS_EXPECT_EQ(WB::clocks(sb).hue_noise_phase, stopped_phase);

  HS_EXPECT_EQ(sb.updateParameter("Hue Noise Speed", 0.000002f),
               ParamSetResult::APPLIED);
  for (int frame = 0; frame < 120; ++frame) {
    sb.draw_frame();
    sb.advance_display();
  }
  HS_EXPECT_NE(WB::clocks(sb).hue_noise_phase, stopped_phase);
}
inline void test_shaderball_surface_noise_geometry_and_composition() {
  using WB = ShaderBallWhiteBox;
  reset_effect_globals();
  WB::SB sb;
  sb.init();
  WB::RequestedConfig config = WB::legacy_config();
  config.slots.surface_lens = WB::SurfaceLens::NONE;
  config.slots.surface_noise = WB::SurfaceNoise::DIRECT;
  config.slots.surface_noise_placement = WB::SurfaceNoisePlacement::AFTER_LENS;
  config.params.surface_noise.basis = WB::NoiseBasis::FBM3;
  config.params.surface_noise.scale = 2.0f;
  config.params.surface_noise.strength = 0.3f;
  config.params.surface_noise.rate = 0.0f;
  WB::request_config(sb, config);
  HS_EXPECT_TRUE(sb.getParameters().find("Surface Noise Placement") != nullptr);
  HS_EXPECT_TRUE(sb.getParameters().find("Surface Noise Direction") != nullptr);
  HS_EXPECT_TRUE(sb.getParameters().find("Lens Mix") == nullptr);

  const std::array<Vector, 8> directions = {
      Vector(1.0f, 0.0f, 0.0f),
      Vector(-1.0f, 0.0f, 0.0f),
      Vector(0.0f, 1.0f, 0.0f),
      Vector(0.0f, -1.0f, 0.0f),
      Vector(0.0f, 0.0f, 1.0f),
      Vector(0.0f, 0.0f, -1.0f),
      Vector(1.0f, 2.0f, 3.0f).normalized(),
      Vector(-2.0f, 1.0f, -0.5f).normalized()};
  // Surface noise geometry is a stage contract, so the sweep prepares each
  // config directly instead of asking the roster to compile it.
  WB::FrameState frame = WB::config_frame(sb, config);
  for (const Vector &v : directions) {
    frame.params.surface_noise.direction = 0.0f;
    const Vector a = WB::surface_noise(v, frame);
    frame.params.surface_noise.direction = 1.0f;
    const Vector b = WB::surface_noise(v, frame);
    HS_EXPECT_NEAR(a.length(), 1.0f, 1e-5f);
    HS_EXPECT_LE(fast_acos(hs::clamp(dot(v, a), -1.0f, 1.0f)), 0.30001f);
    HS_EXPECT_NEAR(a.x, b.x, 1e-5f);
    HS_EXPECT_NEAR(a.y, b.y, 1e-5f);
    HS_EXPECT_NEAR(a.z, b.z, 1e-5f);
    frame.params.surface_noise.strength = 0.0f;
    const Vector identity = WB::surface_noise(v, frame);
    HS_EXPECT_EQ(std::memcmp(&identity, &v, sizeof(Vector)), 0);
    frame.params.surface_noise.strength = 0.3f;
  }

  config.slots.surface_noise = WB::SurfaceNoise::CURL;
  config.params.surface_noise.basis = WB::NoiseBasis::RIDGED3;
  for (WB::SurfaceCurlIntegrator integrator :
       {WB::SurfaceCurlIntegrator::EULER, WB::SurfaceCurlIntegrator::MIDPOINT,
        WB::SurfaceCurlIntegrator::MIDPOINT_2X}) {
    config.params.surface_noise.integrator = integrator;
    frame = WB::config_frame(sb, config);
    for (const Vector &v : directions) {
      const Vector positive = WB::surface_noise(v, frame);
      frame.params.surface_noise.strength = -0.3f;
      const Vector negative = WB::surface_noise(v, frame);
      HS_EXPECT_NEAR(positive.length(), 1.0f, 1e-5f);
      HS_EXPECT_NEAR(negative.length(), 1.0f, 1e-5f);
      HS_EXPECT_LE(fast_acos(hs::clamp(dot(v, positive), -1.0f, 1.0f)),
                   0.3001f);
      frame.params.surface_noise.strength = 0.3f;
    }
  }
  frame = WB::config_frame(sb, config);
  frame.params.surface_noise.integrator = WB::SurfaceCurlIntegrator::EULER;
  const Vector circulation = WB::surface_curl_field(directions.back(), frame);
  frame.params.surface_noise.strength = 0.01f;
  const Vector positive = WB::surface_noise(directions.back(), frame);
  frame.params.surface_noise.strength = -0.01f;
  const Vector negative = WB::surface_noise(directions.back(), frame);
  HS_EXPECT_GT(dot(positive - directions.back(), circulation), 0.0f);
  HS_EXPECT_LT(dot(negative - directions.back(), circulation), 0.0f);

  config.slots.surface_noise = WB::SurfaceNoise::DIRECT;
  config.params.surface_noise.basis = WB::NoiseBasis::FBM3;
  config.params.surface_noise.strength = 0.3f;
  config.slots.hue_shift = WB::HueShiftMode::WARP_DISPLACEMENT;
  config.params.color.hue_shift_amount = 1.0f;
  frame = WB::config_frame(sb, config);
  const WB::ProjectedLookup displaced_projected =
      WB::surface_project(directions.back(), frame);
  const WB::PlanarWarpResult unwarped = WB::warp(displaced_projected, frame);
  const WB::MaterialSample displaced_material =
      WB::material(displaced_projected, unwarped, frame);
  HS_EXPECT_GT(displaced_projected.surface_path_length, 0.0f);
  HS_EXPECT_EQ(displaced_material.path_length,
               displaced_projected.surface_path_length);
  WB::MaterialSample undisplaced_material = displaced_material;
  undisplaced_material.path_length = 0.0f;
  HS_EXPECT_TRUE(WB::colorize(displaced_material, frame).color !=
                 WB::colorize(undisplaced_material, frame).color);

  frame = WB::config_frame(sb, config);
  for (WB::Projection projection :
       {WB::Projection::SINUSOIDAL, WB::Projection::STEREOGRAPHIC,
        WB::Projection::GNOMONIC, WB::Projection::BONNE,
        WB::Projection::PEIRCE_QUINCUNCIAL, WB::Projection::AIROCEAN,
        WB::Projection::EQUIRECTANGULAR})
    for (WB::SurfaceLens lens :
         {WB::SurfaceLens::NONE, WB::SurfaceLens::GLITCH,
          WB::SurfaceLens::TWIST, WB::SurfaceLens::KALEIDOSCOPE,
          WB::SurfaceLens::MOBIUS, WB::SurfaceLens::KALEIDOSCOPE_TETRAHEDRAL,
          WB::SurfaceLens::KALEIDOSCOPE_OCTAHEDRAL,
          WB::SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL,
          WB::SurfaceLens::KALEIDOSCOPE_TRIANGULAR_PRISM,
          WB::SurfaceLens::KALEIDOSCOPE_SQUARE_PRISM,
          WB::SurfaceLens::KALEIDOSCOPE_PENTAGONAL_PRISM,
          WB::SurfaceLens::KALEIDOSCOPE_HEXAGONAL_PRISM,
          WB::SurfaceLens::KALEIDOSCOPE_OCTAGONAL_PRISM}) {
      frame.slots.projection = projection;
      frame.slots.surface_lens = lens;
      const WB::ProjectedLookup projected =
          WB::surface_project(directions.back(), frame);
      HS_EXPECT_TRUE(std::isfinite(projected.coords.re));
      HS_EXPECT_TRUE(std::isfinite(projected.coords.im));
    }
}

/** @brief Sphere contours use v_projection and reject every planar warp. */
inline void test_shaderball_noise_contour_domains() {
  using WB = ShaderBallWhiteBox;
  constexpr std::array PROJECTIONS = {
      WB::Projection::SINUSOIDAL,         WB::Projection::STEREOGRAPHIC,
      WB::Projection::GNOMONIC,           WB::Projection::BONNE,
      WB::Projection::PEIRCE_QUINCUNCIAL, WB::Projection::AIROCEAN,
      WB::Projection::EQUIRECTANGULAR};

  WB::RequestedConfig config = WB::legacy_config();
  config.slots.surface_lens = WB::SurfaceLens::NONE;
  config.slots.warp_program.outer.kind = WB::WarpStageKind::NONE;
  config.slots.warp_program.inner.kind = WB::WarpStageKind::NONE;
  config.slots.function = WB::Function::NOISE_CONTOUR;
  for (WB::Projection projection : PROJECTIONS) {
    config.slots.projection = projection;
    const bool projected_admitted =
        projection == WB::Projection::SINUSOIDAL ||
        projection == WB::Projection::STEREOGRAPHIC ||
        projection == WB::Projection::GNOMONIC ||
        projection == WB::Projection::EQUIRECTANGULAR;
    HS_EXPECT_EQ(WB::valid_config(config), projected_admitted);
  }

  config.slots.function = WB::Function::NOISE_CONTOUR_SPHERE;
  for (WB::Projection projection : PROJECTIONS) {
    config.slots.projection = projection;
    HS_EXPECT_TRUE(WB::valid_config(config));
  }
  config.slots.warp_program.outer.kind = WB::WarpStageKind::AFFINE_FRAME;
  HS_EXPECT_FALSE(WB::valid_config(config));
  config.slots.warp_program.outer.kind = WB::WarpStageKind::NONE;
  config.slots.warp_program.inner.kind = WB::WarpStageKind::VORTEX;
  HS_EXPECT_FALSE(WB::valid_config(config));

  {
    reset_effect_globals();
    WB::SB sb;
    sb.init();
    config.slots.projection = WB::Projection::STEREOGRAPHIC;
    config.slots.warp_program.inner.kind = WB::WarpStageKind::NONE;
    WB::FrameState frame = WB::config_frame(sb, config);
    frame.transforms.projection_conj = make_rotation(Y_AXIS, -0.37f);
    const Vector point = Vector(0.4f, -0.8f, 0.3f).normalized();
    const Vector expected = rotate(point, frame.transforms.projection_conj);
    float sphere_reference = 0.0f;
    for (size_t index = 0; index < PROJECTIONS.size(); ++index) {
      frame.slots.projection = PROJECTIONS[index];
      const WB::ProjectedLookup projected = WB::surface_project(point, frame);
      HS_EXPECT_NEAR(projected.sphere.x, expected.x, 1e-6f);
      HS_EXPECT_NEAR(projected.sphere.y, expected.y, 1e-6f);
      HS_EXPECT_NEAR(projected.sphere.z, expected.z, 1e-6f);
      const float sample = WB::source(projected, frame);
      if (index == 0)
        sphere_reference = sample;
      else
        HS_EXPECT_NEAR(sample, sphere_reference, 1e-6f);
    }

    frame.slots.function = WB::Function::NOISE_CONTOUR;
    frame.slots.projection = WB::Projection::SINUSOIDAL;
    const float sinusoidal =
        WB::source(WB::surface_project(point, frame), frame);
    frame.slots.projection = WB::Projection::STEREOGRAPHIC;
    const float stereographic =
        WB::source(WB::surface_project(point, frame), frame);
    HS_EXPECT_GT(std::fabs(sinusoidal - stereographic), 1e-5f);
  }

  reset_effect_globals();
  WB::SB gui;
  gui.init();
  HS_EXPECT_EQ(
      gui.updateParameter("Planar Warp 1",
                          static_cast<float>(WB::WarpStageKind::MIRROR_TILE)),
      ParamSetResult::APPLIED);
  gui.draw_frame();
  gui.advance_display();
  WB::settle_transition(gui);
  HS_EXPECT_EQ(
      gui.updateParameter(
          "Function", static_cast<float>(WB::Function::NOISE_CONTOUR_SPHERE)),
      ParamSetResult::APPLIED);
  gui.draw_frame();
  gui.advance_display();
  HS_EXPECT_NE(WB::active_slots(gui).function,
               WB::Function::NOISE_CONTOUR_SPHERE);
  const char *warning = WB::parameter_warning(gui, "Function");
  HS_EXPECT_TRUE(warning != nullptr);
  HS_EXPECT_TRUE(std::strstr(warning, "Planar Warp 1 Mirror Tile") != nullptr);
  HS_EXPECT_EQ(gui.updateParameter("Planar Warp 1",
                                   static_cast<float>(WB::WarpStageKind::NONE)),
               ParamSetResult::APPLIED);
  gui.draw_frame();
  gui.advance_display();
  WB::settle_transition(gui);
  HS_EXPECT_EQ(WB::requested_config(gui).slots.function,
               WB::Function::NOISE_CONTOUR_SPHERE);
  HS_EXPECT_EQ(WB::active_slots(gui).function,
               WB::Function::NOISE_CONTOUR_SPHERE);
  warning = WB::parameter_warning(gui, "Function");
  HS_EXPECT_TRUE(warning == nullptr);
  HS_EXPECT_EQ(WB::active_slots(gui).warp_program.outer.kind,
               WB::WarpStageKind::NONE);
  HS_EXPECT_EQ(WB::active_pipeline(gui), WB::InversePipelineId::NONE);
}

/** @brief Module entry point for ShaderBall contract tests. */
inline int run_shaderball_tests() {
  ModuleFixture fixture("shaderball");
  test_shaderball_operator_catalog_census();
  test_shaderball_inverse_pipeline_manifest();
  test_shaderball_inverse_program_equivalence();
  test_shaderball_full_config_snapshot();
  test_shaderball_surface_noise_range_rebind();
  test_shaderball_incompatible_config_snapshot();
  test_shaderball_selector_storage();
  test_shaderball_clocks_wrapped();
  test_shaderball_pause_semantics();
  test_shaderball_manual_edit_timing();
  test_shaderball_pipeline_contract();
  test_shaderball_legacy_spatial_slots();
  test_shaderball_kaleidoscope_reflection_fold();
  test_shaderball_polyhedral_kaleidoscopes();
  test_shaderball_equirectangular_projection();
  test_shaderball_flush_edge_fade();
  test_shaderball_legacy_sources();
  test_shaderball_coupled_source();
  test_shaderball_hue_rotate_lut_gamut();
  test_shaderball_prepared_hue_rotation();
  test_shaderball_prepared_hue_noise();
  test_shaderball_prepared_hue_noise_color();
  test_shaderball_preset_bank();
  test_shaderball_staggered_param_morph();
  test_shaderball_config_admission();
  test_shaderball_deterministic_gui_edits();
  test_shaderball_mode_specific_parameter_warnings();
  test_shaderball_dodecahedral_lattice_edit();
  test_shaderball_polar_gui_repair();
  test_shaderball_structural_admission();
  test_shaderball_strict_seam_admission();
  test_shaderball_additive_delta_precision();
  test_shaderball_profile_presets();
  test_shaderball_manual_preset_navigation();
  test_shaderball_preset_gui_transition();
  test_shaderball_parameter_capacity();
  test_shaderball_gui_catalog();
  test_shaderball_lens_domain_ranges();
  test_shaderball_projection_catalog();
  test_shaderball_fast_peirce_square();
  test_shaderball_projection_and_admission_contracts();
  test_shaderball_planar_warp_animation();
  test_shaderball_dynamic_backend_admission();
  test_shaderball_kernel_catalog();
  test_shaderball_stable_preset_transition();
  test_shaderball_discrete_transition();
  test_shaderball_pause_does_not_hold_through_clear();
  test_shaderball_palette_resources();
  test_shaderball_brightness_envelopes();
  test_shaderball_hue_shift_modes();
  test_shaderball_surface_noise_geometry_and_composition();
  test_shaderball_noise_contour_domains();
  return fixture.result();
}

} // namespace shaderball_tests
} // namespace hs_test
