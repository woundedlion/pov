/*
 * Required Notice: Copyright 2025 Gabriel Levy. All rights reserved.
 * LICENSE: ALL RIGHTS RESERVED. No redistribution or use without explicit
 * permission.
 */
#pragma once

#include <cstdarg>
#include <cstdio>
#include <span>
#include <tuple>
#include <type_traits>
#include <utility>

/**
 * @file ShaderBall.h
 * @brief Typed pullback sphere shader with composable projection and material stages.
 */

#include "core/color/effect_palette_recipes.h"
#include "core/engine/engine.h"
#include "core/math/lenses.h"
#include "core/math/noise_field.h"
#include "core/math/projections.h"
#include "core/math/stereographic.h"
#include "core/render/pullback.h"

namespace hs_test {
namespace shaderball_tests {
struct ShaderBallWhiteBox;
} // namespace shaderball_tests
} // namespace hs_test

#define HS_SHADERBALL_CONFIG_FIELDS(X)                                         \
  X(SLOTS_FUNCTION, slots.function)                                            \
  X(SLOTS_PROJECTION, slots.projection)                                        \
  X(SLOTS_PROJECTION_FRAME, slots.projection_frame)                            \
  X(SLOTS_SURFACE_LENS, slots.surface_lens)                                    \
  X(SLOTS_WARP_OUTER_KIND, slots.warp_program.outer.kind)                      \
  X(SLOTS_WARP_OUTER_BASIS, slots.warp_program.outer.basis)                    \
  X(SLOTS_WARP_OUTER_ENVELOPE, slots.warp_program.outer.envelope)              \
  X(SLOTS_WARP_OUTER_POLAR_MODE, slots.warp_program.outer.polar_mode)          \
  X(SLOTS_WARP_OUTER_CURL_INTEGRATOR,                                          \
    slots.warp_program.outer.curl_integrator)                                  \
  X(SLOTS_WARP_OUTER_POLAR_HARMONIC, slots.warp_program.outer.polar_harmonic)  \
  X(SLOTS_WARP_OUTER_SEED, slots.warp_program.outer.seed)                      \
  X(SLOTS_WARP_OUTER_RESOURCE_ID, slots.warp_program.outer.resource_id)        \
  X(SLOTS_WARP_INNER_KIND, slots.warp_program.inner.kind)                      \
  X(SLOTS_WARP_INNER_BASIS, slots.warp_program.inner.basis)                    \
  X(SLOTS_WARP_INNER_ENVELOPE, slots.warp_program.inner.envelope)              \
  X(SLOTS_WARP_INNER_POLAR_MODE, slots.warp_program.inner.polar_mode)          \
  X(SLOTS_WARP_INNER_CURL_INTEGRATOR,                                          \
    slots.warp_program.inner.curl_integrator)                                  \
  X(SLOTS_WARP_INNER_POLAR_HARMONIC, slots.warp_program.inner.polar_harmonic)  \
  X(SLOTS_WARP_INNER_SEED, slots.warp_program.inner.seed)                      \
  X(SLOTS_WARP_INNER_RESOURCE_ID, slots.warp_program.inner.resource_id)        \
  X(SLOTS_SIGNAL_WEIGHT, slots.signal_weight)                                  \
  X(SLOTS_VALUE_TRANSFER, slots.value_transfer)                                \
  X(SLOTS_COVERAGE, slots.coverage)                                            \
  X(SLOTS_PALETTE, slots.palette)                                              \
  X(SLOTS_PALETTE_MAPPING, slots.palette_mapping)                              \
  X(SLOTS_BRIGHTNESS_ENVELOPE, slots.brightness_envelope)                      \
  X(SLOTS_HUE_SHIFT, slots.hue_shift)                                          \
  X(SLOTS_PEIRCE_LAYOUT, slots.peirce_layout)                                  \
  X(SLOTS_AIROCEAN_LAYOUT, slots.airocean_layout)                              \
  X(SLOTS_BONNE_HEMISPHERE, slots.bonne_hemisphere)                            \
  X(SLOTS_GNOMONIC_HEMISPHERE, slots.gnomonic_hemisphere)                      \
  X(SLOTS_SURFACE_NOISE, slots.surface_noise)                                  \
  X(SLOTS_SURFACE_NOISE_PLACEMENT, slots.surface_noise_placement)              \
  X(SOURCE_PATTERN_FREQ, params.source.pattern_freq)                           \
  X(SOURCE_SPEED, params.source.speed)                                         \
  X(SOURCE_COMPLEXITY, params.source.complexity)                               \
  X(SOURCE_PATTERN_MIX, params.source.pattern_mix)                             \
  X(SOURCE_SECONDARY_RATE, params.source.secondary_rate)                       \
  X(SOURCE_ANGLE_RATE, params.source.angle_rate)                               \
  X(SOURCE_NOISE_SCALE, params.source.noise_scale)                             \
  X(SOURCE_NOISE_CONTRAST, params.source.noise_contrast)                       \
  X(SOURCE_NOISE_RATE, params.source.noise_time_rate)                          \
  X(SOURCE_LATTICE_CELL_SCALE, params.source.lattice_cell_scale)               \
  X(SOURCE_LATTICE_SHAPE_BLEND, params.source.lattice_shape_blend)             \
  X(SOURCE_LATTICE_SOFTNESS, params.source.lattice_softness)                   \
  X(SOURCE_LATTICE_RADIUS, params.source.lattice_radius)                       \
  X(SOURCE_NOISE_BASIS, params.source.noise_basis)                             \
  X(SOURCE_NOISE_SEED, params.source.noise_seed)                               \
  X(SOURCE_NOISE_RESOURCE_ID, params.source.noise_resource_id)                 \
  X(WARP_OUTER_SCALE, params.warp.outer.scale)                                 \
  X(WARP_OUTER_STRENGTH, params.warp.outer.strength)                           \
  X(WARP_OUTER_SPEED, params.warp.outer.speed)                                 \
  X(WARP_OUTER_TRANSLATION_X, params.warp.outer.translation_x)                 \
  X(WARP_OUTER_TRANSLATION_Y, params.warp.outer.translation_y)                 \
  X(WARP_OUTER_ROTATION, params.warp.outer.rotation)                           \
  X(WARP_OUTER_SCALE_X, params.warp.outer.scale_x)                             \
  X(WARP_OUTER_SCALE_Y, params.warp.outer.scale_y)                             \
  X(WARP_OUTER_SHEAR, params.warp.outer.shear)                                 \
  X(WARP_OUTER_FREQUENCY, params.warp.outer.frequency)                         \
  X(WARP_OUTER_FIELD_ANGLE, params.warp.outer.field_angle)                     \
  X(WARP_OUTER_CENTER_X, params.warp.outer.center_x)                           \
  X(WARP_OUTER_CENTER_Y, params.warp.outer.center_y)                           \
  X(WARP_OUTER_RADIUS, params.warp.outer.radius)                               \
  X(WARP_OUTER_TURNS, params.warp.outer.turns)                                 \
  X(WARP_OUTER_CENTER_ORBIT_RADIUS, params.warp.outer.center_orbit_radius)     \
  X(WARP_OUTER_VECTOR_ANGLE, params.warp.outer.vector_angle)                   \
  X(WARP_OUTER_CELL_X, params.warp.outer.cell_x)                               \
  X(WARP_OUTER_CELL_Y, params.warp.outer.cell_y)                               \
  X(WARP_OUTER_OFFSET_X, params.warp.outer.offset_x)                           \
  X(WARP_OUTER_OFFSET_Y, params.warp.outer.offset_y)                           \
  X(WARP_OUTER_RADIAL_SCALE, params.warp.outer.radial_scale)                   \
  X(WARP_OUTER_RADIAL_PHASE, params.warp.outer.radial_phase)                   \
  X(WARP_OUTER_ANGULAR_PHASE, params.warp.outer.angular_phase)                 \
  X(WARP_OUTER_EDGE_WIDTH, params.warp.outer.edge_width)                       \
  X(WARP_INNER_SCALE, params.warp.inner.scale)                                 \
  X(WARP_INNER_STRENGTH, params.warp.inner.strength)                           \
  X(WARP_INNER_SPEED, params.warp.inner.speed)                                 \
  X(WARP_INNER_TRANSLATION_X, params.warp.inner.translation_x)                 \
  X(WARP_INNER_TRANSLATION_Y, params.warp.inner.translation_y)                 \
  X(WARP_INNER_ROTATION, params.warp.inner.rotation)                           \
  X(WARP_INNER_SCALE_X, params.warp.inner.scale_x)                             \
  X(WARP_INNER_SCALE_Y, params.warp.inner.scale_y)                             \
  X(WARP_INNER_SHEAR, params.warp.inner.shear)                                 \
  X(WARP_INNER_FREQUENCY, params.warp.inner.frequency)                         \
  X(WARP_INNER_FIELD_ANGLE, params.warp.inner.field_angle)                     \
  X(WARP_INNER_CENTER_X, params.warp.inner.center_x)                           \
  X(WARP_INNER_CENTER_Y, params.warp.inner.center_y)                           \
  X(WARP_INNER_RADIUS, params.warp.inner.radius)                               \
  X(WARP_INNER_TURNS, params.warp.inner.turns)                                 \
  X(WARP_INNER_CENTER_ORBIT_RADIUS, params.warp.inner.center_orbit_radius)     \
  X(WARP_INNER_VECTOR_ANGLE, params.warp.inner.vector_angle)                   \
  X(WARP_INNER_CELL_X, params.warp.inner.cell_x)                               \
  X(WARP_INNER_CELL_Y, params.warp.inner.cell_y)                               \
  X(WARP_INNER_OFFSET_X, params.warp.inner.offset_x)                           \
  X(WARP_INNER_OFFSET_Y, params.warp.inner.offset_y)                           \
  X(WARP_INNER_RADIAL_SCALE, params.warp.inner.radial_scale)                   \
  X(WARP_INNER_RADIAL_PHASE, params.warp.inner.radial_phase)                   \
  X(WARP_INNER_ANGULAR_PHASE, params.warp.inner.angular_phase)                 \
  X(WARP_INNER_EDGE_WIDTH, params.warp.inner.edge_width)                       \
  X(PROJECTION_POLE_FADE, params.projection.pole_fade)                         \
  X(PROJECTION_SPIN_RATE, params.projection.spin_rate)                         \
  X(PROJECTION_WANDER, params.projection.wander)                               \
  X(PROJECTION_CENTRAL_MERIDIAN, params.projection.central_meridian)           \
  X(PROJECTION_COORDINATE_SCALE, params.projection.coordinate_scale)           \
  X(PROJECTION_BONNE_STANDARD_PARALLEL,                                        \
    params.projection.bonne_standard_parallel)                                 \
  X(PROJECTION_LAYOUT_SCROLL, params.projection.layout_scroll)                 \
  X(LENS_MIX, params.surface_lens.mix)                                         \
  X(LENS_MOBIUS_A_RE, params.surface_lens.mobius.a.re)                         \
  X(LENS_MOBIUS_A_IM, params.surface_lens.mobius.a.im)                         \
  X(LENS_MOBIUS_B_RE, params.surface_lens.mobius.b.re)                         \
  X(LENS_MOBIUS_B_IM, params.surface_lens.mobius.b.im)                         \
  X(LENS_MOBIUS_C_RE, params.surface_lens.mobius.c.re)                         \
  X(LENS_MOBIUS_C_IM, params.surface_lens.mobius.c.im)                         \
  X(LENS_MOBIUS_D_RE, params.surface_lens.mobius.d.re)                         \
  X(LENS_MOBIUS_D_IM, params.surface_lens.mobius.d.im)                         \
  X(VALUE_ISO_LEVEL, params.value.iso_level)                                   \
  X(VALUE_ISO_WIDTH, params.value.iso_width)                                   \
  X(VALUE_BAND_COUNT, params.value.band_count)                                 \
  X(VALUE_BAND_PHASE, params.value.band_phase)                                 \
  X(VALUE_CUTOUT_THRESHOLD, params.value.cutout_threshold)                     \
  X(VALUE_CUTOUT_SOFTNESS, params.value.cutout_softness)                       \
  X(VALUE_EDGE_WIDTH, params.value.edge_width)                                 \
  X(COLOR_HUE_SHIFT_AMOUNT, params.color.hue_shift_amount)                     \
  X(COLOR_HUE_NOISE_SCALE, params.color.hue_noise_scale)                       \
  X(COLOR_HUE_NOISE_SPEED, params.color.hue_noise_speed)                       \
  X(COLOR_PALETTE_CHROMA, params.color.palette_chroma)                         \
  X(COLOR_MAPPING_FREQUENCY, params.color.mapping_frequency)                   \
  X(COLOR_MAPPING_PHASE, params.color.mapping_phase)                           \
  X(COLOR_PHASE_OSCILLATION_DEPTH, params.color.phase_oscillation_depth)       \
  X(COLOR_PHASE_OSCILLATION_SPEED, params.color.phase_oscillation_speed)       \
  X(COLOR_BRIGHTNESS_DEPTH, params.color.brightness_depth)                     \
  X(COLOR_VALUE_OPACITY_LOW, params.color.value_opacity_low)                   \
  X(COLOR_VALUE_OPACITY_HIGH, params.color.value_opacity_high)                 \
  X(CAMERA_WANDER, params.outer_camera.wander)                                 \
  X(SURFACE_NOISE_BASIS, params.surface_noise.basis)                           \
  X(SURFACE_NOISE_INTEGRATOR, params.surface_noise.integrator)                 \
  X(SURFACE_NOISE_SEED, params.surface_noise.seed)                             \
  X(SURFACE_NOISE_SCALE, params.surface_noise.scale)                           \
  X(SURFACE_NOISE_STRENGTH, params.surface_noise.strength)                     \
  X(SURFACE_NOISE_RATE, params.surface_noise.rate)                             \
  X(SURFACE_NOISE_DIRECTION, params.surface_noise.direction)

/**
 * @brief Slot-based sphere shader with an immutable per-frame pullback state.
 * @tparam W Canvas width in pixels.
 * @tparam H Canvas height in pixels.
 */
template <int W, int H> class ShaderBall : public Effect {
private:
  struct FrameState;
  struct LookRuntime;
  struct WalkDeltas;
  using ShadeFunction = Color4 (*)(const Vector &, const FrameState &);
  struct FrameShader {
    using ShadeFunction = ShaderBall::ShadeFunction;

    const FrameState *frame;
    float alpha;
    ShadeFunction shade_function;

    HS_FLASH_MEMBER Color4 operator()(const Vector &view) const {
      Color4 color = shade_function(view, *frame);
      color.alpha *= alpha;
      return color;
    }
  };

public:
  struct Config;
  static constexpr int GAMUT_ANGLE_STEPS = GAMUT_LUT_ANGLE_STEPS;
  static constexpr int GAMUT_L_STEPS = GAMUT_LUT_L_STEPS;

  HS_COLD_MEMBER ShaderBall() : Effect(W, H, {.strobe = true}) {}

protected:
  HS_COLD_MEMBER void
  set_fixed_preset_view(std::span<const uint8_t> source_indices) {
    assert(!source_indices.empty());
#ifndef NDEBUG
    for (uint8_t index : source_indices)
      assert(index < PRESETS.size());
#endif
    preset_view = source_indices;
    fixed_topology = true;
  }

  HS_COLD_MEMBER void start_circular_mobius_animation(float scale,
                                                      int duration) {
    timeline.add_pausable(
        0,
        Animation::MobiusWarpCircular(blend.params.surface_lens.mobius, scale,
                                      duration, true),
        &anims_paused);
  }

public:
  /** @brief Initializes slots, clocks, palette resources, and choreography. */
  HS_COLD_MEMBER void init() override {
    configure_presets(preset_count_for_view());
#if HS_ENABLE_PARAM_GUI_BRIDGE
    set_parameter_updated_hook(&ShaderBall::dispatch_parameter_updated);
#endif
    state = persistent_arena.make<StateBundle>();
    use_parameter_storage(persistent_arena.allocate_n<ParamDef>(PARAM_CAPACITY),
                          PARAM_CAPACITY);
    const Preset &initial = preset_for_view(0);
    active_slots = initial.config.slots;
    active_pipeline = initial.pipeline;
    blend.params = initial.config.params;
    blend.palette_mapping =
        palette_mapping_weights(initial.config.slots.palette_mapping);
#if HS_ENABLE_PARAM_GUI_BRIDGE
    display_config = initial.config;
#endif
    requested_config = initial.config;
    published_config = initial.config;
#if HS_ENABLE_PARAM_GUI_BRIDGE
    accepted_config = initial.config;
#endif
    prepare_resource_union(initial.config, initial.config);

    rebind_parameters();

    timeline.add(0, Animation::RandomWalk<W>(projection_walk, UP,
                                             state->projection_walk_noise));
    timeline.add(
        0, Animation::RandomWalk<W>(outer_walk, UP, state->outer_walk_noise));

    init_gamut_lut(persistent_arena, GAMUT_ANGLE_STEPS, GAMUT_L_STEPS);
    triadic_palette_cycler.init_generated(
        persistent_arena, next_triadic_palette, this, PALETTE_DWELL_FRAMES,
        PALETTE_FADE_FRAMES, ease_in_out_sin);
    complementary_palette_cycler.init_generated(
        persistent_arena, next_complementary_palette, this,
        PALETTE_DWELL_FRAMES, PALETTE_FADE_FRAMES, ease_in_out_sin);
    analogous_palette_cycler.init_generated(
        persistent_arena, next_analogous_palette, this, PALETTE_DWELL_FRAMES,
        PALETTE_FADE_FRAMES, ease_in_out_sin);
    update_palette_chroma(
        preset_for_view(0).config.params.color.palette_chroma);

    enter_preset();
  }

  /** @brief Advances mutable state, snapshots it, and renders one frame. */
  HS_FLASH_MEMBER void draw_frame() override {
    Canvas canvas(*this);
    {
      HS_PROFILE(sb_timeline_step);
      timeline.step(canvas);
    }
    advance_preset_choreography();

    apply_requested_config();
    prepare_param_morph();
    state->render_config.slots = active_slots;
    state->render_config.params = blend.params;
    const WalkDeltas walk_deltas = sample_walk_deltas();
    if (state->transition.active) {
      if (state->transition.elapsed < state->transition.duration / 2)
        advance_runtime(state->transition.from_runtime,
                        state->transition.from_config, walk_deltas);
      advance_runtime(state->transition.to_runtime, state->transition.to_config,
                      walk_deltas);
    } else {
      advance_runtime(runtime, state->render_config, walk_deltas);
    }
    update_palette_chroma(visible_palette_chroma());
    step_generated_palettes(visible_palette_mode());
#if HS_ENABLE_TEST_HOOKS
    ++generated_palette_step_count;
#endif

    if (state->transition.active) {
      draw_through_clear_transition(canvas);
    } else {
      PreparedEndpoint prepared;
      HS_CHECK(prepare_endpoint(state->render_config, runtime, 1.0f,
                                active_pipeline, prepared),
               "ShaderBall active endpoint has no renderer");
      draw_endpoint(canvas, prepared);
    }
    finish_transitions();
    publish_live_config();
  }

protected:
  template <int Program> HS_FLASH_MEMBER void draw_fixed_program_frame() {
    Canvas canvas(*this);
    {
      HS_PROFILE(sb_timeline_step);
      timeline.step(canvas);
    }
    advance_preset_choreography();
    apply_requested_config();
    prepare_param_morph();
    state->render_config.slots = active_slots;
    state->render_config.params = blend.params;
    const WalkDeltas walk_deltas = sample_walk_deltas();
    if (state->transition.active) {
      if (state->transition.elapsed < state->transition.duration / 2)
        advance_runtime(state->transition.from_runtime,
                        state->transition.from_config, walk_deltas);
      advance_runtime(state->transition.to_runtime, state->transition.to_config,
                      walk_deltas);
    } else {
      advance_runtime(runtime, state->render_config, walk_deltas);
    }
    update_palette_chroma(visible_palette_chroma());
    step_generated_palettes(visible_palette_mode());
#if HS_ENABLE_TEST_HOOKS
    ++generated_palette_step_count;
#endif

    if (state->transition.active) {
      const ThroughClearPhase phase = through_clear_phase(
          state->transition.elapsed, state->transition.duration);
      if (!phase.clear) {
        const Config &config = phase.from_endpoint
                                   ? state->transition.from_config
                                   : state->transition.to_config;
        const LookRuntime &look = phase.from_endpoint
                                      ? state->transition.from_runtime
                                      : state->transition.to_runtime;
        draw_fixed_endpoint<Program>(canvas, config, look, phase.alpha);
      }
    } else {
      draw_fixed_endpoint<Program>(canvas, state->render_config, runtime, 1.0f);
    }
    finish_transitions();
    publish_live_config();
  }

private:
  template <typename Pipeline>
  HS_FLASH_MEMBER void
  draw_typed_fixed_endpoint(Canvas &canvas, const Config &config,
                            const LookRuntime &look, float alpha) {
    FrameState &frame = state->frame;
    prepare_frame(config, look, frame);
    HS_CHECK(pipeline_resources_ready(frame),
             "fixed Shader endpoint resources are unavailable");
    draw_fixed_prepared(canvas, frame, alpha, &Pipeline::shade);
  }

  HS_FLASH_MEMBER void draw_fixed_prepared(Canvas &canvas, FrameState &frame,
                                           float alpha, ShadeFunction shade) {
    FrameShader shader{&frame, alpha, shade};
    HS_PROFILE(sb_shader_draw);
    Scan::Shader::draw<W, H, 1>(canvas, shader);
  }

  template <int Program>
  HS_FLASH_MEMBER void draw_fixed_endpoint(Canvas &canvas, const Config &config,
                                           const LookRuntime &look,
                                           float alpha) {
    if constexpr (Program == 0)
      draw_typed_fixed_endpoint<GlitchNoiseGridWaveShearPipeline>(
          canvas, config, look, alpha);
    else if constexpr (Program == 1)
      draw_typed_fixed_endpoint<KaleidoscopeTwinWaveInnerMirrorPipeline>(
          canvas, config, look, alpha);
    else if constexpr (Program == 2)
      draw_typed_fixed_endpoint<GnomonicKaleidoscopeGridMirrorPipeline>(
          canvas, config, look, alpha);
    else if constexpr (Program == 3)
      draw_typed_fixed_endpoint<GnomonicGlitchGridMirrorPipeline>(
          canvas, config, look, alpha);
    else if constexpr (Program == 4)
      draw_typed_fixed_endpoint<PeirceDodecahedralGridPipeline>(canvas, config,
                                                                look, alpha);
    else if constexpr (Program == 5)
      draw_typed_fixed_endpoint<GnomonicDodecahedralGridWaveMirrorPipeline>(
          canvas, config, look, alpha);
    else if constexpr (Program == 6)
      draw_typed_fixed_endpoint<GnomonicAffineLatticeContourPipeline>(
          canvas, config, look, alpha);
    else if constexpr (Program == 8)
      draw_typed_fixed_endpoint<StereographicPrismPolarWaveLatticePipeline>(
          canvas, config, look, alpha);
    else if constexpr (Program == 9)
      draw_typed_fixed_endpoint<GnomonicDodecahedralGridVectorMirrorPipeline>(
          canvas, config, look, alpha);
    else if constexpr (Program == 11)
      draw_typed_fixed_endpoint<
          StereographicHexagonalPrismTwinWaveInnerMirrorPipeline>(
          canvas, config, look, alpha);
    else if constexpr (Program == 12)
      draw_typed_fixed_endpoint<
          EquirectangularDodecahedralGridInnerMirrorPipeline>(canvas, config,
                                                              look, alpha);
    else if constexpr (Program == 13)
      draw_typed_fixed_endpoint<StereographicGlitchGridMirrorPipeline>(
          canvas, config, look, alpha);
    else if constexpr (Program == 14)
      draw_typed_fixed_endpoint<StereographicMobiusTwinWaveInnerMirrorPipeline>(
          canvas, config, look, alpha);
    else
      static_assert(Program >= 0 && Program <= 14,
                    "unsupported promoted Shader program");
  }

public:
#if HS_ENABLE_EFFECT_CONTROL_API
  void profile_select_preset(size_t index) {
    HS_CHECK(index < preset_count_for_view(),
             "ShaderBall profile preset index out of range");
    HS_CHECK(selectPreset(index), "ShaderBall profile preset selection failed");
    hs::log("Profile preset: %u/%u", static_cast<unsigned>(index),
            static_cast<unsigned>(preset_count_for_view()));
  }
#endif

private:
  friend struct ::hs_test::shaderball_tests::ShaderBallWhiteBox;
  using NoiseBasis = ::NoiseBasis;

  HS_COLD_MEMBER bool apply_preset(const PresetChange &change) override {
    const size_t index = change.to;
    if (change.origin == PresetChangeOrigin::AUTOMATIC) {
      const Choreo choreo = preset_choreo(change.from);
      const Preset &to = preset_for_view(index);
      if (!try_apply_config(to.config, choreo.blend_frames, choreo.staggered,
                            true))
        return false;
      requested_config = to.config;
      published_config = to.config;
#if HS_ENABLE_PARAM_GUI_BRIDGE
      accepted_config = to.config;
      pending_edit_count = 0;
#endif
      rebind_parameters();
      return true;
    }

    state->param_morph.active = false;
    state->transition.active = false;
    const SelectedConfig &selected = preset_for_view(index);
    active_slots = selected.config.slots;
    active_pipeline = selected.pipeline;
    blend.params = selected.config.params;
    blend.palette_mapping =
        palette_mapping_weights(selected.config.slots.palette_mapping);
#if HS_ENABLE_PARAM_GUI_BRIDGE
    display_config = selected.config;
#endif
    requested_config = selected.config;
    published_config = selected.config;
#if HS_ENABLE_PARAM_GUI_BRIDGE
    accepted_config = selected.config;
    pending_edit_count = 0;
#endif
    runtime = {};
    HS_CHECK(prepare_resource_union(selected.config, selected.config),
             "ShaderBall preset resources exceed capacity");
    rebind_parameters();
    return true;
  }

  HS_COLD_MEMBER void preset_changed(const PresetChange &) override {
    if (!state->param_morph.active && !state->transition.active)
      enter_preset();
  }

public:
  enum class Function : uint8_t {
    TWIN_WAVE,
    RINGS,
    SPIRAL,
    GRID,
    NOISE_CONTOUR,
    PRIMITIVE_LATTICE,
    NOISE_CONTOUR_SPHERE
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
    KALEIDOSCOPE_TETRAHEDRAL,
    KALEIDOSCOPE_OCTAHEDRAL,
    KALEIDOSCOPE_DODECAHEDRAL,
    KALEIDOSCOPE_TRIANGULAR_PRISM,
    KALEIDOSCOPE_SQUARE_PRISM,
    KALEIDOSCOPE_PENTAGONAL_PRISM,
    KALEIDOSCOPE_HEXAGONAL_PRISM,
    KALEIDOSCOPE_OCTAGONAL_PRISM,
    TANGENT_NOISE = 255
  };
  enum class WarpEnvelope : uint8_t { FLAT, PROJECTION_WEIGHT, EDGE_FADE };
  enum class PolarMode : uint8_t { LINEAR, LOGARITHMIC };
  enum class CurlIntegrator : uint8_t { EULER_1, MIDPOINT_2, MIDPOINT_4 };
  enum class SurfaceCurlIntegrator : uint8_t { EULER, MIDPOINT, MIDPOINT_2X };
  enum class SurfaceNoise : uint8_t { NONE, DIRECT, CURL };
  enum class SurfaceNoisePlacement : uint8_t { BEFORE_LENS, AFTER_LENS };
  enum class WarpStageKind : uint8_t {
    NONE,
    AFFINE_FRAME,
    WAVE_SHEAR,
    VORTEX,
    VECTOR_NOISE,
    CURL_FLOW,
    MIRROR_TILE,
    POLAR_CHART,
    LEGACY_STEREO_NOISE = 255
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

    constexpr bool operator==(const WarpStageSpec &) const = default;
  };
  struct WarpProgram {
    WarpStageSpec outer;
    WarpStageSpec inner;

    constexpr bool operator==(const WarpProgram &) const = default;
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
  enum class PaletteMode : uint8_t { TRIADIC, COMPLEMENTARY, ANALOGOUS };
  enum class PaletteMapping : uint8_t { CUP, BELL, LINEAR, REVERSE };
  using PaletteMappingWeights = Pullback::Color::PaletteMappingWeights;
  enum class BrightnessEnvelope : uint8_t {
    NONE,
    CUP,
    BELL,
    ASCENDING,
    DESCENDING
  };
  enum class HueShiftMode : uint8_t { NONE, NOISE, WARP_DISPLACEMENT };

  struct Slots {
    Function function;
    Projection projection;
    ProjectionFramePolicy projection_frame;
    SurfaceLens surface_lens;
    WarpProgram warp_program;
    SignalWeight signal_weight;
    ValueTransfer value_transfer;
    CoveragePolicy coverage;
    PaletteMode palette;
    PeirceLayout peirce_layout = PeirceLayout::SQUARE;
    AiroceanLayout airocean_layout = AiroceanLayout::VERTICAL;
    BonneHemisphere bonne_hemisphere = BonneHemisphere::NORTH;
    GnomonicHemispherePolicy gnomonic_hemisphere =
        GnomonicHemispherePolicy::FOLDED;
    SurfaceNoise surface_noise = SurfaceNoise::NONE;
    SurfaceNoisePlacement surface_noise_placement =
        SurfaceNoisePlacement::AFTER_LENS;
    HueShiftMode hue_shift = HueShiftMode::NOISE;
    PaletteMapping palette_mapping = PaletteMapping::LINEAR;
    BrightnessEnvelope brightness_envelope = BrightnessEnvelope::NONE;

    constexpr bool operator==(const Slots &) const = default;
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

    HS_COLD_MEMBER constexpr SourceParams() = default;

    constexpr SourceParams(float pattern_freq, float speed, float complexity,
                           float pattern_mix, float secondary_rate,
                           float angle_rate = 0.0f)
        : pattern_freq(pattern_freq), speed(speed), complexity(complexity),
          pattern_mix(pattern_mix), secondary_rate(secondary_rate),
          angle_rate(angle_rate) {}

    HS_COLD_MEMBER bool operator==(const SourceParams &) const = default;

    HS_COLD_MEMBER void lerp(const SourceParams &a, const SourceParams &b,
                             float t) {
      // Trips if the field set changes, so a new field cannot silently go
      // uninterpolated and unsnapped.
      static_assert(sizeof(SourceParams) == 64,
                    "SourceParams field set changed - update lerp");
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
      noise_seed = t < 1.0f ? a.noise_seed : b.noise_seed;
      noise_resource_id = t < 1.0f ? a.noise_resource_id : b.noise_resource_id;
    }
  };

  struct WarpStageParams {
    float scale = 1.0f;
    float strength = 0.0f;
    float speed = 0.0f;
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

    HS_COLD_MEMBER constexpr WarpStageParams() = default;

    constexpr WarpStageParams(float scale, float strength, float speed)
        : scale(scale), strength(strength), speed(speed) {}

    HS_COLD_MEMBER bool operator==(const WarpStageParams &) const = default;

    HS_COLD_MEMBER void lerp(const WarpStageParams &a, const WarpStageParams &b,
                             float t, bool rotation_is_rate = false) {
      static_assert(sizeof(WarpStageParams) == 100,
                    "WarpStageParams field set changed - update lerp");
      scale = hs::lerp(a.scale, b.scale, t);
      strength = hs::lerp(a.strength, b.strength, t);
      speed = hs::lerp(a.speed, b.speed, t);
      translation_x = hs::lerp(a.translation_x, b.translation_x, t);
      translation_y = hs::lerp(a.translation_y, b.translation_y, t);
      rotation = rotation_is_rate ? hs::lerp(a.rotation, b.rotation, t)
                                  : lerp_angle(a.rotation, b.rotation, t);
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

    HS_COLD_MEMBER void lerp(const WarpParams &a, const WarpParams &b, float t,
                             const WarpProgram &program = WarpProgram{}) {
      outer.lerp(a.outer, b.outer, t,
                 program.outer.kind == WarpStageKind::AFFINE_FRAME);
      inner.lerp(a.inner, b.inner, t,
                 program.inner.kind == WarpStageKind::AFFINE_FRAME);
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

    HS_COLD_MEMBER constexpr ProjectionParams() = default;

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
      static_assert(sizeof(ProjectionParams) == 28,
                    "ProjectionParams field set changed - update lerp");
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
    MobiusParams mobius{0.7071067811865475f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,
                        0.7071067811865475f, 0.0f};

    HS_COLD_MEMBER constexpr SurfaceLensParams() = default;

    constexpr SurfaceLensParams(float mix) : mix(mix) {}

    HS_COLD_MEMBER bool operator==(const SurfaceLensParams &other) const {
      return mix == other.mix && mobius.a.re == other.mobius.a.re &&
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
      // Trips if the field set changes, so a new field cannot silently go
      // uninterpolated and unsnapped.
      static_assert(sizeof(SurfaceLensParams) == 36,
                    "SurfaceLensParams field set changed - update lerp");
      mix = hs::lerp(a.mix, b.mix, t);
      mobius = t < 1.0f ? a.mobius : b.mobius;
    }
  };

  struct SurfaceNoiseParams {
    NoiseBasis basis = NoiseBasis::SIMPLEX;
    SurfaceCurlIntegrator integrator = SurfaceCurlIntegrator::EULER;
    int32_t seed = 1337;
    float scale = 1.0f;
    float strength = 0.0f;
    float rate = 0.0f;
    float direction = 0.0f;

    HS_COLD_MEMBER bool operator==(const SurfaceNoiseParams &) const = default;

    HS_COLD_MEMBER void lerp(const SurfaceNoiseParams &a,
                             const SurfaceNoiseParams &b, float t) {
      static_assert(sizeof(SurfaceNoiseParams) == 24,
                    "SurfaceNoiseParams field set changed - update lerp");
      basis = t < 1.0f ? a.basis : b.basis;
      integrator = t < 1.0f ? a.integrator : b.integrator;
      seed = t < 1.0f ? a.seed : b.seed;
      scale = hs::lerp(a.scale, b.scale, t);
      strength = hs::lerp(a.strength, b.strength, t);
      rate = hs::lerp(a.rate, b.rate, t);
      direction = ProjectionParams::lerp_periodic(a.direction, b.direction, t);
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
      static_assert(sizeof(ValueParams) == 28,
                    "ValueParams field set changed - update lerp");
      iso_level = hs::lerp(a.iso_level, b.iso_level, t);
      iso_width = hs::lerp(a.iso_width, b.iso_width, t);
      band_count = t < 1.0f ? a.band_count : b.band_count;
      band_phase = WarpStageParams::lerp_angle(a.band_phase, b.band_phase, t);
      cutout_threshold = hs::lerp(a.cutout_threshold, b.cutout_threshold, t);
      cutout_softness = hs::lerp(a.cutout_softness, b.cutout_softness, t);
      edge_width = hs::lerp(a.edge_width, b.edge_width, t);
    }
  };

  struct ColorParams {
    float hue_shift_amount = 0.0f;
    float hue_noise_scale = 1.0f;
    float hue_noise_speed = 0.0f;
    float palette_chroma = 0.62f;
    float mapping_frequency = 1.0f;
    float mapping_phase = 0.0f;
    float phase_oscillation_depth = 0.0f;
    float phase_oscillation_speed = 0.0f;
    float brightness_depth = 1.0f;
    float value_opacity_low = 1.0f;
    float value_opacity_high = 1.0f;

    HS_COLD_MEMBER constexpr ColorParams() = default;

    constexpr ColorParams(float amount, float scale, float speed)
        : hue_shift_amount(amount), hue_noise_scale(scale),
          hue_noise_speed(speed) {}

    HS_COLD_MEMBER bool operator==(const ColorParams &) const = default;

    HS_COLD_MEMBER void lerp(const ColorParams &a, const ColorParams &b,
                             float t) {
      static_assert(sizeof(ColorParams) == 44,
                    "ColorParams field set changed - update lerp");
      hue_shift_amount = hs::lerp(a.hue_shift_amount, b.hue_shift_amount, t);
      hue_noise_scale = hs::lerp(a.hue_noise_scale, b.hue_noise_scale, t);
      hue_noise_speed = hs::lerp(a.hue_noise_speed, b.hue_noise_speed, t);
      palette_chroma = hs::lerp(a.palette_chroma, b.palette_chroma, t);
      mapping_frequency = hs::lerp(a.mapping_frequency, b.mapping_frequency, t);
      mapping_phase = hs::lerp(a.mapping_phase, b.mapping_phase, t);
      phase_oscillation_depth =
          hs::lerp(a.phase_oscillation_depth, b.phase_oscillation_depth, t);
      phase_oscillation_speed =
          hs::lerp(a.phase_oscillation_speed, b.phase_oscillation_speed, t);
      brightness_depth = hs::lerp(a.brightness_depth, b.brightness_depth, t);
      value_opacity_low = hs::lerp(a.value_opacity_low, b.value_opacity_low, t);
      value_opacity_high =
          hs::lerp(a.value_opacity_high, b.value_opacity_high, t);
    }
  };

  struct OuterCameraParams {
    float wander = 0.0f;

    HS_COLD_MEMBER bool operator==(const OuterCameraParams &) const = default;

    HS_COLD_MEMBER void lerp(const OuterCameraParams &a,
                             const OuterCameraParams &b, float t) {
      static_assert(sizeof(OuterCameraParams) == 4,
                    "OuterCameraParams field set changed - update lerp");
      wander = hs::lerp(a.wander, b.wander, t);
    }
  };

  struct Params {
    SourceParams source;
    WarpParams warp;
    ProjectionParams projection;
    SurfaceLensParams surface_lens;
    ValueParams value;
    ColorParams color;
    OuterCameraParams outer_camera;
    SurfaceNoiseParams surface_noise;

    HS_COLD_MEMBER Params() = default;

    constexpr Params(SourceParams source, WarpParams warp,
                     ProjectionParams projection,
                     SurfaceLensParams surface_lens, ValueParams value,
                     ColorParams color, OuterCameraParams outer_camera,
                     SurfaceNoiseParams surface_noise)
        : source(source), warp(warp), projection(projection),
          surface_lens(surface_lens), value(value), color(color),
          outer_camera(outer_camera), surface_noise(surface_noise) {}

    HS_COLD_MEMBER bool operator==(const Params &) const = default;

    HS_COLD_MEMBER void lerp(const Params &a, const Params &b, float t,
                             const Slots &slots = Slots{}) {
      source.lerp(a.source, b.source, t);
      warp.lerp(a.warp, b.warp, t, slots.warp_program);
      projection.lerp(a.projection, b.projection, t);
      surface_lens.lerp(a.surface_lens, b.surface_lens, t);
      surface_noise.lerp(a.surface_noise, b.surface_noise, t);
      value.lerp(a.value, b.value, t);
      color.lerp(a.color, b.color, t);
      outer_camera.lerp(a.outer_camera, b.outer_camera, t);
    }

    HS_COLD_MEMBER void lerp_staggered(const Params &a, const Params &b,
                                       float t, const Slots &slots = Slots{}) {
      const int phase_count = (a.source != b.source) + (a.warp != b.warp) +
                              (a.projection != b.projection) +
                              (a.surface_lens != b.surface_lens) +
                              (a.value != b.value) + (a.color != b.color) +
                              (a.outer_camera != b.outer_camera) +
                              (a.surface_noise != b.surface_noise);
      int phase = 0;
      source = a.source;
      warp = a.warp;
      projection = a.projection;
      surface_lens = a.surface_lens;
      value = a.value;
      color = a.color;
      outer_camera = a.outer_camera;
      surface_noise = a.surface_noise;
      if (a.source != b.source)
        source.lerp(a.source, b.source, phase_t(t, phase++, phase_count));
      if (a.warp != b.warp)
        warp.lerp(a.warp, b.warp, phase_t(t, phase++, phase_count),
                  slots.warp_program);
      if (a.projection != b.projection)
        projection.lerp(a.projection, b.projection,
                        phase_t(t, phase++, phase_count));
      if (a.surface_lens != b.surface_lens)
        surface_lens.lerp(a.surface_lens, b.surface_lens,
                          phase_t(t, phase++, phase_count));
      if (a.value != b.value)
        value.lerp(a.value, b.value, phase_t(t, phase++, phase_count));
      if (a.color != b.color)
        color.lerp(a.color, b.color, phase_t(t, phase++, phase_count));
      if (a.outer_camera != b.outer_camera)
        outer_camera.lerp(a.outer_camera, b.outer_camera,
                          phase_t(t, phase++, phase_count));
      if (a.surface_noise != b.surface_noise)
        surface_noise.lerp(a.surface_noise, b.surface_noise,
                           phase_t(t, phase, phase_count));
    }

    HS_COLD_MEMBER static float phase_t(float t, int phase, int phase_count) {
      return ease_in_out_sin(hs::clamp(t * phase_count - phase, 0.0f, 1.0f));
    }
  };

  struct Config {
    Slots slots;
    Params params;

    HS_COLD_MEMBER constexpr Config() = default;
    constexpr Config(const Slots &slots, const Params &params)
        : slots(slots), params(params) {}

    HS_COLD_MEMBER bool operator==(const Config &) const = default;
  };
  using RequestedConfig = Config;

  static constexpr uint32_t CONFIG_SCHEMA_VERSION = 7;

  /**
   * @brief Reports whether a persisted snapshot's schema version can be
   *        restored.
   * @param version Snapshot schema version to test.
   * @return true for the current version.
   * @details Single source of truth for the accepted set, so callers that
   *          pre-screen a version (the WASM bridge) cannot drift from what the
   *          effect actually accepts.
   */
  static constexpr bool config_version_supported(uint32_t version) {
    return version == CONFIG_SCHEMA_VERSION;
  }

  enum class ConfigFieldId : uint16_t {
#define HS_SHADERBALL_FIELD_ENUM(name, path) name,
    HS_SHADERBALL_CONFIG_FIELDS(HS_SHADERBALL_FIELD_ENUM)
#undef HS_SHADERBALL_FIELD_ENUM
        COUNT
  };

  static constexpr size_t CONFIG_FIELD_COUNT =
      static_cast<size_t>(ConfigFieldId::COUNT);

  struct ConfigFieldLayout {
    size_t offset;
    size_t size;
  };

  enum class ConfigRestoreResult : uint8_t {
    APPLIED,
    UNSUPPORTED_VERSION,
    INVALID_LENGTH,
    INVALID_VALUE,
    INVALID_ACCEPTED,
    INVALID_PENDING
  };

  enum class RuntimeFieldId : uint8_t {
    SOURCE_PRIMARY,
    SOURCE_SECONDARY,
    SOURCE_ANGLE,
    WARP_OUTER_ROTATION,
    PROJECTION_SPIN,
    HUE_NOISE_PHASE,
    SOURCE_NOISE_PHASE,
    WARP_INNER_ROTATION,
    SURFACE_NOISE_PHASE,
    WARP_OUTER_PHASE,
    WARP_INNER_PHASE,
    PALETTE_OSCILLATION_PHASE,
    COUNT
  };

  static constexpr size_t RUNTIME_FIELD_COUNT =
      static_cast<size_t>(RuntimeFieldId::COUNT);
  using ConfigValues = std::array<uint32_t, CONFIG_FIELD_COUNT>;
  using RuntimeValues = std::array<float, RUNTIME_FIELD_COUNT>;

  struct FullConfigSnapshot {
    uint32_t schema_version = CONFIG_SCHEMA_VERSION;
    ConfigValues accepted{};
    ConfigValues requested{};
    std::array<uint8_t, CONFIG_FIELD_COUNT> pending{};
    bool has_runtime = false;
    RuntimeValues runtime{};
  };

  struct PendingEdit {
    const char *name = nullptr;
    ConfigFieldId id = ConfigFieldId::COUNT;
    size_t offset = 0;
    size_t size = 0;
  };

  /** @brief Stable name for a full-config field ID. */
  static constexpr const char *config_field_name(ConfigFieldId id) {
    switch (id) {
#define HS_SHADERBALL_FIELD_NAME(name, path)                                   \
  case ConfigFieldId::name:                                                    \
    return #path;
      HS_SHADERBALL_CONFIG_FIELDS(HS_SHADERBALL_FIELD_NAME)
#undef HS_SHADERBALL_FIELD_NAME
    case ConfigFieldId::COUNT:
      break;
    }
    return nullptr;
  }

  static ConfigFieldLayout config_field_layout(ConfigFieldId id) {
    Config config{};
    const uintptr_t base = reinterpret_cast<uintptr_t>(&config);
    switch (id) {
#define HS_SHADERBALL_FIELD_LAYOUT(name, path)                                 \
  case ConfigFieldId::name:                                                    \
    return {reinterpret_cast<uintptr_t>(&config.path) - base,                  \
            sizeof(config.path)};
      HS_SHADERBALL_CONFIG_FIELDS(HS_SHADERBALL_FIELD_LAYOUT)
#undef HS_SHADERBALL_FIELD_LAYOUT
    case ConfigFieldId::COUNT:
      break;
    }
    return {sizeof(Config), 0};
  }

#if HS_ENABLE_PARAM_GUI_BRIDGE
  /** @brief Reserved compatibility accessor. */
  const char *config_import_notice() const { return ""; }

  /** @brief Reserved compatibility no-op. */
  void clear_config_import_notice() {}
#endif

private:
  static constexpr float lens_domain_linear_scale(SurfaceLens lens) {
    switch (lens) {
    case SurfaceLens::KALEIDOSCOPE:
      return 1.0f;
    case SurfaceLens::KALEIDOSCOPE_TETRAHEDRAL:
      return 0.20412415f;
    case SurfaceLens::KALEIDOSCOPE_OCTAHEDRAL:
      return 0.14433757f;
    case SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL:
      return 0.09128709f;
    case SurfaceLens::KALEIDOSCOPE_TRIANGULAR_PRISM:
      return 0.28867513f;
    case SurfaceLens::KALEIDOSCOPE_SQUARE_PRISM:
      return 0.25f;
    case SurfaceLens::KALEIDOSCOPE_PENTAGONAL_PRISM:
      return 0.22360680f;
    case SurfaceLens::KALEIDOSCOPE_HEXAGONAL_PRISM:
      return 0.20412415f;
    case SurfaceLens::KALEIDOSCOPE_OCTAGONAL_PRISM:
      return 0.17677670f;
    case SurfaceLens::NONE:
    case SurfaceLens::GLITCH:
    case SurfaceLens::TWIST:
    case SurfaceLens::MOBIUS:
    case SurfaceLens::TANGENT_NOISE:
      return 1.0f;
    }
    __builtin_unreachable();
  }

  static constexpr float domain_scaled_max(float full_domain_max,
                                           float minimum_max,
                                           float domain_scale) {
    const float scaled = full_domain_max * domain_scale;
    return scaled > minimum_max ? scaled : minimum_max;
  }

  HS_COLD_MEMBER void register_clamped_animated_param(const char *name,
                                                      float *target,
                                                      float minimum,
                                                      float maximum) {
    const float clamped = hs::clamp(*target, minimum, maximum);
    registered_range_clamped |= clamped != *target;
    *target = clamped;
    register_animated_param(name, target, minimum, maximum);
  }

  HS_COLD_MEMBER void rebind_parameters() {
    registered_range_clamped = false;
    reset_parameters();
    Slots &slots = requested_config.slots;
    if (!fixed_topology)
      register_animated_param("Function", &slots.function, FUNCTION_OPTIONS,
                              FUNCTION_EXPORT_OPTIONS, NUM_FUNCTIONS);
    const float domain_scale = lens_domain_linear_scale(slots.surface_lens);
    register_source_controls(slots.function, requested_config.params.source,
                             domain_scale);
    if (!fixed_topology)
      register_animated_param("Projection", &slots.projection,
                              PROJECTION_OPTIONS, PROJECTION_EXPORT_OPTIONS,
                              NUM_PROJECTIONS);
    register_projection_controls(slots, requested_config.params);
    if (!fixed_topology)
      register_animated_param(
          "Projection Frame", &slots.projection_frame, PROJECTION_FRAME_OPTIONS,
          PROJECTION_FRAME_EXPORT_OPTIONS, NUM_PROJECTION_FRAMES);
    register_projection_frame_controls(slots.projection_frame,
                                       requested_config.params, domain_scale);
    register_animated_param("Camera Wander",
                            &requested_config.params.outer_camera.wander,
                            WANDER_MIN, WANDER_MAX);
    if (!fixed_topology)
      register_animated_param("Surface Noise", &slots.surface_noise,
                              SURFACE_NOISE_OPTIONS,
                              SURFACE_NOISE_EXPORT_OPTIONS, NUM_SURFACE_NOISE);
    register_surface_noise_controls(
        slots, requested_config.params.surface_noise,
        slots.surface_noise_placement == SurfaceNoisePlacement::AFTER_LENS
            ? domain_scale
            : 1.0f);
    if (!fixed_topology)
      register_animated_param("Lens", &slots.surface_lens, LENS_OPTIONS,
                              LENS_EXPORT_OPTIONS, NUM_LENSES);
    register_lens_controls(slots.surface_lens,
                           requested_config.params.surface_lens);
    if (!fixed_topology) {
      register_animated_param("Planar Warp 1", &slots.warp_program.outer.kind,
                              WARP_OPTIONS, WARP_EXPORT_OPTIONS, NUM_WARPS);
      register_stage_slot_controls(true, slots.warp_program.outer);
    }
    register_active_warp_controls(true, slots.warp_program.outer,
                                  requested_config.params.warp.outer,
                                  domain_scale);
    if (!fixed_topology) {
      register_animated_param("Planar Warp 2", &slots.warp_program.inner.kind,
                              WARP_OPTIONS, WARP_EXPORT_OPTIONS, NUM_WARPS);
      register_stage_slot_controls(false, slots.warp_program.inner);
    }
    register_active_warp_controls(false, slots.warp_program.inner,
                                  requested_config.params.warp.inner,
                                  domain_scale);
    if (!fixed_topology) {
      register_animated_param("Signal Weight", &slots.signal_weight,
                              SIGNAL_OPTIONS, SIGNAL_EXPORT_OPTIONS,
                              NUM_SIGNALS);
      register_animated_param(
          "Value Transfer", &slots.value_transfer, VALUE_TRANSFER_OPTIONS,
          VALUE_TRANSFER_EXPORT_OPTIONS, NUM_VALUE_TRANSFERS);
    }
    register_value_transfer_controls(slots.value_transfer,
                                     requested_config.params.value);
    if (!fixed_topology)
      register_animated_param("Coverage", &slots.coverage, COVERAGE_OPTIONS,
                              COVERAGE_EXPORT_OPTIONS, NUM_COVERAGE_POLICIES);
    register_coverage_controls(slots.coverage, requested_config.params.value);
    if (!fixed_topology)
      register_animated_param("Palette", &slots.palette, PALETTE_OPTIONS,
                              PALETTE_EXPORT_OPTIONS, NUM_PALETTES);
    register_animated_param("Palette Chroma",
                            &requested_config.params.color.palette_chroma,
                            PALETTE_CHROMA_MIN, PALETTE_CHROMA_MAX);
    register_animated_param(
        "Palette Mapping", &slots.palette_mapping, PALETTE_MAPPING_OPTIONS,
        PALETTE_MAPPING_EXPORT_OPTIONS, NUM_PALETTE_MAPPINGS);
    register_animated_param("Mapping Frequency",
                            &requested_config.params.color.mapping_frequency,
                            MAPPING_FREQUENCY_MIN, MAPPING_FREQUENCY_MAX);
    register_animated_param("Mapping Phase",
                            &requested_config.params.color.mapping_phase,
                            MAPPING_PHASE_MIN, MAPPING_PHASE_MAX);
    register_animated_param(
        "Phase Oscillation Depth",
        &requested_config.params.color.phase_oscillation_depth,
        PHASE_OSCILLATION_DEPTH_MIN, PHASE_OSCILLATION_DEPTH_MAX);
    register_animated_param(
        "Phase Oscillation Speed",
        &requested_config.params.color.phase_oscillation_speed,
        -PHASE_OSCILLATION_SPEED_MAX, PHASE_OSCILLATION_SPEED_MAX);
    if (!fixed_topology)
      register_animated_param("Brightness Envelope", &slots.brightness_envelope,
                              BRIGHTNESS_ENVELOPE_OPTIONS,
                              BRIGHTNESS_ENVELOPE_EXPORT_OPTIONS,
                              NUM_BRIGHTNESS_ENVELOPES);
    if (slots.brightness_envelope != BrightnessEnvelope::NONE)
      register_animated_param("Brightness Depth",
                              &requested_config.params.color.brightness_depth,
                              BRIGHTNESS_DEPTH_MIN, BRIGHTNESS_DEPTH_MAX);
    register_animated_param("Value Opacity Low",
                            &requested_config.params.color.value_opacity_low,
                            VALUE_OPACITY_MIN, VALUE_OPACITY_MAX);
    register_animated_param("Value Opacity High",
                            &requested_config.params.color.value_opacity_high,
                            VALUE_OPACITY_MIN, VALUE_OPACITY_MAX);
    if (!fixed_topology)
      register_animated_param("Hue Shift Mode", &slots.hue_shift,
                              HUE_SHIFT_OPTIONS, HUE_SHIFT_EXPORT_OPTIONS,
                              NUM_HUE_SHIFT_MODES);
    register_color_controls(slots.hue_shift, requested_config.params.color,
                            domain_scale);
#if HS_ENABLE_PARAM_GUI_BRIDGE
    const bool post_registration_clamp = clamp_registered_parameter_ranges();
    if (requested_schema_bound &&
        (registered_range_clamped || post_registration_clamp))
      refresh_accepted_config();
    for (size_t index = 0; index < pending_edit_count; ++index) {
      PendingEdit &edit = pending_edits[index];
      edit.name = nullptr;
      const uintptr_t target =
          reinterpret_cast<uintptr_t>(&requested_config) + edit.offset;
      for (const ParamDef &parameter : getParameters()) {
        if (reinterpret_cast<uintptr_t>(parameter.target) == target) {
          edit.name = parameter.name;
          break;
        }
      }
    }
    mirror_parameter_display_state(requested_config, display_config);
    for (size_t index = 0; index < pending_edit_count; ++index)
      if (pending_edits[index].name != nullptr)
        show_requested_parameter_value(pending_edits[index].name);
#endif
    requested_schema_bound = true;
  }

#if HS_ENABLE_PARAM_GUI_BRIDGE
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
    const ConfigFieldId id = config_field_id(offset, size);
    HS_CHECK(id != ConfigFieldId::COUNT,
             "ShaderBall parameter lacks a stable field ID");
    const size_t before_count = pending_edit_count;
    const bool was_pending = pending_edit_at(id) < pending_edit_count;
    remember_pending_edit(name, id, offset, size);
    refresh_accepted_config();
    const bool is_pending = pending_edit_at(id) < pending_edit_count;
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

  size_t pending_edit_at(ConfigFieldId id) const {
    for (size_t index = 0; index < pending_edit_count; ++index)
      if (pending_edits[index].id == id)
        return index;
    return pending_edit_count;
  }

  static ConfigFieldId config_field_id(size_t offset, size_t size) {
    Config config{};
    const uintptr_t base = reinterpret_cast<uintptr_t>(&config);
#define HS_SHADERBALL_FIELD_MATCH(name, path)                                  \
  if (reinterpret_cast<uintptr_t>(&config.path) - base == offset &&            \
      sizeof(config.path) == size)                                             \
    return ConfigFieldId::name;
    HS_SHADERBALL_CONFIG_FIELDS(HS_SHADERBALL_FIELD_MATCH)
#undef HS_SHADERBALL_FIELD_MATCH
    return ConfigFieldId::COUNT;
  }

  void remember_pending_edit(const char *name, ConfigFieldId id, size_t offset,
                             size_t size) {
    const size_t existing = pending_edit_at(id);
    if (existing < pending_edit_count) {
      pending_edits[existing].name = name;
      pending_edits[existing].size = size;
      return;
    }
    HS_CHECK(pending_edit_count < pending_edits.size(),
             "ShaderBall pending edit capacity exceeded");
    pending_edits[pending_edit_count++] = {name, id, offset, size};
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
    if (admissible_config(requested_config)) {
      accepted_config = requested_config;
      pending_edit_count = 0;
      return;
    }

    Config candidate = requested_config;
    for (size_t index = 0; index < pending_edit_count; ++index)
      copy_pending_value(candidate, accepted_config, pending_edits[index]);
    if (admissible_config(candidate))
      accepted_config = candidate;

    bool admitted;
    do {
      admitted = false;
      for (size_t index = 0; index < pending_edit_count;) {
        candidate = accepted_config;
        copy_pending_value(candidate, requested_config, pending_edits[index]);
        if (!admissible_config(candidate)) {
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
    const ParamDef *parameter = getParameters().find(name);
    if (parameter != nullptr && parameter_out_of_range(*parameter))
      return begin_warning(
          "%s %.7g is outside its registered range [%.7g, %.7g]. Set %s "
          "within that range.",
          name, static_cast<double>(parameter->get_requested()),
          static_cast<double>(parameter->min),
          static_cast<double>(parameter->max), name);
    for (size_t index = 0; index < pending_edit_count; ++index) {
      const PendingEdit &edit = pending_edits[index];
      if (edit.name == nullptr || std::strcmp(edit.name, name) != 0)
        continue;
      if (schema_selector(name) && range_repairs_admission())
        return nullptr;
      return admission_warning(requested_config, edit.name);
    }
    return nullptr;
  }

  static bool parameter_out_of_range(const ParamDef &parameter) {
    const float value = parameter.get_requested();
    return value < parameter.min || value > parameter.max;
  }

  bool clamp_registered_parameter_ranges() {
    const uintptr_t requested = reinterpret_cast<uintptr_t>(&requested_config);
    bool clamped = false;
    for (const ParamDef &parameter : getParameters()) {
      if (parameter.is_bool() || parameter.is_enum() ||
          !parameter_out_of_range(parameter))
        continue;
      const uintptr_t target = reinterpret_cast<uintptr_t>(parameter.target);
      const size_t size = parameter_target_size(parameter);
      if (target < requested ||
          target + size > requested + sizeof(requested_config))
        continue;
      ParamDef writable = parameter;
      writable.set(
          hs::clamp(parameter.get_requested(), parameter.min, parameter.max));
      clamped = true;
    }
    return clamped;
  }

  bool range_repairs_admission() const {
    Config candidate = requested_config;
    const uintptr_t requested = reinterpret_cast<uintptr_t>(&requested_config);
    bool repaired = false;
    for (const ParamDef &parameter : getParameters()) {
      if (!parameter_out_of_range(parameter))
        continue;
      const uintptr_t target = reinterpret_cast<uintptr_t>(parameter.target);
      const size_t size = parameter_target_size(parameter);
      if (target < requested ||
          target + size > requested + sizeof(requested_config))
        continue;
      ParamDef candidate_parameter = parameter;
      candidate_parameter.target =
          reinterpret_cast<uint8_t *>(&candidate) + (target - requested);
      candidate_parameter.set(
          hs::clamp(parameter.get_requested(), parameter.min, parameter.max));
      repaired = true;
    }
    return repaired && admissible_config(candidate);
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
           std::strcmp(name, "Surface Noise") == 0 ||
           std::strcmp(name, "Surface Noise Placement") == 0 ||
           std::strcmp(name, "Lens") == 0 ||
           std::strcmp(name, "Planar Warp 1") == 0 ||
           std::strcmp(name, "Planar Warp 2") == 0 ||
           std::strcmp(name, "Value Transfer") == 0 ||
           std::strcmp(name, "Coverage") == 0 ||
           std::strcmp(name, "Palette") == 0 ||
           std::strcmp(name, "Brightness Envelope") == 0 ||
           std::strcmp(name, "Hue Shift Mode") == 0;
  }
#endif

  HS_COLD_MEMBER void register_value_transfer_controls(ValueTransfer transfer,
                                                       ValueParams &params) {
    if (transfer == ValueTransfer::ISO_CONTOUR) {
      register_animated_param("Iso Level", &params.iso_level, 0.0f, 1.0f);
      register_animated_param("Iso Width", &params.iso_width, SOFTNESS_MIN,
                              0.5f);
    } else if (transfer == ValueTransfer::SMOOTH_BANDS) {
      if (!fixed_topology)
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
                                               float domain_scale) {
    if (is_noise_contour(function)) {
      register_clamped_animated_param(
          "Source Noise Scale", &params.noise_scale, SOURCE_NOISE_SCALE_MIN,
          domain_scaled_max(SOURCE_NOISE_SCALE_MAX, 0.5f, domain_scale));
      register_animated_param("Source Noise Contrast", &params.noise_contrast,
                              0.0f, 8.0f);
      register_clamped_animated_param(
          "Source Noise Speed", &params.noise_time_rate,
          -domain_scaled_max(SOURCE_NOISE_RATE_MAX, 1.0f / 4096.0f,
                             domain_scale),
          domain_scaled_max(SOURCE_NOISE_RATE_MAX, 1.0f / 4096.0f,
                            domain_scale));
      if (!fixed_topology)
        register_animated_param("Source Noise Basis", &params.noise_basis,
                                NOISE_BASIS_OPTIONS, NOISE_BASIS_EXPORT_OPTIONS,
                                NUM_NOISE_BASES);
      return;
    }
    if (function == Function::PRIMITIVE_LATTICE) {
      register_clamped_animated_param(
          "Lattice Cell Scale", &params.lattice_cell_scale, CELL_MIN, CELL_MAX);
      register_animated_param("Lattice Shape", &params.lattice_shape_blend,
                              0.0f, 1.0f);
      register_animated_param("Lattice Softness", &params.lattice_softness,
                              SOFTNESS_MIN, 1.0f);
      register_animated_param("Lattice Radius", &params.lattice_radius,
                              1.0f / 64.0f, 0.49f);
      return;
    }
    register_clamped_animated_param("Pattern Freq", &params.pattern_freq,
                                    PATTERN_FREQ_MIN, PATTERN_FREQ_MAX);
    register_clamped_animated_param(
        "Speed", &params.speed, SPEED_MIN,
        domain_scaled_max(SPEED_MAX, 0.5f, domain_scale));
    register_clamped_animated_param(
        "Source Angle Speed", &params.angle_rate, WAVE_SPIN_MIN,
        domain_scaled_max(WAVE_SPIN_MAX, 0.03f, domain_scale));
    if (function == Function::GRID) {
      register_animated_param("Complexity", &params.complexity, COMPLEXITY_MIN,
                              COMPLEXITY_MAX);
      register_animated_param("Pattern Mix", &params.pattern_mix,
                              PATTERN_MIX_MIN, PATTERN_MIX_MAX);
      register_clamped_animated_param(
          "Drift", &params.secondary_rate, PHASE2_RATE_MIN,
          domain_scaled_max(PHASE2_RATE_MAX, 1.25f, domain_scale));
    }
  }

  HS_COLD_MEMBER void register_projection_controls(Slots &slots,
                                                   Params &params) {
    if (!fixed_topology && slots.projection == Projection::PEIRCE_QUINCUNCIAL)
      register_animated_param("Peirce Layout", &slots.peirce_layout,
                              PEIRCE_LAYOUT_OPTIONS,
                              PEIRCE_LAYOUT_EXPORT_OPTIONS, NUM_PEIRCE_LAYOUTS);
    if (!fixed_topology && slots.projection == Projection::AIROCEAN)
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
    if (!fixed_topology && slots.projection == Projection::BONNE)
      register_animated_param(
          "Bonne Hemisphere", &slots.bonne_hemisphere, BONNE_HEMISPHERE_OPTIONS,
          BONNE_HEMISPHERE_EXPORT_OPTIONS, NUM_BONNE_HEMISPHERES);
    if (!fixed_topology && slots.projection == Projection::GNOMONIC)
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
                                     Params &params, float domain_scale) {
    if (frame == ProjectionFramePolicy::SPIN_WANDER) {
      register_clamped_animated_param(
          "Projection Spin Speed", &params.projection.spin_rate, SPIN_RATE_MIN,
          domain_scaled_max(SPIN_RATE_MAX, 0.04f, domain_scale));
      register_animated_param("Projection Wander", &params.projection.wander,
                              WANDER_MIN, WANDER_MAX);
    }
  }

  HS_COLD_MEMBER void register_lens_controls(SurfaceLens lens,
                                             SurfaceLensParams &params) {
    if (lens == SurfaceLens::NONE)
      return;
    if (lens == SurfaceLens::MOBIUS) {
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

  HS_COLD_MEMBER void
  register_surface_noise_controls(Slots &slots, SurfaceNoiseParams &params,
                                  float domain_scale) {
    if (slots.surface_noise == SurfaceNoise::NONE)
      return;
    if (!fixed_topology) {
      register_animated_param(
          "Surface Noise Placement", &slots.surface_noise_placement,
          SURFACE_NOISE_PLACEMENT_OPTIONS,
          SURFACE_NOISE_PLACEMENT_EXPORT_OPTIONS, NUM_SURFACE_NOISE_PLACEMENTS);
      register_animated_param("Surface Noise Basis", &params.basis,
                              NOISE_BASIS_OPTIONS, NOISE_BASIS_EXPORT_OPTIONS,
                              NUM_NOISE_BASES);
    }
    register_clamped_animated_param("Surface Noise Scale", &params.scale,
                                    LENS_NOISE_SCALE_MIN, LENS_NOISE_SCALE_MAX);
    const float strength_min =
        slots.surface_noise == SurfaceNoise::CURL ? -0.5f : 0.0f;
#if HS_ENABLE_PARAM_GUI_BRIDGE
    register_animated_param_preserving_value(
        "Surface Noise Strength", &params.strength, strength_min, 0.5f);
#else
    register_animated_param("Surface Noise Strength", &params.strength,
                            strength_min, 0.5f);
#endif
    const float speed_max =
        domain_scaled_max(NOISE_RATE_MAX, 0.002f, domain_scale);
    register_clamped_animated_param("Surface Noise Speed", &params.rate,
                                    -speed_max, speed_max);
    if (slots.surface_noise == SurfaceNoise::DIRECT)
      register_animated_param("Surface Noise Direction", &params.direction,
                              0.0f, 1.0f);
    else if (!fixed_topology)
      register_animated_param("Surface Noise Integrator", &params.integrator,
                              SURFACE_CURL_INTEGRATOR_OPTIONS,
                              SURFACE_CURL_INTEGRATOR_EXPORT_OPTIONS,
                              NUM_SURFACE_CURL_INTEGRATORS);
  }

  HS_COLD_MEMBER void register_active_warp_controls(bool first,
                                                    const WarpStageSpec &spec,
                                                    WarpStageParams &params,
                                                    float domain_scale) {
    if (spec.kind == WarpStageKind::NONE)
      return;
    const char *const *names =
        first ? FIRST_WARP_PARAM_NAMES : SECOND_WARP_PARAM_NAMES;
    const char *speed_name =
        first ? "Planar Warp 1 Speed" : "Planar Warp 2 Speed";
    auto register_current = [&](const char *name, float *target, float minimum,
                                float maximum) {
      register_clamped_animated_param(name, target, minimum, maximum);
    };
    if (spec.kind == WarpStageKind::WAVE_SHEAR ||
        spec.kind == WarpStageKind::VECTOR_NOISE ||
        spec.kind == WarpStageKind::CURL_FLOW) {
      const char *strength_name =
          first ? "Planar Warp 1 Strength" : "Planar Warp 2 Strength";
      const bool signed_strength = spec.kind == WarpStageKind::WAVE_SHEAR ||
                                   spec.kind == WarpStageKind::CURL_FLOW;
      float strength_max = spec.kind == WarpStageKind::VECTOR_NOISE
                               ? VECTOR_WARP_STRENGTH_MAX
                               : 4.0f;
      if (spec.kind == WarpStageKind::CURL_FLOW)
        strength_max = curl_strength_limit(spec, params);
      register_current(strength_name, &params.strength,
                       signed_strength ? -strength_max : 0.0f, strength_max);
    }
    const float speed_max =
        domain_scaled_max(NOISE_SPEED_MAX, 0.005f, domain_scale);
    register_current(speed_name, &params.speed, -speed_max, speed_max);
    switch (spec.kind) {
    case WarpStageKind::AFFINE_FRAME: {
      const float snapped_x = roundf(params.translation_x);
      const float snapped_y = roundf(params.translation_y);
      registered_range_clamped |= snapped_x != params.translation_x ||
                                  snapped_y != params.translation_y;
      params.translation_x = snapped_x;
      params.translation_y = snapped_y;
      for (int index = 0; index < 6; ++index) {
        float *targets[] = {&params.translation_x, &params.translation_y,
                            &params.rotation,      &params.scale_x,
                            &params.scale_y,       &params.shear};
        const float minimum[] = {-4.0f, -4.0f, -PI_F / 8.0f,
                                 0.25f, 0.25f, -0.75f};
        const float maximum[] = {4.0f, 4.0f, PI_F / 8.0f, 4.0f, 4.0f, 0.75f};
        register_current(names[WARP_NAME_TRANSLATION_X + index], targets[index],
                         minimum[index], maximum[index]);
      }
      break;
    }
    case WarpStageKind::WAVE_SHEAR:
      register_current(names[WARP_NAME_FREQUENCY], &params.frequency, 0.0f,
                       domain_scaled_max(64.0f, 8.0f, domain_scale));
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
                       domain_scaled_max(spec.kind == WarpStageKind::CURL_FLOW
                                             ? CURL_WARP_SCALE_MAX
                                             : VECTOR_WARP_SCALE_MAX,
                                         1.0f, domain_scale));
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

  HS_COLD_MEMBER void register_color_controls(HueShiftMode mode,
                                              ColorParams &params,
                                              float domain_scale) {
    if (mode == HueShiftMode::NONE)
      return;
    register_clamped_animated_param(
        "Hue Shift Amount", &params.hue_shift_amount,
        -hue_shift_amount_max(mode), hue_shift_amount_max(mode));
    if (mode != HueShiftMode::NOISE)
      return;
    register_clamped_animated_param(
        "Hue Noise Scale", &params.hue_noise_scale, HUE_NOISE_SCALE_MIN,
        domain_scaled_max(HUE_NOISE_SCALE_MAX, 2.0f, domain_scale));
    register_clamped_animated_param("Hue Noise Speed", &params.hue_noise_speed,
                                    -HUE_NOISE_SPEED_MAX, HUE_NOISE_SPEED_MAX);
  }

  struct Blend {
    Params params;
    PaletteMappingWeights palette_mapping;
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

  using ProjectedLookup = Pullback::ProjectionSample;

  struct SourceTraits {
    bool y_periodic;
    bool polar_angle_compatible;
  };

  using PlanarWarpResult = Pullback::WarpResult;
  using SurfaceNoiseResult = Pullback::SurfaceResult;
  using PlanarWarpStageResult = Pullback::WarpStepResult;
  using MaterialSample = Pullback::MaterialSample;

  struct ClockState {
    float source_primary = 0.0f;
    float source_secondary = 0.0f;
    float source_angle = 0.0f;
    float projection_spin = 0.0f;
    float hue_noise_phase = 0.0f;
    float source_noise_time = 0.0f;
    float surface_noise_time = 0.0f;
    float warp_outer_phase = 0.0f;
    float warp_inner_phase = 0.0f;
    float warp_outer_rotation = 0.0f;
    float warp_inner_rotation = 0.0f;
    float palette_oscillation_phase = 0.0f;

    HS_COLD_MEMBER constexpr ClockState() = default;
    constexpr ClockState(float source_primary, float source_secondary,
                         float source_angle, float projection_spin,
                         float hue_noise_phase)
        : source_primary(source_primary), source_secondary(source_secondary),
          source_angle(source_angle), projection_spin(projection_spin),
          hue_noise_phase(hue_noise_phase) {}
  };

  struct PreparedTransforms {
    Quaternion projection_conj;
    Quaternion outer_conj;
  };

  struct PreparedAffineFrame {
    float translation_x;
    float translation_y;
    float scale_x;
    float scale_y;
    float shear;
  };

  struct PreparedMirrorTile {
    float offset_x;
    float offset_y;
  };

  struct PreparedVortex {
    float center_x;
    float center_y;
    float radius_sq;
    float angle_numerator;
  };

  struct PreparedNoiseLoop {
    float diagonal;
    float z;
  };

  union PreparedWarpTransform {
    PreparedAffineFrame affine;
    PreparedMirrorTile mirror;
    PreparedVortex vortex;
    PreparedNoiseLoop noise_loop;
  };

  struct PreparedWarpStage {
    float rotation_cos;
    float rotation_sin;
    PreparedWarpTransform transform;
  };

  struct PreparedWarpProgram {
    PreparedWarpStage outer;
    PreparedWarpStage inner;
  };

  struct PreparedSurfaceNoise {
    Vector loop_offset;
    float direction_cos;
    float direction_sin;
  };

  struct PreparedHueRotation {
    static constexpr int VALUE_STEPS = 64;
    static constexpr int HUE_STEPS = 16;
    static constexpr size_t LUT_SIZE = VALUE_STEPS * HUE_STEPS;

    Pixel *lut;
    bool active;
  };

  struct PreparedHueNoise {
    static constexpr int FACE_STEPS = 24;
    static constexpr int FACE_SIZE = FACE_STEPS * FACE_STEPS;
    static constexpr size_t LUT_SIZE = 6 * FACE_SIZE;

    int8_t *lut;
    bool active;
  };

  struct ResourceBindings {
    const FastNoiseLite *outer_warp_noise;
    const FastNoiseLite *inner_warp_noise;
    const FastNoiseLite *source_noise;
    const FastNoiseLite *surface_noise;
    const FastNoiseLite *color_noise;
    const BakedPalette *generated_palette;
  };

  static constexpr size_t MAX_NOISE_RESOURCES = 9;

  struct FrameState {
    Slots slots;
    Params params;
    PaletteMappingWeights palette_mapping;
    ClockState clocks;
    SourceState prepared_source;
    PreparedTransforms transforms;
    PreparedWarpProgram prepared_warp;
    PreparedSurfaceNoise prepared_surface_noise;
    PreparedHueRotation prepared_hue_rotation;
    PreparedHueNoise prepared_hue_noise;
    ResourceBindings resources;
  };

  struct ShaderBallInstrumentation {
#ifdef HS_PROFILE_SHADERBALL_STAGES
    using Token = uint32_t;

    __attribute__((always_inline)) static Token mark() {
      return HS_OS_CYCLES();
    }

    template <Pullback::ProfileEvent Event>
    __attribute__((always_inline)) static void span(Token start) {
      const uint32_t elapsed = HS_OS_CYCLES() - start;
      if constexpr (Event == Pullback::ProfileEvent::LENS)
        hs::g_shaderball_stage_cycles.lens += elapsed;
      else if constexpr (Event == Pullback::ProfileEvent::SURFACE_NOISE)
        hs::g_shaderball_stage_cycles.surface_noise += elapsed;
      else if constexpr (Event == Pullback::ProfileEvent::PROJECTION)
        hs::g_shaderball_stage_cycles.projection += elapsed;
      else if constexpr (Event == Pullback::ProfileEvent::PLANAR_WARP)
        hs::g_shaderball_stage_cycles.planar_warp += elapsed;
      else if constexpr (Event == Pullback::ProfileEvent::MIRROR_TILE)
        hs::g_shaderball_stage_cycles.mirror_tile += elapsed;
      else if constexpr (Event == Pullback::ProfileEvent::SOURCE)
        hs::g_shaderball_stage_cycles.source += elapsed;
      else if constexpr (Event == Pullback::ProfileEvent::MATERIAL)
        hs::g_shaderball_stage_cycles.material += elapsed;
      else
        hs::g_shaderball_stage_cycles.color += elapsed;
    }
#else
    using Token = Pullback::NoInstrumentation::Token;

    __attribute__((always_inline)) static Token mark() { return {}; }

    template <Pullback::ProfileEvent Event>
    __attribute__((always_inline)) static void span(Token) {
      static_cast<void>(Event);
    }
#endif
  };

  struct ShaderBallBinding {
    using FrameState = ShaderBall::FrameState;
    using Instrumentation = ShaderBallInstrumentation;

    template <typename... Stages> struct ExtraValidation {
      using SurfaceStage = std::tuple_element_t<1, std::tuple<Stages...>>;
      using MaterialStage = std::tuple_element_t<4, std::tuple<Stages...>>;
      static constexpr bool value =
          !SurfaceStage::EDGE_DISTANCE_UNCONDITIONAL ||
          MaterialStage::COVERAGE == CoveragePolicy::EDGE_FADE;
    };
  };

  struct OuterCameraStateProvider {
    using Binding = ShaderBallBinding;
    using FrameState = typename Binding::FrameState;

    static const Quaternion &conjugate(const FrameState &frame) {
      return frame.transforms.outer_conj;
    }
  };

  struct ProjectionStateProvider {
    using Binding = ShaderBallBinding;
    using FrameState = typename Binding::FrameState;

    static const Quaternion &conjugate(const FrameState &frame) {
      return frame.transforms.projection_conj;
    }
    static float pole_fade(const FrameState &frame) {
      return frame.params.projection.pole_fade;
    }
    static float central_meridian(const FrameState &frame) {
      return frame.params.projection.central_meridian;
    }
    static float coordinate_scale(const FrameState &frame) {
      return frame.params.projection.coordinate_scale;
    }
    static float standard_parallel(const FrameState &frame) {
      return frame.params.projection.bonne_standard_parallel;
    }
    static float layout_scroll(const FrameState &frame) {
      return frame.params.projection.layout_scroll;
    }
  };

  struct SurfaceStateProvider {
    using Binding = ShaderBallBinding;
    using FrameState = typename Binding::FrameState;

    static const FastNoiseLite &noise(const FrameState &frame) {
      return *frame.resources.surface_noise;
    }
    static float scale(const FrameState &frame) {
      return frame.params.surface_noise.scale;
    }
    static const Vector &loop_offset(const FrameState &frame) {
      return frame.prepared_surface_noise.loop_offset;
    }
    static float strength(const FrameState &frame) {
      return frame.params.surface_noise.strength;
    }
    static float direction_cos(const FrameState &frame) {
      return frame.prepared_surface_noise.direction_cos;
    }
    static float direction_sin(const FrameState &frame) {
      return frame.prepared_surface_noise.direction_sin;
    }
    __attribute__((always_inline)) static bool
    path_length_required(const FrameState &frame) {
      return tracks_displacement(frame);
    }
  };

  struct LensStateProvider {
    using Binding = ShaderBallBinding;
    using FrameState = typename Binding::FrameState;

    static const MobiusParams &params(const FrameState &frame) {
      return frame.params.surface_lens.mobius;
    }
  };

  template <bool Outer> struct WarpStateProvider {
    using Binding = ShaderBallBinding;
    using FrameState = typename Binding::FrameState;

    static const WarpStageParams &params(const FrameState &frame) {
      if constexpr (Outer)
        return frame.params.warp.outer;
      else
        return frame.params.warp.inner;
    }
    static const PreparedWarpStage &prepared(const FrameState &frame) {
      if constexpr (Outer)
        return frame.prepared_warp.outer;
      else
        return frame.prepared_warp.inner;
    }
    static float phase(const FrameState &frame) {
      if constexpr (Outer)
        return frame.clocks.warp_outer_phase;
      else
        return frame.clocks.warp_inner_phase;
    }
    static const FastNoiseLite &noise(const FrameState &frame) {
      if constexpr (Outer)
        return *frame.resources.outer_warp_noise;
      else
        return *frame.resources.inner_warp_noise;
    }
    __attribute__((always_inline)) static bool
    path_length_required(const FrameState &frame) {
      return tracks_displacement(frame);
    }
  };

  struct SourceStateProvider {
    using Binding = ShaderBallBinding;
    using FrameState = typename Binding::FrameState;

    static const SourceParams &params(const FrameState &frame) {
      return frame.params.source;
    }
    static const SourceState &prepared(const FrameState &frame) {
      return frame.prepared_source;
    }
    static const FastNoiseLite &noise(const FrameState &frame) {
      return *frame.resources.source_noise;
    }
    static float noise_scale(const FrameState &frame) {
      return frame.params.source.noise_scale;
    }
    static float noise_time(const FrameState &frame) {
      return frame.clocks.source_noise_time;
    }
    static float noise_contrast(const FrameState &frame) {
      return frame.params.source.noise_contrast;
    }
  };

  struct ValueStateProvider {
    using Binding = ShaderBallBinding;
    using FrameState = typename Binding::FrameState;

    static float iso_level(const FrameState &frame) {
      return frame.params.value.iso_level;
    }
    static float iso_width(const FrameState &frame) {
      return frame.params.value.iso_width;
    }
    static float band_count(const FrameState &frame) {
      return static_cast<float>(frame.params.value.band_count);
    }
    static float band_phase(const FrameState &frame) {
      return frame.params.value.band_phase;
    }
    static float cutout_threshold(const FrameState &frame) {
      return frame.params.value.cutout_threshold;
    }
    static float cutout_width(const FrameState &frame) {
      return frame.params.value.cutout_softness;
    }
    static float edge_width(const FrameState &frame) {
      return frame.params.value.edge_width;
    }
  };

  static_assert(
      static_cast<uint8_t>(PaletteMapping::CUP) ==
          static_cast<uint8_t>(Pullback::Color::PaletteMapping::CUP) &&
      static_cast<uint8_t>(PaletteMapping::BELL) ==
          static_cast<uint8_t>(Pullback::Color::PaletteMapping::BELL) &&
      static_cast<uint8_t>(PaletteMapping::LINEAR) ==
          static_cast<uint8_t>(Pullback::Color::PaletteMapping::LINEAR) &&
      static_cast<uint8_t>(PaletteMapping::REVERSE) ==
          static_cast<uint8_t>(Pullback::Color::PaletteMapping::REVERSE));

  static constexpr PaletteMappingWeights
  palette_mapping_weights(PaletteMapping mapping) {
    return PaletteMappingWeights::single(
        static_cast<Pullback::Color::PaletteMapping>(mapping));
  }
  static_assert(
      static_cast<uint8_t>(BrightnessEnvelope::NONE) ==
          static_cast<uint8_t>(Pullback::Color::BrightnessEnvelope::NONE) &&
      static_cast<uint8_t>(BrightnessEnvelope::CUP) ==
          static_cast<uint8_t>(Pullback::Color::BrightnessEnvelope::CUP) &&
      static_cast<uint8_t>(BrightnessEnvelope::BELL) ==
          static_cast<uint8_t>(Pullback::Color::BrightnessEnvelope::BELL) &&
      static_cast<uint8_t>(BrightnessEnvelope::ASCENDING) ==
          static_cast<uint8_t>(
              Pullback::Color::BrightnessEnvelope::ASCENDING) &&
      static_cast<uint8_t>(BrightnessEnvelope::DESCENDING) ==
          static_cast<uint8_t>(
              Pullback::Color::BrightnessEnvelope::DESCENDING));
  static_assert(
      static_cast<uint8_t>(HueShiftMode::NONE) ==
          static_cast<uint8_t>(Pullback::Color::HueMode::NONE) &&
      static_cast<uint8_t>(HueShiftMode::NOISE) ==
          static_cast<uint8_t>(Pullback::Color::HueMode::NOISE) &&
      static_cast<uint8_t>(HueShiftMode::WARP_DISPLACEMENT) ==
          static_cast<uint8_t>(Pullback::Color::HueMode::PATH_LENGTH));

  struct ColorStateProvider {
    using Binding = ShaderBallBinding;
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
      return frame.clocks.palette_oscillation_phase;
    }
    static Color4 palette(const FrameState &frame, float value) {
      return frame.resources.generated_palette->get(value);
    }
    static Pullback::Color::HueMode hue_mode(const FrameState &frame) {
      return static_cast<Pullback::Color::HueMode>(frame.slots.hue_shift);
    }
    static float hue_shift_amount(const FrameState &frame) {
      return frame.params.color.hue_shift_amount;
    }
    static Pullback::Color::HueRotationLutView
    hue_rotation(const FrameState &frame) {
      return {frame.prepared_hue_rotation.lut,
              frame.prepared_hue_rotation.active};
    }
    static Pullback::Color::HueNoiseLutView hue_noise(const FrameState &frame) {
      return {frame.prepared_hue_noise.lut, frame.prepared_hue_noise.active};
    }
    static Pullback::Color::BrightnessEnvelope
    brightness_envelope(const FrameState &frame) {
      return static_cast<Pullback::Color::BrightnessEnvelope>(
          frame.slots.brightness_envelope);
    }
    static float brightness_depth(const FrameState &frame) {
      return frame.params.color.brightness_depth;
    }
    static float opacity_low(const FrameState &frame) {
      return frame.params.color.value_opacity_low;
    }
    static float opacity_high(const FrameState &frame) {
      return frame.params.color.value_opacity_high;
    }
  };

  using SourceInput = Pullback::SourceInput;
  using MaterialInput = Pullback::MaterialInput;
  using InverseStageKind = Pullback::StageKind;
  using CodeEmission = Pullback::CodeEmission;
  using ApproximationOracleId = Pullback::ApproximationOracleId;
  using ApproximationDomain = Pullback::ApproximationDomain;
  using ApproximationAggregation = Pullback::ApproximationAggregation;
  using ApproximationMetric = Pullback::ApproximationMetric;

  struct TopologyKey {
    Function function{};
    Projection projection{};
    ProjectionFramePolicy projection_frame{};
    SurfaceLens surface_lens{};
    SignalWeight signal_weight{};
    ValueTransfer value_transfer{};
    CoveragePolicy coverage{};
    PeirceLayout peirce_layout{};
    AiroceanLayout airocean_layout{};
    BonneHemisphere bonne_hemisphere{};
    GnomonicHemispherePolicy gnomonic_hemisphere{};
    SurfaceNoise surface_noise{};
    SurfaceNoisePlacement surface_noise_placement{};
    NoiseBasis surface_noise_basis{};
    SurfaceCurlIntegrator surface_curl_integrator{};
    NoiseBasis source_noise_basis{};
    WarpStageKind outer_warp{};
    NoiseBasis outer_warp_basis{};
    WarpEnvelope outer_warp_envelope{};
    PolarMode outer_polar_mode{};
    CurlIntegrator outer_curl_integrator{};
    uint8_t outer_polar_harmonic{};
    WarpStageKind inner_warp{};
    NoiseBasis inner_warp_basis{};
    WarpEnvelope inner_warp_envelope{};
    PolarMode inner_polar_mode{};
    CurlIntegrator inner_curl_integrator{};
    uint8_t inner_polar_harmonic{};

    constexpr bool operator==(const TopologyKey &) const = default;
  };

  enum class InversePipelineId : uint8_t {
    GLITCH_NOISE_GRID_WAVE_SHEAR,
    KALEIDOSCOPE_TWIN_WAVE_INNER_MIRROR,
    GNOMONIC_KALEIDOSCOPE_GRID_MIRROR,
    GNOMONIC_GLITCH_GRID_MIRROR,
    PEIRCE_DODECAHEDRAL_GRID,
    GNOMONIC_DODECAHEDRAL_GRID_WAVE_MIRROR,
    GNOMONIC_AFFINE_LATTICE_CONTOUR,
    SINUSOIDAL_CURL_LATTICE,
    STEREOGRAPHIC_PRISM_POLAR_WAVE_LATTICE,
    GNOMONIC_DODECAHEDRAL_GRID_VECTOR_MIRROR,
    STEREOGRAPHIC_DODECAHEDRAL_GRID_INNER_MIRROR,
    STEREOGRAPHIC_HEXAGONAL_PRISM_TWIN_WAVE_INNER_MIRROR,
    EQUIRECTANGULAR_DODECAHEDRAL_GRID_INNER_MIRROR,
    STEREOGRAPHIC_GLITCH_GRID_MIRROR,
    STEREOGRAPHIC_MOBIUS_TWIN_WAVE_INNER_MIRROR,
    COUNT,
    NONE = 0xff
  };

  struct SelectedConfig {
    Config config;
    InversePipelineId pipeline;
  };
  using Preset = SelectedConfig;

  size_t preset_count_for_view() const {
    return preset_view.empty() ? PRESETS.size() : preset_view.size();
  }

  const Preset &preset_for_view(size_t index) const {
    assert(index < preset_count_for_view());
    return PRESETS[preset_view.empty() ? index : preset_view[index]];
  }

  template <typename Stage>
  using HasInverseStageContract =
      Pullback::HasStageContract<Stage, ShaderBallBinding>;

  template <typename Stage>
  using HasTypedInverseStageContract =
      Pullback::HasTypedStageContract<Stage, ShaderBallBinding>;

  template <bool LevelOneValid, typename... Stages>
  using InversePipelineValidation =
      Pullback::PipelineValidationLevel2<LevelOneValid, ShaderBallBinding,
                                         Stages...>;

  template <typename... Stages>
  using InversePipeline = Pullback::Pipeline<ShaderBallBinding, Stages...>;

#if defined(__IMXRT1062__)
  HS_FLASH_MEMBER
#else
  __attribute__((always_inline))
#endif
  static Vector pullback_outer_camera_lookup(const Vector &input,
                                             const FrameState &frame) {
    return rotate(input, frame.transforms.outer_conj);
  }

  struct OuterCameraStage
      : Pullback::Stage::OuterCamera<ShaderBallBinding,
                                     OuterCameraStateProvider> {
    static constexpr bool implements(const TopologyKey &) { return true; }

    __attribute__((always_inline)) static Vector run(const Vector &input,
                                                     const FrameState &frame) {
      return pullback_outer_camera_lookup(input, frame);
    }
  };

  template <SurfaceLens LensV>
  using LensPolicy = std::conditional_t<
      LensV == SurfaceLens::NONE, Pullback::Lens::Identity,
      std::conditional_t<
          LensV == SurfaceLens::GLITCH, Pullback::Lens::Glitch,
          std::conditional_t<
              LensV == SurfaceLens::KALEIDOSCOPE, Pullback::Lens::Kaleidoscope,
              std::conditional_t<
                  LensV == SurfaceLens::MOBIUS,
                  Pullback::Lens::Mobius<LensStateProvider>,
                  std::conditional_t<
                      LensV == SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL,
                      Pullback::Lens::DodecahedralKaleidoscope,
                      std::conditional_t<
                          LensV == SurfaceLens::KALEIDOSCOPE_HEXAGONAL_PRISM,
                          Pullback::Lens::HexagonalPrismKaleidoscope,
                          Pullback::Lens::TriangularPrismKaleidoscope>>>>>>;

  template <Projection ProjectionV>
  using ProjectionPolicy = std::conditional_t<
      ProjectionV == Projection::STEREOGRAPHIC,
      Pullback::Projection::Stereographic<ProjectionStateProvider>,
      std::conditional_t<
          ProjectionV == Projection::GNOMONIC,
          Pullback::Projection::Gnomonic<
              ProjectionStateProvider,
              Pullback::Projection::GnomonicHemisphere::FOLDED>,
          std::conditional_t<
              ProjectionV == Projection::BONNE,
              Pullback::Projection::Bonne<ProjectionStateProvider, true>,
              std::conditional_t<ProjectionV == Projection::EQUIRECTANGULAR,
                                 Pullback::Projection::Equirectangular<
                                     ProjectionStateProvider>,
                                 Pullback::Projection::PeirceSquare<
                                     ProjectionStateProvider>>>>>;

  template <Projection ProjectionV, SurfaceLens LensV>
  struct SelectedSurfaceProjectStage
      : Pullback::Stage::SurfaceProject<
            ShaderBallBinding, Pullback::Surface::Identity, LensPolicy<LensV>,
            Pullback::Surface::Identity, ProjectionPolicy<ProjectionV>> {
    static constexpr bool implements(const TopologyKey &key) {
      return key.projection == ProjectionV && key.surface_lens == LensV &&
             key.surface_noise == SurfaceNoise::NONE &&
             (ProjectionV != Projection::GNOMONIC ||
              key.gnomonic_hemisphere == GnomonicHemispherePolicy::FOLDED) &&
             (ProjectionV != Projection::PEIRCE_QUINCUNCIAL ||
              key.peirce_layout == PeirceLayout::SQUARE) &&
             (ProjectionV != Projection::BONNE ||
              key.bonne_hemisphere == BonneHemisphere::NORTH);
    }
  };

  using SinusoidalCurlSurfaceImplementation = Pullback::Stage::SurfaceProject<
      ShaderBallBinding,
      Pullback::Surface::CurlNoise<SurfaceStateProvider, NoiseBasis::SIMPLEX,
                                   Pullback::Surface::Euler>,
      Pullback::Lens::Identity, Pullback::Surface::Identity,
      Pullback::Projection::FoldedSinusoidal<ProjectionStateProvider>>;
  struct SinusoidalCurlSurfaceStage
      : Pullback::Stage::Placed<Pullback::CodeEmission::OUT_OF_LINE_FLASH,
                                SinusoidalCurlSurfaceImplementation> {
    static constexpr bool implements(const TopologyKey &key) {
      return key.projection == Projection::SINUSOIDAL &&
             key.surface_lens == SurfaceLens::NONE &&
             key.surface_noise == SurfaceNoise::CURL &&
             key.surface_noise_placement ==
                 SurfaceNoisePlacement::BEFORE_LENS &&
             key.surface_noise_basis == NoiseBasis::SIMPLEX &&
             key.surface_curl_integrator == SurfaceCurlIntegrator::EULER;
    }
  };

  template <WarpStageKind KindV, bool Outer>
  using WarpPolicy = std::conditional_t<
      KindV == WarpStageKind::NONE, Pullback::Warp::Identity,
      std::conditional_t<
          KindV == WarpStageKind::AFFINE_FRAME,
          Pullback::Warp::AffineFrame<WarpStateProvider<Outer>>,
          std::conditional_t<
              KindV == WarpStageKind::WAVE_SHEAR,
              Pullback::Warp::WaveShear<WarpStateProvider<Outer>>,
              std::conditional_t<
                  KindV == WarpStageKind::VECTOR_NOISE,
                  Pullback::Warp::VectorNoise<WarpStateProvider<Outer>,
                                              NoiseBasis::SIMPLEX,
                                              Pullback::Warp::FlatEnvelope>,
                  std::conditional_t<
                      KindV == WarpStageKind::MIRROR_TILE,
                      Pullback::Warp::MirrorTile<WarpStateProvider<Outer>>,
                      Pullback::Warp::PolarChart<WarpStateProvider<Outer>,
                                                 Pullback::Warp::LinearPolar,
                                                 1>>>>>>;

  template <WarpStageKind Outer, WarpStageKind Inner>
  struct SelectedPlanarWarpStage
      : Pullback::Stage::PlanarWarp<ShaderBallBinding, WarpPolicy<Outer, true>,
                                    WarpPolicy<Inner, false>> {
    static constexpr bool implements(const TopologyKey &key) {
      static_assert(Outer == WarpStageKind::NONE ||
                    Outer == WarpStageKind::AFFINE_FRAME ||
                    Outer == WarpStageKind::WAVE_SHEAR ||
                    Outer == WarpStageKind::VECTOR_NOISE ||
                    Outer == WarpStageKind::MIRROR_TILE ||
                    Outer == WarpStageKind::POLAR_CHART);
      static_assert(Inner == WarpStageKind::NONE ||
                    Inner == WarpStageKind::WAVE_SHEAR ||
                    Inner == WarpStageKind::MIRROR_TILE);
      return key.outer_warp == Outer && key.inner_warp == Inner &&
             (Outer != WarpStageKind::WAVE_SHEAR ||
              key.outer_warp_envelope == WarpEnvelope::FLAT) &&
             (Inner != WarpStageKind::WAVE_SHEAR ||
              key.inner_warp_envelope == WarpEnvelope::FLAT);
    }
  };

  template <Function FunctionV>
  using SourcePolicy = std::conditional_t<
      FunctionV == Function::GRID, Pullback::Source::Grid<SourceStateProvider>,
      std::conditional_t<
          FunctionV == Function::PRIMITIVE_LATTICE,
          Pullback::Source::PrimitiveLattice<SourceStateProvider>,
          Pullback::Source::TwinWave<SourceStateProvider>>>;

  template <Function FunctionV>
  struct SourceStage
      : Pullback::Stage::Source<ShaderBallBinding, SourcePolicy<FunctionV>> {
    static constexpr bool implements(const TopologyKey &key) {
      return key.function == FunctionV;
    }
  };

  template <CoveragePolicy CoverageV>
  using CoveragePolicyFor = std::conditional_t<
      CoverageV == CoveragePolicy::OPAQUE, Pullback::Coverage::Opaque,
      std::conditional_t<
          CoverageV == CoveragePolicy::EDGE_FADE,
          Pullback::Coverage::EdgeFade<ValueStateProvider>,
          std::conditional_t<CoverageV == CoveragePolicy::PROJECTION_WEIGHT,
                             Pullback::Coverage::Projection,
                             Pullback::Coverage::ProjectionSquared>>>;

  template <CoveragePolicy CoverageV>
  struct LinearMaterialStage
      : Pullback::Stage::Material<
            ShaderBallBinding, Pullback::Weight::Projection,
            Pullback::Transfer::Linear, CoveragePolicyFor<CoverageV>> {
    static_assert(CoverageV == CoveragePolicy::OPAQUE ||
                  CoverageV == CoveragePolicy::EDGE_FADE ||
                  CoverageV == CoveragePolicy::PROJECTION_WEIGHT_SQUARED ||
                  CoverageV == CoveragePolicy::PROJECTION_WEIGHT);
    static constexpr CoveragePolicy COVERAGE = CoverageV;
    static constexpr bool implements(const TopologyKey &key) {
      return key.signal_weight == SignalWeight::PROJECTION &&
             key.value_transfer == ValueTransfer::LINEAR &&
             key.coverage == CoverageV;
    }
  };

  struct IsoContourProjectionWeightMaterialStage
      : Pullback::Stage::Material<
            ShaderBallBinding, Pullback::Weight::Projection,
            Pullback::Transfer::IsoContour<ValueStateProvider>,
            Pullback::Coverage::Projection> {
    static constexpr CoveragePolicy COVERAGE =
        CoveragePolicy::PROJECTION_WEIGHT;
    static constexpr bool implements(const TopologyKey &key) {
      return key.signal_weight == SignalWeight::PROJECTION &&
             key.value_transfer == ValueTransfer::ISO_CONTOUR &&
             key.coverage == CoveragePolicy::PROJECTION_WEIGHT;
    }
  };

  struct ColorStage
      : Pullback::Stage::Color<
            ShaderBallBinding,
            Pullback::Color::GeneratedPalette<ColorStateProvider>> {
    static constexpr bool implements(const TopologyKey &) { return true; }
  };

  using GlitchNoiseGridWaveShearPipelineBase = InversePipeline<
      OuterCameraStage,
      SelectedSurfaceProjectStage<Projection::STEREOGRAPHIC,
                                  SurfaceLens::GLITCH>,
      SelectedPlanarWarpStage<WarpStageKind::WAVE_SHEAR, WarpStageKind::NONE>,
      SourceStage<Function::GRID>,
      LinearMaterialStage<CoveragePolicy::PROJECTION_WEIGHT_SQUARED>,
      ColorStage>;
  struct GlitchNoiseGridWaveShearPipeline
      : GlitchNoiseGridWaveShearPipelineBase {
    HS_HOT_FLASH_MEMBER static Color4 shade(const Vector &view,
                                            const FrameState &frame) {
      return GlitchNoiseGridWaveShearPipelineBase::template run_stage<0>(view,
                                                                         frame);
    }
  };
  using KaleidoscopeTwinWaveInnerMirrorPipeline = InversePipeline<
      OuterCameraStage,
      SelectedSurfaceProjectStage<Projection::STEREOGRAPHIC,
                                  SurfaceLens::KALEIDOSCOPE>,
      SelectedPlanarWarpStage<WarpStageKind::NONE, WarpStageKind::MIRROR_TILE>,
      SourceStage<Function::TWIN_WAVE>,
      LinearMaterialStage<CoveragePolicy::PROJECTION_WEIGHT_SQUARED>,
      ColorStage>;
  using StereographicMobiusTwinWaveInnerMirrorPipeline = InversePipeline<
      OuterCameraStage,
      SelectedSurfaceProjectStage<Projection::STEREOGRAPHIC,
                                  SurfaceLens::MOBIUS>,
      SelectedPlanarWarpStage<WarpStageKind::NONE, WarpStageKind::MIRROR_TILE>,
      SourceStage<Function::TWIN_WAVE>,
      LinearMaterialStage<CoveragePolicy::PROJECTION_WEIGHT_SQUARED>,
      ColorStage>;
  using StereographicHexagonalPrismTwinWaveInnerMirrorPipeline =
      InversePipeline<
          OuterCameraStage,
          SelectedSurfaceProjectStage<
              Projection::STEREOGRAPHIC,
              SurfaceLens::KALEIDOSCOPE_HEXAGONAL_PRISM>,
          SelectedPlanarWarpStage<WarpStageKind::NONE,
                                  WarpStageKind::MIRROR_TILE>,
          SourceStage<Function::TWIN_WAVE>,
          LinearMaterialStage<CoveragePolicy::PROJECTION_WEIGHT_SQUARED>,
          ColorStage>;
  using GnomonicKaleidoscopeGridMirrorPipeline = InversePipeline<
      OuterCameraStage,
      SelectedSurfaceProjectStage<Projection::GNOMONIC,
                                  SurfaceLens::KALEIDOSCOPE>,
      SelectedPlanarWarpStage<WarpStageKind::MIRROR_TILE, WarpStageKind::NONE>,
      SourceStage<Function::GRID>,
      LinearMaterialStage<CoveragePolicy::EDGE_FADE>, ColorStage>;
  using GnomonicGlitchGridMirrorPipeline = InversePipeline<
      OuterCameraStage,
      SelectedSurfaceProjectStage<Projection::GNOMONIC, SurfaceLens::GLITCH>,
      SelectedPlanarWarpStage<WarpStageKind::MIRROR_TILE, WarpStageKind::NONE>,
      SourceStage<Function::GRID>,
      LinearMaterialStage<CoveragePolicy::EDGE_FADE>, ColorStage>;
  using PeirceDodecahedralGridPipelineBase = InversePipeline<
      OuterCameraStage,
      SelectedSurfaceProjectStage<Projection::PEIRCE_QUINCUNCIAL,
                                  SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL>,
      SelectedPlanarWarpStage<WarpStageKind::NONE, WarpStageKind::NONE>,
      SourceStage<Function::GRID>,
      LinearMaterialStage<CoveragePolicy::EDGE_FADE>, ColorStage>;
  struct PeirceDodecahedralGridPipeline : PeirceDodecahedralGridPipelineBase {
    HS_HOT_FLASH_MEMBER static Color4 shade(const Vector &view,
                                            const FrameState &frame) {
      return PeirceDodecahedralGridPipelineBase::template run_stage<0>(view,
                                                                       frame);
    }
  };
  using GnomonicDodecahedralGridWaveMirrorPipelineBase = InversePipeline<
      OuterCameraStage,
      SelectedSurfaceProjectStage<Projection::GNOMONIC,
                                  SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL>,
      SelectedPlanarWarpStage<WarpStageKind::WAVE_SHEAR,
                              WarpStageKind::MIRROR_TILE>,
      SourceStage<Function::GRID>,
      LinearMaterialStage<CoveragePolicy::PROJECTION_WEIGHT_SQUARED>,
      ColorStage>;
  struct GnomonicDodecahedralGridWaveMirrorPipeline
      : GnomonicDodecahedralGridWaveMirrorPipelineBase {
    HS_HOT_FLASH_MEMBER static Color4 shade(const Vector &view,
                                            const FrameState &frame) {
      return GnomonicDodecahedralGridWaveMirrorPipelineBase::template run_stage<
          0>(view, frame);
    }
  };
  using GnomonicDodecahedralGridVectorMirrorPipelineBase = InversePipeline<
      OuterCameraStage,
      SelectedSurfaceProjectStage<Projection::GNOMONIC,
                                  SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL>,
      SelectedPlanarWarpStage<WarpStageKind::VECTOR_NOISE,
                              WarpStageKind::MIRROR_TILE>,
      SourceStage<Function::GRID>,
      LinearMaterialStage<CoveragePolicy::PROJECTION_WEIGHT_SQUARED>,
      ColorStage>;
  struct GnomonicDodecahedralGridVectorMirrorPipeline
      : GnomonicDodecahedralGridVectorMirrorPipelineBase {
    HS_FLASH_MEMBER __attribute__((noinline)) static Color4
    shade(const Vector &view, const FrameState &frame) {
      return GnomonicDodecahedralGridVectorMirrorPipelineBase::
          template run_stage<0>(view, frame);
    }
  };
  using GnomonicAffineLatticeContourPipeline = InversePipeline<
      OuterCameraStage,
      SelectedSurfaceProjectStage<Projection::GNOMONIC, SurfaceLens::NONE>,
      SelectedPlanarWarpStage<WarpStageKind::AFFINE_FRAME, WarpStageKind::NONE>,
      SourceStage<Function::PRIMITIVE_LATTICE>,
      IsoContourProjectionWeightMaterialStage, ColorStage>;
  using SinusoidalCurlLatticePipeline = InversePipeline<
      OuterCameraStage, SinusoidalCurlSurfaceStage,
      SelectedPlanarWarpStage<WarpStageKind::NONE, WarpStageKind::NONE>,
      SourceStage<Function::PRIMITIVE_LATTICE>,
      LinearMaterialStage<CoveragePolicy::PROJECTION_WEIGHT>, ColorStage>;
  using StereographicPrismPolarWaveLatticePipeline = InversePipeline<
      OuterCameraStage,
      SelectedSurfaceProjectStage<Projection::STEREOGRAPHIC,
                                  SurfaceLens::KALEIDOSCOPE_TRIANGULAR_PRISM>,
      SelectedPlanarWarpStage<WarpStageKind::POLAR_CHART,
                              WarpStageKind::WAVE_SHEAR>,
      SourceStage<Function::PRIMITIVE_LATTICE>,
      LinearMaterialStage<CoveragePolicy::PROJECTION_WEIGHT_SQUARED>,
      ColorStage>;
  using StereographicDodecahedralGridInnerMirrorPipeline = InversePipeline<
      OuterCameraStage,
      SelectedSurfaceProjectStage<Projection::STEREOGRAPHIC,
                                  SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL>,
      SelectedPlanarWarpStage<WarpStageKind::NONE, WarpStageKind::MIRROR_TILE>,
      SourceStage<Function::GRID>,
      LinearMaterialStage<CoveragePolicy::PROJECTION_WEIGHT_SQUARED>,
      ColorStage>;
  using EquirectangularDodecahedralGridInnerMirrorPipeline = InversePipeline<
      OuterCameraStage,
      Pullback::Stage::Placed<
          Pullback::CodeEmission::OUT_OF_LINE_FLASH,
          SelectedSurfaceProjectStage<Projection::EQUIRECTANGULAR,
                                      SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL>>,
      SelectedPlanarWarpStage<WarpStageKind::NONE, WarpStageKind::MIRROR_TILE>,
      SourceStage<Function::GRID>,
      LinearMaterialStage<CoveragePolicy::PROJECTION_WEIGHT_SQUARED>,
      ColorStage>;
  using StereographicGlitchGridMirrorPipeline = InversePipeline<
      OuterCameraStage,
      SelectedSurfaceProjectStage<Projection::STEREOGRAPHIC,
                                  SurfaceLens::GLITCH>,
      SelectedPlanarWarpStage<WarpStageKind::MIRROR_TILE, WarpStageKind::NONE>,
      SourceStage<Function::GRID>,
      LinearMaterialStage<CoveragePolicy::EDGE_FADE>, ColorStage>;

  struct ProgramDescriptor {
    InversePipelineId id;
    TopologyKey key;
    ShadeFunction shade;
    bool (*continuous_parameters_supported)(const Config &);
    bool (*resources_ready)(const FrameState &);
  };

  static constexpr bool source_uses_noise(Function function) {
    return function == Function::NOISE_CONTOUR ||
           function == Function::NOISE_CONTOUR_SPHERE;
  }

  /** @brief Whether the stage scales its amplitude by the warp envelope. */
  static constexpr bool warp_uses_envelope(WarpStageKind kind) {
    return kind == WarpStageKind::WAVE_SHEAR ||
           kind == WarpStageKind::VECTOR_NOISE ||
           kind == WarpStageKind::CURL_FLOW;
  }

  static constexpr void canonicalize_warp_key(WarpStageKind kind,
                                              NoiseBasis &basis,
                                              WarpEnvelope &envelope,
                                              PolarMode &polar_mode,
                                              CurlIntegrator &curl_integrator,
                                              uint8_t &polar_harmonic) {
    if (!warp_uses_noise(kind))
      basis = {};
    if (!warp_uses_envelope(kind))
      envelope = {};
    if (kind != WarpStageKind::CURL_FLOW)
      curl_integrator = {};
    if (kind != WarpStageKind::POLAR_CHART) {
      polar_mode = {};
      polar_harmonic = 0;
    }
  }

  static constexpr TopologyKey make_topology_key(const Config &config) {
    const Slots &slots = config.slots;
    TopologyKey key{
        slots.function,
        slots.projection,
        slots.projection_frame,
        slots.surface_lens,
        slots.signal_weight,
        slots.value_transfer,
        slots.coverage,
        slots.peirce_layout,
        slots.airocean_layout,
        slots.bonne_hemisphere,
        slots.gnomonic_hemisphere,
        slots.surface_noise,
        slots.surface_noise_placement,
        config.params.surface_noise.basis,
        config.params.surface_noise.integrator,
        config.params.source.noise_basis,
        slots.warp_program.outer.kind,
        slots.warp_program.outer.basis,
        slots.warp_program.outer.envelope,
        slots.warp_program.outer.polar_mode,
        slots.warp_program.outer.curl_integrator,
        slots.warp_program.outer.polar_harmonic,
        slots.warp_program.inner.kind,
        slots.warp_program.inner.basis,
        slots.warp_program.inner.envelope,
        slots.warp_program.inner.polar_mode,
        slots.warp_program.inner.curl_integrator,
        slots.warp_program.inner.polar_harmonic,
    };
    if (key.projection != Projection::PEIRCE_QUINCUNCIAL)
      key.peirce_layout = {};
    if (key.projection != Projection::AIROCEAN)
      key.airocean_layout = {};
    if (key.projection != Projection::BONNE)
      key.bonne_hemisphere = {};
    if (key.projection != Projection::GNOMONIC)
      key.gnomonic_hemisphere = {};
    if (key.surface_noise == SurfaceNoise::NONE) {
      key.surface_noise_placement = {};
      key.surface_noise_basis = {};
      key.surface_curl_integrator = {};
    }
    if (!source_uses_noise(key.function))
      key.source_noise_basis = {};
    canonicalize_warp_key(key.outer_warp, key.outer_warp_basis,
                          key.outer_warp_envelope, key.outer_polar_mode,
                          key.outer_curl_integrator, key.outer_polar_harmonic);
    canonicalize_warp_key(key.inner_warp, key.inner_warp_basis,
                          key.inner_warp_envelope, key.inner_polar_mode,
                          key.inner_curl_integrator, key.inner_polar_harmonic);
    return key;
  }

  static constexpr bool all_continuous_parameters_supported(const Config &) {
    return true;
  }

  HS_FLASH_MEMBER static bool
  pipeline_resources_ready(const FrameState &frame) {
    if (warp_uses_noise(frame.slots.warp_program.outer.kind) &&
        frame.resources.outer_warp_noise == nullptr)
      return false;
    if (warp_uses_noise(frame.slots.warp_program.inner.kind) &&
        frame.resources.inner_warp_noise == nullptr)
      return false;
    if (is_noise_contour(frame.slots.function) &&
        frame.resources.source_noise == nullptr)
      return false;
    if (frame.slots.surface_noise != SurfaceNoise::NONE &&
        frame.resources.surface_noise == nullptr)
      return false;
    if (frame.resources.generated_palette == nullptr)
      return false;
    if (frame.slots.hue_shift == HueShiftMode::NOISE &&
        frame.params.color.hue_shift_amount != 0.0f &&
        frame.resources.color_noise == nullptr)
      return false;
    if (frame.prepared_hue_rotation.active &&
        frame.prepared_hue_rotation.lut == nullptr)
      return false;
    return !frame.prepared_hue_noise.active ||
           frame.prepared_hue_noise.lut != nullptr;
  }

  struct PreparedEndpoint {
    FrameState *frame;
    ShadeFunction shade;
    InversePipelineId pipeline;
    float alpha;
#if defined(HS_PROFILE_ENABLE)
    size_t preset;
#endif
  };

  enum class ProfileEndpoint : uint8_t { STEADY, FROM, TO };

  struct LookRuntime {
    ClockState clocks{};
    Quaternion projection_wander;
    Quaternion outer_wander;
    PreparedTransforms transforms;

    HS_COLD_MEMBER LookRuntime() = default;
  };

  template <typename T> static uint32_t encode_field_value(const T &value) {
    static_assert(sizeof(T) <= sizeof(uint32_t));
    uint32_t payload = 0;
    std::memcpy(&payload, &value, sizeof(T));
    return payload;
  }

  template <typename T>
  static bool decode_field_value(uint32_t payload, T &value) {
    static_assert(sizeof(T) <= sizeof(uint32_t));
    if constexpr (sizeof(T) < sizeof(uint32_t)) {
      const uint32_t value_mask = (uint32_t{1} << (sizeof(T) * 8)) - 1;
      if ((payload & ~value_mask) != 0)
        return false;
    }
    if constexpr (std::is_same_v<T, bool>)
      if (payload > 1)
        return false;
    std::memcpy(&value, &payload, sizeof(T));
    return true;
  }

  static ConfigValues encode_config_values(const Config &config) {
    ConfigValues values{};
#define HS_SHADERBALL_ENCODE_FIELD(name, path)                                 \
  values[static_cast<size_t>(ConfigFieldId::name)] =                           \
      encode_field_value(config.path);
    HS_SHADERBALL_CONFIG_FIELDS(HS_SHADERBALL_ENCODE_FIELD)
#undef HS_SHADERBALL_ENCODE_FIELD
    values[static_cast<size_t>(ConfigFieldId::SLOTS_SURFACE_LENS)] =
        surface_lens_storage_id(config.slots.surface_lens);
    values[static_cast<size_t>(ConfigFieldId::SLOTS_WARP_OUTER_KIND)] =
        warp_storage_id(config.slots.warp_program.outer.kind);
    values[static_cast<size_t>(ConfigFieldId::SLOTS_WARP_INNER_KIND)] =
        warp_storage_id(config.slots.warp_program.inner.kind);
    return values;
  }

  static bool decode_config_values(const ConfigValues &values, Config &config) {
    bool valid = true;
#define HS_SHADERBALL_DECODE_FIELD(name, path)                                 \
  valid = decode_field_value(values[static_cast<size_t>(ConfigFieldId::name)], \
                             config.path) &&                                   \
          valid;
    HS_SHADERBALL_CONFIG_FIELDS(HS_SHADERBALL_DECODE_FIELD)
#undef HS_SHADERBALL_DECODE_FIELD
    valid = decode_surface_lens_storage(
                values[static_cast<size_t>(ConfigFieldId::SLOTS_SURFACE_LENS)],
                config.slots.surface_lens) &&
            valid;
    valid =
        decode_warp_storage(
            values[static_cast<size_t>(ConfigFieldId::SLOTS_WARP_OUTER_KIND)],
            config.slots.warp_program.outer.kind) &&
        valid;
    valid =
        decode_warp_storage(
            values[static_cast<size_t>(ConfigFieldId::SLOTS_WARP_INNER_KIND)],
            config.slots.warp_program.inner.kind) &&
        valid;
    return valid;
  }

  static constexpr uint32_t surface_lens_storage_id(SurfaceLens lens) {
    const uint8_t value = static_cast<uint8_t>(lens);
    if (lens == SurfaceLens::TANGENT_NOISE)
      return 5;
    return value < 5 ? value : value + 1;
  }

  static bool decode_surface_lens_storage(uint32_t id, SurfaceLens &lens) {
    if (id <= 4) {
      lens = static_cast<SurfaceLens>(id);
      return true;
    }
    if (id == 5) {
      lens = SurfaceLens::TANGENT_NOISE;
      return true;
    }
    if (id <= 13) {
      lens = static_cast<SurfaceLens>(id - 1);
      return true;
    }
    return false;
  }

  static constexpr uint32_t warp_storage_id(WarpStageKind kind) {
    if (kind == WarpStageKind::LEGACY_STEREO_NOISE)
      return 1;
    const uint8_t value = static_cast<uint8_t>(kind);
    return value == 0 ? 0 : value + 1;
  }

  static bool decode_warp_storage(uint32_t id, WarpStageKind &kind) {
    if (id == 0) {
      kind = WarpStageKind::NONE;
      return true;
    }
    if (id == 1) {
      kind = WarpStageKind::LEGACY_STEREO_NOISE;
      return true;
    }
    if (id <= 8) {
      kind = static_cast<WarpStageKind>(id - 1);
      return true;
    }
    return false;
  }

  /** @brief Tests every slot against the last enumerator its schema admits. */
  static constexpr bool valid_slot_enums(const Slots &slots) {
    return enum_at_most(slots.function, Function::NOISE_CONTOUR_SPHERE) &&
           enum_at_most(slots.projection, Projection::EQUIRECTANGULAR) &&
           enum_at_most(slots.projection_frame,
                        ProjectionFramePolicy::SPIN_WANDER) &&
           enum_at_most(slots.surface_lens,
                        SurfaceLens::KALEIDOSCOPE_OCTAGONAL_PRISM) &&
           enum_at_most(slots.surface_noise, SurfaceNoise::CURL) &&
           enum_at_most(slots.surface_noise_placement,
                        SurfaceNoisePlacement::AFTER_LENS) &&
           enum_at_most(slots.warp_program.outer.kind,
                        WarpStageKind::POLAR_CHART) &&
           enum_at_most(slots.warp_program.inner.kind,
                        WarpStageKind::POLAR_CHART) &&
           valid_warp_spec(slots.warp_program.outer) &&
           valid_warp_spec(slots.warp_program.inner) &&
           enum_at_most(slots.signal_weight, SignalWeight::PROJECTION) &&
           enum_at_most(slots.value_transfer, ValueTransfer::SMOOTH_BANDS) &&
           enum_at_most(slots.coverage, CoveragePolicy::PROJECTION_WEIGHT) &&
           enum_at_most(slots.palette, PaletteMode::ANALOGOUS) &&
           enum_at_most(slots.palette_mapping, PaletteMapping::REVERSE) &&
           enum_at_most(slots.brightness_envelope,
                        BrightnessEnvelope::DESCENDING) &&
           enum_at_most(slots.hue_shift, HueShiftMode::WARP_DISPLACEMENT) &&
           enum_at_most(slots.peirce_layout, PeirceLayout::VERTICAL) &&
           enum_at_most(slots.airocean_layout, AiroceanLayout::HORIZONTAL) &&
           enum_at_most(slots.bonne_hemisphere, BonneHemisphere::SOUTH) &&
           enum_at_most(slots.gnomonic_hemisphere,
                        GnomonicHemispherePolicy::BACK_HEMISPHERE);
  }

  static constexpr bool valid_snapshot_config(const Config &config) {
    return valid_slot_enums(config.slots) && preset_in_ranges(config.params) &&
           hue_shift_amount_in_range(config);
  }

#if HS_ENABLE_PARAM_GUI_BRIDGE
public:
  /** @brief Captures all accepted, requested, pending, and runtime state. */
  HS_COLD_MEMBER FullConfigSnapshot capture_full_config_snapshot() const {
    FullConfigSnapshot snapshot;
    snapshot.accepted = encode_config_values(accepted_config);
    snapshot.requested = encode_config_values(requested_config);
    for (size_t index = 0; index < pending_edit_count; ++index)
      snapshot.pending[static_cast<size_t>(pending_edits[index].id)] = 1;
    snapshot.has_runtime = true;
    const ClockState &clocks = runtime.clocks;
    snapshot.runtime = {
        clocks.source_primary,     clocks.source_secondary,
        clocks.source_angle,       clocks.warp_outer_rotation,
        clocks.projection_spin,    clocks.hue_noise_phase,
        clocks.source_noise_time,  clocks.warp_inner_rotation,
        clocks.surface_noise_time, clocks.warp_outer_phase,
        clocks.warp_inner_phase,   clocks.palette_oscillation_phase};
    return snapshot;
  }

  /**
   * @brief Atomically restores a versioned ShaderBall configuration snapshot.
   * @return APPLIED on success; failures leave the effect and import notice
   * unchanged.
   */
  HS_COLD_MEMBER ConfigRestoreResult
  restore_full_config_snapshot(const FullConfigSnapshot &snapshot) {
    if (!config_version_supported(snapshot.schema_version))
      return ConfigRestoreResult::UNSUPPORTED_VERSION;

    Config next_accepted{};
    Config next_requested{};
    if (!decode_config_values(snapshot.accepted, next_accepted) ||
        !decode_config_values(snapshot.requested, next_requested))
      return ConfigRestoreResult::INVALID_VALUE;
    normalize_config_ranges(next_accepted);
    normalize_config_ranges(next_requested);
    RuntimeValues next_runtime = snapshot.runtime;
    if (!valid_snapshot_config(next_accepted) ||
        !valid_snapshot_config(next_requested))
      return ConfigRestoreResult::INVALID_VALUE;
    if (!admissible_config(next_accepted))
      return ConfigRestoreResult::INVALID_ACCEPTED;
    if (fixed_topology) {
      const InversePipelineId pipeline = preset_for_view(0).pipeline;
      if (resolve_pipeline_id(next_accepted) != pipeline)
        return ConfigRestoreResult::INVALID_ACCEPTED;
      if (resolve_pipeline_id(next_requested) != pipeline)
        return ConfigRestoreResult::INVALID_PENDING;
    }

    const ConfigValues migrated_accepted = encode_config_values(next_accepted);
    const ConfigValues migrated_requested =
        encode_config_values(next_requested);
    size_t next_pending_count = 0;
    for (size_t index = 0; index < CONFIG_FIELD_COUNT; ++index) {
      if (snapshot.pending[index] > 1)
        return ConfigRestoreResult::INVALID_PENDING;
      const bool differs =
          migrated_accepted[index] != migrated_requested[index];
      if ((snapshot.pending[index] != 0) != differs)
        return ConfigRestoreResult::INVALID_PENDING;
      next_pending_count += differs;
    }
    if (next_pending_count > pending_edits.size())
      return ConfigRestoreResult::INVALID_PENDING;

    if (snapshot.has_runtime)
      for (float value : snapshot.runtime)
        if (!std::isfinite(value))
          return ConfigRestoreResult::INVALID_VALUE;
    if (!prepare_resource_union(next_accepted, next_accepted))
      return ConfigRestoreResult::INVALID_ACCEPTED;

    state->param_morph.active = false;
    state->transition.active = false;
    accepted_config = next_accepted;
    requested_config = next_requested;
    published_config = next_accepted;
    active_slots = next_accepted.slots;
    active_pipeline = resolve_pipeline_id(next_accepted);
    blend.params = next_accepted.params;
    blend.palette_mapping =
        palette_mapping_weights(next_accepted.slots.palette_mapping);
    pending_edit_count = 0;
    for (size_t index = 0; index < CONFIG_FIELD_COUNT; ++index) {
      if (migrated_accepted[index] == migrated_requested[index])
        continue;
      const ConfigFieldId id = static_cast<ConfigFieldId>(index);
      const ConfigFieldLayout layout = config_field_layout(id);
      pending_edits[pending_edit_count++] = {nullptr, id, layout.offset,
                                             layout.size};
    }
    display_config = next_requested;
    if (snapshot.has_runtime) {
      ClockState &clocks = runtime.clocks;
      clocks.source_primary = next_runtime[0];
      clocks.source_secondary = next_runtime[1];
      clocks.source_angle = next_runtime[2];
      clocks.warp_outer_rotation = next_runtime[3];
      clocks.projection_spin = next_runtime[4];
      clocks.hue_noise_phase = next_runtime[5];
      clocks.source_noise_time = next_runtime[6];
      clocks.warp_inner_rotation = next_runtime[7];
      clocks.surface_noise_time = next_runtime[8];
      clocks.warp_outer_phase = next_runtime[9];
      clocks.warp_inner_phase = next_runtime[10];
      clocks.palette_oscillation_phase = next_runtime[11];
    }
    rebind_parameters();
    return ConfigRestoreResult::APPLIED;
  }
#endif

private:
  struct WalkDeltas {
    Quaternion projection;
    Quaternion outer;
  };

  struct ParamMorphRuntime {
    Params from;
    Params to;
    PaletteMappingWeights mapping_from;
    PaletteMappingWeights mapping_to;
    PaletteMapping mapping_destination = PaletteMapping::LINEAR;
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
    InversePipelineId from_pipeline = InversePipelineId::NONE;
    InversePipelineId to_pipeline = InversePipelineId::NONE;
  };

  struct StateBundle {
    FrameState frame;
    Config render_config;
    std::array<FastNoiseLite, MAX_NOISE_RESOURCES> noise_resources;
    std::array<NoiseFieldKey, MAX_NOISE_RESOURCES> prepared_noise_keys{};
    std::array<Pixel, PreparedHueRotation::LUT_SIZE> hue_rotation_lut;
    std::array<int8_t, PreparedHueNoise::LUT_SIZE> hue_noise_lut;
    /** @brief Inputs the resident hue-noise table was built from; scale 0 marks
     *         it unbuilt. */
    float hue_noise_lut_scale = 0.0f;
    float hue_noise_lut_phase = 0.0f;
    FastNoiseLite projection_walk_noise;
    FastNoiseLite outer_walk_noise;
    ParamMorphRuntime param_morph;
    TransitionRuntime transition;

    HS_COLD_MEMBER StateBundle() = default;
  };

  struct ThroughClearPhase {
    float alpha;
    bool from_endpoint;
    bool clear;
  };

  HS_COLD_MEMBER static constexpr bool warp_uses_noise(WarpStageKind kind) {
    return kind == WarpStageKind::VECTOR_NOISE ||
           kind == WarpStageKind::CURL_FLOW;
  }

  HS_COLD_MEMBER static constexpr bool seam_sensitive_warp(WarpStageKind kind) {
    return kind == WarpStageKind::VECTOR_NOISE ||
           kind == WarpStageKind::CURL_FLOW;
  }

  HS_COLD_MEMBER static constexpr NoiseFieldKey
  warp_resource_key(const WarpStageSpec &spec) {
    return {NoiseDomain::PROJECTED_2D,
            spec.basis,
            spec.seed,
            spec.kind == WarpStageKind::CURL_FLOW
                ? NoiseChannelLayout::CURL_V1
                : (spec.basis == NoiseBasis::SIMPLEX
                       ? NoiseChannelLayout::DIRECT_VECTOR_V2
                       : NoiseChannelLayout::DIRECT_V1),
            1,
            1,
            static_cast<uint8_t>(spec.kind == WarpStageKind::CURL_FLOW ? 1 : 0),
            FastNoiseLite::NoiseType_OpenSimplex2,
            1.0f};
  }

  HS_COLD_MEMBER static constexpr NoiseFieldKey
  source_resource_key(const Config &config) {
    return {config.slots.function == Function::NOISE_CONTOUR_SPHERE
                ? NoiseDomain::SPHERE_3D
                : NoiseDomain::PROJECTED_2D,
            config.params.source.noise_basis,
            config.params.source.noise_seed,
            NoiseChannelLayout::SCALAR_V1,
            1,
            1,
            0,
            FastNoiseLite::NoiseType_OpenSimplex2,
            1.0f};
  }

  HS_COLD_MEMBER static constexpr NoiseFieldKey
  surface_noise_resource_key(const Config &config) {
    return {NoiseDomain::SPHERE_3D,
            config.params.surface_noise.basis,
            config.params.surface_noise.seed,
            config.slots.surface_noise == SurfaceNoise::CURL
                ? (config.params.surface_noise.basis == NoiseBasis::SIMPLEX
                       ? NoiseChannelLayout::CURL_ANALYTIC_V2
                       : NoiseChannelLayout::CURL_V1)
                : (config.params.surface_noise.basis == NoiseBasis::SIMPLEX
                       ? NoiseChannelLayout::DIRECT_VECTOR_V2
                       : NoiseChannelLayout::DIRECT_V1),
            1,
            1,
            static_cast<uint8_t>(
                config.slots.surface_noise == SurfaceNoise::CURL ? 1 : 0),
            FastNoiseLite::NoiseType_OpenSimplex2,
            1.0f};
  }

  HS_COLD_MEMBER static constexpr NoiseFieldKey color_noise_resource_key() {
    return {NoiseDomain::SPHERE_3D,
            NoiseBasis::SIMPLEX,
            6047,
            NoiseChannelLayout::SCALAR_V1,
            1,
            1,
            0,
            FastNoiseLite::NoiseType_OpenSimplex2,
            1.0f};
  }

  HS_COLD_MEMBER static constexpr bool
  append_resource_key(const NoiseFieldKey &key,
                      std::array<NoiseFieldKey, MAX_NOISE_RESOURCES> &keys,
                      size_t &count) {
    for (size_t index = 0; index < count; ++index)
      if (keys[index] == key)
        return true;
    if (count == keys.size())
      return false;
    keys[count++] = key;
    return true;
  }

  HS_COLD_MEMBER static constexpr bool append_config_resource_keys(
      const Config &config,
      std::array<NoiseFieldKey, MAX_NOISE_RESOURCES> &keys, size_t &count) {
    if (warp_uses_noise(config.slots.warp_program.outer.kind) &&
        !append_resource_key(warp_resource_key(config.slots.warp_program.outer),
                             keys, count))
      return false;
    if (warp_uses_noise(config.slots.warp_program.inner.kind) &&
        !append_resource_key(warp_resource_key(config.slots.warp_program.inner),
                             keys, count))
      return false;
    if (is_noise_contour(config.slots.function) &&
        !append_resource_key(source_resource_key(config), keys, count))
      return false;
    if (config.slots.surface_noise != SurfaceNoise::NONE &&
        !append_resource_key(surface_noise_resource_key(config), keys, count))
      return false;
    if (config.slots.hue_shift == HueShiftMode::NOISE &&
        config.params.color.hue_shift_amount != 0.0f &&
        !append_resource_key(color_noise_resource_key(), keys, count))
      return false;
    return true;
  }

  HS_COLD_MEMBER static constexpr bool resource_union_fits(const Config &from,
                                                           const Config &to) {
    std::array<NoiseFieldKey, MAX_NOISE_RESOURCES> keys{};
    size_t count = 0;
    return append_config_resource_keys(from, keys, count) &&
           append_config_resource_keys(to, keys, count);
  }

  HS_COLD_MEMBER bool prepare_resource_union(const Config &from,
                                             const Config &to) {
    std::array<NoiseFieldKey, MAX_NOISE_RESOURCES> keys{};
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
  resolve_resource(const NoiseFieldKey &key) const {
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
    return is_noise_contour(config.slots.function)
               ? resolve_resource(source_resource_key(config))
               : nullptr;
  }

  HS_COLD_MEMBER const FastNoiseLite *
  resolve_surface_noise_resource(const Config &config) const {
    return config.slots.surface_noise != SurfaceNoise::NONE
               ? resolve_resource(surface_noise_resource_key(config))
               : nullptr;
  }

  HS_COLD_MEMBER const FastNoiseLite *
  resolve_color_noise_resource(const Config &config) const {
    if (config.slots.hue_shift != HueShiftMode::NOISE ||
        config.params.color.hue_shift_amount == 0.0f)
      return nullptr;
    return resolve_resource(color_noise_resource_key());
  }

  HS_COLD_MEMBER const BakedPalette &palette_for(PaletteMode mode) const {
    switch (mode) {
    case PaletteMode::TRIADIC:
      return triadic_palette_cycler.palette();
    case PaletteMode::COMPLEMENTARY:
      return complementary_palette_cycler.palette();
    case PaletteMode::ANALOGOUS:
      return analogous_palette_cycler.palette();
    }
    __builtin_unreachable();
  }

  PaletteMode visible_palette_mode() const {
    if (!state->transition.active)
      return active_slots.palette;
    const ThroughClearPhase phase = through_clear_phase(
        state->transition.elapsed, state->transition.duration);
    return phase.from_endpoint ? state->transition.from_config.slots.palette
                               : state->transition.to_config.slots.palette;
  }

  float visible_palette_chroma() const {
    if (!state->transition.active)
      return blend.params.color.palette_chroma;
    const ThroughClearPhase phase = through_clear_phase(
        state->transition.elapsed, state->transition.duration);
    return phase.from_endpoint
               ? state->transition.from_config.params.color.palette_chroma
               : state->transition.to_config.params.color.palette_chroma;
  }

  void step_generated_palettes(PaletteMode visible) {
    if (visible == PaletteMode::TRIADIC)
      triadic_palette_cycler.step();
    else
      triadic_palette_cycler.advance_without_display();
    if (visible == PaletteMode::COMPLEMENTARY)
      complementary_palette_cycler.step();
    else
      complementary_palette_cycler.advance_without_display();
    if (visible == PaletteMode::ANALOGOUS)
      analogous_palette_cycler.step();
    else
      analogous_palette_cycler.advance_without_display();
  }

  HS_COLD_MEMBER void update_palette_chroma(float chroma) {
    if (chroma == palette_chroma)
      return;
    palette_chroma = chroma;
    triadic_palette_cycler.set_generated_chroma(chroma);
    complementary_palette_cycler.set_generated_chroma(chroma);
    analogous_palette_cycler.set_generated_chroma(chroma);
  }

  HS_FLASH_MEMBER static PreparedWarpStage
  prepare_warp_stage(const WarpStageSpec &spec, const WarpStageParams &params,
                     float stage_phase,
                     const Complex &source_period = Complex(),
                     float affine_rotation = 0.0f) {
    PreparedWarpStage prepared{};
    float rotation = params.rotation;
    if (spec.kind == WarpStageKind::VECTOR_NOISE)
      rotation = params.vector_angle;
    else if (spec.kind == WarpStageKind::WAVE_SHEAR)
      rotation = params.field_angle;
    if (spec.kind == WarpStageKind::AFFINE_FRAME) {
      const float phase = TWO_PI_F * wrap_t(stage_phase);
      const float phase_cos = cosf(phase);
      rotation = affine_rotation;
      prepared.transform.affine = {
          wrap_t(stage_phase) * params.translation_x * source_period.re,
          wrap_t(stage_phase) * params.translation_y * source_period.im,
          powf(params.scale_x, phase_cos),
          powf(params.scale_y, phase_cos),
          params.shear * phase_cos,
      };
    } else if (spec.kind == WarpStageKind::MIRROR_TILE) {
      prepared.transform.mirror = {
          wrap_t(params.offset_x / params.cell_x + stage_phase) * params.cell_x,
          wrap_t(params.offset_y / params.cell_y) * params.cell_y,
      };
    } else if (spec.kind == WarpStageKind::VORTEX) {
      const float orbit_phase = TWO_PI_F * stage_phase;
      prepared.transform.vortex = {
          params.center_x + params.center_orbit_radius * cosf(orbit_phase),
          params.center_y + params.center_orbit_radius * sinf(orbit_phase),
          params.radius * params.radius,
          TWO_PI_F * params.turns,
      };
    } else if (spec.kind == WarpStageKind::VECTOR_NOISE) {
      const float angle = TWO_PI_F * wrap_t(stage_phase);
      prepared.transform.noise_loop = {
          NOISE_LOOP_RADIUS * sinf(angle) * 0.7071067811865475f,
          NOISE_LOOP_RADIUS * cosf(angle),
      };
    }
    prepared.rotation_cos = cosf(rotation);
    prepared.rotation_sin = sinf(rotation);
    return prepared;
  }

  HS_COLD_MEMBER FrameState prepare_frame() const {
    FrameState frame;
    prepare_frame({active_slots, blend.params}, runtime, frame);
    return frame;
  }

  HS_COLD_MEMBER FrameState prepare_frame(const Config &config,
                                          const LookRuntime &look) const {
    FrameState frame;
    prepare_frame(config, look, frame);
    return frame;
  }

  HS_COLD_MEMBER void prepare_frame(const Config &config,
                                    const LookRuntime &look,
                                    FrameState &frame) const {
    const bool animated_projection =
        config.slots.projection_frame == ProjectionFramePolicy::SPIN_WANDER;
    const Complex source_period = source_cartesian_period(config);
    const PreparedWarpProgram prepared_warp{
        prepare_warp_stage(config.slots.warp_program.outer,
                           config.params.warp.outer,
                           look.clocks.warp_outer_phase, source_period,
                           look.clocks.warp_outer_rotation),
        prepare_warp_stage(config.slots.warp_program.inner,
                           config.params.warp.inner,
                           look.clocks.warp_inner_phase, source_period,
                           look.clocks.warp_inner_rotation)};
    const float surface_phase =
        TWO_PI_F * wrap_t(look.clocks.surface_noise_time);
    const float surface_direction =
        TWO_PI_F * config.params.surface_noise.direction;
    const PreparedSurfaceNoise prepared_surface_noise{
        Vector(NOISE_LOOP_RADIUS * cosf(surface_phase),
               NOISE_LOOP_RADIUS * sinf(surface_phase), 0.0f),
        cosf(surface_direction), sinf(surface_direction)};
    const BakedPalette *palette = &palette_for(config.slots.palette);
    PreparedHueRotation prepared_hue_rotation{
        state->hue_rotation_lut.data(),
        config.slots.hue_shift != HueShiftMode::NONE &&
            config.params.color.hue_shift_amount != 0.0f};
    if (prepared_hue_rotation.active)
      prepare_hue_rotation_lut(prepared_hue_rotation, *palette);
    const FastNoiseLite *color_noise = resolve_color_noise_resource(config);
    PreparedHueNoise prepared_hue_noise{
        state->hue_noise_lut.data(),
        config.slots.hue_shift == HueShiftMode::NOISE &&
            config.params.color.hue_shift_amount != 0.0f};
    if (prepared_hue_noise.active &&
        (state->hue_noise_lut_scale != config.params.color.hue_noise_scale ||
         state->hue_noise_lut_phase != look.clocks.hue_noise_phase)) {
      prepare_hue_noise_lut(prepared_hue_noise, *color_noise,
                            config.params.color.hue_noise_scale,
                            look.clocks.hue_noise_phase);
      state->hue_noise_lut_scale = config.params.color.hue_noise_scale;
      state->hue_noise_lut_phase = look.clocks.hue_noise_phase;
    }
    frame.slots = config.slots;
    frame.params = config.params;
    frame.palette_mapping =
        state->param_morph.active
            ? blend.palette_mapping
            : palette_mapping_weights(config.slots.palette_mapping);
    frame.clocks = look.clocks;
    frame.prepared_source = {
        look.clocks.source_primary, look.clocks.source_secondary,
        look.clocks.source_angle, fast_cosf(look.clocks.source_angle),
        fast_sinf(look.clocks.source_angle)};
    frame.transforms = {animated_projection ? look.transforms.projection_conj
                                            : Quaternion(),
                        look.transforms.outer_conj};
    frame.prepared_warp = prepared_warp;
    frame.prepared_surface_noise = prepared_surface_noise;
    frame.prepared_hue_rotation = prepared_hue_rotation;
    frame.prepared_hue_noise = prepared_hue_noise;
    frame.resources = {resolve_warp_resource(config.slots.warp_program.outer),
                       resolve_warp_resource(config.slots.warp_program.inner),
                       resolve_source_resource(config),
                       resolve_surface_noise_resource(config),
                       color_noise,
                       palette};
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

  HS_FLASH_MEMBER void draw_through_clear_transition(Canvas &canvas) {
    const ThroughClearPhase phase = through_clear_phase(
        state->transition.elapsed, state->transition.duration);
    if (phase.clear)
      return;
    PreparedEndpoint prepared;
    const Config &config = phase.from_endpoint ? state->transition.from_config
                                               : state->transition.to_config;
    const LookRuntime &look = phase.from_endpoint
                                  ? state->transition.from_runtime
                                  : state->transition.to_runtime;
    const InversePipelineId pipeline = phase.from_endpoint
                                           ? state->transition.from_pipeline
                                           : state->transition.to_pipeline;
    HS_CHECK(prepare_endpoint(config, look, phase.alpha, pipeline, prepared),
             "ShaderBall transition endpoint has no renderer");
    draw_endpoint(canvas, prepared,
                  phase.from_endpoint ? ProfileEndpoint::FROM
                                      : ProfileEndpoint::TO);
  }

  HS_COLD_MEMBER bool prepare_endpoint(const Config &config,
                                       const LookRuntime &look, float alpha,
                                       InversePipelineId selected,
                                       PreparedEndpoint &prepared) const {
    const ProgramDescriptor *program = get_inverse_program(selected);
    ShadeFunction shade;
    bool (*resources_ready)(const FrameState &);
    if (program != nullptr) {
      if (program->key != make_topology_key(config) ||
          !program->continuous_parameters_supported(config))
        return false;
      shade = program->shade;
      resources_ready = program->resources_ready;
    } else {
#if HS_ENABLE_SHADERBALL_DYNAMIC_BACKEND
      if (selected != InversePipelineId::NONE || !valid_config(config))
        return false;
      shade = &shade_dynamic;
      resources_ready = &pipeline_resources_ready;
#else
      return false;
#endif
    }
    prepared.frame = &state->frame;
    prepare_frame(config, look, *prepared.frame);
    if (!resources_ready(*prepared.frame))
      return false;
    prepared.shade = shade;
    prepared.pipeline = selected;
    prepared.alpha = alpha;
#if defined(HS_PROFILE_ENABLE)
    prepared.preset = selected_preset_index(config, selected);
#endif
    return true;
  }

  HS_FLASH_MEMBER void
  draw_endpoint(Canvas &canvas, PreparedEndpoint &prepared,
                ProfileEndpoint endpoint = ProfileEndpoint::STEADY) {
#if defined(HS_PROFILE_ENABLE)
    emit_pullback_program(prepared, endpoint);
#else
    (void)endpoint;
#endif
    FrameShader shader{prepared.frame, prepared.alpha, prepared.shade};
    HS_PROFILE(sb_shader_draw);
    Scan::Shader::draw<W, H, 1>(canvas, shader);
  }

  static constexpr const char *pipeline_name(InversePipelineId pipeline) {
    switch (pipeline) {
    case InversePipelineId::GLITCH_NOISE_GRID_WAVE_SHEAR:
      return "GLITCH_NOISE_GRID_WAVE_SHEAR";
    case InversePipelineId::KALEIDOSCOPE_TWIN_WAVE_INNER_MIRROR:
      return "KALEIDOSCOPE_TWIN_WAVE_INNER_MIRROR";
    case InversePipelineId::GNOMONIC_KALEIDOSCOPE_GRID_MIRROR:
      return "GNOMONIC_KALEIDOSCOPE_GRID_MIRROR";
    case InversePipelineId::GNOMONIC_GLITCH_GRID_MIRROR:
      return "GNOMONIC_GLITCH_GRID_MIRROR";
    case InversePipelineId::PEIRCE_DODECAHEDRAL_GRID:
      return "PEIRCE_DODECAHEDRAL_GRID";
    case InversePipelineId::GNOMONIC_DODECAHEDRAL_GRID_WAVE_MIRROR:
      return "GNOMONIC_DODECAHEDRAL_GRID_WAVE_MIRROR";
    case InversePipelineId::GNOMONIC_AFFINE_LATTICE_CONTOUR:
      return "GNOMONIC_AFFINE_LATTICE_CONTOUR";
    case InversePipelineId::SINUSOIDAL_CURL_LATTICE:
      return "SINUSOIDAL_CURL_LATTICE";
    case InversePipelineId::STEREOGRAPHIC_PRISM_POLAR_WAVE_LATTICE:
      return "STEREOGRAPHIC_PRISM_POLAR_WAVE_LATTICE";
    case InversePipelineId::GNOMONIC_DODECAHEDRAL_GRID_VECTOR_MIRROR:
      return "GNOMONIC_DODECAHEDRAL_GRID_VECTOR_MIRROR";
    case InversePipelineId::STEREOGRAPHIC_DODECAHEDRAL_GRID_INNER_MIRROR:
      return "STEREOGRAPHIC_DODECAHEDRAL_GRID_INNER_MIRROR";
    case InversePipelineId::
        STEREOGRAPHIC_HEXAGONAL_PRISM_TWIN_WAVE_INNER_MIRROR:
      return "STEREOGRAPHIC_HEXAGONAL_PRISM_TWIN_WAVE_INNER_MIRROR";
    case InversePipelineId::EQUIRECTANGULAR_DODECAHEDRAL_GRID_INNER_MIRROR:
      return "EQUIRECTANGULAR_DODECAHEDRAL_GRID_INNER_MIRROR";
    case InversePipelineId::STEREOGRAPHIC_GLITCH_GRID_MIRROR:
      return "STEREOGRAPHIC_GLITCH_GRID_MIRROR";
    case InversePipelineId::STEREOGRAPHIC_MOBIUS_TWIN_WAVE_INNER_MIRROR:
      return "STEREOGRAPHIC_MOBIUS_TWIN_WAVE_INNER_MIRROR";
    case InversePipelineId::COUNT:
      return "COUNT";
    case InversePipelineId::NONE:
      return "NONE";
    }
    return "NONE";
  }

#if defined(HS_PROFILE_ENABLE)
  static constexpr const char *profile_endpoint_name(ProfileEndpoint endpoint) {
    switch (endpoint) {
    case ProfileEndpoint::STEADY:
      return "steady";
    case ProfileEndpoint::FROM:
      return "from";
    case ProfileEndpoint::TO:
      return "to";
    }
    return "steady";
  }

  size_t selected_preset_index(const Config &config,
                               InversePipelineId pipeline) const {
    for (size_t index = 0; index < preset_count_for_view(); ++index)
      if (preset_for_view(index).pipeline == pipeline &&
          preset_for_view(index).config == config)
        return index;
    return getPresetIndex();
  }

  void emit_pullback_program(const PreparedEndpoint &prepared,
                             ProfileEndpoint endpoint) {
    if (profile_program_valid && profile_program_preset == prepared.preset &&
        profile_program_pipeline == prepared.pipeline &&
        profile_program_endpoint == endpoint)
      return;
    hs::log("Pullback program: preset=%u/%u pipeline=%s endpoint=%s",
            static_cast<unsigned>(prepared.preset),
            static_cast<unsigned>(preset_count_for_view()),
            pipeline_name(prepared.pipeline), profile_endpoint_name(endpoint));
    profile_program_valid = true;
    profile_program_preset = prepared.preset;
    profile_program_pipeline = prepared.pipeline;
    profile_program_endpoint = endpoint;
  }
#endif

#if HS_ENABLE_SHADERBALL_DYNAMIC_BACKEND ||                                    \
    (HS_ENABLE_TEST_HOOKS && HS_ENABLE_TEST_ORACLES)
  /**
   * @brief Shades one sphere sample by pulling it back to a source coordinate.
   * @param view Unit direction of the visible sphere point.
   * @param frame Immutable snapshot of slots, parameters, clocks, transforms,
   *        and palette resources for this frame.
   * @return Straight-alpha colour for the sample.
   * @details Walks outer camera, surface lens, and projection backward. A
   * strict projection whose two lens branches land in different regions cannot
   * be joined in the plane, so the branches are shaded separately and their
   * outputs blended instead.
   */
  static Color4 shade_dynamic(const Vector &view, const FrameState &frame) {
    const Vector outer_local = outer_camera_lookup(view, frame);
    const ProjectedLookup projected =
        surface_lens_project_lookup(outer_local, frame);
    return shade_projected(projected, frame);
  }

  /**
   * @brief Runs the planar half of the pullback: warps, samples, shapes, and
   *        colorizes.
   * @param projected Plane coordinates and seam metadata for the sample.
   * @param frame Frame snapshot.
   * @return Straight-alpha colour for the sample.
   */
  HS_FLASH_MEMBER static Color4
  shade_projected(const ProjectedLookup &projected, const FrameState &frame) {
    HS_SB_STAGE_MARK(stage_start);
    const PlanarWarpResult warped = planar_warp_lookup(projected, frame);
    HS_SB_STAGE_SPAN(planar_warp, stage_start);
    const Complex source_coords = condition_source_coords(warped.coords, frame);
    const float field = sample_source(source_coords, projected, frame);
    HS_SB_STAGE_SPAN(source, stage_start);
    const MaterialSample material =
        shape_material(field, projected, warped, frame);
    HS_SB_STAGE_SPAN(material, stage_start);
    const Color4 color = colorize(material, frame);
    HS_SB_STAGE_SPAN(color, stage_start);
    return color;
  }
#endif

  __attribute__((always_inline)) static ProjectedLookup
  stereographic_lookup(const Vector &local, const FrameState &frame) {
    const Complex coords = stereo(local);
    const float r_sq = coords.re * coords.re + coords.im * coords.im;
    return {coords,
            0,
            0,
            BOUNDARY_SINGULAR,
            std::max(0.0f, 1.0f - local.y),
            pole_attenuation(r_sq, frame.params.projection.pole_fade),
            0,
            0,
            0,
            1.0f,
            local};
  }

  __attribute__((always_inline)) static ProjectedLookup
  selected_stereographic_lookup(const Vector &v, const FrameState &frame) {
    return stereographic_lookup(rotate(v, frame.transforms.projection_conj),
                                frame);
  }

  __attribute__((always_inline)) static ProjectedLookup
  profiled_stereographic_lookup(const Vector &v, const FrameState &frame) {
    HS_SB_STAGE_MARK(stage_start);
    const ProjectedLookup projected = selected_stereographic_lookup(v, frame);
    HS_SB_STAGE_SPAN(projection, stage_start);
    return projected;
  }

  template <SurfaceLens LENS>
  __attribute__((always_inline)) static Vector
  selected_lens_lookup(const Vector &v) {
    static_assert(LENS == SurfaceLens::NONE || LENS == SurfaceLens::GLITCH ||
                  LENS == SurfaceLens::KALEIDOSCOPE ||
                  LENS == SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL ||
                  LENS == SurfaceLens::KALEIDOSCOPE_TRIANGULAR_PRISM);
    HS_SB_STAGE_MARK(stage_start);
    const Vector lensed = [&]() {
      if constexpr (LENS == SurfaceLens::NONE)
        return v;
      else if constexpr (LENS == SurfaceLens::GLITCH)
        return lenses::glitch_lens(v);
      else if constexpr (LENS == SurfaceLens::KALEIDOSCOPE)
        return lenses::kaleidoscope_lens(v);
      else if constexpr (LENS == SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL)
        return lenses::dodecahedral_kaleidoscope_lens(v);
      else
        return lenses::polyhedral_kaleidoscope_lens(
            v, lenses::TRIANGULAR_PRISM_MIRRORS);
    }();
    HS_SB_STAGE_SPAN(lens, stage_start);
    return lensed;
  }

  /** @brief Whether this frame's colorizer reads displacement metadata. */
  __attribute__((always_inline)) static bool
  tracks_displacement(const FrameState &frame) {
    return frame.prepared_hue_rotation.active &&
           frame.slots.hue_shift == HueShiftMode::WARP_DISPLACEMENT;
  }

  __attribute__((always_inline)) static float
  surface_noise_path_length(const Vector &step, const FrameState &frame) {
    if (!tracks_displacement(frame))
      return 0.0f;
    return sqrtf(dot(step, step));
  }

  /**
   * @brief Builds one manifest row.
   * @tparam Pipeline Compiled inverse pipeline the row selects.
   * @tparam Id Stable identifier for the row.
   * @tparam Key Topology the row is matched on.
   * @param continuous Predicate over the continuous parameters @p Pipeline
   *        serves.
   * @details Rejects at compile time a pipeline whose stages hardcode a
   * topology facet that @p Key does not carry.
   */
  template <typename Pipeline, InversePipelineId Id, TopologyKey Key>
  static constexpr ProgramDescriptor
  make_program(bool (*continuous)(const Config &)) {
    static_assert(Pipeline::implements(Key),
                  "inverse pipeline does not implement its topology key");
    return {Id, Key, &Pipeline::shade, continuous, &pipeline_resources_ready};
  }

  HS_COLD_MEMBER static const std::array<ProgramDescriptor, 15> &
  inverse_programs() {
    static constexpr std::array<ProgramDescriptor, 15> PROGRAMS{{
        make_program<GlitchNoiseGridWaveShearPipeline,
                     InversePipelineId::GLITCH_NOISE_GRID_WAVE_SHEAR,
                     make_topology_key(wave_shear_generated_preset())>(
            &all_continuous_parameters_supported),
        make_program<KaleidoscopeTwinWaveInnerMirrorPipeline,
                     InversePipelineId::KALEIDOSCOPE_TWIN_WAVE_INNER_MIRROR,
                     make_topology_key(kaleidoscope_mirror_preset())>(
            &all_continuous_parameters_supported),
        make_program<GnomonicKaleidoscopeGridMirrorPipeline,
                     InversePipelineId::GNOMONIC_KALEIDOSCOPE_GRID_MIRROR,
                     make_topology_key(
                         gnomonic_kaleidoscope_grid_mirror_preset())>(
            &all_continuous_parameters_supported),
        make_program<GnomonicGlitchGridMirrorPipeline,
                     InversePipelineId::GNOMONIC_GLITCH_GRID_MIRROR,
                     make_topology_key(
                         gnomonic_grid_mirror_preset(SurfaceLens::GLITCH))>(
            &all_continuous_parameters_supported),
        make_program<PeirceDodecahedralGridPipeline,
                     InversePipelineId::PEIRCE_DODECAHEDRAL_GRID,
                     make_topology_key(peirce_dodecahedral_generated_preset())>(
            &all_continuous_parameters_supported),
        make_program<GnomonicDodecahedralGridWaveMirrorPipeline,
                     InversePipelineId::GNOMONIC_DODECAHEDRAL_GRID_WAVE_MIRROR,
                     make_topology_key(gnomonic_wave_shear_grid_preset())>(
            &all_continuous_parameters_supported),
        make_program<GnomonicAffineLatticeContourPipeline,
                     InversePipelineId::GNOMONIC_AFFINE_LATTICE_CONTOUR,
                     make_topology_key(
                         gnomonic_affine_lattice_contour_preset())>(
            &all_continuous_parameters_supported),
        make_program<SinusoidalCurlLatticePipeline,
                     InversePipelineId::SINUSOIDAL_CURL_LATTICE,
                     make_topology_key(sinusoidal_lattice_curl_preset(1.0f))>(
            &all_continuous_parameters_supported),
        make_program<StereographicPrismPolarWaveLatticePipeline,
                     InversePipelineId::STEREOGRAPHIC_PRISM_POLAR_WAVE_LATTICE,
                     make_topology_key(
                         stereographic_prism_polar_wave_lattice_preset())>(
            &all_continuous_parameters_supported),
        make_program<
            GnomonicDodecahedralGridVectorMirrorPipeline,
            InversePipelineId::GNOMONIC_DODECAHEDRAL_GRID_VECTOR_MIRROR,
            make_topology_key(
                gnomonic_dodecahedral_vector_mirror_grid_preset())>(
            &all_continuous_parameters_supported),
        make_program<
            StereographicDodecahedralGridInnerMirrorPipeline,
            InversePipelineId::STEREOGRAPHIC_DODECAHEDRAL_GRID_INNER_MIRROR,
            make_topology_key(
                stereographic_dodecahedral_grid_inner_mirror_preset())>(
            &all_continuous_parameters_supported),
        make_program<
            StereographicHexagonalPrismTwinWaveInnerMirrorPipeline,
            InversePipelineId::
                STEREOGRAPHIC_HEXAGONAL_PRISM_TWIN_WAVE_INNER_MIRROR,
            make_topology_key(
                stereographic_hexagonal_prism_twin_wave_mirror_preset())>(
            &all_continuous_parameters_supported),
        make_program<
            EquirectangularDodecahedralGridInnerMirrorPipeline,
            InversePipelineId::EQUIRECTANGULAR_DODECAHEDRAL_GRID_INNER_MIRROR,
            make_topology_key(
                equirectangular_dodecahedral_double_mapping_grid_inner_mirror_preset())>(
            &all_continuous_parameters_supported),
        make_program<StereographicGlitchGridMirrorPipeline,
                     InversePipelineId::STEREOGRAPHIC_GLITCH_GRID_MIRROR,
                     make_topology_key(
                         stereographic_glitch_grid_mirror_preset())>(
            &all_continuous_parameters_supported),
        make_program<
            StereographicMobiusTwinWaveInnerMirrorPipeline,
            InversePipelineId::STEREOGRAPHIC_MOBIUS_TWIN_WAVE_INNER_MIRROR,
            make_topology_key(
                stereographic_mobius_twin_wave_inner_mirror_preset())>(
            &all_continuous_parameters_supported),
    }};
    return PROGRAMS;
  }

  HS_COLD_MEMBER static const ProgramDescriptor *
  find_inverse_program(const Config &config) {
    const TopologyKey key = make_topology_key(config);
    for (const ProgramDescriptor &program : inverse_programs())
      if (program.key == key && program.continuous_parameters_supported(config))
        return &program;
    return nullptr;
  }

  HS_COLD_MEMBER static const ProgramDescriptor *
  get_inverse_program(InversePipelineId id) {
    for (const ProgramDescriptor &program : inverse_programs())
      if (program.id == id)
        return &program;
    return nullptr;
  }

  HS_COLD_MEMBER static InversePipelineId
  resolve_pipeline_id(const Config &config) {
    const ProgramDescriptor *program = find_inverse_program(config);
    return program == nullptr ? InversePipelineId::NONE : program->id;
  }

  HS_COLD_MEMBER static const ProgramDescriptor *
  resolve_inverse_program(const FrameState &frame) {
    const Config config{frame.slots, frame.params};
    const ProgramDescriptor *program = find_inverse_program(config);
    if (program == nullptr || !program->resources_ready(frame))
      return nullptr;
    return program;
  }

#if HS_ENABLE_SHADERBALL_DYNAMIC_BACKEND ||                                    \
    (HS_ENABLE_TEST_HOOKS && HS_ENABLE_TEST_ORACLES)
  HS_COLD_MEMBER static typename FrameShader::ShadeFunction
  resolve_shade_function(const FrameState &frame) {
    const ProgramDescriptor *program = resolve_inverse_program(frame);
    HS_CHECK(program != nullptr,
             "ShaderBall test topology has no compiled inverse pipeline");
    return program->shade;
  }
#endif

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

  __attribute__((always_inline)) static Vector
  outer_camera_lookup(const Vector &view, const FrameState &frame) {
    return rotate(view, frame.transforms.outer_conj);
  }

#if HS_ENABLE_SHADERBALL_DYNAMIC_BACKEND ||                                    \
    (HS_ENABLE_TEST_HOOKS && HS_ENABLE_TEST_ORACLES)
  static ProjectedLookup surface_lens_project_lookup(const Vector &v,
                                                     const FrameState &frame) {
    const Slots &slots = frame.slots;
    Vector pre_lens = v;
    float surface_path_length = 0.0f;
    if (slots.surface_noise != SurfaceNoise::NONE &&
        slots.surface_noise_placement == SurfaceNoisePlacement::BEFORE_LENS) {
      HS_SB_STAGE_MARK(surface_start);
      const SurfaceNoiseResult displaced = apply_surface_noise_result(v, frame);
      pre_lens = displaced.sphere;
      surface_path_length = displaced.path_length;
      HS_SB_STAGE_SPAN(surface_noise, surface_start);
    }
    const Vector lensed = slots.surface_lens == SurfaceLens::NONE
                              ? pre_lens
                              : profiled_apply_lens(pre_lens, frame);
    Vector post_lens = lensed;
    if (slots.surface_noise != SurfaceNoise::NONE &&
        slots.surface_noise_placement == SurfaceNoisePlacement::AFTER_LENS) {
      HS_SB_STAGE_MARK(surface_start);
      const SurfaceNoiseResult displaced =
          apply_surface_noise_result(lensed, frame);
      post_lens = displaced.sphere;
      surface_path_length = displaced.path_length;
      HS_SB_STAGE_SPAN(surface_noise, surface_start);
    }
    ProjectedLookup projected = profiled_project_branch(post_lens, frame);
    projected.surface_path_length = surface_path_length;
    return projected;
  }
#endif

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
    return Pullback::Projection::from_kernel(result, coordinate_scale);
  }

  HS_FLASH_MEMBER static ProjectedLookup
  project_bonne(const Vector &local, const FrameState &frame) {
    return Pullback::Projection::bonne(
        local, frame.params.projection.central_meridian,
        (frame.slots.bonne_hemisphere == BonneHemisphere::NORTH ? 1.0f
                                                                : -1.0f) *
            frame.params.projection.bonne_standard_parallel,
        frame.params.projection.coordinate_scale);
  }

  HS_FLASH_MEMBER static ProjectedLookup
  project_peirce(const Vector &local, const FrameState &frame) {
    if (frame.slots.peirce_layout == PeirceLayout::SQUARE &&
        frame.params.projection.central_meridian == 0.0f &&
        projection_edge_distance_required(frame))
      return Pullback::Projection::peirce_fast_square(
          local, frame.params.projection.coordinate_scale);
    return Pullback::Projection::peirce(
        local, frame.params.projection.central_meridian,
        static_cast<uint8_t>(frame.slots.peirce_layout),
        frame.params.projection.layout_scroll,
        projection_edge_distance_required(frame),
        frame.params.projection.coordinate_scale);
  }

  HS_FLASH_MEMBER static ProjectedLookup
  project_airocean(const Vector &local, const FrameState &frame) {
    return Pullback::Projection::airocean(
        local, frame.params.projection.central_meridian,
        frame.slots.airocean_layout == AiroceanLayout::HORIZONTAL,
        projection_edge_distance_required(frame),
        frame.params.projection.coordinate_scale);
  }

  HS_FLASH_MEMBER static ProjectedLookup
  project_sinusoidal(const Vector &local, const FrameState &frame) {
    return finalize_projection(
        local,
        projections::folded_sinusoidal(
            local, frame.params.projection.central_meridian),
        Projection::SINUSOIDAL, frame.params.projection.pole_fade);
  }

  HS_FLASH_MEMBER static ProjectedLookup
  project_equirectangular(const Vector &local, const FrameState &frame) {
    return Pullback::Projection::equirectangular(
        local, frame.params.projection.central_meridian,
        frame.params.projection.pole_fade);
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
    ProjectedLookup projected = [&]() {
      if (frame.slots.projection != Projection::STEREOGRAPHIC)
        return project_nonstereographic(local, frame);
      return stereographic_lookup(local, frame);
    }();
    projected.sphere = local;
    return projected;
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
    case Projection::STEREOGRAPHIC:
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
    return {
        coords,
        selected->region_id,
        selected->component_id,
        selected->boundary_flags,
        selected->fade_edge_distance,
        pole_attenuation(r_sq, pole_fade),
        selected->flags,
        selected->traits,
        selected->edge_class,
        hs::lerp(direct.domain_coverage, lensed.domain_coverage, mix),
        nlerp_unit(direct.sphere, lensed.sphere, mix),
        hs::lerp(direct.surface_path_length, lensed.surface_path_length, mix)};
  }

  /**
   * @brief Pulls plane coordinates back through both warp stages.
   * @param projected Projection output; supplies the stage input coordinates
   *        and the weight and edge distance the stage envelopes read.
   * @param frame Frame snapshot.
   * @return Source-side coordinates plus the path length accumulated by the
   *         warp program.
   * @details Pullback order is Planar Warp 1 then Planar Warp 2, the reverse
   * of the authored order `source -> Warp 2 -> Warp 1 -> projection`.
   */
#if HS_ENABLE_SHADERBALL_DYNAMIC_BACKEND ||                                    \
    (HS_ENABLE_TEST_HOOKS && HS_ENABLE_TEST_ORACLES)
  HS_FLASH_MEMBER static PlanarWarpResult
  planar_warp_lookup(const ProjectedLookup &projected,
                     const FrameState &frame) {
    const bool path_length_required = tracks_displacement(frame);
    const PlanarWarpStageResult outer = warp_stage_lookup(
        projected.coords, projected, frame.slots.warp_program.outer,
        frame.params.warp.outer, frame.clocks.warp_outer_phase,
        frame.resources.outer_warp_noise, frame.prepared_warp.outer,
        path_length_required);
    const PlanarWarpStageResult inner = warp_stage_lookup(
        outer.coords, projected, frame.slots.warp_program.inner,
        frame.params.warp.inner, frame.clocks.warp_inner_phase,
        frame.resources.inner_warp_noise, frame.prepared_warp.inner,
        path_length_required);
    return {inner.coords, outer.path_length + inner.path_length};
  }
#endif

  HS_FLASH_MEMBER static PlanarWarpStageResult
  finish_closed_form_warp(const Complex &input, const Complex &output,
                          bool path_length_required) {
    return Pullback::Warp::finish_closed_form(input, output,
                                              path_length_required);
  }

  HS_FLASH_MEMBER static PlanarWarpStageResult
  warp_affine_frame(const Complex &input, const WarpStageParams &,
                    const PreparedWarpStage &prepared,
                    bool path_length_required) {
    return Pullback::Warp::affine_frame(input, prepared, path_length_required);
  }

  HS_FLASH_MEMBER static PlanarWarpStageResult
  warp_wave_shear(const Complex &input, const WarpStageParams &params,
                  float stage_phase, float amplitude,
                  const PreparedWarpStage &prepared,
                  bool path_length_required) {
    return Pullback::Warp::wave_shear(input, params, stage_phase, amplitude,
                                      prepared, path_length_required);
  }

  HS_FLASH_MEMBER static PlanarWarpStageResult
  warp_vortex(const Complex &input, const PreparedWarpStage &prepared,
              bool path_length_required) {
    return Pullback::Warp::vortex(input, prepared, path_length_required);
  }

  HS_FLASH_MEMBER static PlanarWarpStageResult
  warp_vector_noise(const Complex &input, const WarpStageSpec &spec,
                    const WarpStageParams &params, float amplitude,
                    const FastNoiseLite &noise,
                    const PreparedWarpStage &prepared,
                    bool path_length_required) {
    return Pullback::Warp::vector_noise(input, params, amplitude, noise,
                                        spec.basis, prepared,
                                        path_length_required);
  }

  HS_FLASH_MEMBER static PlanarWarpStageResult
  warp_curl_flow(const Complex &input, const WarpStageSpec &spec,
                 const WarpStageParams &params, float amplitude,
                 const FastNoiseLite &noise, float stage_phase,
                 bool path_length_required) {
    const uint8_t intervals =
        spec.curl_integrator == CurlIntegrator::EULER_1      ? 1
        : spec.curl_integrator == CurlIntegrator::MIDPOINT_2 ? 2
                                                             : 4;
    return Pullback::Warp::curl_flow(input, noise, spec.basis, intervals,
                                     params.scale, amplitude, stage_phase,
                                     path_length_required);
  }

  HS_FLASH_MEMBER static PlanarWarpStageResult
  warp_polar_chart(const Complex &input, const WarpStageSpec &spec,
                   const WarpStageParams &params, float stage_phase,
                   bool path_length_required) {
    return Pullback::Warp::polar_chart(
        input, params, stage_phase, spec.polar_mode == PolarMode::LOGARITHMIC,
        spec.polar_harmonic, path_length_required);
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
   * @param path_length_required Whether the frame's colorizer reads the
   *        displacement scalar.
   * @return Stage output coordinates, the delta it applied, and the path length
   *         travelled.
   * @details Path length is the direct displacement for the closed-form kinds
   * and the integrated arc length for curl flow. It is zero when
   * @p path_length_required is false.
   */
  HS_FLASH_MEMBER static PlanarWarpStageResult
  warp_stage_lookup(const Complex &input, const ProjectedLookup &projected,
                    const WarpStageSpec &spec, const WarpStageParams &params,
                    float stage_phase, const FastNoiseLite *stage_noise,
                    const PreparedWarpStage &prepared,
                    bool path_length_required) {
    if (spec.kind == WarpStageKind::NONE)
      return {input, Complex(), 0.0f};
    const float envelope =
        warp_envelope(projected, spec.envelope, params.edge_width);
    const float amplitude = params.strength * envelope;
    switch (spec.kind) {
    case WarpStageKind::NONE:
    case WarpStageKind::LEGACY_STEREO_NOISE:
      break;
    case WarpStageKind::AFFINE_FRAME:
      return warp_affine_frame(input, params, prepared, path_length_required);
    case WarpStageKind::WAVE_SHEAR:
      return warp_wave_shear(input, params, stage_phase, amplitude, prepared,
                             path_length_required);
    case WarpStageKind::VORTEX:
      return warp_vortex(input, prepared, path_length_required);
    case WarpStageKind::VECTOR_NOISE:
      if (amplitude == 0.0f)
        return {input, Complex(), 0.0f};
      HS_CHECK(stage_noise != nullptr,
               "ShaderBall vector warp has no noise resource");
      return warp_vector_noise(input, spec, params, amplitude, *stage_noise,
                               prepared, path_length_required);
    case WarpStageKind::CURL_FLOW:
      if (amplitude == 0.0f)
        return {input, Complex(), 0.0f};
      HS_CHECK(stage_noise != nullptr,
               "ShaderBall curl warp has no noise resource");
      return warp_curl_flow(input, spec, params, amplitude, *stage_noise,
                            stage_phase, path_length_required);
    case WarpStageKind::MIRROR_TILE: {
      HS_SB_STAGE_MARK(mirror_start);
      const PlanarWarpStageResult result = finish_closed_form_warp(
          input, mirror_tile(input, params, prepared), path_length_required);
      HS_SB_STAGE_SPAN(mirror_tile, mirror_start);
      return result;
    }
    case WarpStageKind::POLAR_CHART:
      return warp_polar_chart(input, spec, params, stage_phase,
                              path_length_required);
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

  HS_FLASH_MEMBER static Complex
  curl_flow(const Complex &input, const FastNoiseLite &noise,
            const WarpStageSpec &spec, const WarpStageParams &params,
            float distance, float phase, float &path_length) {
    const uint8_t intervals =
        spec.curl_integrator == CurlIntegrator::EULER_1      ? 1
        : spec.curl_integrator == CurlIntegrator::MIDPOINT_2 ? 2
                                                             : 4;
    constexpr bool path_length_required = true;
    const PlanarWarpStageResult result = Pullback::Warp::curl_flow(
        input, noise, spec.basis, intervals, params.scale, distance, phase,
        path_length_required);
    path_length = result.path_length;
    return result.coords;
  }

  HS_FLASH_MEMBER static Complex curl_vector(const Complex &p,
                                             const FastNoiseLite &noise,
                                             NoiseBasis basis, float scale,
                                             float phase) {
    return Pullback::Warp::curl_vector(p, noise, basis, scale, phase);
  }

  __attribute__((always_inline)) static Complex
  mirror_tile(const Complex &input, const WarpStageParams &params,
              const PreparedWarpStage &prepared) {
    return Pullback::Warp::mirror_tile_coords(input, params, prepared);
  }

#if HS_ENABLE_SHADERBALL_DYNAMIC_BACKEND ||                                    \
    (HS_ENABLE_TEST_HOOKS && HS_ENABLE_TEST_ORACLES)
  static Complex condition_source_coords(const Complex &coords,
                                         const FrameState &frame) {
    if (is_noise_contour(frame.slots.function) ||
        frame.slots.function == Function::PRIMITIVE_LATTICE)
      return coords;
    return stereo_pattern_args(coords, frame.params.source.pattern_freq);
  }

  /**
   * @brief Turns the signed source field into a shaped value and a coverage.
   * @param field Signed source sample, nominally in [-1, 1].
   * @param projected Projection output; supplies the signal weight, edge
   *        distance, and out-of-domain coverage.
   * @param warped Warp output from the coordinate stage.
   * @param frame Frame snapshot.
   * @return Value and coverage in [0, 1], plus the surface coordinate.
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
    return {value, coverage, projected.sphere,
            projected.surface_path_length + warped.path_length};
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
      return {value, projected.domain_coverage, projected.sphere,
              projected.surface_path_length + warped.path_length};
    return shape_nontrivial_material(value, projected, warped, frame);
  }
#endif

  /**
   * @brief Samples the selected noise-contour source coordinate.
   * @param q Projected or sphere noise coordinate.
   * @param frame Frame snapshot.
   * @return Signed field value in [-1, 1].
   */
  HS_FLASH_MEMBER static float sample_noise_contour(const Vector &q,
                                                    const FrameState &frame) {
    return Pullback::Source::noise_contour(*frame.resources.source_noise,
                                           frame.params.source.noise_basis, q,
                                           frame.params.source.noise_contrast);
  }

#if HS_ENABLE_SHADERBALL_DYNAMIC_BACKEND ||                                    \
    (HS_ENABLE_TEST_HOOKS && HS_ENABLE_TEST_ORACLES)
  HS_FLASH_MEMBER static float sample_source(const Complex &p,
                                             const ProjectedLookup &projected,
                                             const FrameState &frame) {
    if (frame.slots.function == Function::GRID)
      return grid(p, frame.params.source, frame.prepared_source);
    if (frame.slots.function == Function::NOISE_CONTOUR)
      return sample_noise_contour(
          noise_projected_coordinate(p, frame.params.source.noise_scale,
                                     frame.clocks.source_noise_time),
          frame);
    if (frame.slots.function == Function::NOISE_CONTOUR_SPHERE)
      return sample_noise_contour(
          noise_sphere_coordinate(projected.sphere,
                                  frame.params.source.noise_scale,
                                  frame.clocks.source_noise_time),
          frame);
    if (frame.slots.function == Function::PRIMITIVE_LATTICE)
      return primitive_lattice(p, frame.params.source);
    return sample_function(frame.slots.function, p, frame.prepared_source);
  }
#endif

  HS_FLASH_MEMBER static float primitive_lattice(const Complex &p,
                                                 const SourceParams &params) {
    return Pullback::Source::primitive_lattice(p, params);
  }

  /**
   * @brief Maps a shaped material sample to a palette colour.
   * @param sample Shaped value, coverage, and surface coordinate.
   * @param frame Frame snapshot.
   * @return Colour whose alpha is the palette alpha scaled by the coverage.
   */
  HS_FLASH_MEMBER static Color4 colorize_generated(const MaterialSample &sample,
                                                   const FrameState &frame) {
    return Pullback::Color::GeneratedPalette<ColorStateProvider>::apply(sample,
                                                                        frame);
  }

  __attribute__((always_inline)) static float
  palette_mapping_coordinate(float value, PaletteMapping mapping,
                             float frequency, float offset) {
    return Pullback::Color::palette_mapping_coordinate(
        value, static_cast<Pullback::Color::PaletteMapping>(mapping), frequency,
        offset);
  }

  __attribute__((always_inline)) static float
  brightness_envelope_gain(float value, BrightnessEnvelope envelope,
                           float depth) {
    return Pullback::Color::brightness_envelope_gain(
        value, static_cast<Pullback::Color::BrightnessEnvelope>(envelope),
        depth);
  }

#if HS_ENABLE_SHADERBALL_DYNAMIC_BACKEND ||                                    \
    (HS_ENABLE_TEST_HOOKS && HS_ENABLE_TEST_ORACLES)
  HS_O3_FN static Color4 colorize(const MaterialSample &sample,
                                  const FrameState &frame) {
    return colorize_generated(sample, frame);
  }
#endif

  HS_FLASH_MEMBER static void
  prepare_hue_rotation_lut(PreparedHueRotation &prepared,
                           const BakedPalette &palette) {
    Pullback::Color::prepare_hue_rotation_lut(
        std::span<Pixel, Pullback::Color::HueRotationLutView::SIZE>(
            prepared.lut, PreparedHueRotation::LUT_SIZE),
        palette);
  }

  HS_FLASH_MEMBER static Vector hue_noise_face_direction(int face, float u,
                                                         float v) {
    return Pullback::Color::hue_noise_face_direction(face, u, v);
  }

  HS_FLASH_MEMBER static void prepare_hue_noise_lut(PreparedHueNoise &prepared,
                                                    const FastNoiseLite &noise,
                                                    float scale, float phase) {
    Pullback::Color::prepare_hue_noise_lut(
        std::span<int8_t, Pullback::Color::HueNoiseLutView::SIZE>(
            prepared.lut, PreparedHueNoise::LUT_SIZE),
        noise, scale, phase);
  }

  HS_FLASH_MEMBER static float
  sample_hue_noise_lut(const PreparedHueNoise &prepared, const Vector &v) {
    return Pullback::Color::sample_hue_noise_lut(
        {prepared.lut, prepared.active}, v);
  }

  __attribute__((always_inline)) static Pixel
  sample_hue_rotation_lut(const PreparedHueRotation &prepared, float value,
                          float amount) {
    return Pullback::Color::sample_hue_rotation_lut(
        {prepared.lut, prepared.active}, value, amount);
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
      break;
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

  HS_FLASH_MEMBER static Vector surface_curl_field(const Vector &v,
                                                   const FrameState &frame) {
    return Pullback::Surface::curl_field(
        v, *frame.resources.surface_noise, frame.params.surface_noise.basis,
        frame.params.surface_noise.scale,
        frame.prepared_surface_noise.loop_offset);
  }

  HS_FLASH_MEMBER static SurfaceNoiseResult
  finish_surface_noise_step(const Vector &v, const Vector &step,
                            const FrameState &frame) {
    return Pullback::Surface::finish_step(v, step, tracks_displacement(frame));
  }

  HS_FLASH_MEMBER static SurfaceNoiseResult
  apply_simplex_euler_surface_noise_result(const Vector &v,
                                           const FrameState &frame) {
    return Pullback::Surface::curl_noise(
        v, *frame.resources.surface_noise, NoiseBasis::SIMPLEX,
        Pullback::Surface::Integrator::EULER, frame.params.surface_noise.scale,
        frame.prepared_surface_noise.loop_offset,
        frame.params.surface_noise.strength, tracks_displacement(frame));
  }

  HS_FLASH_MEMBER static SurfaceNoiseResult
  midpoint_surface_curl_step(const Vector &v, float distance,
                             const FrameState &frame) {
    return Pullback::Surface::curl_midpoint_step(
        v, *frame.resources.surface_noise, frame.params.surface_noise.basis,
        frame.params.surface_noise.scale,
        frame.prepared_surface_noise.loop_offset, distance,
        tracks_displacement(frame));
  }

  HS_FLASH_MEMBER static SurfaceNoiseResult
  apply_surface_noise_result(const Vector &v, const FrameState &frame) {
    const SurfaceNoiseParams &params = frame.params.surface_noise;
    const bool path_length_required = tracks_displacement(frame);
    if (frame.slots.surface_noise == SurfaceNoise::DIRECT) {
      return Pullback::Surface::direct_noise(
          v, *frame.resources.surface_noise, params.basis, params.scale,
          frame.prepared_surface_noise.loop_offset, params.strength,
          frame.prepared_surface_noise.direction_cos,
          frame.prepared_surface_noise.direction_sin, path_length_required);
    }
    const Pullback::Surface::Integrator integrator =
        params.integrator == SurfaceCurlIntegrator::EULER
            ? Pullback::Surface::Integrator::EULER
        : params.integrator == SurfaceCurlIntegrator::MIDPOINT
            ? Pullback::Surface::Integrator::MIDPOINT
            : Pullback::Surface::Integrator::MIDPOINT_2X;
    return Pullback::Surface::curl_noise(
        v, *frame.resources.surface_noise, params.basis, integrator,
        params.scale, frame.prepared_surface_noise.loop_offset, params.strength,
        path_length_required);
  }

  HS_FLASH_MEMBER static Vector apply_surface_noise(const Vector &v,
                                                    const FrameState &frame) {
    return apply_surface_noise_result(v, frame).sphere;
  }

  /**
   * @brief Applies a lens whose image depends on the direction alone.
   * @param v Unit direction on the sphere.
   * @param lens Lens to apply; MOBIUS reads FrameState parameters and is
   *        rejected here.
   * @return The lensed direction.
   */
  __attribute__((always_inline)) static Vector
  apply_frame_free_lens(const Vector &v, SurfaceLens lens) {
    switch (lens) {
    case SurfaceLens::NONE:
      return v;
    case SurfaceLens::GLITCH:
      return lenses::glitch_lens(v);
    case SurfaceLens::TWIST:
      return lenses::twist_lens(v);
    case SurfaceLens::KALEIDOSCOPE:
      return lenses::kaleidoscope_lens(v);
    case SurfaceLens::KALEIDOSCOPE_TETRAHEDRAL:
      return lenses::polyhedral_kaleidoscope_lens(v,
                                                  lenses::TETRAHEDRAL_MIRRORS);
    case SurfaceLens::KALEIDOSCOPE_OCTAHEDRAL:
      return lenses::polyhedral_kaleidoscope_lens(v,
                                                  lenses::OCTAHEDRAL_MIRRORS);
    case SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL:
      return lenses::dodecahedral_kaleidoscope_lens(v);
    case SurfaceLens::KALEIDOSCOPE_TRIANGULAR_PRISM:
      return lenses::polyhedral_kaleidoscope_lens(
          v, lenses::TRIANGULAR_PRISM_MIRRORS);
    case SurfaceLens::KALEIDOSCOPE_SQUARE_PRISM:
      return lenses::polyhedral_kaleidoscope_lens(v,
                                                  lenses::SQUARE_PRISM_MIRRORS);
    case SurfaceLens::KALEIDOSCOPE_PENTAGONAL_PRISM:
      return lenses::polyhedral_kaleidoscope_lens(
          v, lenses::PENTAGONAL_PRISM_MIRRORS);
    case SurfaceLens::KALEIDOSCOPE_HEXAGONAL_PRISM:
      return lenses::polyhedral_kaleidoscope_lens(
          v, lenses::HEXAGONAL_PRISM_MIRRORS);
    case SurfaceLens::KALEIDOSCOPE_OCTAGONAL_PRISM:
      return lenses::polyhedral_kaleidoscope_lens(
          v, lenses::OCTAGONAL_PRISM_MIRRORS);
    case SurfaceLens::MOBIUS:
    case SurfaceLens::TANGENT_NOISE:
      HS_CHECK(false, "frame-parameterized lens needs the FrameState overload");
      __builtin_unreachable();
    }
    __builtin_unreachable();
  }

  HS_FLASH_MEMBER static Vector mobius_lens(const Vector &v,
                                            const MobiusParams &params) {
    return Pullback::Lens::mobius(v, params);
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
      return projections::folded_sinusoidal(v);
    case Projection::EQUIRECTANGULAR:
      return projections::equirectangular(v);
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
    case Function::NOISE_CONTOUR_SPHERE:
    case Function::PRIMITIVE_LATTICE:
      break;
    }
    __builtin_unreachable();
  }

  HS_FLASH_MEMBER static float twin_wave(const Complex &p,
                                         const SourceState &source) {
    return Pullback::Source::twin_wave(p, source);
  }

  HS_FLASH_MEMBER static float rings(const Complex &p,
                                     const SourceState &source) {
    return Pullback::Source::rings(p, source);
  }

  HS_FLASH_MEMBER static float spiral(const Complex &p,
                                      const SourceState &source) {
    return Pullback::Source::spiral(p, source);
  }

  HS_O3_FN static float grid(const Complex &p, const SourceParams &params,
                             const SourceState &source) {
    return Pullback::Source::grid(p, params, source);
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
    look.clocks.hue_noise_phase =
        wrap_t(look.clocks.hue_noise_phase + params.color.hue_noise_speed);
    look.clocks.source_noise_time =
        wrap_t(look.clocks.source_noise_time + params.source.noise_time_rate);
    look.clocks.surface_noise_time =
        wrap_t(look.clocks.surface_noise_time + params.surface_noise.rate);
    if (config.slots.warp_program.outer.kind == WarpStageKind::AFFINE_FRAME)
      look.clocks.warp_outer_rotation =
          TWO_PI_F *
          wrap_t((look.clocks.warp_outer_rotation +
                  params.warp.outer.speed * params.warp.outer.rotation) /
                 TWO_PI_F);
    if (config.slots.warp_program.inner.kind == WarpStageKind::AFFINE_FRAME)
      look.clocks.warp_inner_rotation =
          TWO_PI_F *
          wrap_t((look.clocks.warp_inner_rotation +
                  params.warp.inner.speed * params.warp.inner.rotation) /
                 TWO_PI_F);
    look.clocks.warp_outer_phase =
        wrap_t(look.clocks.warp_outer_phase + params.warp.outer.speed);
    look.clocks.warp_inner_phase =
        wrap_t(look.clocks.warp_inner_phase + params.warp.inner.speed);
    look.clocks.palette_oscillation_phase =
        wrap_t(look.clocks.palette_oscillation_phase +
               params.color.phase_oscillation_speed);
    update_spatial_frames(look, config, deltas);
  }

  HS_COLD_MEMBER void prepare_param_morph() {
    if (!state->param_morph.active)
      return;
    const float mix =
        transition_mix(state->param_morph.elapsed, state->param_morph.duration);
    blend.palette_mapping = PaletteMappingWeights::lerp(
        state->param_morph.mapping_from, state->param_morph.mapping_to, mix);
    if (mix == 0.0f)
      blend.params = state->param_morph.from;
    else if (mix == 1.0f)
      blend.params = state->param_morph.to;
    else if (state->param_morph.staggered)
      blend.params.lerp_staggered(state->param_morph.from,
                                  state->param_morph.to, mix, active_slots);
    else
      blend.params.lerp(state->param_morph.from, state->param_morph.to, mix,
                        active_slots);
  }

  static float transition_mix(uint16_t elapsed, uint16_t duration) {
    if (elapsed == 0)
      return 0.0f;
    if (elapsed >= duration)
      return 1.0f;
    return ease_in_out_sin(static_cast<float>(elapsed) / duration);
  }

  __attribute__((noinline)) HS_COLD_MEMBER void apply_requested_config() {
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
    if (!admissible_config(next_config)) {
      reject_requested_config();
      return;
    }
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
    active_pipeline = resolve_pipeline_id(next_config);
    blend.params = next_config.params;
    blend.palette_mapping =
        palette_mapping_weights(next_config.slots.palette_mapping);
#if HS_ENABLE_PARAM_GUI_BRIDGE
    display_config = next_config;
#endif
    published_config = next_config;
#if HS_ENABLE_PARAM_GUI_BRIDGE
    accepted_config = next_config;
#endif
    if (!requested_schema_bound)
      rebind_parameters();
  }

  HS_COLD_MEMBER void reject_requested_config() {
    requested_config = published_config;
#if HS_ENABLE_PARAM_GUI_BRIDGE
    accepted_config = published_config;
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
    if (!admissible_config(candidate) || duration == 0)
      return false;
    if (state->transition.active)
      return false;
    Config &current = state->render_config;
    current.slots = active_slots;
    current.params = blend.params;
    if (!transition_admitted(current, candidate))
      return false;
    if (candidate == current) {
      state->param_morph.active = false;
      blend.palette_mapping =
          palette_mapping_weights(current.slots.palette_mapping);
      return true;
    }
    if (stable_topology(current, candidate)) {
      if (!prepare_resource_union(current, candidate))
        return false;
      state->param_morph = {
          current.params,
          candidate.params,
          blend.palette_mapping,
          palette_mapping_weights(candidate.slots.palette_mapping),
          candidate.slots.palette_mapping,
          0,
          duration,
          staggered,
          continue_choreo,
          true};
      return true;
    }
    if (!prepare_resource_union(current, current))
      return false;
    const uint16_t planned_duration =
        (duration & 1U) != 0 ? duration + 1 : duration;
    state->param_morph.active = false;
    state->transition = {current, candidate,        runtime,         runtime,
                         0,       planned_duration, continue_choreo, true};
    state->transition.from_pipeline = active_pipeline;
    state->transition.to_pipeline = resolve_pipeline_id(candidate);
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
      active_pipeline = state->transition.to_pipeline;
      blend.params = state->transition.to_config.params;
      blend.palette_mapping =
          palette_mapping_weights(active_slots.palette_mapping);
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
    active_slots.palette_mapping = state->param_morph.mapping_destination;
    blend.palette_mapping = state->param_morph.mapping_to;
    state->param_morph.active = false;
    if (continue_choreo)
      enter_preset();
  }

  __attribute__((noinline)) HS_COLD_MEMBER void publish_live_config() {
    if (anims_paused || state->transition.active || state->param_morph.active)
      return;
#if HS_ENABLE_PARAM_GUI_BRIDGE
    if (accepted_config != published_config)
      return;
#endif
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
#if HS_ENABLE_PARAM_GUI_BRIDGE
    accepted_config = published_config;
#endif
  }

#if HS_ENABLE_PARAM_GUI_BRIDGE
  HS_COLD_MEMBER void refresh_parameter_display() override {
    if (state->transition.active) {
      const float mix =
          transition_mix(state->transition.elapsed, state->transition.duration);
      display_config.slots = mix < 0.5f ? state->transition.from_config.slots
                                        : state->transition.to_config.slots;
      display_config.params.lerp(state->transition.from_config.params,
                                 state->transition.to_config.params, mix,
                                 display_config.slots);
      return;
    }
    display_config = {active_slots, blend.params};
  }
#endif

  template <typename Enum>
  HS_COLD_MEMBER static constexpr bool enum_at_most(Enum value, Enum last) {
    return static_cast<uint8_t>(value) <= static_cast<uint8_t>(last);
  }

  HS_COLD_MEMBER static constexpr bool is_noise_contour(Function function) {
    return function == Function::NOISE_CONTOUR ||
           function == Function::NOISE_CONTOUR_SPHERE;
  }

  HS_COLD_MEMBER static constexpr SourceTraits
  source_traits(Function function) {
    switch (function) {
    case Function::GRID:
      return {true, true};
    case Function::PRIMITIVE_LATTICE:
      return {true, true};
    case Function::TWIN_WAVE:
    case Function::RINGS:
    case Function::SPIRAL:
    case Function::NOISE_CONTOUR:
    case Function::NOISE_CONTOUR_SPHERE:
      return {false, false};
    }
    return {false, false};
  }

  HS_COLD_MEMBER static constexpr Complex
  source_cartesian_period(const Config &config) {
    if (config.slots.function != Function::PRIMITIVE_LATTICE ||
        !(config.params.source.lattice_cell_scale > 0.0f))
      return {};
    const float period = 1.0f / config.params.source.lattice_cell_scale;
    return {period, period};
  }

  HS_COLD_MEMBER static constexpr bool
  affine_has_translation(const WarpStageSpec &spec,
                         const WarpStageParams &params) {
    return spec.kind == WarpStageKind::AFFINE_FRAME &&
           (params.translation_x != 0.0f || params.translation_y != 0.0f);
  }

  HS_COLD_MEMBER static constexpr bool whole_affine_winding(float value) {
    if (!(value >= -4.0f && value <= 4.0f))
      return false;
    return value == static_cast<float>(static_cast<int>(value));
  }

  HS_COLD_MEMBER static constexpr float
  hue_shift_amount_max(HueShiftMode mode) {
    return mode == HueShiftMode::NOISE ? HUE_NOISE_AMOUNT_MAX
                                       : HUE_SHIFT_AMOUNT_MAX;
  }

  HS_COLD_MEMBER static constexpr bool
  hue_shift_amount_in_range(const Config &config) {
    return config.params.color.hue_shift_amount >=
               -hue_shift_amount_max(config.slots.hue_shift) &&
           config.params.color.hue_shift_amount <=
               hue_shift_amount_max(config.slots.hue_shift);
  }

  HS_COLD_MEMBER static constexpr bool
  affine_translation_compatible(const Config &config) {
    const bool outer = affine_has_translation(config.slots.warp_program.outer,
                                              config.params.warp.outer);
    const bool inner = affine_has_translation(config.slots.warp_program.inner,
                                              config.params.warp.inner);
    if (!outer && !inner)
      return true;
    if (config.slots.function != Function::PRIMITIVE_LATTICE ||
        (outer &&
         config.slots.warp_program.inner.kind != WarpStageKind::NONE) ||
        (config.slots.hue_shift == HueShiftMode::WARP_DISPLACEMENT &&
         config.params.color.hue_shift_amount != 0.0f))
      return false;
    return (!outer ||
            (whole_affine_winding(config.params.warp.outer.translation_x) &&
             whole_affine_winding(config.params.warp.outer.translation_y))) &&
           (!inner ||
            (whole_affine_winding(config.params.warp.inner.translation_x) &&
             whole_affine_winding(config.params.warp.inner.translation_y)));
  }

  /// Source periods spanned by one angular seam jump of a polar chart.
  /// Primitive Lattice is periodic in its own cell scale and ignores the
  /// pattern frequency, so its seam spans `2*pi*harmonic*cell_scale` cells.
  HS_COLD_MEMBER static constexpr float
  polar_seam_periods(const RequestedConfig &config,
                     const WarpStageSpec &polar) {
    const float harmonic = static_cast<float>(polar.polar_harmonic);
    if (config.slots.function == Function::PRIMITIVE_LATTICE)
      return TWO_PI_F * harmonic * config.params.source.lattice_cell_scale;
    return config.params.source.pattern_freq * harmonic;
  }

  HS_COLD_MEMBER static constexpr bool
  polar_source_compatible(const RequestedConfig &config,
                          const WarpStageSpec &polar) {
    const SourceTraits traits = source_traits(config.slots.function);
    if (!traits.y_periodic || !traits.polar_angle_compatible)
      return false;
    const float periods = polar_seam_periods(config, polar);
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
    if (!valid_slot_enums(slots))
      return false;
    if (slots.function == Function::NOISE_CONTOUR_SPHERE &&
        (slots.warp_program.outer.kind != WarpStageKind::NONE ||
         slots.warp_program.inner.kind != WarpStageKind::NONE))
      return false;
    if (slots.surface_lens == SurfaceLens::TANGENT_NOISE ||
        slots.warp_program.outer.kind == WarpStageKind::LEGACY_STEREO_NOISE ||
        slots.warp_program.inner.kind == WarpStageKind::LEGACY_STEREO_NOISE)
      return false;
    const bool outer_polar =
        slots.warp_program.outer.kind == WarpStageKind::POLAR_CHART;
    const bool inner_polar =
        slots.warp_program.inner.kind == WarpStageKind::POLAR_CHART;
    if ((outer_polar && slots.warp_program.inner.kind != WarpStageKind::NONE &&
         slots.warp_program.inner.kind != WarpStageKind::WAVE_SHEAR) ||
        (inner_polar && slots.warp_program.outer.kind != WarpStageKind::NONE))
      return false;
    if (inner_polar &&
        !polar_source_compatible(candidate, slots.warp_program.inner))
      return false;
    if (outer_polar &&
        slots.warp_program.inner.kind == WarpStageKind::WAVE_SHEAR) {
      const SourceTraits traits = source_traits(slots.function);
      if (!traits.y_periodic || !traits.polar_angle_compatible)
        return false;
    }
    if (outer_polar && slots.warp_program.inner.kind == WarpStageKind::NONE &&
        !polar_source_compatible(candidate, slots.warp_program.outer))
      return false;
    if (!affine_translation_compatible(candidate) ||
        !strict_seam_compatible(candidate) ||
        !preset_in_ranges(candidate.params) ||
        !hue_shift_amount_in_range(candidate))
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
    const SurfaceNoiseParams &surface_noise = candidate.params.surface_noise;
    if (!enum_at_most(surface_noise.basis, NoiseBasis::RIDGED3) ||
        !enum_at_most(surface_noise.integrator,
                      SurfaceCurlIntegrator::MIDPOINT_2X) ||
        surface_noise.scale < LENS_NOISE_SCALE_MIN ||
        surface_noise.scale > LENS_NOISE_SCALE_MAX ||
        surface_noise.strength <
            (slots.surface_noise == SurfaceNoise::CURL ? -0.5f : 0.0f) ||
        surface_noise.strength > 0.5f || surface_noise.rate < NOISE_RATE_MIN ||
        surface_noise.rate > NOISE_RATE_MAX || surface_noise.direction < 0.0f ||
        surface_noise.direction > 1.0f)
      return false;
    return resource_union_fits(candidate, candidate);
  }

  /**
   * @brief Admission test for a requested configuration.
   * @details The simulator accepts every structurally valid configuration;
   * device builds additionally require a compiled inverse pipeline.
   */
  HS_COLD_MEMBER static bool
  admissible_config(const RequestedConfig &candidate) {
    if (!valid_config(candidate))
      return false;
#if HS_ENABLE_SHADERBALL_DYNAMIC_BACKEND
    return true;
#else
    return find_inverse_program(candidate) != nullptr;
#endif
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
    case WarpStageKind::LEGACY_STEREO_NOISE:
      break;
    case WarpStageKind::AFFINE_FRAME:
      append_range_warning("Translate X", params.translation_x, -4.0f, 4.0f);
      append_range_warning("Translate Y", params.translation_y, -4.0f, 4.0f);
      append_range_warning("Rotation", params.rotation, -PI_F / 8.0f,
                           PI_F / 8.0f);
      append_range_warning("Scale X", params.scale_x, 0.25f, 4.0f);
      append_range_warning("Scale Y", params.scale_y, 0.25f, 4.0f);
      append_range_warning("Shear", params.shear, -0.75f, 0.75f);
      break;
    case WarpStageKind::WAVE_SHEAR:
      append_range_warning("Warp Strength", params.strength, -4.0f, 4.0f);
      append_range_warning("Frequency", params.frequency, 0.0f, 64.0f);
      append_range_warning("Warp Speed", params.speed, NOISE_SPEED_MIN,
                           NOISE_SPEED_MAX);
      break;
    case WarpStageKind::VORTEX:
      append_range_warning("Radius", params.radius, 1.0f / 64.0f, 8.0f);
      append_range_warning("Turns", params.turns, -4.0f, 4.0f);
      append_range_warning("Orbit Radius", params.center_orbit_radius, 0.0f,
                           4.0f);
      append_range_warning("Warp Speed", params.speed, NOISE_SPEED_MIN,
                           NOISE_SPEED_MAX);
      break;
    case WarpStageKind::VECTOR_NOISE:
      append_range_warning("Warp Strength", params.strength, 0.0f,
                           VECTOR_WARP_STRENGTH_MAX);
      append_range_warning("Warp Scale", params.scale, 1.0f / 64.0f,
                           VECTOR_WARP_SCALE_MAX);
      append_range_warning("Warp Speed", params.speed, NOISE_SPEED_MIN,
                           NOISE_SPEED_MAX);
      break;
    case WarpStageKind::CURL_FLOW: {
      append_range_warning("Warp Strength", params.strength,
                           -CURL_WARP_STRENGTH_MAX, CURL_WARP_STRENGTH_MAX);
      append_range_warning("Warp Scale", params.scale, 1.0f / 64.0f,
                           CURL_WARP_SCALE_MAX);
      append_range_warning("Warp Speed", params.speed, NOISE_SPEED_MIN,
                           NOISE_SPEED_MAX);
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
      append_range_warning("Rotation", params.rotation, 0.0f, TWO_PI_F);
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
    const Complex source_period = source_cartesian_period(candidate);
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
      bound = stage_coordinate_bound(stages[index], params[index], bound,
                                     source_period);
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
    std::array<NoiseFieldKey, MAX_NOISE_RESOURCES> keys{};
    size_t count = 0;
    auto add = [&](const NoiseFieldKey &key) {
      return !append_resource_key(key, keys, count);
    };
    if (warp_uses_noise(candidate.slots.warp_program.outer.kind) &&
        add(warp_resource_key(candidate.slots.warp_program.outer)))
      return warning_text.data();
    if (warp_uses_noise(candidate.slots.warp_program.inner.kind) &&
        add(warp_resource_key(candidate.slots.warp_program.inner)))
      return warning_text.data();
    if (is_noise_contour(candidate.slots.function) &&
        add(source_resource_key(candidate)))
      return warning_text.data();
    if (candidate.slots.surface_noise != SurfaceNoise::NONE &&
        add(surface_noise_resource_key(candidate)))
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
    if (candidate.slots.function == Function::NOISE_CONTOUR_SPHERE &&
        outer.kind != WarpStageKind::NONE && inner.kind != WarpStageKind::NONE)
      return begin_warning(
          "Noise Contour (Sphere) rejects Planar Warp 1 %s and Planar Warp 2 "
          "%s. Set both warps to None, or select Noise Contour (Projected).",
          WARP_OPTIONS[static_cast<uint8_t>(outer.kind)],
          WARP_OPTIONS[static_cast<uint8_t>(inner.kind)]);
    if (candidate.slots.function == Function::NOISE_CONTOUR_SPHERE &&
        (outer.kind != WarpStageKind::NONE ||
         inner.kind != WarpStageKind::NONE)) {
      const bool outer_active = outer.kind != WarpStageKind::NONE;
      const char *position = outer_active ? "Planar Warp 1" : "Planar Warp 2";
      const WarpStageKind kind = outer_active ? outer.kind : inner.kind;
      return begin_warning(
          "Noise Contour (Sphere) rejects %s %s. Set %s to None, or select "
          "Noise Contour (Projected).",
          position, WARP_OPTIONS[static_cast<uint8_t>(kind)], position);
    }
    if (outer.kind == WarpStageKind::POLAR_CHART &&
        inner.kind != WarpStageKind::NONE &&
        inner.kind != WarpStageKind::WAVE_SHEAR)
      return begin_warning(
          "Planar Warp 1 Polar Chart cannot run while Planar Warp 2 is %s. Set "
          "Planar Warp 2 to None or Wave Shear, or choose a different Planar "
          "Warp 1.",
          WARP_OPTIONS[static_cast<uint8_t>(inner.kind)]);
    if (inner.kind == WarpStageKind::POLAR_CHART &&
        outer.kind != WarpStageKind::NONE)
      return begin_warning(
          "Planar Warp 2 Polar Chart cannot run while Planar Warp 1 is %s. Set "
          "Planar Warp 1 to None or choose a different Planar Warp 2.",
          WARP_OPTIONS[static_cast<uint8_t>(outer.kind)]);
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
      const float periods = polar_seam_periods(candidate, *polar);
      const float nearest_periods = floorf(periods + 0.5f);
      if (candidate.slots.function == Function::PRIMITIVE_LATTICE)
        return begin_warning(
            "%s Polar Chart requires 2*pi x Lattice Cell Scale x Polar "
            "Harmonic to be a whole number. %.7g x %u gives %.7g. Set Lattice "
            "Cell Scale to %.7g or change %s Polar Harmonic.",
            position,
            static_cast<double>(candidate.params.source.lattice_cell_scale),
            static_cast<unsigned>(polar->polar_harmonic),
            static_cast<double>(periods),
            static_cast<double>(
                nearest_periods /
                (TWO_PI_F * static_cast<float>(polar->polar_harmonic))),
            position);
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
    if (!affine_translation_compatible(candidate)) {
      const bool outer_scroll =
          affine_has_translation(outer, candidate.params.warp.outer);
      const char *position = outer_scroll ? "Planar Warp 1" : "Planar Warp 2";
      const WarpStageParams &params = outer_scroll
                                          ? candidate.params.warp.outer
                                          : candidate.params.warp.inner;
      if (!whole_affine_winding(params.translation_x) ||
          !whole_affine_winding(params.translation_y))
        return begin_warning(
            "%s Affine Frame translation must use whole source-cell windings. "
            "Set Translation X and Translation Y to whole numbers.",
            position);
      if (candidate.slots.function != Function::PRIMITIVE_LATTICE)
        return begin_warning(
            "%s Affine Frame translation requires an exactly periodic "
            "Function. Select Primitive Lattice or set both translations to "
            "zero.",
            position);
      if (outer_scroll && inner.kind != WarpStageKind::NONE)
        return begin_warning(
            "Planar Warp 1 Affine Frame translation cannot precede Planar "
            "Warp 2 %s because the later warp breaks its source-period seam. "
            "Set Planar Warp 2 to None or set both translations to zero.",
            WARP_OPTIONS[static_cast<uint8_t>(inner.kind)]);
      return begin_warning(
          "%s Affine Frame translation cannot drive Total Warp Displacement "
          "hue because its path length resets at the source-period seam. "
          "Select Hue Shift None or Noise, or set both translations to zero.",
          position);
    }
    if (!strict_seam_compatible(candidate)) {
      begin_warning(
          "Projection %s requires seam-safe stages.",
          PROJECTION_OPTIONS[static_cast<uint8_t>(candidate.slots.projection)]);
      if (candidate.slots.function == Function::NOISE_CONTOUR)
        append_warning(" Function Noise Contour (Projected) is not seam-safe.");
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
    const SurfaceNoiseParams &surface_noise = candidate.params.surface_noise;
    const float minimum_surface_strength =
        candidate.slots.surface_noise == SurfaceNoise::CURL ? -0.5f : 0.0f;
    if (surface_noise.scale < LENS_NOISE_SCALE_MIN ||
        surface_noise.scale > LENS_NOISE_SCALE_MAX ||
        surface_noise.strength < minimum_surface_strength ||
        surface_noise.strength > 0.5f || surface_noise.rate < NOISE_RATE_MIN ||
        surface_noise.rate > NOISE_RATE_MAX || surface_noise.direction < 0.0f ||
        surface_noise.direction > 1.0f) {
      begin_warning("Surface Noise %s rejected.",
                    SURFACE_NOISE_OPTIONS[static_cast<uint8_t>(
                        candidate.slots.surface_noise)]);
      append_range_warning("Surface Noise Scale", surface_noise.scale,
                           LENS_NOISE_SCALE_MIN, LENS_NOISE_SCALE_MAX);
      append_range_warning("Surface Noise Strength", surface_noise.strength,
                           minimum_surface_strength, 0.5f);
      append_range_warning("Surface Noise Rate", surface_noise.rate,
                           NOISE_RATE_MIN, NOISE_RATE_MAX);
      append_range_warning("Surface Noise Direction", surface_noise.direction,
                           0.0f, 1.0f);
      append_warning(" Set the named Surface Noise control within its range.");
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
    if (!HS_ENABLE_SHADERBALL_DYNAMIC_BACKEND &&
        find_inverse_program(candidate) == nullptr)
      return uncompiled_program_warning(candidate, edited_name);
    return begin_warning(
        "%s was rejected by an unclassified ShaderBall admission rule. Keep "
        "the requested value and report this exact configuration as a bug.",
        edited_name);
  }

  const char *uncompiled_program_warning(const Config &candidate,
                                         const char *edited_name) const {
    const TopologyKey key = make_topology_key(candidate);
    for (const ProgramDescriptor &program : inverse_programs())
      if (program.key == key)
        return begin_warning(
            "%s is outside what the compiled pipeline for this stage "
            "combination supports. Restore %s or change a stage.",
            edited_name, edited_name);
    return begin_warning(
        "This stage combination has no compiled pipeline. Restore %s or "
        "select a combination reachable from a preset.",
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

  /** @brief Maximum component emitted by the bounded curl vector field. */
  HS_COLD_MEMBER static constexpr float
  curl_vector_component_bound(NoiseBasis) {
    return CURL_VECTOR_COMPONENT_MAX;
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
    const float stable_limit =
        0.5f * static_cast<float>(curl_intervals(spec.curl_integrator)) /
        (scale * curl_vector_component_bound(spec.basis));
    return stable_limit < CURL_WARP_STRENGTH_MAX ? stable_limit
                                                 : CURL_WARP_STRENGTH_MAX;
  }

  HS_COLD_MEMBER static constexpr bool
  valid_stage_tuple(const WarpStageSpec &spec, const WarpStageParams &params) {
    switch (spec.kind) {
    case WarpStageKind::NONE:
      return true;
    case WarpStageKind::LEGACY_STEREO_NOISE:
      return false;
    case WarpStageKind::AFFINE_FRAME:
      return params.translation_x >= -4.0f && params.translation_x <= 4.0f &&
             params.translation_y >= -4.0f && params.translation_y <= 4.0f &&
             params.rotation >= -PI_F / 8.0f &&
             params.rotation <= PI_F / 8.0f && params.scale_x >= 0.25f &&
             params.scale_x <= 4.0f && params.scale_y >= 0.25f &&
             params.scale_y <= 4.0f && params.shear >= -0.75f &&
             params.shear <= 0.75f && params.speed >= NOISE_SPEED_MIN &&
             params.speed <= NOISE_SPEED_MAX;
    case WarpStageKind::WAVE_SHEAR:
      return params.strength >= -4.0f && params.strength <= 4.0f &&
             params.frequency >= 0.0f && params.frequency <= 64.0f &&
             params.speed >= NOISE_SPEED_MIN && params.speed <= NOISE_SPEED_MAX;
    case WarpStageKind::VORTEX:
      return params.radius >= 1.0f / 64.0f && params.radius <= 8.0f &&
             params.turns >= -4.0f && params.turns <= 4.0f &&
             params.center_orbit_radius >= 0.0f &&
             params.center_orbit_radius <= 4.0f &&
             params.speed >= NOISE_SPEED_MIN && params.speed <= NOISE_SPEED_MAX;
    case WarpStageKind::VECTOR_NOISE:
      return params.strength >= 0.0f &&
             params.strength <= VECTOR_WARP_STRENGTH_MAX &&
             params.scale >= 1.0f / 64.0f &&
             params.scale <= VECTOR_WARP_SCALE_MAX &&
             params.speed >= NOISE_SPEED_MIN && params.speed <= NOISE_SPEED_MAX;
    case WarpStageKind::CURL_FLOW:
      return params.strength >= -CURL_WARP_STRENGTH_MAX &&
             params.strength <= CURL_WARP_STRENGTH_MAX &&
             params.scale >= 1.0f / 64.0f &&
             params.scale <= CURL_WARP_SCALE_MAX &&
             params.speed >= NOISE_SPEED_MIN &&
             params.speed <= NOISE_SPEED_MAX &&
             params.scale * abs_value(params.strength) *
                     curl_vector_component_bound(spec.basis) /
                     curl_intervals(spec.curl_integrator) <=
                 0.5f;
    case WarpStageKind::MIRROR_TILE:
      return params.rotation >= 0.0f && params.rotation <= TWO_PI_F &&
             params.cell_x >= CELL_MIN && params.cell_x <= CELL_MAX &&
             params.cell_y >= CELL_MIN && params.cell_y <= CELL_MAX &&
             params.speed >= NOISE_SPEED_MIN && params.speed <= NOISE_SPEED_MAX;
    case WarpStageKind::POLAR_CHART:
      return params.radial_scale >= 1.0f / 64.0f &&
             params.radial_scale <= 16.0f && params.speed >= NOISE_SPEED_MIN &&
             params.speed <= NOISE_SPEED_MAX;
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
                         const WarpStageParams &params, float input_bound,
                         const Complex &source_period) {
    switch (spec.kind) {
    case WarpStageKind::NONE:
      return input_bound;
    case WarpStageKind::LEGACY_STEREO_NOISE:
      return WARP_COORD_LIMIT + 1.0f;
    case WarpStageKind::AFFINE_FRAME: {
      const float rotated = 1.414214f * input_bound;
      const float x_bound = rotated / params.scale_x +
                            abs_value(params.shear) * rotated / params.scale_y +
                            abs_value(params.translation_x * source_period.re);
      const float y_bound = rotated / params.scale_y +
                            abs_value(params.translation_y * source_period.im);
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
                               params.scale *
                               curl_vector_component_bound(spec.basis);
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
    const Complex source_period = source_cartesian_period(config);
    const WarpStageSpec stages[] = {config.slots.warp_program.outer,
                                    config.slots.warp_program.inner};
    const WarpStageParams params[] = {config.params.warp.outer,
                                      config.params.warp.inner};
    for (size_t index = 0; index < 2; ++index) {
      if ((stages[index].kind == WarpStageKind::VECTOR_NOISE ||
           stages[index].kind == WarpStageKind::CURL_FLOW) &&
          params[index].scale * (bound + 100.0f) > NOISE_LATTICE_LIMIT)
        return false;
      bound = stage_coordinate_bound(stages[index], params[index], bound,
                                     source_period);
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
    worst.params.source.lattice_cell_scale =
        min_value(from.params.source.lattice_cell_scale,
                  to.params.source.lattice_cell_scale);
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
  affine_winding_pair_stable(const WarpStageSpec &spec,
                             const WarpStageParams &a,
                             const WarpStageParams &b) {
    return spec.kind != WarpStageKind::AFFINE_FRAME ||
           (a.translation_x == b.translation_x &&
            a.translation_y == b.translation_y);
  }

  HS_COLD_MEMBER static constexpr bool
  stable_parameter_path_admitted(const Config &from, const Config &to) {
    return curl_pair_stable(from.slots.warp_program.outer,
                            from.params.warp.outer, to.params.warp.outer) &&
           curl_pair_stable(from.slots.warp_program.inner,
                            from.params.warp.inner, to.params.warp.inner) &&
           affine_winding_pair_stable(from.slots.warp_program.outer,
                                      from.params.warp.outer,
                                      to.params.warp.outer) &&
           affine_winding_pair_stable(from.slots.warp_program.inner,
                                      from.params.warp.inner,
                                      to.params.warp.inner) &&
           polar_pair_stable(from, to) && safe_program_path(from, to);
  }

  HS_COLD_MEMBER static constexpr bool
  same_parameter_topology(const Config &from, const Config &to) {
    Slots from_slots = from.slots;
    Slots to_slots = to.slots;
    from_slots.palette_mapping = PaletteMapping::LINEAR;
    to_slots.palette_mapping = PaletteMapping::LINEAR;
    return from_slots == to_slots &&
           from.params.source.noise_basis == to.params.source.noise_basis &&
           from.params.source.noise_seed == to.params.source.noise_seed &&
           from.params.source.noise_resource_id ==
               to.params.source.noise_resource_id &&
           from.params.value.band_count == to.params.value.band_count &&
           from.params.surface_noise.basis == to.params.surface_noise.basis &&
           from.params.surface_noise.integrator ==
               to.params.surface_noise.integrator &&
           from.params.surface_noise.seed == to.params.surface_noise.seed &&
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
    return scale * distance * curl_vector_component_bound(spec.basis) /
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
              static_cast<unsigned>(preset_count_for_view()));
#endif
    } else {
      preset_dwell_remaining = 1;
      preset_dwell_armed = true;
    }
  }

  static void next_generated_palette(uint32_t &hue, uint32_t sequence,
                                     PaletteHarmony harmony, float chroma,
                                     GenerativePalette &out) {
    if (sequence > 0)
      hue += HUE_STEP;
    out = GenerativePalette{PaletteRecipes::profile(
        PaletteDomain::STRAIGHT, harmony, AxisCurve::ASCENDING,
        PaletteRecipes::hue_turns(hue), chroma)};
  }

  static void next_triadic_palette(void *context, uint32_t sequence,
                                   GenerativePalette &out) {
    ShaderBall &effect = *static_cast<ShaderBall *>(context);
    next_generated_palette(effect.palette_hue, sequence,
                           PaletteHarmony::TRIADIC, effect.palette_chroma, out);
  }

  static void next_complementary_palette(void *context, uint32_t sequence,
                                         GenerativePalette &out) {
    ShaderBall &effect = *static_cast<ShaderBall *>(context);
    next_generated_palette(effect.complementary_hue, sequence,
                           PaletteHarmony::COMPLEMENTARY, effect.palette_chroma,
                           out);
  }

  static void next_analogous_palette(void *context, uint32_t sequence,
                                     GenerativePalette &out) {
    ShaderBall &effect = *static_cast<ShaderBall *>(context);
    next_generated_palette(effect.analogous_hue, sequence,
                           PaletteHarmony::ANALOGOUS, effect.palette_chroma,
                           out);
  }

  static constexpr uint8_t BOUNDARY_CUT = 1U << 0;
  static constexpr uint8_t BOUNDARY_SINGULAR = 1U << 1;
  static constexpr uint8_t PROJECTION_FLAG_FOLDED = 1U << 0;
  static constexpr float GNOMONIC_AXIS_EPS = 1e-3f;
  static constexpr float WARP_COORD_LIMIT = 65536.0f;
  static constexpr float NOISE_LATTICE_LIMIT = 1048576.0f;
  static constexpr float SPIRAL_ARMS = 3.0f;
  static constexpr float ONE_BELOW_UNIT = 0x1.fffffep-1f;
  static constexpr uint32_t HUE_STEP = 159;
  static constexpr int PALETTE_DWELL_FRAMES = 0;
  static constexpr int PALETTE_FADE_FRAMES = 600;
  static constexpr size_t PARAM_CAPACITY = 80;

  static constexpr const char *FUNCTION_OPTIONS[] = {
      "Twin Wave",
      "Rings",
      "Spiral",
      "Grid",
      "Noise Contour (Projected)",
      "Primitive Lattice",
      "Noise Contour (Sphere)"};
  static constexpr const char *FUNCTION_EXPORT_OPTIONS[] = {
      "Function::TWIN_WAVE",
      "Function::RINGS",
      "Function::SPIRAL",
      "Function::GRID",
      "Function::NOISE_CONTOUR",
      "Function::PRIMITIVE_LATTICE",
      "Function::NOISE_CONTOUR_SPHERE"};
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
      "SurfaceLens::KALEIDOSCOPE_TETRAHEDRAL",
      "SurfaceLens::KALEIDOSCOPE_OCTAHEDRAL",
      "SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL",
      "SurfaceLens::KALEIDOSCOPE_TRIANGULAR_PRISM",
      "SurfaceLens::KALEIDOSCOPE_SQUARE_PRISM",
      "SurfaceLens::KALEIDOSCOPE_PENTAGONAL_PRISM",
      "SurfaceLens::KALEIDOSCOPE_HEXAGONAL_PRISM",
      "SurfaceLens::KALEIDOSCOPE_OCTAGONAL_PRISM"};
  static constexpr int NUM_LENSES = std::size(LENS_OPTIONS);
  static constexpr const char *SURFACE_NOISE_OPTIONS[] = {"None", "Direct",
                                                          "Curl"};
  static constexpr const char *SURFACE_NOISE_EXPORT_OPTIONS[] = {
      "SurfaceNoise::NONE", "SurfaceNoise::DIRECT", "SurfaceNoise::CURL"};
  static constexpr int NUM_SURFACE_NOISE = std::size(SURFACE_NOISE_OPTIONS);
  static constexpr const char *SURFACE_NOISE_PLACEMENT_OPTIONS[] = {
      "Before Lens", "After Lens"};
  static constexpr const char *SURFACE_NOISE_PLACEMENT_EXPORT_OPTIONS[] = {
      "SurfaceNoisePlacement::BEFORE_LENS",
      "SurfaceNoisePlacement::AFTER_LENS"};
  static constexpr int NUM_SURFACE_NOISE_PLACEMENTS =
      std::size(SURFACE_NOISE_PLACEMENT_OPTIONS);
  static constexpr const char *SURFACE_CURL_INTEGRATOR_OPTIONS[] = {
      "Euler", "Midpoint", "Midpoint 2x"};
  static constexpr const char *SURFACE_CURL_INTEGRATOR_EXPORT_OPTIONS[] = {
      "SurfaceCurlIntegrator::EULER", "SurfaceCurlIntegrator::MIDPOINT",
      "SurfaceCurlIntegrator::MIDPOINT_2X"};
  static constexpr int NUM_SURFACE_CURL_INTEGRATORS =
      std::size(SURFACE_CURL_INTEGRATOR_OPTIONS);
  static constexpr const char *WARP_OPTIONS[] = {"None",
                                                 "Affine Frame",
                                                 "Wave Shear",
                                                 "Vortex",
                                                 "Projected Vector Noise",
                                                 "Projected Curl Flow",
                                                 "Mirror Tile",
                                                 "Polar Chart"};
  static constexpr const char *WARP_EXPORT_OPTIONS[] = {
      "WarpStageKind::NONE",         "WarpStageKind::AFFINE_FRAME",
      "WarpStageKind::WAVE_SHEAR",   "WarpStageKind::VORTEX",
      "WarpStageKind::VECTOR_NOISE", "WarpStageKind::CURL_FLOW",
      "WarpStageKind::MIRROR_TILE",  "WarpStageKind::POLAR_CHART"};
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
  static constexpr const char *PALETTE_OPTIONS[] = {
      "Generated Triadic", "Generated Complementary", "Generated Analogous"};
  static constexpr const char *PALETTE_EXPORT_OPTIONS[] = {
      "PaletteMode::TRIADIC", "PaletteMode::COMPLEMENTARY",
      "PaletteMode::ANALOGOUS"};
  static constexpr int NUM_PALETTES = std::size(PALETTE_OPTIONS);
  static constexpr const char *PALETTE_MAPPING_OPTIONS[] = {
      "Cup", "Bell", "Linear", "Reverse"};
  static constexpr const char *PALETTE_MAPPING_EXPORT_OPTIONS[] = {
      "PaletteMapping::CUP", "PaletteMapping::BELL", "PaletteMapping::LINEAR",
      "PaletteMapping::REVERSE"};
  static constexpr int NUM_PALETTE_MAPPINGS =
      std::size(PALETTE_MAPPING_OPTIONS);
  static constexpr const char *BRIGHTNESS_ENVELOPE_OPTIONS[] = {
      "None", "Cup", "Bell", "Ascending", "Descending"};
  static constexpr const char *BRIGHTNESS_ENVELOPE_EXPORT_OPTIONS[] = {
      "BrightnessEnvelope::NONE", "BrightnessEnvelope::CUP",
      "BrightnessEnvelope::BELL", "BrightnessEnvelope::ASCENDING",
      "BrightnessEnvelope::DESCENDING"};
  static constexpr int NUM_BRIGHTNESS_ENVELOPES =
      std::size(BRIGHTNESS_ENVELOPE_OPTIONS);
  static constexpr const char *HUE_SHIFT_OPTIONS[] = {
      "None", "Noise", "Total Warp Displacement"};
  static constexpr const char *HUE_SHIFT_EXPORT_OPTIONS[] = {
      "HueShiftMode::NONE", "HueShiftMode::NOISE",
      "HueShiftMode::WARP_DISPLACEMENT"};
  static constexpr int NUM_HUE_SHIFT_MODES = std::size(HUE_SHIFT_OPTIONS);

  static constexpr float WARP_SCALE_MIN = 1.0f / 64.0f;
  static constexpr float WARP_SCALE_MAX = 100.0f;
  static constexpr float WARP_STRENGTH_MIN = -4.0f;
  static constexpr float WARP_STRENGTH_MAX = 30.0f;
  static constexpr float VECTOR_WARP_SCALE_MAX = 4.0f;
  static constexpr float VECTOR_WARP_STRENGTH_MAX = 1.0f;
  static constexpr float CURL_WARP_SCALE_MAX = 2.0f;
  static constexpr float CURL_WARP_STRENGTH_MAX = 1.0f;
  static constexpr float CURL_VECTOR_COMPONENT_MAX = 4.0f;
  static constexpr float WARP_SPEED_MIN = -1.0f / 64.0f;
  static constexpr float WARP_SPEED_MAX = 1.0f;
  static constexpr float PATTERN_FREQ_MIN = 0.1f, PATTERN_FREQ_MAX = 20.0f;
  static constexpr float SPEED_MIN = 0.0f, SPEED_MAX = 5.0f;
  static constexpr float COMPLEXITY_MIN = 0.0f, COMPLEXITY_MAX = 3.0f;
  static constexpr float PATTERN_MIX_MIN = 0.0f, PATTERN_MIX_MAX = 1.0f;
  static constexpr float PHASE2_RATE_MIN = 0.0f, PHASE2_RATE_MAX = 2.0f;
  static constexpr float POLE_FADE_MIN = 1.0f, POLE_FADE_MAX = 20.0f;
  static constexpr float SPIN_RATE_MIN = 0.0f, SPIN_RATE_MAX = 0.05f;
  static constexpr float WANDER_MIN = 0.0f, WANDER_MAX = 1.0f;
  static constexpr float LENS_MIX_MIN = 0.0f, LENS_MIX_MAX = 1.0f;
  static constexpr float HUE_SHIFT_AMOUNT_MAX = 4.0f;
  static constexpr float HUE_NOISE_AMOUNT_MAX = 1.0f;
  static constexpr float HUE_NOISE_SCALE_MIN = 1.0f / 64.0f;
  static constexpr float HUE_NOISE_SCALE_MAX = 8.0f;
  static constexpr float HUE_NOISE_SPEED_MAX = 0.001f;
  static constexpr float PALETTE_CHROMA_MIN = 0.0f;
  static constexpr float PALETTE_CHROMA_MAX = 1.0f;
  static constexpr float MAPPING_FREQUENCY_MIN = 1.0f;
  static constexpr float MAPPING_FREQUENCY_MAX = 32.0f;
  static constexpr float MAPPING_PHASE_MIN = -1.0f;
  static constexpr float MAPPING_PHASE_MAX = 1.0f;
  static constexpr float PHASE_OSCILLATION_DEPTH_MIN = 0.0f;
  static constexpr float PHASE_OSCILLATION_DEPTH_MAX = 1.0f;
  static constexpr float PHASE_OSCILLATION_SPEED_MAX = 0.01f;
  static constexpr float BRIGHTNESS_DEPTH_MIN = 0.0f;
  static constexpr float BRIGHTNESS_DEPTH_MAX = 1.0f;
  static constexpr float VALUE_OPACITY_MIN = 0.0f;
  static constexpr float VALUE_OPACITY_MAX = 1.0f;
  static constexpr float WAVE_SPIN_MIN = 0.0f, WAVE_SPIN_MAX = 0.05f;
  static constexpr float SOURCE_NOISE_SCALE_MIN = 0.0f;
  static constexpr float SOURCE_NOISE_SCALE_MAX = 2.0f;
  static constexpr float SOURCE_NOISE_RATE_MIN = -1.0f / 1024.0f;
  static constexpr float SOURCE_NOISE_RATE_MAX = 1.0f / 1024.0f;
  static constexpr float LENS_NOISE_SCALE_MIN = 1.0f / 64.0f;
  static constexpr float LENS_NOISE_SCALE_MAX = 8.0f;
  static constexpr float NOISE_RATE_MIN = -1.0f / 64.0f;
  static constexpr float NOISE_RATE_MAX = 1.0f / 64.0f;
  static constexpr float NOISE_SPEED_MIN = NOISE_RATE_MIN;
  static constexpr float NOISE_SPEED_MAX = NOISE_RATE_MAX;
  static constexpr float CELL_MIN = 1.0f / 64.0f;
  static constexpr float CELL_MAX = 8.0f;
  static constexpr float SOFTNESS_MIN = 1.0f / 1024.0f;

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
           p.surface_noise.scale >= LENS_NOISE_SCALE_MIN &&
           p.surface_noise.scale <= LENS_NOISE_SCALE_MAX &&
           p.surface_noise.strength >= -0.5f &&
           p.surface_noise.strength <= 0.5f &&
           p.surface_noise.rate >= NOISE_RATE_MIN &&
           p.surface_noise.rate <= NOISE_RATE_MAX &&
           p.surface_noise.direction >= 0.0f &&
           p.surface_noise.direction <= 1.0f &&
           enum_at_most(p.surface_noise.basis, NoiseBasis::RIDGED3) &&
           enum_at_most(p.surface_noise.integrator,
                        SurfaceCurlIntegrator::MIDPOINT_2X) &&
           p.value.iso_level >= 0.0f && p.value.iso_level <= 1.0f &&
           p.value.iso_width >= SOFTNESS_MIN && p.value.iso_width <= 0.5f &&
           p.value.band_count >= 1 && p.value.band_count <= BAND_COUNT_MAX &&
           p.value.band_phase >= 0.0f && p.value.band_phase <= TWO_PI_F &&
           p.value.cutout_threshold >= 0.0f &&
           p.value.cutout_threshold <= 1.0f &&
           p.value.cutout_softness >= SOFTNESS_MIN &&
           p.value.cutout_softness <= 0.5f && p.value.edge_width >= 0.0f &&
           p.value.edge_width <= 0.5f &&
           p.color.hue_shift_amount >= -HUE_SHIFT_AMOUNT_MAX &&
           p.color.hue_shift_amount <= HUE_SHIFT_AMOUNT_MAX &&
           p.color.hue_noise_scale >= HUE_NOISE_SCALE_MIN &&
           p.color.hue_noise_scale <= HUE_NOISE_SCALE_MAX &&
           p.color.hue_noise_speed >= -HUE_NOISE_SPEED_MAX &&
           p.color.hue_noise_speed <= HUE_NOISE_SPEED_MAX &&
           p.color.palette_chroma >= PALETTE_CHROMA_MIN &&
           p.color.palette_chroma <= PALETTE_CHROMA_MAX &&
           p.color.mapping_frequency >= MAPPING_FREQUENCY_MIN &&
           p.color.mapping_frequency <= MAPPING_FREQUENCY_MAX &&
           p.color.mapping_phase >= MAPPING_PHASE_MIN &&
           p.color.mapping_phase <= MAPPING_PHASE_MAX &&
           p.color.phase_oscillation_depth >= PHASE_OSCILLATION_DEPTH_MIN &&
           p.color.phase_oscillation_depth <= PHASE_OSCILLATION_DEPTH_MAX &&
           p.color.phase_oscillation_speed >= -PHASE_OSCILLATION_SPEED_MAX &&
           p.color.phase_oscillation_speed <= PHASE_OSCILLATION_SPEED_MAX &&
           p.color.brightness_depth >= BRIGHTNESS_DEPTH_MIN &&
           p.color.brightness_depth <= BRIGHTNESS_DEPTH_MAX &&
           p.color.value_opacity_low >= VALUE_OPACITY_MIN &&
           p.color.value_opacity_low <= VALUE_OPACITY_MAX &&
           p.color.value_opacity_high >= VALUE_OPACITY_MIN &&
           p.color.value_opacity_high <= VALUE_OPACITY_MAX;
  }

  HS_COLD_MEMBER static constexpr bool
  warp_stage_params_in_ranges(const WarpStageParams &params) {
    static_assert(sizeof(WarpStageParams) == 100,
                  "WarpStageParams field set changed - update the range check");
    return params.scale >= WARP_SCALE_MIN && params.scale <= WARP_SCALE_MAX &&
           params.strength >= WARP_STRENGTH_MIN &&
           params.strength <= WARP_STRENGTH_MAX &&
           params.speed >= WARP_SPEED_MIN && params.speed <= WARP_SPEED_MAX &&
           params.translation_x >= -4.0f && params.translation_x <= 4.0f &&
           params.translation_y >= -4.0f && params.translation_y <= 4.0f &&
           params.rotation >= -PI_F / 8.0f && params.rotation <= TWO_PI_F &&
           params.scale_x >= 0.25f && params.scale_x <= 4.0f &&
           params.scale_y >= 0.25f && params.scale_y <= 4.0f &&
           params.shear >= -0.75f && params.shear <= 0.75f &&
           params.frequency >= 0.0f && params.frequency <= 64.0f &&
           params.field_angle >= 0.0f && params.field_angle <= TWO_PI_F &&
           params.center_x >= -4.0f && params.center_x <= 4.0f &&
           params.center_y >= -4.0f && params.center_y <= 4.0f &&
           params.radius >= 1.0f / 64.0f && params.radius <= 8.0f &&
           params.turns >= -4.0f && params.turns <= 4.0f &&
           params.center_orbit_radius >= 0.0f &&
           params.center_orbit_radius <= 4.0f && params.vector_angle >= 0.0f &&
           params.vector_angle <= TWO_PI_F && params.cell_x >= CELL_MIN &&
           params.cell_x <= CELL_MAX && params.cell_y >= CELL_MIN &&
           params.cell_y <= CELL_MAX && params.offset_x >= -8.0f &&
           params.offset_x <= 8.0f && params.offset_y >= -8.0f &&
           params.offset_y <= 8.0f && params.radial_scale >= 1.0f / 64.0f &&
           params.radial_scale <= 16.0f && params.radial_phase >= 0.0f &&
           params.radial_phase <= TWO_PI_F && params.angular_phase >= 0.0f &&
           params.angular_phase <= TWO_PI_F &&
           params.edge_width >= SOFTNESS_MIN && params.edge_width <= 0.5f;
  }

  static constexpr Slots GENERATED_SURFACE_NOISE_SLOTS{
      Function::GRID,
      Projection::STEREOGRAPHIC,
      ProjectionFramePolicy::SPIN_WANDER,
      SurfaceLens::GLITCH,
      {{WarpStageKind::NONE}, {WarpStageKind::NONE}},
      SignalWeight::PROJECTION,
      ValueTransfer::LINEAR,
      CoveragePolicy::OPAQUE,
      PaletteMode::TRIADIC,
      PeirceLayout::SQUARE,
      AiroceanLayout::VERTICAL,
      BonneHemisphere::NORTH,
      GnomonicHemispherePolicy::FOLDED,
      SurfaceNoise::DIRECT,
      SurfaceNoisePlacement::AFTER_LENS};
  static constexpr Slots KALEIDOSCOPE_GENERATED_SURFACE_NOISE_SLOTS = [] {
    Slots slots = GENERATED_SURFACE_NOISE_SLOTS;
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
                  ColorParams color, OuterCameraParams outer_camera) {
    const WarpStageParams inner_warp{0.1f, 0.0f, 0.0f};
    projection.wander = outer_camera.wander;
    color.hue_noise_speed = hs::clamp(
        color.hue_noise_speed, -HUE_NOISE_SPEED_MAX, HUE_NOISE_SPEED_MAX);
    return {source,       {outer_warp, inner_warp},
            projection,   surface_lens,
            {},           color,
            outer_camera, {}};
  }

  static constexpr void normalize_config_ranges(Config &config) {
    config.params.color.hue_noise_speed =
        hs::clamp(config.params.color.hue_noise_speed, -HUE_NOISE_SPEED_MAX,
                  HUE_NOISE_SPEED_MAX);
    auto snap_affine = [](const WarpStageSpec &spec, WarpStageParams &params) {
      if (spec.kind != WarpStageKind::AFFINE_FRAME)
        return;
      params.translation_x = snap_affine_winding(params.translation_x);
      params.translation_y = snap_affine_winding(params.translation_y);
    };
    snap_affine(config.slots.warp_program.outer, config.params.warp.outer);
    snap_affine(config.slots.warp_program.inner, config.params.warp.inner);
  }

  static constexpr float snap_affine_winding(float value) {
    if (value <= -4.0f)
      return -4.0f;
    if (value >= 4.0f)
      return 4.0f;
    if (value != value)
      return value;
    return value < 0.0f ? static_cast<float>(static_cast<int>(value - 0.5f))
                        : static_cast<float>(static_cast<int>(value + 0.5f));
  }

  static constexpr void
  author_surface_noise(Params &params, const SourceParams &source,
                       const WarpStageParams &legacy_warp) {
    const float scale = legacy_warp.scale * 0.01f;
    const float strength = legacy_warp.strength / 60.0f;
    const float rate = source.speed * legacy_warp.speed / 65536.0f;
    params.surface_noise.scale =
        scale < LENS_NOISE_SCALE_MIN
            ? LENS_NOISE_SCALE_MIN
            : (scale > LENS_NOISE_SCALE_MAX ? LENS_NOISE_SCALE_MAX : scale);
    params.surface_noise.strength =
        strength < 0.0f ? 0.0f : (strength > 0.5f ? 0.5f : strength);
    params.surface_noise.rate =
        rate < NOISE_RATE_MIN ? NOISE_RATE_MIN
                              : (rate > NOISE_RATE_MAX ? NOISE_RATE_MAX : rate);
  }

  static constexpr Params authored_surface_noise_params(
      SourceParams source, WarpStageParams legacy_warp,
      ProjectionParams projection, SurfaceLensParams surface_lens,
      ColorParams color, OuterCameraParams outer_camera) {
    Params params = authored_params(source, legacy_warp, projection,
                                    surface_lens, color, outer_camera);
    author_surface_noise(params, source, legacy_warp);
    return params;
  }

  static constexpr Config wave_shear_generated_preset(
      float pattern_freq = 4.439f, float complexity = 0.5f,
      float warp_strength = 0.5f, float warp_speed = 0.015625f,
      float surface_mix = 1.0f) {
    Slots slots = GENERATED_SURFACE_NOISE_SLOTS;
    slots.warp_program.outer.kind = WarpStageKind::WAVE_SHEAR;
    slots.surface_noise = SurfaceNoise::NONE;
    slots.coverage = CoveragePolicy::PROJECTION_WEIGHT_SQUARED;
    slots.palette_mapping = PaletteMapping::CUP;
    Params params =
        authored_params({pattern_freq, 0.245f, complexity, 0.0f, 0.0f, 0.0f},
                        {1.0f, warp_strength, warp_speed}, {1.0f, 0.0f, 0.0f},
                        {surface_mix}, {0.292f, 0.6304219f, 0.0f}, {0.0f});
    params.color.palette_chroma = 0.788f;
    params.color.mapping_phase = -0.0f;
    return {slots, params};
  }

  static constexpr Config kaleidoscope_mirror_preset() {
    const Slots slots{Function::TWIN_WAVE,
                      Projection::STEREOGRAPHIC,
                      ProjectionFramePolicy::SPIN_WANDER,
                      SurfaceLens::KALEIDOSCOPE,
                      {{WarpStageKind::NONE}, {WarpStageKind::MIRROR_TILE}},
                      SignalWeight::PROJECTION,
                      ValueTransfer::LINEAR,
                      CoveragePolicy::PROJECTION_WEIGHT_SQUARED,
                      PaletteMode::TRIADIC};
    const Params params = authored_params(
        {10.158f, 0.245f, 0.513f, 0.0f, 0.8f, 0.027f}, {0.1f, 0.0f, 0.5f},
        {4.971f, 0.0f, 1.0f}, {1.0f}, {0.0f, 1.0f, 0.0f}, {1.0f});
    return {slots, params};
  }

  static constexpr Config gnomonic_grid_mirror_preset(SurfaceLens lens) {
    Slots slots{Function::GRID,
                Projection::GNOMONIC,
                ProjectionFramePolicy::IDENTITY,
                lens,
                {{WarpStageKind::MIRROR_TILE}, {WarpStageKind::NONE}},
                SignalWeight::PROJECTION,
                ValueTransfer::LINEAR,
                CoveragePolicy::EDGE_FADE,
                PaletteMode::TRIADIC};
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

  static constexpr Config gnomonic_kaleidoscope_grid_mirror_preset() {
    Config config = gnomonic_grid_mirror_preset(SurfaceLens::KALEIDOSCOPE);
    config.params.color.palette_chroma = 0.4f;
    config.params.color.hue_shift_amount = 0.424f;
    config.params.color.hue_noise_scale = 2.2033439f;
    return config;
  }

  static constexpr Config peirce_dodecahedral_generated_preset() {
    Slots slots{Function::GRID,
                Projection::PEIRCE_QUINCUNCIAL,
                ProjectionFramePolicy::SPIN_WANDER,
                SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL,
                {{WarpStageKind::NONE}, {WarpStageKind::NONE}},
                SignalWeight::PROJECTION,
                ValueTransfer::LINEAR,
                CoveragePolicy::EDGE_FADE,
                PaletteMode::TRIADIC};
    slots.peirce_layout = PeirceLayout::SQUARE;
    Params params =
        authored_params({5.0f, 0.1f, 0.5f, 0.0f, 0.8f, 0.0f}, {}, {1.0f, 0.0f},
                        {1.0f}, {0.319f, 1.0f, 0.05f / TWO_PI_F}, {1.0f});
    params.projection.central_meridian = 0.0f;
    params.projection.coordinate_scale = 1.0f;
    params.value.edge_width = 0.1f;
    return {slots, params};
  }

  static constexpr Config gnomonic_wave_shear_grid_preset() {
    Slots slots{Function::GRID,
                Projection::GNOMONIC,
                ProjectionFramePolicy::IDENTITY,
                SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL,
                {{WarpStageKind::WAVE_SHEAR}, {WarpStageKind::MIRROR_TILE}},
                SignalWeight::PROJECTION,
                ValueTransfer::LINEAR,
                CoveragePolicy::PROJECTION_WEIGHT_SQUARED,
                PaletteMode::TRIADIC};
    slots.gnomonic_hemisphere = GnomonicHemispherePolicy::FOLDED;
    WarpStageParams outer_warp;
    outer_warp.strength = -0.176f;
    outer_warp.speed = -0.00325f;
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
        {2.311f, 0.0f}, {1.0f}, {0.721f, 1.0f, 0.0f}, {1.0f});
    params.warp.inner = inner_warp;
    return {slots, params};
  }

  static constexpr Config gnomonic_affine_lattice_contour_preset() {
    const Slots slots{Function::PRIMITIVE_LATTICE,
                      Projection::GNOMONIC,
                      ProjectionFramePolicy::SPIN_WANDER,
                      SurfaceLens::NONE,
                      {{WarpStageKind::AFFINE_FRAME}, {WarpStageKind::NONE}},
                      SignalWeight::PROJECTION,
                      ValueTransfer::ISO_CONTOUR,
                      CoveragePolicy::PROJECTION_WEIGHT,
                      PaletteMode::TRIADIC};
    Params params;
    params.source.pattern_freq = 8.0f;
    params.source.speed = 0.075f;
    params.source.complexity = 0.009122372f;
    params.source.pattern_mix = 1.0f;
    params.source.secondary_rate = 1.146f;
    params.source.lattice_cell_scale = 1.22925f;
    params.source.lattice_shape_blend = 1.0f;
    params.source.lattice_softness = 0.1608203f;
    params.source.lattice_radius = 0.332981884f;
    params.warp.outer.scale = 50.7493f;
    params.warp.outer.strength = 30.0f;
    params.warp.outer.speed = 0.015625f;
    params.warp.outer.translation_x = 4.0f;
    params.warp.outer.translation_y = 4.0f;
    params.warp.outer.shear = -0.0f;
    params.warp.inner.scale = 0.1f;
    params.projection.spin_rate = 0.0208791979f;
    params.projection.wander = 0.00309175253f;
    params.surface_lens.mix = 1.0f;
    params.value.iso_level = 0.138f;
    params.value.iso_width = 0.227034181f;
    params.value.band_count = 19;
    params.value.band_phase = 6.10725641f;
    params.value.cutout_threshold = 0.5f;
    params.value.cutout_softness = 0.05f;
    params.value.edge_width = 0.327f;
    params.color.hue_shift_amount = 0.398f;
    params.color.hue_noise_scale = 0.8300313f;
    params.color.hue_noise_speed = 0.000212000014f;
    params.outer_camera.wander = 1.0f;
    params.surface_noise.scale = 0.507492959f;
    params.surface_noise.strength = 0.5f;
    params.surface_noise.rate = 5.377579e-7f;
    return {slots, params};
  }

  static constexpr Config sinusoidal_lattice_curl_preset(float noise_scale) {
    Slots slots{Function::PRIMITIVE_LATTICE,
                Projection::SINUSOIDAL,
                ProjectionFramePolicy::SPIN_WANDER,
                SurfaceLens::NONE,
                {{WarpStageKind::NONE}, {WarpStageKind::NONE}},
                SignalWeight::PROJECTION,
                ValueTransfer::LINEAR,
                CoveragePolicy::PROJECTION_WEIGHT,
                PaletteMode::TRIADIC};
    slots.palette_mapping = PaletteMapping::CUP;
    slots.brightness_envelope = BrightnessEnvelope::CUP;
    slots.surface_noise = SurfaceNoise::CURL;
    slots.surface_noise_placement = SurfaceNoisePlacement::BEFORE_LENS;
    Params params;
    params.source.pattern_freq = 3.52279997f;
    params.source.speed = 0.1f;
    params.source.complexity = 0.9f;
    params.source.pattern_mix = 1.0f;
    params.source.secondary_rate = 0.8f;
    params.source.lattice_cell_scale = 0.710265636f;
    params.source.lattice_shape_blend = 1.0f;
    params.source.lattice_softness = 0.455532223f;
    params.source.lattice_radius = 0.290762514f;
    params.warp.outer.strength = 1.0f;
    params.warp.outer.speed = 0.000343749998f;
    params.projection.pole_fade = 20.0f;
    params.projection.wander = 1.0f;
    params.surface_lens.mix = 1.0f;
    params.color.hue_shift_amount = 0.268000007f;
    params.color.hue_noise_scale = 2.0f;
    params.color.palette_chroma = 1.0f;
    params.color.mapping_phase = -0.165999994f;
    params.outer_camera.wander = 1.0f;
    params.surface_noise.scale = noise_scale;
    params.surface_noise.strength = 0.0759999976f;
    return {slots, params};
  }

  static constexpr Config stereographic_prism_polar_wave_lattice_preset() {
    Slots slots{Function::PRIMITIVE_LATTICE,
                Projection::STEREOGRAPHIC,
                ProjectionFramePolicy::SPIN_WANDER,
                SurfaceLens::KALEIDOSCOPE_TRIANGULAR_PRISM,
                {{WarpStageKind::POLAR_CHART}, {WarpStageKind::WAVE_SHEAR}},
                SignalWeight::PROJECTION,
                ValueTransfer::LINEAR,
                CoveragePolicy::PROJECTION_WEIGHT_SQUARED,
                PaletteMode::ANALOGOUS};
    slots.palette_mapping = PaletteMapping::CUP;
    slots.surface_noise_placement = SurfaceNoisePlacement::BEFORE_LENS;
    Params params;
    params.source.pattern_freq = 3.52279997f;
    params.source.speed = 0.1f;
    params.source.complexity = 0.9f;
    params.source.pattern_mix = 1.0f;
    params.source.secondary_rate = 0.8f;
    params.source.lattice_cell_scale = 0.774140596f;
    params.source.lattice_shape_blend = 1.0f;
    params.source.lattice_softness = 0.377608389f;
    params.source.lattice_radius = 0.290762514f;
    params.warp.outer.strength = 1.0f;
    params.warp.outer.speed = 0.000343749998f;
    params.warp.outer.translation_x = 4.0f;
    params.warp.inner.speed = 0.000999999931f;
    params.warp.inner.translation_x = -0.0f;
    params.warp.inner.shear = 0.75f;
    params.projection.pole_fade = 2.27300000f;
    params.projection.wander = 1.0f;
    params.surface_lens.mix = 1.0f;
    params.color.hue_shift_amount = 0.268000007f;
    params.color.hue_noise_scale = 2.0f;
    params.color.palette_chroma = 1.0f;
    params.color.mapping_phase = -0.165999994f;
    params.outer_camera.wander = 1.0f;
    params.surface_noise.scale = 3.73634386f;
    params.surface_noise.strength = 0.0759999976f;
    return {slots, params};
  }

  static constexpr Config gnomonic_dodecahedral_vector_mirror_grid_preset() {
    Slots slots{Function::GRID,
                Projection::GNOMONIC,
                ProjectionFramePolicy::IDENTITY,
                SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL,
                {{WarpStageKind::VECTOR_NOISE}, {WarpStageKind::MIRROR_TILE}},
                SignalWeight::PROJECTION,
                ValueTransfer::LINEAR,
                CoveragePolicy::PROJECTION_WEIGHT_SQUARED,
                PaletteMode::TRIADIC};
    slots.palette_mapping = PaletteMapping::CUP;
    slots.brightness_envelope = BrightnessEnvelope::CUP;
    Params params;
    params.source.pattern_freq = 4.9755f;
    params.source.speed = 0.04f;
    params.source.complexity = 1.704f;
    params.source.secondary_rate = 0.8f;
    params.source.angle_rate = 0.027f;
    params.warp.outer.strength = 0.138f;
    params.warp.outer.speed = -0.00005f;
    params.warp.outer.frequency = 1.408f;
    params.warp.outer.field_angle = 2.23053074f;
    params.warp.inner.speed = 0.00327999983f;
    params.projection.pole_fade = 2.311f;
    params.projection.wander = 1.0f;
    params.surface_lens.mix = 1.0f;
    params.color.hue_shift_amount = 0.721f;
    params.color.palette_chroma = 1.0f;
    params.color.brightness_depth = 0.655f;
    params.outer_camera.wander = 1.0f;
    return {slots, params};
  }

  static constexpr Config
  stereographic_dodecahedral_grid_inner_mirror_preset() {
    Slots slots{Function::GRID,
                Projection::STEREOGRAPHIC,
                ProjectionFramePolicy::SPIN_WANDER,
                SurfaceLens::KALEIDOSCOPE_DODECAHEDRAL,
                {{WarpStageKind::NONE}, {WarpStageKind::MIRROR_TILE}},
                SignalWeight::PROJECTION,
                ValueTransfer::LINEAR,
                CoveragePolicy::PROJECTION_WEIGHT_SQUARED,
                PaletteMode::ANALOGOUS};
    slots.palette_mapping = PaletteMapping::CUP;
    Params params;
    params.source.pattern_freq = 2.82629991f;
    params.source.complexity = 0.513f;
    params.source.secondary_rate = 0.8f;
    params.source.angle_rate = 0.0269999988f;
    params.warp.outer.scale = 0.1f;
    params.warp.outer.speed = 0.5f;
    params.warp.inner.scale = 0.1f;
    params.warp.inner.speed = 0.00013f;
    params.warp.inner.cell_y = 0.997703135f;
    params.projection.pole_fade = 3.432f;
    params.projection.wander = 1.0f;
    params.surface_lens.mix = 1.0f;
    params.color.hue_shift_amount = 0.366f;
    params.color.hue_noise_scale = 1.47215629f;
    params.color.palette_chroma = 1.0f;
    params.outer_camera.wander = 1.0f;
    return {slots, params};
  }

  static constexpr Config
  stereographic_dodecahedral_complex_grid_inner_mirror_preset() {
    Config config = stereographic_dodecahedral_grid_inner_mirror_preset();
    config.params.source.complexity = 3.0f;
    config.params.source.pattern_mix = 1.0f;
    return config;
  }

  static constexpr Config
  stereographic_dodecahedral_double_mapping_grid_inner_mirror_preset() {
    Config config =
        stereographic_dodecahedral_complex_grid_inner_mirror_preset();
    config.params.source.pattern_freq = 3.9407f;
    config.params.projection.wander = 0.165f;
    config.params.color.mapping_frequency = 2.0f;
    return config;
  }

  static constexpr Config
  equirectangular_dodecahedral_double_mapping_grid_inner_mirror_preset() {
    Config config =
        stereographic_dodecahedral_double_mapping_grid_inner_mirror_preset();
    config.slots.projection = Projection::EQUIRECTANGULAR;
    config.params.projection.pole_fade = 2.14f;
    return config;
  }

  static constexpr Config
  equirectangular_dodecahedral_grid_inner_mirror_preset() {
    Config config =
        equirectangular_dodecahedral_double_mapping_grid_inner_mirror_preset();
    config.params.color.mapping_frequency = 1.0f;
    return config;
  }

  static constexpr Config
  equirectangular_dodecahedral_fine_grid_inner_mirror_preset() {
    Config config = equirectangular_dodecahedral_grid_inner_mirror_preset();
    config.params.source.pattern_freq = 0.3985f;
    config.params.warp.inner.speed = 0.00058f;
    config.params.warp.inner.cell_y = 0.901890635f;
    config.params.color.mapping_frequency = 21.212f;
    return config;
  }

  static constexpr Config
  stereographic_hexagonal_prism_twin_wave_mirror_preset() {
    Slots slots{Function::TWIN_WAVE,
                Projection::STEREOGRAPHIC,
                ProjectionFramePolicy::SPIN_WANDER,
                SurfaceLens::KALEIDOSCOPE_HEXAGONAL_PRISM,
                {{WarpStageKind::NONE}, {WarpStageKind::MIRROR_TILE}},
                SignalWeight::PROJECTION,
                ValueTransfer::LINEAR,
                CoveragePolicy::PROJECTION_WEIGHT_SQUARED,
                PaletteMode::ANALOGOUS};
    slots.palette_mapping = PaletteMapping::BELL;
    Params params = authored_params(
        {3.881f, 0.128598228f, 0.513f, 0.0f, 0.8f, 0.027f}, {0.1f, 0.0f, 0.5f},
        {4.971f, 0.0f, 1.0f}, {1.0f}, {0.226f, 1.47215629f, 0.000138f}, {1.0f});
    params.color.palette_chroma = 1.0f;
    params.color.mapping_frequency = 1.341f;
    params.color.mapping_phase = -1.0f;
    return {slots, params};
  }

  static constexpr Config stereographic_glitch_grid_mirror_preset() {
    Slots slots{Function::GRID,
                Projection::STEREOGRAPHIC,
                ProjectionFramePolicy::IDENTITY,
                SurfaceLens::GLITCH,
                {{WarpStageKind::MIRROR_TILE}, {WarpStageKind::NONE}},
                SignalWeight::PROJECTION,
                ValueTransfer::LINEAR,
                CoveragePolicy::EDGE_FADE,
                PaletteMode::TRIADIC};
    slots.hue_shift = HueShiftMode::WARP_DISPLACEMENT;
    WarpStageParams outer_warp;
    outer_warp.rotation = 0.295309722f;
    outer_warp.cell_x = 5.381125f;
    outer_warp.offset_x = 1.344f;
    outer_warp.offset_y = -1.456f;
    Params params =
        authored_params({2.5477f, 0.235f, 1.854f, 0.0f, 1.0f, 0.0f}, outer_warp,
                        {1.4f, 0.0f}, {1.0f}, {2.048f, 1.0f, 0.0f}, {1.0f});
    params.value.edge_width = 0.5f;
    params.color.palette_chroma = 0.292f;
    return {slots, params};
  }

  static constexpr Config stereographic_mobius_twin_wave_inner_mirror_preset() {
    Slots slots{Function::TWIN_WAVE,
                Projection::STEREOGRAPHIC,
                ProjectionFramePolicy::SPIN_WANDER,
                SurfaceLens::MOBIUS,
                {{WarpStageKind::NONE}, {WarpStageKind::MIRROR_TILE}},
                SignalWeight::PROJECTION,
                ValueTransfer::LINEAR,
                CoveragePolicy::PROJECTION_WEIGHT_SQUARED,
                PaletteMode::COMPLEMENTARY};
    slots.brightness_envelope = BrightnessEnvelope::CUP;
    slots.hue_shift = HueShiftMode::WARP_DISPLACEMENT;
    Params params = authored_params(
        {10.158f, 0.245f, 0.513f, 0.0f, 0.8f, 0.027f}, {0.1f, 0.0f, 0.5f},
        {2.102f, 0.0f, 1.0f}, {1.0f}, {0.312f, 1.0f, 0.0f}, {1.0f});
    params.surface_lens.mobius = {-1.072f, 0.304f, 0.416f,      0.0f,
                                  0.0f,    0.0f,   0.70710677f, 0.0f};
    params.color.palette_chroma = 0.398f;
    return {slots, params};
  }

  static constexpr Config stereographic_mobius_animated_inner_mirror_preset() {
    Config config = stereographic_mobius_twin_wave_inner_mirror_preset();
    config.params.warp.inner.speed = 0.005875f;
    config.params.warp.inner.cell_x = 0.2791094f;
    config.params.warp.inner.cell_y = 6.810328f;
    return config;
  }

  static constexpr std::array<Preset, 24> PRESETS = {{
      {wave_shear_generated_preset(),
       InversePipelineId::GLITCH_NOISE_GRID_WAVE_SHEAR},
      {kaleidoscope_mirror_preset(),
       InversePipelineId::KALEIDOSCOPE_TWIN_WAVE_INNER_MIRROR},
      {gnomonic_kaleidoscope_grid_mirror_preset(),
       InversePipelineId::GNOMONIC_KALEIDOSCOPE_GRID_MIRROR},
      {gnomonic_grid_mirror_preset(SurfaceLens::GLITCH),
       InversePipelineId::GNOMONIC_GLITCH_GRID_MIRROR},
      {peirce_dodecahedral_generated_preset(),
       InversePipelineId::PEIRCE_DODECAHEDRAL_GRID},
      {gnomonic_wave_shear_grid_preset(),
       InversePipelineId::GNOMONIC_DODECAHEDRAL_GRID_WAVE_MIRROR},
      {gnomonic_affine_lattice_contour_preset(),
       InversePipelineId::GNOMONIC_AFFINE_LATTICE_CONTOUR},
      {sinusoidal_lattice_curl_preset(1.78815627f),
       InversePipelineId::SINUSOIDAL_CURL_LATTICE},
      {sinusoidal_lattice_curl_preset(3.29720306f),
       InversePipelineId::SINUSOIDAL_CURL_LATTICE},
      {stereographic_prism_polar_wave_lattice_preset(),
       InversePipelineId::STEREOGRAPHIC_PRISM_POLAR_WAVE_LATTICE},
      {gnomonic_dodecahedral_vector_mirror_grid_preset(),
       InversePipelineId::GNOMONIC_DODECAHEDRAL_GRID_VECTOR_MIRROR},
      {stereographic_dodecahedral_grid_inner_mirror_preset(),
       InversePipelineId::STEREOGRAPHIC_DODECAHEDRAL_GRID_INNER_MIRROR},
      {stereographic_hexagonal_prism_twin_wave_mirror_preset(),
       InversePipelineId::STEREOGRAPHIC_HEXAGONAL_PRISM_TWIN_WAVE_INNER_MIRROR},
      {stereographic_dodecahedral_complex_grid_inner_mirror_preset(),
       InversePipelineId::STEREOGRAPHIC_DODECAHEDRAL_GRID_INNER_MIRROR},
      {stereographic_dodecahedral_double_mapping_grid_inner_mirror_preset(),
       InversePipelineId::STEREOGRAPHIC_DODECAHEDRAL_GRID_INNER_MIRROR},
      {equirectangular_dodecahedral_double_mapping_grid_inner_mirror_preset(),
       InversePipelineId::EQUIRECTANGULAR_DODECAHEDRAL_GRID_INNER_MIRROR},
      {equirectangular_dodecahedral_grid_inner_mirror_preset(),
       InversePipelineId::EQUIRECTANGULAR_DODECAHEDRAL_GRID_INNER_MIRROR},
      {equirectangular_dodecahedral_fine_grid_inner_mirror_preset(),
       InversePipelineId::EQUIRECTANGULAR_DODECAHEDRAL_GRID_INNER_MIRROR},
      {stereographic_glitch_grid_mirror_preset(),
       InversePipelineId::STEREOGRAPHIC_GLITCH_GRID_MIRROR},
      {stereographic_mobius_twin_wave_inner_mirror_preset(),
       InversePipelineId::STEREOGRAPHIC_MOBIUS_TWIN_WAVE_INNER_MIRROR},
      {stereographic_mobius_animated_inner_mirror_preset(),
       InversePipelineId::STEREOGRAPHIC_MOBIUS_TWIN_WAVE_INNER_MIRROR},
      {wave_shear_generated_preset(3.1447f, 0.5f, 2.72f, 0.00690625f),
       InversePipelineId::GLITCH_NOISE_GRID_WAVE_SHEAR},
      {wave_shear_generated_preset(7.5227f, 1.698f, 0.0f, 0.00690625f),
       InversePipelineId::GLITCH_NOISE_GRID_WAVE_SHEAR},
      {wave_shear_generated_preset(8.8162f, 1.698f, 1.376f, 0.00559375f, 0.0f),
       InversePipelineId::GLITCH_NOISE_GRID_WAVE_SHEAR},
  }};
  static_assert(
      [] {
        for (const Preset &preset : PRESETS)
          if (!valid_config(preset.config) ||
              preset.pipeline == InversePipelineId::NONE)
            return false;
        return true;
      }(),
      "a ShaderBall preset lies outside its registered range");
  static_assert(
      [] {
        for (size_t index = 0; index < PRESETS.size(); ++index)
          if (!transition_admitted(
                  PRESETS[index].config,
                  PRESETS[(index + 1) % PRESETS.size()].config))
            return false;
        return true;
      }(),
      "a ShaderBall preset edge lacks continuous transition admission");

  static constexpr std::array<Choreo, PRESETS.size()> CHOREO = [] {
    std::array<Choreo, PRESETS.size()> choreo;
    for (Choreo &entry : choreo)
      entry = {0, 0, 480, false};
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

  uint32_t palette_hue = 0;
  uint32_t complementary_hue = 0;
  uint32_t analogous_hue = 0;
  float palette_chroma = 0.62f;
  PaletteCycler triadic_palette_cycler;
  PaletteCycler complementary_palette_cycler;
  PaletteCycler analogous_palette_cycler;

  Slots active_slots = PRESETS[0].config.slots;
  InversePipelineId active_pipeline = PRESETS[0].pipeline;
#if HS_ENABLE_PARAM_GUI_BRIDGE
  Config display_config = PRESETS[0].config;
  std::array<PendingEdit, PARAM_CAPACITY> pending_edits{};
  size_t pending_edit_count = 0;
  mutable std::array<char, 1024> warning_text{};
#endif
  RequestedConfig requested_config = PRESETS[0].config;
  Config published_config = PRESETS[0].config;
#if HS_ENABLE_PARAM_GUI_BRIDGE
  Config accepted_config = PRESETS[0].config;
#endif
  bool requested_schema_bound = false;
  bool registered_range_clamped = false;
  bool fixed_topology = false;
  std::span<const uint8_t> preset_view{};
  uint16_t preset_dwell_remaining = 0;
  bool preset_dwell_armed = false;
  Blend blend{PRESETS[0].config.params,
              palette_mapping_weights(PRESETS[0].config.slots.palette_mapping)};
  LookRuntime runtime;
#if HS_ENABLE_TEST_HOOKS
  uint32_t walk_step_count = 0;
  uint32_t generated_palette_step_count = 0;
#endif
#if defined(HS_PROFILE_ENABLE)
  bool profile_program_valid = false;
  size_t profile_program_preset = 0;
  InversePipelineId profile_program_pipeline = InversePipelineId::NONE;
  ProfileEndpoint profile_program_endpoint = ProfileEndpoint::STEADY;
#endif

  static constexpr size_t FOOTPRINT_BYTES =
      gamut_lut_bytes(GAMUT_ANGLE_STEPS, GAMUT_L_STEPS) +
      3 * PaletteCycler::generated_arena_bytes() +
      PARAM_CAPACITY * sizeof(ParamDef) + sizeof(StateBundle) +
      alignof(StateBundle);
  static_assert(
      FOOTPRINT_BYTES <= DEVICE_PERSISTENT_BUDGET,
      "ShaderBall persistent footprint exceeds the default partition");
};
